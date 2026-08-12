# Packages
library(testthat)
library(seqwrap)


# Helpers shared across the tests in this file
sw_test_container <- function(n_genes = 12) {
  dat <- seqwrap::simcounts(
    seed = 1,
    n_samples = 30,
    n_genes = n_genes,
    clusters = 10
  )

  seqwrap::seqwrap_compose(
    modelfun = stats::lm,
    arguments = list(formula = y ~ x),
    data = dat$data,
    metadata = dat$metadata,
    samplename = "sample"
  )
}

sw_test_summary <- function(x) {
  cf <- coef(summary(x))
  data.frame(
    coef = rownames(cf),
    estimate = unname(cf[, 1]),
    se = unname(cf[, 2]),
    stringsAsFactors = FALSE
  )
}

sw_test_eval <- function(x) {
  data.frame(logLik = as.numeric(logLik(x)))
}


test_that("chunk_size does not change results", {
  container <- sw_test_container()

  # 12 targets with return_models = TRUE trips the memory warning by design
  fit <- function(cs) {
    suppressWarnings(seqwrap::seqwrap(
      container,
      summary_fun = sw_test_summary,
      eval_fun = sw_test_eval,
      return_models = TRUE,
      cores = 1,
      verbose = FALSE,
      chunk_size = cs
    ))
  }

  one <- fit(1)
  five <- fit(5)
  auto <- fit(NULL)

  # Target order and content are independent of how work is chunked
  expect_identical(one@targets, five@targets)
  expect_identical(one@targets, auto@targets)
  expect_identical(one@summaries, five@summaries)
  expect_identical(one@summaries, auto@summaries)
  expect_identical(one@evaluations, auto@evaluations)
  expect_equal(
    unname(sapply(one@models, function(m) unname(coef(m)[2]))),
    unname(sapply(five@models, function(m) unname(coef(m)[2])))
  )
})


test_that("retaining many models raises a memory warning", {
  # Twelve targets, above the threshold
  expect_warning(
    seqwrap::seqwrap(
      sw_test_container(n_genes = 12),
      return_models = TRUE,
      cores = 1,
      verbose = FALSE
    ),
    "12 model objects"
  )

  # A small subset, which is what return_models is intended for, stays quiet
  expect_no_warning(
    seqwrap::seqwrap(
      sw_test_container(n_genes = 12),
      return_models = TRUE,
      subset = 1:5,
      cores = 1,
      verbose = FALSE
    )
  )

  # And the default does not retain models at all
  res <- seqwrap::seqwrap(
    sw_test_container(n_genes = 12),
    cores = 1,
    verbose = FALSE
  )
  expect_null(res@models)
})


test_that("invalid chunk_size is rejected", {
  container <- sw_test_container(n_genes = 4)

  expect_error(
    seqwrap::seqwrap(container, chunk_size = 0, cores = 1, verbose = FALSE),
    "positive integer"
  )
  expect_error(
    seqwrap::seqwrap(
      container,
      chunk_size = c(2, 3),
      cores = 1,
      verbose = FALSE
    ),
    "single integer"
  )
})


test_that("cache = 'memory' binds summaries and evaluations", {
  container <- sw_test_container()

  args <- list(
    container,
    summary_fun = sw_test_summary,
    eval_fun = sw_test_eval,
    return_models = FALSE,
    cores = 1,
    verbose = FALSE
  )

  none <- do.call(seqwrap::seqwrap, c(args, list(cache = "none")))
  memory <- do.call(seqwrap::seqwrap, c(args, list(cache = "memory")))

  # The reducing modes return a single data frame per slot
  expect_true(is.list(none@summaries) && !is.data.frame(none@summaries))
  expect_s3_class(memory@summaries, "data.frame")
  expect_s3_class(memory@evaluations, "data.frame")
  expect_null(memory@cache)

  # A target column is added, matching seqwrap_summarise() output
  expect_true("target" %in% names(memory@summaries))
  expect_identical(names(memory@summaries)[1], "target")

  # Combining gives the same answer either way
  expect_identical(
    seqwrap_summarise(none, verbose = FALSE)$summaries,
    seqwrap_summarise(memory, verbose = FALSE)$summaries
  )
  expect_identical(
    seqwrap_summarise(none, verbose = FALSE)$evaluations,
    seqwrap_summarise(memory, verbose = FALSE)$evaluations
  )

  # Conditions are recorded in the same long form regardless of cache mode
  expect_identical(names(memory@errors), names(none@errors))
  expect_identical(
    names(memory@errors),
    c("target", "stage", "type", "class", "message")
  )
})


test_that("cache = 'disk' spills chunks and reads them back", {
  container <- sw_test_container()

  args <- list(
    container,
    summary_fun = sw_test_summary,
    eval_fun = sw_test_eval,
    return_models = FALSE,
    cores = 1,
    verbose = FALSE
  )

  none <- do.call(seqwrap::seqwrap, c(args, list(cache = "none")))
  disk <- do.call(
    seqwrap::seqwrap,
    c(args, list(cache = "disk", chunk_size = 5))
  )

  # Nothing is retained in memory
  expect_null(disk@summaries)
  expect_null(disk@evaluations)
  expect_false(is.null(disk@cache))

  # The default cache location is inside the session temporary directory
  expect_true(
    startsWith(normalizePath(disk@cache$path), normalizePath(tempdir()))
  )
  expect_equal(disk@cache$n_chunks, 3)
  expect_equal(sum(disk@cache$chunks$n), disk@k)
  expect_true(all(file.exists(disk@cache$chunks$summaries)))

  # Reading the cache reproduces the in-memory result exactly
  expect_identical(
    seqwrap_summarise(none, verbose = FALSE)$summaries,
    seqwrap_summarise(disk, verbose = FALSE)$summaries
  )
  expect_identical(
    seqwrap_summarise(none, verbose = FALSE)$evaluations,
    seqwrap_summarise(disk, verbose = FALSE)$evaluations
  )

  # print() copes with the cached representation
  expect_no_error(capture.output(print(disk)))

  # The cache can be removed again
  path <- disk@cache$path
  expect_gt(seqwrap_cache_clear(disk), 0)
  expect_false(dir.exists(path))
})


test_that("a user supplied cache_path is honoured", {
  container <- sw_test_container(n_genes = 6)
  target_dir <- file.path(tempdir(), "seqwrap-user-cache")
  unlink(target_dir, recursive = TRUE)

  res <- seqwrap::seqwrap(
    container,
    summary_fun = sw_test_summary,
    eval_fun = NULL,
    return_models = FALSE,
    cores = 1,
    verbose = FALSE,
    cache = "disk",
    cache_path = target_dir
  )

  expect_equal(res@cache$path, target_dir)
  expect_true(length(list.files(target_dir, pattern = "\\.rds$")) > 0)

  seqwrap_cache_clear(res)
})


test_that("save_models writes one file per target", {
  container <- sw_test_container(n_genes = 8)
  model_dir <- file.path(tempdir(), "seqwrap-test-models")
  unlink(model_dir, recursive = TRUE)
  dir.create(model_dir, recursive = TRUE)

  res <- seqwrap::seqwrap(
    container,
    return_models = FALSE,
    save_models = TRUE,
    model_path = model_dir,
    cores = 1,
    verbose = FALSE
  )

  saved <- list.files(model_dir, pattern = "\\.rds$")

  # One file per target, named after the target rather than a single
  # repeatedly overwritten file
  expect_equal(length(saved), res@k)
  expect_setequal(
    sort(tools::file_path_sans_ext(saved)),
    sort(res@targets)
  )
  expect_s3_class(readRDS(file.path(model_dir, saved[1])), "lm")

  unlink(model_dir, recursive = TRUE)
})


test_that("save_models works together with disk caching", {
  container <- sw_test_container(n_genes = 6)
  model_dir <- file.path(tempdir(), "seqwrap-test-models-cache")
  unlink(model_dir, recursive = TRUE)
  dir.create(model_dir, recursive = TRUE)

  res <- seqwrap::seqwrap(
    container,
    summary_fun = sw_test_summary,
    eval_fun = NULL,
    return_models = FALSE,
    save_models = TRUE,
    model_path = model_dir,
    cache = "disk",
    cores = 1,
    verbose = FALSE
  )

  expect_equal(length(list.files(model_dir, pattern = "\\.rds$")), res@k)
  expect_false(is.null(res@cache))

  seqwrap_cache_clear(res)
  unlink(model_dir, recursive = TRUE)
})


test_that("target-wise data stays aligned when chunked", {
  dat <- seqwrap::simcounts(
    seed = 1,
    n_samples = 30,
    n_genes = 12,
    clusters = 10
  )

  # A distinct offset per target makes misalignment detectable
  targetdata <- data.frame(shift = seq_len(12) * 100)

  container <- seqwrap::seqwrap_compose(
    modelfun = stats::lm,
    arguments = alist(formula = y ~ x, offset = rep(shift, 30)),
    data = dat$data,
    metadata = dat$metadata,
    targetdata = targetdata,
    samplename = "sample"
  )

  fit <- function(cs) {
    res <- seqwrap::seqwrap(
      container,
      summary_fun = sw_test_summary,
      eval_fun = NULL,
      cores = 1,
      verbose = FALSE,
      chunk_size = cs
    )
    combined <- seqwrap_summarise(res, verbose = FALSE)$summaries
    combined$estimate[!duplicated(combined$target)]
  }

  expect_identical(fit(1), fit(5))
  expect_identical(fit(1), fit(NULL))
})
