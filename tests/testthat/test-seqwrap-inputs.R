# Tests for user-facing input handling. Several of these cover paths that were
# previously untested and silently wrong.

library(testthat)
library(seqwrap)


sw_input_data <- function(n_genes = 6) {
  seqwrap::simcounts(
    seed = 1,
    n_samples = 20,
    n_genes = n_genes,
    clusters = 5
  )
}

sw_input_container <- function(dat = sw_input_data(), ...) {
  seqwrap::seqwrap_compose(
    modelfun = stats::lm,
    arguments = list(formula = y ~ x),
    data = dat$data,
    metadata = dat$metadata,
    samplename = "sample",
    ...
  )
}


test_that("rownames = TRUE uses row names as target identifiers", {
  dat <- sw_input_data()

  # The same data expressed two ways: identifiers in the first column, and
  # identifiers as row names
  reference <- seqwrap::seqwrap(
    sw_input_container(dat),
    cores = 1,
    verbose = FALSE
  )

  by_rowname <- dat$data[, -1]
  rownames(by_rowname) <- dat$data[[1]]

  rn <- seqwrap::seqwrap(
    seqwrap::seqwrap_compose(
      modelfun = stats::lm,
      arguments = list(formula = y ~ x),
      data = by_rowname,
      rownames = TRUE,
      metadata = dat$metadata,
      samplename = "sample"
    ),
    cores = 1,
    verbose = FALSE
  )

  # Previously the first sample column was used as the identifier, which
  # silently produced the wrong targets and dropped a sample from the model
  expect_identical(rn@targets, reference@targets)
  expect_equal(rn@k, reference@k)
  expect_equal(rn@summaries$estimate, reference@summaries$estimate)
})


test_that("objects in `exported` are available to summary and eval functions", {
  scale_factor <- 7

  summary_uses_exported <- function(x) {
    data.frame(scaled = unname(coef(x)[2]) * scale_factor)
  }
  eval_uses_exported <- function(x) {
    data.frame(seen = scale_factor)
  }

  res <- seqwrap::seqwrap(
    sw_input_container(),
    summary_fun = summary_uses_exported,
    eval_fun = eval_uses_exported,
    exported = list(scale_factor = scale_factor),
    cores = 2,
    verbose = FALSE
  )

  # The objects have to exist on the worker; a function defined in the global
  # environment loses that environment when it is sent there
  expect_equal(nrow(seqwrap_errors(res, type = "error")), 0)
  expect_false(is.null(res@summaries))
  expect_equal(nrow(res@summaries), res@k)
  expect_true(all(res@evaluations$seen == scale_factor))
})


test_that("unnamed elements in `exported` are rejected", {
  expect_error(
    seqwrap::seqwrap(
      sw_input_container(),
      summary_fun = function(x) data.frame(a = 1),
      eval_fun = NULL,
      exported = list(1, 2),
      cores = 1,
      verbose = FALSE
    ),
    "must be named"
  )
})


test_that("cores accepts NULL and rejects malformed values", {
  container <- sw_input_container(sw_input_data(n_genes = 3))

  # NULL means sequential
  expect_no_error(
    seqwrap::seqwrap(container, cores = NULL, verbose = FALSE)
  )

  # Validation happens before any cluster is created
  expect_error(
    seqwrap::seqwrap(container, cores = "two", verbose = FALSE),
    "must be NULL, a single integer"
  )
  expect_error(
    seqwrap::seqwrap(container, cores = c(1, 2), verbose = FALSE),
    "must be NULL, a single integer"
  )
})


test_that("missing or blank target identifiers are rejected clearly", {
  dat <- sw_input_data()

  blank <- dat$data
  blank[[1]][2] <- ""
  expect_error(
    seqwrap::seqwrap(sw_input_container(list(data = blank,
                                             metadata = dat$metadata)),
                     cores = 1, verbose = FALSE),
    "missing or blank"
  )

  missing <- dat$data
  missing[[1]][3] <- NA_character_
  expect_error(
    seqwrap::seqwrap(sw_input_container(list(data = missing,
                                             metadata = dat$metadata)),
                     cores = 1, verbose = FALSE),
    "missing or blank"
  )
})


test_that("an empty subset returns an empty result rather than failing", {
  res <- seqwrap::seqwrap(
    sw_input_container(),
    subset = integer(0),
    cores = 1,
    verbose = FALSE
  )

  expect_equal(res@k, 0)
  expect_identical(res@targets, character(0))
  expect_equal(nrow(res@errors), 0)
  expect_no_error(capture.output(print(res)))
  expect_no_error(seqwrap_summarise(res, verbose = FALSE))
})


test_that("target identifiers unsafe as file names are sanitised", {
  dat <- sw_input_data()
  dat$data[[1]] <- c("a/b", "c:d", "CON", "e|f", "normal", "with space")

  model_dir <- file.path(tempdir(), "seqwrap-unsafe-names")
  unlink(model_dir, recursive = TRUE)
  dir.create(model_dir, recursive = TRUE)

  res <- seqwrap::seqwrap(
    sw_input_container(dat),
    save_models = TRUE,
    model_path = model_dir,
    cores = 1,
    verbose = FALSE
  )

  saved <- list.files(model_dir, pattern = "\\.rds$")

  # One file per target, none overwritten, and nothing with a path separator
  # or a Windows reserved device name
  expect_equal(length(saved), res@k)
  expect_false(any(grepl("[/:|]", saved)))
  expect_false(any(toupper(tools::file_path_sans_ext(saved)) == "CON"))

  unlink(model_dir, recursive = TRUE)
})


test_that("results do not depend on the number of cores", {
  container <- sw_input_container(sw_input_data(n_genes = 8))

  one <- seqwrap::seqwrap(container, cores = 1, verbose = FALSE)
  two <- seqwrap::seqwrap(container, cores = 2, verbose = FALSE)

  expect_identical(one@targets, two@targets)
  expect_equal(one@summaries, two@summaries)
  expect_equal(one@evaluations, two@evaluations)
})


test_that("a named list of data frames is accepted as data", {
  dat <- sw_input_data()

  res <- seqwrap::seqwrap(
    seqwrap::seqwrap_compose(
      modelfun = stats::lm,
      arguments = list(formula = y ~ x),
      data = list(y = dat$data, w = dat$data),
      metadata = dat$metadata,
      samplename = "sample"
    ),
    cores = 1,
    verbose = FALSE
  )

  expect_equal(res@k, nrow(dat$data))
  expect_equal(length(res@targets), nrow(dat$data))
})


test_that("seqwrap_compose updates an existing container", {
  container <- sw_input_container()

  updated <- seqwrap::seqwrap_compose(
    x = container,
    update = list(samplename = "sample", additional_vars = "cluster")
  )

  expect_true(S7::S7_inherits(updated, seqwrap::swcontainer))
  expect_identical(updated@additional_vars, "cluster")

  # An empty update returns the container unchanged
  expect_identical(
    seqwrap::seqwrap_compose(x = container, update = list())@samplename,
    container@samplename
  )
})


test_that("DGEList input is rejected with a pointer to the alternative", {
  skip_if_not_installed("edgeR")

  dat <- sw_input_data()
  counts <- as.matrix(dat$data[, -1])
  rownames(counts) <- dat$data[[1]]
  dge <- edgeR::DGEList(counts = counts)

  expect_error(
    seqwrap::seqwrap_compose(x = dge),
    "DGEList input is no longer accepted"
  )
  expect_error(
    seqwrap::seqwrap(dge, cores = 1, verbose = FALSE),
    "DGEList input is no longer accepted"
  )
})


test_that("input that is neither a container nor NULL is rejected", {
  expect_error(
    seqwrap::seqwrap(data.frame(a = 1), cores = 1, verbose = FALSE),
    "must be a swcontainer object or NULL"
  )
})


test_that("additional_vars keeps extra metadata available to the model", {
  dat <- sw_input_data()

  expect_no_error(
    seqwrap::seqwrap(
      sw_input_container(dat, additional_vars = "cluster"),
      cores = 1,
      verbose = FALSE
    )
  )
})
