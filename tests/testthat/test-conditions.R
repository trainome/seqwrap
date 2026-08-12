# Packages
library(testthat)
library(seqwrap)


sw_cond_container <- function(modelfun, n_genes = 6) {
  dat <- seqwrap::simcounts(
    seed = 1,
    n_samples = 30,
    n_genes = n_genes,
    clusters = 10
  )

  seqwrap::seqwrap_compose(
    modelfun = modelfun,
    arguments = list(formula = y ~ x),
    data = dat$data,
    metadata = dat$metadata,
    samplename = "sample"
  )
}

# Warns for roughly half the targets, and always returns a usable model. The
# decision is driven by the data so that it survives being sent to a worker.
sw_half_warn <- function(formula, data) {
  m <- stats::lm(formula, data = data)
  if (as.integer(round(sum(data$y))) %% 2 == 0) {
    warning("simulated convergence warning")
  }
  m
}

sw_always_warn <- function(formula, data) {
  m <- stats::lm(formula, data = data)
  warning("simulated convergence warning")
  m
}


test_that("a warning during fitting does not discard the model", {
  res <- suppressWarnings(seqwrap::seqwrap(
    sw_cond_container(sw_always_warn),
    cores = 1,
    verbose = FALSE,
    return_models = TRUE
  ))

  # Every target warned, and every target must still have its model
  expect_equal(nrow(seqwrap_errors(res, "fit", "warning")), res@k)
  expect_equal(sum(sapply(res@models, is.null)), 0)
  expect_s3_class(res@models[[1]], "lm")

  # The warning is still recorded
  expect_true(all(seqwrap_errors(res, "fit", "warning")$class == "simpleWarning"))
  expect_match(
    seqwrap_errors(res, "fit", "warning")$message[1],
    "simulated convergence warning"
  )
})


test_that("targets that warn still contribute summaries and evaluations", {
  res <- suppressWarnings(seqwrap::seqwrap(
    sw_cond_container(sw_half_warn),
    cores = 1,
    verbose = FALSE,
    return_models = FALSE
  ))

  n_warn <- nrow(seqwrap_errors(res, "fit", "warning"))
  expect_gt(n_warn, 0)
  expect_lt(n_warn, res@k)

  # No target should be missing a summary or an evaluation merely because the
  # fitting algorithm warned
  expect_equal(length(unique(res@summaries$target)), res@k)
  expect_equal(length(unique(res@evaluations$target)), res@k)

  combined <- seqwrap_summarise(res, verbose = FALSE)
  expect_equal(length(unique(combined$summaries$target)), res@k)
  expect_equal(length(unique(combined$evaluations$target)), res@k)
})


test_that("a warning in the summary function does not discard the summary", {
  warn_summary <- function(x) {
    warning("simulated summary warning")
    cf <- coef(summary(x))
    data.frame(coef = rownames(cf), estimate = unname(cf[, 1]))
  }

  res <- suppressWarnings(seqwrap::seqwrap(
    sw_cond_container(stats::lm),
    summary_fun = warn_summary,
    eval_fun = NULL,
    cores = 1,
    verbose = FALSE,
    return_models = FALSE
  ))

  expect_equal(nrow(seqwrap_errors(res, "summary", "warning")), res@k)
  expect_equal(length(unique(res@summaries$target)), res@k)
  expect_equal(
    nrow(seqwrap_summarise(res, verbose = FALSE)$summaries),
    res@k * 2
  )
})


test_that("an error during fitting still yields NULL results and is recorded", {
  always_fail <- function(formula, data) stop("simulated fitting failure")

  res <- seqwrap::seqwrap(
    sw_cond_container(always_fail),
    cores = 1,
    verbose = FALSE,
    return_models = TRUE
  )

  expect_equal(nrow(seqwrap_errors(res, "fit", "error")), res@k)
  expect_equal(sum(sapply(res@models, is.null)), res@k)

  # Summary and evaluation functions are not called for a failed fit, so they
  # report neither a result nor a spurious error of their own
  expect_null(res@summaries)
  expect_null(res@evaluations)
  expect_equal(nrow(seqwrap_errors(res, "summary", "error")), 0)
})


test_that("degenerate zero column results are dropped rather than counted", {
  # broom.mixed::tidy(NULL) returns a 0 x 0 tibble, which is neither NULL nor an
  # error. Adding a target column to one would invent an all-missing row.
  degenerate <- function(x) tibble::tibble()

  res <- seqwrap::seqwrap(
    sw_cond_container(stats::lm),
    summary_fun = degenerate,
    eval_fun = NULL,
    cores = 1,
    verbose = FALSE,
    return_models = FALSE
  )

  combined <- seqwrap_summarise(res, verbose = FALSE)
  expect_null(combined$summaries)
})


test_that("summaries with differing columns are bound rather than erroring", {
  # A defensive summary function that returns a different shape for some
  # targets should not bring down the whole run
  varying <- function(x) {
    cf <- coef(summary(x))
    out <- data.frame(coef = rownames(cf), estimate = unname(cf[, 1]))
    if (unname(cf[1, 1]) > 0) out$extra <- "present"
    out
  }

  for (cache in c("none", "memory", "disk")) {
    res <- seqwrap::seqwrap(
      sw_cond_container(stats::lm),
      summary_fun = varying,
      eval_fun = NULL,
      cores = 1,
      verbose = FALSE,
      return_models = FALSE,
      cache = cache
    )

    combined <- seqwrap_summarise(res, verbose = FALSE)
    expect_true("extra" %in% names(combined$summaries))
    expect_equal(length(unique(combined$summaries$target)), res@k)

    if (!is.null(res@cache)) seqwrap_cache_clear(res)
  }
})


test_that("warnings are handled the same way in parallel and when chunked", {
  reference <- suppressWarnings(seqwrap::seqwrap(
    sw_cond_container(sw_half_warn),
    cores = 1,
    verbose = FALSE,
    return_models = FALSE,
    chunk_size = 1
  ))

  for (cores in c(1, 2)) {
    for (cs in list(2, NULL)) {
      res <- suppressWarnings(seqwrap::seqwrap(
        sw_cond_container(sw_half_warn),
        cores = cores,
        verbose = FALSE,
        return_models = FALSE,
        chunk_size = cs
      ))
      expect_identical(res@summaries, reference@summaries)
      expect_equal(
        nrow(seqwrap_errors(res, "fit", "warning")),
        nrow(seqwrap_errors(reference, "fit", "warning"))
      )
    }
  }
})
