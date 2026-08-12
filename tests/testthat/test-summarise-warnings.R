# Packages
library(testthat)
library(seqwrap)


# Warns for roughly half the targets while always returning a usable model.
# The decision is driven by the data so it survives being sent to a worker.
sw_warn_results <- function(cache = "none") {
  dat <- seqwrap::simcounts(
    seed = 1,
    n_samples = 30,
    n_genes = 6,
    clusters = 10
  )

  half_warn <- function(formula, data) {
    m <- stats::lm(formula, data = data)
    if (as.integer(round(sum(data$y))) %% 2 == 0) {
      warning("simulated singular fit")
    }
    m
  }

  suppressWarnings(seqwrap::seqwrap(
    seqwrap::seqwrap_compose(
      modelfun = half_warn,
      arguments = list(formula = y ~ x),
      data = dat$data,
      metadata = dat$metadata,
      samplename = "sample"
    ),
    cores = 1,
    verbose = FALSE,
    return_models = FALSE,
    cache = cache
  ))
}


test_that("targets that warned are kept by default", {
  res <- sw_warn_results()
  n_warn <- nrow(seqwrap_errors(res, "fit", "warning"))
  expect_gt(n_warn, 0)

  combined <- seqwrap_summarise(res, verbose = FALSE)

  # Every target is reported, warning or not
  expect_equal(length(unique(combined$summaries$target)), res@k)
  expect_equal(length(unique(combined$evaluations$target)), res@k)
  expect_null(combined$dropped)
})


test_that("the errors element reports one row per condition", {
  res <- sw_warn_results()
  n_warn <- nrow(seqwrap_errors(res, "fit", "warning"))

  combined <- seqwrap_summarise(res, verbose = FALSE)

  expect_false(is.null(combined$errors))
  expect_identical(
    names(combined$errors),
    c("target", "stage", "type", "class", "message")
  )
  expect_equal(nrow(combined$errors), n_warn)
  expect_true(all(combined$errors$stage == "fit"))
  expect_true(all(combined$errors$type == "warning"))
  expect_true(all(grepl("simulated singular fit", combined$errors$message)))

  # It can be switched off
  expect_null(seqwrap_summarise(res, errors = FALSE, verbose = FALSE)$errors)
})


test_that("drop_warnings removes warned targets and records which", {
  res <- sw_warn_results()
  n_warn <- nrow(seqwrap_errors(res, "fit", "warning"))

  combined <- seqwrap_summarise(res, drop_warnings = TRUE, verbose = FALSE)

  expect_equal(
    length(unique(combined$summaries$target)),
    res@k - n_warn
  )
  expect_equal(
    length(unique(combined$evaluations$target)),
    res@k - n_warn
  )

  # The removal is recorded rather than silent
  expect_equal(length(combined$dropped), n_warn)
  expect_setequal(
    combined$dropped,
    seqwrap_errors(res, "fit", "warning")$target
  )

  # The full record of conditions survives the drop
  expect_equal(nrow(combined$errors), n_warn)
})


test_that("drop_warnings accepts a selection of stages", {
  res <- sw_warn_results()

  by_stage <- seqwrap_summarise(
    res,
    drop_warnings = "fit",
    verbose = FALSE
  )
  by_all <- seqwrap_summarise(res, drop_warnings = TRUE, verbose = FALSE)
  expect_setequal(by_stage$dropped, by_all$dropped)

  # Only fitting warnings were raised, so selecting another stage drops nothing
  other <- seqwrap_summarise(
    res,
    drop_warnings = c("summary", "evaluation"),
    verbose = FALSE
  )
  expect_null(other$dropped)
  expect_equal(length(unique(other$summaries$target)), res@k)

  expect_error(
    seqwrap_summarise(res, drop_warnings = "nonsense", verbose = FALSE),
    "drop_warnings"
  )
})


test_that("dropping is announced rather than silent", {
  res <- sw_warn_results()

  # seqwrap_summarise() emits several messages, so collect them all rather than
  # depending on which one arrives first
  messages <- character(0)
  withCallingHandlers(
    invisible(seqwrap_summarise(res, drop_warnings = TRUE, verbose = TRUE)),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  expect_true(any(grepl("Dropped", messages)))

  # and stays quiet when the user asked for quiet
  quiet <- character(0)
  withCallingHandlers(
    invisible(seqwrap_summarise(res, drop_warnings = TRUE, verbose = FALSE)),
    message = function(m) {
      quiet <<- c(quiet, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_equal(length(quiet), 0)
})


test_that("condition reporting works across cache modes", {
  reference <- seqwrap_summarise(sw_warn_results("none"), verbose = FALSE)

  for (cache in c("memory", "disk")) {
    res <- sw_warn_results(cache)
    combined <- seqwrap_summarise(res, verbose = FALSE)

    expect_equal(nrow(combined$errors), nrow(reference$errors))
    expect_setequal(combined$errors$target, reference$errors$target)

    dropped <- seqwrap_summarise(res, drop_warnings = TRUE, verbose = FALSE)
    expect_setequal(dropped$dropped, reference$errors$target)

    if (!is.null(res@cache)) seqwrap_cache_clear(res)
  }
})
