# Guards on the shapes that user code depends on. These exist so that an
# accidental rename or reordering during refactoring fails here rather than in
# somebody's analysis script.

library(testthat)
library(seqwrap)


sw_contract_results <- function(...) {
  dat <- seqwrap::simcounts(
    seed = 1,
    n_samples = 20,
    n_genes = 6,
    clusters = 5
  )

  seqwrap::seqwrap(
    seqwrap::seqwrap_compose(
      modelfun = stats::lm,
      arguments = list(formula = y ~ x),
      data = dat$data,
      metadata = dat$metadata,
      samplename = "sample"
    ),
    cores = 1,
    verbose = FALSE,
    ...
  )
}


test_that("seqwrapResults exposes the documented properties", {
  res <- sw_contract_results()

  expect_true(all(
    c(
      "models", "summaries", "evaluations", "errors", "targets",
      "cache", "n", "k", "call_arguments", "call_engine", "elapsed_time"
    ) %in% names(S7::props(res))
  ))

  expect_type(res@targets, "character")
  expect_s3_class(res@errors, "data.frame")
  expect_equal(res@k, length(res@targets))
})


test_that("the errors table keeps its documented columns", {
  res <- sw_contract_results()

  expect_identical(
    names(res@errors),
    c("target", "stage", "type", "class", "message")
  )

  # Same columns whether or not anything was raised
  expect_identical(names(seqwrap_errors(res)), names(res@errors))

  # And the filters only ever return the documented levels
  all_conditions <- seqwrap_errors(res)
  expect_true(all(all_conditions$stage %in%
    c("fit", "summary", "evaluation")))
  expect_true(all(all_conditions$type %in% c("error", "warning")))
})


test_that("summaries and evaluations lead with the target column", {
  res <- sw_contract_results()

  expect_identical(names(res@summaries)[1], "target")
  expect_identical(names(res@evaluations)[1], "target")
  expect_setequal(unique(res@summaries$target), res@targets)
})


test_that("seqwrap_summarise returns the documented elements", {
  res <- sw_contract_results()
  combined <- seqwrap_summarise(res, verbose = FALSE)

  expect_type(combined, "list")
  expect_true(all(names(combined) %in%
    c("summaries", "evaluations", "errors", "dropped")))
  expect_s3_class(combined$summaries, "data.frame")
  expect_identical(names(combined$summaries)[1], "target")

  # Switching a component off removes it rather than returning an empty one
  only_summaries <- seqwrap_summarise(
    res,
    evaluations = FALSE,
    errors = FALSE,
    verbose = FALSE
  )
  expect_null(only_summaries$evaluations)
  expect_null(only_summaries$errors)
})


test_that("generic_evaluation keeps its documented columns", {
  res <- sw_contract_results()

  expect_true(all(
    c("engine", "converged", "code", "message", "singular", "iterations") %in%
      names(res@evaluations)
  ))
})


test_that("every cache mode yields the same combined results", {
  reference <- seqwrap_summarise(
    sw_contract_results(cache = "none"),
    verbose = FALSE
  )

  for (mode in c("memory", "disk")) {
    res <- sw_contract_results(cache = mode)
    combined <- seqwrap_summarise(res, verbose = FALSE)

    expect_equal(combined$summaries, reference$summaries)
    expect_equal(combined$evaluations, reference$evaluations)

    if (!is.null(res@cache)) seqwrap_cache_clear(res)
  }
})


test_that("seqwrap_cache_clear is a no-op without a cache", {
  res <- sw_contract_results()

  expect_equal(seqwrap_cache_clear(res), 0L)
  expect_error(seqwrap_cache_clear("not a results object"), "seqwrapResults")
})


test_that("seqwrap_errors validates its filters", {
  res <- sw_contract_results()

  expect_error(seqwrap_errors(res, stage = "nonsense"), "stage")
  expect_error(seqwrap_errors(res, type = "nonsense"), "type")
  expect_error(seqwrap_errors("not a results object"), "seqwrapResults")
})


test_that("summarise refuses input that is not a seqwrapResults object", {
  expect_error(seqwrap_summarise(list()), "seqwrapResults")
  expect_error(seqwrap_summarise(data.frame()), "seqwrapResults")
})
