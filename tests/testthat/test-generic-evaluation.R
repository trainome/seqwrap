# Packages
library(testthat)
library(seqwrap)


# The columns every method must return, so that results from any engine can be
# bound together by seqwrap_summarise()
eval_schema <- c(
  "engine",
  "converged",
  "code",
  "message",
  "singular",
  "iterations"
)

sw_eval_testdata <- function() {
  set.seed(1)
  n <- 120
  d <- data.frame(
    x = rep(c("a", "b"), each = n / 2),
    g = factor(rep(1:20, length.out = n)),
    t = stats::runif(n)
  )
  d$y <- stats::rpois(n, exp(2 + 0.4 * (d$x == "b")))
  d
}


test_that("generic_evaluation returns a stable schema for base engines", {
  d <- sw_eval_testdata()

  m_lm <- stats::lm(t ~ x, data = d)
  m_glm <- stats::glm(y ~ x, data = d, family = stats::poisson())

  for (m in list(m_lm, m_glm)) {
    out <- generic_evaluation(m)
    expect_s3_class(out, "data.frame")
    expect_identical(names(out), eval_schema)
    expect_equal(nrow(out), 1)
    expect_type(out$engine, "character")
    expect_type(out$converged, "logical")
    expect_type(out$code, "integer")
    expect_type(out$message, "character")
    expect_type(out$singular, "logical")
    expect_type(out$iterations, "integer")
  }

  expect_identical(generic_evaluation(m_lm)$engine, "lm")
  expect_identical(generic_evaluation(m_glm)$engine, "glm")
  expect_true(generic_evaluation(m_glm)$converged)
})


test_that("rows from different engines bind together", {
  d <- sw_eval_testdata()

  rows <- list(
    generic_evaluation(stats::lm(t ~ x, data = d)),
    generic_evaluation(stats::glm(y ~ x, data = d, family = stats::poisson()))
  )

  bound <- do.call(rbind, rows)
  expect_equal(nrow(bound), 2)
  expect_identical(names(bound), eval_schema)
})


test_that("generic_evaluation returns NULL for a failed fit", {
  expect_null(generic_evaluation(NULL))
})


test_that("unknown engines fall through to the default method", {
  weird <- structure(list(converged = FALSE), class = "not_a_real_engine")

  out <- generic_evaluation(weird)
  expect_identical(names(out), eval_schema)
  expect_identical(out$engine, "not_a_real_engine")
  expect_false(out$converged)

  # An object with no convergence information at all still returns a row
  bare <- structure(list(), class = "engine_without_diagnostics")
  expect_true(is.na(generic_evaluation(bare)$converged))
})


test_that("rank deficient fits are flagged as singular", {
  d <- sw_eval_testdata()
  d$dup <- as.numeric(d$x == "b")

  # x and dup are collinear, so one coefficient is aliased
  out <- generic_evaluation(stats::lm(t ~ x + dup, data = d))

  expect_true(out$singular)
  expect_match(out$message, "rank deficient")
})


test_that("glmmTMB diagnostics separate convergence from a degenerate fit", {
  skip_if_not_installed("glmmTMB")

  d <- sw_eval_testdata()
  m <- suppressWarnings(glmmTMB::glmmTMB(
    y ~ x + (1 | g),
    data = d,
    family = glmmTMB::nbinom2
  ))

  out <- generic_evaluation(m)
  expect_identical(names(out), eval_schema)
  expect_identical(out$engine, "glmmTMB")
  expect_identical(out$code, as.integer(m$fit$convergence))

  # singular tracks the positive definiteness of the Hessian, which glmmTMB
  # reports independently of the optimiser's own convergence code
  expect_identical(out$singular, !isTRUE(m$sdr$pdHess))
})


test_that("lme4 diagnostics report singularity alongside code 0", {
  skip_if_not_installed("lme4")

  d <- sw_eval_testdata()
  # A response that is constant within groups produces a singular fit
  d$flat <- rep(1:6, length.out = nrow(d)) + 0.0
  m <- suppressMessages(lme4::lmer(flat ~ x + (1 | g), data = d))

  out <- generic_evaluation(m)
  expect_identical(names(out), eval_schema)
  expect_identical(out$engine, "lmerMod")

  # The optimiser is satisfied, but the fit is singular. This is precisely why
  # the two are reported separately.
  expect_identical(out$code, 0L)
  expect_true(out$singular)
  expect_false(is.na(out$message))
})


test_that("a single merMod method covers lmer and glmer", {
  skip_if_not_installed("lme4")

  d <- sw_eval_testdata()
  m_lmer <- suppressMessages(lme4::lmer(t ~ x + (1 | g), data = d))
  m_glmer <- suppressMessages(lme4::glmer(
    y ~ x + (1 | g),
    data = d,
    family = stats::poisson()
  ))

  expect_identical(generic_evaluation(m_lmer)$engine, "lmerMod")
  expect_identical(generic_evaluation(m_glmer)$engine, "glmerMod")
  expect_identical(names(generic_evaluation(m_glmer)), eval_schema)
})


test_that("nlme diagnostics are reported", {
  skip_if_not_installed("nlme")

  d <- sw_eval_testdata()
  m_lme <- nlme::lme(t ~ x, random = ~ 1 | g, data = d)
  m_gls <- nlme::gls(t ~ x, data = d)

  expect_identical(generic_evaluation(m_lme)$engine, "lme")
  expect_identical(generic_evaluation(m_gls)$engine, "gls")
  expect_identical(names(generic_evaluation(m_lme)), eval_schema)
})


test_that("glm.nb reports theta estimation problems", {
  skip_if_not_installed("MASS")

  d <- sw_eval_testdata()
  m <- suppressWarnings(MASS::glm.nb(y ~ x, data = d))

  out <- generic_evaluation(m)
  expect_identical(out$engine, "negbin")
  expect_identical(names(out), eval_schema)

  # theta is estimated in a loop separate from the IRLS step and can fail on
  # its own, which has to be reflected in the reported convergence
  if (!is.null(m$th.warn)) {
    expect_false(out$converged)
    expect_match(out$message, "theta")
  }
})


test_that("mgcv diagnostics cover both fitting stages", {
  skip_if_not_installed("mgcv")

  d <- sw_eval_testdata()
  d$z <- d$t + stats::rnorm(nrow(d))
  m <- mgcv::gam(t ~ s(z), data = d)

  out <- generic_evaluation(m)
  expect_identical(names(out), eval_schema)
  expect_identical(out$engine, "gam")

  # converged combines the penalised IRLS step with the smoothing parameter
  # optimisation, and singular reports the Hessian
  expect_true(out$converged)
  expect_false(out$singular)
})


test_that("evaluations are produced through seqwrap, including in parallel", {
  dat <- seqwrap::simcounts(
    seed = 1,
    n_samples = 30,
    n_genes = 8,
    clusters = 10
  )

  container <- seqwrap::seqwrap_compose(
    modelfun = stats::lm,
    arguments = list(formula = y ~ x),
    data = dat$data,
    metadata = dat$metadata,
    samplename = "sample"
  )

  # S3 dispatch has to resolve on the workers, not only in this process
  for (cores in c(1, 2)) {
    res <- seqwrap::seqwrap(
      container,
      cores = cores,
      return_models = FALSE,
      verbose = FALSE
    )
    evaluations <- seqwrap_summarise(res, verbose = FALSE)$evaluations

    expect_equal(nrow(evaluations), 8)
    expect_true(all(evaluations$engine == "lm"))
    expect_true(all(c("target", eval_schema) %in% names(evaluations)))
  }
})


test_that("residual_diagnostics reproduces the previous DHARMa behaviour", {
  skip_if_not_installed("DHARMa")

  d <- sw_eval_testdata()
  m <- stats::glm(y ~ x, data = d, family = stats::poisson())

  out <- suppressMessages(residual_diagnostics(m, n = 50))

  expect_identical(names(out), c("uniformity", "dispersion", "outliers"))
  expect_equal(nrow(out), 1)
  expect_null(residual_diagnostics(NULL))
})
