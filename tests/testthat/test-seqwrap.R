# Packages
library(testthat)
library(seqwrap)


test_that("seqwrap returns a list of models in the model
          slot when asked to return models", {

    # Simulate data
    dat <- seqwrap::simcounts(seed = 1,
                               n_samples = 40,
                               n_genes = 10,
                               clusters = 20)

    swobject <- seqwrap::seqwrap_compose(
      modelfun = glmmTMB::glmmTMB,
      arguments = list(
        formula = y ~
          x +
          (1 | cluster),
        family = glmmTMB::nbinom1
      ),
      data = dat$data,
      metadata = dat$metadata,
      samplename = "sample",
      additional_vars = NULL
    )





    test_glmmtmb <- seqwrap::seqwrap(swobject,
                                     return_models = TRUE,
                                     cores = 1)



    expect_s3_class(test_glmmtmb@models[[1]], "glmmTMB")

  test_glmnb <- seqwrap::seqwrap(
    swobject,
    arguments = list(formula = y ~ x),
    modelfun = MASS::glm.nb,
    return_models = TRUE,
    cores = 1)

  expect_s3_class(test_glmnb@models[[1]], c("glm", "lm", "negbin"))







  test_lm <- seqwrap::seqwrap(
    swobject,
    modelfun = stats::lm,
    arguments = list(formula = y ~ x),
    return_models = TRUE,
    cores = 1
  )

  expect_s3_class(test_lm@models[[1]], "lm")

  test_glm <- seqwrap::seqwrap(
    swobject,
    modelfun = stats::glm,
    arguments = alist(
      formula = y ~ x,
      family = poisson(link = "log")
    ),
    return_models = TRUE,
    cores = 1
  )

  expect_s3_class(test_glm@models[[1]], "glm")



  test_lme <- seqwrap::seqwrap(
    swobject,
    modelfun = nlme::lme,
    arguments = alist(
      fixed = y ~ x,
      random = ~ 1 | cluster
    ),
    return_models = TRUE,
    cores = 1
  )


  expect_s3_class(test_lme@models[[1]], "lme")


          })


test_that("Model summaries and evaluations returns expected results", {
  ## Create a model summary function
  summaryfun_glmmtmb <- function(x) {
    # Extract conditional effects and store as a tibble
    results <- tibble::as_tibble(coef(summary(x))$cond, rownames = "coef") |>
      dplyr::select(
        coef,
        estimate = Estimate,
        se = "Std. Error",
        z = "z value",
        p = "Pr(>|z|)"
      )

    ## Return results
    return(results)
  }

  evalfun_glmmtmb <- function(x) {
    simresid <- DHARMa::simulateResiduals(
      fittedModel = x,
      plot = FALSE,
      n = 1000
    )

    unif <- DHARMa::testUniformity(simresid, plot = FALSE)
    disp <- DHARMa::testDispersion(simresid, plot = FALSE)
    pdhess <- x$sdr$pdHess

    return(data.frame(
      unif = unif$p.value,
      disp = disp$p.value,
      pdhess = pdhess
    ))
  }


  # Simulate data
  dat <- seqwrap:::simcounts(seed = 1,
                             n_samples = 40,
                             n_genes = 10,
                             clusters = 20)

  swobject <- seqwrap::seqwrap_compose(
    modelfun = glmmTMB::glmmTMB,
    arguments = list(
      formula = y ~
        x +
        (1 | cluster),
      family = glmmTMB::nbinom1
    ),
    data = dat$data,
    metadata = dat$metadata,
    samplename = "sample",
    additional_vars = NULL
  )



  testsummary_glmmtmb <- seqwrap::seqwrap(
    swobject,
    summary_fun = summaryfun_glmmtmb,
    eval_fun = evalfun_glmmtmb,
    cores = 1
  )

  ## Expected output from summary/evaluation functions
  expect_s3_class(testsummary_glmmtmb@summaries[[1]], "tbl")
  expect_s3_class(testsummary_glmmtmb@evaluations[[1]], "data.frame")

  ## Expect no errors in the error data
  expect_null(testsummary_glmmtmb@errors$err_sum[[1]])
  expect_null(testsummary_glmmtmb@errors$err_eval[[1]])

  # TODO add tests for the other model types
})


test_that("seqwrap is silent when setting verbose = FALSE", {
  # Simulate data
  dat <- seqwrap::simcounts(seed = 1,
                            n_samples = 40,
                            n_genes = 10,
                            clusters = 20)

  swobject <- seqwrap::seqwrap_compose(
    modelfun = glmmTMB::glmmTMB,
    arguments = list(
      formula = y ~
        x +
        (1 | cluster),
      family = glmmTMB::nbinom1
    ),
    data = dat$data,
    metadata = dat$metadata,
    samplename = "sample",
    additional_vars = NULL
  )


  expect_no_message({
    test_glmmtmb <- seqwrap::seqwrap(swobject,
                                     return_models = TRUE,
                                     verbose = FALSE,
                                     cores = 1)
    })


  expect_no_message({
    temp <-  seqwrap_summarise(test_glmmtmb, verbose = FALSE)
    })


          })





