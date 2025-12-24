
library(testthat)
library(seqwrap)


test_that("seqwrap accepts NULL (the default) as target data", {

  # Simulate data
  dat <- seqwrap::simcounts(seed = 1,
                            n_samples = 40,
                            n_genes = 10,
                            clusters = 20)


 expect_no_error(
   seqwrap::seqwrap_compose(
   modelfun = glmmTMB::glmmTMB,
   # No targetdata -- Explicit NULL
   targetdata = NULL,
   arguments = alist(
     formula = y ~
       x +
       (1 | cluster),
     family = glmmTMB::nbinom2,
     priors = data.frame(prior = Prior,
                         class = Class,
                         coef = Coef)
   ),
   data = dat$data,
   metadata = dat$metadata,
   samplename = "sample",
   additional_vars = NULL
 ))


 expect_no_error(
   seqwrap::seqwrap_compose(
     modelfun = glmmTMB::glmmTMB,
  # No target data -- Implicit null
  # targetdata = NULL,
     arguments = alist(
       formula = y ~
         x +
         (1 | cluster),
       family = glmmTMB::nbinom2,
       priors = data.frame(prior = Prior,
                           class = Class,
                           coef = Coef)
     ),
     data = dat$data,
     metadata = dat$metadata,
     samplename = "sample",
     additional_vars = NULL
   ))



})








test_that("seqwrap returns a list of models when target data (data frame) is used to pass
          variables used in arguments.", {

            # Simulate data
            dat <- seqwrap::simcounts(seed = 1,
                                       n_samples = 40,
                                       n_genes = 10,
                                       clusters = 20)

            # Target data created with variables
            tardat <- data.frame(Prior = paste0("normal(",
                                                rnorm(10, 2, 0.1),
                                                ",",
                                                rnorm(10, 0.3, 0.001),
                                                ")"),
                                 Class = rep("fixef_disp", 10),
                                 Coef = rep(1, 10))


            swobject <- seqwrap::seqwrap_compose(
              modelfun = glmmTMB::glmmTMB,
              # Target data with multiple variables
              targetdata = tardat,
              arguments = alist(
                formula = y ~
                  x +
                  (1 | cluster),
                family = glmmTMB::nbinom2,
                priors = data.frame(prior = Prior,
                                   class = Class,
                                   coef = Coef)
              ),
              data = dat$data,
              metadata = dat$metadata,
              samplename = "sample",
              additional_vars = NULL
            )





            test_glmmtmb <- seqwrap::seqwrap(swobject,
                                             return_models = TRUE,
                                             cores = 1)

            # Subset test
            test_glmmtmb_subset <- seqwrap::seqwrap(swobject,
                                             return_models = TRUE,
                                             subset = 1:5,
                                             cores = 1)

            expect_s3_class(test_glmmtmb@models[[1]], "glmmTMB")
            expect_s3_class(test_glmmtmb_subset@models[[1]], "glmmTMB")

            # Check that glmmTMB has used priors
            expect_no_error(test_glmmtmb@models[[1]]$modelInfo$priors)
            expect_no_error(test_glmmtmb_subset@models[[1]]$modelInfo$priors)



          })


test_that("seqwrap returns a list of models when target data (list) is used
to pass variables used in arguments.", {

            # Simulate data
            dat <- seqwrap::simcounts(seed = 1,
                                      n_samples = 40,
                                      n_genes = 10,
                                      clusters = 20)

            # Target data created with variables
            tardat <- data.frame(prior = c("normal(0, 1)","normal(2, 1)"),
                                 class = c("beta", "fixef_disp"),
                                 coef = c("x", "1"))

            tardat_list <- list()

            for(i in 1:10) {
              tardat_list[[i]] <- tardat
            }




            swobject <- seqwrap::seqwrap_compose(
              modelfun = glmmTMB::glmmTMB,
              # Target data with multiple variables
              targetdata = tardat_list,
              arguments = alist(
                formula = y ~
                  x +
                  (1 | cluster),
                family = glmmTMB::nbinom2,
                priors = data.frame(prior = prior,
                                    class = class,
                                    coef = coef)
              ),
              data = dat$data,
              metadata = dat$metadata,
              samplename = "sample",
              additional_vars = NULL
            )





            test_glmmtmb <- seqwrap::seqwrap(swobject,
                                             return_models = TRUE,
                                             cores = 1)

            # Subset test
            test_glmmtmb_subset <- seqwrap::seqwrap(swobject,
                                                    return_models = TRUE,
                                                    subset = 1:5,
                                                    cores = 1)

            expect_s3_class(test_glmmtmb@models[[1]], "glmmTMB")
            expect_s3_class(test_glmmtmb_subset@models[[1]], "glmmTMB")

            # Check that glmmTMB has used priors
            expect_no_error(test_glmmtmb@models[[1]]$modelInfo$priors)
            expect_no_error(test_glmmtmb_subset@models[[1]]$modelInfo$priors)


          })







