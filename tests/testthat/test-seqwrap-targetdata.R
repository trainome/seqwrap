
library(testthat)
library(seqwrap)

test_that("seqwrap returns a list of models when target data is used to pass
          variables used in arguments.", {

            # Simulate data
            dat <- seqwrap:::simcounts(seed = 1,
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

          })

