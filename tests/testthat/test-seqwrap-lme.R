# Use this for specific tests of lme functionality
#
#
#
library(testthat)
library(seqwrap)

test_that("seqwrap lme and gls models from the nlme package", {

            # Simulate data
            dat <- seqwrap::simcounts(seed = 1,
                                       n_samples = 40,
                                       n_genes = 10,
                                       clusters = 20)

            swobject <- seqwrap::seqwrap_compose(
              modelfun = nlme::lme,
              arguments = alist(
                fixed = y ~ x,
                random = ~ 1|cluster
              ),
              data = dat$data,
              metadata = dat$metadata,
              samplename = "sample",
              additional_vars = NULL
            )


            test_lme <- seqwrap::seqwrap(
              swobject,
              return_models = TRUE,
              cores = 1
            )

            expect_s3_class(test_lme@models[[1]], "lme")



            test_gls <- seqwrap::seqwrap(
              swobject,
              modelfun = nlme::gls,
              arguments = alist(
                model = y ~ x,
                correlation = nlme::corCompSymm(form = ~ 1|cluster) ),
              additional_vars = "cluster",
              return_models = TRUE,
              cores = 1
            )


            expect_s3_class(test_gls@models[[1]], "gls")




          })
