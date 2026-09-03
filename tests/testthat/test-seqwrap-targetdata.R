
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





            expect_warning(
              test_glmmtmb <- seqwrap::seqwrap(swobject, return_models = TRUE,
                                               cores = 1),
              "matched to targets by position")

            # Subset test
            expect_warning(
              test_glmmtmb_subset <- seqwrap::seqwrap(swobject, return_models = TRUE,
                                                      subset = 1:5, cores = 1),
              "matched to targets by position")

            expect_s3_class(test_glmmtmb@models[[1]], "glmmTMB")
            expect_s3_class(test_glmmtmb_subset@models[[1]], "glmmTMB")

            # Check that glmmTMB has used priors
            expect_no_error(test_glmmtmb@models[[1]]$modelInfo$priors)
            expect_no_error(test_glmmtmb_subset@models[[1]]$modelInfo$priors)


          })









test_that("target data is paired with targets by data row, not sorted order", {

  dat <- seqwrap::simcounts(seed = 1, n_genes = 6, n_samples = 8)
  counts <- dat$data
  # Identifiers whose row order differs from their sorted order
  counts[[1]] <- c("g3", "g1", "g6", "g2", "g5", "g4")

  # A fitting function that keeps the target-wise value it was given
  tagged_lm <- function(formula, data, tag) {
    m <- stats::lm(formula, data = data)
    m$tag <- tag
    m
  }

  container <- seqwrap::seqwrap_compose(
    modelfun = tagged_lm,
    arguments = alist(formula = y ~ x, tag = tag),
    data = counts,
    metadata = dat$metadata,
    samplename = "sample",
    targetdata = data.frame(tag = counts[[1]]),
    summary_fun = function(m) data.frame(tag = m$tag),
    eval_fun = function(m) NULL
  )

  # Data frame target data
  res <- seqwrap::seqwrap(container, cores = 1, verbose = FALSE)
  expect_identical(res@summaries$tag, res@summaries$target)

  # An unnamed list is matched by position, with a warning
  as_list <- lapply(counts[[1]], function(z) data.frame(tag = z))
  expect_warning(
    res_list <- seqwrap::seqwrap(
      seqwrap::seqwrap_compose(container, update = list(targetdata = as_list)),
      cores = 1, verbose = FALSE
    ),
    "matched to targets by position"
  )
  expect_identical(res_list@summaries$tag, res_list@summaries$target)

  # A list named by target identifier is matched by name, whatever its order
  named <- stats::setNames(as_list, counts[[1]])[c(6, 2, 4, 1, 5, 3)]
  expect_no_warning(
    res_named <- seqwrap::seqwrap(
      seqwrap::seqwrap_compose(container, update = list(targetdata = named)),
      cores = 1, verbose = FALSE
    )
  )
  expect_identical(res_named@summaries$tag, res_named@summaries$target)

  # A named list may hold more targets than the data, and is subset by name
  extra <- c(named, list(g9 = data.frame(tag = "g9")))
  res_extra <- seqwrap::seqwrap(
    seqwrap::seqwrap_compose(container, update = list(targetdata = extra)),
    subset = c(2, 4), cores = 1, verbose = FALSE
  )
  expect_identical(res_extra@summaries$tag, res_extra@summaries$target)
  expect_setequal(res_extra@summaries$target, c("g1", "g2"))

  # A named list missing a target is rejected
  expect_error(
    seqwrap::seqwrap_compose(container, update = list(targetdata = named[-1])),
    "no element"
  )
  expect_error(
    seqwrap::seqwrap(
      seqwrap::seqwrap_compose(container, update = list(targetdata = named[-1])),
      cores = 1, verbose = FALSE
    ),
    "no element"
  )

  # Chunked runs and subsets keep the pairing
  res_chunk <- seqwrap::seqwrap(container, chunk_size = 2, cores = 1,
                                verbose = FALSE)
  expect_identical(res_chunk@summaries$tag, res_chunk@summaries$target)

  res_sub <- seqwrap::seqwrap(container, subset = c(1, 3, 5), cores = 1,
                              verbose = FALSE)
  expect_identical(res_sub@summaries$tag, res_sub@summaries$target)
  expect_setequal(res_sub@summaries$target, c("g3", "g6", "g5"))
})
