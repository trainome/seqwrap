
library(testthat)



test_that("simcounts2 returns the documented structure with defaults", {
  skip_on_cran() # simulation thus no need to be checked by CRAN
  out <- simcounts2()

  # The documented return is a list with (1) counts, (2) eta, (3) phi, (4) meta
  expect_type(out, "list")
  expect_true(all(c("counts", "eta", "phi", "meta") %in% names(out)))

  # counts: matrix of non-negative integers; dimensions: genes x samples
  expect_true(is.matrix(out$counts))
  expect_true(is.numeric(out$counts))
  expect_true(all(is.finite(out$counts)))
  expect_true(all(out$counts >= 0))


  # allow >=0 and integer-valued
  expect_true(all(out$counts == floor(out$counts)))

  # eta: matrix;
  expect_true(is.matrix(out$eta))
  expect_equal(dim(out$eta), dim(out$counts))
  expect_true(is.numeric(out$eta))
  expect_true(all(is.finite(out$eta)))

  # phi: numeric vector
  expect_true(is.numeric(out$phi))
  expect_true(all(is.finite(out$phi)))
  expect_gt(length(out$phi), 0)

  # meta: data.frame with one row per sample
  expect_s3_class(out$meta, "data.frame")
  # expected samples: (n1 + n2) * 3 from defaults (5+5)*3 = 30
  expect_equal(nrow(out$meta), (5 + 5) * 3)

  # The number of columns in counts should match number of samples
  expect_equal(ncol(out$counts), nrow(out$meta))

  # The number of genes should match length(beta0) (default beta0 length = 2)
  expect_equal(nrow(out$counts), length(c(2.3, 3)))

  # If library sizes are included in metadata, sanity-check them
  if ("library_size" %in% names(out$meta)) {
    expect_true(is.numeric(out$meta$library_size))
    expect_true(all(out$meta$library_size > 0))
    expect_equal(length(out$meta$library_size), nrow(out$meta))
  }
})
