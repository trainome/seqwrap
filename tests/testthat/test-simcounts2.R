
library(testthat)
library(seqwrap)



test_that("simcounts2 returns the documented structure with defaults", {
  skip_on_cran() # simulation thus no need to be checked by CRAN

  #### Use variables so that we know nrows, n1 n2 etc?


  out <- simcounts2()

  # The documented return is a list with (1) counts, (2) eta, (3) phi, (4) meta
  expect_type(out, "list")
  expect_true(all(c("counts", "eta", "phi", "metadata") %in% names(out)))

  # counts: matrix of non-negative integers; dimensions: genes x samples
  expect_true(is.data.frame(out$counts)) #### Re-write with  expect_s3_class
  expect_true(is.numeric(out$counts[,2]))     ########### Re-write with apply
  expect_true(all(is.finite(out$counts))) ### Re-write with apply
  expect_true(all(out$counts >= 0))


  # allow >=0 and integer-valued
  expect_true(all(as.matrix(out$counts[,-1]) == floor(as.matrix(out$counts[,-1]))))

  # eta: matrix;
  expect_true(is.data.frame(out$eta))
  # Check that eta and counts have the same dimensions
  expect_equal(dim(out$eta), dim(out$counts))
  # Check that all eta values are numeric
  expect_true(all(apply(out$eta[,-1], 2, is.numeric)))
  # Check tha no values are missing or inf
  expect_true(all(apply(out$eta[,-1], 2, is.finite)))


  # phi: numeric vector (second column) in a data frame,
  expect_true(is.numeric(out$phi[,2]))
  expect_true(all(is.finite(out$phi[,2])))
  # Check that number of columns are greater than 0
  expect_gt(length(out$phi), 0)
  # Check number of rows
  expect_gt(nrow(out$phi), 0)


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
