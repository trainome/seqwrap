
library(testthat)
library(seqwrap)


test_that("simcounts returns the documented structure with defaults", {
  skip_on_cran()
  out <- simcounts()

  # Return should be a list with data, parameters, metadata
  expect_type(out, "list")
  expect_true(all(c("data", "parameters", "metadata") %in% names(out)))

  # data: data.frame with gene rows and sample columns (plus 'targetid')
  expect_s3_class(out$data, "data.frame")
  expect_true("targetid" %in% names(out$data))
  expect_equal(nrow(out$data), 10)           # default n_genes
  expect_equal(ncol(out$data), 1 + 16)       # 1 targetid + default n_samples

  # parameters: one row per gene with beta_0, beta_1, overdispersion
  expect_s3_class(out$parameters, "data.frame")
  expect_true(all(c("targetid", "beta_0", "beta_1", "overdispersion") %in% names(out$parameters)))
  expect_equal(nrow(out$parameters), 10)
  expect_true(is.numeric(out$parameters$beta_0))
  expect_true(is.numeric(out$parameters$beta_1))
  expect_true(is.numeric(out$parameters$overdispersion))
  expect_true(all(is.finite(out$parameters$overdispersion)))
  expect_true(all(out$parameters$overdispersion >= 1 & out$parameters$overdispersion <= 10))

  # metadata: one row per sample with sample, cluster, x
  expect_s3_class(out$metadata, "data.frame")
  expect_true(all(c("sample", "cluster", "x") %in% names(out$metadata)))
  expect_equal(nrow(out$metadata), 16)

  # Column names of data (excluding targetid) should match metadata$sample
  expect_identical(colnames(out$data)[-1], out$metadata$sample)

  # Counts should be non-negative integer-like
  counts <- as.matrix(out$data[,-1])
  expect_true(is.numeric(counts))
  expect_true(all(counts >= 0))
  expect_true(all(counts == floor(counts)))
})

test_that("dimensions match specified n_genes and n_samples", {
  n_genes <- 6
  n_samples <- 12
  clusters <- 6  # must divide n_samples

  out <- simcounts(
    n_genes = n_genes,
    n_samples = n_samples,
    clusters = clusters,
    seed = 1
  )

  # data dims: genes x samples (+ targetid column)
  expect_equal(nrow(out$data), n_genes)
  expect_equal(ncol(out$data), 1 + n_samples)

  # parameters: one row per gene
  expect_equal(nrow(out$parameters), n_genes)

  # metadata: one row per sample
  expect_equal(nrow(out$metadata), n_samples)
})

