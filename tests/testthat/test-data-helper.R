test_that("data_helper returns a list when given a formatted data frame or list.", {
  local({
    # Set up data frames and lists for data_helper tests

   expect_type(seqwrap:::data_helper(dat = data.frame(
     gene_id = c("aaa", "bbb"),
     sample1 = c(21, 26),
     sample2 = c(11, 22)
   )), "list")

   expect_type(seqwrap:::data_helper(dat = data.frame(
     row.names = c("aaa", "bbb"),
     sample1 = c(21, 26),
     sample2 = c(11, 22)
   ), rownames = TRUE), "list")


    expect_type(seqwrap:::data_helper(dat = list(
        y1 = data.frame(
          row.names = c("aaa", "bbb"),
          sample1 = c(21, 26),
          sample2 = c(11, 22)
          ),
        y2 = data.frame(
          row.names = c("aaa", "bbb"),
          sample1 = c(21, 26),
          sample2 = c(11, 22)
        )), rownames = TRUE), "list")
  })
})

test_that("data_helper throws error on duplicate targets in data frame", {
  # Create a data frame with duplicate target IDs in the first column
  dat_with_duplicates <- data.frame(
    gene_id = c("aaa", "bbb", "aaa", "ccc"),  # "aaa" appears twice
    sample1 = c(21, 26, 30, 15),
    sample2 = c(11, 22, 33, 44)
  )

  # Expect error with informative message about duplicate target
  expect_error(
    seqwrap:::data_helper(dat = dat_with_duplicates),
    "Duplicate target found: 'aaa' appears 2 times in the first column"
  )

  # Test with multiple duplicates
  dat_multiple_duplicates <- data.frame(
    gene_id = c("aaa", "bbb", "aaa", "bbb", "aaa"),
    sample1 = c(1, 2, 3, 4, 5),
    sample2 = c(10, 20, 30, 40, 50)
  )

  # Should error on first duplicate encountered (either "aaa" or "bbb")
  expect_error(
    seqwrap:::data_helper(dat = dat_multiple_duplicates),
    "Duplicate target found"
  )
})
