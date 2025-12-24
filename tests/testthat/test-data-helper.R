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
