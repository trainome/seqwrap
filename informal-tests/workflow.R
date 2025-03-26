library(seqwrap)



# Populate a seqwrap container
test <- seqwrap_compose(modelfun = stats::lm,
                        arguments = list(formula = y ~ x),
                        data = data.frame(targetid = paste0("t", 1:10),
                                          x = rpois(10,5),
                                          y = rpois(10,5),
                                          z = rpois(10,5),
                                          k = rpois(10,5),
                                          l = rpois(10,5),
                                          m = rpois(10,5)),
                        targetdata = NULL,
                        metadata = data.frame(seq_sample_id = c("x",
                                                                "y",
                                                                "z",
                                                                "k",
                                                                "l",
                                                                "m"),
                                              x = c(0,0,0,1,1,1)),
                        samplename = "seq_sample_id")

# Run the seqwrap container
seqwrap(y = test)

