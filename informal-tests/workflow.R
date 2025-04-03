library(seqwrap)


## A simple summary function

sumfun <- function(x) {

  df <- data.frame(coef(summary(x)))

  df$coef <- rownames(df)
  rownames(df) <- NULL
  colnames(df) <- c("estimate", "se", "tval", "pval", "coef")


  df <- df[,c("coef","estimate", "se", "tval", "pval")]



  return(df)

}






# Populate a seqwrap container
test <- seqwrap_compose(modelfun = stats::glm,
                        arguments = list(formula = y ~ x,
                                         family = poisson(link = "log")),
                        data = data.frame(targetid = paste0("t", 1:10),
                                          x = rpois(10,5),
                                          y = rpois(10,5),
                                          z = rpois(10,5),
                                          k = rpois(10,5),
                                          l = rpois(10,5),
                                          m = rpois(10,5)),
                        targetdata = NULL,
                        summary_fun = sumfun,
                        metadata = data.frame(seq_sample_id = c("x",
                                                                "y",
                                                                "z",
                                                                "k",
                                                                "l",
                                                                "m"),
                                              x = c(0,0,0,1,1,1)),
                        samplename = "seq_sample_id")

# Run the seqwrap container


ls_names <- names(test@arguments)



substitute(test@arguments)

testres <- seqwrap(y = test, return_models = FALSE)

testres@ca


testres@summaries

results <- seqwrap_summarise(testres)



results$summaries













ls <- list(disp = c(1, 2, 3),
           unf = c(3, 4, 5))

subset <- 1:2

lapply(ls, function(x) x[subset])



temp <- seqwrap:::data_helper(
  dat = dat2[1:10,],
  targetdat = dispersions,
  rownames = FALSE
)


target.wise <- temp[[1]][[2]]

list2env(target.wise, envir = .GlobalEnv)




rm(disp)
target.wise[[disp]]



























