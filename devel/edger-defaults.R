

# Download data and do salmon alignment
library(GEOquery)


# Get the GEO dataset information
gse <- getGEO("GSE172421", GSEMatrix = FALSE)




library(edgeR)
library(seqwarp)
library(tidyverse)


seqwrap::rna_seq_sample

meta <- seqwrap::rna_seq_metadata |>
  filter(!is.na(seq_sample_id),
         leg == "R",
         time == "w0") |>
  print()


data <- seqwrap::rna_seq_sample |>
  select(transcript_id, all_of(meta$seq_sample_id)) |>
  mutate(across(all_of(meta$seq_sample_id), as.integer)) |>
  print()



y <- DGEList(counts = data)

y$samples$sex <- ifelse(meta$sex == "male", 0, 1)

y$samples$group <- ifelse(y$samples$sex == 0, 1, 2)

keep <- filterByExpr(y, group = y$samples$group)

y <- y[keep,, keep.lib.sizes = FALSE]


y <- normLibSizes(y)

y$samples

design <- model.matrix(~y$samples$sex)

y <- estimateDisp(y, design)



fit <- glmQLFit(y, design)


fit$coefficients[,2]

test <- glmQLFTest(fit, coef =2)

topTags(test)




