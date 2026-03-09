library(GEOquery)
library(tidyverse)

# The data set is available at GEO accession number GSE195707
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE195nnn/GSE195707/suppl/GSE195707%5Fcount%5Fmatrix%2Etxt%2Egz"

# Download the count file
temp <- tempfile(fileext = ".txt.gz")
download.file(url, temp, mode = "wb")
counts <- read.table(temp, header = TRUE, sep = "\t")
unlink(temp)

# Process the counts
counts <- counts |>
  # Change estimated counts to integers
  mutate(across(where(is.numeric), ~ as.integer(round(.)))) |>
  # Filter rows with min less than 10
  filter(apply(across(where(is.numeric)), 1, min) >= 10)


# Create a meta data df
gse <- getGEO("GSE195707", GSEMatrix = TRUE)


metadata <- data.frame(gse$GSE195707_series_matrix.txt.gz, row.names = NULL) |>
  dplyr::select(
    seq_sample_id = description,
    treatment = characteristics_ch1.1,
    surgery = characteristics_ch1.2
  ) |>
  mutate(
    treatment = gsub("treatment: ", "", treatment),
    surgery = gsub("surgery: ", "", surgery),
    treatment = factor(treatment, levels = c("Vehicle", "Senolytic")),
    surgery = factor(surgery, levels = c("Sham", "Overload"))
  )


# Save the data as one file
dungan_counts <- list(metadata = metadata, countdata = counts)

usethis::use_data(dungan_counts, compress = "xz", overwrite = TRUE)
