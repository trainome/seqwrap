# Preparing example counts for the package tutorials etc.
#
# The data used in presentation of the package comes from
# Pillon et al. 2022. This script downloads data and stores it in
# two data frames pillon_count and pillon_metadata.

library(GEOquery)

# Get the gene counts
counts <- getRNASeqData("GSE202295")

# Sort row annotations and convert to symbols
annotation <- SummarizedExperiment::rowData(counts) |>
  data.frame() |>
  dplyr::select(GeneID, Symbol) |>
  dplyr::mutate(GeneID = as.character(GeneID)) |>
  data.frame(row.names = NULL)

# Get the count matrix and add symbols
counts_mat <- counts@assays@data@listData$counts

countdata <- counts_mat |>
  data.frame() |>
  tibble::rownames_to_column("GeneID") |>
  dplyr::inner_join(annotation) |>
  dplyr::select(geneid = Symbol, tidyselect::starts_with("GSM"))

# Remove all zero row counts
countdata <- countdata[rowSums(countdata[, -1]) != 0, ]


# Sort meta data
metadata <- counts@colData |>
  data.frame() |>
  dplyr::select(seq_sample_id = geo_accession, group_id_time = title) |>
  tidyr::separate(group_id_time, into = c("group", "id", "time")) |>
  dplyr::mutate(time = factor(time, levels = c("basal", "post", "rec"))) |>
  data.frame(row.names = NULL)


# Save the data as one file
pillon_counts <- list(metadata = metadata, countdata = countdata)

usethis::use_data(pillon_counts, compress = "xz", overwrite = TRUE)
