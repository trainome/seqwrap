#' Download format gene expression data from Gene Expression Omnibus
#'
#' Currently, a data set from Shetty et al. is made available for the user
#' after downloading from the Gene Expression Omnibus. The dataset can be used
#' for replicating analyses in the Vignette and upcoming paper.
#'
#' @param dir Directory for downloading temporary files. Defaults to "data-raw"
#' @param gse GSE accession number
#' @references Shetty AC, Sivinski J, Cornell J, McCracken C et al. Peripheral
#' blood transcriptomic profiling of molecular mechanisms commonly
#' regulated by binge drinking and placebo effects.
#' \emph{Sci Rep} 2024 May 10;14(1):10733. PMID: 38730024
#' @export
geo_data <- function(dir = "./data-raw"){

  if (!dir.exists(dir)) dir.create(dir)


  # Raw count matrix and annotation from GEO
  geo_count <- "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE232408&format=file&file=GSE232408_raw_counts_GRCh38.p13_NCBI.tsv.gz"
  annotation <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz"

  download.file(geo_count, "./raw-data/GSE23240_count.tsv.gz", mode = "wb")
  download.file(annotation, "./raw-data/GSE23240_annotation.tsv.gz", mode = "wb")

  dat <- read.table("./raw-data/GSE23240_count.tsv.gz", header = TRUE)
  annotation <- read.delim2("./raw-data/GSE23240_annotation.tsv.gz",
                            header = TRUE)

  # Compile count data
  GSE23240_countdata <- dat |>
    dplyr::inner_join(select(annotation, GeneID, Symbol)) |>
    dplyr::select(geneid = Symbol, tidyr::starts_with("GSM"))


  # Download and compile meta data

  gset <- getGEO("GSE232408",
                 destdir = "./raw-data",
                 getGPL = FALSE)

  # Create metadata
  GSE23240_metadata <- data.frame(participant = gset$GSE232408_series_matrix.txt.gz$title,
                         tissue = gset$GSE232408_series_matrix.txt.gz$`tissue:ch1`,
                         gender = gset$GSE232408_series_matrix.txt.gz$`gender:ch1`,
                         dose = gset$GSE232408_series_matrix.txt.gz$`dose:ch1`,
                         age = gset$GSE232408_series_matrix.txt.gz$`age (yrs):ch1`,
                         study_site = gset$GSE232408_series_matrix.txt.gz$`study site:ch1`,
                         seq_sample_id = gset$GSE232408_series_matrix.txt.gz$geo_accession) |>
    separate(participant, into = c("participant", "experiment")) |>
    separate(experiment, into = c("period", "time"), sep = "D") |>
    mutate(period = gsub("V", "", period))


  return(list(metadata = GSE23240_metadata,
              countdata = GSE23240_countdata))



}
