#' Dungan et al 2022 counts
#'
#' Raw count data from Dungan et al. 2022 was downloaded from NCBI Gene
#' Expression Omnibus (GEO) and prepared for use in
#' the seqwrap package. Preparation included adding gene symbols to the first
#' column of the count data frame (counts), and sorting relevant variables to
#' metadata data frame (metadata).
#'
#' @format A list with two data frames:
#' \describe{
#'   \item{counts}{Raw gene expression count matrix (genes x samples).
#'     \describe{
#'       \item{genesymbol}{HGNC gene symbol}
#'       \item{OS.. OV..}{Raw counts, one column for each sample containing integer
#'       read counts, named by experimental id.}
#'     }
#'   }
#'   \item{metadata}{Sample-level metadata.
#'     \describe{
#'       \item{seq_sample_id}{Experimental unit id (mouse identifier)}
#'       \item{treatment}{Senolytic or control (Vehicle) treatment}
#'       \item{surgery}{Surgery inducing overload (synergist ablation) or sham surgury}
#'     }
#'   }
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195707}
#' @references Dungan, C. M. et al. Senolytic treatment rescues blunted
#' muscle hypertrophy in old mice. \emph{GeroScience} 44, 1925–1940 (2022).
"dungan_counts"
