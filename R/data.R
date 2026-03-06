#' Pillon et al 2022 counts
#'
#' Raw count data from Pillon et al. 2022 was downloaded from NCBI Gene
#' Expression Omnibus (GEO) using the GEOquery package and prepared for use in
#' the seqwrap package. Preparation included adding gene symbols to the first
#' column of the count data frame (counts), and sorting relevant variables to
#' metadata data frame (metadata).
#'
#' @format A list with two data frames:
#' \describe{
#'   \item{counts}{Raw gene expression count matrix (genes x samples).
#'     \describe{
#'       \item{genesymbol}{HGNC gene symbol}
#'       \item{GSM...}{Raw counts, one column for each sample containing integer
#'       read counts, named by GEO accession.}
#'     }
#'   }
#'   \item{metadata}{Sample-level metadata.
#'     \describe{
#'       \item{seq_sample_id}{GEO sample accession (GSM identifier)}
#'       \item{group}{Disease group. Type 2 Diabetes, T2D;
#'       NGT, normal glucose tolerance}
#'       \item{id}{Participant id}
#'       \item{time}{Sampling time-point. Before exercise, basal; after
#'       exercise, post; after recovery (3 h), rec.}
#'     }
#'   }
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202295}
#' @references Pillon NJ, Smith JAB, Alm PS, Chibalin AV et al.
#' Distinctive exercise-induced inflammatory response and exerkine
#'  induction in skeletal muscle of people with type 2 diabetes.
#' \emph{Sci Adv} 2022 Sep 9;8(36):eabo3192. PMID: 36070371
"pillon_counts"
