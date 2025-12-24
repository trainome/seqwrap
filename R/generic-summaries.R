#' Generic summaries for model parameter estimates.
#'
#' When no summary function is provided to seqwrap (summary_fun), this function
#' uses broom(.mixed)::tidy to give a table of model parameter estimates.
#'
#' @param x A model fitted in seqwrap.
#'
#' @returns A tidy data frame possible to bind using seqwrap_summarise
#' @export
generic_summary <- function(x){
 out <- broom.mixed::tidy(x)
 return(out)
}
#' Generic evaluation for model fits
#'
#' When no evaluation function is provided to seqwrap (eval_fun), this function
#' uses DHARMa to simulate residuals and test for uniformity, dispersion and
#' outliers. These tests could guide troubleshooting and detection of model
#' misspecification.
#'
#' @param x A model fitted in seqwrap.
#'
#' @returns A tidy data frame possible to combine using seqwrap_summarise
#' @export
generic_evaluation <- function(x){
 # Simulate residuals using DHARMa
 sim_resid <- DHARMa::simulateResiduals(x, n = 250, plot = FALSE)
 # Combine potential metrics
 out <- tibble::tibble(
   uniformity = DHARMa::testUniformity(sim_resid, plot = FALSE)$p.value,
   dispersion = DHARMa::testDispersion(sim_resid, plot = FALSE)$p.value,
   outliers = DHARMa::testOutliers(sim_resid, plot = FALSE)$p.value)
 return(out)
}
