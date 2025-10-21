#' Simulate counts from a parallel groups design with three time-points
#'
#'
#' Simulate gene counts from a Negative Binomial distribution conditional
#' on a two-condition, repeated measures experiment with three time points. The
#' function is used for internal testing and tutorials.
#'
#' @details
#' The function simulation of counts from a conditional Negative Binomial
#' distribution where the \deqn{\mu = \text{exp}(\eta)}, and
#' \deqn{\eta = \beta_0 + \beta1 \text{conditionB} + \beta_2 \text{timet2} +
#' \beta_3 \text{timet3} + \beta_4 \text{conditionB:timet2} +\beta_5
#' \text{conditionB:timet3} + \text{offset}(\text{library size}) +
#' \text{b}_0 + \text{b}_1 + \text{b}_2}.
#' Varying effects are added as
#' \deqn{\text{b}_l \sim \text{Normal}(0, \sigma_l)}. The library size is
#' used as an offset for simulating data. Library sizes are simulated from a
#' log-normal distribution or provided in library_size. In the simulation of
#' counts, the library size is entered as offset after scaling to the median
#' library size. Raw library sizes are included in the meta data.
#' The dispersion parameter can be simulated from a model (lm or loess)
#' provided to the function where the predictor variable should be
#' \deqn{\text{log}(\mu)} and
#' the outcome variable (dispersion) should be \deqn{\text{log}(\phi)}.
#' If no model is provided, values are simulated from
#' a log-normal distribution based on hard coded parameter values
#' from an observed mean-dispersion relationship.
#'
#'
#
#' @param n1 Number clusters in group 1 (control)
#' @param n2 Number of clusters in group 2 (treatment)
#' @param beta0 The intercept parameter on the log scale.
#' @param conditionB Baseline differences between conditions
#' @param timet2 Changes from baseline to time-point 2 in the reference group
#' @param timet3 Changes from baseline to time-point 3 in the reference group
#' @param conditionB_timet2 Difference from reference group in changes from
#' baseline to time-point 2
#' @param conditionB_timet3 Difference from reference group in changes from
#' baseline to time-point 3
#' @param b0 Vector of SD of the between cluster variation in intercept
#' @param b1 Vector of SD of the between cluster variation in timet1 effects
#' @param b2 Vector of SD of the between cluster variation in timet2 effects
#' @param phi_model A model (lm or loess) of a dispersion ~ log_mu relationship
#' @param library_size A vector of library sizes with the length (n1 + n2) * 3
#' or NULL if library sizes are to be simulated
#' @param lib_size_mean Mean of the distribution of library sizes.
#' @param lib_size_cv Coefficient of variation for the distribution of library
#' sizes.
#' @param seed Random seed.
#' @importFrom stats rnorm runif rnbinom
#' @return A list of (1) simulated counts, (2) the eta, and (3) phi values used
#' for simulating data, and (4) the meta data data frame.
#'
#'
#' @export
simcounts2 <- function(n1 = 5,
                       n2 = 5,
                       beta0 = c(2.3, 3),
                       conditionB = c(0.2, 0.1),
                       timet2 = c(0, 0),
                       timet3 = c(0.5, 0.25),
                       conditionB_timet2 = c(0.1, 0.2) ,
                       conditionB_timet3 = c(0.5, 0.6),
                       b0 = c(1, 1),
                       b1 = c(0.5, 0.5),
                       b2 = c(0.2, 0.2),
                       phi_model = NULL,
                       library_size = NULL,
                       lib_size_mean = 1e6,
                       lib_size_cv = 0.3,
                       seed = 1
                       ) {

  # Generate design data

  len_beta0 <- length(beta0)
  len_b1 <- length(conditionB)
  len_b2 <- length(timet2)
  len_b3 <- length(timet3)
  len_b4 <- length(conditionB_timet2)
  len_b5 <- length(conditionB_timet3)

  # Check that all are similar length
 if(!all(sapply(list(len_b1, len_b2, len_b3, len_b4, len_b5),
                FUN = identical, len_beta0))) {
   stop("All parameter values (conditionB, timet2, ...)
        must be of the same length.")
 }

  # Generate library sizes
  if(is.null(library_size)) {
    n_samples <- (n1 + n2) * 3  # total samples

    # Simulate realistic library sizes (log-normal distribution)
    sigma_log <- sqrt(log(1 + lib_size_cv^2))
    mu_log <- log(lib_size_mean) - sigma_log^2/2

    library_size <- round(exp(rnorm(n_samples, mu_log, sigma_log)))
  }





  # A function to predict phi from model parameters
  # or use hard-coded values from approximation
  phi_predict <- function(phi_model, log_mu) {

    if(is.null(phi_model)) {

      phi <- 0.9273216 +
        (-0.3868206) * log_mu +
        0.3816686 * log_mu^2 +
        (-0.0549217) * log_mu^3 +
        0.0021181 * log_mu^4

    } else {

      phi <- predict(phi_model, newdata = data.frame(log_mu = log_mu))
      phi # This is on the log scale (noted in details)
    }

    return(phi)

  }

  # A function to extract the sigma from a phi-prediction model
  phi_sigma <- function(phi_model) {

    if(is.null(phi_model)) {
      phi_sigma <- 4.25
    }

    if(inherits(phi_model, "loess")) {
      phi_sigma <- phi_model$s
    }

    if(inherits(phi_model, "lm")) {
      phi_sigma <- sigma(phi_model)
    }

    return(phi_sigma)
  }





  # The coefficient matrix contain all parameter values collected
  # row-wise.
  coef_mat <- matrix(c(beta0,
                       conditionB,
                       timet2,
                       timet3,
                       conditionB_timet2,
                       conditionB_timet3),
                     ncol = 6)
  colnames(coef_mat) <- c("Intercept", "conditionB", "timet2",
                          "timet3", "conditionB:timet2", "conditionB:timet3")


# Create the predictor data frame

 design <-  rbind(

   expand.grid(id = paste0("A", 1:n1),
               time = c("t1", "t2", "t3"),
               condition = c("A")),

    expand.grid(id = paste0("B", 1:n2),
              time = c("t1", "t2", "t3"),
              condition = c("B"))

    )

 # Adding library sizes to the design
 design$library_size <- library_size

 # Create a model prediction matrix
 mod_mat <- model.matrix(~ condition * time, data = design)

 # Calculate mean responses. The model matrix is multiplied with
 # the parameter values for each gene.
 # The parameter values are determined by the linear predictor eta and
 # includes the varying effects (eta_full).

 # From the loop below we will store the phi, eta_full, parameter values
 # and simulated counts

 # Counts
 counts_list <- list()
 # eta
 eta_list <- list()
 # phi
 phi_list <- list()

 for(i in seq_along(beta0)) {

  # Combine fixed effects to get the eta
  eta_fixed <-  mod_mat %*% coef_mat[i,]

  # Add library size offset
  offset <- log(library_size / median(library_size))


  # Adding varying effects (non-correlated)
  b_0 <- c(rep(rnorm(n1, 0, b0[i]), 3), rep(rnorm(n2, 0, b0[i]), 3))
  b_1 <- c(rep(0, n1), rnorm(n1, 0, b1[i]), rep(0, n1),
           rep(0, n2), rnorm(n2, 0, b1[i]), rep(0, n2))
  b_2 <- c(rep(0, n1), rep(0, n1), rnorm(n1, 0, b2[i]),
           rep(0, n2), rep(0, n2), rnorm(n2, 0, b2[i]))

  # Combining all parameters into the linear predictor eta.
  eta_full <- eta_fixed + offset + b_0 + b_1 + b_2

  # Get the phi from predictive model and combine with phi_sigma
  log_phi <- phi_predict(phi_model, log(mean(exp(eta_fixed))))
  log_phi_sd <- phi_sigma(phi_model)

  phi_g <- exp(rnorm(1, log_phi, log_phi_sd))

  # Simulate from the Negative binomial distribution
  counts <- rnbinom(n = length(eta_full), size = phi_g, mu = exp(eta_full))
  names(counts) <- paste0(design$id,"_",design$time)
  names(eta_full) <- paste0(design$id,"_",design$time)
  # Collect variables
  counts_list[[i]] <- cbind(data.frame(gene = i), t(counts))
  eta_list[[i]] <- cbind(data.frame(gene = i), t(eta_full))
  phi_list[[i]] <- cbind(data.frame(gene = i), phi_g)


 }

 # Combine all list and return values
 counts <- do.call(rbind, counts_list)
 eta <- do.call(rbind, eta_list)
 phi <- do.call(rbind, phi_list)
 design$seq_sample_id <- paste0(design$id,"_",design$time)

 return(list(counts = counts,
             eta = eta,
             phi = phi,
             metadata = design))



}

