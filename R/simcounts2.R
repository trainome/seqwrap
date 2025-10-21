#' Simulate counts from a parallell groups design with three time-points
#' and two groups.


# Create example data

#' @param n1 Number clusters in group 1 (control)
#' @param n2 Number of clusters in group 2 (treatment)
#' @param b0 Vector sd of the
#' @param b1
#' @param b2



simcounts2 <- function(n1 = 5,
                       n2 = 5,
                       beta0 = c(2.3, 3),
                       conditioncon = c(0.2, 0.1),
                       timetime2 = c(0, 0),
                       timetime3 = c(0.5, 0.25),
                       conditioncon_timetime2 = c(0.1, 0.2) ,
                       conditioncon_timetime3 = c(0.5, 0.6)
                       ) {

  # Generate design data

  len_beta0 <- length(beta0)
  len_b1 <- length(conditionA)
  len_b2 <- length(timetime2)
  len_b3 <- length(timetime3)
  len_b4 <- length(conditioncon_timetime2)
  len_b5 <- length(conditioncon_timetime3)

  # Check that all are similar length
 if(!all(sapply(list(len_b1, len_b2, len_b3, len_b4, len_b5),
                FUN = identical, len_beta0))) {
   stop("All parameter values (conditioncon, timetime2, ...)
        must be of the same length.")
 }

  coef_mat <- matrix(c(beta0,
                       conditioncon,
                       timetime2,
                       timetime3,
                       conditioncon_timetime2,
                       conditioncon_timetime2),
                     ncol = 6)


# Treatment

 design <-  rbind(
    expand.grid(id = paste0("B", 1:n2),
              time = c("t1", "t2", "t3"),
              condition = c("B")),

    expand.grid(id = paste0("A", 1:n1),
                time = c("t1", "t2", "t3"),
                condition = c("A"))
  )


 mod_mat <- model.matrix(~ condition * time, data = design)


 for(i in seq_along(len_b1)) {


   mod_mat %*% coef_mat[2,]

 }






}




mod <- simcounts2()

#' Simulate counts from a simple experiment (paired or unpaierd) using
#' Poisson or negative binomial distributions.
#'
#'
#' @param n_genes Number of genes to simulate.
#' @param n_samples Number of samples to simulate.
#' @param beta_0 Gene intercepts.
#' @param sigma_0 SD of gene intercepts
#' @param beta_1 Gene slopes.
#' @param sigma_1 SD of gene slopes
#' @param b_0 Cluster-specific intercept. If NULL, non-paired data is simulated.
#' @param clusters Number of clusters. If NULL, non-paired data is simulated.
#' @param sample_sd SD of sample-specific effects.
#' @param overdispersion_min_max Range of overdispersion parameter. If NULL,
#' the Poisson distrubution is used to simulate data.
#' @param seed Random seed for reproducibility.
#' @importFrom stats rnorm runif rpois rnbinom
#'
#'
#' @keywords internal
simcounts <- function(n_genes = 10,
                      n_samples = 16,
                      beta_0 = 1,
                      sigma_0 = 0.5,
                      beta_1 = 1,
                      sigma_1 = 0.5,
                      b_0 = 0.5,
                      clusters = 8,
                      sample_sd = 0.5,
                      overdispersion_min_max = c(1,10),
                      seed = 123) {


  set.seed(seed)
  # Draw gene-specific effects
  beta_0 <- rnorm(n_genes, mean = beta_0, sd = sigma_0)
  beta_1 <- rnorm(n_genes, mean = beta_1, sd = sigma_1)
  # Draw sample-specific effects
  sample_effects <- rnorm(n_samples, mean = 0, sd = sample_sd)
  # Draw cluster-specific effects
  if (!is.null(b_0)) {
    if (is.null(clusters)) {
      stop("Clusters must be provided when b_0 is not NULL.")
    }
    cluster_effects <- rnorm(clusters, mean = 0, sd = b_0)
    cluster_effects <- rep(cluster_effects, n_samples/clusters)
  } else {
    cluster_effects <- rep(0, n_samples)
  }

  # Expression values
  expr_mat <- matrix(0, nrow = n_genes, ncol = n_samples)

  # Design matrix
  x <- rep(c(0,1), each = n_samples/2)
  if (!is.null(b_0)) {
    cluster <- rep(1:clusters, length = n_samples)
  }

  # If overdispersion calculate gene-wise overdispersion parameter
  # Draw overdispersion parameter from uniform distribution
  if (!is.null(overdispersion_min_max)) {
    overdispersion <- runif(n_genes, min = overdispersion_min_max[1],
                            max = overdispersion_min_max[2])
  }






  metadata <- data.frame(
    sample = paste0("samp", seq_len(n_samples)),
    cluster = paste0("c",cluster),
    x = x
  )


    for(gene in seq_len(n_genes)) {
      for (sample in seq_len(n_samples)) {
        lambda <- exp(beta_0[gene] +
                        beta_1[gene] * x[sample] +
                        cluster_effects[cluster[sample]] +
                        sample_effects[sample])

        if (is.null(overdispersion_min_max)) {
          # Poisson distribution
          expr_mat[gene, sample] <- rpois(1, lambda)
        } else {
          # Negative binomial distribution

          expr_mat[gene, sample] <- rnbinom(1,
                                            size = overdispersion[gene],
                                            mu = lambda)
        }



      }
    }




  data <- data.frame(targetid = paste0("gene", seq_len(n_genes)),
                     expr_mat)

  colnames(data)[-1] <- metadata$sample


  parameters <- data.frame(
    targetid = paste0("gene", seq_len(n_genes)),
    beta_0 = beta_0,
    beta_1 = beta_1,
    overdispersion = overdispersion
  )


  # Return list with true parameters,  expression matrix and metadata
  return(list(data = data,
              parameters = parameters,
              metadata = metadata))

}




#' Simulate counts from a simple experiment (paired or unpaired) using
#' Poisson or negative binomial distributions.
#'
#'
#' @param n_genes Number of genes to simulate.
#' @param n_samples Number of samples to simulate.
#' @param beta_0 Gene intercepts.
#' @param sigma_0 SD of gene intercepts
#' @param beta_1 Gene slopes.
#' @param sigma_1 SD of gene slopes
#' @param b_0 Cluster-specific intercept. If NULL, non-paired data is simulated.
#' @param clusters Number of clusters. If NULL, non-paired data is simulated.
#' @param sample_sd SD of sample-specific effects.
#' @param overdispersion_min_max Range of overdispersion parameter. If NULL,
#' the Poisson distrubution is used to simulate data.
#' @param seed Random seed for reproducibility.
#' @importFrom stats rnorm runif rpois rnbinom
#'
#'
#' @keywords internal
simcounts <- function(n_genes = 10,
                      n_samples = 16,
                      beta_0 = 1,
                      sigma_0 = 0.5,
                      beta_1 = 1,
                      sigma_1 = 0.5,
                      b_0 = 0.5,
                      clusters = 8,
                      sample_sd = 0.5,
                      overdispersion_min_max = c(1,10),
                      seed = 123) {


  set.seed(seed)
  # Draw gene-specific effects
  beta_0 <- rnorm(n_genes, mean = beta_0, sd = sigma_0)
  beta_1 <- rnorm(n_genes, mean = beta_1, sd = sigma_1)
  # Draw sample-specific effects
  sample_effects <- rnorm(n_samples, mean = 0, sd = sample_sd)
  # Draw cluster-specific effects
  if (!is.null(b_0)) {
    if (is.null(clusters)) {
      stop("Clusters must be provided when b_0 is not NULL.")
    }
    cluster_effects <- rnorm(clusters, mean = 0, sd = b_0)
    cluster_effects <- rep(cluster_effects, n_samples/clusters)
  } else {
    cluster_effects <- rep(0, n_samples)
  }

  # Expression values
  expr_mat <- matrix(0, nrow = n_genes, ncol = n_samples)

  # Design matrix
  x <- rep(c(0,1), each = n_samples/2)
  if (!is.null(b_0)) {
    cluster <- rep(1:clusters, length = n_samples)
  }

  # If overdispersion calculate gene-wise overdispersion parameter
  # Draw overdispersion parameter from uniform distribution
  if (!is.null(overdispersion_min_max)) {
    overdispersion <- runif(n_genes, min = overdispersion_min_max[1],
                            max = overdispersion_min_max[2])
  }






  metadata <- data.frame(
    sample = paste0("samp", seq_len(n_samples)),
    cluster = paste0("c",cluster),
    x = x
  )


  for(gene in seq_len(n_genes)) {
    for (sample in seq_len(n_samples)) {
      lambda <- exp(beta_0[gene] +
                      beta_1[gene] * x[sample] +
                      cluster_effects[cluster[sample]] +
                      sample_effects[sample])

      if (is.null(overdispersion_min_max)) {
        # Poisson distribution
        expr_mat[gene, sample] <- rpois(1, lambda)
      } else {
        # Negative binomial distribution

        expr_mat[gene, sample] <- rnbinom(1,
                                          size = overdispersion[gene],
                                          mu = lambda)
      }



    }
  }




  data <- data.frame(targetid = paste0("gene", seq_len(n_genes)),
                     expr_mat)

  colnames(data)[-1] <- metadata$sample


  parameters <- data.frame(
    targetid = paste0("gene", seq_len(n_genes)),
    beta_0 = beta_0,
    beta_1 = beta_1,
    overdispersion = overdispersion
  )


  # Return list with true parameters,  expression matrix and metadata
  return(list(data = data,
              parameters = parameters,
              metadata = metadata))

}













