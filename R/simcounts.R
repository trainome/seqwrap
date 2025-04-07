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









