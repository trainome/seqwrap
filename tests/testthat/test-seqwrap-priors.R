library(testthat)
library(seqwrap)

# A first run without priors, shared by the tests below. Target identifiers
# are shuffled so that data row order differs from sorted order, which is what
# a per-target prior list has to survive.
prior_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) {
      return(cache)
    }
    skip_if_not_installed("glmmTMB")

    dat <- seqwrap::simcounts(
      seed = 11,
      n_genes = 40,
      n_samples = 30,
      clusters = 15
    )
    counts <- dat$data
    set.seed(3)
    counts <- counts[sample(nrow(counts)), ]
    rownames(counts) <- NULL
    metadata <- dat$metadata
    metadata$ln_libsize <- log(colSums(counts[, -1]))

    container <- seqwrap::seqwrap_compose(
      modelfun = glmmTMB::glmmTMB,
      arguments = list(
        formula = y ~ x + (1 | cluster) + offset(ln_libsize),
        family = glmmTMB::nbinom2
      ),
      data = counts,
      metadata = metadata,
      samplename = "sample",
      eval_fun = seqwrap::dispersion_evaluation
    )

    results <- seqwrap::seqwrap(container, cores = 1, verbose = FALSE)

    cache <<- list(counts = counts, metadata = metadata, results = results)
    cache
  }
})


test_that("dispersion_evaluation reads the log dispersion and mean count", {
  skip_if_not_installed("glmmTMB")

  dat <- seqwrap::simcounts(seed = 5, n_genes = 3, n_samples = 30)
  metadata <- dat$metadata
  metadata$y <- as.integer(dat$data[1, -1])

  m <- glmmTMB::glmmTMB(y ~ x, data = metadata, family = glmmTMB::nbinom2)

  out <- dispersion_evaluation(m)

  expect_s3_class(out, "data.frame")
  expect_identical(nrow(out), 1L)
  expect_named(out, c("dispersion", "dispersion.se", "log_mu"))
  expect_equal(out$dispersion, log(stats::sigma(m)), tolerance = 1e-6)
  expect_true(is.finite(out$dispersion.se) && out$dispersion.se > 0)
  expect_equal(out$log_mu, log(mean(metadata$y)))

  expect_null(dispersion_evaluation(NULL))
  expect_error(
    dispersion_evaluation(stats::lm(mpg ~ wt, data = mtcars)),
    "glmmTMB"
  )
})


test_that("seqwrap_priors returns one prior data frame per data row", {
  fx <- prior_fixture()

  priors <- seqwrap_priors(fx$results, data = fx$counts)

  expect_s3_class(priors, "seqwrap_priors")
  expect_true(is.list(priors))
  expect_length(priors, nrow(fx$counts))
  # Same order as the data, not sorted order
  expect_identical(names(priors), as.character(fx$counts[[1]]))

  one <- priors[[1]]
  expect_s3_class(one, "data.frame")
  expect_named(one, c("prior", "class", "coef"))
  expect_true(all(vapply(one, is.character, logical(1))))

  # A fixed effect on x, a random effect SD on cluster and a dispersion prior
  expect_identical(one$class, c("fixef", "ranef_sd", "fixef_disp"))
  expect_identical(one$coef, c("x", "cluster", "1"))
  expect_match(one$prior[1], "^normal\\(0, [0-9.]+\\)$")
  expect_match(one$prior[2], "^gamma\\([0-9.]+, 2\\)$")
  expect_match(one$prior[3], "^normal\\(-?[0-9.]+, [0-9.]+\\)$")

  # Shared priors are the same for every target, the dispersion prior is not
  shared <- unique(lapply(priors, function(p) p$prior[1:2]))
  expect_length(shared, 1L)
  disp <- vapply(priors, function(p) p$prior[3], character(1))
  expect_gt(length(unique(disp)), 1L)

  # Attributes describing the construction
  common <- attr(priors, "common")
  expect_s3_class(common, "data.frame")
  expect_named(common, c("prior", "class", "coef", "location", "scale", "n"))
  expect_identical(common$n, rep(nrow(fx$counts), 2L))

  disp_table <- attr(priors, "dispersion")
  expect_identical(disp_table$target, names(priors))
  expect_equal(
    disp_table$log_mu,
    log(rowMeans(fx$counts[, -1]))
  )
  expect_s3_class(attr(priors, "trend"), "loess")
  expect_identical(attr(priors, "n"), nrow(fx$counts))

  # The estimates behind each prior are kept for plotting
  pd <- attr(priors, "data")
  expect_named(pd, c("estimates", "dispersion"))
  est <- pd$estimates
  expect_named(est, c("target", "class", "coef", "estimate"))
  expect_identical(sort(unique(est$class)), c("fixef", "ranef_sd"))
  fixef_x <- est[est$class == "fixef" & est$coef == "x", ]
  expect_identical(nrow(fixef_x), nrow(fx$counts))
  expect_setequal(fixef_x$target, names(priors))
  # The prior scale is the SD of those estimates
  expect_equal(common$scale[common$coef == "x"], sd(fixef_x$estimate))
  ranef_cl <- est[est$class == "ranef_sd", ]
  expect_identical(unique(ranef_cl$coef), "cluster")
  expect_equal(common$location[common$coef == "cluster"], mean(ranef_cl$estimate))
  obs <- pd$dispersion
  expect_named(obs, c("target", "dispersion", "log_mu"))
  expect_identical(nrow(obs), nrow(fx$counts))
  expect_equal(
    obs$log_mu,
    disp_table$log_mu[match(obs$target, disp_table$target)]
  )

  # Printing dispatches to the seqwrap_priors method and returns the object
  expect_message(out <- print(priors), "seqwrap priors")
  expect_identical(out, priors)
})


test_that("seqwrap_priors options change the priors as documented", {
  fx <- prior_fixture()
  summaries <- seqwrap_summarise(fx$results, verbose = FALSE)

  # Centered on the mean across targets
  centered <- seqwrap_priors(fx$results, data = fx$counts, center = "mean")
  common <- attr(centered, "common")
  expect_equal(
    common$location[common$coef == "x"],
    mean(summaries$summaries$estimate[summaries$summaries$term == "x"])
  )
  expect_false(grepl("^normal\\(0,", centered[[1]]$prior[1]))

  # Naming terms, including the intercept
  with_int <- seqwrap_priors(
    fx$results,
    data = fx$counts,
    terms = c("(Intercept)", "x"),
    center = "mean"
  )
  expect_identical(with_int[[1]]$coef, c("(Intercept)", "x", "cluster", "1"))
  expect_error(
    seqwrap_priors(fx$results, data = fx$counts, terms = "nope"),
    "not fixed-effect terms"
  )

  # No random-effect or dispersion priors
  bare <- seqwrap_priors(
    fx$results,
    data = fx$counts,
    ranef = FALSE,
    dispersion = NULL
  )
  expect_identical(bare[[1]]$class, "fixef")
  expect_null(attr(bare, "dispersion"))
  expect_null(attr(bare, "trend"))

  # A constant dispersion prior is shared by every target
  constant <- seqwrap_priors(fx$results, data = fx$counts, trend = "constant")
  disp <- vapply(constant, function(p) p$prior[3], character(1))
  expect_length(unique(disp), 1L)
  expect_null(attr(constant, "trend"))

  # Gamma shape
  shaped <- seqwrap_priors(fx$results, data = fx$counts, shape = 4)
  expect_match(shaped[[1]]$prior[2], ", 4\\)$")

  # Robust summaries use the median and MAD. Random-effect standard
  # deviations are often at the boundary, in which case the prior is skipped
  # with a warning, so that warning is tolerated here.
  robust <- suppressWarnings(
    seqwrap_priors(fx$results, data = fx$counts, robust = TRUE)
  )
  common_r <- attr(robust, "common")
  expect_equal(
    common_r$scale[common_r$coef == "x"],
    stats::mad(summaries$summaries$estimate[summaries$summaries$term == "x"])
  )

  # The combined list from seqwrap_summarise() is accepted in place of results
  from_list <- seqwrap_priors(summaries, data = fx$counts)
  expect_identical(
    attr(from_list, "common")$prior,
    attr(seqwrap_priors(fx$results, data = fx$counts), "common")$prior
  )
})


test_that("seqwrap_priors reports missing inputs clearly", {
  fx <- prior_fixture()

  expect_error(seqwrap_priors(list(a = 1), data = fx$counts), "seqwrapResults")
  expect_error(
    seqwrap_priors(fx$results, data = list(fx$counts)),
    "single data frame"
  )
  expect_error(
    seqwrap_priors(fx$results, data = fx$counts, shape = -1),
    "positive"
  )

  # No dispersion column when the default evaluation was used
  container <- seqwrap::seqwrap_compose(
    modelfun = glmmTMB::glmmTMB,
    arguments = list(
      formula = y ~ x + (1 | cluster) + offset(ln_libsize),
      family = glmmTMB::nbinom2
    ),
    data = fx$counts,
    metadata = fx$metadata,
    samplename = "sample"
  )
  plain <- seqwrap::seqwrap(container, subset = 1:12, cores = 1,
                            verbose = FALSE)

  expect_error(
    seqwrap_priors(plain, data = fx$counts),
    "dispersion_evaluation"
  )
  expect_no_error(
    seqwrap_priors(plain, data = fx$counts, dispersion = NULL)
  )
})


test_that("seqwrap_priors output fits as targetdata and reaches each target", {
  fx <- prior_fixture()

  priors <- seqwrap_priors(fx$results, data = fx$counts)

  container <- seqwrap::seqwrap_compose(
    modelfun = glmmTMB::glmmTMB,
    arguments = alist(
      formula = y ~ x + (1 | cluster) + offset(ln_libsize),
      family = glmmTMB::nbinom2,
      priors = data.frame(prior = prior, class = class, coef = coef)
    ),
    data = fx$counts,
    metadata = fx$metadata,
    samplename = "sample",
    targetdata = priors
  )

  # Named by target, so matched by name without a positional warning
  expect_no_warning(results <- seqwrap::seqwrap(
    container,
    subset = 1:6,
    return_models = TRUE,
    cores = 1,
    verbose = FALSE
  ))

  expect_identical(nrow(seqwrap_errors(results, type = "error")), 0L)

  # Every model received the prior built for its own target, even though the
  # data rows are not in sorted order
  for (nm in names(results@models)) {
    used <- results@models[[nm]]$modelInfo$priors
    expect_identical(as.character(used$prior), priors[[nm]]$prior)
    expect_identical(as.character(used$coef), priors[[nm]]$coef)
  }
})


test_that("seqwrap_priors works from custom summaries and skips degenerate SDs", {
  counts <- data.frame(gene = paste0("g", 1:6), s1 = 1:6, s2 = 2:7)
  set.seed(4)
  sm <- data.frame(
    target = rep(counts$gene, each = 3),
    term = rep(c("(Intercept)", "x", "sd__(Intercept)"), 6),
    estimate = as.vector(rbind(rnorm(6, 5), rnorm(6, 0.5, 0.2), rep(1e-5, 6))),
    stringsAsFactors = FALSE
  )

  # Without effect and group columns terms are classified by name, and no
  # random-effect prior can be addressed to a grouping variable
  expect_message(
    p <- seqwrap_priors(list(summaries = sm), data = counts,
                        dispersion = NULL),
    "group"
  )
  expect_identical(p[[1]]$coef, "x")
  expect_identical(names(p), counts$gene)

  # Standard deviations at the boundary give no usable gamma prior
  sm$effect <- ifelse(grepl("^sd__", sm$term), "ran_pars", "fixed")
  sm$group <- ifelse(sm$effect == "ran_pars", "id", NA_character_)
  expect_warning(
    p2 <- seqwrap_priors(list(summaries = sm), data = counts,
                         dispersion = NULL),
    "essentially zero"
  )
  expect_identical(p2[[1]]$class, "fixef")

  # A positive mean gives the gamma prior on the grouping variable
  sm$estimate[sm$term == "sd__(Intercept)"] <- 0.5
  p3 <- seqwrap_priors(list(summaries = sm), data = counts, dispersion = NULL)
  expect_identical(p3[[1]]$prior[2], "gamma(0.5, 2)")
  expect_identical(p3[[1]]$coef[2], "id")

  # Centering on the mean uses the average estimate
  p4 <- seqwrap_priors(list(summaries = sm), data = counts,
                       dispersion = NULL, center = "mean", digits = 2)
  expect_identical(
    p4[[1]]$prior[1],
    sprintf("normal(%s, %s)",
            round(mean(sm$estimate[sm$term == "x"]), 2),
            round(sd(sm$estimate[sm$term == "x"]), 2))
  )
})


test_that("seqwrap_priors can be plotted from the retained estimates", {
  fx <- prior_fixture()
  priors <- seqwrap_priors(fx$results, data = fx$counts)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  # Draws without conditions and returns the plotting data invisibly
  expect_no_condition(out <- plot(priors))
  expect_identical(out, attr(priors, "data"))
  expect_no_condition(plot(priors, which = "dispersion", layout = FALSE))
  expect_no_condition(plot(priors, which = c("fixef", "ranef"), trim = 1))

  # Constant trend draws a flat prior
  flat <- seqwrap_priors(fx$results, data = fx$counts, trend = "constant")
  expect_no_condition(plot(flat, which = "dispersion"))

  # Asking for panels the object does not hold is reported, not an error
  no_ranef <- seqwrap_priors(fx$results, data = fx$counts, ranef = FALSE)
  expect_false("ranef_sd" %in% attr(no_ranef, "data")$estimates$class)
  expect_message(plot(no_ranef, which = "ranef"), "Nothing to plot")
  no_disp <- seqwrap_priors(fx$results, data = fx$counts, dispersion = NULL)
  expect_null(attr(no_disp, "data")$dispersion)
  expect_message(plot(no_disp, which = "dispersion"), "Nothing to plot")

  # plot = TRUE draws while building
  expect_no_condition(
    built <- seqwrap_priors(fx$results, data = fx$counts, plot = TRUE)
  )
  expect_identical(built[[1]], priors[[1]])

  # Argument checks
  expect_error(plot(priors, which = "nope"), "should be one of")
  expect_error(plot(priors, trim = 0), "trim")
  expect_error(plot(priors, trim = 1.5), "trim")
  expect_error(
    seqwrap_priors(fx$results, data = fx$counts, plot = NA),
    "plot"
  )
  expect_error(plot(structure(list(), class = "seqwrap_priors")), "data")
})
