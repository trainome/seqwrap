## Empirical Bayes priors for glmmTMB --------------------------------------
##
## glmmTMB accepts a data frame of priors, which makes it possible to
## regularise item-by-item models with information pooled across targets: fit
## every target once without priors, summarise the distribution of each
## parameter over targets, and feed those summaries back as priors in a second
## run. The functions here build that prior data frame, one per target, in the
## shape seqwrap_compose() expects for `targetdata`.


#' Dispersion and mean count from a glmmTMB model
#'
#' An evaluation function for `eval_fun` in [seqwrap()] and
#' [seqwrap_compose()] that records the estimated dispersion of a negative
#' binomial `glmmTMB` model together with the mean of the observed counts. It
#' provides the per-target quantities that [seqwrap_priors()] uses to build a
#' dispersion prior that depends on the expression level of each target.
#'
#' @param x A model fitted with `glmmTMB::glmmTMB`, or NULL when fitting
#'   failed.
#' @param ... Currently unused.
#'
#' @return A one row tibble, or NULL when `x` is NULL, with columns:
#' \describe{
#'   \item{dispersion}{The dispersion parameter on the log scale, as estimated
#'     by `glmmTMB`. For `nbinom2` this is `log(theta)`, which equals
#'     `log(sigma(x))`.}
#'   \item{dispersion.se}{The standard error of `dispersion`, NA when
#'     `glmmTMB` did not compute one.}
#'   \item{log_mu}{The natural logarithm of the mean of the observed response
#'     values used to fit the model.}
#' }
#'
#' @details
#' The dispersion is read from the `betadisp` element of the model's `sdreport`.
#' When a dispersion formula contains covariates only its first element, the
#' intercept, is reported. The convergence diagnostics returned by
#' [generic_evaluation()] are not included; to keep both, wrap the two in a
#' function that binds their results.
#'
#' @seealso [seqwrap_priors()], which consumes the columns returned here.
#'
#' @examples
#' library(seqwrap)
#'
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'   library(glmmTMB)
#'
#'   dat <- simcounts(n_genes = 10, n_samples = 30)
#'
#'   counts <- dat$data
#'   metadata <- dat$metadata
#'   metadata$ln_libsize <- log(colSums(counts[, -1]))
#'   metadata$y <- as.integer(counts[1, -1])
#'
#'   m <- glmmTMB(y ~ x + offset(ln_libsize),
#'                data = metadata,
#'                family = nbinom2)
#'
#'   dispersion_evaluation(m)
#'
#'   # Equivalent to the dispersion reported by sigma() on the log scale
#'   log(sigma(m))
#' }
#'
#' @export
dispersion_evaluation <- function(x, ...) {
  if (is.null(x)) {
    return(NULL)
  }

  if (!inherits(x, "glmmTMB")) {
    cli::cli_abort(c(
      "{.fn dispersion_evaluation} expects a model fitted with
       {.fn glmmTMB::glmmTMB}.",
      "x" = "Got an object of class {.cls {class(x)}}."
    ))
  }

  estimate <- NA_real_
  se <- NA_real_

  sdr <- x$sdr
  if (!is.null(sdr)) {
    par <- sdr$par.fixed
    idx <- which(names(par) == "betadisp")
    if (length(idx) > 0L) {
      estimate <- unname(par[idx[1]])
      if (!is.null(sdr$cov.fixed)) {
        se <- unname(sqrt(diag(sdr$cov.fixed))[idx[1]])
      }
    }
  }

  # Without an sdreport (se = FALSE in glmmTMBControl) fall back on sigma(),
  # which returns the dispersion on the natural scale.
  if (is.na(estimate)) {
    estimate <- tryCatch(
      log(stats::sigma(x)),
      error = function(e) NA_real_
    )
  }

  # glmmTMB records which column of the model frame holds the response.
  resp_col <- x$modelInfo$respCol
  y <- if (!is.null(resp_col)) x$frame[[resp_col]] else x$frame[["y"]]

  tibble::tibble(
    dispersion = estimate,
    dispersion.se = se,
    log_mu = log(mean(as.numeric(y), na.rm = TRUE))
  )
}


#' Build empirical Bayes priors for glmmTMB from a seqwrap run
#'
#' Summarises the parameter estimates from a completed [seqwrap()] run and
#' turns them into target-specific prior specifications for
#' `glmmTMB::glmmTMB`. The result is a list with one prior data frame per
#' target, ready to be passed as `targetdata` to [seqwrap_compose()] for a
#' second, regularised run. This is a simple empirical Bayes strategy: every
#' target is first fitted on its own, the spread of each parameter across
#' targets then becomes the prior for that parameter.
#'
#' @param x A seqwrapResults object from a run without priors, or the list
#'   returned by [seqwrap_summarise()] for that run.
#' @param data The target data frame given to [seqwrap_compose()], with target
#'   identifiers in the first column (or in the row names when
#'   `rownames = TRUE`) and one column per sample. It determines the targets
#'   for which priors are built, their order, and the mean count of each target.
#'   When [seqwrap_compose()] was given a list of data frames, pass the element
#'   holding the counts.
#' @param rownames Logical, are target identifiers held in the row names of
#'   `data`? Defaults to FALSE, as in [seqwrap_compose()].
#' @param terms Character vector naming the fixed-effect terms to place priors
#'   on. NULL (the default) uses every fixed-effect term in the summaries except
#'   the intercept. The intercept can be included by naming it, `"(Intercept)"`,
#'   which is sensible only with `center = "mean"`.
#' @param center Where to center the fixed-effect priors. `"zero"` (the
#'   default) shrinks estimates towards no effect, using `normal(0, s)` with `s`
#'   the spread of the estimates across targets. `"mean"` centers each prior on
#'   the average estimate across targets, `normal(m, s)`.
#' @param robust Logical. When TRUE, the median and median absolute deviation
#'   replace the mean and standard deviation when summarising estimates, and the
#'   dispersion trend is fitted with a robust loss. This limits the influence
#'   of targets whose fits went astray. Defaults to FALSE.
#' @param ranef Logical, should priors be placed on random-effect standard
#'   deviations? Defaults to TRUE. Ignored when the summaries contain no random
#'   effect terms. A grouping variable whose standard deviation estimates are
#'   essentially zero across targets is skipped with a warning, since the
#'   gamma prior needs a positive mean.
#' @param shape The shape of the gamma prior placed on random-effect standard
#'   deviations. `glmmTMB` parameterises this prior by its mean and shape, and
#'   the mean is taken from the estimates across targets. Defaults to 2.
#' @param dispersion The name of the column in the evaluations that holds the
#'   log dispersion of each target, as produced by [dispersion_evaluation()].
#'   Defaults to `"dispersion"`. Set to NULL to build no dispersion prior.
#' @param trend How the dispersion prior depends on expression level.
#'   `"loess"` (the default) fits a loess curve of log dispersion against log
#'   mean count and gives every target a prior centered on the curve at its own
#'   mean count, with the residual standard error of the curve as prior
#'   standard deviation. `"constant"` gives every target the same prior,
#'   centered on the average log dispersion across targets.
#' @param span The span of the loess fit when `trend = "loess"`. Defaults to
#'   0.75.
#' @param drop_warnings Passed to [seqwrap_summarise()] when `x` is a
#'   seqwrapResults object, to exclude targets that raised warnings from the
#'   prior estimation. Defaults to FALSE. Ignored when `x` is already a list of
#'   combined results.
#' @param digits The number of decimals kept in the prior specifications.
#'   Defaults to 3.
#' @param plot Logical, should the priors be plotted over the estimates they
#'   were built from before the object is returned? Defaults to FALSE. See
#'   [plot.seqwrap_priors()], which can also be called on the result later.
#'
#' @return A named list of class `seqwrap_priors` with one element per row of
#'   `data`, in the same order and named by target identifier. Each element is
#'   a data frame with the columns `prior`, `class` and `coef` that
#'   `glmmTMB::glmmTMB` expects in its `priors` argument. The list carries
#'   attributes that describe how it was built: `common`, a data frame of the
#'   priors shared by all targets with their location and scale;
#'   `dispersion`, a data frame with the log mean count and prior location of
#'   each target, or NULL; `trend`, the fitted loess model, or NULL; `n`,
#'   the number of targets that contributed to the estimation; and `data`,
#'   the per-target estimates the priors were built from, see below.
#'
#' @section Estimates used to build the priors:
#' The `data` attribute holds the estimates that went into each prior, so that
#' the priors can be inspected and plotted with other tools than
#' [plot.seqwrap_priors()]. It is a list with two elements.
#'
#' * `estimates`, a data frame with one row per target and prior, with the
#'   columns `target`, `class` (`"fixef"` or `"ranef_sd"`), `coef` and
#'   `estimate`. For fixed effects `estimate` is the coefficient of the term
#'   in the target's model. For random effects it is the estimated standard
#'   deviation, one row per standard deviation term within the grouping
#'   variable.
#' * `dispersion`, a data frame with the columns `target`, `log_mu` and
#'   `dispersion` holding the log mean count and the estimated log dispersion
#'   of every target used to estimate the dispersion trend, or NULL when no
#'   dispersion prior was built.
#'
#' The location and scale of each prior is found in the `common` and
#' `dispersion` attributes, and the loess fit in `trend`.
#'
#' @details
#' Three kinds of prior are built.
#'
#' * Fixed effects (class `"fixef"`). For each term the estimates across
#'   targets are summarised by their location and spread, giving a normal prior
#'   `normal(location, spread)`. Terms are identified from the `effect` and
#'   `component` columns of a [generic_summary()] result. With a custom summary
#'   function that lacks these columns, every term other than the intercept and
#'   terms starting with `sd__` or `cor__` is treated as a fixed effect.
#' * Random-effect standard deviations (class `"ranef_sd"`). For each grouping
#'   variable the standard deviation estimates across targets are averaged and
#'   used as the mean of a gamma prior, `gamma(mean, shape)`, on the standard
#'   deviation scale. A grouping variable with several standard deviation terms
#'   receives one prior for all of them.
#' * Dispersion (class `"fixef_disp"`). The log dispersion of each target
#'   depends strongly on its expression level, so a trend of log dispersion
#'   against log mean count is estimated and each target receives a normal
#'   prior centered on the trend at its own mean count. Targets whose mean count
#'   lies outside the range seen during estimation, including targets with zero
#'   counts, receive the prior at the nearest end of that range.
#'
#' The spread of estimates across targets includes the estimation error of
#' each estimate, so the resulting priors are somewhat wider than the true
#' between-target variation, which makes them conservative.
#'
#' The container for the second run must ask for the priors by name, using
#' `alist()` so that the data frame is built on the worker for each target:
#'
#' ```
#' arguments = alist(
#'   formula = ...,
#'   family = glmmTMB::nbinom2,
#'   priors = data.frame(prior = prior, class = class, coef = coef)
#' )
#' ```
#'
#' Priors are experimental in `glmmTMB`; see its `priors` help page and
#' vignette for the parameterisation of each distribution.
#'
#' @seealso [dispersion_evaluation()] to record the dispersion during the first
#'   run, and [seqwrap_compose()] for the `targetdata` argument.
#'
#' @examples
#' library(seqwrap)
#'
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'   dat <- simcounts(n_genes = 20, n_samples = 30, clusters = 15)
#'   counts <- dat$data
#'   metadata <- dat$metadata
#'   metadata$ln_libsize <- log(colSums(counts[, -1]))
#'
#'   # First run, without priors, recording the dispersion of each target
#'   container <- seqwrap_compose(
#'     modelfun = glmmTMB::glmmTMB,
#'     arguments = list(
#'       formula = y ~ x + (1 | cluster) + offset(ln_libsize),
#'       family = glmmTMB::nbinom2
#'     ),
#'     data = counts,
#'     metadata = metadata,
#'     samplename = "sample",
#'     eval_fun = dispersion_evaluation
#'   )
#'
#'   results <- seqwrap(container, cores = 1, verbose = FALSE)
#'
#'   # Priors pooled across targets, one data frame per target
#'   priors <- seqwrap_priors(results, data = counts)
#'   priors
#'   priors[[1]]
#'
#'   \donttest{
#'   # Second run, with the priors passed as target-specific data
#'   container_prior <- seqwrap_compose(
#'     modelfun = glmmTMB::glmmTMB,
#'     arguments = alist(
#'       formula = y ~ x + (1 | cluster) + offset(ln_libsize),
#'       family = glmmTMB::nbinom2,
#'       priors = data.frame(prior = prior, class = class, coef = coef)
#'     ),
#'     data = counts,
#'     metadata = metadata,
#'     samplename = "sample",
#'     targetdata = priors
#'   )
#'
#'   results_prior <- seqwrap(container_prior, subset = 1:5,
#'                            return_models = TRUE, cores = 1,
#'                            verbose = FALSE)
#'
#'   # The fitted model carries the priors it was given
#'   results_prior@models[[1]]$modelInfo$priors
#'   }
#' }
#'
#' @export
seqwrap_priors <- function(
  x,
  data,
  rownames = FALSE,
  terms = NULL,
  center = c("zero", "mean"),
  robust = FALSE,
  ranef = TRUE,
  shape = 2,
  dispersion = "dispersion",
  trend = c("loess", "constant"),
  span = 0.75,
  drop_warnings = FALSE,
  digits = 3,
  plot = FALSE
) {
  center <- match.arg(center)
  trend <- match.arg(trend)

  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    cli::cli_abort("{.arg plot} must be TRUE or FALSE.")
  }

  if (!is.numeric(shape) || length(shape) != 1L || !is.finite(shape) ||
        shape <= 0) {
    cli::cli_abort("{.arg shape} must be a single positive number.")
  }
  if (!is.numeric(span) || length(span) != 1L || !is.finite(span) ||
        span <= 0) {
    cli::cli_abort("{.arg span} must be a single positive number.")
  }
  if (!is.numeric(digits) || length(digits) != 1L || digits < 0) {
    cli::cli_abort("{.arg digits} must be a single non-negative number.")
  }

  # Combined summaries and evaluations, whichever form the run is in
  if (S7::S7_inherits(x, seqwrapResults)) {
    combined <- seqwrap_summarise(
      x,
      errors = FALSE,
      drop_warnings = drop_warnings,
      verbose = FALSE
    )
  } else if (is.list(x) && !is.null(x[["summaries"]])) {
    combined <- x
  } else {
    cli::cli_abort(
      "{.arg x} must be a seqwrapResults object or the list returned by
       {.fn seqwrap_summarise}."
    )
  }

  summaries <- combined[["summaries"]]
  evaluations <- combined[["evaluations"]]

  if (is.null(summaries) || !all(c("target", "term", "estimate") %in%
                                   colnames(summaries))) {
    cli::cli_abort(c(
      "No model summaries with {.field target}, {.field term} and
       {.field estimate} columns are available in {.arg x}.",
      "i" = "Priors are estimated from the summaries produced by
             {.fn generic_summary} or a summary function returning the same
             columns."
    ))
  }

  target_info <- sw_prior_targets(data, rownames)

  # Fixed effects, shared by all targets
  fixef <- sw_prior_fixef(summaries, terms, center, robust, digits)

  # Random-effect standard deviations, shared by all targets
  ranef_priors <- if (isTRUE(ranef)) {
    sw_prior_ranef(summaries, shape, robust, digits)
  } else {
    NULL
  }
  if (is.null(ranef_priors)) {
    ranef_priors <- list(table = NULL, estimates = NULL)
  }

  # Dispersion, specific to each target
  disp <- if (!is.null(dispersion)) {
    sw_prior_dispersion(
      evaluations,
      column = dispersion,
      ids = target_info$ids,
      log_mu = target_info$log_mu,
      trend = trend,
      span = span,
      robust = robust,
      digits = digits
    )
  } else {
    NULL
  }

  common <- rbind(fixef$table, ranef_priors$table)
  rownames(common) <- NULL

  # The estimates behind each prior are kept so that a user can inspect or
  # plot them in other ways than plot.seqwrap_priors() offers.
  estimates <- rbind(fixef$estimates, ranef_priors$estimates)
  rownames(estimates) <- NULL
  plot_data <- list(
    estimates = estimates,
    dispersion = if (is.null(disp)) NULL else disp$estimates
  )

  if (nrow(common) == 0L && is.null(disp)) {
    cli::cli_abort(
      "No priors could be built: no fixed-effect terms were found and no
       dispersion prior was requested."
    )
  }

  shared <- common[, c("prior", "class", "coef"), drop = FALSE]

  priors <- lapply(seq_along(target_info$ids), function(i) {
    out <- shared
    if (!is.null(disp)) {
      out <- rbind(
        out,
        data.frame(
          prior = disp$table$prior[i],
          class = "fixef_disp",
          coef = "1",
          stringsAsFactors = FALSE
        )
      )
    }
    rownames(out) <- NULL
    out
  })
  names(priors) <- target_info$ids

  out <- structure(
    priors,
    class = c("seqwrap_priors", "list"),
    common = common,
    dispersion = if (is.null(disp)) NULL else disp$table,
    trend = if (is.null(disp)) NULL else disp$model,
    n = length(unique(summaries$target)),
    data = plot_data
  )

  if (plot) {
    plot.seqwrap_priors(out)
  }

  out
}


#' Print method for seqwrap_priors objects
#'
#' @description
#' Printing a `seqwrap_priors` object lists the priors it holds and how they
#' were built.
#'
#' @details
#' The method is called as `print(x, ...)`, where `x` is the object returned by
#' [seqwrap_priors()]. Further arguments (`...`) are accepted for
#' compatibility with the [print()] generic but are not used.
#'
#' @return Invisibly returns the object.
#'
#' @usage NULL
#'
#' @method print seqwrap_priors
#' @name print.seqwrap_priors
# Registered in .onLoad() rather than through an S3method() directive. The S7
# print method for seqwrapResults leaves a `print` binding in this namespace,
# and R then files S3method(print, ...) registrations from NAMESPACE under the
# namespace's own methods table, where base::print never looks.
sw_print_priors <- function(x, ...) {
  common <- attr(x, "common")
  disp <- attr(x, "dispersion")
  trend <- attr(x, "trend")
  n <- attr(x, "n")

  cli::cli_h1("seqwrap priors")
  cli::cli_inform(
    "Priors for {length(x)} target{?s}, estimated from {n} fitted target{?s}."
  )

  for (cls in c("fixef", "ranef_sd")) {
    rows <- common[common$class == cls, , drop = FALSE]
    if (nrow(rows) == 0L) next
    label <- if (cls == "fixef") {
      "Fixed effects (class \"fixef\")"
    } else {
      "Random-effect standard deviations (class \"ranef_sd\")"
    }
    cli::cli_h3(label)
    bullets <- sprintf("%s: %s", rows$coef, rows$prior)
    cli::cli_inform(stats::setNames(bullets, rep("*", length(bullets))))
  }

  if (!is.null(disp)) {
    cli::cli_h3("Dispersion (class \"fixef_disp\")")
    if (is.null(trend)) {
      cli::cli_inform(c("*" = "All targets: {disp$prior[1]}"))
    } else {
      scale <- sw_prior_num(unique(disp$scale)[1], 3)
      cli::cli_inform(c(
        "*" = "normal(trend, {scale}) with the location read from a loess
               curve of log dispersion against log mean count"
      ))
    }
  }

  cli::cli_alert_info(
    "Pass this object as {.arg targetdata} to {.fn seqwrap_compose} and refer
     to {.code prior}, {.code class} and {.code coef} in {.arg arguments}."
  )

  invisible(x)
}


#' Plot priors over the estimates they were built from
#'
#' @description
#' Draws each prior held by a `seqwrap_priors` object over the per-target
#' estimates it summarises, as a visual check that the priors are reasonable
#' before they are used in a second run. Fixed effects and random-effect
#' standard deviations are shown as a histogram of the estimates across
#' targets with the prior density drawn on top. The dispersion prior is shown
#' as the log dispersion of each target against its log mean count, with the
#' fitted trend and a band of one prior standard deviation around it.
#'
#' @details
#' The plots are drawn with base graphics. The estimates behind them are kept
#' in the `data` attribute of the object (see [seqwrap_priors()]), so that a
#' different display can be built from the same values, for example with
#' `ggplot2`. The method returns that attribute invisibly.
#'
#' Random-effect priors are gamma distributions parameterised by their mean
#' and shape, as in `glmmTMB`, so the density drawn has shape `shape` and rate
#' `shape / mean`.
#'
#' @param x A `seqwrap_priors` object from [seqwrap_priors()].
#' @param which Character vector selecting the panels to draw, any of
#'   `"fixef"`, `"ranef"` and `"dispersion"`. Defaults to all three. Panels
#'   for prior classes the object does not hold are skipped.
#' @param trim The central fraction of the estimates shown in each histogram.
#'   Defaults to 0.99, so that a few targets whose fits went astray do not
#'   compress the histogram; estimates outside the shown range are counted in
#'   the panel subtitle rather than drawn. Set to 1 to show every estimate.
#' @param layout Logical, should the panels be arranged in a grid on the
#'   current device? Defaults to TRUE, in which case the graphical parameters
#'   are restored on exit. Set to FALSE to draw the panels one after another
#'   under a layout of your own, set with [graphics::par()] or
#'   [graphics::layout()].
#' @param ... Accepted for compatibility with the [plot()] generic, not used.
#'
#' @return Invisibly, the `data` attribute of `x`: a list with the data frames
#'   `estimates` and `dispersion` that the panels were drawn from.
#'
#' @seealso [seqwrap_priors()], whose `plot` argument calls this method.
#'
#' @examples
#' library(seqwrap)
#'
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'   dat <- simcounts(n_genes = 30, n_samples = 30, clusters = 15)
#'   counts <- dat$data
#'   metadata <- dat$metadata
#'   metadata$ln_libsize <- log(colSums(counts[, -1]))
#'
#'   first <- seqwrap_compose(
#'     modelfun = glmmTMB::glmmTMB,
#'     arguments = list(
#'       formula = y ~ x + (1 | cluster) + offset(ln_libsize),
#'       family = glmmTMB::nbinom2
#'     ),
#'     data = counts,
#'     metadata = metadata,
#'     samplename = "sample",
#'     eval_fun = dispersion_evaluation
#'   )
#'   first_results <- seqwrap(first, cores = 1, verbose = FALSE)
#'
#'   priors <- seqwrap_priors(first_results, data = counts)
#'   plot(priors)
#'
#'   # Only the dispersion trend
#'   plot(priors, which = "dispersion")
#'
#'   # The estimates behind the panels, for a display of your own
#'   str(attr(priors, "data"))
#' }
#'
#' @method plot seqwrap_priors
#' @export
plot.seqwrap_priors <- function(
  x,
  which = c("fixef", "ranef", "dispersion"),
  trim = 0.99,
  layout = TRUE,
  ...
) {
  which <- match.arg(which, several.ok = TRUE)
  if (!is.numeric(trim) || length(trim) != 1L || !is.finite(trim) ||
        trim <= 0 || trim > 1) {
    cli::cli_abort("{.arg trim} must be a single number in (0, 1].")
  }
  if (!is.logical(layout) || length(layout) != 1L || is.na(layout)) {
    cli::cli_abort("{.arg layout} must be TRUE or FALSE.")
  }

  common <- attr(x, "common")
  plot_data <- attr(x, "data")
  disp <- attr(x, "dispersion")

  if (is.null(plot_data)) {
    cli::cli_abort(
      "{.arg x} carries no {.field data} attribute; it was not built by
       {.fn seqwrap_priors}."
    )
  }

  # Collect the panels to draw before drawing, so that the layout can be
  # sized to the number of panels.
  panels <- list()
  if ("fixef" %in% which) {
    rows <- common[common$class == "fixef", , drop = FALSE]
    for (i in seq_len(nrow(rows))) {
      panels[[length(panels) + 1L]] <- list(class = "fixef", row = rows[i, ])
    }
  }
  if ("ranef" %in% which) {
    rows <- common[common$class == "ranef_sd", , drop = FALSE]
    for (i in seq_len(nrow(rows))) {
      panels[[length(panels) + 1L]] <- list(class = "ranef_sd", row = rows[i, ])
    }
  }
  if ("dispersion" %in% which && !is.null(disp)) {
    panels[[length(panels) + 1L]] <- list(class = "fixef_disp")
  }

  if (length(panels) == 0L) {
    cli::cli_inform(
      "Nothing to plot: the object holds no priors of the requested kind."
    )
    return(invisible(plot_data))
  }

  if (layout && length(panels) > 1L) {
    nc <- min(length(panels), 3L)
    nr <- ceiling(length(panels) / nc)
    old <- graphics::par(mfrow = c(nr, nc))
    on.exit(graphics::par(old), add = TRUE)
  }

  for (panel in panels) {
    if (panel$class == "fixef_disp") {
      sw_plot_prior_dispersion(plot_data$dispersion, disp, attr(x, "trend"))
    } else {
      sw_plot_prior_hist(panel$row, plot_data$estimates, trim)
    }
  }

  invisible(plot_data)
}


#' Histogram of estimates with a prior density on top
#'
#' @param row One row of the `common` attribute of a `seqwrap_priors` object.
#' @param estimates The `estimates` element of the `data` attribute.
#' @param trim Central fraction of the estimates to show.
#' @return NULL, invisibly. Called for its side effect.
#' @keywords internal
#' @noRd
sw_plot_prior_hist <- function(row, estimates, trim) {
  est <- estimates$estimate[
    estimates$class == row$class & estimates$coef == row$coef
  ]
  est <- est[is.finite(est)]

  if (row$class == "fixef") {
    label <- "Fixed effect"
    xlab <- sprintf("Estimate of %s across targets", row$coef)
    density <- function(v) stats::dnorm(v, row$location, row$scale)
  } else {
    label <- "Random-effect SD"
    xlab <- sprintf("SD estimate for %s across targets", row$coef)
    # glmmTMB parameterises gamma(mean, shape)
    density <- function(v) {
      stats::dgamma(v, shape = row$scale, rate = row$scale / row$location)
    }
  }
  main <- sprintf("%s %s: %s", label, row$coef, row$prior)

  if (length(est) == 0L) {
    graphics::plot.new()
    graphics::title(main = main, cex.main = 0.95)
    graphics::text(0.5, 0.5, "No finite estimates")
    return(invisible(NULL))
  }

  # Show the central part of the estimates and count the rest, so that a few
  # runaway fits do not flatten the histogram.
  tail <- (1 - trim) / 2
  lim <- stats::quantile(est, c(tail, 1 - tail), names = FALSE)
  # Always show the prior itself, at least to two standard deviations
  if (row$class == "fixef") {
    lim <- range(lim, row$location + c(-2, 2) * row$scale)
  } else {
    lim <- range(0, lim, row$location * 2)
  }
  shown <- est >= lim[1] & est <= lim[2]
  n_out <- sum(!shown)

  h <- graphics::hist(est[shown], plot = FALSE)
  grid <- seq(lim[1], lim[2], length.out = 200)
  dens <- density(grid)
  dens[!is.finite(dens)] <- NA
  ylim <- c(0, max(c(h$density, dens), na.rm = TRUE) * 1.05)

  graphics::hist(
    est[shown],
    breaks = h$breaks,
    freq = FALSE,
    col = "grey88",
    border = "grey60",
    xlim = lim,
    ylim = ylim,
    main = main,
    cex.main = 0.95,
    xlab = xlab,
    ylab = "Density"
  )
  graphics::lines(grid, dens, col = "firebrick", lwd = 2)
  graphics::abline(v = row$location, col = "firebrick", lty = 2)

  sub <- sprintf("%d target%s", length(est), if (length(est) == 1L) "" else "s")
  if (n_out > 0L) {
    sub <- sprintf("%s, %d outside the shown range", sub, n_out)
  }
  graphics::mtext(sub, side = 3, line = 0.2, cex = 0.75, col = "grey30")

  invisible(NULL)
}


#' Log dispersion against log mean count with the prior trend
#'
#' @param observed The `dispersion` element of the `data` attribute.
#' @param table The `dispersion` attribute, one row per target.
#' @param model The loess fit, or NULL for a constant prior.
#' @return NULL, invisibly. Called for its side effect.
#' @keywords internal
#' @noRd
sw_plot_prior_dispersion <- function(observed, table, model) {
  scale <- table$scale[1]
  xlim <- range(observed$log_mu)
  ylim <- range(observed$dispersion)

  if (is.null(model)) {
    grid <- xlim
    fit <- rep(table$location[1], 2L)
    main <- sprintf("Dispersion: %s", table$prior[1])
  } else {
    grid <- seq(xlim[1], xlim[2], length.out = 200)
    fit <- as.numeric(stats::predict(model, newdata = data.frame(log_mu = grid)))
    main <- sprintf("Dispersion: normal(trend, %s)", sw_prior_num(scale, 3))
  }
  ylim <- range(ylim, fit - scale, fit + scale)

  graphics::plot(
    observed$log_mu,
    observed$dispersion,
    pch = 16,
    cex = 0.6,
    col = grDevices::adjustcolor("grey30", alpha.f = 0.5),
    xlim = xlim,
    ylim = ylim,
    main = main,
    cex.main = 0.95,
    xlab = "log mean count",
    ylab = "log dispersion"
  )
  graphics::polygon(
    c(grid, rev(grid)),
    c(fit - scale, rev(fit + scale)),
    col = grDevices::adjustcolor("firebrick", alpha.f = 0.15),
    border = NA
  )
  graphics::lines(grid, fit, col = "firebrick", lwd = 2)

  sub <- sprintf(
    "%d target%s, band is one prior SD",
    nrow(observed),
    if (nrow(observed) == 1L) "" else "s"
  )
  graphics::mtext(sub, side = 3, line = 0.2, cex = 0.75, col = "grey30")

  invisible(NULL)
}


#' Target identifiers and log mean counts from a target data frame
#'
#' @param data A data frame of targets by samples.
#' @param rownames Are identifiers held in the row names?
#' @return A list with `ids` (character) and `log_mu` (numeric).
#' @keywords internal
#' @noRd
sw_prior_targets <- function(data, rownames = FALSE) {
  if (is.list(data) && !is.data.frame(data)) {
    cli::cli_abort(c(
      "{.arg data} must be a single data frame.",
      "i" = "When {.fn seqwrap_compose} was given a list of data frames, pass
             the element holding the counts."
    ))
  }
  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame.")
  }

  data <- as.data.frame(data)

  if (isTRUE(rownames)) {
    ids <- rownames(data)
    values <- data
  } else {
    ids <- as.character(data[[1]])
    values <- data[, -1, drop = FALSE]
  }

  if (anyNA(ids) || any(!nzchar(trimws(ids)))) {
    cli::cli_abort("Target identifiers in {.arg data} must not be missing.")
  }
  if (anyDuplicated(ids)) {
    cli::cli_abort("Target identifiers in {.arg data} must be unique.")
  }

  values <- suppressWarnings(
    matrix(as.numeric(as.matrix(values)), nrow = nrow(values))
  )
  log_mu <- log(rowMeans(values, na.rm = TRUE))

  list(ids = ids, log_mu = log_mu)
}


#' Format numbers for use inside a prior specification
#'
#' @param v A numeric vector.
#' @param digits Number of decimals to keep.
#' @return A character vector without scientific notation.
#' @keywords internal
#' @noRd
sw_prior_num <- function(v, digits) {
  vapply(
    v,
    function(z) {
      format(
        round(z, digits),
        scientific = FALSE,
        trim = TRUE,
        drop0trailing = TRUE
      )
    },
    character(1)
  )
}


#' Location and scale of a set of estimates
#'
#' @param est A numeric vector of estimates, may contain NA.
#' @param robust Use median and MAD rather than mean and SD?
#' @return A list with `location`, `scale` and `n`.
#' @keywords internal
#' @noRd
sw_prior_moments <- function(est, robust = FALSE) {
  est <- est[is.finite(est)]
  n <- length(est)
  if (n < 2L) {
    return(list(location = NA_real_, scale = NA_real_, n = n))
  }
  if (robust) {
    list(location = stats::median(est), scale = stats::mad(est), n = n)
  } else {
    list(location = mean(est), scale = stats::sd(est), n = n)
  }
}


#' Fixed-effect priors from combined summaries
#'
#' @param summaries The combined summaries data frame.
#' @param terms Terms to include, NULL for all except the intercept.
#' @param center "zero" or "mean".
#' @param robust Use robust location and scale?
#' @param digits Decimals kept in the prior strings.
#' @return A list with `table`, a data frame with columns prior, class, coef,
#'   location, scale, n, and `estimates`, a data frame with the per-target
#'   estimates (columns target, class, coef, estimate) that the table
#'   summarises.
#' @keywords internal
#' @noRd
sw_prior_fixef <- function(summaries, terms, center, robust, digits) {
  term <- as.character(summaries$term)
  keep <- !is.na(term)

  # generic_summary() (broom.mixed) labels effects; a custom summary may not
  if ("effect" %in% colnames(summaries)) {
    keep <- keep & summaries$effect %in% "fixed"
  } else {
    keep <- keep & !grepl("^(sd|cor)__", term)
  }
  if ("component" %in% colnames(summaries)) {
    component <- summaries$component
    keep <- keep & (is.na(component) | component %in% "cond")
  }

  if (is.null(terms)) {
    keep <- keep & term != "(Intercept)"
  } else {
    missing <- setdiff(terms, unique(term[keep]))
    if (length(missing) > 0L) {
      cli::cli_abort(c(
        "Some {.arg terms} are not fixed-effect terms in the summaries.",
        "x" = "Not found: {.val {missing}}.",
        "i" = "Available: {.val {unique(term[keep])}}."
      ))
    }
    keep <- keep & term %in% terms
  }

  use <- summaries[keep, , drop = FALSE]
  if (nrow(use) == 0L) {
    return(list(table = sw_prior_empty(), estimates = sw_prior_empty_est()))
  }

  # Keep the order in which terms appear in the summaries
  term_levels <- unique(as.character(use$term))
  if (!is.null(terms)) term_levels <- terms

  use <- use[as.character(use$term) %in% term_levels, , drop = FALSE]
  estimates <- data.frame(
    target = as.character(use$target),
    class = "fixef",
    coef = as.character(use$term),
    estimate = suppressWarnings(as.numeric(use$estimate)),
    stringsAsFactors = FALSE
  )

  rows <- lapply(term_levels, function(tm) {
    m <- sw_prior_moments(use$estimate[use$term == tm], robust = robust)
    location <- if (center == "zero") 0 else m$location
    data.frame(
      class = "fixef",
      coef = tm,
      location = location,
      scale = m$scale,
      n = m$n,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)

  sw_prior_check_scale(out, "fixed-effect term", digits)

  out$prior <- sprintf(
    "normal(%s, %s)",
    sw_prior_num(out$location, digits),
    sw_prior_num(out$scale, digits)
  )
  list(
    table = out[, c("prior", "class", "coef", "location", "scale", "n")],
    estimates = estimates
  )
}


#' Random-effect standard deviation priors from combined summaries
#'
#' @inheritParams sw_prior_fixef
#' @param shape Shape of the gamma prior.
#' @return A list as from `sw_prior_fixef()`, or NULL when the summaries hold
#'   no random-effect standard deviations.
#' @keywords internal
#' @noRd
sw_prior_ranef <- function(summaries, shape, robust, digits) {
  term <- as.character(summaries$term)
  keep <- !is.na(term) & grepl("^sd__", term)
  if ("effect" %in% colnames(summaries)) {
    keep <- keep & summaries$effect %in% "ran_pars"
  }
  if ("component" %in% colnames(summaries)) {
    component <- summaries$component
    keep <- keep & (is.na(component) | component %in% "cond")
  }

  if (!any(keep)) {
    return(NULL)
  }

  # The prior is addressed by grouping variable, which only a summary with a
  # group column can provide. Residual standard deviations are not random
  # effects.
  if (!"group" %in% colnames(summaries)) {
    cli::cli_inform(c(
      "i" = "Random-effect terms were found but the summaries have no
             {.field group} column, so no random-effect priors are built."
    ))
    return(NULL)
  }
  keep <- keep & !summaries$group %in% "Residual" & !is.na(summaries$group)
  if (!any(keep)) {
    return(NULL)
  }

  use <- summaries[keep, , drop = FALSE]
  groups <- unique(as.character(use$group))

  estimates <- data.frame(
    target = as.character(use$target),
    class = "ranef_sd",
    coef = as.character(use$group),
    estimate = suppressWarnings(as.numeric(use$estimate)),
    stringsAsFactors = FALSE
  )

  rows <- lapply(groups, function(g) {
    m <- sw_prior_moments(use$estimate[use$group == g], robust = robust)
    data.frame(
      class = "ranef_sd",
      coef = g,
      location = m$location,
      scale = shape,
      n = m$n,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)

  # A gamma prior needs a positive mean. When nearly every target has a
  # random-effect standard deviation at the boundary the location rounds to
  # zero, which glmmTMB rejects; skipping the prior keeps the workflow going.
  bad <- !is.finite(out$location) | round(out$location, digits) <= 0
  if (any(bad)) {
    bad_coef <- out$coef[bad]
    bad_location <- signif(out$location[bad], 3)
    cli::cli_warn(c(
      "No random-effect prior is built for {.val {bad_coef}}: the standard
       deviation estimates across targets are essentially zero
       (location {.val {bad_location}}).",
      "i" = "Set {.code ranef = FALSE} to skip random-effect priors without
             this warning."
    ))
    out <- out[!bad, , drop = FALSE]
    estimates <- estimates[!estimates$coef %in% bad_coef, , drop = FALSE]
    if (nrow(out) == 0L) {
      return(NULL)
    }
  }

  out$prior <- sprintf(
    "gamma(%s, %s)",
    sw_prior_num(out$location, digits),
    sw_prior_num(shape, digits)
  )
  list(
    table = out[, c("prior", "class", "coef", "location", "scale", "n")],
    estimates = estimates
  )
}


#' Dispersion priors, one per target
#'
#' @param evaluations The combined evaluations data frame, or NULL.
#' @param column Name of the log dispersion column.
#' @param ids Target identifiers from the target data.
#' @param log_mu Log mean count of each target, aligned with `ids`.
#' @param trend "loess" or "constant".
#' @param span Loess span.
#' @param robust Use a robust loess fit and robust moments?
#' @param digits Decimals kept in the prior strings.
#' @return A list with `table`, a data frame with one row per target,
#'   `model`, the loess fit or NULL, and `estimates`, a data frame with the
#'   target, log mean count and log dispersion of every target used to
#'   estimate the trend.
#' @keywords internal
#' @noRd
sw_prior_dispersion <- function(
  evaluations,
  column,
  ids,
  log_mu,
  trend,
  span,
  robust,
  digits
) {
  if (is.null(evaluations) || !column %in% colnames(evaluations)) {
    cli::cli_abort(c(
      "No {.field {column}} column is available in the evaluations.",
      "i" = "Run {.fn seqwrap} with {.code eval_fun = dispersion_evaluation}
             to record the dispersion of each target, or set
             {.code dispersion = NULL} to build no dispersion prior."
    ))
  }

  fitted <- data.frame(
    target = as.character(evaluations$target),
    dispersion = suppressWarnings(as.numeric(evaluations[[column]])),
    stringsAsFactors = FALSE
  )
  fitted$log_mu <- log_mu[match(fitted$target, ids)]

  if (!any(fitted$target %in% ids)) {
    cli::cli_abort(
      "None of the targets in the evaluations are found in {.arg data}."
    )
  }

  fitted <- fitted[is.finite(fitted$dispersion) & is.finite(fitted$log_mu), ]
  fitted <- fitted[!duplicated(fitted$target), ]

  min_n <- if (trend == "loess") 10L else 2L
  if (nrow(fitted) < min_n) {
    cli::cli_abort(c(
      "At least {min_n} targets with a finite dispersion and mean count are
       needed to estimate a {.val {trend}} dispersion prior.",
      "x" = "Only {nrow(fitted)} {?is/are} available."
    ))
  }

  if (trend == "constant") {
    m <- sw_prior_moments(fitted$dispersion, robust = robust)
    location <- rep(m$location, length(ids))
    scale <- m$scale
    model <- NULL
  } else {
    model <- stats::loess(
      dispersion ~ log_mu,
      data = fitted,
      span = span,
      family = if (robust) "symmetric" else "gaussian"
    )

    # loess does not extrapolate, so targets outside the range seen during
    # estimation (including all-zero targets with log_mu = -Inf) take the
    # prior at the nearest end of the range.
    lo <- min(fitted$log_mu)
    hi <- max(fitted$log_mu)
    at <- log_mu
    at[!is.finite(at)] <- lo
    at <- pmin(pmax(at, lo), hi)

    location <- as.numeric(stats::predict(model, newdata = data.frame(log_mu = at)))
    scale <- model$s
  }

  if (!is.finite(scale) || round(scale, digits) <= 0) {
    cli::cli_abort(c(
      "The dispersion prior needs a positive standard deviation.",
      "x" = "The estimated value is {.val {scale}}.",
      "i" = "Increase {.arg digits} if the value is small but positive."
    ))
  }

  table <- data.frame(
    target = ids,
    log_mu = log_mu,
    location = location,
    scale = scale,
    prior = sprintf(
      "normal(%s, %s)",
      sw_prior_num(location, digits),
      sw_prior_num(scale, digits)
    ),
    stringsAsFactors = FALSE
  )
  rownames(table) <- NULL

  rownames(fitted) <- NULL
  list(table = table, model = model, estimates = fitted)
}


#' Abort when a prior scale is missing or rounds to zero
#'
#' @param out A data frame with coef, scale and n columns.
#' @param what Label used in the message.
#' @param digits Decimals kept in the prior strings.
#' @keywords internal
#' @noRd
sw_prior_check_scale <- function(out, what, digits) {
  few <- out$n < 2L
  if (any(few)) {
    cli::cli_abort(c(
      "At least two targets with an estimate are needed for each {what}.",
      "x" = "{.val {out$coef[few]}} {?has/have} {out$n[few]} estimate{?s}."
    ))
  }
  zero <- !is.finite(out$scale) | round(out$scale, digits) <= 0
  if (any(zero)) {
    cli::cli_abort(c(
      "The prior for {what} {.val {out$coef[zero]}} has no positive scale.",
      "i" = "Estimates may be identical across targets, or the value may be
             smaller than {.arg digits} allows; increase {.arg digits}."
    ))
  }
  invisible(TRUE)
}


#' An empty prior table
#'
#' @return A zero row data frame with the prior table columns.
#' @keywords internal
#' @noRd
sw_prior_empty <- function() {
  data.frame(
    prior = character(0),
    class = character(0),
    coef = character(0),
    location = numeric(0),
    scale = numeric(0),
    n = integer(0),
    stringsAsFactors = FALSE
  )
}


#' An empty estimates data frame for priors
#'
#' @return A zero-row data frame with the columns of the `estimates` element
#'   of the `data` attribute of a `seqwrap_priors` object.
#' @keywords internal
#' @noRd
sw_prior_empty_est <- function() {
  data.frame(
    target = character(0),
    class = character(0),
    coef = character(0),
    estimate = numeric(0),
    stringsAsFactors = FALSE
  )
}
