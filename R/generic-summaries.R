#' Generic summaries for model parameter estimates.
#'
#' When no summary function is provided to seqwrap (summary_fun), this function
#' uses broom(.mixed)::tidy to give a table of model parameter estimates.
#'
#' @param x A model fitted in seqwrap.
#'
#' @returns A tidy data frame possible to bind using seqwrap_summarise
#' @examples
#' # The generic summary works on model objects
#' library(seqwrap)
#'
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'   library(glmmTMB)
#'
#'   dat <- simcounts(n_genes = 1000,
#'                    n_samples = 30,
#'                    beta_0 = 5,
#'                    overdispersion_min_max = c(1, 10))
#'
#'   counts <- dat$data
#'   metadata <- dat$metadata
#'   metadata$ln_libsize <- log(colSums(counts[,-1]))
#'
#'   # Save counts for gene i in the same data frame
#'   metadata$y <- as.integer(counts[1, -1])
#'
#'   m <- glmmTMB(y ~ as.factor(x) + offset(ln_libsize),
#'                data = metadata,
#'                family = nbinom2)
#'
#'   generic_summary(m)
#'   generic_evaluation(m)
#' }
#'
#' @export
generic_summary <- function(x) {
  out <- broom.mixed::tidy(x)
  return(out)
}


#' Build one row of convergence diagnostics
#'
#' All `generic_evaluation()` methods return the same columns so that results
#' from any engine can be bound together by `seqwrap_summarise()`. Arguments
#' accept NULL, which becomes the appropriate missing value, so methods can pass
#' model components straight through without checking whether they exist.
#'
#' @param engine Character, the model class.
#' @param converged Logical, did the fitting algorithm report convergence?
#' @param code Integer convergence code from the optimiser.
#' @param message Character message from the optimiser.
#' @param singular Logical, is the fit degenerate?
#' @param iterations Integer iteration or function evaluation count.
#' @return A one row tibble.
#' @keywords internal
#' @noRd
sw_eval_row <- function(
  engine,
  converged = NA,
  code = NA_integer_,
  message = NA_character_,
  singular = NA,
  iterations = NA_integer_
) {
  # Collapse multi-element messages, several engines report more than one.
  if (length(message) > 1L) {
    message <- paste(message, collapse = "; ")
  }

  tibble::tibble(
    engine = as.character(engine)[1],
    converged = as.logical(converged)[1],
    code = suppressWarnings(as.integer(code)[1]),
    message = as.character(message)[1],
    singular = as.logical(singular)[1],
    iterations = suppressWarnings(as.integer(iterations)[1])
  )
}


#' Generic evaluation for model fits
#'
#' When no evaluation function is provided to seqwrap (eval_fun), this function
#' collects the convergence diagnostics reported by the fitting algorithm. It is
#' an S3 generic, so the diagnostics are read from whichever model object the
#' engine produced. Methods are supplied for `glmmTMB`, `lme4` (`merMod`),
#' `nlme` (`lme` and `gls`), `mgcv` (`gam` and `bam`), `MASS::glm.nb`
#' (`negbin`), and the base `glm` and `lm` classes. Unrecognised model classes
#' fall through to a default method.
#'
#' @param x A model fitted in seqwrap, or NULL when fitting failed.
#' @param ... Currently unused, present so that methods can be extended.
#'
#' @returns A one row tibble that can be combined using `seqwrap_summarise()`,
#' or NULL when `x` is NULL. The columns are:
#'
#' \describe{
#'   \item{engine}{The model class the diagnostics were read from.}
#'   \item{converged}{Whether the fitting algorithm reported successful
#'     convergence. NA when the engine reports no convergence signal.}
#'   \item{code}{The numeric convergence code, where the engine supplies one.
#'     Zero conventionally indicates success.}
#'   \item{message}{Any convergence message, NA when the fit was clean.}
#'   \item{singular}{Whether the fit is degenerate: on the boundary of the
#'     parameter space, rank deficient, or with a covariance matrix that is not
#'     positive definite.}
#'   \item{iterations}{Iterations or function evaluations used, where reported.}
#' }
#'
#' @details
#' `converged` and `singular` are deliberately separate. A fit can satisfy the
#' optimiser's convergence criterion and still be untrustworthy: `glmmTMB`
#' regularly returns convergence code 0 alongside a non-positive-definite
#' Hessian, and `lme4` returns code 0 for singular random effect structures.
#' Screening a large set of fits therefore usually means requiring
#' `converged & !singular` rather than `converged` alone.
#'
#' Earlier versions of this function ran DHARMa residual simulations. That is
#' several orders of magnitude more expensive per target, which is prohibitive
#' at the target counts this package is built for, and DHARMa is now an optional
#' dependency. Use `residual_diagnostics()` to get the previous behaviour.
#'
#' @examples
#' # The generic evaluation works on model objects
#' library(seqwrap)
#'
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'   library(glmmTMB)
#'
#'   dat <- simcounts(n_genes = 1000,
#'                    n_samples = 30,
#'                    beta_0 = 5,
#'                    overdispersion_min_max = c(1, 10))
#'
#'   counts <- dat$data
#'   metadata <- dat$metadata
#'   metadata$ln_libsize <- log(colSums(counts[,-1]))
#'
#'   # Save counts for gene i in the same data frame
#'   metadata$y <- as.integer(counts[1, -1])
#'
#'   m <- glmmTMB(y ~ as.factor(x) + offset(ln_libsize),
#'                data = metadata,
#'                family = nbinom2)
#'
#'   generic_evaluation(m)
#' }
#'
#' # Base engines are supported as well
#' fit <- stats::glm(mpg ~ wt, data = mtcars, family = stats::gaussian())
#' generic_evaluation(fit)
#'
#' @export
generic_evaluation <- function(x, ...) {
  UseMethod("generic_evaluation")
}


#' @rdname generic_evaluation
#' @export
generic_evaluation.default <- function(x, ...) {
  # Several engines follow the base R convention of a `converged` element.
  converged <- tryCatch(x[["converged"]], error = function(e) NULL)

  sw_eval_row(
    engine = class(x)[1],
    converged = if (is.null(converged)) NA else isTRUE(converged)
  )
}


#' @rdname generic_evaluation
#' @export
`generic_evaluation.NULL` <- function(x, ...) {
  # Fitting failed for this target; the error itself is recorded in the errors
  # slot of the seqwrap results.
  NULL
}


#' @rdname generic_evaluation
#' @export
generic_evaluation.glmmTMB <- function(x, ...) {
  fit <- x$fit
  pd_hess <- x$sdr$pdHess

  sw_eval_row(
    engine = "glmmTMB",
    converged = if (is.null(fit$convergence)) NA else fit$convergence == 0,
    code = fit$convergence,
    message = fit$message,
    # A non-positive-definite Hessian is glmmTMB's signal that the fit is
    # degenerate, and it occurs independently of the optimiser's own code.
    singular = if (is.null(pd_hess)) NA else !isTRUE(pd_hess),
    iterations = fit$iterations
  )
}


#' @rdname generic_evaluation
#' @export
generic_evaluation.merMod <- function(x, ...) {
  optinfo <- x@optinfo
  messages <- optinfo$conv$lme4$messages

  singular <- tryCatch(
    if (requireNamespace("lme4", quietly = TRUE)) lme4::isSingular(x) else NA,
    error = function(e) NA
  )

  sw_eval_row(
    engine = class(x)[1],
    converged = if (is.null(optinfo$conv$opt)) NA else optinfo$conv$opt == 0,
    code = optinfo$conv$opt,
    message = if (length(messages) > 0L) messages else NA_character_,
    singular = singular,
    # lme4 reports function evaluations rather than outer iterations.
    iterations = optinfo$feval
  )
}


#' @rdname generic_evaluation
#' @export
generic_evaluation.negbin <- function(x, ...) {
  out <- generic_evaluation.glm(x, ...)
  out$engine <- "negbin"

  # glm.nb estimates theta in a separate loop that can fail on its own while
  # the IRLS step reports success.
  if (!is.null(x$th.warn)) {
    theta_msg <- paste0("theta estimation: ", paste(x$th.warn, collapse = "; "))
    out$message <- if (is.na(out$message)) {
      theta_msg
    } else {
      paste(out$message, theta_msg, sep = "; ")
    }
    out$converged <- FALSE
  }

  out
}


#' @rdname generic_evaluation
#' @export
generic_evaluation.glm <- function(x, ...) {
  boundary <- isTRUE(x$boundary)
  aliased <- anyNA(stats::coef(x))

  sw_eval_row(
    engine = class(x)[1],
    converged = if (is.null(x$converged)) NA else isTRUE(x$converged),
    code = if (isTRUE(x$converged)) 0L else 1L,
    message = if (boundary) {
      "fit reached the boundary of the parameter space"
    } else if (aliased) {
      "rank deficient fit, some coefficients are aliased"
    } else {
      NA_character_
    },
    singular = boundary || aliased,
    iterations = x$iter
  )
}


#' @rdname generic_evaluation
#' @export
generic_evaluation.lm <- function(x, ...) {
  aliased <- anyNA(stats::coef(x))

  sw_eval_row(
    engine = "lm",
    # lm is solved directly, so there is no iterative convergence to report.
    converged = TRUE,
    code = 0L,
    message = if (aliased) {
      "rank deficient fit, some coefficients are aliased"
    } else {
      NA_character_
    },
    singular = aliased,
    iterations = NA_integer_
  )
}


#' @rdname generic_evaluation
#' @export
generic_evaluation.lme <- function(x, ...) {
  # nlme stores a character string in apVar when the approximate
  # variance-covariance matrix could not be computed.
  ap_var <- x$apVar
  degenerate <- is.character(ap_var)

  sw_eval_row(
    engine = class(x)[1],
    # lme() signals failure by raising an error rather than by returning a
    # model carrying a convergence code.
    converged = TRUE,
    code = 0L,
    message = if (degenerate) as.character(ap_var) else NA_character_,
    singular = degenerate,
    iterations = x$numIter
  )
}


#' @rdname generic_evaluation
#' @export
generic_evaluation.gls <- function(x, ...) {
  ap_var <- x$apVar
  degenerate <- is.character(ap_var)

  sw_eval_row(
    engine = class(x)[1],
    converged = TRUE,
    code = 0L,
    message = if (degenerate) as.character(ap_var) else NA_character_,
    singular = degenerate,
    iterations = x$numIter
  )
}


#' @rdname generic_evaluation
#' @export
generic_evaluation.gam <- function(x, ...) {
  conv <- x$mgcv.conv

  # x$converged covers the penalised IRLS step, mgcv.conv the smoothing
  # parameter optimisation. Both have to succeed.
  irls_ok <- if (is.null(x$converged)) NA else isTRUE(x$converged)
  outer_ok <- if (is.null(conv$fully.converged)) {
    NA
  } else {
    isTRUE(conv$fully.converged)
  }
  converged <- if (is.na(irls_ok) && is.na(outer_ok)) {
    NA
  } else {
    isTRUE(irls_ok) && !isFALSE(outer_ok)
  }

  sw_eval_row(
    engine = class(x)[1],
    converged = converged,
    code = if (isTRUE(converged)) 0L else 1L,
    message = if (isFALSE(outer_ok)) {
      "smoothing parameter optimisation did not fully converge"
    } else {
      NA_character_
    },
    singular = if (is.null(conv$hess.pos.def)) NA else !isTRUE(conv$hess.pos.def),
    iterations = if (!is.null(conv$iter)) conv$iter else x$iter
  )
}


#' Residual diagnostics using DHARMa
#'
#' Simulates scaled residuals with DHARMa and tests them for uniformity,
#' dispersion and outliers. This was the behaviour of `generic_evaluation()` in
#' seqwrap 0.7.0 and earlier, and is kept here as an opt-in evaluation function.
#'
#' @param x A model fitted in seqwrap, or NULL when fitting failed.
#' @param n Number of simulations passed to `DHARMa::simulateResiduals()`.
#'
#' @returns A one row tibble of p-values, or NULL when `x` is NULL.
#'
#' @details
#' Simulating residuals is far more expensive than reading a convergence code:
#' each target is refit-equivalent `n` times. On data sets with many thousands
#' of targets this dominates the total run time, which is why it is no longer
#' the default. Consider applying it to a subset of targets, or lowering `n`.
#'
#' DHARMa is an optional dependency, so it has to be installed before this
#' function can be used.
#'
#' @examples
#' library(seqwrap)
#'
#' if (requireNamespace("DHARMa", quietly = TRUE)) {
#'   fit <- stats::glm(mpg ~ wt, data = mtcars, family = stats::gaussian())
#'   residual_diagnostics(fit, n = 50)
#' }
#'
#' @export
residual_diagnostics <- function(x, n = 250) {
  if (is.null(x)) {
    return(NULL)
  }

  if (!requireNamespace("DHARMa", quietly = TRUE)) {
    stop(
      "Package 'DHARMa' is needed for residual_diagnostics(). ",
      "Install it with install.packages('DHARMa').",
      call. = FALSE
    )
  }

  sim_resid <- DHARMa::simulateResiduals(x, n = n, plot = FALSE)

  tibble::tibble(
    uniformity = DHARMa::testUniformity(sim_resid, plot = FALSE)$p.value,
    dispersion = DHARMa::testDispersion(sim_resid, plot = FALSE)$p.value,
    outliers = DHARMa::testOutliers(sim_resid, plot = FALSE)$p.value
  )
}
