#' A function to fit models with a chosen fitting algorithm, used in seqwrap.
#' @param fun Name of a fitting function, like glmmTMB::glmmTMB
#' @param arg A list of arguments that can be evaluated by the fitting function
#' @param vars A list of variables to be used in the fitting function
#' @return A fitted model object returned by `fun`.
#' @keywords internal
#' @noRd
fit_fun <- function(fun, arg, vars = NULL) {
  if (!is.null(vars)) {
    # Convert data frame to list if needed
    if (is.data.frame(vars) && nrow(vars) == 1) {
      vars <- as.list(vars)
    }

    # Create an evaluation environment with the variables
    eval_env <- list2env(vars, parent = environment())

    # Helper function to check if x is a formula or unevaluated formula expression
    is_formula_expr <- function(x) {
      # Check for actual formula object
      if (inherits(x, "formula")) return(TRUE)
      # Check for unevaluated formula expression (call with ~ operator)
      if (is.call(x) && length(x) >= 1 && identical(x[[1]], as.symbol("~"))) return(TRUE)
      return(FALSE)
    }

    # Recursive function to replace and evaluate
    replace_and_eval <- function(x, env) {
      if (is_formula_expr(x)) {
        # Don't replace variables inside formulas - they reference data columns
        return(x)
      } else if (is.symbol(x)) {
        # If it's a symbol in our vars, return its value
        var_name <- as.character(x)
        if (var_name %in% names(vars)) {
          return(vars[[var_name]])
        }
        return(x)
      } else if (is.call(x)) {
        # Process each element of the call recursively first
        for (i in seq_along(x)) {
          x[[i]] <- replace_and_eval(x[[i]], env)
        }

        # After replacement, try to evaluate the call
        tryCatch({
          return(eval(x, envir = env))
        }, error = function(e) {
          # If evaluation fails, return the modified call
          return(x)
        })
      } else if (is.pairlist(x)) {
        # Handle pairlists (arguments in alist)
        return(as.pairlist(lapply(x, function(elem) replace_and_eval(elem, env))))
      } else if (is.list(x)) {
        # Handle regular lists
        return(lapply(x, function(elem) replace_and_eval(elem, env)))
      } else {
        return(x)
      }
    }

    # Process the argument list
    arg <- replace_and_eval(arg, eval_env)
  }

  fittedmodel <- do.call(fun, arg)

  return(fittedmodel)
}

#' Special fitting function for lme
#' @param fun The fitting function, nlme::lme
#' @param arg A list of arguments that can be evaluated by the fitting function
#' @param vars A list of variables to be used in the fitting function
#' @return A fitted `lme`/`gls` model object returned by `fun`.
#' @keywords internal
#' @noRd
fit_fun_lme <- function(fun, arg, vars = NULL) {

  # Create an environment for the evaluation
  eval_env <- new.env(parent = environment())

  # If we have target-specific variables, handle replacement and evaluation
  if (!is.null(vars)) {
    # Convert data frame to list if needed
    if (is.data.frame(vars) && nrow(vars) == 1) {
      vars <- as.list(vars)
    }

    # Add variables to the evaluation environment
    for (var_name in names(vars)) {
      eval_env[[var_name]] <- vars[[var_name]]
    }

    # Helper function to check if x is a formula or unevaluated formula expression
    is_formula_expr <- function(x) {
      # Check for actual formula object
      if (inherits(x, "formula")) return(TRUE)
      # Check for unevaluated formula expression (call with ~ operator)
      if (is.call(x) && length(x) >= 1 && identical(x[[1]], as.symbol("~"))) return(TRUE)
      return(FALSE)
    }

    # Recursive function to replace and evaluate
    replace_and_eval <- function(x, env) {
      if (is.symbol(x)) {
        # If it's a symbol in our vars, return its value
        var_name <- as.character(x)
        if (var_name %in% names(vars)) {
          return(vars[[var_name]])
        }
        return(x)
      } else if (is_formula_expr(x)) {
        # Don't replace variables inside formulas - they reference data columns
        return(x)
      } else if (is.call(x)) {
        # Process each element of the call recursively first
        for (i in seq_along(x)) {
          x[[i]] <- replace_and_eval(x[[i]], env)
        }

        # After replacement, try to evaluate the call
        tryCatch({
          return(eval(x, envir = env))
        }, error = function(e) {
          # If evaluation fails, return the modified call
          return(x)
        })
      } else if (is.pairlist(x)) {
        # Handle pairlists (arguments in alist)
        return(as.pairlist(lapply(x, function(elem) replace_and_eval(elem, env))))
      } else if (is.list(x)) {
        # Handle regular lists
        return(lapply(x, function(elem) replace_and_eval(elem, env)))
      } else {
        return(x)
      }
    }

    # Process the argument list
    arg <- replace_and_eval(arg, eval_env)
  }

  # Add data to the environment (after replacement)
  if ("data" %in% names(arg)) {
    data_obj <- arg$data
    for (col_name in names(data_obj)) {
      eval_env[[col_name]] <- data_obj[[col_name]]
    }
  }

  # Call the function with explicit environment
  tryCatch(
    {
      # Make a copy of the arguments
      arg_copy <- arg

      # Ensure formulas have the right environment
      if ("fixed" %in% names(arg_copy) && inherits(arg_copy$fixed, "formula")) {
        environment(arg_copy$fixed) <- eval_env
      }

      if ("random" %in% names(arg_copy)) {
        if (inherits(arg_copy$random, "formula")) {
          environment(arg_copy$random) <- eval_env
        } else if (is.list(arg_copy$random)) {
          for (i in seq_along(arg_copy$random)) {
            if (inherits(arg_copy$random[[i]], "formula")) {
              environment(arg_copy$random[[i]]) <- eval_env
            }
          }
        }
      }

      fittedmodel <- do.call(fun, arg_copy)
      return(fittedmodel)
    },
    error = function(e) {
      # If that fails, try using do.call with the environment parameter
      tryCatch(
        {
          fittedmodel <- do.call(fun, arg, envir = eval_env)
          return(fittedmodel)
        },
        error = function(e2) {
          # If that still fails, try a direct call to lme
          if (
            identical(fun, nlme::lme) ||
            (is.character(fun) && fun == "nlme::lme")
          ) {
            # Extract common parameters
            fixed_form <- arg$fixed
            random_form <- arg$random
            data_obj <- arg$data

            # Try with a minimal set of arguments
            return(nlme::lme(
              fixed = fixed_form,
              random = random_form,
              data = data_obj
            ))
          } else {
            # If all else fails, re-throw the original error
            stop(e)
          }
        }
      )
    }
  )
}

#' Transforma and merge data and fit models. The function is used inside
#' seq_wrapper to combine metadata with target-level data and perform the
#' model fitting.
#' The function is used in a call to pbapply::pblapply.
#' @param x A data frame of target-specific quantities, in seqwrap a list
#' from create_list is used iteratively.
#' @param samp_name Sample names from the upper level function
#' @param metdat Metadata from the upper level function
#' @param add_vars Additional variables to keep from the metadata
#' @param arg_list Arguments from the upper level function
#' @param mt_summary_fun Summary function from the upper level function
#' @param mt_eval_fun Evaluation function from the upper level function
#' @param return_mod Logical, should the models be returned as part of the
#'   results?
#' @param save_mods Logical, should the models be saved?
#' @param mod_path Path to save the models
#' @param target_name The target identifier, used to name the saved model file.
#'   Falls back to `names(x)` when not supplied.
#' @param is_lme Logical, is `ffun` `nlme::lme`? Computed from `ffun` when not
#'   supplied. Passed in by `seqwrap()` so it is determined once per run rather
#'   than once per target.
#' @param ffun the fitting function from the upper level function
#' @return A list with the fitted model (if `return_mod = TRUE`), summary
#'   and evaluation outputs, and any errors/warnings captured during fitting.
#' @importFrom stats as.formula
#' @keywords internal
#' @noRd
seqwrap_mtf <- function(
  x,
  samp_name,
  metdat,
  arg_list,
  add_vars,
  mt_summary_fun,
  mt_eval_fun,
  ffun,
  return_mod,
  save_mods,
  mod_path,
  target_name = NULL,
  is_lme = NULL
) {
  # Extracting the specific target-specific data and transposing
  transposed <- data.frame(t(x[[1]]))
  transposed$temp <- rownames(transposed)
  colnames(transposed)[ncol(transposed)] <- samp_name
  rownames(transposed) <- NULL

  # Merging target-specific data with meta data
  df <- merge(transposed, data.frame(metdat), by = samp_name)

  # Keep only data needed for fitting
  # Improved formula variable extraction for lme compatibility
  parsed <- character(0)

  # Handle standard formula parameter (lm, glm, etc.)
  if ("formula" %in% names(arg_list)) {
    parsed <- c(parsed, all.vars(as.formula(arg_list$formula)))
  }

  # Handle model parameter (various functions)
  if ("model" %in% names(arg_list)) {
    parsed <- c(parsed, all.vars(as.formula(arg_list$model)))
  }

  # Handle fixed/random parameters (for lme, lmer, etc.)
  if ("fixed" %in% names(arg_list)) {
    parsed <- c(parsed, all.vars(as.formula(arg_list$fixed)))
  }

  # Specialized handling for nlme::lme random effects
  if ("random" %in% names(arg_list)) {
    # lme can have random as a formula
    if (inherits(arg_list$random, "formula")) {
      parsed <- c(parsed, all.vars(arg_list$random))
    } # lme can have random as a list of formulas
    else if (is.list(arg_list$random)) {
      for (i in seq_along(arg_list$random)) {
        if (inherits(arg_list$random[[i]], "formula")) {
          parsed <- c(parsed, all.vars(arg_list$random[[i]]))
        } else if (is.call(arg_list$random[[i]])) {
          # Extract variable names from calls (e.g., pdDiag(~time))
          parsed <- c(parsed, all.vars(arg_list$random[[i]]))
        }
      }
    } # Handle the special case where random is a call (typically for nlme)
    else if (is.call(arg_list$random)) {
      parsed <- c(parsed, all.vars(arg_list$random))
    }
  }

  # Remove duplicates
  parsed <- unique(parsed)

  # Filter for columns that actually exist in the data
  parsed <- parsed[parsed %in% colnames(df)]

  # Keep also additional variables that exists in the meta data data set
  if (!is.null(add_vars)) parsed <- c(parsed, add_vars)

  df <- df[, parsed, drop = FALSE]

  # Remove attributes from the list of arguments
  # (this solves an issue when using glmmTMB)
  # Create final arguments list, handling special case for lme
  arguments_final <- append(arg_list, list(data = df))

  # Determine if we're using lme from nlme. The answer depends only on the
  # fitting function, which is fixed for the whole run, so seqwrap() computes it
  # once and passes it in; deparsing the fitting function for every target is
  # pure waste. The fallback keeps this function usable on its own.
  if (is.null(is_lme)) {
    is_lme <- sw_is_lme(ffun)
  }

  # Handle formula environments to ensure proper evaluation
  for (i in seq_along(arguments_final)) {
    if (inherits(arguments_final[[i]], "formula")) {
      # Set formula environment to current environment for all functions
      # This is needed for glmmTMB with complex random effects (e.g., ||)
      # and for lme compatibility
      environment(arguments_final[[i]]) <- environment()
    } else if (is.list(arguments_final[[i]])) {
      # Handle nested formulas in lists (like in random effects)
      for (j in seq_along(arguments_final[[i]])) {
        if (inherits(arguments_final[[i]][[j]], "formula")) {
          environment(arguments_final[[i]][[j]]) <- environment()
        }
      }
    }
  }

  # Making target-wise data available
  # for the model fitting algorithm.
  if (!is.null(x[[2]])) {
    target.wise <- x[[2]]
    # NEW: Convert single-row data frame to named list
    if (is.data.frame(target.wise) && nrow(target.wise) == 1) {
      target.wise <- as.list(target.wise)
    }


  } else {
    target.wise <- NULL
  }
  # Add warning/errors to outputs
  warn <- NULL
  err <- NULL
  warn_sum <- NULL
  warn_eval <- NULL
  err_sum <- NULL
  err_eval <- NULL

  # Adding null values to model outputs
  mod <- NULL
  # Adding null values to outputs from summaries and evaluations
  mod_sum <- NULL
  mod_eval <- NULL

  ## Fit the model, recording conditions without discarding the result.
  ##
  ## A warning handler passed to tryCatch() unwinds the stack: the handler runs
  ## and tryCatch() returns, so the expression never finishes. Model fitting
  ## routinely warns while still producing a usable fit (glmmTMB on a
  ## convergence warning, lme4 on a singular fit, glm.nb when theta hits the
  ## iteration limit), and those fits were being thrown away. withCallingHandlers
  ## records the warning and resumes, so only genuine errors lose the model.
  mod <- tryCatch(
    withCallingHandlers(
      {
        if (is_lme) {
          # For lme, use a more cautious approach
          fit_fun_lme(ffun, arguments_final, target.wise)
        } else {
          # For other functions, use the original approach
          fit_fun(ffun, arguments_final, vars = target.wise)
        }
      },
      warning = function(w) {
        # Keep the first warning, matching what the previous tryCatch() saw
        if (is.null(warn)) warn <<- w
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      err <<- e
      NULL
    }
  )

  ## Do summarize function if it exists. A failed fit has nothing to summarise,
  ## and passing NULL to a summary function tends to produce a degenerate
  ## zero-column result rather than an honest NULL.
  if (!is.null(mt_summary_fun) && !is.null(mod)) {
    mod_sum <- tryCatch(
      withCallingHandlers(
        do.call("mt_summary_fun", list(mod)),
        warning = function(w) {
          if (is.null(warn_sum)) warn_sum <<- w
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        err_sum <<- e
        NULL
      }
    )
  }

  ## Do evaluation function if it exists
  if (!is.null(mt_eval_fun) && !is.null(mod)) {
    mod_eval <- tryCatch(
      withCallingHandlers(
        do.call("mt_eval_fun", list(mod)),
        warning = function(w) {
          if (is.null(warn_eval)) warn_eval <<- w
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        err_eval <<- e
        NULL
      }
    )
  }

  # Save the model if requested. The target identifier has to be passed in:
  # lapply()/Map() hand over list elements without their names, so names(x) is
  # NULL here and every model would otherwise be written to the same file.
  if (save_mods) {
    id <- if (!is.null(target_name)) target_name else names(x)
    saveRDS(
      mod,
      file = file.path(mod_path, paste0(sw_safe_filename(id), ".rds"))
    )
  }

  # Return the model if requested
  if (return_mod) {
    return(list(
      model = mod,
      summaries = mod_sum,
      evaluation = mod_eval,
      warn = warn,
      err = err,
      warn_sum = warn_sum,
      warn_eval = warn_eval,
      err_sum = err_sum,
      err_eval = err_eval
    ))
  } else {
    return(list(
      summaries = mod_sum,
      evaluation = mod_eval,
      warn = warn,
      err = err,
      warn_sum = warn_sum,
      warn_eval = warn_eval,
      err_sum = err_sum,
      err_eval = err_eval
    ))
  }
}



