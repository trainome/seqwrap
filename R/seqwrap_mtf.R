#' A function to fit models with a chosen fitting algorithm, used in seqwrap.
#' @param fun Name of a fitting function, like glmmTMB::glmmTMB
#' @param arg A list of arguments that can be evaluated by the fitting function
#' @param vars A list of variables to be used in the fitting function
#' @keywords internal
fit_fun <- function(fun, arg, vars = NULL) {

  if (!is.null(vars)) {
    # Un-pack the list of variables
    replace_param_value <- function(nested_list, replacement_list) {
      # Define a recursive helper function
      replace_recursive <- function(x) {
        if (is.call(x) || is.expression(x) || is.name(x)) {
          # If it's a name that matches one in our replacement list
          if (is.name(x) && as.character(x) %in% names(replacement_list)) {
            return(replacement_list[[as.character(x)]])
          }

          # If it's a call, process each argument recursively
          if (is.call(x)) {
            for (i in seq_along(x)) {
              x[[i]] <- replace_recursive(x[[i]])
            }
          }
          return(x)
        } else if (is.list(x)) {
          # Process each list element recursively
          for (i in seq_along(x)) {
            x[[i]] <- replace_recursive(x[[i]])
          }
          return(x)
        } else {
          # Return other data types unchanged
          return(x)
        }
      }

      # Call the helper function on the nested list
      return(replace_recursive(nested_list))
    }

    arg <- replace_param_value(arg, vars)

    }


  fittedmodel <- do.call(fun, arg)

  # return the model
  return(fittedmodel)
}

#' Special fitting function for lme
#' @param fun The fitting function, nlme::lme
#' @param arg A list of arguments that can be evaluated by the fitting function
#' @param vars A list of variables to be used in the fitting function
#' @keywords internal
fit_fun_lme <- function(fun, arg, vars = NULL) {
  # For lme models, we need to be more careful with formula environments
  # and how we handle the random effects structure
  
  # Create an environment for the evaluation
  eval_env <- new.env(parent = environment())
  
  # If we have target-specific variables, add them to the evaluation environment
  if (!is.null(vars)) {
    for (var_name in names(vars)) {
      eval_env[[var_name]] <- vars[[var_name]]
    }
  }
  
  # Add data to the environment
  if ("data" %in% names(arg)) {
    data_obj <- arg$data
    for (col_name in names(data_obj)) {
      eval_env[[col_name]] <- data_obj[[col_name]]
    }
  }
  
  # Call the function with explicit environment
  tryCatch({
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
  }, error = function(e) {
    # If that fails, try using do.call with the environment parameter
    tryCatch({
      fittedmodel <- do.call(fun, arg, envir = eval_env)
      return(fittedmodel)
    }, error = function(e2) {
      # If that still fails, try a direct call to lme
      if (identical(fun, nlme::lme) || 
          (is.character(fun) && fun == "nlme::lme")) {
        # Extract common parameters
        fixed_form <- arg$fixed
        random_form <- arg$random
        data_obj <- arg$data
        
        # Try with a minimal set of arguments
        return(nlme::lme(fixed = fixed_form, 
                         random = random_form, 
                         data = data_obj))
      } else {
        # If all else fails, re-throw the original error
        stop(e)
      }
    })
  })
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
#' @return_mod Logical, should the models be returned as part of the results?
#' @param save_mods Logical, should the models be saved?
#' @param mod_path Path to save the models
#' @param ffun the fitting function from the upper level function
#' @importFrom stats as.formula
#' @keywords internal
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
    mod_path
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
    } 
    # lme can have random as a list of formulas
    else if (is.list(arg_list$random)) {
     for (i in seq_along(arg_list$random)) {
        if (inherits(arg_list$random[[i]], "formula")) {
          parsed <- c(parsed, all.vars(arg_list$random[[i]]))
        } else if (is.call(arg_list$random[[i]])) {
         # Extract variable names from calls (e.g., pdDiag(~time))
         parsed <- c(parsed, all.vars(arg_list$random[[i]]))
        }
      }
    }
    # Handle the special case where random is a call (typically for nlme)
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

# Determine if we're using lme from nlme
is_lme <- identical(ffun, nlme::lme) || 
          (is.character(ffun) && ffun == "nlme::lme") ||
          any(grepl("lme$", deparse(ffun)))

# Handle formula environments differently based on function
for (i in seq_along(arguments_final)) {
  if (inherits(arguments_final[[i]], "formula")) {
    if (!is_lme) {
      # For non-lme functions, remove environment as before
      environment(arguments_final[[i]]) <- NULL
    } else {
      # For lme, ensure the environment is set correctly
      environment(arguments_final[[i]]) <- environment()
    }
  } else if (is.list(arguments_final[[i]]) && is_lme) {
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


  ## Fit the model with improved error handling
  tryCatch({
    if (is_lme) {
      # For lme, use a more cautious approach
      mod <- fit_fun_lme(ffun, arguments_final, target.wise)
    } else {
      # For other functions, use the original approach
      mod <- fit_fun(ffun, arguments_final, vars = target.wise)
    }
  }, warning = function(w) warn <<- w,
     error = function(e) err <<- e
  )

  ## Do summarize function if it exists
  if (!is.null(mt_summary_fun)) {
    tryCatch(
      mod_sum <- do.call("mt_summary_fun", list(mod)),
      warning = function(w) warn_sum <<- w,
      error = function(e) err_sum <<- e
    )
  }

  ## Do evaluation function if it exists
  if (!is.null(mt_eval_fun)) {
    tryCatch(
      mod_eval <- do.call("mt_eval_fun", list(mod)),
      warning = function(w) warn_eval <<- w,
      error = function(e) err_eval <<- e
    )
  }

  # Save the model if requested
  mod_path <- if (!is.null(mod_path)) {
    mod_path
  } else {
    paste0(getwd(), "/seqwrap-output")
  }

  if (save_mods) saveRDS(mod, file = paste0(mod_path, "/", names(x), ".rds"))

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
