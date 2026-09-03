# Defines custom property types for use in class constructors #

# Define a modelfun property that correctly updates model_print
modelfun_prop <- S7::new_property(
  validator = function(value) {
    if (!(is.null(value) || is.function(value))) "must be NULL or a function"
  }
)

# Accepts either NULL or a list, this will also validate data frames.
null_or_list <- S7::new_property(
  # The validator function checks if the value is NULL or a list
  validator = function(value) {
    if (!(is.null(value) || is.list(value))) "must be NULL or a list"
  }
)

# Define a custom property type that accepts either NULL or a list
null_or_df <- S7::new_property(
  # The validator function checks if the value is NULL or a list
  validator = function(value) {
    if (!(is.null(value) || is.data.frame(value)))
      "must be NULL or a data.frame"
  }
)


# Define a custom property type that accepts either data frame or a list
df_or_list <- S7::new_property(
  validator = function(value) {
    if (!(is.null(value) || is.data.frame(value) || is.list(value))) {
      "must be NULL, a data frame, or a list"
    }
  },
  default = NULL # Change default to NULL
)

# Define a call property for use in the swcontainer function
call_prop <- S7::new_property(
  validator = function(value) {
    if (!(is.call(value) || is.null(value))) {
      return("must be a language object representation of a call")
    }
    return(NULL) # Explicitly return NULL for valid values
  }
)

# Define a custom property type that accepts either NULL or function
null_or_function <- S7::new_property(
  validator = function(value) {
    if (!(is.null(value) || is.function(value))) "must be NULL or a function"
  }
)


#' seqwrapResults class
#'
#' seqwrapResults constructor function.
#'
#' @param .data A list containing initial data (inherited from S7::class_list)
#' @param models A list of fitted model objects
#' @param summaries A list of model summaries
#' @param evaluations A list of model evaluations
#' @param errors A data frame of errors and warnings with one row per condition
#' raised, holding the columns `target`, `stage`, `type`, `class` and `message`.
#' Targets that raised nothing do not appear. Defaults to an empty data frame.
#' @param targets A character vector of every target identifier, in the order
#' the targets were fitted
#' @param cache A list describing an on-disk cache of summaries and evaluations,
#' or NULL when these are held in memory. See `seqwrap()`.
#' @param n Number of samples
#' @param k Number of targets
#' @param call_arguments Character string of function arguments used
#' @param call_engine Character string of modeling engine used
#' @param elapsed_time An proc.time() object giving the time needed to complete
#' iterative model fitting.
#'
#' @return An S7 object of class `seqwrap_results` storing fitted models,
#' summaries, evaluations, and diagnostic information from a `seqwrap()` run.
#'
#' @examples
#' # seqwrapResults is the S7 class returned by seqwrap(). End users do
#' # not normally call the constructor directly -- it is invoked
#' # internally once iterative model fitting has finished.
#'
#' library(seqwrap)
#'
#' dat <- simcounts(n_genes = 5)
#'
#' container <- seqwrap_compose(
#'   modelfun = stats::lm,
#'   arguments = list(formula = y ~ x),
#'   data = dat$data,
#'   metadata = dat$metadata,
#'   samplename = "sample"
#' )
#'
#' results <- seqwrap(container, subset = 1:2, cores = 1)
#'
#' # The object returned by seqwrap() is a seqwrapResults instance
#' S7::S7_inherits(results, seqwrapResults)
#'
#' @export
seqwrapResults <- S7::new_class(
  name = "seqwrap_results",
  parent = S7::class_list,
  properties = list(
    models = null_or_list,
    summaries = null_or_list,
    evaluations = null_or_list,
    errors = S7::new_property(
      S7::class_data.frame,
      default = quote(data.frame())
    ),
    targets = S7::new_property(S7::class_character, default = character(0)),
    cache = null_or_list,
    n = S7::new_property(S7::class_numeric, default = integer(0)),
    k = S7::new_property(S7::class_numeric, default = integer(0)),
    call_arguments = S7::new_property(
      S7::class_character,
      default = character(0)
    ),
    call_engine = S7::new_property(S7::class_character, default = character(0)),
    elapsed_time = S7::new_property(S7::class_numeric, default = numeric(0))
  )
)


#' seqwrap container class
#'
#' The seqwrap container is used to store and validate data input
#' to the to seqwrap.
#'
#' @param .data A list containing initial data (inherited from S7::class_list)
#' @param modelfun A function used to fit models
#' @param arguments A list of arguments for the fitting function
#' @param data A data frame or list of data frames with target data
#' @param rownames Logical, should row names be used as target IDs?
#' @param metadata A data frame with sample information
#' @param targetdata A data frame with target-wise information
#' @param samplename Character for sample name identification
#' @param additional_vars Character vector of additional variables
#' @param summary_fun A function for summarizing models
#' @param eval_fun A function for evaluating models
#' @param exported A list of objects to export to workers
#' @param model_print Character representation of model function
#' @param arguments_print Character representation of arguments
#'
#' @return An S7 object of class `swcontainer` bundling data, metadata,
#' and the modelling function plus arguments to be consumed by `seqwrap()`.
#'
#' @examples
#' # swcontainer is the S7 class produced by seqwrap_compose(). End users
#' # normally build one via seqwrap_compose() rather than calling this
#' # constructor directly.
#'
#' library(seqwrap)
#'
#' dat <- simcounts(n_genes = 5)
#'
#' container <- seqwrap_compose(
#'   modelfun = stats::lm,
#'   arguments = list(formula = y ~ x),
#'   data = dat$data,
#'   metadata = dat$metadata,
#'   samplename = "sample"
#' )
#'
#' # The object returned by seqwrap_compose() is a swcontainer instance
#' S7::S7_inherits(container, swcontainer)
#'
#' @export
swcontainer <- S7::new_class(
  name = "swcontainer",
  parent = S7::class_list,
  properties = list(
    modelfun = modelfun_prop,
    arguments = null_or_list,
    data = null_or_list,
    rownames = S7::new_property(S7::class_logical, default = logical(0)),
    metadata = S7::new_property(
      S7::class_data.frame,
      default = quote(data.frame())
    ),
    targetdata = df_or_list,
    samplename = S7::new_property(S7::class_character, default = character(0)),
    additional_vars = S7::new_property(
      S7::class_character,
      default = character(0)
    ),
    summary_fun = null_or_function,
    eval_fun = null_or_function,
    exported = S7::new_property(S7::class_list, default = list()),
    model_print = S7::new_property(S7::class_character, default = character(0)),
    arguments_print = S7::new_property(
      S7::class_character,
      default = character(0)
    )
  )
)


#' Print method for objects of class seqwrapResults
#'
#' @description
#' Invoking the print method on seqwrapResults gives a summary
#' of the fitted objects.
#'
#' @details
#' The method is registered for the S7 class `seqwrapResults` and is called
#' as `print(x, ...)`, where `x` is the object returned by [seqwrap()].
#' Further arguments (`...`) are accepted for compatibility with the
#' [print()] generic but are not used.
#'
#' @return Invisibly returns the object
#'
#' @usage NULL
#'
#' @examples
#' # Load packages and prepare data for examples --------------------------------
#'
#' library(seqwrap)
#'
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#' library(glmmTMB)
#'
#' dat <- simcounts2()
#'
#' # Save simulated data as separate objects
#' counts <- dat$counts
#' metadata <- dat$metadata
#'
#' # Prepare library size for use as offset
#' metadata$ln_libsize <- log(metadata$library_size)
#'
#'
#' # A mixed effects negative binomial model of RNA-seq counts ------------------
#'
#' # Populate the seqwrap container
#' container <- seqwrap_compose(
#'   modelfun = glmmTMB::glmmTMB,
#'   arguments = list(
#'     formula = y ~ time * condition + (1|id) + offset(ln_libsize),
#'     family = glmmTMB::nbinom2()
#'   ),
#'   data = counts,
#'   metadata = metadata,
#'   samplename = "seq_sample_id"
#' )
#'
#' # Run seqwrap using the container
#' results <- seqwrap(container,
#'                    cores = 1)
#' print(results)
#' }
#'
#' @method print seqwrapResults
#' @name print.seqwrapResults
S7::method(print, seqwrapResults) <- function(x, ...) {
  cli::cli_h1("seqwrap")
  cli::cli_inform(
    "A total of {x@n} sample{?s} and {x@k} target{?s} where
                  used in {.code {x@call_engine}} with arguments
                  {.code {x@call_arguments}}"
  )

  if (!is.null(x@cache)) {
    cli::cli_alert_info(
      "Summaries and evaluations are cached on disk in
      {.path {x@cache$path}} ({x@cache$n_chunks} chunk file{?s}).
      Use {.code seqwrap_summarise()} to combine them."
    )
  }

  if (NROW(x@errors) > 0L) {
    cli::cli_alert_info("Some targets had associated errors or warnings")
    cli::cli_inform(sw_condition_bullets(x@errors, x@k))
  } else cli::cli_alert_info("No targets had associated errors or warnings")

  invisible(x)
}


#' Format per stage condition counts for printing
#'
#' @param errors A long form condition table.
#' @param k The number of targets, used for percentages.
#' @return A named character vector of cli bullets.
#' @keywords internal
#' @noRd
sw_condition_bullets <- function(errors, k) {
  labels <- c(
    fit = "Fitting algorithm",
    summary = "Summary function",
    evaluation = "Evaluation function"
  )

  bullets <- character(0)
  for (stage in .sw_stages) {
    for (type in .sw_types) {
      n <- sum(errors$stage == stage & errors$type == type)
      percent <- if (length(k) == 1L && k > 0) round(100 * n / k) else 0
      bullets <- c(
        bullets,
        sprintf("%s (%ss): n = %d (%d%%)", labels[[stage]], type, n, percent)
      )
    }
  }

  stats::setNames(bullets, rep("*", length(bullets)))
}


#' Compose a swcontainer object for use in the seqwrap function.
#'
#' This function makes it possible to compose and  run checks on combined
#' data sets (meta data and target data) and fitting functions to avoid issues
#' in iterative modelling. See examples and vignettes for details.
#'
#' @param x An optional swcontainer object to modify. When supplied, the
#' properties named in `update` are changed and the container is returned; all
#' other arguments are ignored.
#' @param modelfun A model fitting function like stats::lm,
#' glmmTMB::glmmTMB or lme4::lmer
#' @param arguments An alist or list of arguments to be passed to the fitting
#'  function, this should not contain data. Note that the formula must have
#'  y as the dependent variable.
#' @param data A data frame or a list of data frames with targets (e.g. genes,
#' transcripts) as rows and sample names as columns.
#' If rownames = FALSE (default), each data frame should have target
#' identifications as the first column in the data frame(s). If rownames = TRUE
#' row names will be converted to target identifications. If data is provided as
#' a list, each element of the list should be named. The corresponding names
#' be available as variables for the fitting function.
#' @param rownames should row names in data be used as target identifications?
#' Defaults to FALSE.
#' @param metadata A data frame with sample names (corresponding to column
#' names in the target matrix)
#' and design variables.
#' @param targetdata A data frame or a list with target-wise values
#' (e.g. dispersion or start values) for each target. This data is made
#' available for the model fitting function and can be used to specify target
#' specific data in each iteration of seqwrap. When a data frame is provided
#' each row corresponds to the target specific value and each column is
#' available by name. When a list is provided, each element is a data frame
#' whose columns are available by name. A list whose elements are named by
#' target identifier is matched to targets by name and may hold more targets
#' than `data`. A data frame, or a list without names, is matched by position,
#' row or element `i` belonging to row `i` of `data`; `seqwrap()` warns when a
#' list is matched this way. See `seqwrap_priors()` for building
#' target-specific priors for `glmmTMB`.
#' @param samplename A character value indicating the variable by which
#' metadata can merge with the target data. This defaults to "seq_sample_id"
#' as this is used in the trainomeMetaData package.
#' @param additional_vars A vector of additional variables that is contained
#' in the metadata data set that is needed to fit the model. By default the
#' metadata is reduced to variables contained in the slots
#' formula/model/fixed/random in additional arguments.
#' More variables may be needed for offsets, weights etc.
#' @param summary_fun A custom (user-created) function for
#' evaluating/summarizing models. If NULL, `generic_summary()` is used, which
#' tidies model parameters with `broom.mixed::tidy()`.
#' @param eval_fun A custom (user-created) function for model
#' diagnostics/evaluation. If NULL, `generic_evaluation()` is used, which
#' reports the convergence diagnostics supplied by the fitting algorithm. Pass
#' `residual_diagnostics` for DHARMa based residual checks instead, bearing in
#' mind that these are far more expensive per target.
#' @param exported A list of functions, values etc. to be passed to
#' summary_fun and eval_fun. This list must contain any functions that
#' should be used in model summarise or evaluations.
#' @param update A list of named parameters to update a swcontainer object.
#' @return A swcontainer object for direct use in seqwrap.
#' @examples
#' # Load packages and prepare data for examples -------------------------------
#'
#' library(seqwrap)
#'
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#' library(glmmTMB)
#'
#' dat <- simcounts2()
#'
#' # Save simulated data as separate objects
#' counts <- dat$counts
#' metadata <- dat$metadata
#'
#' # Prepare library size for use as offset
#' metadata$ln_libsize <- log(metadata$library_size)
#'
#'
#' # A mixed effects negative binomial model of RNA-seq counts -----------------
#'
#' # Populate the seqwrap container
#' container <- seqwrap_compose(
#'   modelfun = glmmTMB::glmmTMB,
#'   arguments = list(
#'     formula = y ~ time * condition + (1|id) + offset(ln_libsize),
#'     family = glmmTMB::nbinom2()
#'   ),
#'   data = counts,
#'   metadata = metadata,
#'   samplename = "seq_sample_id"
#' )
#'
#' # Run seqwrap using the container. Models are returned here only so that
#' # they can be compared with the prior-based fits further down; a full run
#' # would normally leave return_models at its default of FALSE.
#' results <- seqwrap(container,
#'                    return_models = TRUE,
#'                    cores = 1)
#'
#'
#'
#'
#' # Summarise results
#' summaries <- seqwrap_summarise(results)
#'
#'
#' # Including target-specific data --------------------------------------------
#'
#' \donttest{
#'
#' # Target specific data can be supplied to seqwrap_compose to enable, e.g.,
#' # the use of priors for empirical Bayes shrinkage. In this example we are
#' # setting a dummy-prior on the `condition` parameter and target-specific
#' # prior on the dispersion parameter.
#'
#' fixef_priors <- data.frame(
#'   prior = "normal(0, 1)",
#'   class = "fixef",
#'   coef = "conditionB"
#' )
#'
#' dips_priors <- data.frame(
#'   prior = c(
#'     "normal(0.5, 0.25)",
#'     "normal(1, 0.25)"),
#'   class = "fixef_disp",
#'   coef = 1
#' )
#'
#'
#' # Combine information in a target-specific list, named by target
#' prior_list <- list()
#' for(i in 1:2) {
#'
#'   prior_list[[counts[i, 1]]] <- rbind(
#'     fixef_priors,
#'     dips_priors[i,]
#'   )
#'
#'
#' }
#'
#'
#'
#' container <- seqwrap_compose(
#'   modelfun = glmmTMB::glmmTMB,
#'   # NOTE: The use of `alist` prevents evaluation of list components
#'   arguments = alist(
#'     formula = y ~ time * condition + (1|id) + offset(ln_libsize),
#'     family = glmmTMB::nbinom2(),
#'     priors = data.frame(
#'       prior = prior,
#'       class = class,
#'       coef = coef
#'     )
#'   ),
#'   data = counts,
#'   metadata = metadata,
#'   targetdata = prior_list,
#'   samplename = "seq_sample_id"
#' )
#'
#' # Run seqwrap using the container
#' results_prior <- seqwrap(container,
#'                          # Return models to confirm use of prior information
#'                          # The use of prior information is not recommended
#'                          return_models = TRUE,
#'                          cores = 1)
#'
#' # Confirm prior information
#' summary(results_prior@models[[1]])
#' # Compare to naive model
#' summary(results@models[[1]])
#' }
#' }
#'
#' @export
seqwrap_compose <- function(
  x = NULL,
  modelfun,
  arguments,
  data,
  rownames = FALSE,
  metadata,
  targetdata = NULL,
  samplename = "seq_sample_id",
  additional_vars = NULL,
  summary_fun = NULL,
  eval_fun = NULL,
  exported = list(),
  update = list()
) {
  # Check if x is a swcontainer for updating
  if (S7::S7_inherits(x, swcontainer)) {
    # Pass updates as a list for the named parameter
    return(seqwrap_update(container = x, update))
  }

  # Checked before any of the modelling arguments are touched, so that a
  # DGEList reports the real problem rather than a missing argument
  if (inherits(x, "DGEList")) {
    cli::cli_abort(c(
      "DGEList input is no longer accepted.",
      "i" = "Extract the counts and sample information and pass them as
             {.arg data} and {.arg metadata}, for example
             {.code data = data.frame(target = rownames(dge), dge$counts)}."
    ))
  }

  # Extract the modelling algorithm for printing
  if (!is.null(modelfun)) {
    call_str <- match.call()
  }

  # Initialize an empty container and add objects if they exist
  container <- swcontainer()

  # Set other properties
  if (!is.null(modelfun)) container@modelfun <- modelfun
  if (!is.null(modelfun)) container@model_print <- deparse(call_str$modelfun)
  if (!is.null(arguments)) container@arguments <- arguments
  if (!is.null(arguments))
    container@arguments_print <- deparse(call_str$arguments)
  if (!is.null(data)) container@data <- data
  if (!is.null(rownames)) container@rownames <- rownames
  if (!is.null(metadata)) container@metadata <- metadata
  if (!is.null(targetdata)) container@targetdata <- targetdata
  if (!is.null(samplename)) container@samplename <- samplename
  if (!is.null(additional_vars)) container@additional_vars <- additional_vars
  if (!is.null(summary_fun)) container@summary_fun <- summary_fun
  if (!is.null(eval_fun)) container@eval_fun <- eval_fun
  if (!is.null(exported)) container@exported <- exported

  sw_validate_targetdata(container)

  # Return the populated container
  return(container)
}


# Define a generic with named parameters only
seqwrap_update <- S7::new_generic("seqwrap_update", "container")

# Implement method for swcontainer class
S7::method(seqwrap_update, swcontainer) <- function(
  container,
  update = list()
) {
  if (!S7::S7_inherits(container, swcontainer)) {
    stop("First argument must be a swcontainer object")
  }

  if (length(update) == 0) {
    return(container) # Nothing to update
  }

  # Get property names
  name <- names(update)

  # Special handling for modelfun to update model_print
  if ("modelfun" %in% name) {
    # The model print information is not updated
    # due to difficulties getting the call from the
    # method. This is a temporary solution.
    update$model_print <- "The model function has been updated"
  }

  # Update the container object
  S7::props(container) <- update

  # Updated data or target-wise data have to agree, as in a fresh container
  if (any(c("data", "targetdata", "rownames") %in% name)) {
    sw_validate_targetdata(container)
  }

  return(container)
}


#' Check that target-wise data agrees with the data
#'
#' A data frame or an unnamed list is paired by position and must have one
#' entry per target; a list named by target identifier is paired by name and
#' must cover every target. Nothing is done when either slot is empty.
#'
#' @param container A swcontainer object.
#' @return TRUE (invisibly), or an error.
#' @keywords internal
#' @noRd
sw_validate_targetdata <- function(container) {
  if (is.null(container@data) || is.null(container@targetdata)) {
    return(invisible(TRUE))
  }

  target_ids <- sw_target_ids(container@data, isTRUE(container@rownames))
  n_data_rows <- length(target_ids)

  if (is.data.frame(container@targetdata)) {
    n_target_rows <- nrow(container@targetdata)
    if (n_data_rows != n_target_rows) {
      cli::cli_abort(
        "targetdata must have the same number of rows ({n_target_rows})
        as data has rows ({n_data_rows})"
      )
    }
  } else {
    sw_align_targetdata(container@targetdata, target_ids, warn = FALSE)
  }

  invisible(TRUE)
}


#' Check swcontainer objects for seqwrap
#'
#' This function performs verbose checks/diagnostics of a swcontainer object.
#'
#' @param x A swcontainer object
#' @param verbose Logical, should the function print diagnostics? Default TRUE
#' @return `TRUE` (invisibly) when the container is consistent; otherwise
#' calls `cli::cli_abort()` with a descriptive message.
#' @keywords internal
#' @noRd
seqwrap_check <- function(x, verbose = TRUE) {
  # Check if the container is a swcontainer object
  if (!S7::S7_inherits(x, swcontainer)) {
    stop("The container must be a swcontainer object")
  }

  # checks if arguments for the provided fitting function matches arguments
  if (!all(names(x@arguments) %in% names(formals(x@modelfun)))) {
    cli::cli_abort(
      "Arguments do not match named arguments of the selected
                   \nmodel fitting function ('modelfun')."
    )
  }

  # Check if the data has a character vector for first column
  # (indicating transcript id).
  if (!is.character(x@data[, 1])) {
    cli::cli_abort(
      "The first column of the data is not character or
            factor,\ncheck if this column indicate target
            identifications."
    )
  }

  # Check if the sample name is present in the meta data
  if (!x@samplename %in% colnames(x@metadata)) {
    cli::cli_abort(
      "The samplename does not exist in the metadata,\nno
            variable for matching metadata and
            data."
    )
  }

  # Check that data is formatted correctly
  if (!all(x@metadata[, x@samplename] %in% names(x@data[, -1]))) {
    cli::cli_abort(
      "The sample names in the metadata does
         not match the sample column names in the seqdata."
    )
  }

  # Count the number of unique targets in the data
  if (x@rownames) {
    n_targets <- length(unique(rownames(x@data)))
  } else {
    n_targets <- length(unique(x@data[, 1]))
  }

  # If the target unique indicator is not the same as number of rows abort
  if (!n_targets == nrow(x@data)) {
    cli::cli_abort(
      "The number of unique target identifications does not
        match the number of rows in the data. You might have duplicate id's."
    )
  }

  # Count the number of samples
  n_samples <- length(unique((x@metadata[, x@samplename])))

  # If the unique sample indicator is not the same as number of rows abort
  if (!n_samples == nrow(x@metadata)) {
    cli::cli_abort(
      "The number of unique sample identifications does not
        match the number of rows in the meta data.
      You might have duplicate sample id's."
    )
  }

  # Check that sample id match between data and meta data
  meta_data_sample_id <- unique((x@metadata[, x@samplename]))

  if (x@rownames) {
    data_sample_id <- colnames(x@data)
  } else {
    data_sample_id <- colnames(x@data[, -1])
  }

  # Save sample id for printing
  if (length(data_sample_id) > 3) {
    sample_id_print <- paste(data_sample_id[1:3], collapse = ", ")
  } else {
    sample_id_print <- paste0(paste(data_sample_id, collapse = ", "), ", ...")
  }

  if (
    !all(meta_data_sample_id %in% data_sample_id) &&
      all(data_sample_id %in% meta_data_sample_id)
  ) {
    cli::cli_abort(
      "The sample names in the meta data does not match the sample column names
    in the data. Check if the sample names are correctly formatted and that
    `samplename` looks for the correct variable in the meta data."
    )
  }

  # If verbose is TRUE, print the diagnostics
  if (verbose) {
    printfun <- function() {
      cli::cli_h1("seqwrap diagnostics")

      cli::cli_inform(
        "The swcontainer object has been checked for
                      consistency with the following diagnostics:\n"
      )

      cli::cli_ul(
        "The number of unique target identifications is {n_targets}."
      )

      cli::cli_ul(
        "The number of unique sample identifications is {n_samples}."
      )

      cli::cli_ul(
        "Sample identifications ({sample_id_print}) match between data
        and meta data."
      )

      cli::cli_ul(
        "The modelling function to be used is {x@model_print}."
      )
    }
    printfun()
  }

  # If the swcontainer is a fully populated data container return TRUE
  return(TRUE)
}


# #' Print method for objects of class swcontainer
# #'
# #' @description
# #' Invoking the print method on swcontainer gives a summary
# #' of data container.
# #'
# #' @param x A swcontainer object
# #'
# #' @return Invisibly returns the object
# #'
# #' @examples
# #' seqwrap.dat <- seqwrap_compose(...)
# #' print(seqwrap.dat)
# #'
# #' @method print swcontainer
# #' @name print.swcontainer
# S7::method(print, swcontainer) <- function(x) {
#   cli::cli_h1("seqwrap data container")
#
#   if (!is.null(x@fitting_fun)) {
#
#   }
#
#
#  # invisible(x)
# }

#' A flexible upper-level wrapper for iterative modelling using any available
#' fitting algorithm
#'
#'
#' @inheritParams seqwrap_compose
#' @param y An swcontainer object created with `seqwrap_compose()`, or NULL to
#' build one from the arguments given here.
#' @param return_models Logical, should models be returned as part of the
#' output? Defaults to FALSE, as retaining every fitted model is the largest
#' memory cost a run can incur. Set it to TRUE while developing a model on a
#' subset of targets, where inspecting the fitted objects is useful. A warning
#' is raised when more than ten models would be retained. To keep models from a
#' full run without holding them in memory, use `save_models` and `model_path`.
#' @param save_models Logical, should models be saved? Models may be saved
#' on disk to save working memory.
#' @param model_path A character. The path to saved models.
#' @param subset A sequence, random samples or integers to indicate which
#' rows to keep in data. This is useful if you want to test the model in a
#' subset of targets. If left to the default (NULL), all rows will be used.
#' @param chunk_size An integer giving the number of targets handed to a worker
#' as a single unit of work. Larger chunks reduce the per-task overhead of
#' sending data to and from workers, which matters when the number of targets is
#' large. If NULL (default) a chunk size is chosen automatically so that each
#' worker receives several chunks. Chunking does not change the results.
#' @param cache One of `"none"`, `"memory"` or `"disk"`, controlling how
#' summaries and evaluations are accumulated. See Details.
#' @param cache_path A character path to the directory used when
#' `cache = "disk"`. If NULL (default) a directory inside the session temporary
#' directory is created and removed when the session ends.
#' @param cores An integer indicating the number of cores to be used in parallel
#'  computations. If NULL, a sequential for loop is used. If "max", all
#'  available cores are used.
#' @param verbose Logical, should the function print diagnostics after checking
#' the data container?
#' @return A nested list with three upper levels slots: models, a list of
#' fitted objects; summaries, a list of summaries created from the summary_fun
#' function; evaluations, a list of diagnostics created from eval_fun.
#' @details This function provides a flexible wrapper to fit, summarize and
#' evaluate statistical models fitted to high dimensional omics-type data.
#' Models are fitted and passed to user defined functions to summarize and
#' evaluate models.
#'
#' ## Caching summaries and evaluations
#'
#' Summary and evaluation functions return one small data frame per target.
#' Holding these as separate objects is convenient but costs roughly ten times
#' more memory than the same rows bound into a single data frame, because the
#' per-data-frame overhead is paid once per target. At high target counts
#' (hundreds of thousands, as in array-scale data) that overhead dominates. The
#' `cache` argument controls the trade-off:
#'
#' * `"memory"` (default) binds summaries and evaluations into a single data
#'   frame per slot while still on the worker, adding a `target` column. This
#'   requires that `summary_fun` and `eval_fun` return data frames.
#' * `"none"` keeps one data frame per target, so `@summaries` and
#'   `@evaluations` are named lists indexed by target. Use this when a summary
#'   or evaluation function returns something that cannot be row-bound, or when
#'   per-target access is more convenient than filtering a combined table.
#' * `"disk"` additionally writes each chunk's bound data frames to
#'   `cache_path` and leaves `@summaries` and `@evaluations` empty. Parent
#'   memory then stays flat regardless of the number of targets. Use
#'   `seqwrap_summarise()` to read and combine the cached chunks, and
#'   `seqwrap_cache_clear()` to remove them.
#'
#' `seqwrap_summarise()` returns the same combined data frames whichever mode
#' was used, so only code reading `@summaries` or `@evaluations` directly is
#' affected by the choice.
#'
#' Caching is independent of `save_models` and `model_path`, which continue to
#' write one file per fitted model.
#'
#' Note that `cache = "disk"` caches summaries and evaluations only. With
#' `return_models = TRUE` the fitted models are still collected in memory, so
#' the two options are normally combined as `return_models = FALSE`.
#' @examples
#' library(seqwrap)
#'
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#' library(glmmTMB)
#'
#' # Simulate n targets
#' dat <- simcounts(n_genes = 10000)
#'
#' # Save simulated data as separate objects
#' counts <- dat$data
#' metadata <- dat$metadata
#'
#' # Prepare library size for use as offset
#' metadata$ln_libsize <- log(colSums(counts[,-1]))
#'
#'
#' # A mixed effects negative binomial model of RNA-seq counts -----------------
#'
#' # Populate the seqwrap container
#' container <- seqwrap_compose(
#'   modelfun = glmmTMB::glmmTMB,
#'   arguments = list(
#'     formula = y ~ x + (1|cluster) + offset(ln_libsize),
#'     family = glmmTMB::nbinom2()
#'   ),
#'   data = counts,
#'   metadata = metadata,
#'   samplename = "sample"
#' )
#'
#' # Run seqwrap using the container on a subset of targets
#' results <- seqwrap(container,
#'                    subset = 1:5,
#'                    cores = 1)
#'
#'
#' # Summarise results only contains the subset
#' summaries <- seqwrap_summarise(results)
#' }
#'
#' @export
seqwrap <- function(
  y = NULL,
  modelfun = NULL,
  arguments = NULL,
  data = NULL,
  metadata = NULL,
  samplename = NULL,
  additional_vars = NULL,
  summary_fun = NULL,
  eval_fun = NULL,
  exported = list(),
  return_models = FALSE,
  save_models = FALSE,
  model_path = NULL,
  subset = NULL,
  chunk_size = NULL,
  cache = c("memory", "none", "disk"),
  cache_path = NULL,
  cores = 1,
  verbose = TRUE
) {
  cache <- match.arg(cache)

  # Elapsed time
  start_time <- proc.time()

  # If the input is a swcontainer object, use the object
  if (S7::S7_inherits(y, seqwrap::swcontainer)) {
    container <- y

    # If variables are to be updated, update the container
    updates <- list()

    if (!is.null(modelfun)) updates$modelfun <- modelfun
    if (!is.null(arguments)) updates$arguments <- arguments
    if (!is.null(data)) updates$data <- data
    if (!is.null(metadata)) updates$metadata <- metadata
    if (!is.null(samplename)) updates$samplename <- samplename
    if (!is.null(additional_vars)) updates$additional_vars <- additional_vars
    if (!is.null(summary_fun)) updates$summary_fun <- summary_fun
    if (!is.null(eval_fun)) updates$eval_fun <- eval_fun
    if (!is.null(exported)) updates$exported <- exported

    # Update the container
    if (length(updates) > 0)
      container <- seqwrap_compose(container, update = updates)
  } else if (is.null(y)) {
    # If the input is NULL, compose a new container
    container <- seqwrap_compose(
      modelfun = modelfun,
      arguments = arguments,
      data = data,
      metadata = metadata,
      samplename = samplename,
      additional_vars = additional_vars,
      summary_fun = summary_fun,
      eval_fun = eval_fun,
      exported = exported
    )
  } else if (inherits(y, "DGEList")) {
    cli::cli_abort(c(
      "DGEList input is no longer accepted.",
      "i" = "Extract the counts and sample information and pass them as
             {.arg data} and {.arg metadata}, for example
             {.code data = data.frame(target = rownames(dge), dge$counts)}."
    ))
  } else {
    stop("The input must be a swcontainer object or NULL")
  }

  # If eval_fun or summary_fun are NULL, supply generic functions
  if (is.null(container@summary_fun)) container@summary_fun <- generic_summary
  if (is.null(container@eval_fun)) container@eval_fun <- generic_evaluation

  # If subset is provided, subset the data
  if (!is.null(subset)) {
    # Subset the data
    container@data <- container@data[subset, , drop = FALSE]

    # Subset the targetdata if it is not NULL. A list named by target
    # identifier is aligned by name below, after the data has been subset.
    if (!is.null(container@targetdata)) {
      if (is.data.frame(container@targetdata)) {
        container@targetdata <- container@targetdata[subset, , drop = FALSE]
      } else if (sw_is_unnamed(container@targetdata)) {
        container@targetdata <- container@targetdata[subset]
      }
    }
  }

  # Validate model_path when save_models is TRUE
  if (save_models && is.null(model_path)) {
    stop("'model_path' must be specified when 'save_models' is TRUE")
  }

  # Check the container for consistency
  # seqwrap_check(container, verbose = verbose)

  # Validate cache_path when cache is "disk"
  if (cache == "disk") {
    if (is.null(cache_path)) {
      # CRAN policy allows packages to write to the session temporary
      # directory only, unless the user names a location explicitly.
      cache_path <- tempfile("seqwrap-cache-")
    }
    dir.create(cache_path, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(cache_path)) {
      cli::cli_abort("Could not create the cache directory {.path {cache_path}}")
    }
  }

  # Get the number of targets and samples
  k <- if (is.data.frame(container@data)) {
    nrow(container@data)
  } else {
    nrow(container@data[[1]])
  }
  n <- nrow(container@metadata)

  # Target identifiers name every result, and are used for model file names, so
  # a missing or blank one fails deep inside a worker with an opaque message.
  target_ids <- if (isTRUE(container@rownames)) {
    rownames(container@data)
  } else if (is.data.frame(container@data)) {
    container@data[[1]]
  } else {
    container@data[[1]][[1]]
  }
  bad_ids <- is.na(target_ids) | !nzchar(trimws(as.character(target_ids)))
  if (any(bad_ids)) {
    # Bind each quantity to its own name: a cli string carrying two quantities
    # and a plural marker cannot tell which one to agree with.
    n_bad <- sum(bad_ids)
    n_total <- length(target_ids)
    first_rows <- utils::head(which(bad_ids), 5)
    cli::cli_abort(c(
      "Target identifiers must not be missing or blank.",
      "x" = "{n_bad} of {n_total} targets have a missing or empty identifier.",
      "i" = "Affected row{?s}: {first_rows}."
    ))
  }

  # Pair target-wise data with targets: by name when the list is named by
  # target identifier, otherwise by position with a warning
  container@targetdata <- sw_align_targetdata(
    container@targetdata,
    as.character(target_ids)
  )

  # Retaining fitted models is the largest memory cost a run can incur, so
  # flag it whenever more than a handful would be kept. Prototyping on a small
  # subset, which is what return_models is meant for, stays quiet.
  if (return_models && k > 10) {
    cli::cli_warn(
      "You are saving {k} model objects which could be very memory intensive.
      If this was intended, ignore this warning."
    )
  }

  # Determine the number of cores
  if (is.null(cores)) {
    num_cores <- 1L
  } else if (is.character(cores) && identical(cores, "max")) {
    num_cores <- parallel::detectCores()
  } else if (is.numeric(cores) && length(cores) == 1L) {
    num_cores <- min(as.integer(cores), parallel::detectCores())
  } else {
    stop("'cores' must be NULL, a single integer, or the string \"max\"")
  }

  # Build chunks of targets. Each chunk carries a slice of the raw target data
  # and is expanded into per-target data frames on the worker, so the full
  # per-target expansion never exists in this process.
  if (is.null(chunk_size)) {
    chunk_size <- sw_default_chunk_size(k, num_cores)
  } else if (!is.numeric(chunk_size) || length(chunk_size) != 1L) {
    stop("'chunk_size' must be NULL or a single integer")
  } else if (chunk_size < 1) {
    stop("'chunk_size' must be a positive integer")
  }

  # data_helper() orders targets by identifier for data frame input and by row
  # for list input. Chunking follows that order so that a chunked run returns
  # targets in the same sequence as an unchunked one.
  use_rownames <- isTRUE(container@rownames)
  target_order <- sw_target_order(container@data, rownames = use_rownames)
  chunk_positions <- sw_chunk_index(k, chunk_size)

  chunks <- lapply(seq_along(chunk_positions), function(i) {
    pos <- chunk_positions[[i]]
    list(
      id = i,
      data = sw_slice_data(container@data, target_order[pos]),
      # Target-wise data is paired with targets by their row in the data, so it
      # follows the same reordering as the data itself.
      targetdata = sw_slice_targetdata(container@targetdata, target_order[pos])
    )
  })

  # Catch the function calls for printing
  funcall <- match.call()
  ## Arguments to string
  call_arguments <- deparse(funcall$arguments)
  call_arguments <- sub("^.*?\\((.*)\\).*$", "\\1", call_arguments)
  ## modelfun to string
  call_engine <- deparse(funcall$modelfun)

  # Print pre-fit information

  if (verbose) {
    cli::cli_h1("seqwrap")
    cli::cli_inform(
      "Initiating clusters for parallel processing with {num_cores} core{?s}"
    )
    cli::cli_inform(
      "Fitting {k} target{?s} in {length(chunks)} chunk{?s}
      of up to {chunk_size} target{?s}"
    )
    # The progress bar advances once per batch of chunks, one chunk per core,
    # because workers cannot report progress from inside a chunk. Say so, since
    # a bar that moves in a few large steps otherwise looks stalled.
    n_updates <- max(1L, ceiling(length(chunks) / max(1L, num_cores)))
    cli::cli_inform(
      "The progress bar updates {n_updates} time{?s}
      (once per {num_cores} chunk{?s} completed), in steps of about
      {round(100 / n_updates)}%; use a smaller {.arg chunk_size}
      for more frequent updates"
    )
    if (cache == "disk") {
      cli::cli_inform("Caching summaries and evaluations in {.path {cache_path}}")
      if (return_models) {
        cli::cli_alert_warning(
          "{.code return_models = TRUE} keeps every fitted model in memory;
          consider {.code return_models = FALSE} when caching to disk."
        )
      }
    }
  }

  ## Applying the model function in parallel ##

  # Making variables from the container available in the parallel environment
  metadata <- container@metadata
  arguments <- container@arguments
  modelfun <- container@modelfun
  samplename <- container@samplename
  additional_vars <- container@additional_vars
  summary_fun <- container@summary_fun
  eval_fun <- container@eval_fun
  exported <- container@exported

  # Create a cluster using the number of cores specified
  cl <- parallel::makeCluster(num_cores)

  # Load the seqwrap namespace (and its Imports) on each worker
  parallel::clusterEvalQ(cl, requireNamespace("seqwrap"))

  ## Export data to clusters
  parallel::clusterExport(
    cl,
    varlist = c(
      "metadata",
      "arguments",
      "modelfun",
      "samplename",
      "additional_vars",
      "save_models",
      "exported",
      "model_path",
      "return_models",
      "summary_fun",
      "eval_fun"
    ),
    envir = environment()
  )

  # Place the contents of `exported` on the workers, not just the list itself.
  # Summary and evaluation functions refer to these objects by name, and a
  # function defined in the global environment loses that environment when it is
  # sent to a worker, so the objects have to exist there independently.
  if (length(exported) > 0) {
    if (is.null(names(exported)) || any(!nzchar(names(exported)))) {
      cli::cli_abort("Every element of {.arg exported} must be named.")
    }
    parallel::clusterExport(
      cl,
      varlist = names(exported),
      envir = list2env(exported)
    )
  }

  # Parallel execution of the fitting process
  if (verbose) {
    cli::cli_inform("Merging and modelling data")
  }

  chunk_results <- pbapply::pblapply(
    cl = cl,
    X = chunks,
    FUN = seqwrap_chunk,
    samp_name = samplename,
    metdat = metadata,
    arg_list = arguments,
    add_vars = additional_vars,
    mt_summary_fun = summary_fun,
    mt_eval_fun = eval_fun,
    ffun = modelfun,
    return_mod = return_models,
    save_mods = save_models,
    mod_path = model_path,
    reduce = cache != "none",
    cache_path = if (cache == "disk") cache_path else NULL,
    # Depends only on the fitting function, so resolve it once for the run
    # rather than deparsing the fitting function for every target
    is_lme = sw_is_lme(modelfun),
    rownames = use_rownames
  )

  parallel::stopCluster(cl)

  # Combine results

  models <- NULL
  summaries <- NULL
  evaluations <- NULL
  errors <- NULL
  targets <- character(0)
  cache_manifest <- NULL

  if (cache == "none") {
    # Flatten the chunks back into one entry per target. The result is
    # identical to mapping over targets individually.
    results <- unlist(
      lapply(chunk_results, `[[`, "fits"),
      recursive = FALSE,
      use.names = TRUE
    )

    # Collect models
    if (return_models) models <- lapply(results, `[[`, "model")

    if (!is.null(summary_fun)) summaries <- lapply(results, `[[`, "summaries")
    if (!is.null(eval_fun)) evaluations <- lapply(results, `[[`, "evaluation")

    # names() is NULL for an empty result set, e.g. subset = integer(0)
    targets <- names(results)
    if (is.null(targets)) targets <- character(0)

    ## One row per error or warning actually raised
    errors <- sw_conditions_long(results, targets)
  } else {
    if (return_models) {
      models <- unlist(
        lapply(chunk_results, `[[`, "models"),
        recursive = FALSE,
        use.names = TRUE
      )
    }

    # Chunks already bound their summaries and evaluations, so only the much
    # smaller per-chunk data frames are combined here.
    if (cache == "memory") {
      summaries <- sw_read_chunks(chunk_results, "summaries")
      evaluations <- sw_read_chunks(chunk_results, "evaluations")
    }

    targets <- unlist(
      lapply(chunk_results, `[[`, "targets"),
      use.names = FALSE
    )
    if (is.null(targets)) targets <- character(0)

    errors <- do.call(rbind, lapply(chunk_results, `[[`, "errors"))
    if (is.null(errors)) errors <- sw_empty_conditions()

    if (cache == "disk") {
      cache_manifest <- list(
        path = cache_path,
        format = "rds",
        n_chunks = length(chunk_results),
        chunks = data.frame(
          chunk = seq_along(chunk_results),
          n = vapply(
            chunk_results,
            function(z) length(z$targets),
            integer(1)
          ),
          summaries = vapply(
            chunk_results,
            `[[`,
            character(1),
            "file_summaries"
          ),
          evaluations = vapply(
            chunk_results,
            `[[`,
            character(1),
            "file_evaluations"
          ),
          stringsAsFactors = FALSE
        )
      )
    }
  }

  # Track the elapsed time
  elapsed_time <- proc.time() - start_time
  sys_time <- Sys.time()

  ## Evaluate errors for the resulting print function
  if (verbose) {
    cli::cli_inform({
      paste0(
        format(sys_time, "%b %d %Y %X"),
        ": Completed model fitting and evaluation. Elapsed time was ",
        round(elapsed_time[[3]] / 60, 2),
        " minutes"
      )
    })

    if (NROW(errors) > 0L) {
      cli::cli_alert_info("Some targets had associated errors or warnings")
      cli::cli_inform(sw_condition_bullets(errors, k))
    }
  }

  ## Combine the results into a seqwrapResults
  comb_results <- seqwrapResults(
    models = models,
    summaries = summaries,
    evaluations = evaluations,
    errors = errors,
    targets = targets,
    cache = cache_manifest,
    n = n,
    k = k,
    call_arguments = container@arguments_print,
    call_engine = container@model_print,
    elapsed_time = elapsed_time
  )

  return(comb_results)
}


#' Summarise seqwrapResults objects
#'
#' @param x A seqwrapResults object
#' @param summaries Logical, should summaries be combined?
#' @param evaluations Logical, should evaluations be combined?
#' @param verbose, Logical should progress be printed? Default TRUE
#' @param errors Logical, should errors and warnings be combined? When TRUE the
#' returned list gains an `errors` element with one row per condition raised.
#' @param drop_warnings Should targets whose fitting, summary or evaluation
#' raised a warning be removed from the combined results? FALSE (the default)
#' keeps every target. TRUE removes targets that warned at any stage. A
#' character vector selects stages, one or more of `"fit"`, `"summary"` and
#' `"evaluation"`. See Details.
#'
#' @details
#' This functions attempts to summarize results from the summary and evaluation
#' functions applied in each iteration during modelling. The function expects
#' that the summary and evaluation functions return data frames.
#'
#' ## Targets that raised warnings
#'
#' `drop_warnings` defaults to FALSE, so a target that produced a warning is
#' still reported. This is deliberate. Warnings from fitting algorithms differ
#' widely in severity: a singular fit from `lme4` means a variance component has
#' reached the boundary of the parameter space, which often leaves inference on
#' the fixed effects intact, whereas a non-positive-definite Hessian from
#' `glmmTMB` means the standard errors cannot be trusted. Treating both as
#' grounds for removal discards usable results.
#'
#' Removing them silently is also a form of selection. Targets that warn are not
#' a random subset: low count and low variance targets produce singular fits far
#' more often than others, so dropping them biases the remaining set and changes
#' the number of tests carried into any multiplicity correction.
#'
#' For a graded alternative, `generic_evaluation()` reports `converged` and
#' `singular` per target, which can be joined onto the summaries and filtered
#' explicitly. When `drop_warnings` is used, the identifiers of the removed
#' targets are returned in the `dropped` element so that the choice stays
#' visible and reversible.
#'
#' @return A list (invisibly) with up to four elements: `summaries`, combined
#' parameter summaries from each model; `evaluations`, combined diagnostics from
#' each model; `errors`, one row per error or warning with columns `target`,
#' `stage`, `type` and `message`; and `dropped`, the identifiers of any targets
#' removed by `drop_warnings`. Entries are omitted when the corresponding slot
#' of `x` is empty or the user disables them via the `summaries`,
#' `evaluations` and `errors` arguments.
#' @examples
#' # Load packages and prepare data for examples -------------------------------
#'
#' library(seqwrap)
#'
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#' library(glmmTMB)
#'
#' # Simulate n targets
#' dat <- simcounts(n_genes = 10000)
#'
#' # Save simulated data as separate objects
#' counts <- dat$data
#' metadata <- dat$metadata
#'
#' # Prepare library size for use as offset
#' metadata$ln_libsize <- log(colSums(counts[,-1]))
#'
#'
#' # A mixed effects negative binomial model of RNA-seq counts -----------------
#'
#' # Populate the seqwrap container
#' container <- seqwrap_compose(
#'   modelfun = glmmTMB::glmmTMB,
#'   arguments = list(
#'     formula = y ~ x + (1|cluster) + offset(ln_libsize),
#'     family = glmmTMB::nbinom2()
#'   ),
#'   data = counts,
#'   metadata = metadata,
#'   samplename = "sample"
#' )
#'
#' # Run seqwrap using the container on a subset of targets
#' results <- seqwrap(container,
#'                    subset = 1:5,
#'                    cores = 1)
#'
#'
#' # Summarise results only contains the subset
#' summaries <- seqwrap_summarise(results)
#'
#' # Get summaries
#' summaries$summaries
#'
#' # Get model evaluations
#' summaries$evaluations
#' }
#'
#'
#' @export
seqwrap_summarise <- function(
  x,
  summaries = TRUE,
  evaluations = TRUE,
  errors = TRUE,
  drop_warnings = FALSE,
  verbose = TRUE
) {
  # Check if the input is a seqwrapResults object
  if (!S7::S7_inherits(x, seqwrapResults)) {
    stop("The input must be a seqwrapResults object")
  }

  drop_stages <- sw_drop_stages(drop_warnings)
  dropped_targets <- sw_warned_targets(x@errors, drop_stages)

  # Removes the dropped targets from a combined data frame
  apply_drop <- function(df) {
    if (is.null(df) || length(dropped_targets) == 0L) {
      return(df)
    }
    out <- df[!df[["target"]] %in% dropped_targets, , drop = FALSE]
    rownames(out) <- NULL
    if (nrow(out) == 0L) NULL else out
  }

  ## Print information from the seqwrapResults object
  print_info <- function() {
    cli::cli_h1("seqwrap summarise")

    cli::cli_li(
      "A total of {x@n} sample{?s} and {x@k} target{?s} where
                  used in {.code {x@call_engine}} with arguments
                  {.code {x@call_arguments}}"
    )
    if (summaries) {
      cli::cli_li(
        "Attempting to combine results from
                  the provided summary function."
      )
    }

    if (evaluations) {
      cli::cli_li(
        "Attempting to combine results from
                  the provided evaluations function."
      )
    }

    if (errors) {
      cli::cli_li(
        "Attempting to summarise errors and
                  warnings from the fitting process."
      )
    }
  }

  if (verbose) print_info()

  ## Initialize results variables ##
  summarised_results_final <- NULL
  evaluated_results_final <- NULL

  # print(x)

  ## Extract summarises
  if (summaries) {
    summarised_results_final <- apply_drop(sw_collect(x, "summaries"))

    if (is.null(summarised_results_final)) {
      print_summary <- function() {
        cli::cli_h1("Model summaries")
        cli::cli_alert_info("No summaries available")
      }
      if (verbose) print_summary()
    } else {
      n_summaries <- length(unique(summarised_results_final[["target"]]))

      print_summary <- function() {
        cli::cli_h1("Model summaries")
        cli::cli_li("{n_summaries} targets have associated summaries")
      }

      if (verbose) print_summary()
    }
  } # End if (summaries)

  ## Extract evaluations
  if (evaluations) {
    evaluated_results_final <- apply_drop(sw_collect(x, "evaluations"))

    if (is.null(evaluated_results_final)) {
      print_evaluation <- function() {
        cli::cli_h1("Model evaluations")
        cli::cli_alert_info("No evaluations available")
      }
      if (verbose) print_evaluation()
    } else {
      n_evals <- length(unique(evaluated_results_final[["target"]]))

      print_evaluation <- function() {
        cli::cli_h1("Model evaluations")
        cli::cli_li("{n_evals} targets have associated evaluation results")
      }

      if (verbose) print_evaluation()
    }
  } # End if (evaluations)

  ## Combine results
  # Create empty results list
  results <- list()

  if (!is.null(summarised_results_final)) {
    results$summaries <- summarised_results_final
  }

  if (!is.null(evaluated_results_final)) {
    results$evaluations <- evaluated_results_final
  }

  ## Extract errors and warnings, one row per condition raised
  if (errors) {
    condition_table <- if (NROW(x@errors) > 0L) x@errors else NULL

    if (is.null(condition_table)) {
      if (verbose) {
        cli::cli_h1("Errors and warnings")
        cli::cli_alert_info("No errors or warnings were raised")
      }
    } else {
      results$errors <- condition_table

      if (verbose) {
        n_targets <- length(unique(condition_table$target))
        n_errors <- sum(condition_table$type == "error")
        n_warnings <- sum(condition_table$type == "warning")

        cli::cli_h1("Errors and warnings")
        cli::cli_li(
          "{n_targets} target{?s} raised {n_errors} error{?s} and
          {n_warnings} warning{?s}"
        )
      }
    }
  }

  ## Record which targets were removed, so that the choice stays visible
  if (length(dropped_targets) > 0L) {
    results$dropped <- dropped_targets

    if (verbose) {
      cli::cli_alert_warning(
        "Dropped {length(dropped_targets)} of {x@k} target{?s} that raised a
        warning during {.val {drop_stages}}. Their identifiers are returned in
        the {.field dropped} element."
      )
    }
  }

  # Check if the final list is empty
  if (length(results) == 0) {
    cli::cli_alert_warning("No results were generated, check your input.")
  }

  # Check if the final list is empty
  if (length(results) != 0) {
    if (verbose) {
      cli::cli_h1("Combined results")
      cli::cli_alert_info(
        "Combined results have been generated and were
                        silently returned."
      )
    }
  }

  # Return the final list object
  return(invisible(results))
}
