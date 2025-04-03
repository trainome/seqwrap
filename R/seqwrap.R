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
  # The validator function checks if the value is NULL or a list
  validator = function(value) {
    if (!(is.data.frame(value) || is.list(value))) {
      "must be a data frame or a list"
    }
  }
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







#' seqwrapResults class
#'
#' seqwrapResults constructor function.
#'
#' @slot models description
#' @slot summaries description
#' @slot evaluations description
#' @slot errors description
#' @slot n description
#' @slot k description
#' @slot call_arguments description
#' @slot call_engine description
#'
#'
#' @export
seqwrapResults <- S7::new_class(
  name = "seqwrap_results",
  parent = S7::class_list,
  properties = list(
    models = null_or_list,
    summaries = null_or_list,
    evaluations = null_or_list,
    errors = S7::class_data.frame,
    n = S7::class_numeric,
    k = S7::class_numeric,
    call_arguments = S7::class_character,
    call_engine = S7::class_character
  )
)

#' seqwrap container class
#'
#' The seqwrap container is used to store and validate data input
#' to the to seqwrap.
#'
#' @slot modelfun A function used to fit models using data and meta data
#' as input.
#' @slot arguments A list of arguments provided to the fitting function.
#' @slot data A data frame or a list of data frames
#' @slot rownames A logical indicating if row names should be used as target id
#' @slot metadata A data frame with sample information
#' @slot targetdata A data frame with sample information
#' @slot samplename A character for sample name identification
#' @slot additional_vars A character vector of additional variables exported to
#' @slot summary_fun A function
#' @slot eval_fun A function
#' @slot exported A list
#'
#' @export
swcontainer <- S7::new_class(
  name = "swcontainer",
  parent = S7::class_list,
  properties = list(
    modelfun = modelfun_prop,
    arguments = null_or_list,
    data = null_or_list,
    rownames = S7::class_logical,
    metadata = S7::class_data.frame,
    targetdata = null_or_list,
    samplename = S7::class_character,
    additional_vars = S7::class_character,
    summary_fun = S7::class_function,
    eval_fun = S7::class_function,
    exported = S7::class_list,
    model_print = S7::class_character,
    arguments_print = S7::class_character
  )
)


#' Print method for objects of class seqwrapResults
#'
#' @description
#' Invoking the print method on seqwrapResults gives a summary
#' of the fitted objects.
#'
#' @param x A seqwrapResults object
#'
#' @return Invisibly returns the object
#'
#' @examples
#' results <- seqwrap(...)
#' print(results)
#'
#' @method print seqwrapResults
#' @name print.seqwrapResults
S7::method(print, seqwrapResults) <- function(x) {
  cli::cli_h1("seqwrap")
  cli::cli_inform(
    "A total of {x@n} sample{?s} and {x@k} target{?s} where
                  used in {.code {x@call_engine}} with arguments
                  {.code {x@call_arguments}}"
  )

  # Count non-null values in the error/warning data frame
  errors_sum <- sapply(x@errors, function(x) sum(!sapply(x, is.null)))

  if (any(errors_sum[-1] > 0)) {
    cli::cli_alert_info("Some targets had associated errors or warnings")

    cli::cli_inform(c(
      "*" = "Fitting algorithm (errors): n = {errors_sum[2]}
      ({100 * (errors_sum[2]/k)}%)",
      "*" = "Fitting algorithm (warnings): n = {errors_sum[3]}
      ({100 * (errors_sum[3]/k)}%)",
      "*" = "Summary function (errors): n = {errors_sum[4]}
      ({100 * (errors_sum[4]/k)}%)",
      "*" = "Summary function (warnings): n = {errors_sum[5]}
      ({100 * (errors_sum[5]/k)}%)",
      "*" = "Evaluation function (errors): n = {errors_sum[6]}
      ({100 * (errors_sum[6]/k)}%)",
      "*" = "Evaluation function (warnings): n = {errors_sum[7]}
      ({100 * (errors_sum[7]/k)}%)"
    ))
  } else cli::cli_alert_info("No targets had associated errors or warnings")

  invisible(x)
}


#' Compose a swcontainer object for use in the seqwrap function.
#'
#' This function makes it possible to run checks on combined data sets
#' (meta data and target data) and fitting functions to avoid issues
#' in iterative modelling.
#'
#' @param x An optional named list or DGEList object (DESeqDataSet not yet
#' implemented).
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
#' @param samplename A character value indicating the variable by which
#' metadata can merge with the target data. This defaults to "seq_sample_id"
#' as this is used in the trainomeMetaData package.
#' @param additional_vars A vector of additional variables that is contained
#' in the metadata data set that is needed to fit the model. By default the
#' metadata is reduced to variables contained in the slots
#' formula/model/fixed/random in additional arguments.
#' More variables may be needed for offsets, weights etc.
#' @param summary_fun A custom (user-created) function for
#' evaluating/summarizing models. If NULL, a list of fitted models are returned
#' @param eval_fun A custom (user-created) function for model
#' diagnostics/evaluation. If NULL, no evaluation/diagnostics of models are
#' returned
#' @param exported A list of functions, values etc. to be passed to
#' summary_fun and eval_fun. This list must contain any functions that
#' should be used in model summarise or evaluations.
#' @param update A list of named parameters to update a swcontainer object.
#' @return A swcontainer object for direct use in seqwrap.
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


  # Extract the modelling algorithm for printing
  if (!is.null(modelfun)) {
    call_str <- match.call()
  }

  if (class(x) == "DGEList") {
    # Prepare count data from DGEList

    # If gene annotations are not available, use rownames
    if (is.null(x$genes)) {
      # If row names are numbers (1:nrow) add "target" for readability
      if (identical(rownames(x), as.character(1:nrow(x$counts)))) {
        xdata <- data.frame(
          target = paste0("target", 1:nrow(x$counts)),
          x$count
        )
        # Change name on sample names if provided
        if (!is.null(samplename)) colnames(xdata)[1] <- samplename
      } else {
        xdata <- data.frame(target = rownames(x$counts), x$count)
        # Change name on sample names if provided
        if (!is.null(samplename)) colnames(xdata)[1] <- samplename
      }
    }
    # If gene annotations are available, use these as identifiers
    if (!is.null(x$genes)) {
      xdata <- data.frame(target = x$genes[, 1], x$count)
      if (!is.null(samplename)) colnames(xdata)[1] <- samplename
    }

    # Prepare metadata from DGEList

    xmetadata <- data.frame(
      samplename = rownames(x$samples),
      x$samples,
      row.names = NULL
    )
    if (!is.null(samplename)) colnames(xmetadata)[1] <- samplename

    # Prepare gene-wise dispersion data

    container <- swcontainer(
      data = data_helper(data.frame(x$counts))
    )
  }

  # Initialize an empty container and add objects if they exist
  container <- swcontainer()



  # Set other properties
  if (!is.null(modelfun)) container@modelfun <- modelfun
  if (!is.null(modelfun)) container@model_print <- deparse(call_str$modelfun)
  if (!is.null(arguments)) container@arguments <- arguments
  if (!is.null(arguments)) container@arguments_print <- deparse(call_str$arguments)
  if (!is.null(data)) container@data <- data
  if (!is.null(rownames)) container@rownames <- rownames
  if (!is.null(metadata)) container@metadata <- metadata
  if (!is.null(targetdata)) container@targetdata <- targetdata
  if (!is.null(samplename)) container@samplename <- samplename
  if (!is.null(additional_vars)) container@additional_vars <- additional_vars
  if (!is.null(summary_fun)) container@summary_fun <- summary_fun
  if (!is.null(eval_fun)) container@eval_fun <- eval_fun
  if (!is.null(exported)) container@exported <- exported

  # Return the populated container
  return(container)
}


# Define a generic with named parameters only
seqwrap_update <- S7::new_generic("seqwrap_update", "container")

# Implement method for swcontainer class
S7::method(seqwrap_update, swcontainer) <- function(container, update = list()) {

  if (!S7::S7_inherits(container, swcontainer)) {
    stop("First argument must be a swcontainer object")
  }


  if (length(update) == 0) {
    return(container) # Nothing to update
  }


  # Get property names
 name <- names(update)


 # Special handling for modelfun to update model_print
  if ("modelfun" %in% name)  {

    # The model print information is not updated
    # due to difficulties getting the call from the
    # method. This is a temporary solution.
    update$model_print <- "The model function has been updated"

     }

  # Update the container object
  S7::props(container) <- update

 return(container)
}


#' Check swcontainer objects for seqwrap
#'
#' This function performs verbose checks/diagnostics of a swcontainer object
#'
#' @param x A swcontainer object
#' @return A list of diagnostics
seqwrap_check <- function(x, verbose = TRUE) {
  # Check if the container is a swcontainer object
  if (!S7::S7_inherits(container, swcontainer)) {
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
    n_targets <- length(unique(x@data[,1]))
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


  if (!all(meta_data_sample_id %in% data_sample_id) &&
      all(data_sample_id %in% meta_data_sample_id)) {
    cli::cli_abort(
    "The sample names in the meta data does not match the sample column names
    in the data. Check if the sample names are correctly formatted and that
    `samplename` looks for the correct variable in the meta data."
    )
  }

  # If verbose is TRUE, print the diagnostics
  if(verbose) {

    printfun <- function() {

      cli::cli_h1("seqwrap diagnostics")

      cli::cli_inform("The swcontainer object has been checked for
                      consistency with the following diagnostics:\n")

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
#' fitting algorithm.
#'
#'
#' @inheritParams seqwrap_compose
#' @param y, An swcontainer object created with seqwrap_compose, a named list
#' or a DGEList object.
#' @param return_models Logical, should models be returned as part of the
#' output? Save models during development on subsets of the data.
#' If used on large data sets, this will result in large memory usage.
#' @param save_models Logical, should models be saved? Models may be saved
#' on disk to save working memory.
#' @param model_path A character. The path to saved models.
#' @param subset A sequence, random samples or integers to indicate which
#' rows to keep in data. This is useful if you want to test the model in a
#' subset of targets. If keft to the default (NULL), all rows will be used.
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
#' @export
seqwrap <- function(
  y = NULL,
  modelfun = NULL,
  arguments = NULL,
  data = NULL,
  metadata = NULL,
  samplename = "seq_sample_id",
  additional_vars = NULL,
  summary_fun = NULL,
  eval_fun = NULL,
  exported = list(),
  return_models = TRUE,
  save_models = FALSE,
  model_path = NULL,
  subset = NULL,
  cores = 1,
  verbose = TRUE
) {
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
    if (length(updates) > 0) container <- seqwrap_compose(container,
                                                          update = updates)


  } else if (class(y) == "DGEList") {
    # If the input is a DGEList object, compose a new container
    container <- seqwrap_compose(x = y)

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
  } else {
    stop("The input must be a swcontainer object, a DGEList object or NULL")
  }

  # If subset is provided, subset the data
  if (!is.null(subset)) {
    if (is.numeric(subset)) {
      container@data <- container@data[subset, ]
      if (!is.null(container@targetdata)) {
        container@targetdata <- lapply(container@targetdata, function(x) {
          x[subset]
        })
      }
    } else {
      stop("Subset must be a numeric vector")
    }
  }


  # Check the container for consistency
  # seqwrap_check(container, verbose = verbose)

  # Get the number of targets and samples
  k <- nrow(container@data)
  n <- nrow(container@metadata)



  # data_helper function. Combine data into a list of data frames
  # containing variables, y in case of user-provided data frame;
  # names of variables from list names in case of user-provided
  # named list.

  if (is.null(container@targetdata)) {
    dfs <- data_helper(container@data)
  } else {
    dfs <- data_helper(container@data, container@targetdata)
  }

  # Determine the number of cores
  if (is.null(cores)) num_cores <- 1
  if (cores >= parallel::detectCores()) num_cores <- parallel::detectCores()
  if (cores <= parallel::detectCores()) num_cores <- cores

  # Catch the function calls for printing
  funcall <- match.call()
  ## Arguments to string
  call_arguments <- deparse(funcall$arguments)
  call_arguments <- sub("^.*?\\((.*)\\).*$", "\\1", call_arguments)
  ## modelfun to string
  call_engine <- deparse(funcall$modelfun)

  # Print pre-fit information
  cli::cli_h1("seqwrap")
  cli::cli_inform(
    "Initiating clusters for parallel processing with {num_cores} core{?s}"
  )

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

  # Parallel execution of the fitting process
  cli::cli_inform("Merging and modelling data")
  results <- pbapply::pblapply(
    cl = cl,
    X = dfs,
    FUN = seqwrap_mtf,
    samp_name = samplename,
    metdat = metadata,
    arg_list = arguments,
    add_vars = additional_vars,
    mt_summary_fun = summary_fun,
    mt_eval_fun = eval_fun,
    ffun = modelfun,
    return_mod = return_models,
    save_mods = save_models,
    mod_path = model_path
  )

  parallel::stopCluster(cl)

  # Combine results

  models <- NULL
  summaries <- NULL
  evaluations <- NULL
  errors <- NULL

  # Collect models
  if (return_models) models <- lapply(results, `[[`, "model")

  if (!is.null(summary_fun)) summaries <- lapply(results, `[[`, "summaries")
  if (!is.null(eval_fun)) evaluations <- lapply(results, `[[`, "evaluation")

  ## Create a data frame of all errors/warnings
  errors <- tibble::tibble(
    target = names(results),
    errors_fit = lapply(results, `[[`, "err"),
    warnings_fit = lapply(results, `[[`, "warn"),
    err_sum = lapply(results, `[[`, "err_sum"),
    warn_sum = lapply(results, `[[`, "warn_sum"),
    err_eval = lapply(results, `[[`, "err_eval"),
    warn_eval = lapply(results, `[[`, "warn_eval")
  )

  # Count non-null values in the error/warning data frame
  errors_sum <- sapply(errors, function(x) sum(!sapply(x, is.null)))

  ## Evaluate errors for the resulting print function
  cli::cli_inform("Completed model fitting and evaluation")

  if (any(errors_sum[-1] > 0)) {
    cli::cli_alert_info("Some targets had associated errors or warnings")

    cli::cli_inform(c(
      "*" = "Fitting algorithm (errors): n = {errors_sum[2]}
      ({100 * (errors_sum[2]/k)}%)",
      "*" = "Fitting algorithm (warnings): n = {errors_sum[3]}
      ({100 * (errors_sum[3]/k)}%)",
      "*" = "Summary function (errors): n = {errors_sum[4]}
      ({100 * (errors_sum[4]/k)}%)",
      "*" = "Summary function (warnings): n = {errors_sum[5]}
      ({100 * (errors_sum[5]/k)}%)",
      "*" = "Evaluation function (errors): n = {errors_sum[6]}
      ({100 * (errors_sum[6]/k)}%)",
      "*" = "Evaluation function (warnings): n = {errors_sum[7]}
      ({100 * (errors_sum[7]/k)}%)"
    ))
  }

  ## Combine the results into a seqwrapResults
  comb_results <- seqwrapResults(
    models = models,
    summaries = summaries,
    evaluations = evaluations,
    errors = errors,
    n = n,
    k = k,
    call_arguments = container@arguments_print,
    call_engine = container@model_print
  )

  return(comb_results)
}



#' Summarise seqwrapResults objects
#'
#' @param x A seqwrapResults object
#' @param summaries Logical, should summaries be combined?
#' @param evaluations Logical, should evaluations be combined?
#' @param errors Logical, should errors be combined?
#'
#' @details
#' This functions attempts to summarise results from the summary and evaluation
#' functions applied in each iteration during modelling. The function expects
#' that the summary and evaluation functions return data frames.
#'
#'
#'
#' @export
seqwrap_summarise <- function(x,
                              summaries = TRUE,
                              evaluations = TRUE,
                              errors = TRUE,
                              verbose = TRUE) {

  # Check if the input is a seqwrapResults object
  if (!S7::S7_inherits(x, seqwrapResults)) {
    stop("The input must be a seqwrapResults object")
  }


  ## Print information from the seqwrapResults object
  print_info <-function() {
    cli::cli_h1("seqwrap summarise")

    cli::cli_li(
      "A total of {x@n} sample{?s} and {x@k} target{?s} where
                  used in {.code {x@call_engine}} with arguments
                  {.code {x@call_arguments}}"
    )
    if (summaries) {
      cli::cli_li("Attempting to combine results from
                  the provided summary function.")
    }

    if (evaluations) {
      cli::cli_li("Attempting to combine results from
                  the provided evaluations function.")
    }

    if (errors) {
      cli::cli_li("Attempting to summarise errors and
                  warnings from the fitting process.")
    }

  }

  if (verbose) print_info()



  ## Initialize results variables ##
  summarised_results_final <- NULL
  evaluated_results_final <- NULL

  # print(x)

  ## Extract summarises
  if (summaries) {


    # Check for NULL or empty list
    if (is.null(x@summaries) || all(sapply(x@summaries, is.null)))  {

      print_summary <- function() {
      cli::cli_h1("Model summaries")
      cli::cli_alert_info("No summaries available")
      }
      if (verbose) print_summary()

      } else {

      # Count the number of non-null elements in the list of summaries
      n_summaries <- sum(!sapply(x@summaries, is.null))



      print_summary <- function() {
        cli::cli_h1("Model summaries")
        cli::cli_li("{n_summaries} targets have associated summaries")

      }

      if (verbose) print_summary()

      # Filter out NULL elements from the list of summaries
      x@summaries <- x@summaries[!sapply(x@summaries, is.null)]

      # Use Map to put names in the target column
      temp_list_summaries <- Map(function(df, df_name) {
        df[["target"]] <- df_name
        # Move the target column to the first column
        df <- df[, c("target", setdiff(colnames(df), "target")),
                 drop = FALSE]
        # Return data frames
        df
      }, x@summaries, names(x@summaries))

      # Bind the list of data frames
      summarised_results_final <- do.call(rbind, temp_list_summaries)
      # NOW set rownames to NULL on the combined data frame
      rownames(summarised_results_final) <- NULL

    }
  } # End if (summaries)

  ## Extract evaluations
  if (evaluations) {

    # Check for NULL or empty list
    if (is.null(x@evaluations) || all(sapply(x@evaluations, is.null))) {

      print_evaluation <- function() {
        cli::cli_h1("Model evaluations")
        cli::cli_alert_info("No evaluations available")
      }
      if (verbose) print_evaluation()
      # evaluated_results_final remains NULL
    } else {

      # Count the number of non-null elements in the list of evaluations
      n_evals <- sum(!sapply(x@evaluations, is.null))



      print_evaluation <- function() {
        cli::cli_h1("Model evaluations")
        cli::cli_li("{n_evals} targets have associated evaluation results")

      }

      if (verbose) print_evaluation()

      # Filter out NULL elements from the list of evaluations
      x@evaluations <- x@evaluations[!sapply(x@evaluations, is.null)]

      # Use Map to process each data frame in the list
      temp_list_evaluations <- Map(function(df, df_name) {
        df[["target"]] <- df_name
        # Move the target column to the first column
        df <- df[, c("target", setdiff(colnames(df), "target")),
                 drop = FALSE]

        df
      }, x@evaluations, names(x@evaluations))

      # Bind the list of data frames
      evaluated_results_final <- do.call(rbind, temp_list_evaluations)
      rownames(evaluated_results_final) <- NULL

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

  # Add error handling results here if implemented

  # Check if the final list is empty
  if (length(results) == 0) {
    cli::cli_alert_warning("No results were generated, check your input.")
  }

  # Check if the final list is empty
  if (length(results) != 0) {
  if (verbose) {
    cli::cli_h1("Combined results")
    cli::cli_alert_info("Combined results have been generated and were
                        silently returned.")
  }

  }

  # Return the final list object
  return(invisible(results))

}















