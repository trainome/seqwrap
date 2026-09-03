## Chunk-level helpers used by seqwrap() ------------------------------------
##
## The unit of work handed to a worker is a *chunk* of targets rather than a
## single target. This keeps the number of cluster round trips proportional to
## the number of chunks rather than the number of targets, lets summary and
## evaluation data frames be reduced (rbind-ed) while still on the worker, and
## provides a natural granularity for spilling results to disk.
##
## Per-target fitting itself is unchanged: seqwrap_chunk() calls seqwrap_mtf()
## once per target, so model saving (save_models/model_path) behaves exactly as
## it did when seqwrap_mtf() was mapped over targets directly.


#' Target order used by data_helper()
#'
#' `data_helper()` orders targets differently depending on its input: the data
#' frame branch groups with `split()`, which sorts target identifiers, while the
#' list branch walks rows in their original order. Chunking must reproduce that
#' ordering so that a chunked run returns targets in the same sequence as an
#' unchunked one.
#'
#' @param data The `data` slot of a swcontainer.
#' @param rownames Logical, are row names used as target identifiers?
#' @return An integer vector of row positions in `data_helper()` output order.
#' @keywords internal
#' @noRd
sw_target_order <- function(data, rownames = FALSE) {
  if (is.data.frame(data)) {
    # data_helper() groups with split(), which sorts on whichever column holds
    # the target identifier
    if (isTRUE(rownames)) order(rownames(data)) else order(data[[1]])
  } else {
    seq_len(nrow(data[[1]]))
  }
}


#' Split target positions into contiguous chunks
#'
#' @param k Number of targets.
#' @param chunk_size Number of targets per chunk.
#' @return A list of integer vectors of positions.
#' @keywords internal
#' @noRd
sw_chunk_index <- function(k, chunk_size) {
  if (k < 1L) {
    return(list())
  }
  chunk_size <- max(1L, as.integer(chunk_size))
  starts <- seq.int(1L, k, by = chunk_size)
  lapply(starts, function(s) seq.int(s, min(s + chunk_size - 1L, k)))
}


#' Choose a chunk size
#'
#' Aims for several chunks per worker so that load balancing still works, while
#' keeping chunks large enough to amortise the cost of a cluster round trip.
#'
#' @param k Number of targets.
#' @param num_cores Number of workers.
#' @return A single integer.
#' @keywords internal
#' @noRd
sw_default_chunk_size <- function(k, num_cores) {
  if (k < 1L) {
    return(1L)
  }
  max(1L, min(2000L, as.integer(ceiling(k / (max(1L, num_cores) * 4L)))))
}


#' Subset container data by row position
#'
#' @param data A data frame, or a named list of data frames.
#' @param rows Integer row positions.
#' @return An object of the same shape as `data`, restricted to `rows`.
#' @keywords internal
#' @noRd
sw_slice_data <- function(data, rows) {
  if (is.data.frame(data)) {
    data[rows, , drop = FALSE]
  } else {
    lapply(data, function(d) d[rows, , drop = FALSE])
  }
}


#' Subset target-wise data by position
#'
#' @param targetdata A data frame, a list, or NULL.
#' @param pos Integer positions.
#' @return An object of the same shape as `targetdata`, or NULL.
#' @keywords internal
#' @noRd
sw_slice_targetdata <- function(targetdata, pos) {
  if (is.null(targetdata)) {
    return(NULL)
  }
  if (is.data.frame(targetdata)) {
    targetdata[pos, , drop = FALSE]
  } else {
    targetdata[pos]
  }
}


#' Make a target identifier safe to use as a file name
#'
#' Target identifiers routinely contain characters that are not valid in file
#' names (`/`, `:`, `|`), and Windows additionally reserves a handful of device
#' names. Unsafe characters are replaced with an underscore.
#'
#' @param x A character vector of target identifiers.
#' @return A character vector safe for use as a file name.
#' @keywords internal
#' @noRd
sw_safe_filename <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return("target")
  }
  out <- gsub("[^A-Za-z0-9._-]", "_", as.character(x))
  # Reserved device names on Windows, matched without their extension.
  reserved <- c(
    "CON", "PRN", "AUX", "NUL",
    paste0("COM", 1:9), paste0("LPT", 1:9)
  )
  bad <- is.na(out) | out == "" | toupper(out) %in% reserved
  out[bad] <- paste0("target_", seq_along(out))[bad]
  out
}


#' Bind a named list of per-target data frames into a single data frame
#'
#' NULL entries (targets whose summary or evaluation function failed) are
#' dropped. The target identifier is added as the first column, matching the
#' layout produced by `seqwrap_summarise()`.
#'
#' @param x A named list of data frames, possibly containing NULLs.
#' @return A single data frame, or NULL when every entry is NULL.
#' @keywords internal
#' @noRd
sw_bind_targets <- function(x) {
  if (length(x) == 0L) {
    return(NULL)
  }

  # Drop entries that carry no information. Besides NULL, a summary function
  # handed a failed fit often returns a zero-column data frame rather than
  # signalling an error: broom.mixed::tidy(NULL) returns a 0 x 0 tibble. Adding
  # a target column to one of those would silently invent a row of all-missing
  # values, so they are discarded here.
  keep <- vapply(
    x,
    function(df) !is.null(df) && NCOL(df) > 0L && NROW(df) > 0L,
    logical(1)
  )
  if (!any(keep)) {
    return(NULL)
  }
  x <- x[keep]

  # Binding requires data frames. Without this check the failure would surface
  # as an opaque subsetting error from deep inside a worker.
  not_a_frame <- !vapply(x, is.data.frame, logical(1))
  if (any(not_a_frame)) {
    first <- which(not_a_frame)[1]
    cli::cli_abort(c(
      "Summary and evaluation functions must return a data frame when
       {.arg cache} combines results.",
      "x" = "Target {.val {names(x)[first]}} returned
             {.cls {class(x[[first]])}}.",
      "i" = "Use {.code cache = \"none\"} to keep one result per target."
    ))
  }

  bound <- Map(
    function(df, nm) {
      df[["target"]] <- nm
      df[, c("target", setdiff(colnames(df), "target")), drop = FALSE]
    },
    x,
    names(x)
  )

  sw_rbind_fill(bound)
}


#' Bind data frames that may not share the same columns
#'
#' `rbind()` fails outright when its arguments have different columns, which
#' would discard a whole run's results because a handful of targets returned a
#' differently shaped summary. Missing columns are filled with NA instead.
#'
#' @param parts A list of data frames.
#' @return A single data frame, or NULL when `parts` is empty.
#' @keywords internal
#' @noRd
sw_rbind_fill <- function(parts) {
  if (length(parts) == 0L) {
    return(NULL)
  }

  columns <- lapply(parts, colnames)
  all_columns <- unique(unlist(columns, use.names = FALSE))

  # The common case: every part has the same columns in the same order
  consistent <- all(vapply(
    columns,
    function(nm) identical(nm, all_columns),
    logical(1)
  ))

  if (!consistent) {
    parts <- lapply(parts, function(df) {
      missing <- setdiff(all_columns, colnames(df))
      for (nm in missing) df[[nm]] <- NA
      df[, all_columns, drop = FALSE]
    })
  }

  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}


# Maps the condition fields returned by seqwrap_mtf() onto the stage of the
# run and the kind of condition they represent.
.sw_condition_map <- list(
  list(field = "err", stage = "fit", type = "error"),
  list(field = "warn", stage = "fit", type = "warning"),
  list(field = "err_sum", stage = "summary", type = "error"),
  list(field = "warn_sum", stage = "summary", type = "warning"),
  list(field = "err_eval", stage = "evaluation", type = "error"),
  list(field = "warn_eval", stage = "evaluation", type = "warning")
)

# The stages and condition types a user can filter on
.sw_stages <- c("fit", "summary", "evaluation")
.sw_types <- c("error", "warning")


#' Is the fitting function nlme::lme?
#'
#' `nlme::lme` needs a more cautious calling convention than the other engines,
#' so it has to be recognised before fitting. The answer depends only on the
#' fitting function, which is fixed for a whole run, so this is computed once in
#' `seqwrap()` and passed down rather than recomputed for every target. The
#' `deparse()` fallback is the expensive part: between 0.25 and 0.8 ms per call
#' depending on the size of the function.
#'
#' @param ffun The fitting function.
#' @return TRUE when `ffun` is, or wraps, `nlme::lme`.
#' @keywords internal
#' @noRd
sw_is_lme <- function(ffun) {
  identical(ffun, nlme::lme) ||
    (is.character(ffun) && ffun == "nlme::lme") ||
    any(grepl("lme$", deparse(ffun)))
}


#' An empty condition table
#'
#' @return A zero row tibble with the condition table columns.
#' @keywords internal
#' @noRd
sw_empty_conditions <- function() {
  tibble::tibble(
    target = character(0),
    stage = character(0),
    type = character(0),
    class = character(0),
    message = character(0)
  )
}


#' Collect conditions into one row per condition raised
#'
#' Storing one entry per target for each of the six stage and type combinations
#' costs memory proportional to the number of targets even when nothing went
#' wrong, and condition objects themselves are around 1.6 kB each because they
#' carry their originating call. Recording only the conditions that were
#' actually raised, as message text plus the condition class, is far smaller and
#' much easier to filter.
#'
#' @param fits A named list of `seqwrap_mtf()` return values.
#' @param targets Target identifiers, in the same order as `fits`.
#' @return A tibble with columns target, stage, type, class and message.
#' @keywords internal
#' @noRd
sw_conditions_long <- function(fits, targets) {
  parts <- lapply(.sw_condition_map, function(entry) {
    conditions <- lapply(fits, `[[`, entry$field)
    hit <- which(!vapply(conditions, is.null, logical(1)))

    if (length(hit) == 0L) {
      return(NULL)
    }
    conditions <- conditions[hit]

    tibble::tibble(
      target = targets[hit],
      stage = entry$stage,
      type = entry$type,
      class = vapply(conditions, function(z) class(z)[1], character(1)),
      message = vapply(conditions, conditionMessage, character(1))
    )
  })

  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) == 0L) {
    return(sw_empty_conditions())
  }

  out <- do.call(rbind, parts)
  # Report in the order targets were fitted rather than grouped by stage
  out <- out[order(match(out$target, targets), out$stage, out$type), ]
  rownames(out) <- NULL
  out
}


#' Fit every target in one chunk
#'
#' Expands a chunk of raw target data into per-target data frames on the worker,
#' fits each target with `seqwrap_mtf()`, and optionally reduces and spills the
#' resulting summary and evaluation data frames.
#'
#' @param chunk A list with elements `data`, `targetdata` and `id`.
#' @param samp_name Sample name variable.
#' @param metdat Metadata.
#' @param arg_list Arguments passed to the fitting function.
#' @param add_vars Additional metadata variables to retain.
#' @param mt_summary_fun Summary function.
#' @param mt_eval_fun Evaluation function.
#' @param ffun Fitting function.
#' @param return_mod Should models be returned?
#' @param save_mods Should models be written to `mod_path`?
#' @param mod_path Directory for saved models.
#' @param reduce Should summaries and evaluations be bound within the chunk?
#' @param cache_path Directory for spilled chunk results, or NULL.
#' @param is_lme Logical, is `ffun` `nlme::lme`? Computed from `ffun` when not
#'   supplied.
#' @param rownames Logical, should row names of `chunk$data` be used as target
#'   identifiers?
#' @return A list carrying either the raw per-target fits (`reduce = FALSE`) or
#'   the reduced chunk results, plus the file paths written when spilling.
#' @keywords internal
#' @noRd
seqwrap_chunk <- function(
  chunk,
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
  reduce = FALSE,
  cache_path = NULL,
  is_lme = NULL,
  rownames = FALSE
) {
  # Determined once per run by seqwrap(); resolved here if called directly
  if (is.null(is_lme)) {
    is_lme <- sw_is_lme(ffun)
  }

  # Expand this chunk into per-target data frames on the worker. Building these
  # here rather than in the parent is what keeps the parent's memory flat.
  dfs <- data_helper(chunk$data, chunk$targetdata, rownames = rownames)

  # Map() rather than lapply() so that each target's identifier reaches
  # seqwrap_mtf(), which needs it to name saved model files.
  fits <- Map(
    function(one, nm) {
      seqwrap_mtf(
        one,
        samp_name = samp_name,
        metdat = metdat,
        arg_list = arg_list,
        add_vars = add_vars,
        mt_summary_fun = mt_summary_fun,
        mt_eval_fun = mt_eval_fun,
        ffun = ffun,
        return_mod = return_mod,
        save_mods = save_mods,
        mod_path = mod_path,
        target_name = nm,
        is_lme = is_lme
      )
    },
    dfs,
    names(dfs)
  )
  names(fits) <- names(dfs)

  # Without reduction the chunk is transparent: the caller sees exactly the
  # per-target structure it would have seen from an unchunked run.
  if (!reduce) {
    return(list(fits = fits))
  }

  summaries <- sw_bind_targets(lapply(fits, `[[`, "summaries"))
  evaluations <- sw_bind_targets(lapply(fits, `[[`, "evaluation"))

  errors <- sw_conditions_long(fits, names(fits))

  models <- if (return_mod) lapply(fits, `[[`, "model") else NULL

  file_summaries <- NA_character_
  file_evaluations <- NA_character_

  if (!is.null(cache_path)) {
    if (!is.null(summaries)) {
      file_summaries <- file.path(
        cache_path,
        sprintf("summaries-%06d.rds", chunk$id)
      )
      # Compression dominates write time for objects this small and the files
      # are transient, so it is not worth paying for.
      saveRDS(summaries, file = file_summaries, compress = FALSE)
    }
    if (!is.null(evaluations)) {
      file_evaluations <- file.path(
        cache_path,
        sprintf("evaluations-%06d.rds", chunk$id)
      )
      saveRDS(evaluations, file = file_evaluations, compress = FALSE)
    }
    summaries <- NULL
    evaluations <- NULL
  }

  list(
    targets = names(fits),
    models = models,
    summaries = summaries,
    evaluations = evaluations,
    errors = errors,
    file_summaries = file_summaries,
    file_evaluations = file_evaluations
  )
}


#' Bind the reduced data frames carried back by each chunk
#'
#' @param chunk_results The list returned by the chunk-level `pblapply()`.
#' @param field Either `"summaries"` or `"evaluations"`.
#' @return A single data frame, or NULL when every chunk was empty.
#' @keywords internal
#' @noRd
sw_read_chunks <- function(chunk_results, field) {
  parts <- lapply(chunk_results, `[[`, field)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  sw_rbind_fill(parts)
}


#' Read a set of cached chunk files and bind them
#'
#' @param files A character vector of file paths, NA entries are skipped.
#' @return A single data frame, or NULL when nothing could be read.
#' @keywords internal
#' @noRd
sw_read_cache <- function(files) {
  files <- files[!is.na(files)]
  files <- files[file.exists(files)]
  if (length(files) == 0L) {
    return(NULL)
  }

  sw_rbind_fill(lapply(files, readRDS))
}


#' Normalise the drop_warnings argument to a set of stages
#'
#' @param drop_warnings FALSE, TRUE, or a character vector of stages.
#' @return A character vector of stages, possibly empty.
#' @keywords internal
#' @noRd
sw_drop_stages <- function(drop_warnings) {
  if (is.null(drop_warnings) || isFALSE(drop_warnings)) {
    return(character(0))
  }
  if (isTRUE(drop_warnings)) {
    return(.sw_stages)
  }

  stages <- as.character(drop_warnings)
  invalid <- setdiff(stages, .sw_stages)
  if (length(invalid) > 0L) {
    # Bound to a local name: cli reads an interpolated name starting with a dot
    # as a cli literal rather than as a variable.
    valid <- .sw_stages
    cli::cli_abort(
      "{.arg drop_warnings} must be TRUE, FALSE, or a subset of
      {.val {valid}}; got {.val {invalid}}."
    )
  }

  stages
}


#' Identify targets that raised a warning at the given stages
#'
#' @param errors The errors slot of a seqwrapResults object, in long form.
#' @param stages A character vector of stages.
#' @return A character vector of target identifiers.
#' @keywords internal
#' @noRd
sw_warned_targets <- function(errors, stages) {
  if (length(stages) == 0L || NROW(errors) == 0L) {
    return(character(0))
  }

  unique(errors$target[errors$type == "warning" & errors$stage %in% stages])
}


#' Collect summaries or evaluations from a seqwrapResults object
#'
#' Handles the three representations produced by the `cache` argument of
#' `seqwrap()`: a named list of per-target data frames (`"none"`), a single
#' bound data frame (`"memory"`), or a set of chunk files on disk (`"disk"`).
#'
#' @param x A seqwrapResults object.
#' @param slot Either `"summaries"` or `"evaluations"`.
#' @return A single data frame with a `target` column, or NULL.
#' @keywords internal
#' @noRd
sw_collect <- function(x, slot) {
  if (!is.null(x@cache)) {
    return(sw_read_cache(x@cache$chunks[[slot]]))
  }

  value <- if (identical(slot, "summaries")) x@summaries else x@evaluations

  if (is.null(value)) {
    return(NULL)
  }
  if (is.data.frame(value)) {
    return(value)
  }

  sw_bind_targets(value)
}


#' Errors and warnings from a seqwrap run
#'
#' Returns the errors and warnings raised during a `seqwrap()` run, with one row
#' per condition and optional filtering. Targets that completed cleanly do not
#' appear, so the result is proportional to the number of problems rather than
#' to the number of targets.
#'
#' @param x A seqwrapResults object.
#' @param stage Optional. One or more of `"fit"`, `"summary"` and
#'   `"evaluation"`, to restrict the result to conditions raised at those
#'   stages.
#' @param type Optional. One or both of `"error"` and `"warning"`.
#' @param target Optional. Target identifiers to restrict the result to.
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{target}{The target identifier.}
#'   \item{stage}{Where the condition arose: `"fit"`, `"summary"` or
#'     `"evaluation"`.}
#'   \item{type}{Either `"error"` or `"warning"`.}
#'   \item{class}{The class of the condition object, useful for separating
#'     specific warning types from generic ones.}
#'   \item{message}{The condition message.}
#' }
#'
#' @details
#' To identify targets that completed without any condition, compare against the
#' `targets` property, which lists every target that was fitted:
#' `setdiff(x@targets, seqwrap_errors(x)$target)`.
#'
#' @examples
#' library(seqwrap)
#'
#' dat <- simcounts(n_genes = 6)
#'
#' container <- seqwrap_compose(
#'   modelfun = stats::lm,
#'   arguments = list(formula = y ~ x),
#'   data = dat$data,
#'   metadata = dat$metadata,
#'   samplename = "sample"
#' )
#'
#' results <- seqwrap(container, cores = 1, verbose = FALSE)
#'
#' # Every condition raised during the run
#' seqwrap_errors(results)
#'
#' # Only warnings from the fitting algorithm
#' seqwrap_errors(results, stage = "fit", type = "warning")
#'
#' # Targets that completed cleanly
#' setdiff(results@targets, seqwrap_errors(results)$target)
#'
#' @export
seqwrap_errors <- function(x, stage = NULL, type = NULL, target = NULL) {
  if (!S7::S7_inherits(x, seqwrapResults)) {
    stop("The input must be a seqwrapResults object")
  }

  # Validate before the empty shortcut, so that a mistyped filter is reported
  # whether or not the run happened to raise any conditions
  if (!is.null(stage)) {
    invalid <- setdiff(stage, .sw_stages)
    if (length(invalid) > 0L) {
      valid <- .sw_stages
      cli::cli_abort(
        "{.arg stage} must be a subset of {.val {valid}};
        got {.val {invalid}}."
      )
    }
  }
  if (!is.null(type)) {
    invalid <- setdiff(type, .sw_types)
    if (length(invalid) > 0L) {
      valid <- .sw_types
      cli::cli_abort(
        "{.arg type} must be a subset of {.val {valid}};
        got {.val {invalid}}."
      )
    }
  }

  out <- x@errors
  if (NROW(out) == 0L) {
    return(sw_empty_conditions())
  }

  if (!is.null(stage)) {
    out <- out[out$stage %in% stage, , drop = FALSE]
  }
  if (!is.null(type)) {
    out <- out[out$type %in% type, , drop = FALSE]
  }

  if (!is.null(target)) {
    out <- out[out$target %in% target, , drop = FALSE]
  }

  rownames(out) <- NULL
  out
}


#' Remove a seqwrap disk cache
#'
#' Deletes the chunk files written when `seqwrap()` was called with
#' `cache = "disk"`. A cache placed in the session temporary directory (the
#' default) is removed when the R session ends, so calling this function is only
#' necessary to reclaim space earlier, or when a user-supplied `cache_path` was
#' used.
#'
#' @param x A seqwrapResults object created with `cache = "disk"`.
#' @param recursive Logical. Should the cache directory itself be removed in
#'   addition to the chunk files? Defaults to TRUE.
#'
#' @return The number of files removed, invisibly.
#'
#' @examples
#' library(seqwrap)
#'
#' dat <- simcounts(n_genes = 6)
#'
#' container <- seqwrap_compose(
#'   modelfun = stats::lm,
#'   arguments = list(formula = y ~ x),
#'   data = dat$data,
#'   metadata = dat$metadata,
#'   samplename = "sample"
#' )
#'
#' results <- seqwrap(
#'   container,
#'   eval_fun = function(x) NULL,
#'   cache = "disk",
#'   cores = 1,
#'   verbose = FALSE
#' )
#'
#' # Combine the cached chunks before clearing them
#' summaries <- seqwrap_summarise(results, verbose = FALSE)
#'
#' seqwrap_cache_clear(results)
#'
#' @export
seqwrap_cache_clear <- function(x, recursive = TRUE) {
  if (!S7::S7_inherits(x, seqwrapResults)) {
    stop("The input must be a seqwrapResults object")
  }

  if (is.null(x@cache)) {
    cli::cli_alert_info("No disk cache is associated with this object.")
    return(invisible(0L))
  }

  files <- c(x@cache$chunks$summaries, x@cache$chunks$evaluations)
  files <- files[!is.na(files)]
  files <- files[file.exists(files)]

  removed <- sum(file.remove(files))

  if (recursive && dir.exists(x@cache$path)) {
    unlink(x@cache$path, recursive = TRUE)
  }

  invisible(removed)
}


#' Target identifiers of a data slot, in row order
#'
#' @param data A data frame, or a named list of data frames.
#' @param rownames Logical, are row names used as target identifiers?
#' @return A character vector.
#' @keywords internal
#' @noRd
sw_target_ids <- function(data, rownames = FALSE) {
  ids <- if (isTRUE(rownames)) {
    if (is.data.frame(data)) rownames(data) else rownames(data[[1]])
  } else if (is.data.frame(data)) {
    data[[1]]
  } else {
    data[[1]][[1]]
  }
  as.character(ids)
}


#' Align target-wise data with the targets in the data
#'
#' A data frame or an unnamed list is paired with targets by position, element
#' `i` belonging to row `i` of the data. A list named by target identifier is
#' paired by name, which also allows it to hold more targets than the data,
#' as when priors built from a full run are applied to a subset.
#'
#' @param targetdata A data frame, a list, or NULL.
#' @param target_ids Target identifiers in data row order.
#' @param warn Logical, warn when an unnamed list is matched by position?
#' @return `targetdata` reordered to `target_ids` when matched by name,
#'   otherwise unchanged.
#' @keywords internal
#' @noRd
sw_align_targetdata <- function(targetdata, target_ids, warn = TRUE) {
  if (is.null(targetdata) || is.data.frame(targetdata)) {
    return(targetdata)
  }

  nms <- names(targetdata)
  unnamed <- sw_is_unnamed(targetdata)

  if (unnamed) {
    if (length(targetdata) != length(target_ids)) {
      cli::cli_abort(
        "{.arg targetdata} must have the same number of elements
         ({length(targetdata)}) as {.arg data} has rows ({length(target_ids)})."
      )
    }
    if (warn) {
      cli::cli_warn(c(
        "{.arg targetdata} is an unnamed list and is matched to targets by
         position.",
        "i" = "Make sure that element {.code i} of the list belongs to row
               {.code i} of {.arg data}, or name the elements by target
               identifier to match them by name."
      ))
    }
    return(targetdata)
  }

  if (any(is.na(nms) | !nzchar(nms))) {
    cli::cli_abort(
      "Either every element of {.arg targetdata} must be named by target
       identifier, or none of them."
    )
  }
  if (anyDuplicated(nms)) {
    dups <- unique(nms[duplicated(nms)])
    cli::cli_abort(c(
      "Element names of {.arg targetdata} must be unique.",
      "x" = "Duplicated: {.val {utils::head(dups, 5)}}."
    ))
  }

  missing <- setdiff(target_ids, nms)
  if (length(missing) > 0L) {
    n_missing <- length(missing)
    cli::cli_abort(c(
      "{.arg targetdata} is matched to targets by name, but {n_missing}
       target{?s} in {.arg data} {?has/have} no element.",
      "x" = "Missing: {.val {utils::head(missing, 5)}}."
    ))
  }

  targetdata[match(target_ids, nms)]
}


#' Does a list carry no usable names?
#'
#' @param x A list.
#' @return TRUE when `x` has no names, or only empty ones.
#' @keywords internal
#' @noRd
sw_is_unnamed <- function(x) {
  nms <- names(x)
  is.null(nms) || all(is.na(nms) | !nzchar(nms))
}
