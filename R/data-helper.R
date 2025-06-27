#' Data helper function.
#'
#' The function takes a named list of data frames
#' and combines them into target-wise data frames
#' for use in subsequent modelling.
#'
#' If the input is a data frame, the function returns a
#' a list of data frames with a single row (y in subsequent modelling)
#'
#' @param dat Sample/target-wise data, either a data frame or a named list.
#' @param targetdat, target-wise data, e.g. targe-wise dispersion
#' @returns A list of lists, each element of the list has a target specific
#' data frame over all samples and a target-specific data frame with values
#' specific to target.
#' @keywords internal
data_helper <- function(dat, targetdat = NULL, rownames = FALSE) {

  # Convert tibble to data.frame. This assumes that data comes in
  # list or data frame, not nested in tibbles.
  if (any(c("tbl_df", "tbl") %in% class(dat))) dat <- data.frame(dat)

  # Check if input is a list or a data frame
  if (!is.list(dat) && !is.data.frame(dat)) {
    stop("Input must be either a list of data frames or a single data frame")
  }

  # If the input is a data frame
  if (inherits(dat, "data.frame")) {

    if (rownames) {
      cli::cli_inform("Using row names as target identification")

      col_names <- colnames(dat)
      dat[["target_id"]] <- rownames(dat)
      # Move this column to be the first column
      dat <- dat[, c("target_id", col_names)]
      # Reset row names to numeric
      rownames(dat) <- NULL
    }

    # Split the data into a list of data frames for each target
    dfs <- split(dat[, seq_len(ncol(dat))[-1]], dat[, 1])

    dfs <- lapply(dfs, function(x) {
      rownames(x) <- "y"
      return(x)
    })


  }

  if (inherits(dat, "list")) {
    if (rownames) {
      cli::cli_inform("Using row names as target identification")
      dat <- lapply(dat, function(x) {
        col_names <- colnames(x)
        x[["target_id"]] <- rownames(x)
        # Move this column to be the first column
        x <- x[, c("target_id", col_names)]
        # Reset row names to numeric
        rownames(x) <- NULL
        return(x)
      })
    }

    # Check if all elements are data frames
    if (!all(sapply(dat, is.data.frame))) {
      stop("All elements in the list must be data frames")
    }

    # Get the names of the input data frames
    df_names <- names(dat)
    if (is.null(df_names)) {
      stop("The list of data frames must be named")
    }

    # Find the number of rows in each data frame
    row_counts <- sapply(dat, nrow)

    # Check if all data frames have the same number of rows
    if (length(unique(row_counts)) != 1) {
      stop("All data frames must have the same number of rows")
    }

    # Check if all data frames have at least one column
    if (any(sapply(dat, ncol) < 1)) {
      stop("All data frames must have at least one column")
    }

    # Create a list to store the new data frames
    dfs <- list()

    # For each row index
    for (i in 1:row_counts[1]) {
      # Create a new data frame for this row
      new_df <- data.frame(row.names = df_names)

      # Get the name for this result element (from first column of first data frame)
      first_df <- dat[[1]]
      first_col_name <- colnames(first_df)[1]
      result_name <- as.character(first_df[i, 1])

      # For each original data frame
      for (j in 1:length(dat)) {
        df_name <- df_names[j]
        current_df <- dat[[j]]

        # Extract the row from the current data frame
        row_data <- current_df[i, , drop = FALSE]

        # For each column in the current data frame (skip the first column)
        for (col_name in colnames(current_df)[-1]) {
          # Add the value to the new data frame
          new_df[df_name, col_name] <- row_data[[col_name]]
        }
      }

      # Add the new data frame to the result list with appropriate name
      dfs[[result_name]] <- new_df
    }


  }

  # Include targetwise data if it is not null
  if(!is.null(targetdat)) {
    # Target-wise data, common to all samples
    # Check if input is a data frame or a vector
    if (!is.vector(targetdat) && !is.data.frame(targetdat)) {
      stop("Target-wise data must be either a vector of data frames or
         a single data frame")
    }

    # If the target-wise data is a vector, adding "target.wise"
    # as heading
    if (is.vector(targetdat)) {
      target.wise.dat <- data.frame(target.wise = targetdat)

      dfs.targetwise <- split(target.wise.dat, seq(nrow(target.wise.dat)))
      names(dfs.targetwise) <- rownames(target.wise.dat)


    }
    # If the target-wise data is a data frame, columns of the data available
    # data frame will be as is
    if (is.data.frame(targetdat)) {
      dfs.targetwise <- split(targetdat, seq(nrow(targetdat)))
      names(dfs.targetwise) <- rownames(targetdat)
    }

    # Combine into a list: target/sample data and target-wise data
    combined.dfs <- mapply(function(x, y) list(x, y),
                           dfs,
                           dfs.targetwise,
                           SIMPLIFY = FALSE)



  }
  # If targetdat is not provided, keep structure of the list anyway
  if (is.null(targetdat)) {
    combined.dfs <- lapply(dfs, function(x) list(x, NULL))
  }

  return(combined.dfs)

}
