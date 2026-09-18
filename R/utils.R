# ============================================================================
# Shared utility functions for fastsae package
# ============================================================================

#' Extract variable from data frame or use as-is
#'
#' @description This helper function accepts a variable as a vector, column name,
#' or one-sided formula and returns the corresponding data column.
#'
#' @param data A data frame or data frame extension.
#' @param variable Either a vector (used as-is), a character column name,
#'   or a one-sided formula referencing a column.
#'
#' @return The variable as a vector.
#'
#' @keywords internal
.get_variable <- function(data, variable) {
  if (length(variable) == nrow(data)) {
    return(variable)
  } else if (methods::is(variable, "character")) {
    if (variable %in% colnames(data)) {
      variable <- data[[variable]]
    } else {
      cli::cli_abort('variable "{variable}" is not found in the data')
    }
  } else if (methods::is(variable, "formula")) {
    # extract column name (class character) from formula
    variable <- data[[all.vars(variable)]]
  } else {
    cli::cli_abort('variable "{variable}" is not found in the data')
  }
  return(variable)
}

#' Extract or generate domain identifier
#'
#' @description This helper function attempts to find a domain identifier column
#' in the data frame. If none is found, it generates sequential indices.
#'
#' @param data A data frame or data frame extension.
#'
#' @return A vector of domain identifiers.
#'
#' @keywords internal
.get_domain_id <- function(data) {
  # Try common domain column names
  domain_cols <- c("area", "domain", "id", "region", "kabupaten", "kota")
  for (col in domain_cols) {
    if (col %in% colnames(data)) {
      return(data[[col]])
    }
  }
  # If not found, generate index
  return(seq_len(nrow(data)))
}
