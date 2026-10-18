# ============================================================================
# Official Statistical Publication Export and Reliability Flagging for fastsae
# ============================================================================

#' Add Official Statistical Reliability Flags
#'
#' @description
#' Categorizes small area estimates into publication reliability tiers based on their
#' relative standard error (RSE, \%). The default thresholds follow the criteria used by
#' BPS-Statistics Indonesia:
#' \itemize{
#'   \item \strong{Reliable} (\eqn{\text{RSE} \le 25\%}): suitable for publication.
#'   \item \strong{Use with caution} (\eqn{25\% < \text{RSE} \le 50\%}): publish with a
#'     cautionary note on sampling variability.
#'   \item \strong{Unreliable} (\eqn{\text{RSE} > 50\%}): suppress or aggregate to a higher level.
#' }
#' Other agencies use other cut-offs (for example 25\% and 50\% at the Australian Bureau of
#' Statistics); set \code{thresholds} accordingly.
#'
#' @param object A \code{fastsae} model object or a data frame containing an RSE column.
#' @param rse_col Character string specifying the name of the RSE column. Default is \code{"rse"}.
#'   If it is missing, the RSE is computed from the \code{mse} and estimate columns.
#' @param thresholds Numeric vector of length 2 with the RSE (\%) cut-offs. Default is \code{c(25, 50)}.
#'
#' @return The estimation table (a data frame) with an additional factor column \code{reliability}.
#'
#' @export
#' @examples
#' library(fastsae)
#' data(mys)
#'
#' fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)
#' df_flagged <- add_reliability_flags(fit_fh)
#' table(df_flagged$reliability)
add_reliability_flags <- function(object, rse_col = "rse", thresholds = c(25, 50)) {
  .check_thresholds(thresholds)
  df <- if (inherits(object, "fastsae")) .sae_table(object) else as.data.frame(object)
  if (is.null(df)) {
    cli::cli_abort("Could not extract estimation data from {.arg object}.")
  }

  if (!rse_col %in% names(df)) {
    est_col <- .sae_est_col(df)
    if (!"mse" %in% names(df) || is.na(est_col)) {
      cli::cli_abort("Column {.val {rse_col}} was not found and cannot be computed (needs {.val mse} and an estimate column).")
    }
    df[[rse_col]] <- (sqrt(df$mse) / abs(df[[est_col]])) * 100
  }

  df$reliability <- .reliability_tier(df[[rse_col]], thresholds)
  df
}

.has_pkg <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

#' Export Small Area Estimation Results to Excel or CSV
#'
#' @description
#' Exports comprehensive estimation tables and model summaries to multi-sheet
#' Excel workbooks (\code{.xlsx}) or comma-separated files (\code{.csv}). Formatted
#' to meet official statistics dissemination standards.
#'
#' The generated Excel workbook contains up to four sheets:
#' \itemize{
#'   \item \strong{Estimates}: Domain IDs, Direct survey estimates, SAE model predictions,
#'     Standard Errors (SE), MSE, RSE (\%), 95\% confidence / credible intervals, and
#'     official Reliability Flags.
#'   \item \strong{Model_Summary}: Model formula, method, convergence status, regression
#'     coefficients (\eqn{\hat{\beta}}, SE, p-values), variance components (\eqn{\sigma_u^2},
#'     spatial \eqn{\rho}, temporal \eqn{\rho_t}, \eqn{\phi}), and goodness-of-fit metrics.
#'   \item \strong{Benchmarked} (Optional): Pre- vs post-calibration values, adjustments,
#'     and percentage shifts when benchmarking calibration is applied.
#'   \item \strong{Benchmarked_Groups} (Optional): the group-level (stage 1) calibration of a
#'     hierarchical benchmark.
#' }
#' With a \code{.csv} file, the estimates are written to \code{file} and the model summary to
#' \code{<file>_summary.csv}; the benchmark sheets are only written to Excel workbooks.
#'
#' @param object A fitted \code{fastsae} object or a \code{fastsae_benchmark} object.
#' @param file Character string specifying the target file path (must end with \code{.xlsx} or \code{.csv}).
#' @param benchmark Optional \code{fastsae_benchmark} object to include calibration results alongside the model.
#' @param thresholds Numeric vector of length 2 with the RSE (\%) cut-offs of the reliability flags
#'   (see \code{\link{add_reliability_flags}}). Default is \code{c(25, 50)}.
#' @param overwrite Logical indicating whether to overwrite an existing file. Default is \code{TRUE}.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns a named list of data frames corresponding to the exported sheets.
#'
#' @export
#' @examples
#' library(fastsae)
#' data(mys)
#'
#' fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)
#' tmp_file <- tempfile(fileext = ".xlsx")
#' export_sae(fit_fh, file = tmp_file)
#' unlink(tmp_file)
export_sae <- function(
  object,
  file,
  benchmark = NULL,
  thresholds = c(25, 50),
  overwrite = TRUE,
  ...
) {
  if (missing(file) || !is.character(file) || length(file) != 1) {
    cli::cli_abort("Please specify a valid file path in {.arg file} (e.g. {.path 'estimates.xlsx'}).")
  }

  .check_thresholds(thresholds)
  ext <- tolower(tools::file_ext(file))
  if (!ext %in% c("xlsx", "csv")) {
    cli::cli_abort("Target file extension must be either {.val xlsx} or {.val csv}.")
  }

  if (file.exists(file) && !overwrite) {
    cli::cli_abort("File {.path {file}} already exists and {.arg overwrite = FALSE}.")
  }

  # Build Sheets
  sheets <- list()

  # 1. Sheet: Estimates
  sheets$Estimates <- .build_estimates_sheet(object, thresholds = thresholds)

  # 2. Sheet: Model Summary
  sheets$Model_Summary <- .build_model_summary_sheet(object)

  # 3. Sheet: Benchmarked (if provided)
  bm_obj <- benchmark %||% (if (inherits(object, "fastsae_benchmark")) object else NULL)
  if (!is.null(bm_obj) && inherits(bm_obj, "fastsae_benchmark")) {
    sheets$Benchmarked <- as.data.frame(bm_obj)
    s1 <- attr(bm_obj, "stage1_summary")
    if (!is.null(s1)) {
      sheets$Benchmarked_Groups <- s1
    }
  }

  # Write output
  if (ext == "xlsx") {
    if (.has_pkg("writexl")) {
      writexl::write_xlsx(sheets, path = file)
      cli::cli_alert_success("Successfully exported SAE results to Excel: {.path {file}}")
    } else if (.has_pkg("openxlsx")) {
      openxlsx::write.xlsx(sheets, file = file)
      cli::cli_alert_success("Successfully exported SAE results to Excel: {.path {file}}")
    } else {
      cli::cli_warn("Neither {.pkg writexl} nor {.pkg openxlsx} is installed. Falling back to CSV.")
      csv_file <- paste0(tools::file_path_sans_ext(file), ".csv")
      utils::write.csv(sheets$Estimates, file = csv_file, row.names = FALSE)
      cli::cli_alert_success("Exported primary estimates sheet to CSV: {.path {csv_file}}")
    }
  } else if (ext == "csv") {
    utils::write.csv(sheets$Estimates, file = file, row.names = FALSE)
    cli::cli_alert_success("Successfully exported SAE estimates to CSV: {.path {file}}")
    if (length(sheets) > 1) {
      summary_file <- paste0(tools::file_path_sans_ext(file), "_summary.csv")
      utils::write.csv(sheets$Model_Summary, file = summary_file, row.names = FALSE)
      cli::cli_alert_info("Exported model summary sheet to: {.path {summary_file}}")
    }
  }

  invisible(sheets)
}

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

#' Construct Estimates Sheet
#' @noRd
.build_estimates_sheet <- function(object, thresholds = c(25, 50)) {
  if (inherits(object, "fastsae_benchmark")) {
    bm_df <- as.data.frame(object)
    adj_pct <- bm_df$adjustment_pct %||% bm_df$rel_adjustment_pct
    res <- data.frame(
      Domain = bm_df$domain,
      Original_Estimate = bm_df$original,
      Benchmarked_Estimate = bm_df$benchmarked,
      Adjustment = bm_df$adjustment,
      Adjustment_Pct = adj_pct,
      stringsAsFactors = FALSE
    )
    return(res)
  }

  df <- if (is.data.frame(object)) object else .sae_table(object)
  if (is.null(df)) cli::cli_abort("Cannot extract estimation table from object.")

  dom_col <- intersect(c("domain", "area", "id"), names(df))[1]
  dom_vals <- df[[dom_col]] %||% seq_len(nrow(df))

  est_col <- .sae_est_col(df)
  if (is.na(est_col)) cli::cli_abort("Cannot find an estimate column in the estimation table.")
  est_val <- df[[est_col]]

  y_val <- if ("y" %in% names(df)) df$y else rep(NA_real_, nrow(df))
  mse_val <- df$mse %||% (if (!is.null(df$sd)) df$sd^2 else rep(NA_real_, nrow(df)))
  se_val <- sqrt(mse_val)

  rse_val <- df$rse
  if (is.null(rse_val) && !all(is.na(mse_val))) {
    rse_val <- (se_val / abs(est_val)) * 100
  }

  ci_l <- df$ci_lower %||% (est_val - 1.96 * se_val)
  ci_u <- df$ci_upper %||% (est_val + 1.96 * se_val)

  rel_flag <- .reliability_tier(rse_val, thresholds)

  res <- data.frame(
    Domain = as.character(dom_vals),
    Direct_Estimate = round(as.numeric(y_val), 5),
    Model_Estimate = round(as.numeric(est_val), 5),
    Standard_Error = round(as.numeric(se_val), 5),
    MSE = round(as.numeric(mse_val), 6),
    RSE_Pct = round(as.numeric(rse_val), 2),
    CI_Lower_95 = round(as.numeric(ci_l), 5),
    CI_Upper_95 = round(as.numeric(ci_u), 5),
    Reliability_Flag = as.character(rel_flag),
    stringsAsFactors = FALSE
  )

  return(res)
}

#' Construct Model Summary Sheet
#' @noRd
.build_model_summary_sheet <- function(object) {
  if (!inherits(object, "fastsae")) {
    return(data.frame(Parameter = "Info", Value = "Custom data frame export", stringsAsFactors = FALSE))
  }

  summary_rows <- list()
  add_row <- function(param, val) {
    summary_rows <<- append(summary_rows, list(data.frame(Parameter = as.character(param), Value = as.character(val), stringsAsFactors = FALSE)))
  }

  add_row("Package", "fastsae")
  add_row("Model Type", object$model %||% "SAE")
  if (!is.null(object$family)) add_row("Distribution Family", object$family)
  if (!is.null(object$method)) add_row("Fitting Method", object$method)
  if (!is.null(object$call)) add_row("Call Formula", deparse(object$call))
  if (!is.null(object$convergence)) add_row("Convergence", ifelse(isTRUE(object$convergence), "Yes", "No"))

  # Variance Components (two values for two-fold models: sigma2_v, sigma2_u)
  rev <- object$random_effect_var
  if (!is.null(rev)) {
    rev_names <- names(rev) %||% if (length(rev) == 1) "sigma2_u" else paste0("component_", seq_along(rev))
    for (i in seq_along(rev)) add_row(paste0("Random Effect Variance (", rev_names[i], ")"), sprintf("%.6f", rev[[i]]))
  }
  if (!is.null(object$rho)) add_row("Spatial Autocorrelation (rho)", sprintf("%.4f", object$rho))
  if (!is.null(object$rho_time)) add_row("Temporal Autocorrelation (rho_t)", sprintf("%.4f", object$rho_time))
  if (!is.null(object$phi)) add_row("Spatial Mixing Fraction (phi)", sprintf("%.4f", object$phi))

  # Goodness of Fit
  if (!is.null(object$goodness)) {
    for (nm in names(object$goodness)) {
      add_row(paste0("Fit Metric (", toupper(nm), ")"), sprintf("%.4f", object$goodness[[nm]]))
    }
  }

  # Coefficients
  if (!is.null(object$estcoef)) {
    for (rn in rownames(object$estcoef)) {
      b_val <- object$estcoef[rn, "beta"]
      se_val <- object$estcoef[rn, "std.error"]
      p_val <- object$estcoef[rn, "pvalue"]
      add_row(paste0("Beta: ", rn), sprintf("%.4f (SE: %.4f, p: %.4e)", b_val, se_val, p_val))
    }
  }

  do.call(rbind, summary_rows)
}
