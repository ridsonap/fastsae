# ============================================================================
# Official Statistical Publication Export and Reliability Flagging for fastsae
# ============================================================================

#' Add Official Statistical Reliability Flags
#'
#' @description
#' Categorizes Small Area Estimates into standard statistical publication reliability tiers
#' based on their Relative Standard Error (RSE / CV \%), following guidelines from
#' official statistical agencies (such as BPS-Statistics Indonesia, Eurostat, and US Census Bureau).
#'
#' Default classification thresholds:
#' \itemize{
#'   \item \strong{Reliable} (\eqn{\text{RSE} < 20\%}): Suitable for unconditional official publication.
#'   \item \strong{Use with Caution} (\eqn{20\% \le \text{RSE} < 30\%}): Usable with cautionary footnotes regarding sampling variability.
#'   \item \strong{Unreliable} (\eqn{\text{RSE} \ge 30\%}): High sampling error; suppression or aggregation to higher geographic levels recommended.
#' }
#'
#' @param object A \code{fastsae} model object or a data frame containing an RSE column.
#' @param rse_col Character string specifying the name of the RSE column. Default is \code{"rse"}.
#' @param thresholds Numeric vector of length 2 defining the RSE (\%) cutoffs. Default is \code{c(20, 30)}.
#'
#' @return The input object or data frame with an additional factor column \code{reliability}.
#'
#' @export
#' @examples
#' library(fastsae)
#' data(mys)
#'
#' fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)
#' df_flagged <- add_reliability_flags(fit_fh)
#' table(df_flagged$reliability)
add_reliability_flags <- function(object, rse_col = "rse", thresholds = c(20, 30)) {
  if (length(thresholds) != 2 || thresholds[1] >= thresholds[2]) {
    cli::cli_abort("{.arg thresholds} must be a numeric vector of length 2 with increasing values.")
  }

  is_model <- inherits(object, "fastsae")
  df <- if (is_model) (object$df_hb %||% object$df_ebp %||% object$df_eblup) else as.data.frame(object)

  if (is.null(df)) {
    cli::cli_abort("Could not extract estimation data from {.arg object}.")
  }

  if (!rse_col %in% names(df)) {
    if ("mse" %in% names(df)) {
      est_col <- intersect(c("hb", "ebp", "eblup", "estimate", "est"), names(df))[1]
      if (is.null(est_col) || is.na(est_col)) cli::cli_abort("Cannot calculate RSE: missing estimate column.")
      df[[rse_col]] <- (sqrt(df$mse) / abs(df[[est_col]])) * 100
    } else {
      cli::cli_abort("Column {.val {rse_col}} was not found in data.")
    }
  }

  labels_tier <- c(
    paste0("Reliable (< ", thresholds[1], "%)"),
    paste0("Use with Caution (", thresholds[1], "-", thresholds[2], "%)"),
    paste0("Unreliable (\u2265 ", thresholds[2], "%)")
  )

  df$reliability <- cut(
    df[[rse_col]],
    breaks = c(-Inf, thresholds[1], thresholds[2], Inf),
    labels = labels_tier,
    right = FALSE
  )

  if (is_model) {
    if (!is.null(object$df_hb)) object$df_hb <- df
    else if (!is.null(object$df_ebp)) object$df_ebp <- df
    else object$df_eblup <- df
    return(df)
  }

  return(df)
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
#' The generated Excel workbook contains up to three dedicated sheets:
#' \itemize{
#'   \item \strong{Estimates}: Domain IDs, Direct survey estimates, SAE model predictions,
#'     Standard Errors (SE), MSE, RSE (\%), 95\% confidence / credible intervals, and
#'     official Reliability Flags.
#'   \item \strong{Model_Summary}: Model formula, method, convergence status, regression
#'     coefficients (\eqn{\hat{\beta}}, SE, p-values), variance components (\eqn{\sigma_u^2},
#'     spatial \eqn{\rho}, temporal \eqn{\rho_t}, \eqn{\phi}), and goodness-of-fit metrics.
#'   \item \strong{Benchmarked} (Optional): Pre- vs post-calibration values, adjustments,
#'     and percentage shifts when benchmarking calibration is applied.
#' }
#'
#' @param object A fitted \code{fastsae} object or a \code{fastsae_benchmark} object.
#' @param file Character string specifying the target file path (must end with \code{.xlsx} or \code{.csv}).
#' @param benchmark Optional \code{fastsae_benchmark} object to include calibration results alongside the model.
#' @param thresholds Numeric vector of length 2 defining the RSE (\%) thresholds for reliability flags. Default is \code{c(20, 30)}.
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
  thresholds = c(20, 30),
  overwrite = TRUE,
  ...
) {
  if (missing(file) || !is.character(file) || length(file) != 1) {
    cli::cli_abort("Please specify a valid file path in {.arg file} (e.g. {.path 'estimates.xlsx'}).")
  }

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
.build_estimates_sheet <- function(object, thresholds = c(20, 30)) {
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

  df <- object$df_hb %||% object$df_ebp %||% object$df_eblup
  if (is.null(df) && is.data.frame(object)) df <- object
  if (is.null(df)) cli::cli_abort("Cannot extract estimation table from object.")

  dom_col <- intersect(c("domain", "area", "id"), names(df))[1]
  dom_vals <- df[[dom_col]] %||% seq_len(nrow(df))

  est_col <- intersect(c("hb", "ebp", "eblup", "estimate", "est"), names(df))[1]
  est_val <- df[[est_col]]

  y_val <- if ("y" %in% names(df)) df$y else rep(NA_real_, nrow(df))
  mse_val <- df$mse %||% (if (!is.null(df$sd)) df$sd^2 else rep(NA_real_, nrow(df)))
  se_val  <- sqrt(mse_val)

  rse_val <- df$rse
  if (is.null(rse_val) && !all(is.na(mse_val))) {
    rse_val <- (se_val / abs(est_val)) * 100
  }

  ci_l <- df$ci_lower %||% (est_val - 1.96 * se_val)
  ci_u <- df$ci_upper %||% (est_val + 1.96 * se_val)

  labels_tier <- c(
    paste0("Reliable (< ", thresholds[1], "%)"),
    paste0("Use with Caution (", thresholds[1], "-", thresholds[2], "%)"),
    paste0("Unreliable (\u2265 ", thresholds[2], "%)")
  )
  rel_flag <- cut(
    rse_val,
    breaks = c(-Inf, thresholds[1], thresholds[2], Inf),
    labels = labels_tier,
    right = FALSE
  )

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
  if (!is.null(object$n_domains)) add_row("Number of Domains", object$n_domains)

  # Variance Components
  if (!is.null(object$random_effect_var)) add_row("Random Effect Variance (sigma2_u)", sprintf("%.6f", object$random_effect_var))
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
      p_val  <- object$estcoef[rn, "pvalue"]
      add_row(paste0("Beta: ", rn), sprintf("%.4f (SE: %.4f, p: %.4e)", b_val, se_val, p_val))
    }
  }

  do.call(rbind, summary_rows)
}
