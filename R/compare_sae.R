# ============================================================================
# Model Comparison and Concordance Evaluation for fastsae (compare_sae)
# ============================================================================

#' Compare Two Small Area Estimation Models
#'
#' @description
#' Provides a comprehensive evaluation and concordance comparison between two
#' fitted Small Area Estimation (SAE) models (e.g. comparing Classical Fay-Herriot
#' vs Spatial Fay-Herriot, Area-level vs Unit-level BHF, or EBLUP vs EBP).
#'
#' Computes empirical agreement and efficiency metrics including:
#' \itemize{
#'   \item Pearson linear correlation (\eqn{r}) and Spearman rank correlation (\eqn{\rho}).
#'   \item Mean Absolute Difference (MAE) and Root Mean Squared Difference (RMSD).
#'   \item Relative efficiency metrics: Mean MSE ratio (\eqn{\text{MSE}_1 / \text{MSE}_2}),
#'     and the proportion of domains where Model 2 achieves greater precision.
#'   \item Average Relative Standard Error (RSE \%) across models and precision gains.
#' }
#'
#' @param model1 A fitted \code{fastsae} model object, or a list containing two \code{fastsae} models.
#' @param model2 A second fitted \code{fastsae} model object. Ignored if \code{model1} is a list.
#' @param names Optional character vector of length 2 specifying descriptive labels
#'   for the models. Default is inferred from model types or call arguments.
#' @param thresholds Numeric vector of length 2 defining the RSE (\%) thresholds for reliability.
#'   Default is \code{c(20, 30)}.
#' @param x An object of class \code{fastsae_comparison} (for \code{print} and \code{plot} methods).
#' @param object An object of class \code{fastsae_comparison} (for \code{summary} and \code{autoplot} methods).
#' @param y Ignored argument for compatibility with the generic \code{plot} method.
#' @param type Character string indicating comparison plot type: \code{"scatter"},
#'   \code{"comparison"}, \code{"difference"}, \code{"mse"}, or \code{"rse"}. Default is \code{"scatter"}.
#' @param title Optional character string specifying a custom plot title.
#' @param ... Additional arguments.
#'
#' @return An S3 object of class \code{fastsae_comparison} containing:
#'   \itemize{
#'     \item \code{metrics}: A data frame of overall concordance and efficiency statistics.
#'     \item \code{data}: A domain-level data frame with aligned estimates, MSE, and RSE values.
#'     \item \code{names}: Vector of the two model names.
#'     \item \code{thresholds}: The RSE thresholds used.
#'   }
#'
#' @seealso \code{\link{autoplot.fastsae_comparison}}, \code{\link{map_sae}}
#'
#' @export
#' @examples
#' library(fastsae)
#' data(mys)
#'
#' # Fit two models
#' fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)
#' fit_sfh <- eblup_sfh(y ~ x1 + x2, vardir = "vardir", data = mys, W = mys_proxmat)
#'
#' # Compare models
#' comp <- compare_sae(fit_fh, fit_sfh, names = c("FH", "Spatial FH"))
#' print(comp)
#' summary(comp)
#'
#' # Plot comparison
#' autoplot(comp, type = "scatter")
#' autoplot(comp, type = "difference")
compare_sae <- function(model1, model2 = NULL, names = NULL, thresholds = c(20, 30), ...) {
  UseMethod("compare_sae")
}

#' @rdname compare_sae
#' @export
compare_sae.default <- function(model1, model2 = NULL, names = NULL, thresholds = c(20, 30), ...) {
  # Support passing a list of two models
  if (is.list(model1) && !inherits(model1, "fastsae") && length(model1) >= 2) {
    if (is.null(names)) {
      names <- names(model1)[1:2]
    }
    model2 <- model1[[2]]
    model1 <- model1[[1]]
  }

  if (is.null(model2)) {
    cli::cli_abort("Please provide two models to compare via {.arg model1} and {.arg model2}.")
  }

  # Infer names if not supplied
  if (is.null(names) || length(names) < 2) {
    name1 <- .get_model_short_name(model1, fallback = "Model 1")
    name2 <- .get_model_short_name(model2, fallback = "Model 2")
    if (name1 == name2) {
      name1 <- paste0(name1, " (1)")
      name2 <- paste0(name2, " (2)")
    }
    names <- c(name1, name2)
  }

  # Extract domain-level data
  df1 <- .extract_domain_data(model1)
  df2 <- .extract_domain_data(model2)

  # Merge on domain ID
  merged <- merge(df1, df2, by = "domain", suffixes = c("_1", "_2"), all = FALSE)

  if (nrow(merged) == 0) {
    cli::cli_abort("No overlapping domains found between {.arg model1} and {.arg model2}.")
  }

  # Compute domain differences
  merged$diff <- merged$estimate_2 - merged$estimate_1
  merged$abs_diff <- abs(merged$diff)

  # MSE ratio (efficiency gain of Model 2 over Model 1)
  has_mse <- !all(is.na(merged$mse_1)) && !all(is.na(merged$mse_2))
  if (has_mse) {
    merged$mse_ratio <- merged$mse_1 / merged$mse_2
  } else {
    merged$mse_ratio <- NA_real_
  }

  # RSE difference
  has_rse <- !all(is.na(merged$rse_1)) && !all(is.na(merged$rse_2))
  if (has_rse) {
    merged$rse_diff <- merged$rse_1 - merged$rse_2
  } else {
    merged$rse_diff <- NA_real_
  }

  # Summary metrics
  cor_p <- stats::cor(merged$estimate_1, merged$estimate_2, use = "complete.obs")
  cor_s <- stats::cor(merged$estimate_1, merged$estimate_2, method = "spearman", use = "complete.obs")
  mae_val <- mean(merged$abs_diff, na.rm = TRUE)
  rmsd_val <- sqrt(mean(merged$diff^2, na.rm = TRUE))

  mean_ratio <- if (has_mse) mean(merged$mse_ratio, na.rm = TRUE) else NA_real_
  med_ratio  <- if (has_mse) stats::median(merged$mse_ratio, na.rm = TRUE) else NA_real_
  n_eff2     <- if (has_mse) sum(merged$mse_2 < merged$mse_1, na.rm = TRUE) else NA_integer_
  pct_eff2   <- if (has_mse) (n_eff2 / nrow(merged)) * 100 else NA_real_

  mean_rse1 <- if (has_rse) mean(merged$rse_1, na.rm = TRUE) else NA_real_
  mean_rse2 <- if (has_rse) mean(merged$rse_2, na.rm = TRUE) else NA_real_
  rse_gain  <- if (has_rse) mean_rse1 - mean_rse2 else NA_real_

  metrics_df <- data.frame(
    Metric = c(
      "Matched Domains",
      "Pearson Correlation (r)",
      "Spearman Correlation (rho)",
      "Mean Absolute Difference (MAE)",
      "Root Mean Squared Difference (RMSD)",
      paste0("Mean MSE Ratio (", names[1], " / ", names[2], ")"),
      paste0("Median MSE Ratio (", names[1], " / ", names[2], ")"),
      paste0("Domains where ", names[2], " has lower MSE"),
      paste0("Percentage of domains where ", names[2], " is more efficient (%)"),
      paste0("Mean RSE % (", names[1], ")"),
      paste0("Mean RSE % (", names[2], ")"),
      paste0("Average RSE Reduction (%)")
    ),
    Value = c(
      as.character(nrow(merged)),
      sprintf("%.4f", cor_p),
      sprintf("%.4f", cor_s),
      sprintf("%.4f", mae_val),
      sprintf("%.4f", rmsd_val),
      if (!is.na(mean_ratio)) sprintf("%.4f", mean_ratio) else "N/A",
      if (!is.na(med_ratio)) sprintf("%.4f", med_ratio) else "N/A",
      if (!is.na(n_eff2)) sprintf("%d / %d", n_eff2, nrow(merged)) else "N/A",
      if (!is.na(pct_eff2)) sprintf("%.1f%%", pct_eff2) else "N/A",
      if (!is.na(mean_rse1)) sprintf("%.2f%%", mean_rse1) else "N/A",
      if (!is.na(mean_rse2)) sprintf("%.2f%%", mean_rse2) else "N/A",
      if (!is.na(rse_gain)) sprintf("%.2f%%", rse_gain) else "N/A"
    ),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      metrics = metrics_df,
      data = merged,
      names = names,
      thresholds = thresholds,
      model1 = model1,
      model2 = model2
    ),
    class = "fastsae_comparison"
  )
}

#' @rdname compare_sae
#' @export
print.fastsae_comparison <- function(x, ...) {
  cli::cli_h1("Small Area Estimation Model Concordance & Comparison")
  cli::cli_text("{.strong Model 1}: {x$names[1]}")
  cli::cli_text("{.strong Model 2}: {x$names[2]}")
  cli::cli_text("{.strong Matched Domains}: {nrow(x$data)}")
  cli::cli_rule()

  cat("\nKey Concordance & Efficiency Metrics:\n")
  print(x$metrics, row.names = FALSE)

  cat("\n")
  invisible(x)
}

#' @rdname compare_sae
#' @export
summary.fastsae_comparison <- function(object, ...) {
  print(object, ...)

  cat("Distribution of Domain Point Estimates:\n")
  est_mat <- cbind(
    stats::quantile(object$data$estimate_1, na.rm = TRUE),
    stats::quantile(object$data$estimate_2, na.rm = TRUE),
    stats::quantile(object$data$diff, na.rm = TRUE)
  )
  colnames(est_mat) <- c(object$names[1], object$names[2], "Difference")
  print(round(est_mat, 4))

  if (!all(is.na(object$data$rse_1)) && !all(is.na(object$data$rse_2))) {
    cat("\nDistribution of Relative Standard Error (RSE %):\n")
    rse_mat <- cbind(
      stats::quantile(object$data$rse_1, na.rm = TRUE),
      stats::quantile(object$data$rse_2, na.rm = TRUE),
      stats::quantile(object$data$rse_diff, na.rm = TRUE)
    )
    colnames(rse_mat) <- c(object$names[1], object$names[2], "RSE Diff (1 - 2)")
    print(round(rse_mat, 2))
  }

  cat("\n")
  invisible(object)
}

#' @rdname compare_sae
#' @export
plot.fastsae_comparison <- function(x, y = NULL, ...) {
  ggplot2::autoplot(x, ...)
}

#' @rdname compare_sae
#' @export
autoplot.fastsae_comparison <- function(
  object,
  type = c("scatter", "comparison", "difference", "mse", "rse"),
  title = NULL,
  ...
) {
  type <- match.arg(type)
  df <- object$data
  n1 <- object$names[1]
  n2 <- object$names[2]

  if (type == "scatter") {
    r_val <- stats::cor(df$estimate_1, df$estimate_2, use = "complete.obs")
    lims <- range(c(df$estimate_1, df$estimate_2), na.rm = TRUE)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$estimate_1, y = .data$estimate_2)) +
      ggplot2::geom_point(ggplot2::aes(color = .data$diff), size = 2.5, alpha = 0.8) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#d7191c", linewidth = 0.8) +
      ggplot2::scale_color_gradient2(
        low = "#2b83ba", mid = "#ffffbf", high = "#d7191c", midpoint = 0,
        name = "Difference"
      ) +
      ggplot2::coord_fixed(xlim = lims, ylim = lims) +
      ggplot2::labs(
        title = title %||% paste("Model Concordance:", n1, "vs", n2),
        subtitle = paste0("Pearson r = ", round(r_val, 4), " (Dashed red line indicates 1:1 identity)"),
        x = paste(n1, "Estimate"),
        y = paste(n2, "Estimate")
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

    return(p)

  } else if (type == "difference") {
    df$domain_fac <- factor(df$domain, levels = df$domain[order(df$diff)])

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$domain_fac, y = .data$diff)) +
      ggplot2::geom_segment(
        ggplot2::aes(x = .data$domain_fac, xend = .data$domain_fac, y = 0, yend = .data$diff),
        color = "gray70"
      ) +
      ggplot2::geom_point(ggplot2::aes(color = .data$diff > 0), size = 2.5) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      ggplot2::scale_color_manual(
        values = c("TRUE" = "#2b83ba", "FALSE" = "#d7191c"),
        labels = c("TRUE" = paste(n2, ">", n1), "FALSE" = paste(n2, "<", n1)),
        name = "Shift"
      ) +
      ggplot2::labs(
        title = title %||% paste("Domain Estimate Shift (", n2, "minus", n1, ")"),
        subtitle = "Ordered by magnitude of adjustment",
        x = "Domain",
        y = "Difference"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "bottom"
      )

    return(p)

  } else if (type == "mse") {
    long_mse <- rbind(
      data.frame(domain = df$domain, mse = df$mse_1, Model = n1),
      data.frame(domain = df$domain, mse = df$mse_2, Model = n2)
    )
    long_mse$domain <- factor(long_mse$domain, levels = unique(df$domain))

    p <- ggplot2::ggplot(long_mse, ggplot2::aes(x = .data$domain, y = .data$mse, fill = .data$Model)) +
      ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
      ggplot2::scale_fill_brewer(palette = "Set1") +
      ggplot2::labs(
        title = title %||% "Mean Squared Error (MSE) Comparison",
        subtitle = "Lower MSE denotes superior estimation precision",
        x = "Domain",
        y = "MSE"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "bottom"
      )

    return(p)

  } else if (type == "rse") {
    long_rse <- rbind(
      data.frame(domain = df$domain, rse = df$rse_1, Model = n1),
      data.frame(domain = df$domain, rse = df$rse_2, Model = n2)
    )
    long_rse$domain <- factor(long_rse$domain, levels = unique(df$domain))

    p <- ggplot2::ggplot(long_rse, ggplot2::aes(x = .data$domain, y = .data$rse, color = .data$Model, group = .data$Model)) +
      ggplot2::geom_line(alpha = 0.7) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_hline(yintercept = object$thresholds[1], linetype = "dashed", color = "#2A9D8F") +
      ggplot2::geom_hline(yintercept = object$thresholds[2], linetype = "dashed", color = "#E76F51") +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::labs(
        title = title %||% "Relative Standard Error (RSE %) Comparison",
        subtitle = paste0("Official statistics reliability thresholds: < ", object$thresholds[1], "% (Reliable) and \u2265 ", object$thresholds[2], "% (Unreliable)"),
        x = "Domain",
        y = "RSE (%)"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "bottom"
      )

    return(p)

  } else {
    # Default comparison: Overlaid point estimates with lines
    long_est <- rbind(
      data.frame(domain = df$domain, est = df$estimate_1, Model = n1),
      data.frame(domain = df$domain, est = df$estimate_2, Model = n2)
    )
    long_est$domain <- factor(long_est$domain, levels = unique(df$domain))

    p <- ggplot2::ggplot(long_est, ggplot2::aes(x = .data$domain, y = .data$est, color = .data$Model, group = .data$Model)) +
      ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.4), size = 2) +
      ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.4), alpha = 0.6) +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::labs(
        title = title %||% "Domain Estimation Profile Across Models",
        x = "Domain",
        y = "Point Estimate"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "bottom"
      )

    return(p)
  }
}

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

#' Extract standardized domain-level data frame
#' @noRd
.extract_domain_data <- function(model) {
  if (inherits(model, "fastsae")) {
    df <- model$df_ebp %||% model$df_eblup
    if (is.null(df)) cli::cli_abort("Model does not contain estimation data frame.")
    dom_col <- intersect(c("domain", "area", "id"), names(df))[1]
    dom_vals <- df[[dom_col]] %||% seq_len(nrow(df))
    est_val <- df$ebp %||% df$eblup %||% df$est
    mse_val <- df$mse %||% (if (!is.null(df$sd)) df$sd^2 else NA_real_)
    rse_val <- df$rse %||% (if (!all(is.na(mse_val))) (sqrt(mse_val) / abs(est_val)) * 100 else NA_real_)

    return(data.frame(
      domain = as.character(dom_vals),
      estimate = as.numeric(est_val),
      mse = as.numeric(mse_val),
      rse = as.numeric(rse_val),
      stringsAsFactors = FALSE
    ))
  } else if (is.data.frame(model)) {
    dom_col <- intersect(c("domain", "area", "id", "code"), names(model))[1]
    est_col <- intersect(c("estimate", "ebp", "eblup", "est", "y_hat"), names(model))[1]
    if (is.null(est_col)) cli::cli_abort("Data frame must contain an estimate column.")

    dom_vals <- if (!is.null(dom_col)) model[[dom_col]] else seq_len(nrow(model))
    est_val <- model[[est_col]]
    mse_val <- model$mse %||% rep(NA_real_, nrow(model))
    rse_val <- model$rse %||% (if (!all(is.na(mse_val))) (sqrt(mse_val) / abs(est_val)) * 100 else rep(NA_real_, nrow(model)))

    return(data.frame(
      domain = as.character(dom_vals),
      estimate = as.numeric(est_val),
      mse = as.numeric(mse_val),
      rse = as.numeric(rse_val),
      stringsAsFactors = FALSE
    ))
  } else {
    cli::cli_abort("Unsupported model object class: {.cls {class(model)}}.")
  }
}

#' Infer short model name
#' @noRd
.get_model_short_name <- function(model, fallback = "Model") {
  if (inherits(model, "fastsae")) {
    if (!is.null(model$model)) {
      m <- switch(model$model,
        "FH" = "FH",
        "SFH" = "Spatial FH",
        "ST" = "Spatio-Temporal FH",
        "BHF" = "BHF",
        "EBP" = paste0("EBP-", toupper(model$family %||% "Beta")),
        model$model
      )
      return(m)
    }
  }
  return(fallback)
}
