#' @import ggplot2
#' @importFrom rlang .data
NULL

#' Re-export ggplot2's autoplot generic
#'
#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

#' Autoplot Method for fastsae Objects
#'
#' @description
#' Creates diagnostic and comparison plots for Small Area Estimation (SAE)
#' models fitted with \pkg{fastsae}. This extends the generic
#' \code{\link[ggplot2]{autoplot}} from \pkg{ggplot2}.
#'
#' @param object An object of class \code{fastsae}, or a (named) list of
#'   \code{fastsae} objects.
#' @param type Type of plot to create.
#'   \itemize{
#'     \item For a single \code{fastsae} object: \code{"comparison"},
#'       \code{"mse"}, \code{"estimates"}, or \code{"scatter"}.
#'     \item For a list of \code{fastsae} objects: \code{"comparison"},
#'       \code{"mse"}, or \code{"scatter"}.
#'   }
#' @param ... Additional arguments passed to internal plotting helpers
#'   (e.g. \code{title}) or to \pkg{ggplot2} layers.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' library(fastsae)
#'
#' # Single model plot
#' fit_fh <- eblup_fh(y ~ x1 + x2 + x3, data = mys, vardir = "vardir")
#' autoplot(fit_fh, type = "estimates")
#'
#' # Compare two models
#' fit_sfh <- eblup_sfh(y ~ x1 + x2 + x3, data = mys, vardir = "vardir", W = mys_proxmat)
#' autoplot(list("Fay-Herriot" = fit_fh, "Spatial FH" = fit_sfh), type = "comparison")
#'
#' @importFrom ggplot2 autoplot
#' @export
autoplot.fastsae <- function(object, type = c("comparison", "mse", "estimates", "scatter"), ...) {
  type <- match.arg(type)
  switch(type,
    comparison = .autoplot_single_comparison(object, ...),
    mse        = .autoplot_mse(object, ...),
    estimates  = .autoplot_estimates(object, ...),
    scatter    = .autoplot_scatter_single(object, ...)
  )
}

#' @rdname autoplot.fastsae
#' @export
autoplot.list <- function(object, type = c("comparison", "mse", "scatter"), ...) {
  if (!all(vapply(object, inherits, logical(1), "fastsae"))) {
    # bukan list of fastsae -> lempar ke method/default lain, hentikan eksekusi di sini
    return(NextMethod())
  }
  type <- match.arg(type)
  switch(type,
    comparison = .autoplot_multi_comparison(object, ...),
    mse        = .autoplot_multi_mse(object, ...),
    scatter    = .autoplot_multi_scatter(object, ...)
  )
}

# ------------------------------------------------------------------------------
# Single Model Plots
# ------------------------------------------------------------------------------

#' @noRd
.autoplot_single_comparison <- function(x, title = NULL, ...) {
  df <- x$df_ebp %||% x$df_eblup
  est_col <- if ("ebp" %in% names(df)) "ebp" else "eblup"
  est_label <- if (est_col == "ebp") "EBP" else "EBLUP"

  has_ci <- all(c("ci_lower", "ci_upper") %in% names(df)) && !all(is.na(df$ci_lower))
  has_mse <- "mse" %in% names(df) && !all(is.na(df$mse))
  if (!has_ci && has_mse) {
    df$ci_lower <- df[[est_col]] - 1.96 * sqrt(df$mse)
    df$ci_upper <- df[[est_col]] + 1.96 * sqrt(df$mse)
    has_ci <- TRUE
  }

  if (is.null(title)) {
    title <- paste(est_label, "Estimates with 95% Confidence / Credible Bands")
  }

  p <- ggplot(df, aes(x = .data$domain, y = .data[[est_col]])) +
    geom_point(color = "#2E86AB", size = 2)

  if (has_ci) {
    p <- p + geom_ribbon(aes(ymin = .data$ci_lower, ymax = .data$ci_upper, group = 1),
      fill = "#2E86AB", alpha = 0.2
    )
  }

  p + geom_line(aes(group = 1), color = "#2E86AB", alpha = 0.6) +
    labs(
      title = title,
      x = "Domain",
      y = paste(est_label, "Estimate")
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

#' @noRd
.autoplot_estimates <- function(x, title = NULL, ...) {
  df <- x$df_ebp %||% x$df_eblup

  if (!"y" %in% names(df)) {
    cli::cli_abort(c(
      "Plot type 'estimates' requires direct estimates 'y' in estimation table.",
      "i" = "This plot type is only applicable for area-level models."
    ))
  }

  # Filter out NA values for direct estimates
  df <- df[!is.na(df$y), ]
  est_col <- if ("ebp" %in% names(df)) "ebp" else "eblup"
  est_label <- if (est_col == "ebp") "EBP" else "EBLUP"

  if (is.null(title)) {
    title <- paste(est_label, "Estimates vs Direct Estimates")
  }

  # Determine range for 45-degree line
  min_val <- min(c(df$y, df[[est_col]]), na.rm = TRUE)
  max_val <- max(c(df$y, df[[est_col]]), na.rm = TRUE)

  ggplot(df, aes(x = .data$y, y = .data[[est_col]])) +
    geom_point(color = "#2E86AB", size = 2.5, alpha = 0.7) +
    geom_abline(
      intercept = 0, slope = 1, linetype = "dashed",
      color = "#E94F37", linewidth = 1
    ) +
    geom_hline(
      yintercept = mean(df[[est_col]], na.rm = TRUE),
      linetype = "dotted", color = "gray50"
    ) +
    geom_vline(
      xintercept = mean(df$y, na.rm = TRUE),
      linetype = "dotted", color = "gray50"
    ) +
    expand_limits(x = c(min_val, max_val), y = c(min_val, max_val)) +
    labs(
      title = title,
      x = "Direct Estimate (y)",
      y = paste(est_label, "Estimate")
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
}

#' @noRd
.autoplot_scatter_single <- function(x, title = NULL, ...) {
  cli::cli_warn("Scatter plot requires at least two models. Use a named list of models.")
  .autoplot_estimates(x, title = title, ...)
}

# ------------------------------------------------------------------------------
# Multiple Model Plots
# ------------------------------------------------------------------------------

#' @noRd
.autoplot_multi_comparison <- function(x, title = NULL, ...) {
  if (length(x) < 2) {
    cli::cli_abort("Comparison plot requires at least two models")
  }

  # Combine data from all models
  plot_data <- lapply(names(x), function(name) {
    df <- x[[name]]$df_ebp %||% x[[name]]$df_eblup
    est_col <- if ("ebp" %in% names(df)) "ebp" else "eblup"
    has_ci <- all(c("ci_lower", "ci_upper") %in% names(df)) && !all(is.na(df$ci_lower))
    has_mse <- "mse" %in% names(df) && !all(is.na(df$mse))
    ci_l <- if (has_ci) df$ci_lower else if (has_mse) df[[est_col]] - 1.96 * sqrt(df$mse) else df[[est_col]]
    ci_u <- if (has_ci) df$ci_upper else if (has_mse) df[[est_col]] + 1.96 * sqrt(df$mse) else df[[est_col]]
    data.frame(
      domain = df$domain,
      estimate = df[[est_col]],
      ci_lower = ci_l,
      ci_upper = ci_u,
      model = name,
      stringsAsFactors = FALSE
    )
  })
  plot_data <- do.call(rbind, plot_data)

  if (is.null(title)) {
    title <- "Comparison of Small Area Estimates Across Models"
  }

  # Check if domain names are unique or need model prefix
  plot_data$domain_label <- as.character(plot_data$domain)

  ggplot(plot_data, aes(
    x = .data$domain_label, y = .data$estimate,
    color = .data$model, group = .data$model
  )) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
      position = position_dodge(width = 0.5), width = 0.3, alpha = 0.6
    ) +
    geom_line(position = position_dodge(width = 0.5), alpha = 0.5) +
    labs(
      title = title,
      x = "Domain",
      y = "Small Area Estimate",
      color = "Model"
    ) +
    scale_color_brewer(palette = "Set1") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
}

#' @noRd
.autoplot_multi_mse <- function(x, title = NULL, ...) {
  if (length(x) < 1) {
    cli::cli_abort("MSE plot requires at least one model")
  }

  # Combine MSE data from all models
  plot_data <- lapply(names(x), function(name) {
    df <- x[[name]]$df_ebp %||% x[[name]]$df_eblup
    data.frame(
      domain = as.character(df$domain),
      mse = df$mse,
      model = name,
      stringsAsFactors = FALSE
    )
  })
  plot_data <- do.call(rbind, plot_data)

  if (is.null(title)) {
    title <- "Mean Squared Error Comparison Across Domains"
  }

  # Check if all models have same domains
  same_domains <- all(sapply(x, function(m) {
    df_m <- m$df_ebp %||% m$df_eblup
    df_1 <- x[[1]]$df_ebp %||% x[[1]]$df_eblup
    identical(df_m$domain, df_1$domain)
  }))

  if (same_domains && length(x) > 1) {
    # Faceted plot for same domains
    ggplot(plot_data, aes(x = .data$domain, y = .data$mse, fill = .data$model)) +
      geom_bar(position = position_dodge(preserve = "single"), stat = "identity", alpha = 0.8) +
      labs(
        title = title,
        x = "Domain",
        y = "MSE",
        fill = "Model"
      ) +
      scale_fill_brewer(palette = "Set1") +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
      )
  } else {
    # Separate facets for different domains
    ggplot(plot_data, aes(x = .data$domain, y = .data$mse)) +
      geom_bar(stat = "identity", fill = "#2E86AB", alpha = 0.7) +
      facet_wrap(~model, scales = "free_x", ncol = 1) +
      labs(
        title = title,
        x = "Domain",
        y = "MSE"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
  }
}

#' @noRd
.autoplot_multi_scatter <- function(x, title = NULL, ...) {
  if (length(x) != 2) {
    cli::cli_abort("Scatter plot for model comparison requires exactly two models")
  }

  model_names <- names(x)
  df1_raw <- x[[1]]$df_ebp %||% x[[1]]$df_eblup
  df2_raw <- x[[2]]$df_ebp %||% x[[2]]$df_eblup

  col1_name <- if ("ebp" %in% names(df1_raw)) "ebp" else "eblup"
  col2_name <- if ("ebp" %in% names(df2_raw)) "ebp" else "eblup"
  label1 <- if (col1_name == "ebp") "EBP" else "EBLUP"
  label2 <- if (col2_name == "ebp") "EBP" else "EBLUP"

  df1 <- data.frame(domain = df1_raw$domain, est1 = df1_raw[[col1_name]])
  df2 <- data.frame(domain = df2_raw$domain, est2 = df2_raw[[col2_name]])

  names(df1) <- c("domain", paste0("est_", model_names[1]))
  names(df2) <- c("domain", paste0("est_", model_names[2]))

  plot_data <- merge(df1, df2, by = "domain", all.x = TRUE, all.y = TRUE)
  plot_data <- plot_data[stats::complete.cases(plot_data), ]

  if (nrow(plot_data) == 0) {
    cli::cli_abort("No common domains found between the two models")
  }

  if (is.null(title)) {
    title <- paste("Comparison:", model_names[1], "vs", model_names[2])
  }

  # Calculate correlation
  corr <- stats::cor(plot_data[[2]], plot_data[[3]], use = "complete.obs")

  # Determine range for 45-degree line
  all_vals <- c(plot_data[[2]], plot_data[[3]])
  min_val <- min(all_vals, na.rm = TRUE)
  max_val <- max(all_vals, na.rm = TRUE)

  col1 <- names(plot_data)[2]
  col2 <- names(plot_data)[3]

  ggplot(plot_data, aes(x = .data[[col1]], y = .data[[col2]])) +
    geom_point(color = "#2E86AB", size = 2.5, alpha = 0.7) +
    geom_abline(
      intercept = 0, slope = 1, linetype = "dashed",
      color = "#E94F37", linewidth = 1
    ) +
    expand_limits(x = c(min_val, max_val), y = c(min_val, max_val)) +
    labs(
      title = paste0(title, " (r = ", round(corr, 4), ")"),
      x = paste(label1, "-", model_names[1]),
      y = paste(label2, "-", model_names[2])
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

# ------------------------------------------------------------------------------
# MSE plot for single model (alias for first MSE function)
# ------------------------------------------------------------------------------

#' @noRd
.autoplot_mse <- function(x, title = NULL, ...) {
  df <- x$df_ebp %||% x$df_eblup

  if (!"mse" %in% names(df) || all(is.na(df$mse))) {
    cli::cli_abort(c(
      "MSE estimates are not available in this object.",
      "i" = "Fit the model with MSE computation enabled (e.g. compute_mse = TRUE)."
    ))
  }

  df$domain <- as.character(df$domain)

  if (is.null(title)) {
    title <- "Mean Squared Error by Domain"
  }

  ggplot(df, aes(x = .data$domain, y = .data$mse)) +
    geom_bar(stat = "identity", fill = "#2E86AB", alpha = 0.7) +
    labs(
      title = title,
      x = "Domain",
      y = "MSE"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

#' Autoplot Method for fastsae_diagnose Objects
#'
#' @description
#' Generates diagnostic visual inspections for small area estimation models evaluated
#' via \code{\link{diagnose}}.
#'
#' @param object An object of class \code{"fastsae_diagnose"} returned by \code{\link{diagnose}}.
#' @param type Character string indicating the diagnostic plot type:
#'   \itemize{
#'     \item \code{"all"}: Combined multi-metric inspection (Calibration and RSE reduction).
#'     \item \code{"calibration"}: Direct estimates vs SAE predictions with 1:1 identity line.
#'     \item \code{"rse"}: RSE comparison between direct estimator and model-based predictions.
#'     \item \code{"residuals"}: Standardized residuals versus fitted values.
#'     \item \code{"qq"}: Normal Q-Q plot of standardized residuals.
#'   }
#' @param ... Additional arguments passed to \pkg{ggplot2} layers.
#'
#' @return A \code{ggplot} object.
#'
#' @export
autoplot.fastsae_diagnose <- function(object,
                                      type = c("all", "calibration", "rse", "residuals", "qq"),
                                      ...) {
  type <- match.arg(type)
  df <- object$df_diag
  df_sampled <- df[!is.na(df$direct_y), , drop = FALSE]

  if (type == "calibration") {
    # 1. Calibration Plot (Direct vs SAE Prediction)
    p <- ggplot(df_sampled, aes(x = .data$sae_pred, y = .data$direct_y)) +
      geom_point(color = "#2E86AB", size = 2.5, alpha = 0.8) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#D90429", linewidth = 0.9) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "#2B2D42", linetype = "solid", linewidth = 0.7) +
      labs(
        title = "Bias & Calibration Diagnostic (Brown et al., 2001)",
        subtitle = "Points should fluctuate symmetrically around the red 1:1 line",
        x = "SAE Prediction",
        y = "Direct Estimate"
      ) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold"))
    return(p)

  } else if (type == "rse") {
    # 2. RSE Reduction Plot
    valid_rse <- df_sampled[!is.na(df_sampled$direct_rse) & !is.na(df_sampled$sae_rse), ]
    p <- ggplot(valid_rse, aes(x = .data$direct_rse, y = .data$sae_rse)) +
      geom_point(color = "#3A86FF", size = 2.5, alpha = 0.8) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#E63946", linewidth = 0.9) +
      geom_hline(yintercept = object$precision$rse_threshold, linetype = "dotted", color = "#6C757D") +
      labs(
        title = "Precision & Efficiency Gain (RSE Comparison)",
        subtitle = paste0("Points below the red line demonstrate variance reduction (Threshold: ", object$precision$rse_threshold, "%)"),
        x = "Direct Estimator RSE (%)",
        y = "SAE Model RSE (%)"
      ) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold"))
    return(p)

  } else if (type == "residuals") {
    # 3. Residuals vs Fitted
    p <- ggplot(df_sampled, aes(x = .data$sae_pred, y = .data$std_residuals)) +
      geom_point(color = "#1D3557", size = 2.5, alpha = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "#E63946") +
      geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "#6C757D") +
      labs(
        title = "Residual Diagnostic: Standardized Residuals vs Fitted",
        subtitle = "Check for homoscedasticity and random scatter around zero",
        x = "Fitted SAE Prediction",
        y = "Standardized Residual"
      ) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold"))
    return(p)

  } else if (type == "qq") {
    # 4. Normal Q-Q Plot
    res_clean <- df_sampled$std_residuals[!is.na(df_sampled$std_residuals)]
    qq_df <- data.frame(std_residuals = res_clean)
    p <- ggplot(qq_df, aes(sample = .data$std_residuals)) +
      ggplot2::stat_qq(color = "#457B9D", size = 2.5, alpha = 0.8) +
      ggplot2::stat_qq_line(color = "#E63946", linetype = "dashed", linewidth = 0.8) +
      labs(
        title = "Normal Q-Q Plot of Standardized Residuals",
        subtitle = "Points should adhere closely to the theoretical line",
        x = "Theoretical Quantiles",
        y = "Sample Quantiles"
      ) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold"))
    return(p)

  } else {
    # 5. Combined Dual Inspection (Calibration & RSE side by side)
    df_long <- data.frame(
      Domain = as.character(df_sampled$domain),
      Direct = df_sampled$direct_y,
      SAE = df_sampled$sae_pred,
      RSE_Direct = df_sampled$direct_rse,
      RSE_SAE = df_sampled$sae_rse
    )

    p <- ggplot(df_sampled, aes(x = .data$sae_pred, y = .data$direct_y)) +
      geom_point(aes(color = .data$sae_rse), size = 3, alpha = 0.85) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#D90429", linewidth = 0.9) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "#2B2D42", linewidth = 0.7) +
      ggplot2::scale_color_gradient(low = "#2A9D8F", high = "#E76F51", name = "SAE RSE (%)") +
      labs(
        title = "Diagnostic Overview: Calibration & Estimation Precision",
        subtitle = "Direct vs SAE with 1:1 line (dashed red); color represents SAE RSE",
        x = "SAE Prediction",
        y = "Direct Estimate"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "right"
      )
    return(p)
  }
}
