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
#'     \item For a single \code{fastsae} object: \code{"comparison"}, \code{"ribbon"},
#'       \code{"mse"}, \code{"estimates"}, \code{"scatter"}, \code{"rse"}, or \code{"map"};
#'       for a fit with \code{\link{jax_engine}()} also \code{"trace"} and \code{"rank"},
#'       the MCMC convergence plots (see Details).
#'     \item For a list of \code{fastsae} objects: \code{"comparison"},
#'       \code{"mse"}, \code{"scatter"}, or \code{"map"}.
#'   }
#' @param ... Additional arguments passed to internal plotting helpers
#'   (e.g. \code{title}) or to \pkg{ggplot2} layers. For \code{"trace"} and \code{"rank"}:
#'   \code{pars}, the parameters to show (names in \code{fit$mcmc_diagnostics$parameter},
#'   such as \code{"sigma2_u"}, a coefficient name or \code{"hb[5]"} for row 5 of
#'   \code{df_hb}), and for \code{"rank"} \code{bins} (default 20).
#'
#' @details
#' \strong{MCMC convergence plots.} \code{"trace"} draws each chain's sampled values against
#' the iteration; \code{"rank"} the rank plot of Vehtari et al. (2021): the draws of all chains
#' are ranked together and each chain's ranks are counted in \code{bins} bins. If the chains
#' sample the same distribution each chain's counts lie near the dashed line (an equal share);
#' a chain stuck in one region has too many low or high ranks. The facet labels give the
#' rank-normalized R-hat and bulk ESS. By default the plots show the variance and other
#' hyperparameters, up to eight coefficients and the area or domain with the largest R-hat;
#' all R-hats and ESSs are in \code{fit$mcmc_diagnostics}. The numbers, not the plots, are the
#' check (R-hat < 1.01 and ESS > 100 per chain); the plots show why a check fails.
#'
#' @references Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & B\enc{ü}{u}rkner, P.-C.
#'   (2021). Rank-normalization, folding, and localization: An improved R-hat for assessing
#'   convergence of MCMC. \emph{Bayesian Analysis}, 16(2), 667-718.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' library(fastsae)
#'
#' # Single model plot
#' fit_fh <- eblup_fh(y ~ x1 + x2 + x3, data = mys, vardir = "vardir")
#' autoplot(fit_fh, type = "estimates")
#' autoplot(fit_fh, type = "rse")
#'
#' # Compare two models
#' fit_sfh <- eblup_sfh(y ~ x1 + x2 + x3, data = mys, vardir = "vardir", W = mys_proxmat)
#' autoplot(list("Fay-Herriot" = fit_fh, "Spatial FH" = fit_sfh), type = "comparison")
#'
#' @importFrom ggplot2 autoplot
#' @export
autoplot.fastsae <- function(
  object,
  type = c("comparison", "ribbon", "mse", "estimates", "scatter", "rse", "map", "trace", "rank"),
  ...
) {
  type <- match.arg(type)
  if (type == "map") {
    return(map_sae(object, ...))
  }
  if (type %in% c("trace", "rank")) {
    return(.autoplot_mcmc(object, type, ...))
  }
  switch(type,
    comparison = .autoplot_single_comparison(object, ...),
    ribbon     = .autoplot_single_comparison(object, ...),
    mse        = .autoplot_mse(object, ...),
    estimates  = .autoplot_estimates(object, ...),
    scatter    = .autoplot_scatter_single(object, ...),
    rse        = .autoplot_rse(object, ...)
  )
}

#' @rdname autoplot.fastsae
#' @export
autoplot.list <- function(object, type = c("comparison", "mse", "scatter", "map"), ...) {
  if (!all(vapply(object, inherits, logical(1), "fastsae"))) {
    # bukan list of fastsae -> lempar ke method/default lain, hentikan eksekusi di sini
    return(NextMethod())
  }
  type <- match.arg(type)
  if (length(object) > 0 && (is.null(names(object)) || any(!nzchar(names(object))))) {
    names(object) <- paste("Model", seq_along(object))
  }
  if (type == "map") {
    return(map_sae(object, ...))
  }
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
  df <- .sae_table(x)
  est_col <- .sae_est_col(df)
  est_label <- if (est_col == "hb") "HB" else "EBLUP"

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
  df <- .sae_table(x)

  if (!"y" %in% names(df)) {
    cli::cli_abort(c(
      "Plot type 'estimates' requires direct estimates 'y' in estimation table.",
      "i" = "This plot type is only applicable for area-level models."
    ))
  }

  # Filter out NA values for direct estimates
  df <- df[!is.na(df$y), ]
  est_col <- .sae_est_col(df)
  est_label <- if (est_col == "hb") "HB" else "EBLUP"

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

#' @noRd
.autoplot_rse <- function(x, title = NULL, thresholds = c(25, 50), ...) {
  .check_thresholds(thresholds)
  df <- .sae_table(x)
  est_col <- .sae_est_col(df)

  # Calculate SAE RSE if not already present
  if (!"rse" %in% names(df)) {
    mse_val <- df$mse %||% (if (!is.null(df$sd)) df$sd^2 else NA_real_)
    df$rse <- (sqrt(mse_val) / abs(df[[est_col]])) * 100
  }

  # Check if direct vardir/rse exists
  has_direct_rse <- FALSE
  if ("vardir" %in% names(df) && "y" %in% names(df)) {
    df$direct_rse <- (sqrt(df$vardir) / abs(df$y)) * 100
    has_direct_rse <- TRUE
  } else if ("direct_rse" %in% names(df)) {
    has_direct_rse <- TRUE
  }

  lab <- .row_labels(df)
  df$domain_char <- factor(lab, levels = lab)

  if (has_direct_rse) {
    long_rse <- rbind(
      data.frame(domain = df$domain_char, rse = df$direct_rse, Type = "Direct Survey"),
      data.frame(domain = df$domain_char, rse = df$rse, Type = "SAE Model")
    )
    p <- ggplot(long_rse, aes(x = .data$domain, y = .data$rse, color = .data$Type, group = .data$Type)) +
      geom_point(size = 2) +
      geom_line(alpha = 0.7) +
      ggplot2::scale_color_manual(values = c("Direct Survey" = "#E63946", "SAE Model" = "#1D3557"))
  } else {
    p <- ggplot(df, aes(x = .data$domain_char, y = .data$rse, group = 1)) +
      geom_point(color = "#1D3557", size = 2) +
      geom_line(color = "#1D3557", alpha = 0.7)
  }

  p <- p +
    geom_hline(yintercept = thresholds[1], linetype = "dashed", color = "#2A9D8F", linewidth = 0.8) +
    geom_hline(yintercept = thresholds[2], linetype = "dashed", color = "#E76F51", linewidth = 0.8) +
    ggplot2::annotate("text",
      x = 1, y = thresholds[1] - 1.5, label = paste0("Reliable (\u2264 ", thresholds[1], "%)"),
      hjust = 0, size = 3, color = "#2A9D8F", fontface = "italic"
    ) +
    ggplot2::annotate("text",
      x = 1, y = thresholds[2] + 1.5, label = paste0("Unreliable (> ", thresholds[2], "%)"),
      hjust = 0, size = 3, color = "#E76F51", fontface = "italic"
    ) +
    labs(
      title = title %||% "Relative Standard Error (RSE %) across Domains",
      subtitle = paste0("Reliability thresholds: RSE \u2264 ", thresholds[1], "% (reliable) and > ", thresholds[2], "% (unreliable)"),
      x = "Domain",
      y = "RSE (%)",
      color = "Method"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  return(p)
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
    df <- .sae_table(x[[name]])
    est_col <- .sae_est_col(df)
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
    df <- .sae_table(x[[name]])
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
    df_m <- .sae_table(m)
    df_1 <- .sae_table(x[[1]])
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
  df1_raw <- .sae_table(x[[1]])
  df2_raw <- .sae_table(x[[2]])

  col1_name <- .sae_est_col(df1_raw)
  col2_name <- .sae_est_col(df2_raw)
  label1 <- if (col1_name == "hb") "HB" else "EBLUP"
  label2 <- if (col2_name == "hb") "HB" else "EBLUP"

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
  df <- .sae_table(x)

  if (!"mse" %in% names(df) || all(is.na(df$mse))) {
    cli::cli_abort(c(
      "MSE estimates are not available in this object.",
      "i" = "Fit the model with an MSE estimator (an {.arg mse_type} other than \"none\")."
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
                                      type = c("all", "calibration", "rse", "residuals", "qq", "pit", "cpo"),
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
  } else if (type == "pit") {
    # 5. PIT Histogram (Probability Integral Transform)
    if (!"pit" %in% names(df) || all(is.na(df$pit))) {
      cli::cli_abort("PIT (Probability Integral Transform) values not available in diagnostic object.")
    }
    pit_clean <- df$pit[!is.na(df$pit) & df$pit >= 0 & df$pit <= 1]
    pit_df <- data.frame(pit = pit_clean)
    p <- ggplot(pit_df, aes(x = .data$pit)) +
      ggplot2::geom_histogram(aes(y = ggplot2::after_stat(density)), bins = 10, fill = "#457B9D", color = "white", alpha = 0.8) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "#D90429", linewidth = 0.9) +
      labs(
        title = "Bayesian Calibration: Probability Integral Transform (PIT)",
        subtitle = "A well-calibrated predictive distribution adheres to Uniform(0,1) density = 1 (dashed red line)",
        x = "PIT Value",
        y = "Density"
      ) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold"))
    return(p)
  } else if (type == "cpo") {
    # 6. CPO Index Plot
    if (!"cpo" %in% names(df) || all(is.na(df$cpo))) {
      cli::cli_abort("CPO (Conditional Predictive Ordinate) values not available in diagnostic object.")
    }
    df_cpo <- df[!is.na(df$cpo), , drop = FALSE]
    lab <- .row_labels(df_cpo)
    df_cpo$domain_factor <- factor(lab, levels = lab)
    p <- ggplot(df_cpo, aes(x = .data$domain_factor, y = .data$cpo)) +
      ggplot2::geom_segment(aes(x = .data$domain_factor, xend = .data$domain_factor, y = 0, yend = .data$cpo), color = "#A8DADC") +
      ggplot2::geom_point(color = "#2E86AB", size = 2.5) +
      labs(
        title = "Leave-One-Out Cross-Validation: CPO by Domain",
        subtitle = "Lower CPO values identify potential area outliers / model fit deficiencies",
        x = "Domain",
        y = "CPO Value"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(face = "bold")
      )
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

#' Per-chain draws for the MCMC plots: iterations x chains x parameters
#'
#' The stored trace of the key parameters, plus rows of `draws` (chain-major,
#' as NumPyro concatenates the chains) for parameters `hb[j]` not in it.
#' @noRd
.mcmc_trace <- function(object, pars) {
  trace <- object$trace
  if (is.null(trace)) {
    cli::cli_abort(c(
      "Trace and rank plots need the MCMC draws of a fit with {.code engine = jax_engine()}.",
      "i" = "This fit used the {.val {object$engine %||% 'frequentist'}} engine."
    ))
  }
  if (is.null(pars)) {
    return(trace)
  }
  avail <- dimnames(trace)$parameter
  extra <- setdiff(pars, avail)
  rows <- suppressWarnings(as.integer(sub("^hb\\[([0-9]+)\\]$", "\\1", extra)))
  bad <- extra[is.na(rows) | !grepl("^hb\\[", extra)]
  if (length(bad) > 0) {
    cli::cli_abort(c(
      "Unknown parameter{?s} {.val {bad}}.",
      "i" = "Available: {.val {avail}}, and {.code hb[j]} for row j of {.field df_hb}."
    ))
  }
  if (length(extra) > 0) {
    draws <- object$draws
    if (is.null(draws) || any(rows > nrow(draws))) {
      cli::cli_abort("{.val {extra}} need{?s/} the draws of the fit ({.code jax_engine(keep_draws = TRUE)}) and a valid row of {.field df_hb}.")
    }
    n_iter <- dim(trace)[1]
    n_chain <- dim(trace)[2]
    add <- array(t(draws[rows, , drop = FALSE]), dim = c(n_iter, n_chain, length(rows)))
    trace <- array(c(trace, add),
      dim = c(n_iter, n_chain, dim(trace)[3] + length(rows)),
      dimnames = list(iteration = NULL, chain = seq_len(n_chain), parameter = c(avail, extra))
    )
  }
  trace[, , pars, drop = FALSE]
}

#' Trace and rank plots of a jax_engine() fit
#' @noRd
.autoplot_mcmc <- function(object, type, pars = NULL, bins = 20, title = NULL, ...) {
  trace <- .mcmc_trace(object, pars)
  dims <- dim(trace)
  n_iter <- dims[1]
  n_chain <- dims[2]
  par_names <- dimnames(trace)$parameter
  df <- data.frame(
    iteration = rep(seq_len(n_iter), times = n_chain * dims[3]),
    chain = factor(rep(rep(seq_len(n_chain), each = n_iter), times = dims[3])),
    parameter = rep(par_names, each = n_iter * n_chain),
    value = as.vector(trace)
  )
  # facet labels with the diagnostics of each parameter
  tab <- object$mcmc_diagnostics
  lab <- par_names
  if (!is.null(tab)) {
    i <- match(par_names, tab$parameter)
    lab <- ifelse(is.na(i), par_names,
      sprintf("%s\nR-hat %.3f, ESS %.0f", par_names, tab$rhat[i], tab$ess_bulk[i])
    )
  }
  df$parameter <- factor(df$parameter, levels = par_names, labels = lab)

  if (type == "trace") {
    p <- ggplot(df, aes(x = .data$iteration, y = .data$value, colour = .data$chain)) +
      geom_line(linewidth = 0.25, alpha = 0.8) +
      facet_wrap(~parameter, scales = "free_y") +
      labs(
        x = "Iteration (after warmup)", y = NULL, colour = "Chain",
        title = title %||% "MCMC trace plot"
      )
  } else {
    df$rank <- stats::ave(df$value, df$parameter, FUN = function(v) rank(v, ties.method = "average"))
    n_total <- n_iter * n_chain
    df$bin <- pmin(floor((df$rank - 1) / n_total * bins), bins - 1)
    counts <- stats::aggregate(list(count = df$rank), df[c("parameter", "chain", "bin")], length)
    full <- expand.grid(
      parameter = levels(df$parameter), chain = levels(df$chain), bin = 0:(bins - 1),
      stringsAsFactors = FALSE
    )
    counts$parameter <- as.character(counts$parameter)
    counts$chain <- as.character(counts$chain)
    full <- merge(full, counts, all.x = TRUE)
    full$count[is.na(full$count)] <- 0
    full$parameter <- factor(full$parameter, levels = levels(df$parameter))
    full$chain <- factor(full$chain, levels = levels(df$chain))
    # step from the left edge of each bin; the last bin is closed on the right
    last <- full[full$bin == bins - 1, ]
    last$bin <- bins
    full <- rbind(full, last)
    p <- ggplot(full, aes(x = .data$bin / bins, y = .data$count, colour = .data$chain)) +
      geom_hline(yintercept = n_iter / bins, linetype = "dashed", colour = "grey40") +
      geom_step(linewidth = 0.4) +
      facet_wrap(~parameter, scales = "free_y") +
      scale_x_continuous(labels = function(x) paste0(round(100 * x), "%")) +
      labs(
        x = "Rank (share of all draws)", y = "Draws per bin", colour = "Chain",
        title = title %||% "MCMC rank plot",
        subtitle = "Chains that sample the same distribution stay near the dashed line"
      )
  }
  p + theme_minimal(base_size = 11) + theme(legend.position = "bottom")
}
