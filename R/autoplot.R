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
  df <- x$df_eblup
  has_mse <- "mse" %in% names(df) && !all(is.na(df$mse))
  if (has_mse) {
    df$ci_lower <- df$eblup - 1.96 * sqrt(df$mse)
    df$ci_upper <- df$eblup + 1.96 * sqrt(df$mse)
  }

  if (is.null(title)) {
    title <- "EBLUP Estimates with 95% Confidence Bands"
  }

  p <- ggplot(df, aes(x = .data$domain, y = .data$eblup)) +
    geom_point(color = "#2E86AB", size = 2)

  if (has_mse) {
    p <- p + geom_ribbon(aes(ymin = .data$ci_lower, ymax = .data$ci_upper, group = 1),
      fill = "#2E86AB", alpha = 0.2
    )
  }

  p + geom_line(aes(group = 1), color = "#2E86AB", alpha = 0.6) +
    labs(
      title = title,
      x = "Domain",
      y = "EBLUP Estimate"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

#' @noRd
.autoplot_estimates <- function(x, title = NULL, ...) {
  df <- x$df_eblup

  if (!"y" %in% names(df)) {
    cli::cli_abort(c(
      "Plot type 'estimates' requires direct estimates 'y' in df_eblup.",
      "i" = "This plot type is only applicable for area-level models (FH, SFH, STFH)."
    ))
  }

  # Filter out NA values for direct estimates
  df <- df[!is.na(df$y), ]

  if (is.null(title)) {
    title <- "EBLUP Estimates vs Direct Estimates"
  }

  # Determine range for 45-degree line
  min_val <- min(c(df$y, df$eblup), na.rm = TRUE)
  max_val <- max(c(df$y, df$eblup), na.rm = TRUE)

  ggplot(df, aes(x = .data$y, y = .data$eblup)) +
    geom_point(color = "#2E86AB", size = 2.5, alpha = 0.7) +
    geom_abline(
      intercept = 0, slope = 1, linetype = "dashed",
      color = "#E94F37", linewidth = 1
    ) +
    geom_hline(
      yintercept = mean(df$eblup, na.rm = TRUE),
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
      y = "EBLUP Estimate"
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
    df <- x[[name]]$df_eblup
    has_mse <- "mse" %in% names(df) && !all(is.na(df$mse))
    ci_l <- if (has_mse) df$eblup - 1.96 * sqrt(df$mse) else df$eblup
    ci_u <- if (has_mse) df$eblup + 1.96 * sqrt(df$mse) else df$eblup
    data.frame(
      domain = df$domain,
      eblup = df$eblup,
      ci_lower = ci_l,
      ci_upper = ci_u,
      model = name,
      stringsAsFactors = FALSE
    )
  })
  plot_data <- do.call(rbind, plot_data)

  if (is.null(title)) {
    title <- "Comparison of EBLUP Estimates Across Models"
  }

  # Check if domain names are unique or need model prefix
  plot_data$domain_label <- as.character(plot_data$domain)

  ggplot(plot_data, aes(
    x = .data$domain_label, y = .data$eblup,
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
      y = "EBLUP Estimate",
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
    df <- x[[name]]$df_eblup
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
  same_domains <- all(sapply(x, function(m) identical(m$df_eblup$domain, x[[1]]$df_eblup$domain)))

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
  df1 <- x[[1]]$df_eblup[, c("domain", "eblup")]
  df2 <- x[[2]]$df_eblup[, c("domain", "eblup")]

  names(df1) <- c("domain", paste0("eblup_", model_names[1]))
  names(df2) <- c("domain", paste0("eblup_", model_names[2]))

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
      x = paste("EBLUP -", model_names[1]),
      y = paste("EBLUP -", model_names[2])
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
  df <- x$df_eblup

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
