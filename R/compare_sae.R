# ============================================================================
# Model Comparison and Concordance Evaluation for fastsae (compare_sae)
# ============================================================================

#' Compare Two Small Area Estimation Models
#'
#' @description
#' Provides a comprehensive evaluation and concordance comparison between two
#' fitted Small Area Estimation (SAE) models (e.g. comparing Classical Fay-Herriot
#' vs Spatial Fay-Herriot, Area-level vs Unit-level BHF, or EBLUP vs HB).
#'
#' Computes agreement, ranking and efficiency metrics:
#' \itemize{
#'   \item Agreement of the estimates (Bland and Altman 1986): the mean difference
#'     \eqn{\bar d} of Model 2 minus Model 1, the 95\% limits of agreement
#'     \eqn{\bar d \pm 1.96\, s_d}, the mean absolute difference (MAE) and the root mean
#'     squared difference (RMSD), and Lin's (1989) concordance correlation coefficient
#'     \eqn{\rho_c = 2 s_{12} / (s_1^2 + s_2^2 + (\bar\theta_1 - \bar\theta_2)^2)} (moments
#'     with divisor \eqn{n}), which is 1 only when the estimates lie on the identity line.
#'   \item Differences relative to the uncertainty: \eqn{z_d = (\hat\theta_{2d} -
#'     \hat\theta_{1d}) / \sqrt{\text{MSE}_{1d}}}, summarised by its median absolute value and
#'     the share of domains with \eqn{|z_d| > 1}. This is descriptive, not a test: both
#'     estimators use the same data, so their difference has a smaller variance than either.
#'   \item Ranking of the domains: Spearman's rank correlation, the mean absolute change of
#'     rank, and the share of the \code{top} highest and lowest domains of Model 1 that are
#'     also among the highest and lowest of Model 2 (the domains a programme would target).
#'   \item Efficiency: mean and median MSE ratio (\eqn{\text{MSE}_1 / \text{MSE}_2}), the
#'     share of domains where Model 2 has the lower MSE, and the mean RSE (\%) of each model.
#' }
#'
#' Pearson's correlation is not reported: two reasonable models of the same data correlate
#' almost perfectly whatever their differences, and a correlation measures association, not
#' agreement (a model that is 10\% higher everywhere still has \eqn{r = 1}; Bland and
#' Altman 1986).
#'
#' @param model1 A fitted \code{fastsae} model object, or a list containing two \code{fastsae} models.
#' @param model2 A second fitted \code{fastsae} model object. Ignored if \code{model1} is a list.
#' @param names Optional character vector of length 2 specifying descriptive labels
#'   for the models. Default is inferred from model types or call arguments.
#' @param thresholds Numeric vector of length 2 with the RSE (\%) cut-offs drawn in
#'   \code{autoplot(type = "rse")} (see \code{\link{add_reliability_flags}}). Default \code{c(25, 50)}.
#' @param top Share (below 1) or number (1 or more) of the highest and of the lowest domains
#'   whose overlap between the models is reported. Default \code{0.1} (10\%, at least one domain).
#' @param x An object of class \code{fastsae_comparison} (for \code{print} and \code{plot} methods).
#' @param object An object of class \code{fastsae_comparison} (for \code{summary} and \code{autoplot} methods).
#' @param y Ignored argument for compatibility with the generic \code{plot} method.
#' @param type Character string indicating comparison plot type: \code{"scatter"},
#'   \code{"comparison"}, \code{"difference"} (Bland-Altman plot: the difference against the
#'   mean of the two estimates, with the mean difference and the 95\% limits of agreement),
#'   \code{"mse"}, or \code{"rse"}. Default is \code{"scatter"}.
#' @param title Optional character string specifying a custom plot title.
#' @param ... Additional arguments.
#'
#' @return An S3 object of class \code{fastsae_comparison} containing:
#'   \itemize{
#'     \item \code{metrics}: A data frame of the metrics above (columns \code{Group}, \code{Metric},
#'       \code{Value}), formatted for printing.
#'     \item \code{agreement}: The same agreement, uncertainty and ranking metrics as a named
#'       numeric vector.
#'     \item \code{data}: A domain-level data frame with aligned estimates, MSE and RSE values,
#'       their differences, ranks (\code{rank_1}, \code{rank_2}, \code{rank_change}; 1 is the
#'       highest estimate) and \code{z_diff}.
#'     \item \code{names}: Vector of the two model names.
#'     \item \code{thresholds}: The RSE thresholds used.
#'   }
#'
#' @references
#' Bland, J. M. and Altman, D. G. (1986). Statistical methods for assessing agreement
#' between two methods of clinical measurement. \emph{The Lancet}, 327(8476), 307-310.
#'
#' Lin, L. I.-K. (1989). A concordance correlation coefficient to evaluate reproducibility.
#' \emph{Biometrics}, 45(1), 255-268.
#'
#' @seealso \code{\link{autoplot.fastsae_comparison}}, \code{\link{map_sae}},
#'   \code{\link{diagnose}} (model estimates against the direct estimates)
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
compare_sae <- function(model1, model2 = NULL, names = NULL, thresholds = c(25, 50), top = 0.1, ...) {
  .check_thresholds(thresholds)
  if (!is.numeric(top) || length(top) != 1L || is.na(top) || top <= 0) {
    cli::cli_abort("{.arg top} must be a positive number: a share below 1 or a number of domains.")
  }
  UseMethod("compare_sae")
}

#' @rdname compare_sae
#' @export
compare_sae.default <- function(model1, model2 = NULL, names = NULL, thresholds = c(25, 50), top = 0.1, ...) {
  # Support passing a list of two models
  if (is.list(model1) && !is.data.frame(model1) && !inherits(model1, "fastsae") && length(model1) >= 2) {
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

  # Merge on domain ID; for panel (domain x time) data also merge on time
  # to avoid a Cartesian product
  by_cols <- "domain"
  if ("subarea" %in% names(df1) && "subarea" %in% names(df2)) {
    by_cols <- c("domain", "subarea")
  } else if ("time" %in% names(df1) && "time" %in% names(df2)) {
    by_cols <- c("domain", "time")
  }
  merged <- merge(df1, df2, by = by_cols, suffixes = c("_1", "_2"), all = FALSE)

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

  # Ranks (1 = highest estimate) and differences in units of the RMSE of Model 1
  merged$rank_1 <- rank(-merged$estimate_1, na.last = "keep", ties.method = "average")
  merged$rank_2 <- rank(-merged$estimate_2, na.last = "keep", ties.method = "average")
  merged$rank_change <- merged$rank_2 - merged$rank_1
  merged$z_diff <- if (has_mse) merged$diff / sqrt(merged$mse_1) else NA_real_

  ag <- .agreement_stats(merged$estimate_1, merged$estimate_2)
  n_dom <- sum(stats::complete.cases(merged$estimate_1, merged$estimate_2))
  k_top <- max(1L, if (top < 1) as.integer(round(top * n_dom)) else as.integer(top))
  k_top <- min(k_top, n_dom)
  top_hi <- .top_overlap(merged$estimate_1, merged$estimate_2, k_top)
  top_lo <- .top_overlap(-merged$estimate_1, -merged$estimate_2, k_top)
  cor_s <- stats::cor(merged$estimate_1, merged$estimate_2, method = "spearman", use = "complete.obs")
  mean_rank_change <- mean(abs(merged$rank_change), na.rm = TRUE)

  med_z <- if (has_mse) stats::median(abs(merged$z_diff), na.rm = TRUE) else NA_real_
  pct_z1 <- if (has_mse) 100 * mean(abs(merged$z_diff) > 1, na.rm = TRUE) else NA_real_

  mean_ratio <- if (has_mse) mean(merged$mse_ratio, na.rm = TRUE) else NA_real_
  med_ratio <- if (has_mse) stats::median(merged$mse_ratio, na.rm = TRUE) else NA_real_
  n_eff2 <- if (has_mse) sum(merged$mse_2 < merged$mse_1, na.rm = TRUE) else NA_integer_
  pct_eff2 <- if (has_mse) (n_eff2 / nrow(merged)) * 100 else NA_real_

  mean_rse1 <- if (has_rse) mean(merged$rse_1, na.rm = TRUE) else NA_real_
  mean_rse2 <- if (has_rse) mean(merged$rse_2, na.rm = TRUE) else NA_real_

  # A posterior variance and a frequentist MSE are not the same quantity, so a
  # smaller one is not evidence of a more efficient estimator
  mixed <- xor(inherits(model1, "fastsae_hb"), inherits(model2, "fastsae_hb"))
  lower_label <- if (mixed) "reports lower uncertainty" else "has lower MSE"
  fmt <- function(v, f = "%.4f") if (is.na(v)) "N/A" else sprintf(f, v)

  metrics_df <- data.frame(
    Group = c(
      "Domains", rep("Agreement", 5), rep("Uncertainty", 2), rep("Ranking", 4),
      rep("Efficiency", 5)
    ),
    Metric = c(
      "Matched Domains",
      paste0("Mean Difference (", names[2], " - ", names[1], ")"),
      "95% Limits of Agreement",
      "Mean Absolute Difference (MAE)",
      "Root Mean Squared Difference (RMSD)",
      "Lin's Concordance Correlation (CCC)",
      paste0("Median |Difference| / RMSE of ", names[1]),
      paste0("% of domains with |Difference| > RMSE of ", names[1]),
      "Spearman Rank Correlation (rho)",
      "Mean Absolute Rank Change",
      sprintf("Shared Top %d (highest estimates)", k_top),
      sprintf("Shared Bottom %d (lowest estimates)", k_top),
      paste0("Mean MSE Ratio (", names[1], " / ", names[2], ")"),
      paste0("Median MSE Ratio (", names[1], " / ", names[2], ")"),
      if (mixed) {
        paste0("% of domains where ", names[2], " ", lower_label)
      } else {
        paste0("% of domains where ", names[2], " is more efficient")
      },
      paste0("Mean RSE % (", names[1], ")"),
      paste0("Mean RSE % (", names[2], ")")
    ),
    Value = c(
      as.character(nrow(merged)),
      fmt(ag$mean_diff),
      if (is.na(ag$loa[1])) "N/A" else sprintf("[%.4f, %.4f]", ag$loa[1], ag$loa[2]),
      fmt(ag$mae),
      fmt(ag$rmsd),
      fmt(ag$ccc),
      fmt(med_z, "%.3f"),
      fmt(pct_z1, "%.1f%%"),
      fmt(cor_s),
      fmt(mean_rank_change, "%.2f"),
      sprintf("%d / %d", top_hi, k_top),
      sprintf("%d / %d", top_lo, k_top),
      fmt(mean_ratio),
      fmt(med_ratio),
      fmt(pct_eff2, "%.1f%%"),
      fmt(mean_rse1, "%.2f%%"),
      fmt(mean_rse2, "%.2f%%")
    ),
    stringsAsFactors = FALSE
  )
  agreement <- c(
    mean_diff = ag$mean_diff, loa_lower = ag$loa[1], loa_upper = ag$loa[2],
    mae = ag$mae, rmsd = ag$rmsd, ccc = ag$ccc, spearman = cor_s,
    mean_abs_rank_change = mean_rank_change, top_shared = top_hi, bottom_shared = top_lo,
    top_k = k_top, median_abs_z = med_z, pct_abs_z_gt1 = pct_z1
  )

  structure(
    list(
      metrics = metrics_df,
      agreement = agreement,
      data = merged,
      names = names,
      thresholds = thresholds,
      mixed_inference = mixed,
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

  cat("\nAgreement, Ranking & Efficiency Metrics:\n")
  print(x$metrics, row.names = FALSE, right = FALSE)
  if (isTRUE(x$mixed_inference)) {
    cat("\n")
    cli::cli_alert_info(paste(
      "One model is Hierarchical Bayes: its posterior variance and the frequentist MSE of the",
      "other are not the same quantity, so a smaller value does not show a more accurate estimator."
    ))
  }

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
  # one label per row (domain, or domain/sub-area or domain:time for nested/panel data)
  df$label <- .row_labels(df)
  n1 <- object$names[1]
  n2 <- object$names[2]

  if (type == "scatter") {
    ag <- object$agreement
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
        subtitle = sprintf(
          "Lin's CCC = %.4f, Spearman rho = %.4f (dashed red line: identity)",
          ag[["ccc"]], ag[["spearman"]]
        ),
        x = paste(n1, "Estimate"),
        y = paste(n2, "Estimate")
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

    return(p)
  } else if (type == "difference") {
    # Bland-Altman plot: difference against the mean of the two estimates
    ag <- object$agreement
    df$mean_est <- (df$estimate_1 + df$estimate_2) / 2
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$mean_est, y = .data$diff)) +
      ggplot2::geom_hline(yintercept = 0, color = "gray60") +
      ggplot2::geom_hline(yintercept = ag[["mean_diff"]], color = "#2b83ba", linewidth = 0.8) +
      ggplot2::geom_hline(
        yintercept = c(ag[["loa_lower"]], ag[["loa_upper"]]),
        linetype = "dashed", color = "#d7191c", linewidth = 0.7
      ) +
      ggplot2::geom_point(size = 2.2, alpha = 0.8) +
      ggplot2::labs(
        title = title %||% paste("Agreement of", n2, "and", n1, "(Bland-Altman)"),
        subtitle = sprintf(
          "Mean difference %.4f (solid); 95%% limits of agreement [%.4f, %.4f] (dashed)",
          ag[["mean_diff"]], ag[["loa_lower"]], ag[["loa_upper"]]
        ),
        x = "Mean of the two estimates",
        y = paste(n2, "minus", n1)
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

    return(p)
  } else if (type == "mse") {
    long_mse <- rbind(
      data.frame(domain = df$label, mse = df$mse_1, Model = n1),
      data.frame(domain = df$label, mse = df$mse_2, Model = n2)
    )
    long_mse$domain <- factor(long_mse$domain, levels = df$label)

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
      data.frame(domain = df$label, rse = df$rse_1, Model = n1),
      data.frame(domain = df$label, rse = df$rse_2, Model = n2)
    )
    long_rse$domain <- factor(long_rse$domain, levels = df$label)

    p <- ggplot2::ggplot(long_rse, ggplot2::aes(x = .data$domain, y = .data$rse, color = .data$Model, group = .data$Model)) +
      ggplot2::geom_line(alpha = 0.7) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_hline(yintercept = object$thresholds[1], linetype = "dashed", color = "#2A9D8F") +
      ggplot2::geom_hline(yintercept = object$thresholds[2], linetype = "dashed", color = "#E76F51") +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::labs(
        title = title %||% "Relative Standard Error (RSE %) Comparison",
        subtitle = paste0("Reliability thresholds: RSE \u2264 ", object$thresholds[1], "% (reliable) and > ", object$thresholds[2], "% (unreliable)"),
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
      data.frame(domain = df$label, est = df$estimate_1, Model = n1),
      data.frame(domain = df$label, est = df$estimate_2, Model = n2)
    )
    long_est$domain <- factor(long_est$domain, levels = df$label)

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

#' Agreement of two sets of estimates: Bland-Altman mean difference and 95% limits of
#' agreement (difference est2 - est1), MAE, RMSD and Lin's concordance correlation
#' coefficient (moments with divisor n, as Lin 1989)
#' @noRd
.agreement_stats <- function(est1, est2) {
  ok <- stats::complete.cases(est1, est2)
  x <- est1[ok]
  y <- est2[ok]
  d <- y - x
  n <- length(d)
  sd_d <- if (n > 1) stats::sd(d) else NA_real_
  m1 <- mean(x)
  m2 <- mean(y)
  s12 <- mean((x - m1) * (y - m2))
  denom <- mean((x - m1)^2) + mean((y - m2)^2) + (m1 - m2)^2
  list(
    mean_diff = mean(d),
    loa = mean(d) + c(-1, 1) * 1.96 * sd_d,
    mae = mean(abs(d)),
    rmsd = sqrt(mean(d^2)),
    ccc = if (denom > 0) 2 * s12 / denom else NA_real_
  )
}

#' Number of the k domains with the highest est1 that are also among the k highest of est2
#' @noRd
.top_overlap <- function(est1, est2, k) {
  ok <- which(stats::complete.cases(est1, est2))
  top1 <- ok[order(est1[ok], decreasing = TRUE)][seq_len(k)]
  top2 <- ok[order(est2[ok], decreasing = TRUE)][seq_len(k)]
  length(intersect(top1, top2))
}

#' Extract standardized domain-level data frame
#' @noRd
.extract_domain_data <- function(model) {
  if (inherits(model, "fastsae")) {
    df <- .sae_table(model)
    if (is.null(df)) cli::cli_abort("Model does not contain estimation data frame.")
    dom_col <- intersect(c("domain", "area", "id"), names(df))[1]
    dom_vals <- df[[dom_col]] %||% seq_len(nrow(df))
    est_val <- df[[.sae_est_col(df)]]
    mse_val <- df$mse %||% (if (!is.null(df$sd)) df$sd^2 else NA_real_)
    rse_val <- df$rse %||% (if (!all(is.na(mse_val))) (sqrt(mse_val) / abs(est_val)) * 100 else NA_real_)

    out <- data.frame(
      domain = as.character(dom_vals),
      estimate = as.numeric(est_val),
      mse = as.numeric(mse_val),
      rse = as.numeric(rse_val),
      stringsAsFactors = FALSE
    )
    # Preserve time and subarea index to avoid Cartesian join
    time_col <- intersect(c("time", "period", "year"), names(df))
    if (length(time_col) > 0) out$time <- df[[time_col[1]]]
    subarea_col <- intersect(c("subarea", "sub_id", "subdomain"), names(df))
    if (length(subarea_col) > 0) out$subarea <- as.character(df[[subarea_col[1]]])
    return(out)
  } else if (is.data.frame(model)) {
    dom_col <- intersect(c("domain", "area", "id", "code"), names(model))
    dom_vals <- if (length(dom_col) > 0) model[[dom_col[1]]] else seq_len(nrow(model))
    est_col <- intersect(c("estimate", "hb", "eblup", "est", "y_hat"), names(model))
    if (length(est_col) == 0) cli::cli_abort("Data frame must contain an estimate column.")

    est_val <- model[[est_col[1]]]
    mse_val <- model$mse %||% rep(NA_real_, nrow(model))
    rse_val <- model$rse %||% (if (!all(is.na(mse_val))) (sqrt(mse_val) / abs(est_val)) * 100 else rep(NA_real_, nrow(model)))

    out <- data.frame(
      domain = as.character(dom_vals),
      estimate = as.numeric(est_val),
      mse = as.numeric(mse_val),
      rse = as.numeric(rse_val),
      stringsAsFactors = FALSE
    )
    time_col <- intersect(c("time", "period", "year"), names(model))
    if (length(time_col) > 0) out$time <- model[[time_col[1]]]
    subarea_col <- intersect(c("subarea", "sub_id", "subdomain"), names(model))
    if (length(subarea_col) > 0) out$subarea <- as.character(model[[subarea_col[1]]])
    return(out)
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
        "TWOFOLD" = "Two-fold FH",
        model$model
      )
      return(m)
    }
  }
  return(fallback)
}
