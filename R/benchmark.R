# ============================================================================
# Benchmarking Methods for Small Area Estimation (fastsae)
# ============================================================================

#' Benchmark Small Area Estimation Predictions to Aggregate Targets
#'
#' @description
#' Calibrates small area model predictions (from \code{ebp_area}, \code{eblup_fh},
#' \code{eblup_sfh}, \code{eblup_stfh}, or \code{eblup_bhf}) so that their weighted
#' aggregate matches a known or direct benchmark at an overarching level (e.g. provincial
#' or national totals/means), as required in official statistics production.
#'
#' Supports Ratio (multiplicative), Difference (additive), Optimal (quadratic loss / MSE-weighted),
#' and Logit (bounded rate) benchmarking methods based on established SAE literature
#' (Rao and Molina, 2015; Datta et al., 2011; Steorts et al., 2014; Berg and Fuller, 2014).
#'
#' @param object A fitted \code{fastsae} model object (e.g., from \code{ebp_area},
#'   \code{eblup_fh}, etc.) or a numeric vector of model predictions.
#' @param target Numeric value or named vector specifying the aggregate benchmark
#'   target \eqn{T_g}. If \code{group} is specified, \code{target} can be a single
#'   numeric (if identical target applies to all groups), a named numeric vector
#'   matching group levels, or a data frame with columns \code{group} and \code{target}.
#' @param weight Optional numeric vector or character string naming the domain
#'   benchmark weight column in \code{object$data} (e.g., population sizes \eqn{N_d}
#'   or population shares \eqn{N_d / \sum N_j}). If \code{NULL}, attempts to detect
#'   population/weight columns automatically, or defaults to equal weights.
#' @param method Character string specifying the benchmarking calibration method:
#'   \itemize{
#'     \item \code{"ratio"}: Proportional / multiplicative adjustment (\eqn{\hat{\theta}_d^{\text{BM}} = \hat{\theta}_d \times (T_g / \sum_j w_j \hat{\theta}_j)}).
#'       Preserves non-negativity and area proportions (You and Rao, 2002; Rao and Molina, 2015).
#'     \item \code{"difference"}: Uniform additive adjustment (\eqn{\hat{\theta}_d^{\text{BM}} = \hat{\theta}_d + (T_g - \sum_j w_j \hat{\theta}_j) / \sum_j w_j}).
#'       Preserves absolute differences between domain estimates.
#'     \item \code{"optimal"}: Variance-weighted quadratic loss benchmarking (Datta et al., 2011;
#'       Steorts et al., 2014; Rao and Molina, 2015, Section 10.3). Adjusts domains proportional
#'       to their uncertainty (\eqn{\text{MSE}_d / w_d}), so domains with higher estimation error
#'       absorb larger adjustments while highly precise domains remain stable.
#'     \item \code{"logit"}: Logit-scale additive shift for bounded indicators \eqn{\theta_d \in (0, 1)}
#'       (Berg and Fuller, 2014). Guarantees that benchmarked rates strictly remain in the unit interval
#'       via 1D root-finding.
#'   }
#' @param group Optional vector or character string naming the grouping/stratum column
#'   (e.g., province or region) for multi-level hierarchical benchmarking.
#' @param type Character string: \code{"mean"} (default, benchmark target represents
#'   weighted average/rate, \eqn{\sum_{d \in g} w_d \hat{\theta}_d = T_g} with normalized weights)
#'   or \code{"total"} (benchmark target represents population total, \eqn{\sum_{d \in g} w_d \hat{\theta}_d = T_g}).
#' @param ... Additional arguments passed to methods.
#'
#' @return An object of class \code{c("fastsae_benchmark", "data.frame")} containing:
#' \itemize{
#'   \item \code{domain}: Domain identifier.
#'   \item \code{group}: Stratum / overarching group identifier (if specified).
#'   \item \code{weight}: Benchmark weight \eqn{w_d}.
#'   \item \code{original}: Model-based prediction before benchmarking (\eqn{\hat{\theta}_d}).
#'   \item \code{benchmarked}: Calibrated estimate satisfying the benchmark constraint (\eqn{\hat{\theta}_d^{\text{BM}}}).
#'   \item \code{adjustment}: Absolute adjustment (\eqn{\hat{\theta}_d^{\text{BM}} - \hat{\theta}_d}).
#'   \item \code{rel_adjustment}: Relative adjustment percentage (\%).
#'   \item \code{target}: Corresponding benchmark target \eqn{T_g}.
#' }
#'
#' @references
#' \enumerate{
#'   \item Rao, J. N. K., and Molina, I. (2015). \emph{Small Area Estimation} (2nd ed.).
#'     John Wiley & Sons. Chapter 10: "Benchmarking and Other Issues", pp. 297-315.
#'   \item You, Y., and Rao, J. N. K. (2002). A pseudo-empirical best linear unbiased
#'     prediction approach to small area estimation using survey weights.
#'     \emph{The Canadian Journal of Statistics}, 30(3), 431-439.
#'   \item Datta, G. S., Ghosh, M., Steorts, R., and Maples, J. (2011). Bayesian
#'     benchmarking with applications to small area estimation. \emph{Test}, 20(3), 574-588.
#'   \item Steorts, R. C., Hall, P., and Ghosh, M. (2014). General benchmarking under
#'     quadratic loss with applications to small area estimation.
#'     \emph{Journal of Survey Statistics and Methodology}, 2(2), 173-193.
#'   \item Berg, E., and Fuller, W. A. (2014). Small area prediction of proportions
#'     with a constrained multinomial logit model. \emph{Journal of Survey Statistics and Methodology}, 2(3), 256-283.
#' }
#'
#' @export
#' @examples
#' library(fastsae)
#' data(mys)
#'
#' # 1. Fit Fay-Herriot model
#' fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)
#'
#' # 2. Ratio benchmarking to state/national average target (e.g. target = 6.5)
#' bm_ratio <- benchmark(fit_fh, target = 6.5, method = "ratio")
#' head(bm_ratio)
#'
#' # 3. Optimal (MSE-weighted) benchmarking
#' bm_opt <- benchmark(fit_fh, target = 6.5, method = "optimal")
#' head(bm_opt)
#'
#' # 4. Logit benchmarking for bounded rates (proportions in [0, 1])
#' mys$prop <- mys$y / 100
#' mys$var_prop <- mys$vardir / 10000
#' fit_beta <- ebp_area(prop ~ x1, vardir = "var_prop", data = mys, family = "beta")
#' bm_logit <- benchmark(fit_beta, target = 0.065, method = "logit")
#' head(bm_logit)
benchmark <- function(object, ...) {
  UseMethod("benchmark")
}

#' @rdname benchmark
#' @export
benchmark_sae <- function(object, ...) {
  benchmark(object, ...)
}

#' @rdname benchmark
#' @export
benchmark.fastsae <- function(
  object,
  target,
  weight = NULL,
  method = c("ratio", "difference", "optimal", "logit"),
  group = NULL,
  type = c("mean", "total"),
  ...
) {
  method <- match.arg(method)
  type <- match.arg(type)

  # 1. Extract estimates and domain identifiers
  df_est <- object$df_ebp %||% object$df_eblup
  if (is.null(df_est)) {
    cli::cli_abort("The {.cls fastsae} object does not contain fitted area estimates ({.code df_ebp} or {.code df_eblup}).")
  }

  domain_vec <- df_est$domain %||% df_est$area %||% seq_len(nrow(df_est))
  y_hat <- df_est$ebp %||% df_est$eblup %||% df_est$est
  if (is.null(y_hat)) {
    cli::cli_abort("Could not find estimation column ({.code ebp}, {.code eblup}, or {.code est}) in model predictions.")
  }

  mse_vec <- df_est$mse %||% (if (!is.null(df_est$sd)) df_est$sd^2 else NULL)

  # 2. Extract weights
  w_vec <- NULL
  if (is.character(weight) && length(weight) == 1) {
    if (!is.null(object$data) && weight %in% names(object$data)) {
      w_vec <- object$data[[weight]]
    } else if (weight %in% names(df_est)) {
      w_vec <- df_est[[weight]]
    } else {
      cli::cli_abort("Weight column {.val {weight}} not found in object data.")
    }
  } else if (is.numeric(weight)) {
    if (length(weight) != length(y_hat)) {
      cli::cli_abort("Length of {.arg weight} ({length(weight)}) does not match number of domains ({length(y_hat)}).")
    }
    w_vec <- as.numeric(weight)
  } else if (is.null(weight)) {
    # Default to equal weights across domains
    w_vec <- rep(1, length(y_hat))
  }

  # 3. Extract grouping
  group_vec <- NULL
  if (is.character(group) && length(group) == 1) {
    if (!is.null(object$data) && group %in% names(object$data)) {
      group_vec <- as.character(object$data[[group]])
    } else if (group %in% names(df_est)) {
      group_vec <- as.character(df_est[[group]])
    } else {
      cli::cli_abort("Grouping column {.val {group}} not found in object data.")
    }
  } else if (!is.null(group)) {
    if (length(group) != length(y_hat)) {
      cli::cli_abort("Length of {.arg group} ({length(group)}) does not match number of domains ({length(y_hat)}).")
    }
    group_vec <- as.character(group)
  }

  # Call internal benchmarking worker
  .benchmark_worker(
    domain = domain_vec,
    y_hat = y_hat,
    target = target,
    weight = w_vec,
    method = method,
    group = group_vec,
    type = type,
    mse = mse_vec,
    call = match.call(),
    model_family = object$family %||% "gaussian"
  )
}

#' @rdname benchmark
#' @export
benchmark.default <- function(
  object,
  target,
  weight = NULL,
  method = c("ratio", "difference", "optimal", "logit"),
  group = NULL,
  type = c("mean", "total"),
  ...
) {
  method <- match.arg(method)
  type <- match.arg(type)

  if (!is.numeric(object)) {
    cli::cli_abort("{.arg object} must be a numeric vector of predictions or a {.cls fastsae} object.")
  }

  y_hat <- as.numeric(object)
  n <- length(y_hat)
  domain_vec <- names(object) %||% paste0("domain_", seq_len(n))

  w_vec <- if (is.null(weight)) rep(1, n) else as.numeric(weight)
  if (length(w_vec) != n) {
    cli::cli_abort("Length of {.arg weight} ({length(w_vec)}) does not match length of {.arg object} ({n}).")
  }

  group_vec <- if (!is.null(group)) as.character(group) else NULL
  if (!is.null(group_vec) && length(group_vec) != n) {
    cli::cli_abort("Length of {.arg group} ({length(group_vec)}) does not match length of {.arg object} ({n}).")
  }

  .benchmark_worker(
    domain = domain_vec,
    y_hat = y_hat,
    target = target,
    weight = w_vec,
    method = method,
    group = group_vec,
    type = type,
    mse = NULL,
    call = match.call(),
    model_family = "generic"
  )
}

#' Internal worker to compute small area benchmarking calibrations
#' @noRd
.benchmark_worker <- function(
  domain,
  y_hat,
  target,
  weight,
  method,
  group = NULL,
  type = "mean",
  mse = NULL,
  call = NULL,
  model_family = "gaussian"
) {
  n <- length(y_hat)

  if (any(weight <= 0, na.rm = TRUE)) {
    cli::cli_abort("Benchmark {.arg weight} must be strictly positive for all domains.")
  }

  # Build data frame
  df_work <- data.frame(
    domain = domain,
    y_hat = y_hat,
    weight = weight,
    mse = if (!is.null(mse) && length(mse) == n) mse else rep(NA_real_, n),
    stringsAsFactors = FALSE
  )
  if (!is.null(group)) {
    df_work$group <- as.character(group)
  } else {
    df_work$group <- "All"
  }

  unique_groups <- unique(df_work$group)
  n_groups <- length(unique_groups)

  # Validate and parse target per group
  target_map <- numeric(n_groups)
  names(target_map) <- unique_groups

  if (is.data.frame(target)) {
    if (!all(c("group", "target") %in% names(target))) {
      cli::cli_abort("Target data frame must contain columns {.val group} and {.val target}.")
    }
    for (g in unique_groups) {
      val <- target$target[target$group == g]
      if (length(val) == 0) {
        cli::cli_abort("Target value for group {.val {g}} not found in target data frame.")
      }
      target_map[g] <- as.numeric(val[1])
    }
  } else if (is.numeric(target)) {
    if (length(target) == 1 && n_groups == 1) {
      target_map[1] <- as.numeric(target)
    } else if (length(target) == 1 && n_groups > 1) {
      # Same target for each group (e.g. constant mean rate)
      target_map[] <- as.numeric(target)
    } else if (!is.null(names(target))) {
      for (g in unique_groups) {
        if (!g %in% names(target)) {
          cli::cli_abort("Named target vector does not contain group {.val {g}}.")
        }
        target_map[g] <- as.numeric(target[g])
      }
    } else if (length(target) == n_groups) {
      target_map[] <- as.numeric(target)
    } else {
      cli::cli_abort("Length of {.arg target} ({length(target)}) does not match number of groups ({n_groups}).")
    }
  } else {
    cli::cli_abort("{.arg target} must be a numeric value, named numeric vector, or data frame.")
  }

  df_work$target <- target_map[df_work$group]
  df_work$y_bm <- NA_real_

  # Compute benchmarking calibration per group
  for (g in unique_groups) {
    idx_g <- which(df_work$group == g)
    y_g <- df_work$y_hat[idx_g]
    w_g <- df_work$weight[idx_g]
    T_g <- target_map[g]

    # Normalize weights if type == "mean"
    w_calc <- if (type == "mean") w_g / sum(w_g) else w_g
    agg_initial <- sum(w_calc * y_g)

    if (method == "ratio") {
      # Ratio / Multiplicative Benchmarking (You and Rao, 2002; Rao and Molina, 2015)
      if (abs(agg_initial) < .Machine$double.eps) {
        cli::cli_abort("Initial weighted aggregation for group {.val {g}} is virtually zero; ratio benchmarking is undefined.")
      }
      ratio_adj <- T_g / agg_initial
      df_work$y_bm[idx_g] <- y_g * ratio_adj

    } else if (method == "difference") {
      # Difference / Additive Benchmarking (Rao and Molina, 2015, Sec 10.3)
      diff_val <- T_g - agg_initial
      # Uniform addition across domains: diff_val / sum(w_calc)
      df_work$y_bm[idx_g] <- y_g + (diff_val / sum(w_calc))

    } else if (method == "optimal") {
      # Optimal Quadratic Loss / Variance-Weighted Benchmarking (Datta et al., 2011; Steorts et al., 2014)
      # Minimizes sum_d (y_bm_d - y_d)^2 / (MSE_d / w_d^2) subject to sum w_d y_bm_d = T_g
      mse_g <- df_work$mse[idx_g]
      diff_val <- T_g - agg_initial

      if (any(is.na(mse_g)) || all(mse_g <= 0)) {
        # Fallback to difference benchmarking if MSE is not available
        cli::cli_warn("MSE not available or non-positive for group {.val {g}}; falling back to difference benchmarking.")
        df_work$y_bm[idx_g] <- y_g + (diff_val / sum(w_calc))
      } else {
        # Optimal formula: y_bm = y + (MSE / w) * [diff_val / sum(MSE)]
        # When w_calc is used: c_d = mse_g / (w_calc^2), so adjustment is c_d * w_calc * (diff_val / sum(c_j * w_calc^2))
        inv_prec <- pmax(mse_g, 1e-8)
        denom <- sum(inv_prec)
        lambda <- diff_val / denom
        df_work$y_bm[idx_g] <- y_g + (inv_prec / w_calc) * lambda
      }

    } else if (method == "logit") {
      # Logit-scale additive shift for bounded indicators in (0, 1) (Berg and Fuller, 2014)
      if (any(y_g <= 0 | y_g >= 1)) {
        cli::cli_abort(c(
          "Logit benchmarking requires model estimates to be strictly within (0, 1).",
          "x" = "Got estimates outside (0, 1) in group {.val {g}}."
        ))
      }
      if (type == "mean" && (T_g <= 0 || T_g >= 1)) {
        cli::cli_abort("For logit benchmarking with type = 'mean', target must be strictly in (0, 1).")
      }

      logit_y <- stats::qlogis(y_g)

      f_obj <- function(alpha) {
        sum(w_calc * stats::plogis(logit_y + alpha)) - T_g
      }

      # Solve for alpha using uniroot
      # Find bracket:
      low_b <- -20
      upp_b <- 20
      f_low <- f_obj(low_b)
      f_upp <- f_obj(upp_b)

      if (f_low * f_upp > 0) {
        # Expand bracket
        low_b <- -50
        upp_b <- 50
      }

      sol <- stats::uniroot(f_obj, interval = c(low_b, upp_b), tol = 1e-9)
      alpha_opt <- sol$root
      df_work$y_bm[idx_g] <- stats::plogis(logit_y + alpha_opt)
    }
  }

  # Build return table
  res_df <- data.frame(
    domain = df_work$domain,
    stringsAsFactors = FALSE
  )
  if (!is.null(group)) {
    res_df$group <- df_work$group
  }
  res_df$weight <- df_work$weight
  res_df$original <- df_work$y_hat
  res_df$benchmarked <- df_work$y_bm
  res_df$adjustment <- df_work$y_bm - df_work$y_hat
  res_df$rel_adjustment_pct <- ifelse(abs(df_work$y_hat) > .Machine$double.eps,
                                      (res_df$adjustment / df_work$y_hat) * 100, NA_real_)
  res_df$target <- df_work$target

  # Calculate verification aggregation
  verification <- list()
  for (g in unique_groups) {
    idx_g <- which(df_work$group == g)
    w_calc <- if (type == "mean") df_work$weight[idx_g] / sum(df_work$weight[idx_g]) else df_work$weight[idx_g]
    orig_sum <- sum(w_calc * df_work$y_hat[idx_g])
    bm_sum <- sum(w_calc * df_work$y_bm[idx_g])
    t_val <- target_map[g]
    verification[[g]] <- list(
      target = t_val,
      original_sum = orig_sum,
      benchmarked_sum = bm_sum,
      discrepancy = abs(bm_sum - t_val)
    )
  }

  attr(res_df, "method") <- method
  attr(res_df, "type") <- type
  attr(res_df, "verification") <- verification
  attr(res_df, "call") <- call

  class(res_df) <- c("fastsae_benchmark", "data.frame")
  return(res_df)
}

#' Print a fastsae_benchmark object
#'
#' @param x An object of class \code{fastsae_benchmark}.
#' @param ... Additional arguments.
#'
#' @return The original object invisibly.
#' @export
print.fastsae_benchmark <- function(x, ...) {
  cat("=== fastsae Small Area Benchmark Calibration ===\n")

  method_name <- switch(attr(x, "method") %||% "ratio",
    "ratio"      = "Ratio (Multiplicative Adjustment)",
    "difference" = "Difference (Additive Adjustment)",
    "optimal"    = "Optimal (Variance-Weighted Quadratic Loss)",
    "logit"      = "Logit (Bounded Unit Interval Shift)",
    attr(x, "method")
  )
  type_str <- if (identical(attr(x, "type"), "total")) "Population Total" else "Weighted Mean / Rate"

  cat("Method:", method_name, "\n")
  cat("Target Type:", type_str, "\n")
  cat("Total Domains:", nrow(x), "\n")

  # Print verification
  verif <- attr(x, "verification")
  if (!is.null(verif) && length(verif) > 0) {
    cat("\n-- Aggregate Consistency Check --\n")
    for (g in names(verif)) {
      v <- verif[[g]]
      g_label <- if (g == "All") "Target" else paste0("Group '", g, "'")
      disc_ok <- v$discrepancy < 1e-5
      status_icon <- if (disc_ok) "[CONSISTENT]" else "[APPROX]"
      cat(sprintf("%s: Target = %.5f | Original = %.5f -> Calibrated = %.5f %s\n",
                  g_label, v$target, v$original_sum, v$benchmarked_sum, status_icon))
    }
  }

  cat("\nFirst 6 benchmarked domains:\n")
  print(utils::head(as.data.frame(x), 6), ...)
  if (nrow(x) > 6) {
    cat("... and", nrow(x) - 6, "more domains.\n")
  }

  cat("\n")
  invisible(x)
}

#' Summary of a fastsae_benchmark object
#'
#' @param object An object of class \code{fastsae_benchmark}.
#' @param ... Additional arguments.
#'
#' @return A summary table invisibly.
#' @export
summary.fastsae_benchmark <- function(object, ...) {
  print(object, ...)
  cat("-- Adjustment Statistics --\n")
  cat("Absolute Adjustment Summary:\n")
  print(summary(object$adjustment))
  cat("\nRelative Adjustment (%) Summary:\n")
  print(summary(object$rel_adjustment_pct))
  invisible(object)
}

#' Plot Method for fastsae_benchmark Objects
#'
#' Visualizes original model predictions versus benchmarked calibrated estimates.
#'
#' @param x An object of class \code{fastsae_benchmark}.
#' @param ... Additional arguments passed to plotting methods.
#'
#' @return A \code{ggplot2} plot object.
#' @export
autoplot.fastsae_benchmark <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggplot2} is required to plot benchmarked objects.")
  }

  df_plot <- as.data.frame(x)
  has_group <- "group" %in% names(df_plot) && length(unique(df_plot$group)) > 1

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = .data$original, y = .data$benchmarked)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50")

  if (has_group) {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = .data$group), alpha = 0.85, size = 2) +
      ggplot2::labs(color = "Group")
  } else {
    p <- p + ggplot2::geom_point(color = "#0d6efd", alpha = 0.85, size = 2)
  }

  method_name <- switch(attr(x, "method") %||% "ratio",
    "ratio"      = "Ratio (Multiplicative)",
    "difference" = "Difference (Additive)",
    "optimal"    = "Optimal (Variance-Weighted)",
    "logit"      = "Logit (Bounded (0, 1))",
    attr(x, "method")
  )

  p <- p +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::labs(
      title = paste0("Small Area Benchmark Calibration (", method_name, ")"),
      subtitle = "Original Model Prediction vs Benchmarked Target-Calibrated Estimate",
      x = "Original Model Prediction",
      y = "Benchmarked Estimate"
    )

  return(p)
}

#' @rdname autoplot.fastsae_benchmark
#' @export
plot.fastsae_benchmark <- function(x, ...) {
  p <- autoplot(x, ...)
  print(p)
  invisible(p)
}
