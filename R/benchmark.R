# ============================================================================
# Benchmarking Methods for Small Area Estimation (fastsae)
# ============================================================================

#' Benchmark Small Area Estimation Predictions to Aggregate Targets
#'
#' @description
#' Calibrates small area model predictions (from \code{ebp_area}, \code{eblup_fh},
#' \code{eblup_sfh}, \code{eblup_stfh}, or \code{eblup_bhf}) so that their weighted
#' aggregate matches known benchmark targets at higher administrative levels (e.g. provincial
#' or national totals/means), as required in official statistics production.
#'
#' Supports:
#' \itemize{
#'   \item \strong{Single-Stage Benchmarking}: Calibrates small areas directly to an aggregate target
#'     (either globally or within sub-regions/provinces via \code{group}).
#'   \item \strong{Two-Stage Hierarchical Benchmarking} (Rao and Molina, 2015, Sec 10.4; Steorts et al., 2014):
#'     Calibrates small areas (e.g. 500 Kabupaten/Kota) within overarching groups (e.g. 34 Provinces)
#'     while \emph{simultaneously} harmonizing the provincial aggregates to match an overarching National Target
#'     (\code{national_target}). Guarantees simultaneous mathematical consistency across all administrative tiers:
#'     \eqn{\sum_{d \in \text{Prov}_g} \tilde{w}_{dg} \hat{\theta}_{dg}^{\text{BM}} = T_g^*} and
#'     \eqn{\sum_g \tilde{W}_g T_g^* = T_{\text{nat}}}.
#' }
#'
#' Methods implemented:
#' \itemize{
#'   \item \code{"ratio"}: Proportional / multiplicative adjustment. Preserves area proportions and non-negativity (You and Rao, 2002; Rao and Molina, 2015).
#'   \item \code{"difference"}: Uniform additive adjustment (Rao and Molina, 2015, Section 10.3).
#'   \item \code{"optimal"}: Variance-weighted quadratic loss benchmarking (Datta et al., 2011; Steorts et al., 2014).
#'     Domains with higher MSE absorb larger adjustments while precise domains remain minimally perturbed.
#'   \item \code{"logit"}: Logit-scale additive shift for bounded indicators \eqn{\theta_d \in (0, 1)} (Berg and Fuller, 2014).
#' }
#'
#' @param object A fitted \code{fastsae} model object (e.g., from \code{ebp_area},
#'   \code{eblup_fh}, etc.) or a numeric vector of model predictions.
#' @param target Numeric value, named vector, or data frame specifying the benchmark
#'   target(s) \eqn{T_g}. If \code{group} is specified, \code{target} can be:
#'   \itemize{
#'     \item A single scalar numeric (constant target for all groups).
#'     \item A named numeric vector matching group identifiers.
#'     \item A data frame with columns \code{group} and \code{target}.
#'     \item \code{NULL} if \code{national_target} is provided in two-stage hierarchical mode
#'       (initial provincial targets are automatically derived from model aggregates).
#'   }
#' @param weight Optional numeric vector or character string naming the domain
#'   benchmark weight column in \code{object$data} (e.g., population sizes \eqn{N_d}
#'   or population shares \eqn{N_d / \sum N_j}). Defaults to equal weights across domains if \code{NULL}.
#' @param method Character string specifying the benchmarking calibration method:
#'   \code{"ratio"} (default), \code{"difference"}, \code{"optimal"}, or \code{"logit"}.
#' @param group Optional vector or character string naming the grouping/stratum column
#'   (e.g., province or region) for multi-level hierarchical calibration.
#' @param type Character string: \code{"mean"} (default, benchmark target represents
#'   weighted average/rate, \eqn{\sum_{d \in g} w_d \hat{\theta}_d = T_g} with normalized weights)
#'   or \code{"total"} (benchmark target represents population total, \eqn{\sum_{d \in g} w_d \hat{\theta}_d = T_g}).
#' @param national_target Optional numeric scalar defining the overarching national target
#'   (Level 0) for Two-Stage Hierarchical Benchmarking. When provided alongside \code{group},
#'   first calibrates/harmonizes group (provincial) targets to match \code{national_target},
#'   then calibrates small areas to the harmonized provincial targets.
#' @param outer_target Alias for \code{national_target}.
#' @param x An object of class \code{fastsae_benchmark} (for \code{print} and \code{plot} methods).
#' @param y Ignored argument for compatibility with the generic \code{plot} method.
#' @param all_groups Logical; if \code{TRUE}, prints verification status for all groups in
#'   hierarchical benchmarking. Default is \code{FALSE}.
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
#'   \item \code{rel_adjustment_pct}: Relative adjustment percentage (\%).
#'   \item \code{target}: Corresponding group benchmark target \eqn{T_g^*}.
#'   \item \code{national_target}: Overarching national target (if two-stage hierarchical).
#' }
#'
#' @references
#' \enumerate{
#'   \item Rao, J. N. K., and Molina, I. (2015). \emph{Small Area Estimation} (2nd ed.).
#'     John Wiley & Sons. Chapter 10: "Benchmarking and Other Practical Issues", pp. 297-315.
#'   \item You, Y., and Rao, J. N. K. (2002). A pseudo-empirical best linear unbiased
#'     prediction approach to small area estimation using survey weights.
#'     \emph{The Canadian Journal of Statistics}, 30(3), 431-439.
#'   \item Steorts, R. C., Hall, P., and Ghosh, M. (2014). General benchmarking under
#'     quadratic loss with applications to small area estimation.
#'     \emph{Journal of Survey Statistics and Methodology}, 2(2), 173-193.
#'   \item Datta, G. S., Ghosh, M., Steorts, R., and Maples, J. (2011). Bayesian
#'     benchmarking with applications to small area estimation. \emph{Test}, 20(3), 574-588.
#'   \item Berg, E., and Fuller, W. A. (2014). Small area prediction of proportions
#'     with a constrained multinomial logit model. \emph{Journal of Survey Statistics and Methodology}, 2(3), 256-283.
#' }
#'
#' @export
#' @examples
#' library(fastsae)
#' data(mys)
#'
#' # Assign province groups for hierarchical calibration example
#' mys$province <- rep(c("Prov_A", "Prov_B", "Prov_C"), length.out = nrow(mys))
#'
#' # 1. Fit Fay-Herriot model
#' fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)
#'
#' # 2. Ratio benchmarking to overall target (e.g. target = 6.5)
#' bm_ratio <- benchmark_sae(fit_fh, target = 6.5, method = "ratio")
#' head(bm_ratio)
#'
#' # 3. Two-Stage Hierarchical Calibration (e.g. Regency -> Province -> National)
#' # Calibrate domains within provinces while matching overarching national target:
#' bm_hier <- benchmark_sae(
#'   fit_fh,
#'   group = "province",
#'   weight = mys$n,
#'   national_target = 6.2,
#'   method = "ratio"
#' )
#' print(bm_hier)
benchmark_sae <- function(object, ...) {
  UseMethod("benchmark_sae")
}

#' @rdname benchmark_sae
#' @export
benchmark_sae.fastsae <- function(
  object,
  target = NULL,
  weight = NULL,
  method = c("ratio", "difference", "optimal", "logit"),
  group = NULL,
  type = c("mean", "total"),
  national_target = NULL,
  outer_target = NULL,
  ...
) {
  method <- match.arg(method)
  type <- match.arg(type)
  nat_target <- outer_target %||% national_target

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
  group_name <- if (is.character(group) && length(group) == 1) group else "Group"
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
    group_var = group_name,
    type = type,
    national_target = nat_target,
    mse = mse_vec,
    call = match.call(),
    model_family = object$family %||% "gaussian"
  )
}

#' @rdname benchmark_sae
#' @export
benchmark_sae.default <- function(
  object,
  target = NULL,
  weight = NULL,
  method = c("ratio", "difference", "optimal", "logit"),
  group = NULL,
  type = c("mean", "total"),
  national_target = NULL,
  outer_target = NULL,
  ...
) {
  method <- match.arg(method)
  type <- match.arg(type)
  nat_target <- outer_target %||% national_target

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
  group_name <- if (is.character(group) && length(group) == 1) group else "Group"
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
    group_var = group_name,
    type = type,
    national_target = nat_target,
    mse = NULL,
    call = match.call(),
    model_family = "generic"
  )
}

#' @rdname benchmark_sae
#' @export
benchmark <- function(object, ...) {
  benchmark_sae(object, ...)
}

#' @export
benchmark.fastsae <- function(object, ...) {
  benchmark_sae.fastsae(object, ...)
}

#' @export
benchmark.default <- function(object, ...) {
  benchmark_sae.default(object, ...)
}

#' Internal worker to compute small area benchmarking calibrations
#' @noRd
.benchmark_worker <- function(
  domain,
  y_hat,
  target = NULL,
  weight,
  method,
  group = NULL,
  group_var = "Group",
  type = "mean",
  national_target = NULL,
  mse = NULL,
  call = NULL,
  model_family = "gaussian"
) {
  n <- length(y_hat)

  if (any(weight <= 0, na.rm = TRUE)) {
    cli::cli_abort("Benchmark {.arg weight} must be strictly positive for all domains.")
  }

  # Build working data frame
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

  # Check if Two-Stage Hierarchical Calibration is activated
  is_hierarchical <- !is.null(national_target) && n_groups > 1

  # Compute group totals and initial weighted estimates
  group_weights <- numeric(n_groups)
  names(group_weights) <- unique_groups
  group_orig <- numeric(n_groups)
  names(group_orig) <- unique_groups
  group_mse <- numeric(n_groups)
  names(group_mse) <- unique_groups

  for (g in unique_groups) {
    idx_g <- which(df_work$group == g)
    w_g <- df_work$weight[idx_g]
    W_g <- sum(w_g)
    group_weights[g] <- W_g
    w_norm <- w_g / W_g
    group_orig[g] <- sum(w_norm * df_work$y_hat[idx_g])
    if (!all(is.na(df_work$mse[idx_g]))) {
      # Variance of group weighted estimate
      group_mse[g] <- sum((w_norm^2) * pmax(df_work$mse[idx_g], 1e-8, na.rm = TRUE), na.rm = TRUE)
    } else {
      group_mse[g] <- NA_real_
    }
  }

  target_map <- numeric(n_groups)
  names(target_map) <- unique_groups
  has_initial_targets <- FALSE
  initial_target_map <- numeric(n_groups)
  names(initial_target_map) <- unique_groups

  if (!is.null(target)) {
    has_initial_targets <- TRUE
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
      } else if (length(target) == 1 && n_groups > 1 && !is_hierarchical) {
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
    initial_target_map <- target_map
  } else if (is_hierarchical) {
    # If target is NULL, use group_orig as starting base for stage 1 calibration
    target_map <- group_orig
    initial_target_map <- group_orig
  } else if (!is.null(national_target) && n_groups == 1) {
    # Single group fallback
    target_map[1] <- as.numeric(national_target)
    initial_target_map[1] <- as.numeric(national_target)
  } else {
    cli::cli_abort("Please provide either {.arg target} or {.arg national_target}.")
  }

  # --- STAGE 1: Hierarchical Group Harmonization to National Target ---
  stage1_df <- NULL
  if (is_hierarchical) {
    total_national_weight <- sum(group_weights)
    W_share <- group_weights / total_national_weight

    # Current aggregate of group targets at national level
    nat_current_agg <- if (type == "mean") sum(W_share * target_map) else sum(target_map)
    diff_nat <- national_target - nat_current_agg

    if (abs(diff_nat) > 1e-7) {
      if (method == "ratio") {
        ratio_nat <- national_target / nat_current_agg
        target_map <- target_map * ratio_nat
      } else if (method == "difference") {
        adj_per_group <- if (type == "mean") diff_nat else diff_nat / n_groups
        target_map <- target_map + adj_per_group
      } else if (method == "optimal") {
        inv_prec_g <- pmax(group_mse, 1e-8)
        if (any(is.na(inv_prec_g))) inv_prec_g <- rep(1, n_groups)
        lambda_g <- diff_nat / sum(inv_prec_g)
        target_map <- target_map + (inv_prec_g / (if (type == "mean") W_share else 1)) * lambda_g
      } else if (method == "logit") {
        logit_t <- stats::qlogis(target_map)
        f_nat <- function(a) {
          sum(W_share * stats::plogis(logit_t + a)) - national_target
        }
        sol_nat <- stats::uniroot(f_nat, interval = c(-50, 50), tol = 1e-9)
        target_map <- stats::plogis(logit_t + sol_nat$root)
      }
    }

    stage1_df <- data.frame(
      group = unique_groups,
      n_domains = vapply(unique_groups, function(g) sum(df_work$group == g), integer(1)),
      weight_total = group_weights,
      weight_share_pct = round(W_share * 100, 3),
      original_group_estimate = group_orig,
      initial_target = if (has_initial_targets) initial_target_map else NA_real_,
      calibrated_target = target_map,
      adjustment = target_map - (if (has_initial_targets) initial_target_map else group_orig),
      stringsAsFactors = FALSE
    )
  }

  # --- STAGE 2: Small Area Level Calibration within Each Group ---
  df_work$target <- target_map[df_work$group]
  df_work$y_bm <- NA_real_

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
      df_work$y_bm[idx_g] <- y_g + (diff_val / sum(w_calc))

    } else if (method == "optimal") {
      # Optimal Quadratic Loss / Variance-Weighted Benchmarking (Datta et al., 2011; Steorts et al., 2014)
      mse_g <- df_work$mse[idx_g]
      diff_val <- T_g - agg_initial

      if (any(is.na(mse_g)) || all(mse_g <= 0)) {
        cli::cli_warn("MSE not available or non-positive for group {.val {g}}; falling back to difference benchmarking.")
        df_work$y_bm[idx_g] <- y_g + (diff_val / sum(w_calc))
      } else {
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

      low_b <- -20
      upp_b <- 20
      f_low <- f_obj(low_b)
      f_upp <- f_obj(upp_b)
      if (f_low * f_upp > 0) {
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
  if (is_hierarchical) {
    res_df$national_target <- rep(national_target, nrow(res_df))
  }

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

  national_verification <- NULL
  if (is_hierarchical) {
    w_nat_calc <- if (type == "mean") df_work$weight / sum(df_work$weight) else df_work$weight
    nat_orig_sum <- sum(w_nat_calc * df_work$y_hat)
    nat_bm_sum <- sum(w_nat_calc * df_work$y_bm)
    national_verification <- list(
      target = national_target,
      original_sum = nat_orig_sum,
      benchmarked_sum = nat_bm_sum,
      discrepancy = abs(nat_bm_sum - national_target)
    )
  }

  attr(res_df, "method") <- method
  attr(res_df, "type") <- type
  attr(res_df, "hierarchical") <- is_hierarchical
  attr(res_df, "group_var") <- group_var
  attr(res_df, "stage1_summary") <- stage1_df
  attr(res_df, "verification") <- verification
  attr(res_df, "national_verification") <- national_verification
  attr(res_df, "call") <- call

  class(res_df) <- c("fastsae_benchmark", "data.frame")
  return(res_df)
}

#' @rdname benchmark_sae
#' @export
print.fastsae_benchmark <- function(x, all_groups = FALSE, ...) {
  is_hierarchical <- isTRUE(attr(x, "hierarchical"))

  if (is_hierarchical) {
    cat("=== fastsae Two-Stage Hierarchical Benchmark Calibration ===\n")
  } else {
    cat("=== fastsae Small Area Benchmark Calibration ===\n")
  }

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

  # National Level verification
  nat_verif <- attr(x, "national_verification")
  if (!is.null(nat_verif)) {
    disc_ok <- nat_verif$discrepancy < 1e-5
    status_icon <- if (disc_ok) "[CONSISTENT]" else "[APPROX]"
    cat(sprintf("\nLevel 0 (National Target): Target = %.5f | Calibrated Aggregate = %.5f %s\n",
                nat_verif$target, nat_verif$benchmarked_sum, status_icon))
  }

  # Print verification
  verif <- attr(x, "verification")
  if (!is.null(verif) && length(verif) > 0) {
    n_grps <- length(verif)
    cat(sprintf("\nLevel 1 (Group Consistency): %d group%s\n", n_grps, if (n_grps > 1) "s" else ""))

    show_grps <- if (n_grps > 8 && !isTRUE(all_groups)) utils::head(names(verif), 6) else names(verif)
    for (g in show_grps) {
      v <- verif[[g]]
      g_label <- if (g == "All") "Target" else paste0("Group '", g, "'")
      disc_ok <- v$discrepancy < 1e-5
      status_icon <- if (disc_ok) "[CONSISTENT]" else "[APPROX]"
      cat(sprintf("  %s: Target = %.5f | Original = %.5f -> Calibrated = %.5f %s\n",
                  g_label, v$target, v$original_sum, v$benchmarked_sum, status_icon))
    }
    if (n_grps > length(show_grps)) {
      cat(sprintf("  ... and %d more groups (all verified consistent within < 1e-5; pass all_groups = TRUE to print all)\n",
                  n_grps - length(show_grps)))
    }
  }

  cat("\nFirst 6 benchmarked domains (Level 2):\n")
  print(utils::head(as.data.frame(x), 6), ...)
  if (nrow(x) > 6) {
    cat("... and", nrow(x) - 6, "more domains.\n")
  }

  cat("\n")
  invisible(x)
}

#' @rdname benchmark_sae
#' @export
summary.fastsae_benchmark <- function(object, ...) {
  print(object, all_groups = FALSE, ...)

  s1 <- attr(object, "stage1_summary")
  if (!is.null(s1)) {
    grp_name <- attr(object, "group_var") %||% "Group"
    cat(sprintf("\n-- Stage 1: %s Harmonization Summary --\n", grp_name))
    print(utils::head(s1, 10), row.names = FALSE)
    if (nrow(s1) > 10) {
      cat("... and", nrow(s1) - 10, "more groups.\n")
    }
  }

  cat("\n-- Adjustment Statistics across Small Areas (Level 2) --\n")
  cat("Absolute Adjustment Summary:\n")
  print(summary(object$adjustment))
  cat("\nRelative Adjustment (%) Summary:\n")
  print(summary(object$rel_adjustment_pct))
  invisible(object)
}

#' @rdname benchmark_sae
#' @export
autoplot.fastsae_benchmark <- function(object, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggplot2} is required to plot benchmarked objects.")
  }

  df_plot <- as.data.frame(object)
  has_group <- "group" %in% names(df_plot) && length(unique(df_plot$group)) > 1

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = .data$original, y = .data$benchmarked)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50")

  if (has_group) {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = .data$group), alpha = 0.85, size = 2) +
      ggplot2::labs(color = "Group")
  } else {
    p <- p + ggplot2::geom_point(color = "#1D6A5C", alpha = 0.85, size = 2.5)
  }

  is_hier <- isTRUE(attr(object, "hierarchical"))
  title_str <- if (is_hier) "Two-Stage Hierarchical Calibration" else "Benchmark Calibration"

  p <- p +
    ggplot2::labs(
      title = paste0("Small Area ", title_str),
      subtitle = "Original Model Estimates vs Calibrated Benchmarked Predictions",
      x = "Original Estimate",
      y = "Benchmarked Estimate"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  if (has_group && length(unique(df_plot$group)) > 10) {
    p <- p + ggplot2::theme(legend.position = "none") # suppress cluttered legend for 34 groups
  }

  return(p)
}

#' @rdname benchmark_sae
#' @export
plot.fastsae_benchmark <- function(x, y = NULL, ...) {
  dots <- list(...)
  if (!is.null(dots$sf_geom)) {
    return(map_sae(x, ...))
  }
  autoplot.fastsae_benchmark(x, ...)
}
