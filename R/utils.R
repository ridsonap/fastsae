# ============================================================================
# Shared utility functions for fastsae package
# ============================================================================

utils::globalVariables(c("density", ".data"))

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
#' @noRd
.get_variable <- function(data, variable) {
  if (is.character(variable) && length(variable) == 1) {
    if (variable %in% colnames(data)) {
      return(data[[variable]])
    } else {
      cli::cli_abort('variable "{variable}" is not found in the data')
    }
  } else if (inherits(variable, "formula")) {
    v_names <- all.vars(variable)
    if (length(v_names) == 1 && v_names %in% colnames(data)) {
      return(data[[v_names]])
    } else {
      cli::cli_abort("formula does not reference a valid single column in data")
    }
  } else if (length(variable) == nrow(data)) {
    return(variable)
  } else {
    cli::cli_abort("variable is not valid or length does not match data ({length(variable)} vs {nrow(data)})")
  }
}

#' Setup self-benchmarking data structures via Wang-Fuller-Qu (2008) augmented model
#'
#' @param data A data frame containing the model variables.
#' @param formula The model formula.
#' @param y Direct estimates response vector.
#' @param vardir Sampling variances vector.
#' @param weight Benchmark weight (vector, column name, or formula).
#' @param target Benchmark target (numeric scalar, named vector, or NULL).
#' @param group Group identifier (vector, column name, or formula, or NULL).
#'
#' @return A list with data_aug, formula_aug, z_names, w_norm, group_vec, unique_groups, target_vec, shift, original_y.
#' @noRd
.setup_self_benchmark <- function(
  data,
  formula,
  y,
  vardir = NULL,
  weight = NULL,
  target = NULL,
  group = NULL
) {
  n_obs <- nrow(data)
  response_var <- all.vars(formula)[1]

  # 1. Weights extraction & validation
  w_raw <- if (!is.null(weight)) .get_variable(data, weight) else rep(1, n_obs)
  if (length(w_raw) != n_obs) {
    cli::cli_abort("Length of {.arg benchmark_weight} ({length(w_raw)}) must match number of rows in {.arg data} ({n_obs}).")
  }
  if (any(w_raw < 0, na.rm = TRUE)) {
    cli::cli_abort("{.arg benchmark_weight} must be non-negative.")
  }
  w_raw[is.na(w_raw)] <- 0

  # 2. Group extraction & validation
  if (!is.null(group)) {
    group_vec <- as.character(.get_variable(data, group))
  } else {
    group_vec <- rep("All", n_obs)
  }
  unique_groups <- unique(group_vec)

  # 3. Normalize weights per group over sampled domains (!is.na(y))
  w_norm <- rep(0, n_obs)
  target_vec <- numeric(length(unique_groups))
  names(target_vec) <- unique_groups
  direct_agg_vec <- numeric(length(unique_groups))
  names(direct_agg_vec) <- unique_groups
  shift <- rep(0, n_obs)

  for (i in seq_along(unique_groups)) {
    g <- unique_groups[i]
    idx_g <- which(group_vec == g & !is.na(y))
    if (length(idx_g) == 0) {
      next
    }
    sum_w_g <- sum(w_raw[idx_g])
    if (sum_w_g <= 0) {
      cli::cli_abort("Sum of {.arg benchmark_weight} for sampled domains in group {.val {g}} must be positive.")
    }
    w_norm[idx_g] <- w_raw[idx_g] / sum_w_g
    dir_agg <- sum(w_norm[idx_g] * y[idx_g])
    direct_agg_vec[g] <- dir_agg

    # Determine target for group g
    if (is.null(target)) {
      tgt <- dir_agg
    } else if (is.numeric(target)) {
      if (!is.null(names(target)) && g %in% names(target)) {
        tgt <- target[g]
      } else if (length(target) == 1) {
        tgt <- target
      } else if (length(target) == length(unique_groups)) {
        tgt <- target[i]
      } else {
        cli::cli_abort("Dimensions of {.arg benchmark_target} do not match groups.")
      }
    } else {
      cli::cli_abort("{.arg benchmark_target} must be numeric or NULL.")
    }
    target_vec[g] <- tgt
    shift[idx_g] <- tgt - dir_agg
  }

  # 4. Construct augmented data and shifted response if target != dir_agg
  data_aug <- data
  if (any(abs(shift) > 1e-12)) {
    data_aug[[response_var]] <- data_aug[[response_var]] + shift
  }

  # 5. Augmented covariate construction (Wang, Fuller, and Qu, 2008)
  # z_i = w_i * D_i
  vd <- if (!is.null(vardir)) vardir else rep(1, n_obs)
  vd_clean <- ifelse(is.na(vd) | !is.finite(vd), 0, vd)

  z_names <- character()
  if (length(unique_groups) > 1) {
    for (i in seq_along(unique_groups)) {
      g <- unique_groups[i]
      col_z <- paste0("bench_aug_", make.names(g))
      data_aug[[col_z]] <- ifelse(group_vec == g, w_norm * vd_clean, 0)
      z_names <- c(z_names, col_z)
    }
  } else {
    col_z <- "bench_aug"
    data_aug[[col_z]] <- w_norm * vd_clean
    z_names <- col_z
  }

  # 6. Updated formula (ponytail: preserve terms env via reformulate; avoids dropping offset/intercept)
  formula_terms <- attr(stats::terms(formula), "term.labels")
  # Detect if original had intercept
  has_intercept <- attr(stats::terms(formula), "intercept") == 1
  int_str <- if (has_intercept) NULL else "-1"
  all_fixed <- c(formula_terms, z_names)
  # Use reformulate to keep environment attached
  rhs_parts <- c(all_fixed, int_str)
  rhs_parts <- rhs_parts[!sapply(rhs_parts, is.null)]
  if (length(rhs_parts) == 0) rhs_parts <- "1"
  formula_aug <- stats::reformulate(rhs_parts, response = response_var, env = environment(formula) %||% parent.frame())
  # Ensure intercept handling
  if (!has_intercept) {
    # reformulate with intercept=FALSE not directly; adjust via update if needed
    formula_aug <- stats::update(formula_aug, ~ . -1)
  }

  list(
    data_aug = data_aug,
    formula_aug = formula_aug,
    z_names = z_names,
    w_norm = w_norm,
    group_vec = group_vec,
    unique_groups = unique_groups,
    target_vec = target_vec,
    direct_agg_vec = direct_agg_vec,
    shift = shift,
    original_y = y
  )
}

#' Compute benchmark summary table
#' @noRd
.compute_benchmark_summary <- function(
  unique_groups,
  group_vec,
  w_norm,
  y,
  target_vec,
  estimates
) {
  res_list <- vector("list", length(unique_groups))
  for (i in seq_along(unique_groups)) {
    g <- unique_groups[i]
    idx_g <- which(group_vec == g & !is.na(y))
    dir_agg <- if (length(idx_g) > 0) sum(w_norm[idx_g] * y[idx_g]) else NA_real_
    mod_agg <- if (length(idx_g) > 0) sum(w_norm[idx_g] * estimates[idx_g]) else NA_real_
    tgt <- target_vec[g]
    res_list[[i]] <- data.frame(
      group = g,
      target = tgt,
      direct_aggregate = dir_agg,
      model_aggregate = mod_agg,
      discrepancy = mod_agg - tgt,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, res_list)
}
