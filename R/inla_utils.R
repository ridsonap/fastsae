# ============================================================================
# Utility functions for INLA integration in fastsae
# ============================================================================

#' Check if INLA package is available
#'
#' @noRd
.check_inla_installed <- function() {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg INLA} is required for this estimation method.",
      "i" = "Install it using: {.code install.packages('INLA', repos = c(getOption('repos'), INLA = 'https://inla.r-inla-download.org/R/stable'))}"
    ))
  }
}

#' Convert and validate spatial proximity / weight matrix for INLA models
#'
#' @param W Matrix, Matrix object, `spdep` `nb`, or `listw` object.
#' @param n_domains Expected number of unique domains.
#' @param spatial Spatial model type ("bym2", "bym", "besag", "generic1", "slm", "none").
#' @param domain_names Optional vector of unique domain names for alignment.
#'
#' @return A list with `graph` (symmetric adjacency matrix or inla graph object)
#'   and `W_mat` (numeric matrix for SLM/Leroux).
#' @noRd
.convert_spatial_weights <- function(W, n_domains, spatial = "bym2", domain_names = NULL) {
  if (is.null(W)) {
    cli::cli_abort("Spatial weight/proximity matrix {.arg W} must be provided when {.code spatial != 'none'}.")
  }

  # If W is an spdep nb or listw object
  if (inherits(W, "listw")) {
    if (requireNamespace("spdep", quietly = TRUE)) {
      W_mat <- spdep::listw2mat(W)
    } else {
      cli::cli_abort("Package {.pkg spdep} is required to convert a {.cls listw} object.")
    }
  } else if (inherits(W, "nb")) {
    if (requireNamespace("spdep", quietly = TRUE)) {
      W_mat <- spdep::nb2mat(W, style = "B", zero.policy = TRUE)
    } else {
      cli::cli_abort("Package {.pkg spdep} is required to convert an {.cls nb} object.")
    }
  } else if (is.matrix(W) || inherits(W, "Matrix")) {
    W_mat <- as.matrix(W)
  } else if (inherits(W, "inla.graph")) {
    return(list(graph = W, W_mat = NULL))
  } else {
    cli::cli_abort("Unsupported spatial object type for {.arg W}. Expected matrix, Matrix, nb, or listw.")
  }

  # If domain_names provided and W has row/colnames, reorder W to match domain order
  if (!is.null(domain_names) && !is.null(rownames(W_mat))) {
    dom_chr <- as.character(domain_names)
    if (all(dom_chr %in% rownames(W_mat))) {
      W_mat <- W_mat[dom_chr, dom_chr, drop = FALSE]
    }
  }

  # Validate dimensions
  if (nrow(W_mat) != n_domains || ncol(W_mat) != n_domains) {
    cli::cli_abort(c(
      "{.arg W} must be a square matrix with dimensions equal to total number of domains ({n_domains} x {n_domains}).",
      "x" = "Got a {nrow(W_mat)} x {ncol(W_mat)} matrix."
    ))
  }

  # For graph-based spatial models (BYM2, BYM, Besag, Leroux/generic1)
  # INLA requires a symmetric binary adjacency matrix (or symmetric weights)
  # with zero on the diagonal.
  adj_mat <- (W_mat > 0 | t(W_mat) > 0) * 1
  diag(adj_mat) <- 0

  # Check if disconnected components or isolated nodes exist
  degree <- rowSums(adj_mat)
  if (any(degree == 0)) {
    isolated_idx <- which(degree == 0)
    cli::cli_warn(
      "{length(isolated_idx)} domain(s) have no spatial neighbors in {.arg W} (isolated nodes): {paste(isolated_idx, collapse = ', ')}."
    )
  }

  return(list(
    graph = adj_mat,
    W_mat = W_mat
  ))
}

#' Extract and format results from an INLA model fit
#'
#' @noRd
.extract_inla_results <- function(
  fit,
  data,
  y,
  domain,
  time = NULL,
  family,
  spatial,
  temporal = "none",
  st_interaction = "none",
  vardir = NULL,
  trials = NULL,
  exposure = NULL,
  X_mat = NULL,
  rho_range = NULL,
  domain_id = NULL,
  time_id = NULL,
  unique_domains = NULL,
  unique_times = NULL,
  call = NULL
) {
  n_obs <- length(y)
  n_domains <- if (!is.null(unique_domains)) length(unique_domains) else length(unique(domain))
  n_times <- if (!is.null(unique_times)) length(unique_times) else (if (!is.null(time)) length(unique(time)) else 1)

  # 1. Fixed effects coefficients
  if (spatial == "slm" && !is.null(X_mat)) {
    p_cov <- ncol(X_mat)
    colnames_X <- colnames(X_mat)
    rf_df <- as.data.frame(fit$summary.random[[1]])
    coef_rows <- rf_df[(n_domains + 1):(n_domains + p_cov), ]
    estcoef <- data.frame(
      beta = coef_rows$mean,
      std.error = coef_rows$sd,
      zvalue = coef_rows$mean / ifelse(coef_rows$sd > 0, coef_rows$sd, NA_real_),
      pvalue = 2 * stats::pnorm(abs(coef_rows$mean / ifelse(coef_rows$sd > 0, coef_rows$sd, NA_real_)), lower.tail = FALSE),
      ci_lower = coef_rows[["0.025quant"]],
      ci_upper = coef_rows[["0.975quant"]],
      row.names = colnames_X
    )
  } else {
    fixed_summary <- as.data.frame(fit$summary.fixed)
    estcoef <- data.frame(
      beta = fixed_summary$mean,
      std.error = fixed_summary$sd,
      zvalue = fixed_summary$mean / ifelse(fixed_summary$sd > 0, fixed_summary$sd, NA_real_),
      pvalue = 2 * stats::pnorm(abs(fixed_summary$mean / ifelse(fixed_summary$sd > 0, fixed_summary$sd, NA_real_)), lower.tail = FALSE),
      ci_lower = fixed_summary[["0.025quant"]],
      ci_upper = fixed_summary[["0.975quant"]],
      row.names = rownames(fixed_summary)
    )
  }

  # 2. Hyperparameters
  hyper_summary <- if (!is.null(fit$summary.hyperpar)) as.data.frame(fit$summary.hyperpar) else NULL

  # Extract variance component(s) and spatial/temporal correlation parameters
  random_effect_var <- NULL
  random_effect_var_time <- NULL
  phi <- NULL
  rho <- NULL
  rho_time <- NULL

  if (!is.null(hyper_summary) && nrow(hyper_summary) > 0) {
    # Spatial / Area random effect precision
    prec_rows <- grep("Precision for (\\.\\.domain_id\\.\\.|domain)", rownames(hyper_summary), value = TRUE)
    if (length(prec_rows) == 0 && spatial == "none" && temporal == "none") {
      # Fallback for simple IID
      prec_rows <- grep("Precision for", rownames(hyper_summary), value = TRUE)
    }
    if (length(prec_rows) > 0) {
      prec_est <- hyper_summary[prec_rows[1], "mean"]
      if (!is.na(prec_est) && prec_est > 0) {
        random_effect_var <- 1 / prec_est
      }
    }

    # Temporal precision
    prec_time_rows <- grep("Precision for (\\.\\.time_id\\.\\.|time)", rownames(hyper_summary), value = TRUE)
    if (length(prec_time_rows) > 0) {
      prec_time_est <- hyper_summary[prec_time_rows[1], "mean"]
      if (!is.na(prec_time_est) && prec_time_est > 0) {
        random_effect_var_time <- 1 / prec_time_est
      }
    }

    # Look for Phi (spatial proportion in BYM2), Beta (generic1 / Leroux), or Rho (slm)
    phi_rows <- grep("Phi for", rownames(hyper_summary), value = TRUE)
    beta_rows <- grep("Beta for", rownames(hyper_summary), value = TRUE)
    rho_rows <- grep("Rho for (\\.\\.domain_id\\.\\.|domain)", rownames(hyper_summary), value = TRUE)
    rho_time_rows <- grep("Rho for (\\.\\.time_id\\.\\.|time)", rownames(hyper_summary), value = TRUE)

    if (length(phi_rows) > 0) {
      phi <- hyper_summary[phi_rows[1], "mean"]
      rho <- phi
    } else if (length(beta_rows) > 0) {
      # In generic1 / Leroux CAR, Beta represents spatial autocorrelation parameter rho
      phi <- hyper_summary[beta_rows[1], "mean"]
      rho <- phi
    } else if (length(rho_rows) > 0) {
      # In SLM, Rho represents spatial autocorrelation
      rho_raw <- hyper_summary[rho_rows[1], "mean"]
      if (!is.null(rho_range)) {
        rho <- rho_range[1] + rho_raw * (rho_range[2] - rho_range[1])
      } else {
        rho <- rho_raw
      }
      phi <- rho
    }

    if (length(rho_time_rows) > 0) {
      rho_time <- hyper_summary[rho_time_rows[1], "mean"]
    }
  }

  # 3. Random effects extraction per observation
  rand_eff <- NULL
  rand_eff_spatial <- NULL
  rand_eff_temporal <- NULL

  if (!is.null(fit$summary.random) && length(fit$summary.random) > 0) {
    rand_names <- names(fit$summary.random)

    # Spatial random effect
    if ("..domain_id.." %in% rand_names) {
      rf_df_spat <- fit$summary.random[["..domain_id.."]]
      if ("mean" %in% names(rf_df_spat)) {
        spat_vals <- rf_df_spat$mean[seq_len(n_domains)]
        if (!is.null(domain_id)) {
          rand_eff_spatial <- spat_vals[domain_id]
        } else if (length(spat_vals) == n_obs) {
          rand_eff_spatial <- spat_vals
        }
      }
    }

    # Temporal random effect
    if ("..time_id.." %in% rand_names) {
      rf_df_time <- fit$summary.random[["..time_id.."]]
      if ("mean" %in% names(rf_df_time)) {
        if (st_interaction == "domain-specific" || length(rf_df_time$mean) == n_domains * n_times) {
          # Domain-specific temporal dynamics (indexed by (time, domain) or group)
          # In INLA, group = domain_id with time_id creates n_times * n_domains rows
          # INLA group ordering: for each domain, all time periods
          if (!is.null(domain_id) && !is.null(time_id)) {
            # INLA groups time_id by domain_id: index is (domain_id - 1) * n_times + time_id
            group_idx <- (domain_id - 1) * n_times + time_id
            rand_eff_temporal <- rf_df_time$mean[group_idx]
          }
        } else {
          # Main temporal trend (length n_times)
          time_vals <- rf_df_time$mean[seq_len(n_times)]
          if (!is.null(time_id)) {
            rand_eff_temporal <- time_vals[time_id]
          }
        }
      }
    }

    # Total random effect
    if (!is.null(rand_eff_spatial) && !is.null(rand_eff_temporal)) {
      rand_eff <- rand_eff_spatial + rand_eff_temporal
    } else if (!is.null(rand_eff_spatial)) {
      rand_eff <- rand_eff_spatial
    } else if (!is.null(rand_eff_temporal)) {
      rand_eff <- rand_eff_temporal
    } else {
      # Fallback to first random effect
      rf_first <- fit$summary.random[[1]]$mean
      if (length(rf_first) >= n_obs) {
        rand_eff <- rf_first[seq_len(n_obs)]
      } else if (!is.null(domain_id) && length(rf_first) >= n_domains) {
        rand_eff <- rf_first[domain_id]
      }
    }
  }

  if (is.null(rand_eff) || length(rand_eff) != n_obs) {
    rand_eff <- rep(NA_real_, n_obs)
  }

  # 4. Fitted values (Predictions on original response scale)
  fitted_summary <- as.data.frame(fit$summary.fitted.values[seq_len(n_obs), , drop = FALSE])
  linpred_summary <- as.data.frame(fit$summary.linear.predictor[seq_len(n_obs), , drop = FALSE])

  ebp_est <- fitted_summary$mean
  ebp_sd  <- fitted_summary$sd
  ebp_mse <- ebp_sd^2
  ebp_rse <- ifelse(abs(ebp_est) < .Machine$double.eps, NA_real_, (ebp_sd / abs(ebp_est)) * 100)
  ci_lower <- fitted_summary[["0.025quant"]]
  ci_upper <- fitted_summary[["0.975quant"]]
  linpred  <- linpred_summary$mean

  # Assemble df_ebp
  if (!is.null(time)) {
    df_ebp <- data.frame(
      domain = domain,
      time = time,
      y = y,
      ebp = ebp_est,
      linear_pred = linpred,
      sd = ebp_sd,
      mse = ebp_mse,
      rse = ebp_rse,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      random_effect = rand_eff,
      stringsAsFactors = FALSE
    )
    if (!is.null(rand_eff_spatial)) df_ebp$random_effect_spatial <- rand_eff_spatial
    if (!is.null(rand_eff_temporal)) df_ebp$random_effect_temporal <- rand_eff_temporal
  } else {
    df_ebp <- data.frame(
      domain = domain,
      y = y,
      ebp = ebp_est,
      linear_pred = linpred,
      sd = ebp_sd,
      mse = ebp_mse,
      rse = ebp_rse,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      random_effect = rand_eff,
      stringsAsFactors = FALSE
    )
  }

  # If vardir provided (Gaussian FH, Beta, or Gamma)
  if (!is.null(vardir)) {
    df_ebp$vardir <- vardir
    if (family == "beta") {
      phi_dir <- (y * (1 - y) / vardir) - 1
      df_ebp$precision <- pmax(phi_dir, 1, na.rm = TRUE)
    } else if (family == "gamma") {
      s_dir <- (y^2) / vardir
      df_ebp$precision <- s_dir
      df_ebp$cv_dir <- sqrt(vardir) / y
    }
  }

  # If trials provided (Binomial or Beta)
  if (!is.null(trials)) {
    df_ebp$trials <- trials
    if (family == "binomial") {
      df_ebp$estimated_total <- ebp_est * trials
    } else if (family == "beta" && is.null(vardir)) {
      df_ebp$precision <- pmax(trials - 1, 1, na.rm = TRUE)
    }
  }

  # If exposure provided (Poisson / NegBinom)
  if (!is.null(exposure)) {
    df_ebp$exposure <- exposure
    df_ebp$rate <- ebp_est
    df_ebp$estimated_count <- ebp_est * exposure
  }

  # 5. Goodness of fit measures
  goodness <- c(
    DIC = if (!is.null(fit$dic$dic)) fit$dic$dic else NA_real_,
    pD = if (!is.null(fit$dic$p.eff)) fit$dic$p.eff else NA_real_,
    WAIC = if (!is.null(fit$waic$waic)) fit$waic$waic else NA_real_,
    pWAIC = if (!is.null(fit$waic$p.eff)) fit$waic$p.eff else NA_real_,
    Marginal_LogLik = if (!is.null(fit$mlik)) fit$mlik[1, 1] else NA_real_
  )

  # Format model description
  model_label <- paste0("EBP-", toupper(family))
  if (temporal != "none" || spatial != "none") {
    comps <- c()
    if (spatial != "none") comps <- c(comps, toupper(spatial))
    if (temporal != "none") comps <- c(comps, toupper(temporal))
    if (st_interaction != "none") comps <- c(comps, paste0("ST:", toupper(st_interaction)))
    model_label <- paste0(model_label, " (", paste(comps, collapse = " + "), ")")
  } else {
    model_label <- paste0(model_label, " (Non-spatial)")
  }

  out <- list(
    df_ebp = df_ebp,
    ebp = df_ebp, # backward compatibility with fastsae conventions
    df_eblup = df_ebp, # compatibility with eblup methods
    estcoef = estcoef,
    hyperpar = hyper_summary,
    random_effect_var = random_effect_var,
    random_effect_var_time = random_effect_var_time,
    phi = phi,
    rho = rho,
    rho_time = rho_time,
    goodness = goodness,
    family = family,
    spatial = spatial,
    temporal = temporal,
    st_interaction = st_interaction,
    level = "area",
    model = model_label,
    convergence = TRUE,
    fit = fit,
    call = call
  )

  class(out) <- c("fastsae_ebp_area", "fastsae")
  return(out)
}

