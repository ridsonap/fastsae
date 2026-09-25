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
#' @param n_domains Expected number of domains.
#' @param spatial Spatial model type ("bym2", "bym", "besag", "generic1", "slm", "none").
#'
#' @return A list with `graph` (symmetric adjacency matrix or inla graph object)
#'   and `W_mat` (numeric matrix for SLM/Leroux).
#' @noRd
.convert_spatial_weights <- function(W, n_domains, spatial = "bym2") {
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
  family,
  spatial,
  vardir = NULL,
  trials = NULL,
  exposure = NULL,
  X_mat = NULL,
  rho_range = NULL,
  call = NULL
) {
  n_domains <- length(domain)

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
  
  # Extract variance component(s) and spatial correlation parameters
  random_effect_var <- NULL
  phi <- NULL
  rho <- NULL
  if (!is.null(hyper_summary) && nrow(hyper_summary) > 0) {
    # Look for domain precision
    prec_rows <- grep("Precision for", rownames(hyper_summary), value = TRUE)
    if (length(prec_rows) > 0) {
      prec_est <- hyper_summary[prec_rows[1], "mean"]
      if (!is.na(prec_est) && prec_est > 0) {
        random_effect_var <- 1 / prec_est
      }
    }
    # Look for Phi (spatial proportion in BYM2), Beta (generic1 / Leroux), or Rho (slm)
    phi_rows <- grep("Phi for", rownames(hyper_summary), value = TRUE)
    beta_rows <- grep("Beta for", rownames(hyper_summary), value = TRUE)
    rho_rows <- grep("Rho for", rownames(hyper_summary), value = TRUE)

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
  }

  # 3. Random effects per domain
  rand_eff <- NULL
  if (!is.null(fit$summary.random) && length(fit$summary.random) > 0) {
    rf_df <- fit$summary.random[[1]]
    if ("mean" %in% names(rf_df)) {
      # In BYM2, INLA returns 2 * n_domains rows (first n is marginal effect, second is spatial)
      # In SLM, INLA returns n_domains + p rows (first n is spatial effect, rest are beta coefficients)
      rand_eff <- rf_df$mean[seq_len(min(nrow(rf_df), n_domains))]
    }
  }
  if (is.null(rand_eff) || length(rand_eff) != n_domains) {
    rand_eff <- rep(NA_real_, n_domains)
  }

  # 4. Fitted values (Predictions on original response scale)
  fitted_summary <- as.data.frame(fit$summary.fitted.values[seq_len(n_domains), , drop = FALSE])
  linpred_summary <- as.data.frame(fit$summary.linear.predictor[seq_len(n_domains), , drop = FALSE])

  ebp_est <- fitted_summary$mean
  ebp_sd  <- fitted_summary$sd
  ebp_mse <- ebp_sd^2
  ebp_rse <- ifelse(abs(ebp_est) < .Machine$double.eps, NA_real_, (ebp_sd / abs(ebp_est)) * 100)
  ci_lower <- fitted_summary[["0.025quant"]]
  ci_upper <- fitted_summary[["0.975quant"]]
  linpred  <- linpred_summary$mean

  # Assemble df_ebp
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

  out <- list(
    df_ebp = df_ebp,
    ebp = df_ebp, # backward compatibility with fastsae conventions
    df_eblup = df_ebp, # compatibility with eblup methods
    estcoef = estcoef,
    hyperpar = hyper_summary,
    random_effect_var = random_effect_var,
    phi = phi,
    rho = rho,
    goodness = goodness,
    family = family,
    spatial = spatial,
    level = "area",
    model = paste0("EBP-", toupper(family), if (spatial != "none") paste0(" (", toupper(spatial), ")") else " (Non-spatial)"),
    convergence = TRUE,
    fit = fit,
    call = call
  )

  class(out) <- c("fastsae_ebp_area", "fastsae")
  return(out)
}
