# ============================================================================
# Utility functions for INLA integration in fastsae
# ============================================================================

# Is the package 'sn' (needed by INLA::inla.posterior.sample) available?
.has_sn <- function() requireNamespace("sn", quietly = TRUE)

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

#' Posterior mean of a variance from the marginal of its precision
#'
#' E[1 / tau | y] (the posterior mean of the variance), not 1 / E[tau | y],
#' which is smaller by Jensen's inequality. Falls back to 1 / E[tau | y] when
#' the marginal is not available.
#' @noRd
.post_mean_variance <- function(fit, row_name, prec_mean) {
  marg <- fit$marginals.hyperpar[[row_name]]
  v <- if (!is.null(marg)) {
    tryCatch(INLA::inla.emarginal(function(x) 1 / x, marg), error = function(e) NA_real_)
  } else {
    NA_real_
  }
  if (!is.finite(v) && !is.na(prec_mean) && prec_mean > 0) v <- 1 / prec_mean
  v
}


#' Number of hyperparameters of a Gaussian area model with a Knorr-Held interaction
#'
#' Mirrors the terms that hb_area() builds for st_interaction = "type1" to
#' "type4": the area effect (iid without a spatial effect, one precision for
#' besag, two for bym, bym2, Leroux and SLM), the temporal main effect (a
#' precision, plus the correlation of ar1) and the interaction (a precision,
#' plus the correlation of ar1 in type II and IV), and the residual precision
#' when the sampling variances are not given.
#' @noRd
.n_hyper_st <- function(spatial, temporal, st_interaction, vardir_known) {
  n_space <- if (spatial %in% c("none", "besag")) 1L else 2L
  n_time <- if (temporal == "ar1") 2L else 1L
  n_int <- if (temporal == "ar1" && st_interaction %in% c("type2", "type4")) 2L else 1L
  n_space + n_time + n_int + if (vardir_known) 0L else 1L
}

#' Knorr-Held space-time interaction with all identifiability constraints
#'
#' With a random walk in time the interaction is intrinsic in the direction "one
#' level per domain" (rw1) and also "one linear trend per domain" (rw2). INLA's
#' grouped models constrain only within each group, so those directions kept a
#' flat prior and a domain without direct estimates had a posterior SD of 36 or
#' more; extraconstr of a grouped model is also applied per group. The
#' interaction is therefore a generic0 model on the full domain x time index
#' with the scaled Kronecker structure of Knorr-Held (2000) and all constraints:
#'   type4: R_t (x) R_s, element (t - 1) * D + d; sum_d delta_dt = 0 per time and
#'          connected component of the graph, sum_t delta_dt = 0 (and
#'          sum_t (t - mean(t)) delta_dt = 0 for rw2) per domain, one domain per
#'          component dropped (those rows are implied by the others);
#'   type2: I_D (x) R_t, element (d - 1) * T + t; per-domain sum (and trend).
#' Each factor is scaled to a generalised variance of 1 as with scale.model:
#' R_t with its random-walk constraints, R_s per connected component. A domain
#' without neighbours gets a spatial factor of 1 (an iid term).
#' @return NULL when the grouped INLA model is already identified, else a list
#'   with Cmatrix, rankdef, extraconstr and the index of every observation.
#' @noRd
.st_generic0 <- function(st_interaction, temporal, domain_id, time_id, n_domains, n_times,
                         graph = NULL, scale_model = TRUE) {
  if (!temporal %in% c("rw1", "rw2")) {
    return(NULL)
  }
  if (!(st_interaction == "type4" || (st_interaction == "type2" && temporal == "rw2"))) {
    return(NULL)
  }
  D <- n_domains
  Tt <- n_times
  k_t <- if (temporal == "rw1") 1L else 2L
  tc <- seq_len(Tt) - (Tt + 1) / 2
  Dmat <- diff(diag(Tt), differences = k_t)
  R_t <- crossprod(Dmat)
  if (isTRUE(scale_model)) {
    A_t <- if (k_t == 1L) matrix(1, 1, Tt) else rbind(rep(1, Tt), tc)
    R_t <- as.matrix(INLA::inla.scale.model(R_t, constr = list(A = A_t, e = rep(0, nrow(A_t)))))
  }
  row_dt <- function(d, w) {
    a <- numeric(D * Tt)
    a[if (st_interaction == "type4") (seq_len(Tt) - 1) * D + d else (d - 1) * Tt + seq_len(Tt)] <- w
    a
  }
  if (st_interaction == "type4") {
    adj <- (as.matrix(graph) > 0) * 1
    diag(adj) <- 0
    R_s <- diag(rowSums(adj)) - adj
    comp <- .graph_components(adj)
    sizes <- table(comp)
    multi <- as.integer(names(sizes)[sizes > 1])
    for (k in unique(comp)) {
      m <- which(comp == k)
      if (length(m) == 1) {
        R_s[m, m] <- 1 # a domain without neighbours: iid term
      } else if (isTRUE(scale_model)) {
        R_s[m, m] <- as.matrix(INLA::inla.scale.model(R_s[m, m, drop = FALSE],
          constr = list(A = matrix(1, 1, length(m)), e = 0)
        ))
      }
    }
    Cmat <- kronecker(R_t, R_s)
    rows_time <- unlist(lapply(seq_len(Tt), function(t) {
      lapply(multi, function(k) {
        a <- numeric(D * Tt)
        a[(t - 1) * D + which(comp == k)] <- 1
        a
      })
    }), recursive = FALSE)
    doms <- setdiff(seq_len(D), vapply(multi, function(k) which(comp == k)[1], integer(1)))
    rows_dom <- lapply(doms, row_dt, w = 1)
    if (k_t == 2L) rows_dom <- c(rows_dom, lapply(doms, row_dt, w = tc))
    rows <- c(rows_time, rows_dom)
    rankdef <- k_t * D + Tt * length(multi) - k_t * length(multi)
    index <- (time_id - 1L) * D + domain_id
  } else {
    Cmat <- kronecker(diag(D), R_t)
    rows <- c(lapply(seq_len(D), row_dt, w = 1), lapply(seq_len(D), row_dt, w = tc))
    rankdef <- 2L * D
    index <- (domain_id - 1L) * Tt + time_id
  }
  A <- do.call(rbind, rows)
  list(
    Cmatrix = Matrix::Matrix(Cmat, sparse = TRUE), rankdef = rankdef,
    extraconstr = list(A = A, e = rep(0, nrow(A))), index = index
  )
}

#' Connected components of an adjacency matrix (breadth-first search)
#' @noRd
.graph_components <- function(adj) {
  adj <- as.matrix(adj) > 0
  n <- nrow(adj)
  comp <- integer(n)
  k <- 0L
  for (i in seq_len(n)) {
    if (comp[i] > 0) next
    k <- k + 1L
    queue <- i
    comp[i] <- k
    while (length(queue) > 0) {
      v <- queue[1]
      queue <- queue[-1]
      nb <- which(adj[v, ] & comp == 0)
      comp[nb] <- k
      queue <- c(queue, nb)
    }
  }
  comp
}

#' Posterior draws of an INLA fit with a reproducible seed
#'
#' INLA chooses the hyperparameter configurations of the draws with R's RNG,
#' so a seed given only to INLA::inla.posterior.sample() does not make them
#' reproducible: the seed is also set on the R side, and the session's random
#' number stream is restored afterwards. seed = 0 is INLA's random seed.
#' @noRd
.inla_draws <- function(fit, n, seed = 1L) {
  seed <- as.integer(seed)
  .with_seed(
    if (seed > 0) seed else -1L,
    INLA::inla.posterior.sample(n = n, result = fit, seed = seed)
  )
}

#' Posterior draws of the linear predictor from INLA's hyperparameter configurations
#'
#' With a seed, INLA::inla.posterior.sample() samples serially and runs the inla
#' program once per configuration (~0.3 s each; 36-143 configurations in
#' hb_twofold()). This picks the configurations as INLA does (R's RNG, by their
#' posterior weight) and draws the latent field of each in R from its Gaussian
#' approximation N(improved mean, Q^-1), conditioned on the linear constraints,
#' then maps it to the predictor (A x in INLA's compact mode, the Predictor rows
#' in classic mode). INLA's skewness correction of the marginals is not applied:
#' on hb_twofold() fits it moved the area aggregates less than the Monte Carlo
#' error of 20000 draws. Returns the n_pred x n matrix, or NULL when the fit has
#' no configurations (control.compute = list(config = TRUE)).
#' @noRd
.inla_predictor_draws <- function(fit, n, n_pred, seed = 1L) {
  cs <- fit$misc$configs
  if (is.null(cs) || length(cs$config) == 0L) {
    return(NULL)
  }
  A <- cs$A
  pred_rows <- NULL
  if (is.null(A) || nrow(A) == 0L) {
    i <- match("Predictor", cs$contents$tag)
    if (is.na(i)) {
      return(NULL)
    }
    pred_rows <- cs$contents$start[i] + seq_len(n_pred) - 1L
  } else if (nrow(A) != n_pred) {
    return(NULL)
  }
  seed <- as.integer(seed)
  .with_seed(if (seed > 0) seed else -1L, {
    lp <- vapply(cs$config, function(cfg) cfg$log.posterior, numeric(1))
    k_of <- sample.int(length(lp), n, replace = TRUE, prob = exp(lp - max(lp)))
    out <- matrix(NA_real_, n_pred, n)
    for (k in sort(unique(k_of))) {
      j <- which(k_of == k)
      cfg <- cs$config[[k]]
      # INLA stores the upper triangle of Q
      L <- Matrix::Cholesky(Matrix::forceSymmetric(cfg$Q, uplo = "U"), perm = TRUE, LDL = FALSE)
      z <- matrix(stats::rnorm(nrow(cfg$Q) * length(j)), nrow(cfg$Q))
      x <- Matrix::solve(L, Matrix::solve(L, z, system = "Lt"), system = "Pt") +
        (cfg$improved.mean %||% cfg$mean)
      if (!is.null(cs$constr)) {
        Ac <- cs$constr$A
        QiA <- Matrix::solve(L, t(Ac))
        x <- x - QiA %*% Matrix::solve(Ac %*% QiA, Ac %*% x - cs$constr$e)
      }
      out[, j] <- as.matrix(if (is.null(pred_rows)) A %*% x else x[pred_rows, , drop = FALSE])
    }
    out
  })
}

#' Convert and validate spatial proximity / weight matrix for INLA models
#'
#' @param W Matrix, Matrix object, `spdep` `nb`, or `listw` object.
#' @param n_domains Expected number of unique domains.
#' @param spatial Spatial model type ("bym2", "bym", "besag", "generic1", "sar", "slm", "none").
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

  # CAR models (besag, bym, bym2, leroux/generic1) use only the pattern W > 0. A
  # dense matrix such as distance-decay weights becomes a (nearly) complete
  # graph, on which the "spatial" effect is almost exchangeable.
  n_adj <- nrow(adj_mat)
  density <- if (n_adj > 1) sum(adj_mat) / (n_adj * (n_adj - 1)) else 0
  if (!spatial %in% c("sar", "slm") && n_adj >= 4 && density > 0.5) {
    cli::cli_warn(c(
      "{.arg W} links {round(100 * density)}% of all pairs of domains. The {.val {spatial}} model uses only the pattern {.code W > 0}, not the weights, so a dense matrix (e.g. distance-decay weights) gives an almost unstructured spatial effect.",
      "i" = "Use a neighbourhood graph: {.code create_weights(method = \"queen\")} for polygons or {.code create_weights(coords, method = \"knn\")} for centroids."
    ))
  }

  # Check if disconnected components or isolated nodes exist
  degree <- rowSums(adj_mat)
  if (any(degree == 0)) {
    isolated_idx <- which(degree == 0)
    cli::cli_warn(
      "{length(isolated_idx)} domain(s) have no spatial neighbors in {.arg W} (isolated nodes): {paste(isolated_idx, collapse = ', ')}."
    )
  }

  # SAR: the weights of eblup_sfh(), which expects a row-standardised W
  if (spatial == "sar") {
    rs <- rowSums(W_mat)
    if (any(rs > 0 & abs(rs - 1) > 1e-8)) {
      cli::cli_inform("{.code spatial = \"sar\"}: {.arg W} is row-standardised.")
      W_mat <- W_mat / ifelse(rs > 0, rs, 1)
    }
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
  call = NULL,
  beta_phi = NULL
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
    hyper_names <- rownames(hyper_summary)
    # Interaction terms on copied indices (..domain_id2.., ..time_id2..) are not main effects
    main_rows <- function(pattern) {
      rows <- grep(pattern, hyper_names, value = TRUE)
      rows[!grepl("_id2\\.\\.", rows)]
    }

    # Spatial / Area random effect precision
    prec_rows <- main_rows("Precision for (\\.\\.domain_id\\.\\.|domain)")
    # bym: report the structured (spatial) component, not the iid one
    spat_comp <- grep("spatial component", prec_rows, value = TRUE)
    if (length(spat_comp) > 0) prec_rows <- spat_comp
    if (length(prec_rows) == 0 && spatial == "none" && temporal == "none") {
      # Fallback for simple IID; never the precision of the observations
      prec_rows <- main_rows("^Precision for")
      prec_rows <- prec_rows[!grepl("observations", prec_rows)]
    }
    if (length(prec_rows) > 0) {
      prec_est <- hyper_summary[prec_rows[1], "mean"]
      if (!is.na(prec_est) && prec_est > 0) {
        random_effect_var <- .post_mean_variance(fit, prec_rows[1], prec_est)
      }
    }

    # Temporal precision
    prec_time_rows <- main_rows("Precision for (\\.\\.time_id\\.\\.|time)")
    if (length(prec_time_rows) > 0) {
      prec_time_est <- hyper_summary[prec_time_rows[1], "mean"]
      if (!is.na(prec_time_est) && prec_time_est > 0) {
        random_effect_var_time <- .post_mean_variance(fit, prec_time_rows[1], prec_time_est)
      }
    }

    # Look for Phi (spatial proportion in BYM2), Beta (generic1 / Leroux), or Rho (sar, slm)
    phi_rows <- main_rows("Phi for")
    beta_rows <- main_rows("^(Beta|Lambda) for")
    rho_rows <- main_rows("^Rho for (\\.\\.domain_id\\.\\.|domain)")
    rho_time_rows <- main_rows("^Rho for (\\.\\.time_id\\.\\.|time)")
    # separable: the AR(1) of the spatial field over time is a group correlation
    if (length(rho_time_rows) == 0) rho_time_rows <- main_rows("^GroupRho for")

    if (length(phi_rows) > 0) {
      phi <- hyper_summary[phi_rows[1], "mean"]
      rho <- phi
    } else if (length(beta_rows) > 0) {
      # In generic1 / Leroux CAR, Beta represents spatial autocorrelation parameter rho
      phi <- hyper_summary[beta_rows[1], "mean"]
      rho <- phi
    } else if (length(rho_rows) > 0) {
      # In SAR and SLM, Rho is the spatial autocorrelation, reported by INLA on (0, 1)
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

  hb_est <- fitted_summary$mean
  hb_sd <- fitted_summary$sd
  hb_mse <- hb_sd^2
  hb_rse <- ifelse(abs(hb_est) < .Machine$double.eps, NA_real_, (hb_sd / abs(hb_est)) * 100)
  ci_lower <- fitted_summary[["0.025quant"]]
  ci_upper <- fitted_summary[["0.975quant"]]
  linpred <- linpred_summary$mean

  # Assemble df_hb
  if (!is.null(time)) {
    df_hb <- data.frame(
      domain = domain,
      time = time,
      y = y,
      hb = hb_est,
      linear_pred = linpred,
      sd = hb_sd,
      mse = hb_mse,
      rse = hb_rse,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      random_effect = rand_eff,
      stringsAsFactors = FALSE
    )
    if (!is.null(rand_eff_spatial)) df_hb$random_effect_spatial <- rand_eff_spatial
    if (!is.null(rand_eff_temporal)) df_hb$random_effect_temporal <- rand_eff_temporal
  } else {
    df_hb <- data.frame(
      domain = domain,
      y = y,
      hb = hb_est,
      linear_pred = linpred,
      sd = hb_sd,
      mse = hb_mse,
      rse = hb_rse,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      random_effect = rand_eff,
      stringsAsFactors = FALSE
    )
  }

  df_hb <- .hb_aux_columns(df_hb, family, y, vardir, trials, exposure, beta_phi = beta_phi)

  # 5. Goodness of fit measures
  goodness <- c(
    DIC = if (!is.null(fit$dic$dic)) fit$dic$dic else NA_real_,
    pD = if (!is.null(fit$dic$p.eff)) fit$dic$p.eff else NA_real_,
    WAIC = if (!is.null(fit$waic$waic)) fit$waic$waic else NA_real_,
    pWAIC = if (!is.null(fit$waic$p.eff)) fit$waic$p.eff else NA_real_,
    Marginal_LogLik = if (!is.null(fit$mlik)) fit$mlik[1, 1] else NA_real_
  )

  model_label <- .hb_model_label(family, spatial, temporal, st_interaction)

  out <- list(
    df_hb = df_hb,
    hb = df_hb,
    df_eblup = df_hb, # compatibility with eblup methods
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
    convergence = .inla_convergence(fit),
    inla_mode_status = fit$mode$mode.status,
    fit = fit,
    call = call
  )

  class(out) <- c("fastsae_hb_area", "fastsae_hb", "fastsae")
  return(out)
}


#' Posterior (Bayesian) self-benchmarking of an INLA area-level fit
#'
#' Draws from the joint posterior of the linear predictor, maps them to the
#' response scale, adjusts every draw so that the benchmark constraint holds
#' exactly (see `.benchmark_adjust()`), and summarises the adjusted draws.
#' Because the constraint holds draw by draw it also holds for the posterior
#' mean, and the reported SD and credible interval account for it (Bayesian
#' benchmarking, cf. Datta et al., 2011).
#'
#' @param res Result of `.extract_inla_results()` (needs `res$fit`, `res$df_hb`).
#' @param sb_info List from `.setup_self_benchmark()`.
#' @param family Response family.
#' @param n_draws Number of posterior draws.
#' @param seed Seed passed to `INLA::inla.posterior.sample()` (0 = random).
#' @param link Link function (`NULL`: canonical link of `family`).
#' @return `res` with benchmarked `df_hb` columns; the unbenchmarked estimate
#'   is kept in `hb_unbenchmarked`.
#' @noRd
.posterior_benchmark <- function(res, sb_info, family, n_draws = 1000, seed = 0L, link = NULL) {
  n <- nrow(res$df_hb)
  smp <- .inla_draws(res$fit, n_draws, seed)
  rn <- rownames(smp[[1]]$latent)
  pred_rows <- grep("^Predictor:", rn)[seq_len(n)]
  eta <- vapply(smp, function(s) s$latent[pred_rows, 1], numeric(n))
  theta <- .inla_inv_link(eta, family, link)
  theta_bm <- .benchmark_adjust(theta, sb_info, .benchmark_type_for_family(family))
  df <- res$df_hb
  df$hb_unbenchmarked <- df$hb
  df$hb <- rowMeans(theta_bm)
  df$sd <- apply(theta_bm, 1, stats::sd)
  df$mse <- df$sd^2
  df$rse <- ifelse(abs(df$hb) < .Machine$double.eps, NA_real_, df$sd / abs(df$hb) * 100)
  df$ci_lower <- apply(theta_bm, 1, stats::quantile, probs = 0.025, names = FALSE)
  df$ci_upper <- apply(theta_bm, 1, stats::quantile, probs = 0.975, names = FALSE)
  res$df_hb <- df
  res$hb <- df
  if (!is.null(res$df_eblup)) res$df_eblup <- df
  res$benchmark_draws <- n_draws
  res
}

#' Model description shared by the hb_area() engines
#' @noRd
.hb_model_label <- function(family, spatial, temporal, st_interaction) {
  model_label <- paste0("HB-", toupper(family))
  if (temporal != "none" || spatial != "none") {
    comps <- c()
    if (spatial != "none") comps <- c(comps, toupper(spatial))
    if (temporal != "none") comps <- c(comps, toupper(temporal))
    if (st_interaction != "none") comps <- c(comps, paste0("ST:", toupper(st_interaction)))
    paste0(model_label, " (", paste(comps, collapse = " + "), ")")
  } else {
    paste0(model_label, " (Non-spatial)")
  }
}

#' Auxiliary df_hb columns shared by the hb_area() engines
#'
#' Sampling variances / precisions, binomial totals and count rates, computed
#' from the HB estimate already in `df_hb$hb`. `beta_phi` is the beta precision
#' the model was fitted with (second stage of .beta_phi_from_vardir()).
#' @noRd
.hb_aux_columns <- function(df_hb, family, y, vardir = NULL, trials = NULL, exposure = NULL,
                            beta_phi = NULL) {
  hb_est <- df_hb$hb
  # If vardir provided (Gaussian FH, Beta, or Gamma)
  if (!is.null(vardir)) {
    df_hb$vardir <- vardir
    if (family == "beta") {
      df_hb$precision <- if (!is.null(beta_phi)) beta_phi else pmax(y * (1 - y) / vardir - 1, 1, na.rm = TRUE)
    } else if (family == "gamma") {
      s_dir <- (y^2) / vardir
      df_hb$precision <- s_dir
      df_hb$cv_dir <- sqrt(vardir) / y
    }
  }

  # If trials provided (Binomial or Beta)
  if (!is.null(trials)) {
    df_hb$trials <- trials
    if (family == "binomial") {
      df_hb$estimated_total <- hb_est * trials
    } else if (family == "beta" && is.null(vardir)) {
      df_hb$precision <- pmax(trials - 1, 1, na.rm = TRUE)
    }
  }

  # If exposure provided (Poisson / NegBinom)
  if (!is.null(exposure)) {
    df_hb$exposure <- exposure
    df_hb$rate <- hb_est
    df_hb$estimated_count <- hb_est * exposure
  }
  df_hb
}

#' Link functions supported per response family (INLA names)
#'
#' The first entry is the canonical link used when `link = NULL`.
#' @noRd
.inla_links <- function(family) {
  switch(family,
    gaussian = "identity",
    binomial = ,
    beta = c("logit", "probit", "cloglog", "loglog"),
    "log" # poisson, nbinomial, gamma
  )
}

#' Validate `inla_engine(link = )` against the response family
#' @return The link to use (canonical when `link` is `NULL`).
#' @noRd
.check_inla_link <- function(link, family) {
  ok <- .inla_links(family)
  if (is.null(link)) {
    return(ok[1])
  }
  link <- tolower(link)
  if (!link %in% ok) {
    cli::cli_abort(c(
      "{.code link = \"{link}\"} is not available for {.code family = \"{family}\"}.",
      "i" = "Use {.or {.val {ok}}}."
    ))
  }
  link
}

#' Inverse link: linear predictor to the mean scale
#' @noRd
.inla_inv_link <- function(eta, family, link = NULL) {
  link <- link %||% .inla_links(family)[1]
  switch(link,
    identity = eta,
    logit = stats::plogis(eta),
    probit = stats::pnorm(eta),
    cloglog = 1 - exp(-exp(eta)),
    loglog = exp(-exp(-eta)),
    log = exp(eta)
  )
}

#' Convergence of an INLA fit
#'
#' INLA returns a fit, not an error, when the search for the posterior mode of
#' the hyperparameters fails: it then reports \code{ok = FALSE} or a non-zero
#' \code{mode$mode.status}. A failure is reported as a warning.
#' @noRd
.inla_convergence <- function(fit) {
  status <- fit$mode$mode.status
  ok <- !isFALSE(fit$ok) && (is.null(status) || isTRUE(status == 0))
  if (!ok) {
    cli::cli_warn(c(
      "INLA did not locate the posterior mode of the hyperparameters reliably{if (!is.null(status)) paste0(' (mode status ', status, ')') else ''}.",
      "i" = "The estimates are returned but are not reliable; try a more informative prior, rescaled covariates or {.code inla_engine(strategy = \"laplace\")}."
    ))
  }
  ok
}
