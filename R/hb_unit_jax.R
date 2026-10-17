#' Fit hb_unit() by MCMC (NumPyro NUTS on JAX)
#'
#' Receives the data prepared by `hb_unit()`: the sampled units followed by one
#' prediction row per domain (response NA) at the covariate mean of the
#' non-sampled units (`popsize_var` given) or of the population. The backend fits
#' the nested error model of Battese, Harter and Fuller (1988), y_dj =
#' x_dj' beta + u_d + e_dj, or its logistic / Poisson version, and returns the
#' draws of the prediction rows only. Returns the same
#' `c("fastsae_hb_unit", "fastsae_hb", "fastsae")` structure as the INLA engine.
#' @noRd
.fit_jax_unit <- function(
  formula,
  term_labels,
  combined_df,
  n_sample,
  y_s,
  dom_s,
  all_domains,
  Xpop_domain,
  popsize_var,
  family,
  spatial,
  adj_graph,
  unsampled_domains,
  engine,
  call
) {
  n_domains <- length(all_domains)
  rhs <- stats::reformulate(term_labels)
  X_mat <- stats::model.matrix(rhs, stats::model.frame(rhs, combined_df, na.action = stats::na.pass))
  if (nrow(X_mat) != nrow(combined_df) || anyNA(X_mat)) {
    cli::cli_abort("The covariates of the sample and of {.arg Xpop} must be complete.")
  }
  pred_rows <- n_sample + seq_len(n_domains)

  spec <- NULL
  if (spatial != "none") {
    if (!is.matrix(adj_graph)) {
      cli::cli_abort("{.code engine = \"jax\"} needs {.arg W} as a matrix, {.cls Matrix}, {.cls nb} or {.cls listw} object, not an INLA graph.")
    }
    spec <- .icar_spectral(adj_graph)
  }

  run <- .run_numpyro(
    engine,
    X = X_mat,
    y = as.numeric(combined_df$y),
    domain_idx = combined_df$..dom_id.. - 1L,
    trials = if (family == "binomial") rep(1, nrow(combined_df)) else NULL,
    D = n_domains,
    T = 1L,
    family = family,
    spatial = spatial,
    temporal = "none",
    st_interaction = "none",
    W_adj = if (!is.null(spec)) adj_graph else NULL,
    scale_factor = if (!is.null(spec)) spec$scale_factor else 1,
    eig_values = spec$eig_values,
    eig_vectors = spec$eig_vectors,
    draw_rows = pred_rows - 1L
  )
  fit_py <- run$fit
  hyper <- fit_py$hyperparameters
  sigma2_e <- hyper$sigma2_e %||% NA_real_
  # Draws of the domain mean of the model, g^{-1}(xbar' beta + u_d), D x S
  mu_draws <- run$draws

  dom_factor <- factor(dom_s, levels = all_domains)
  samp_size <- as.integer(table(dom_factor))
  sample_mean_y <- as.numeric(tapply(y_s, dom_factor, mean, na.rm = TRUE))

  # With a finite population the domain mean is Ybar_d = f_d ybar_{d,s} +
  # (1 - f_d) Ybar_{d,r}. For the Gaussian model Ybar_{d,r} = mu_{d,r} + ebar_{d,r}
  # with ebar_{d,r} ~ N(0, sigma_e^2 / (N_d - n_d)) given the parameters, drawn
  # here with each draw of sigma_e^2 (the INLA engine adds E[sigma_e^2] / (N_d - n_d)
  # to the variance). As there, binomial and Poisson leave out ebar_{d,r}.
  theta <- mu_draws
  pop_size <- rep(NA_real_, n_domains)
  if (!is.null(popsize_var)) {
    pop_size <- as.numeric(Xpop_domain[[popsize_var]])
    f_d <- samp_size / pop_size
    f_d[is.na(f_d) | f_d < 0] <- 0
    f_d[f_d > 1] <- 1
    s2e_draws <- as.numeric(fit_py$sigma2_e_draws)
    if (family == "gaussian" && length(s2e_draws) == ncol(theta)) {
      n_r <- pop_size - samp_size
      inv_r <- ifelse(!is.na(n_r) & n_r > 0, 1 / n_r, 0)
      z <- .with_seed(engine$seed, matrix(stats::rnorm(length(theta)), nrow(theta)))
      theta <- theta + z * sqrt(outer(inv_r, s2e_draws))
    }
    ybar <- ifelse(samp_size > 0, sample_mean_y, 0)
    theta <- f_d * ybar + (1 - f_d) * theta
  }

  hb <- rowMeans(theta)
  hb_sd <- apply(theta, 1, stats::sd)
  q <- apply(theta, 1, stats::quantile, probs = c(0.025, 0.975), names = FALSE)

  df_hb <- data.frame(
    domain = all_domains,
    y = sample_mean_y,
    hb = hb,
    linear_pred = as.numeric(fit_py$linpred_mean)[pred_rows],
    vardir = .hb_unit_vardir(family, sigma2_e, samp_size, sample_mean_y),
    sd = hb_sd,
    mse = hb_sd^2,
    rse = ifelse(abs(hb) < .Machine$double.eps, NA_real_, (hb_sd / abs(hb)) * 100),
    ci_lower = q[1, ],
    ci_upper = q[2, ],
    samp_size = samp_size,
    sample_mean = sample_mean_y,
    random_effect = as.numeric(fit_py$rand_eff_mean)[pred_rows],
    stringsAsFactors = FALSE
  )
  if (!is.null(popsize_var)) {
    df_hb$pop_size <- pop_size
    df_hb$estimated_total <- hb * pop_size
  }

  hyperpar <- if (length(hyper) > 0) {
    data.frame(Parameter = names(hyper), Estimate = as.numeric(unlist(hyper)), row.names = NULL)
  } else {
    data.frame(Parameter = character(0), Estimate = numeric(0))
  }

  rhat_max <- run$rhat_max
  result <- list(
    df_hb = df_hb,
    estcoef = run$estcoef,
    hyperpar = hyperpar,
    random_effect_var = hyper$sigma2_u %||% NA_real_,
    residual_var = sigma2_e,
    phi = hyper$phi %||% NA_real_,
    goodness = list(
      dic = as.numeric(fit_py$dic),
      p_eff_dic = as.numeric(fit_py$p_dic),
      waic = as.numeric(fit_py$waic),
      p_eff_waic = as.numeric(fit_py$p_waic),
      mlik = NA_real_
    ),
    family = family,
    spatial = spatial,
    unsampled_domains = unsampled_domains,
    model = paste0(
      "HB-UNIT-", toupper(family),
      if (spatial != "none") paste0(" (", toupper(spatial), ")") else " (Non-spatial)"
    ),
    method = paste0("NumPyro NUTS (", engine$chains, " x ", engine$samples, " draws, ", fit_py$device_used, ")"),
    convergence = if (!is.null(rhat_max)) rhat_max < 1.05 else NULL,
    rhat_max = rhat_max,
    device = fit_py$device_used,
    elapsed_seconds = run$elapsed,
    level = "unit",
    formula = formula,
    fit = fit_py,
    engine = "jax",
    engine_spec = engine,
    call = call
  )
  if (isTRUE(engine$keep_draws)) result$draws <- theta

  class(result) <- c("fastsae_hb_unit", "fastsae_hb", "fastsae")
  result
}
