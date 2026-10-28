#' Fit hb_twofold() by MCMC (NumPyro NUTS on JAX)
#'
#' Receives the inputs checked by `hb_twofold()` and returns the same
#' `c("fastsae_hb_twofold", "fastsae_hb", "fastsae")` structure as the INLA
#' engine. The backend's nested model has the area effect `u` (iid, Besag or
#' BYM2) and an iid sub-area effect, i.e. v_i and u_ij of Torabi and Rao (2014);
#' the area aggregates come from the MCMC draws through `.twofold_area_table()`.
#' @noRd
.fit_jax_twofold <- function(
  formula,
  mf,
  y,
  domain_vec,
  subarea_vec,
  unique_domains,
  family,
  spatial,
  adj_graph,
  vardir,
  trials,
  exposure,
  weight,
  compute_area,
  engine,
  call
) {
  X_mat <- stats::model.matrix(attr(mf, "terms"), mf)
  y_fit <- as.numeric(y)
  n <- length(y_fit)
  obs <- !is.na(y_fit)

  # binomial cases may be non-integer (effective sample sizes); hb_area() and
  # hb_twofold() have checked them against trials
  if (family == "poisson" && any(y_fit[obs] < 0 | y_fit[obs] %% 1 != 0)) {
    cli::cli_abort("For {.code family = \"{family}\"}, the response must be non-negative integer counts.")
  }
  trials_num <- if (!is.null(trials)) as.numeric(trials) else NULL
  if (family == "binomial" && any(!is.finite(trials_num[obs]) | trials_num[obs] <= 0)) {
    cli::cli_abort("{.arg trials} must be finite and strictly positive for sampled subareas.")
  }
  exposure_num <- if (!is.null(exposure)) as.numeric(exposure) else NULL
  if (!is.null(exposure_num) && any(!is.finite(exposure_num[obs]) | exposure_num[obs] <= 0)) {
    cli::cli_abort("{.arg exposure} must be finite and strictly positive for sampled subareas.")
  }

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
    y = y_fit,
    domain_idx = match(domain_vec, unique_domains) - 1L,
    subarea_idx = seq_len(n) - 1L,
    num_subareas = n,
    vardir = if (family == "gaussian") as.numeric(vardir) else NULL,
    trials = trials_num,
    exposure = exposure_num,
    D = length(unique_domains),
    T = 1L,
    family = family,
    spatial = spatial,
    temporal = "none",
    st_interaction = "none",
    W_adj = if (!is.null(spec)) adj_graph else NULL,
    scale_factor = if (!is.null(spec)) spec$scale_factor else 1,
    eig_values = spec$eig_values,
    eig_vectors = spec$eig_vectors
  )
  fit_py <- run$fit
  draws <- run$draws

  hb_pred <- as.numeric(fit_py$hb_mean)
  hb_sd <- as.numeric(fit_py$hb_sd)
  or_na <- function(x) if (length(x) == n) as.numeric(x) else rep(NA_real_, n)
  df_hb <- data.frame(
    domain = domain_vec,
    subarea = subarea_vec,
    y = y,
    hb = hb_pred,
    linear_pred = as.numeric(fit_py$linpred_mean),
    vardir = if (!is.null(vardir)) vardir else rep(NA_real_, n),
    sd = hb_sd,
    mse = hb_sd^2,
    rse = ifelse(abs(hb_pred) < .Machine$double.eps, NA_real_, (hb_sd / abs(hb_pred)) * 100),
    ci_lower = as.numeric(fit_py$hb_ci_lower),
    ci_upper = as.numeric(fit_py$hb_ci_upper),
    random_effect_area = or_na(fit_py$rand_eff_area_mean),
    random_effect_subarea = or_na(fit_py$rand_eff_subarea_mean),
    stringsAsFactors = FALSE
  )

  df_area <- NULL
  if (compute_area) {
    norm_weight <- .twofold_weights(weight, domain_vec, unique_domains, subarea_vec)
    df_area <- .twofold_area_table(draws, norm_weight, domain_vec, unique_domains, hb_pred, hb_sd)
  }

  hyper <- fit_py$hyperparameters
  hyperpar <- if (length(hyper) > 0) {
    data.frame(Parameter = names(hyper), Estimate = as.numeric(unlist(hyper)), row.names = NULL)
  } else {
    data.frame(Parameter = character(0), Estimate = numeric(0))
  }
  # The backend calls the area variance sigma2_u and the sub-area one sigma2_subarea;
  # hb_twofold() reports them as sigma2_v (v_i) and sigma2_u (u_ij)
  random_effect_var <- c(sigma2_v = hyper$sigma2_u %||% NA_real_, sigma2_u = hyper$sigma2_subarea %||% NA_real_)

  rhat_max <- run$rhat_max
  result <- list(
    df_hb = df_hb,
    df_subarea = df_hb,
    df_area = df_area,
    estcoef = run$estcoef,
    hyperpar = hyperpar,
    random_effect_var = random_effect_var,
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
    model = paste0(
      "HB-TWOFOLD-", toupper(family),
      if (spatial != "none") paste0(" (", toupper(spatial), ")") else " (Non-spatial)"
    ),
    method = paste0("NumPyro NUTS (", engine$chains, " x ", engine$samples, " draws, ", fit_py$device_used, ")"),
    convergence = if (!is.null(rhat_max)) rhat_max < 1.05 else NULL,
    rhat_max = rhat_max,
    device = fit_py$device_used,
    elapsed_seconds = run$elapsed,
    level = "subarea",
    formula = formula,
    fit = fit_py,
    engine = "jax",
    engine_spec = engine,
    call = call
  )
  if (isTRUE(engine$keep_draws)) result$draws <- draws
  result <- c(result, .mcmc_fields(run))

  class(result) <- c("fastsae_hb_twofold", "fastsae_hb", "fastsae")
  result
}
