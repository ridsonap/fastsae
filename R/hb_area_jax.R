#' Check that engine = "jax" supports the requested model
#'
#' Runs before any data processing so that an unsupported combination fails fast.
#' @noRd
.check_jax_model <- function(spatial, temporal, family, self_benchmark, benchmark_method) {
  if (spatial == "slm") {
    cli::cli_abort("{.code spatial = \"slm\"} is not available with {.fn jax_engine}.")
  }
  if (temporal == "rw2") {
    cli::cli_abort("{.code temporal = \"rw2\"} is not available with {.fn jax_engine}.")
  }
  if (isTRUE(self_benchmark) && benchmark_method == "wfq" && family != "gaussian") {
    cli::cli_abort(c(
      "{.code benchmark_method = \"wfq\"} with {.fn jax_engine} is available for {.code family = \"gaussian\"} only.",
      "i" = "Use {.code benchmark_method = \"posterior\"}."
    ))
  }
  invisible(NULL)
}

#' ICAR scaling factor and spectral decomposition of a graph Laplacian
#'
#' Scaling factor: geometric mean of the marginal variances of the ICAR field
#' (Riebler et al., 2016). The eigendecomposition is reused by the Python
#' backend for the ICAR and Leroux bases.
#' @noRd
.icar_spectral <- function(adj_mat) {
  L <- diag(rowSums(adj_mat)) - adj_mat
  eig <- eigen(L, symmetric = TRUE)
  pos_idx <- which(eig$values > 1e-6)
  if (length(pos_idx) > 0) {
    V <- eig$vectors[, pos_idx, drop = FALSE]
    q_inv_diag <- rowSums((V^2) * rep(1 / eig$values[pos_idx], each = nrow(V)))
    scale_factor <- exp(mean(log(pmax(q_inv_diag, 1e-8))))
  } else {
    scale_factor <- 1
  }
  list(scale_factor = scale_factor, eig_values = eig$values, eig_vectors = eig$vectors)
}

#' Run the NumPyro backend with the settings of a jax_engine()
#'
#' Sets the JAX platform for the device, calls `fit_numpyro_hb()` with the
#' model arguments in `...` and returns the backend output (`fit`), the n x S
#' draws of the area parameters (`draws`), the coefficient table (`estcoef`),
#' the largest rank-normalized split R-hat (`rhat_max`), the convergence
#' diagnostics (`mcmc`, see `.mcmc_summary()`) and the sampling time (`elapsed`).
#' Shared by `.fit_jax_area()` and `.fit_jax_twofold()`.
#' @noRd
.run_numpyro <- function(engine, X, ...) {
  # JAX reads JAX_PLATFORMS when it initializes, i.e. at the first fit of the session
  old_plat <- Sys.getenv("JAX_PLATFORMS", unset = NA)
  on.exit(
    {
      if (is.na(old_plat) || !nzchar(old_plat)) {
        Sys.unsetenv("JAX_PLATFORMS")
      } else {
        Sys.setenv(JAX_PLATFORMS = old_plat)
      }
    },
    add = TRUE
  )
  if (!nzchar(Sys.getenv("JAX_PLATFORMS"))) {
    if (engine$device == "cpu" || (engine$device == "auto" && Sys.info()["sysname"] == "Darwin")) {
      Sys.setenv(JAX_PLATFORMS = "cpu")
    } else if (engine$device == "cuda") {
      Sys.setenv(JAX_PLATFORMS = "cuda,cpu")
    }
  }

  backend <- .get_numpyro_backend(device = engine$device)
  cli::cli_alert_info("Running NumPyro NUTS ({engine$chains} chain{?s} x {engine$samples} draws, {engine$warmup} warmup)...")
  t0 <- Sys.time()
  fit_py <- backend$fit_numpyro_hb(
    X = X,
    ...,
    num_warmup = engine$warmup,
    num_samples = engine$samples,
    num_chains = engine$chains,
    device = engine$device,
    seed = engine$seed,
    prior_beta = engine$prior_beta,
    priors = engine$priors,
    pc_alpha = .prior_pc_alpha,
    robust = engine$robust,
    return_draws = TRUE,
    progress_bar = engine$progress
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cli::cli_alert_success("Sampling completed in {round(elapsed, 2)} seconds.")
  .warn_device_fallback(engine$device, fit_py$device_used)

  # Draws of the area parameters, n x S
  draws <- t(as.matrix(fit_py$hb_draws))
  fit_py$hb_draws <- NULL

  beta_mean <- as.numeric(fit_py$beta_mean)
  beta_sd <- as.numeric(fit_py$beta_sd)
  zval <- beta_mean / pmax(beta_sd, 1e-8)
  estcoef <- data.frame(
    beta = beta_mean,
    std.error = beta_sd,
    zvalue = zval,
    pvalue = 2 * stats::pnorm(abs(zval), lower.tail = FALSE),
    ci_lower = as.numeric(fit_py$beta_ci_lower),
    ci_upper = as.numeric(fit_py$beta_ci_upper),
    row.names = colnames(X)
  )
  if (engine$prior_beta == "horseshoe" && length(fit_py$shrinkage_weights) == ncol(X)) {
    estcoef$shrinkage_factor <- round(as.numeric(fit_py$shrinkage_weights), 4)
  }

  rhat_max <- if (!is.null(fit_py$rhat_max)) as.numeric(fit_py$rhat_max) else NULL
  mcmc <- .mcmc_summary(fit_py, colnames(X), engine$chains)
  fit_py$diagnostics <- NULL
  fit_py$trace <- NULL
  .warn_mcmc(rhat_max, mcmc, engine$chains)

  list(fit = fit_py, draws = draws, estcoef = estcoef, rhat_max = rhat_max, mcmc = mcmc, elapsed = elapsed)
}

#' Convergence diagnostics of a NumPyro fit
#'
#' The backend's table of rank-normalized split R-hat and bulk/tail ESS
#' (Vehtari et al., 2021) per hyperparameter, coefficient (named after the
#' columns of X) and area parameter (`hb[j]`, row j of df_hb), the per-chain
#' draws of the key parameters as an iterations x chains x parameters array,
#' and the numbers of divergent transitions and of iterations at the maximum
#' tree depth. NULL entries when the backend gives none.
#' @noRd
.mcmc_summary <- function(fit_py, coef_names, chains) {
  tab <- fit_py$diagnostics
  if (is.null(tab)) {
    return(list(
      table = NULL, trace = NULL, ess_min = NULL, divergences = fit_py$num_divergent,
      treedepth_hits = fit_py$treedepth_hits
    ))
  }
  rename <- function(p) {
    k <- suppressWarnings(as.integer(sub("^beta\\[([0-9]+)\\]$", "\\1", p)))
    ifelse(grepl("^beta\\[", p) & !is.na(k) & k <= length(coef_names), coef_names[pmax(k, 1)], p)
  }
  table <- data.frame(
    parameter = rename(as.character(unlist(tab$parameter))),
    mean = as.numeric(tab$mean),
    sd = as.numeric(tab$sd),
    rhat = as.numeric(tab$rhat),
    ess_bulk = as.numeric(tab$ess_bulk),
    ess_tail = as.numeric(tab$ess_tail),
    stringsAsFactors = FALSE
  )
  trace <- NULL
  if (!is.null(fit_py$trace)) {
    tr <- as.array(fit_py$trace$draws)
    trace <- aperm(array(tr, dim = dim(tr)), c(2, 1, 3))
    dimnames(trace) <- list(
      iteration = NULL, chain = seq_len(dim(trace)[2]),
      parameter = rename(as.character(unlist(fit_py$trace$names)))
    )
  }
  ess <- c(table$ess_bulk, table$ess_tail)
  list(
    table = table,
    trace = trace,
    ess_min = if (any(is.finite(ess))) min(ess, na.rm = TRUE) else NA_real_,
    divergences = fit_py$num_divergent,
    treedepth_hits = fit_py$treedepth_hits
  )
}

#' Warn about MCMC convergence (thresholds of Vehtari et al., 2021)
#'
#' R-hat >= 1.05 marks the fit as not converged (as before); 1.01 <= R-hat < 1.05
#' only informs. Bulk or tail ESS below 100 per chain and divergent transitions
#' warn; iterations at the maximum tree depth (an efficiency issue) inform.
#' The warnings have class `fastsae_mcmc_warning`.
#' @noRd
.warn_mcmc <- function(rhat_max, mcmc, chains) {
  tab <- mcmc$table
  worst <- function(col, fun) {
    if (is.null(tab) || !any(is.finite(tab[[col]]))) {
      return(NULL)
    }
    tab$parameter[fun(tab[[col]])]
  }
  if (!is.null(rhat_max) && !(rhat_max < 1.05)) {
    par <- worst("rhat", which.max)
    cli::cli_warn(c(
      "The MCMC chains have not converged (max R-hat = {round(rhat_max, 3)} >= 1.05{if (!is.null(par)) paste0(', ', par)}).",
      "i" = "Increase {.arg warmup} and {.arg samples} in {.fn jax_engine}; see {.code autoplot(fit, type = \"trace\")}."
    ), class = "fastsae_mcmc_warning")
  } else if (!is.null(rhat_max) && rhat_max >= 1.01) {
    cli::cli_alert_warning("Max R-hat = {round(rhat_max, 3)} ({worst('rhat', which.max)}); values below 1.01 are recommended.")
  }
  if (!is.null(mcmc$ess_min) && is.finite(mcmc$ess_min) && mcmc$ess_min < 100 * chains) {
    par <- tab$parameter[which.min(pmin(tab$ess_bulk, tab$ess_tail))]
    cli::cli_warn(c(
      "Low effective sample size: min bulk/tail ESS = {round(mcmc$ess_min)} ({par}), below 100 per chain ({100 * chains}).",
      "i" = "The posterior SDs and intervals are unreliable; increase {.arg samples} in {.fn jax_engine}."
    ), class = "fastsae_mcmc_warning")
  }
  if (!is.null(mcmc$divergences) && mcmc$divergences > 0) {
    cli::cli_warn(c(
      "{mcmc$divergences} divergent transition{?s} after warmup: the sampler may miss parts of the posterior.",
      "i" = "Increase {.arg warmup}, or simplify the model (e.g. a weakly identified variance)."
    ), class = "fastsae_mcmc_warning")
  }
  if (!is.null(mcmc$treedepth_hits) && mcmc$treedepth_hits > 0) {
    cli::cli_alert_warning("{mcmc$treedepth_hits} iteration{?s} reached the maximum tree depth (slow sampling, not a validity problem).")
  }
  invisible(NULL)
}

#' Convergence fields of a jax_engine() result
#' @noRd
.mcmc_fields <- function(run) {
  list(
    mcmc_diagnostics = run$mcmc$table,
    ess_min = run$mcmc$ess_min,
    divergences = run$mcmc$divergences,
    treedepth_hits = run$mcmc$treedepth_hits,
    trace = run$mcmc$trace
  )
}

#' Fit hb_area() by MCMC (NumPyro NUTS on JAX)
#'
#' Receives the data prepared by `hb_area()` and returns the same
#' `c("fastsae_hb_area", "fastsae_hb", "fastsae")` structure as the INLA engine.
#' Self-benchmarking adjusts the posterior draws with `.benchmark_adjust()`,
#' exactly as `.posterior_benchmark()` does for INLA.
#' @noRd
.fit_jax_area <- function(
  formula,
  data,
  y,
  domain,
  time = NULL,
  family,
  spatial = "none",
  temporal = "none",
  st_interaction = "none",
  W_obj = NULL,
  vardir = NULL,
  trials = NULL,
  exposure = NULL,
  engine = jax_engine(),
  sb_info = NULL,
  benchmark_weight = NULL,
  benchmark_target = NULL,
  benchmark_group = NULL,
  call = NULL
) {
  ctrl <- engine
  prior_beta <- engine$prior_beta
  robust <- engine$robust

  # Response and design (after WFQ augmentation, when used)
  mf <- stats::model.frame(formula, data, na.action = stats::na.pass)
  y_fit <- as.numeric(stats::model.response(mf))
  X_mat <- stats::model.matrix(attr(mf, "terms"), mf)
  if (anyNA(X_mat) || any(!is.finite(X_mat))) {
    cli::cli_abort("Covariates must not contain missing or non-finite values with {.code engine = \"jax\"}.")
  }
  obs <- !is.na(y_fit)

  vardir_num <- if (!is.null(vardir)) as.numeric(vardir) else NULL
  if (!is.null(vardir_num) && family %in% c("gaussian", "beta")) {
    bad <- obs & !(is.finite(vardir_num) & vardir_num > 0)
    if (any(bad)) {
      cli::cli_abort("{.arg vardir} must be finite and strictly positive for sampled domains ({sum(bad)} invalid).")
    }
  }
  if (family == "gamma" && !is.null(vardir_num)) {
    cli::cli_alert_info("{.code engine = \"jax\"} estimates a common Gamma shape; {.arg vardir} is not used for {.code family = \"gamma\"}.")
    vardir_num <- NULL
  }
  if (family %in% c("binomial", "poisson", "nbinomial") && any(y_fit[obs] < 0 | y_fit[obs] %% 1 != 0)) {
    cli::cli_abort("For {.code family = \"{family}\"}, the response must be non-negative integer counts.")
  }
  if (family == "beta" && any(y_fit[obs] <= 0 | y_fit[obs] >= 1)) {
    cli::cli_abort("For {.code family = \"beta\"}, the response must lie strictly in (0, 1).")
  }
  trials_num <- if (!is.null(trials)) as.numeric(trials) else NULL
  if (family == "binomial" && any(!is.finite(trials_num[obs]) | trials_num[obs] <= 0)) {
    cli::cli_abort("{.arg trials} must be finite and strictly positive for sampled domains.")
  }
  exposure_num <- if (!is.null(exposure)) as.numeric(exposure) else NULL
  if (!is.null(exposure_num) && any(!is.finite(exposure_num[obs]) | exposure_num[obs] <= 0)) {
    cli::cli_abort("{.arg exposure} must be finite and strictly positive for sampled domains.")
  }

  unique_domains <- unique(domain)
  domain_idx <- as.integer(factor(domain, levels = unique_domains)) - 1L
  unique_times <- if (!is.null(time)) sort(unique(time)) else NULL
  time_idx <- if (!is.null(time)) as.integer(factor(time, levels = unique_times)) - 1L else NULL

  spec <- NULL
  if (spatial != "none") {
    if (!is.matrix(W_obj$graph)) {
      cli::cli_abort("{.code engine = \"jax\"} needs {.arg W} as a matrix, {.cls Matrix}, {.cls nb} or {.cls listw} object, not an INLA graph.")
    }
    if (spatial != "sar") spec <- .icar_spectral(W_obj$graph)
  }

  run <- .run_numpyro(
    engine,
    X = X_mat,
    y = y_fit,
    domain_idx = domain_idx,
    time_idx = time_idx,
    vardir = vardir_num,
    trials = trials_num,
    exposure = exposure_num,
    D = length(unique_domains),
    T = if (!is.null(unique_times)) length(unique_times) else 1L,
    family = family,
    spatial = if (spatial == "generic1") "leroux" else spatial,
    temporal = temporal,
    st_interaction = st_interaction,
    W_adj = if (!is.null(spec)) W_obj$graph else NULL,
    scale_factor = if (!is.null(spec)) spec$scale_factor else 1,
    eig_values = spec$eig_values,
    eig_vectors = spec$eig_vectors,
    W_sar = if (spatial == "sar") W_obj$W_mat else NULL,
    benchmark = FALSE
  )
  fit_py <- run$fit
  draws <- run$draws
  estcoef <- run$estcoef
  elapsed <- run$elapsed

  hb_est <- as.numeric(fit_py$hb_mean)
  hb_sd <- as.numeric(fit_py$hb_sd)
  df_hb <- data.frame(domain = domain, stringsAsFactors = FALSE)
  if (!is.null(time)) df_hb$time <- time
  df_hb$y <- y
  df_hb$hb <- hb_est
  df_hb$linear_pred <- as.numeric(fit_py$linpred_mean)
  df_hb$sd <- hb_sd
  df_hb$mse <- hb_sd^2
  df_hb$rse <- ifelse(abs(hb_est) < .Machine$double.eps, NA_real_, (hb_sd / abs(hb_est)) * 100)
  df_hb$ci_lower <- as.numeric(fit_py$hb_ci_lower)
  df_hb$ci_upper <- as.numeric(fit_py$hb_ci_upper)
  df_hb$random_effect <- as.numeric(fit_py$rand_eff_mean)
  df_hb <- .hb_aux_columns(df_hb, family, y, vardir, trials, exposure)

  hyper <- fit_py$hyperparameters
  hyperpar <- if (length(hyper) > 0) {
    data.frame(Parameter = names(hyper), Estimate = as.numeric(unlist(hyper)), row.names = NULL)
  } else {
    data.frame(Parameter = character(0), Estimate = numeric(0))
  }

  model_label <- .hb_model_label(family, spatial, temporal, st_interaction)
  extras <- c(if (prior_beta == "horseshoe") "Horseshoe", if (isTRUE(robust)) "Student-t")
  if (length(extras) > 0) model_label <- paste0(model_label, " [", paste(extras, collapse = " + "), "]")

  rhat_max <- run$rhat_max

  out <- list(
    df_hb = df_hb,
    hb = df_hb,
    df_eblup = df_hb,
    estcoef = estcoef,
    hyperpar = hyperpar,
    random_effect_var = hyper$sigma2_u,
    random_effect_var_time = hyper$sigma2_t,
    phi = hyper$phi,
    rho = hyper$rho_spatial,
    rho_time = hyper$rho_t,
    goodness = c(
      DIC = as.numeric(fit_py$dic), pD = as.numeric(fit_py$p_dic),
      WAIC = as.numeric(fit_py$waic), pWAIC = as.numeric(fit_py$p_waic)
    ),
    family = family,
    spatial = spatial,
    temporal = temporal,
    st_interaction = st_interaction,
    level = "area",
    model = model_label,
    convergence = if (!is.null(rhat_max)) rhat_max < 1.05 else NULL,
    rhat_max = rhat_max,
    method = paste0("NumPyro NUTS (", ctrl$chains, " x ", ctrl$samples, " draws, ", fit_py$device_used, ")"),
    device = fit_py$device_used,
    elapsed_seconds = elapsed,
    prior_beta = prior_beta,
    robust = robust,
    engine_spec = engine,
    fit = fit_py,
    call = call
  )

  if (!is.null(sb_info)) {
    if (identical(sb_info$method, "posterior")) {
      draws <- .benchmark_adjust(draws, sb_info, .benchmark_type_for_family(family))
      df <- out$df_hb
      df$hb_unbenchmarked <- df$hb
      df$hb <- rowMeans(draws)
      df$sd <- apply(draws, 1, stats::sd)
      df$mse <- df$sd^2
      df$rse <- ifelse(abs(df$hb) < .Machine$double.eps, NA_real_, df$sd / abs(df$hb) * 100)
      df$ci_lower <- apply(draws, 1, stats::quantile, probs = 0.025, names = FALSE)
      df$ci_upper <- apply(draws, 1, stats::quantile, probs = 0.975, names = FALSE)
      if (!is.null(df$estimated_total)) df$estimated_total <- df$hb * trials
      if (!is.null(df$rate)) df$rate <- df$hb
      if (!is.null(df$estimated_count)) df$estimated_count <- df$hb * exposure
      out$df_hb <- out$hb <- out$df_eblup <- df
    }
    out$self_benchmark <- TRUE
    out$benchmark_method <- sb_info$method
    out$benchmark_weight <- benchmark_weight
    out$benchmark_target <- benchmark_target
    out$benchmark_group <- benchmark_group
    out$benchmark_summary <- .compute_benchmark_summary(
      unique_groups = sb_info$unique_groups,
      group_vec = sb_info$group_vec,
      w_norm = sb_info$w_norm,
      y = sb_info$original_y,
      target_vec = sb_info$target_vec,
      estimates = out$df_hb$hb
    )
    if (identical(sb_info$method, "wfq")) {
      out$df_hb$y <- out$hb$y <- out$df_eblup$y <- sb_info$original_y
    }
  }
  if (isTRUE(ctrl$keep_draws)) out$draws <- draws
  out <- c(out, .mcmc_fields(run))

  if (inherits(data, "sf") && requireNamespace("sf", quietly = TRUE) && nrow(out$df_hb) == nrow(data)) {
    out$df_hb <- sf::st_set_geometry(out$df_hb, sf::st_geometry(data))
  }

  class(out) <- c("fastsae_hb_area", "fastsae_hb", "fastsae")
  out
}
