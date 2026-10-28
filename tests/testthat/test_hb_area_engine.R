# hb_area(engine = ...): engine selection, validation, and the NumPyro/JAX engine.
# The R side of the JAX engine is tested against a fake Python backend, so these
# tests run without Python; the end-to-end fits at the bottom need jax + numpyro.

sim_engine_df <- function(D = 12, seed = 1) {
  set.seed(seed)
  x <- rnorm(D)
  data.frame(
    area = paste0("a", seq_len(D)),
    x = x,
    y = 1 + 0.5 * x + rnorm(D, sd = 0.4),
    vardir = 0.1,
    w = seq_len(D),
    prov = rep(c("P1", "P2", "P3"), length.out = D)
  )
}

chain_W <- function(D) {
  W <- matrix(0, D, D)
  W[cbind(1:(D - 1), 2:D)] <- 1
  W + t(W)
}

# Stand-in for hb_area_numpyro.fit_numpyro_hb(): draws centred on 1 + 0.5 * x
fake_backend <- function(env) {
  list(fit_numpyro_hb = function(y, X, domain_idx, ...) {
    args <- list(...)
    env$args <- c(list(y = y, X = X, domain_idx = domain_idx), args)
    n <- nrow(X)
    S <- 200
    set.seed(3)
    centre <- if (args$family %in% c("binomial", "beta")) stats::plogis(drop(X %*% rep(0.3, ncol(X)))) else 1 + 0.5 * X[, ncol(X)]
    draws <- matrix(rep(centre, each = S), S, n) * exp(matrix(stats::rnorm(S * n, sd = 0.05), S, n))
    list(
      beta_mean = rep(0.5, ncol(X)), beta_sd = rep(0.1, ncol(X)),
      beta_ci_lower = rep(0.3, ncol(X)), beta_ci_upper = rep(0.7, ncol(X)),
      hb_mean = colMeans(draws), hb_sd = apply(draws, 2, stats::sd),
      hb_ci_lower = apply(draws, 2, stats::quantile, 0.025, names = FALSE),
      hb_ci_upper = apply(draws, 2, stats::quantile, 0.975, names = FALSE),
      linpred_mean = colMeans(draws), rand_eff_mean = rep(0, n),
      hyperparameters = list(sigma2_u = 0.2, phi = 0.4),
      dic = 10, p_dic = 3, waic = 11, p_waic = 3.5,
      device_used = "cpu", rhat_max = 1.01, shrinkage_weights = NULL,
      hb_draws = draws
    )
  })
}

local_fake_backend <- function(env, frame = parent.frame()) {
  local_mocked_bindings(.get_numpyro_backend = function(...) fake_backend(env), .env = frame)
}

# print() writes through cli (messages) and cat (stdout)
printed <- function(x) {
  out <- NULL
  msgs <- capture_messages(out <- utils::capture.output(print(x)))
  paste(c(msgs, out), collapse = "\n")
}

skip_if_no_numpyro <- function() {
  skip_on_cran()
  skip_if_not_installed("reticulate")
  if (!suppressMessages(check_numpyro_available())) skip("NumPyro / JAX not available")
}

# ---------------------------------------------------------------------------
# Engine selection and validation (no Python needed)
# ---------------------------------------------------------------------------

test_that("engine constructors validate their settings", {
  expect_error(jax_engine(chains = 0), "whole number")
  expect_error(jax_engine(warmup = 1.5), "whole number")
  expect_error(jax_engine(seed = NA), "whole number")
  expect_error(jax_engine(keep_draws = "yes"), "TRUE")
  expect_error(jax_engine(robust = NA), "TRUE")
  expect_error(jax_engine(device = "tpu"))
  expect_error(jax_engine(prior_beta = "lasso"))
  eng <- jax_engine(samples = 10, device = "CPU", prior_beta = "horseshoe")
  expect_s3_class(eng, c("fastsae_engine_jax", "fastsae_engine"))
  expect_identical(eng$samples, 10L)
  expect_identical(eng$device, "cpu")
  expect_identical(eng$prior_beta, "horseshoe")

  expect_error(inla_engine(strategy = "exact"))
  expect_error(inla_engine(link = 1), "single string")
  expect_error(inla_engine(prior_prec = "pc.prec"), "INLA prior syntax")
  expect_error(inla_engine(draws = -1), "whole number")
  eng <- inla_engine(strategy = "simplified.laplace", control.compute = list(cpo = TRUE))
  expect_s3_class(eng, c("fastsae_engine_inla", "fastsae_engine"))
  expect_identical(eng$strategy, "simplified.laplace")
  expect_identical(eng$inla_args, list(control.compute = list(cpo = TRUE)))

  expect_match(printed(jax_engine()), "JAX engine")
  expect_match(printed(eng), "control.compute")
})

test_that("engine strings map to default engine objects", {
  expect_identical(.as_hb_engine("inla"), inla_engine())
  expect_identical(.as_hb_engine("JAX"), jax_engine())
  eng <- jax_engine(samples = 5)
  expect_identical(.as_hb_engine(eng), eng)
  expect_error(.as_hb_engine("laplace"), "removed")
  expect_error(.as_hb_engine("stan"), "must be")
  expect_error(.as_hb_engine(list(name = "inla")), "must be")
  df <- sim_engine_df()
  expect_error(hb_area(y ~ x, data = df, vardir = "vardir", engine = "laplace", print_result = FALSE), "removed")
})

test_that("jax_engine() rejects unsupported models and INLA arguments before fitting", {
  df <- sim_engine_df()
  fit_jax <- function(...) hb_area(y ~ x, data = df, vardir = "vardir", engine = "jax", print_result = FALSE, ...)
  expect_error(fit_jax(spatial = "slm", W = chain_W(12)), "not available")
  expect_error(fit_jax(time = rep(1:2, 6), domain = rep(1:6, 2), temporal = "rw2"), "not available")
  expect_error(fit_jax(strategy = "gaussian"), "cannot be used with")
  expect_error(fit_jax(control.compute = list(cpo = TRUE)), "cannot be used with")
  df$yc <- rpois(12, 5)
  expect_error(
    hb_area(yc ~ x,
      data = df, family = "poisson", engine = "jax", self_benchmark = TRUE,
      benchmark_method = "wfq", print_result = FALSE
    ),
    "gaussian"
  )
})

test_that("MCMC settings given to hb_area() point to jax_engine()", {
  df <- sim_engine_df()
  for (bad in list(list(robust = TRUE), list(prior_beta = "horseshoe"), list(control = list(samples = 10)))) {
    expect_error(
      do.call(hb_area, c(list(y ~ x, data = df, vardir = "vardir", print_result = FALSE), bad)),
      "jax_engine"
    )
  }
})

test_that("deprecated method and INLA arguments of hb_area() map onto inla_engine()", {
  expect_warning(eng <- .legacy_engine_args("inla", TRUE, "laplace", list()), "deprecated")
  expect_identical(eng, inla_engine(strategy = "laplace"))
  expect_warning(eng <- .legacy_engine_args("inla", TRUE, "inla", list()), "deprecated")
  expect_identical(eng, inla_engine())
  expect_error(.legacy_engine_args("jax", FALSE, "laplace", list()), "deprecated")

  expect_warning(
    eng <- .legacy_engine_args("inla", TRUE, NULL, list(
      strategy = "gaussian", benchmark_draws = 50, benchmark_seed = 3,
      prior_rho = list(prior = "normal", param = c(0, 1)), num.threads = "1:1"
    )),
    "deprecated"
  )
  expect_identical(eng, inla_engine(
    strategy = "gaussian", draws = 50, seed = 3,
    prior_rho = list(prior = "normal", param = c(0, 1)), num.threads = "1:1"
  ))

  # merged into an explicit inla_engine(), unless the same setting is given twice
  expect_warning(
    eng <- .legacy_engine_args(inla_engine(draws = 10), FALSE, NULL, list(strategy = "gaussian")),
    "deprecated"
  )
  expect_identical(eng, inla_engine(draws = 10, strategy = "gaussian"))
  expect_error(
    suppressWarnings(.legacy_engine_args(inla_engine(strategy = "gaussian"), FALSE, NULL, list(strategy = "laplace"))),
    "both"
  )
  expect_error(
    suppressWarnings(.legacy_engine_args(inla_engine(num.threads = 1), FALSE, NULL, list(num.threads = 2))),
    "both"
  )
  expect_identical(.legacy_engine_args(jax_engine(), FALSE, NULL, list()), jax_engine())
})

test_that("hb_unit() and hb_twofold() take inla_engine() with their own defaults", {
  pc <- list(prior = "pc.prec", param = c(1, 0.01))
  unit_defaults <- list(
    prior_prec = list(prior = "loggamma", param = c(0.01, 0.01)),
    prior_phi = list(prior = "pc", param = c(0.5, 0.5))
  )
  eng <- .resolve_inla_engine("inla", TRUE, list(), "hb_unit",
    defaults = unit_defaults,
    used = c("scale_model", "prior_prec", "prior_phi")
  )
  expect_identical(eng$prior_prec, unit_defaults$prior_prec)
  eng <- .resolve_inla_engine(inla_engine(prior_prec = pc), FALSE, list(), "hb_unit",
    defaults = unit_defaults, used = c("scale_model", "prior_prec", "prior_phi")
  )
  expect_identical(eng$prior_prec, pc)

  # settings the function does not use are ignored with a warning
  expect_warning(
    .resolve_inla_engine(inla_engine(prior_rho_time = pc, link = "log"), FALSE, list(), "hb_unit",
      defaults = unit_defaults, used = c("scale_model", "prior_prec", "prior_phi")
    ),
    "not used by"
  )

  # the old hb_twofold() arguments are renamed
  expect_warning(
    eng <- .resolve_inla_engine("inla", TRUE, list(n_samples = 50, prior_prec_area = pc, strategy = "gaussian"),
      "hb_twofold",
      defaults = list(draws = 200L),
      used = c("scale_model", "prior_prec", "prior_prec_subarea", "prior_phi", "draws")
    ),
    "deprecated"
  )
  expect_identical(eng$draws, 50L)
  expect_identical(eng$prior_prec, pc)
  expect_identical(eng$strategy, "gaussian")
  eng <- .resolve_inla_engine("inla", TRUE, list(), "hb_twofold",
    defaults = list(draws = 200L),
    used = c("scale_model", "prior_prec", "prior_prec_subarea", "prior_phi", "draws")
  )
  expect_identical(eng$draws, 200L)

  # hb_unit() and hb_twofold() also take jax_engine(); see test_hb_unit_engine.R and
  # test_hb_twofold_engine.R
  expect_error(hb_unit(y ~ x, data.frame(), data.frame(), "d", engine = "laplace"), "removed")
})

test_that(".icar_spectral matches the Riebler scaling of a chain graph", {
  W <- chain_W(5)
  sp <- .icar_spectral(W)
  Q <- diag(rowSums(W)) - W
  Q_inv <- MASS::ginv(Q)
  expect_equal(sp$scale_factor, exp(mean(log(diag(Q_inv)))), tolerance = 1e-8)
  expect_length(sp$eig_values, 5)
  expect_equal(.icar_spectral(matrix(0, 2, 2))$scale_factor, 1)
})

# ---------------------------------------------------------------------------
# R side of engine = "jax" with a fake backend
# ---------------------------------------------------------------------------

test_that("engine = 'jax' maps arguments to the backend and returns the fastsae structure", {
  env <- new.env()
  local_fake_backend(env)
  df <- sim_engine_df()
  df$y[3] <- NA

  fit <- hb_area(y ~ x,
    data = df, domain = "area", vardir = "vardir", spatial = "leroux",
    W = chain_W(12), engine = jax_engine(samples = 50, device = "cpu", prior_beta = "horseshoe"),
    print_result = FALSE
  )

  a <- env$args
  expect_equal(a$spatial, "leroux")
  expect_equal(a$domain_idx, 0:11)
  expect_true(is.na(a$y[3]))
  expect_equal(a$num_samples, 50L)
  expect_equal(a$device, "cpu")
  expect_equal(a$prior_beta, "horseshoe")
  expect_false(a$robust)
  expect_true(a$return_draws)
  expect_false(a$benchmark)
  expect_equal(a$scale_factor, .icar_spectral(chain_W(12))$scale_factor)
  expect_equal(dim(a$eig_vectors), c(12, 12))

  expect_s3_class(fit, c("fastsae_hb_area", "fastsae_hb", "fastsae"))
  expect_equal(fit$engine, "jax")
  expect_identical(fit$engine_spec, jax_engine(samples = 50, device = "cpu", prior_beta = "horseshoe"))
  expect_equal(fit$spatial, "generic1")
  expect_equal(
    names(fit$df_hb),
    c("domain", "y", "hb", "linear_pred", "sd", "mse", "rse", "ci_lower", "ci_upper", "random_effect", "vardir")
  )
  expect_equal(fit$df_hb$domain, df$area)
  expect_equal(rownames(fit$estcoef), c("(Intercept)", "x"))
  expect_equal(fit$random_effect_var, 0.2)
  expect_equal(fit$phi, 0.4)
  expect_true(fit$convergence)
  expect_equal(dim(fit$draws), c(12, 200))
  expect_named(fit$goodness, c("DIC", "pD", "WAIC", "pWAIC"))

  expect_match(printed(fit), "Engine: jax")
  expect_match(printed(fit), "max R-hat = 1.010")
  expect_match(printed(summary(fit)), "Engine: jax")
  expect_length(coef(fit), 2)
})

test_that("engine = 'jax' posterior self-benchmarking holds per group", {
  env <- new.env()
  local_fake_backend(env)
  df <- sim_engine_df()
  fit <- hb_area(y ~ x,
    data = df, domain = "area", vardir = "vardir", engine = "jax",
    self_benchmark = TRUE, benchmark_weight = "w", benchmark_group = "prov",
    print_result = FALSE
  )
  expect_true(fit$self_benchmark)
  expect_equal(fit$benchmark_method, "posterior")
  expect_true("hb_unbenchmarked" %in% names(fit$df_hb))
  expect_equal(fit$benchmark_summary$discrepancy, rep(0, 3), tolerance = 1e-8)
  # every draw satisfies the constraint
  for (g in unique(df$prov)) {
    idx <- df$prov == g
    w <- df$w[idx] / sum(df$w[idx])
    target <- sum(w * df$y[idx])
    expect_equal(colSums(w * fit$draws[idx, ]), rep(target, ncol(fit$draws)), tolerance = 1e-8)
  }
  expect_match(printed(fit), "benchmarked posterior draws")
})

test_that("engine = 'jax' converts binomial proportions and keeps draws optional", {
  env <- new.env()
  local_fake_backend(env)
  set.seed(1)
  df <- data.frame(id = 1:15, x = rnorm(15), n = 40)
  df$p <- rbinom(15, df$n, 0.3) / df$n
  expect_message(
    fit <- hb_area(p ~ x,
      data = df, domain = "id", family = "binomial", trials = "n",
      engine = jax_engine(keep_draws = FALSE), print_result = FALSE
    ),
    "proportions"
  )
  expect_equal(env$args$y, round(df$p * df$n))
  expect_equal(env$args$family, "binomial")
  expect_null(fit$draws)
  expect_true(all(c("trials", "estimated_total") %in% names(fit$df_hb)))
})

test_that("engine = 'jax' takes non-integer effective numbers of cases", {
  env <- new.env()
  local_fake_backend(env)
  set.seed(2)
  n <- sample(30:120, 15, replace = TRUE)
  df <- data.frame(id = 1:15, x = rnorm(15), n_eff = n / runif(15, 1.2, 2.5))
  df$p <- rbinom(15, n, 0.3) / n
  expect_message(
    hb_area(p ~ x, data = df, domain = "id", family = "binomial", trials = "n_eff", engine = "jax", print_result = FALSE),
    "not rounded"
  )
  expect_equal(env$args$y, df$p * df$n_eff)
  expect_equal(env$args$trials, df$n_eff)
  df$cnt <- round(df$p * df$n_eff)
  df$cnt[1] <- df$n_eff[1] + 2
  expect_error(
    hb_area(cnt ~ x, data = df, domain = "id", family = "binomial", trials = "n_eff", engine = "jax", print_result = FALSE),
    "between 0 and"
  )
})

test_that("benchmark_sae() uses the stored draws of engine = 'jax'", {
  env <- new.env()
  local_fake_backend(env)
  df <- sim_engine_df()
  fit <- hb_area(y ~ x, data = df, domain = "area", vardir = "vardir", engine = "jax", print_result = FALSE)
  M <- .mse_matrix(fit)
  expect_equal(dim(M), c(12, 12))
  expect_equal(diag(M), fit$df_hb$mse, tolerance = 1e-6)
  bm <- benchmark_sae(fit, target = 2, weight = df$w, method = "ratio")
  expect_s3_class(bm, "fastsae_benchmark")
})

test_that("engine = 'jax' validates sampling variances and counts", {
  env <- new.env()
  local_fake_backend(env)
  df <- sim_engine_df()
  df$vardir[2] <- NA
  expect_error(
    hb_area(y ~ x, data = df, vardir = "vardir", engine = "jax", print_result = FALSE),
    "strictly positive"
  )
  df <- sim_engine_df()
  expect_error(
    hb_area(y ~ x, data = df, family = "poisson", engine = "jax", print_result = FALSE),
    "non-negative integer"
  )
})

# ---------------------------------------------------------------------------
# End-to-end fits (need Python with jax and numpyro)
# ---------------------------------------------------------------------------

jax_eng <- jax_engine(warmup = 300, samples = 300, chains = 2, device = "cpu", progress = FALSE)

test_that("the Metal replacement of NUTS' checkpoint index matches NumPyro", {
  skip_if_no_numpyro()
  backend <- fastsae:::.get_numpyro_backend(device = "cpu")
  hmc_util <- reticulate::import("numpyro.infer.hmc_util")
  jnp <- reticulate::import("jax.numpy")
  n <- jnp$arange(0L, 4096L, dtype = jnp$int32)
  ours <- backend$`_leaf_idx_to_ckpt_idxs`(n)
  ref <- hmc_util$`_leaf_idx_to_ckpt_idxs`(n)
  expect_identical(
    as.integer(reticulate::py_to_r(jnp$asarray(ours[[1]])$tolist())),
    as.integer(reticulate::py_to_r(jnp$asarray(ref[[1]])$tolist()))
  )
  expect_identical(
    as.integer(reticulate::py_to_r(jnp$asarray(ours[[2]])$tolist())),
    as.integer(reticulate::py_to_r(jnp$asarray(ref[[2]])$tolist()))
  )
})

test_that("engine = 'jax' Fay-Herriot agrees with INLA", {
  skip_if_no_numpyro()
  skip_if_not_installed("INLA")
  data("mys", package = "fastsae")
  fit_jax <- quiet_mcmc(hb_area(y ~ x1 + x2 + x3,
    data = mys, vardir = "vardir",
    engine = jax_eng, print_result = FALSE
  ))
  fit_inla <- hb_area(y ~ x1 + x2 + x3, data = mys, vardir = "vardir", print_result = FALSE)
  expect_identical(names(fit_jax$df_hb), names(fit_inla$df_hb))
  expect_gt(cor(fit_jax$df_hb$hb, fit_inla$df_hb$hb), 0.99)
  expect_lt(median(abs(fit_jax$df_hb$hb / fit_inla$df_hb$hb - 1)), 0.05)
  expect_true(fit_jax$convergence)
  expect_equal(dim(fit_jax$draws), c(nrow(mys), 600))
})

test_that("engine = 'jax' WAIC of a Gaussian model uses the residual SD of each draw", {
  skip_if_no_numpyro()
  data("mys", package = "fastsae")
  # one observation per area leaves the area and residual variances weakly identified (the
  # chains may not mix), which does not matter for the WAIC identity checked here
  fit <- suppressWarnings(suppressMessages(
    hb_area(y ~ x1 + x2 + x3, data = mys, engine = jax_eng, print_result = FALSE)
  ))
  s2e <- as.numeric(fit$fit$sigma2_e_draws)
  expect_length(s2e, 600)
  expect_equal(fit$residual_var %||% mean(s2e), mean(s2e))
  # lppd and p_WAIC from the pointwise log-likelihood, S x n (Gelman et al., 2014)
  sd_e <- sqrt(s2e)
  obs <- !is.na(mys$y)
  ll <- stats::dnorm(matrix(mys$y[obs], 600, sum(obs), byrow = TRUE), t(fit$draws[obs, ]), sd_e, log = TRUE)
  lppd <- sum(log(colMeans(exp(ll))))
  p_waic <- sum(apply(ll, 2, stats::var))
  # the backend works in single precision
  expect_equal(unname(fit$goodness["WAIC"]), -2 * (lppd - p_waic), tolerance = 1e-4)
})

test_that("engine = 'jax' fits count, spatial, and temporal models end to end", {
  skip_if_no_numpyro()
  data("mys", package = "fastsae")
  data("sim_panel", package = "fastsae")
  data("mys_proxmat", package = "fastsae")

  fit_bym2 <- quiet_mcmc(hb_area(y ~ x1 + x2 + x3,
    data = mys, vardir = "vardir", spatial = "bym2",
    W = mys_knn, engine = jax_eng, print_result = FALSE
  ))
  expect_false(is.null(fit_bym2$phi))
  expect_true(all(is.finite(fit_bym2$df_hb$hb)))

  fit_pois <- quiet_mcmc(hb_area(y_poisson ~ x1 + x2,
    data = sim_panel, domain = "area", time = "year",
    family = "poisson", exposure = "exposure", temporal = "ar1",
    engine = jax_eng, print_result = FALSE
  ))
  expect_true(all(c("time", "rate", "estimated_count") %in% names(fit_pois$df_hb)))
  expect_true(all(is.finite(fit_pois$df_hb$hb) & fit_pois$df_hb$hb > 0))
  expect_false(is.null(fit_pois$rho_time))

  fit_bm <- quiet_mcmc(hb_area(y ~ x1 + x2 + x3,
    data = mys, vardir = "vardir", engine = jax_eng,
    self_benchmark = TRUE, benchmark_weight = "n",
    print_result = FALSE
  ))
  expect_equal(fit_bm$benchmark_summary$discrepancy, 0, tolerance = 1e-8)
})

test_that("engine = 'jax' samples an AR(1) time effect without divergent transitions", {
  skip_if_no_numpyro()
  data("sim_panel", package = "fastsae")
  # With 5 periods the level of a stationary AR(1) effect trades off against the
  # intercept along a ridge whose width depends on rho: NUTS used to report 30-200
  # divergent transitions of 4000 here (R-hat up to 1.04, 0.2 posterior SD from INLA)
  eng <- jax_engine(warmup = 500, samples = 1000, chains = 4, seed = 42, device = "cpu", progress = FALSE)
  fit <- quiet_mcmc(hb_area(y_gaussian ~ x1 + x2,
    vardir = "vardir", data = sim_panel,
    domain = "area", time = "year", temporal = "ar1",
    engine = eng, print_result = FALSE
  ))
  expect_lte(fit$divergences, 20)
  # 1.01 is within Monte Carlo noise of the intercept at ESS ~ 500 (as for temporal = "iid")
  expect_lt(fit$rhat_max, 1.02)

  skip_if_not_installed("INLA")
  ref <- suppressWarnings(hb_area(y_gaussian ~ x1 + x2,
    vardir = "vardir", data = sim_panel,
    domain = "area", time = "year", temporal = "ar1",
    print_result = FALSE
  ))
  expect_lt(max(abs(fit$df_hb$hb - ref$df_hb$hb) / ref$df_hb$sd), 0.15)
})

test_that("engine = 'jax' passes vardir of a gamma model to the backend", {
  env <- new.env()
  local_fake_backend(env)
  df <- sim_engine_df()
  df$y <- exp(df$y)
  df$v <- (0.2 * df$y)^2
  fit <- hb_area(y ~ x,
    data = df, domain = "area", vardir = "v", family = "gamma",
    engine = jax_engine(samples = 50, device = "cpu"), print_result = FALSE
  )
  expect_equal(env$args$vardir, df$v)
})

test_that("engine = 'jax' reports the beta precision of the second stage", {
  env <- new.env()
  local_fake_backend(env)
  df <- sim_engine_df()
  df$y <- stats::plogis(df$y - 1)
  df$v <- 0.002
  fit <- hb_area(y ~ x,
    data = df, domain = "area", vardir = "v", family = "beta",
    engine = jax_engine(samples = 50, device = "cpu"), print_result = FALSE
  )
  expect_equal(fit$df_hb$precision, env$args$beta_phi)
})

test_that("engine = 'jax' gamma model without vardir warns as INLA does", {
  env <- new.env()
  local_fake_backend(env)
  df <- sim_engine_df()
  df$y <- exp(df$y)
  expect_warning(
    hb_area(y ~ x,
      data = df, domain = "area", family = "gamma",
      engine = jax_engine(samples = 50, device = "cpu"), print_result = FALSE
    ),
    "weakly identified"
  )
})

test_that("engine = 'jax' gamma model with vardir agrees with INLA", {
  skip_if_no_numpyro()
  skip_if_not_installed("INLA")
  set.seed(11)
  m <- 40
  d <- data.frame(x = rnorm(m))
  mu <- exp(0.5 + 0.4 * d$x + rnorm(m, sd = 0.5))
  d$y <- rgamma(m, shape = 20, rate = 20 / mu)
  d$v <- mu^2 / 20
  d$y[1:4] <- NA
  fj <- quiet_mcmc(hb_area(y ~ x,
    data = d, vardir = "v", family = "gamma",
    engine = jax_engine(warmup = 1000, samples = 1000, chains = 4, device = "cpu", seed = 1),
    print_result = FALSE
  ))
  fi <- hb_area(y ~ x, data = d, vardir = "v", family = "gamma", print_result = FALSE)
  z <- abs(fj$df_hb$hb - fi$df_hb$hb) / fi$df_hb$sd
  expect_lt(max(z), 0.15)
  expect_lt(max(abs(fj$df_hb$sd / fi$df_hb$sd - 1)), 0.15)
})

sim_nb_vardir <- function(m = 40, seed = 12) {
  set.seed(seed)
  d <- data.frame(x = rnorm(m), E = runif(m, 10, 100))
  mu <- d$E * exp(-0.5 + 0.4 * d$x + rnorm(m, sd = 0.5))
  d$y <- rnbinom(m, mu = mu, size = 10)
  d$v <- mu + mu^2 / 10
  d$y[1:4] <- NA
  d
}

test_that("a negative binomial model with vardir fixes the size of each area", {
  skip_if_not_installed("INLA")
  d <- sim_nb_vardir()
  w <- character(0)
  fit <- withCallingHandlers(
    hb_area(y ~ x, data = d, family = "nbinomial", exposure = "E", vardir = "v", print_result = FALSE),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("weakly identified", w)))
  expect_false(any(grepl("size", rownames(fit$fit$summary.hyperpar))))
  expect_warning(
    hb_area(y ~ x, data = d, family = "nbinomial", exposure = "E", print_result = FALSE),
    "weakly identified"
  )
})

test_that("engine = 'jax' passes vardir of a negative binomial model to the backend", {
  env <- new.env()
  local_fake_backend(env)
  d <- sim_nb_vardir(m = 12)
  d$y[1:4] <- 5
  d$area <- paste0("a", 1:12)
  hb_area(y ~ x,
    data = d, domain = "area", vardir = "v", exposure = "E", family = "nbinomial",
    engine = jax_engine(samples = 50, device = "cpu"), print_result = FALSE
  )
  expect_equal(env$args$vardir, d$v)
})

test_that("engine = 'jax' negative binomial model with vardir agrees with INLA", {
  skip_if_no_numpyro()
  skip_if_not_installed("INLA")
  d <- sim_nb_vardir()
  fj <- quiet_mcmc(hb_area(y ~ x,
    data = d, vardir = "v", exposure = "E", family = "nbinomial",
    engine = jax_engine(warmup = 1000, samples = 1000, chains = 4, device = "cpu", seed = 1),
    print_result = FALSE
  ))
  fi <- hb_area(y ~ x, data = d, vardir = "v", exposure = "E", family = "nbinomial", print_result = FALSE)
  # The first stage (one estimated size) already differs by up to 0.25 posterior SD:
  # INLA's approximation of the weakly identified NB model (see hb_area()); locally the
  # two stages differed by at most 0.17 SD, with a median of 0.03
  z <- abs(fj$df_hb$hb - fi$df_hb$hb) / fi$df_hb$sd
  expect_lt(max(z), 0.25)
  expect_lt(median(z), 0.1)
  expect_lt(max(abs(fj$df_hb$sd / fi$df_hb$sd - 1)), 0.15)
})

test_that("the INLA negative binomial size gets the half-normal(10) prior of jax_engine()", {
  skip_if_not_installed("INLA")
  d <- sim_nb_vardir()
  fit <- suppressWarnings(hb_area(y ~ x, data = d, family = "nbinomial", exposure = "E", print_result = FALSE))
  prior <- fit$fit$.args$control.family[[1]]$hyper$theta$prior
  expect_match(prior, "^expression: s = 10; size = exp\\(theta\\)")
  # a prior given by the user is kept
  own <- list(prior = "pc.mgamma", param = 7)
  fit_own <- suppressWarnings(hb_area(y ~ x,
    data = d, family = "nbinomial", exposure = "E", print_result = FALSE,
    engine = inla_engine(control.family = list(hyper = list(theta = own)))
  ))
  expect_match(fit_own$fit$.args$control.family[[1]]$hyper$theta$prior, "^pc.mgamma$")
})

sim_beta_vardir <- function(m = 30, b0 = -2.5) {
  set.seed(12)
  d <- data.frame(x = rnorm(m), n = sample(10:80, m, replace = TRUE))
  p <- plogis(b0 + 0.5 * d$x + rnorm(m, sd = 0.5))
  d$y <- pmin(pmax(rbeta(m, p * (d$n - 1), (1 - p) * (d$n - 1)), 1e-4), 1 - 1e-4)
  d$v <- p * (1 - p) / d$n
  d$y[1:3] <- NA
  d$v[1:3] <- NA
  d
}

test_that("INLA beta model with vardir takes the precision at the model mean (two stages)", {
  skip_if_not_installed("INLA")
  d <- sim_beta_vardir()
  fit <- hb_area(y ~ x, data = d, vardir = "v", family = "beta", print_result = FALSE)
  obs <- !is.na(d$y)
  phi <- fit$fit$.args$scale[obs]
  # not the plug-in of the direct estimate, which gives the model variance
  # V mu (1 - mu) / (y (1 - y)) instead of V (81-84% coverage for p near 0.08)
  plug_in <- pmax(d$y * (1 - d$y) / d$v - 1, 1)[obs]
  expect_gt(max(abs(phi / plug_in - 1)), 0.2)
  # but the precision at a model mean (the first-stage posterior mean): the beta variance
  # mu (1 - mu) / (phi + 1) is then vardir
  mu <- c(0.05, 0.3, 0.5, 0.9)
  v <- c(0.001, 0.01, 0.02, 0.002)
  expect_equal(mu * (1 - mu) / (.beta_phi_from_vardir(mu, v) + 1), v)
  expect_equal(.beta_phi_from_vardir(c(0.05, NA), c(0.04, 0.01)), c(1, 1))
})

test_that("engine = 'jax' beta model with vardir passes the two-stage precision", {
  env <- new.env()
  local_fake_backend(env)
  d <- sim_beta_vardir()
  d$area <- seq_len(nrow(d))
  hb_area(y ~ x,
    data = d, domain = "area", vardir = "v", family = "beta",
    engine = jax_engine(samples = 50, device = "cpu"), print_result = FALSE
  )
  expect_length(env$args$beta_phi, nrow(d))
  expect_true(all(env$args$beta_phi >= 1))
})

test_that("engine = 'jax' WAIC of a beta model with trials takes the precision n - 1", {
  skip_if_no_numpyro()
  d <- sim_beta_vardir()
  fit <- quiet_mcmc(hb_area(y ~ x,
    data = d, trials = "n", family = "beta", engine = jax_eng, print_result = FALSE
  ))
  obs <- !is.na(d$y)
  # the likelihood's precision phi_d = n_d - 1 (not a fixed 10)
  phi <- matrix(pmax(d$n[obs] - 1, 1), 600, sum(obs), byrow = TRUE)
  p <- pmin(pmax(t(fit$draws[obs, ]), 1e-5), 1 - 1e-5)
  y <- matrix(d$y[obs], 600, sum(obs), byrow = TRUE)
  ll <- stats::dbeta(y, p * phi, (1 - p) * phi, log = TRUE)
  lppd <- sum(log(colMeans(exp(ll))))
  p_waic <- sum(apply(ll, 2, stats::var))
  expect_equal(unname(fit$goodness["WAIC"]), -2 * (lppd - p_waic), tolerance = 1e-4)
})

test_that("engine = 'jax' beta model with vardir agrees with INLA", {
  skip_if_no_numpyro()
  skip_if_not_installed("INLA")
  d <- sim_beta_vardir()
  fj <- quiet_mcmc(hb_area(y ~ x,
    data = d, vardir = "v", family = "beta",
    engine = jax_engine(warmup = 1000, samples = 1000, chains = 4, device = "cpu", seed = 1),
    print_result = FALSE
  ))
  fi <- hb_area(y ~ x, data = d, vardir = "v", family = "beta", print_result = FALSE)
  z <- abs(fj$df_hb$hb - fi$df_hb$hb) / fi$df_hb$sd
  expect_lt(max(z), 0.25)
  expect_lt(max(abs(fj$df_hb$sd / fi$df_hb$sd - 1)), 0.2)
  expect_lt(fj$rhat_max, 1.05)
})
