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
    hb_area(yc ~ x, data = df, family = "poisson", engine = "jax", self_benchmark = TRUE,
            benchmark_method = "wfq", print_result = FALSE),
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
  unit_defaults <- list(prior_prec = list(prior = "loggamma", param = c(0.01, 0.01)),
                        prior_phi = list(prior = "pc", param = c(0.5, 0.5)))
  eng <- .resolve_inla_engine("inla", TRUE, list(), "hb_unit", defaults = unit_defaults,
                              used = c("scale_model", "prior_prec", "prior_phi"))
  expect_identical(eng$prior_prec, unit_defaults$prior_prec)
  eng <- .resolve_inla_engine(inla_engine(prior_prec = pc), FALSE, list(), "hb_unit",
                              defaults = unit_defaults, used = c("scale_model", "prior_prec", "prior_phi"))
  expect_identical(eng$prior_prec, pc)

  # settings the function does not use are ignored with a warning
  expect_warning(
    .resolve_inla_engine(inla_engine(prior_rho_time = pc, link = "log"), FALSE, list(), "hb_unit",
                         defaults = unit_defaults, used = c("scale_model", "prior_prec", "prior_phi")),
    "not used by"
  )

  # the old hb_twofold() arguments are renamed
  expect_warning(
    eng <- .resolve_inla_engine("inla", TRUE, list(n_samples = 50, prior_prec_area = pc, strategy = "gaussian"),
                                "hb_twofold", defaults = list(draws = 200L),
                                used = c("scale_model", "prior_prec", "prior_prec_subarea", "prior_phi", "draws")),
    "deprecated"
  )
  expect_identical(eng$draws, 50L)
  expect_identical(eng$prior_prec, pc)
  expect_identical(eng$strategy, "gaussian")
  eng <- .resolve_inla_engine("inla", TRUE, list(), "hb_twofold", defaults = list(draws = 200L),
                              used = c("scale_model", "prior_prec", "prior_prec_subarea", "prior_phi", "draws"))
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

  fit <- hb_area(y ~ x, data = df, domain = "area", vardir = "vardir", spatial = "leroux",
                 W = chain_W(12), engine = jax_engine(samples = 50, device = "cpu", prior_beta = "horseshoe"),
                 print_result = FALSE)

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
  fit <- hb_area(y ~ x, data = df, domain = "area", vardir = "vardir", engine = "jax",
                 self_benchmark = TRUE, benchmark_weight = "w", benchmark_group = "prov",
                 print_result = FALSE)
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
    fit <- hb_area(p ~ x, data = df, domain = "id", family = "binomial", trials = "n",
                   engine = jax_engine(keep_draws = FALSE), print_result = FALSE),
    "proportions"
  )
  expect_equal(env$args$y, round(df$p * df$n))
  expect_equal(env$args$family, "binomial")
  expect_null(fit$draws)
  expect_true(all(c("trials", "estimated_total") %in% names(fit$df_hb)))
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
  expect_identical(as.integer(reticulate::py_to_r(jnp$asarray(ours[[1]])$tolist())),
                   as.integer(reticulate::py_to_r(jnp$asarray(ref[[1]])$tolist())))
  expect_identical(as.integer(reticulate::py_to_r(jnp$asarray(ours[[2]])$tolist())),
                   as.integer(reticulate::py_to_r(jnp$asarray(ref[[2]])$tolist())))
})

test_that("engine = 'jax' Fay-Herriot agrees with INLA", {
  skip_if_no_numpyro()
  skip_if_not_installed("INLA")
  data("mys", package = "fastsae")
  fit_jax <- suppressMessages(hb_area(y ~ x1 + x2 + x3, data = mys, vardir = "vardir",
                                      engine = jax_eng, print_result = FALSE))
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

  fit_bym2 <- suppressMessages(hb_area(y ~ x1 + x2 + x3, data = mys, vardir = "vardir", spatial = "bym2",
                                       W = mys_knn, engine = jax_eng, print_result = FALSE))
  expect_false(is.null(fit_bym2$phi))
  expect_true(all(is.finite(fit_bym2$df_hb$hb)))

  fit_pois <- suppressMessages(hb_area(y_poisson ~ x1 + x2, data = sim_panel, domain = "area", time = "year",
                                       family = "poisson", exposure = "exposure", temporal = "ar1",
                                       engine = jax_eng, print_result = FALSE))
  expect_true(all(c("time", "rate", "estimated_count") %in% names(fit_pois$df_hb)))
  expect_true(all(is.finite(fit_pois$df_hb$hb) & fit_pois$df_hb$hb > 0))
  expect_false(is.null(fit_pois$rho_time))

  fit_bm <- suppressMessages(hb_area(y ~ x1 + x2 + x3, data = mys, vardir = "vardir", engine = jax_eng,
                                     self_benchmark = TRUE, benchmark_weight = "n",
                                     print_result = FALSE))
  expect_equal(fit_bm$benchmark_summary$discrepancy, 0, tolerance = 1e-8)
})
