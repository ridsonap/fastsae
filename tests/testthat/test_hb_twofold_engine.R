# hb_twofold(engine = jax_engine()): the R side is tested against a fake Python
# backend, so these tests run without Python; the end-to-end fits at the bottom
# need jax + numpyro.

tf_engine_data <- function(m = 12, seed = 1) {
  set.seed(seed)
  dat <- do.call(rbind, lapply(seq_len(m), function(d) {
    nd <- sample(2:4, 1)
    data.frame(area = d, subarea = paste0(d, "-", seq_len(nd)), x1 = rnorm(nd),
               vardir = runif(nd, 0.2, 0.8), pop = runif(nd, 50, 500))
  }))
  dat$y <- 1 + dat$x1 + rnorm(m, sd = 0.8)[dat$area] + rnorm(nrow(dat), sd = 0.6) +
    rnorm(nrow(dat), sd = sqrt(dat$vardir))
  dat
}

tf_chain_W <- function(D) {
  W <- matrix(0, D, D)
  W[cbind(1:(D - 1), 2:D)] <- 1
  W + t(W)
}

# Stand-in for hb_area_numpyro.fit_numpyro_hb() with the nested sub-area effect
fake_twofold_backend <- function(env) {
  list(fit_numpyro_hb = function(y, X, domain_idx, ...) {
    args <- list(...)
    env$args <- c(list(y = y, X = X, domain_idx = domain_idx), args)
    n <- nrow(X)
    S <- 300
    set.seed(3)
    eta <- 0.2 + 0.5 * X[, ncol(X)]
    noise <- matrix(stats::rnorm(S * n, sd = 0.1), S, n)
    draws <- switch(args$family,
      binomial = stats::plogis(matrix(rep(eta, each = S), S, n) + noise),
      poisson = exp(matrix(rep(eta, each = S), S, n) + noise),
      matrix(rep(1 + X[, ncol(X)], each = S), S, n) + noise
    )
    list(
      beta_mean = rep(0.5, ncol(X)), beta_sd = rep(0.1, ncol(X)),
      beta_ci_lower = rep(0.3, ncol(X)), beta_ci_upper = rep(0.7, ncol(X)),
      hb_mean = colMeans(draws), hb_sd = apply(draws, 2, stats::sd),
      hb_ci_lower = apply(draws, 2, stats::quantile, 0.025, names = FALSE),
      hb_ci_upper = apply(draws, 2, stats::quantile, 0.975, names = FALSE),
      linpred_mean = colMeans(draws), rand_eff_mean = rep(0, n),
      rand_eff_area_mean = seq_len(n) / 10, rand_eff_subarea_mean = -seq_len(n) / 10,
      hyperparameters = list(sigma2_u = 0.3, sigma2_subarea = 0.1, phi = 0.4),
      dic = 10, p_dic = 3, waic = 11, p_waic = 3.5,
      device_used = "cpu", rhat_max = 1.01, shrinkage_weights = NULL,
      hb_draws = draws
    )
  })
}

local_fake_twofold_backend <- function(env, frame = parent.frame()) {
  local_mocked_bindings(.get_numpyro_backend = function(...) fake_twofold_backend(env), .env = frame)
}

tf_printed <- function(x) {
  out <- NULL
  msgs <- capture_messages(out <- utils::capture.output(print(x)))
  paste(c(msgs, out), collapse = "\n")
}

skip_if_no_numpyro_tf <- function() {
  skip_on_cran()
  skip_if_not_installed("reticulate")
  if (!suppressMessages(check_numpyro_available())) skip("NumPyro / JAX not available")
}

twofold_hb_columns <- c("domain", "subarea", "y", "hb", "linear_pred", "vardir", "sd", "mse", "rse",
                        "ci_lower", "ci_upper", "random_effect_area", "random_effect_subarea")

# ---------------------------------------------------------------------------
# R side of engine = "jax" with a fake backend
# ---------------------------------------------------------------------------

test_that("hb_twofold(engine = 'jax') maps the nested model to the backend", {
  env <- new.env()
  local_fake_twofold_backend(env)
  dat <- tf_engine_data()
  dat$y[2] <- NA
  eng <- jax_engine(samples = 50, device = "cpu")

  fit <- suppressMessages(hb_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                                     data = dat, weight = "pop", engine = eng, print_result = FALSE))
  n <- nrow(dat)
  a <- env$args
  expect_equal(a$domain_idx, dat$area - 1L)
  expect_equal(a$subarea_idx, 0:(n - 1))
  expect_equal(a$num_subareas, n)
  expect_equal(a$D, 12L)
  expect_null(a$time_idx)
  expect_equal(a$temporal, "none")
  expect_equal(a$st_interaction, "none")
  expect_equal(a$spatial, "none")
  expect_equal(a$family, "gaussian")
  expect_equal(a$vardir, dat$vardir)
  expect_true(is.na(a$y[2]))
  expect_null(a$twofold_weights)
  expect_true(a$return_draws)
  expect_equal(a$num_samples, 50L)

  expect_s3_class(fit, c("fastsae_hb_twofold", "fastsae_hb", "fastsae"))
  expect_equal(fit$engine, "jax")
  expect_identical(fit$engine_spec, eng)
  expect_named(fit$df_hb, twofold_hb_columns)
  expect_identical(fit$df_subarea, fit$df_hb)
  expect_equal(fit$df_hb$subarea, dat$subarea)
  expect_equal(fit$df_hb$random_effect_area, seq_len(n) / 10)
  expect_equal(fit$df_hb$random_effect_subarea, -seq_len(n) / 10)
  # area (v_i) and sub-area (u_ij) variances in the notation of Torabi and Rao (2014)
  expect_equal(fit$random_effect_var, c(sigma2_v = 0.3, sigma2_u = 0.1))
  expect_true(fit$convergence)
  expect_equal(dim(fit$draws), c(n, 300))
  expect_equal(fit$level, "subarea")

  # area aggregates from the joint draws with population weights normalised per area
  w <- ave(dat$pop, dat$area, FUN = function(v) v / sum(v))
  area_draws <- rowsum(w * fit$draws, dat$area)
  expect_equal(fit$df_area$domain, 1:12)
  expect_equal(fit$df_area$hb_area, unname(rowMeans(area_draws)), tolerance = 1e-10)
  expect_equal(fit$df_area$sd_area, unname(apply(area_draws, 1, stats::sd)), tolerance = 1e-10)
  expect_equal(fit$df_area$ci_upper_area,
               unname(apply(area_draws, 1, stats::quantile, 0.975)), tolerance = 1e-10)
  expect_equal(fit$df_area$n_subareas, as.integer(table(dat$area)))

  expect_match(tf_printed(fit), "Engine: jax")
  expect_match(tf_printed(fit), "max R-hat = 1.010")
})

test_that("hb_twofold(engine = 'jax') passes counts, spatial structure and MCMC options", {
  env <- new.env()
  local_fake_twofold_backend(env)
  dat <- tf_engine_data()
  dat$n <- 40
  dat$cnt <- stats::rbinom(nrow(dat), dat$n, 0.3)

  fit <- suppressMessages(hb_twofold(cnt ~ x1, domain = "area", data = dat, family = "binomial",
                                     trials = "n", spatial = "bym2", W = tf_chain_W(12),
                                     engine = jax_engine(robust = TRUE, prior_beta = "horseshoe",
                                                         keep_draws = FALSE),
                                     print_result = FALSE))
  a <- env$args
  expect_equal(a$family, "binomial")
  expect_equal(a$trials, dat$n)
  expect_equal(a$spatial, "bym2")
  expect_equal(a$scale_factor, .icar_spectral(tf_chain_W(12))$scale_factor)
  expect_true(a$robust)
  expect_equal(a$prior_beta, "horseshoe")
  expect_null(fit$draws)
  expect_equal(fit$phi, 0.4)
  expect_false(anyNA(fit$df_area$sd_area))
  expect_true(all(fit$df_hb$hb > 0 & fit$df_hb$hb < 1))

  # compute_area = FALSE leaves the area table out
  fit2 <- suppressMessages(hb_twofold(cnt ~ x1, domain = "area", data = dat, family = "binomial",
                                      trials = "n", engine = "jax", compute_area = FALSE,
                                      print_result = FALSE))
  expect_null(fit2$df_area)
})

test_that("hb_twofold(engine = 'jax') rejects INLA settings and keeps the input checks", {
  dat <- tf_engine_data()
  fit_jax <- function(...) hb_twofold(y ~ x1, domain = "area", data = dat, engine = "jax",
                                      print_result = FALSE, ...)
  expect_error(fit_jax(vardir = "vardir", strategy = "gaussian"), "cannot be used with")
  expect_error(fit_jax(vardir = "vardir", n_samples = 50), "cannot be used with")
  expect_error(fit_jax(), "vardir")
  expect_error(fit_jax(vardir = "vardir", family = "binomial"), "trials")
})

# ---------------------------------------------------------------------------
# End-to-end against INLA (needs jax + numpyro)
# ---------------------------------------------------------------------------

test_that("hb_twofold(engine = 'jax') agrees with INLA", {
  skip_if_no_numpyro_tf()
  skip_if_not_installed("INLA")
  dat <- tf_engine_data(m = 20, seed = 2)
  dat$n <- 60
  dat$cnt_b <- stats::rbinom(nrow(dat), dat$n, stats::plogis(-0.5 + 0.5 * dat$x1))
  dat$E <- 40
  dat$cnt_p <- stats::rpois(nrow(dat), dat$E * exp(-1 + 0.4 * dat$x1))
  eng <- jax_engine(warmup = 500, samples = 500, chains = 2, device = "cpu", progress = FALSE)
  # With 4 x 1000 draws and m = 20: correlations >= 0.9995, median |difference| of
  # the area means 0.02-0.07 posterior SD, area SD ratios 0.85-1.14 (priors differ)
  cfgs <- list(
    list(y ~ x1, vardir = "vardir"),
    list(cnt_b ~ x1, family = "binomial", trials = "n"),
    list(cnt_p ~ x1, family = "poisson", exposure = "E")
  )
  for (cfg in cfgs) {
    args <- c(cfg, list(domain = "area", subarea = "subarea", data = dat, weight = "pop", print_result = FALSE))
    fit_jax <- suppressMessages(do.call(hb_twofold, c(args, list(engine = eng))))
    fit_inla <- do.call(hb_twofold, args)
    expect_named(fit_jax$df_hb, names(fit_inla$df_hb))
    expect_named(fit_jax$df_area, names(fit_inla$df_area))
    expect_true(fit_jax$convergence)
    expect_gt(cor(fit_jax$df_hb$hb, fit_inla$df_hb$hb), 0.99)
    expect_gt(cor(fit_jax$df_area$hb_area, fit_inla$df_area$hb_area), 0.99)
    expect_lt(median(abs(fit_jax$df_area$hb_area - fit_inla$df_area$hb_area) / fit_inla$df_area$sd_area), 0.25)
    sd_ratio <- fit_jax$df_area$sd_area / fit_inla$df_area$sd_area
    expect_true(all(sd_ratio > 0.7 & sd_ratio < 1.4))
  }
})
