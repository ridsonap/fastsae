# hb_unit(engine = jax_engine()): the R side is tested against a fake Python
# backend, so these tests run without Python; the end-to-end fits at the bottom
# need jax + numpyro.

# Sample of 6 domains (domain 6 unsampled) and their population means
unit_engine_data <- function(seed = 1) {
  set.seed(seed)
  n_d <- c(4, 2, 5, 3, 1, 0)
  s <- data.frame(dom = rep(seq_along(n_d), n_d))
  s$x1 <- rnorm(nrow(s), s$dom / 3)
  s$y <- 1 + 2 * s$x1 + rnorm(6)[s$dom] + rnorm(nrow(s))
  s$yb <- rbinom(nrow(s), 1, 0.4)
  pop <- data.frame(dom = seq_along(n_d), x1 = seq_along(n_d) / 3, N = c(40, 10, 60, 30, 20, 50))
  list(sample = s, pop = pop)
}

# Stand-in for hb_area_numpyro.fit_numpyro_hb(): returns the draws of the rows in draw_rows
fake_unit_backend <- function(env) {
  list(fit_numpyro_hb = function(y, X, domain_idx, ...) {
    args <- list(...)
    env$args <- c(list(y = y, X = X, domain_idx = domain_idx), args)
    n <- nrow(X)
    S <- 200
    set.seed(3)
    eta <- 0.2 + 0.5 * X[, ncol(X)]
    noise <- matrix(stats::rnorm(S * n, sd = 0.1), S, n)
    draws <- switch(args$family,
      binomial = stats::plogis(matrix(rep(eta, each = S), S, n) + noise),
      matrix(rep(1 + X[, ncol(X)], each = S), S, n) + noise
    )
    rows <- args$draw_rows + 1L
    gauss <- args$family == "gaussian"
    list(
      beta_mean = rep(0.5, ncol(X)), beta_sd = rep(0.1, ncol(X)),
      beta_ci_lower = rep(0.3, ncol(X)), beta_ci_upper = rep(0.7, ncol(X)),
      hb_mean = colMeans(draws), hb_sd = apply(draws, 2, stats::sd),
      hb_ci_lower = apply(draws, 2, stats::quantile, 0.025, names = FALSE),
      hb_ci_upper = apply(draws, 2, stats::quantile, 0.975, names = FALSE),
      linpred_mean = colMeans(draws), rand_eff_mean = seq_len(n) / 10,
      hyperparameters = c(
        list(sigma2_u = 0.3), if (gauss) list(sigma2_e = 1.5),
        if (args$spatial == "bym2") list(phi = 0.4)
      ),
      dic = 10, p_dic = 3, waic = 11, p_waic = 3.5,
      device_used = "cpu", rhat_max = 1.01, shrinkage_weights = NULL,
      hb_draws = draws[, rows, drop = FALSE],
      sigma2_e_draws = if (gauss) seq(1, 2, length.out = S)
    )
  })
}

local_fake_unit_backend <- function(env, frame = parent.frame()) {
  local_mocked_bindings(.get_numpyro_backend = function(...) fake_unit_backend(env), .env = frame)
}

skip_if_no_numpyro_unit <- function() {
  skip_on_cran()
  skip_if_not_installed("reticulate")
  if (!suppressMessages(check_numpyro_available())) skip("NumPyro / JAX not available")
}

unit_hb_columns <- c(
  "domain", "y", "hb", "linear_pred", "vardir", "sd", "mse", "rse", "ci_lower",
  "ci_upper", "samp_size", "sample_mean", "random_effect"
)

# ---------------------------------------------------------------------------
# R side of engine = "jax" with a fake backend
# ---------------------------------------------------------------------------

test_that("hb_unit(engine = 'jax') fits the sampled units and predicts each domain", {
  env <- new.env()
  local_fake_unit_backend(env)
  d <- unit_engine_data()
  n <- nrow(d$sample)
  eng <- jax_engine(samples = 50, device = "cpu")

  fit <- suppressMessages(hb_unit(y ~ x1,
    unit_data = d$sample, Xpop = d$pop, domain_var = "dom",
    engine = eng, print_result = FALSE
  ))
  a <- env$args
  # the sampled units, then one prediction row per domain at the population mean
  expect_equal(a$y, c(d$sample$y, rep(NA, 6)))
  expect_equal(unname(a$X[, "x1"]), c(d$sample$x1, d$pop$x1))
  expect_equal(a$domain_idx, c(d$sample$dom, 1:6) - 1L)
  expect_equal(a$draw_rows, n + 0:5)
  expect_equal(a$D, 6L)
  expect_null(a$trials)
  expect_null(a$vardir)
  expect_null(a$subarea_idx)
  expect_equal(a$spatial, "none")
  expect_equal(a$temporal, "none")
  expect_equal(a$num_samples, 50L)

  expect_s3_class(fit, c("fastsae_hb_unit", "fastsae_hb", "fastsae"))
  expect_equal(fit$engine, "jax")
  expect_identical(fit$engine_spec, eng)
  expect_named(fit$df_hb, unit_hb_columns)
  # without popsize_var the domain means are the model means at the prediction rows
  expect_equal(fit$df_hb$hb, unname(rowMeans(fit$draws)))
  expect_equal(fit$df_hb$random_effect, (n + 1:6) / 10)
  expect_equal(fit$df_hb$samp_size, c(4L, 2L, 5L, 3L, 1L, 0L))
  expect_equal(fit$df_hb$vardir, c(1.5 / c(4, 2, 5, 3, 1), NA))
  expect_equal(fit$random_effect_var, 0.3)
  expect_equal(fit$residual_var, 1.5)
  expect_equal(fit$unsampled_domains, 6L)
  expect_equal(dim(fit$draws), c(6, 200))
  expect_true(fit$convergence)
  expect_true(is.na(fit$goodness$mlik))

  msgs <- capture_messages(utils::capture.output(print(fit)))
  expect_match(paste(msgs, collapse = "\n"), "Engine: jax")
})

test_that("hb_unit(engine = 'jax') draws the finite population mean", {
  env <- new.env()
  local_fake_unit_backend(env)
  d <- unit_engine_data()
  n <- nrow(d$sample)
  fit <- suppressMessages(hb_unit(y ~ x1,
    unit_data = d$sample, Xpop = d$pop, domain_var = "dom",
    popsize_var = "N", engine = jax_engine(seed = 9), print_result = FALSE
  ))

  # prediction rows at the mean covariates of the non-sampled units
  n_d <- c(4, 2, 5, 3, 1, 0)
  xbar_s <- unname(c(tapply(d$sample$x1, d$sample$dom, mean), 0))
  expect_equal(unname(env$args$X[n + 1:6, "x1"]), (d$pop$N * d$pop$x1 - n_d * xbar_s) / (d$pop$N - n_d))

  # Ybar_d = f_d ybar_s + (1 - f_d) (mu_r + ebar_r), ebar_r ~ N(0, sigma_e^2 / (N_d - n_d))
  mu <- t(fake_unit_backend(new.env())$fit_numpyro_hb(
    env$args$y, env$args$X, env$args$domain_idx,
    family = "gaussian", spatial = "none",
    draw_rows = n + 0:5
  )$hb_draws)
  s2e <- seq(1, 2, length.out = 200)
  z <- .with_seed(9, matrix(rnorm(6 * 200), 6))
  f <- n_d / d$pop$N
  ybar <- unname(c(tapply(d$sample$y, d$sample$dom, mean), 0))
  theta <- f * ybar + (1 - f) * (mu + z * sqrt(outer(1 / (d$pop$N - n_d), s2e)))
  expect_equal(fit$draws, theta)
  expect_equal(fit$df_hb$hb, rowMeans(theta))
  expect_equal(fit$df_hb$sd, apply(theta, 1, sd))
  expect_equal(fit$df_hb$ci_upper, apply(theta, 1, quantile, 0.975, names = FALSE))
  expect_equal(fit$df_hb$pop_size, d$pop$N)
  expect_equal(fit$df_hb$estimated_total, fit$df_hb$hb * d$pop$N)
})

test_that("hb_unit(engine = 'jax') passes Bernoulli trials and the spatial structure", {
  env <- new.env()
  local_fake_unit_backend(env)
  d <- unit_engine_data()
  W <- matrix(0, 6, 6)
  W[cbind(1:5, 2:6)] <- 1
  W <- W + t(W)
  fit <- suppressWarnings(suppressMessages(hb_unit(
    yb ~ x1,
    unit_data = d$sample, Xpop = d$pop, domain_var = "dom", family = "binomial",
    spatial = "bym2", W = W, popsize_var = "N", engine = jax_engine(keep_draws = FALSE),
    print_result = FALSE
  )))
  a <- env$args
  expect_equal(a$trials, rep(1, nrow(d$sample) + 6))
  expect_equal(a$spatial, "bym2")
  expect_equal(a$scale_factor, .icar_spectral(W)$scale_factor)
  expect_null(fit$draws)
  expect_equal(fit$phi, 0.4)
  expect_true(is.na(fit$residual_var))
  expect_true(all(fit$df_hb$hb > 0 & fit$df_hb$hb < 1))
})

test_that("hb_unit() checks the unit responses of the binomial and Poisson models", {
  d <- unit_engine_data()
  d$sample$cnt <- 2
  expect_error(suppressWarnings(hb_unit(cnt ~ x1,
    unit_data = d$sample, Xpop = d$pop, domain_var = "dom",
    family = "binomial", engine = "jax"
  )), "0 or 1")
  d$sample$cnt <- 1.5
  expect_error(suppressWarnings(hb_unit(cnt ~ x1,
    unit_data = d$sample, Xpop = d$pop, domain_var = "dom",
    family = "poisson", engine = "jax"
  )), "integer counts")
  expect_error(hb_unit(y ~ x1,
    unit_data = d$sample, Xpop = d$pop, domain_var = "dom",
    engine = "jax", strategy = "gaussian"
  ), "cannot be used with")
})

# ---------------------------------------------------------------------------
# End-to-end (needs jax + numpyro)
# ---------------------------------------------------------------------------

test_that("hb_unit(engine = 'jax') matches Stan for the same model and priors", {
  skip_if_no_numpyro_unit()
  # fixtures/make-hbunit-stan.R: 4 x 5000 Stan draws, R-hat < 1.01, ESS > 6000.
  # With 4 x 2000 NumPyro draws the domain means differed by at most 0.035 posterior
  # SD, the SDs by 2-5% and sigma2_u by 0-2%: Monte Carlo error.
  fx <- readRDS(test_path("fixtures", "hbunit-stan.rds"))
  d <- fx$data
  eng <- jax_engine(warmup = 500, samples = 1000, chains = 2, device = "cpu", progress = FALSE)
  cfgs <- list(
    gaussian = list(y_gauss ~ x1),
    gaussian_bym2 = list(y_gauss ~ x1, spatial = "bym2", W = d$W),
    binomial = list(y_bin ~ x1, family = "binomial"),
    poisson = list(y_pois ~ x1, family = "poisson")
  )
  for (k in names(cfgs)) {
    fit <- quiet_mcmc(suppressWarnings(do.call(hb_unit, c(cfgs[[k]], list(
      unit_data = d$sample, Xpop = d$pop, domain_var = "dom", engine = eng, print_result = FALSE
    )))))
    ref <- fx$ref[[k]]
    z <- abs(fit$df_hb$hb - ref$mean) / ref$sd
    expect_true(fit$convergence, label = k)
    expect_lt(median(z), 0.1, label = k)
    expect_lt(max(z), 0.3, label = k)
    sd_ratio <- fit$df_hb$sd / ref$sd
    expect_true(all(sd_ratio > 0.8 & sd_ratio < 1.25), label = k)
    expect_equal(fit$random_effect_var, ref$sigma2_u, tolerance = 0.15, label = k)
    expect_equal(fit$estcoef$beta, unname(ref$beta), tolerance = 0.05, label = k)
    if (startsWith(k, "gaussian")) expect_equal(fit$residual_var, ref$sigma2_e, tolerance = 0.05, label = k)
  }
})

test_that("hb_unit(engine = 'jax') agrees with hbsae and INLA", {
  skip_if_no_numpyro_unit()
  skip_if_not_installed("hbsae")
  skip_if_not_installed("INLA")
  fx <- readRDS(test_path("fixtures", "hbunit-stan.rds"))
  d <- fx$data
  s <- d$sample
  eng <- jax_engine(warmup = 500, samples = 1000, chains = 2, device = "cpu", progress = FALSE)
  # 30 domains, three unsampled. With 4 x 2000 draws the domain means differed from
  # hbsae's (numerical integration over the variance ratio, flat priors) by at most
  # 0.094 posterior SD and from INLA's by at most 0.051, SD ratios 0.94-1.08.
  Xp <- cbind(1, d$pop$x1)
  rownames(Xp) <- d$pop$dom
  h <- hbsae::fSAE.Unit(s$y_gauss, cbind(1, s$x1), factor(s$dom),
    Narea = setNames(d$pop$N, d$pop$dom),
    Xpop = Xp, fpc = TRUE, method = "HB", silent = TRUE
  )
  fit <- suppressMessages(hb_unit(y_gauss ~ x1,
    unit_data = s, Xpop = d$pop, domain_var = "dom",
    popsize_var = "N", engine = eng, print_result = FALSE
  ))
  ref_sd <- sqrt(as.numeric(h$mse[rownames(Xp)]))
  z <- abs(fit$df_hb$hb - as.numeric(h$est[rownames(Xp)])) / ref_sd
  expect_lt(median(z), 0.1)
  expect_lt(max(z), 0.3)
  expect_true(all(fit$df_hb$sd / ref_sd > 0.85 & fit$df_hb$sd / ref_sd < 1.2))

  for (k in c("binomial", "poisson")) {
    args <- list(as.formula(paste(if (k == "binomial") "y_bin" else "y_pois", "~ x1")),
      unit_data = s,
      Xpop = d$pop, domain_var = "dom", family = k, print_result = FALSE
    )
    fit_jax <- quiet_mcmc(suppressWarnings(do.call(hb_unit, c(args, list(engine = eng)))))
    fit_inla <- suppressWarnings(do.call(hb_unit, args))
    expect_named(fit_jax$df_hb, names(fit_inla$df_hb))
    z <- abs(fit_jax$df_hb$hb - fit_inla$df_hb$hb) / fit_inla$df_hb$sd
    expect_lt(median(z), 0.1, label = k)
    expect_lt(max(z), 0.3, label = k)
    sd_ratio <- fit_jax$df_hb$sd / fit_inla$df_hb$sd
    expect_true(all(sd_ratio > 0.85 & sd_ratio < 1.2), label = k)
  }
})
