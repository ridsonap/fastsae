# hb_ratio(engine = ...): the JAX engine of the sub-proportion model. The R side is
# tested against a fake Python backend (no Python needed); the end-to-end fits at the
# bottom need jax + numpyro.

sim_ratio_engine_df <- function(D = 40, seed = 7) {
  set.seed(seed)
  d <- data.frame(area = sprintf("a%02d", D:1), x = rnorm(D), pop = rpois(D, 1e4), prov = rep(1:4, length.out = D))
  d$p0 <- plogis(-2 + 0.3 * rnorm(D))
  d$p0_se <- 0.15 * d$p0
  d$n_poor <- rpois(D, 30) + 1
  d$s_dir <- rbinom(D, d$n_poor, plogis(-1.5 + 0.5 * d$x + rnorm(D, 0, 0.3))) / d$n_poor
  d
}

# Stand-in for hb_area_numpyro.fit_numpyro_hb(): 2 chains x 50 draws of s_d around
# plogis(-1 + 0.3 x), with the convergence diagnostics of three parameters
fake_ratio_backend <- function(env) {
  list(fit_numpyro_hb = function(y, X, domain_idx, ...) {
    env$args <- c(list(y = y, X = X, domain_idx = domain_idx), list(...))
    n <- nrow(X)
    S <- 100
    set.seed(3)
    centre <- stats::plogis(-1 + 0.3 * X[, ncol(X)])
    draws <- matrix(rep(centre, each = S), S, n) * exp(matrix(stats::rnorm(S * n, sd = 0.05), S, n))
    list(
      beta_mean = rep(0.5, ncol(X)), beta_sd = rep(0.1, ncol(X)),
      beta_ci_lower = rep(0.3, ncol(X)), beta_ci_upper = rep(0.7, ncol(X)),
      hb_mean = colMeans(draws), hb_sd = apply(draws, 2, stats::sd),
      hb_ci_lower = apply(draws, 2, stats::quantile, 0.025, names = FALSE),
      hb_ci_upper = apply(draws, 2, stats::quantile, 0.975, names = FALSE),
      linpred_mean = colMeans(stats::qlogis(draws)), rand_eff_mean = rep(0, n),
      hyperparameters = list(sigma2_u = 0.2),
      dic = 10, p_dic = 3, waic = 11, p_waic = 3.5,
      device_used = "cpu", rhat_max = 1.003, shrinkage_weights = NULL,
      diagnostics = list(
        parameter = list("sigma2_u", "beta[2]", "hb[1]"), mean = c(0.2, 0.5, 0.3), sd = c(0.1, 0.1, 0.02),
        rhat = c(1.002, 1.001, 1.003), ess_bulk = c(400, 450, 500), ess_tail = c(410, 460, 510)
      ),
      trace = list(names = list("sigma2_u", "beta[2]", "hb[1]"), draws = array(stats::rnorm(2 * 50 * 3), c(2, 50, 3))),
      num_divergent = 0L, treedepth_hits = 0L, ebfmi = c(0.9, 1.1),
      hb_draws = draws
    )
  })
}

local_fake_ratio_backend <- function(env, frame = parent.frame()) {
  local_mocked_bindings(.get_numpyro_backend = function(...) fake_ratio_backend(env), .env = frame)
}

skip_if_no_numpyro_ratio <- function() {
  skip_on_cran()
  skip_if_not_installed("reticulate")
  if (!suppressMessages(check_numpyro_available())) skip("NumPyro / JAX not available")
}

# ---------------------------------------------------------------------------
# R side of the JAX engine (no Python needed)
# ---------------------------------------------------------------------------

test_that("hb_ratio(engine = 'jax') uses the MCMC draws of the sub-proportion", {
  env <- new.env()
  local_fake_ratio_backend(env)
  d <- sim_ratio_engine_df()
  d$s_dir[3] <- NA # unsampled area
  eng <- jax_engine(samples = 50, chains = 2, device = "cpu", progress = FALSE)
  fit <- suppressMessages(hb_ratio(s_dir ~ x,
    data = d, trials = "n_poor", parent = "p0", domain = "area", engine = eng, print_result = FALSE
  ))
  df <- fit$df_hb
  expect_s3_class(fit, c("fastsae_hb_ratio", "fastsae_hb", "fastsae"))
  expect_identical(fit$engine, "jax")
  expect_equal(env$args$family, "binomial")
  expect_equal(env$args$trials, d$n_poor)
  expect_true(is.na(env$args$y[3]))
  # the draws are those of the sampler: all of them, and exactly the sub-proportion fit
  expect_equal(dim(fit$draws), c(nrow(d), 100L))
  expect_identical(fit$ratio_fit$draws * d$p0, fit$draws)
  expect_equal(df$ratio, fit$ratio_fit$df_hb$hb)
  expect_equal(df$ratio_sd, fit$ratio_fit$df_hb$sd)
  expect_equal(df$hb, df$ratio * d$p0)
  expect_true(all(df$hb <= d$p0) && all(df$ci_upper <= d$p0))
  expect_equal(df$rse, df$rse_cond)
  expect_true(is.na(df$y[3]) && !is.na(df$hb[3]))
  expect_equal(df$y[-3], d$s_dir[-3] * d$p0[-3])
  # convergence of the sub-proportion model
  expect_true(fit$convergence)
  expect_equal(fit$rhat_max, 1.003)
  expect_equal(fit$ess_min, 400)
  expect_equal(fit$divergences, 0L)
  expect_equal(fit$ebfmi, c(0.9, 1.1))
  expect_identical(fit$mcmc_diagnostics, fit$ratio_fit$mcmc_diagnostics)
  expect_equal(fit$mcmc_diagnostics$parameter, c("sigma2_u", "x", "hb[1]"))
  expect_match(fit$method, "NumPyro NUTS")
  expect_null(fit$inla_mode_status)
  out <- paste(capture_messages(utils::capture.output(print(fit))), collapse = "\n")
  expect_match(out, "Ratio model")
  expect_match(out, "max R-hat = 1.003")
  expect_equal(summary(fit)$rhat_max, 1.003)
  # trace and rank plots of the sub-proportion model; hb[j] is s_d of row j
  p <- autoplot(fit, type = "trace", pars = c("sigma2_u", "hb[4]"))
  expect_s3_class(p, "ggplot")
  hb4 <- p$data[startsWith(as.character(p$data$parameter), "hb[4]"), ]
  expect_equal(hb4$value[hb4$chain == "2"], fit$ratio_fit$draws[4, 51:100])
  expect_s3_class(autoplot(fit, type = "rank", bins = 10), "ggplot")
})

test_that("hb_ratio(engine = 'jax') benchmarks the draws and adds the parent's uncertainty", {
  env <- new.env()
  local_fake_ratio_backend(env)
  d <- sim_ratio_engine_df()
  fit <- suppressMessages(hb_ratio(s_dir ~ x,
    data = d, trials = "n_poor", parent = "p0", parent_se = "p0_se", domain = "area",
    engine = "jax", benchmark_target = 0.03, benchmark_weight = "pop", print_result = FALSE
  ))
  df <- fit$df_hb
  expect_equal(sum(d$pop * df$hb) / sum(d$pop), 0.03, tolerance = 1e-8)
  expect_equal(unname(fit$benchmark$achieved), 0.03, tolerance = 1e-8)
  expect_true(all(df$hb <= d$p0))
  expect_true(all(df$rse > df$rse_cond))
  # the benchmark moves s_d away from the sub-proportion fit, which stays as fitted
  expect_false(isTRUE(all.equal(df$ratio, fit$ratio_fit$df_hb$hb)))
  # the parent draws follow `seed`, the sampler jax_engine(seed = )
  fit2 <- suppressMessages(hb_ratio(s_dir ~ x,
    data = d, trials = "n_poor", parent = "p0", parent_se = "p0_se", domain = "area",
    engine = "jax", benchmark_target = 0.03, benchmark_weight = "pop", seed = 5, print_result = FALSE
  ))
  expect_identical(fit2$df_hb$hb, df$hb)
  expect_false(identical(fit2$df_hb$sd, df$sd))
  tg <- c("1" = 0.02, "2" = 0.03, "3" = 0.025, "4" = 0.035)
  fg <- suppressMessages(hb_ratio(s_dir ~ x,
    data = d, trials = "n_poor", parent = "p0", domain = "area", engine = "jax",
    benchmark_target = tg, benchmark_weight = "pop", benchmark_group = "prov", print_result = FALSE
  ))
  ach <- tapply(d$pop * fg$df_hb$hb, d$prov, sum) / tapply(d$pop, d$prov, sum)
  expect_equal(as.numeric(ach), as.numeric(tg), tolerance = 1e-8)
})

test_that("hb_ratio(engine = 'jax') ignores n_draws and keeps its draws without keep_draws", {
  env <- new.env()
  local_fake_ratio_backend(env)
  d <- sim_ratio_engine_df(D = 12)
  args <- list(formula = s_dir ~ x, data = d, trials = "n_poor", parent = "p0", print_result = FALSE)
  run <- function(...) suppressMessages(do.call(hb_ratio, utils::modifyList(args, list(...))))
  expect_warning(fit <- run(engine = "jax", n_draws = 500), "n_draws")
  expect_equal(ncol(fit$draws), 100L)
  expect_no_warning(run(engine = "jax"))
  # keep_draws = FALSE drops the draws of s_d from the sub-proportion fit only
  eng <- jax_engine(device = "cpu", keep_draws = FALSE, progress = FALSE)
  fit <- run(engine = eng)
  expect_equal(dim(fit$draws), c(12L, 100L))
  expect_null(fit$ratio_fit$draws)
  expect_identical(fit$ratio_fit$engine_spec, eng)
  expect_error(autoplot(fit, type = "trace", pars = "hb[4]"), "keep_draws")
  expect_s3_class(autoplot(fit, type = "trace"), "ggplot")
  # models and arguments the JAX engine does not have fail as in hb_area()
  expect_error(run(engine = "jax", spatial = "slm", W = diag(12)), "slm")
  expect_error(run(engine = "jax", strategy = "gaussian"), "INLA argument")
})

test_that("an INLA hb_ratio fit has no MCMC diagnostics", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("sn")
  d <- sim_ratio_engine_df(D = 20)
  fit <- hb_ratio(s_dir ~ x,
    data = d, trials = "n_poor", parent = "p0", engine = inla_engine(strategy = "gaussian"),
    n_draws = 50, print_result = FALSE
  )
  expect_identical(fit$engine, "inla")
  expect_identical(fit$ratio_fit$engine_spec$strategy, "gaussian")
  expect_equal(ncol(fit$draws), 50L)
  expect_null(fit$rhat_max)
  expect_false("mcmc_diagnostics" %in% names(fit))
  expect_error(autoplot(fit, type = "trace"), "jax_engine")
})

# ---------------------------------------------------------------------------
# End to end (needs jax + numpyro)
# ---------------------------------------------------------------------------

test_that("hb_ratio(engine = 'jax') agrees with the INLA engine", {
  skip_if_no_numpyro_ratio()
  skip_if_not_installed("INLA")
  skip_if_not_installed("sn")
  d <- sim_ratio_engine_df()
  d$s_dir[3] <- NA # unsampled area
  eng <- jax_engine(warmup = 1000, samples = 1000, chains = 4, device = "cpu", seed = 1, progress = FALSE)
  fit_ratio <- function(engine, ...) {
    suppressWarnings(suppressMessages(hb_ratio(s_dir ~ x,
      data = d, trials = "n_poor", parent = "p0", parent_se = "p0_se", domain = "area",
      engine = engine, print_result = FALSE, ...
    )), classes = "fastsae_mcmc_warning")
  }
  fj <- fit_ratio(eng)
  fi <- fit_ratio("inla", n_draws = 4000)
  expect_identical(names(fj$df_hb), names(fi$df_hb))
  expect_equal(dim(fj$draws), c(nrow(d), 4000L))
  expect_lt(fj$rhat_max, 1.05)
  expect_true(fj$convergence)
  expect_true(all(fj$df_hb$hb <= d$p0))
  expect_equal(fj$df_hb$ratio, fj$ratio_fit$df_hb$hb, tolerance = 1e-6)
  expect_lt(max(abs(fj$df_hb$hb - fi$df_hb$hb) / fi$df_hb$sd), 0.25)
  expect_lt(max(abs(fj$df_hb$sd / fi$df_hb$sd - 1)), 0.2)
  expect_lt(max(abs(fj$df_hb$ratio_sd / fi$df_hb$ratio_sd - 1)), 0.2)

  # benchmarked: both engines reach the target, and still agree
  bj <- fit_ratio(eng, benchmark_target = 0.03, benchmark_weight = "pop")
  bi <- fit_ratio("inla", n_draws = 4000, benchmark_target = 0.03, benchmark_weight = "pop")
  expect_equal(sum(d$pop * bj$df_hb$hb) / sum(d$pop), 0.03, tolerance = 1e-8)
  expect_true(all(bj$df_hb$hb <= d$p0))
  expect_lt(max(abs(bj$df_hb$hb - bi$df_hb$hb) / bi$df_hb$sd), 0.25)
  expect_lt(max(abs(bj$df_hb$sd / bi$df_hb$sd - 1)), 0.2)
})
