# Default priors of the HB models (.default_prior_prec), shared by inla_engine() and
# jax_engine(): half-normal random-effect SDs, scale sd(y) for Gaussian models and 1 otherwise.

# scale s of the half-normal default prior, read back from its INLA expression
hn_scale <- function(prior) as.numeric(sub("^expression: s = ([^;]+);.*$", "\\1", prior$prior))

test_that("the default SD prior is half-normal(sd y) for Gaussian models, half-normal(1) otherwise", {
  y <- c(10, 14, NA, 22, 30)
  expect_equal(hn_scale(.default_prior_prec(y, "gaussian")), sd(y, na.rm = TRUE))
  expect_equal(hn_scale(.default_prior_prec(y, "binomial")), 1)
  expect_equal(hn_scale(.default_prior_prec(NULL, "gaussian")), 1)
  expect_equal(hn_scale(.default_prior_prec(c(3, 3, NA), "gaussian")), 1)
  # the expression is the half-normal density of sigma = exp(-theta / 2) on theta
  s <- 2.5
  logdens <- function(theta) {
    sigma <- exp(-theta / 2)
    log(2) - log(s) - 0.5 * log(2 * pi) - 0.5 * (sigma / s)^2 + log(sigma / 2)
  }
  expect_equal(integrate(function(t) exp(logdens(t)), -40, 40)$value, 1, tolerance = 1e-6)
  expect_equal(
    integrate(function(t) exp(logdens(t)) * (exp(-t / 2) > s), -40, 40)$value,
    2 * stats::pnorm(-1),
    tolerance = 1e-4
  )
})

test_that("hb_area(), hb_unit() and hb_twofold() use the scaled prior unless one is given", {
  skip_if_not_installed("INLA")
  data("mys", package = "fastsae")
  s_y <- sd(mys$y, na.rm = TRUE)
  fit <- hb_area(y ~ x1 + x2 + x3, data = mys, vardir = "vardir", print_result = FALSE)
  expect_equal(hn_scale(fit$engine_spec$prior_prec), s_y)
  expect_equal(hn_scale(fit$engine_spec$prior_prec_time), s_y)
  own <- list(prior = "pc.prec", param = c(1, 0.01))
  fit_own <- hb_area(y ~ x1 + x2 + x3,
    data = mys, vardir = "vardir", print_result = FALSE,
    engine = inla_engine(prior_prec = own)
  )
  expect_identical(fit_own$engine_spec$prior_prec, own)

  set.seed(1)
  sub <- data.frame(
    area = rep(1:8, each = 3), sub = 1:24, x = rnorm(24),
    vardir = runif(24, 0.2, 0.5)
  )
  sub$y <- 50 + 10 * sub$x + rnorm(8, sd = 5)[sub$area] + rnorm(24, sd = 3)
  tf <- hb_twofold(y ~ x,
    vardir = "vardir", domain = "area", subarea = "sub", data = sub,
    print_result = FALSE
  )
  expect_equal(hn_scale(tf$engine_spec$prior_prec), sd(sub$y))
  expect_equal(hn_scale(tf$engine_spec$prior_prec_subarea), sd(sub$y))
})

test_that("Gaussian hb_area() estimates do not depend on the unit of the response", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  data("mys", package = "fastsae")
  d10 <- transform(mys, y = 10 * y, vardir = 100 * vardir)
  f1 <- hb_area(y ~ x1 + x2 + x3, data = mys, vardir = "vardir", print_result = FALSE)
  f10 <- hb_area(y ~ x1 + x2 + x3, data = d10, vardir = "vardir", print_result = FALSE)
  # With the unscaled PC(1, 0.01) prior the area SD shrank in the larger unit: the
  # estimates moved by up to 4.1 posterior SD (median 0.74) and the SDs fell to 0.17-0.93
  expect_equal(f10$df_hb$hb / 10, f1$df_hb$hb, tolerance = 1e-3)
  expect_equal(f10$df_hb$sd / 10, f1$df_hb$sd, tolerance = 1e-2)
})

test_that("jax_engine() passes the prior set to the backend", {
  expect_equal(jax_engine()$priors, "inla")
  expect_equal(jax_engine(priors = "weakly_informative")$priors, "weakly_informative")
  expect_error(jax_engine(priors = "flat"), "should be one of")

  env <- new.env()
  local_mocked_bindings(.get_numpyro_backend = function(...) {
    list(fit_numpyro_hb = function(y, X, domain_idx, ...) {
      env$args <- list(...)
      n <- nrow(X)
      draws <- matrix(rnorm(100 * n, 1, 0.1), 100, n)
      list(
        beta_mean = rep(0.5, ncol(X)), beta_sd = rep(0.1, ncol(X)),
        beta_ci_lower = rep(0.3, ncol(X)), beta_ci_upper = rep(0.7, ncol(X)),
        hb_mean = colMeans(draws), hb_sd = apply(draws, 2, sd),
        hb_ci_lower = apply(draws, 2, quantile, 0.025, names = FALSE),
        hb_ci_upper = apply(draws, 2, quantile, 0.975, names = FALSE),
        linpred_mean = colMeans(draws), rand_eff_mean = rep(0, n),
        hb_draws = draws, hyperparameters = list(sigma2_u = 0.3),
        rhat_max = 1.001, dic = 1, p_dic = 1, waic = 1, p_waic = 1, device_used = "cpu",
        shrinkage_weights = NULL
      )
    })
  })
  data("mys", package = "fastsae")
  hb_area(y ~ x1, data = mys, vardir = "vardir", engine = jax_engine(samples = 100), print_result = FALSE)
  expect_equal(env$args$priors, "inla")
  hb_area(y ~ x1,
    data = mys, vardir = "vardir", print_result = FALSE,
    engine = jax_engine(samples = 100, priors = "weakly_informative")
  )
  expect_equal(env$args$priors, "weakly_informative")
})

skip_if_no_numpyro_backend <- function() {
  skip_on_cran()
  skip_if_not_installed("reticulate")
  if (!suppressMessages(check_numpyro_available())) skip("NumPyro / JAX not available")
}

test_that("the JAX PC prior of the BYM2 phi has the shape of INLA's", {
  skip_if_no_numpyro_backend()
  skip_if_not_installed("INLA")
  W <- matrix(0, 12, 12)
  W[cbind(1:11, 2:12)] <- 1
  W <- W + t(W)
  W[1, 5] <- W[5, 1] <- 1
  sp <- .icar_spectral(W)
  mod <- .get_numpyro_backend(device = "cpu")
  st <- mod$pc_bym2_phi_setup(sp$eig_values, sp$scale_factor, 0.5, 0.5)
  phi <- c(0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99)
  ours <- vapply(phi, function(p) mod$pc_bym2_phi_logdens(p, st[[1]], st[[2]]), numeric(1))
  ref <- INLA:::inla.pc.bym.phi(graph = INLA::inla.read.graph(W), u = 0.5, alpha = 0.5)(phi)
  # equal up to INLA's normalisation over its truncated support
  expect_lt(diff(range(ours - ref)), 0.01)
})

test_that("with the default priors JAX and INLA fit the same model", {
  skip_if_no_numpyro_backend()
  skip_if_not_installed("INLA")
  data("mys", package = "fastsae")
  eng <- jax_engine(warmup = 1000, samples = 1000, chains = 4, device = "cpu", seed = 1)
  for (sp in c("none", "bym2")) {
    W <- if (sp == "bym2") mys_proxmat else NULL
    fj <- suppressMessages(hb_area(y ~ x1 + x2 + x3,
      data = mys, vardir = "vardir", spatial = sp, W = W,
      engine = eng, print_result = FALSE
    ))
    fi <- hb_area(y ~ x1 + x2 + x3,
      data = mys, vardir = "vardir", spatial = sp, W = W,
      print_result = FALSE
    )
    z <- abs(fj$df_hb$hb - fi$df_hb$hb) / fi$df_hb$sd
    expect_lt(max(z), 0.15)
    expect_lt(max(abs(fj$df_hb$sd / fi$df_hb$sd - 1)), 0.1)
    expect_lt(fj$rhat_max, 1.02)
  }
})
