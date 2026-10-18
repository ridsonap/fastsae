# hb_area(spatial = "sar"): the SAR random effect of eblup_sfh(),
# theta = X beta + (I - rho W)^{-1} eps.

# Simulated SAR data on a symmetric 4-nearest-neighbour graph (row-standardised W);
# `n_miss` domains have no direct estimate.
sim_sar <- function(D = 60, rho = 0.6, sigma = 0.7, seed = 1, n_miss = 5) {
  set.seed(seed)
  xy <- cbind(runif(D), runif(D))
  dd <- as.matrix(stats::dist(xy))
  diag(dd) <- Inf
  A <- matrix(0, D, D)
  for (i in seq_len(D)) A[i, order(dd[i, ])[1:4]] <- 1
  A <- pmax(A, t(A))
  W <- A / rowSums(A)
  x <- rnorm(D)
  u <- solve(diag(D) - rho * W, rnorm(D, sd = sigma))
  vardir <- runif(D, 0.1, 0.4)
  dat <- data.frame(area = seq_len(D), x = x, vardir = vardir)
  dat$theta <- 1 + 0.8 * x + u
  dat$y <- dat$theta + rnorm(D, sd = sqrt(vardir))
  dat$y[sample(D, n_miss)] <- NA
  list(data = dat, W = W, A = A)
}

# Exact posterior of theta under the INLA SAR model with known vardir: a grid over
# (log sigma, rho) with beta ~ N(0, prec_beta^-1 I) integrated analytically;
# sigma ~ Exponential(-log(0.01)) (= pc.prec(1, 0.01)), rho ~ U(-0.999, 0.999).
exact_sar_posterior <- function(dat, W, prec_beta = 1e-4, n_grid = 70) {
  D <- nrow(dat)
  X <- cbind(1, dat$x)
  s <- !is.na(dat$y)
  ys <- dat$y[s]
  XSX <- X %*% t(X) / prec_beta
  lam <- -log(0.01)
  lsig <- seq(log(0.05), log(3), length.out = n_grid)
  rhos <- seq(-0.995, 0.995, length.out = n_grid)
  out <- list()
  lw <- numeric(0)
  for (ls in lsig) {
    for (r in rhos) {
      Ainv <- solve(diag(D) - r * W)
      Ct <- XSX + exp(2 * ls) * Ainv %*% t(Ainv) # Cov(theta)
      Vy <- Ct[s, s] + diag(dat$vardir[s])
      R <- chol(Vy)
      a <- backsolve(R, ys, transpose = TRUE)
      ll <- -sum(log(diag(R))) - 0.5 * sum(a^2)
      lw <- c(lw, ll + log(lam) - lam * exp(ls) + ls) # prior on log sigma
      K <- Ct[, s] %*% chol2inv(R)
      out[[length(out) + 1]] <- list(m = drop(K %*% ys), v = diag(Ct - K %*% Ct[s, ]), rho = r)
    }
  }
  w <- exp(lw - max(lw))
  w <- w / sum(w)
  m <- Reduce(`+`, Map(function(o, wi) wi * o$m, out, w))
  m2 <- Reduce(`+`, Map(function(o, wi) wi * (o$v + o$m^2), out, w))
  list(mean = m, sd = sqrt(m2 - m^2), rho = sum(w * vapply(out, `[[`, 0, "rho")))
}

test_that("spatial = 'sar' validates its inputs", {
  sim <- sim_sar(D = 12)
  fit <- function(...) {
    hb_area(y ~ x,
      data = sim$data, domain = "area", vardir = "vardir", spatial = "sar",
      print_result = FALSE, ...
    )
  }
  expect_error(fit(), "must be provided")
  panel <- rbind(transform(sim$data, year = 1), transform(sim$data, year = 2))
  for (st in c("type3", "type4", "separable")) {
    expect_error(
      hb_area(y ~ x,
        data = panel, domain = "area", time = "year", vardir = "vardir",
        spatial = "sar", W = sim$W, temporal = "ar1", st_interaction = st, print_result = FALSE
      ),
      "sar"
    )
  }
})

test_that("spatial = 'sar' passes the row-standardised W to the JAX backend", {
  env <- new.env()
  local_mocked_bindings(
    .get_numpyro_backend = function(...) {
      list(fit_numpyro_hb = function(y, X, domain_idx, ...) {
        args <- list(...)
        env$args <- args
        n <- nrow(X)
        draws <- matrix(rnorm(100 * n, 1, 0.1), 100, n)
        list(
          beta_mean = rep(0.5, ncol(X)), beta_sd = rep(0.1, ncol(X)),
          beta_ci_lower = rep(0.3, ncol(X)), beta_ci_upper = rep(0.7, ncol(X)),
          hb_mean = colMeans(draws), hb_sd = apply(draws, 2, sd),
          hb_ci_lower = apply(draws, 2, quantile, 0.025, names = FALSE),
          hb_ci_upper = apply(draws, 2, quantile, 0.975, names = FALSE),
          linpred_mean = colMeans(draws), rand_eff_mean = rep(0, n),
          hb_draws = draws, hyperparameters = list(sigma2_u = 0.3, rho_spatial = 0.55),
          rhat_max = 1.001, dic = 1, p_dic = 1, waic = 1, p_waic = 1, device_used = "cpu",
          shrinkage_weights = NULL
        )
      })
    }
  )
  sim <- sim_sar(D = 12)
  fit <- hb_area(y ~ x,
    data = sim$data, domain = "area", vardir = "vardir", spatial = "sar",
    W = sim$A, engine = jax_engine(samples = 100, device = "cpu"), print_result = FALSE
  )
  expect_equal(env$args$spatial, "sar")
  expect_equal(unname(env$args$W_sar), unname(sim$W)) # binary A row-standardised
  expect_equal(fit$spatial, "sar")
  expect_equal(fit$rho, 0.55)
  expect_equal(fit$random_effect_var, 0.3)
})

test_that("spatial = 'sar' with INLA matches the exact posterior", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  sim <- sim_sar(D = 40)
  fit <- hb_area(y ~ x,
    data = sim$data, domain = "area", vardir = "vardir", spatial = "sar", W = sim$W,
    engine = inla_engine(control.fixed = list(
      mean = 0, prec = 1e-4, mean.intercept = 0, prec.intercept = 1e-4
    )),
    print_result = FALSE
  )
  ex <- exact_sar_posterior(sim$data, sim$W)
  z <- abs(fit$df_hb$hb - ex$mean) / ex$sd
  expect_lt(max(z), 0.1)
  expect_gt(min(fit$df_hb$sd / ex$sd), 0.9)
  expect_lt(max(fit$df_hb$sd / ex$sd), 1.1)
  expect_lt(abs(fit$rho - ex$rho), 0.05)
})

test_that("spatial = 'sar' is close to eblup_sfh() with many domains", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  sim <- sim_sar(D = 150, seed = 2, n_miss = 10)
  hb <- hb_area(y ~ x,
    data = sim$data, domain = "area", vardir = "vardir", spatial = "sar",
    W = sim$W, print_result = FALSE
  )
  sfh <- eblup_sfh(y ~ x, vardir = ~vardir, data = sim$data, W = sim$W, print_result = FALSE)
  z <- abs(hb$df_hb$hb - sfh$df_eblup$eblup) / hb$df_hb$sd
  expect_lt(median(z), 0.1)
  expect_lt(max(z), 0.3)
  r <- hb$df_hb$sd / sqrt(sfh$df_eblup$mse)
  expect_gt(median(r), 0.85)
  expect_lt(median(r), 1.2)
  # unsampled domains: kriging from the neighbours in both
  miss <- is.na(sim$data$y)
  expect_lt(max(z[miss]), 0.3)
})

test_that("spatial = 'sar' with JAX agrees with INLA", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("reticulate")
  if (!suppressMessages(check_numpyro_available())) skip("NumPyro / JAX not available")
  sim <- sim_sar(D = 60, seed = 3)
  fj <- suppressMessages(hb_area(y ~ x,
    data = sim$data, domain = "area", vardir = "vardir", spatial = "sar", W = sim$W,
    engine = jax_engine(warmup = 500, samples = 1000, chains = 4, device = "cpu", seed = 1),
    print_result = FALSE
  ))
  fi <- hb_area(y ~ x,
    data = sim$data, domain = "area", vardir = "vardir", spatial = "sar", W = sim$W,
    print_result = FALSE
  )
  expect_gt(cor(fj$df_hb$hb, fi$df_hb$hb), 0.99)
  # the default priors differ (half-normal SD and N(0, 10^2) coefficients on the
  # standardised scale against PC priors): 0.22 SD at most locally
  expect_lt(max(abs(fj$df_hb$hb - fi$df_hb$hb) / fi$df_hb$sd), 0.3)
  expect_lt(fj$rhat_max, 1.02)
  expect_lte(fj$divergences, 10)
})
