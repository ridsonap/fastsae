# hb_area(var_x = ): covariates measured with error (structural and functional models).
# The data come from sim_me_df() (helper-measurement-error.R); the JAX fits end to end
# are in test_hb_area_engine.R, which the JAX job of the CI runs.

fit_me <- function(formula, data, ...) {
  hb_area(formula, data = data, domain = "area", vardir = "vardir", print_result = FALSE, ...)
}

# print() writes through cli (messages) and cat (stdout)
printed_me <- function(x) {
  out <- NULL
  msgs <- capture_messages(out <- utils::capture.output(print(x)))
  paste(c(msgs, out), collapse = "\n")
}

# Stand-in for hb_area_numpyro.fit_numpyro_hb() that records its arguments
fake_me_backend <- function(env) {
  list(fit_numpyro_hb = function(y, X, domain_idx, ...) {
    env$args <- c(list(y = y, X = X, domain_idx = domain_idx), list(...))
    n <- nrow(X)
    set.seed(3)
    draws <- matrix(stats::rnorm(100 * n, 1, 0.1), 100, n)
    list(
      beta_mean = rep(0.5, ncol(X)), beta_sd = rep(0.1, ncol(X)),
      beta_ci_lower = rep(0.3, ncol(X)), beta_ci_upper = rep(0.7, ncol(X)),
      hb_mean = colMeans(draws), hb_sd = apply(draws, 2, stats::sd),
      hb_ci_lower = apply(draws, 2, stats::quantile, 0.025, names = FALSE),
      hb_ci_upper = apply(draws, 2, stats::quantile, 0.975, names = FALSE),
      linpred_mean = colMeans(draws), rand_eff_mean = rep(0, n),
      hyperparameters = stats::setNames(list(0.2, 5, 4), c("sigma2_u", "mu_x[1]", "sigma2_x[1]")),
      dic = 10, p_dic = 3, waic = 11, p_waic = 3.5,
      device_used = "cpu", rhat_max = 1.01, shrinkage_weights = NULL, hb_draws = draws
    )
  })
}

test_that("var_x is checked before fitting", {
  d <- sim_me_df()
  expect_error(fit_me(y ~ w, d, var_x = "vx"), "named")
  expect_error(fit_me(y ~ w, d, var_x = c(w = "vx", w = "vx")), "named")
  expect_error(fit_me(y ~ w, d, var_x = 1:3), "named")
  expect_error(fit_me(y ~ w, d, var_x = c(z = "vx")), "not a main effect")
  expect_error(fit_me(y ~ log(w), d, var_x = c(w = "vx")), "not a main effect")
  expect_error(fit_me(y ~ w + w:z, d, var_x = c(w = "vx")), "also enters")
  expect_error(fit_me(y ~ w + I(w^2), d, var_x = c(w = "vx")), "also enters")
  d$wf <- factor(d$w > 5)
  expect_error(fit_me(y ~ wf, d, var_x = c(wf = "vx")), "numeric")
  expect_error(fit_me(y ~ w, d, var_x = c(w = "nope")), "not valid")
  expect_error(fit_me(y ~ w, d, var_x = list(w = 1:3)), "not valid")
  expect_error(fit_me(y ~ w, d, var_x = c(w = "area")), "must be numeric")
  d$vbad <- d$vx
  d$vbad[3] <- -1
  expect_error(fit_me(y ~ w, d, var_x = c(w = "vbad")), "non-negative")
  d$vbad[3] <- NA
  expect_error(fit_me(y ~ w, d, var_x = c(w = "vbad")), "non-negative")
  expect_error(fit_me(y ~ w, d, var_x = c(w = "vx"), me_model = "berkson"), "should be one of")
  expect_error(fit_me(y ~ w, d, var_x = c(w = "vx"), spatial = "slm", W = diag(nrow(d))), "slm")
  # without an intercept a constant covariate passes the design check
  d0 <- transform(d, w = 2)
  for (m in c("structural", "functional")) {
    expect_error(fit_me(y ~ 0 + w, d0, var_x = c(w = "vx"), me_model = m), "same value in every row")
  }
  # me_model alone does nothing; the input error below stops before fitting
  expect_warning(
    expect_error(hb_area(y ~ w, data = "not a data frame", me_model = "functional"), "data frame"),
    "ignored without"
  )
  expect_no_warning(expect_error(hb_area(y ~ w, data = "not a data frame"), "data frame"))
})

test_that(".me_spec() shares a true value only within a domain", {
  X <- cbind("(Intercept)" = 1, w = c(1, 1, 2, 1, 1))
  spec <- .me_spec(
    list(w = c(0.1, 0.1, 0.1, 0.1, 0.2)), "functional", y ~ w, data.frame(y = 1:5), X,
    domain = c("a", "b", "a", "a", "a")
  )
  # rows 1 and 4: same domain, value and variance; row 5 has another variance
  expect_equal(unname(spec$group[, "w"]), c(1L, 2L, 3L, 1L, 4L))
  expect_equal(spec$model, "functional")
  expect_equal(unname(spec$var[, "w"]), c(0.1, 0.1, 0.1, 0.1, 0.2))
})

test_that(".me_jitter_ties() makes repeated values distinct and keeps the others", {
  x <- c(0, 1, 1, 1e7, 1e7, 2, 1)
  j <- .me_jitter_ties(x)
  expect_equal(anyDuplicated(j), 0L)
  expect_equal(j, x, tolerance = 1e-12)
  expect_identical(j[c(1, 2, 4, 6)], x[c(1, 2, 4, 6)])
})

test_that(".me_inla_terms() sets up the mec terms", {
  me <- list(
    model = "structural", cov = "w", w = cbind(w = c(1, 2, 2, 4)),
    var = cbind(w = c(0, 0.1, 0.1, 0.2)), group = cbind(w = 1:4)
  )
  X <- cbind("(Intercept)" = 1, w = me$w[, 1])
  tm <- .me_inla_terms(me, y = c(1, 1, 1, NA), X = X, family = "gaussian")
  t1 <- tm$terms[[1]]
  # areas 2 and 3 have the same observed value but their own true values
  expect_equal(anyDuplicated(t1$values), 0L)
  expect_equal(t1$x, t1$values)
  # a variance of 0 is floored at 1e-8 var(w)
  expect_equal(max(t1$scale), 1e8)
  expect_true(t1$hyper$mean.x$fixed && t1$hyper$prec.x$fixed)
  # y does not depend on w: the initial coefficient is kept away from 0
  expect_gt(abs(t1$hyper$beta$initial), 0)
  expect_equal(tm$coef_order, c("(Intercept)", "w"))
  expect_equal(t1$sigma2_x, unname(.me_x_mode(t1$values, 1 / t1$scale)["sigma2"]) * t1$s^2, tolerance = 1e-6)
})

test_that("measurement error helpers handle every family and form of prior", {
  X <- cbind("(Intercept)" = 1, w = 1:6)
  y <- c(2, 4.1, 5.9, 8.2, 9.9, NA)
  b <- .me_initial_beta(X, y, "gaussian")
  expect_equal(unname(b["w"]), unname(stats::coef(stats::lm(y ~ X[, "w"]))[2]))
  expect_gt(attr(b, "scale"), 0)
  counts <- c(1, 3, 6, 10, 20, 31)
  expect_true(all(is.finite(.me_initial_beta(X, c(0, 2, 5, 9, 14, 20), "binomial", trials = rep(20, 6)))))
  expect_true(all(is.finite(.me_initial_beta(X, counts, "poisson", exposure = rep(2, 6)))))
  expect_true(all(is.finite(.me_initial_beta(X, counts, "nbinomial"))))
  expect_true(all(is.finite(.me_initial_beta(X, c(0.1, 0.2, 0.3, 0.5, 0.7, 0.9), "beta"))))
  expect_true(all(is.finite(.me_initial_beta(X, c(1, 2, 4, 8, 16, 32), "gamma"))))

  expect_equal(.me_beta_prior(NULL, "w"), c(mean = 0, prec = 0.001))
  expect_equal(.me_beta_prior(list(mean = 1, prec = list(w = 0.5, default = 0.01)), "w"), c(mean = 1, prec = 0.5))
  expect_equal(.me_beta_prior(list(prec = list(default = 0.01)), "w"), c(mean = 0, prec = 0.01))

  eta <- c(-1, 0.5)
  for (lk in c("logit", "probit", "cloglog", "loglog")) {
    expect_equal(.inla_link(.inla_inv_link(eta, "binomial", lk), "binomial", lk), eta)
  }
  for (fam in c("gaussian", "poisson")) {
    expect_equal(.inla_link(.inla_inv_link(eta, fam), fam), eta)
  }
})

test_that(".me_rename_mcmc() names the true-value parameters after the covariates", {
  tab <- data.frame(parameter = c("sigma2_u", "mu_x[1]", "sigma2_x[2]", "mu_x[3]"))
  tr <- array(0, c(2, 1, 2), dimnames = list(iteration = NULL, chain = "1", parameter = c("mu_x[1]", "beta[1]")))
  out <- .me_rename_mcmc(list(table = tab, trace = tr), c("a", "b"))
  expect_equal(out$table$parameter, c("sigma2_u", "mu_x[a]", "sigma2_x[b]", "mu_x[3]"))
  expect_equal(dimnames(out$trace)$parameter, c("mu_x[a]", "beta[1]"))
  expect_null(.me_rename_mcmc(list(table = NULL, trace = NULL), "a")$table)
})

test_that("measurement error hyperparameters count against the grid of interaction models", {
  expect_equal(.me_n_hyper(NULL), 0L)
  expect_equal(.me_n_hyper(list(cov = c("a", "b"), model = "structural")), 2L)
  expect_equal(.me_n_hyper(list(cov = "a", model = "functional")), 1L)
  expect_equal(.n_hyper_st("none", "rw1", "type1", vardir_known = TRUE), 3L)
  expect_equal(.n_hyper_st("none", "rw1", "type1", vardir_known = TRUE, n_me = 1L), 4L)
})

test_that("zero error variances reproduce the fit with the covariate as a fixed effect", {
  skip_if_not_installed("INLA")
  d <- sim_me_df(n_unsampled = 5)
  d$v0 <- 0
  fit_fe <- fit_me(y ~ w + z, d)
  for (model in c("structural", "functional")) {
    fit_0 <- fit_me(y ~ w + z, d, var_x = c(w = "v0"), me_model = model)
    expect_lt(max(abs(fit_0$df_hb$hb - fit_fe$df_hb$hb) / fit_fe$df_hb$sd), 0.05)
    expect_lt(max(abs(fit_0$df_hb$sd / fit_fe$df_hb$sd - 1)), 0.02)
    expect_equal(fit_0$estcoef$beta, fit_fe$estcoef$beta, tolerance = 0.02)
    expect_equal(fit_0$estcoef$std.error, fit_fe$estcoef$std.error, tolerance = 0.03)
  }
})

test_that("a fit with measurement error keeps the structure of hb_area() results", {
  skip_if_not_installed("INLA")
  d <- sim_me_df(n_unsampled = 5)
  fit <- fit_me(y ~ z + w, d, var_x = c(w = "vx"))
  fit_fe <- fit_me(y ~ z + w, d)

  expect_s3_class(fit, c("fastsae_hb_area", "fastsae_hb", "fastsae"))
  expect_equal(rownames(fit$estcoef), c("(Intercept)", "z", "w"))
  expect_named(coef(fit), c("(Intercept)", "z", "w"))
  expect_identical(names(fit$estcoef), names(fit_fe$estcoef))
  expect_identical(names(fit$df_hb), names(fit_fe$df_hb))
  expect_false(anyNA(fit$df_hb$hb))
  expect_true(fit$convergence)

  me <- fit$measurement_error
  expect_equal(me$model, "structural")
  expect_equal(me$covariates, "w")
  # the mode given the observed values: mean, and variance net of the errors
  expect_equal(unname(me$mu_x), mean(d$w), tolerance = 0.05)
  expect_equal(unname(me$sigma2_x), var(d$w) - mean(d$vx), tolerance = 0.2)
  expect_null(fit_fe$measurement_error)
  expect_true(fit$fit$.args$inla.mode == "classic")

  # the coefficient is a hyperparameter of the term, in the unit of the covariate
  expect_true("MEC beta for w" %in% rownames(fit$hyperpar))
  expect_false(any(grepl("mean_x|prec_x|\\.\\.me_", rownames(fit$hyperpar))))
  expect_equal(fit$hyperpar["MEC beta for w", "mean"], fit$estcoef["w", "beta"])

  expect_match(fit$model, "[ME: w]", fixed = TRUE)
  expect_match(printed_me(fit), "Measurement error (structural): w", fixed = TRUE)
  expect_match(printed_me(summary(fit)), "Measurement error (structural): w", fixed = TRUE)

  fit_f <- fit_me(y ~ z + w, d, var_x = list(w = ~vx), me_model = "functional")
  expect_equal(fit_f$measurement_error$model, "functional")
  expect_null(fit_f$measurement_error$mu_x)
  expect_false(any(grepl("mean_x|prec_x", rownames(fit_f$hyperpar))))
})

test_that("the posterior SD of a domain without a direct estimate grows with its error variance", {
  skip_if_not_installed("INLA")
  d <- sim_me_df(n_unsampled = 6)
  d$vx[1:3] <- 0.05
  d$vx[4:6] <- 3
  fit <- fit_me(y ~ w, d, var_x = c(w = "vx"))
  fit_fe <- fit_me(y ~ w, d)
  expect_gt(mean(fit$df_hb$sd[4:6]) / mean(fit$df_hb$sd[1:3]), 1.2)
  # without the measurement error model the variance does not enter
  expect_lt(mean(fit_fe$df_hb$sd[4:6]) / mean(fit_fe$df_hb$sd[1:3]), 1.1)
})

test_that("the structural model recovers the coefficient that measurement error attenuates", {
  skip_if_not_installed("INLA")
  d <- sim_me_df(D = 150, seed = 3)
  fit_s <- fit_me(y ~ w, d, var_x = c(w = "vx"))
  fit_f <- fit_me(y ~ w, d, var_x = c(w = "vx"), me_model = "functional")
  fit_n <- fit_me(y ~ w, d)
  b <- function(fit) fit$estcoef["w", "beta"]
  expect_lt(abs(b(fit_s) - 0.8), abs(b(fit_n) - 0.8))
  expect_lt(fit_s$estcoef["w", "ci_lower"], 0.8)
  expect_gt(fit_s$estcoef["w", "ci_upper"], 0.8)
  # the flat prior on the true values gives about the naive coefficient
  expect_lt(abs(b(fit_f) - b(fit_n)), abs(b(fit_s) - b(fit_n)))
  expect_equal(unname(fit_s$measurement_error$sigma2_x), 4, tolerance = 0.35)
})

test_that("the functional model agrees with the Ybarra-Lohr EBLUP of emdi", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("emdi")
  d <- sim_me_df(D = 80)
  fit <- fit_me(y ~ w, d, var_x = c(w = "vx"), me_model = "functional")
  Ci <- array(0, c(2, 2, nrow(d)))
  Ci[2, 2, ] <- d$vx
  yl <- suppressMessages(suppressWarnings(emdi::fh(
    fixed = y ~ w, vardir = "vardir", combined_data = d, domains = "area",
    method = "me", Ci = Ci, MSE = TRUE, mse_type = "jackknife"
  )))
  i <- match(d$area, as.character(yl$ind$Domain))
  expect_gt(cor(fit$df_hb$hb, yl$ind$FH[i]), 0.99)
  expect_gt(cor(fit$df_hb$mse, yl$MSE$FH[i]), 0.95)
})

test_that("domains have their own true values, shared across periods only for a reused estimate", {
  skip_if_not_installed("INLA")
  d <- sim_me_df(D = 30)
  d$w[2] <- d$w[1]
  d$vx[2] <- d$vx[1]
  fit <- fit_me(y ~ w, d, var_x = c(w = "vx"))
  expect_equal(nrow(fit$fit$summary.random$..me_1..), 30)

  dp <- rbind(transform(d, year = 1), transform(d, year = 2, y = y + 0.2))
  dp$w[dp$year == 2 & dp$area == "a3"] <- d$w[3] + 0.5
  fit_p <- fit_me(y ~ w, dp, time = "year", temporal = "iid", var_x = c(w = "vx"))
  expect_equal(nrow(fit_p$fit$summary.random$..me_1..), 31)
})

test_that("covariates in large units and covariates without an effect are fitted", {
  skip_if_not_installed("INLA")
  d <- sim_me_df()
  d$w_big <- d$w * 1e5
  d$vx_big <- d$vx * 1e10
  fit <- fit_me(y ~ w, d, var_x = c(w = "vx"))
  fit_big <- fit_me(y ~ w_big, d, var_x = c(w_big = "vx_big"))
  expect_equal(fit_big$estcoef["w_big", "beta"] * 1e5, fit$estcoef["w", "beta"], tolerance = 0.02)
  expect_equal(fit_big$df_hb$hb, fit$df_hb$hb, tolerance = 0.01)

  d0 <- sim_me_df(beta = 0)
  fit_0 <- fit_me(y ~ w, d0, var_x = c(w = "vx"))
  expect_true(fit_0$convergence)
  expect_lt(fit_0$estcoef["w", "ci_lower"], 0)
  expect_gt(fit_0$estcoef["w", "ci_upper"], 0)
})

test_that("measurement error works with a binomial likelihood", {
  skip_if_not_installed("INLA")
  set.seed(5)
  D <- 50
  x <- rnorm(D)
  vx <- runif(D, 0.1, 0.5)
  n <- rpois(D, 40) + 5
  d <- data.frame(
    area = seq_len(D), n = n, vx = vx, w = x + rnorm(D, 0, sqrt(vx)),
    y = rbinom(D, n, plogis(-1 + 0.7 * x + rnorm(D, 0, 0.3)))
  )
  d$y[1:5] <- NA
  fit <- hb_area(y ~ w,
    data = d, domain = "area", family = "binomial", trials = "n",
    var_x = c(w = "vx"), print_result = FALSE
  )
  expect_true(all(fit$df_hb$hb > 0 & fit$df_hb$hb < 1))
  expect_equal(rownames(fit$estcoef), c("(Intercept)", "w"))
  expect_gt(fit$estcoef["w", "beta"], 0)
})

test_that("self-benchmarking holds with measurement error", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("sn")
  d <- sim_me_df()
  d$pop <- seq_len(nrow(d))
  fit <- fit_me(y ~ w, d,
    var_x = c(w = "vx"), self_benchmark = TRUE, benchmark_weight = "pop",
    engine = inla_engine(draws = 200)
  )
  expect_equal(sum(fit$df_hb$hb * d$pop) / sum(d$pop), sum(d$y * d$pop) / sum(d$pop), tolerance = 1e-8)
  fit_wfq <- fit_me(y ~ w, d,
    var_x = c(w = "vx"), self_benchmark = TRUE, benchmark_weight = "pop", benchmark_method = "wfq"
  )
  expect_equal(rownames(fit_wfq$estcoef), c("(Intercept)", "w", "bench_aug"))
})

test_that("Gaussian interaction models with measurement error leave INLA's integration", {
  skip_if_not_installed("INLA")
  d <- sim_me_df(D = 20)
  dp <- rbind(transform(d, year = 1), transform(d, year = 2), transform(d, year = 3))
  dp$y <- dp$y + rnorm(nrow(dp), 0, 0.3)
  args <- list(time = "year", temporal = "rw1", st_interaction = "type1")
  fit_fe <- do.call(fit_me, c(list(y ~ w, dp), args))
  fit <- do.call(fit_me, c(list(y ~ w, dp), args, list(var_x = c(w = "vx"))))
  expect_identical(fit_fe$fit$.args$control.inla$int.strategy, "grid")
  expect_false(identical(fit$fit$.args$control.inla$int.strategy, "grid"))
})

test_that("diagnose(), benchmark_sae() and hb_ratio() work with measurement error", {
  skip_if_not_installed("INLA")
  d <- sim_me_df()
  fit <- fit_me(y ~ w, d, var_x = c(w = "vx"))
  expect_s3_class(suppressWarnings(diagnose(fit)), "fastsae_diagnose")
  bm <- suppressMessages(benchmark_sae(fit, target = 6, weight = seq_len(nrow(d))))
  expect_s3_class(bm, "fastsae_benchmark")
  expect_equal(sum(bm$benchmarked * bm$weight) / sum(bm$weight), 6)

  set.seed(2)
  D <- 40
  r <- data.frame(area = seq_len(D), x = rnorm(D), vx = 0.1, p0 = 0.2, n_poor = rpois(D, 40) + 5)
  r$s_dir <- rbinom(D, r$n_poor, plogis(-1 + 0.5 * r$x)) / r$n_poor
  fit_r <- hb_ratio(s_dir ~ x,
    data = r, trials = "n_poor", parent = "p0", domain = "area",
    var_x = c(x = "vx"), n_draws = 200, print_result = FALSE
  )
  expect_equal(fit_r$ratio_fit$measurement_error$covariates, "x")
})

test_that("engine = 'jax' passes the covariates measured with error to the backend", {
  env <- new.env()
  local_mocked_bindings(.get_numpyro_backend = function(...) fake_me_backend(env))
  d <- sim_me_df(D = 12)
  d$vx[2] <- 0
  dp <- rbind(transform(d, year = 1), transform(d, year = 2))
  fit <- hb_area(y ~ z + w,
    data = dp, domain = "area", time = "year", temporal = "iid", vardir = "vardir",
    var_x = c(w = "vx"), engine = jax_engine(samples = 10, device = "cpu"), print_result = FALSE
  )
  a <- env$args
  expect_equal(a$me_model, "structural")
  expect_equal(a$me_cols, list(2L))
  # one true value per domain: the same estimate in both periods
  expect_equal(a$me_group[[1]], rep(0:11, 2))
  expect_equal(a$me_w[[1]], d$w)
  expect_equal(a$me_var[[1]][-2], d$vx[-2])
  expect_equal(a$me_var[[1]][2], 1e-8 * stats::var(dp$w))
  expect_equal(
    fit$measurement_error,
    list(model = "structural", covariates = "w", mu_x = c(w = 5), sigma2_x = c(w = 4))
  )
  expect_equal(fit$hyperpar$Parameter, c("sigma2_u", "mu_x[w]", "sigma2_x[w]"))
  expect_match(fit$model, "[ME: w]", fixed = TRUE)

  fit_f <- hb_area(y ~ z + w,
    data = d, domain = "area", vardir = "vardir", var_x = c(w = "vx"), me_model = "functional",
    engine = jax_engine(samples = 10, device = "cpu"), print_result = FALSE
  )
  expect_equal(env$args$me_model, "functional")
  expect_equal(env$args$me_group[[1]], 0:11)
  expect_null(fit_f$measurement_error$mu_x)
})

test_that("step_sae() can drop a covariate measured with error", {
  skip_if_not_installed("INLA")
  d <- sim_me_df(D = 30)
  # step_sae() evaluates the stored call again: no wrapper
  fit <- hb_area(y ~ w + z,
    data = d, domain = "area", vardir = "vardir", var_x = c(w = "vx"),
    me_model = "functional", print_result = FALSE
  )
  # a criterion that prefers every model without w
  local_mocked_bindings(.sae_criterion = function(fit, criterion) if ("w" %in% rownames(fit$estcoef)) 1 else 0)
  expect_no_warning(sel <- step_sae(fit, direction = "backward", trace = FALSE))
  expect_equal(deparse(sel$final), "y ~ z")
  expect_null(sel$model$measurement_error)
  expect_null(sel$model$call$var_x)
})
