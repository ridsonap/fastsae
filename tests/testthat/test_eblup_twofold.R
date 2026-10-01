test_that("eblup_twofold recovers two-fold parameters from simulated data", {
  set.seed(42)
  m <- 40
  s2v_true <- 1.0  # area variance (paper notation)
  s2u_true <- 0.5  # subarea variance (paper notation)
  beta_true <- c(2, 1.5)

  dat <- do.call(rbind, lapply(seq_len(m), function(d) {
    nd <- sample(2:6, 1)
    data.frame(
      area = d,
      subarea = paste0(d, "-", seq_len(nd)),
      x1 = rnorm(nd),
      vardir = runif(nd, 0.2, 1.0)
    )
  }))
  v <- rnorm(m, sd = sqrt(s2v_true))[dat$area]  # area effect
  u <- rnorm(nrow(dat), sd = sqrt(s2u_true))    # subarea effect
  e <- rnorm(nrow(dat), sd = sqrt(dat$vardir))
  dat$y <- beta_true[1] + beta_true[2] * dat$x1 + v + u + e
  theta <- beta_true[1] + beta_true[2] * dat$x1 + v + u

  fit <- eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                       data = dat, print_result = FALSE)

  expect_true(fit$convergence)
  expect_equal(class(fit), "fastsae")
  expect_true(fit$model %in% c("TWOFOLD", "TFH"))

  s2 <- fit$random_effect_var
  expect_named(s2, c("sigma2_v", "sigma2_u"))
  # loose tolerance: variance components are noisy with m = 40
  expect_lt(abs(s2[["sigma2_v"]] - s2v_true), 0.7)
  expect_lt(abs(s2[["sigma2_u"]] - s2u_true), 0.35)
  expect_lt(max(abs(fit$estcoef$beta - beta_true)), 0.3)

  eb <- fit$df_eblup$eblup
  expect_gt(cor(eb, theta), 0.9)
  expect_true(all(fit$df_eblup$mse > 0))
  expect_true(all(is.finite(fit$df_eblup$mse)))
  # analytical MSE roughly calibrated: mean squared error / mean mse ~ O(1)
  ratio <- mean((eb - theta)^2) / mean(fit$df_eblup$mse)
  expect_gt(ratio, 0.2)
  expect_lt(ratio, 5)
})

test_that("eblup_twofold EBLUP matches a pure-R reference given the same parameters", {
  set.seed(7)
  m <- 10
  dat <- do.call(rbind, lapply(seq_len(m), function(d) {
    nd <- sample(2:4, 1)
    data.frame(area = d, subarea = seq_len(nd),
               x1 = rnorm(nd), vardir = runif(nd, 0.3, 0.8))
  }))
  dat$y <- 1 + dat$x1 + rnorm(m)[dat$area] + rnorm(nrow(dat), sd = 0.4) +
    rnorm(nrow(dat), sd = sqrt(dat$vardir))

  fit <- eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                       data = dat, print_result = FALSE)

  s2v <- fit$random_effect_var[["sigma2_v"]]  # area
  s2u <- fit$random_effect_var[["sigma2_u"]]  # subarea
  beta <- fit$estcoef$beta
  X <- cbind(1, dat$x1)

  # plain-R BLUP per area with dense solve (no Woodbury)
  ref <- numeric(nrow(dat))
  for (d in unique(dat$area)) {
    ii <- which(dat$area == d)
    nd <- length(ii)
    Vd <- s2v * matrix(1, nd, nd) + diag(s2u + dat$vardir[ii])
    Vi <- solve(Vd)
    Bd <- (s2v * matrix(1, nd, nd) + s2u * diag(nd)) %*% Vi
    ref[ii] <- X[ii, , drop = FALSE] %*% beta + Bd %*% (dat$y[ii] - X[ii, , drop = FALSE] %*% beta)
  }
  expect_equal(unname(fit$df_eblup$eblup), unname(ref), tolerance = 1e-6)
})

test_that("eblup_twofold handles non-sampled subareas and degenerate variances", {
  set.seed(11)
  m <- 15
  dat <- do.call(rbind, lapply(seq_len(m), function(d) {
    nd <- sample(1:4, 1)
    data.frame(area = d, subarea = seq_len(nd),
               x1 = rnorm(nd), vardir = runif(nd, 0.3, 0.8))
  }))
  # no subarea random effect in the truth -> sigma2_v should be ~ 0
  dat$y <- 1 + dat$x1 + rnorm(m, sd = 1)[dat$area] +
    rnorm(nrow(dat), sd = sqrt(dat$vardir))
  dat$y[c(3, 10)] <- NA_real_

  fit <- eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                       data = dat, print_result = FALSE)

  expect_true(fit$convergence)
  expect_lt(fit$random_effect_var[["sigma2_u"]], 0.2)
  expect_true(all(is.finite(fit$df_eblup$eblup)))
  expect_true(all(is.finite(fit$df_eblup$mse)))
  # non-sampled rows get synthetic predictions
  expect_true(all(!is.na(fit$df_eblup$eblup[c(3, 10)])))
})

test_that("eblup_twofold bootstrap MSE matches analytical MSE", {
  set.seed(42)
  m <- 20
  dat <- do.call(rbind, lapply(seq_len(m), function(d) {
    nd <- sample(2:4, 1)
    data.frame(area = d, subarea = paste0(d, "-", seq_len(nd)),
               x1 = rnorm(nd), vardir = runif(nd, 0.3, 1.0))
  }))
  dat$y <- 1 + dat$x1 + rnorm(m, sd = 1)[dat$area] +
    rnorm(nrow(dat), sd = 0.7) + rnorm(nrow(dat), sd = sqrt(dat$vardir))

  fit_a <- eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                         data = dat, mse = "analytical", print_result = FALSE)
  fit_b <- eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                         data = dat, mse = "bootstrap", B = 100, seed = 123,
                         print_result = FALSE)

  # EBLUP identical regardless of MSE method
  expect_equal(fit_a$df_eblup$eblup, fit_b$df_eblup$eblup)
  # MSE positively correlated and same order of magnitude
  expect_gt(cor(fit_a$df_eblup$mse, fit_b$df_eblup$mse), 0.7)
  ratio <- mean(fit_b$df_eblup$mse) / mean(fit_a$df_eblup$mse)
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2)
  expect_true(all(fit_b$df_eblup$mse > 0))
})

test_that("eblup_twofold requires a domain identifier", {
  expect_error(
    eblup_twofold(y ~ x1, vardir = "vardir", data = data.frame(y = 1, x1 = 1, vardir = 1)),
    "domain"
  )
})

test_that("eblup_twofold REML variance components match pure-R profile REML optimization", {
  set.seed(42)
  m <- 20
  dat <- do.call(rbind, lapply(seq_len(m), function(d) {
    nd <- sample(2:4, 1)
    data.frame(
      area = d, subarea = paste0(d, "-", seq_len(nd)),
      x1 = rnorm(nd), vardir = runif(nd, 0.2, 0.8)
    )
  }))
  dat$y <- 1 + dat$x1 + rnorm(m)[dat$area] + rnorm(nrow(dat), sd = 0.5) + rnorm(nrow(dat), sd = sqrt(dat$vardir))

  fit_tfh <- eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                           data = dat, method = "REML", print_result = FALSE)

  # Exact profile REML log-likelihood
  areas <- unique(dat$area)
  X_list <- lapply(areas, function(d) cbind(1, dat$x1[dat$area == d]))
  y_list <- lapply(areas, function(d) dat$y[dat$area == d])
  psi_list <- lapply(areas, function(d) dat$vardir[dat$area == d])

  neg_loglik_reml <- function(theta) {
    s2v <- exp(theta[1])
    s2u <- exp(theta[2])

    sum_Xt_Vi_X <- 0
    sum_Xt_Vi_y <- 0
    log_det_V <- 0

    Vi_list <- list()
    for (i in seq_along(areas)) {
      ni <- length(y_list[[i]])
      Vi <- s2v * matrix(1, ni, ni) + diag(s2u + psi_list[[i]], ni, ni)
      chol_Vi <- chol(Vi)
      log_det_V <- log_det_V + 2 * sum(log(diag(chol_Vi)))
      inv_Vi <- chol2inv(chol_Vi)
      Vi_list[[i]] <- inv_Vi

      Xi <- X_list[[i]]
      yi <- y_list[[i]]
      sum_Xt_Vi_X <- sum_Xt_Vi_X + t(Xi) %*% inv_Vi %*% Xi
      sum_Xt_Vi_y <- sum_Xt_Vi_y + t(Xi) %*% inv_Vi %*% yi
    }

    chol_Xt_Vi_X <- chol(sum_Xt_Vi_X)
    log_det_Xt_Vi_X <- 2 * sum(log(diag(chol_Xt_Vi_X)))
    beta_hat <- chol2inv(chol_Xt_Vi_X) %*% sum_Xt_Vi_y

    quad_form <- 0
    for (i in seq_along(areas)) {
      ri <- y_list[[i]] - X_list[[i]] %*% beta_hat
      quad_form <- quad_form + t(ri) %*% Vi_list[[i]] %*% ri
    }

    val <- 0.5 * (log_det_V + log_det_Xt_Vi_X + quad_form)
    return(as.numeric(val))
  }

  opt <- stats::optim(c(0, 0), neg_loglik_reml, method = "BFGS")
  s2v_opt <- exp(opt$par[1])
  s2u_opt <- exp(opt$par[2])

  expect_equal(unname(fit_tfh$random_effect_var["sigma2_v"]), s2v_opt, tolerance = 1e-4)
  expect_equal(unname(fit_tfh$random_effect_var["sigma2_u"]), s2u_opt, tolerance = 1e-4)
})

test_that("eblup_twofold matches Bayesian two-fold hierarchical model via INLA", {
  skip_if_not_installed("INLA")

  set.seed(42)
  m <- 20
  dat <- do.call(rbind, lapply(seq_len(m), function(d) {
    nd <- sample(2:4, 1)
    data.frame(
      area = d, subarea = paste0(d, "-", seq_len(nd)),
      x1 = rnorm(nd), vardir = runif(nd, 0.2, 0.8)
    )
  }))
  dat$y <- 1 + dat$x1 + rnorm(m)[dat$area] + rnorm(nrow(dat), sd = 0.5) + rnorm(nrow(dat), sd = sqrt(dat$vardir))
  dat$sub_id <- seq_len(nrow(dat))

  fit_tfh <- eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                           data = dat, print_result = FALSE)

  form_inla <- y ~ x1 + f(area, model = "iid") + f(sub_id, model = "iid")
  prec_obs <- 1 / dat$vardir
  fit_inla <- INLA::inla(
    form_inla,
    data = dat,
    family = "gaussian",
    control.family = list(hyper = list(prec = list(initial = 0, fixed = TRUE))),
    scale = prec_obs,
    control.predictor = list(compute = TRUE),
    num.threads = 1
  )

  eb_pred <- fit_tfh$df_eblup$eblup
  inla_pred <- fit_inla$summary.fitted.values[seq_len(nrow(dat)), "mean"]

  expect_gt(stats::cor(eb_pred, inla_pred), 0.95)
})

test_that("eblup_twofold is consistent with standard Fay-Herriot on single-subarea boundary case", {
  set.seed(123)
  n <- 30
  dat <- data.frame(
    area = seq_len(n),
    subarea = paste0(seq_len(n), "-1"),
    x1 = rnorm(n),
    vardir = runif(n, 0.2, 0.8)
  )
  dat$y <- 2 + 1.5 * dat$x1 + rnorm(n, sd = 0.6) + rnorm(n, sd = sqrt(dat$vardir))

  fit_fh <- eblup_fh(y ~ x1, vardir = "vardir", data = dat, print_result = FALSE)
  fit_tfh <- eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                           data = dat, print_result = FALSE)

  expect_gt(stats::cor(fit_fh$df_eblup$eblup, fit_tfh$df_eblup$eblup), 0.97)
  expect_lt(max(abs(fit_fh$estcoef$beta - fit_tfh$estcoef$beta)), 0.05)
})

test_that("eblup_tfh alias is identical to eblup_twofold", {
  expect_identical(eblup_tfh, eblup_twofold)
})


