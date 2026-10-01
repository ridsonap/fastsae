test_that("eblup_tfh recovers two-fold parameters from simulated data", {
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

  fit <- eblup_tfh(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                   data = dat, print_result = FALSE)

  expect_true(fit$convergence)
  expect_equal(class(fit), "fastsae")
  expect_equal(fit$model, "TFH")

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

test_that("eblup_tfh EBLUP matches a pure-R reference given the same parameters", {
  set.seed(7)
  m <- 10
  dat <- do.call(rbind, lapply(seq_len(m), function(d) {
    nd <- sample(2:4, 1)
    data.frame(area = d, subarea = seq_len(nd),
               x1 = rnorm(nd), vardir = runif(nd, 0.3, 0.8))
  }))
  dat$y <- 1 + dat$x1 + rnorm(m)[dat$area] + rnorm(nrow(dat), sd = 0.4) +
    rnorm(nrow(dat), sd = sqrt(dat$vardir))

  fit <- eblup_tfh(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
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

test_that("eblup_tfh handles non-sampled subareas and degenerate variances", {
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

  fit <- eblup_tfh(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                   data = dat, print_result = FALSE)

  expect_true(fit$convergence)
  expect_lt(fit$random_effect_var[["sigma2_u"]], 0.2)
  expect_true(all(is.finite(fit$df_eblup$eblup)))
  expect_true(all(is.finite(fit$df_eblup$mse)))
  # non-sampled rows get synthetic predictions
  expect_true(all(!is.na(fit$df_eblup$eblup[c(3, 10)])))
})

test_that("eblup_tfh requires a domain identifier", {
  expect_error(
    eblup_tfh(y ~ x1, vardir = "vardir", data = data.frame(y = 1, x1 = 1, vardir = 1)),
    "domain"
  )
})
