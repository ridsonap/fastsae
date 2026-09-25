test_that("ebp_area works for Gaussian Fay-Herriot model with INLA", {
  skip_if_not_installed("INLA")

  data("mys", package = "fastsae")

  # Fit standard Gaussian FH
  fit_norm <- ebp_area(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = "vardir",
    family = "gaussian",
    spatial = "none",
    print_result = FALSE
  )

  expect_s3_class(fit_norm, "fastsae_ebp_area")
  expect_s3_class(fit_norm, "fastsae")

  # Verify df_ebp columns
  expect_true(all(c("domain", "y", "ebp", "sd", "mse", "rse", "ci_lower", "ci_upper") %in% names(fit_norm$df_ebp)))
  expect_equal(nrow(fit_norm$df_ebp), nrow(mys))

  # Verify unsampled domains get predictions
  unsampled_idx <- which(is.na(mys$y))
  expect_true(length(unsampled_idx) > 0)
  expect_false(any(is.na(fit_norm$df_ebp$ebp[unsampled_idx])))
  expect_false(any(is.na(fit_norm$df_ebp$sd[unsampled_idx])))

  # Verify coefficients
  expect_true(nrow(fit_norm$estcoef) == 4)
  expect_true(all(c("beta", "std.error", "zvalue", "pvalue") %in% names(fit_norm$estcoef)))

  # Verify methods
  expect_equal(length(fitted(fit_norm)), nrow(mys))
  expect_equal(length(coef(fit_norm)), 4)
  expect_output(print(fit_norm))
  expect_output(print(summary(fit_norm)))
})

test_that("ebp_area works with BYM2 and Besag spatial priors", {
  skip_if_not_installed("INLA")

  data("mys", package = "fastsae")
  data("mys_proxmat", package = "fastsae")

  # 1. BYM2 via ebp_area
  fit_bym2 <- ebp_area(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = "vardir",
    family = "gaussian",
    spatial = "bym2",
    W = mys_proxmat,
    print_result = FALSE
  )

  expect_s3_class(fit_bym2, "fastsae_ebp_area")
  expect_false(is.null(fit_bym2$phi))
  expect_true(fit_bym2$phi >= 0 && fit_bym2$phi <= 1)

  # Check unsampled domains in BYM2
  unsampled_idx <- which(is.na(mys$y))
  expect_false(any(is.na(fit_bym2$df_ebp$ebp[unsampled_idx])))

  # 2. BYM2 with spatial weights W
  fit_spatial <- ebp_area(
    y ~ x1 + x2 + x3,
    data = mys,
    W = mys_proxmat,
    vardir = "vardir",
    family = "gaussian",
    spatial = "bym2",
    print_result = FALSE
  )
  expect_s3_class(fit_spatial, "fastsae_ebp_area")
  expect_equal(length(fit_spatial$df_ebp$ebp), nrow(mys))

  # 3. Besag ICAR model
  fit_besag <- ebp_area(
    y ~ x1 + x2 + x3,
    data = mys,
    W = mys_proxmat,
    vardir = "vardir",
    family = "gaussian",
    spatial = "besag",
    print_result = FALSE
  )
  expect_s3_class(fit_besag, "fastsae_ebp_area")
  expect_false(any(is.na(fit_besag$df_ebp$ebp)))
})

test_that("ebp_area works for Binomial response", {
  skip_if_not_installed("INLA")

  set.seed(42)
  D <- 25
  n_trials <- sample(30:100, D, replace = TRUE)
  x <- rnorm(D)
  prob_true <- plogis(-0.5 + 1.2 * x)
  y <- rbinom(D, size = n_trials, prob = prob_true)
  # Unsampled domains
  y[c(5, 12)] <- NA

  df_bin <- data.frame(id = 1:D, y = y, x = x, n = n_trials)

  # Non-spatial Binomial Logit
  fit_bin <- ebp_area(
    y ~ x,
    data = df_bin,
    domain = "id",
    family = "binomial",
    trials = "n",
    print_result = FALSE
  )

  expect_s3_class(fit_bin, "fastsae_ebp_area")
  # Probabilities should be strictly in (0, 1)
  expect_true(all(fit_bin$df_ebp$ebp >= 0 & fit_bin$df_ebp$ebp <= 1))
  # Unsampled domains predicted
  expect_false(any(is.na(fit_bin$df_ebp$ebp[c(5, 12)])))
  # Check estimated totals column
  expect_true("estimated_total" %in% names(fit_bin$df_ebp))
})

test_that("ebp_area works for Poisson response with exposure", {
  skip_if_not_installed("INLA")

  set.seed(42)
  D <- 20
  exposure <- runif(D, 100, 500)
  x <- rnorm(D)
  lambda_true <- exp(-1 + 0.5 * x)
  y <- rpois(D, lambda = exposure * lambda_true)
  y[3] <- NA

  df_pois <- data.frame(id = 1:D, y = y, x = x, E = exposure)

  # Non-spatial Poisson Log
  fit_pois <- ebp_area(
    y ~ x,
    data = df_pois,
    domain = "id",
    family = "poisson",
    exposure = "E",
    print_result = FALSE
  )

  expect_s3_class(fit_pois, "fastsae_ebp_area")
  # Rates should be positive
  expect_true(all(fit_pois$df_ebp$ebp > 0))
  expect_false(is.na(fit_pois$df_ebp$ebp[3]))
  expect_true("estimated_count" %in% names(fit_pois$df_ebp))
})

test_that("ebp_area works with method = 'laplace' (Frequentist GLMM)", {
  set.seed(42)
  D <- 20
  x <- rnorm(D)
  n <- rep(50, D)
  p <- plogis(-0.5 + 1.2 * x)
  y <- rbinom(D, size = n, prob = p)
  df_bin <- data.frame(id = factor(1:D), y = y, x = x, n = n)

  fit_lap <- ebp_area(
    y ~ x,
    data = df_bin,
    domain = "id",
    family = "binomial",
    trials = "n",
    method = "laplace",
    print_result = FALSE
  )

  expect_s3_class(fit_lap, "fastsae_ebp_area")
  expect_equal(length(fit_lap$df_ebp$ebp), D)
  expect_equal(nrow(fit_lap$estcoef), 2)
  expect_true(all(fit_lap$df_ebp$ebp >= 0 & fit_lap$df_ebp$ebp <= 1))
})

test_that("ebp_area input validations and error handling work", {
  data("mys", package = "fastsae")

  # 1. Error if spatial != "none" but W is NULL
  expect_error(
    ebp_area(y ~ x1, data = mys, vardir = "vardir", spatial = "bym2", W = NULL),
    "Spatial weight/proximity matrix"
  )

  # 2. Error if W dimensions do not match data
  W_bad <- matrix(0, 10, 10)
  expect_error(
    ebp_area(y ~ x1, data = mys, vardir = "vardir", spatial = "bym2", W = W_bad),
    "must be a square matrix with dimensions equal to total number of domains"
  )

  # 3. Error if binomial without trials
  df_test <- data.frame(y = c(1, 2), x = c(1, 2))
  expect_error(
    ebp_area(y ~ x, data = df_test, family = "binomial"),
    "trials"
  )
})

test_that("autoplot works dynamically with ebp and eblup objects", {
  skip_if_not_installed("INLA")
  data("mys", package = "fastsae")
  data("mys_proxmat", package = "fastsae")

  fit_ebp0 <- ebp_area(y ~ x1 + x2 + x3, data = mys, vardir = "vardir", family = "gaussian", print_result = FALSE)
  fit_ebp_sp <- ebp_area(y ~ x1 + x2 + x3, data = mys, W = mys_proxmat, vardir = "vardir", family = "gaussian", spatial = "bym2", print_result = FALSE)

  # Check that eblup is not in df_ebp (no redundancy)
  expect_false("eblup" %in% names(fit_ebp0$df_ebp))
  expect_true("ebp" %in% names(fit_ebp0$df_ebp))

  # Single model autoplots
  p1 <- autoplot(fit_ebp0, type = "estimates")
  expect_s3_class(p1, "ggplot")
  p2 <- autoplot(fit_ebp0, type = "comparison")
  expect_s3_class(p2, "ggplot")
  p3 <- autoplot(fit_ebp0, type = "mse")
  expect_s3_class(p3, "ggplot")

  # Multi-model autoplots
  p_multi <- autoplot(list("Non-Spatial" = fit_ebp0, "BYM2" = fit_ebp_sp), type = "comparison")
  expect_s3_class(p_multi, "ggplot")

  p_scatter <- autoplot(list("Non-Spatial" = fit_ebp0, "BYM2" = fit_ebp_sp), type = "scatter")
  expect_s3_class(p_scatter, "ggplot")
})

test_that("ebp_area works for Gamma response with known sampling variances", {
  skip_if_not_installed("INLA")
  data("mys", package = "fastsae")

  set.seed(123)
  mys_gamma <- mys
  mys_gamma$y_pos <- exp(0.5 + 0.3 * mys$x1 + rnorm(nrow(mys), 0, 0.2))
  mys_gamma$vardir_pos <- (0.15 * mys_gamma$y_pos)^2 # direct CV ~ 15%
  mys_gamma$y_pos[c(3, 7)] <- NA

  fit_gamma <- ebp_area(
    y_pos ~ x1,
    data = mys_gamma,
    vardir = "vardir_pos",
    family = "gamma",
    spatial = "none",
    print_result = FALSE
  )

  expect_s3_class(fit_gamma, "fastsae_ebp_area")
  expect_true(all(c("vardir", "precision", "cv_dir") %in% names(fit_gamma$df_ebp)))
  expect_true(all(fit_gamma$df_ebp$ebp > 0))
  # Unsampled domains predicted
  expect_false(any(is.na(fit_gamma$df_ebp$ebp[c(3, 7)])))
  # Check CV direct is roughly 0.15 for sampled
  sampled_idx <- which(!is.na(mys_gamma$y_pos))
  expect_equal(round(fit_gamma$df_ebp$cv_dir[sampled_idx[1]], 2), 0.15)
})

test_that("ebp_area automatically handles Binomial response with proportions", {
  skip_if_not_installed("INLA")
  data("mys", package = "fastsae")

  set.seed(42)
  mys_bin <- mys
  mys_bin$p_direct <- plogis(-0.2 + 0.4 * mys$x1)
  mys_bin$n_samp <- sample(30:80, nrow(mys), replace = TRUE)
  mys_bin$p_direct[c(2, 5)] <- NA

  # Should issue an informative message about proportions conversion
  expect_message(
    fit_bin <- ebp_area(
      p_direct ~ x1,
      data = mys_bin,
      family = "binomial",
      trials = "n_samp",
      print_result = FALSE
    ),
    "converting to integer counts"
  )

  expect_s3_class(fit_bin, "fastsae_ebp_area")
  expect_true("estimated_total" %in% names(fit_bin$df_ebp))
  expect_true(all(fit_bin$df_ebp$ebp >= 0 & fit_bin$df_ebp$ebp <= 1))
  expect_false(any(is.na(fit_bin$df_ebp$ebp[c(2, 5)])))
})

test_that("ebp_area works for Leroux CAR (spatial = generic1)", {
  skip_if_not_installed("INLA")
  data("mys", package = "fastsae")
  data("mys_proxmat", package = "fastsae")

  fit_leroux <- ebp_area(
    y ~ x1 + x2 + x3,
    data = mys,
    W = mys_proxmat,
    vardir = "vardir",
    family = "gaussian",
    spatial = "generic1",
    print_result = FALSE
  )

  expect_s3_class(fit_leroux, "fastsae_ebp_area")
  expect_false(is.null(fit_leroux$rho))
  expect_true(fit_leroux$rho >= 0 && fit_leroux$rho <= 1)
  expect_false(any(is.na(fit_leroux$df_ebp$ebp)))
})

test_that("ebp_area works for Spatial Lag Model (spatial = slm)", {
  skip_if_not_installed("INLA")
  data("mys", package = "fastsae")
  data("mys_proxmat", package = "fastsae")

  fit_slm <- ebp_area(
    y ~ x1 + x2 + x3,
    data = mys,
    W = mys_proxmat,
    vardir = "vardir",
    family = "gaussian",
    spatial = "slm",
    print_result = FALSE
  )

  expect_s3_class(fit_slm, "fastsae_ebp_area")
  expect_false(is.null(fit_slm$rho))
  expect_equal(nrow(fit_slm$estcoef), 4)
  expect_true(all(c("beta", "std.error", "zvalue", "pvalue") %in% names(fit_slm$estcoef)))
  expect_false(any(is.na(fit_slm$df_ebp$ebp)))
})
