library(testthat)
library(fastsae)

test_that("hb_area works for Gaussian Fay-Herriot model with INLA", {
  skip_if_not_installed("INLA")

  data("mys", package = "fastsae")

  # Fit standard Gaussian FH
  fit_norm <- hb_area(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = "vardir",
    family = "gaussian",
    spatial = "none",
    print_result = FALSE
  )

  expect_s3_class(fit_norm, "fastsae_hb_area")
  expect_s3_class(fit_norm, "fastsae")

  # Verify df_hb columns
  expect_true(all(c("domain", "y", "hb", "sd", "mse", "rse", "ci_lower", "ci_upper") %in% names(fit_norm$df_hb)))
  expect_equal(nrow(fit_norm$df_hb), nrow(mys))

  # Verify unsampled domains get predictions
  unsampled_idx <- which(is.na(mys$y))
  expect_true(length(unsampled_idx) > 0)
  expect_false(any(is.na(fit_norm$df_hb$hb[unsampled_idx])))
  expect_false(any(is.na(fit_norm$df_hb$sd[unsampled_idx])))

  # Verify coefficients
  expect_true(nrow(fit_norm$estcoef) == 4)
  expect_true(all(c("beta", "std.error", "zvalue", "pvalue") %in% names(fit_norm$estcoef)))

  # Verify methods
  expect_equal(length(fitted(fit_norm)), nrow(mys))
  expect_equal(length(coef(fit_norm)), 4)
  expect_output(print(fit_norm))
  expect_output(print(summary(fit_norm)))
})

test_that("hb_area works with BYM2 and Besag spatial priors", {
  skip_if_not_installed("INLA")

  data("mys", package = "fastsae")
  data("mys_proxmat", package = "fastsae")

  # 1. BYM2 via hb_area
  fit_bym2 <- hb_area(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = "vardir",
    family = "gaussian",
    spatial = "bym2",
    W = mys_proxmat,
    print_result = FALSE
  )

  expect_s3_class(fit_bym2, "fastsae_hb_area")
  expect_false(is.null(fit_bym2$phi))
  expect_true(fit_bym2$phi >= 0 && fit_bym2$phi <= 1)

  # Check unsampled domains in BYM2
  unsampled_idx <- which(is.na(mys$y))
  expect_false(any(is.na(fit_bym2$df_hb$hb[unsampled_idx])))

  # 2. BYM2 with spatial weights W
  fit_spatial <- hb_area(
    y ~ x1 + x2 + x3,
    data = mys,
    W = mys_proxmat,
    vardir = "vardir",
    family = "gaussian",
    spatial = "bym2",
    print_result = FALSE
  )
  expect_s3_class(fit_spatial, "fastsae_hb_area")
  expect_equal(length(fit_spatial$df_hb$hb), nrow(mys))

  # 3. Besag ICAR model
  fit_besag <- hb_area(
    y ~ x1 + x2 + x3,
    data = mys,
    W = mys_proxmat,
    vardir = "vardir",
    family = "gaussian",
    spatial = "besag",
    print_result = FALSE
  )
  expect_s3_class(fit_besag, "fastsae_hb_area")
  expect_false(any(is.na(fit_besag$df_hb$hb)))
})

test_that("hb_area works for Binomial response", {
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
  fit_bin <- hb_area(
    y ~ x,
    data = df_bin,
    domain = "id",
    family = "binomial",
    trials = "n",
    print_result = FALSE
  )

  expect_s3_class(fit_bin, "fastsae_hb_area")
  # Probabilities should be strictly in (0, 1)
  expect_true(all(fit_bin$df_hb$hb >= 0 & fit_bin$df_hb$hb <= 1))
  # Unsampled domains predicted
  expect_false(any(is.na(fit_bin$df_hb$hb[c(5, 12)])))
  # Check estimated totals column
  expect_true("estimated_total" %in% names(fit_bin$df_hb))
})

test_that("hb_area works for Poisson response with exposure", {
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
  fit_pois <- hb_area(
    y ~ x,
    data = df_pois,
    domain = "id",
    family = "poisson",
    exposure = "E",
    print_result = FALSE
  )

  expect_s3_class(fit_pois, "fastsae_hb_area")
  # Rates should be positive
  expect_true(all(fit_pois$df_hb$hb > 0))
  expect_false(is.na(fit_pois$df_hb$hb[3]))
  expect_true("estimated_count" %in% names(fit_pois$df_hb))
})

test_that("hb_area works with method = 'laplace' (Frequentist GLMM)", {
  set.seed(42)
  D <- 20
  x <- rnorm(D)
  n <- rep(50, D)
  p <- plogis(-0.5 + 1.2 * x)
  y <- rbinom(D, size = n, prob = p)
  df_bin <- data.frame(id = factor(1:D), y = y, x = x, n = n)

  fit_lap <- hb_area(
    y ~ x,
    data = df_bin,
    domain = "id",
    family = "binomial",
    trials = "n",
    method = "laplace",
    print_result = FALSE
  )

  expect_s3_class(fit_lap, "fastsae_hb_area")
  expect_equal(length(fit_lap$df_hb$hb), D)
  expect_equal(nrow(fit_lap$estcoef), 2)
  expect_true(all(fit_lap$df_hb$hb >= 0 & fit_lap$df_hb$hb <= 1))
})

test_that("hb_area input validations and error handling work", {
  data("mys", package = "fastsae")

  # 1. Error if spatial != "none" but W is NULL
  expect_error(
    hb_area(y ~ x1, data = mys, vardir = "vardir", spatial = "bym2", W = NULL),
    "Spatial weight/proximity matrix"
  )

  # 2. Error if W dimensions do not match data
  W_bad <- matrix(0, 10, 10)
  expect_error(
    hb_area(y ~ x1, data = mys, vardir = "vardir", spatial = "bym2", W = W_bad),
    "must be a square matrix with dimensions equal to total number of domains"
  )

  # 3. Error if binomial without trials
  df_test <- data.frame(y = c(1, 2), x = c(1, 2))
  expect_error(
    hb_area(y ~ x, data = df_test, family = "binomial"),
    "trials"
  )
})

test_that("autoplot works dynamically with hb and eblup objects", {
  skip_if_not_installed("INLA")
  data("mys", package = "fastsae")
  data("mys_proxmat", package = "fastsae")

  fit_hb0 <- hb_area(y ~ x1 + x2 + x3, data = mys, vardir = "vardir", family = "gaussian", print_result = FALSE)
  fit_hb_sp <- hb_area(y ~ x1 + x2 + x3, data = mys, W = mys_proxmat, vardir = "vardir", family = "gaussian", spatial = "bym2", print_result = FALSE)

  # Check that eblup is not in df_hb (no redundancy)
  expect_false("eblup" %in% names(fit_hb0$df_hb))
  expect_true("hb" %in% names(fit_hb0$df_hb))

  # Single model autoplots
  p1 <- autoplot(fit_hb0, type = "estimates")
  expect_s3_class(p1, "ggplot")
  p2 <- autoplot(fit_hb0, type = "comparison")
  expect_s3_class(p2, "ggplot")
  p3 <- autoplot(fit_hb0, type = "mse")
  expect_s3_class(p3, "ggplot")

  # Multi-model autoplots
  p_multi <- autoplot(list("Non-Spatial" = fit_hb0, "BYM2" = fit_hb_sp), type = "comparison")
  expect_s3_class(p_multi, "ggplot")

  p_scatter <- autoplot(list("Non-Spatial" = fit_hb0, "BYM2" = fit_hb_sp), type = "scatter")
  expect_s3_class(p_scatter, "ggplot")
})

test_that("hb_area works for Gamma response with known sampling variances", {
  skip_if_not_installed("INLA")
  data("mys", package = "fastsae")

  set.seed(123)
  mys_gamma <- mys
  mys_gamma$y_pos <- exp(0.5 + 0.3 * mys$x1 + rnorm(nrow(mys), 0, 0.2))
  mys_gamma$vardir_pos <- (0.15 * mys_gamma$y_pos)^2 # direct CV ~ 15%
  mys_gamma$y_pos[c(3, 7)] <- NA

  fit_gamma <- hb_area(
    y_pos ~ x1,
    data = mys_gamma,
    vardir = "vardir_pos",
    family = "gamma",
    spatial = "none",
    print_result = FALSE
  )

  expect_s3_class(fit_gamma, "fastsae_hb_area")
  expect_true(all(c("vardir", "precision", "cv_dir") %in% names(fit_gamma$df_hb)))
  expect_true(all(fit_gamma$df_hb$hb > 0))
  # Unsampled domains predicted
  expect_false(any(is.na(fit_gamma$df_hb$hb[c(3, 7)])))
  # Check CV direct is roughly 0.15 for sampled
  sampled_idx <- which(!is.na(mys_gamma$y_pos))
  expect_equal(round(fit_gamma$df_hb$cv_dir[sampled_idx[1]], 2), 0.15)
})

test_that("hb_area automatically handles Binomial response with proportions", {
  skip_if_not_installed("INLA")
  data("mys", package = "fastsae")

  set.seed(42)
  mys_bin <- mys
  mys_bin$p_direct <- plogis(-0.2 + 0.4 * mys$x1)
  mys_bin$n_samp <- sample(30:80, nrow(mys), replace = TRUE)
  mys_bin$p_direct[c(2, 5)] <- NA

  # Should issue an informative message about proportions conversion
  expect_message(
    fit_bin <- hb_area(
      p_direct ~ x1,
      data = mys_bin,
      family = "binomial",
      trials = "n_samp",
      print_result = FALSE
    ),
    "converting to integer counts"
  )

  expect_s3_class(fit_bin, "fastsae_hb_area")
  expect_true("estimated_total" %in% names(fit_bin$df_hb))
  expect_true(all(fit_bin$df_hb$hb >= 0 & fit_bin$df_hb$hb <= 1))
  expect_false(any(is.na(fit_bin$df_hb$hb[c(2, 5)])))
})

test_that("hb_area works for Leroux CAR (spatial = generic1)", {
  skip_if_not_installed("INLA")
  data("mys", package = "fastsae")
  data("mys_proxmat", package = "fastsae")

  fit_leroux <- hb_area(
    y ~ x1 + x2 + x3,
    data = mys,
    W = mys_proxmat,
    vardir = "vardir",
    family = "gaussian",
    spatial = "generic1",
    print_result = FALSE
  )

  expect_s3_class(fit_leroux, "fastsae_hb_area")
  expect_false(is.null(fit_leroux$rho))
  expect_true(fit_leroux$rho >= 0 && fit_leroux$rho <= 1)
  expect_false(any(is.na(fit_leroux$df_hb$hb)))
})

test_that("hb_area works for Spatial Lag Model (spatial = slm)", {
  skip_if_not_installed("INLA")
  data("mys", package = "fastsae")
  data("mys_proxmat", package = "fastsae")

  fit_slm <- hb_area(
    y ~ x1 + x2 + x3,
    data = mys,
    W = mys_proxmat,
    vardir = "vardir",
    family = "gaussian",
    spatial = "slm",
    print_result = FALSE
  )

  expect_s3_class(fit_slm, "fastsae_hb_area")
  expect_false(is.null(fit_slm$rho))
  expect_equal(nrow(fit_slm$estcoef), 4)
  expect_true(all(c("beta", "std.error", "zvalue", "pvalue") %in% names(fit_slm$estcoef)))
  expect_false(any(is.na(fit_slm$df_hb$hb)))
})

test_that("hb_area works for method = laplace via lme4 (binomial & poisson) and fallback messages", {
  set.seed(42)
  D <- 20
  df_glmm <- data.frame(
    area = 1:D,
    x = rnorm(D),
    n = sample(20:50, D, replace = TRUE),
    exposure = sample(50:100, D, replace = TRUE)
  )
  df_glmm$y_bin <- rbinom(D, size = df_glmm$n, prob = plogis(-0.2 + 0.5 * df_glmm$x))
  df_glmm$y_pois <- rpois(D, lambda = df_glmm$exposure * exp(-1.5 + 0.3 * df_glmm$x))

  # 1. method = "laplace" for binomial
  fit_lap_bin <- hb_area(
    y_bin ~ x,
    data = df_glmm,
    domain = "area",
    family = "binomial",
    trials = "n",
    method = "laplace",
    print_result = FALSE
  )
  expect_s3_class(fit_lap_bin, "fastsae_hb_area")
  expect_equal(fit_lap_bin$model, "HB-BINOMIAL (Laplace GLMM)")

  # 2. method = "laplace" for poisson with exposure
  fit_lap_pois <- hb_area(
    y_pois ~ x,
    data = df_glmm,
    domain = "area",
    family = "poisson",
    exposure = "exposure",
    method = "laplace",
    print_result = FALSE
  )
  expect_s3_class(fit_lap_pois, "fastsae_hb_area")
  expect_equal(fit_lap_pois$model, "HB-POISSON (Laplace GLMM)")

  # 3. method = "laplace" for gaussian triggers message and INLA laplace strategy
  data("mys", package = "fastsae")
  expect_message(
    fit_lap_gauss <- hb_area(
      y ~ x1,
      data = mys,
      vardir = "vardir",
      family = "gaussian",
      method = "laplace",
      print_result = FALSE
    ),
    "utilizing INLA"
  )
  expect_s3_class(fit_lap_gauss, "fastsae_hb_area")

  # 4. Gaussian model with vardir = NULL (estimated sampling variance)
  fit_no_vardir <- hb_area(
    y ~ x1 + x2,
    data = mys,
    vardir = NULL,
    family = "gaussian",
    print_result = FALSE
  )
  expect_s3_class(fit_no_vardir, "fastsae_hb_area")
  expect_true(!is.null(fit_no_vardir$random_effect_var))

  # 5. method = "laplace" with spatial != "none" and non-gaussian triggers message
  expect_message(
    hb_area(
      y_bin ~ x,
      data = df_glmm,
      domain = "area",
      family = "binomial",
      trials = "n",
      spatial = "bym2",
      W = mys_proxmat[1:D, 1:D],
      method = "laplace",
      print_result = FALSE
    ),
    "Spatial and temporal models require INLA"
  )

  # 6. method = "laplace" without exposure (Poisson)
  fit_lap_pois2 <- hb_area(y_pois ~ x, data = df_glmm, domain = "area", family = "poisson", method = "laplace", print_result = FALSE)
  expect_s3_class(fit_lap_pois2, "fastsae_hb_area")
})

test_that("hb_area Beta likelihood and vardir validation checks work", {
  skip_if_not_installed("INLA")
  set.seed(42)
  D <- 15
  df_beta <- data.frame(
    area = 1:D,
    x = rnorm(D),
    y = runif(D, 0.1, 0.8),
    vardir = runif(D, 0.001, 0.01),
    n = sample(30:60, D, replace = TRUE)
  )

  # 1. Beta with vardir (Janicki 2020)
  fit_beta_v <- hb_area(y ~ x, data = df_beta, domain = "area", family = "beta", vardir = "vardir", print_result = FALSE)
  expect_s3_class(fit_beta_v, "fastsae_hb_area")
  expect_true("precision" %in% names(fit_beta_v$df_hb))

  # 2. Beta with trials (no vardir)
  fit_beta_t <- hb_area(y ~ x, data = df_beta, domain = "area", family = "beta", trials = "n", print_result = FALSE)
  expect_s3_class(fit_beta_t, "fastsae_hb_area")
  expect_true("precision" %in% names(fit_beta_t$df_hb))

  # 3. Beta with non-positive vardir errors
  df_beta_bad <- df_beta
  df_beta_bad$vardir[1] <- 0
  expect_error(hb_area(y ~ x, data = df_beta_bad, domain = "area", family = "beta", vardir = "vardir", print_result = FALSE), "strictly positive")

  # 4. Gamma with non-positive vardir errors
  df_gamma_bad <- df_beta
  df_gamma_bad$vardir[1] <- -1
  expect_error(hb_area(y ~ x, data = df_gamma_bad, domain = "area", family = "gamma", vardir = "vardir", print_result = FALSE), "strictly positive")

  # 5. Respect CRAN limit check
  old_env <- Sys.getenv("_R_CHECK_PACKAGE_NAME_")
  Sys.setenv("_R_CHECK_PACKAGE_NAME_" = "fastsae")
  fit_cran <- hb_area(y ~ x, data = df_beta, domain = "area", family = "beta", trials = "n", print_result = FALSE)
  expect_s3_class(fit_cran, "fastsae_hb_area")
  if (nzchar(old_env)) Sys.setenv("_R_CHECK_PACKAGE_NAME_" = old_env) else Sys.unsetenv("_R_CHECK_PACKAGE_NAME_")
})

test_that("hb_area handles ST interaction types", {
  skip_if_not_installed("INLA")
  set.seed(42)
  D <- 6
  T_per <- 2
  N <- D * T_per

  df_st <- expand.grid(time = 1:T_per, domain = 1:D)
  df_st$x <- rnorm(N)
  df_st$y <- 5 + 0.5 * df_st$x + rnorm(N, 0, 0.2)
  df_st$vardir <- rep(0.1, N)

  W_st <- matrix(0, D, D)
  for (i in 1:(D - 1)) {
    W_st[i, i + 1] <- W_st[i + 1, i] <- 1
  }

  # separable with bym2
  fit_sep_bym2 <- hb_area(y ~ x, data = df_st, domain = "domain", time = "time", vardir = "vardir",
                          spatial = "bym2", temporal = "ar1", st_interaction = "separable", W = W_st, print_result = FALSE)
  expect_s3_class(fit_sep_bym2, "fastsae_hb_area")

  # separable with bym
  fit_sep_bym <- hb_area(y ~ x, data = df_st, domain = "domain", time = "time", vardir = "vardir",
                         spatial = "bym", temporal = "ar1", st_interaction = "separable", W = W_st, print_result = FALSE)
  expect_s3_class(fit_sep_bym, "fastsae_hb_area")

  # separable with generic1 (crashes inside INLA GMRFLib)
  expect_error(
    hb_area(y ~ x, data = df_st, domain = "domain", time = "time", vardir = "vardir",
            spatial = "generic1", temporal = "ar1", st_interaction = "separable", W = W_st, print_result = FALSE),
    "Fitting model with INLA failed"
  )

  # separable with spatial = none
  fit_sep_none <- hb_area(y ~ x, data = df_st, domain = "domain", time = "time", vardir = "vardir",
                          spatial = "none", temporal = "ar1", st_interaction = "separable", print_result = FALSE)
  expect_s3_class(fit_sep_none, "fastsae_hb_area")

  # type1
  fit_t1 <- hb_area(y ~ x, data = df_st, domain = "domain", time = "time", vardir = "vardir",
                    spatial = "none", temporal = "ar1", st_interaction = "type1", print_result = FALSE)
  expect_s3_class(fit_t1, "fastsae_hb_area")

  # type2 (INLA failure)
  expect_error(
    hb_area(y ~ x, data = df_st, domain = "domain", time = "time", vardir = "vardir",
            spatial = "none", temporal = "ar1", st_interaction = "type2", print_result = FALSE),
    "Fitting model with INLA failed"
  )

  # type3 (INLA failure)
  expect_error(
    hb_area(y ~ x, data = df_st, domain = "domain", time = "time", vardir = "vardir",
            spatial = "besag", temporal = "ar1", st_interaction = "type3", W = W_st, print_result = FALSE),
    "Fitting model with INLA failed"
  )

  # type4 (INLA failure)
  expect_error(
    hb_area(y ~ x, data = df_st, domain = "domain", time = "time", vardir = "vardir",
            spatial = "besag", temporal = "ar1", st_interaction = "type4", W = W_st, print_result = FALSE),
    "Fitting model with INLA failed"
  )
})

test_that(".fit_glmm_laplace and hb_area prior/spatial branches work", {
  skip_if_not_installed("lme4")
  # 1. .fit_glmm_laplace gaussian and binomial without trials
  df_test <- data.frame(y = rnorm(10), x = rnorm(10), d = rep(1:5, 2))
  fit_g <- fastsae:::.fit_glmm_laplace(y ~ x, data = df_test, domain = df_test$d, family = "gaussian")
  expect_s3_class(fit_g, "fastsae")

  df_bin <- data.frame(y = rbinom(10, 1, 0.5), x = rnorm(10), d = rep(1:5, 2))
  fit_b <- fastsae:::.fit_glmm_laplace(y ~ x, data = df_bin, domain = df_bin$d, family = "binomial")
  expect_s3_class(fit_b, "fastsae")

  # 2. vardir <= 0 error for gaussian
  df_bad_var <- data.frame(y = 1:5, x = 1:5, vardir = c(1, 0, 1, 1, 1))
  expect_error(hb_area(y ~ x, data = df_bad_var, vardir = "vardir", family = "gaussian"), "strictly positive")

  skip_if_not_installed("INLA")
  data("mys", package = "fastsae")
  data("mys_proxmat", package = "fastsae")

  # 3. spatial = "bym"
  fit_bym <- hb_area(y ~ x1, data = mys, vardir = "vardir", spatial = "bym", W = mys_proxmat, print_result = FALSE)
  expect_s3_class(fit_bym, "fastsae_hb_area")

  # 4. temporal = "ar1" with prior_rho_time and domain-specific st_interaction
  panel_data <- mys[rep(1:10, each = 3), ]
  panel_data$time <- rep(1:3, 10)
  panel_data$area <- rep(1:10, each = 3)
  fit_ds <- hb_area(
    y ~ x1,
    data = panel_data,
    domain = "area",
    time = "time",
    vardir = "vardir",
    temporal = "ar1",
    prior_rho_time = list(prior = "normal", param = c(0, 0.15)),
    st_interaction = "domain-specific",
    print_result = FALSE
  )
  expect_s3_class(fit_ds, "fastsae_hb_area")

  # 5. generic1 and slm with prior_rho
  fit_gen <- hb_area(
    y ~ x1,
    data = mys,
    vardir = "vardir",
    spatial = "generic1",
    W = mys_proxmat,
    prior_rho = list(prior = "normal", param = c(0, 1)),
    print_result = FALSE
  )
  expect_s3_class(fit_gen, "fastsae_hb_area")

  fit_slm <- hb_area(
    y ~ x1,
    data = mys,
    vardir = "vardir",
    spatial = "slm",
    W = mys_proxmat,
    prior_rho = list(prior = "normal", param = c(0, 1)),
    print_result = FALSE
  )
  expect_s3_class(fit_slm, "fastsae_hb_area")
})


