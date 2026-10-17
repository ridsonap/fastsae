library(testthat)
library(fastsae)

test_that("hb_area matches tipsae::fit_sae on spatio-temporal Beta SAE", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("tipsae")
  skip_if_not_installed("spdep")

  data("emilia", package = "tipsae")
  data("emilia_shp", package = "tipsae")

  # 1. Convert shapefile to spatial adjacency matrix W
  nb <- spdep::poly2nb(emilia_shp)
  W <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
  rownames(W) <- colnames(W) <- emilia_shp$NAME_DISTRICT

  # 2. Converged tipsae (Stan) posterior, stored because Stan draws differ
  # across platforms; rebuild with fixtures/make-tipsae-emilia.R
  ref <- readRDS(test_path("fixtures", "tipsae-emilia.rds"))
  expect_lt(attr(ref, "info")$max_rhat_theta, 1.01)

  # 3. Fit fastsae INLA model with equivalent domain-specific RW1 + Besag structure
  fit_inla <- hb_area(
    hcr ~ x,
    data = emilia,
    domain = "id",
    time = "year",
    vardir = "vars",
    family = "beta",
    spatial = "besag",
    temporal = "rw1",
    st_interaction = "domain-specific",
    W = W,
    print_result = FALSE
  )

  expect_s3_class(fit_inla, "fastsae_hb_area")
  expect_s3_class(fit_inla, "fastsae")

  # 4. Check output structure
  df_hb <- fit_inla$df_hb
  expect_true(all(c("domain", "time", "y", "hb", "sd", "mse", "rse", "ci_lower", "ci_upper") %in% names(df_hb)))
  expect_true(all(c("random_effect", "random_effect_spatial", "random_effect_temporal") %in% names(df_hb)))
  expect_equal(nrow(df_hb), nrow(emilia))

  # 5. Merge and compare estimates
  df_comp <- merge(
    data.frame(id = df_hb$domain, year = df_hb$time, inla_est = df_hb$hb, inla_sd = df_hb$sd),
    ref,
    by = c("id", "year")
  )
  expect_equal(nrow(df_comp), nrow(emilia))

  cor_val <- stats::cor(df_comp$inla_est, df_comp$hb)
  mae_val <- mean(abs(df_comp$inla_est - df_comp$hb))
  sd_ratio <- df_comp$inla_sd / df_comp$sd

  # Against the converged tipsae fit: r = 0.999, MAE = 0.0017 and SD ratios
  # 0.81 to 1.22 (the priors differ)
  expect_gt(cor_val, 0.99)
  expect_lt(mae_val, 0.005)
  expect_true(all(sd_ratio > 0.7 & sd_ratio < 1.4))

  # Check hyperparameters
  expect_false(is.null(fit_inla$random_effect_var))
  expect_false(is.null(fit_inla$random_effect_var_time))
  expect_gt(fit_inla$random_effect_var, 0)
  expect_gt(fit_inla$random_effect_var_time, 0)
})

test_that("hb_area works for Gaussian Spatio-Temporal Fay-Herriot on panel data", {
  skip_if_not_installed("INLA")

  data("mys_panel", package = "fastsae")
  data("mys_proxmat", package = "fastsae")

  # Subset to 3 years for fast test execution
  panel_sub <- mys_panel[mys_panel$year %in% c(2022, 2023, 2024), ]

  # 1. Spatio-temporal with BYM2 + AR(1)
  fit_st <- hb_area(
    y ~ x1 + x2 + x3,
    data = panel_sub,
    domain = "area",
    time = "year",
    vardir = "vardir",
    family = "gaussian",
    spatial = "bym2",
    temporal = "ar1",
    st_interaction = "none",
    W = mys_knn,
    print_result = FALSE
  )

  expect_s3_class(fit_st, "fastsae_hb_area")
  expect_equal(nrow(fit_st$df_hb), nrow(panel_sub))
  expect_true("time" %in% names(fit_st$df_hb))
  expect_false(any(is.na(fit_st$df_hb$hb)))
  expect_false(any(is.na(fit_st$df_hb$sd)))

  # Check temporal hyperparameters
  expect_false(is.null(fit_st$random_effect_var_time))
  expect_false(is.null(fit_st$rho_time))
  expect_true(fit_st$rho_time >= -1 && fit_st$rho_time <= 1)

  # 2. Test separable space-time interaction (ST-FH)
  fit_sep <- hb_area(
    y ~ x1 + x2,
    data = panel_sub,
    domain = "area",
    time = "year",
    vardir = "vardir",
    family = "gaussian",
    spatial = "besag",
    temporal = "ar1",
    st_interaction = "separable",
    W = mys_knn,
    print_result = FALSE
  )
  expect_s3_class(fit_sep, "fastsae_hb_area")
  expect_false(any(is.na(fit_sep$df_hb$hb)))
})

test_that("hb_area works for Poisson and Binomial Spatio-Temporal models", {
  skip_if_not_installed("INLA")

  set.seed(123)
  D <- 15
  T_periods <- 3
  N <- D * T_periods

  df_panel <- expand.grid(time = 2021:(2021 + T_periods - 1), domain = 1:D)
  df_panel$x <- rnorm(N)
  df_panel$exposure <- runif(N, 50, 200)
  df_panel$trials <- sample(40:120, N, replace = TRUE)

  # Simulate Poisson counts
  lambda_true <- df_panel$exposure * exp(-1.5 + 0.4 * df_panel$x)
  df_panel$y_pois <- rpois(N, lambda = lambda_true)

  # Simulate Binomial counts
  prob_true <- plogis(-0.8 + 0.5 * df_panel$x)
  df_panel$y_bin <- rbinom(N, size = df_panel$trials, prob = prob_true)

  # Synthetic W matrix
  W_syn <- matrix(0, D, D)
  for (i in 1:(D - 1)) {
    W_syn[i, i + 1] <- 1
    W_syn[i + 1, i] <- 1
  }

  # 1. Poisson Spatio-Temporal
  fit_pois <- hb_area(
    y_pois ~ x,
    data = df_panel,
    domain = "domain",
    time = "time",
    exposure = "exposure",
    family = "poisson",
    spatial = "besag",
    temporal = "rw1",
    W = W_syn,
    print_result = FALSE
  )
  expect_s3_class(fit_pois, "fastsae_hb_area")
  expect_true(all(fit_pois$df_hb$hb > 0))
  expect_true("rate" %in% names(fit_pois$df_hb))
  expect_true("estimated_count" %in% names(fit_pois$df_hb))

  # 2. Binomial Spatio-Temporal
  fit_bin <- hb_area(
    y_bin ~ x,
    data = df_panel,
    domain = "domain",
    time = "time",
    trials = "trials",
    family = "binomial",
    spatial = "bym2",
    temporal = "ar1",
    W = W_syn,
    print_result = FALSE
  )
  expect_s3_class(fit_bin, "fastsae_hb_area")
  expect_true(all(fit_bin$df_hb$hb >= 0 & fit_bin$df_hb$hb <= 1))
  expect_true("estimated_total" %in% names(fit_bin$df_hb))
})

test_that("hb_area works for Negative Binomial and Gamma Spatio-Temporal models", {
  skip_if_not_installed("INLA")

  set.seed(456)
  D <- 12
  T_periods <- 3
  N <- D * T_periods

  df_panel <- expand.grid(time = 1:T_periods, domain = 1:D)
  df_panel$x <- rnorm(N)
  df_panel$exposure <- runif(N, 100, 300)
  df_panel$y_nbin <- rnbinom(N, size = 5, mu = df_panel$exposure * exp(-2 + 0.3 * df_panel$x))
  df_panel$y_gamma <- rgamma(N, shape = 4, scale = exp(1 + 0.2 * df_panel$x) / 4)
  df_panel$vardir <- (df_panel$y_gamma^2) / 4

  W_syn <- matrix(0, D, D)
  for (i in 1:(D - 1)) {
    W_syn[i, i + 1] <- 1
    W_syn[i + 1, i] <- 1
  }

  # 1. Negative Binomial ST
  fit_nbin <- hb_area(
    y_nbin ~ x,
    data = df_panel,
    domain = "domain",
    time = "time",
    exposure = "exposure",
    family = "nbinomial",
    spatial = "bym2",
    temporal = "rw1",
    W = W_syn,
    print_result = FALSE
  )
  expect_s3_class(fit_nbin, "fastsae_hb_area")
  expect_true(all(fit_nbin$df_hb$hb > 0))

  # 2. Gamma ST with known vardir
  fit_gamma <- hb_area(
    y_gamma ~ x,
    data = df_panel,
    domain = "domain",
    time = "time",
    vardir = "vardir",
    family = "gamma",
    spatial = "besag",
    temporal = "rw1",
    W = W_syn,
    print_result = FALSE
  )
  expect_s3_class(fit_gamma, "fastsae_hb_area")
  expect_true(all(fit_gamma$df_hb$hb > 0))
  expect_true("precision" %in% names(fit_gamma$df_hb))
})

test_that("hb_area temporal input validation and error handling work", {
  skip_if_not_installed("INLA")
  data("mys", package = "fastsae")

  # 1. Error when temporal != "none" but time is NULL
  expect_error(
    hb_area(y ~ x1, data = mys, vardir = "vardir", temporal = "rw1"),
    "time"
  )

  # 2. Error when time is specified without domain
  expect_error(
    hb_area(y ~ x1, data = mys, vardir = "vardir", time = "year", temporal = "rw1"),
    "domain.*must also be specified"
  )

  # 3. Error when time has < 2 unique periods
  df_single_time <- mys
  df_single_time$area_id <- seq_len(nrow(mys))
  df_single_time$year <- 2024
  expect_error(
    hb_area(y ~ x1, data = df_single_time, domain = "area_id", vardir = "vardir", time = "year", temporal = "rw1"),
    "at least 2 distinct time periods"
  )

  # 4. Print and summary work with temporal model
  df_temp <- data.frame(
    y = rnorm(20, mean = 5),
    vardir = runif(20, 0.1, 0.5),
    area = rep(1:10, 2),
    year = rep(c(2023, 2024), each = 10)
  )
  fit_pure_temp <- hb_area(
    y ~ 1,
    data = df_temp,
    domain = "area",
    time = "year",
    vardir = "vardir",
    family = "gaussian",
    temporal = "rw1",
    print_result = FALSE
  )
  expect_output(print(fit_pure_temp))
  expect_output(print(summary(fit_pure_temp)))
})
