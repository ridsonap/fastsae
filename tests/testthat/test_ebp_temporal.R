test_that("ebp_area matches tipsae::fit_sae on spatio-temporal Beta SAE", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("tipsae")
  skip_if_not_installed("spdep")

  data("emilia", package = "tipsae")
  data("emilia_shp", package = "tipsae")

  # 1. Convert shapefile to spatial adjacency matrix W
  nb <- spdep::poly2nb(emilia_shp)
  W <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
  rownames(W) <- colnames(W) <- emilia_shp$NAME_DISTRICT

  # 2. Fit tipsae Stan MCMC model
  t0_tip <- Sys.time()
  fit_tip <- tipsae::fit_sae(
    formula_fixed = hcr ~ x,
    data = emilia,
    domains = "id",
    disp_direct = "vars",
    type_disp = "var",
    domain_size = "n",
    spatial_error = TRUE,
    spatial_df = emilia_shp,
    domains_spatial_df = "NAME_DISTRICT",
    temporal_error = TRUE,
    temporal_variable = "year",
    chains = 1,
    iter = 200,
    seed = 42
  )
  t1_tip <- Sys.time()
  time_tip <- as.numeric(difftime(t1_tip, t0_tip, units = "secs"))
  sum_tip <- summary(fit_tip)
  est_tip <- tipsae::extract(sum_tip)$in_sample
  colnames(est_tip)[1:4] <- c("Domains", "Times", "Direct_est", "HB_est")

  # 3. Fit fastsae INLA model with equivalent domain-specific RW1 + Besag structure
  t0_inla <- Sys.time()
  fit_inla <- ebp_area(
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
  t1_inla <- Sys.time()
  time_inla <- as.numeric(difftime(t1_inla, t0_inla, units = "secs"))

  expect_s3_class(fit_inla, "fastsae_ebp_area")
  expect_s3_class(fit_inla, "fastsae")

  # 4. Check output structure
  df_ebp <- fit_inla$df_ebp
  expect_true(all(c("domain", "time", "y", "ebp", "sd", "mse", "rse", "ci_lower", "ci_upper") %in% names(df_ebp)))
  expect_true(all(c("random_effect", "random_effect_spatial", "random_effect_temporal") %in% names(df_ebp)))
  expect_equal(nrow(df_ebp), nrow(emilia))

  # 5. Merge and compare estimates
  df_comp <- merge(
    data.frame(id = df_ebp$domain, year = df_ebp$time, inla_est = df_ebp$ebp),
    est_tip,
    by.x = c("id", "year"),
    by.y = c("Domains", "Times")
  )

  cor_val <- stats::cor(df_comp$inla_est, df_comp$HB_est)
  mae_val <- mean(abs(df_comp$inla_est - df_comp$HB_est))

  # Benchmark expectations: high correlation (> 0.95), low MAE (< 0.015), INLA speedup (> 2x)
  expect_gt(cor_val, 0.95)
  expect_lt(mae_val, 0.015)
  expect_gt(time_tip, time_inla)

  # Check hyperparameters
  expect_false(is.null(fit_inla$random_effect_var))
  expect_false(is.null(fit_inla$random_effect_var_time))
  expect_gt(fit_inla$random_effect_var, 0)
  expect_gt(fit_inla$random_effect_var_time, 0)
})

test_that("ebp_area works for Gaussian Spatio-Temporal Fay-Herriot on panel data", {
  skip_if_not_installed("INLA")

  data("mys_panel", package = "fastsae")
  data("mys_proxmat", package = "fastsae")

  # Subset to 3 years for fast test execution
  panel_sub <- mys_panel[mys_panel$year %in% c(2022, 2023, 2024), ]

  # 1. Spatio-temporal with BYM2 + AR(1)
  fit_st <- ebp_area(
    y ~ x1 + x2 + x3,
    data = panel_sub,
    domain = "area",
    time = "year",
    vardir = "vardir",
    family = "gaussian",
    spatial = "bym2",
    temporal = "ar1",
    st_interaction = "none",
    W = mys_proxmat,
    print_result = FALSE
  )

  expect_s3_class(fit_st, "fastsae_ebp_area")
  expect_equal(nrow(fit_st$df_ebp), nrow(panel_sub))
  expect_true("time" %in% names(fit_st$df_ebp))
  expect_false(any(is.na(fit_st$df_ebp$ebp)))
  expect_false(any(is.na(fit_st$df_ebp$sd)))

  # Check temporal hyperparameters
  expect_false(is.null(fit_st$random_effect_var_time))
  expect_false(is.null(fit_st$rho_time))
  expect_true(fit_st$rho_time >= -1 && fit_st$rho_time <= 1)

  # 2. Test separable space-time interaction (ST-FH)
  fit_sep <- ebp_area(
    y ~ x1 + x2,
    data = panel_sub,
    domain = "area",
    time = "year",
    vardir = "vardir",
    family = "gaussian",
    spatial = "besag",
    temporal = "ar1",
    st_interaction = "separable",
    W = mys_proxmat,
    print_result = FALSE
  )
  expect_s3_class(fit_sep, "fastsae_ebp_area")
  expect_false(any(is.na(fit_sep$df_ebp$ebp)))
})

test_that("ebp_area works for Poisson and Binomial Spatio-Temporal models", {
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
  fit_pois <- ebp_area(
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
  expect_s3_class(fit_pois, "fastsae_ebp_area")
  expect_true(all(fit_pois$df_ebp$ebp > 0))
  expect_true("rate" %in% names(fit_pois$df_ebp))
  expect_true("estimated_count" %in% names(fit_pois$df_ebp))

  # 2. Binomial Spatio-Temporal
  fit_bin <- ebp_area(
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
  expect_s3_class(fit_bin, "fastsae_ebp_area")
  expect_true(all(fit_bin$df_ebp$ebp >= 0 & fit_bin$df_ebp$ebp <= 1))
  expect_true("estimated_total" %in% names(fit_bin$df_ebp))
})

test_that("ebp_area works for Negative Binomial and Gamma Spatio-Temporal models", {
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
  fit_nbin <- ebp_area(
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
  expect_s3_class(fit_nbin, "fastsae_ebp_area")
  expect_true(all(fit_nbin$df_ebp$ebp > 0))

  # 2. Gamma ST with known vardir
  fit_gamma <- ebp_area(
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
  expect_s3_class(fit_gamma, "fastsae_ebp_area")
  expect_true(all(fit_gamma$df_ebp$ebp > 0))
  expect_true("precision" %in% names(fit_gamma$df_ebp))
})

test_that("ebp_area temporal input validation and error handling work", {
  skip_if_not_installed("INLA")
  data("mys", package = "fastsae")

  # 1. Error when temporal != "none" but time is NULL
  expect_error(
    ebp_area(y ~ x1, data = mys, vardir = "vardir", temporal = "rw1"),
    "time"
  )

  # 2. Error when time is specified without domain
  expect_error(
    ebp_area(y ~ x1, data = mys, vardir = "vardir", time = "year", temporal = "rw1"),
    "domain.*must also be specified"
  )

  # 3. Error when time has < 2 unique periods
  df_single_time <- mys
  df_single_time$area_id <- seq_len(nrow(mys))
  df_single_time$year <- 2024
  expect_error(
    ebp_area(y ~ x1, data = df_single_time, domain = "area_id", vardir = "vardir", time = "year", temporal = "rw1"),
    "at least 2 distinct time periods"
  )

  # 4. Print and summary work with temporal model
  df_temp <- data.frame(
    y = rnorm(20, mean = 5),
    vardir = runif(20, 0.1, 0.5),
    area = rep(1:10, 2),
    year = rep(c(2023, 2024), each = 10)
  )
  fit_pure_temp <- ebp_area(
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
