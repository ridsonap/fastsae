library(testthat)
library(fastsae)

test_that("hb_area Gaussian estimates match sae::mseFH benchmark (r > 0.99)", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("sae")

  data("mys", package = "fastsae")
  mysnona <- mys[!is.na(mys[["y"]]), ]

  fit_sae <- sae::mseFH(
    mysnona[["y"]] ~ mysnona[["x1"]] + mysnona[["x2"]] + mysnona[["x3"]],
    vardir = mysnona[["vardir"]],
    method = "REML"
  )

  fit_hb <- hb_area(
    y ~ x1 + x2 + x3,
    data = mysnona,
    vardir = "vardir",
    domain = "area",
    family = "gaussian",
    spatial = "none",
    print_result = FALSE
  )

  sae_eblup <- as.numeric(fit_sae[["est"]][["eblup"]])
  hb_pred <- as.numeric(fit_hb[["df_hb"]][["hb"]])

  # Correlation between Frequentist EBLUP (REML) and Bayesian Posterior Mean (INLA)
  r_est <- cor(sae_eblup, hb_pred)
  expect_gt(r_est, 0.99)

  # Average difference should be small (< 0.25 on this scale)
  expect_lt(mean(abs(sae_eblup - hb_pred)), 0.25)

  # MSE / Posterior Variance correlation
  sae_mse <- as.numeric(fit_sae[["mse"]])
  hb_mse <- as.numeric(fit_hb[["df_hb"]][["mse"]])
  r_mse <- cor(sae_mse, hb_mse)
  expect_gt(r_mse, 0.98)
})

test_that("hb_area achieves parameter recovery and nominal coverage rate on synthetic data", {
  skip_if_not_installed("INLA")
  sim <- sim_area_data(
    D = 40,
    beta = c(10.0, 1.0, 0.5),
    sigma_u = 0.5,
    n_unsampled = 0,
    seed = 2026
  )

  fit_sim <- hb_area(
    y_gaussian ~ x1 + x2,
    data = sim$data,
    vardir = "vardir",
    family = "gaussian",
    spatial = "none",
    print_result = FALSE
  )

  diag_sim <- diagnose(fit_sim, truth = sim$data$truth_gaussian)
  expect_s3_class(diag_sim, "fastsae_diagnose")

  # Parameter recovery: True fixed effects within 95% Credible Interval
  b0_ci <- c(fit_sim$estcoef["(Intercept)", "ci_lower"], fit_sim$estcoef["(Intercept)", "ci_upper"])
  expect_true(10.0 >= b0_ci[1] && 10.0 <= b0_ci[2])

  b1_ci <- c(fit_sim$estcoef["x1", "ci_lower"], fit_sim$estcoef["x1", "ci_upper"])
  expect_true(1.0 >= b1_ci[1] && 1.0 <= b1_ci[2])

  b2_ci <- c(fit_sim$estcoef["x2", "ci_lower"], fit_sim$estcoef["x2", "ci_upper"])
  expect_true(0.5 >= b2_ci[1] && 0.5 <= b2_ci[2])

  # Low relative bias across areas (< 10%)
  expect_lt(diag_sim$simulation_metrics$mean_abs_rb, 10)

  # 95% Credible Interval Coverage Rate should be at least 85%
  expect_gte(diag_sim$simulation_metrics$coverage_rate, 85)
})

test_that("diagnose extracts Bayesian metrics and autoplot supports pit and cpo", {
  skip_if_not_installed("INLA")
  data("mys", package = "fastsae")

  fit_hb <- hb_area(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = ~vardir,
    domain = ~area,
    family = "gaussian",
    print_result = FALSE
  )

  d_hb <- diagnose(fit_hb)
  expect_s3_class(d_hb, "fastsae_diagnose")

  # Bayesian metrics slot
  expect_true(!is.null(d_hb$bayesian_metrics))
  expect_true(is.numeric(d_hb$bayesian_metrics$waic))
  expect_true(is.numeric(d_hb$bayesian_metrics$dic))
  expect_true(is.numeric(d_hb$bayesian_metrics$log_mlik))
  expect_true(!is.null(d_hb$bayesian_metrics$cpo_summary))
  expect_true(!is.null(d_hb$bayesian_metrics$pit_test))

  # Dataframe includes cpo and pit columns
  expect_true("cpo" %in% names(d_hb$df_diag))
  expect_true("pit" %in% names(d_hb$df_diag))

  # Autoplot methods for Bayesian calibration
  p_pit <- autoplot(d_hb, type = "pit")
  expect_s3_class(p_pit, "ggplot")

  p_cpo <- autoplot(d_hb, type = "cpo")
  expect_s3_class(p_cpo, "ggplot")

  # Print report contains Bayesian section
  expect_message(print(d_hb), "Bayesian Information Criteria")
})

test_that("hb_area input validation and error checking", {
  df_dummy <- data.frame(y = 1:5, x = 1:5, area = 1:5, year = c(1, 1, 2, 2, 3), vardir = rep(0.5, 5))

  # 1. Non-dataframe data
  expect_error(hb_area(y ~ x, data = "not_df"), "must be a data frame")

  # 2. time specified without domain
  expect_error(hb_area(y ~ x, data = df_dummy, time = "year"), "domain.*must also be specified")

  # 3. temporal != none without time
  expect_error(hb_area(y ~ x, data = df_dummy, domain = "area", temporal = "ar1"), "time.*must be specified")

  # 4. temporal != none with < 2 unique times
  df_one_time <- df_dummy
  df_one_time$year <- 1
  expect_error(hb_area(y ~ x, data = df_one_time, domain = "area", time = "year", temporal = "ar1"), "at least 2 distinct time periods")

  # 5. Binomial without trials
  expect_error(hb_area(y ~ x, data = df_dummy, domain = "area", family = "binomial"), "trials.*must be specified")

  # 6. Gamma with non-positive y
  df_gamma_bad <- df_dummy
  df_gamma_bad$y[1] <- 0
  expect_error(hb_area(y ~ x, data = df_gamma_bad, domain = "area", family = "gamma"), "strictly positive")

  # 7. Spatial != none without W
  expect_error(hb_area(y ~ x, data = df_dummy, domain = "area", spatial = "bym2"), "must be provided")
})

test_that("inla_utils .convert_spatial_weights checks", {
  # 1. Unsupported spatial object type
  expect_error(fastsae::: .convert_spatial_weights("bad_W", n_domains = 5), "Unsupported spatial object type")

  # 2. Dimension mismatch
  bad_mat <- matrix(0, nrow = 3, ncol = 3)
  expect_error(fastsae::: .convert_spatial_weights(bad_mat, n_domains = 5), "must be a square matrix with dimensions equal to total number of domains")

  # 3. Isolated nodes warning
  iso_mat <- matrix(0, nrow = 4, ncol = 4)
  iso_mat[1, 2] <- iso_mat[2, 1] <- 1
  expect_warning(fastsae::: .convert_spatial_weights(iso_mat, n_domains = 4), "isolated nodes")

  # 4. inla.graph object
  g_obj <- structure(list(), class = "inla.graph")
  res_g <- fastsae::: .convert_spatial_weights(g_obj, n_domains = 5)
  expect_identical(res_g$graph, g_obj)

  # 5. listw and nb objects from spdep
  if (requireNamespace("spdep", quietly = TRUE)) {
    W_mat <- matrix(0, 4, 4)
    W_mat[1, 2] <- W_mat[2, 1] <- W_mat[2, 3] <- W_mat[3, 2] <- W_mat[3, 4] <- W_mat[4, 3] <- 1
    lw <- spdep::mat2listw(W_mat, style = "B")
    res_lw <- fastsae::: .convert_spatial_weights(lw, n_domains = 4)
    expect_equal(unname(as.matrix(res_lw$W_mat)), W_mat, ignore_attr = TRUE)

    res_nb <- fastsae::: .convert_spatial_weights(lw$neighbours, n_domains = 4)
    expect_equal(unname(as.matrix(res_nb$W_mat)), W_mat, ignore_attr = TRUE)
  }
})

test_that("hb_twofold input validation and error checking", {
  df_sub <- data.frame(y = 1:5, x = 1:5, domain = c(1, 1, 2, 2, 3), subarea = 1:5, vardir = rep(0.5, 5))

  # 1. Non-dataframe data
  expect_error(hb_twofold(y ~ x, data = "not_df", domain = "domain", vardir = "vardir"), "must be a data frame")

  # 2. Empty data
  expect_error(hb_twofold(y ~ x, data = df_sub[0, ], domain = "domain", vardir = "vardir"), "contains 0 rows")

  # 3. Missing domain
  expect_error(hb_twofold(y ~ x, data = df_sub, vardir = "vardir"), "domain.*must be specified")

  # 4. Gaussian without vardir
  expect_error(hb_twofold(y ~ x, data = df_sub, domain = "domain", family = "gaussian"), "vardir.*must be provided")

  # 5. Binomial without trials
  expect_error(hb_twofold(y ~ x, data = df_sub, domain = "domain", family = "binomial"), "trials.*must be provided")
})

test_that("hb_unit input validation and error checking", {
  df_unit <- data.frame(y = 1:10, x = 1:10, d = rep(1:2, each = 5))
  df_pop <- data.frame(d = 1:2, x = c(5, 6), N = c(100, 100))

  # 1. Non-formula
  expect_error(hb_unit("y ~ x", unit_data = df_unit, Xpop = df_pop, domain_var = "d"), "valid formula object")

  # 2. Non-dataframe unit_data
  expect_error(hb_unit(y ~ x, unit_data = "bad", Xpop = df_pop, domain_var = "d"), "must be a data frame")

  # 3. Non-dataframe Xpop
  expect_error(hb_unit(y ~ x, unit_data = df_unit, Xpop = "bad", domain_var = "d"), "must be a data frame")

  # 4. domain_var not a single character string
  expect_error(hb_unit(y ~ x, unit_data = df_unit, Xpop = df_pop, domain_var = 123), "single character string")

  # 5. domain_var missing in unit_data
  expect_error(hb_unit(y ~ x, unit_data = df_unit, Xpop = df_pop, domain_var = "nonexistent"), "not found in.*unit_data")

  # 6. domain_var missing in Xpop
  df_pop_bad <- df_pop
  names(df_pop_bad)[1] <- "other_d"
  expect_error(hb_unit(y ~ x, unit_data = df_unit, Xpop = df_pop_bad, domain_var = "d"), "not found in.*Xpop")

  # 7. popsize_var missing in Xpop
  expect_error(hb_unit(y ~ x, unit_data = df_unit, Xpop = df_pop, domain_var = "d", popsize_var = "wrong_pop"), "not found in.*Xpop")

  # 8. Domains in unit_data missing from Xpop
  df_pop_missing <- df_pop[1, , drop = FALSE]
  expect_error(hb_unit(y ~ x, unit_data = df_unit, Xpop = df_pop_missing, domain_var = "d"), "missing from `Xpop`")
})

test_that(".extract_inla_results handles all model edge cases and parameter combinations", {
  mock_fit <- list(
    summary.fixed = data.frame(
      mean = c(1.0, 2.0),
      sd = c(0.1, 0.2),
      "0.025quant" = c(0.8, 1.6),
      "0.975quant" = c(1.2, 2.4),
      row.names = c("(Intercept)", "x1"),
      check.names = FALSE
    ),
    summary.hyperpar = data.frame(
      mean = c(2.0, 0.5, 0.4),
      row.names = c("Precision for ..domain_id..", "Phi for ..domain_id..", "Beta for ..domain_id..")
    ),
    summary.random = list(
      `..domain_id..` = data.frame(mean = c(0.1, -0.1, 0.2))
    ),
    summary.fitted.values = data.frame(
      mean = c(1.1, 2.1, 3.1),
      sd = c(0.2, 0.2, 0.3),
      "0.025quant" = c(0.7, 1.7, 2.5),
      "0.975quant" = c(1.5, 2.5, 3.7),
      check.names = FALSE
    ),
    summary.linear.predictor = data.frame(
      mean = c(1.1, 2.1, 3.1)
    ),
    dic = list(dic = 100, p.eff = 2),
    waic = list(waic = 105, p.eff = 2.1),
    mlik = matrix(c(-50, -50), nrow = 1)
  )

  # 1. domain_id = NULL and length(spat_vals) == n_obs
  res1 <- fastsae:::.extract_inla_results(
    fit = mock_fit,
    data = data.frame(y = 1:3),
    y = 1:3,
    domain = 1:3,
    family = "gaussian",
    spatial = "bym2",
    domain_id = NULL
  )
  expect_s3_class(res1, "fastsae_hb_area")
  expect_equal(res1$phi, 0.5)

  # 2. st_interaction == "domain-specific" with domain_id and time_id
  mock_fit_st <- mock_fit
  mock_fit_st$summary.random[["..time_id.."]] <- data.frame(mean = rep(0.05, 6))
  mock_fit_st$summary.fitted.values <- data.frame(
    mean = rep(1.0, 6), sd = rep(0.1, 6), "0.025quant" = rep(0.8, 6), "0.975quant" = rep(1.2, 6), check.names = FALSE
  )
  mock_fit_st$summary.linear.predictor <- data.frame(mean = rep(1.0, 6))
  res2 <- fastsae:::.extract_inla_results(
    fit = mock_fit_st,
    data = data.frame(y = 1:6),
    y = 1:6,
    domain = rep(1:3, each = 2),
    time = rep(1:2, 3),
    family = "gaussian",
    spatial = "bym2",
    temporal = "ar1",
    st_interaction = "domain-specific",
    domain_id = rep(1:3, each = 2),
    time_id = rep(1:2, 3),
    unique_domains = 1:3,
    unique_times = 1:2
  )
  expect_s3_class(res2, "fastsae_hb_area")
  expect_true("random_effect_temporal" %in% names(res2$df_hb))

  # 3. Fallback to first random effect
  mock_fit_fb <- mock_fit
  mock_fit_fb$summary.random <- list(other_rf = data.frame(mean = c(0.1, 0.2, 0.3)))
  res3 <- fastsae:::.extract_inla_results(
    fit = mock_fit_fb,
    data = data.frame(y = 1:3),
    y = 1:3,
    domain = 1:3,
    family = "gaussian",
    spatial = "none",
    domain_id = 1:3
  )
  expect_s3_class(res3, "fastsae_hb_area")

  # 4. SLM spatial model with X_mat and rho_range
  mock_fit_slm <- list(
    summary.fixed = NULL,
    summary.hyperpar = data.frame(
      mean = c(2.0, 0.4),
      row.names = c("Precision for ..domain_id..", "Rho for ..domain_id..")
    ),
    summary.random = list(
      slm_effect = data.frame(
        mean = c(0.1, 0.2, 0.3, 1.5, 2.5),
        sd = c(0.1, 0.1, 0.1, 0.2, 0.3),
        "0.025quant" = c(0, 0, 0, 1.1, 1.9),
        "0.975quant" = c(0, 0, 0, 1.9, 3.1),
        check.names = FALSE
      )
    ),
    summary.fitted.values = data.frame(
      mean = c(1, 2, 3), sd = c(0.1, 0.1, 0.1), "0.025quant" = c(0.8, 1.8, 2.8), "0.975quant" = c(1.2, 2.2, 3.2), check.names = FALSE
    ),
    summary.linear.predictor = data.frame(mean = c(1, 2, 3))
  )
  X_mat <- matrix(1:6, nrow = 3, dimnames = list(NULL, c("Intercept", "x1")))
  res_slm <- fastsae:::.extract_inla_results(
    fit = mock_fit_slm,
    data = data.frame(y = 1:3),
    y = 1:3,
    domain = 1:3,
    family = "gaussian",
    spatial = "slm",
    X_mat = X_mat,
    rho_range = c(-1, 1)
  )
  expect_s3_class(res_slm, "fastsae_hb_area")
  expect_equal(nrow(res_slm$estcoef), 2)

  # 5. Beta family with trials (vardir = NULL) and Gamma with vardir
  res_beta <- fastsae:::.extract_inla_results(
    fit = mock_fit,
    data = data.frame(y = c(0.2, 0.4, 0.6)),
    y = c(0.2, 0.4, 0.6),
    domain = 1:3,
    family = "beta",
    spatial = "none",
    trials = c(10, 20, 30)
  )
  expect_true("precision" %in% names(res_beta$df_hb))

  res_gamma <- fastsae:::.extract_inla_results(
    fit = mock_fit,
    data = data.frame(y = 1:3),
    y = 1:3,
    domain = 1:3,
    family = "gamma",
    spatial = "none",
    vardir = c(0.1, 0.2, 0.3)
  )
  expect_true("cv_dir" %in% names(res_gamma$df_hb))

  # 6. Exposure provided for Poisson
  res_pois <- fastsae:::.extract_inla_results(
    fit = mock_fit,
    data = data.frame(y = 1:3),
    y = 1:3,
    domain = 1:3,
    family = "poisson",
    spatial = "none",
    exposure = c(100, 200, 300)
  )
  expect_true("estimated_count" %in% names(res_pois$df_hb))
})

test_that("inla_utils helper checks and random effect extraction fallbacks work", {
  # 1. .check_inla_installed when INLA is missing
  testthat::with_mocked_bindings(
    requireNamespace = function(package, ...) {
      if (package == "INLA") return(FALSE)
      TRUE
    },
    .package = "base",
    {
      expect_error(fastsae:::.check_inla_installed(), "Package INLA is required")
    }
  )

  # 2. .convert_spatial_weights when spdep is missing for listw / nb
  skip_if_not_installed("spdep")
  W_mat <- matrix(c(0,1,0,1,0,1,0,1,0), 3, 3)
  listw_obj <- spdep::mat2listw(W_mat, style = "W")
  nb_obj <- spdep::mat2listw(W_mat, style = "B")$neighbours
  testthat::with_mocked_bindings(
    requireNamespace = function(package, ...) {
      if (package == "spdep") return(FALSE)
      TRUE
    },
    .package = "base",
    {
      expect_error(fastsae:::.convert_spatial_weights(listw_obj, 3, spatial = "bym2"), "spdep is required to convert")
      expect_error(fastsae:::.convert_spatial_weights(nb_obj, 3, spatial = "bym2"), "spdep is required to convert")
    }
  )

  # 3. .extract_inla_results with temporal random effect only (line 258)
  mock_fit_temp <- list(
    summary.fixed = data.frame(mean = 1, sd = 0.1, "0.025quant" = 0.8, "0.975quant" = 1.2, row.names = "(Intercept)", check.names = FALSE),
    summary.hyperpar = data.frame(mean = 2, row.names = "Precision for ..time_id.."),
    summary.random = list(`..time_id..` = data.frame(mean = c(0.1, 0.2, 0.3))),
    summary.fitted.values = data.frame(mean = c(1, 2, 3), sd = c(0.1, 0.1, 0.1), "0.025quant" = c(0.8, 1.8, 2.8), "0.975quant" = c(1.2, 2.2, 3.2), check.names = FALSE),
    summary.linear.predictor = data.frame(mean = c(1, 2, 3), sd = c(0.1, 0.1, 0.1), "0.025quant" = c(0.8, 1.8, 2.8), "0.975quant" = c(1.2, 2.2, 3.2), check.names = FALSE)
  )
  res_t <- fastsae:::.extract_inla_results(
    fit = mock_fit_temp, data = data.frame(y = 1:3), y = 1:3, domain = 1:3,
    family = "gaussian", spatial = "none", temporal = "rw1", time = 1:3
  )
  expect_equal(res_t$df_hb$random_effect, c(0.1, 0.2, 0.3))

  # 4. Fallback for simple IID with custom precision name (line 161) and domain_id mapping (line 265)
  mock_fit_custom <- list(
    summary.fixed = data.frame(mean = 1, sd = 0.1, "0.025quant" = 0.8, "0.975quant" = 1.2, row.names = "(Intercept)", check.names = FALSE),
    summary.hyperpar = data.frame(mean = 4, row.names = "Precision for custom_effect"),
    summary.random = list(custom_effect = data.frame(mean = c(0.5, 0.6))),
    summary.fitted.values = data.frame(mean = 1:4, sd = rep(0.1, 4), "0.025quant" = rep(0.8, 4), "0.975quant" = rep(1.2, 4), check.names = FALSE),
    summary.linear.predictor = data.frame(mean = 1:4, sd = rep(0.1, 4), "0.025quant" = rep(0.8, 4), "0.975quant" = rep(1.2, 4), check.names = FALSE)
  )
  res_c <- fastsae:::.extract_inla_results(
    fit = mock_fit_custom, data = data.frame(y = 1:4), y = 1:4, domain = c(1, 1, 2, 2),
    family = "gaussian", spatial = "none", temporal = "none",
    domain_id = c(1, 1, 2, 2)
  )
  expect_equal(res_c$random_effect_var, 1 / 4)
  expect_equal(res_c$df_hb$random_effect, c(0.5, 0.5, 0.6, 0.6))

  # 5. Invalid length fallback to NA (line 271)
  mock_fit_na <- mock_fit_custom
  mock_fit_na$summary.random <- list(custom_effect = data.frame(mean = 0.5))
  res_na <- fastsae:::.extract_inla_results(
    fit = mock_fit_na, data = data.frame(y = 1:4), y = 1:4, domain = c(1, 1, 2, 2),
    family = "gaussian", spatial = "none", temporal = "none"
  )
  expect_true(all(is.na(res_na$df_hb$random_effect)))
})


