library(testthat)
library(fastsae)

test_that("ebp_area Gaussian estimates match sae::mseFH benchmark (r > 0.99)", {
  skip_if_not_installed("sae")

  data("mys", package = "fastsae")
  mysnona <- mys[!is.na(mys[["y"]]), ]

  fit_sae <- sae::mseFH(
    mysnona[["y"]] ~ mysnona[["x1"]] + mysnona[["x2"]] + mysnona[["x3"]],
    vardir = mysnona[["vardir"]],
    method = "REML"
  )

  fit_ebp <- ebp_area(
    y ~ x1 + x2 + x3,
    data = mysnona,
    vardir = "vardir",
    domain = "area",
    family = "gaussian",
    spatial = "none",
    print_result = FALSE
  )

  sae_eblup <- as.numeric(fit_sae[["est"]][["eblup"]])
  ebp_pred <- as.numeric(fit_ebp[["df_ebp"]][["ebp"]])

  # Correlation between Frequentist EBLUP (REML) and Bayesian Posterior Mean (INLA)
  r_est <- cor(sae_eblup, ebp_pred)
  expect_gt(r_est, 0.99)

  # Average difference should be small (< 0.25 on this scale)
  expect_lt(mean(abs(sae_eblup - ebp_pred)), 0.25)

  # MSE / Posterior Variance correlation
  sae_mse <- as.numeric(fit_sae[["mse"]])
  ebp_mse <- as.numeric(fit_ebp[["df_ebp"]][["mse"]])
  r_mse <- cor(sae_mse, ebp_mse)
  expect_gt(r_mse, 0.98)
})

test_that("ebp_area achieves parameter recovery and nominal coverage rate on synthetic data", {
  sim <- sim_area_data(
    D = 40,
    beta = c(10.0, 1.0, 0.5),
    sigma_u = 0.5,
    n_unsampled = 0,
    seed = 2026
  )

  fit_sim <- ebp_area(
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
  data("mys", package = "fastsae")

  fit_ebp <- ebp_area(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = ~vardir,
    domain = ~area,
    family = "gaussian",
    print_result = FALSE
  )

  d_ebp <- diagnose(fit_ebp)
  expect_s3_class(d_ebp, "fastsae_diagnose")

  # Bayesian metrics slot
  expect_true(!is.null(d_ebp$bayesian_metrics))
  expect_true(is.numeric(d_ebp$bayesian_metrics$waic))
  expect_true(is.numeric(d_ebp$bayesian_metrics$dic))
  expect_true(is.numeric(d_ebp$bayesian_metrics$log_mlik))
  expect_true(!is.null(d_ebp$bayesian_metrics$cpo_summary))
  expect_true(!is.null(d_ebp$bayesian_metrics$pit_test))

  # Dataframe includes cpo and pit columns
  expect_true("cpo" %in% names(d_ebp$df_diag))
  expect_true("pit" %in% names(d_ebp$df_diag))

  # Autoplot methods for Bayesian calibration
  p_pit <- autoplot(d_ebp, type = "pit")
  expect_s3_class(p_pit, "ggplot")

  p_cpo <- autoplot(d_ebp, type = "cpo")
  expect_s3_class(p_cpo, "ggplot")

  # Print report contains Bayesian section
  expect_message(print(d_ebp), "Bayesian Information Criteria")
})
