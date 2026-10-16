library(testthat)
library(fastsae)

test_that("diagnose works for eblup_fh objects", {
  data("mys", package = "fastsae")

  fit_fh <- eblup_fh(y ~ x1 + x2 + x3, data = mys, vardir = ~vardir, domain = ~area, print_result = FALSE)
  d_fh <- diagnose(fit_fh, rse_threshold = 25)

  expect_s3_class(d_fh, "fastsae_diagnose")
  expect_named(d_fh, c("model_info", "precision", "brown_test", "spatial_test",
                       "bayesian_metrics", "simulation_metrics", "df_diag", "status"))

  # Model info
  expect_equal(d_fh$model_info$N_total, nrow(mys))
  expect_equal(d_fh$model_info$N_sampled, sum(!is.na(mys$y)))

  # Precision
  expect_true(d_fh$precision$prop_reliable >= 0 && d_fh$precision$prop_reliable <= 100)
  expect_true(!is.na(d_fh$precision$mean_sae_rse))
  expect_true(d_fh$precision$mean_sae_rse < d_fh$precision$mean_direct_rse)
  expect_true(!is.null(d_fh$precision$eff_ratio_summary))

  # Brown diagnostics: descriptive bias regression + upper-tail W test
  expect_true(!is.null(d_fh$brown_test$bias_regression))
  expect_true(!is.null(d_fh$brown_test$goodness_of_fit))
  expect_true(is.numeric(d_fh$brown_test$bias_regression$beta))
  expect_null(d_fh$brown_test$bias_regression$p_value)
  gof <- d_fh$brown_test$goodness_of_fit
  expect_equal(gof$p_value, pchisq(gof$w_stat, gof$df, lower.tail = FALSE))
  expect_identical(gof$is_good_fit, gof$p_value >= 0.05)

  # Print method
  expect_message(print(d_fh), "Diagnostic Report")
})

test_that("diagnose works for eblup_sfh with spatial autocorrelation test", {
  data("mys", package = "fastsae")
  data("mys_proxmat", package = "fastsae")

  fit_sfh <- eblup_sfh(y ~ x1 + x2 + x3, data = mys, vardir = ~vardir,
                       domain = ~area, W = mys_proxmat, print_result = FALSE)

  # diagnose should automatically detect W from fit_sfh
  d_sfh <- diagnose(fit_sfh)
  expect_s3_class(d_sfh, "fastsae_diagnose")

  # Spatial Moran test on residuals
  expect_false(is.null(d_sfh$spatial_test))
  expect_true(is.numeric(d_sfh$spatial_test$moran_I))
  expect_true(is.numeric(d_sfh$spatial_test$p_value))
  expect_true(d_sfh$spatial_test$p_value >= 0 && d_sfh$spatial_test$p_value <= 1)
})

test_that("diagnose works for hb_area and handles simulation truth", {
  skip_if_not(requireNamespace("INLA", quietly = TRUE), "INLA not installed")

  data("sim_area", package = "fastsae")
  data("mys_proxmat", package = "fastsae")

  fit_hb <- hb_area(y_gaussian ~ x1 + x2, data = sim_area,
                      vardir = "vardir", spatial = "bym2", W = mys_knn,
                      print_result = FALSE)

  # Supply dummy ground truth for simulation verification
  mock_truth <- fit_hb$df_hb$hb * runif(nrow(sim_area), 0.98, 1.02)
  d_hb <- diagnose(fit_hb, truth = mock_truth)

  expect_s3_class(d_hb, "fastsae_diagnose")
  expect_false(is.null(d_hb$simulation_metrics))
  expect_true(is.numeric(d_hb$simulation_metrics$mean_abs_rb))
  expect_true(is.numeric(d_hb$simulation_metrics$rrmse))
  err <- fit_hb$df_hb$hb - mock_truth
  expect_equal(d_hb$simulation_metrics$rrmse, 100 * sqrt(mean(err^2)) / mean(abs(mock_truth)))
})

test_that("autoplot.fastsae_diagnose generates valid ggplot objects", {
  data("mys", package = "fastsae")
  fit_fh <- eblup_fh(y ~ x1 + x2 + x3, data = mys, vardir = ~vardir, domain = ~area, print_result = FALSE)
  d_fh <- diagnose(fit_fh)

  # Test each plot type
  p_all <- autoplot(d_fh, type = "all")
  expect_s3_class(p_all, "ggplot")

  p_calib <- autoplot(d_fh, type = "calibration")
  expect_s3_class(p_calib, "ggplot")

  p_rse <- autoplot(d_fh, type = "rse")
  expect_s3_class(p_rse, "ggplot")

  p_res <- autoplot(d_fh, type = "residuals")
  expect_s3_class(p_res, "ggplot")

  p_qq <- autoplot(d_fh, type = "qq")
  expect_s3_class(p_qq, "ggplot")
})

test_that("diagnose error handling and input validation", {
  # 1. Non-fastsae object
  expect_error(diagnose(list(a = 1)), "must be a fitted model of class")

  # 2. fastsae object missing estimation tables
  dummy_obj <- structure(list(model = "FH"), class = "fastsae")
  expect_error(diagnose(dummy_obj), "Could not find estimation results table")

  # 3. Table missing prediction column
  dummy_obj2 <- structure(list(df_eblup = data.frame(a = 1:5)), class = "fastsae")
  expect_error(diagnose(dummy_obj2), "Prediction column")
})

test_that("diagnose handles custom W, truth without CI, and status branches", {
  data("mys", package = "fastsae")
  data("mys_proxmat", package = "fastsae")

  fit_fh <- eblup_fh(y ~ x1 + x2 + x3, data = mys, vardir = ~vardir, domain = ~area, print_result = FALSE)

  # Supply W explicitly
  d_w <- diagnose(fit_fh, W = mys_proxmat)
  expect_false(is.null(d_w$spatial_test))

  # Supply truth without CI in df_eblup
  fit_no_ci <- fit_fh
  fit_no_ci$df_eblup$ci_lower <- NULL
  fit_no_ci$df_eblup$ci_upper <- NULL
  truth_vals <- fit_fh$df_eblup$eblup + rnorm(nrow(mys), sd = 0.1)
  d_truth <- diagnose(fit_no_ci, truth = truth_vals)
  expect_true(is.na(d_truth$simulation_metrics$coverage_rate))
  expect_true(is.numeric(d_truth$simulation_metrics$mean_rb))

  # Status caution and warning branches
  # 1. Force low precision: set rse_threshold very low so prop_reliable < 80
  d_low_prec <- diagnose(fit_fh, rse_threshold = 1)
  expect_match(d_low_prec$status, "WARNING|CAUTION")

  # 2. Check print output for low precision
  expect_message(print(d_low_prec), "Caution: low precision")
})

test_that("diagnose and print.fastsae_diagnose cover all alert branches and Bayesian components", {
  data("mys", package = "fastsae")
  fit_fh <- eblup_fh(y ~ x1 + x2 + x3, data = mys, vardir = ~vardir, domain = ~area, print_result = FALSE)

  # Build a diagnostic object with specific components to exercise all print paths
  d_full <- diagnose(fit_fh)

  # Branch 1: prop_reliable < 80, diff_rse <= 0, prop_gain < 80
  d_full$precision$prop_reliable <- 50
  d_full$precision$mean_direct_rse <- 10
  d_full$precision$mean_sae_rse <- 15
  d_full$precision$prop_gain <- 60
  d_full$precision$eff_ratio_summary <- c(Min = 0.5, Q1 = 0.7, Median = 0.9, Mean = 0.9, Q3 = 1.1, Max = 1.5)

  # Goodness of fit poor
  d_full$brown_test$goodness_of_fit$is_good_fit <- FALSE

  # Spatial residual autocorrelation present
  d_full$spatial_test <- list(
    moran_I = 0.45,
    expected_I = -0.05,
    sd_I = 0.1,
    z_stat = 5.0,
    p_value = 0.0001,
    no_residual_autocorrelation = FALSE
  )

  # Simulation metrics with high bias and low coverage
  d_full$simulation_metrics <- list(
    mean_rb = 12.5,
    mean_abs_rb = 12.5,
    rrmse = 18.2,
    coverage_rate = 75.0
  )

  # Bayesian metrics with failure > 0 and uncalibrated pit
  d_full$bayesian_metrics <- list(
    waic = 210.5,
    pWAIC = 12.1,
    dic = 208.3,
    pD = 11.4,
    log_mlik = -105.2,
    cpo_summary = list(mean_cpo = 0.12, min_cpo = 0.001, failures = 3),
    pit_test = list(d_stat = 0.35, p_value = 0.002, is_calibrated = FALSE)
  )

  # Status with CAUTION
  d_full$status <- "CAUTION: Model estimates differ from the direct estimates more than their uncertainty allows (Brown W test)."
  expect_message(print(d_full), "Differs from direct estimates")
  expect_message(print(d_full), "Bias regression")
  expect_message(print(d_full), "Residual Spatial Pattern Remains")
  expect_message(print(d_full), "Below nominal 90%")
  expect_message(print(d_full), "Potential dispersion/skewness deviation")
  expect_message(print(d_full), "flagged with numerical caution")
  expect_message(print(d_full), "CAUTION")

  # Status with WARNING
  d_full$status <- "WARNING: High proportions of domains with RSE >= threshold."
  expect_message(print(d_full), "WARNING")

  # Bias regression is NULL
  d_full$brown_test$bias_regression <- NULL
  expect_message(print(d_full), "insufficient sampled observations")
})

test_that("autoplot.fastsae_diagnose pit and cpo plots work and error when missing", {
  data("mys", package = "fastsae")
  fit_fh <- eblup_fh(y ~ x1 + x2 + x3, data = mys, vardir = ~vardir, domain = ~area, print_result = FALSE)
  d_fh <- diagnose(fit_fh)

  # When pit and cpo are missing in df_diag
  expect_error(autoplot(d_fh, type = "pit"), "PIT .* not available")
  expect_error(autoplot(d_fh, type = "cpo"), "CPO .* not available")

  # Add mock pit and cpo to df_diag
  d_fh$df_diag$pit <- runif(nrow(d_fh$df_diag), 0.05, 0.95)
  d_fh$df_diag$cpo <- runif(nrow(d_fh$df_diag), 0.1, 0.8)

  p_pit <- autoplot(d_fh, type = "pit")
  expect_s3_class(p_pit, "ggplot")

  p_cpo <- autoplot(d_fh, type = "cpo")
  expect_s3_class(p_cpo, "ggplot")
})

test_that("diagnose handles a minimal object and fully passing metrics print", {
  # 1. Minimal object
  obj_ebp <- structure(
    list(
      model = "EBP Model",
      df_eblup = data.frame(domain = 1:5, y = c(1, 2, 3, 4, 5), eblup = c(1.1, 1.9, 3.1, 3.9, 5.0), vardir = rep(0.1, 5), mse = rep(0.01, 5))
    ),
    class = "fastsae"
  )
  d_ebp <- diagnose(obj_ebp, truth = 1:5)
  expect_s3_class(d_ebp, "fastsae_diagnose")

  # 2. Perfect passing status & print
  d_pass <- d_ebp
  d_pass$status <- "PASS: Model estimates are consistent with the direct estimates (Brown W test) and reliable."
  d_pass$precision$prop_reliable <- 95
  d_pass$precision$prop_gain <- 85
  d_pass$precision$mean_direct_rse <- 25
  d_pass$precision$mean_sae_rse <- 10
  d_pass$brown_test$goodness_of_fit <- list(is_good_fit = TRUE, w_stat = 3.2, df = 5, p_value = 0.6)
  d_pass$spatial_test <- list(no_residual_autocorrelation = TRUE, moran_I = 0.01, expected_I = -0.05, p_value = 0.7)
  d_pass$simulation_metrics <- list(mean_rb = 1.2, mean_abs_rb = 3.1, rrmse = 4.5, coverage_rate = 94.5)

  msgs <- testthat::capture_messages(print(d_pass))
  combined <- paste(msgs, collapse = " ")
  expect_true(grepl("Final Assessment.*PASS", combined))
  expect_true(grepl("High precision", combined))
  expect_true(grepl("Consistent with direct estimates", combined))
  expect_true(grepl("No Residual Spatial Autocorrelation", combined))
  expect_true(grepl("Low bias < 5%", combined))
  expect_true(grepl("Nominal target met", combined))
})

test_that("diagnose Brown W test keeps its size and is not driven by shrinkage", {
  # Correct FH model without spatial effects: the upper-tail W test must not
  # reject, and the status must not report bias (the old F test of alpha = 0,
  # beta = 1 rejected ~93% of such fits; the old two-sided W test 100%).
  set.seed(11)
  D <- 80
  rej <- vapply(1:40, function(r) {
    x1 <- rnorm(D)
    psi <- runif(D, 0.3, 1.5)
    dat <- data.frame(y = 1 + x1 + rnorm(D, sd = sqrt(0.6)) + rnorm(D, sd = sqrt(psi)), x1 = x1, psi = psi)
    fit <- eblup_fh(y ~ x1, vardir = "psi", data = dat, print_result = FALSE)
    d <- diagnose(fit)
    c(gof = !d$brown_test$goodness_of_fit$is_good_fit, caution = grepl("^CAUTION", d$status),
      beta = d$brown_test$bias_regression$beta)
  }, numeric(3))
  expect_lte(mean(rej["gof", ]), 0.1)
  expect_equal(mean(rej["caution", ]), mean(rej["gof", ]))
  # the descriptive slope is above 1 under shrinkage, as documented
  expect_gt(mean(rej["beta", ]), 1)
})

test_that("diagnose reads Bayesian criteria of hb_unit-style goodness lists", {
  obj <- structure(
    list(
      model = "HB-UNIT",
      df_hb = data.frame(domain = 1:6, y = c(1, 2, 3, 4, 5, 6), hb = c(1.1, 1.9, 3.1, 3.9, 5.1, 6.0),
                         vardir = rep(0.2, 6), mse = rep(0.05, 6)),
      goodness = list(dic = 101.5, p_eff_dic = 4.2, waic = 102.3, p_eff_waic = 4.8, mlik = -60.1),
      fit = list(cpo = list(cpo = runif(20), pit = runif(20)))  # unit rows: not per domain
    ),
    class = c("fastsae_hb_unit", "fastsae_hb", "fastsae")
  )
  d <- diagnose(obj)
  expect_equal(d$bayesian_metrics$dic, 101.5)
  expect_equal(d$bayesian_metrics$waic, 102.3)
  expect_equal(d$bayesian_metrics$pWAIC, 4.8)
  expect_equal(d$bayesian_metrics$log_mlik, -60.1)
  expect_null(d$bayesian_metrics$cpo_summary)   # CPO rows do not match the domains
  expect_null(d$df_diag$cpo)
})

test_that("diagnose standardises residuals by sqrt(vardir - MSE)", {
  set.seed(21)
  D <- 300
  x1 <- rnorm(D); psi <- runif(D, 0.3, 1.5)
  dat <- data.frame(y = 1 + x1 + rnorm(D, sd = sqrt(0.6)) + rnorm(D, sd = sqrt(psi)), x1 = x1, psi = psi)
  fit <- eblup_fh(y ~ x1, vardir = "psi", data = dat, print_result = FALSE)
  d <- diagnose(fit)$df_diag
  expect_equal(d$std_residuals, d$residuals / sqrt(d$vardir - d$mse))
  # approximately N(0, 1) under the correct model (was ~0.76 with sqrt(vardir))
  expect_gt(sd(d$std_residuals), 0.85)
  expect_lt(sd(d$std_residuals), 1.15)
})
