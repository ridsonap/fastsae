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

  # Brown calibration test
  expect_true(!is.null(d_fh$brown_test$wald_test))
  expect_true(!is.null(d_fh$brown_test$goodness_of_fit))
  expect_true(is.numeric(d_fh$brown_test$wald_test$p_value))

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

test_that("diagnose works for ebp_area and handles simulation truth", {
  skip_if_not(requireNamespace("INLA", quietly = TRUE), "INLA not installed")

  data("sae_area_multi", package = "fastsae")
  data("mys_proxmat", package = "fastsae")

  fit_ebp <- ebp_area(y_gaussian ~ x1 + x2, data = sae_area_multi,
                      vardir = "vardir", spatial = "bym2", W = mys_proxmat,
                      print_result = FALSE)

  # Supply dummy ground truth for simulation verification
  mock_truth <- fit_ebp$df_ebp$ebp * runif(nrow(sae_area_multi), 0.98, 1.02)
  d_ebp <- diagnose(fit_ebp, truth = mock_truth)

  expect_s3_class(d_ebp, "fastsae_diagnose")
  expect_false(is.null(d_ebp$simulation_metrics))
  expect_true(is.numeric(d_ebp$simulation_metrics$mean_abs_rb))
  expect_true(is.numeric(d_ebp$simulation_metrics$mean_rrmse))
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
