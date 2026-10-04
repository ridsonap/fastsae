library(fastsae)

test_that("gvf_smooth works with log_linear method", {
  res_ll <- gvf_smooth(
    vardir = "vardir",
    y = "y",
    n = "n",
    data = mys,
    method = "log_linear",
    bias_correction = TRUE
  )

  sampled_idx <- which(!is.na(mys$y) & mys$vardir > 0)
  expect_s3_class(res_ll, "gvf_smooth")
  expect_length(res_ll$smooth_vardir, nrow(mys))
  expect_true(all(res_ll$smooth_vardir[sampled_idx] > 0))
  expect_true(res_ll$r_squared >= 0 && res_ll$r_squared <= 1)
  expect_equal(res_ll$method, "log_linear")
  expect_true("smooth_vardir" %in% names(res_ll$df))

  # Test without bias correction
  res_nobc <- gvf_smooth(
    vardir = mys$vardir,
    y = mys$y,
    n = mys$n,
    method = "log_linear",
    bias_correction = FALSE
  )
  expect_true(all(res_nobc$smooth_vardir[sampled_idx] < res_ll$smooth_vardir[sampled_idx]))
})

test_that("gvf_smooth works with gamma_glm method", {
  res_gamma <- gvf_smooth(
    vardir = "vardir",
    y = "y",
    n = "n",
    data = mys,
    method = "gamma_glm"
  )

  sampled_idx <- which(!is.na(mys$y) & mys$vardir > 0)
  expect_s3_class(res_gamma, "gvf_smooth")
  expect_length(res_gamma$smooth_vardir, nrow(mys))
  expect_true(all(res_gamma$smooth_vardir[sampled_idx] > 0))
  expect_equal(res_gamma$method, "gamma_glm")
  expect_s3_class(res_gamma$fitted_model, "glm")
})

test_that("gvf_smooth works with cv_power method", {
  res_cv <- gvf_smooth(
    vardir = "vardir",
    y = "y",
    n = "n",
    data = mys,
    method = "cv_power"
  )

  sampled_idx <- which(!is.na(mys$y) & mys$vardir > 0)
  expect_s3_class(res_cv, "gvf_smooth")
  expect_length(res_cv$smooth_vardir, nrow(mys))
  expect_true(all(res_cv$smooth_vardir[sampled_idx] > 0))
  expect_equal(res_cv$method, "cv_power")
})

test_that("gvf_smooth integrates seamlessly with eblup_fh", {
  gvf_fit <- gvf_smooth(
    vardir = "vardir",
    y = "y",
    n = "n",
    data = mys,
    method = "log_linear"
  )

  # Pass gvf_fit object directly to eblup_fh
  fit_fh <- eblup_fh(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = gvf_fit
  )

  expect_s3_class(fit_fh, "fastsae")
  expect_true(fit_fh$convergence)
  expect_equal(fit_fh$df_eblup$vardir, gvf_fit$smooth_vardir)
  expect_true(all(!is.na(fit_fh$df_eblup$eblup)))
  expect_true(all(fit_fh$df_eblup$mse > 0))
})

test_that("gvf_smooth S3 methods and error handling work", {
  gvf_fit <- gvf_smooth(
    vardir = "vardir",
    y = "y",
    data = mys,
    method = "log_linear"
  )

  # S3 output tests
  expect_output(print(gvf_fit), "Generalized Variance Function")
  expect_output(summary(gvf_fit), "Model Coefficients")

  # Error handling: non-positive vardir
  expect_error(
    gvf_smooth(vardir = c(1, -2, 3), y = c(1, 2, 3)),
    "strictly positive"
  )

  # Error handling: missing both y and n
  expect_error(
    gvf_smooth(vardir = c(1, 2, 3)),
    "At least one explanatory variable"
  )

  # Error handling: length mismatch
  expect_error(
    gvf_smooth(vardir = c(1, 2, 3), y = c(1, 2)),
    "must match length of 'vardir'"
  )
})
