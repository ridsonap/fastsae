library(testthat)
library(fastsae)

skip_if_not_installed("sae")

# ------------------------------------------------------------------
# Shared objects
# ------------------------------------------------------------------

idx <- !is.na(mys$y)

mysnona <- mys[idx, ]
mys_proxmat_nona <- mys_proxmat[idx, idx]

fit_fast <- eblup_sfh(
  y ~ x1 + x2 + x3,
  vardir = "vardir",
  domain = "area",
  method = "REML",
  data = mysnona,
  W = mys_proxmat_nona,
  print_result = FALSE
)

fit_sae <- sae::mseSFH(
  mysnona$y ~ mysnona$x1 + mysnona$x2 + mysnona$x3,
  vardir = mysnona$vardir,
  proxmat = mys_proxmat_nona,
  method = "REML"
)


fit_fast_ml <- eblup_sfh(
  y ~ x1 + x2 + x3,
  vardir = "vardir",
  domain = "area",
  method = "ML",
  data = mysnona,
  W = mys_proxmat_nona,
  print_result = FALSE
)

fit_sae_ml <- sae::mseSFH(
  mysnona$y ~ mysnona$x1 + mysnona$x2 + mysnona$x3,
  vardir = mysnona$vardir,
  proxmat = mys_proxmat_nona,
  method = "ML"
)

tol <- 1e-6


# -----------------------------------------------------------------
# Structure
# ------------------------------------------------------------------

test_that("eblup_sfh returns valid structure", {
  expect_true(is.list(fit_fast))

  # Standardized fields
  expect_true("estcoef" %in% names(fit_fast))
  expect_true("random_effect_var" %in% names(fit_fast))
  expect_true("rho" %in% names(fit_fast))
  expect_true("estvarcomp" %in% names(fit_fast))
  expect_true("goodness" %in% names(fit_fast))
  expect_true("df_eblup" %in% names(fit_fast))
  expect_true("model" %in% names(fit_fast))
  expect_true("level" %in% names(fit_fast))
  expect_true("n_iter" %in% names(fit_fast))
  expect_true("convergence" %in% names(fit_fast))
  expect_true("method" %in% names(fit_fast))

  # Check estvarcomp structure
  expect_true("parameter" %in% names(fit_fast$estvarcomp))
  expect_true("estimate" %in% names(fit_fast$estvarcomp))
  expect_true("std.error" %in% names(fit_fast$estvarcomp))

  # Check values
  expect_equal(fit_fast$model, "SFH")
  expect_equal(fit_fast$level, "area")
  expect_true(fit_fast$convergence)
  expect_true(fit_fast$n_iter > 0)

  expect_length(
    fit_fast$df_eblup$eblup,
    nrow(mysnona)
  )

  expect_length(
    fit_fast$df_eblup$mse,
    nrow(mysnona)
  )
})

# ------------------------------------------------------------------
# EBLUP
# ------------------------------------------------------------------

test_that("EBLUP agrees with sae::mseSFH", {
  expect_equal(
    fit_fast$df_eblup$eblup,
    as.numeric(fit_sae$est$eblup),
    tolerance = tol
  )

  expect_equal(
    fit_fast_ml$df_eblup$eblup,
    as.numeric(fit_sae_ml$est$eblup),
    tolerance = tol
  )
})

# ------------------------------------------------------------------
# MSE
# ------------------------------------------------------------------

test_that("MSE agrees with sae::mseSFH", {
  expect_equal(
    fit_fast$df_eblup$mse,
    fit_sae$mse,
    tolerance = tol
  )

  expect_equal(
    fit_fast_ml$df_eblup$mse,
    fit_sae_ml$mse,
    tolerance = tol
  )
})

# ------------------------------------------------------------------
# Variance component
# ------------------------------------------------------------------

test_that("random effect variance agrees with sae::mseSFH", {
  expect_equal(
    fit_fast$random_effect_var,
    fit_sae$est$fit$refvar,
    tolerance = tol
  )

  expect_equal(
    fit_fast_ml$random_effect_var,
    fit_sae_ml$est$fit$refvar,
    tolerance = tol
  )
})

test_that("rho parameter is in reasonable range", {
  # Spatial autocorrelation should be in (-1, 1)
  expect_true(fit_fast$rho >= -1 && fit_fast$rho <= 1)
})

# ------------------------------------------------------------------
# Goodness of fit
# ------------------------------------------------------------------

test_that("goodness statistics agree with sae::mseSFH", {
  expect_equal(
    as.numeric(fit_fast$goodness),
    as.numeric(fit_sae$est$fit$goodness[-4]),
    tolerance = tol
  )

  expect_equal(
    as.numeric(fit_fast_ml$goodness),
    as.numeric(fit_sae_ml$est$fit$goodness[-4]),
    tolerance = tol
  )
})


# ------------------------------------------------------------------
# Regression coefficients
# ------------------------------------------------------------------

test_that("beta estimates agree with sae::mseSFH", {
  expect_identical(
    names(fit_fast$estcoef$beta),
    names(fit_sae$est$fit$estcoef$beta)
  )

  expect_equal(
    fit_fast$estcoef$beta,
    fit_sae$est$fit$estcoef$beta,
    tolerance = tol
  )

  expect_identical(
    names(fit_fast_ml$estcoef$beta),
    names(fit_sae_ml$est$fit$estcoef$beta)
  )

  expect_equal(
    fit_fast_ml$estcoef$beta,
    fit_sae_ml$est$fit$estcoef$beta,
    tolerance = tol
  )
})

# ------------------------------------------------------------------
# Standard errors
# ------------------------------------------------------------------

test_that("beta standard errors agree with sae::mseSFH", {
  expect_equal(
    fit_fast$estcoef$stderr_beta,
    fit_sae$est$fit$estcoef$std.error,
    tolerance = tol
  )

  expect_equal(
    fit_fast_ml$estcoef$stderr_beta,
    fit_sae_ml$est$fit$estcoef$std.error,
    tolerance = tol
  )
})

# ------------------------------------------------------------------
# Error handling
# ------------------------------------------------------------------

test_that("non-square proximity matrix throws error", {
  W_bad <- mys_proxmat_nona[-1, ]

  expect_error(
    eblup_sfh(
      y ~ x1 + x2 + x3,
      vardir = "vardir",
      data = mysnona,
      W = W_bad,
      print_result = FALSE
    )
  )
})

test_that("negative sampling variance throws error", {
  dat_bad <- mysnona
  dat_bad$vardir[1] <- -1

  expect_error(
    eblup_sfh(
      y ~ x1 + x2 + x3,
      vardir = "vardir",
      data = dat_bad,
      W = mys_proxmat_nona,
      print_result = FALSE
    )
  )
})
