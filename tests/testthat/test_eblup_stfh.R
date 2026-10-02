library(testthat)
library(fastsae)

skip_if_not_installed("sae")

# ------------------------------------------------------------------
# Shared objects
# ------------------------------------------------------------------

library(dplyr)
mys_panel_nona <- mys_panel |>
  filter(!is.na(y))
mys_proxmat_nona <- mys_proxmat[-c(21, 25), -c(21, 25)]


fit_fast <- eblup_stfh(
  y ~ x1 + x2 + x3,
  data = mys_panel_nona,
  vardir = ~vardir,
  domain = ~area,
  time = ~year,
  W = mys_proxmat_nona,
  model = "ST",
  print_result = FALSE
)

fit_sae <- sae::eblupSTFH(
  mys_panel_nona$y ~ mys_panel_nona$x1 + mys_panel_nona$x2 + mys_panel_nona$x3,
  vardir = mys_panel_nona$vardir,
  D = length(unique(mys_panel_nona$area)),
  T = length(unique(mys_panel_nona$year)),
  proxmat = mys_proxmat_nona
)


tol <- 1e-6


# -----------------------------------------------------------------
# Structure
# ------------------------------------------------------------------

test_that("eblup_stfh returns valid structure matching seblup_area", {
  expect_true(is.list(fit_fast))

  # Top-level fields (matching seblup_area structure)
  expect_true("estcoef" %in% names(fit_fast))
  expect_true("goodness" %in% names(fit_fast))
  expect_true("df_eblup" %in% names(fit_fast))
  expect_true("convergence" %in% names(fit_fast))
  expect_true("n_iter" %in% names(fit_fast))

  # Additional for spatio-temporal
  expect_true("estvarcomp" %in% names(fit_fast))
  expect_true("model" %in% names(fit_fast))
  expect_true("formula" %in% names(fit_fast))

  # df_eblup has eblup column
  expect_true("eblup" %in% names(fit_fast$df_eblup))
  expect_true("random_effect_u1" %in% names(fit_fast$df_eblup))
  expect_true("random_effect_u2" %in% names(fit_fast$df_eblup))

  expect_length(fit_fast$df_eblup$eblup, nrow(mys_panel_nona))
})

# ------------------------------------------------------------------
# EBLUP
# ------------------------------------------------------------------

test_that("EBLUP agrees with sae::eblupSTFH", {
  expect_equal(
    fit_fast$df_eblup$eblup,
    as.numeric(fit_sae$eblup),
    tolerance = tol
  )
})

# ------------------------------------------------------------------
# Variance component (estvarcomp - spatio-temporal specific)
# ------------------------------------------------------------------

test_that("estvarcomp agrees with sae::eblupSTFH", {
  expect_equal(
    fit_fast$estvarcomp$estimate,
    as.numeric(fit_sae$fit$estvarcomp$estimate),
    tolerance = tol
  )
  expect_equal(
    fit_fast$estvarcomp$std.error,
    as.numeric(fit_sae$fit$estvarcomp$std.error),
    tolerance = tol
  )
})

# ------------------------------------------------------------------
# Goodness of fit
# ------------------------------------------------------------------

test_that("goodness statistics agree with sae::eblupSTFH", {
  expect_equal(
    as.numeric(fit_fast$goodness),
    as.numeric(fit_sae$fit$goodness),
    tolerance = tol
  )
})

# ------------------------------------------------------------------
# Regression coefficients
# ------------------------------------------------------------------

test_that("beta estimates agree with sae::eblupSTFH", {
  expect_identical(
    names(fit_fast$estcoef$beta),
    names(fit_sae$fit$estcoef$beta)
  )

  expect_equal(
    fit_fast$estcoef$beta,
    fit_sae$fit$estcoef$beta,
    tolerance = tol
  )
})

# ------------------------------------------------------------------
# Beta Standard errors
# ------------------------------------------------------------------

test_that("beta standard errors agree with sae::eblupSTFH", {
  expect_equal(
    fit_fast$estcoef$std.error,
    fit_sae$fit$estcoef$std.error,
    tolerance = tol
  )
})


# ------------------------------------------------------------------
# Error handling
# ------------------------------------------------------------------

test_that("non-square proximity matrix throws error", {
  W_bad <- mys_proxmat_nona[-1, ]

  expect_error(
    eblup_stfh(
      y ~ x1 + x2 + x3,
      data = mys_panel_nona,
      vardir = ~vardir,
      domain = ~area,
      time = ~year,
      W = W_bad,
      model = "ST"
    )
  )
})

test_that("negative sampling variance throws error", {
  dat_bad <- mys_panel_nona
  dat_bad$vardir[1] <- -1

  expect_error(
    eblup_stfh(
      y ~ x1 + x2 + x3,
      data = dat_bad,
      vardir = ~vardir,
      domain = ~area,
      time = ~year,
      W = mys_proxmat_nona,
      model = "ST"
    )
  )
})

test_that("eblup_stfh validation and edge cases", {
  # Unsampled areas (NA in response)
  expect_error(
    eblup_stfh(
      y ~ x1 + x2 + x3,
      data = mys_panel,
      vardir = ~vardir,
      domain = ~area,
      time = ~year,
      W = mys_proxmat,
      model = "ST",
      print_result = FALSE
    ),
    "does not support unsampled areas"
  )

  # W contains NA
  W_na <- mys_proxmat_nona
  W_na[1, 2] <- NA
  expect_error(
    eblup_stfh(
      y ~ x1 + x2 + x3,
      data = mys_panel_nona,
      vardir = ~vardir,
      domain = ~area,
      time = ~year,
      W = W_na,
      model = "ST",
      print_result = FALSE
    ),
    "Argument W contains NA"
  )

  # Dimension mismatch
  bad_panel <- mys_panel_nona[-1, ]
  expect_error(
    eblup_stfh(
      y ~ x1 + x2 + x3,
      data = bad_panel,
      vardir = ~vardir,
      domain = ~area,
      time = ~year,
      W = mys_proxmat_nona,
      model = "ST",
      print_result = FALSE
    ),
    "Dimensions mismatch"
  )

  # Custom starting values
  fit_custom_start <- eblup_stfh(
    y ~ x1 + x2 + x3,
    data = mys_panel_nona,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat_nona,
    model = "ST",
    sigma21_start = 1.0,
    sigma22_start = 0.5,
    rho1_start = 0.3,
    rho2_start = 0.4,
    print_result = FALSE
  )
  expect_s3_class(fit_custom_start, "fastsae")

  # Non-convergence with maxiter = 1
  fit_st_nc <- eblup_stfh(
    y ~ x1 + x2 + x3,
    data = mys_panel_nona,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat_nona,
    model = "ST",
    maxiter = 1,
    print_result = FALSE
  )
  expect_false(fit_st_nc$convergence)

  # print_result = TRUE with model = "S"
  expect_output(
    eblup_stfh(
      y ~ x1 + x2 + x3,
      data = mys_panel_nona,
      vardir = ~vardir,
      domain = ~area,
      time = ~year,
      W = mys_proxmat_nona,
      model = "S",
      print_result = TRUE
    )
  )

  # Parameter bounds errors
  expect_error(
    eblup_stfh(y ~ x1 + x2, data = mys_panel_nona, vardir = ~vardir, domain = ~area, time = ~year,
               W = mys_proxmat_nona, sigma21_start = -1),
    "sigma21_start must be >= 0"
  )
  expect_error(
    eblup_stfh(y ~ x1 + x2, data = mys_panel_nona, vardir = ~vardir, domain = ~area, time = ~year,
               W = mys_proxmat_nona, sigma22_start = -1),
    "sigma22_start must be >= 0"
  )
  expect_error(
    eblup_stfh(y ~ x1 + x2, data = mys_panel_nona, vardir = ~vardir, domain = ~area, time = ~year,
               W = mys_proxmat_nona, rho1_start = -1.5),
    "rho1_start must be in the interval"
  )
  expect_error(
    eblup_stfh(y ~ x1 + x2, data = mys_panel_nona, vardir = ~vardir, domain = ~area, time = ~year,
               W = mys_proxmat_nona, model = "ST", rho2_start = 1.5),
    "rho2_start must be in the interval"
  )

  # Parametric bootstrap MSE
  expect_output(
    fit_pb <- eblup_stfh(
      y ~ x1 + x2,
      data = mys_panel_nona,
      vardir = ~vardir,
      domain = ~area,
      time = ~year,
      W = mys_proxmat_nona,
      model = "ST",
      compute_mse = TRUE,
      B = 5,
      seed = 42,
      print_result = TRUE
    ),
    "Fixed Effects Coefficients"
  )
  expect_true("mse_pb" %in% names(fit_pb$df_eblup))
  expect_true(all(fit_pb$df_eblup$mse_pb >= 0))
  expect_equal(fit_pb$B, 5)

  # Direct nrow(mf) != length(vardir) check
  testthat::with_mocked_bindings(
    .get_variable = function(data, variable) c(1, 2),
    .package = "fastsae",
    {
      expect_error(
        eblup_stfh(
          y ~ x1 + x2,
          data = mys_panel_nona,
          vardir = ~vardir,
          domain = ~area,
          time = ~year,
          W = mys_proxmat_nona,
          print_result = FALSE
        ),
        "Length of 'vardir' must equal number of observations"
      )
    }
  )
})

