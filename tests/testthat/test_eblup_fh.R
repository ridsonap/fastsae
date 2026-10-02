library(testthat)
library(fastsae)

skip_if_not_installed("sae")

quiet <- function(expr) {
  result <- NULL
  invisible(capture.output(result <- suppressWarnings(suppressMessages(expr))))
  result
}

# ------------------------------------------------------------------
# Shared objects
# ------------------------------------------------------------------

mysnona <- mys[!is.na(mys$y), ]

quiet(
  expect_output(
    fit_fast <- eblup_fh(
      y ~ x1 + x2 + x3,
      vardir = "vardir",
      domain = "area",
      method = "REML",
      data = mysnona,
      print_result = TRUE
    )
  )
)

fit_sae <- sae::mseFH(
  mysnona$y ~ mysnona$x1 + mysnona$x2 + mysnona$x3,
  vardir = mysnona$vardir,
  method = "REML"
)


fit_fast_ml <- eblup_fh(
  y ~ x1 + x2 + x3,
  vardir = "vardir",
  domain = "area",
  method = "ML",
  data = mysnona,
  print_result = FALSE
)


fit_sae_ml <- sae::mseFH(
  mysnona$y ~ mysnona$x1 + mysnona$x2 + mysnona$x3,
  vardir = mysnona$vardir,
  method = "ML"
)

tol <- 1e-5


# Structure ---------------------------------------------------------------
test_that("eblup_fh returns valid structure", {
  expect_true(is.list(fit_fast))

  expect_true("df_eblup" %in% names(fit_fast))
  expect_true("random_effect_var" %in% names(fit_fast))
  expect_true("estcoef" %in% names(fit_fast))

  # output length
  expect_length(
    fit_fast$df_eblup$eblup,
    nrow(mysnona)
  )

  # mse output length
  expect_length(
    fit_fast$df_eblup$mse,
    nrow(mysnona)
  )
})



# EBLUP est ---------------------------------------------------------------
test_that("EBLUP agrees with sae::mseFH", {
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



# MSE ---------------------------------------------------------------------
test_that("MSE agrees with sae::mseFH", {
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


# Random effect variance --------------------------------------------------
test_that("random effect variance agrees with sae::mseFH", {
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


# Goodness of fit ---------------------------------------------------------
test_that("goodness statistics agree with sae::mseFH", {
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


# Regression coefficients and Standard errors -----------------------------
test_that("beta estimates agree with sae::mseFH", {
  expect_equal(
    fit_fast$estcoef$beta,
    fit_sae$est$fit$estcoef$beta,
    tolerance = tol
  )

  expect_equal(
    fit_fast$estcoef$stderr_beta,
    fit_sae$est$fit$estcoef$std.error,
    tolerance = tol
  )

  expect_equal(
    fit_fast_ml$estcoef$beta,
    fit_sae_ml$est$fit$estcoef$beta,
    tolerance = tol
  )

  expect_equal(
    fit_fast_ml$estcoef$stderr_beta,
    fit_sae_ml$est$fit$estcoef$std.error,
    tolerance = tol
  )
})


# ------------------------------------------------------------------
# Input validation
# ------------------------------------------------------------------

test_that("negative or zero vardir throws error", {
  dat_bad <- mysnona
  dat_bad$vardir[1] <- -1
  dat_bad$vardir[2] <- 0

  expect_error(
    eblup_fh(
      y ~ x1 + x2 + x3,
      vardir = "vardir",
      data = dat_bad,
      print_result = FALSE
    )
  )
})


# ------------------------------------------------------------------
# Method argument
# ------------------------------------------------------------------

test_that("invalid method throws error", {
  expect_error(
    eblup_fh(
      y ~ x1 + x2 + x3,
      vardir = "vardir",
      method = "INVALID",
      data = mysnona,
      print_result = FALSE
    )
  )

  # na auxiliary variable
  mysnona_bad <- mysnona
  mysnona_bad$x1[1] <- NA
  expect_error(
    eblup_fh(
      y ~ x1 + x2 + x3,
      vardir = "vardir",
      method = "ML",
      data = mysnona_bad,
      print_result = FALSE
    )
  )

  # vardir length
  expect_error(
    eblup_fh(
      y ~ x1 + x2 + x3,
      vardir = 1:10,
      method = "ML",
      data = mysnona,
      print_result = FALSE
    )
  )

  quiet(
    # maxiter no convergence
    fit_nc <- eblup_fh(
      y ~ x1 + x2 + x3,
      vardir = 'vardir',
      method = "ML",
      maxiter = 1,
      data = mysnona,
      print_result = FALSE
    )
  )
  expect_false(fit_nc$convergence)
})



test_that(".eblup_core direct checks and eblup_fh edge cases", {
  # 1. No sampled areas
  expect_error(fastsae:::.eblup_core(matrix(1, 2, 1), c(NA_real_, NA_real_), c(1, 1)), "No sampled areas found")

  # 2. Invalid method in .eblup_core
  expect_error(fastsae:::.eblup_core(matrix(1, 2, 1), c(1, 2), c(1, 1), method = "BAD"), "method must be 'ML' or 'REML'")

  # 3. Non-positive vardir
  expect_error(fastsae:::.eblup_core(matrix(1, 2, 1), c(1, 2), c(0, 1)), "strictly positive")

  # 4. domain = NULL and print_result = TRUE
  fit_default_dom <- eblup_fh(
    y ~ x1 + x2,
    data = mys,
    vardir = "vardir",
    domain = NULL,
    method = "ML",
    print_result = FALSE
  )

  expect_s3_class(fit_default_dom, "fastsae")
  expect_equal(fit_default_dom$df_eblup$domain, seq_len(nrow(mys)))
})


test_that("eblup_fh aborts when vardir length does not match data", {
  data("mys", package = "fastsae")
  testthat::with_mocked_bindings(
    .get_variable = function(data, variable) c(1, 2),
    .package = "fastsae",
    {
      expect_error(
        eblup_fh(y ~ x1 + x2, data = mys, vardir = "vardir"),
        "Length of 'vardir' must equal number of observations"
      )
    }
  )
})
