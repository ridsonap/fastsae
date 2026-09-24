library(testthat)
library(fastsae)

skip_if_not_installed("INLA")

# Prepare test data
mysnona <- mys[!is.na(mys$y), ]

# ------------------------------------------------------------------
# Basic structure tests
# ------------------------------------------------------------------

test_that("ebp_area returns valid fastsae structure", {
  fit <- ebp_area(
    y ~ x1 + x2 + x3,
    data = mysnona,
    family = "gaussian",
    spatial = "none",
    print_result = FALSE
  )

  expect_s3_class(fit, "fastsae")
  expect_true("df_eblup" %in% names(fit))
  expect_true("estcoef" %in% names(fit))
  expect_true("random_effect_var" %in% names(fit))
  expect_true("goodness" %in% names(fit))
  expect_equal(fit$model, "INLA")
  expect_equal(fit$family, "gaussian")
  expect_equal(fit$spatial, "none")
})

# ------------------------------------------------------------------
# Gaussian non-spatial tests
# ------------------------------------------------------------------

test_that("ebp_area Gaussian returns reasonable estimates", {
  # Skip if sae package not available
  skip_if_not_installed("sae")

  # Fit with INLA
  fit_inla <- ebp_area(
    y ~ x1 + x2 + x3,
    data = mysnona,
    family = "gaussian",
    spatial = "none",
    print_result = FALSE
  )

  # INLA should return estimates in a reasonable range
  # (compared to the original y values)
  expect_true(all(fit_inla$df_eblup$eblup > 0))
  expect_true(all(fit_inla$df_eblup$mse >= 0))
  expect_true(all(fit_inla$df_eblup$rse >= 0, na.rm = TRUE))
})

test_that("ebp_area returns correct number of estimates", {
  n_domains <- nrow(mysnona)

  fit <- ebp_area(
    y ~ x1 + x2 + x3,
    data = mysnona,
    family = "gaussian",
    spatial = "none",
    print_result = FALSE
  )

  expect_equal(nrow(fit$df_eblup), n_domains)
  expect_equal(length(fit$df_eblup$eblup), n_domains)
  expect_equal(length(fit$df_eblup$mse), n_domains)
})

# ------------------------------------------------------------------
# Spatial model tests
# ------------------------------------------------------------------

test_that("ebp_area BYM spatial model works", {
  # Use correct W dimensions
  W_sub <- mys_proxmat[1:32, 1:32]
  fit_bym <- ebp_area(
    y ~ x1 + x2 + x3,
    data = mysnona,
    family = "gaussian",
    spatial = "bym",
    W = W_sub,
    print_result = FALSE
  )

  expect_s3_class(fit_bym, "fastsae")
  expect_equal(fit_bym$spatial, "bym")
  expect_false(is.null(fit_bym$rho))
})

test_that("ebp_area BYM2 spatial model works", {
  W_sub <- mys_proxmat[1:32, 1:32]
  fit_bym2 <- ebp_area(
    y ~ x1 + x2 + x3,
    data = mysnona,
    family = "gaussian",
    spatial = "bym2",
    W = W_sub,
    print_result = FALSE
  )

  expect_s3_class(fit_bym2, "fastsae")
  expect_equal(fit_bym2$spatial, "bym2")
})

# ------------------------------------------------------------------
# Besagproper spatial model tests
# ------------------------------------------------------------------

test_that("ebp_area Besagproper model works", {
  skip_if_not_installed("INLA")

  W_sub <- mys_proxmat[1:32, 1:32]
  fit_besag <- ebp_area(
    y ~ x1 + x2 + x3,
    data = mysnona,
    family = "gaussian",
    spatial = "besagproper",
    W = W_sub,
    print_result = FALSE
  )

  expect_s3_class(fit_besag, "fastsae")
  expect_equal(fit_besag$spatial, "besagproper")
  expect_true(all(fit_besag$df_eblup$mse >= 0))
})

test_that("ebp_area Besagproper returns reasonable estimates", {
  skip_if_not_installed("INLA")

  W_sub <- mys_proxmat[1:32, 1:32]
  fit_besag <- ebp_area(
    y ~ x1 + x2 + x3,
    data = mysnona,
    family = "gaussian",
    spatial = "besagproper",
    W = W_sub,
    print_result = FALSE
  )

  expect_true(all(fit_besag$df_eblup$eblup > 0))
  expect_true(all(fit_besag$df_eblup$rse >= 0, na.rm = TRUE))
})

# ------------------------------------------------------------------
# Random effect variance tests for all spatial models
# ------------------------------------------------------------------

test_that("ebp_area returns random effect variance for all spatial models", {
  skip_if_not_installed("INLA")

  W_sub <- mys_proxmat[1:32, 1:32]

  # Test iid (none)
  fit_none <- ebp_area(
    y ~ x1 + x2 + x3,
    data = mysnona,
    family = "gaussian",
    spatial = "none",
    print_result = FALSE
  )
  expect_true(is.numeric(fit_none$random_effect_var))

  # Test Besagproper
  fit_besag <- ebp_area(
    y ~ x1 + x2 + x3,
    data = mysnona,
    family = "gaussian",
    spatial = "besagproper",
    W = W_sub,
    print_result = FALSE
  )
  expect_true(is.numeric(fit_besag$random_effect_var))
  expect_true(fit_besag$random_effect_var >= 0)
})

test_that("ebp_area errors without W for spatial models", {
  expect_error(
    ebp_area(
      y ~ x1 + x2 + x3,
      data = mysnona,
      family = "gaussian",
      spatial = "bym",
      print_result = FALSE
    ),
    "W.*required"
  )

  expect_error(
    ebp_area(
      y ~ x1 + x2 + x3,
      data = mysnona,
      family = "gaussian",
      spatial = "bym2",
      print_result = FALSE
    ),
    "W.*required"
  )

  expect_error(
    ebp_area(
      y ~ x1 + x2 + x3,
      data = mysnona,
      family = "gaussian",
      spatial = "besagproper",
      print_result = FALSE
    ),
    "W.*required"
  )
})

test_that("ebp_area errors with wrong W dimensions", {
  # Wrong dimension W
  wrong_W <- matrix(0.25, 10, 10)

  expect_error(
    ebp_area(
      y ~ x1 + x2 + x3,
      data = mysnona,
      family = "gaussian",
      spatial = "bym",
      W = wrong_W,
      print_result = FALSE
    ),
    "dimension"
  )
})

# ------------------------------------------------------------------
# Domain handling tests
# ------------------------------------------------------------------

test_that("ebp_area handles domain argument", {
  fit <- ebp_area(
    y ~ x1 + x2 + x3,
    data = mysnona,
    family = "gaussian",
    spatial = "none",
    domain = "area",
    print_result = FALSE
  )

  expect_true("domain" %in% colnames(fit$df_eblup))
  expect_equal(fit$df_eblup$domain, mysnona$area)
})

test_that("ebp_area auto-generates domain if NULL", {
  fit <- ebp_area(
    y ~ x1 + x2 + x3,
    data = mysnona,
    family = "gaussian",
    spatial = "none",
    domain = NULL,
    print_result = FALSE
  )

  expect_true("domain" %in% colnames(fit$df_eblup))
  expect_equal(length(unique(fit$df_eblup$domain)), nrow(mysnona))
})

# ------------------------------------------------------------------
# Poisson family tests
# ------------------------------------------------------------------

test_that("ebp_area Poisson model works", {
  # Create count data
  set.seed(123)
  n <- 30
  count_data <- data.frame(
    y = rpois(n, lambda = exp(0.5 + 0.8 * rnorm(n))),
    x1 = rnorm(n),
    x2 = rnorm(n),
    area = paste0("area_", 1:n)
  )

  fit_poisson <- ebp_area(
    y ~ x1 + x2,
    data = count_data,
    family = "poisson",
    spatial = "none",
    print_result = FALSE
  )

  expect_s3_class(fit_poisson, "fastsae")
  expect_equal(fit_poisson$family, "poisson")
})

# ------------------------------------------------------------------
# Binomial family tests
# ------------------------------------------------------------------

test_that("ebp_area Binomial model works", {
  # Create binary data (simpler case)
  set.seed(456)
  n <- 30
  p <- plogis(0.5 + 0.6 * rnorm(n))
  y <- rbinom(n, size = 1, prob = p)  # Binary outcome

  binom_data <- data.frame(
    y = y,
    x1 = rnorm(n),
    x2 = rnorm(n),
    area = paste0("area_", 1:n)
  )

  fit_binom <- ebp_area(
    y ~ x1 + x2,
    data = binom_data,
    family = "binomial",
    spatial = "none",
    print_result = FALSE
  )

  expect_s3_class(fit_binom, "fastsae")
  expect_equal(fit_binom$family, "binomial")
})

# ------------------------------------------------------------------
# Model comparison metrics tests
# ------------------------------------------------------------------

test_that("ebp_area returns goodness of fit metrics", {
  fit <- ebp_area(
    y ~ x1 + x2 + x3,
    data = mysnona,
    family = "gaussian",
    spatial = "none",
    print_result = FALSE
  )

  expect_true("goodness" %in% names(fit))
  expect_true("dic" %in% names(fit$goodness))
  expect_true("waic" %in% names(fit$goodness))
})

# ------------------------------------------------------------------
# Input validation tests
# ------------------------------------------------------------------

test_that("ebp_area validates family argument", {
  expect_error(
    ebp_area(
      y ~ x1 + x2 + x3,
      data = mysnona,
      family = "invalid_family",
      spatial = "none",
      print_result = FALSE
    )
  )
})

test_that("ebp_area validates spatial argument", {
  expect_error(
    ebp_area(
      y ~ x1 + x2 + x3,
      data = mysnona,
      family = "gaussian",
      spatial = "invalid_spatial",
      print_result = FALSE
    )
  )
})

test_that("ebp_area checks for NA in auxiliary variables", {
  bad_data <- mysnona
  bad_data$x1[1] <- NA

  expect_error(
    ebp_area(
      y ~ x1 + x2 + x3,
      data = bad_data,
      family = "gaussian",
      spatial = "none",
      print_result = FALSE
    ),
    "NA"
  )
})

# ------------------------------------------------------------------
# MSE and RSE tests
# ------------------------------------------------------------------

test_that("ebp_area returns valid MSE and RSE", {
  fit <- ebp_area(
    y ~ x1 + x2 + x3,
    data = mysnona,
    family = "gaussian",
    spatial = "none",
    print_result = FALSE
  )

  # MSE should be non-negative
  expect_true(all(fit$df_eblup$mse >= 0))

  # RSE should be non-negative
  expect_true(all(fit$df_eblup$rse >= 0, na.rm = TRUE))

  # RSE should be defined for non-zero estimates
  nonzero_eblup <- fit$df_eblup$eblup != 0
  expect_true(all(is.finite(fit$df_eblup$rse[nonzero_eblup])))
})
