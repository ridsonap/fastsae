# Tests for autoplot.fastsae methods
library(testthat)

# Load fastsae functions
load_all_path <- system.file(package = "fastsae")
if (load_all_path == "") {
  # Not installed as package, try devtools
  devtools::load_all()
} else {
  library(fastsae)
}

skip_if_not_installed("sae")

test_that("autoplot works for single model", {
  skip_if_not_installed("sae")

  # Test Fay-Herriot model
  fit_fh <- eblup_fh(y ~ x1 + x2 + x3, data = mys, vardir = "vardir", print_result = FALSE)

  # Verify domain column exists
  expect_true("domain" %in% names(fit_fh$df_eblup))
  expect_equal(length(fit_fh$df_eblup$domain), nrow(mys))

  # Test estimates plot
  p1 <- autoplot(fit_fh, type = "estimates")
  expect_s3_class(p1, "ggplot")
  # Verify the plot has layers (not empty)
  expect_gte(length(p1$layers), 1)

  # Test mse plot
  p2 <- autoplot(fit_fh, type = "mse")
  expect_s3_class(p2, "ggplot")

  # Test comparison plot (single model)
  p3 <- autoplot(fit_fh, type = "comparison")
  expect_s3_class(p3, "ggplot")
})

test_that("autoplot.list works for multiple models", {
  skip_if_not_installed("sae")

  fit_fh <- eblup_fh(y ~ x1 + x2 + x3, data = mys, vardir = "vardir", print_result = FALSE)
  fit_sfh <- eblup_sfh(y ~ x1 + x2 + x3,
    data = mys, vardir = "vardir",
    W = mys_proxmat, print_result = FALSE, mse_method = "analytical"
  )

  # Create list of models
  fits_list <- list("FH" = fit_fh, "SFH" = fit_sfh)

  # Test comparison plot using the method directly
  p1 <- autoplot.list(fits_list, type = "comparison")
  expect_s3_class(p1, "ggplot")

  # Test mse plot
  p2 <- autoplot.list(fits_list, type = "mse")
  expect_s3_class(p2, "ggplot")

  # Test scatter plot (exactly 2 models)
  p3 <- autoplot.list(fits_list, type = "scatter")
  expect_s3_class(p3, "ggplot")
})

test_that("autoplot.list scatter requires exactly 2 models", {
  skip_if_not_installed("sae")

  fit_fh <- eblup_fh(y ~ x1 + x2 + x3, data = mys, vardir = "vardir", print_result = FALSE)

  # Single model should error
  expect_error(
    autoplot.list(list("FH" = fit_fh), type = "scatter"),
    "exactly two models"
  )

  # Three models should error
  fits_three <- list("FH1" = fit_fh, "FH2" = fit_fh, "FH3" = fit_fh)
  expect_error(
    autoplot.list(fits_three, type = "scatter"),
    "exactly two models"
  )
})

test_that("df_eblup has required columns", {
  skip_if_not_installed("sae")

  fit_fh <- eblup_fh(y ~ x1 + x2 + x3, data = mys, vardir = "vardir", print_result = FALSE)

  required_cols <- c("y", "eblup", "vardir", "mse", "rse", "domain")
  for (col in required_cols) {
    expect_true(col %in% names(fit_fh$df_eblup),
      info = paste("Column", col, "should be in df_eblup")
    )
  }
})

test_that("autoplot works with NA in y (unsampled domains)", {
  skip_if_not_installed("sae")

  # Create test data with NA values
  test_data <- mys
  test_data$y[1:5] <- NA

  fit_fh <- eblup_fh(y ~ x1 + x2 + x3, data = test_data, vardir = "vardir", print_result = FALSE)

  # Should still produce plot
  p1 <- autoplot(fit_fh, type = "estimates")
  expect_s3_class(p1, "ggplot")

  p2 <- autoplot(fit_fh, type = "mse")
  expect_s3_class(p2, "ggplot")
})

test_that("autoplot type argument validation works", {
  skip_if_not_installed("sae")

  fit_fh <- eblup_fh(y ~ x1 + x2 + x3, data = mys, vardir = "vardir", print_result = FALSE)

  # Invalid type should error
  expect_error(
    autoplot(fit_fh, type = "invalid_type"),
    NULL
  ) # match.arg throws simple error
})
