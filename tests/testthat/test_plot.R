# Tests for autoplot.fastsae methods
library(testthat)

library(fastsae)

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

  # Test comparison plot using S3 generic
  p1 <- autoplot(fits_list, type = "comparison")
  expect_s3_class(p1, "ggplot")

  # Test mse plot
  p2 <- autoplot(fits_list, type = "mse")
  expect_s3_class(p2, "ggplot")

  # Test scatter plot (exactly 2 models)
  p3 <- autoplot(fits_list, type = "scatter")
  expect_s3_class(p3, "ggplot")
})

test_that("autoplot.list scatter requires exactly 2 models", {
  skip_if_not_installed("sae")

  fit_fh <- eblup_fh(y ~ x1 + x2 + x3, data = mys, vardir = "vardir", print_result = FALSE)

  # Single model should error
  expect_error(
    autoplot(list("FH" = fit_fh), type = "scatter"),
    "exactly two models"
  )

  # Three models should error
  fits_three <- list("FH1" = fit_fh, "FH2" = fit_fh, "FH3" = fit_fh)
  expect_error(
    autoplot(fits_three, type = "scatter"),
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

test_that("autoplot covers ribbon, rse, single scatter warning, map, and error paths", {
  fit_fh <- eblup_fh(y ~ x1 + x2 + x3, data = mys, vardir = "vardir", print_result = FALSE)

  # 1. Ribbon plot
  p_ribbon <- autoplot(fit_fh, type = "ribbon")
  expect_s3_class(p_ribbon, "ggplot")

  # 2. RSE plot with direct RSE
  p_rse <- autoplot(fit_fh, type = "rse")
  expect_s3_class(p_rse, "ggplot")

  # 3. RSE plot without vardir/direct RSE
  fit_no_vardir <- fit_fh
  fit_no_vardir$df_eblup$vardir <- NULL
  fit_no_vardir$df_eblup$direct_rse <- NULL
  p_rse2 <- autoplot(fit_no_vardir, type = "rse")
  expect_s3_class(p_rse2, "ggplot")

  # 4. Scatter plot for single model (warns and defaults to estimates)
  expect_warning(p_scat <- autoplot(fit_fh, type = "scatter"), "requires at least two models")
  expect_s3_class(p_scat, "ggplot")

  # 5. Estimates plot without y errors
  fit_no_y <- fit_fh
  fit_no_y$df_eblup$y <- NULL
  expect_error(autoplot(fit_no_y, type = "estimates"), "requires direct estimates 'y'")

  # 6. MSE plot without mse errors
  fit_no_mse <- fit_fh
  fit_no_mse$df_eblup$mse <- NA_real_
  expect_error(autoplot(fit_no_mse, type = "mse"), "MSE estimates are not available")

  # 7. Map type via autoplot
  if (requireNamespace("sf", quietly = TRUE)) {
    grid_sf <- sf::st_make_grid(
      sf::st_polygon(list(matrix(c(0,0, 6,0, 6,7, 0,7, 0,0), ncol = 2, byrow = TRUE))),
      cellsize = c(1, 1),
      what = "polygons"
    )[seq_len(nrow(mys))]
    mys_sf <- sf::st_sf(area = mys$area, geometry = grid_sf)

    p_map <- autoplot(fit_fh, type = "map", sf_geom = mys_sf)
    expect_s3_class(p_map, "ggplot")

    p_map_list <- autoplot(list("FH1" = fit_fh, "FH2" = fit_fh), type = "map", sf_geom = mys_sf)
    expect_s3_class(p_map_list, "ggplot")
  }
})

test_that("autoplot.list covers edge cases and multi-model branches", {
  fit_fh <- eblup_fh(y ~ x1 + x2 + x3, data = mys, vardir = "vardir", print_result = FALSE)

  # 1. Non-fastsae list calls NextMethod
  expect_error(autoplot(list(1, 2)))

  # 2. Comparison plot with < 2 models errors
  expect_error(autoplot(list(fit_fh), type = "comparison"), "requires at least two models")

  # 3. Multi-model MSE with 1 model succeeds
  p_mse1 <- autoplot(list("FH" = fit_fh), type = "mse")
  expect_s3_class(p_mse1, "ggplot")

  # 4. Multi-model MSE with empty list errors
  empty_list <- list()
  class(empty_list) <- "list"
  # autoplot.list checks vapply inherits fastsae on empty list -> TRUE -> checks length < 1
  fake_empty <- structure(list(), class = "list")
  # A list of fastsae with length 0
  # Note: all(vapply(list(), inherits, logical(1), "fastsae")) is TRUE in R!
  expect_error(autoplot(structure(list(), class = "list"), type = "mse"), "requires at least one model")

  # 5. Multi-model MSE with different domains (faceted free_x)
  fit_diff_dom <- fit_fh
  fit_diff_dom$df_eblup$domain <- paste0("Area_", seq_len(nrow(fit_diff_dom$df_eblup)))
  p_diff_mse <- autoplot(list("M1" = fit_fh, "M2" = fit_diff_dom), type = "mse")
  expect_s3_class(p_diff_mse, "ggplot")

  # 6. Multi-model scatter with no common domains errors
  fit_no_overlap <- fit_fh
  fit_no_overlap$df_eblup$domain <- paste0("Z_", seq_len(nrow(fit_no_overlap$df_eblup)))
  expect_error(autoplot(list("M1" = fit_fh, "M2" = fit_no_overlap), type = "scatter"), "No common domains found")

  # 7. Multi-model comparison without ci columns uses mse to construct bands
  fit_no_ci <- fit_fh
  fit_no_ci$df_eblup$ci_lower <- NULL
  fit_no_ci$df_eblup$ci_upper <- NULL
  p_cmp_no_ci <- autoplot(list("M1" = fit_no_ci, "M2" = fit_no_ci), type = "comparison")
  expect_s3_class(p_cmp_no_ci, "ggplot")

  # 8. RSE plot when rse is missing and direct_rse is present without vardir
  fit_calc_rse <- fit_fh
  fit_calc_rse$df_eblup$rse <- NULL
  fit_calc_rse$df_eblup$vardir <- NULL
  fit_calc_rse$df_eblup$direct_rse <- 15
  p_calc_rse <- autoplot(fit_calc_rse, type = "rse")
  expect_s3_class(p_calc_rse, "ggplot")
})

