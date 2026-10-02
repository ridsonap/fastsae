library(testthat)
library(fastsae)

test_that("compare_sae works for two fitted models", {
  data(mys)
  fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys, print_result = FALSE)
  fit_sfh <- eblup_sfh(y ~ x1 + x2, vardir = "vardir", data = mys, W = mys_proxmat, print_result = FALSE)

  comp <- compare_sae(fit_fh, fit_sfh, names = c("FH", "Spatial FH"))

  expect_s3_class(comp, "fastsae_comparison")
  expect_equal(comp$names, c("FH", "Spatial FH"))
  expect_equal(nrow(comp$data), nrow(mys))
  expect_true("Metric" %in% names(comp$metrics))
  expect_true("Value" %in% names(comp$metrics))

  # Test print and summary
  expect_output(print(comp), "Key Concordance & Efficiency Metrics")
  expect_output(summary(comp), "Distribution of Domain Point Estimates")

  # Test autoplots
  p_scat <- autoplot(comp, type = "scatter")
  expect_s3_class(p_scat, "ggplot")

  p_diff <- autoplot(comp, type = "difference")
  expect_s3_class(p_diff, "ggplot")

  p_mse <- autoplot(comp, type = "mse")
  expect_s3_class(p_mse, "ggplot")

  p_rse <- autoplot(comp, type = "rse")
  expect_s3_class(p_rse, "ggplot")

  p_cmp <- autoplot(comp, type = "comparison")
  expect_s3_class(p_cmp, "ggplot")

  # Test plot.fastsae_comparison
  p_plot <- plot(comp, type = "scatter")
  expect_s3_class(p_plot, "ggplot")

  # Spatial map plot dispatch
  if (requireNamespace("sf", quietly = TRUE)) {
    grid_sf <- sf::st_make_grid(
      sf::st_polygon(list(matrix(c(0,0, 6,0, 6,7, 0,7, 0,0), ncol = 2, byrow = TRUE))),
      cellsize = c(1, 1),
      what = "polygons"
    )[seq_len(nrow(mys))]
    mys_sf <- sf::st_sf(area = mys$area, geometry = grid_sf)
    p_map <- plot(comp, sf_geom = mys_sf)
    expect_s3_class(p_map, "ggplot")
  }
})

test_that("compare_sae works with auto names, same names, and list inputs", {
  data(mys)
  fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys, print_result = FALSE)
  fit_sfh <- eblup_sfh(y ~ x1 + x2, vardir = "vardir", data = mys, W = mys_proxmat, print_result = FALSE)

  # 1. Auto-detected names
  comp_auto <- compare_sae(fit_fh, fit_sfh)
  expect_equal(comp_auto$names, c("FH", "Spatial FH"))

  # 2. Same models (duplicate short names)
  comp_same <- compare_sae(fit_fh, fit_fh)
  expect_equal(comp_same$names, c("FH (1)", "FH (2)"))

  # 3. Unnamed list
  comp_unnamed <- compare_sae(list(fit_fh, fit_sfh))
  expect_equal(comp_unnamed$names, c("FH", "Spatial FH"))

  # 4. Named list
  comp_list <- compare_sae(list("FH_Model" = fit_fh, "SFH_Model" = fit_sfh))
  expect_equal(comp_list$names, c("FH_Model", "SFH_Model"))
})

test_that("compare_sae works with custom data frames, time, and subarea indices", {
  # 1. Custom data frames
  df1 <- data.frame(domain = 1:5, estimate = c(1, 2, 3, 4, 5), mse = rep(0.5, 5))
  df2 <- data.frame(domain = 1:5, estimate = c(1.1, 2.1, 2.9, 4.2, 4.8), mse = rep(0.4, 5))
  comp_df <- compare_sae(df1, df2)
  expect_s3_class(comp_df, "fastsae_comparison")
  expect_equal(nrow(comp_df$data), 5)

  # 2. Data frames without MSE
  df1_nomse <- data.frame(domain = 1:3, estimate = 1:3)
  df2_nomse <- data.frame(domain = 1:3, estimate = 2:4)
  comp_nomse <- compare_sae(df1_nomse, df2_nomse)
  expect_true(all(is.na(comp_nomse$data$mse_ratio)))

  # 3. Data frames with time column
  dft1 <- data.frame(domain = c(1, 1, 2, 2), time = c(2020, 2021, 2020, 2021), estimate = 1:4)
  dft2 <- data.frame(domain = c(1, 1, 2, 2), time = c(2020, 2021, 2020, 2021), estimate = 1.5:4.5)
  comp_time <- compare_sae(dft1, dft2)
  expect_equal(nrow(comp_time$data), 4)

  # 4. Data frames with subarea column
  dfs1 <- data.frame(domain = c(1, 1), subarea = c("A", "B"), estimate = 1:2)
  dfs2 <- data.frame(domain = c(1, 1), subarea = c("A", "B"), estimate = 2:3)
  comp_sub <- compare_sae(dfs1, dfs2)
  expect_equal(nrow(comp_sub$data), 2)
})

test_that("compare_sae error handling", {
  data(mys)
  fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys, print_result = FALSE)

  # Missing second model
  expect_error(compare_sae(fit_fh), "Please provide two models")

  # No common domains
  df_a <- data.frame(domain = 1:3, estimate = 1:3)
  df_b <- data.frame(domain = 4:6, estimate = 1:3)
  expect_error(compare_sae(df_a, df_b), "No overlapping domains found")

  # Unsupported object class
  expect_error(compare_sae(123, 456), "Unsupported model object class")

  # Data frame missing estimate column
  expect_error(compare_sae(data.frame(domain = 1:3), df_b), "must contain an estimate column")

  # Model missing estimation data frame
  dummy_empty <- structure(list(), class = "fastsae")
  expect_error(compare_sae(dummy_empty, fit_fh), "does not contain estimation data frame")

  # Invalid autoplot type
  comp <- compare_sae(df_a, df_a)
  expect_error(autoplot(comp, type = "invalid"), "should be one of")
})

test_that("compare_sae preserves time column when present in model estimation table", {
  df1 <- structure(
    list(df_eblup = data.frame(domain = 1:5, eblup = 1:5, mse = 0.1, time = 1:5)),
    class = "fastsae"
  )
  df2 <- structure(
    list(df_eblup = data.frame(domain = 1:5, eblup = 1.1:5.1, mse = 0.2, time = 1:5)),
    class = "fastsae"
  )
  comp <- compare_sae(df1, df2)
  expect_true("time" %in% names(comp$data))
  expect_equal(comp$data$time, 1:5)
})
