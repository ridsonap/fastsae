test_that("spatial_test works with standard data and returns expected structure", {
  data(mys)
  data(mys_proxmat)

  st <- spatial_test(
    formula = y ~ x1 + x2,
    vardir = "vardir",
    data = mys,
    W = mys_proxmat
  )

  expect_s3_class(st, "fastsae_spatial_test")
  expect_true("response_test" %in% names(st))
  expect_true("residual_test" %in% names(st))
  expect_true("model_comparison" %in% names(st))
  expect_true("spatial_recommended" %in% names(st))
  expect_type(st$spatial_recommended, "logical")
  expect_type(st$reason, "character")

  # 1. Direct response Moran's I
  expect_true(!is.na(st$response_test$I))
  expect_true(!is.na(st$response_test$p_value))
  expect_equal(st$response_test$expected, -1 / (st$n - 1))

  # 2. Residual Moran's I and LM tests
  expect_true(!is.na(st$residual_test$moran$I))
  expect_true(!is.na(st$residual_test$lm_error$statistic))
  expect_true(!is.na(st$residual_test$lm_error$p_value))
  expect_true(!is.na(st$residual_test$lm_lag$statistic))

  # 3. Model comparison
  expect_true(!is.null(st$model_comparison))
  expect_true(!is.na(st$model_comparison$lrt_stat))
  expect_true(!is.na(st$model_comparison$delta_aic))
  expect_true(!is.na(st$model_comparison$rho))

  # 4. Print method
  expect_message(print(st), "Spatial")
  expect_message(print(st), "Direct Response Variable")
})

test_that("spatial_test handles data without vardir (ESDA only)", {
  data(mys)
  data(mys_proxmat)

  # Run spatial test without vardir
  st_no_var <- spatial_test(
    formula = y ~ x1 + x2,
    data = mys,
    W = mys_proxmat
  )

  expect_s3_class(st_no_var, "fastsae_spatial_test")
  expect_true(!is.na(st_no_var$response_test$I))
  expect_true(!is.na(st_no_var$residual_test$moran$I))
  expect_null(st_no_var$model_comparison)
  expect_message(print(st_no_var), "Spatial")
})

test_that("spatial_test auto-constructs W when data is sf and W is NULL", {
  skip_if_not_installed("sf")

  # Create small 9-polygon grid
  grid_geom <- sf::st_make_grid(
    sf::st_polygon(list(matrix(c(0,0, 3,0, 3,3, 0,3, 0,0), ncol = 2, byrow = TRUE))),
    cellsize = c(1, 1)
  )
  set.seed(42)
  df_sf <- sf::st_sf(
    domain = 1:9,
    y = rnorm(9, 10, 2),
    x1 = runif(9, 1, 5),
    vardir = rep(0.5, 9),
    geometry = grid_geom
  )

  # Call spatial_test with W = NULL
  st_sf <- spatial_test(
    formula = y ~ x1,
    vardir = "vardir",
    data = df_sf
  )

  expect_s3_class(st_sf, "fastsae_spatial_test")
  expect_true(!is.null(st_sf$W))
  expect_message(print(st_sf), "Spatial")
})
