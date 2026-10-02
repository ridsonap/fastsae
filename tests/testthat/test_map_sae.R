library(testthat)
library(fastsae)

test_that("map_sae works for single model eblup_fh with auto and explicit keys", {
  skip_if_not_installed("sf")

  data(mys)
  fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)

  # Synthetic sf grid
  grid_sf <- sf::st_make_grid(
    sf::st_polygon(list(matrix(c(0,0, 6,0, 6,7, 0,7, 0,0), ncol = 2, byrow = TRUE))),
    cellsize = c(1, 1),
    what = "polygons"
  )[seq_len(nrow(mys))]
  mys_sf <- sf::st_sf(area = mys$area, geometry = grid_sf)

  # 1. Default estimate map with auto-detected key
  p_est <- map_sae(fit_fh, sf_geom = mys_sf)
  expect_s3_class(p_est, "ggplot")

  # 2. Explicit key
  p_key <- map_sae(fit_fh, sf_geom = mys_sf, key = "area")
  expect_s3_class(p_key, "ggplot")

  # 3. Named key
  p_named <- map_sae(fit_fh, sf_geom = mys_sf, key = c("domain" = "area"))
  expect_s3_class(p_named, "ggplot")

  # 4. RSE map
  p_rse <- map_sae(fit_fh, sf_geom = mys_sf, type = "rse")
  expect_s3_class(p_rse, "ggplot")

  # 5. Reliability flag map (BPS traffic light)
  p_rel <- map_sae(fit_fh, sf_geom = mys_sf, type = "reliability")
  expect_s3_class(p_rel, "ggplot")

  # 6. Direct vs SAE comparison map
  p_cmp <- map_sae(fit_fh, sf_geom = mys_sf, type = "comparison")
  expect_s3_class(p_cmp, "ggplot")

  # 7. Difference map
  p_dif <- map_sae(fit_fh, sf_geom = mys_sf, type = "difference")
  expect_s3_class(p_dif, "ggplot")
})

test_that("map_sae works with zero-friction when model data is sf", {
  skip_if_not_installed("sf")

  data(mys)
  grid_sf <- sf::st_make_grid(
    sf::st_polygon(list(matrix(c(0,0, 6,0, 6,7, 0,7, 0,0), ncol = 2, byrow = TRUE))),
    cellsize = c(1, 1),
    what = "polygons"
  )[seq_len(nrow(mys))]
  mys_sf_data <- sf::st_sf(mys, geometry = grid_sf)

  fit_sf_data <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys_sf_data)

  # No sf_geom needed!
  p <- map_sae(fit_sf_data)
  expect_s3_class(p, "ggplot")
})

test_that("map_sae works for benchmarked objects", {
  skip_if_not_installed("sf")

  data(mys)
  fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)
  bm <- benchmark_sae(fit_fh, target = 6.5, method = "ratio")

  grid_sf <- sf::st_make_grid(
    sf::st_polygon(list(matrix(c(0,0, 6,0, 6,7, 0,7, 0,0), ncol = 2, byrow = TRUE))),
    cellsize = c(1, 1),
    what = "polygons"
  )[seq_len(nrow(mys))]
  mys_sf <- sf::st_sf(area = mys$area, geometry = grid_sf)

  p_bm_cmp <- map_sae(bm, sf_geom = mys_sf, type = "comparison")
  expect_s3_class(p_bm_cmp, "ggplot")

  p_bm_dif <- map_sae(bm, sf_geom = mys_sf, type = "difference")
  expect_s3_class(p_bm_dif, "ggplot")

  p_bm_est <- map_sae(bm, sf_geom = mys_sf, type = "estimate")
  expect_s3_class(p_bm_est, "ggplot")
})

test_that("map_sae works for two-model comparison", {
  skip_if_not_installed("sf")

  data(mys)
  fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)
  fit_sfh <- eblup_sfh(y ~ x1 + x2, vardir = "vardir", data = mys, W = mys_proxmat)

  grid_sf <- sf::st_make_grid(
    sf::st_polygon(list(matrix(c(0,0, 6,0, 6,7, 0,7, 0,0), ncol = 2, byrow = TRUE))),
    cellsize = c(1, 1),
    what = "polygons"
  )[seq_len(nrow(mys))]
  mys_sf <- sf::st_sf(area = mys$area, geometry = grid_sf)

  # map_sae with model2
  p_two <- map_sae(fit_fh, model2 = fit_sfh, sf_geom = mys_sf)
  expect_s3_class(p_two, "ggplot")

  # map_sae difference
  p_diff <- map_sae(fit_fh, model2 = fit_sfh, sf_geom = mys_sf, type = "difference")
  expect_s3_class(p_diff, "ggplot")

  # map_sae with list
  p_list <- map_sae(list("FH" = fit_fh, "SFH" = fit_sfh), sf_geom = mys_sf)
  expect_s3_class(p_list, "ggplot")
})

test_that("autoplot and plot dispatch correctly to map_sae", {
  skip_if_not_installed("sf")

  data(mys)
  fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)

  grid_sf <- sf::st_make_grid(
    sf::st_polygon(list(matrix(c(0,0, 6,0, 6,7, 0,7, 0,0), ncol = 2, byrow = TRUE))),
    cellsize = c(1, 1),
    what = "polygons"
  )[seq_len(nrow(mys))]
  mys_sf <- sf::st_sf(area = mys$area, geometry = grid_sf)

  # autoplot with type = "map"
  p_ap <- autoplot(fit_fh, type = "map", sf_geom = mys_sf)
  expect_s3_class(p_ap, "ggplot")

  # autoplot with type = "rse"
  p_rse <- autoplot(fit_fh, type = "rse")
  expect_s3_class(p_rse, "ggplot")

  # plot.fastsae with sf_geom
  p_plot <- plot(fit_fh, sf_geom = mys_sf)
  expect_s3_class(p_plot, "ggplot")
})

test_that("map_sae handles errors cleanly", {
  skip_if_not_installed("sf")

  data(mys)
  fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)

  # Missing sf_geom
  expect_error(map_sae(fit_fh), "must be provided")

  # Non-existent key
  grid_sf <- sf::st_make_grid(
    sf::st_polygon(list(matrix(c(0,0, 6,0, 6,7, 0,7, 0,0), ncol = 2, byrow = TRUE))),
    cellsize = c(1, 1),
    what = "polygons"
  )[seq_len(nrow(mys))]
  mys_sf <- sf::st_sf(area = mys$area, geometry = grid_sf)

  expect_error(map_sae(fit_fh, sf_geom = mys_sf, key = "non_existent_column"), "was not found")
})

test_that("map_sae.default works with data frame and checks errors", {
  skip_if_not_installed("sf")
  data(mys)
  grid_sf <- sf::st_make_grid(
    sf::st_polygon(list(matrix(c(0,0, 6,0, 6,7, 0,7, 0,0), ncol = 2, byrow = TRUE))),
    cellsize = c(1, 1),
    what = "polygons"
  )[seq_len(nrow(mys))]
  mys_sf <- sf::st_sf(area = mys$area, geometry = grid_sf)

  df_est <- data.frame(area = mys$area, eblup = runif(nrow(mys), 5, 10), mse = runif(nrow(mys), 0.1, 0.5))

  # 1. Plain data frame
  p_df <- map_sae(df_est, sf_geom = mys_sf)
  expect_s3_class(p_df, "ggplot")

  # 2. RSE map from data frame (rse computed from mse)
  p_df_rse <- map_sae(df_est, sf_geom = mys_sf, type = "rse")
  expect_s3_class(p_df_rse, "ggplot")

  # 3. Reliability map from data frame
  p_df_rel <- map_sae(df_est, sf_geom = mys_sf, type = "reliability")
  expect_s3_class(p_df_rel, "ggplot")

  # 4. Error if object is not a data frame
  expect_error(map_sae(12345), "must be a <fastsae> object, a benchmark object, or a data frame")

  # 5. Error if data frame has no estimate column
  expect_error(map_sae(data.frame(area = 1:5, val = 1:5), sf_geom = mys_sf), "Data frame must contain an estimate column")

  # 6. Error if no sf_geom provided to data frame
  expect_error(map_sae(df_est), "Please provide spatial polygons")
})

test_that("map_sae.list edge cases and two-model errors", {
  skip_if_not_installed("sf")
  data(mys)
  fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)

  grid_sf <- sf::st_make_grid(
    sf::st_polygon(list(matrix(c(0,0, 6,0, 6,7, 0,7, 0,0), ncol = 2, byrow = TRUE))),
    cellsize = c(1, 1),
    what = "polygons"
  )[seq_len(nrow(mys))]
  mys_sf <- sf::st_sf(area = mys$area, geometry = grid_sf)

  # 1. Empty list error
  expect_error(map_sae(list()), "Empty list provided")

  # 2. Single element list delegates to single map_sae
  p_single_list <- map_sae(list(fit_fh), sf_geom = mys_sf)
  expect_s3_class(p_single_list, "ggplot")

  # 3. Two models without sf_geom errors
  expect_error(map_sae(fit_fh, model2 = fit_fh), "Please provide spatial polygons")

  # 4. Benchmarked object without sf_geom errors
  bm <- benchmark_sae(fit_fh, target = 6.5, method = "ratio")
  expect_error(map_sae(bm), "Please provide spatial polygons")
})

test_that("map_sae smart key matching errors and NA direct branches", {
  skip_if_not_installed("sf")
  data(mys)
  fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)

  grid_sf <- sf::st_make_grid(
    sf::st_polygon(list(matrix(c(0,0, 6,0, 6,7, 0,7, 0,0), ncol = 2, byrow = TRUE))),
    cellsize = c(1, 1),
    what = "polygons"
  )[seq_len(nrow(mys))]
  mys_sf <- sf::st_sf(area = mys$area, geometry = grid_sf)

  # 1. Named key with missing model column
  expect_error(map_sae(fit_fh, sf_geom = mys_sf, key = c("wrong_col" = "area")), "Model domain column .* was not found")

  # 2. Named key with missing spatial column
  expect_error(map_sae(fit_fh, sf_geom = mys_sf, key = c("domain" = "wrong_col")), "Spatial key column .* was not found")

  # 3. Auto key match fails when columns have no standard name and < 40% overlap
  mys_sf_mismatch <- sf::st_sf(unrelated_col = paste0("UNKNOWN_", 1:nrow(mys)), geometry = grid_sf)
  expect_error(map_sae(fit_fh, sf_geom = mys_sf_mismatch), "Could not automatically match model domain IDs")

  # 4. Model with all NA direct estimates
  fit_na_direct <- fit_fh
  fit_na_direct$df_eblup$y <- NA_real_

  # Comparison warns and falls back to estimate map
  expect_warning(p_na_cmp <- map_sae(fit_na_direct, sf_geom = mys_sf, type = "comparison"), "Direct estimates are not available")
  expect_s3_class(p_na_cmp, "ggplot")

  # Difference aborts with error
  expect_error(map_sae(fit_na_direct, sf_geom = mys_sf, type = "difference"), "Direct estimates are required")
})

