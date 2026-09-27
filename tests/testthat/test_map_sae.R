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
  bm <- benchmark(fit_fh, target = 6.5, method = "ratio")

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
