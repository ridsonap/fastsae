test_that("create_weights works correctly with coordinate matrix", {
  set.seed(42)
  coords <- matrix(runif(40), ncol = 2)
  rownames(coords) <- paste0("Area_", 1:20)

  # 1. k-NN method with row-standardization
  W_knn <- create_weights(coords, method = "knn", k = 4, style = "W")
  expect_s3_class(W_knn, "fastsae_weights")
  expect_equal(dim(W_knn), c(20, 20))
  expect_equal(rownames(W_knn), rownames(coords))
  expect_equal(colnames(W_knn), rownames(coords))
  expect_equal(as.numeric(diag(W_knn)), rep(0, 20))
  expect_equal(unname(rowSums(W_knn)), rep(1, 20), tolerance = 1e-7)

  # 2. k-NN method with binary adjacency
  W_bin <- create_weights(coords, method = "knn", k = 3, style = "B")
  expect_true(all(W_bin %in% c(0, 1)))
  expect_equal(unname(rowSums(W_bin)), rep(3, 20))
  expect_equal(as.numeric(diag(W_bin)), rep(0, 20))

  # 3. Distance band method
  W_dist <- create_weights(coords, method = "distance", style = "W")
  expect_s3_class(W_dist, "fastsae_weights")
  expect_true(all(rowSums(W_dist) > 0))

  # 4. IDW method
  W_idw <- create_weights(coords, method = "idw", alpha = 1, style = "W")
  expect_s3_class(W_idw, "fastsae_weights")
  expect_equal(unname(rowSums(W_idw)), rep(1, 20), tolerance = 1e-7)
  expect_equal(as.numeric(diag(W_idw)), rep(0, 20))

  # 5. S3 methods print and summary
  expect_message(print(W_knn), "Spatial Weights Matrix")
  expect_output(print(W_knn), "Matrix preview")
  expect_message(summary(W_knn), "Spatial Neighbor Distribution")
})

test_that("create_weights handles isolated islands with fallback", {
  # Create coordinates where 1 area is far away (isolated island)
  coords_clustered <- matrix(runif(18, min = 0, max = 1), ncol = 2)
  coords_island <- matrix(c(100, 100), ncol = 2)
  coords_all <- rbind(coords_clustered, coords_island)
  rownames(coords_all) <- paste0("Dom_", 1:10)

  # With distance threshold that leaves the island isolated
  W_fallback <- create_weights(
    coords_all,
    method = "distance",
    d_max = 2.0,
    island_fallback = TRUE,
    k_fallback = 1,
    style = "W"
  )

  expect_s3_class(W_fallback, "fastsae_weights")
  # Island should have at least 1 neighbor connected via fallback
  expect_gt(rowSums(W_fallback)[10], 0)
  expect_equal(attr(W_fallback, "n_islands"), 1)
  expect_equal(attr(W_fallback, "island_domains"), "Dom_10")
})

test_that("create_weights works with sf polygon data and Queen/Rook", {
  skip_if_not_installed("sf")

  # Create a simple 3x3 regular grid of polygons
  grid_geom <- sf::st_make_grid(
    sf::st_polygon(list(matrix(c(0,0, 3,0, 3,3, 0,3, 0,0), ncol = 2, byrow = TRUE))),
    cellsize = c(1, 1)
  )
  sf_grid <- sf::st_sf(id = paste0("G", 1:9), geometry = grid_geom)

  # 1. Queen contiguity
  W_queen <- create_weights(sf_grid, method = "queen", style = "W")
  expect_s3_class(W_queen, "fastsae_weights")
  expect_equal(dim(W_queen), c(9, 9))
  expect_equal(unname(rowSums(W_queen)), rep(1, 9), tolerance = 1e-7)
  # Center cell (index 5) in a 3x3 grid has 8 Queen neighbors
  expect_equal(sum(W_queen[5, ] > 0), 8)

  # 2. Rook contiguity
  W_rook <- create_weights(sf_grid, method = "rook", style = "W")
  expect_s3_class(W_rook, "fastsae_weights")
  # Center cell in a 3x3 grid has 4 Rook neighbors (only edges, no corners)
  expect_equal(sum(W_rook[5, ] > 0), 4)

  # Queen must have >= neighbors than Rook for all cells
  expect_true(all(rowSums(W_queen > 0) >= rowSums(W_rook > 0)))

  # 3. Binary adjacency for BYM2 (style = "B")
  W_bym2 <- create_weights(sf_grid, method = "queen", style = "B")
  expect_true(all(W_bym2 %in% c(0, 1)))
  # Must be symmetric for INLA
  expect_equal(W_bym2, t(W_bym2))
})

test_that("eblup_sfh automatically constructs weights when data is sf and W is NULL", {
  skip_if_not_installed("sf")

  # Create 9-polygon dataset with simulated variables (3x3 grid)
  grid_geom <- sf::st_make_grid(
    sf::st_polygon(list(matrix(c(0,0, 3,0, 3,3, 0,3, 0,0), ncol = 2, byrow = TRUE))),
    cellsize = c(1, 1)
  )
  set.seed(123)
  df_sf <- sf::st_sf(
    domain = 1:9,
    y = c(10.2, 12.5, 9.8, 14.1, 11.0, 13.2, 8.9, 15.0, 12.1),
    x = c(1.1, 2.0, 0.9, 2.5, 1.8, 2.2, 0.8, 2.9, 1.9),
    vardir = rep(0.5, 9),
    geometry = grid_geom
  )

  # Run eblup_sfh without passing W
  fit_auto_w <- eblup_sfh(
    formula = y ~ x,
    vardir = "vardir",
    data = df_sf,
    print_result = FALSE
  )

  expect_s3_class(fit_auto_w, "fastsae")
  expect_true(!is.null(fit_auto_w$W))
  expect_equal(dim(fit_auto_w$W), c(9, 9))
  # Check that geometry is preserved on df_eblup
  expect_s3_class(fit_auto_w$df_eblup, "sf")
})
