library(testthat)
library(fastsae)

test_that("benchmark.default works for ratio, difference, optimal, and logit methods", {
  set.seed(42)
  n <- 20
  y_hat <- runif(n, 0.1, 0.4)
  w <- runif(n, 50, 500)
  target_mean <- 0.28

  # 1. Ratio benchmarking
  bm_ratio <- benchmark_sae(y_hat, target = target_mean, weight = w, method = "ratio")
  expect_s3_class(bm_ratio, "fastsae_benchmark")
  expect_equal(sum((w / sum(w)) * bm_ratio$benchmarked), target_mean, tolerance = 1e-8)
  expect_true(all(bm_ratio$benchmarked > 0))

  # 2. Difference benchmarking
  bm_diff <- benchmark_sae(y_hat, target = target_mean, weight = w, method = "difference")
  expect_s3_class(bm_diff, "fastsae_benchmark")
  expect_equal(sum((w / sum(w)) * bm_diff$benchmarked), target_mean, tolerance = 1e-8)

  # 3. Optimal benchmarking
  bm_opt <- benchmark_sae(y_hat, target = target_mean, weight = w, method = "optimal")
  expect_s3_class(bm_opt, "fastsae_benchmark")
  expect_equal(sum((w / sum(w)) * bm_opt$benchmarked), target_mean, tolerance = 1e-8)

  # 4. Logit benchmarking (bounded in (0, 1))
  bm_logit <- benchmark_sae(y_hat, target = target_mean, weight = w, method = "logit")
  expect_s3_class(bm_logit, "fastsae_benchmark")
  expect_equal(sum((w / sum(w)) * bm_logit$benchmarked), target_mean, tolerance = 1e-6)
  expect_true(all(bm_logit$benchmarked > 0 & bm_logit$benchmarked < 1))

  # 5. Type = 'total'
  target_total <- 2500
  bm_total <- benchmark_sae(y_hat, target = target_total, weight = w, method = "ratio", type = "total")
  expect_equal(sum(w * bm_total$benchmarked), target_total, tolerance = 1e-8)
})

test_that("benchmark works on fastsae eblup_fh models", {
  data(mys, package = "fastsae")
  fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)

  target_val <- 7.5
  bm <- benchmark_sae(fit_fh, target = target_val, method = "ratio")
  expect_s3_class(bm, "fastsae_benchmark")
  expect_equal(mean(bm$benchmarked), target_val, tolerance = 1e-8)

  # Optimal benchmarking with MSE
  bm_opt <- benchmark_sae(fit_fh, target = target_val, method = "optimal")
  expect_s3_class(bm_opt, "fastsae_benchmark")
  expect_equal(mean(bm_opt$benchmarked), target_val, tolerance = 1e-8)

  # Custom weight
  bm_wt <- benchmark_sae(fit_fh, target = target_val, weight = mys$n, method = "difference")
  expect_equal(sum((mys$n / sum(mys$n)) * bm_wt$benchmarked), target_val, tolerance = 1e-8)
})

test_that("benchmark works with grouping (hierarchical calibration)", {
  set.seed(123)
  n <- 30
  y_hat <- runif(n, 5, 15)
  groups <- rep(c("Region_A", "Region_B", "Region_C"), each = 10)
  targets <- c(Region_A = 8.0, Region_B = 11.5, Region_C = 14.0)

  bm_grp <- benchmark_sae(y_hat, target = targets, group = groups, method = "ratio")
  expect_s3_class(bm_grp, "fastsae_benchmark")

  # Check consistency for each group
  for (g in names(targets)) {
    sub_bm <- bm_grp[bm_grp$group == g, ]
    expect_equal(mean(sub_bm$benchmarked), unname(targets[g]), tolerance = 1e-8)
  }
})

test_that("benchmark handles input errors gracefully", {
  y_hat <- c(0.2, 0.4, 0.6)
  expect_error(benchmark_sae(y_hat, target = 0.5, weight = c(1, -1, 1)), "strictly positive")
  expect_error(benchmark_sae(y_hat, target = 0.5, weight = c(1, 1)), "does not match")
  expect_error(benchmark_sae(c(-0.1, 0.4, 0.6), target = 0.5, method = "logit"), "strictly within")
  expect_error(benchmark_sae(y_hat, target = 1.5, method = "logit"), "strictly in")
})

test_that("benchmark print, summary, and autoplot methods work", {
  data(mys, package = "fastsae")
  fit_fh <- eblup_fh(y ~ x1, vardir = "vardir", data = mys)
  bm <- benchmark_sae(fit_fh, target = 7.0, method = "ratio")

  expect_output(print(bm), "fastsae Small Area Benchmark Calibration")
  expect_output(summary(bm), "Adjustment Statistics")

  p <- autoplot(bm)
  expect_s3_class(p, "ggplot")

  # plot.fastsae_benchmark dispatch
  p_plot <- plot(bm)
  expect_s3_class(p_plot, "ggplot")
})

test_that("benchmark works as alias for benchmark_sae", {
  y_hat <- c(0.2, 0.4, 0.6)
  bm1 <- benchmark_sae(y_hat, target = 0.5, method = "ratio")
  bm2 <- benchmark(y_hat, target = 0.5, method = "ratio")
  expect_equal(bm1$benchmarked, bm2$benchmarked)
  expect_s3_class(bm2, "fastsae_benchmark")
})

test_that("benchmark covers all validation, targets, and methods branches", {
  data(mys, package = "fastsae")
  fit_fh <- eblup_fh(y ~ x1, vardir = "vardir", data = mys, print_result = FALSE)

  # 1. Target as data frame with group and target columns
  groups_df <- rep(c("G1", "G2"), each = 21)
  target_df <- data.frame(group = c("G1", "G2"), target = c(6.5, 7.5))
  bm_df_target <- benchmark_sae(fit_fh, target = target_df, group = groups_df)
  expect_s3_class(bm_df_target, "fastsae_benchmark")

  # Target data frame missing columns
  expect_error(benchmark_sae(fit_fh, target = data.frame(a = 1), group = groups_df), "must contain columns")
  # Target data frame missing a group
  expect_error(benchmark_sae(fit_fh, target = data.frame(group = "G1", target = 6.5), group = groups_df), "not found in target data frame")

  # 2. Named vector target missing a group
  expect_error(benchmark_sae(fit_fh, target = c(G1 = 6.5, G3 = 7.5), group = groups_df), "does not contain group")
  # Length mismatch
  expect_error(benchmark_sae(fit_fh, target = c(6.5, 7.5, 8.5), group = groups_df), "does not match number of groups")
  # Invalid target class
  expect_error(benchmark_sae(fit_fh, target = "invalid", group = groups_df), "must be a numeric value")

  # 3. Neither target nor national target provided
  expect_error(benchmark_sae(fit_fh), "Please provide either")

  # 4. Weight as column name in object data and df_eblup
  bm_col_wt <- benchmark_sae(fit_fh, target = 7.0, weight = "n")
  expect_s3_class(bm_col_wt, "fastsae_benchmark")
  expect_error(benchmark_sae(fit_fh, target = 7.0, weight = "nonexistent_col"), "not found in object data")
  expect_error(benchmark_sae(fit_fh, target = 7.0, weight = 1:5), "does not match number of domains")

  # 5. Group as column name and error handling
  expect_error(benchmark_sae(fit_fh, target = 7.0, group = "nonexistent_grp"), "not found in object data")
  expect_error(benchmark_sae(fit_fh, target = 7.0, group = 1:5), "does not match number of domains")

  # 6. Invalid fastsae object
  expect_error(benchmark_sae(structure(list(), class = "fastsae"), target = 7.0), "does not contain fitted area estimates")
  expect_error(benchmark_sae(structure(list(df_eblup = data.frame(a = 1)), class = "fastsae"), target = 7.0), "Could not find estimation column")

  # 7. Non-numeric object in benchmark_sae.default
  expect_error(benchmark_sae("not_numeric", target = 7.0), "must be a numeric vector")
  expect_error(benchmark_sae(1:5, target = 7.0, weight = 1:3), "Length of `weight`")
  expect_error(benchmark_sae(1:5, target = 7.0, group = 1:3), "Length of `group`")

  # 8. Hierarchical Stage 1 with difference, optimal, logit methods
  y_test <- runif(20, 0.1, 0.9)
  grp_test <- rep(c("A", "B"), each = 10)
  bm_h_diff <- benchmark_sae(y_test, national_target = 0.5, group = grp_test, method = "difference")
  expect_s3_class(bm_h_diff, "fastsae_benchmark")

  bm_h_opt <- benchmark_sae(y_test, national_target = 0.5, group = grp_test, method = "optimal")
  expect_s3_class(bm_h_opt, "fastsae_benchmark")

  bm_h_logit <- benchmark_sae(y_test, national_target = 0.5, group = grp_test, method = "logit")
  expect_s3_class(bm_h_logit, "fastsae_benchmark")

  # 9. Optimal method fallback to difference when MSE contains NA or non-positive
  expect_warning(benchmark_sae(y_test, target = 0.5, method = "optimal"), "MSE not available or non-positive")

  # 10. Aggregation virtually zero for ratio
  expect_error(benchmark_sae(c(0, 0, 0), target = 5.0, method = "ratio"), "virtually zero")

  # 11. Print all_groups = TRUE and summary.fastsae_benchmark
  expect_output(print(bm_h_diff, all_groups = TRUE))
  expect_output(summary(bm_h_diff), "Stage 1: Group Harmonization Summary")

  # 12. Benchmark with > 10 groups to test truncation in print and autoplot legend suppression
  y_12 <- runif(24, 0.2, 0.8)
  grp_12 <- rep(paste0("G", 1:12), each = 2)
  bm_12 <- benchmark_sae(y_12, national_target = 0.5, group = grp_12, method = "ratio")
  expect_output(print(bm_12, all_groups = FALSE), "and 6 more groups")
  expect_output(summary(bm_12), "and 2 more groups")

  p_12 <- autoplot(bm_12)
  expect_s3_class(p_12, "ggplot")

  # 13. plot with sf_geom
  if (requireNamespace("sf", quietly = TRUE)) {
    grid_sf <- sf::st_make_grid(
      sf::st_polygon(list(matrix(c(0,0, 6,0, 6,7, 0,7, 0,0), ncol = 2, byrow = TRUE))),
      cellsize = c(1, 1),
      what = "polygons"
    )[seq_len(nrow(mys))]
    mys_sf <- sf::st_sf(area = mys$area, geometry = grid_sf)
    p_sf <- plot(bm_col_wt, sf_geom = mys_sf)
    expect_s3_class(p_sf, "ggplot")
  }
})

