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
})

test_that("benchmark works as alias for benchmark_sae", {
  y_hat <- c(0.2, 0.4, 0.6)
  bm1 <- benchmark_sae(y_hat, target = 0.5, method = "ratio")
  bm2 <- benchmark(y_hat, target = 0.5, method = "ratio")
  expect_equal(bm1$benchmarked, bm2$benchmarked)
  expect_s3_class(bm2, "fastsae_benchmark")
})
