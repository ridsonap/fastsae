test_that("two-stage hierarchical benchmarking works for D = 500 and 34 provinces (national target only)", {
  set.seed(42)
  d <- 500
  n_prov <- 34
  prov_names <- paste0("Prov_", sprintf("%02d", 1:n_prov))
  domain_prov <- sample(prov_names, d, replace = TRUE)

  x1 <- rnorm(d)
  x2 <- runif(d)
  vardir <- runif(d, 0.1, 0.8)
  y <- 5 + 0.8 * x1 - 0.5 * x2 + rnorm(d, 0, sqrt(vardir))
  pop <- round(runif(d, 50000, 500000))

  df_indo <- data.frame(
    domain = paste0("Kab_", sprintf("%03d", 1:d)),
    provinsi = domain_prov,
    y = y,
    vardir = vardir,
    x1 = x1,
    x2 = x2,
    pop = pop,
    stringsAsFactors = FALSE
  )

  fit <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = df_indo, print_result = FALSE)

  nat_target <- 6.5

  # Two-stage hierarchical calibration with ratio method
  bm <- benchmark_sae(
    object = fit,
    national_target = nat_target,
    group = "provinsi",
    weight = "pop",
    method = "ratio",
    type = "mean"
  )

  expect_s3_class(bm, "fastsae_benchmark")
  expect_true(isTRUE(attr(bm, "hierarchical")))
  expect_equal(nrow(bm), d)
  expect_true("national_target" %in% names(bm))
  expect_equal(unique(bm$national_target), nat_target)

  # 1. Verify Level 0 (National Aggregate)
  w_norm_nat <- df_indo$pop / sum(df_indo$pop)
  nat_calib_sum <- sum(w_norm_nat * bm$benchmarked)
  expect_equal(nat_calib_sum, nat_target, tolerance = 1e-7)

  # 2. Verify Level 1 (All 34 Provincial Groups)
  s1_sum <- attr(bm, "stage1_summary")
  expect_equal(nrow(s1_sum), n_prov)

  for (pr in prov_names) {
    sub_bm <- bm[bm$group == pr, ]
    sub_w <- df_indo$pop[df_indo$provinsi == pr]
    w_norm_prov <- sub_w / sum(sub_w)
    prov_calib_sum <- sum(w_norm_prov * sub_bm$benchmarked)

    expected_target <- s1_sum$calibrated_target[s1_sum$group == pr]
    expect_equal(prov_calib_sum, expected_target, tolerance = 1e-7)
  }

  # Test print and summary with 34 groups
  expect_output(print(bm), "Two-Stage Hierarchical Benchmark Calibration")
  expect_output(print(bm), "Level 0 \\(National Target\\)")
  expect_output(print(bm), "34 groups")
  expect_output(summary(bm), "Stage 1: provinsi Harmonization Summary")
})

test_that("two-stage hierarchical benchmarking works with initial provincial targets", {
  set.seed(123)
  d <- 300
  n_prov <- 10
  prov_names <- paste0("Region_", 1:n_prov)
  domain_prov <- sample(prov_names, d, replace = TRUE)

  w <- runif(d, 100, 1000)
  y <- runif(d, 4, 12)

  # Initial user targets for 10 provinces (which don't sum to national target)
  init_targets <- stats::setNames(runif(n_prov, 6.0, 9.0), prov_names)
  nat_target <- 7.25

  # Two-stage calibration with difference method
  bm_diff <- benchmark_sae(
    y,
    target = init_targets,
    national_target = nat_target,
    group = domain_prov,
    weight = w,
    method = "difference",
    type = "mean"
  )

  expect_s3_class(bm_diff, "fastsae_benchmark")
  expect_true(isTRUE(attr(bm_diff, "hierarchical")))

  # Check national consistency
  nat_agg <- sum((w / sum(w)) * bm_diff$benchmarked)
  expect_equal(nat_agg, nat_target, tolerance = 1e-7)

  # Check each province consistency
  s1 <- attr(bm_diff, "stage1_summary")
  for (pr in prov_names) {
    idx <- which(domain_prov == pr)
    w_p <- w[idx] / sum(w[idx])
    p_agg <- sum(w_p * bm_diff$benchmarked[idx])
    expect_equal(p_agg, s1$calibrated_target[s1$group == pr], tolerance = 1e-7)
  }
})
