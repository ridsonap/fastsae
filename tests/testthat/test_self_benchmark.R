library(testthat)
library(fastsae)

test_that("hb_area self-benchmarking works for Gaussian Fay-Herriot with default direct target", {
  skip_if_not_installed("INLA")

  data("mys", package = "fastsae")

  # Fit self-benchmarked model (Wang, Fuller, and Qu, 2008)
  fit_sb <- hb_area(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = "vardir",
    family = "gaussian",
    self_benchmark = TRUE,
    benchmark_weight = "n",
    print_result = FALSE
  )

  expect_s3_class(fit_sb, "fastsae_hb_area")
  expect_true(fit_sb$self_benchmark)
  expect_true(!is.null(fit_sb$benchmark_summary))

  # Verify benchmarking constraint: weighted sum of hb equals direct target
  idx_samp <- which(!is.na(mys$y))
  w <- mys$n[idx_samp] / sum(mys$n[idx_samp])
  direct_target <- sum(w * mys$y[idx_samp])
  model_aggregate <- sum(w * fit_sb$df_hb$hb[idx_samp])

  expect_equal(model_aggregate, direct_target, tolerance = 1e-4)
  expect_true(abs(fit_sb$benchmark_summary$discrepancy[1]) < 1e-4)

  # Check that original y is preserved in df_hb
  expect_equal(fit_sb$df_hb$y, mys$y)

  # Verify print and summary methods work
  expect_output(print(fit_sb))
  expect_output(print(summary(fit_sb)))
})

test_that("hb_area self-benchmarking works with external fixed target", {
  skip_if_not_installed("INLA")

  data("mys", package = "fastsae")

  ext_target <- 6.5
  fit_sb_ext <- hb_area(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = "vardir",
    family = "gaussian",
    self_benchmark = TRUE,
    benchmark_weight = "n",
    benchmark_target = ext_target,
    print_result = FALSE
  )

  expect_true(fit_sb_ext$self_benchmark)
  expect_equal(fit_sb_ext$benchmark_summary$target[1], ext_target)

  idx_samp <- which(!is.na(mys$y))
  w <- mys$n[idx_samp] / sum(mys$n[idx_samp])
  model_aggregate <- sum(w * fit_sb_ext$df_hb$hb[idx_samp])

  expect_equal(model_aggregate, ext_target, tolerance = 1e-4)
})

test_that("hb_area self-benchmarking works with groups (e.g. provinces)", {
  skip_if_not_installed("INLA")

  data("mys", package = "fastsae")
  mys$prov <- rep(c("Region_A", "Region_B", "Region_C"), length.out = nrow(mys))

  fit_sb_grp <- hb_area(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = "vardir",
    family = "gaussian",
    self_benchmark = TRUE,
    benchmark_weight = "n",
    benchmark_group = "prov",
    print_result = FALSE
  )

  expect_true(fit_sb_grp$self_benchmark)
  expect_equal(nrow(fit_sb_grp$benchmark_summary), 3)

  # Verify each group matches its target
  for (i in 1:nrow(fit_sb_grp$benchmark_summary)) {
    expect_true(abs(fit_sb_grp$benchmark_summary$discrepancy[i]) < 1e-4)
  }
})

test_that("hb_area self-benchmarking works with spatial BYM2 model", {
  skip_if_not_installed("INLA")

  data("mys", package = "fastsae")
  data("mys_proxmat", package = "fastsae")

  fit_sb_spatial <- hb_area(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = "vardir",
    family = "gaussian",
    spatial = "bym2",
    W = mys_proxmat,
    self_benchmark = TRUE,
    benchmark_weight = "n",
    print_result = FALSE
  )

  expect_true(fit_sb_spatial$self_benchmark)
  expect_false(is.null(fit_sb_spatial$phi))

  idx_samp <- which(!is.na(mys$y))
  w <- mys$n[idx_samp] / sum(mys$n[idx_samp])
  direct_target <- sum(w * mys$y[idx_samp])
  model_aggregate <- sum(w * fit_sb_spatial$df_hb$hb[idx_samp])

  expect_equal(model_aggregate, direct_target, tolerance = 1e-4)
})

test_that("eblup_fh self-benchmarking works (Frequentist Wang-Fuller-Qu)", {
  data("mys", package = "fastsae")

  # 1. Global self-benchmarking
  fit_fh_sb <- eblup_fh(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = "vardir",
    self_benchmark = TRUE,
    benchmark_weight = "n",
    print_result = FALSE
  )

  expect_s3_class(fit_fh_sb, "fastsae")
  expect_true(fit_fh_sb$self_benchmark)

  idx_samp <- which(!is.na(mys$y))
  w <- mys$n[idx_samp] / sum(mys$n[idx_samp])
  direct_target <- sum(w * mys$y[idx_samp])
  model_aggregate <- sum(w * fit_fh_sb$df_eblup$eblup[idx_samp])

  # In frequentist GLS, Wang-Fuller-Qu satisfies the constraint to machine precision
  expect_equal(model_aggregate, direct_target, tolerance = 1e-10)
  expect_true(abs(fit_fh_sb$benchmark_summary$discrepancy[1]) < 1e-10)

  # 2. Grouped self-benchmarking in eblup_fh
  mys$prov <- rep(c("P1", "P2", "P3"), length.out = nrow(mys))
  fit_fh_grp <- eblup_fh(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = "vardir",
    self_benchmark = TRUE,
    benchmark_weight = "n",
    benchmark_group = "prov",
    print_result = FALSE
  )

  expect_equal(nrow(fit_fh_grp$benchmark_summary), 3)
  for (i in 1:nrow(fit_fh_grp$benchmark_summary)) {
    expect_true(abs(fit_fh_grp$benchmark_summary$discrepancy[i]) < 1e-10)
  }
})

test_that("self-benchmarking handles input validation properly", {
  data("mys", package = "fastsae")

  # Negative weights error
  bad_mys <- mys
  bad_mys$n[1] <- -10
  expect_error(
    eblup_fh(y ~ x1 + x2, data = bad_mys, vardir = "vardir", self_benchmark = TRUE, benchmark_weight = "n"),
    "must be non-negative"
  )

  # Non-matching weight length
  expect_error(
    eblup_fh(y ~ x1 + x2, data = mys, vardir = "vardir", self_benchmark = TRUE, benchmark_weight = 1:5),
    "length does not match"
  )
})

test_that("eblup_sfh self-benchmarking works for Spatial Fay-Herriot", {
  data("mys", package = "fastsae")
  data("mys_proxmat", package = "fastsae")

  # 1. Global self-benchmarking with analytical MSE
  fit_sfh_sb <- eblup_sfh(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = "vardir",
    W = mys_proxmat,
    self_benchmark = TRUE,
    benchmark_weight = "n",
    print_result = FALSE
  )

  expect_s3_class(fit_sfh_sb, "fastsae")
  expect_true(fit_sfh_sb$self_benchmark)
  expect_false(is.null(fit_sfh_sb$benchmark_summary))

  idx_samp <- which(!is.na(mys$y))
  w <- mys$n[idx_samp] / sum(mys$n[idx_samp])
  direct_target <- sum(w * mys$y[idx_samp])
  model_aggregate <- sum(w * fit_sfh_sb$df_eblup$eblup[idx_samp])

  expect_equal(model_aggregate, direct_target, tolerance = 1e-10)
  expect_true(abs(fit_sfh_sb$benchmark_summary$discrepancy[1]) < 1e-10)

  # 2. External target in eblup_sfh
  ext_target <- 6.25
  fit_sfh_ext <- eblup_sfh(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = "vardir",
    W = mys_proxmat,
    self_benchmark = TRUE,
    benchmark_weight = "n",
    benchmark_target = ext_target,
    print_result = FALSE
  )

  expect_equal(fit_sfh_ext$benchmark_summary$target[1], ext_target)
  model_agg_ext <- sum(w * fit_sfh_ext$df_eblup$eblup[idx_samp])
  expect_equal(model_agg_ext, ext_target, tolerance = 1e-10)

  # 3. Grouped self-benchmarking in eblup_sfh
  mys$reg <- rep(c("North", "South"), length.out = nrow(mys))
  fit_sfh_grp <- eblup_sfh(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = "vardir",
    W = mys_proxmat,
    self_benchmark = TRUE,
    benchmark_weight = "n",
    benchmark_group = "reg",
    print_result = FALSE
  )

  expect_equal(nrow(fit_sfh_grp$benchmark_summary), 2)
  for (i in 1:nrow(fit_sfh_grp$benchmark_summary)) {
    expect_true(abs(fit_sfh_grp$benchmark_summary$discrepancy[i]) < 1e-10)
  }

  # 4. Print and summary methods
  expect_output(print(fit_sfh_sb))
  expect_output(print(summary(fit_sfh_sb)))
})

test_that("eblup_twofold self-benchmarking works for subarea models", {
  set.seed(42)
  m <- 8
  dat <- do.call(rbind, lapply(1:m, function(d) {
    nd <- 4
    data.frame(
      area = paste0("Kab_", d),
      subarea = paste0("Kec_", d, "_", seq_len(nd)),
      x1 = rnorm(nd),
      vardir = runif(nd, 0.1, 0.4),
      pop = sample(100:500, nd),
      stringsAsFactors = FALSE
    )
  }))
  dat$y <- 2 + 1.2 * dat$x1 + rnorm(nrow(dat), 0, 0.5)

  # 1. Global self-benchmarking
  fit_tf_global <- eblup_twofold(
    y ~ x1,
    vardir = "vardir",
    domain = "area",
    subarea = "subarea",
    data = dat,
    self_benchmark = TRUE,
    benchmark_weight = "pop",
    print_result = FALSE
  )

  expect_s3_class(fit_tf_global, "fastsae")
  expect_true(fit_tf_global$self_benchmark)
  expect_equal(nrow(fit_tf_global$benchmark_summary), 1)
  expect_true(abs(fit_tf_global$benchmark_summary$discrepancy[1]) < 1e-10)

  # 2. Area-level (Subarea-to-Area) self-benchmarking
  fit_tf_area <- eblup_twofold(
    y ~ x1,
    vardir = "vardir",
    domain = "area",
    subarea = "subarea",
    data = dat,
    self_benchmark = TRUE,
    benchmark_weight = "pop",
    benchmark_group = "area",
    print_result = FALSE
  )

  expect_equal(nrow(fit_tf_area$benchmark_summary), m)
  for (i in seq_len(m)) {
    expect_true(abs(fit_tf_area$benchmark_summary$discrepancy[i]) < 1e-10)
  }
})

test_that("eblup_stfh self-benchmarking works for Spatio-Temporal models", {
  skip_if_not_installed("sae")
  library(dplyr)
  data("mys_panel", package = "fastsae")
  data("mys_proxmat", package = "fastsae")

  mys_panel_nona <- mys_panel |> filter(!is.na(y))
  mys_proxmat_nona <- mys_proxmat[-c(21, 25), -c(21, 25)]

  # 1. Global self-benchmarking
  fit_st_global <- eblup_stfh(
    y ~ x1 + x2 + x3,
    data = mys_panel_nona,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat_nona,
    model = "ST",
    self_benchmark = TRUE,
    print_result = FALSE
  )

  expect_s3_class(fit_st_global, "fastsae")
  expect_true(fit_st_global$self_benchmark)
  expect_true(abs(fit_st_global$benchmark_summary$discrepancy[1]) < 1e-10)

  # 2. Per-year self-benchmarking
  fit_st_year <- eblup_stfh(
    y ~ x1 + x2 + x3,
    data = mys_panel_nona,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat_nona,
    model = "ST",
    self_benchmark = TRUE,
    benchmark_group = "year",
    print_result = FALSE
  )

  n_years <- length(unique(mys_panel_nona$year))
  expect_equal(nrow(fit_st_year$benchmark_summary), n_years)
  for (i in seq_len(n_years)) {
    expect_true(abs(fit_st_year$benchmark_summary$discrepancy[i]) < 1e-10)
  }
})


