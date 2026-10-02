library(testthat)
library(fastsae)

test_that("add_reliability_flags adds correct tiers and handles all edge cases", {
  data(mys)
  fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys, print_result = FALSE)

  df_flagged <- add_reliability_flags(fit_fh)
  expect_true("reliability" %in% names(df_flagged))
  expect_s3_class(df_flagged$reliability, "factor")
  expect_equal(
    levels(df_flagged$reliability),
    c("Reliable (< 20%)", "Use with Caution (20-30%)", "Unreliable (\u2265 30%)")
  )

  # Threshold validation errors
  expect_error(add_reliability_flags(fit_fh, thresholds = c(20)), "must be a numeric vector of length 2")
  expect_error(add_reliability_flags(fit_fh, thresholds = c(30, 20)), "with increasing values")

  # Null estimation data
  expect_error(add_reliability_flags(structure(list(x = 1), class = "fastsae")), "Could not extract estimation data")

  # Data frame without rse_col, but has mse and estimate
  df_mse <- data.frame(estimate = c(10, 20), mse = c(1, 4))
  res_mse <- add_reliability_flags(df_mse)
  expect_true("rse" %in% names(res_mse))
  expect_true("reliability" %in% names(res_mse))

  # Data frame with mse but no estimate column
  df_no_est <- data.frame(mse = c(1, 4))
  expect_error(add_reliability_flags(df_no_est), "missing estimate column")

  # Data frame without rse and without mse
  df_no_rse <- data.frame(val = c(1, 2))
  expect_error(add_reliability_flags(df_no_rse), "was not found in data")

  # Model with df_hb
  mod_hb <- structure(list(df_hb = data.frame(hb = 10, mse = 1)), class = "fastsae")
  res_hb <- add_reliability_flags(mod_hb)
  expect_true("reliability" %in% names(res_hb))

  # Model with df_ebp
  mod_ebp <- structure(list(df_ebp = data.frame(ebp = 10, mse = 1)), class = "fastsae")
  res_ebp <- add_reliability_flags(mod_ebp)
  expect_true("reliability" %in% names(res_ebp))
})

test_that("export_sae exports to Excel and CSV correctly with all options", {
  data(mys)
  fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys, print_result = FALSE)
  bm <- benchmark_sae(fit_fh, target = 6.5, method = "ratio")

  # 1. Export to Excel (.xlsx) if writexl or openxlsx installed
  if (requireNamespace("writexl", quietly = TRUE) || requireNamespace("openxlsx", quietly = TRUE)) {
    xlsx_file <- tempfile(fileext = ".xlsx")
    # attach stage1_summary to benchmark
    attr(bm, "stage1_summary") <- data.frame(Group = "Total", Shift = 0.05)
    res_xlsx <- export_sae(fit_fh, file = xlsx_file, benchmark = bm)

    expect_true(file.exists(xlsx_file))
    expect_true("Estimates" %in% names(res_xlsx))
    expect_true("Model_Summary" %in% names(res_xlsx))
    expect_true("Benchmarked" %in% names(res_xlsx))
    expect_true("Benchmarked_Groups" %in% names(res_xlsx))
    expect_equal(nrow(res_xlsx$Estimates), nrow(mys))
    expect_true("Reliability_Flag" %in% names(res_xlsx$Estimates))
    unlink(xlsx_file)
  }

  # 2. Export to CSV (.csv)
  csv_file <- tempfile(fileext = ".csv")
  res_csv <- export_sae(fit_fh, file = csv_file)

  expect_true(file.exists(csv_file))
  expect_true("Estimates" %in% names(res_csv))
  summary_csv <- paste0(tools::file_path_sans_ext(csv_file), "_summary.csv")
  expect_true(file.exists(summary_csv))
  unlink(csv_file)
  unlink(summary_csv)

  # 3. Export benchmark object directly
  bm_csv <- tempfile(fileext = ".csv")
  res_bm <- export_sae(bm, file = bm_csv)
  expect_true(file.exists(bm_csv))
  unlink(bm_csv)
  summary_bm <- paste0(tools::file_path_sans_ext(bm_csv), "_summary.csv")
  if (file.exists(summary_bm)) unlink(summary_bm)

  # 4. Export custom data frame with sd instead of mse
  custom_df <- data.frame(domain = c("A", "B"), estimate = c(5, 10), sd = c(0.5, 1.0))
  custom_csv <- tempfile(fileext = ".csv")
  res_custom <- export_sae(custom_df, file = custom_csv)
  expect_true(file.exists(custom_csv))
  unlink(custom_csv)

  # 5. Export comprehensive model summary with various attributes
  comp_model <- structure(
    list(
      model = "Spatial-Temporal",
      family = "gaussian",
      method = "REML",
      call = quote(test_call()),
      convergence = FALSE,
      n_domains = 10,
      random_effect_var = 1.23,
      rho = 0.45,
      rho_time = 0.67,
      phi = 0.89,
      goodness = c(AIC = 120.5, BIC = 135.2),
      estcoef = data.frame(beta = c(1, 2), std.error = c(0.1, 0.2), pvalue = c(0.01, 0.05), row.names = c("(Intercept)", "x1")),
      df_eblup = data.frame(domain = 1:5, y = 1:5, eblup = 1:5, mse = rep(0.5, 5), rse = rep(10, 5))
    ),
    class = "fastsae"
  )
  comp_csv <- tempfile(fileext = ".csv")
  res_comp <- export_sae(comp_model, file = comp_csv)
  expect_true(file.exists(comp_csv))
  unlink(comp_csv)
  unlink(paste0(tools::file_path_sans_ext(comp_csv), "_summary.csv"))
})

test_that("export_sae error handling", {
  data(mys)
  fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys, print_result = FALSE)

  # Invalid file path / missing
  expect_error(export_sae(fit_fh), "Please specify a valid file path")
  expect_error(export_sae(fit_fh, file = c("a.csv", "b.csv")), "Please specify a valid file path")
  expect_error(export_sae(fit_fh, file = 123), "Please specify a valid file path")

  # Invalid extension
  expect_error(export_sae(fit_fh, file = "invalid.txt"), "Target file extension must be either")

  # Overwrite = FALSE when file exists
  tmp <- tempfile(fileext = ".csv")
  writeLines("test", tmp)
  expect_error(export_sae(fit_fh, file = tmp, overwrite = FALSE), "already exists")
  unlink(tmp)

  # Cannot extract estimation table from invalid object
  expect_error(export_sae(list(a = 1), file = tempfile(fileext = ".csv")), "Cannot extract estimation table")
})

test_that("export_sae handles openxlsx branch and package missing fallback", {
  data(mys)
  fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys, print_result = FALSE)

  # 1. writexl missing, openxlsx available
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    tmp_xlsx <- tempfile(fileext = ".xlsx")
    testthat::with_mocked_bindings(
      {
        res <- export_sae(fit_fh, file = tmp_xlsx)
      },
      .has_pkg = function(pkg) {
        if (pkg == "writexl") return(FALSE)
        requireNamespace(pkg, quietly = TRUE)
      },
      .package = "fastsae"
    )
    expect_true(file.exists(tmp_xlsx))
    unlink(tmp_xlsx)
  }

  # 2. Both writexl and openxlsx missing -> fallback to CSV with warning
  tmp_xlsx2 <- tempfile(fileext = ".xlsx")
  expect_warning(
    testthat::with_mocked_bindings(
      {
        res <- export_sae(fit_fh, file = tmp_xlsx2)
      },
      .has_pkg = function(pkg) FALSE,
      .package = "fastsae"
    ),
    "Neither writexl nor openxlsx is installed"
  )
  fallback_csv <- paste0(tools::file_path_sans_ext(tmp_xlsx2), ".csv")
  expect_true(file.exists(fallback_csv))
  unlink(fallback_csv)
  if (file.exists(tmp_xlsx2)) unlink(tmp_xlsx2)
})

