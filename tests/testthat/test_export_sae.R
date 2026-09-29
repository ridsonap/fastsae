test_that("add_reliability_flags adds correct tiers", {
  data(mys)
  fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)

  df_flagged <- add_reliability_flags(fit_fh)
  expect_true("reliability" %in% names(df_flagged))
  expect_s3_class(df_flagged$reliability, "factor")
  expect_equal(
    levels(df_flagged$reliability),
    c("Reliable (< 20%)", "Use with Caution (20-30%)", "Unreliable (\u2265 30%)")
  )
})

test_that("export_sae exports to Excel and CSV correctly", {
  skip_if_not_installed("writexl")

  data(mys)
  fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)
  bm <- benchmark_sae(fit_fh, target = 6.5, method = "ratio")

  # 1. Export to Excel (.xlsx)
  xlsx_file <- tempfile(fileext = ".xlsx")
  res_xlsx <- export_sae(fit_fh, file = xlsx_file, benchmark = bm)

  expect_true(file.exists(xlsx_file))
  expect_true("Estimates" %in% names(res_xlsx))
  expect_true("Model_Summary" %in% names(res_xlsx))
  expect_true("Benchmarked" %in% names(res_xlsx))
  expect_equal(nrow(res_xlsx$Estimates), nrow(mys))
  expect_true("Reliability_Flag" %in% names(res_xlsx$Estimates))
  unlink(xlsx_file)

  # 2. Export to CSV (.csv)
  csv_file <- tempfile(fileext = ".csv")
  res_csv <- export_sae(fit_fh, file = csv_file)

  expect_true(file.exists(csv_file))
  expect_true("Estimates" %in% names(res_csv))
  unlink(csv_file)
})

test_that("export_sae error handling", {
  data(mys)
  fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)

  expect_error(export_sae(fit_fh, file = "invalid.txt"), "Target file extension must be either")
  expect_error(export_sae(fit_fh), "Please specify a valid file path")
})
