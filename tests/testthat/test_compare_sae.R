test_that("compare_sae works for two fitted models", {
  data(mys)
  fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)
  fit_sfh <- eblup_sfh(y ~ x1 + x2, vardir = "vardir", data = mys, W = mys_proxmat)

  comp <- compare_sae(fit_fh, fit_sfh, names = c("FH", "Spatial FH"))

  expect_s3_class(comp, "fastsae_comparison")
  expect_equal(comp$names, c("FH", "Spatial FH"))
  expect_equal(nrow(comp$data), nrow(mys))
  expect_true("Metric" %in% names(comp$metrics))
  expect_true("Value" %in% names(comp$metrics))

  # Test print and summary
  expect_output(print(comp), "Key Concordance & Efficiency Metrics")
  expect_output(summary(comp), "Distribution of Domain Point Estimates")

  # Test autoplots
  p_scat <- autoplot(comp, type = "scatter")
  expect_s3_class(p_scat, "ggplot")

  p_diff <- autoplot(comp, type = "difference")
  expect_s3_class(p_diff, "ggplot")

  p_mse <- autoplot(comp, type = "mse")
  expect_s3_class(p_mse, "ggplot")

  p_rse <- autoplot(comp, type = "rse")
  expect_s3_class(p_rse, "ggplot")

  p_cmp <- autoplot(comp, type = "comparison")
  expect_s3_class(p_cmp, "ggplot")

  # Test plot.fastsae_comparison
  p_plot <- plot(comp, type = "scatter")
  expect_s3_class(p_plot, "ggplot")
})

test_that("compare_sae works with a list of two models", {
  data(mys)
  fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)
  fit_sfh <- eblup_sfh(y ~ x1 + x2, vardir = "vardir", data = mys, W = mys_proxmat)

  comp <- compare_sae(list("FH_Model" = fit_fh, "SFH_Model" = fit_sfh))
  expect_s3_class(comp, "fastsae_comparison")
  expect_equal(comp$names, c("FH_Model", "SFH_Model"))
})

test_that("compare_sae error handling", {
  data(mys)
  fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)

  expect_error(compare_sae(fit_fh), "Please provide two models")
})
