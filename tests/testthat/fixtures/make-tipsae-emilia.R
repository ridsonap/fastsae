# Builds tipsae-emilia.rds: the tipsae (Stan) posterior for the spatio-temporal
# Beta model on tipsae's `emilia` data, the reference for test_hb_temporal.R.
# Stan draws differ across platforms and one run can fail to mix, so the
# reference is computed once, checked for convergence and stored.
# Run from the package root: Rscript tests/testthat/fixtures/make-tipsae-emilia.R
# (about 20 minutes on 4 cores). Set TIPSAE_FIT_CACHE to an .rds path to keep
# the Stan fit between runs.

data("emilia", package = "tipsae")
data("emilia_shp", package = "tipsae")

cache <- Sys.getenv("TIPSAE_FIT_CACHE")
if (nzchar(cache) && file.exists(cache)) {
  fit <- readRDS(cache)
} else {
  fit <- tipsae::fit_sae(
    formula_fixed = hcr ~ x,
    data = emilia,
    domains = "id",
    disp_direct = "vars",
    type_disp = "var",
    domain_size = "n",
    spatial_error = TRUE,
    spatial_df = emilia_shp,
    domains_spatial_df = "NAME_DISTRICT",
    temporal_error = TRUE,
    temporal_variable = "year",
    chains = 4,
    iter = 10000,
    warmup = 2000,
    cores = 4,
    seed = 42
  )
  if (nzchar(cache)) saveRDS(fit, cache)
}

# tipsae imposes sum-to-zero on the spatial and temporal effects with a tight
# soft constraint, so NUTS hits the maximum treedepth and the raw effects and
# variances mix slowly. The stored quantity is theta (the area means), which is
# identified; require it to have converged.
smry <- rstan::summary(fit$stanfit)$summary
group <- sub("\\[.*$", "", rownames(smry))
print(t(sapply(split(seq_along(group), group), function(i) {
  c(max_rhat = max(smry[i, "Rhat"], na.rm = TRUE), min_ess = min(smry[i, "n_eff"], na.rm = TRUE))
})))
is_theta <- group == "theta"
rhat <- max(smry[is_theta, "Rhat"])
ess <- min(smry[is_theta, "n_eff"])
if (rhat > 1.01 || ess < 400) {
  stop(sprintf("theta did not converge (max R-hat %.3f, min ESS %.0f)", rhat, ess))
}

est <- tipsae::extract(summary(fit))$in_sample
stopifnot(c("Domains", "Times", "HB est.", "sd") %in% names(est))
ref <- data.frame(
  id = est[["Domains"]],
  year = est[["Times"]],
  hb = est[["HB est."]],
  sd = est[["sd"]]
)
attr(ref, "info") <- list(
  tipsae = as.character(utils::packageVersion("tipsae")),
  rstan = as.character(utils::packageVersion("rstan")),
  chains = 4, iter = 10000, warmup = 2000, seed = 42, max_rhat_theta = rhat, min_ess_theta = ess,
  max_rhat_all = max(smry[, "Rhat"], na.rm = TRUE),
  created = Sys.Date()
)
saveRDS(ref, "tests/testthat/fixtures/tipsae-emilia.rds", version = 2)
