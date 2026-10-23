#!/usr/bin/env Rscript

# ==============================================================================
# Concordance: fastsae (hb_area) vs tipsae (fit_sae) for a Beta model, D = 1000
#
# Same data as the D = 1000 row of bench_beta_tipsae.R: sim_area_data(D = 1000,
# spatial_type = "knn", seed = 1042), with 6 unsampled areas. The Beta sampling
# variance is computed from the direct estimate, y (1 - y) / (kappa + 1) with the
# precision kappa = 25 of the simulator. The true area mean is
# mu = plogis(eta - 1), truncated to [0.02, 0.98] as in sim_area_data(), with
# eta = truth_gaussian.
#
# Both packages run with their defaults: hb_area() with INLA, the two-stage
# precision and its default priors; tipsae with type_disp = "var" and Stan's
# defaults (4 chains x 2000 iterations, 1000 warmup). The priors differ
# (half-normal(1) vs half-normal(2.5) on the area SD; N(0, 1000) vs N(0, 2.5^2)
# on the slopes), as does the precision (fixed at the posterior mean of a first
# fit vs a function of the model mean inside the model), so the comparison is
# between what a user gets from each package, not between two
# approximations of one posterior.
#
# Writes inst/extdata/beta_concordance_n1000.rds with
#   $est      domain, sampled, y, truth and, per package, mean, sd, lower, upper
#   $metrics  per package: RMSE and coverage of the 95% intervals against the
#             truth, mean interval width, and the agreement of the two packages
#   $meta     date, versions and settings
#
# Usage: Rscript benchmarks/bench_beta_concordance.R [--D=1000]
# ==============================================================================

suppressPackageStartupMessages({
  library(fastsae)
  library(tipsae)
})

args <- commandArgs(trailingOnly = TRUE)
D_arg <- grep("^--D=", args, value = TRUE)
D <- if (length(D_arg) > 0) as.integer(sub("^--D=", "", D_arg[1])) else 1000L

sim <- sim_area_data(D = D, spatial_type = "knn", seed = 42 + D)
df <- sim$data
df$domain <- as.character(df$domain)
df$vardir_beta <- df$y_beta * (1 - df$y_beta) / 26
truth <- pmin(pmax(stats::plogis(df$truth_gaussian - 1), 0.02), 0.98)

t_fastsae <- system.time(
  fit_f <- hb_area(y_beta ~ x1 + x2, data = df, domain = "domain", vardir = "vardir_beta",
                   family = "beta", print_result = FALSE)
)[["elapsed"]]
est_f <- fit_f$df_hb[match(df$domain, fit_f$df_hb$domain), ]

t_tipsae <- system.time(
  fit_t <- tipsae::fit_sae(y_beta ~ x1 + x2, data = df, domains = "domain",
                           disp_direct = "vardir_beta", type_disp = "var", likelihood = "beta",
                           spatial_error = FALSE, temporal_error = FALSE, refresh = 0, seed = 42)
)[["elapsed"]]
s_t <- summary(fit_t)
est_t <- rbind(s_t$model_estimates, s_t$model_estimates_oos)
est_t <- est_t[match(df$domain, as.character(est_t$Domains)), ]
rhat_tipsae <- max(rstan::summary(fit_t$stanfit)$summary[, "Rhat"], na.rm = TRUE)

est <- data.frame(
  domain = df$domain, sampled = !is.na(df$y_beta), y = df$y_beta, truth = truth,
  fastsae_mean = est_f$hb, fastsae_sd = est_f$sd,
  fastsae_lower = est_f$ci_lower, fastsae_upper = est_f$ci_upper,
  tipsae_mean = est_t$mean, tipsae_sd = est_t$sd,
  tipsae_lower = est_t[["2.5%"]], tipsae_upper = est_t[["97.5%"]]
)

pkg_metrics <- function(p) {
  m <- est[[paste0(p, "_mean")]]
  lo <- est[[paste0(p, "_lower")]]
  hi <- est[[paste0(p, "_upper")]]
  data.frame(
    package = p,
    rmse = sqrt(mean((m - est$truth)^2)),
    coverage = mean(lo <= est$truth & est$truth <= hi),
    width = mean(hi - lo)
  )
}
metrics <- do.call(rbind, lapply(c("fastsae", "tipsae"), pkg_metrics))
metrics$time <- c(t_fastsae, t_tipsae)
agreement <- data.frame(
  cor_mean = stats::cor(est$fastsae_mean, est$tipsae_mean),
  mad_mean = mean(abs(est$fastsae_mean - est$tipsae_mean)),
  ratio_sd = stats::median(est$fastsae_sd / est$tipsae_sd),
  ratio_sd_q05 = unname(stats::quantile(est$fastsae_sd / est$tipsae_sd, 0.05)),
  ratio_sd_q95 = unname(stats::quantile(est$fastsae_sd / est$tipsae_sd, 0.95)),
  rhat_tipsae = rhat_tipsae
)

meta <- list(
  date = as.character(Sys.Date()), D = D,
  r_version = R.version.string,
  fastsae = as.character(utils::packageVersion("fastsae")),
  inla = as.character(utils::packageVersion("INLA")),
  tipsae = as.character(utils::packageVersion("tipsae")),
  cpu = tryCatch(system2("sysctl", c("-n", "machdep.cpu.brand_string"), stdout = TRUE),
                 error = function(e) NA_character_)
)

dest <- file.path("inst", "extdata", "beta_concordance_n1000.rds")
saveRDS(list(est = est, metrics = metrics, agreement = agreement, meta = meta), dest)
print(metrics)
print(agreement)
cat("Saved", dest, "\n")
