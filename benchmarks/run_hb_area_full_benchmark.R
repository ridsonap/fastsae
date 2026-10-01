suppressPackageStartupMessages({
  library(fastsae)
  library(fastsaehb)
  library(tipsae)
  library(rstan)
  library(bench)
  library(tibble)
  library(dplyr)
  library(sf)
  library(spdep)
})

cat("==============================================================================\n")
cat("       BENCHMARK LENGKAP: fastsae::hb_area VS PACKAGE SAE HB LAINNYA          \n")
cat("==============================================================================\n\n")

# ------------------------------------------------------------------------------
# BAGIAN 1: PERBANDINGAN HASIL ESTIMASI & PARAMETER (MODEL NORMAL / FAY-HERRIOT)
# ------------------------------------------------------------------------------
cat("==============================================================================\n")
cat("BAGIAN 1: MODEL NORMAL (FAY-HERRIOT AREA LEVEL)\n")
cat("Perbandingan: fastsae (hb_area) vs fastsaehb (hb_normal Gibbs C++)\n")
cat("==============================================================================\n\n")

data("data_normal_sae", package = "fastsaehb")
df_norm <- data_normal_sae
D_norm <- nrow(df_norm)

cat(sprintf("Jumlah area (D): %d\n", D_norm))

# Fit fastsaehb (Gibbs MCMC)
set.seed(123)
t0 <- Sys.time()
fit_norm_mcmc <- fastsaehb::hb_normal(
  y ~ x1 + x2,
  vardir = "vardir",
  data = df_norm,
  n.iter = 10000,
  n.burnin = 2000,
  n.chains = 3,
  seed = 123
)
t_norm_mcmc <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# Fit fastsae (INLA)
t0 <- Sys.time()
fit_norm_fast <- fastsae::hb_area(
  y ~ x1 + x2,
  vardir = "vardir",
  data = df_norm,
  family = "gaussian",
  print_result = FALSE
)
t_norm_fast <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

est_norm_mcmc <- fit_norm_mcmc$summary_theta$HB_Mean
est_norm_fast <- fit_norm_fast$df_hb$hb
sd_norm_mcmc  <- fit_norm_mcmc$summary_theta$HB_SD
sd_norm_fast  <- fit_norm_fast$df_hb$sd

cor_norm <- cor(est_norm_mcmc, est_norm_fast)
mad_norm <- mean(abs(est_norm_mcmc - est_norm_fast))
max_norm <- max(abs(est_norm_mcmc - est_norm_fast))
rmse_norm <- sqrt(mean((est_norm_mcmc - est_norm_fast)^2))
cor_sd_norm <- cor(sd_norm_mcmc, sd_norm_fast)

cat(sprintf("Waktu eksekusi:\n  - fastsaehb (MCMC Gibbs 3 rantai, 10000 iter): %.3f detik\n  - fastsae (INLA Laplace): %.3f detik\n\n",
            t_norm_mcmc, t_norm_fast))

cat("Metrik Kesamaan Estimasi Domain (HB Area Estimates):\n")
cat(sprintf("  - Korelasi Pearson (r):       %.7f\n", cor_norm))
cat(sprintf("  - Mean Absolute Diff (MAD):   %.6f\n", mad_norm))
cat(sprintf("  - Max Absolute Diff:          %.6f\n", max_norm))
cat(sprintf("  - Root Mean Sq Diff (RMSE):   %.6f\n", rmse_norm))
cat(sprintf("  - Korelasi Posterior SD:      %.7f\n\n", cor_sd_norm))

cat("Perbandingan Koefisien Regresi (Beta):\n")
coef_comp_norm <- data.frame(
  Parameter = c("Intercept", "x1", "x2"),
  fastsae_est = fit_norm_fast$estcoef$beta,
  fastsae_se  = fit_norm_fast$estcoef$std.error,
  fastsaehb_est = fit_norm_mcmc$summary_beta$Mean,
  fastsaehb_se  = fit_norm_mcmc$summary_beta$SD
)
print(coef_comp_norm, row.names = FALSE)

cat("\nPerbandingan Varians Random Effect (sigma^2_u):\n")
cat(sprintf("  - fastsae (mode/mean posterior): %.5f\n", fit_norm_fast$random_effect_var))
cat(sprintf("  - fastsaehb (mean posterior):    %.5f\n\n", fit_norm_mcmc$summary_sigma2_u$Mean[1]))


# ------------------------------------------------------------------------------
# BAGIAN 2: MODEL BETA SAE NON-SPASIAL (DATA REAL EMILIA-ROMAGNA)
# ------------------------------------------------------------------------------
cat("==============================================================================\n")
cat("BAGIAN 2: MODEL BETA SAE NON-SPASIAL (DATA EMILIA-ROMAGNA D=38)\n")
cat("Perbandingan: fastsae (hb_area) vs tipsae (fit_sae Stan HMC) vs fastsaehb (hb_beta C++)\n")
cat("==============================================================================\n\n")

data("emilia_cs", package = "tipsae")
df_beta <- emilia_cs

# 1. tipsae (Stan HMC)
set.seed(42)
t0 <- Sys.time()
fit_beta_tipsae <- tipsae::fit_sae(
  formula_fixed = hcr ~ x,
  data = df_beta,
  domains = "id",
  disp_direct = "vars",
  type_disp = "var",
  likelihood = "beta",
  spatial_error = FALSE,
  temporal_error = FALSE,
  chains = 2,
  iter = 2000,
  refresh = 0,
  seed = 42
)
t_beta_tipsae <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# 2. fastsae (INLA)
t0 <- Sys.time()
fit_beta_fast <- fastsae::hb_area(
  formula = hcr ~ x,
  data = df_beta,
  domain = "id",
  vardir = "vars",
  family = "beta",
  print_result = FALSE
)
t_beta_fast <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# 3. fastsaehb (MCMC Gibbs C++)
t0 <- Sys.time()
fit_beta_fastsaehb <- fastsaehb::hb_beta(
  formula = hcr ~ x,
  vardir = "vars",
  data = df_beta,
  n.iter = 5000,
  n.burnin = 1000,
  seed = 42
)
t_beta_fastsaehb <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cat(sprintf("Waktu eksekusi:\n  - fastsae (INLA):     %.3f detik\n  - tipsae (Stan HMC):  %.3f detik\n  - fastsaehb (Gibbs):  %.3f detik\n\n",
            t_beta_fast, t_beta_tipsae, t_beta_fastsaehb))

# Ekstraksi dan penyelarasan estimasi
tip_res <- tipsae::extract(summary(fit_beta_tipsae))$in_sample
fast_res <- fit_beta_fast$df_hb
fastsaehb_res <- fit_beta_fastsaehb$summary_theta

merged_beta <- data.frame(
  domain = fast_res$domain,
  fastsae_hb = fast_res$hb,
  fastsae_sd = fast_res$sd,
  fastsaehb_hb = fastsaehb_res$HB_Mean,
  fastsaehb_sd = fastsaehb_res$HB_SD,
  stringsAsFactors = FALSE
)
merged_beta <- merge(merged_beta, tip_res[, c("Domains", "HB est.", "sd")], by.x = "domain", by.y = "Domains")
names(merged_beta)[names(merged_beta) == "HB est."] <- "tipsae_hb"
names(merged_beta)[names(merged_beta) == "sd"] <- "tipsae_sd"

cat("Perbandingan Estimasi Domain (Point Estimates HB):\n")
cat(sprintf("  - Korelasi fastsae vs tipsae:      %.7f\n", cor(merged_beta$fastsae_hb, merged_beta$tipsae_hb)))
cat(sprintf("  - MAD fastsae vs tipsae:           %.6f\n", mean(abs(merged_beta$fastsae_hb - merged_beta$tipsae_hb))))
cat(sprintf("  - Max Diff fastsae vs tipsae:      %.6f\n", max(abs(merged_beta$fastsae_hb - merged_beta$tipsae_hb))))
cat(sprintf("  - RMSE fastsae vs tipsae:          %.6f\n", sqrt(mean((merged_beta$fastsae_hb - merged_beta$tipsae_hb)^2))))
cat(sprintf("  - Korelasi fastsae vs fastsaehb:   %.7f\n", cor(merged_beta$fastsae_hb, merged_beta$fastsaehb_hb)))
cat(sprintf("  - MAD fastsae vs fastsaehb:        %.6f\n", mean(abs(merged_beta$fastsae_hb - merged_beta$fastsaehb_hb))))
cat(sprintf("  - Korelasi tipsae vs fastsaehb:    %.7f\n\n", cor(merged_beta$tipsae_hb, merged_beta$fastsaehb_hb)))

cat("Perbandingan Posterior SD:\n")
cat(sprintf("  - Korelasi SD fastsae vs tipsae:    %.7f\n", cor(merged_beta$fastsae_sd, merged_beta$tipsae_sd)))
cat(sprintf("  - MAD SD fastsae vs tipsae:         %.6f\n", mean(abs(merged_beta$fastsae_sd - merged_beta$tipsae_sd))))
cat(sprintf("  - Korelasi SD fastsae vs fastsaehb: %.7f\n\n", cor(merged_beta$fastsae_sd, merged_beta$fastsaehb_sd)))

cat("Perbandingan Koefisien Regresi (Beta):\n")
tip_stan_sum <- rstan::summary(fit_beta_tipsae$stanfit, pars = c("beta0", "beta"))$summary
coef_comp_beta <- data.frame(
  Parameter = c("Intercept", "x"),
  fastsae = fit_beta_fast$estcoef$beta,
  tipsae  = tip_stan_sum[, "mean"],
  fastsaehb = fit_beta_fastsaehb$summary_beta$Mean
)
print(coef_comp_beta, row.names = FALSE)


# ------------------------------------------------------------------------------
# BAGIAN 3: MODEL SPASIAL BESAG (CAR) PADA DATA BETA
# ------------------------------------------------------------------------------
cat("\n==============================================================================\n")
cat("BAGIAN 3: MODEL BETA SPASIAL (BESAG / CAR) PADA GRID D=60\n")
cat("Perbandingan: fastsae (hb_area besag) vs tipsae (fit_sae spatial_error=TRUE)\n")
cat("==============================================================================\n\n")

# Buat grid spasial D=60
D_sp <- 60L
set.seed(42)
poly <- sf::st_polygon(list(matrix(c(0,0, 1,0, 1,1, 0,1, 0,0), ncol = 2, byrow = TRUE)))
grid <- sf::st_make_grid(poly, n = c(6, 10))
domain_ids <- sprintf("area_%02d", seq_len(D_sp))
grid_sf <- sf::st_sf(domain = domain_ids, geometry = grid[1:D_sp])
nb <- spdep::poly2nb(grid_sf)
W <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
dimnames(W) <- list(domain_ids, domain_ids)

x_sp <- rnorm(D_sp, 1.0, 0.5)
vardir_sp <- runif(D_sp, 0.001, 0.003)
y_sp <- pmin(pmax(plogis(-1.2 + 0.5 * x_sp + rnorm(D_sp, 0, 0.2)), 0.01), 0.99)

df_sp <- data.frame(
  domain = domain_ids,
  y = y_sp,
  vardir = vardir_sp,
  x = x_sp,
  stringsAsFactors = FALSE
)

# 1. tipsae spatial
t0 <- Sys.time()
fit_sp_tipsae <- tipsae::fit_sae(
  formula_fixed = y ~ x,
  data = df_sp,
  domains = "domain",
  disp_direct = "vardir",
  type_disp = "var",
  likelihood = "beta",
  spatial_error = TRUE,
  spatial_df = grid_sf,
  domains_spatial_df = "domain",
  temporal_error = FALSE,
  chains = 2,
  iter = 1500,
  refresh = 0,
  seed = 42
)
t_sp_tipsae <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# 2. fastsae spatial
t0 <- Sys.time()
fit_sp_fast <- fastsae::hb_area(
  formula = y ~ x,
  data = df_sp,
  domain = "domain",
  vardir = "vardir",
  family = "beta",
  spatial = "besag",
  W = W,
  print_result = FALSE
)
t_sp_fast <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cat(sprintf("Waktu eksekusi Spasial (D=60):\n  - fastsae (INLA Besag):  %.3f detik\n  - tipsae (Stan Spasial): %.3f detik\n  - Speedup fastsae:       %.1fx lebih cepat\n\n",
            t_sp_fast, t_sp_tipsae, t_sp_tipsae / t_sp_fast))

tip_sp_res <- tipsae::extract(summary(fit_sp_tipsae))$in_sample
fast_sp_res <- fit_sp_fast$df_hb

merged_sp <- merge(fast_sp_res, tip_sp_res, by.x = "domain", by.y = "Domains")
cor_sp <- cor(merged_sp$hb, merged_sp[["HB est."]])
mad_sp <- mean(abs(merged_sp$hb - merged_sp[["HB est."]]))
cat(sprintf("Metrik Kesamaan Estimasi Spasial:\n  - Korelasi HB: %.7f\n  - MAD:          %.6f\n  - Korelasi SD:  %.7f\n\n",
            cor_sp, mad_sp, cor(merged_sp$sd.x, merged_sp$sd.y)))


# ------------------------------------------------------------------------------
# BAGIAN 4: BENCHMARK FORMAL WAKTU & MEMORI SKALABILITAS
# ------------------------------------------------------------------------------
cat("==============================================================================\n")
cat("BAGIAN 4: BENCHMARK FORMAL WAKTU & MEMORI (SCALING across Domain D)\n")
cat("==============================================================================\n\n")

if (file.exists("inst/extdata/beta_benchmark.rds")) {
  b_beta <- readRDS("inst/extdata/beta_benchmark.rds")
  cat("Tabel Benchmark Waktu & Memori (Beta SAE Cross-Sectional):\n")
  print(as.data.frame(b_beta))
}

if (file.exists("inst/extdata/spatial_beta_benchmark.rds")) {
  b_sp <- readRDS("inst/extdata/spatial_beta_benchmark.rds")
  cat("\nTabel Benchmark Waktu & Memori (Beta SAE Spasial):\n")
  print(as.data.frame(b_sp))
}
