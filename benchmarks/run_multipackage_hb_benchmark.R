suppressPackageStartupMessages({
  library(fastsae)
  library(hbsae)
  library(CARBayes)
  library(fastsaehb)
  library(bench)
  library(tibble)
  library(dplyr)
  library(sf)
  library(spdep)
})

cat("==============================================================================\n")
cat("    BENCHMARK KOMPREHENSIF MULTI-PACKAGE: fastsae::hb_area VS PACKAGES SAE HB\n")
cat("    Paket Pembanding: hbsae, CARBayes, fastsaehb, dan tipsae                 \n")
cat("==============================================================================\n\n")

# ==============================================================================
# 1. MODUL FAY-HERRIOT GAUSSIAN: fastsae vs hbsae vs fastsaehb
# ==============================================================================
cat("------------------------------------------------------------------------------\n")
cat("MODUL 1: GAUSSIAN FAY-HERRIOT (hb_area vs hbsae::fSAE.Area vs fastsaehb)\n")
cat("------------------------------------------------------------------------------\n")

data("data_normal_sae", package = "fastsaehb")
df_norm <- data_normal_sae
D_norm <- nrow(df_norm)
X_mat <- model.matrix(~ x1 + x2, data = df_norm)

# 1. fastsae (INLA)
t0 <- Sys.time()
fit_fs_norm <- fastsae::hb_area(
  y ~ x1 + x2,
  vardir = "vardir",
  data = df_norm,
  family = "gaussian",
  print_result = FALSE
)
t_fs_norm <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# 2. hbsae (Boonstra - Statistics Netherlands, Numerical Integration / Gibbs)
t0 <- Sys.time()
fit_hb_norm <- hbsae::fSAE.Area(
  est.init = df_norm$y,
  var.init = df_norm$vardir,
  X = X_mat,
  method = "HB"
)
t_hb_norm <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# 3. fastsaehb (C++ Gibbs Sampler MCMC)
t0 <- Sys.time()
fit_mcmc_norm <- fastsaehb::hb_normal(
  y ~ x1 + x2,
  vardir = "vardir",
  data = df_norm,
  n.iter = 10000,
  n.burnin = 2000,
  seed = 42
)
t_mcmc_norm <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

est_fs_norm <- fit_fs_norm$df_hb$hb
est_hb_norm <- as.numeric(hbsae::EST(fit_hb_norm))
est_mc_norm <- fit_mcmc_norm$summary_theta$HB_Mean

sd_fs_norm <- fit_fs_norm$df_hb$sd
sd_hb_norm <- as.numeric(hbsae::RMSE(fit_hb_norm))
sd_mc_norm <- fit_mcmc_norm$summary_theta$HB_SD

cat(sprintf("Runtime (D=%d):\n  - fastsae (INLA):     %.4f detik\n  - hbsae (Boonstra):   %.4f detik\n  - fastsaehb (Gibbs):  %.4f detik\n\n",
            D_norm, t_fs_norm, t_hb_norm, t_mcmc_norm))

cat("Metrik Kesamaan Estimasi Domain HB (Point Estimates):\n")
cat(sprintf("  - Korelasi fastsae vs hbsae:     %.7f\n", cor(est_fs_norm, est_hb_norm)))
cat(sprintf("  - MAD fastsae vs hbsae:          %.6f\n", mean(abs(est_fs_norm - est_hb_norm))))
cat(sprintf("  - RMSE fastsae vs hbsae:         %.6f\n", sqrt(mean((est_fs_norm - est_hb_norm)^2))))
cat(sprintf("  - Korelasi fastsae vs fastsaehb: %.7f\n", cor(est_fs_norm, est_mc_norm)))
cat(sprintf("  - MAD fastsae vs fastsaehb:      %.6f\n", mean(abs(est_fs_norm - est_mc_norm))))
cat(sprintf("  - Korelasi hbsae vs fastsaehb:   %.7f\n\n", cor(est_hb_norm, est_mc_norm)))

cat("Metrik Kesamaan Posterior Standard Error (SD/RMSE):\n")
cat(sprintf("  - Korelasi SD fastsae vs hbsae:     %.7f\n", cor(sd_fs_norm, sd_hb_norm)))
cat(sprintf("  - MAD SD fastsae vs hbsae:          %.6f\n", mean(abs(sd_fs_norm - sd_hb_norm))))
cat(sprintf("  - Korelasi SD fastsae vs fastsaehb: %.7f\n\n", cor(sd_fs_norm, sd_mc_norm)))

cat("Perbandingan Estimasi Koefisien Regresi (Beta):\n")
coef_table_norm <- data.frame(
  Parameter = c("Intercept", "x1", "x2"),
  fastsae_beta = fit_fs_norm$estcoef$beta,
  hbsae_beta   = as.numeric(fit_hb_norm$beta),
  fastsaehb_beta = fit_mcmc_norm$summary_beta$Mean
)
print(coef_table_norm, row.names = FALSE)


# ==============================================================================
# 2. MODUL POISSON SPASIAL (DISEASE MAPPING / COUNT SAE): fastsae vs CARBayes
# ==============================================================================
cat("\n------------------------------------------------------------------------------\n")
cat("MODUL 2: SPATIAL POISSON (hb_area vs CARBayes::S.CARbym)\n")
cat("------------------------------------------------------------------------------\n")

set.seed(123)
D_poi <- 50L
poly_poi <- sf::st_polygon(list(matrix(c(0,0, 1,0, 1,1, 0,1, 0,0), ncol = 2, byrow = TRUE)))
grid_poi <- sf::st_make_grid(poly_poi, n = c(5, 10))
grid_sf_poi <- sf::st_sf(domain = paste0("area_", seq_len(D_poi)), geometry = grid_poi[1:D_poi])
nb_poi <- spdep::poly2nb(grid_sf_poi)
W_poi <- spdep::nb2mat(nb_poi, style = "B", zero.policy = TRUE)

E_poi <- runif(D_poi, 80, 200)
x_poi <- rnorm(D_poi, 0, 1)
eta_true_poi <- -0.6 + 0.35 * x_poi + rnorm(D_poi, 0, 0.15)
y_poi <- rpois(D_poi, lambda = E_poi * exp(eta_true_poi))
df_poi <- data.frame(
  domain = paste0("area_", seq_len(D_poi)),
  y = y_poi,
  E = E_poi,
  x = x_poi,
  stringsAsFactors = FALSE
)

# 1. fastsae Poisson BYM
t0 <- Sys.time()
fit_fs_poi <- fastsae::hb_area(
  formula = y ~ x,
  data = df_poi,
  family = "poisson",
  exposure = "E",
  spatial = "bym",
  W = W_poi,
  print_result = FALSE
)
t_fs_poi <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# 2. CARBayes Poisson BYM (MCMC 20000 iterasi)
t0 <- Sys.time()
fit_cb_poi <- CARBayes::S.CARbym(
  formula = y ~ x + offset(log(E)),
  family = "poisson",
  data = df_poi,
  W = W_poi,
  burnin = 5000,
  n.sample = 20000,
  verbose = FALSE
)
t_cb_poi <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

rr_fs <- fit_fs_poi$df_hb$hb
rr_cb <- as.numeric(fit_cb_poi$fitted.values / df_poi$E)

cat(sprintf("Runtime Spatial Poisson (D=%d):\n  - fastsae (INLA BYM):    %.3f detik\n  - CARBayes (MCMC 20k):   %.3f detik\n  - Speedup fastsae:       %.1fx lebih cepat\n\n",
            D_poi, t_fs_poi, t_cb_poi, t_cb_poi / t_fs_poi))

cat("Metrik Kesamaan Relative Risk (RR = exp(eta)):\n")
cat(sprintf("  - Korelasi Pearson (r):     %.7f\n", cor(rr_fs, rr_cb)))
cat(sprintf("  - Mean Absolute Diff (MAD): %.6f\n", mean(abs(rr_fs - rr_cb))))
cat(sprintf("  - Max Absolute Diff:        %.6f\n", max(abs(rr_fs - rr_cb))))
cat(sprintf("  - RMSE:                     %.6f\n\n", sqrt(mean((rr_fs - rr_cb)^2))))

cat("Perbandingan Koefisien Regresi (Beta):\n")
coef_table_poi <- data.frame(
  Parameter = c("Intercept", "x"),
  fastsae = fit_fs_poi$estcoef$beta,
  CARBayes = fit_cb_poi$summary.results[c("(Intercept)", "x"), "Mean"]
)
print(coef_table_poi, row.names = FALSE)


# ==============================================================================
# 3. MODUL BINOMIAL SPASIAL (PROPORTION SAE): fastsae vs CARBayes
# ==============================================================================
cat("\n------------------------------------------------------------------------------\n")
cat("MODUL 3: SPATIAL BINOMIAL (hb_area vs CARBayes::S.CARbym)\n")
cat("------------------------------------------------------------------------------\n")

n_bin <- round(runif(D_poi, 100, 300))
prob_true_bin <- plogis(-0.9 + 0.4 * x_poi + rnorm(D_poi, 0, 0.2))
y_bin <- rbinom(D_poi, size = n_bin, prob = prob_true_bin)
df_bin <- data.frame(
  domain = paste0("area_", seq_len(D_poi)),
  y = y_bin,
  n = n_bin,
  x = x_poi,
  stringsAsFactors = FALSE
)

# 1. fastsae Binomial BYM
t0 <- Sys.time()
fit_fs_bin <- fastsae::hb_area(
  formula = y ~ x,
  data = df_bin,
  family = "binomial",
  trials = "n",
  spatial = "bym",
  W = W_poi,
  print_result = FALSE
)
t_fs_bin <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# 2. CARBayes Binomial BYM
t0 <- Sys.time()
fit_cb_bin <- CARBayes::S.CARbym(
  formula = y ~ x,
  family = "binomial",
  trials = df_bin$n,
  data = df_bin,
  W = W_poi,
  burnin = 5000,
  n.sample = 20000,
  verbose = FALSE
)
t_cb_bin <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

p_fs <- fit_fs_bin$df_hb$hb
p_cb <- as.numeric(fit_cb_bin$fitted.values / df_bin$n)

cat(sprintf("Runtime Spatial Binomial (D=%d):\n  - fastsae (INLA BYM):    %.3f detik\n  - CARBayes (MCMC 20k):   %.3f detik\n  - Speedup fastsae:       %.1fx lebih cepat\n\n",
            D_poi, t_fs_bin, t_cb_bin, t_cb_bin / t_fs_bin))

cat("Metrik Kesamaan Estimasi Proporsi (p = y/n):\n")
cat(sprintf("  - Korelasi Pearson (r):     %.7f\n", cor(p_fs, p_cb)))
cat(sprintf("  - Mean Absolute Diff (MAD): %.6f\n", mean(abs(p_fs - p_cb))))
cat(sprintf("  - Max Absolute Diff:        %.6f\n", max(abs(p_fs - p_cb))))
cat(sprintf("  - RMSE:                     %.6f\n\n", sqrt(mean((p_fs - p_cb)^2))))

cat("Perbandingan Koefisien Regresi (Beta):\n")
coef_table_bin <- data.frame(
  Parameter = c("Intercept", "x"),
  fastsae = fit_fs_bin$estcoef$beta,
  CARBayes = fit_cb_bin$summary.results[c("(Intercept)", "x"), "Mean"]
)
print(coef_table_bin, row.names = FALSE)


# ==============================================================================
# 4. BENCHMARK SKALABILITAS WAKTU & MEMORI: GAUSSIAN FH across D = 100, 250, 500
# ==============================================================================
cat("\n==============================================================================\n")
cat("MODUL 4: FORMAL BENCHMARK WAKTU & MEMORI (Gaussian FH across Domain sizes)\n")
cat("==============================================================================\n\n")

d_sizes <- c(100, 250, 500)
bench_results <- list()

for (D_val in d_sizes) {
  set.seed(42 + D_val)
  x1_v <- rnorm(D_val)
  x2_v <- runif(D_val, 10, 20)
  vardir_v <- runif(D_val, 0.05, 0.25)
  u_v <- rnorm(D_val, 0, sqrt(0.2))
  y_v <- 1.5 + 0.8 * x1_v - 0.1 * x2_v + u_v + rnorm(D_val, 0, sqrt(vardir_v))
  df_v <- data.frame(y = y_v, vardir = vardir_v, x1 = x1_v, x2 = x2_v, domain = paste0("area_", seq_len(D_val)))
  X_v <- model.matrix(~ x1 + x2, data = df_v)
  
  cat(sprintf("Running bench::mark for D = %d (5 iterasi each)...\n", D_val))
  res <- bench::mark(
    fastsae = fastsae::hb_area(y ~ x1 + x2, vardir = "vardir", data = df_v, family = "gaussian", print_result = FALSE),
    hbsae   = hbsae::fSAE.Area(est.init = df_v$y, var.init = df_v$vardir, X = X_v, method = "HB"),
    fastsaehb = fastsaehb::hb_normal(y ~ x1 + x2, vardir = "vardir", data = df_v, n.iter = 5000, n.burnin = 1000),
    check = FALSE,
    iterations = 5,
    memory = TRUE,
    time_unit = "s"
  )
  
  bench_df <- tibble::tibble(
    D = D_val,
    Method = as.character(res$expression),
    `Min (s)` = as.numeric(res$min),
    `Median (s)` = as.numeric(res$median),
    `Memory (MB)` = as.numeric(res$mem_alloc) / (1024^2)
  )
  bench_results[[as.character(D_val)]] <- bench_df
}

final_bench_df <- do.call(rbind, bench_results)
print(as.data.frame(final_bench_df))

cat("\n==============================================================================\n")
cat("BENCHMARK MULTI-PACKAGE SELESAI!\n")
cat("==============================================================================\n")
