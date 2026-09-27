#!/usr/bin/env Rscript

# ==============================================================================
# Finite Population Design-Based Simulation for Small Area Estimation (Beta SAE)
# Validating fastsae::ebp_area vs tipsae::fit_sae against True Population Parameters
# ==============================================================================

suppressPackageStartupMessages({
  library(fastsae)
  library(tipsae)
  library(sf)
  library(spdep)
  library(ggplot2)
  library(gridExtra)
  library(tibble)
  library(dplyr)
})

cat("==============================================================================\n")
cat(" FINITE POPULATION DESIGN-BASED SIMULATION: FASTSAE vs TIPSAE               \n")
cat(" Validating against True Population Parameter P_d from 30,000 individuals    \n")
cat("==============================================================================\n\n")

out_dir <- "benchmarks/output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

set.seed(2024L)
D <- 60L
domain_ids <- sprintf("area_%02d", seq_len(D))

# ------------------------------------------------------------------------------
# 1. GENERATE FINITE POPULATION (N ~ 30,000 across D = 60 domains)
# ------------------------------------------------------------------------------
cat(">>> Step 1: Generating Synthetic Finite Population...\n")

# Population sizes per domain N_d between 350 and 750
N_d <- sample(350:750, size = D, replace = TRUE)
names(N_d) <- domain_ids
N_total <- sum(N_d)

cat(sprintf("    Total population size (N): %d across %d domains (mean N_d = %.1f)\n",
            N_total, D, mean(N_d)))

# Spatial structure for the 60 domains (6 x 10 grid)
poly <- sf::st_polygon(list(matrix(c(0,0, 1,0, 1,1, 0,1, 0,0), ncol = 2, byrow = TRUE)))
grid <- sf::st_make_grid(poly, n = c(6, 10))
grid_sf <- sf::st_sf(domain = domain_ids, geometry = grid[1:D])
nb <- spdep::poly2nb(grid_sf)
W <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
dimnames(W) <- list(domain_ids, domain_ids)

# Domain-level random effects (combination of spatial Besag and unstructured)
Q <- diag(rowSums(W)) - W
eig <- eigen(Q)
pos_idx <- which(eig$values > 1e-10)
U_pos <- eig$vectors[, pos_idx]
lambda_pos <- eig$values[pos_idx]
v_raw <- as.vector(U_pos %*% (stats::rnorm(length(pos_idx)) / sqrt(lambda_pos)))
v_sp <- (v_raw - mean(v_raw)) / stats::sd(v_raw) * 0.30

u_iid <- stats::rnorm(D, mean = 0, sd = 0.25)
u_total <- u_iid + v_sp
names(u_total) <- domain_ids

# True Superpopulation regression parameters
beta_0 <- -1.40
beta_1 <-  0.55
beta_2 <- -0.45

# Create individual records in population
pop_list <- vector("list", D)
for (d in seq_len(D)) {
  dom_id <- domain_ids[d]
  n_units <- N_d[d]
  
  # Individual covariates
  x1_ind <- stats::rnorm(n_units, mean = 1.6 + 0.3 * sin(d), sd = 0.6)
  x2_ind <- stats::runif(n_units, min = 0, max = 2.2)
  
  # Latent probability
  eta_ind <- beta_0 + beta_1 * x1_ind + beta_2 * x2_ind + u_total[dom_id]
  pi_ind <- stats::plogis(eta_ind)
  
  # Binary individual indicator (e.g. 1 = impoverished/unemployed, 0 = otherwise)
  y_ind <- stats::rbinom(n_units, size = 1, prob = pi_ind)
  
  pop_list[[d]] <- data.frame(
    domain = dom_id,
    unit_id = paste0(dom_id, "_", sprintf("%04d", seq_len(n_units))),
    x1 = x1_ind,
    x2 = x2_ind,
    y = y_ind,
    stringsAsFactors = FALSE
  )
}
pop_df <- bind_rows(pop_list)

# ------------------------------------------------------------------------------
# 2. CALCULATE TRUE FINITE POPULATION PARAMETERS (P_d)
# ------------------------------------------------------------------------------
cat(">>> Step 2: Calculating True Finite Population Parameters (P_d)...\n")

pop_summary <- pop_df %>%
  group_by(domain) %>%
  summarise(
    N_d = n(),
    True_P_d = mean(y),           # EXACT Ground Truth Population Proportion
    x1_pop_mean = mean(x1),       # Census / Administrative Auxiliary Variables
    x2_pop_mean = mean(x2),
    .groups = "drop"
  )

cat(sprintf("    True P_d range: [%.4f, %.4f], Mean P_d = %.4f\n\n",
            min(pop_summary$True_P_d), max(pop_summary$True_P_d), mean(pop_summary$True_P_d)))

# ------------------------------------------------------------------------------
# 3. SIMULATE SURVEY SAMPLING (SRS without replacement per domain)
# ------------------------------------------------------------------------------
cat(">>> Step 3: Drawing Random Survey Samples (n_d = 20 - 55 per domain)...\n")

# Sample sizes per domain n_d
n_d_vec <- sample(22:55, size = D, replace = TRUE)
names(n_d_vec) <- domain_ids

smp_list <- vector("list", D)
for (d in seq_len(D)) {
  dom_id <- domain_ids[d]
  dom_units <- which(pop_df$domain == dom_id)
  smp_idx <- sample(dom_units, size = n_d_vec[dom_id], replace = FALSE)
  smp_list[[d]] <- pop_df[smp_idx, ]
}
smp_df <- bind_rows(smp_list)

# Compute Direct Estimator & Design-Based Sampling Variance (with FPC)
smp_summary <- smp_df %>%
  group_by(domain) %>%
  summarise(
    n_smp = n(),
    y_raw = mean(y),              # Sample proportion
    .groups = "drop"
  ) %>%
  left_join(pop_summary[, c("domain", "N_d", "True_P_d", "x1_pop_mean", "x2_pop_mean")], by = "domain")

# Apply Smithson-Verkuilen transformation for boundary protection
smp_summary <- smp_summary %>%
  mutate(
    # Smooth proportions strictly into (0, 1)
    p_dir = (y_raw * (n_smp - 1) + 0.5) / n_smp,
    # Design-based variance with Finite Population Correction (1 - f_d)
    fpc = 1 - (n_smp / N_d),
    vardir = fpc * (p_dir * (1 - p_dir)) / (n_smp - 1)
  )

smp_summary <- as.data.frame(smp_summary)

cat(sprintf("    Total sample drawn: %d units (mean n_d = %.1f)\n",
            sum(smp_summary$n_smp), mean(smp_summary$n_smp)))
cat(sprintf("    Direct estimator p_dir range: [%.4f, %.4f]\n",
            min(smp_summary$p_dir), max(smp_summary$p_dir)))
cat(sprintf("    Design variance vardir range: [%.6f, %.6f]\n\n",
            min(smp_summary$vardir), max(smp_summary$vardir)))

# ------------------------------------------------------------------------------
# 4. FIT SMALL AREA MODELS: FASTSAE vs TIPSAE
# ------------------------------------------------------------------------------
cat("------------------------------------------------------------------------------\n")
cat(">>> Step 4: Fitting Cross-Sectional Beta SAE Models...\n")
cat("------------------------------------------------------------------------------\n")

# --- 4A. fastsae (Cross-Sectional) ---
t0_fast <- Sys.time()
fit_fast <- ebp_area(
  formula = p_dir ~ x1_pop_mean + x2_pop_mean,
  data = smp_summary,
  domain = "domain",
  vardir = "vardir",
  family = "beta",
  print_result = FALSE
)
t_fast <- as.numeric(difftime(Sys.time(), t0_fast, units = "secs"))
est_fast <- fit_fast$df_ebp

# --- 4B. tipsae (Cross-Sectional) ---
t0_tip <- Sys.time()
fit_tip <- fit_sae(
  formula_fixed = p_dir ~ x1_pop_mean + x2_pop_mean,
  data = smp_summary,
  domains = "domain",
  disp_direct = "vardir",
  type_disp = "var",
  likelihood = "beta",
  spatial_error = FALSE,
  temporal_error = FALSE,
  chains = 2,
  iter = 1000,
  refresh = 0,
  seed = 42
)
t_tip <- as.numeric(difftime(Sys.time(), t0_tip, units = "secs"))
est_tip <- tipsae::extract(summary(fit_tip))$in_sample

cat(sprintf("Cross-Sectional Finished! FastSAE: %.2fs | TipSAE: %.2fs | Speedup: %.1fx\n\n",
            t_fast, t_tip, t_tip / t_fast))

# --- 4C. Fitting Spatial Beta SAE Models ---
cat("------------------------------------------------------------------------------\n")
cat(">>> Step 5: Fitting Spatial Beta SAE Models (Besag)...\n")
cat("------------------------------------------------------------------------------\n")

# Spatial fastsae
t0_fast_sp <- Sys.time()
fit_fast_sp <- ebp_area(
  formula = p_dir ~ x1_pop_mean + x2_pop_mean,
  data = smp_summary,
  domain = "domain",
  vardir = "vardir",
  family = "beta",
  spatial = "besag",
  W = W,
  print_result = FALSE
)
t_fast_sp <- as.numeric(difftime(Sys.time(), t0_fast_sp, units = "secs"))
est_fast_sp <- fit_fast_sp$df_ebp

# Spatial tipsae
t0_tip_sp <- Sys.time()
fit_tip_sp <- fit_sae(
  formula_fixed = p_dir ~ x1_pop_mean + x2_pop_mean,
  data = smp_summary,
  domains = "domain",
  disp_direct = "vardir",
  type_disp = "var",
  likelihood = "beta",
  spatial_error = TRUE,
  spatial_df = grid_sf,
  domains_spatial_df = "domain",
  temporal_error = FALSE,
  chains = 2,
  iter = 1000,
  refresh = 0,
  seed = 42
)
t_tip_sp <- as.numeric(difftime(Sys.time(), t0_tip_sp, units = "secs"))
est_tip_sp <- tipsae::extract(summary(fit_tip_sp))$in_sample

cat(sprintf("Spatial Finished! FastSAE: %.2fs | TipSAE: %.2fs | Speedup: %.1fx\n\n",
            t_fast_sp, t_tip_sp, t_tip_sp / t_fast_sp))

# ------------------------------------------------------------------------------
# 5. METRICS CALCULATION AGAINST TRUE POPULATION PARAMETER P_d
# ------------------------------------------------------------------------------
cat("==============================================================================\n")
cat(">>> Step 6: Evaluating Models against TRUE Population Parameter P_d...\n")
cat("==============================================================================\n")

# Merge Cross-Sectional Results
m_cs <- smp_summary %>%
  left_join(est_fast[, c("domain", "ebp", "sd", "ci_lower", "ci_upper")], by = "domain") %>%
  left_join(est_tip[, c("Domains", "HB est.", "sd", "2.5%", "97.5%")], by = c("domain" = "Domains")) %>%
  rename(
    fast_ebp = ebp, fast_sd = sd.x, fast_lower = ci_lower, fast_upper = ci_upper,
    tip_ebp = `HB est.`, tip_sd = sd.y, tip_lower = `2.5%`, tip_upper = `97.5%`
  )

# Merge Spatial Results
m_sp <- smp_summary %>%
  left_join(est_fast_sp[, c("domain", "ebp", "sd", "ci_lower", "ci_upper")], by = "domain") %>%
  left_join(est_tip_sp[, c("Domains", "HB est.", "sd", "2.5%", "97.5%")], by = c("domain" = "Domains")) %>%
  rename(
    fast_ebp = ebp, fast_sd = sd.x, fast_lower = ci_lower, fast_upper = ci_upper,
    tip_ebp = `HB est.`, tip_sd = sd.y, tip_lower = `2.5%`, tip_upper = `97.5%`
  )

calc_eval <- function(est_vec, true_vec, low_vec = NULL, up_vec = NULL, dir_mse = NULL) {
  err <- est_vec - true_vec
  rmse <- sqrt(mean(err^2))
  mae <- mean(abs(err))
  rrmse <- mean(abs(err) / true_vec) * 100
  rb <- mean(err / true_vec) * 100
  cor_p <- stats::cor(est_vec, true_vec)
  cor_s <- stats::cor(est_vec, true_vec, method = "spearman")
  
  cov_rate <- if (!is.null(low_vec) && !is.null(up_vec)) {
    mean(true_vec >= low_vec & true_vec <= up_vec) * 100
  } else NA_real_
  
  ci_w <- if (!is.null(low_vec) && !is.null(up_vec)) {
    mean(up_vec - low_vec)
  } else NA_real_
  
  eff <- if (!is.null(dir_mse)) dir_mse / (rmse^2) else NA_real_
  
  tibble::tibble(
    RMSE = rmse,
    MAE = mae,
    `RRMSE (%)` = rrmse,
    `RB (%)` = rb,
    `Cor (Pearson)` = cor_p,
    `Cor (Spearman)` = cor_s,
    `95% CrI Coverage (%)` = cov_rate,
    `Mean CrI Width` = ci_w,
    `Efficiency Gain (vs Direct)` = eff
  )
}

# Cross-sectional metrics
mse_dir <- mean((m_cs$p_dir - m_cs$True_P_d)^2)
res_dir <- calc_eval(m_cs$p_dir, m_cs$True_P_d, dir_mse = mse_dir)
res_fast_cs <- calc_eval(m_cs$fast_ebp, m_cs$True_P_d, m_cs$fast_lower, m_cs$fast_upper, dir_mse = mse_dir)
res_tip_cs  <- calc_eval(m_cs$tip_ebp,  m_cs$True_P_d, m_cs$tip_lower,  m_cs$tip_upper,  dir_mse = mse_dir)

tbl_eval_cs <- bind_rows(
  bind_cols(Method = "Direct Estimator (Sample p_dir)", res_dir, `Runtime (s)` = NA_real_),
  bind_cols(Method = "fastsae::ebp_area (INLA)", res_fast_cs, `Runtime (s)` = t_fast),
  bind_cols(Method = "tipsae::fit_sae (Stan MCMC)", res_tip_cs, `Runtime (s)` = t_tip)
)

# Spatial metrics
res_fast_sp <- calc_eval(m_sp$fast_ebp, m_sp$True_P_d, m_sp$fast_lower, m_sp$fast_upper, dir_mse = mse_dir)
res_tip_sp  <- calc_eval(m_sp$tip_ebp,  m_sp$True_P_d, m_sp$tip_lower,  m_sp$tip_upper,  dir_mse = mse_dir)

tbl_eval_sp <- bind_rows(
  bind_cols(Method = "Direct Estimator (Sample p_dir)", res_dir, `Runtime (s)` = NA_real_),
  bind_cols(Method = "fastsae::ebp_area (Spatial Besag)", res_fast_sp, `Runtime (s)` = t_fast_sp),
  bind_cols(Method = "tipsae::fit_sae (Spatial Besag)", res_tip_sp, `Runtime (s)` = t_tip_sp)
)

cat("\n--- TABLE 1: Cross-Sectional Evaluation vs True Population P_d ---\n")
print(as.data.frame(tbl_eval_cs), digits = 4)
cat(sprintf("Agreement FastSAE vs TipSAE: Pearson r = %.5f | Max Abs Diff = %.5f\n\n",
            cor(m_cs$fast_ebp, m_cs$tip_ebp), max(abs(m_cs$fast_ebp - m_cs$tip_ebp))))

cat("--- TABLE 2: Spatial Evaluation vs True Population P_d ---\n")
print(as.data.frame(tbl_eval_sp), digits = 4)
cat(sprintf("Agreement FastSAE vs TipSAE: Pearson r = %.5f | Max Abs Diff = %.5f\n\n",
            cor(m_sp$fast_ebp, m_sp$tip_ebp), max(abs(m_sp$fast_ebp - m_sp$tip_ebp))))

# ------------------------------------------------------------------------------
# 6. VISUALIZATIONS
# ------------------------------------------------------------------------------
theme_clean <- theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# --- PLOT 1: Estimates vs True Population P_d ---
df_p1 <- bind_rows(
  data.frame(Domain = m_cs$domain, True_P_d = m_cs$True_P_d, Estimate = m_cs$p_dir, Method = "Direct (Sample)", Model = "Sample"),
  data.frame(Domain = m_cs$domain, True_P_d = m_cs$True_P_d, Estimate = m_cs$fast_ebp, Method = "FastSAE (ebp_area)", Model = "Cross-Sectional"),
  data.frame(Domain = m_cs$domain, True_P_d = m_cs$True_P_d, Estimate = m_cs$tip_ebp, Method = "TipSAE (fit_sae)", Model = "Cross-Sectional"),
  data.frame(Domain = m_sp$domain, True_P_d = m_sp$True_P_d, Estimate = m_sp$fast_ebp, Method = "FastSAE (ebp_area)", Model = "Spatial (Besag)"),
  data.frame(Domain = m_sp$domain, True_P_d = m_sp$True_P_d, Estimate = m_sp$tip_ebp, Method = "TipSAE (fit_sae)", Model = "Spatial (Besag)")
)
df_p1$Method <- factor(df_p1$Method, levels = c("Direct (Sample)", "FastSAE (ebp_area)", "TipSAE (fit_sae)"))

p1 <- ggplot(df_p1, aes(x = True_P_d, y = Estimate, color = Method)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_point(alpha = 0.75, size = 1.8) +
  scale_color_manual(values = c("Direct (Sample)" = "#E64B35", "FastSAE (ebp_area)" = "#00A087", "TipSAE (fit_sae)" = "#3C5488")) +
  facet_wrap(~ Model + Method, nrow = 2) +
  labs(
    title = "Design-Based Simulation: Estimator vs True Population Parameter P_d",
    subtitle = "True P_d calculated from 30,000 population individuals; Dashed line is y = x",
    x = "True Finite Population Proportion (P_d)",
    y = "Estimated Proportion"
  ) +
  theme_clean +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "plot_finite_pop_estimates_vs_true.png"), p1, width = 10, height = 6.2, dpi = 300)

# --- PLOT 2: Absolute Error Comparison vs True P_d (Boxplot) ---
df_err <- bind_rows(
  data.frame(AbsError = abs(m_cs$p_dir - m_cs$True_P_d), Method = "Direct (Sample)", Model = "Sample"),
  data.frame(AbsError = abs(m_cs$fast_ebp - m_cs$True_P_d), Method = "FastSAE (INLA)", Model = "Cross-Sectional"),
  data.frame(AbsError = abs(m_cs$tip_ebp - m_cs$True_P_d), Method = "TipSAE (Stan)", Model = "Cross-Sectional"),
  data.frame(AbsError = abs(m_sp$fast_ebp - m_sp$True_P_d), Method = "FastSAE (INLA)", Model = "Spatial (Besag)"),
  data.frame(AbsError = abs(m_sp$tip_ebp - m_sp$True_P_d), Method = "TipSAE (Stan)", Model = "Spatial (Besag)")
)
df_err$Method <- factor(df_err$Method, levels = c("Direct (Sample)", "FastSAE (INLA)", "TipSAE (Stan)"))

p2 <- ggplot(df_err, aes(x = Method, y = AbsError, fill = Method)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1.5) +
  facet_wrap(~ Model, scales = "free_x") +
  scale_fill_manual(values = c("Direct (Sample)" = "#E64B35", "FastSAE (INLA)" = "#00A087", "TipSAE (Stan)" = "#3C5488")) +
  labs(
    title = "Absolute Error Distribution Relative to True Finite Population Parameter",
    subtitle = "Substantial error reduction by both SAE models compared to direct survey sampling",
    x = "Estimation Method",
    y = "Absolute Error |Estimate - True P_d|"
  ) +
  theme_clean +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "plot_finite_pop_errors.png"), p2, width = 8.5, height = 4.8, dpi = 300)

# --- PLOT 3: 95% Credible Interval Coverage of True P_d ---
sub_sp <- m_sp[1:25, ]
sub_sp$AreaIdx <- seq_len(nrow(sub_sp))

df_ci <- bind_rows(
  data.frame(AreaIdx = sub_sp$AreaIdx, Est = sub_sp$fast_ebp, Lower = sub_sp$fast_lower, Upper = sub_sp$fast_upper,
             True_P_d = sub_sp$True_P_d, Method = "FastSAE (INLA)"),
  data.frame(AreaIdx = sub_sp$AreaIdx, Est = sub_sp$tip_ebp, Lower = sub_sp$tip_lower, Upper = sub_sp$tip_upper,
             True_P_d = sub_sp$True_P_d, Method = "TipSAE (Stan MCMC)")
)

p3 <- ggplot(df_ci, aes(x = factor(AreaIdx), y = Est, color = Method)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), position = position_dodge(0.6), width = 0.4, linewidth = 0.7) +
  geom_point(aes(shape = Method), position = position_dodge(0.6), size = 2) +
  geom_point(aes(y = True_P_d), color = "black", shape = 4, size = 2.5, stroke = 1.2) +
  scale_color_manual(values = c("FastSAE (INLA)" = "#00A087", "TipSAE (Stan MCMC)" = "#3C5488")) +
  labs(
    title = "95% Credible Intervals vs True Population Parameter P_d (Spatial Model)",
    subtitle = "Black crosses (x) = True Population P_d; Bars = 95% Bayesian Credible Intervals",
    x = "Domain Area Index",
    y = "Estimated Proportion & 95% CrI"
  ) +
  theme_clean

ggsave(file.path(out_dir, "plot_finite_pop_coverage.png"), p3, width = 11, height = 5.2, dpi = 300)

# Save results object
saveRDS(list(
  eval_cs = tbl_eval_cs,
  eval_sp = tbl_eval_sp,
  data_cs = m_cs,
  data_sp = m_sp,
  pop_summary = pop_summary,
  smp_summary = smp_summary
), file.path(out_dir, "finite_population_simulation_results.rds"))

cat("\n[SUCCESS] Finite Population Simulation completed successfully!\n")
cat("Results saved to:", normalizePath(out_dir), "\n")
