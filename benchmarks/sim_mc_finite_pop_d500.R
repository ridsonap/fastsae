#!/usr/bin/env Rscript

# ==============================================================================
# Large-Scale National Monte Carlo Simulation: D = 500 Domains, N = 250,000 Units
# R = 25 Replications comparing fastsae::ebp_area vs tipsae::fit_sae against True P_d
# ==============================================================================

suppressPackageStartupMessages({
  library(fastsae)
  library(tipsae)
  library(sf)
  library(spdep)
  library(parallel)
  library(ggplot2)
  library(tibble)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
R_val <- 25L
r_match <- grep("^--R=", args, value = TRUE)
if (length(r_match) > 0) {
  R_val <- as.integer(sub("^--R=", "", r_match[1]))
}

cores_val <- min(5L, max(1L, parallel::detectCores() - 2L))
core_match <- grep("^--cores=", args, value = TRUE)
if (length(core_match) > 0) {
  cores_val <- as.integer(sub("^--cores=", "", core_match[1]))
}

out_dir <- "benchmarks/output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("==============================================================================\n")
cat(sprintf(" NATIONAL-SCALE MONTE CARLO SIMULATION: D = 500 Domains (R = %d Replications)\n", R_val))
cat(sprintf(" Parallel Workers: %d cores | Total Population: N ~ 250,000 Units\n", cores_val))
cat("==============================================================================\n\n")

# ------------------------------------------------------------------------------
# 1. FIXED SYNTHETIC NATIONAL POPULATION (D = 500, N ~ 250,000)
# ------------------------------------------------------------------------------
set.seed(999L)
D <- 500L
domain_ids <- sprintf("area_%03d", seq_len(D))

# Population size per domain N_d between 400 and 600 (Mean = 500)
N_d <- sample(400:600, size = D, replace = TRUE)
names(N_d) <- domain_ids
N_total <- sum(N_d)

# Spatial grid: 20 x 25 = 500 contiguous polygons
poly <- sf::st_polygon(list(matrix(c(0,0, 1,0, 1,1, 0,1, 0,0), ncol = 2, byrow = TRUE)))
grid <- sf::st_make_grid(poly, n = c(20, 25))
grid_sf <- sf::st_sf(domain = domain_ids, geometry = grid[1:D])
nb <- spdep::poly2nb(grid_sf)
W <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
dimnames(W) <- list(domain_ids, domain_ids)

# Spatial Besag + IID Random Effects across 500 domains
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

beta_0 <- -1.40
beta_1 <-  0.55
beta_2 <- -0.45

cat(">>> Generating 250,000 unit records across 500 domains...\n")
pop_list <- vector("list", D)
for (d in seq_len(D)) {
  dom_id <- domain_ids[d]
  n_units <- N_d[d]
  x1_ind <- stats::rnorm(n_units, mean = 1.6 + 0.2 * sin(d / 10), sd = 0.6)
  x2_ind <- stats::runif(n_units, min = 0, max = 2.2)
  eta_ind <- beta_0 + beta_1 * x1_ind + beta_2 * x2_ind + u_total[dom_id]
  pi_ind <- stats::plogis(eta_ind)
  y_ind <- stats::rbinom(n_units, size = 1, prob = pi_ind)
  
  pop_list[[d]] <- data.frame(
    domain = dom_id,
    x1 = x1_ind,
    x2 = x2_ind,
    y = y_ind,
    stringsAsFactors = FALSE
  )
}
pop_df <- bind_rows(pop_list)

# True Finite Population Quantities (Fixed Ground Truth)
pop_summary <- pop_df %>%
  group_by(domain) %>%
  summarise(
    N_d = n(),
    True_P_d = mean(y),
    x1_pop_mean = mean(x1),
    x2_pop_mean = mean(x2),
    .groups = "drop"
  ) %>%
  as.data.frame()
rownames(pop_summary) <- pop_summary$domain

True_P <- pop_summary$True_P_d
names(True_P) <- pop_summary$domain

cat(sprintf("    Total population size (N): %d across %d domains\n", N_total, D))
cat(sprintf("    True P_d Mean: %.4f (range [%.4f, %.4f])\n\n",
            mean(True_P), min(True_P), max(True_P)))

# ------------------------------------------------------------------------------
# 2. REPLICATION WORKER FUNCTION (D = 500)
# ------------------------------------------------------------------------------
run_one_rep_d500 <- function(r) {
  rep_seed <- 8000L + r
  set.seed(rep_seed)
  
  # Sample sizes n_d in [20, 50] per domain (total n ~ 17,500)
  n_d_vec <- sample(20:50, size = D, replace = TRUE)
  names(n_d_vec) <- domain_ids
  
  smp_list <- vector("list", D)
  for (d in seq_len(D)) {
    dom_id <- domain_ids[d]
    dom_units <- which(pop_df$domain == dom_id)
    smp_idx <- sample(dom_units, size = n_d_vec[dom_id], replace = FALSE)
    smp_list[[d]] <- pop_df[smp_idx, ]
  }
  smp_df <- bind_rows(smp_list)
  
  smp_summary <- smp_df %>%
    group_by(domain) %>%
    summarise(
      n_smp = n(),
      y_raw = mean(y),
      .groups = "drop"
    ) %>%
    left_join(pop_summary[, c("domain", "N_d", "True_P_d", "x1_pop_mean", "x2_pop_mean")], by = "domain") %>%
    mutate(
      p_dir = (y_raw * (n_smp - 1) + 0.5) / n_smp,
      fpc = 1 - (n_smp / N_d),
      vardir = fpc * (p_dir * (1 - p_dir)) / (n_smp - 1)
    )
  smp_summary <- as.data.frame(smp_summary)
  rownames(smp_summary) <- smp_summary$domain
  
  # A. FastSAE (Spatial Besag on D = 500)
  t0_fast <- Sys.time()
  fit_fast <- fastsae::ebp_area(
    formula = p_dir ~ x1_pop_mean + x2_pop_mean,
    data = smp_summary,
    domain = "domain",
    vardir = "vardir",
    family = "beta",
    spatial = "besag",
    W = W,
    print_result = FALSE
  )
  time_fast <- as.numeric(difftime(Sys.time(), t0_fast, units = "secs"))
  df_fast <- fit_fast$df_ebp
  rownames(df_fast) <- df_fast$domain
  df_fast <- df_fast[domain_ids, ]
  
  # B. TipSAE (Spatial Besag on D = 500)
  t0_tip <- Sys.time()
  fit_tip <- tipsae::fit_sae(
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
    chains = 1,
    iter = 500,
    warmup = 150,
    refresh = 0,
    seed = rep_seed
  )
  time_tip <- as.numeric(difftime(Sys.time(), t0_tip, units = "secs"))
  df_tip <- tipsae::extract(summary(fit_tip))$in_sample
  rownames(df_tip) <- df_tip$Domains
  df_tip <- df_tip[domain_ids, ]
  
  list(
    r = r,
    time_fast = time_fast,
    time_tip = time_tip,
    p_dir = smp_summary[domain_ids, "p_dir"],
    fast_est = df_fast$ebp,
    fast_low = df_fast$ci_lower,
    fast_up  = df_fast$ci_upper,
    tip_est  = df_tip$`HB est.`,
    tip_low  = df_tip$`2.5%`,
    tip_up   = df_tip$`97.5%`
  )
}

# ------------------------------------------------------------------------------
# 3. RUN PARALLEL MONTE CARLO (D = 500, R = 25)
# ------------------------------------------------------------------------------
cat(sprintf(">>> Starting Monte Carlo Replications (R = %d) for D = 500 domains...\n", R_val))
t_mc_start <- Sys.time()

mc_results <- parallel::mclapply(
  seq_len(R_val),
  run_one_rep_d500,
  mc.cores = cores_val,
  mc.preschedule = FALSE
)

t_mc_elapsed <- as.numeric(difftime(Sys.time(), t_mc_start, units = "secs"))
cat(sprintf("\n>>> Monte Carlo Completed in %.1f seconds (%.2f minutes)!\n\n",
            t_mc_elapsed, t_mc_elapsed / 60))

# ------------------------------------------------------------------------------
# 4. AGGREGATE EMPIRICAL PROPERTIES (500 DOMAINS x 25 REPLICATIONS)
# ------------------------------------------------------------------------------
cat(">>> Aggregating Empirical Properties across 500 domains...\n")

mat_dir  <- matrix(NA_real_, nrow = D, ncol = R_val, dimnames = list(domain_ids, NULL))
mat_fast <- matrix(NA_real_, nrow = D, ncol = R_val, dimnames = list(domain_ids, NULL))
mat_tip  <- matrix(NA_real_, nrow = D, ncol = R_val, dimnames = list(domain_ids, NULL))

cov_fast <- matrix(NA_real_, nrow = D, ncol = R_val, dimnames = list(domain_ids, NULL))
cov_tip  <- matrix(NA_real_, nrow = D, ncol = R_val, dimnames = list(domain_ids, NULL))

width_fast <- matrix(NA_real_, nrow = D, ncol = R_val, dimnames = list(domain_ids, NULL))
width_tip  <- matrix(NA_real_, nrow = D, ncol = R_val, dimnames = list(domain_ids, NULL))

times_fast <- numeric(R_val)
times_tip  <- numeric(R_val)

for (r in seq_len(R_val)) {
  res_r <- mc_results[[r]]
  if (is.null(res_r) || inherits(res_r, "try-error")) next
  
  times_fast[r] <- res_r$time_fast
  times_tip[r]  <- res_r$time_tip
  
  mat_dir[, r]  <- res_r$p_dir
  mat_fast[, r] <- res_r$fast_est
  mat_tip[, r]  <- res_r$tip_est
  
  cov_fast[, r] <- as.numeric(True_P >= res_r$fast_low & True_P <= res_r$fast_up)
  cov_tip[, r]  <- as.numeric(True_P >= res_r$tip_low  & True_P <= res_r$tip_up)
  
  width_fast[, r] <- res_r$fast_up - res_r$fast_low
  width_tip[, r]  <- res_r$tip_up - res_r$tip_low
}

calc_dom_metrics <- function(est_mat, cov_mat = NULL, width_mat = NULL) {
  emp_mean <- rowMeans(est_mat, na.rm = TRUE)
  emp_bias <- emp_mean - True_P
  emp_rb <- (emp_bias / True_P) * 100
  emp_marb <- rowMeans(abs(sweep(est_mat, 1, True_P, "-")) / True_P, na.rm = TRUE) * 100
  emp_mse <- rowMeans((sweep(est_mat, 1, True_P, "-"))^2, na.rm = TRUE)
  emp_rmse <- sqrt(emp_mse)
  emp_rrmse <- (emp_rmse / True_P) * 100
  emp_cr <- if (!is.null(cov_mat)) rowMeans(cov_mat, na.rm = TRUE) * 100 else rep(NA_real_, D)
  emp_w <- if (!is.null(width_mat)) rowMeans(width_mat, na.rm = TRUE) else rep(NA_real_, D)
  
  data.frame(
    domain = domain_ids,
    True_P = True_P,
    Emp_Mean = emp_mean,
    Emp_Bias = emp_bias,
    Emp_RB = emp_rb,
    Emp_MARB = emp_marb,
    Emp_MSE = emp_mse,
    Emp_RMSE = emp_rmse,
    Emp_RRMSE = emp_rrmse,
    Emp_CR = emp_cr,
    Emp_Width = emp_w,
    stringsAsFactors = FALSE
  )
}

dom_dir  <- calc_dom_metrics(mat_dir)
dom_fast <- calc_dom_metrics(mat_fast, cov_fast, width_fast)
dom_tip  <- calc_dom_metrics(mat_tip, cov_tip, width_tip)

summary_table_d500 <- tibble::tibble(
  Method = c("Direct Estimator (Sample p_dir)", "fastsae::ebp_area (Spatial Besag)", "tipsae::fit_sae (Spatial Besag)"),
  `Empirical RMSE` = c(mean(dom_dir$Emp_RMSE), mean(dom_fast$Emp_RMSE), mean(dom_tip$Emp_RMSE)),
  `Empirical RRMSE (%)` = c(mean(dom_dir$Emp_RRMSE), mean(dom_fast$Emp_RRMSE), mean(dom_tip$Emp_RRMSE)),
  `Empirical MARB (%)` = c(mean(dom_dir$Emp_MARB), mean(dom_fast$Emp_MARB), mean(dom_tip$Emp_MARB)),
  `Empirical Relative Bias (%)` = c(mean(dom_dir$Emp_RB), mean(dom_fast$Emp_RB), mean(dom_tip$Emp_RB)),
  `Empirical 95% CrI Coverage (%)` = c(NA_real_, mean(dom_fast$Emp_CR), mean(dom_tip$Emp_CR)),
  `Empirical Mean CrI Width` = c(NA_real_, mean(dom_fast$Emp_Width), mean(dom_tip$Emp_Width)),
  `Empirical Efficiency Gain` = c(1.000, mean(dom_dir$Emp_MSE) / mean(dom_fast$Emp_MSE), mean(dom_dir$Emp_MSE) / mean(dom_tip$Emp_MSE)),
  `Avg Runtime (s/fit)` = c(NA_real_, mean(times_fast), mean(times_tip))
)

cat("\n==============================================================================\n")
cat("NATIONAL-SCALE SIMULATION SUMMARY TABLE (D = 500 DOMAINS, R = 25 REPLICATIONS)\n")
cat("==============================================================================\n")
print(as.data.frame(summary_table_d500), digits = 4)
cat(sprintf("\nComputational Speedup at D = 500: FastSAE is %.1fx faster per fit than TipSAE!\n",
            mean(times_tip) / mean(times_fast)))

# ------------------------------------------------------------------------------
# 5. VISUALIZATIONS FOR D = 500
# ------------------------------------------------------------------------------
theme_clean <- theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Plot A: Domain Empirical RMSE (Ordered across 500 domains)
ord_idx <- order(dom_dir$Emp_RMSE)
df_rmse_d500 <- bind_rows(
  data.frame(Rank = 1:D, RMSE = dom_dir$Emp_RMSE[ord_idx], Method = "Direct Estimator"),
  data.frame(Rank = 1:D, RMSE = dom_fast$Emp_RMSE[ord_idx], Method = "FastSAE (INLA)"),
  data.frame(Rank = 1:D, RMSE = dom_tip$Emp_RMSE[ord_idx], Method = "TipSAE (Stan)")
)
df_rmse_d500$Method <- factor(df_rmse_d500$Method, levels = c("Direct Estimator", "FastSAE (INLA)", "TipSAE (Stan)"))

pa <- ggplot(df_rmse_d500, aes(x = Rank, y = RMSE, color = Method)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  scale_color_manual(values = c("Direct Estimator" = "#E64B35", "FastSAE (INLA)" = "#00A087", "TipSAE (Stan)" = "#3C5488")) +
  labs(
    title = "National Scale: Empirical RMSE across D = 500 Domains (R = 25 Replications)",
    subtitle = "Domains ordered by direct survey error; dramatic variance reduction across all 500 areas",
    x = "Small Area Domain Rank (1 to 500)",
    y = "Empirical RMSE vs True P_d"
  ) +
  theme_clean

ggsave(file.path(out_dir, "plot_d500_mc_empirical_rmse.png"), pa, width = 10, height = 4.8, dpi = 300)

# Plot B: Empirical RRMSE Boxplot across 500 domains
df_box_d500 <- bind_rows(
  data.frame(RRMSE = dom_dir$Emp_RRMSE, Method = "Direct Estimator"),
  data.frame(RRMSE = dom_fast$Emp_RRMSE, Method = "FastSAE (INLA)"),
  data.frame(RRMSE = dom_tip$Emp_RRMSE, Method = "TipSAE (Stan)")
)
df_box_d500$Method <- factor(df_box_d500$Method, levels = c("Direct Estimator", "FastSAE (INLA)", "TipSAE (Stan)"))

pb <- ggplot(df_box_d500, aes(x = Method, y = RRMSE, fill = Method)) +
  geom_boxplot(alpha = 0.8, width = 0.5, outlier.size = 1) +
  scale_fill_manual(values = c("Direct Estimator" = "#E64B35", "FastSAE (INLA)" = "#00A087", "TipSAE (Stan)" = "#3C5488")) +
  labs(
    title = "National Scale: Distribution of Empirical RRMSE across D = 500 Domains",
    subtitle = "Relative Root Mean Squared Error (%) relative to True Population Parameter P_d",
    x = "Estimator",
    y = "Empirical RRMSE (%)"
  ) +
  theme_clean +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "plot_d500_mc_rrmse_box.png"), pb, width = 7.5, height = 4.8, dpi = 300)

# Plot C: Agreement Scatter across 500 domains
df_agree_d500 <- data.frame(
  FastSAE = dom_fast$Emp_Mean,
  TipSAE = dom_tip$Emp_Mean,
  True_P = True_P
)

pc <- ggplot(df_agree_d500, aes(x = FastSAE, y = TipSAE)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "firebrick", linewidth = 0.8) +
  geom_point(color = "#00A087", alpha = 0.6, size = 1.6) +
  labs(
    title = sprintf("National Scale Agreement: FastSAE vs TipSAE across D = 500 Domains (r = %.4f)",
                    cor(dom_fast$Emp_Mean, dom_tip$Emp_Mean)),
    subtitle = "Points represent empirical mean estimates across R = 25 replications; line is y = x",
    x = "fastsae::ebp_area Empirical Mean",
    y = "tipsae::fit_sae Empirical Mean"
  ) +
  theme_clean

ggsave(file.path(out_dir, "plot_d500_mc_agreement.png"), pc, width = 6.5, height = 5.5, dpi = 300)

# Save RDS
results_d500 <- list(
  R = R_val,
  D = D,
  summary_table = summary_table_d500,
  dom_dir = dom_dir,
  dom_fast = dom_fast,
  dom_tip = dom_tip,
  True_P = True_P,
  pop_summary = pop_summary,
  mc_elapsed_sec = t_mc_elapsed
)
saveRDS(results_d500, file.path(out_dir, "mc_finite_population_results_D500_R25.rds"))

cat("\n[SUCCESS] National Scale (D = 500) Simulation Results saved to:", normalizePath(out_dir), "\n")
