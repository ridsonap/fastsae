#!/usr/bin/env Rscript

# ==============================================================================
# Design-based Monte Carlo simulation: D = 500 domains, N ~ 250,000 units.
#
# A fixed finite population (binary outcome, logistic in two covariates with ICAR +
# iid area effects) is sampled R times by stratified simple random sampling (20-50
# units per domain). Each sample is fitted by
#   - fastsae::eblup_fh()   Fay-Herriot EBLUP, analytical MSE
#   - fastsae::eblup_sfh()  spatial Fay-Herriot EBLUP (row-standardised W), analytical MSE
#   - fastsae::hb_area()    beta model with an ICAR effect and vardir (defaults)
#   - tipsae::fit_sae()     beta model with an ICAR effect, type_disp = "var"
#                           (Stan, 2 chains of 1,000 iterations after 500 warm-up)
# and evaluated against the true domain proportions: RMSE, RRMSE, MARB, coverage of
# the 95% intervals (EBLUP +- 1.96 sqrt(mse) for the EBLUPs, credible intervals for
# the HB fits) with its Monte Carlo standard error across replicates, interval width,
# and the ratio of the mean estimated MSE (posterior variance for HB) to the empirical
# MSE, by domain.
#
# Usage (from the package root):
#   Rscript benchmarks/sim_mc_finite_pop_d500.R [--R=200] [--cores=5]
# Writes benchmarks/output/mc_design_based_D500_R<R>.rds.
# ==============================================================================

suppressPackageStartupMessages({
  library(fastsae)
  library(tipsae)
  library(sf)
  library(spdep)
  library(parallel)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
arg_int <- function(name, default) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m)) as.integer(sub(paste0("^--", name, "="), "", m[1])) else default
}
R_val <- arg_int("R", 200L)
# each tipsae fit runs its 2 chains in parallel, so 5 workers use 10 cores
cores_val <- arg_int("cores", min(5L, max(1L, parallel::detectCores() %/% 2L)))
tip_chains <- 2L
tip_iter <- 1500L
tip_warmup <- 500L

out_dir <- "benchmarks/output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
cat(sprintf("Design-based simulation: D = 500, R = %d, %d workers\n", R_val, cores_val))

# ------------------------------------------------------------------------------
# 1. Fixed finite population (D = 500, N ~ 250,000)
# ------------------------------------------------------------------------------
set.seed(999L)
D <- 500L
domain_ids <- sprintf("area_%03d", seq_len(D))

N_d <- sample(400:600, size = D, replace = TRUE)
names(N_d) <- domain_ids

# 20 x 25 grid of contiguous squares
poly <- sf::st_polygon(list(matrix(c(0, 0, 1, 0, 1, 1, 0, 1, 0, 0), ncol = 2, byrow = TRUE)))
grid <- sf::st_make_grid(poly, n = c(20, 25))
grid_sf <- sf::st_sf(domain = domain_ids, geometry = grid[1:D])
nb <- spdep::poly2nb(grid_sf)
W <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
dimnames(W) <- list(domain_ids, domain_ids)
W_std <- W / rowSums(W)

# ICAR + iid area effects
Q <- diag(rowSums(W)) - W
eig <- eigen(Q)
pos_idx <- which(eig$values > 1e-10)
v_raw <- as.vector(eig$vectors[, pos_idx] %*% (stats::rnorm(length(pos_idx)) / sqrt(eig$values[pos_idx])))
v_sp <- (v_raw - mean(v_raw)) / stats::sd(v_raw) * 0.30
u_iid <- stats::rnorm(D, mean = 0, sd = 0.25)
u_total <- u_iid + v_sp
names(u_total) <- domain_ids

beta_0 <- -1.40
beta_1 <- 0.55
beta_2 <- -0.45

pop_list <- vector("list", D)
for (d in seq_len(D)) {
  dom_id <- domain_ids[d]
  n_units <- N_d[d]
  x1_ind <- stats::rnorm(n_units, mean = 1.6 + 0.2 * sin(d / 10), sd = 0.6)
  x2_ind <- stats::runif(n_units, min = 0, max = 2.2)
  eta_ind <- beta_0 + beta_1 * x1_ind + beta_2 * x2_ind + u_total[dom_id]
  y_ind <- stats::rbinom(n_units, size = 1, prob = stats::plogis(eta_ind))
  pop_list[[d]] <- data.frame(domain = dom_id, x1 = x1_ind, x2 = x2_ind, y = y_ind,
                              stringsAsFactors = FALSE)
}
pop_df <- bind_rows(pop_list)

pop_summary <- pop_df %>%
  group_by(domain) %>%
  summarise(N_d = n(), True_P_d = mean(y), x1_pop_mean = mean(x1), x2_pop_mean = mean(x2),
            .groups = "drop") %>%
  as.data.frame()
rownames(pop_summary) <- pop_summary$domain
True_P <- pop_summary$True_P_d
names(True_P) <- pop_summary$domain
cat(sprintf("Population: N = %d, true proportions %.3f to %.3f\n", sum(N_d), min(True_P), max(True_P)))

estimators <- c("Direct", "FH", "SFH", "fastsae", "tipsae")

# ------------------------------------------------------------------------------
# 2. One replicate
# ------------------------------------------------------------------------------
run_one_rep <- function(r) {
  rep_seed <- 8000L + r
  set.seed(rep_seed)

  n_d_vec <- sample(20:50, size = D, replace = TRUE)
  names(n_d_vec) <- domain_ids
  smp_df <- bind_rows(lapply(domain_ids, function(dom_id) {
    dom_units <- which(pop_df$domain == dom_id)
    pop_df[sample(dom_units, size = n_d_vec[dom_id], replace = FALSE), ]
  }))

  # Direct proportions with a continuity adjustment (keeps them in (0, 1), so that
  # the sampling variances are positive) and their SRS variances with the fpc
  smp <- smp_df %>%
    group_by(domain) %>%
    summarise(n_smp = n(), y_raw = mean(y), .groups = "drop") %>%
    left_join(pop_summary[, c("domain", "N_d", "x1_pop_mean", "x2_pop_mean")], by = "domain") %>%
    mutate(
      p_dir = (y_raw * (n_smp - 1) + 0.5) / n_smp,
      vardir = (1 - n_smp / N_d) * p_dir * (1 - p_dir) / (n_smp - 1)
    ) %>%
    as.data.frame()
  rownames(smp) <- smp$domain
  smp <- smp[domain_ids, ]
  fml <- p_dir ~ x1_pop_mean + x2_pop_mean

  out <- list(r = r, est = list(), low = list(), up = list(), unc = list(), time = list())
  put <- function(name, est, low, up, unc, time) {
    out$est[[name]] <<- est
    out$low[[name]] <<- low
    out$up[[name]] <<- up
    out$unc[[name]] <<- unc
    out$time[[name]] <<- time
  }
  put("Direct", smp$p_dir, smp$p_dir - 1.96 * sqrt(smp$vardir), smp$p_dir + 1.96 * sqrt(smp$vardir),
      smp$vardir, NA_real_)

  # A. Fay-Herriot and spatial Fay-Herriot EBLUPs with analytical MSE
  eblup_put <- function(name, call) {
    t0 <- Sys.time()
    fit <- tryCatch(suppressMessages(suppressWarnings(call())), error = function(e) NULL)
    tm <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (is.null(fit)) return(put(name, rep(NA_real_, D), rep(NA_real_, D), rep(NA_real_, D), rep(NA_real_, D), tm))
    df <- fit$df_eblup[match(domain_ids, fit$df_eblup$domain), ]
    se <- sqrt(df$mse)
    put(name, df$eblup, df$eblup - 1.96 * se, df$eblup + 1.96 * se, df$mse, tm)
  }
  eblup_put("FH", function() eblup_fh(fml, vardir = "vardir", domain = "domain", data = smp,
                                      print_result = FALSE))
  eblup_put("SFH", function() eblup_sfh(fml, vardir = "vardir", domain = "domain", data = smp,
                                        W = W_std, print_result = FALSE))

  # B. fastsae::hb_area, beta model with an ICAR effect
  t0 <- Sys.time()
  fit_fast <- fastsae::hb_area(formula = fml, data = smp, domain = "domain", vardir = "vardir",
                               family = "beta", spatial = "besag", W = W, print_result = FALSE)
  tm <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  df <- fit_fast$df_hb[match(domain_ids, fit_fast$df_hb$domain), ]
  put("fastsae", df$hb, df$ci_lower, df$ci_upper, df$sd^2, tm)

  # C. tipsae::fit_sae, beta model with an ICAR effect
  t0 <- Sys.time()
  fit_tip <- tipsae::fit_sae(
    formula_fixed = fml, data = smp, domains = "domain", disp_direct = "vardir",
    type_disp = "var", likelihood = "beta", spatial_error = TRUE, spatial_df = grid_sf,
    domains_spatial_df = "domain", temporal_error = FALSE,
    chains = tip_chains, cores = tip_chains, iter = tip_iter, warmup = tip_warmup,
    refresh = 0, seed = rep_seed
  )
  tm <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  df_tip <- tipsae::extract(summary(fit_tip))$in_sample
  df_tip <- df_tip[match(domain_ids, df_tip$Domains), ]
  put("tipsae", df_tip$`HB est.`, df_tip$`2.5%`, df_tip$`97.5%`, df_tip$sd^2, tm)
  s <- rstan::summary(fit_tip$stanfit)$summary
  out$rhat_tip <- max(s[, "Rhat"], na.rm = TRUE)
  out$rhat_tip_theta <- max(s[grep("theta[", rownames(s), fixed = TRUE), "Rhat"], na.rm = TRUE)
  out
}

# ------------------------------------------------------------------------------
# 3. Run the replicates
# ------------------------------------------------------------------------------
t_start <- Sys.time()
mc_results <- parallel::mclapply(seq_len(R_val), function(r) {
  tryCatch(run_one_rep(r), error = function(e) structure(conditionMessage(e), class = "try-error"))
}, mc.cores = cores_val, mc.preschedule = FALSE)
mc_elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
ok <- !vapply(mc_results, inherits, logical(1), "try-error")
cat(sprintf("%d of %d replicates completed in %.1f min\n", sum(ok), R_val, mc_elapsed / 60))
if (any(!ok)) print(unique(unlist(mc_results[!ok])))
mc_results <- mc_results[ok]
R_ok <- length(mc_results)

# ------------------------------------------------------------------------------
# 4. Performance measures
# ------------------------------------------------------------------------------
mat <- function(field, name) vapply(mc_results, function(x) x[[field]][[name]], numeric(D))

summ <- lapply(estimators, function(m) {
  est <- mat("est", m)
  err <- est - True_P
  covered <- (True_P >= mat("low", m)) & (True_P <= mat("up", m))
  width <- mat("up", m) - mat("low", m)
  unc <- mat("unc", m)
  mse_d <- rowMeans(err^2, na.rm = TRUE)
  cov_rep <- colMeans(covered, na.rm = TRUE)
  dom <- data.frame(
    domain = domain_ids, True_P = True_P,
    Emp_Mean = rowMeans(est, na.rm = TRUE),
    Emp_RB = (rowMeans(est, na.rm = TRUE) - True_P) / True_P * 100,
    Emp_MARB = rowMeans(abs(err) / True_P, na.rm = TRUE) * 100,
    Emp_MSE = mse_d, Emp_RMSE = sqrt(mse_d), Emp_RRMSE = sqrt(mse_d) / True_P * 100,
    Emp_CR = rowMeans(covered, na.rm = TRUE) * 100,
    Emp_Width = rowMeans(width, na.rm = TRUE),
    MSE_ratio = rowMeans(unc, na.rm = TRUE) / mse_d
  )
  list(dom = dom, summary = data.frame(
    Estimator = m,
    RMSE = mean(dom$Emp_RMSE), RRMSE = mean(dom$Emp_RRMSE), MARB = mean(dom$Emp_MARB),
    RB = mean(dom$Emp_RB),
    Coverage = mean(cov_rep) * 100, Coverage_mcse = stats::sd(cov_rep) / sqrt(length(cov_rep)) * 100,
    Width = mean(dom$Emp_Width),
    MSE_ratio = stats::median(dom$MSE_ratio),
    Efficiency = NA_real_,
    Time = mean(vapply(mc_results, function(x) x$time[[m]], numeric(1)), na.rm = TRUE),
    Failed = sum(colSums(is.na(est)) > 0)
  ))
})
names(summ) <- estimators
summary_table <- do.call(rbind, lapply(summ, `[[`, "summary"))
mse_direct <- mean(summ$Direct$dom$Emp_MSE)
summary_table$Efficiency <- mse_direct / vapply(summ, function(s) mean(s$dom$Emp_MSE), numeric(1))
rownames(summary_table) <- NULL
rhat <- vapply(mc_results, `[[`, numeric(1), "rhat_tip")
rhat_theta <- vapply(mc_results, `[[`, numeric(1), "rhat_tip_theta")

print(summary_table, digits = 4)
cat(sprintf("tipsae: largest R-hat above 1.01 in %d and above 1.05 in %d of %d fits\n",
            sum(rhat > 1.01), sum(rhat > 1.05), R_ok))

saveRDS(list(
  R = R_ok, R_planned = R_val, D = D, N = sum(N_d),
  summary_table = summary_table,
  dom = lapply(summ, `[[`, "dom"),
  rhat_tip = rhat, rhat_tip_theta = rhat_theta,
  tipsae_settings = list(chains = tip_chains, iter = tip_iter, warmup = tip_warmup),
  True_P = True_P, pop_summary = pop_summary,
  mc_elapsed_sec = mc_elapsed, workers = cores_val,
  versions = vapply(c("fastsae", "tipsae", "INLA"), function(p) as.character(packageVersion(p)), "")
), file.path(out_dir, sprintf("mc_design_based_D500_R%d.rds", R_val)))
