#!/usr/bin/env Rscript

# ==============================================================================
# Benchmark: fastsae frequentist EBLUP + analytical MSE (FH, SFH, STFH)
#
# Re-measures ONLY the fastsae rows of inst/extdata/{eblup,seblup,steblup}_benchmark.rds
# and keeps the existing sae/emdi rows. The original generator of those files is
# not available, so the data below follow the documented design (D domains,
# 5 auxiliary covariates, T = 5 periods for STFH) but are not the same draws as
# the ones used for the sae/emdi rows.
#
# Usage: Rscript benchmarks/bench_frequentist_fastsae.R [D ...] [--iter=20]
# ==============================================================================

suppressPackageStartupMessages({
  library(fastsae)
  library(bench)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
n_args <- suppressWarnings(as.integer(args[!grepl("^--", args)]))
n_args <- n_args[!is.na(n_args)]
iter_arg <- 20L
iter_match <- grep("^--iter=", args, value = TRUE)
if (length(iter_match) > 0) iter_arg <- as.integer(sub("^--iter=", "", iter_match[1]))
D_values <- if (length(n_args) > 0) n_args else c(30, 50, 100, 250, 500, 1000)

T_periods <- 5L
p_cov <- 5L

# Area-level data: y = X beta + v + e, v ~ SAR(rho) (rho = 0 for FH)
sim_fh <- function(D, W = NULL, rho = 0, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(D * p_cov), D, p_cov, dimnames = list(NULL, paste0("x", 1:p_cov)))
  beta <- c(1, seq(0.8, -0.8, length.out = p_cov))
  u <- rnorm(D, sd = 0.7)
  v <- if (is.null(W)) u else solve(diag(D) - rho * W, u)
  vardir <- runif(D, 0.3, 1.2)
  data.frame(y = drop(cbind(1, X) %*% beta) + v + rnorm(D, sd = sqrt(vardir)), X, vardir = vardir)
}

# Spatio-temporal data (domain-major): v1 ~ SAR(rho1) over domains, v2 ~ AR(1)(rho2) within domain
sim_stfh <- function(D, W, rho1 = 0.5, rho2 = 0.6, seed = 1) {
  set.seed(seed)
  M <- D * T_periods
  X <- matrix(rnorm(M * p_cov), M, p_cov, dimnames = list(NULL, paste0("x", 1:p_cov)))
  beta <- c(1, seq(0.8, -0.8, length.out = p_cov))
  v1 <- solve(diag(D) - rho1 * W, rnorm(D, sd = 0.6))
  v2 <- as.vector(sapply(seq_len(D), function(d) {
    e <- numeric(T_periods)
    e[1] <- rnorm(1, sd = 0.4 / sqrt(1 - rho2^2))
    for (t in 2:T_periods) e[t] <- rho2 * e[t - 1] + rnorm(1, sd = 0.4)
    e
  }))
  vardir <- runif(M, 0.3, 1.2)
  data.frame(
    area = rep(seq_len(D), each = T_periods), year = rep(seq_len(T_periods), D),
    y = drop(cbind(1, X) %*% beta) + rep(v1, each = T_periods) + v2 + rnorm(M, sd = sqrt(vardir)),
    X, vardir = vardir
  )
}

fml <- as.formula(paste("y ~", paste0("x", 1:p_cov, collapse = " + ")))

to_rows <- function(res, n_val) {
  tibble::tibble(
    Method = "fastsae",
    `Min (s)` = as.numeric(res$min),
    `Median (s)` = as.numeric(res$median),
    `Iterations/sec` = as.numeric(res$`itr/sec`),
    `Memory (MB)` = as.numeric(res$mem_alloc) / (1024^2),
    n_iter = as.integer(res$n_itr),
    n = as.numeric(n_val)
  )
}

update_file <- function(file, new_rows) {
  old <- as.data.frame(readRDS(file))
  keep <- old[!(old$Method == "fastsae" & old$n %in% new_rows$n), ]
  out <- rbind(keep, as.data.frame(new_rows))
  out <- out[order(out$n, match(out$Method, c("fastsae", "emdi", "sae"))), ]
  out <- tibble::as_tibble(out)
  class(out) <- c("bench_mark", "tbl_df", "tbl", "data.frame")
  saveRDS(out, file)
  invisible(out)
}

rows <- list(FH = list(), SFH = list(), STFH = list())
for (n_val in D_values) {
  cat(sprintf(">>> D = %d (%d iterations)\n", n_val, iter_arg))
  W <- as.matrix(sim_spatial_weights(D = n_val, type = "knn", style = "W", k = 4, seed = n_val))
  d_fh <- sim_fh(n_val, seed = n_val)
  d_sfh <- sim_fh(n_val, W = W, rho = 0.5, seed = n_val)
  d_st <- sim_stfh(n_val, W, seed = n_val)

  r_fh <- bench::mark(
    eblup_fh(fml, vardir = "vardir", data = d_fh, print_result = FALSE),
    iterations = iter_arg, check = FALSE, memory = TRUE
  )
  r_sfh <- bench::mark(
    eblup_sfh(fml, vardir = "vardir", data = d_sfh, W = W, print_result = FALSE),
    iterations = iter_arg, check = FALSE, memory = TRUE
  )
  r_st <- bench::mark(
    suppressMessages(eblup_stfh(fml, data = d_st, vardir = ~vardir, domain = ~area, time = ~year,
                                W = W, compute_mse = TRUE, print_result = FALSE)),
    iterations = iter_arg, check = FALSE, memory = TRUE
  )
  rows$FH[[length(rows$FH) + 1]] <- to_rows(r_fh, n_val)
  rows$SFH[[length(rows$SFH) + 1]] <- to_rows(r_sfh, n_val)
  rows$STFH[[length(rows$STFH) + 1]] <- to_rows(r_st, n_val)
  cat(sprintf("    FH %.4f s | SFH %.4f s | STFH %.4f s (medians)\n",
              as.numeric(r_fh$median), as.numeric(r_sfh$median), as.numeric(r_st$median)))
}

files <- c(FH = "eblup_benchmark.rds", SFH = "seblup_benchmark.rds", STFH = "steblup_benchmark.rds")
for (m in names(files)) {
  out <- update_file(file.path("inst/extdata", files[[m]]), do.call(rbind, rows[[m]]))
  cat("\n==", files[[m]], "\n")
  print(as.data.frame(out)[, c("Method", "Median (s)", "Memory (MB)", "n_iter", "n")], digits = 4)
}
