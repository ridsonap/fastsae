#!/usr/bin/env Rscript

# ==============================================================================
# One cell of the frequentist benchmark (called by bench_frequentist.R).
#
# Usage: Rscript benchmarks/bench_frequentist_child.R <model> <method> <D> <iter> <out.rds>
#   model  : FH, SFH or STFH
#   method : fastsae, fastsae_mse (STFH only), sae, emdi, or baseline (loads the packages and builds the
#            data without fitting, to measure the resident memory of the process
#            before any fit)
#
# Every method fits the same simulated data set (seed = D) with point estimates and
# analytical MSE, except for STFH (see below). The parent runs this script under /usr/bin/time to record the peak
# resident set size of the process.
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
model <- args[1]
method <- args[2]
D <- as.integer(args[3])
iter <- as.integer(args[4])
out_file <- args[5]

suppressPackageStartupMessages({
  library(fastsae)
  library(sae)
  library(emdi)
  library(bench)
})

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
  data.frame(area = seq_len(D), y = drop(cbind(1, X) %*% beta) + v + rnorm(D, sd = sqrt(vardir)),
             X, vardir = vardir)
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
W <- as.matrix(fastsae:::sim_spatial_weights(D = D, type = "knn", style = "W", k = 4, seed = D))
dat <- switch(model,
  FH = sim_fh(D, seed = D),
  SFH = sim_fh(D, W = W, rho = 0.5, seed = D),
  STFH = sim_stfh(D, W, seed = D)
)

fit <- switch(paste(model, method),
  "FH fastsae" = function() eblup_fh(fml, vardir = "vardir", data = dat, print_result = FALSE),
  "FH sae" = function() sae::mseFH(fml, vardir, method = "REML", data = dat),
  "FH emdi" = function() emdi::fh(fixed = fml, vardir = "vardir", combined_data = dat,
                                  domains = "area", method = "reml", MSE = TRUE),
  "SFH fastsae" = function() eblup_sfh(fml, vardir = "vardir", data = dat, W = W, print_result = FALSE),
  "SFH sae" = function() sae::mseSFH(fml, vardir, W, method = "REML", data = dat),
  "SFH emdi" = function() emdi::fh(fixed = fml, vardir = "vardir", combined_data = dat,
                                   domains = "area", method = "reml", MSE = TRUE,
                                   correlation = "spatial", corMatrix = W),
  # sae has no analytical MSE for the spatio-temporal model (only pbmseSTFH()), so
  # the comparison is of point estimates; "fastsae_mse" adds the analytical MSE
  "STFH fastsae" = function() suppressMessages(eblup_stfh(fml, data = dat, vardir = ~vardir,
                                                          domain = ~area, time = ~year, W = W,
                                                          mse_type = "none",
                                                          print_result = FALSE)),
  "STFH fastsae_mse" = function() suppressMessages(eblup_stfh(fml, data = dat, vardir = ~vardir,
                                                              domain = ~area, time = ~year, W = W,
                                                              mse_type = "analytical",
                                                              print_result = FALSE)),
  "STFH sae" = function() sae::eblupSTFH(fml, D = D, T = T_periods, vardir = vardir, proxmat = W,
                                         data = dat),
  NULL
)

if (method == "baseline") {
  saveRDS(list(model = model, method = method, D = D), out_file)
  quit(save = "no")
}
if (is.null(fit)) stop("no implementation for ", model, " / ", method)

# bench::mark() does not keep the result when check = FALSE; keep the last one
last <- NULL
res <- bench::mark({ last <- fit(); NULL }, iterations = iter, check = FALSE, memory = TRUE)

# Convergence of the fitting algorithm (a non-converged fit runs to its iteration limit)
converged <- function(x) {
  v <- switch(method,
    fastsae = , fastsae_mse = x$convergence,
    sae = if (model == "STFH") x$fit$convergence else x$est$fit$convergence,
    emdi = if (is.data.frame(x$model$variance)) x$model$variance$convergence
  )
  if (is.null(v)) NA else as.logical(v)[1]
}

saveRDS(list(
  model = model, method = method, D = D,
  min = as.numeric(res$min), median = as.numeric(res$median),
  itr_sec = as.numeric(res$`itr/sec`), mem_alloc_mb = as.numeric(res$mem_alloc) / 1024^2,
  n_iter = as.integer(res$n_itr), converged = converged(last)
), out_file)
