#!/usr/bin/env Rscript

# ==============================================================================
# Benchmark: fastsae hb_area() (INLA) for the Bayesian benchmark files
#
# Re-measures ONLY the fastsae rows of inst/extdata/{beta,spatial_beta,
# spatio_temporal_beta,hb_normal,hb}_benchmark.rds and keeps the tipsae/hbsae
# rows. Data, model and seeds follow the original generators:
#   beta                 bench_beta_tipsae.R            (Beta, iid)
#   spatial_beta         bench_spatial_beta.R           (Beta, Besag ICAR on a grid)
#   spatio_temporal_beta bench_spatio_temporal_beta.R   (Beta, Besag + RW1, T = 5)
#   hb_normal / hb       bench_hb_normal.R              (Gaussian, iid)
# so the fastsae rows are measured on the same draws as the comparator rows.
#
# Usage: Rscript benchmarks/bench_bayes_fastsae.R [beta|spatial_beta|st_beta|hb ...] [--iter=10]
# ==============================================================================

suppressPackageStartupMessages({
  library(fastsae)
  library(bench)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
iter_arg <- 10L
iter_match <- grep("^--iter=", args, value = TRUE)
if (length(iter_match) > 0) iter_arg <- as.integer(sub("^--iter=", "", iter_match[1]))
which_args <- args[!grepl("^--", args)]
designs <- if (length(which_args) > 0) which_args else c("beta", "spatial_beta", "st_beta", "hb")

T_periods <- 5L

# Regular contiguous grid of n cells (as in the original scripts)
make_grid_sf <- function(n_val) {
  factors <- which(n_val %% 1:n_val == 0)
  nx <- factors[which.min(abs(factors - sqrt(n_val)))]
  ny <- n_val / nx
  if (nx > ny) {
    tmp <- nx; nx <- ny; ny <- tmp
  }
  poly <- sf::st_polygon(list(matrix(c(0, 0, 1, 0, 1, 1, 0, 1, 0, 0), ncol = 2, byrow = TRUE)))
  grid <- sf::st_make_grid(poly, n = c(nx, ny))
  sf::st_sf(domain = as.character(1:n_val), geometry = grid[1:n_val])
}

grid_W <- function(n_val) {
  grid_sf <- make_grid_sf(n_val)
  W <- spdep::nb2mat(spdep::poly2nb(grid_sf), style = "B", zero.policy = TRUE)
  dimnames(W) <- list(grid_sf[["domain"]], grid_sf[["domain"]])
  W
}

# One fastsae call per design; returns the expression to time for domain size n
fit_call <- list(
  beta = function(n_val) {
    df <- sim_area_data(D = n_val, spatial_type = "knn", seed = 42 + n_val)$data
    df$domain <- as.character(df$domain)
    function() hb_area(y_beta ~ x1 + x2, data = df, domain = "domain", vardir = "vardir",
                       family = "beta", print_result = FALSE)
  },
  spatial_beta = function(n_val) {
    W <- grid_W(n_val)
    df <- sim_area_data(D = n_val, W = W, spatial_type = "knn", seed = 42 + n_val)$data
    df$domain <- as.character(df$domain)
    function() hb_area(y_beta ~ x1 + x2, data = df, domain = "domain", vardir = "vardir",
                       family = "beta", spatial = "besag", W = W, print_result = FALSE)
  },
  st_beta = function(n_val) {
    W <- grid_W(n_val)
    df <- sim_series_data(D = n_val, T = T_periods, W = W, n_unsampled = 0,
                          prop_intermittent = 0, seed = 42 + n_val)$data
    df$domain <- as.character(df$area)
    function() hb_area(y_beta ~ x1 + x2, data = df, domain = "domain", time = "year",
                       vardir = "vardir", family = "beta", spatial = "besag",
                       temporal = "rw1", st_interaction = "domain-specific", W = W,
                       print_result = FALSE)
  },
  hb = function(n_val) {
    df <- sim_area_data(D = n_val, n_unsampled = 0, seed = 42 + n_val)$data
    function() hb_area(y_gaussian ~ x1 + x2, vardir = "vardir", data = df,
                       family = "gaussian", print_result = FALSE)
  }
)

files <- list(
  beta = "beta_benchmark.rds",
  spatial_beta = "spatial_beta_benchmark.rds",
  st_beta = "spatio_temporal_beta_benchmark.rds",
  hb = c("hb_normal_benchmark.rds", "hb_benchmark.rds")
)

for (design in designs) {
  paths <- file.path("inst/extdata", files[[design]])
  old <- as.data.frame(readRDS(paths[1]))
  n_values <- sort(unique(old$n))
  cat(sprintf(">>> %s: n = %s (%d iterations)\n", design, paste(n_values, collapse = ", "), iter_arg))

  new_rows <- lapply(n_values, function(n_val) {
    f <- fit_call[[design]](n_val)
    f() # warm-up: loads INLA and the graph code outside the timing
    res <- bench::mark(fastsae = f(), check = FALSE, iterations = iter_arg,
                       filter_gc = FALSE, memory = TRUE, time_unit = "s")
    cat(sprintf("    n = %4d: median %.3f s\n", n_val, as.numeric(res$median)))
    data.frame(
      Method = "fastsae",
      `Min (s)` = as.numeric(res$min),
      `Median (s)` = as.numeric(res$median),
      `Iterations/sec` = as.numeric(res$`itr/sec`),
      `Memory (MB)` = as.numeric(res$mem_alloc) / (1024^2),
      n_iter = as.integer(res$n_itr),
      n = as.numeric(n_val),
      check.names = FALSE
    )
  })

  out <- rbind(old[old$Method != "fastsae", ], do.call(rbind, new_rows))
  out <- out[order(out$n, out$Method != "fastsae"), ]
  rownames(out) <- NULL
  out <- tibble::as_tibble(out)
  class(out) <- c("bench_mark", "tbl_df", "tbl", "data.frame")
  for (p in paths) saveRDS(out, p)
  cat(sprintf("    [Saved] %s\n", paste(paths, collapse = ", ")))
}
