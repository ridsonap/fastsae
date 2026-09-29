#!/usr/bin/env Rscript

# ==============================================================================
# Benchmark: fastsae (hb_area) vs tipsae (fit_sae) for Spatio-Temporal Beta SAE
# Supports incremental appending across domain sizes n
# ==============================================================================

suppressPackageStartupMessages({
  library(fastsae)
  library(tipsae)
  library(sf)
  library(spdep)
  library(bench)
  library(tibble)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

# Parse domain values, iterations, and time periods
n_args <- suppressWarnings(as.integer(args[!grepl("^--", args)]))
n_args <- n_args[!is.na(n_args)]

iter_arg <- 5L
iter_match <- grep("^--iter=", args, value = TRUE)
if (length(iter_match) > 0) {
  iter_arg <- as.integer(sub("^--iter=", "", iter_match[1]))
}

T_arg <- 5L
T_match <- grep("^--T=", args, value = TRUE)
if (length(T_match) > 0) {
  T_arg <- as.integer(sub("^--T=", "", T_match[1]))
}

D_values <- if (length(n_args) > 0) n_args else c(30, 50, 100, 250, 500, 1000)

cat("==============================================================================\n")
cat("Starting Benchmark: fastsae (hb_area) vs tipsae (fit_sae) [Spatio-Temporal Beta]\n")
cat(sprintf("Domain sizes: %s | Time periods (T): %d | Iterations: %d each\n",
            paste(D_values, collapse = ", "), T_arg, iter_arg))
cat("==============================================================================\n\n")

# Setup output directories
out_dirs <- c("inst/exdata", "inst/extdata")
for (d in out_dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

save_all_locations <- function(obj) {
  filenames <- c(
    "spatio_temporal_beta_benchmark.rds"
  )
  for (d in out_dirs) {
    for (f in filenames) {
      saveRDS(obj, file.path(d, f))
    }
  }
}

file_primary <- file.path("inst/exdata", "spatio_temporal_beta_benchmark.rds")

# Load existing benchmark data if present
df_existing <- NULL
if (file.exists(file_primary)) {
  df_existing <- readRDS(file_primary)
  cat(sprintf("Loaded existing benchmark data with %d rows across n = %s\n\n",
              nrow(df_existing), paste(unique(df_existing$n), collapse = ", ")))
}

# Function to generate spatial contiguous grid
make_grid_sf <- function(n_val) {
  factors <- which(n_val %% 1:n_val == 0)
  nx <- factors[which.min(abs(factors - sqrt(n_val)))]
  ny <- n_val / nx
  if (nx > ny) {
    tmp <- nx; nx <- ny; ny <- tmp
  }
  poly <- sf::st_polygon(list(matrix(c(0,0, 1,0, 1,1, 0,1, 0,0), ncol = 2, byrow = TRUE)))
  grid <- sf::st_make_grid(poly, n = c(nx, ny))
  grid_sf <- sf::st_sf(domain = as.character(1:n_val), geometry = grid[1:n_val])
  return(grid_sf)
}

for (n_val in D_values) {
  current_iter <- iter_arg
  cat(sprintf(">>> Running benchmark for n = %d, T = %d (%d iterations each)...\n",
              n_val, T_arg, current_iter))
  
  grid_sf <- make_grid_sf(n_val)
  nb <- spdep::poly2nb(grid_sf)
  W <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
  dimnames(W) <- list(grid_sf[["domain"]], grid_sf[["domain"]])
  
  sim <- sim_series_data(
    D = n_val,
    T = T_arg,
    W = W,
    n_unsampled = 0,
    prop_intermittent = 0,
    seed = 42 + n_val
  )
  df <- sim[["data"]]
  df[["domain"]] <- as.character(df[["area"]])
  
  t_start <- Sys.time()
  
  res <- bench::mark(
    fastsae = fastsae::hb_area(
      formula = y_beta ~ x1 + x2,
      data = df,
      domain = "domain",
      time = "year",
      vardir = "vardir",
      family = "beta",
      spatial = "besag",
      temporal = "rw1",
      st_interaction = "domain-specific",
      W = W,
      print_result = FALSE
    ),
    tipsae = tipsae::fit_sae(
      formula_fixed = y_beta ~ x1 + x2,
      data = df,
      domains = "domain",
      disp_direct = "vardir",
      type_disp = "var",
      likelihood = "beta",
      spatial_error = TRUE,
      spatial_df = grid_sf,
      domains_spatial_df = "domain",
      temporal_error = TRUE,
      temporal_variable = "year",
      refresh = 0,
      seed = 42
    ),
    check = FALSE,
    iterations = current_iter,
    filter_gc = FALSE,
    memory = TRUE,
    time_unit = "s"
  )
  
  t_elapsed <- round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1)
  cat(sprintf("    Finished n = %d in %s seconds.\n", n_val, t_elapsed))
  
  df_n <- tibble::tibble(
    Method = as.character(res$expression),
    `Min (s)` = as.numeric(res$min),
    `Median (s)` = as.numeric(res$median),
    `Iterations/sec` = as.numeric(res$`itr/sec`),
    `Memory (MB)` = as.numeric(res$mem_alloc) / (1024^2),
    n_iter = as.integer(res$n_itr),
    n = as.numeric(n_val)
  )
  
  # Remove existing rows for this n if already present, then append cleanly
  if (!is.null(df_existing)) {
    df_existing_df <- as.data.frame(df_existing)
    df_existing_clean <- df_existing_df[df_existing_df$n != n_val, ]
    df_combined <- rbind(df_existing_clean, as.data.frame(df_n))
  } else {
    df_combined <- as.data.frame(df_n)
  }
  
  # Sort by n
  df_combined <- df_combined[order(df_combined$n), ]
  df_combined <- tibble::as_tibble(df_combined)
  class(df_combined) <- c("bench_mark", "tbl_df", "tbl", "data.frame")
  df_existing <- df_combined
  
  # Save immediately after each n completes
  save_all_locations(df_existing)
  cat(sprintf("    [Saved] Updated spatio_temporal_beta_benchmark.rds in inst/exdata and inst/extdata (Total rows: %d)\n\n",
              nrow(df_existing)))
}

cat("==============================================================================\n")
cat("Spatio-Temporal Beta Benchmark Complete! Final Summary Table:\n")
cat("==============================================================================\n\n")

print(df_existing)

cat("\nSummary of Median Runtime & Speedup:\n")
all_n <- sort(unique(df_existing$n))
for (n_val in all_n) {
  sub <- df_existing[df_existing$n == n_val, ]
  med_fast <- sub$`Median (s)`[sub$Method == "fastsae"]
  med_tip  <- sub$`Median (s)`[sub$Method == "tipsae"]
  mem_fast <- sub$`Memory (MB)`[sub$Method == "fastsae"]
  mem_tip  <- sub$`Memory (MB)`[sub$Method == "tipsae"]
  speedup  <- med_tip / med_fast
  mem_ratio <- mem_tip / mem_fast
  cat(sprintf("  n = %4d: fastsae = %7.4f s (%5.1f MB) | tipsae = %8.4f s (%7.1f MB) | Speedup = %5.2fx | RAM = %5.1fx less\n",
              n_val, med_fast, mem_fast, med_tip, mem_tip, speedup, mem_ratio))
}
