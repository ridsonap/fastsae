#!/usr/bin/env Rscript

# ==============================================================================
# Benchmark: fastsae (ebp_area) vs tipsae (fit_sae) for Beta SAE
# Supports incremental appending across domain sizes n
# ==============================================================================

suppressPackageStartupMessages({
  library(fastsae)
  library(tipsae)
  library(bench)
  library(tibble)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

# Parse domain values and iterations
n_args <- suppressWarnings(as.integer(args[!grepl("^--", args)]))
n_args <- n_args[!is.na(n_args)]

iter_arg <- 20L
iter_match <- grep("^--iter=", args, value = TRUE)
if (length(iter_match) > 0) {
  iter_arg <- as.integer(sub("^--iter=", "", iter_match[1]))
}

D_values <- if (length(n_args) > 0) n_args else c(250, 500)

cat("==============================================================================\n")
cat("Starting Benchmark: fastsae (ebp_area) vs tipsae (fit_sae) [Beta Family]\n")
cat(sprintf("Domain sizes to run: %s | Iterations: %d each\n", paste(D_values, collapse = ", "), iter_arg))
cat("==============================================================================\n\n")

out_dir <- "inst/exdata"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_primary <- file.path(out_dir, "beta_benchmark.rds")
file_alias   <- file.path(out_dir, "benchmark_beta.rds")

# Load existing benchmark data if present
df_existing <- NULL
if (file.exists(file_primary)) {
  df_existing <- readRDS(file_primary)
  cat(sprintf("Loaded existing benchmark data with %d rows across n = %s\n\n",
              nrow(df_existing), paste(unique(df_existing$n), collapse = ", ")))
} else if (file.exists(file_alias)) {
  df_existing <- readRDS(file_alias)
}

for (n_val in D_values) {
  current_iter <- iter_arg
  cat(sprintf(">>> Running benchmark for n = %d (%d iterations each)...\n", n_val, current_iter))
  
  sim <- sim_area_data(D = n_val, spatial_type = "knn", seed = 42 + n_val)
  df <- sim$data
  df$domain <- as.character(df$domain)
  
  t_start <- Sys.time()
  
  res <- bench::mark(
    fastsae = fastsae::ebp_area(
      formula = y_beta ~ x1 + x2,
      data = df,
      domain = "domain",
      vardir = "vardir",
      family = "beta",
      print_result = FALSE
    ),
    tipsae = tipsae::fit_sae(
      formula_fixed = y_beta ~ x1 + x2,
      data = df,
      domains = "domain",
      disp_direct = "vardir",
      type_disp = "var",
      likelihood = "beta",
      spatial_error = FALSE,
      temporal_error = FALSE,
      refresh = 0,
      seed = 42
    ),
    check = FALSE,
    iterations = current_iter,
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
  
  df_combined <- tibble::as_tibble(df_combined)
  class(df_combined) <- c("bench_mark", "tbl_df", "tbl", "data.frame")
  df_existing <- df_combined
  
  # Save immediately to both filenames after each n completes
  saveRDS(df_existing, file_primary)
  saveRDS(df_existing, file_alias)
  cat(sprintf("    [Saved] Updated %s and %s (Total rows: %d)\n\n", file_primary, file_alias, nrow(df_existing)))
}

cat("==============================================================================\n")
cat("Benchmark Complete! Final Summary Table:\n")
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
