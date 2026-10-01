#!/usr/bin/env Rscript

# ==============================================================================
# Benchmark: fastsae (hb_area) vs hbsae (fSAE.Area) vs fastsaehb (hb_normal)
# Fay-Herriot Normal / Gaussian Hierarchical Bayes Model
# Across Domain Sizes n = 30, 50, 100, 250
# ==============================================================================

suppressPackageStartupMessages({
  library(fastsae)
  library(hbsae)
  library(fastsaehb)
  library(bench)
  library(tibble)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

# Parse domain values and iterations
n_args <- suppressWarnings(as.integer(args[!grepl("^--", args)]))
n_args <- n_args[!is.na(n_args)]

iter_arg <- 10L
iter_match <- grep("^--iter=", args, value = TRUE)
if (length(iter_match) > 0) {
  iter_arg <- as.integer(sub("^--iter=", "", iter_match[1]))
}

D_values <- if (length(n_args) > 0) n_args else c(30, 50, 100, 250)

cat("==============================================================================\n")
cat("BENCHMARK HIERARCHICAL BAYES NORMAL (FAY-HERRIOT AREA LEVEL)\n")
cat("Metode yang diuji: fastsae (hb_area) vs hbsae (fSAE.Area) vs fastsaehb (hb_normal)\n")
cat(sprintf("Domain sizes (n): %s | Iterations: %d each\n", paste(D_values, collapse = ", "), iter_arg))
cat("==============================================================================\n\n")

out_dirs <- c("inst/extdata", "inst/extdata")
for (d in out_dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

file_primary <- file.path("inst/extdata", "hb_normal_benchmark.rds")

df_results <- list()

for (n_val in D_values) {
  current_iter <- iter_arg
  cat(sprintf(">>> Running benchmark for n = %d (%d iterations each)...\n", n_val, current_iter))
  
  sim <- fastsae::sim_area_data(D = n_val, n_unsampled = 0, seed = 42 + n_val)
  df <- sim$data
  X <- stats::model.matrix(~ x1 + x2, data = df)
  
  t_start <- Sys.time()
  
  res <- bench::mark(
    fastsae = fastsae::hb_area(
      formula = y_gaussian ~ x1 + x2,
      vardir = "vardir",
      data = df,
      family = "gaussian",
      print_result = FALSE
    ),
    hbsae = hbsae::fSAE.Area(
      est.init = df$y_gaussian,
      var.init = df$vardir,
      X = X,
      method = "HB"
    ),
    fastsaehb = fastsaehb::hb_normal(
      formula = y_gaussian ~ x1 + x2,
      vardir = "vardir",
      data = df,
      n.iter = 5000,
      n.burnin = 1000
    ),
    check = FALSE,
    iterations = current_iter,
    filter_gc = FALSE,
    memory = TRUE,
    time_unit = "s"
  )
  
  t_elapsed <- round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 2)
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
  
  df_results[[as.character(n_val)]] <- df_n
}

df_all <- do.call(rbind, df_results)
df_all <- tibble::as_tibble(df_all)
class(df_all) <- c("bench_mark", "tbl_df", "tbl", "data.frame")

# Save RDS in target directories
for (d in out_dirs) {
  saveRDS(df_all, file.path(d, "hb_normal_benchmark.rds"))
  saveRDS(df_all, file.path(d, "hb_benchmark.rds"))
  cat(sprintf("    [Saved] %s (Total rows: %d)\n", file.path(d, "hb_normal_benchmark.rds"), nrow(df_all)))
}

cat("\n==============================================================================\n")
cat("BENCHMARK SELESAI! TABEL HASIL LENGKAP:\n")
cat("==============================================================================\n\n")

print(df_all)

cat("\nRingkasan Metrik (Waktu Median & Alokasi RAM):\n")
for (n_val in D_values) {
  sub <- df_all[df_all$n == n_val, ]
  cat(sprintf("\n--- Domain n = %d ---\n", n_val))
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("  %-10s | Median: %8.4f s | Memori: %8.2f MB | Itr/sec: %7.2f\n",
                sub$Method[i], sub$`Median (s)`[i], sub$`Memory (MB)`[i], sub$`Iterations/sec`[i]))
  }
}
