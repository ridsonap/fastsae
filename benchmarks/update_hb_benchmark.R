suppressPackageStartupMessages({
  library(fastsae)
  library(hbsae)
  library(bench)
  library(tibble)
  library(dplyr)
})

cat("==============================================================================\n")
cat("UPDATE BENCHMARK: HAPUS FASTSAEHB & TAMBAHKAN hbsae vs fastsae (n = 500, 1000)\n")
cat("==============================================================================\n\n")

file_primary <- "inst/extdata/hb_normal_benchmark.rds"

if (!file.exists(file_primary)) {
  stop("File benchmark tidak ditemukan: ", file_primary)
}

df_existing <- readRDS(file_primary)
cat(sprintf("Data awal memiliki %d baris:\n", nrow(df_existing)))
print(df_existing)

# 1. Hapus baris fastsaehb dan ubah ke data.frame biasa
df_existing_df <- as.data.frame(df_existing)
df_filtered <- df_existing_df[df_existing_df$Method != "fastsaehb", ]

cat(sprintf("\nSetelah menghapus fastsaehb: %d baris tersisa (n = %s)\n\n",
            nrow(df_filtered), paste(unique(df_filtered$n), collapse = ", ")))

# 2. Jalankan benchmark untuk n = 500 dan n = 1000
new_D <- c(500, 1000)
new_results <- list()

for (n_val in new_D) {
  cat(sprintf(">>> Menjalankan benchmark hbsae vs fastsae untuk n = %d (10 iterasi)...\n", n_val))
  
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
    check = FALSE,
    iterations = 10,
    filter_gc = FALSE,
    memory = TRUE,
    time_unit = "s"
  )
  
  t_elapsed <- round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 2)
  cat(sprintf("    Selesai n = %d dalam %s detik.\n", n_val, t_elapsed))
  
  df_n <- data.frame(
    Method = as.character(res$expression),
    `Min (s)` = as.numeric(res$min),
    `Median (s)` = as.numeric(res$median),
    `Iterations/sec` = as.numeric(res$`itr/sec`),
    `Memory (MB)` = as.numeric(res$mem_alloc) / (1024^2),
    n_iter = as.integer(res$n_itr),
    n = as.numeric(n_val),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  new_results[[as.character(n_val)]] <- df_n
}

df_new <- do.call(rbind, new_results)

# 3. Gabungkan dan urutkan berdasarkan n
df_final_df <- rbind(df_filtered, df_new)
df_final_df <- df_final_df[order(df_final_df$n), ]
df_final <- tibble::as_tibble(df_final_df)
class(df_final) <- c("bench_mark", "tbl_df", "tbl", "data.frame")

# 4. Simpan ke inst/extdata dan inst/extdata
out_dirs <- c("inst/extdata", "inst/extdata")
for (d in out_dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  saveRDS(df_final, file.path(d, "hb_normal_benchmark.rds"))
  saveRDS(df_final, file.path(d, "hb_benchmark.rds"))
  cat(sprintf("[Disimpan] %s (Total: %d baris)\n", file.path(d, "hb_normal_benchmark.rds"), nrow(df_final)))
}

cat("\n==============================================================================\n")
cat("HASIL AKHIR BENCHMARK: hbsae vs fastsae (n = 30, 50, 100, 250, 500, 1000)\n")
cat("==============================================================================\n\n")

print(df_final)

cat("\nRingkasan Perbandingan (Waktu & Memori):\n")
for (n_val in sort(unique(df_final$n))) {
  sub <- df_final[df_final$n == n_val, ]
  t_fast <- sub$`Median (s)`[sub$Method == "fastsae"]
  t_hb   <- sub$`Median (s)`[sub$Method == "hbsae"]
  m_fast <- sub$`Memory (MB)`[sub$Method == "fastsae"]
  m_hb   <- sub$`Memory (MB)`[sub$Method == "hbsae"]
  cat(sprintf("n = %4d | Waktu: fastsae = %6.4f s, hbsae = %6.4f s | Memori: fastsae = %6.2f MB, hbsae = %6.2f MB\n",
              n_val, t_fast, t_hb, m_fast, m_hb))
}
