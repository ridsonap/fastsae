skip_if_not_installed("emdi")
skip_if_not_installed("MASS")

suppressMessages({
  library(emdi)
  library(MASS) # untuk mvrnorm jika diperlukan
})

# helper ------------------------------------------------------------------
quiet <- function(expr) {
  result <- NULL
  # Tangkap output ke memori (bukan file/nullfile()) -- beberapa lingkungan
  # sandbox memblokir penulisan ke /dev/null (nullfile() di Unix), jadi ini
  # lebih portabel sekaligus tetap cross-platform.
  invisible(capture.output(result <- suppressWarnings(suppressMessages(expr))))
  result
}


# ------------------------- 1. Desain simulasi -------------------------------
# Populasi "aman" (lihat helper-simulate_safe_sfh.R): W dari kontiguitas
# grid row-standardized & rasio sinyal/noise dijaga supaya rho tidak mentok
# batas +-1 di hampir semua replikasi (beda dari desain lama berbasis
# mys/mys_proxmat, yang mentok batas di ~19 dari 20 replikasi). rho_true
# dipertahankan di 0.7 (sudah divalidasi: 0/100 mentok batas, 0/100 gagal
# konvergen dengan desain ini).
pop <- simulate_safe_sfh(rho_true = 0.7)
m <- pop$m
W <- pop$W
X <- pop$X
beta_true <- pop$beta_true
sigma2_true <- pop$sigma2_true
rho_true <- pop$rho_true
vardir <- pop$vardir
u_true <- pop$u_true
theta_true <- pop$theta_true

# ------------------------- 2. Pengaturan simulasi ---------------------------
nsim <- 100
B <- 100

# Struktur penyimpanan hasil: array [replikasi, area]
eblup_fast <- eblup_emdi <- matrix(NA_real_, nsim, m)
mse_fast <- mse_emdi <- matrix(NA_real_, nsim, m)


# ------------------------- 3. Loop simulasi Monte Carlo ---------------------
for (s in seq_len(nsim)) {
  # message(sprintf("Replikasi %d / %d", s, nsim))
  e_s <- rnorm(m, mean = 0, sd = sqrt(vardir))
  y_s <- theta_true + e_s

  dat_s <- data.frame(
    y      = y_s,
    x1     = X[, "x1"],
    x2     = X[, "x2"],
    x3     = X[, "x3"],
    vardir = vardir
  )


  seed_s <- 1000 + s

  set.seed(seed_s)
  fit_fast_s <- tryCatch(
    quiet(
      eblup_sfh(
        y ~ x1 + x2 + x3,
        vardir = "vardir",
        method = "REML",
        data = dat_s,
        W = W,
        mse_method = "pbmse",
        print_result = FALSE,
        B = B
      )
    ),
    error = function(e) NULL
  )

  set.seed(seed_s)
  fit_emdi_s <- tryCatch(
    quiet(
      emdi::fh(
        y ~ x1 + x2 + x3,
        vardir = "vardir",
        method = "reml",
        combined_data = dat_s,
        correlation = "spatial",
        corMatrix = W,
        MSE = TRUE,
        mse_type = "spatialparbootbc",
        B = B
      )
    ),
    error = function(e) NULL
  )

  # simpan hasil
  if (!is.null(fit_fast_s) && !is.null(fit_emdi_s)) {
    eblup_fast[s, ] <- fit_fast_s$df_eblup$eblup
    mse_fast[s, ]   <- fit_fast_s$df_eblup$mse_pbbc
    eblup_emdi[s, ] <- fit_emdi_s$ind$FH
    mse_emdi[s, ]   <- fit_emdi_s$MSE$FH
  }
}


# ------------------------- 4. True MSE empiris per area ---------------------
true_mse <- function(eblup_mat, theta_true) {
  colMeans((sweep(eblup_mat, 2, theta_true, "-"))^2, na.rm = TRUE)
}

true_mse_fast <- true_mse(eblup_fast, theta_true)
true_mse_emdi <- true_mse(eblup_emdi, theta_true)

# pertama
test_that("Unbias Monte Carlo", {
  expect_equal(
    true_mse_fast,
    true_mse_emdi,
    tolerance = 1e-6
  )
})

test_that("Parametric Bootstrap MSE agrees with emdi::fh", {
  mse_domain_sae  <- colMeans(mse_emdi, na.rm = TRUE)
  mse_domain_fast <- colMeans(mse_fast, na.rm = TRUE)

  # 1. Selisih relatif rata-rata MSE per area harus kecil (< 5%)
  expect_equal(
    mse_domain_sae,
    mse_domain_fast,
    tolerance = 0.05
  )

  # 2. Pola/ranking MSE antar area harus konsisten (m = 20-50 area -> ambang 0.9 bermakna)
  cor_mse <- cor(mse_domain_sae, mse_domain_fast, method = "spearman")
  expect_gt(cor_mse, 0.9, label = paste("Spearman correlation of MSE:", cor_mse))

  # 3. Safety net longgar: guard kalau tolerance di atas suatu saat dilonggarkan,
  #    order of magnitude tetap tidak boleh menyimpang jauh.
  #    (Catatan: dengan tolerance 5% di atas, cek ini praktis selalu lolos
  #    kalau expect_equal juga lolos -- dipertahankan sebagai dokumentasi niat,
  #    bukan lapisan proteksi independen.)
  ratio <- mean(mse_domain_sae) / mean(mse_domain_fast)
  expect_gt(ratio, 0.5, label = paste("Ratio of mean MSE:", ratio))
  expect_lt(ratio, 2.0, label = paste("Ratio of mean MSE:", ratio))
})



