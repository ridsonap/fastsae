# Walkthrough: Adaptasi Metode Bawaan `fastsae` untuk Membaca `ebp` / `eblup` Secara Fleksibel

Kita telah menyelesaikan restrukturisasi metode bawaan [`fastsae`](file:///Volumes/work/_MainR/fastsae-inla) (`autoplot`, `print`, `summary`, `fitted`, `residuals`) agar secara cerdas dan fleksibel membaca kolom **`ebp`** maupun **`eblup`**, serta membersihkan tabel hasil `ebp_area` dari kolom duplikat.

---

## 1. Perubahan yang Dilakukan (Changes Made)

### A. Pembersihan Tabel Hasil `ebp_area`
* **[`R/inla_utils.R`](file:///Volumes/work/_MainR/fastsae-inla/R/inla_utils.R)**:
  - Menghapus kolom duplikat `eblup = ebp_est` dari pembentukan `df_ebp`.
  - Tabel hasil `df_ebp` kini murni dan bersih hanya memuat kolom:
    `domain`, `y`, `ebp`, `linear_pred`, `sd`, `mse`, `rse`, `ci_lower`, `ci_upper`, `random_effect` (serta kolom opsional `vardir`, `trials`, `exposure`).
* **[`R/ebp_area.R`](file:///Volumes/work/_MainR/fastsae-inla/R/ebp_area.R)**:
  - Menghapus kolom duplikat `eblup = preds` pada fitting `.fit_glmm_laplace`.

### B. Adaptasi Fleksibel pada `autoplot()`
* **[`R/autoplot.R`](file:///Volumes/work/_MainR/fastsae-inla/R/autoplot.R)**:
  - **Single Model Plot (`comparison`)**:
    Secara dinamis mendeteksi keberadaan kolom `ebp` atau `eblup`. Jika data memiliki `ci_lower` dan `ci_upper` (seperti pada EBP INLA), langsung dimanfaatkan sebagai pita ketidakpastian. Sumbu Y dan judul otomatis beradaptasi menjadi `"EBP Estimate"` atau `"EBLUP Estimate"`.
  - **Estimates Plot (`estimates`)**:
    Menampilkan scatter plot nilai langsung $y$ vs estimasi (`ebp` atau `eblup`) dengan judul adaptif.
  - **Multi-Model Comparison (`comparison`)**:
    Mendukung perbandingan campuran antar model (misal membandingkan model klasik `eblup_fh` dengan model INLA `ebp_area`), mengekstrak estimasi masing-masing secara independen.
  - **Multi-Model Scatter (`scatter`)**:
    Mengekstrak estimasi model 1 dan model 2 secara dinamis dengan label sumbu yang akurat.
  - **MSE Plot (`mse`)**:
    Mendukung visualisasi MSE baik untuk objek EBP maupun EBLUP.

### C. Pengujian Komprehensif
* **[`tests/testthat/test_ebp_area.R`](file:///Volumes/work/_MainR/fastsae-inla/tests/testthat/test_ebp_area.R)**:
  - Menambahkan verifikasi eksplisit bahwa kolom `eblup` tidak ada di `df_ebp` (`expect_false("eblup" %in% names(df_ebp))`).
  - Menambahkan pengujian `autoplot(type = "estimates")`, `autoplot(type = "comparison")`, `autoplot(type = "mse")`, serta perbandingan multi-model `autoplot(list(m0, m1))`.

---

## 2. Hasil Verifikasi & Testing (Validation Results)

1. **Unit Test EBP (`test_ebp_area.R`)**:
   Seluruh **43 pengujian** berhasil lulus 100%:
   ```
   ══ Testing test_ebp_area.R ═════════════════════════════════════════════════════
   [ FAIL 0 | WARN 0 | SKIP 0 | PASS 43 ] Done!
   ```
2. **Backward Compatibility Test (`test_plot.R`)**:
   Seluruh **20 pengujian** autoplot model lama (`eblup_fh`, `eblup_sfh`, `eblup_bhf`) tetap berjalan dan lulus 100%:
   ```
   ══ Testing test_plot.R ═════════════════════════════════════════════════════════
   [ FAIL 0 | WARN 0 | SKIP 0 | PASS 20 ] Done!
   ```
3. **BHF Test (`test_eblup_bhf.R`)**:
   Lulus 100% tanpa kendala.

---

## 3. Contoh Penggunaan Visualisasi

```r
library(devtools)
load_all()

# 1. Model Non-Spasial INLA
m0 <- ebp(y ~ x1 + x2 + x3, data = mys, vardir = "vardir", family = "gaussian")

# Tabel hasil bersih:
head(m0$df_ebp)
# Hanya ada kolom 'ebp' dan 'linear_pred', tidak ada lagi kolom duplikat 'eblup'!

# 2. Model Spasial BYM2 INLA
m_bym2 <- ebp_spatial(y ~ x1 + x2 + x3, data = mys, W = mys_proxmat, vardir = "vardir", spatial = "bym2")

# 3. Model Klasik Fay-Herriot (EBLUP)
m_fh <- eblup_fh(y ~ x1 + x2 + x3, data = mys, vardir = "vardir")

# 4. Visualisasi Tunggal (otomatis berlabel "EBP Estimate")
autoplot(m_bym2, type = "estimates")
autoplot(m_bym2, type = "comparison")

# 5. Visualisasi Perbandingan Campuran (EBLUP klasik vs EBP INLA)
autoplot(list("EBLUP Klasik" = m_fh, "EBP Spatial INLA" = m_bym2), type = "comparison")
```
