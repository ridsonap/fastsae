# Rencana Implementasi: Adaptasi Metode Bawaan fastsae untuk Membaca `ebp` / `eblup`

## 1. Deskripsi Tujuan (Goal Description)

Saat ini, pada tabel hasil `ebp_area()` dan `ebp_spatial_area()`, terdapat kolom `ebp` dan kolom duplikat `eblup` dengan nilai yang sama persis demi kompatibilitas fungsi visualisasi bawaan (`autoplot()`) dan fungsi S3.

Sesuai arahan pengguna (*"rancang saja metode bawaan fastsae agar bisa membaca ebp / eblup"*), kita akan merestrukturisasi metode bawaan `fastsae` agar secara dinamis dan cerdas mampu membaca objek estimasi yang memiliki kolom **`ebp`** maupun **`eblup`**, serta menyesuaikan label grafik dan ringkasan secara kontekstual:
1. **Pembersihan Tabel Output**: Menghapus kolom duplikat `eblup` pada tabel hasil `ebp_area()` sehingga tabel menjadi bersih, ringkas, dan murni menggunakan kolom `ebp`.
2. **Fleksibilitas Metode Bawaan**:
   - `autoplot()`: Mampu memplot model EBLUP (`eblup_fh`, `eblup_sfh`, `eblup_stfh`, `eblup_bhf`) maupun model EBP (`ebp_area`, `ebp_spatial_area`), termasuk perbandingan multi-model antara model EBLUP dan EBP!
   - `print()` & `summary()`: Menampilkan label "EBP" atau "EBLUP" sesuai tipe model dan kolom yang tersedia.
   - `fitted()` & `residuals()`: Otomatis mengekstrak nilai estimasi (`ebp` atau `eblup`) tanpa perlu duplikasi kolom.

---

## 2. User Review Required

> [!NOTE]
> **Penghapusan Kolom Redundan pada Tabel `df_ebp`**:
> Pada fungsi `ebp_area()` dan `ebp_spatial_area()`, kolom `eblup` akan dihapus dari data frame `df_ebp`.
> Tabel `df_ebp` hanya akan memuat:
> `domain`, `y`, `ebp`, `linear_pred`, `sd`, `mse`, `rse`, `ci_lower`, `ci_upper`, `random_effect` (serta kolom opsional `vardir`, `trials`, `exposure`).
> Metode bawaan `autoplot()`, `summary()`, `fitted()`, dan `residuals()` akan otomatis membaca `ebp` jika model bertipe EBP, atau `eblup` jika model bertipe EBLUP.

---

## 3. Rencana Perubahan Detail (Proposed Changes)

### Komponen 1: Pembersihan Tabel Hasil `ebp_area`
#### [MODIFY] [`R/inla_utils.R`](file:///Volumes/work/_MainR/fastsae-inla/R/inla_utils.R)
- Hapus kolom `eblup = ebp_est` dari `df_ebp` di baris 166.
- Tabel hasil murni hanya memuat `ebp`.

#### [MODIFY] [`R/ebp_area.R`](file:///Volumes/work/_MainR/fastsae-inla/R/ebp_area.R)
- Hapus kolom `eblup = preds` dari `.fit_glmm_laplace` di baris 438.

---

### Komponen 2: Adaptasi `autoplot` untuk Mendukung `ebp` dan `eblup`
#### [MODIFY] [`R/autoplot.R`](file:///Volumes/work/_MainR/fastsae-inla/R/autoplot.R)
1. **`.autoplot_single_comparison(x, title = NULL, ...)`**:
   - Deteksi data frame: `df <- x$df_ebp %||% x$df_eblup`.
   - Deteksi kolom estimasi: `est_col <- if ("ebp" %in% names(df)) "ebp" else "eblup"`.
   - Deteksi label: `est_label <- if (est_col == "ebp") "EBP" else "EBLUP"`.
   - Manfaatkan `ci_lower` dan `ci_upper` jika sudah ada (seperti pada EBP INLA), atau hitung dari `mse` jika belum ada.
   - Sumbu Y dan judul grafik otomatis menampilkan `"EBP Estimate"` atau `"EBLUP Estimate"`.
2. **`.autoplot_estimates(x, title = NULL, ...)`**:
   - Deteksi kolom estimasi `df[[est_col]]`.
   - Label judul otomatis: `"EBP Estimates vs Direct Estimates"` atau `"EBLUP Estimates vs Direct Estimates"`.
3. **`.autoplot_multi_comparison(x, title = NULL, ...)`**:
   - Mendukung perbandingan campuran antar model (misal membandingkan model klasik `eblup_fh` dengan model INLA `ebp_area`):
     Masing-masing model diekstrak nilai estimasinya secara dinamis (`m$df_ebp %||% m$df_eblup`).
4. **`.autoplot_multi_mse(x, title = NULL, ...)`**:
   - Mendukung perbandingan MSE/varians baik dari `df_ebp` maupun `df_eblup`.
5. **`.autoplot_multi_scatter(x, title = NULL, ...)`**:
   - Scatter plot perbandingan dua model membaca kolom estimasi dinamis dari model 1 dan model 2.
6. **`.autoplot_mse(x, title = NULL, ...)`**:
   - Membaca `df <- x$df_ebp %||% x$df_eblup`.

---

### Komponen 3: Pengujian & Verifikasi
#### [MODIFY] [`tests/testthat/test_ebp_area.R`](file:///Volumes/work/_MainR/fastsae-inla/tests/testthat/test_ebp_area.R)
- Verifikasi bahwa `df_ebp` tidak lagi memuat kolom duplikat `eblup`.
- Uji coba `autoplot(fit_ebp, type = "estimates")`.
- Uji coba `autoplot(fit_ebp, type = "comparison")`.
- Uji coba perbandingan silang: `autoplot(list("FH-Linear" = fit_fh, "FH-INLA" = fit_ebp), type = "comparison")`.

---

## 4. Verification Plan

### Automated Tests
```bash
# 1. Jalankan unit test spesifik
Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test_ebp_area.R')"

# 2. Jalankan test autoplot yang sudah ada untuk memastikan backward compatibility
Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test_plot.R')"
```

### Manual Verification
1. Uji visualisasi model EBP tunggal:
   ```r
   m <- ebp(y ~ x1 + x2, data = mys, vardir = "vardir", family = "gaussian")
   autoplot(m, type = "estimates")
   autoplot(m, type = "comparison")
   ```
2. Uji visualisasi perbandingan model EBLUP klasik (`eblup_fh`) vs model EBP INLA (`ebp_spatial`):
   ```r
   m_fh <- eblup_fh(y ~ x1 + x2, data = mys, vardir = "vardir")
   m_sp <- ebp_spatial(y ~ x1 + x2, data = mys, vardir = "vardir", W = mys_proxmat)
   autoplot(list("EBLUP FH" = m_fh, "EBP Spatial INLA" = m_sp), type = "comparison")
   ```
