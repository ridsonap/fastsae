# Ringkasan Implementasi: Fitur Uji Diagnostik Universal (`diagnose`)

## 1. Fitur yang Telah Diimplementasikan

### A. Fungsi Evaluasi Diagnostik Universal (`diagnose`)
File: [`R/diagnose.R`](file:///Volumes/work/_MainR/fastsae-inla/R/diagnose.R)
Fungsi `diagnose()` bekerja secara universal untuk semua model di `fastsae` (`ebp_area`, `eblup_fh`, `eblup_sfh`, `eblup_stfh`). Fungsi ini mengevaluasi 4 pilar kualitas model SAE:

1. **Evaluasi Presisi & Peningkatan Efisiensi**:
   - Persentase domain yang memenuhi batas toleransi RSE (default $< 25\%$).
   - Perbandingan rata-rata dan median RSE Direct vs RSE Model.
   - Rasio Penurunan Varians (*Efficiency Gain* $\text{vardir} / \text{MSE}$) (Min, Q1, Median, Mean, Q3, Max) serta persentase area yang mengalami reduksi varians.

2. **Uji Kalibrasi & Bias Eksternal (Brown et al., 2001)**:
   - **Bias Regression**: Mengestimasi regresi OLS $y_d^{dir} = \alpha + \beta \hat{\theta}_d^{SAE}$ dan menguji hipotesis simultan $H_0: \alpha = 0, \beta = 1$ via Wald $F$-test.
   - **Goodness-of-Fit Chi-Square ($W$-statistic)**: Menguji apakah deviasi terstandar sesuai dengan varians sampling ($W = \sum \frac{(y_d - \hat{\theta}_d)^2}{\text{vardir}_d + \text{MSE}_d} \sim \chi^2(D_s)$).

3. **Uji Sisa Autokorelasi Spasial (Residual Moran's $I$)**:
   - Menghitung indeks Moran's $I$, nilai ekspektasi, varians teoritis, $z$-score, dan $p$-value dua arah pada residual model.
   - Otomatis mendeteksi matriks spasial $W$ dari objek model spasial (`eblup_sfh` atau `ebp_area` spasial).

4. **Metrik Evaluasi Simulasi (*Ground Truth*)**:
   - Jika argumen `truth` diberikan (misal dari data simulasi `sim_area_data` atau `sim_series_data`), otomatis menghitung:
     - Mean Relative Bias (MRB)
     - Mean Absolute Relative Bias (MARB)
     - Mean Relative RMSE (RRMSE)
     - Coverage Rate (CR) interval kredibel 95%.

### B. Visualisasi Diagnostik Terpadu via `autoplot(diag_obj)`
File: [`R/autoplot.R`](file:///Volumes/work/_MainR/fastsae-inla/R/autoplot.R)
Mendukung 5 jenis visualisasi diagnostik berbasis `ggplot2`:
* `type = "all"`: Overview grafis komposit yang memadukan titik kalibrasi 1:1 dan gradien warna presisi RSE SAE.
* `type = "calibration"`: Scatter plot $y^{dir}$ vs $\hat{\theta}^{SAE}$ dengan garis identitas $y=x$ (merah putus-putus) dan garis regresi (hitam).
* `type = "rse"`: Scatter plot RSE Model vs RSE Direct membuktikan penurunan varians (titik di bawah garis diagonal).
* `type = "residuals"`: Scatter plot residual terstandar vs nilai dugaan fitted values.
* `type = "qq"`: Normal Q-Q plot dari residual terstandar untuk memeriksa asumsi normalitas.

---

## 2. Contoh Penggunaan Cepat

```r
library(fastsae)
data(mys)
data(mys_proxmat)

# 1. Fit Model EBLUP Spasial (atau ebp_area)
fit_sfh <- eblup_sfh(y ~ x1 + x2, data = mys, vardir = ~vardir, domain = ~area, W = mys_proxmat)

# 2. Jalankan Diagnostik Lengkap
diag_res <- diagnose(fit_sfh)
print(diag_res)

# 3. Buat Grafik Diagnostik
autoplot(diag_res, type = "calibration")
autoplot(diag_res, type = "rse")
autoplot(diag_res, type = "residuals")
autoplot(diag_res, type = "qq")
```

---

## 3. Tampilan Modern Berbasis Terminal (`cli`)

Metode cetak `print.fastsae_diagnose()` menggunakan paket `cli` dengan indikator visual modern:
- `cli::cli_alert_success` (`✔`) untuk uji yang lolos (unbiased, precision gain, tidak ada autokorelasi sisa).
- `cli::cli_alert_warning` (`!`) untuk peringatan deviasi moderat atau presisi rendah.
- `cli::cli_alert_danger` (`✖`) untuk uji yang gagal secara signifikan.

```text
── fastsae Small Area Estimation Diagnostic Report ─────────────────────────────
Model: "FH" | Domains: 42 (Sampled: 32, Unsampled: 10)

── 1. Precision & Efficiency Gain ──

! Domains with RSE < 25%: 52.4% (Caution: low precision)
✔ Average RSE reduction: Direct 29.49% -> SAE 26.27% (Gain: 3.22%)
✔ Variance reduction in 100% of areas (MSE ratio median: 1.29, max: 3.91)

── 2. Brown et al. (2001) Calibration Tests ──

✔ Bias Regression Test (H0: alpha = 0, beta = 1): F = 3.144, p-value = 0.0576 [Statistically Unbiased]
Estimated parameters: alpha = -0.7317, beta = 1.1557
! Goodness-of-Fit Statistic W (Chi-Square): W = 4.53 (df = 32), p-value = 1 [Deviation from Survey Variance]
────────────────────────────────────────────────────────────────────────────────
! Final Assessment: CAUTION: Model goodness-of-fit indicates notable deviation from survey variance.
────────────────────────────────────────────────────────────────────────────────
```

---

## 4. Hasil Pengujian Unit

Seluruh 26 pengujian pada [`tests/testthat/test_diagnose.R`](file:///Volumes/work/_MainR/fastsae-inla/tests/testthat/test_diagnose.R) **lulus 100%**:
- Validasi diagnostik pada model Fay-Herriot (`eblup_fh`).
- Validasi diagnostik dan uji Moran's $I$ pada model Spatial Fay-Herriot (`eblup_sfh`).
- Validasi diagnostik dan metrik ground truth pada Bayesian EBP (`ebp_area`).
- Validasi pembuatan grafik `ggplot2` untuk seluruh 5 variasi plot `autoplot(diag)`.
