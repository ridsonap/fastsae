# Ringkasan Eksekusi: Validasi Ilmiah Model `ebp_area` & Diagnostik Bayesian

## 1. Fitur & Peningkatan yang Telah Diimplementasikan

### A. Diagnostik Bayesian Internal pada `diagnose()`
File: [`R/diagnose.R`](file:///Volumes/work/_MainR/fastsae-inla/R/diagnose.R)
Fungsi `diagnose()` kini secara otomatis mengekstrak metrik Bayesian jika model yang dievaluasi merupakan model Bayesian (`fastsae_ebp_area`):
1. **Kriteria Informasi**:
   - **WAIC** (Watanabe-Akaike Information Criterion) beserta parameter efektif $p_{\text{WAIC}}$.
   - **DIC** (Deviance Information Criterion) beserta parameter efektif $p_D$.
   - **Marginal Log-Likelihood** ($\log \text{ML}$) untuk perbandingan Bayes Factor antar struktur spasial.
2. **Uji Kalibrasi Prediktif (Leave-One-Out CPO & PIT)**:
   - **Probability Integral Transform (PIT)**: Menjalankan uji Kolmogorov-Smirnov terhadap distribusi teoritis $\text{Uniform}(0, 1)$ (`stats::ks.test`). Model terkalibrasi baik jika $p \ge 0.05$.
   - **Conditional Predictive Ordinate (CPO)**: Menghitung mean, minimum CPO, serta mendeteksi ada/tidaknya peringatan kegagalan numerik pada observasi pencilan (*outliers*).
   - Menyimpan kolom `cpo` dan `pit` langsung ke dalam tabel `df_diag` pada objek diagnostik.

### B. Visualisasi Diagnostik Bayesian pada `autoplot()`
File: [`R/autoplot.R`](file:///Volumes/work/_MainR/fastsae-inla/R/autoplot.R)
Dua varian visualisasi baru ditambahkan ke `autoplot(diag_obj)`:
- `type = "pit"`: Histogram densitas nilai PIT dengan garis referensi horizontal merah putus-putus pada $\text{density} = 1.0$. Memeriksa kalibrasi sebaran residual (mendeteksi *overdispersion* atau *underdispersion*).
- `type = "cpo"`: Index plot nilai CPO per domain (dengan garis segmen dan titik) untuk mengidentifikasi area yang memiliki kecocokan model rendah (*potential outliers*).

### C. Suite Pengujian Validasi Ilmiah Otomatis
File: [`tests/testthat/test_ebp_validation.R`](file:///Volumes/work/_MainR/fastsae-inla/tests/testthat/test_ebp_validation.R)
Mencakup 3 pengujian validasi komprehensif:
1. **Benchmark Equivalence vs `sae::mseFH`**:
   - Menguji estimasi titik `ebp_area` vs `sae::mseFH` ($r > 0.99$, hasil aktual $r = 0.9974$).
   - Menguji estimasi ketidakpastian / varians ($r > 0.98$, hasil aktual $r = 0.9958$).
   - Rata-rata selisih absolut $< 0.25$ ($0.117$).
2. **Parameter Recovery & Nominal Coverage Rate**:
   - Mensimulasikan data dengan ground truth ($\boldsymbol{\beta} = [10.0, 1.0, 0.5]$).
   - Memverifikasi bahwa seluruh $\beta$ asli masuk ke dalam 95% Credible Interval.
   - Mean Absolute Relative Bias $< 10\%$ (aktual: $1.62\%$).
   - Cakupan interval kredibel 95% memenuhi target nominal $\ge 85\%$ (aktual: $92.5\%$).
3. **Bayesian Diagnostics Verification**:
   - Memverifikasi ekstraksi WAIC, DIC, Marginal LogLik, CPO, dan uji PIT.
   - Memverifikasi pembentukan plot `type = "pit"` dan `type = "cpo"`.

---

## 2. Contoh Tampilan Terminal Diagnostik Bayesian

```text
── fastsae Small Area Estimation Diagnostic Report ─────────────────────────────
Model: "EBP-GAUSSIAN (Non-spatial)" | Domains: 42 (Sampled: 32, Unsampled: 10)

── 1. Precision & Efficiency Gain ──

! Domains with RSE < 25%: 59.5% (Caution: low precision)
✔ Average RSE reduction: Direct 29.49% -> SAE 23.52% (Gain: 5.97%)
✔ Variance reduction in 100% of areas (MSE ratio median: 1.49, max: 5.16)

── 2. Brown et al. (2001) Calibration Tests ──

✔ Bias Regression Test (H0: alpha = 0, beta = 1): F = 2.518, p-value = 0.0975 [Statistically Unbiased]
Estimated parameters: alpha = -0.8053, beta = 1.1734
! Goodness-of-Fit Statistic W (Chi-Square): W = 7.42 (df = 32), p-value = 1 [Deviation from Survey Variance]

── 5. Bayesian Information Criteria & Predictive Diagnostics ──

ℹ WAIC: 118.63 (p_eff: 14.38) | DIC: 119.62 (p_eff: 20.46)
✔ PIT Calibration Test vs Uniform(0,1): D = 0.096, p-value = 0.9035 [Well-calibrated predictive distribution]
✔ Leave-One-Out CPO: No numerical approximation issues (min CPO = 0.0091)
────────────────────────────────────────────────────────────────────────────────
! Final Assessment: CAUTION: Model goodness-of-fit indicates notable deviation from survey variance.
────────────────────────────────────────────────────────────────────────────────
```

---

## 3. Hasil Pengujian Unit

- `test_ebp_validation.R`: **21/21 PASS (100%)**
- `test_diagnose.R`: **26/26 PASS (100%)**
- `test_sim_data.R`: **51/51 PASS (100%)**
- `test_sim_series.R`: **35/35 PASS (100%)**
- `test_ebp_area.R`: **43/43 PASS (100%)**
