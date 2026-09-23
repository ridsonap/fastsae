# Rencana Metodologis & Implementasi: Uji Kesesuaian & Diagnostik Estimasi SAE

Dokumen ini memaparkan **kerangka kerja komprehensif (metodologis & teknis)** untuk menguji apakah hasil estimasi Small Area Estimation (EBP / EBLUP) sudah tepat, valid, dan reliabel, serta merancang fitur diagnostik otomatis di paket [`fastsae`](file:///Volumes/work/_MainR/fastsae-inla).

---

## 1. Kerangka Metodologi Evaluasi Hasil Estimasi SAE

Dalam literatur Small Area Estimation (Rao & Molina 2015; Pfeffermann 2013; Brown et al. 2001; Tzavidis et al. 2018), pengujian kualitas estimasi model-based terbagi menjadi **4 Pilar Utama**:

```
                       ┌──────────────────────────────────────────────┐
                       │        EVALUASI KESESUAIAN ESTIMASI          │
                       │           SMALL AREA ESTIMATION              │
                       └──────────────────────┬───────────────────────┘
                                              │
         ┌──────────────────┬─────────────────┴────────────────┬──────────────────┐
         ▼                  ▼                                  ▼                  ▼
┌──────────────────┐ ┌───────────────┐               ┌──────────────────┐ ┌────────────────┐
│  Pilar 1: Presisi│ │Pilar 2: Uji   │               │ Pilar 3: Asumsi  │ │Pilar 4: Evalu- │
│  & Efisiensi     │ │Kalibrasi Brown│               │ Model & Residual │ │asi Simulasi    │
│                  │ │(Bias & Fit)   │               │                  │ │(Jika Ada Truth)│
│- RSE / CV < 25%  │ │- Regresi H0:  │               │- Residual Norm.  │ │- Relative Bias │
│- Gain in MSE     │ │  α = 0, β = 1 │               │- Ranef Normality │ │- RRMSE         │
│- Rata-rata MSE   │ │- Statistik W  │               │- Moran's I resid │ │- Coverage Rate │
│  turun signifikan│ │  Chi-Square   │               │- DIC / WAIC / AIC│ │  (95% CI)      │
└──────────────────┘ └───────────────┘               └──────────────────┘ └────────────────┘
```

---

### Pilar 1: Evaluasi Presisi & Peningkatan Efisiensi (*Gain in Efficiency*)

Tujuan utama SAE adalah meminjam kekuatan (*borrowing strength*) untuk menekan varians sampling yang besar pada penduga langsung (*direct estimator*).

1. **Relative Standard Error (RSE / Coefficient of Variation)**:
   $$\text{RSE}(\hat{\theta}_d) = \frac{\sqrt{\text{MSE}(\hat{\theta}_d)}}{\hat{\theta}_d} \times 100\%$$
   *Standar Acuan (BPS, Eurostat, US Census Bureau)*:
   - **RSE < 25%** (atau < 20%): Estimasi presisi tinggi, sangat reliabel untuk rilis data resmi.
   - **25% $\le$ RSE < 50%**: Estimasi cukup reliabel, perlu diberi catatan kehati-hatian (*use with caution*).
   - **RSE $\ge$ 50%**: Estimasi tidak reliabel untuk interpretasi individual area.

2. **Rasio Efisiensi (*Efficiency Gain*)**:
   Membandingkan varians penduga langsung ($\text{vardir}_d$) dengan MSE dari EBP/EBLUP:
   $$\text{Gain}_d = \frac{\text{vardir}_d - \text{MSE}_d^{\text{EBP}}}{\text{vardir}_d} \times 100\% \quad \text{atau} \quad \text{Ratio}_d = \frac{\text{vardir}_d}{\text{MSE}_d^{\text{EBP}}}$$
   - **Kriteria Kesesuaian**: Model dinyatakan berhasil jika $\text{MSE}_d^{\text{EBP}} < \text{vardir}_d$ (Ratio > 1) pada mayoritas area tersampel. Nilai efisiensi gain biasanya berkisar antara **30% hingga 80%**.

---

### Pilar 2: Uji Validasi & Kalibrasi Eksternal (*Brown et al., 2001*)

Uji diagnostik baku yang dirumuskan oleh Brown, Chambers, Heady, & Heasman (2001) untuk menguji apakah estimasi model konsisten dan tidak menyimpang (*unbiased*) terhadap estimasi langsung:

1. **Bias Diagnostic (Regresi Wald Test)**:
   Melakukan regresi antara penduga langsung $y_d^{dir}$ terhadap penduga model $\hat{\theta}_d^{EBP}$:
   $$y_d^{dir} = \alpha + \beta \hat{\theta}_d^{EBP} + \varepsilon_d$$
   - **Hipotesis**: $H_0: \alpha = 0 \text{ dan } \beta = 1$.
   - **Uji**: Wald $F$-test simultan.
   - **Kriteria**: Jika $p$-value $> 0.05$, $H_0$ **gagal ditolak**, membuktikan bahwa estimasi EBP tidak bias secara sistematis (*statistically unbiased*).

2. **Goodness-of-Fit Diagnostic (Statistik $W$ Chi-Square)**:
   Menguji apakah selisih kuadrat terstandar antara estimasi langsung dan model konsisten dengan varians sampling:
   $$W = \sum_{d=1}^{D_s} \frac{(y_d^{dir} - \hat{\theta}_d^{EBP})^2}{\text{vardir}_d + \text{MSE}_d^{EBP}} \sim \chi^2(D_s)$$
   - **Kriteria**: Jika nilai $W$ berada di dalam selang kritis $\left[\chi^2_{0.025, D_s}, \, \chi^2_{0.975, D_s}\right]$ (atau $p$-value $> 0.05$), model dinyatakan cocok (*good fit*) dengan data.

3. **Plot Kalibrasi (Direct vs EBP Scatter Plot)**:
   - Titik-titik $(\hat{\theta}_d^{EBP}, y_d^{dir})$ harus menyebar simetris di sekitar garis diagonal $y = x$.
   - Tidak boleh ada deviasi lengkung (*non-linearity*) atau penyimpangan kemiringan (*slope deflection*).

---

### Pilar 3: Asumsi Model & Diagnostik Residual

1. **Uji Normalitas & Homoskedastisitas**:
   - **Q-Q Plot Residual**: Memeriksa apakah residual terstandar berdistribusi normal $N(0, 1)$.
   - **Q-Q Plot Random Effects**: Memeriksa apakah efek acak domain $\hat{u}_d$ berdistribusi normal $N(0, \sigma_u^2)$.
   - **Residuals vs Fitted**: Memeriksa apakah residual bebas dari pola heteroskedastisitas.

2. **Uji Sisa Autokorelasi Spasial (*Residual Moran's I*)**:
   - Menghitung indeks Moran's $I$ pada residual model dengan matriks ketetanggaan $W$.
   - **Kriteria**: Nilai $p$-value Moran's $I$ harus $> 0.05$. Ini membuktikan bahwa seluruh dependensi spasial telah berhasil diserap oleh komponen spasial model (misal BYM2 atau Besag), sehingga residualnya bersifat *white noise* independen.

3. **Kriteria Informasi Model (Model Selection)**:
   - **DIC & WAIC** (Bayesian INLA): Semakin kecil nilai DIC dan WAIC, semakin baik keseimbangan antara kecocokan model (*fit*) dan kompleksitas parameter (*parsimony*).
   - Penurunan $\Delta \text{DIC} > 5$ mengindikasikan model spasial secara nyata lebih unggul daripada non-spasial.

---

### Pilar 4: Validasi Berbasis Data Simulasi (*Ground Truth*)

Khusus ketika menggunakan data simulasi (seperti dari fungsi [`sim_area_data`](file:///Volumes/work/_MainR/fastsae-inla/R/sim_area_data.R) atau [`sim_series_data`](file:///Volumes/work/_MainR/fastsae-inla/R/sim_series_data.R)), nilai parameter sejati $\theta_d$ diketahui secara pasti:
1. **Relative Bias (RB)**:
   $$\text{RB}_d = \frac{\hat{\theta}_d^{EBP} - \theta_d}{\theta_d} \times 100\% \quad \Rightarrow \quad \overline{\text{ARB}} = \frac{1}{D} \sum_{d=1}^D |\text{RB}_d| < 5\%$$
2. **Relative Root Mean Squared Error (RRMSE)**:
   $$\text{RRMSE}_d = \frac{\sqrt{(\hat{\theta}_d^{EBP} - \theta_d)^2}}{\theta_d} \times 100\%$$
3. **Coverage Rate (CR)**:
   Persentase domain di mana interval kredibel $95\%$ $(\text{CI}_{lower}, \text{CI}_{upper})$ memuat nilai parameter sejati $\theta_d$ (target ideal: $\approx 95\%$).

---

## 2. Rencana Implementasi: Fitur Uji Diagnostik Otomatis (`sae_diagnostic`)

Untuk memudahkan pengguna menjalankan seluruh pengujian di atas tanpa perlu menghitung rumus manual satu per satu, kita dapat mengimplementasikan fungsi baru:

### Nama Fungsi: `sae_diagnostic()`
File: `R/diagnostic.R`

```r
sae_diagnostic(
  object,              # Objek hasil estimasi fastsae (ebp_area, eblup_fh, dll.)
  W = NULL,            # Matriks spasial untuk uji Moran's I residual (opsional)
  truth = NULL,        # Vektor nilai sejati jika data simulasi (opsional)
  rse_threshold = 25   # Ambang batas RSE yang dianggap reliabel (default 25%)
)
```

### Nilai Kembalian (`fastsae_diagnostic`):
Objek list yang memuat:
* **`$precision`**: Persentase domain dengan RSE < 25%, ringkasan efisiensi gain (Min, Median, Mean, Max penurunan varians).
* **`$brown_test`**:
  - `wald_test`: Nilai $F$-statistic, p-value uji hipotesis $\alpha = 0, \beta = 1$.
  - `goodness_of_fit`: Nilai statistik $W$, derajat bebas $D_s$, nilai kritis $\chi^2$, dan $p$-value.
* **`$spatial_test`**: Nilai Moran's $I$ pada residual dan p-value (jika $W$ tersedia).
* **`$simulation_metrics`**: Nilai Mean Relative Bias (MRB), Mean RRMSE, dan Coverage Rate (jika `truth` diberikan).
* **`$status`**: Kesimpulan status evaluasi ringkas (*"Model Pas & Reliabel"*, *"Model Mengandung Bias"*, atau *"Presisi Rendah"*).

### Visualisasi Diagnostik Terpadu via `autoplot(diag_obj)`:
Menghasilkan grafik 4-panel elegan menggunakan `ggplot2`:
1. **Panel 1: Direct vs EBP Calibration** (Scatter plot dengan garis identitas $y = x$ dan garis regresi estimasi).
2. **Panel 2: RSE Comparison** (Scatter plot RSE EBP vs RSE Direct membuktikan penurunan varians).
3. **Panel 3: Residuals vs Fitted** (Pengecekan heteroskedastisitas).
4. **Panel 4: Residual Normal Q-Q Plot** (Pengecekan normalitas).

---

## 3. Contoh Tampilan Output di Konsol

```
========================================================================
             fastsae Small Area Estimation Diagnostic Report           
========================================================================
Model: EBP-GAUSSIAN (BYM2) | Domains: 42 (Sampled: 32, Unsampled: 10)

1. PRECISION & EFFICIENCY GAIN:
  - Domains with RSE < 25%:     32 / 32 (100.0%) [PASS]
  - Average Direct RSE:         17.42%
  - Average EBP RSE:            8.15%  (Average Gain: 53.2%)
  - MSE Ratio (vardir / MSE):   Min: 1.45 | Median: 2.15 | Max: 4.80

2. BROWN ET AL. (2001) CALIBRATION TESTS:
  - Bias Diagnostic (H0: alpha = 0, beta = 1):
    F-statistic: 0.412, p-value: 0.6657 [PASS - No Systemic Bias]
  - Goodness-of-Fit Diagnostic (Chi-Square W):
    W-statistic: 28.64 (df = 32), p-value: 0.6371 [PASS - Good Fit]

3. RESIDUAL SPATIAL AUTOCORRELATION:
  - Residual Moran's I: 0.021, p-value: 0.384 [PASS - No Residual Spatial Autocorrelation]

========================================================================
OVERALL ASSESSMENT: MODEL VALID & RELIABLE FOR PUBLICATION
========================================================================
```

---

## 4. User Review Required

> [!IMPORTANT]
> **Pilihan Fitur Uji Diagnostik**:
> Apakah Anda ingin kita langsung menambahkan fungsi `sae_diagnostic()` dan `autoplot.fastsae_diagnostic()` ke dalam paket [`fastsae`](file:///Volumes/work/_MainR/fastsae-inla)?
> 
> Fungsi ini akan otomatis dapat dipakai untuk semua model:
> - Model EBP baru: `fit <- ebp_area(...)` $\rightarrow$ `sae_diagnostic(fit)`
> - Model EBLUP yang sudah ada: `fit <- eblup_fh(...)` / `eblup_sfh(...)` $\rightarrow$ `sae_diagnostic(fit)`
