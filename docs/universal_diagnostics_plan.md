# Rencana Arsitektur & Metodologi: Uji Diagnostik Universal untuk Seluruh Model `fastsae`

Dokumen ini menjawab pertanyaan: **"Apakah kerangka pengujian diagnostik ini berlaku untuk `eblup_fh`, `eblup_sfh`, `eblup_stfh`, dan model lainnya juga?"**

> [!IMPORTANT]
> **Jawaban Tegas: YA, 100%!**
> Kerangka kerja evaluasi presisi, efisiensi gain, uji kalibrasi bias Brown et al. (2001), goodness-of-fit $\chi^2$, uji residual Moran's $I$, dan plot diagnostik adalah **metodologi universal dalam literatur Small Area Estimation** (Rao & Molina 2015).
> Metrik-metrik ini berlaku sama persis untuk model Frequentist (`eblup_fh`, `eblup_sfh`, `eblup_stfh`) maupun Bayesian EBP (`ebp_area`).

---

## 1. Matriks Keselarasan Diagnostik Antar Model di `fastsae`

Tabel berikut merangkum uji diagnostik yang berlaku untuk masing-masing fungsi model di [`fastsae`](file:///Volumes/work/_MainR/fastsae-inla):

| Uji Diagnostik / Evaluasi | `eblup_fh()` | `eblup_sfh()` | `eblup_stfh()` | `ebp_area()` |
| :--- | :---: | :---: | :---: | :---: |
| **Pilar 1: Presisi (RSE < 25%)** | ✅ Ya | ✅ Ya | ✅ Ya (jika MSE dihitung) | ✅ Ya |
| **Pilar 1: Efficiency Gain (`vardir / MSE`)** | ✅ Ya | ✅ Ya | ✅ Ya | ✅ Ya (Gaussian FH) |
| **Pilar 2: Brown Bias Regression ($H_0: \alpha=0, \beta=1$)** | ✅ Ya | ✅ Ya | ✅ Ya | ✅ Ya |
| **Pilar 2: Brown Goodness-of-Fit ($\chi^2$ test $W$)** | ✅ Ya | ✅ Ya | ✅ Ya | ✅ Ya |
| **Pilar 2: Calibration Plot ($y^{dir}$ vs $\hat{\theta}^{SAE}$)** | ✅ Ya | ✅ Ya | ✅ Ya | ✅ Ya |
| **Pilar 3: Q-Q Plot Residual & Normalitas** | ✅ Ya | ✅ Ya | ✅ Ya | ✅ Ya (Gaussian/Laplace) |
| **Pilar 3: Uji Residual Moran's $I$ (Spasial)** | ❌ (Non-spasial) | ✅ Ya (Otomatis via $W$) | ✅ Ya (via $W$) | ✅ Ya (jika ada $W$) |
| **Pilar 4: Evaluasi Simulasi (RB, RRMSE, CR)** | ✅ Ya | ✅ Ya | ✅ Ya | ✅ Ya |

---

## 2. Mengapa Metodologi Ini Universal?

1. **Inti dari Small Area Estimation adalah Mengungguli *Direct Estimator***:
   - Apapun modelnya (apakah Fay-Herriot standar, Spatial FH, Spatio-Temporal FH, atau Bayesian INLA), tujuan praktis statistikawan adalah membuktikan:
     a. **Varians turun**: $\text{MSE}(\hat{\theta}_d) < \text{Var}(y_d^{dir})$.
     b. **Estimator tidak bias**: Garis regresi $y_d^{dir} = \alpha + \beta \hat{\theta}_d$ harus memiliki $\alpha = 0$ dan $\beta = 1$.
     c. **Estimator andal**: Persentase wilayah dengan $\text{RSE} < 25\%$ meningkat drastis.

2. **Konsistensi Struktur Objek di `fastsae`**:
   Semua fungsi model di paket `fastsae` mengembalikan objek S3 berkelas `fastsae` yang memuat tabel hasil yang kompatibel:
   - Pada `eblup_fh`, `eblup_sfh`, `eblup_stfh`: kolom tersimpan di `$df_eblup` (`eblup`, `y`, `vardir`, `mse`, `rse`, `random_effect`).
   - Pada `ebp_area`: kolom tersimpan di `$df_ebp` (`ebp`, `y`, `vardir`, `mse`, `rse`, `random_effect`).

---

## 3. Desain API Polimorfik: `sae_diagnostic()`

Fungsi `sae_diagnostic()` dirancang bersifat **polimorfik** (otomatis mendeteksi tipe model):

```r
sae_diagnostic(
  object,              # Objek model: eblup_fh, eblup_sfh, eblup_stfh, atau ebp_area
  W = NULL,            # Matriks spasial (otomatis diambil jika ada di model spasial)
  truth = NULL,        # Vektor nilai sejati jika data simulasi
  rse_threshold = 25,  # Ambang batas RSE reliabel (default 25%)
  alpha_level = 0.05   # Tingkat signifikansi uji hipotesis
)
```

### Alur Kerja Deteksi Otomatis:
```
sae_diagnostic(object)
│
├── 1. Deteksi Kolom Estimasi & Data:
│   - Jika ada df_ebp: ambil object$df_ebp$ebp
│   - Jika ada df_eblup: ambil object$df_eblup$eblup
│   - Ambil y_direct = df$y dan vardir = df$vardir
│
├── 2. Ekstraksi Matriks Spasial W:
│   - Jika W diberikan secara eksplisit: gunakan W
│   - Jika W = NULL: periksa apakah object$W ada (pada eblup_sfh / ebp_area spasial)
│
├── 3. Eksekusi Uji Diagnostik:
│   ├── (a) RSE summary (< 25%, Mean, Median) & Efficiency Gain Ratio
│   ├── (b) Uji Regresi Brown et al. (H0: alpha = 0, beta = 1) via Wald Test
│   ├── (c) Uji Goodness-of-Fit Brown (Statistik W Chi-Square)
│   ├── (d) Uji Moran's I pada residual (jika W tersedia)
│   └── (e) Evaluasi Truth (Relative Bias & RRMSE) jika truth != NULL
│
└── 4. Kembalikan Objek S3 `fastsae_diagnostic`
    - Method `print()`: Ringkasan rapi di konsol dengan indikator [PASS / CAUTION / FAIL]
    - Method `autoplot()`: Grafik diagnostik 4-panel elegan ggplot2
```

---

## 4. Contoh Penggunaan Lintas Model

### Contoh 1: Diagnostik pada Model Standar Fay-Herriot (`eblup_fh`)
```r
library(fastsae)
data(mys)

# Fit Fay-Herriot
fit_fh <- eblup_fh(y ~ x1 + x2, data = mys, vardir = ~vardir, domain = ~area)

# Jalankan diagnostik otomatis
diag_fh <- sae_diagnostic(fit_fh)
print(diag_fh)
autoplot(diag_fh)
```

### Contoh 2: Diagnostik pada Model Spatial Fay-Herriot (`eblup_sfh`)
```r
data(mys_proxmat)

# Fit Spatial Fay-Herriot
fit_sfh <- eblup_sfh(y ~ x1 + x2, data = mys, vardir = ~vardir, domain = ~area, W = mys_proxmat)

# Otomatis menjalankan uji Moran's I residual karena W tersimpan di dalam model
diag_sfh <- sae_diagnostic(fit_sfh)
print(diag_sfh)
autoplot(diag_sfh)
```

### Contoh 3: Diagnostik pada Model Bayesian INLA (`ebp_area`)
```r
# Fit Bayesian EBP
fit_ebp <- ebp_area(y ~ x1 + x2, data = mys, vardir = "vardir", spatial = "bym2", W = mys_proxmat)

# Diagnostik model Bayesian
diag_ebp <- sae_diagnostic(fit_ebp)
print(diag_ebp)
autoplot(diag_ebp)
```

---

## 5. Rencana Perubahan Berkas (Proposed Changes)

1. **[NEW] [`R/sae_diagnostic.R`](file:///Volumes/work/_MainR/fastsae-inla/R/sae_diagnostic.R)**:
   - Implementasi `sae_diagnostic()`
   - Method `print.fastsae_diagnostic()`
   - Perhitungan uji Wald $\alpha=0, \beta=1$, statistik Chi-Square Brown, dan Moran's $I$.
2. **[MODIFY] [`R/autoplot.R`](file:///Volumes/work/_MainR/fastsae-inla/R/autoplot.R)**:
   - Implementasi `autoplot.fastsae_diagnostic()` (4-panel diagnostic plot).
3. **[MODIFY] [`DESCRIPTION`](file:///Volumes/work/_MainR/fastsae-inla/DESCRIPTION)** & **[`NAMESPACE`](file:///Volumes/work/_MainR/fastsae-inla/NAMESPACE)**:
   - Ekspor fungsi `sae_diagnostic`.
4. **[NEW] [`tests/testthat/test_diagnostic.R`](file:///Volumes/work/_MainR/fastsae-inla/tests/testthat/test_diagnostic.R)**:
   - Unit test untuk `eblup_fh`, `eblup_sfh`, dan `ebp_area`.

---

## 6. Verification Plan

### Automated Tests
```bash
# 1. Update dokumentasi
Rscript -e "devtools::document()"

# 2. Uji diagnostik lintas model
Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test_diagnostic.R')"
```

### Manual Verification
Jalankan `sae_diagnostic()` pada `eblup_fh`, `eblup_sfh`, dan `ebp_area` menggunakan dataset `mys`. Periksa bahwa laporan diagnostik tercetak dengan indikator kelayakan yang akurat.
