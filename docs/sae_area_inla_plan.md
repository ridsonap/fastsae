# Rencana Implementasi: Area-Level Small Area Estimation dengan INLA dan Laplace Approximation

## 1. Deskripsi Tujuan (Goal Description)

Implementasi model **Area-Level Small Area Estimation** berbasis **INLA (Integrated Nested Laplace Approximations)** dan **Laplace Approximation** ke dalam paket [`fastsae`](file:///Volumes/work/_MainR/fastsae-inla).

### Arsitektur Desain Fungsi (Unified & Modular API)

Sesuai masukan terbaik dari pengguna: \* **Fungsi Utama**: **`ebp_area()`** sebagai *single unified entry-point*. - Memiliki argumen `spatial = c("none", "bym2", "bym", "besag", "generic1", "slm")`. - Default `spatial = "none"` untuk model non-spasial (IID random intercept / Generalized Fay-Herriot). - Ketika `spatial = "bym2"` (atau `"bym"`, `"besag"`, dll.), pengguna cukup menambahkan matriks spasial `W`. \* **Fungsi Khusus / Wrapper**: **`ebp_spatial_area()`** - Disediakan sebagai shortcut khusus model spasial (dengan argumen `W` wajib), konsisten dengan pasangan `eblup_fh()` dan `eblup_sfh()` yang sudah ada di `fastsae`. - Mengarahkan komputasi langsung ke `ebp_area(..., spatial = spatial, W = W)`.

------------------------------------------------------------------------

## 2. Keuntungan Pendekatan Unified `ebp_area(..., spatial = "none")`

1.  **Perbandingan Model Sangat Mudah**: Pengguna dapat dengan mudah membandingkan model non-spasial vs spasial pada data yang sama hanya dengan mengubah 1 parameter:

    ``` r
    # 1. Model Non-Spasial
    m0 <- ebp_area(y ~ x1 + x2, data = mys, vardir = "vardir")

    # 2. Model Spasial BYM2 (cukup tambahkan spatial & W)
    m1 <- ebp_area(y ~ x1 + x2, data = mys, vardir = "vardir", spatial = "bym2", W = mys_proxmat)

    # Evaluasi model selection (DIC & WAIC)
    m0$goodness
    m1$goodness
    ```

2.  **Fleksibilitas Penuh**: Pengguna yang menyukai fungsi terpisah tetap dapat memanggil `ebp_spatial_area(..., W = mys_proxmat)`.

------------------------------------------------------------------------

## 3. Cakupan Model & Likelihood

```         
ebp_area(formula, data, family = ..., spatial = ...)
│
├── 1. Pilihan Distribusi / Likelihood (family)
│   ├── "gaussian":   Fay-Herriot klasik (vardir diketahui) atau homoscedastic error (vardir = NULL)
│   ├── "binomial":   Proporsi / binary count, link logit/probit, argumen `trials = n_d`
│   ├── "poisson":    Cacahan kasus / angka langka, link log, argumen `exposure = E_d`
│   ├── "nbinomial":  Cacahan terdispersi (overdispersion), link log
│   ├── "beta":       Proporsi kontinu (0, 1), link logit
│   └── "gamma":      Nilai kontinu positif condong (skewed positive), link log
│
└── 2. Pilihan Struktur Spasial (spatial)
    ├── "none":       IID random intercept (non-spasial)
    ├── "bym2":       Besag-York-Mollié 2 (scaled, memisahkan varians spasial & iid via parameter mixing φ)
    ├── "bym":        Besag-York-Mollié klasik
    ├── "besag":      Intrinsic Conditional Autoregressive (ICAR) murni
    ├── "generic1":   Leroux model (Q = τ [ρ (D - W) + (1 - ρ) I])
    └── "slm":        Spatial Lag Model (SAR)
```

------------------------------------------------------------------------

## 4. Tanda Tangan Fungsi (API Design)

``` r
ebp_area <- function(
  formula,
  data,
  domain = NULL,
  family = c("gaussian", "binomial", "poisson", "nbinomial", "beta", "gamma"),
  spatial = c("none", "bym2", "bym", "besag", "generic1", "slm"),
  W = NULL,                # Matriks ketetanggaan/bobot spasial, listw, atau nb
  vardir = NULL,           # Khusus Gaussian: varians sampling per area
  trials = NULL,           # Khusus Binomial: jumlah sampel/trials (n_d)
  exposure = NULL,         # Khusus Poisson/NegBinom: offset/populasi (E_d)
  method = c("inla", "laplace"),
  strategy = c("simplified.laplace", "laplace", "gaussian"),
  link = NULL,             # Default: identity (gaussian), logit (binomial), log (poisson)
  scale_model = TRUE,      # Penyekalaan graf untuk interpretasi varians stabil
  prior_prec = list(prior = "pc.prec", param = c(1, 0.01)),
  prior_phi = list(prior = "pc", param = c(0.5, 0.5)), # Khusus BYM2: porsi spasial
  print_result = TRUE,
  ...
)

ebp_spatial_area <- function(
  formula,
  data,
  W,
  spatial = c("bym2", "bym", "besag", "generic1", "slm"),
  ...
) {
  spatial <- match.arg(spatial)
  ebp_area(formula = formula, data = data, W = W, spatial = spatial, ...)
}
```

------------------------------------------------------------------------

## 5. Output Objek S3 (`class = c("fastsae_ebp_area", "fastsae")`)

Objek hasil estimasi memuat: \* `df_ebp`: Data frame lengkap hasil estimasi per domain: - `domain`: Nama / ID domain - `y`: Nilai observasi asli (bisa `NA` untuk area yang tidak tersampel) - `ebp`: Estimasi titik pada skala asli (misal: mean, proporsi $p_d$, laju $\lambda_d$) - `linear_pred`: Prediktor linier $\eta_d = \mathbf{x}_d^\top \boldsymbol{\beta} + \psi_d$ - `sd`: Standar deviasi posterior - `mse`: Mean Squared Error ($\approx \text{sd}^2$) - `rse`: Relative Standard Error / CV (%) = $(\text{sd} / |\text{ebp}|) \times 100\%$ - `ci_lower`, `ci_upper`: 95% Credible Interval (persentil 2.5% dan 97.5%) - `random_effect`: Efek acak area $\psi_d$ \* `estcoef`: Koefisien regresi tetap ($\boldsymbol{\beta}$) dengan SD, z-value, p-value / 95% CI. \* `hyperpar`: Estimasi hiperparameter (varians acak $\sigma_u^2$, presisi $\tau$, porsi spasial $\phi$). \* `goodness`: DIC, WAIC, Marginal Log-Likelihood, CPO. \* `fit`: Objek mentah INLA / model fit.

------------------------------------------------------------------------

## 6. Proposed Changes (Daftar File)

### Komponen 1: Konfigurasi & Utilitas Helper

- **[MODIFY] [`DESCRIPTION`](file:///Volumes/work/_MainR/fastsae-inla/DESCRIPTION)**: Daftarkan `INLA` dan `spdep` di `Suggests:`.
- **[NEW] [`R/inla_utils.R`](file:///Volumes/work/_MainR/fastsae-inla/R/inla_utils.R)**:
  - Pengecekan paket INLA (`.check_inla_installed()`).
  - Helper konversi matriks spasial $W$ (matrix, sparse Matrix, `nb`, `listw`) ke graf INLA biner simetris.
  - Helper ekstraksi hasil fitting INLA ke data frame standar `df_ebp`.

### Komponen 2: Implementasi Fungsi Inti

- **[NEW] [`R/ebp_area.R`](file:///Volumes/work/_MainR/fastsae-inla/R/ebp_area.R)**: Fungsi `ebp_area()` (menyatukan model non-spasial dan spasial).
- **[NEW] [`R/ebp_spatial_area.R`](file:///Volumes/work/_MainR/fastsae-inla/R/ebp_spatial_area.R)**: Wrapper `ebp_spatial_area()`.

### Komponen 3: Metode S3 & Visualisasi

- **[MODIFY] [`R/methods.R`](file:///Volumes/work/_MainR/fastsae-inla/R/methods.R)**: S3 methods: `print.fastsae_ebp_area`, `summary.fastsae_ebp_area`, `coef`, `fitted`, `residuals`.
- **[MODIFY] [`R/autoplot.R`](file:///Volumes/work/_MainR/fastsae-inla/R/autoplot.R)**: Visualisasi ggplot2: perbandingan direct vs ebp, ranking domain dengan 95% CI, dan distribusi RSE.

### Komponen 4: Pengujian & Dokumentasi

- **[NEW] [`tests/testthat/test_ebp_area.R`](file:///Volumes/work/_MainR/fastsae-inla/tests/testthat/test_ebp_area.R)**: Test suite komprehensif untuk Gaussian, Binomial, Poisson, NegBinom, Beta, Gamma, non-spasial vs spasial BYM2/Besag, serta area unsampled.
- **[MODIFY] [`main.R`](file:///Volumes/work/_MainR/fastsae-inla/main.R)**: Skrip contoh interaktif menjalankan `ebp_area()` dan `ebp_spatial_area()`.

------------------------------------------------------------------------

## 7. Verification Plan

### Automated Tests

``` bash
# 1. Dokumentasi & NAMESPACE
Rscript -e "devtools::document()"

# 2. Jalankan test suite
Rscript -e "testthat::test_file('tests/testthat/test_ebp_area.R')"
Rscript -e "devtools::test()"
```

### Manual Verification

1.  Verifikasi `ebp_area(spatial = "none")` dan `ebp_area(spatial = "bym2")` pada dataset `mys`.
2.  Verifikasi kesesuaian estimasi area unsampled (`NA`).
3.  Verifikasi `autoplot(fit)` menghasilkan visualisasi diagnostik yang rapi.
