# Rencana Implementasi: Fitur Pembangkitan Data Deret Waktu Spatio-Temporal (`sim_series_data`)

## 1. Deskripsi Tujuan (Goal Description)

Menambahkan fitur pembangkitan data simulasi **deret waktu / spatio-temporal panel** (`sim_series_data`) ke dalam paket [`fastsae`](file:///Volumes/work/_MainR/fastsae-inla). Fitur ini memperluas fungsionalitas simulasi cross-sectional yang sudah ada (`sim_area_data`) menjadi data multi-dimensi (Domain $\times$ Waktu, $D \times T$).

Fitur ini dirancang untuk:
1. **Mensimulasikan dinamika temporal AR(1)** untuk setiap domain $d$ sepanjang periode $t = 1, \dots, T$, dipadukan dengan autokorelasi spasial antar domain $W$ (SAR / GMRF).
2. **Menghasilkan respon multi-distribusi simultan** sepanjang waktu:
   - Kontinu Gaussian (`y_gaussian`, `vardir`) yang langsung kompatibel dengan [`eblup_stfh()`](file:///Volumes/work/_MainR/fastsae-inla/R/eblup_stfh.R) (Marhuenda et al., 2013).
   - Cacahan Poisson (`y_poisson`, `exposure`) dengan offset populasi bervariasi terhadap waktu.
   - Cacahan/Proporsi Binomial (`y_binomial`, `trials`) dari survei berulang (panel Susenas/Sakernas).
   - Proporsi Kontinu Beta (`y_beta`) $\in (0, 1)$.
   - Cacahan Terdispersi Negatif Binomial (`y_nbinomial`).
   - Kontinu Miring Gamma (`y_gamma`).
3. **Mendukung pola observasi tak tersampel realistis**:
   - Area tidak tersampel persisten (*persistent unsampled domains*).
   - Pengamatan hilang berkala (*intermittent non-sampling* pada tahun tertentu).
4. **Menyediakan dataset bawaan siap pakai (`sae_panel_multi`)**:
   - Panel 42 domain $\times$ 5 periode waktu (210 baris observasi) berbasis matriks spasial [`mys_proxmat`](file:///Volumes/work/_MainR/fastsae-inla/data/mys_proxmat.rda).

---

## 2. Formulasi Matematis Model Spatio-Temporal

Model spatio-temporal area-level mengikuti formulasi baku Marhuenda, Molina & Morales (2013) dan Rao & Molina (2015):

Untuk domain $d = 1, \dots, D$ dan waktu $t = 1, \dots, T$:
$$\eta_{d, t} = \mathbf{x}_{d, t}^\top \boldsymbol{\beta} + \gamma \cdot t + u_{1, d} + u_{2, d, t}$$

1. **Efek Kovariat & Tren Waktu**:
   - $\mathbf{x}_{d, t} = (x_{1, d, t}, x_{2, d, t})^\top$: Kovariat dinamis yang bervariasi menurut area dan waktu.
   - $\gamma \cdot t$: Komponen tren deterministik (bisa linear, random walk, atau stationary).
2. **Efek Acak Spasial Area ($u_{1, d}$)**:
   - Terstruktur spasial via proses SAR: $\mathbf{u}_1 = \sigma_1 (I_D - \rho_1 W_{std})^{-1} \boldsymbol{\epsilon}_1$, di mana $\boldsymbol{\epsilon}_1 \sim \mathcal{N}(0, I_D)$.
3. **Efek Acak Spatio-Temporal Dinamis ($u_{2, d, t}$)**:
   - Mengikuti proses autoregresif AR(1) independen untuk tiap domain:
     $$u_{2, d, t} = \rho_2 u_{2, d, t-1} + \epsilon_{2, d, t}, \quad \epsilon_{2, d, t} \sim \mathcal{N}\left(0, \sigma_2^2 (1 - \rho_2^2)\right)$$
     (dengan varians marjinal stasioner $\text{Var}(u_{2, d, t}) = \sigma_2^2$).
4. **Prediktor Linear Gabungan**:
   $$\eta_{d, t} = \beta_0 + \beta_1 x_{1, d, t} + \beta_2 x_{2, d, t} + \text{trend}_t + u_{1, d} + u_{2, d, t}$$

---

## 3. Desain API `sim_series_data()`

```r
sim_series_data(
  D = 42,
  T = 5,
  time_start = 2022,
  W = NULL,
  spatial_type = c("knn", "grid", "ring"),
  rho_s = 0.5,             # Autokorelasi spasial area (rho1 pada eblup_stfh)
  rho_t = 0.6,             # Autokorelasi temporal AR(1) (rho2 pada eblup_stfh)
  sigma_s = 0.4,           # Standar deviasi efek spasial area (u1)
  sigma_t = 0.3,           # Standar deviasi efek spatio-temporal AR(1) (u2)
  trend = c("linear", "random_walk", "ar1", "none"),
  trend_slope = 0.05,      # Kemiringan tren waktu
  beta = c(1.0, 0.8, -0.5),# Koefisien regresi fixed effects
  n_unsampled = 4,         # Domain tak tersampel di semua tahun
  prop_intermittent = 0.05,# Proporsi observasi hilang acak pada tahun tertentu
  sort_order = c("domain-major", "time-major"), # Format urutan baris
  seed = NULL
)
```

### Format Nilai Kembalian (`fastsae_sim_series`):
* **`$data`**: Data frame $(D \cdot T) \times 15$ dengan kolom:
  - `area` (atau `domain`): ID domain
  - `year` (atau `time`): Periode waktu ($2022, 2023, \dots$)
  - `x1`, `x2`: Kovariat dinamis
  - `y_gaussian`, `vardir`: Respon kontinu & varians sampling
  - `y_poisson`, `exposure`: Respon cacahan & exposure offset
  - `y_binomial`, `trials`: Respon binomial & jumlah trials
  - `y_beta`: Respon proporsi $(0, 1)$
  - `y_nbinomial`: Respon overdispersed count
  - `y_gamma`: Respon kontinu positif
  - `x_coord`, `y_coord`: Koordinat centroid domain
* **`$W`**: Matriks ketetanggaan spasial $D \times D$.
* **`$W_std`**: Matriks spasial row-standardized $D \times D$.
* **`$u_spatial`**: Vektor efek spasial sejati $u_{1, d}$.
* **`$u_temporal`**: Matriks $D \times T$ efek waktu sejati $u_{2, d, t}$.
* **`$parameters`**: List nilai parameter sejati ($\rho_s, \rho_t, \sigma_s, \sigma_t, \boldsymbol{\beta}, \gamma$).

---

## 4. User Review Required

> [!IMPORTANT]
> **Pilihan Urutan Data (Sorting Order)**:
> Fungsi `eblup_stfh()` pada `fastsae` (dan paket `sae::eblupSTFH`) secara ketat mewajibkan urutan **domain-major order** (seluruh tahun untuk domain 1 lebih dahulu, lalu domain 2, dst.).
> Secara *default*, `sim_series_data()` akan menyajikan data dalam urutan **`domain-major`** agar pengguna dapat langsung memanggil:
> ```r
> eblup_stfh(y_gaussian ~ x1 + x2, data = sim$data, domain = ~area, time = ~year, vardir = ~vardir, W = sim$W_std)
> ```
> Opsi `sort_order = "time-major"` juga disediakan jika pengguna ingin menganalisis per panel tahunan (misal `year` 2022 untuk semua area, lalu 2023, dst.).

---

## 5. Proposed Changes (Daftar File)

### Komponen 1: Fungsi Pembangkit Spatio-Temporal
* **[NEW] [`R/sim_series_data.R`](file:///Volumes/work/_MainR/fastsae-inla/R/sim_series_data.R)**:
  - Implementasi fungsi `sim_series_data()` lengkap dengan parameter spasial dan AR(1) temporal.
  - S3 method `print.fastsae_sim_series()` untuk mencetak ringkasan dimensi $D \times T$, missing values, parameter autokorelasi, dan variabel yang tersedia.

### Komponen 2: Dataset Bawaan Spatio-Temporal
* **[NEW] [`data/sae_panel_multi.rda`](file:///Volumes/work/_MainR/fastsae-inla/data/sae_panel_multi.rda)**:
  - Dataset panel sintetis 42 domain $\times$ 5 periode (210 baris) menggunakan struktur spasial [`mys_proxmat`](file:///Volumes/work/_MainR/fastsae-inla/data/mys_proxmat.rda).
* **[MODIFY] [`R/data.R`](file:///Volumes/work/_MainR/fastsae-inla/R/data.R)**:
  - Dokumentasi roxygen2 untuk `sae_panel_multi`.

### Komponen 3: Konfigurasi & Ekspor Paket
* **[MODIFY] [`DESCRIPTION`](file:///Volumes/work/_MainR/fastsae-inla/DESCRIPTION)**:
  - Menambahkan `sim_series_data.R` pada `Collate:`.
* **[MODIFY] [`NAMESPACE`](file:///Volumes/work/_MainR/fastsae-inla/NAMESPACE)**:
  - Mengekspor `sim_series_data` dan mendaftarkan S3 method `print.fastsae_sim_series`.

### Komponen 4: Pengujian & Contoh
* **[NEW] [`tests/testthat/test_sim_series.R`](file:///Volumes/work/_MainR/fastsae-inla/tests/testthat/test_sim_series.R)**:
  - Uji validasi struktur dimensi ($D \cdot T$), urutan domain-major vs time-major.
  - Uji sifat AR(1) temporal dan korelasi spasial.
  - Uji penanganan missing data (persisten dan intermittent).
  - Uji integrasi langsung dengan `eblup_stfh()`.
  - Uji dataset bawaan `sae_panel_multi`.
* **[MODIFY] [`main.R`](file:///Volumes/work/_MainR/fastsae-inla/main.R)**:
  - Menambahkan contoh praktis simulasi deret waktu dan estimasi model spatio-temporal.

---

## 6. Verification Plan

### Automated Tests
```bash
# 1. Generate dokumentasi dan verifikasi roxygen2
Rscript -e "devtools::document()"

# 2. Jalankan test suite pembangkitan deret waktu spatio-temporal
Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test_sim_series.R')"

# 3. Jalankan pengujian estimasi spatio-temporal yang ada
Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test_eblup_stfh.R')"
```

### Manual Verification
1. Bangkitkan panel 30 domain sepanjang 4 tahun:
   ```r
   sim <- sim_series_data(D = 30, T = 4, time_start = 2021, seed = 123)
   ```
2. Jalankan estimasi `eblup_stfh`:
   ```r
   fit_stfh <- eblup_stfh(
     y_gaussian ~ x1 + x2,
     data = sim$data,
     domain = ~area,
     time = ~year,
     vardir = ~vardir,
     W = sim$W_std
   )
   ```
3. Pastikan algoritma Fisher-scoring konvergen dan mengembalikan estimasi varians/autokorelasi $\rho_1$ dan $\rho_2$ yang mendekati parameter simulasi.
