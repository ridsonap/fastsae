# Ringkasan Implementasi: Pembangkitan Data Deret Waktu Spatio-Temporal (`sim_series_data`)

## 1. Fitur yang Telah Diimplementasikan

### A. Pembangkit Data Deret Waktu Spatio-Temporal (`sim_series_data`)
File: [`R/sim_series_data.R`](file:///Volumes/work/_MainR/fastsae-inla/R/sim_series_data.R)
Fungsi ini mensimulasikan panel data berdimensi $D \times T$ yang memadukan autokorelasi spasial area (SAR) dan dinamika temporal AR(1):
* **Struktur Acak Spatio-Temporal**:
  - Efek spasial area $u_{1, d}$ terstruktur عبر matriks ketetanggaan spasial $W$ ($SAR$ process dengan parameter $\rho_s$).
  - Efek waktu dinamis $u_{2, d, t}$ mengikuti proses autoregresif AR(1) stasioner ($u_{2, d, t} = \rho_t u_{2, d, t-1} + \epsilon_t$).
  - Komponen tren deterministik: `trend = c("linear", "random_walk", "ar1", "none")`.
* **Respon Multi-Distribusi Dinamis**:
  - `y_gaussian` & `vardir`: Fay-Herriot spatio-temporal.
  - `y_poisson` & `exposure`: Pemetaan kejadian/penyakit dengan offset populasi bertumbuh.
  - `y_binomial` & `trials`: Proporsi survei berulang (panel Susenas/Sakernas).
  - `y_beta`: Proporsi kontinu $(0, 1)$ untuk regresi Beta.
  - `y_nbinomial`: Cacahan terdispersi (*overdispersed count*).
  - `y_gamma`: Variabel kontinu positif miring.
* **Urutan Baris (*Sorting Order*)**:
  - Mendukung `sort_order = "domain-major"` (standar wajib untuk [`eblup_stfh()`](file:///Volumes/work/_MainR/fastsae-inla/R/eblup_stfh.R) dan `sae::eblupSTFH`).
  - Mendukung `sort_order = "time-major"` untuk analisis per gelombang tahun.
* **Pola Observasi Hilang / Unsampled**:
  - `n_unsampled`: Domain yang persisten tak tersampel di semua tahun (`NA`).
  - `prop_intermittent`: Observasi yang hilang berkala pada tahun-tahun tertentu.
* **Objek Kembalian**: Berkelas `fastsae_sim_series` yang dilengkapi dengan `print()` method yang informatif, menyimpan data frame `$data`, matriks spasial `$W` & `$W_std`, serta *ground-truth* efek acak `$u_spatial` dan `$u_temporal`.

### B. Dataset Bawaan: `sae_panel_multi`
File: [`data/sae_panel_multi.rda`](file:///Volumes/work/_MainR/fastsae-inla/data/sae_panel_multi.rda) dan dokumentasi di [`R/data.R`](file:///Volumes/work/_MainR/fastsae-inla/R/data.R)
* Dataset panel 42 domain $\times$ 5 periode tahun (2022-2026, total 210 baris observasi).
* Kompatibel langsung dengan matriks spasial bawaan [`mys_proxmat`](file:///Volumes/work/_MainR/fastsae-inla/data/mys_proxmat.rda).

---

## 2. Contoh Penggunaan Cepat

```r
library(fastsae)

# 1. Bangkitkan panel 30 domain sepanjang 4 tahun
sim_panel <- sim_series_data(
  D = 30,
  T = 4,
  time_start = 2021,
  rho_s = 0.5,
  rho_t = 0.6,
  n_unsampled = 0,
  prop_intermittent = 0,
  seed = 123
)
print(sim_panel)
head(sim_panel$data)

# 2. Estimasi Model Spatio-Temporal Fay-Herriot (eblup_stfh)
fit_stfh <- eblup_stfh(
  y_gaussian ~ x1 + x2,
  data = sim_panel$data,
  domain = ~area,
  time = ~year,
  vardir = ~vardir,
  W = sim_panel$W_std
)
summary(fit_stfh)

# 3. Menggunakan dataset bawaan sae_panel_multi
data(sae_panel_multi)
head(sae_panel_multi)
```

---

## 3. Hasil Pengujian Unit

Seluruh 35 pengujian pada [`tests/testthat/test_sim_series.R`](file:///Volumes/work/_MainR/fastsae-inla/tests/testthat/test_sim_series.R) **lulus 100%**:
- Validasi dimensi $(D \cdot T)$, urutan domain-major vs time-major.
- Validasi sifat batasan nilai sebaran (Beta $\in (0, 1)$, Poisson $\ge 0$, Binomial $\le$ trials).
- Validasi domain tak tersampel persisten dan intermittent.
- Validasi fitting model `eblup_stfh()` pada data panel hasil generate dan dataset bawaan `sae_panel_multi`.
