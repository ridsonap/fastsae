# Ringkasan Implementasi: Pembangkitan Data Multi-Distribusi & Matriks Spasial $W$

## 1. Fitur yang Telah Diimplementasikan

### A. Pembangkit Matriks Ketetanggaan / Bobot Spasial (`sim_spatial_weights`)
File: [`R/sim_spatial_weights.R`](file:///Volumes/work/_MainR/fastsae-inla/R/sim_spatial_weights.R)
Fungsi ini menghasilkan matriks spasial sintetis dengan topologi yang realistis:
* **`type = "knn"`** (*k-Nearest Neighbors*): Menghasilkan koordinat acak 2D dan menghubungkan $k$-tetangga terdekat (mutual adjacency simetris). Sangat cocok untuk mensimulasikan batas wilayah administratif riil.
* **`type = "grid"`** (*Regular Lattice*): Menata domain dalam kisi 2D ($r \times c$) menggunakan *Rook contiguity* (berbagi sisi bersama).
* **`type = "ring"`** (*1D Circular Lattice*): Domain tersambung dalam cincin melingkar ($i-1$ dan $i+1$).
* **`style = "B"`**: Matriks biner $0/1$ simetris tanpa diagonal (standar untuk model graf INLA: **BYM2**, **BYM**, **Besag**).
* **`style = "W"`**: Matriks *row-standardized* ($\sum_j W_{ij} = 1$, standar untuk model autoregresif **SAR**).

### B. Pembangkit Data Tabular Multi-Respons (`sim_area_data`)
File: [`R/sim_area_data.R`](file:///Volumes/work/_MainR/fastsae-inla/R/sim_area_data.R)
Membangkitkan dataset simulasi *area-level* dengan kovariat bersama dan autokorelasi spasial gabungan (SAR + BYM2):
* **`y_gaussian`** & **`vardir`**: Model Fay-Herriot standar kontinu dengan sampling error heterogen.
* **`y_poisson`** & **`exposure`**: Model cacahan dengan *offset exposure* populasi.
* **`y_binomial`** & **`trials`**: Model binomial proporsi/cacahan sukses dari ukuran sampel domain.
* **`y_beta`**: Model proporsi kontinu strictly in $(0, 1)$ untuk regresi Beta.
* **`y_nbinomial`**: Model cacahan dengan dispersi berlebih (*overdispersed count*).
* **`y_gamma`**: Model kontinu miring positif (*positive skewed continuous*).
* **Penanganan Domain Tak Tersampel (`n_unsampled`)**: Secara otomatis menandai nilai respon sebagai `NA` sementara variabel pendukung tetap tersedia (sesuai kerangka sensus riil).
* Mengembalikan objek kelas `fastsae_sim_data` dengan fungsi `print()` informatif, memuat komponen `$data`, `$W`, `$W_std`, `$coords`, serta *ground-truth* efek acak `$u_spatial` dan `$u_iid`.

### C. Dataset Bawaan: `sae_area_multi`
File: [`data/sae_area_multi.rda`](file:///Volumes/work/_MainR/fastsae-inla/data/sae_area_multi.rda) dan dokumentasi di [`R/data.R`](file:///Volumes/work/_MainR/fastsae-inla/R/data.R)
* Dataset 42 domain yang kompatibel 100% dengan matriks spasial bawaan [`mys_proxmat`](file:///Volumes/work/_MainR/fastsae-inla/data/mys_proxmat.rda).
* Pengguna dapat langsung memuat via `data(sae_area_multi)` tanpa konfigurasi tambahan.

---

## 2. Contoh Penggunaan Cepat

```r
library(fastsae)

# 1. Bangkitkan matriks spasial KNN 30 domain
W_knn <- sim_spatial_weights(D = 30, type = "knn", k = 4, style = "B", seed = 123)

# 2. Bangkitkan data multi-distribusi dengan efek spasial
sim <- sim_area_data(W = W_knn, rho = 0.5, phi = 0.6, n_unsampled = 4, seed = 2026)
print(sim)
head(sim$data)

# 3. Model Poisson Spasial (BYM2)
fit_pois <- ebp_spatial(
  y_poisson ~ x1 + x2,
  data = sim$data,
  exposure = "exposure",
  family = "poisson",
  W = sim$W,
  spatial = "bym2"
)
summary(fit_pois)

# 4. Menggunakan dataset bawaan sae_area_multi dan mys_proxmat
data(sae_area_multi)
data(mys_proxmat)

fit_bin <- ebp_spatial(
  y_binomial ~ x1 + x2,
  data = sae_area_multi,
  trials = "trials",
  family = "binomial",
  W = mys_proxmat
)
summary(fit_bin)
```

---

## 3. Hasil Pengujian Unit

Seluruh 51 pengujian pada [`tests/testthat/test_sim_data.R`](file:///Volumes/work/_MainR/fastsae-inla/tests/testthat/test_sim_data.R) lulus 100%:
* `sim_spatial_weights`: validasi simetri, diagonal 0, rowSums untuk style W, struktur grid contiguity, ring contiguity, dan reprodusibilitas seed.
* `sim_area_data`: validasi semua kolom respons, batasan nilai (Beta $\in (0, 1)$, Poisson $\ge 0$, Binomial $\le$ trials), dan domain tak tersampel (`NA`).
* Integrasi `ebp` dan `ebp_spatial` pada data hasil generate dan dataset `sae_area_multi`.
