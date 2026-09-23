# Rencana Implementasi: Dataset Hasil Generate dengan Multi-Distribusi (Poisson, Beta, Binomial, dll.)

## 1. Deskripsi Tujuan (Goal Description)

Menambahkan dataset simulasi area-level multi-respons (**`sae_area_multi`**) ke dalam paket [`fastsae`](file:///Volumes/work/_MainR/fastsae-inla). Dataset ini dirancang khusus untuk memfasilitasi pengujian, demonstrasi, dan analisis Small Area Estimation lintas berbagai keluarga distribusi (Poisson, Beta, Binomial, Gaussian, Negative Binomial, Gamma) dan struktur spasial (menggunakan matriks ketetanggaan `mys_proxmat`).

Selain menyediakan dataset bawaan yang siap pakai (`data/sae_area_multi.rda`), kita juga menyediakan fungsi generator **`sim_area_data()`** agar pengguna dapat membangkitkan data sintetis baru kapan saja dengan jumlah domain ($D$), matriks spasial, dan seed acak yang fleksibel.

---

## 2. Struktur Dataset yang Dibangkitkan (`sae_area_multi`)

Dataset terdiri dari **42 domain** (sesuai dimensi matriks spasial bawaan `mys_proxmat` $42 \times 42$), memuat 6 domain tak tersampel (*unsampled* dengan nilai `NA`), dengan rincian kolom:

| Nama Kolom | Tipe / Distribusi | Keterangan & Kegunaan SAE |
|:---|:---|:---|
| **`domain`** | Karakter / Faktor | Identitas unik domain (`"Area_01"` s/d `"Area_42"`) |
| **`x1`**, **`x2`** | Kontinu | Variabel kovariat / indikator pendukung (auxiliary variables) |
| **`y_gaussian`** | Gaussian | Respon kontinu (Fay-Herriot klasik) |
| **`vardir`** | Positif | Varians sampling estimator langsung untuk `y_gaussian` |
| **`y_poisson`** | Integer / Cacahan | Respon cacahan (angka kasus penyakit / kemiskinan ekstrim) |
| **`exposure`** | Integer Positif | Offset populasi / expected cases ($E_d$) untuk Poisson & NegBinom |
| **`y_binomial`** | Integer / Sukses | Jumlah sukses / kejadian (misal jumlah ruta miskin dalam sampel) |
| **`trials`** | Integer Positif | Jumlah ukuran sampel per area ($n_d$, $y_{bin} \le n_d$) |
| **`y_nbinomial`**| Integer / Cacahan | Data cacahan dengan overdispersi |
| **`y_beta`** | Kontinu $(0, 1)$ | Respon proporsi/indeks kontinu strictly in $(0, 1)$ |
| **`y_gamma`** | Kontinu Positif | Respon kontinu bernilai positif dan miring (pengeluaran agregat) |

> [!NOTE]
> **Area Tidak Tersampel (*Unsampled Domains*)**:
> Domain pada indeks `c(5, 12, 19, 27, 35, 41)` memiliki nilai `NA` pada seluruh kolom respon `y_...`. Hal ini memungkinkan pengguna menguji kemampuan prediksi sintetis dan *spatial kriging* INLA pada seluruh jenis distribusi.

---

## 3. Rencana Perubahan Detail (Proposed Changes)

### Komponen 1: Pembangkitan Data & Penyimpanan ke `data/`
#### [NEW] [`data/sae_area_multi.rda`](file:///Volumes/work/_MainR/fastsae-inla/data/sae_area_multi.rda)
- File binary dataset R (.rda) berukuran ringan, otomatis dimuat via `LazyData: true` (sehingga pengguna cukup memanggil `data(sae_area_multi)` atau langsung `sae_area_multi`).

#### [NEW] [`R/sim_area_data.R`](file:///Volumes/work/_MainR/fastsae-inla/R/sim_area_data.R)
- Fungsi generator fleksibel:
  ```r
  sim_area_data <- function(
    D = 42,
    W = NULL,
    seed = NULL,
    n_unsampled = 6
  )
  ```
  Menghasilkan data frame dengan kolom yang sama persis, berguna untuk uji simulasi berulang (*Monte Carlo simulation study*).

---

### Komponen 2: Dokumentasi Dataset
#### [MODIFY] [`R/data.R`](file:///Volumes/work/_MainR/fastsae-inla/R/data.R)
- Menambahkan dokumentasi roxygen2 lengkap untuk `sae_area_multi`:
  - Format data frame, deskripsi tiap kolom, dan contoh pemodelan untuk masing-masing distribusi (Gaussian, Binomial, Poisson, Beta, Gamma).

---

### Komponen 3: Unit Testing & Skrip Contoh
#### [NEW] [`tests/testthat/test_sae_area_multi.R`](file:///Volumes/work/_MainR/fastsae-inla/tests/testthat/test_sae_area_multi.R)
- Menguji integritas dataset `sae_area_multi` (dimensi, tipe data, ketiadaan nilai invalid).
- Menguji fitting `ebp()` dan `ebp_spatial()` pada masing-masing kolom respon (`y_poisson`, `y_beta`, `y_binomial`, `y_gaussian`, `y_gamma`, `y_nbinomial`).
- Menguji fungsi `sim_area_data()`.

#### [MODIFY] [`main.R`](file:///Volumes/work/_MainR/fastsae-inla/main.R)
- Menambahkan blok contoh interaktif menjalankan model EBP untuk semua jenis distribusi menggunakan data `sae_area_multi`.

---

## 4. Verification Plan

### Automated Tests
```bash
# 1. Generate dokumentasi paket (man/sae_area_multi.Rd)
Rscript -e "devtools::document()"

# 2. Jalankan unit test dataset baru
Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test_sae_area_multi.R')"

# 3. Jalankan test suite ebp_area yang sudah ada
Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test_ebp_area.R')"
```

### Manual Verification
Jalankan skrip di `main.R`:
```r
# Cek pemodelan Poisson Spasial (BYM2)
m_pois <- ebp_spatial(y_poisson ~ x1 + x2, data = sae_area_multi, exposure = "exposure", family = "poisson", spatial = "bym2", W = mys_proxmat)
summary(m_pois)

# Cek pemodelan Beta (Proporsi)
m_beta <- ebp(y_beta ~ x1 + x2, data = sae_area_multi, family = "beta")
summary(m_beta)

# Cek pemodelan Binomial Spasial
m_bin <- ebp_spatial(y_binomial ~ x1 + x2, data = sae_area_multi, trials = "trials", family = "binomial", spatial = "bym2", W = mys_proxmat)
summary(m_bin)
```
