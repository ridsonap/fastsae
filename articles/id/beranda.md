# fastsae: Estimasi Area Kecil Frekuentis dan Bayesian Berkinerja Tinggi di R

[![CRAN
status](https://www.r-pkg.org/badges/version/fastsae)](https://CRAN.R-project.org/package=fastsae)
[![R-CMD-check](https://github.com/ridsonap/fastsae/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ridsonap/fastsae/actions/workflows/R-CMD-check.yaml)
[![License:
GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**fastsae** menyediakan model-model Estimasi Area Kecil (*Small Area
Estimation*, SAE) yang lazim dipakai untuk data survei. Estimasi
frekuentis dijalankan dalam kode **C++** terkompilasi (`RcppArmadillo`)
dengan *multi-threading* **OpenMP**. Estimasi Bayesian secara bawaan
memakai **INLA** (*Integrated Nested Laplace Approximations*), dan model
yang sama (kecuali efek *spatial lag* dan RW2) juga dapat diestimasi
dengan MCMC penuh melalui **NumPyro di atas JAX**, baik di CPU maupun
GPU. Untuk ratusan hingga ribuan domain, fastsae jauh lebih cepat
daripada implementasi R yang sudah ada (lihat [Ringkasan
Kinerja](#ringkasan-kinerja); untuk beberapa puluh domain, `sae` kurang
lebih sama cepatnya), dan hasilnya tetap sesuai secara numerik dengan
paket rujukan seperti `sae`.

------------------------------------------------------------------------

## Mengapa fastsae?

- ⚡ **Kecepatan**: dengan 1.000 domain, estimasi *Fisher scoring* di
  C++ beserta MSE analitiknya 29 hingga 1.250 kali lebih cepat daripada
  `sae` (model spatio-temporal: 40 detik, bukan 19 menit). Untuk
  beberapa puluh domain, selisihnya hanya sepersekian detik.
- 🎯 **Kesesuaian numerik**: pada data yang sama, koefisien, komponen
  varians, EBLUP, dan MSE analitiknya sama dengan `sae` (Molina &
  Marhuenda) hingga sekitar $`10^{-14}`$. Pengecualian yang
  terdokumentasi hanya dua: korelasi spasial yang berada di batas ruang
  parameter dan MSE ML pada model spasial.
- 🔀 **Frekuentis dan Bayesian dalam satu paket**: model EBLUP di C++,
  serta model Bayesian hierarkis melalui INLA atau, dengan
  `engine = jax_engine()`, melalui MCMC (NUTS NumPyro di atas JAX, di
  CPU atau GPU) lengkap dengan diagnostik konvergensinya.
- 🧵 **Bootstrap paralel**: estimasi MSE dengan bootstrap parametrik
  maupun nonparametrik, berjalan di banyak inti melalui OpenMP.
- 🗺️ **Area tak tersampel**: area, subarea, dan periode tanpa estimasi
  langsung (`y = NA`) tetap disertakan dalam model dan diprediksi
  beserta MSE-nya.
- 💾 **Hemat memori**: alokasi $`O(m^2)`$ / $`O(m^3)`$ yang tidak perlu
  dihindari.
- 📊 **Antarmuka R yang standar**: metode S3 (`summary`, `coef`,
  `fitted`, `residuals`) dan
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  yang grafiknya siap dipublikasikan.

------------------------------------------------------------------------

## Mengenal INLA: Inferensi Bayesian yang Cepat tanpa MCMC

**INLA** (*Integrated Nested Laplace Approximations*; Rue dkk., 2009)
adalah alternatif deterministik bagi *Markov Chain Monte Carlo* (MCMC)
untuk inferensi Bayesian pada model Gaussian laten. Alih-alih menarik
ribuan sampel, INLA:

- **mengaproksimasi** distribusi marginal posterior secara analitik
  dengan aproksimasi Laplace;
- **mengintegrasikan** observasi non-Gaussian secara analitik;
- **memberikan** hasil dalam hitungan detik, bukan jam.

`fastsae` memakai INLA melalui paket R
[`INLA`](https://www.r-inla-download.org) untuk mengestimasi model area
kecil Bayesian hierarkis dengan:

- **6 keluarga likelihood**: Gaussian, Beta, Binomial, Poisson, Binomial
  Negatif, Gamma;
- **efek acak spasial**: Besag ICAR, BYM, BYM2, Leroux, SAR, *spatial
  lag*;
- **efek acak temporal**: RW1, RW2, AR(1), iid;
- **interaksi spatio-temporal**: tren khusus per domain, medan
  ruang-waktu terpisah (*separable*), dan Knorr-Held tipe I–IV.

**Dibandingkan dengan MCMC**: tidak ada rantai yang harus dijalankan dan
dipantau, hasilnya deterministik, dan estimasinya selesai dalam hitungan
detik. INLA tetaplah sebuah aproksimasi. Karena itu, pengaturan INLA di
fastsae dipilih dengan membandingkannya terhadap posterior eksak dan
hasil Stan yang dijalankan panjang. Bila perlu pemeriksaan,
`engine = jax_engine()` mengestimasi model yang sama dengan MCMC.

------------------------------------------------------------------------

## Model Subarea Dua Tingkat (*Two-Fold*): Struktur Hierarkis Bersarang

Model subarea dua tingkat (Torabi & Rao, 2014) memperluas model
Fay-Herriot standar untuk data hierarkis yang **subareanya tersarang di
dalam area**:

``` math
y_{ij} = x_{ij}^\top \beta + v_i + u_{ij} + e_{ij}
```

dengan:

- $`v_i \sim N(0, \sigma_v^2)`$: **efek acak tingkat area**;
- $`u_{ij} \sim N(0, \sigma_u^2)`$: **efek acak tingkat subarea**;
- $`e_{ij} \sim N(0, D_{ij})`$: galat sampling dengan varians yang
  diketahui.

**Implementasi yang tersedia:**

- [`eblup_twofold()`](https://ridsonap.is-a.dev/fastsae/reference/eblup_twofold.md):
  EBLUP frekuentis melalui *Fisher scoring* REML atau ML (C++), dengan
  MSE analitik atau bootstrap serta rata-rata tingkat area;
- [`hb_twofold()`](https://ridsonap.is-a.dev/fastsae/reference/hb_twofold.md):
  model Bayesian hierarkis (INLA atau MCMC) untuk respons Gaussian,
  binomial, dan Poisson, dengan efek area spasial sebagai pilihan.

------------------------------------------------------------------------

## Perbandingan dengan Paket Lain

| Fitur | `fastsae` | `sae` (Molina & Rao) | `emdi` (Kreutzmann dkk.) |
|:---|:--:|:--:|:--:|
| **Komputasi inti** | **C++ (RcppArmadillo)** | R murni | R / lme4 / nlme |
| ***Multi-threading*** | **OpenMP bawaan** (`n_threads`) | Satu *thread* | Opsional lewat foreach/parallel |
| **Kesesuaian dengan `sae` pada data yang sama** | **Sekitar $`10^{-14}`$** | Rujukan | Sekitar $`10^{-6}`$ (toleransi konvergensi) |
| **Dukungan area tak tersampel** | **Otomatis (sintetis; *kriging* untuk model spasial)** | Data harus dipilah manual | Sintetis (domain di luar sampel) |
| **MSE bootstrap** | **C++ paralel (OpenMP)** | R | R |
| **Dukungan metode S3** | `print`, `summary`, `coef`, `fitted`, `residuals`, `autoplot` | List khusus | S3 standar |

------------------------------------------------------------------------

## Model yang Didukung

### Model Frekuentis (EBLUP, C++)

| Model | Fungsi | Struktur efek acak | Estimasi MSE |
|:---|:---|:---|:---|
| **Fay-Herriot** (tingkat area) | [`eblup_fh()`](https://ridsonap.is-a.dev/fastsae/reference/eblup_fh.md) | Efek area saling bebas (transformasi log sebagai pilihan) | Analitik (Prasad-Rao / Datta-Lahiri) |
| **Fay-Herriot Spasial** | [`eblup_sfh()`](https://ridsonap.is-a.dev/fastsae/reference/eblup_sfh.md) | Autokorelasi spasial SAR(1) | Analitik, bootstrap parametrik dan nonparametrik |
| **Fay-Herriot Spatio-Temporal** | [`eblup_stfh()`](https://ridsonap.is-a.dev/fastsae/reference/eblup_stfh.md) | SAR(1) + AR(1) temporal | Analitik, bootstrap parametrik |
| **Battese-Harter-Fuller** (tingkat unit) | [`eblup_bhf()`](https://ridsonap.is-a.dev/fastsae/reference/eblup_bhf.md) | *Intercept* acak yang tersarang dalam domain | Bootstrap |
| **Subarea Dua Tingkat** | [`eblup_twofold()`](https://ridsonap.is-a.dev/fastsae/reference/eblup_twofold.md) | Efek area + subarea yang tersarang | Analitik, bootstrap parametrik |

### Model Bayesian (INLA atau MCMC)

| Model | Fungsi | Distribusi | Spasial/Temporal |
|:---|:---|:---|:---|
| **HB tingkat area** | [`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md) | Gaussian, Beta, Binomial, Poisson, NB, Gamma | ✓ Dukungan penuh |
| **HB tingkat unit** | [`hb_unit()`](https://ridsonap.is-a.dev/fastsae/reference/hb_unit.md) | Gaussian, Binomial, Poisson | ✓ BYM2/Besag |
| **HB dua tingkat** | [`hb_twofold()`](https://ridsonap.is-a.dev/fastsae/reference/hb_twofold.md) | Gaussian, Binomial, Poisson | ✓ BYM2/Besag |
| **Subindikator dari indikator induk** | [`hb_ratio()`](https://ridsonap.is-a.dev/fastsae/reference/hb_ratio.md) | Pangsa binomial dari indikator induk, tidak pernah melebihinya | Melalui [`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md) |

Model Bayesian diestimasi dengan `engine = inla_engine()` (bawaan).
Keempatnya juga dapat diestimasi dengan MCMC di Python (NumPyro di atas
JAX, di CPU atau GPU) melalui `engine = jax_engine()`; hasilnya
melaporkan R-hat, ukuran sampel efektif, transisi divergen, dan E-BFMI.
[`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
juga menerima kovariat yang mengandung galat pengukuran (`var_x`). Kedua
*engine* memakai prior bawaan yang sama, yang dijelaskan di [artikel
prior
bawaan](https://ridsonap.is-a.dev/fastsae/articles/id/hb-priors.md).

------------------------------------------------------------------------

## Perangkat di Sekitar Model

| Tahap | Fungsi |
|:---|:---|
| Sebelum estimasi | [`create_weights()`](https://ridsonap.is-a.dev/fastsae/reference/create_weights.md) (matriks pembobot spasial), [`spatial_test()`](https://ridsonap.is-a.dev/fastsae/reference/spatial_test.md) (Moran’s I, uji LM), [`gvf_smooth()`](https://ridsonap.is-a.dev/fastsae/reference/gvf_smooth.md) (penghalusan varians sampling), [`effective_trials()`](https://ridsonap.is-a.dev/fastsae/reference/effective_trials.md) (ukuran sampel efektif untuk model binomial), [`step_sae()`](https://ridsonap.is-a.dev/fastsae/reference/step_sae.md) (pemilihan kovariat) |
| Memeriksa hasil estimasi | [`diagnose()`](https://ridsonap.is-a.dev/fastsae/reference/diagnose.md) (peningkatan presisi, kecocokan model, kalibrasi, autokorelasi spasial residual, PIT dan CPO), [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html), [`compare_sae()`](https://ridsonap.is-a.dev/fastsae/reference/compare_sae.md) (kesesuaian antara dua model) |
| Benchmarking | [`benchmark_sae()`](https://ridsonap.is-a.dev/fastsae/reference/benchmark_sae.md) setelah estimasi (rasio, selisih, optimal, logit), atau `self_benchmark = TRUE` di dalam fungsi model |
| Pelaporan | [`add_reliability_flags()`](https://ridsonap.is-a.dev/fastsae/reference/add_reliability_flags.md) (tingkat keandalan RSE, atau kelas beserta peluangnya), [`map_sae()`](https://ridsonap.is-a.dev/fastsae/reference/map_sae.md) (peta koroplet), [`export_sae()`](https://ridsonap.is-a.dev/fastsae/reference/export_sae.md) (lembar kerja) |
| Simulasi | [`sim_area_data()`](https://ridsonap.is-a.dev/fastsae/reference/sim_area_data.md), [`sim_series_data()`](https://ridsonap.is-a.dev/fastsae/reference/sim_series_data.md) |

------------------------------------------------------------------------

## Ringkasan Kinerja

#### 1. EBLUP frekuentis dengan MSE analitik (waktu dan memori puncak)

| Model ($`n = 1.000`$ domain, 5 kovariat) | `fastsae` | `sae` | `emdi` |
|:---|:--:|:--:|:--:|
| **Fay-Herriot** | **0,000639 detik** (1,91 MB) | 0,799 detik (39,0 MB; 1.250x lebih lambat) | 13,5 detik (222 MB; 21.071x lebih lambat) |
| **Fay-Herriot Spasial** | **0,518 detik** (78,7 MB) | 24,1 detik (239 MB; 46x lebih lambat) | 25,2 detik (310 MB; 49x lebih lambat) |
| **FH Spatio-Temporal ($`T = 5`$), estimasi titik** | **39,8 detik** (495 MB) | 1160 detik (3,60 GB; 29x lebih lambat) | \- |

#### 2. Model Beta Bayesian: INLA vs MCMC Stan ($`n`$ terbesar yang diuji)

| Model | $`n`$ | [`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md) (INLA) | [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html) (Stan) | Rasio waktu |
|:---|:--:|:--:|:--:|:--:|
| **Beta** | 1.000 | **3,73 detik** | 180 detik | 48x |
| **Beta Spasial (Besag ICAR)** | 1.000 | **4,28 detik** | 362 detik | 85x |
| **Beta Spatio-Temporal (Besag + RW1, $`T = 5`$)** | 50 | **3,80 detik** | 716 detik | 189x |

Waktu yang dilaporkan adalah median dari `bench::mark()`. Baris
frekuentis diukur pada 2026-10-26 dengan fastsae 1.0.0, sae 1.3, dan
emdi 2.2.3 pada himpunan data simulasi yang sama (Apple M4, R version
4.6.1 (2026-06-24); `benchmarks/bench_frequentist.R`); setiap sel
dijalankan dalam proses R tersendiri. Memori adalah kenaikan memori
residen puncak proses tersebut selama estimasi, termasuk memori yang
dialokasikan oleh kode terkompilasi (kenaikan beberapa puluh MB
berubah-ubah antarproses sehingga tidak layak dibandingkan). `sae` tidak
memiliki MSE analitik untuk model spatio-temporal, jadi baris itu
membandingkan estimasi titik saja; dengan MSE analitiknya, fastsae
memerlukan 39,3 detik. Baris `fastsae` Bayesian diukur ulang pada
2026-10-14 (`benchmarks/bench_bayes_fastsae.R`), sedangkan baris
`tipsae` (4 rantai, masing-masing 2.000 iterasi) berasal dari pengukuran
sebelumnya pada data yang sama. Hasil untuk semua ukuran ada di [artikel
tolok
ukur](https://ridsonap.is-a.dev/fastsae/articles/id/benchmarks.md).

------------------------------------------------------------------------

## Instalasi

CRAN baru memuat fastsae 0.1.0, versi lama yang belum memiliki model
Bayesian hierarkis dan sebagian besar perangkat di atas. Versi 1.0.0
yang dijelaskan di situs ini dipasang dari GitHub:

``` r

# install.packages("remotes")
remotes::install_github("ridsonap/fastsae")
```

Model Bayesian memerlukan **INLA**, yang tidak tersedia di CRAN:

``` r

install.packages(
  "INLA",
  repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
  dep = TRUE
)
```

*Engine* MCMC yang bersifat opsional
([`jax_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md))
memerlukan Python dengan JAX dan NumPyro:
[`setup_numpyro_env()`](https://ridsonap.is-a.dev/fastsae/reference/setup_numpyro_env.md)
memasangnya melalui `reticulate`, dan
`setup_numpyro_env(check_only = TRUE)` memeriksa pemasangannya.

------------------------------------------------------------------------

## Dokumentasi dan Panduan Pengguna

Panduan rinci, contoh lengkap, dan analisis kinerja:

- [**Memulai
  fastsae**](https://ridsonap.is-a.dev/fastsae/articles/id/getting-started.md):
  langkah demi langkah mengestimasi model, mendiagnosisnya, dan
  mengambil hasil dengan metode S3.
- [**Alur Kerja SAE yang
  Optimal**](https://ridsonap.is-a.dev/fastsae/articles/id/optimal-workflow.md):
  urutan analisis yang lengkap, dari keandalan estimasi langsung hingga
  benchmarking dan pelaporan.
- [**Model Spasial & Spatio-Temporal
  (EBLUP)**](https://ridsonap.is-a.dev/fastsae/articles/id/spatial-temporal.md):
  bekerja dengan matriks pembobot spasial, *kriging* spasial, dan data
  panel spatio-temporal.
- [**Model Bayesian Tingkat Area
  (INLA)**](https://ridsonap.is-a.dev/fastsae/articles/id/hb-area-inla.md):
  SAE Bayesian hierarkis untuk 6 distribusi serta efek spasial dan
  temporal, pemeriksaan konvergensi MCMC, kovariat dengan galat
  pengukuran, dan
  [`hb_ratio()`](https://ridsonap.is-a.dev/fastsae/reference/hb_ratio.md).
- [**Prior Bawaan Model
  HB**](https://ridsonap.is-a.dev/fastsae/articles/id/hb-priors.md):
  prior yang dipakai secara bawaan dan alasan pemilihannya.
- [**Diagnostik & Kalibrasi
  Model**](https://ridsonap.is-a.dev/fastsae/articles/id/model-diagnostics.md):
  diagnostik model dan benchmarking ke target agregat.
- [**Estimasi Tingkat Unit (BHF dan
  HB)**](https://ridsonap.is-a.dev/fastsae/articles/id/unit-level-bhf.md):
  mengestimasi model Battese-Harter-Fuller dengan data survei tingkat
  individu.
- [**Model Subarea Dua
  Tingkat**](https://ridsonap.is-a.dev/fastsae/articles/id/twofold-subarea.md):
  estimasi subarea yang tersarang di dalam area.
- [**Tolok Ukur
  Kinerja**](https://ridsonap.is-a.dev/fastsae/articles/id/benchmarks.md):
  grafik interaktif dan perbandingan percepatan untuk berbagai ukuran
  sampel.
- [**Engine
  hb_area()**](https://ridsonap.is-a.dev/fastsae/articles/id/hb-area-engines.md):
  INLA dibandingkan dengan MCMC JAX di CPU dan GPU, serta dengan tipsae.
- [**Validasi Populasi
  Terhingga**](https://ridsonap.is-a.dev/fastsae/articles/id/finite-population-validation.md):
  simulasi berbasis desain terhadap nilai populasi yang diketahui.
- [**Apakah fastsae Lebih
  Baik?**](https://ridsonap.is-a.dev/fastsae/articles/id/fastsae-evaluation.md):
  akurasi, kecepatan, dan memori dibandingkan sae, emdi, dan tipsae.
- [**Daftar
  Istilah**](https://ridsonap.is-a.dev/fastsae/articles/id/glosarium.md):
  padanan istilah Inggris–Indonesia yang dipakai di situs ini.
- [**Referensi
  Fungsi**](https://ridsonap.is-a.dev/fastsae/reference/index.md):
  dokumentasi lengkap API dan argumennya (dalam bahasa Inggris).

------------------------------------------------------------------------

## Daftar Pustaka

- Fay, R. E., & Herriot, R. A. (1979). Estimates of income for small
  places: An application of James-Stein procedures to Census data.
  *Journal of the American Statistical Association*, 74(366), 269–277.
- Battese, G. E., Harter, R. M., & Fuller, W. A. (1988). An
  error-components model for prediction of county crop areas using
  survey and satellite data. *Journal of the American Statistical
  Association*, 83(401), 28–36.
- Marhuenda, Y., Molina, I., & Morales, D. (2013). Small area estimation
  with spatio-temporal Fay-Herriot models. *Computational Statistics &
  Data Analysis*, 58, 308–325.
- Pratesi, M., & Salvati, N. (2008). Small area estimation: The EBLUP
  estimator based on spatially correlated random area effects.
  *Statistical Methods & Applications*, 17(1), 113–141.
- Torabi, M., & Rao, J. N. K. (2014). On small area estimation under a
  sub-area level model. *Journal of Multivariate Analysis*, 127, 36–55.
- Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian
  inference for latent Gaussian models using integrated nested Laplace
  approximations (with discussion). *Journal of the Royal Statistical
  Society: Series B*, 71(2), 319–392.
- Rao, J. N. K., & Molina, I. (2015). *Small Area Estimation* (2nd ed.).
  John Wiley & Sons.
