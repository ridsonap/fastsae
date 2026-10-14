# CODE_MAP — Peta kode paket `fastsae`

Sumber: `NAMESPACE`, `R/*.R`, `src/*.cpp`, `src/Makevars`, `DESCRIPTION`.
Nomor baris mengacu pada keadaan repo per 14 Okt 2026 (`main` setelah `4a66caa`;
versi pertama: commit `0252055`). Dibuat dengan membaca kode, bukan dari
deskripsi luar.

---

## 1. Fungsi yang diekspor (`NAMESPACE`)

### 1.1 Model estimasi

| Fungsi ekspor | File R | Model | Metode estimasi | Metode MSE | Bagian C++ | OpenMP |
|---|---|---|---|---|---|---|
| `eblup_fh` | `R/eblup_fh.R:98` | Fay–Herriot area-level; `transform = "log"` (Log-FH) | Fisher scoring ML/REML (C++, `.eblup_core`) | Analitik $g_1+g_2+2g_3$ (REML) / dikurangi bias ML (C++); Log-FH: $\hat\theta^2\,\mathrm{MSE}_{\log}$ | `src/eblup_fh.cpp:8` | tidak |
| `eblup_bhf` | `R/eblup_bhf.R:62` | Battese–Harter–Fuller unit-level | REML/ML terprofil atas $\gamma=\sigma_u^2/\sigma_e^2$ (grid log + Brent, C++ `.bhf_fit_cpp`); prediktor `.eblup_bhf_cpp` | Bootstrap parametrik (C++ `.pbmse_bhf_cpp`) | `src/eblup_bhf_core.cpp:221,246`, `src/eblup_unit.cpp:8` | ya: `eblup_bhf_core.cpp:283` |
| `eblup_twofold` | `R/eblup_twofold.R:76` | Two-fold subarea | Fisher scoring 2 parameter $(\sigma_v^2,\sigma_u^2)$ (C++), informasi REML eksak | Analitik $g_1+g_2+2g_3$ (− bias ML) + bootstrap parametrik | `src/eblup_twofold.cpp:303` | ya: `eblup_twofold.cpp:558` |
| `eblup_sfh` | `R/eblup_sfh.R:131` | Spatial FH (SAR pada semua $D$ domain) | Fisher scoring $(\sigma^2,\rho)$ + fallback terkendala; faktorisasi sparse `SarEngine` | Analitik $g_1+g_2+2g_3-g_4$ (+ kriging orde kedua); `pbmse`; `npbmse` | `src/eblup_sfh_core.cpp:393`, `eblup_sfh_pbmse.cpp:15`, `eblup_sfh_npbmse.cpp:15`, `sar_engine.h` | ya: `eblup_sfh_pbmse.cpp:116`, `eblup_sfh_npbmse.cpp:188` |
| `eblup_stfh` | `R/eblup_stfh.R:164` | Spatio-temporal FH (efek spasial konstan-waktu + AR(1)) pada grid $D\times T$; sel `y = NA` (periode hilang, domain tanpa sampel) diprediksi dengan BLUP | REML 3/4 parameter terstruktur pada sel tersampel (Woodbury, blok per domain, tanpa matriks $n\times n$) + fallback terkendala | Tersampel: analitik $g_1+g_2+2g_3-g_4$; tanpa sampel: $g_1+g_2+g_3$; bootstrap parametrik untuk semua sel | `src/eblup_stfh_core.cpp:826`, `eblup_stfh_pbmse.cpp:62` | ya: `eblup_stfh_pbmse.cpp:193` |
| `hb_area` | `R/hb_area.R:231` | HB area-level (6 family; iid/BYM/BYM2/Besag/generic1(Leroux)/SLM + temporal + KH 1–4, separable) | `inla_engine()` (INLA) atau `jax_engine()` (NumPyro NUTS, `R/hb_area_jax.R`) | $\mathrm{sd}^2$ posterior; self-benchmarking posterior via `inla.posterior.sample` | — | — |
| `hb_unit` | `R/hb_unit.R:123` | HB unit-level (gaussian/binomial/poisson; iid/BYM2/Besag) | INLA, prior `loggamma(0.01,0.01)`; Gaussian non-spasial: start mode REML | $\mathrm{sd}^2$ posterior; koreksi populasi hingga | — | — |
| `hb_twofold` | `R/hb_twofold.R:136` | HB two-fold subarea | INLA | $\mathrm{sd}^2$; tingkat area via sampling posterior | — | — |
| `inla_engine`, `jax_engine` | `R/engine.R:83,155` | Objek engine untuk `hb_*` | — | — | — | — |
| `setup_numpyro_env`, `check_numpyro_available` | `R/python_env.R:81,16` | Lingkungan Python untuk `jax_engine` | — | — | — | — |

### 1.2 Utilitas

| Fungsi ekspor | File | Fungsi | Metode |
|---|---|---|---|
| `benchmark_sae` | `R/benchmark.R:120` | Kalibrasi agregat | Ratio, difference, optimal (KKT, $\delta_i\propto w_i\mathrm{mse}_i$), logit (`uniroot`); target nasional dua tahap; `mse_benchmarked` dengan matriks MSE penuh (`.mse_operator`) |
| `diagnose` | `R/diagnose.R:75` | Diagnostik model & indikator keandalan | Uji Brown (regresi bias deskriptif; $W$ vs $\chi^2_m$ ekor atas), residual terstandar $(y-\hat\theta)/\sqrt{\psi-\mathrm{MSE}}$, Moran's I + varian Cliff–Ord, flag RSE |
| `compare_sae` | `R/compare_sae.R:62` | Perbandingan dua model | Pearson/Spearman, MAE, RMSD, rasio MSE, gain RSE |
| `export_sae`, `add_reliability_flags` | `R/export_sae.R:100,35` | Ekspor tabel; flag RSE | CSV/XLSX; ambang BPS `c(25, 50)` |
| `map_sae` | `R/map_sae.R:95` | Peta choropleth | ggplot2/sf |
| `create_weights` | `R/create_weights.R:69` | Matriks bobot spasial | queen/rook/knn/distance/idw; gaya W/B/C/U/minmax; fallback pulau |
| `spatial_test` | `R/spatial_test.R:62` | Uji pra-pemodelan spasial | Moran's I (Cliff–Ord), LM-error/LM-lag (Anselin), LRT ML untuk $\rho=0$; satu uji keputusan (LRT bila `vardir` ada, selain itu LM-error) |
| `gvf_smooth` | `R/gvf.R:67` | Penghalusan `vardir` | GVF log-linear, Gamma GLM, CV power |
| `autoplot`, `plot`, `print`, `summary`, `coef`, `fitted`, `residuals` | `R/methods.R`, `R/autoplot.R` | Metode S3 | — |
| `sim_spatial_weights`, `sim_area_data`, `sim_series_data` | `R/sim_*.R:47,47,66` | Data simulasi | DGP area-level dan STFH |

Data: `mys`, `mys_panel`, `mys_proxmat`, `cornsoybean`, `cornsoybeanmeans` (`R/data.R`, `LazyData: true`).

---

## 2. Lapisan C++ (`src/`)

| Berkas | Fungsi Rcpp | Peran | OpenMP |
|---|---|---|---|
| `eblup_fh.cpp` | `.eblup_core` | Fisher scoring + EBLUP + MSE FH | – |
| `eblup_unit.cpp` | `.eblup_bhf_cpp` | Prediktor area BHF ($f_d$) | – |
| `eblup_bhf_core.cpp` | `.bhf_fit_cpp`, `.pbmse_bhf_cpp` | REML/ML BHF terprofil + bootstrap | line 283 |
| `eblup_twofold.cpp` | `.eblup_twofold_core` | Fit two-fold + MSE analitik + bootstrap | line 558 |
| `eblup_sfh_core.cpp` + `sar_engine.h` | `.seblup_core` | Fisher scoring spatial FH (sparse) + MSE analitik + kriging | – |
| `eblup_sfh_pbmse.cpp` | `.seblup_pbmse` | Bootstrap parametrik spatial FH | line 116 |
| `eblup_sfh_npbmse.cpp` | `.seblup_npbmse` | Bootstrap nonparametrik spatial FH | line 188 |
| `eblup_stfh_core.cpp` | `.eblup_stfh_core` | Fit STFH terstruktur + MSE analitik | – |
| `eblup_stfh_pbmse.cpp` | `.pbmse_stfh` | Bootstrap parametrik STFH (estimasi ulang $\theta$; semua sel grid) | line 193 |
| `RcppExports.cpp` / `R/RcppExports.R` | stub registrasi | useDynLib(registration=TRUE) | – |

`src/Makevars`: `CXX_STD = CXX17`, OpenMP lewat `SHLIB_OPENMP_CXXFLAGS`,
LAPACK/BLAS; `LinkingTo`: `Rcpp`, `RcppArmadillo`, `RcppEigen`.

## 2b. Lapisan Python (`inst/python/`, engine JAX)

| Berkas / fungsi | Peran |
|---|---|
| `hb_area_numpyro.py:76` `init_jax_environment` | Memilih perangkat JAX (`cpu`/`cuda`/`metal`) lewat `JAX_PLATFORMS`, mengunci `jax_default_device`; memanggil `_patch_numpyro_for_metal` bila perangkatnya Metal (baris 107) |
| `hb_area_numpyro.py:25,36` `_popcount32`, `_leaf_idx_to_ckpt_idxs` | Pengganti `jnp.bitwise_count` (`popcnt`, tidak ada di jax-metal) untuk indeks checkpoint U-turn NUTS; identik dengan NumPyro untuk 0–4095 (tes `test_hb_area_engine.R`) |
| `hb_area_numpyro.py:42` `_patch_numpyro_for_metal` | Menambal NumPyro di sesi Metal: indeks checkpoint di atas, dan `lax.cond` yang keluarannya dilewatkan `lax.optimization_barrier` (jax-metal gagal menserialisasi `cond` yang semua cabangnya mengeluarkan konstanta) |
| `hb_area_numpyro.py:452-508` | Di Metal: rantai sekuensial dan dijalankan per langkah dari Python (jalur progress-bar NumPyro; loop rantai yang digabung tidak pernah selesai di Metal), bar disenyapkan bila `progress = FALSE`; fallback ke CPU bila kompilasi gagal |
| `R/hb_area_jax.R:27` `.icar_spectral` | Basis eigen $D\times D$ Laplacian graf untuk efek ICAR/BYM2 di JAX: $O(D^3)$ sekali, $O(D^2)$ per evaluasi gradien |
| `R/python_env.R:145` `.warn_device_fallback` | Peringatan bila GPU yang diminta tidak dipakai |

Benchmark engine: `benchmarks/bench_hb_area_engines.R` (INLA, JAX CPU/Metal/CUDA, `tipsae`;
satu proses per engine; hasil `inst/extdata/hb_engine_benchmark*.rds`), notebook CUDA di
`benchmarks/kaggle/` dan `benchmarks/colab/`.

Ketergantungan (`Imports`): `cli`, `ggplot2`, `Matrix`, `methods`, `Rcpp`,
`rlang`, `stats`, `utils`. `INLA`, `reticulate`, `sf`, `spdep`, `sn`, `lme4`
(hanya pembanding di tes), `sae`, `emdi`, `tipsae`, dll. di `Suggests`.

---

## 3. Pokok per bagian laporan

- **CB1** Pendahuluan — SAE, direct vs model-based, arsitektur.
- **CB2** Desain perangkat lunak — inti C++, Fisher scoring, blok/Woodbury, memori, bootstrap paralel.
- **CB3** Benchmark — metodologi + hasil dari `inst/extdata`.
- **CB4** Latar statistik — notasi, LMM, GLS/BLUP/EBLUP, ML/REML, Fisher scoring, $g_1g_2g_3$, bootstrap, unsampled, benchmarking.
- **CB5** Per model — FH, SFH, STFH, BHF, two-fold, HB area/unit/two-fold + utilitas.
- **CB6** Bab Bayes — INLA, GMRF, CAR/ICAR/proper CAR/SAR, BYM/BYM2, PC prior, likelihood/link, posteros, TikZ.
- **CB7** Diagnostik, perbandingan, ekspor, pemetaan.
- **CB8** Lampiran — derivasi, glosarium, reproduktibilitas.
