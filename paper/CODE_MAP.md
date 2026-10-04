# CODE_MAP — Peta kode paket `fastsae`

Sumber: `NAMESPACE`, `R/*.R`, `src/*.cpp`, `src/Makevars`, `DESCRIPTION`.
Semua nomor baris mengacu pada keadaan repo `main` (commit `0252055`).
Dibuat dengan membaca kode, bukan dari deskripsi luar.

---

## 1. Fungsi yang diekspor (`NAMESPACE`)

### 1.1 Model estimasi

| Fungsi ekspor | File R | Model | Metode estimasi | Metode MSE | Bagian C++ | OpenMP |
|---|---|---|---|---|---|---|
| `eblup_fh` | `R/eblup_fh.R:49` | Fay–Herriot area-level, $y_d=x_d'\beta+u_d+e_d$ | Fisher scoring ML/REML (C++, `.eblup_core`) | Analitik $g_1+g_2+2g_3$ (REML) / $g_1+g_2-g_3$-bentuk ML (C++) | `src/eblup_fh.cpp:9` (`Rcpp::export(.eblup_core)`) | tidak |
| `eblup_bhf` | `R/eblup_bhf.R:44` | Battese–Harter–Fuller unit-level | `lme4::lmer` (ML/REML) untuk $\beta,\sigma^2$; prediktor area di C++ (`.eblup_bhf_cpp`) | Analitik `VarPred` (R) + bootstrap parametrik `.pbmse_unit` (R, serial) | `src/eblup_unit.cpp:91` | tidak |
| `eblup_twofold` | `R/eblup_twofold.R:64` | Two-fold subarea: $y_{ijk}=x_{ij}'\beta+u_i+v_{ij}+e_{ijk}$ | Fisher scoring 2 parameter $(\sigma_v^2,\sigma_u^2)$ di R, fit inti di C++ | Analitik $g_1{+}g_2{+}g_3{+}g_4$ + bootstrap parametrik di C++ | `src/eblup_twofold.cpp:273` (`.eblup_twofold_core`) | ya: `src/eblup_twofold.cpp:481` (`parallel for schedule(static)`, replikat bootstrap) |
| `eblup_sfh` | `R/eblup_sfh.R:85` | Spatial FH (error-components/SAR, $A=(I-\rho W')(I-\rho W)$) | Fisher scoring 2 parameter $(\sigma^2,\rho)$ di C++ | Analitik $g_1{+}g_2{+}2g_3{-}g_4$; `pbmse`; `npbmse` | `src/eblup_sfh_core.cpp:209` (`.seblup_core`), `eblup_sfh_pbmse.cpp:15` (`.seblup_pbmse`), `eblup_sfh_npbmse.cpp:15` (`.seblup_npbmse`) | ya: `eblup_sfh_pbmse.cpp:108`, `eblup_sfh_npbmse.cpp:159` (`parallel for schedule(dynamic)`) |
| `eblup_stfh` | `R/eblup_stfh.R:103` | Spatio-temporal FH (efek spasial konstan-waktu + AR(1) waktu per domain) | Fisher scoring/REML 4 parameter $(\sigma_1^2,\rho_1,\sigma_2^2,\rho_2)$ di C++ | Analitik + bootstrap parametrik | `src/eblup_stfh_core.cpp:160` (`.eblup_stfh_core`), `eblup_stfh_pbmse.cpp:232` (`.pbmse_stfh`) | ya: `eblup_stfh_pbmse.cpp:352` (`parallel for schedule(dynamic)`) |
| `hb_area` | `R/hb_area.R:148` | Hierarki Bayes area-level (Gaussian, beta, binomial, poisson, nbinomial, gamma; struktur iid/BYM/BYM2/BESAG/GENERIC1/SLM + temporal + interaksi space-time KH type 1–4, separable) | INLA (nested Laplace); fallback `glmm` → `.fit_glmm_laplace` | Rumus posterior: $\mathrm{sd}^2$ / logit-scale (frequentis), eksponensial-momen (overdispersi), $1/\mathrm{sd}^2$ (normal) | — (hanya memanggil INLA) | — |
| `hb_unit` | `R/hb_unit.R:169` | Hierarki Bayes unit-level (linear normal; $u_i$ via `f(iid)` atau `f()` on-the-fly CAR) | INLA (`strategy` MCMC default), prior `loggamma(0.01,0.01)` | $\mathrm{sd}^2$ analitik dari `summary.fitted.values` | — | — |
| `hb_twofold` | `R/hb_twofold.R:138` | Hierarki Bayes two-fold subarea | INLA; interaksi spacetime KH type 1–4; `f(zz, model=slm)` | $\mathrm{sd}^2$ analitik | — | — |

### 1.2 Utilitas

| Fungsi ekspor | File | Fungsi | Metode |
|---|---|---|---|
| `benchmark_sae` / `benchmark` | `R/benchmark.R:120,264` | Kalibrasi agregat | Ratio, difference, optimal (pemberat MSE), logit (`uniroot`); target nasional + hierarki bertingkat `outer_target` |
| `diagnose` | `R/diagnose.R:...` | Diagnostik model & indikator keandalan | Uji Brown (F, df (2,n−2)), uji kecocokan $W=\sum(y-\hat y)^2/(v_d+\mathrm{mse})\sim\chi^2$, Moran's I + varian Cliff–Ord, flag RSE |
| `compare_sae` | `R/compare_sae.R` | Perbandingan dua model | Pearson/Spearman, MAE, RMSD, rasio MSE, gain RSE |
| `export_sae` | `R/export_sae.R` | Ekspor tabel hasil ke spreadsheet | Tulis CSV/XLSX |
| `map_sae` | `R/map_sae.R` | Peta choropleth (S3: `fastsae`, `fastsae_benchmark`, `list`, `default`) | ggplot2/sf |
| `autoplot`, `plot`, `print`, `summary`, `coef`, `fitted`, `residuals` | `R/methods.R`, `R/autoplot.R` | Metode S3 untuk kelas `fastsae`, `fastsae_diagnose`, `fastsae_benchmark`, `fastsae_comparison` | — |
| `sim_spatial_weights` | `R/sim_spatial_weights.R` | Matriks tetangga rook/queen + link skala untuk grid | Deterministik (indeks `expand.grid`) |
| `sim_area_data` | `R/sim_area_data.R` | Data simulasi area-level | DGP Fay–Herriot (deterministik dgn `seed`) |
| `sim_series_data` | `R/sim_series_data.R` | Data simulasi deret spasio-temporal | DGP STFH |
| `add_reliability_flags` | `R/diagnose.R` (bagian) | Flag keandalan berbasis RSE | Threshold default 25% |

Data: `mys`, `mys_panel`, `mys_proxmat`, `cornsoybean`, `cornsoybeanmeans` (`R/data.R`, `LazyData: true`).

---

## 2. Lapisan C++ (`src/`)

| Berkas | Fungsi Rcpp | Peran | OpenMP |
|---|---|---|---|
| `eblup_fh.cpp` | `.eblup_core` | Fisher scoring + EBLUP + MSE FH | – |
| `eblup_unit.cpp` | `.eblup_bhf_cpp` | Prediktor area BHF (rata-rata sampel & populasi, $f_d$) | – |
| `eblup_twofold.cpp` | `.eblup_twofold_core` | Fit two-fold + MSE analitik + bootstrap | line 481 |
| `eblup_sfh_core.cpp` | `.seblup_core` | Fisher scoring spatial FH + MSE analitik | – |
| `eblup_sfh_pbmse.cpp` | `.seblup_pbmse` | Bootstrap parametrik spatial FH | line 108 |
| `eblup_sfh_npbmse.cpp` | `.seblup_npbmse` | Bootstrap nonparametrik spatial FH | line 159 |
| `eblup_stfh_core.cpp` | `.eblup_stfh_core` | Fit STFH + MSE analitik | – |
| `eblup_stfh_pbmse.cpp` | `.pbmse_stfh` | Bootstrap parametrik STFH | line 352 |
| `RcppExports.cpp` / `R/RcppExports.R` | stub registrasi | useDynLib(registration=TRUE) | – |

`src/Makevars`: `CXX_STD = CXX17`, `PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS)`,
`PKG_LIBS = $(SHLIB_OPENMP_CXXFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)` —
Armado/LAPACK-BC + OpenMP; wajib lewat `R CMD SHLIB`/`Rcpp::compileAttributes()`.

Ketergantungan: `Rcpp`, `RcppArmadillo` (LinkingTo), `lme4` (BHF), `INLA` (Bayes),
`cli`, `ggplot2`, `scales`, `utils` (Suggests: `sf`, `testthat`, `bench`).

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
