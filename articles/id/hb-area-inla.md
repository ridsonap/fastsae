# Bayesian Area-Level Small Area Estimation with INLA

*Terjemahan dari [versi bahasa
Inggris](https://ridsonap.is-a.dev/fastsae/articles/hb-area-inla.md).
Penjelasan istilah ada di [Daftar
Istilah](https://ridsonap.is-a.dev/fastsae/articles/id/glosarium.md).*

## Pendahuluan

Model Small Area Estimation (SAE) klasik, seperti model Fay-Herriot
(1979), secara tradisional mengasumsikan distribusi Gaussian untuk
direct estimator survei dengan varians sampling error yang diketahui.
Namun, statistik resmi di dunia nyata sering berhadapan dengan indikator
non-Gaussian:

- **Angka kemiskinan, tingkat pengangguran, dan rasio prevalensi**:
  terbatas ketat pada interval terbuka $`(0, 1)`$.
- **Jumlah kasus penyakit, kejadian kriminal, dan peristiwa langka**:
  data cacah diskret yang mengalami overdispersi Poisson atau Negative
  Binomial dengan eksposur populasi.
- **Agregat survei biner**: proporsi yang berasal dari ukuran sampel
  kecil (percobaan Binomial).
- **Pengeluaran ekonomi yang menceng**: hasil kontinu yang positif ketat
  (Gamma).

Selain itu, bila data dikumpulkan di berbagai area geografis dan diulang
dari waktu ke waktu (survei panel/longitudinal), hasilnya menunjukkan
ketergantungan spasial sekaligus autokorelasi temporal.

Paket **fastsae** menghadirkan
[`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md),
kerangka hierarchical Bayesian terpadu yang ditenagai oleh **Integrated
Nested Laplace Approximations (INLA)**. Metode ini memberi aproksimasi
deterministik yang cepat untuk marginal posterior, tanpa chain, warm-up,
dan tuning seperti pada Markov chain Monte Carlo (MCMC); dibandingkan
dengan eksekusi Stan yang panjang dengan model dan prior yang sama,
posterior mean dan SD-nya sangat sesuai (lihat bagian Engines pada
[`?hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)).

------------------------------------------------------------------------

## Family Distribusi yang Didukung

[`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
mendukung 6 distribusi probabilitas untuk respons langsung:

| Family | Kasus Penggunaan SAE Tipikal | Input Sampling yang Diperlukan | Parameter Skala Model |
|:---|:---|:---|:---|
| **`"gaussian"`** | Indikator kontinu (mis. rata-rata pendapatan, rata-rata lama sekolah) | `vardir` (sampling variance $`D_d`$) | Kebalikan sampling variance $`\text{scale} = 1/D_d`$ |
| **`"beta"`** | Rate dan proporsi yang ketat di $`(0, 1)`$ | `vardir` atau `trials` | Presisi $`\phi_d = \hat\mu_d(1 - \hat\mu_d)/V_d - 1`$ dari fit pertama (Janicki, 2020), atau $`n_d - 1`$ dari `trials` |
| **`"poisson"`** | Data cacah dan rate peristiwa langka | `exposure` (expected count / offset populasi) | Offset log-link $`\log(E_d)`$ |
| **`"nbinomial"`** | Data cacah overdispersed dan peristiwa berkelompok | `exposure`; `vardir` opsional | Offset log-link; size $`\hat\mu_d^2/(V_d - \hat\mu_d)`$ dari fit pertama bila `vardir` diberikan |
| **`"binomial"`** | Proporsi biner teragregasi | `trials` (ukuran sampel domain $`n_d`$) | Jumlah trial per domain |
| **`"gamma"`** | Indikator positif menceng ke kanan (mis. pengeluaran) | `vardir` | Presisi (shape) $`y_d^2/V_d`$ |

------------------------------------------------------------------------

## Formulasi Matematis

Model hierarchical Bayesian yang mendasari
[`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
memperluas kerangka Fay-Herriot klasik untuk menangani respons
non-Gaussian dan random effect terstruktur. Model ini terdiri dari tiga
level:

### Level 1: Sampling Model (Likelihood)

Untuk domain $`d`$ dengan direct estimator $`y_d`$ dan sampling variance
$`D_d`$ yang diketahui:

``` math
\pi(y_d \mid \eta_d) = \text{ExpFamily}(y_d \mid \eta_d, D_d)
```

Fungsi link menghubungkan linear predictor $`\eta_d`$ dengan rata-rata
$`\mu_d = E(y_d)`$:

| Family | Fungsi Link | Linear Predictor |
|:---|:---|:---|
| Gaussian | Identitas: $`\eta_d = \mu_d`$ | $`\eta_d = x_d^\top \beta + u_d`$ |
| Beta | Logit: $`\eta_d = \log(\mu_d / (1-\mu_d))`$ | $`\eta_d = x_d^\top \beta + u_d`$ |
| Poisson | Log: $`\eta_d = \log(\mu_d)`$ | $`\eta_d = x_d^\top \beta + u_d + \log(E_d)`$ |
| Binomial | Logit: $`\eta_d = \log(\mu_d / (1-\mu_d))`$ | $`\eta_d = x_d^\top \beta + u_d`$ |
| Negative Binomial | Log: $`\eta_d = \log(\mu_d)`$ | $`\eta_d = x_d^\top \beta + u_d + \log(E_d)`$ |
| Gamma | Log: $`\eta_d = \log(\mu_d)`$ | $`\eta_d = x_d^\top \beta + u_d`$ |

### Level 2: Prior untuk Fixed Effect

``` math
\beta_0 \propto 1, \qquad \beta_j \sim N(0, 1000), \quad j = 1, \ldots, p
```

Intercept memakai prior flat INLA dan koefisien lainnya memakai prior
normal dengan varians 1000, yang hampir flat pada skala kebanyakan
covariate.

### Level 3: Prior untuk Random Effect

``` math
u \mid \sigma_u \sim N(0, \sigma_u^2 R^{-}), \qquad \sigma_u \sim \text{half-normal}(s)
```

Random effect $`u = (u_1, \ldots, u_m)^\top`$ mengikuti Gaussian Markov
random field berrata-rata nol yang matriks strukturnya $`R`$ mengodekan
ketergantungan spasial atau temporal yang dipilih ($`R = I`$ untuk efek
iid; Besag, BYM2, RW1, dan seterusnya untuk selainnya). Skala prior
half-normal pada standar deviasinya adalah $`s = \mathrm{sd}(y)`$, yaitu
standar deviasi direct estimate, untuk model Gaussian, dan $`s = 1`$
pada skala logit atau log untuk family lainnya. Hyperparameter lainnya
(parameter mixing BYM2, autokorelasi, parameter dispersi) beserta alasan
pemilihan default ini tercantum di
[`vignette("hb-priors")`](https://ridsonap.is-a.dev/fastsae/articles/hb-priors.md);
[`inla_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)
dapat menggantikan salah satunya.

### Aproksimasi INLA

Alih-alih sampling MCMC, INLA mengaproksimasi marginal posterior
hyperparameter $`\theta`$ dan latent field $`x = (\beta, u)`$ (Rue dkk.,
2009):

``` math
\tilde{\pi}(\theta \mid y) \propto \left. \frac{\pi(x, \theta, y)}{\tilde{\pi}_G(x \mid \theta, y)} \right|_{x = x^*(\theta)}
```

dengan $`\tilde{\pi}_G`$ adalah aproksimasi Gaussian dari full
conditional $`x`$ dan $`x^*(\theta)`$ adalah modusnya, serta

``` math
\tilde{\pi}(x_d \mid y) = \sum_k \tilde{\pi}(x_d \mid \theta_k, y) \, \tilde{\pi}(\theta_k \mid y) \, \Delta_k
```

dengan integral dievaluasi memakai **Laplace approximation** dan
**integrasi numerik**. Tiga strategi komputasi tersedia melalui
`hb_area(..., engine = inla_engine(strategy = ...))`, yang juga menerima
prior dan argumen
[`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html) lainnya:

- **`strategy = "gaussian"`**: Tercepat, terbaik untuk masalah yang
  berperilaku baik
- **`strategy = "laplace"`**: Default, keseimbangan akurasi/kecepatan
  yang baik
- **`strategy = "simplified.laplace"`**: Laplace yang sedikit lebih
  cepat, baik untuk model besar

Model yang sama dapat diestimasi dengan MCMC di CPU atau GPU melalui
`engine = jax_engine()` (NumPyro NUTS di atas JAX; lihat
[`?jax_engine`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)).
Hasil MCMC baru dapat dipakai bila chain-nya telah konvergen. Fungsi
melaporkan rank-normalized R-hat dan effective sample size (ESS) untuk
setiap parameter (Vehtari dkk., 2021) dan memberi peringatan bila R-hat
paling sedikit 1.05, ESS di bawah 100 per chain, atau terdapat divergent
transition; plot trace dan rank menunjukkan parameter mana yang
bermasalah:

``` r

fit_jax <- hb_area(y ~ x1 + x2 + x3,
  data = mys, vardir = "vardir",
  engine = jax_engine(chains = 4, samples = 1000)
)
fit_jax$rhat_max # below 1.01 is recommended
fit_jax$ess_min # above 100 per chain
head(fit_jax$mcmc_diagnostics) # R-hat and ESS of each parameter
autoplot(fit_jax, type = "rank") # chains near the dashed line: mixed
autoplot(fit_jax, type = "trace", pars = c("sigma2_u", "hb[5]"))
```

------------------------------------------------------------------------

## Struktur Spatial dan Spatio-Temporal

### Prior Spatial

1.  **BYM2 (`spatial = "bym2"`)**: Model Besag-York-Mollié yang di-scale
    (Riebler dkk., 2016) menguraikan random effect area menjadi komponen
    terstruktur spasial dan komponen IID tak terstruktur:
    ``` math
    u_d = \frac{1}{\sqrt{\tau_u}} \left( \sqrt{1 - \phi} \, v_d + \sqrt{\phi} \, u^*_d \right)
    ```
    dengan $`\phi \in [0, 1]`$ mengukur proporsi varians yang dijelaskan
    oleh struktur spasial dan memiliki penalised complexity (PC) prior
    $`P(\phi < 0.5) = 0.5`$ (Simpson dkk., 2017).
2.  **Besag ICAR (`spatial = "besag"`)**: Model spasial Intrinsic
    Conditional Autoregressive.
3.  **Efek area SAR (`spatial = "sar"`)**:
    $`u = (I - \rho W)^{-1}\varepsilon`$, random effect dari model
    Spatial Fay-Herriot pada
    [`eblup_sfh()`](https://ridsonap.is-a.dev/fastsae/reference/eblup_sfh.md).
4.  **Spatial Lag Model (`spatial = "slm"`)**:
    $`\eta = (I - \rho W)^{-1}(X\beta + \varepsilon)`$; covariate ikut
    difilter, sehingga ini bukan model
    [`eblup_sfh()`](https://ridsonap.is-a.dev/fastsae/reference/eblup_sfh.md).

### Model Temporal

Bila data diamati secara longitudinal pada $`T \ge 2`$ periode waktu,
tentukan `time`: - **`temporal = "rw1"`**: Random Walk orde pertama
($`\Delta t_\tau \sim N(0, \sigma_t^2)`$). - **`temporal = "rw2"`**:
Random Walk orde kedua untuk tren temporal nonlinear yang mulus. -
**`temporal = "ar1"`**: Proses autoregresif orde pertama dengan
autokorelasi $`|\rho_t| < 1`$. - **`temporal = "iid"`**: Random effect
waktu yang tak terstruktur.

### Struktur Interaksi Ruang-Waktu (`st_interaction`)

- **`"none"` (Aditif)**: Efek utama spasial dan tren temporal yang
  terpisah:
  ``` math
  \eta_{dt} = x_{dt}^\top \beta + s_d + \gamma_t
  ```
- **`"domain-specific"`** (sesuai `tipsae`): random walk / proses AR(1)
  spesifik domain dengan presisi bersama:
  ``` math
  \eta_{dt} = x_{dt}^\top \beta + s_d + t_{d, t}, \quad t_{d, \cdot} \sim \text{RW1}(\sigma_t^2)
  ```
- **`"separable"`** (sesuai `eblup_stfh` klasik): spatial random field
  yang berevolusi dinamis sepanjang waktu dengan korelasi temporal
  AR(1).
- **`"type1"` sampai `"type4"`**: interaksi Knorr-Held (2000)
  $`\delta_{dt}`$ yang ditambahkan ke efek utama,
  $`\eta_{dt} = x_{dt}^\top \beta + s_d + \gamma_t + \delta_{dt}`$:
  - `"type1"`: $`\delta_{dt}`$ saling bebas (ruang tak terstruktur
    $`\times`$ waktu tak terstruktur);
  - `"type2"`: proses `temporal` yang saling bebas untuk setiap domain
    (ruang tak terstruktur $`\times`$ waktu terstruktur);
  - `"type3"`: field ICAR yang saling bebas untuk setiap periode (ruang
    terstruktur $`\times`$ waktu tak terstruktur);
  - `"type4"`: field ICAR yang berevolusi mengikuti struktur `temporal`
    (ruang terstruktur $`\times`$ waktu terstruktur).

  `"type3"` dan `"type4"` memerlukan efek `spatial` terstruktur
  (`"bym2"`, `"bym"`, `"besag"`, atau `"leroux"`).

------------------------------------------------------------------------

## Studi Kasus: SAE Beta Spatio-Temporal pada Kemiskinan Emilia-Romagna

Kita mendemonstrasikan
[`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
dengan data kemiskinan Italia `emilia` dari paket `tipsae` (38 distrik
kesehatan di Emilia-Romagna selama 5 tahun, $`N = 190`$). Tujuannya
mengestimasi Head Count Ratio (`hcr`), yaitu proporsi rumah tangga di
bawah garis kemiskinan.

### 1. Memuat Data dan Menyusun Matriks Adjacency Spasial

``` r

library(fastsae)
library(tipsae)
library(spdep)

data("emilia")
data("emilia_shp")

# Construct binary spatial adjacency matrix W from district polygons
nb <- poly2nb(emilia_shp)
W <- nb2mat(nb, style = "B", zero.policy = TRUE)
rownames(W) <- colnames(W) <- as.character(emilia_shp$NAME_DISTRICT)

head(emilia[, c("id", "year", "hcr", "vars", "x")])
#>                    id year    hcr         vars       x
#> 1 CASALECCHIO DI RENO 2014 0.0404 9.090478e-05 -0.2624
#> 2   CITTA' DI BOLOGNA 2014 0.0825 6.404001e-05 -0.0008
#> 3               IMOLA 2014 0.1033 3.120275e-04 -0.0522
#> 4         PIANURA EST 2014 0.0633 1.025764e-04 -0.4007
#> 5       PIANURA OVEST 2014 0.0625 1.562500e-04 -0.2277
#> 6      PORRETTA TERME 2014 0.1276 6.643609e-04 -0.4434
```

### 2. Mengestimasi Model Beta Spatio-Temporal

Kita memodelkan proporsi kemiskinan dengan likelihood Beta, random
effect spasial Besag ICAR, dinamika temporal RW(1) spesifik domain, dan
sampling variance `vars` yang diketahui:

``` r

fit_beta_st <- hb_area(
  formula = hcr ~ x,
  data = emilia,
  domain = "id",
  time = "year",
  vardir = "vars",
  family = "beta",
  spatial = "besag",
  temporal = "rw1",
  st_interaction = "domain-specific",
  W = W,
  print_result = FALSE
)

summary(fit_beta_st)
#> 
#> Variance Components:
#> sigma2_u: 0.074071 
#> sigma2_t (temporal): 0.008529 
#> 
#> Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept) -2.2696e+00  1.5370e-02 -1.4766e+02  0.0000e+00 -2.2998e+00
#> x            4.1691e-01  6.6161e-02  6.3015e+00  2.9476e-10  2.8685e-01
#>             ci_upper
#> (Intercept)  -2.2395
#> x             0.5465
#> 
#> Hyperparameters:
#>                                  mean        sd 0.025quant 0.5quant 0.975quant
#> Precision for ..domain_id..  15.25279  5.508663   7.224529  14.3383   28.60678
#> Precision for ..time_id..   136.59658 57.258577  60.484581 125.0454  281.21531
#>                                  mode
#> Precision for ..domain_id..  12.67085
#> Precision for ..time_id..   104.92205
#> 
#> Goodness of Fit:
#>                                                   DIC 
#>                                            -971.94532 
#>                                                    pD 
#>                                              60.69149 
#>                                                  WAIC 
#>                                            -979.78129 
#>                                                 pWAIC 
#>                                              42.46186 
#> Marginal_LogLik.log marginal-likelihood (integration) 
#>                                             432.98607 
#> 
#> HB Summary Statistics:
#>        hb           linear_pred           sd                mse           
#>  Min.   :0.04905   Min.   :-2.970   Min.   :0.004907   Min.   :2.408e-05  
#>  1st Qu.:0.08142   1st Qu.:-2.428   1st Qu.:0.007907   1st Qu.:6.253e-05  
#>  Median :0.09444   Median :-2.266   Median :0.009723   Median :9.454e-05  
#>  Mean   :0.09673   Mean   :-2.268   Mean   :0.009911   Mean   :1.057e-04  
#>  3rd Qu.:0.10925   3rd Qu.:-2.103   3rd Qu.:0.011574   3rd Qu.:1.339e-04  
#>  Max.   :0.15774   Max.   :-1.678   Max.   :0.018014   Max.   :3.245e-04  
#>       rse        
#>  Min.   : 6.545  
#>  1st Qu.: 9.167  
#>  Median :10.224  
#>  Mean   :10.312  
#>  3rd Qu.:11.594  
#>  Max.   :14.895
```

### 3. Memeriksa Prediksi Posterior dan Ketidakpastiannya

Objek keluaran memuat prediksi lengkap per domain-periode, standard
error posterior, dan credible interval 95%:

``` r

head(fit_beta_st$df_hb[, c("domain", "time", "y", "hb", "sd", "rse", "ci_lower", "ci_upper")])
#>                domain time      y         hb          sd       rse   ci_lower
#> 1 CASALECCHIO DI RENO 2014 0.0404 0.05102824 0.005596856 10.968154 0.04062429
#> 2   CITTA' DI BOLOGNA 2014 0.0825 0.08182496 0.006281352  7.676572 0.07003562
#> 3               IMOLA 2014 0.1033 0.09310380 0.010890467 11.697124 0.07366716
#> 4         PIANURA EST 2014 0.0633 0.06365559 0.006408108 10.066843 0.05185445
#> 5       PIANURA OVEST 2014 0.0625 0.06937249 0.007927316 11.427176 0.05493781
#> 6      PORRETTA TERME 2014 0.1276 0.09446583 0.014070632 14.894944 0.06968690
#>     ci_upper
#> 1 0.06256395
#> 2 0.09467501
#> 3 0.11635698
#> 4 0.07700998
#> 5 0.08602892
#> 6 0.12479138
```

------------------------------------------------------------------------

## Perbandingan dan Validasi terhadap `tipsae` (Stan MCMC)

Pada
`tipsae::fit_sae(..., spatial_error = TRUE, temporal_error = TRUE)`,
model diestimasi dengan Hamiltonian Monte Carlo (HMC) di Stan.
[`fastsae::hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
mengestimasi model berstruktur sama dengan INLA; prior kedua paket
berbeda.

| Karakteristik | [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html) (Stan MCMC) | [`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md) (INLA) |
|:---|:---|:---|
| **Model** | Beta, Besag ICAR + RW(1) spesifik domain | `spatial = "besag"`, `temporal = "rw1"`, `st_interaction = "domain-specific"` |
| **Estimasi** | MCMC (NUTS) | Integrated nested Laplace approximation (mode `"compact"` bawaan INLA untuk Beta) |
| **Waktu eksekusi** | 22 menit (4 chain x 10.000 iterasi, $`\hat R`$ terbesar 1.007); 4 chain x 2.000 iterasi (120 s) belum konvergen | 1.1 s |
| **Posterior mean vs `tipsae` yang konvergen** | Acuan | $`r = 0.999`$, selisih absolut rata-rata 0.0017 |

Diukur pada 16-10-2026 (acuan `tipsae` diestimasi ulang pada 17-10-2026)
pada satu mesin (Apple silicon); waktu bergantung pada perangkat keras.
Akurasi
[`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
juga diperiksa terhadap Stan dengan model dan prior yang identik; lihat
bagian Engines pada
[`?hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md).

### Penjelajah Skalabilitas Interaktif ($`n = 30`$ sampai $`n = 1.000`$)

Gunakan kontrol interaktif di bawah untuk membandingkan runtime,
penggunaan memori, dan faktor speedup antara
[`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
dan [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html)
pada model **Standard Beta**, **Spatial Beta (Besag ICAR)**, dan
**Spatio-Temporal Beta**:

Standard Beta

Spatial Beta (ICAR)

Spatio-Temporal Beta

Execution Time

Memory Allocated

Throughput (iter/s)

Execution time (median, seconds) — Standard Beta · log scale

Speedup Factor & Efficiency Multiplier (fastsae vs tipsae) — Standard
Beta

#### Kesesuaian dengan Stan MCMC (model Beta, $`n = 1.000`$)

Untuk memastikan kecepatan
[`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
tidak mengorbankan kualitas estimasinya, kedua paket diestimasi dengan
default masing-masing pada data $`n = 1.000`$ dari benchmark Beta (efek
area independen, sampling variance dari direct estimate;
`benchmarks/bench_beta_concordance.R`). Kedua paket tidak mengestimasi
model yang identik (prior-nya berbeda, dan fastsae menetapkan presisi
pada estimasi tahap pertama sedangkan tipsae menjadikannya fungsi dari
rata-rata model), sehingga perbandingan ini adalah antara apa yang
diperoleh pengguna dari masing-masing paket:

| Metrik | [`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md) (INLA) | [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html) (Stan) |
|:---|:--:|:--:|
| RMSE terhadap rata-rata area sebenarnya | 0.0678 | 0.0684 |
| Coverage interval 95% | 93.3% | 91.1% |
| Lebar rata-rata interval 95% | 0.252 | 0.259 |
| Waktu | 4.3 s | 184 s |

Posterior mean kedua paket berbeda 0.0096 secara rata-rata (selisih
absolut), dan posterior SD fastsae adalah 0.99 kali SD tipsae pada area
median (90% area antara 0.88 dan 1.17). $`\hat{R}`$ terbesar pada fit
tipsae adalah 1.008.

------------------------------------------------------------------------

## Data Cacah Cross-Sectional: Poisson dan Negative Binomial

Untuk data cacah seperti kejadian penyakit atau peristiwa kriminal,
[`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
mendukung distribusi Poisson dan Negative Binomial dengan offset
populasi:

``` r

# Simulate 42-domain count data with known spatial structure
data("sim_area", package = "fastsae")

# The BYM2 model uses a neighbourhood graph: the 4 nearest centroids of each domain
W_knn <- create_weights(cbind(sim_area$x_coord, sim_area$y_coord), method = "knn", k = 4)

# Fit Poisson Spatial Model with population exposure offset
fit_pois <- hb_area(
  formula = y_poisson ~ x1 + x2,
  data = sim_area,
  exposure = "exposure",
  family = "poisson",
  spatial = "bym2",
  W = W_knn,
  print_result = FALSE
)

summary(fit_pois)
#> 
#> Variance Components:
#> sigma2_u: 0.069376 
#> rho (spatial): 0.1595 
#> phi (spatial fraction): 0.1595 
#> 
#> Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  1.4161e-01  1.3846e-01  1.0228e+00  3.0643e-01 -1.3129e-01
#> x1           3.8646e-01  5.0100e-02  7.7137e+00  1.2220e-14  2.8774e-01
#> x2          -2.9875e-01  3.3665e-02 -8.8741e+00  7.0481e-19 -3.6547e-01
#>             ci_upper
#> (Intercept)   0.4149
#> x1            0.4853
#> x2           -0.2327
#> 
#> Hyperparameters:
#>                                   mean        sd  0.025quant   0.5quant
#> Precision for ..domain_id.. 15.6289975 4.4994465 8.488621940 15.0603840
#> Phi for ..domain_id..        0.1595245 0.1562726 0.007220142  0.1055163
#>                             0.975quant        mode
#> Precision for ..domain_id..  26.045092 14.01273894
#> Phi for ..domain_id..         0.593418  0.01699927
#> 
#> Goodness of Fit:
#>                                                   DIC 
#>                                             317.74616 
#>                                                    pD 
#>                                              31.88700 
#>                                                  WAIC 
#>                                             309.67484 
#>                                                 pWAIC 
#>                                              17.15906 
#> Marginal_LogLik.log marginal-likelihood (integration) 
#>                                            -185.87987 
#> 
#> HB Summary Statistics:
#>        hb          linear_pred            sd               mse           
#>  Min.   :0.3061   Min.   :-1.1959   Min.   :0.02962   Min.   :0.0008776  
#>  1st Qu.:0.7648   1st Qu.:-0.2791   1st Qu.:0.05813   1st Qu.:0.0033800  
#>  Median :1.3850   Median : 0.3240   Median :0.07685   Median :0.0059062  
#>  Mean   :1.5294   Mean   : 0.2381   Mean   :0.16304   Mean   :0.0690516  
#>  3rd Qu.:2.3016   3rd Qu.: 0.8023   3rd Qu.:0.12189   3rd Qu.:0.0148681  
#>  Max.   :3.6016   Max.   : 1.2810   Max.   :0.76549   Max.   :0.5859684  
#>       rse        
#>  Min.   : 2.791  
#>  1st Qu.: 5.164  
#>  Median : 6.265  
#>  Mean   :10.478  
#>  3rd Qu.: 9.915  
#>  Max.   :30.137
```

------------------------------------------------------------------------

## Daftar Pustaka

- De Nicolò, S., & Gardini, A. (2024). The R package tipsae: Tools for
  mapping proportions and indicators on the unit interval. *Journal of
  Statistical Software*, 108(1), 1–36.
- Janicki, R. (2020). Properties of the beta regression model for small
  area estimation of proportions and application to estimation of
  poverty rates. *Communications in Statistics - Theory and Methods*,
  49(9), 2264–2284.
- Knorr-Held, L. (2000). Bayesian modelling of inhomogeneous spatial and
  temporal variation in rates. *Statistics in Medicine*, 19(17-18),
  2555–2567.
- Riebler, A., Sørbye, S. H., Simpson, D., & Rue, H. (2016). An
  intuitive Bayesian spatial model for disease mapping that accounts for
  scaling. *Statistical Methods in Medical Research*, 25(4), 1145–1165.
- Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian
  inference for latent Gaussian models by using integrated nested
  Laplace approximations. *Journal of the Royal Statistical Society:
  Series B*, 71(2), 319–392.
- Simpson, D., Rue, H., Riebler, A., Martins, T. G., & Sørbye, S. H.
  (2017). Penalising model component complexity: A principled, practical
  approach to constructing priors. *Statistical Science*, 32(1), 1–28.
- Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C.
  (2021). Rank-normalization, folding, and localization: An improved
  R-hat for assessing convergence of MCMC. *Bayesian Analysis*, 16(2),
  667–718.
