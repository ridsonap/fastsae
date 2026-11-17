# Estimasi Area Kecil Bayesian Tingkat Area dengan INLA

## Pendahuluan

Model Estimasi Area Kecil (SAE) klasik, seperti model Fay-Herriot
(1979), umumnya mengasumsikan bahwa estimator survei langsung
berdistribusi Gaussian dengan varians galat sampling yang diketahui.
Padahal, dalam statistik resmi sehari-hari, banyak indikator yang tidak
berdistribusi Gaussian:

- **Tingkat kemiskinan, tingkat pengangguran, dan rasio prevalensi**:
  nilainya selalu berada dalam selang terbuka $`(0, 1)`$.
- **Jumlah kasus penyakit, kejadian kriminal, dan peristiwa langka**:
  cacahan diskret yang mengikuti Poisson atau Binomial Negatif (bila
  terjadi overdispersi), dengan paparan berupa jumlah penduduk.
- **Agregat biner hasil survei**: proporsi yang dihitung dari ukuran
  sampel kecil (percobaan Binomial).
- **Pengeluaran ekonomi yang menceng**: respons kontinu yang selalu
  positif (Gamma).

Selain itu, bila data dikumpulkan di berbagai wilayah geografis dan
diulang dari waktu ke waktu (survei panel atau longitudinal), responsnya
menunjukkan ketergantungan spasial sekaligus autokorelasi temporal.

Paket **fastsae** menyediakan
[`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md),
kerangka Bayesian hierarkis terpadu yang ditenagai **INLA** (*Integrated
Nested Laplace Approximations*). INLA memberikan aproksimasi
deterministik yang cepat untuk distribusi marginal posterior, tanpa
rantai, *warm-up*, dan penyetelan seperti pada *Markov chain Monte
Carlo* (MCMC). Dibandingkan dengan Stan yang dijalankan panjang dengan
model dan prior yang sama, rata-rata dan simpangan baku posteriornya
sangat berdekatan (lihat bagian *Engines* di
[`?hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)).

------------------------------------------------------------------------

## Keluarga Distribusi yang Didukung

[`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
mendukung 6 distribusi peluang untuk respons langsung:

| Keluarga | Kegunaan umum dalam SAE | Masukan sampling yang diperlukan | Parameter skala model |
|:---|:---|:---|:---|
| **`"gaussian"`** | Indikator kontinu (misalnya rata-rata pendapatan, rata-rata lama sekolah) | `vardir` (varians sampling $`D_d`$) | Kebalikan varians sampling, $`\text{scale} = 1/D_d`$ |
| **`"beta"`** | Tingkat dan proporsi yang nilainya selalu dalam $`(0, 1)`$ | `vardir` atau `trials` | Presisi $`\phi_d = \hat\mu_d(1 - \hat\mu_d)/V_d - 1`$ dari estimasi tahap pertama (Janicki, 2020), atau $`n_d - 1`$ dari `trials` |
| **`"poisson"`** | Cacahan dan tingkat kejadian langka | `exposure` (cacahan harapan / *offset* jumlah penduduk) | *Offset* pada link log, $`\log(E_d)`$ |
| **`"nbinomial"`** | Cacahan yang overdispersi dan kejadian yang mengelompok | `exposure`; `vardir` bersifat opsional | *Offset* pada link log; ukuran $`\hat\mu_d^2/(V_d - \hat\mu_d)`$ dari estimasi tahap pertama bila `vardir` diberikan |
| **`"binomial"`** | Proporsi biner teragregasi | `trials` (ukuran sampel domain $`n_d`$, atau ukuran sampel efektif $`n_d/\mathrm{deff}_d`$ untuk proporsi berbobot survei; tidak dibulatkan) | Banyaknya percobaan per domain |
| **`"gamma"`** | Indikator positif yang menceng ke kanan (misalnya pengeluaran) | `vardir` | Presisi (*shape*) $`y_d^2/V_d`$ |

------------------------------------------------------------------------

## Formulasi Matematis

Model Bayesian hierarkis yang mendasari
[`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
memperluas kerangka Fay-Herriot klasik agar dapat menangani respons
non-Gaussian dan efek acak yang terstruktur. Model ini terdiri atas tiga
tingkat.

### Tingkat 1: Model Sampling (Likelihood)

Untuk domain $`d`$ dengan estimator langsung $`y_d`$ dan varians
sampling yang diketahui $`D_d`$:

``` math
\pi(y_d \mid \eta_d) = \text{ExpFamily}(y_d \mid \eta_d, D_d)
```

Fungsi link menghubungkan prediktor linear $`\eta_d`$ dengan rata-rata
$`\mu_d = E(y_d)`$:

| Keluarga | Fungsi link | Prediktor linear |
|:---|:---|:---|
| Gaussian | Identitas: $`\eta_d = \mu_d`$ | $`\eta_d = x_d^\top \beta + u_d`$ |
| Beta | Logit: $`\eta_d = \log(\mu_d / (1-\mu_d))`$ | $`\eta_d = x_d^\top \beta + u_d`$ |
| Poisson | Log: $`\eta_d = \log(\mu_d)`$ | $`\eta_d = x_d^\top \beta + u_d + \log(E_d)`$ |
| Binomial | Logit: $`\eta_d = \log(\mu_d / (1-\mu_d))`$ | $`\eta_d = x_d^\top \beta + u_d`$ |
| Binomial Negatif | Log: $`\eta_d = \log(\mu_d)`$ | $`\eta_d = x_d^\top \beta + u_d + \log(E_d)`$ |
| Gamma | Log: $`\eta_d = \log(\mu_d)`$ | $`\eta_d = x_d^\top \beta + u_d`$ |

### Tingkat 2: Prior untuk Efek Tetap

``` math
\beta_0 \propto 1, \qquad \beta_j \sim N(0, 1000), \quad j = 1, \ldots, p
```

*Intercept* memakai prior datar bawaan INLA, sedangkan koefisien lainnya
memakai prior normal dengan varians 1000, yang hampir datar pada skala
sebagian besar kovariat.

### Tingkat 3: Prior untuk Efek Acak

``` math
u \mid \sigma_u \sim N(0, \sigma_u^2 R^{-}), \qquad \sigma_u \sim \text{half-normal}(s)
```

Efek acak $`u = (u_1, \ldots, u_m)^\top`$ mengikuti medan acak Markov
Gaussian (*Gaussian Markov random field*) dengan rata-rata nol, dan
matriks strukturnya $`R`$ menyatakan ketergantungan spasial atau
temporal yang dipilih ($`R = I`$ untuk efek iid; selain itu Besag, BYM2,
RW1, dan seterusnya). Skala prior *half-normal* untuk simpangan bakunya
adalah $`s = \mathrm{sd}(y)`$, yaitu simpangan baku estimasi langsung,
untuk model Gaussian, dan $`s = 1`$ pada skala logit atau log untuk
keluarga lainnya. Hiperparameter lainnya (parameter pencampur BYM2,
autokorelasi, parameter dispersi) beserta alasan pemilihan nilai
bawaannya dijelaskan di [artikel prior
bawaan](https://ridsonap.is-a.dev/fastsae/articles/id/hb-priors.md);
semuanya dapat diganti melalui
[`inla_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md).

### Aproksimasi INLA

Alih-alih penarikan sampel MCMC, INLA mengaproksimasi distribusi
marginal posterior dari hiperparameter $`\theta`$ dan dari medan laten
$`x = (\beta, u)`$ (Rue dkk., 2009):

``` math
\tilde{\pi}(\theta \mid y) \propto \left. \frac{\pi(x, \theta, y)}{\tilde{\pi}_G(x \mid \theta, y)} \right|_{x = x^*(\theta)}
```

dengan $`\tilde{\pi}_G`$ sebagai aproksimasi Gaussian bagi distribusi
bersyarat penuh dari $`x`$ dan $`x^*(\theta)`$ sebagai modusnya, serta

``` math
\tilde{\pi}(x_d \mid y) = \sum_k \tilde{\pi}(x_d \mid \theta_k, y) \, \tilde{\pi}(\theta_k \mid y) \, \Delta_k
```

dengan integral yang dihitung melalui **aproksimasi Laplace** dan
**integrasi numerik**. Ada tiga strategi komputasi yang dapat dipilih
melalui `hb_area(..., engine = inla_engine(strategy = ...))`; fungsi ini
juga menerima prior dan argumen
[`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html) lainnya:

- **`strategy = "gaussian"`**: paling cepat, paling cocok untuk masalah
  yang berperilaku baik;
- **`strategy = "laplace"`**: bawaan, keseimbangan yang baik antara
  akurasi dan kecepatan;
- **`strategy = "simplified.laplace"`**: Laplace yang sedikit lebih
  cepat, cocok untuk model besar.

Model yang sama dapat diestimasi dengan MCMC di CPU atau GPU melalui
`engine = jax_engine()` (NUTS NumPyro di atas JAX; lihat
[`?jax_engine`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)).
Hasil MCMC baru dapat dipakai bila rantainya sudah konvergen. Hasil
estimasi melaporkan R-hat ternormalisasi peringkat (*rank-normalized*)
dan ukuran sampel efektif (ESS) setiap parameter (Vehtari dkk., 2021)
serta E-BFMI setiap rantai (Betancourt, 2016). Peringatan akan muncul
bila R-hat paling sedikit 1,05, ESS kurang dari 100 per rantai, terjadi
transisi divergen, atau E-BFMI suatu rantai di bawah 0,3. Grafik *trace*
dan grafik peringkat memperlihatkan parameter mana yang bermasalah:

``` r

fit_jax <- hb_area(y ~ x1 + x2 + x3,
  data = mys, vardir = "vardir",
  engine = jax_engine(chains = 4, samples = 1000)
)
fit_jax$rhat_max # disarankan di bawah 1,01
fit_jax$ess_min # di atas 100 per rantai
fit_jax$ebfmi # satu nilai per rantai; di bawah 0,3 memicu peringatan
head(fit_jax$mcmc_diagnostics) # R-hat dan ESS setiap parameter
autoplot(fit_jax, type = "rank") # rantai di sekitar garis putus-putus: sudah tercampur baik
autoplot(fit_jax, type = "trace", pars = c("sigma2_u", "hb[5]"))
```

Banyaknya transisi divergen tercatat di `fit_jax$divergences`, dan
`fit_jax$divergent_draws` menandai tarikan mana yang divergen (iterasi x
rantai). `type = "divergence"` memplot tarikan dua parameter dengan
tarikan divergen berwarna merah (Gabry dkk., 2019). Secara bawaan,
grafik ini menampilkan area dengan R-hat terbesar terhadap logaritma
varians area, yaitu bentuk corong (*funnel*) yang khas pada model
hierarkis; parameter lain dapat dipilih dengan `pars` dan `log`.
Divergensi yang berkumpul di satu sudut, biasanya di tempat variansnya
kecil, menunjukkan bagian posterior yang tidak dapat diikuti oleh
*sampler*. `target_accept` yang lebih besar membuat langkahnya lebih
kecil, yang biasanya menghilangkan divergensi tersebut meskipun
penarikan sampel menjadi lebih lama:

``` r

fit_jax$divergences
autoplot(fit_jax, type = "divergence")
autoplot(fit_jax, type = "divergence", pars = c("x1", "sigma2_u"))

# bawaannya 0,85; nilai yang dipakai tersimpan di fit_jax2$fit$target_accept_prob
fit_jax2 <- hb_area(y ~ x1 + x2 + x3,
  data = mys, vardir = "vardir",
  engine = jax_engine(chains = 4, target_accept = 0.95)
)
fit_jax2$divergences
```

------------------------------------------------------------------------

## Struktur Spasial dan Spatio-Temporal

### Prior Spasial

1.  **BYM2 (`spatial = "bym2"`)**: model Besag-York-Mollié terskala
    (Riebler dkk., 2016) menguraikan efek acak area menjadi komponen
    spasial yang terstruktur dan komponen IID yang tidak terstruktur:
    ``` math
    u_d = \frac{1}{\sqrt{\tau_u}} \left( \sqrt{1 - \phi} \, v_d + \sqrt{\phi} \, u^*_d \right)
    ```
    dengan $`\phi \in [0, 1]`$ sebagai proporsi varians yang dijelaskan
    oleh struktur spasial, dengan prior *penalised complexity* (PC)
    $`P(\phi < 0{,}5) = 0{,}5`$ (Simpson dkk., 2017).
2.  **Besag ICAR (`spatial = "besag"`)**: model spasial *Intrinsic
    Conditional Autoregressive*.
3.  **BYM (`spatial = "bym"`)**: model Besag-York-Mollié klasik, yaitu
    efek iid dan efek Besag yang masing-masing memiliki variansnya
    sendiri.
4.  **Leroux (`spatial = "leroux"`, di INLA disebut `"generic1"`)**:
    presisi $`\sigma_u^{-2}\{(1 - \lambda) I + \lambda R\}`$, campuran
    struktur iid dan ICAR.
5.  **Efek area SAR (`spatial = "sar"`)**:
    $`u = (I - \rho W)^{-1}\varepsilon`$, yaitu efek acak model
    Fay-Herriot spasial pada
    [`eblup_sfh()`](https://ridsonap.is-a.dev/fastsae/reference/eblup_sfh.md).
6.  **Model *Spatial Lag* (`spatial = "slm"`)**:
    $`\eta = (I - \rho W)^{-1}(X\beta + \varepsilon)`$; kovariatnya ikut
    tersaring, sehingga model ini bukan model
    [`eblup_sfh()`](https://ridsonap.is-a.dev/fastsae/reference/eblup_sfh.md).

### Model Temporal

Bila data diamati secara longitudinal selama $`T \ge 2`$ periode,
tentukan argumen `time`:

- **`temporal = "rw1"`**: *random walk* orde pertama
  ($`\Delta t_\tau \sim N(0, \sigma_t^2)`$);
- **`temporal = "rw2"`**: *random walk* orde kedua untuk tren temporal
  nonlinear yang mulus;
- **`temporal = "ar1"`**: proses autoregresif orde pertama dengan
  autokorelasi $`|\rho_t| < 1`$;
- **`temporal = "iid"`**: efek acak waktu yang tidak terstruktur.

### Struktur Interaksi Ruang-Waktu (`st_interaction`)

- **`"none"` (aditif)**: efek utama spasial dan tren temporal yang
  terpisah:
  ``` math
  \eta_{dt} = x_{dt}^\top \beta + s_d + \gamma_t
  ```
- **`"domain-specific"`** (setara dengan `tipsae`): *random walk* /
  proses AR(1) khusus per domain dengan presisi yang sama:
  ``` math
  \eta_{dt} = x_{dt}^\top \beta + s_d + t_{d, t}, \quad t_{d, \cdot} \sim \text{RW1}(\sigma_t^2)
  ```
- **`"separable"`** (setara dengan `eblup_stfh` klasik): medan acak
  spasial yang berkembang dari waktu ke waktu dengan korelasi temporal
  AR(1).
- **`"type1"` hingga `"type4"`**: interaksi Knorr-Held (2000)
  $`\delta_{dt}`$ yang ditambahkan pada efek utama,
  $`\eta_{dt} = x_{dt}^\top \beta + s_d + \gamma_t + \delta_{dt}`$:
  - `"type1"`: $`\delta_{dt}`$ saling bebas (ruang tak terstruktur
    $`\times`$ waktu tak terstruktur);
  - `"type2"`: proses `temporal` yang saling bebas untuk setiap domain
    (ruang tak terstruktur $`\times`$ waktu terstruktur);
  - `"type3"`: medan ICAR yang saling bebas untuk setiap periode (ruang
    terstruktur $`\times`$ waktu tak terstruktur);
  - `"type4"`: medan ICAR yang berkembang mengikuti struktur `temporal`
    (ruang terstruktur $`\times`$ waktu terstruktur).

  `"type3"` dan `"type4"` memerlukan efek `spatial` yang terstruktur
  (`"bym2"`, `"bym"`, `"besag"`, atau `"leroux"`).

------------------------------------------------------------------------

## Studi Kasus: SAE Beta Spatio-Temporal untuk Kemiskinan di Emilia-Romagna

[`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
diperagakan dengan data kemiskinan Italia `emilia` dari paket `tipsae`
(38 distrik kesehatan di Emilia-Romagna selama 5 tahun, $`N = 190`$).
Tujuannya adalah mengestimasi *Head Count Ratio* (`hcr`), yaitu proporsi
rumah tangga yang berada di bawah garis kemiskinan.

### 1. Memuat Data dan Menyusun Matriks Ketetanggaan Spasial

``` r

library(fastsae)
library(tipsae)
library(spdep)

data("emilia")
data("emilia_shp")

# Susun matriks ketetanggaan spasial biner W dari poligon distrik
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

Proporsi kemiskinan dimodelkan dengan likelihood Beta, efek acak spasial
Besag ICAR, dinamika temporal RW(1) khusus per domain, dan varians
sampling `vars` yang diketahui:

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
#> sigma2_u: 0.074072 
#> sigma2_t (temporal): 0.008529 
#> 
#> Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept) -2.2696e+00  1.5370e-02 -1.4766e+02  0.0000e+00 -2.2998e+00
#> x            4.1691e-01  6.6161e-02  6.3015e+00  2.9475e-10  2.8685e-01
#>             ci_upper
#> (Intercept)  -2.2395
#> x             0.5465
#> 
#> Hyperparameters:
#>                                  mean        sd 0.025quant  0.5quant 0.975quant
#> Precision for ..domain_id..  15.25259  5.508617   7.224329  14.33813    28.6064
#> Precision for ..time_id..   136.59180 57.254571  60.482863 125.04200   281.1988
#>                                  mode
#> Precision for ..domain_id..  12.67073
#> Precision for ..time_id..   104.92080
#> 
#> Goodness of Fit:
#>                                                   DIC 
#>                                            -971.94509 
#>                                                    pD 
#>                                              60.69169 
#>                                                  WAIC 
#>                                            -979.78129 
#>                                                 pWAIC 
#>                                              42.46184 
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

Objek hasil berisi prediksi untuk setiap pasangan domain-periode,
simpangan baku posterior, dan selang kredibel 95%:

``` r

head(fit_beta_st$df_hb[, c("domain", "time", "y", "hb", "sd", "rse", "ci_lower", "ci_upper")])
#>                domain time      y         hb          sd       rse   ci_lower
#> 1 CASALECCHIO DI RENO 2014 0.0404 0.05102822 0.005596847 10.968141 0.04062430
#> 2   CITTA' DI BOLOGNA 2014 0.0825 0.08182494 0.006281352  7.676575 0.07003561
#> 3               IMOLA 2014 0.1033 0.09310383 0.010890468 11.697121 0.07366716
#> 4         PIANURA EST 2014 0.0633 0.06365559 0.006408114 10.066852 0.05185445
#> 5       PIANURA OVEST 2014 0.0625 0.06937249 0.007927324 11.427186 0.05493781
#> 6      PORRETTA TERME 2014 0.1276 0.09446585 0.014070678 14.894990 0.06968683
#>     ci_upper
#> 1 0.06256394
#> 2 0.09467500
#> 3 0.11635699
#> 4 0.07700999
#> 5 0.08602894
#> 6 0.12479143
```

------------------------------------------------------------------------

## Perbandingan dan Validasi terhadap `tipsae` (MCMC Stan)

Pada
`tipsae::fit_sae(..., spatial_error = TRUE, temporal_error = TRUE)`,
model diestimasi dengan *Hamiltonian Monte Carlo* (HMC) di Stan.
[`fastsae::hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
mengestimasi model dengan struktur yang sama melalui INLA; prior kedua
paket berbeda.

| Karakteristik | [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html) (MCMC Stan) | [`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md) (INLA) |
|:---|:---|:---|
| **Model** | Beta, Besag ICAR + RW(1) khusus per domain | `spatial = "besag"`, `temporal = "rw1"`, `st_interaction = "domain-specific"` |
| **Estimasi** | MCMC (NUTS) | *Integrated nested Laplace approximation* (mode `"compact"` bawaan INLA untuk Beta) |
| **Waktu eksekusi** | 22 menit (4 rantai x 10.000 iterasi, $`\hat R`$ terbesar 1,007); dengan 4 rantai x 2.000 iterasi (120 detik), rantainya belum konvergen | 1,1 detik |
| **Rata-rata posterior vs `tipsae` yang sudah konvergen** | Rujukan | $`r = 0{,}999`$, rata-rata selisih absolut 0,0017 |

Diukur pada 2026-10-16 (rujukan `tipsae` diestimasi ulang pada
2026-10-17) di satu mesin (Apple silicon); waktu eksekusi bergantung
pada perangkat kerasnya. Akurasi
[`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
juga sudah diperiksa terhadap Stan dengan model dan prior yang identik;
lihat bagian *Engines* di
[`?hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md).

### Penjelajah Skalabilitas Interaktif ($`n = 30`$ hingga $`n = 1.000`$)

Gunakan tombol-tombol di bawah untuk membandingkan waktu eksekusi,
kebutuhan memori, dan faktor percepatan antara
[`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
dan [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html),
baik untuk model **Beta Standar** maupun **Beta Spasial (Besag ICAR)**:

Beta Standar

Beta Spasial (ICAR)

Beta Spatio-Temporal

Waktu Eksekusi

Memori yang Dialokasikan

Laju (iterasi/detik)

Waktu eksekusi (median, detik) — Beta Standar · skala log

Faktor Percepatan & Pengali Efisiensi (fastsae vs tipsae) — Beta Standar

#### Kesesuaian dengan MCMC Stan (model Beta, $`n = 1.000`$)

Untuk memastikan bahwa kecepatan
[`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
tidak mengorbankan kualitas estimasinya, kedua paket dijalankan dengan
pengaturan bawaannya pada himpunan data $`n = 1.000`$ dari tolok ukur
Beta (efek area saling bebas, varians sampling dari estimasi langsung;
`benchmarks/bench_beta_concordance.R`). Model yang diestimasi kedua
paket tidak identik (prior keduanya berbeda; fastsae menetapkan presisi
pada estimasi tahap pertama, sedangkan tipsae menjadikannya fungsi dari
rata-rata model), sehingga yang dibandingkan adalah hasil yang diperoleh
pengguna dari masing-masing paket:

| Metrik | [`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md) (INLA) | [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html) (Stan) |
|:---|:--:|:--:|
| RMSE terhadap rata-rata area yang sebenarnya | 0,0678 | 0,0684 |
| Cakupan selang 95% | 93,3% | 91,1% |
| Rata-rata lebar selang 95% | 0,252 | 0,259 |
| Waktu | 4,3 detik | 184 detik |

Rata-rata posterior kedua paket rata-rata berselisih 0,0096 (selisih
absolut), dan simpangan baku posterior fastsae sebesar 0,99 kali
simpangan baku tipsae pada area median (90% area berada antara 0,88 dan
1,17). $`\hat{R}`$ terbesar pada hasil tipsae adalah 1,008.

------------------------------------------------------------------------

## Data Cacahan Lintas Area: Poisson & Binomial Negatif

Untuk data cacahan seperti insidensi penyakit atau kejadian kriminal,
[`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
mendukung distribusi Poisson dan Binomial Negatif dengan *offset* jumlah
penduduk:

``` r

# Data cacahan simulasi untuk 42 domain dengan struktur spasial yang diketahui
data("sim_area", package = "fastsae")

# Model BYM2 memakai graf ketetanggaan: 4 sentroid terdekat dari setiap domain
W_knn <- create_weights(cbind(sim_area$x_coord, sim_area$y_coord), method = "knn", k = 4)

# Estimasi model spasial Poisson dengan offset paparan jumlah penduduk
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
#> x2          -2.9875e-01  3.3665e-02 -8.8741e+00  7.0478e-19 -3.6547e-01
#>             ci_upper
#> (Intercept)   0.4149
#> x1            0.4853
#> x2           -0.2327
#> 
#> Hyperparameters:
#>                                   mean        sd  0.025quant   0.5quant
#> Precision for ..domain_id.. 15.6289890 4.4994215 8.488628544 15.0603865
#> Phi for ..domain_id..        0.1595236 0.1562695 0.007220675  0.1055172
#>                             0.975quant        mode
#> Precision for ..domain_id.. 26.0450024 14.01276535
#> Phi for ..domain_id..        0.5934074  0.01700112
#> 
#> Goodness of Fit:
#>                                                   DIC 
#>                                             317.74616 
#>                                                    pD 
#>                                              31.88700 
#>                                                  WAIC 
#>                                             309.67485 
#>                                                 pWAIC 
#>                                              17.15907 
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
#>  Max.   :3.6016   Max.   : 1.2810   Max.   :0.76549   Max.   :0.5859676  
#>       rse        
#>  Min.   : 2.791  
#>  1st Qu.: 5.164  
#>  Median : 6.265  
#>  Mean   :10.478  
#>  3rd Qu.: 9.915  
#>  Max.   :30.137
```

------------------------------------------------------------------------

## Kovariat dengan Galat Pengukuran

Kovariat tingkat area sering kali juga merupakan hasil estimasi,
misalnya rata-rata dari survei lain yang memiliki varians samplingnya
sendiri. Bila dipakai seolah-olah nilainya pasti, koefisien kovariat
tersebut mengalami atenuasi, semua area disusutkan dengan kadar yang
sama betapa pun tidak presisinya kovariat masing-masing, dan estimasi
berbasis model bisa lebih buruk daripada estimasi langsungnya (Ybarra
dan Lohr, 2008). `var_x` memuat varians galat-galat tersebut: nilai
sebenarnya dari kovariat lalu menggantikan nilai teramati di dalam
model, dengan distribusi normal (`me_model = "structural"`, bawaan; Bell
dkk., 2019) atau prior datar (`me_model = "functional"`, padanan
Bayesian dari model Ybarra–Lohr; Arima, Datta, dan Liseo, 2015). Varians
0 menandai nilai yang tidak mengandung galat.

``` r

set.seed(2)
D <- 60
x_true <- rnorm(D, 5, 2)                              # kovariat yang sebenarnya, tidak teramati
var_x1 <- sample(c(0.05, 0.5, 4), D, replace = TRUE)  # varians estimasinya
d_me <- data.frame(
  area = seq_len(D),
  x1 = x_true + rnorm(D, 0, sqrt(var_x1)),            # hasil estimasinya
  var_x1 = var_x1,
  vardir = runif(D, 0.5, 1.5)
)
d_me$y <- 1 + 0.8 * x_true + rnorm(D, 0, 0.7) + rnorm(D, 0, sqrt(d_me$vardir))
d_me$y[1:6] <- NA                                     # area tanpa estimasi langsung

fit_naive <- hb_area(y ~ x1, data = d_me, domain = "area", vardir = "vardir", print_result = FALSE)
fit_me <- hb_area(y ~ x1,
  data = d_me, domain = "area", vardir = "vardir",
  var_x = c(x1 = "var_x1"), print_result = FALSE
)

# Koefisien kovariat: 0,8 dalam simulasi
rbind(
  naive = fit_naive$estcoef["x1", c("beta", "ci_lower", "ci_upper")],
  structural = fit_me$estcoef["x1", c("beta", "ci_lower", "ci_upper")]
)
#>                 beta  ci_lower  ci_upper
#> naive      0.5855003 0.4338812 0.7374280
#> structural 0.7276754 0.5812973 0.8782197

# Rata-rata simpangan baku posterior area menurut varians galat kovariatnya
sds <- data.frame(var_x1 = d_me$var_x1, naive = fit_naive$df_hb$sd, structural = fit_me$df_hb$sd)
aggregate(cbind(naive, structural) ~ var_x1, data = sds, FUN = mean)
#>   var_x1     naive structural
#> 1   0.05 0.8049362  0.6122791
#> 2   0.50 0.7694051  0.6625098
#> 3   4.00 0.7674870  0.7990784
```

Selang koefisien dari model naif tidak mencakup nilai sebenarnya, dan
model naif memberikan simpangan baku posterior yang hampir sama untuk
area dengan kovariat yang presisi maupun yang tidak. Model galat
pengukuran mempersempit selang untuk area yang kovariatnya presisi dan
memperlebar selang untuk area yang kovariatnya tidak presisi. Dengan
`me_model = "functional"`, estimasi area sesuai dengan EBLUP Ybarra–Lohr
(`emdi::fh(method = "me")`), tetapi koefisiennya tetap dekat dengan
koefisien model naif: mengintegralkan nilai sebenarnya menyisakan
regresi pada nilai teramati, hanya dengan varians yang lebih besar.
Rinciannya, termasuk validasi terhadap Stan, ada di bagian *Covariates
measured with error* pada
[`?hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md).

------------------------------------------------------------------------

## Subindikator dari Indikator Induk: `hb_ratio()`

Ada indikator yang merupakan bagian dari indikator lain: kemiskinan
ekstrem di dalam kemiskinan, atau *stunting* berat di dalam *stunting*.
Bila indikator induk $`P_d`$ sudah diestimasi (misalnya tingkat
kemiskinan yang sudah dirilis secara resmi), memodelkan subindikator
secara terpisah dapat menghasilkan estimasi yang melebihi induknya.
[`hb_ratio()`](https://ridsonap.is-a.dev/fastsae/reference/hb_ratio.md)
memodelkan pangsa $`s_d`$, yaitu bagian populasi induk yang memiliki
karakteristik turunan tersebut, dengan model
[`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
binomial yang banyaknya percobaan adalah kasus induk yang tersampel,
lalu mengestimasi $`\theta_d = s_d P_d`$. Karena $`0 \le s_d \le 1`$,
tidak ada estimasi yang melebihi induknya.

Model binomial memerlukan banyaknya percobaan *efektif*,
$`n_d / \mathrm{deff}_d`$, agar varians binomialnya sama dengan varians
desain.
[`effective_trials()`](https://ridsonap.is-a.dev/fastsae/reference/effective_trials.md)
menghitungnya dari proporsi langsung, varians desainnya, dan ukuran
sampel. Efek desain untuk satu area saja tidak stabil, dan tidak
terdefinisi bila proporsi langsungnya 0 atau 1, sehingga secara bawaan
efek desain digabungkan di dalam kelompok (di sini provinsi):

``` r

set.seed(1)
D <- 60
d_r <- data.frame(area = 1:D, prov = rep(1:6, each = 10), x = rnorm(D), pop = rpois(D, 5e4))
d_r$p0 <- plogis(-2 + 0.3 * rnorm(D)) # tingkat indikator induk (kemiskinan) yang sudah dirilis
d_r$p0_se <- 0.15 * d_r$p0 # beserta galat bakunya
d_r$n_poor <- rpois(D, 40) # rumah tangga miskin yang tersampel
s_true <- plogis(-1.5 + 0.5 * d_r$x + rnorm(D, 0, 0.3))
d_r$s_dir <- rbinom(D, d_r$n_poor, s_true) / d_r$n_poor # pangsa langsung penduduk miskin ekstrem
d_r$v_s <- 1.5 * d_r$s_dir * (1 - d_r$s_dir) / d_r$n_poor # varians desainnya (deff 1,5)

eff <- effective_trials(~s_dir, vardir = ~v_s, n = ~n_poor, group = ~prov, data = d_r)
head(eff)
#>    n deff deff_source   trials
#> 1 45  1.5       group 30.00000
#> 2 32  1.5       group 21.33333
#> 3 33  1.5       group 22.00000
#> 4 33  1.5       group 22.00000
#> 5 36  1.5       group 24.00000
#> 6 42  1.5       group 28.00000
d_r$n_eff <- eff$trials
```

`parent_se` menambahkan ketidakpastian indikator induk ke `sd`, `mse`,
dan `rse` (`rse_cond` adalah RSE bersyarat pada indikator induk).
Benchmarking pada skala logit dari $`s_d`$ membuat rata-rata berbobot
penduduk sama dengan target (di sini tingkat nasional 2,5%) dan tetap
menjaga estimasi di bawah induknya:

``` r

fit_ratio <- hb_ratio(s_dir ~ x,
  data = d_r, trials = "n_eff", parent = "p0", parent_se = "p0_se", domain = "area",
  benchmark_target = 0.025, benchmark_weight = "pop", print_result = FALSE
)
head(fit_ratio$df_hb[, c("domain", "hb", "rse", "rse_cond", "ratio", "parent")])
#>   domain         hb      rse  rse_cond     ratio     parent
#> 1      1 0.01082259 19.87851 12.837866 0.1369326 0.07903595
#> 2      2 0.02623409 18.05624 10.474397 0.2030470 0.12920204
#> 3      3 0.01018307 20.10055 13.716438 0.1294855 0.07864257
#> 4      4 0.04154269 18.28158 10.153290 0.3775102 0.11004389
#> 5      5 0.02305273 17.98405  9.478757 0.2226433 0.10354109
#> 6      6 0.01270608 20.53951 13.435543 0.1268784 0.10014374
all(fit_ratio$df_hb$hb <= d_r$p0) # tidak pernah melebihi induknya
#> [1] TRUE
fit_ratio$benchmark$achieved # rata-rata berbobot setelah benchmarking
#>   all 
#> 0.025
```

`fit_ratio$ratio_fit` adalah hasil
[`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
untuk pangsa tersebut, yang dapat diperiksa dengan
[`diagnose()`](https://ridsonap.is-a.dev/fastsae/reference/diagnose.md).
`fit_ratio$draws` berisi tarikan posterior $`\theta_d`$ (area x
tarikan). Dengan tarikan ini,
[`add_reliability_flags()`](https://ridsonap.is-a.dev/fastsae/reference/add_reliability_flags.md)
memberikan kelas setiap estimasi beserta peluang posterior bahwa nilai
sebenarnya berada di kelas itu. Untuk indikator yang langka, estimasi
yang mendekati 0 memiliki RSE besar meskipun galat absolutnya kecil;
dalam situasi ini, yang dipublikasikan dapat berupa kelas beserta
peluangnya:

``` r

flags <- add_reliability_flags(fit_ratio, classes = c(0.01, 0.02, 0.03), draws = fit_ratio$draws)
head(flags[, c("domain", "hb", "rse", "est_class", "p_class")])
#>   domain         hb      rse    est_class p_class
#> 1      1 0.01082259 19.87851 0.01 - <0.02  0.6590
#> 2      2 0.02623409 18.05624 0.02 - <0.03  0.7020
#> 3      3 0.01018307 20.10055 0.01 - <0.02  0.5055
#> 4      4 0.04154269 18.28158      >= 0.03  0.9450
#> 5      5 0.02305273 17.98405 0.02 - <0.03  0.7135
#> 6      6 0.01270608 20.53951 0.01 - <0.02  0.8420
```

Model pangsa juga dapat diestimasi dengan MCMC. Hasilnya lalu memuat
diagnostik konvergensi model tersebut, dan grafik *trace*, peringkat,
serta divergensi menampilkan rantainya, dengan `hb[j]` sebagai $`s_d`$
untuk baris `j`:

``` r

fit_ratio_jax <- hb_ratio(s_dir ~ x,
  data = d_r, trials = "n_eff", parent = "p0", domain = "area",
  engine = jax_engine(chains = 4)
)
fit_ratio_jax$rhat_max
fit_ratio_jax$ebfmi
autoplot(fit_ratio_jax, type = "trace")
```

------------------------------------------------------------------------

## Daftar Pustaka

- Arima, S., Datta, G. S., & Liseo, B. (2015). Bayesian estimators for
  small area models when auxiliary information is measured with error.
  *Scandinavian Journal of Statistics*, 42(2), 518–529.
- Bell, W. R., Chung, H. C., Datta, G. S., & Franco, C. (2019).
  Measurement error in small area estimation: Functional versus
  structural versus naive models. *Survey Methodology*, 45(1), 61–80.
- Betancourt, M. (2016). Diagnosing suboptimal cotangent disintegrations
  in Hamiltonian Monte Carlo. *arXiv:1604.00695*.
- De Nicolò, S., & Gardini, A. (2024). The R package tipsae: Tools for
  mapping proportions and indicators on the unit interval. *Journal of
  Statistical Software*, 108(1), 1–36.
- Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., & Gelman, A.
  (2019). Visualization in Bayesian workflow. *Journal of the Royal
  Statistical Society: Series A*, 182(2), 389–402.
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
- Ybarra, L. M. R., & Lohr, S. L. (2008). Small area estimation when
  auxiliary information is measured with error. *Biometrika*, 95(4),
  919–931.
