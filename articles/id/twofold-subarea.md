# Sub-Area Level Estimation with Two-Fold Models

*Terjemahan dari [versi bahasa
Inggris](https://ridsonap.is-a.dev/fastsae/articles/twofold-subarea.md).
Penjelasan istilah ada di [Daftar
Istilah](https://ridsonap.is-a.dev/fastsae/articles/id/glosarium.md).*

## Gambaran Umum

Survei sering dirancang untuk area yang besar (provinsi, wilayah),
padahal estimasinya dibutuhkan untuk sub-area di dalamnya (kabupaten,
kota). Sub-area dari satu area memiliki karakteristik area tersebut,
sehingga estimasinya semestinya meminjam informasi tidak hanya dari
regresi pada auxiliary variable, tetapi juga dari sub-area lain dalam
area yang sama.

Sub-area level model two-fold dari Torabi dan Rao (2014) melakukannya
dengan dua random effect bersarang. **fastsae** mengestimasinya dengan
dua cara:

| Fungsi | Pendekatan | Ketidakpastian | Agregat area |
|:---|:---|:---|:---|
| [`eblup_twofold()`](https://ridsonap.is-a.dev/fastsae/reference/eblup_twofold.md) | EBLUP, REML atau ML | MSE analitik (tipe Prasad-Rao) atau parametric bootstrap | `df_area`, dibobot dengan `weight` |
| [`hb_twofold()`](https://ridsonap.is-a.dev/fastsae/reference/hb_twofold.md) | Hierarchical Bayes dengan INLA | Posterior SD dan credible interval | `df_area`, dibobot dengan `weight` |

[`hb_twofold()`](https://ridsonap.is-a.dev/fastsae/reference/hb_twofold.md)
juga mengestimasi respons binomial dan Poisson serta efek area spasial
(BYM2, Besag).

------------------------------------------------------------------------

## Model

Untuk sub-area $`j`$ ($`j = 1, \dots, N_i`$) dari area $`i`$
($`i = 1, \dots, m`$), direct estimator $`y_{ij}`$ dengan sampling
variance $`\psi_{ij}`$ yang diketahui mengikuti

``` math
y_{ij} = \theta_{ij} + e_{ij}, \qquad \theta_{ij} = x_{ij}^\top \beta + v_i + u_{ij},
```

dengan $`v_i \sim N(0, \sigma_v^2)`$ adalah efek area,
$`u_{ij} \sim N(0, \sigma_u^2)`$ adalah efek sub-area, dan
$`e_{ij} \sim N(0, \psi_{ij})`$ adalah sampling error, semuanya saling
bebas.

Dibandingkan dengan model Fay-Herriot, yang memperlakukan setiap
sub-area sebagai tidak berhubungan, efek area $`v_i`$ diestimasi dari
semua sub-area tersampel di area $`i`$. Hal ini paling berarti untuk
sub-area yang tidak tersampel: mereka memperoleh
$`x_{ij}^\top \hat\beta + \hat v_i`$, bukan hanya synthetic estimator
$`x_{ij}^\top \hat\beta`$. Hanya sub-area dari area yang sama sekali
tidak memiliki sampel yang diprediksi secara sintetik.

Keluaran mengikuti notasi di atas: `random_effect_var` memuat `sigma2_v`
(area) dan `sigma2_u` (sub-area), sedangkan `df_eblup` (atau `df_hb`)
memuat `random_effect_area` dan `random_effect_subarea`.

------------------------------------------------------------------------

## Contoh dengan Data Simulasi

Kita mensimulasikan 30 area dengan masing-masing 3 hingga 6 sub-area,
sehingga rata-rata sub-area sebenarnya $`\theta_{ij}`$ diketahui. Lima
belas sub-area tidak memiliki sampel (direct estimate-nya `NA`), dan
semua sub-area dari area `A30` tidak tersampel.

``` r

library(fastsae)

set.seed(2024)
m <- 30
dat <- do.call(rbind, lapply(seq_len(m), function(i) {
  ni <- sample(3:6, 1)
  data.frame(
    area = sprintf("A%02d", i),
    subarea = sprintf("A%02d-%d", i, seq_len(ni)),
    x1 = rnorm(ni, mean = 2),
    pop = sample(200:2000, ni)
  )
}))
v <- rnorm(m, sd = 1)
names(v) <- unique(dat$area)
dat$theta <- 10 + 2 * dat$x1 + v[dat$area] + rnorm(nrow(dat), sd = sqrt(0.5))
dat$vardir <- runif(nrow(dat), 0.5, 2)
dat$y <- dat$theta + rnorm(nrow(dat), sd = sqrt(dat$vardir))

# Unsampled sub-areas: 15 at random, plus every sub-area of A30
unsampled <- union(sample(which(dat$area != "A30"), 15), which(dat$area == "A30"))
dat$y[unsampled] <- NA

c(subareas = nrow(dat), sampled = sum(!is.na(dat$y)))
#> subareas  sampled 
#>      131      111
head(dat)
#>   area subarea        x1  pop    theta   vardir        y
#> 1  A01   A01-1 1.5347257 1857 13.78128 1.663801 12.20118
#> 2  A01   A01-2 2.5191534 1961 15.03099 1.105026 14.48113
#> 3  A01   A01-3 2.5284898 1222 14.79953 1.301555 15.01471
#> 4  A01   A01-4 1.4848072 1252 11.84648 1.406429 11.36094
#> 5  A02   A02-1 1.8729653 1255 13.61192 1.986447 12.02304
#> 6  A02   A02-2 0.7751246  924 12.23727 1.268523       NA
```

Setiap baris adalah satu sub-area. Sub-area yang tidak tersampel tetap
ada di data dengan `y = NA` beserta covariate-nya, sehingga diprediksi.

### Mengestimasi EBLUP

``` r

fit_tf <- eblup_twofold(
  y ~ x1,
  vardir = "vardir",
  domain = "area",
  subarea = "subarea",
  data = dat,
  weight = "pop",
  print_result = FALSE
)

fit_tf$random_effect_var
#>  sigma2_v  sigma2_u 
#> 1.1067031 0.1486881
fit_tf$estcoef
#>                 beta std.error stderr_beta   zvalue        pvalue
#> (Intercept) 9.605159 0.3223129   0.3223129 29.80073 3.822665e-195
#> x1          2.104314 0.1144109   0.1144109 18.39260  1.505751e-75
head(fit_tf$df_eblup)
#>   domain subarea        y    eblup   vardir random_effect_area
#> 1    A01   A01-1 12.20118 12.39678 1.663801         -0.4204473
#> 2    A01   A01-2 14.48113 14.48525 1.105026         -0.4204473
#> 3    A01   A01-3 15.01471 14.55766 1.301555         -0.4204473
#> 4    A01   A01-4 11.36094 12.21855 1.406429         -0.4204473
#> 5    A02   A02-1 12.02304 12.66335 1.986447         -0.8351927
#> 6    A02   A02-2       NA 10.40107 1.268523         -0.8351927
#>   random_effect_subarea       mse      rse
#> 1         -0.0174796591 0.4205889 5.231427
#> 2         -0.0005543296 0.3941804 4.334326
#> 3          0.0522123369 0.4063272 4.378712
#> 4         -0.0906664052 0.4122688 5.254976
#> 5         -0.0479275040 0.4428309 5.254973
#> 6          0.0000000000 0.4782345 6.648787
```

Bawaan `mse_type = "analytical"` menghasilkan MSE tipe Prasad-Rao dari
model REML (untuk ML ditambahkan koreksi bias). Sub-area `A01-3` tidak
tersampel: ia tetap memperoleh efek area dari `A01`, sedangkan efek
sub-area-nya 0.

`weight` (di sini populasi setiap sub-area) mengagregasi EBLUP sub-area
menjadi rata-rata area,
$`\tilde\theta_i = \sum_j w_{ij}\tilde\theta_{ij} / \sum_j w_{ij}`$, di
`df_area`. MSE-nya memperhitungkan error yang dimiliki bersama oleh
sub-area dalam satu area (efek area yang diestimasi dan $`\hat\beta`$),
sehingga lebih besar daripada yang disarankan oleh
$`\sum_j w_{ij}^2\,\mathrm{mse}_{ij}`$.

``` r

head(fit_tf$df_area)
#>   domain eblup_area  mse_area rse_area n_subareas n_sampled
#> 1    A01   13.43189 0.2748473 3.903090          4         4
#> 2    A02   11.19056 0.3146821 5.012844          5         4
#> 3    A03   14.07817 0.2208185 3.337885          3         3
#> 4    A04   16.64890 0.3174408 3.384120          6         4
#> 5    A05   13.25981 0.4443843 5.027383          4         3
#> 6    A06   12.79993 0.2838182 4.162101          3         3
```

### Seberapa Besar Manfaat Efek Area?

Karena datanya simulasi, kita dapat membandingkan setiap estimator
dengan $`\theta_{ij}`$ sebenarnya. Model Fay-Herriot, yang diestimasi
dengan
[`eblup_fh()`](https://ridsonap.is-a.dev/fastsae/reference/eblup_fh.md)
pada baris yang sama, mengabaikan area.

``` r

fit_fh <- eblup_fh(y ~ x1, vardir = "vardir", data = dat, print_result = FALSE)

sampled <- !is.na(dat$y)
in_sampled_area <- !sampled & dat$area != "A30"
emp_mse <- function(est, rows) mean((est[rows] - dat$theta[rows])^2)

data.frame(
  rows = c("sampled sub-areas", "unsampled sub-areas of sampled areas"),
  direct = c(emp_mse(dat$y, sampled), NA),
  fay_herriot = c(
    emp_mse(fit_fh$df_eblup$eblup, sampled),
    emp_mse(fit_fh$df_eblup$eblup, in_sampled_area)
  ),
  two_fold = c(
    emp_mse(fit_tf$df_eblup$eblup, sampled),
    emp_mse(fit_tf$df_eblup$eblup, in_sampled_area)
  )
)
#>                                   rows   direct fay_herriot  two_fold
#> 1                    sampled sub-areas 1.211463   0.5087563 0.5291281
#> 2 unsampled sub-areas of sampled areas       NA   1.5247133 0.7748664
```

Kedua model memperbaiki direct estimate pada sub-area tersampel, dan
pada data ini keduanya kurang lebih sama akuratnya di sana: sub-area
tersampel sudah punya direct estimate sendiri, sehingga efek area hanya
menambah sedikit. Perbedaannya ada pada sub-area tidak tersampel di area
yang tersampel. Model Fay-Herriot hanya dapat memberinya synthetic
estimate, sedangkan model two-fold menambahkan estimasi efek area dan
menurunkan error sekitar separuhnya. Ini hanya satu data simulasi, bukan
studi simulasi: besar manfaatnya bergantung pada seberapa besar
$`\sigma_v^2`$ relatif terhadap $`\sigma_u^2`$.

### Bootstrap MSE

Parametric bootstrap mengestimasi ulang model pada `B` data yang
dibangkitkan dari model hasil estimasi, lalu memprediksi dengan
estimator yang sama seperti estimasi titik. `seed` membuatnya dapat
direproduksi tanpa mengubah aliran bilangan acak sesi.

``` r

fit_boot <- eblup_twofold(
  y ~ x1,
  vardir = "vardir",
  domain = "area",
  subarea = "subarea",
  data = dat,
  mse_type = "bootstrap",
  B = 100,
  seed = 1,
  print_result = FALSE
)

summary(fit_boot$df_eblup$mse / fit_tf$df_eblup$mse)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.6322  0.8257  0.9056  0.9205  0.9971  1.4684
```

Kedua estimator MSE sesuai secara rata-rata; dengan `B = 100`, nilai
bootstrap tunggal masih bising. Gunakan `mse_type = "none"` bila hanya
estimasi titik yang dibutuhkan.

### Diagnostik dan Plot

Perangkat bawaan paket dapat dipakai pada hasilnya:

``` r

autoplot(fit_tf, type = "comparison")
```

![](twofold-subarea_files/figure-html/diagnostics-1.png)

``` r

diagnose(fit_tf)
#> ── fastsae Small Area Estimation Diagnostic Report ─────────────────────────────
#> Model: "TWOFOLD" | Domains: 131 (Sampled: 111, Unsampled: 20)
#> 
#> ── 1. Precision & Efficiency Gain ──
#> 
#> ✔ Domains with RSE < 25%: 100% (High precision)
#> ✔ Average RSE reduction: Direct "8.54%" -> SAE "4.97%" (Gain: 3.57%)
#> ✔ Variance reduction in 100% of areas (MSE ratio median: 3.2, max: 5.34)
#> 
#> ── 2. Brown et al. (2001) Diagnostics ──
#> 
#> ✔ Goodness-of-Fit Statistic W: W = 60.62 (df = 111), p-value = 1 [Consistent with direct estimates]
#> ℹ Bias regression (direct ~ model, descriptive): alpha = -0.7789 (SE 0.5383), beta = 1.0544 (SE 0.0388), R^2 = 0.872
#> A slope above 1 is expected for shrinkage estimators; see ?diagnose.
#> ────────────────────────────────────────────────────────────────────────────────
#> ✔ Final Assessment: PASS: Model estimates are consistent with the direct estimates (Brown W test) and reliable.
#> ────────────────────────────────────────────────────────────────────────────────
```

------------------------------------------------------------------------

## Hierarchical Bayes dengan `hb_twofold()`

[`hb_twofold()`](https://ridsonap.is-a.dev/fastsae/reference/hb_twofold.md)
mengestimasi model yang sama dengan INLA. Fungsi ini membutuhkan paket
**INLA**, yang tidak ada di CRAN; lihat
<https://www.r-inla.org/download/>. Chunk di bawah hanya dijalankan bila
INLA terpasang.

Argumen `weight` (di sini populasi setiap sub-area) mengagregasi
posterior sub-area menjadi estimasi area,
$`\hat\theta_i = \sum_j w_{ij} \theta_{ij} / \sum_j w_{ij}`$, beserta
ketidakpastian posterior-nya di `df_area`. Pada model Gaussian, ini
adalah kombinasi linear dari prediktor sub-area, yang posteriornya
dihitung INLA secara eksak; agregat binomial dan Poisson berasal dari
posterior draw (`inla_engine(draws = )`).

``` r

fit_hb <- hb_twofold(
  y ~ x1,
  vardir = "vardir",
  domain = "area",
  subarea = "subarea",
  weight = "pop",
  data = dat,
  print_result = FALSE
)

fit_hb$random_effect_var
#>  sigma2_v  sigma2_u 
#> 1.2370025 0.1859788
head(fit_hb$df_hb)
#>   domain subarea        y       hb linear_pred   vardir        sd       mse
#> 1    A01   A01-1 12.20118 12.39290    12.39290 1.663801 0.6267376 0.3928000
#> 2    A01   A01-2 14.48113 14.48306    14.48306 1.105026 0.6022376 0.3626901
#> 3    A01   A01-3 15.01471 14.56061    14.56061 1.301555 0.6149109 0.3781155
#> 4    A01   A01-4 11.36094 12.20725    12.20725 1.406429 0.6221921 0.3871231
#> 5    A02   A02-1 12.02304 12.65286    12.65286 1.986447 0.6498927 0.4223605
#> 6    A02   A02-2       NA 10.39838    10.39838 1.268523 0.7105738 0.5049151
#>        rse  ci_lower ci_upper random_effect_area random_effect_subarea
#> 1 5.057232 11.151455 13.62611         -0.4211314         -2.072824e-02
#> 2 4.158220 13.296659 15.66950         -0.4211314         -1.332579e-03
#> 3 4.223113 13.363453 15.78882         -0.4211314          5.683072e-02
#> 4 5.096906 10.953783 13.40860         -0.4211314         -9.938372e-02
#> 5 5.136332 11.353954 13.91986         -0.8380897         -5.535978e-02
#> 6 6.833501  8.993467 11.81479         -0.8380897          6.810325e-11
head(fit_hb$df_area)
#>   domain  hb_area   sd_area  mse_area rse_area ci_lower_area ci_upper_area
#> 1    A01 13.42839 0.5161961 0.2664584 3.844065      12.41291      14.43824
#> 2    A02 11.18134 0.5548966 0.3079103 4.962701      10.08431      12.26192
#> 3    A03 14.08162 0.4637612 0.2150745 3.293380      13.17467      14.99399
#> 4    A04 16.66043 0.5773106 0.3332875 3.465160      15.53871      17.80395
#> 5    A05 13.24769 0.6585018 0.4336246 4.970692      11.94804      14.53387
#> 6    A06 12.79910 0.5224836 0.2729891 4.082192      11.77397      13.82418
#>   n_subareas
#> 1          4
#> 2          5
#> 3          3
#> 4          6
#> 5          4
#> 6          3
```

Dengan prior bawaan, posterior mean sangat dekat dengan EBLUP (selisih
rata-rata, limits of agreement, dan Lin’s concordance correlation di
bawah):

``` r

compare_sae(fit_tf, fit_hb, names = c("EBLUP two-fold", "HB two-fold"))
#> 
#> ── Small Area Estimation Model Concordance & Comparison ────────────────────────
#> Model 1: EBLUP two-fold
#> Model 2: HB two-fold
#> Matched Domains: 131
#> ────────────────────────────────────────────────────────────────────────────────
#> 
#> Agreement, Ranking & Efficiency Metrics:
#>  Group       Metric                                                  
#>  Domains     Matched Domains                                         
#>  Agreement   Mean Difference (HB two-fold - EBLUP two-fold)          
#>  Agreement   95% Limits of Agreement                                 
#>  Agreement   Mean Absolute Difference (MAE)                          
#>  Agreement   Root Mean Squared Difference (RMSD)                     
#>  Agreement   Lin's Concordance Correlation (CCC)                     
#>  Uncertainty Median |Difference| / RMSE of EBLUP two-fold            
#>  Uncertainty % of domains with |Difference| > RMSE of EBLUP two-fold 
#>  Ranking     Spearman Rank Correlation (rho)                         
#>  Ranking     Mean Absolute Rank Change                               
#>  Ranking     Shared Top 13 (highest estimates)                       
#>  Ranking     Shared Bottom 13 (lowest estimates)                     
#>  Efficiency  Mean MSE Ratio (EBLUP two-fold / HB two-fold)           
#>  Efficiency  Median MSE Ratio (EBLUP two-fold / HB two-fold)         
#>  Efficiency  % of domains where HB two-fold reports lower uncertainty
#>  Efficiency  Mean RSE % (EBLUP two-fold)                             
#>  Efficiency  Mean RSE % (HB two-fold)                                
#>  Value            
#>  131              
#>  -0.0004          
#>  [-0.0209, 0.0201]
#>  0.0077           
#>  0.0104           
#>  1.0000           
#>  0.008            
#>  0.0%             
#>  0.9999           
#>  0.18             
#>  13 / 13          
#>  13 / 13          
#>  1.0438           
#>  1.0526           
#>  78.6%            
#>  4.97%            
#>  4.89%
#> ℹ One model is Hierarchical Bayes: its posterior variance and the frequentist MSE of the other are not the same quantity, so a smaller value does not show a more accurate estimator.
```

[`compare_sae()`](https://ridsonap.is-a.dev/fastsae/reference/compare_sae.md)
melaporkan varians posterior yang lebih kecil daripada MSE EBLUP.
Keduanya bukan besaran yang sama, dan di sini posterior mean
$`\sigma_u^2`$ juga lebih kecil daripada estimasi REML (keduanya di
bawah nilai sebenarnya 0,5, yang sulit dipisahkan dari sampling variance
sebesar 0,5 hingga 2). Karena itu, ketidakpastian yang dilaporkan lebih
kecil bukanlah bukti estimator yang lebih akurat. Pada data simulasi,
error terhadap nilai sebenarnya yang menentukan; di sini keduanya
praktis sama:

``` r

c(
  eblup_two_fold = emp_mse(fit_tf$df_eblup$eblup, sampled | in_sampled_area),
  hb_two_fold = emp_mse(fit_hb$df_hb$hb, sampled | in_sampled_area)
)
#> eblup_two_fold    hb_two_fold 
#>      0.5583827      0.5570989
```

Untuk data cacah atau proporsi, gunakan `family = "poisson"` dengan
`exposure` atau `family = "binomial"` dengan `trials`; dengan
`spatial = "bym2"` dan neighbourhood graph `W` antar-area (misalnya dari
[`create_weights()`](https://ridsonap.is-a.dev/fastsae/reference/create_weights.md)),
efek area berkorelasi spasial.

### MCMC dengan `engine = jax_engine()`

Model yang sama dapat disampling dengan MCMC (No-U-Turn sampler NumPyro
di atas JAX) sebagai pengganti INLA. Ini membutuhkan Python dengan `jax`
dan `numpyro` (lihat
[`?setup_numpyro_env`](https://ridsonap.is-a.dev/fastsae/reference/setup_numpyro_env.md)),
sehingga chunk ini tidak dijalankan di sini. Hasilnya berstruktur sama,
ditambah posterior draw di `draws` dan split R-hat terbesar di
`rhat_max`; agregat area di `df_area` berasal dari draw MCMC.

``` r

fit_jax <- hb_twofold(
  y ~ x1,
  vardir = "vardir",
  domain = "area",
  subarea = "subarea",
  weight = "pop",
  data = dat,
  engine = jax_engine(warmup = 1000, samples = 1000, chains = 4),
  print_result = FALSE
)
fit_jax$rhat_max
head(fit_jax$df_area)
```

Kedua engine memakai prior yang sama secara bawaan
(`jax_engine(priors = "inla")`), sehingga estimasinya hanya berbeda
karena aproksimasi INLA dan Monte Carlo error.

------------------------------------------------------------------------

## Catatan Praktis

- Berikan satu baris per sub-area. Kolom `domain` mengidentifikasi area
  dan wajib ada.
- `vardir` harus berhingga dan positif untuk setiap sub-area tersampel;
  sub-area tidak tersampel boleh `NA`. Covariate harus lengkap untuk
  semua baris, tersampel maupun tidak.
- Minimal dua area harus memiliki sub-area tersampel, jika tidak
  $`\sigma_v^2`$ tidak dapat diestimasi.
- [`eblup_twofold()`](https://ridsonap.is-a.dev/fastsae/reference/eblup_twofold.md)
  memiliki `self_benchmark = TRUE` agar estimasi menjumlah ke total yang
  diketahui (lihat
  [`?eblup_twofold`](https://ridsonap.is-a.dev/fastsae/reference/eblup_twofold.md)).

## Daftar Pustaka

- Torabi, M., & Rao, J. N. K. (2014). On small area estimation under a
  sub-area level model. *Journal of Multivariate Analysis*, 127, 36–55.
- Rao, J. N. K., & Molina, I. (2015). *Small Area Estimation* (2nd ed.).
  John Wiley & Sons.
- Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian
  inference for latent Gaussian models by using integrated nested
  Laplace approximations. *Journal of the Royal Statistical Society:
  Series B*, 71(2), 319–392.
