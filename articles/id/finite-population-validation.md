# Simulasi Populasi Terhingga: Validasi fastsae terhadap Nilai Sebenarnya

## Pendahuluan: Mengapa Simulasi Populasi Terhingga?

Dalam Estimasi Area Kecil (SAE), studi simulasi adalah tolok ukur utama
untuk memvalidasi model statistik sebelum dipakai dalam statistik resmi
(misalnya pemetaan kemiskinan atau indikator kesehatan daerah). Ada dua
paradigma pokok:

1.  **Simulasi berbasis model**: Data dibangkitkan langsung dari
    distribusi parametrik yang diasumsikan model (misalnya menarik
    $`y_d \sim \text{Beta}(\dots)`$). Cara ini berguna untuk memastikan
    algoritmenya benar dan parameternya dapat dipulihkan, tetapi model
    hanya dinilai di bawah asumsi idealnya sendiri.
2.  **Simulasi berbasis desain pada populasi terhingga (*standar
    emas*)**: Dibangkitkan sebuah populasi terhingga yang realistis
    berupa individu atau rumah tangga diskret ($`N = 33.511`$ di
    $`D = 60`$ domain). Besaran tingkat area yang sebenarnya, $`P_d`$,
    dihitung langsung dengan mencacah seluruh unit di setiap domain:
    ``` math
    P_d = \frac{1}{N_d} \sum_{i=1}^{N_d} y_{di} \quad (\text{nilai sebenarnya yang eksak dan nonparametrik})
    ```
    Sampel survei berukuran kecil ($`n_d \ll N_d`$, misalnya
    $`n_d = 20 - 55`$) ditarik berulang kali. Estimator langsung dan
    varians sampling berbasis desain dihitung dari sampel survei, lalu
    dimasukkan ke model area kecil.

Simulasi berbasis desain ini menguji secara ketat apakah:

- model Beta kontinu tingkat area (`hb_area` dan `tipsae`) tetap tangguh
  bila data dasarnya berupa respons biner individual yang diskret
  ($`y_{di} \in \{0, 1\}`$);
- varians sampling berbasis desain $`\widehat{\text{Var}}(\hat{p}_d)`$
  dengan koreksi populasi terhingga ditangani dengan tepat oleh model
  penghubungnya (*linking model*);
- selang kredibel model tetap mendekati cakupan nominal 95% di bawah
  penarikan sampel survei berulang.

------------------------------------------------------------------------

## Rancangan Simulasi

Rancangan simulasi ini meniru survei rumah tangga nasional (seperti
Susenas atau EU-SILC):

    +-------------------------------------------------------------------------------+
    | 1. POPULASI TERHINGGA SINTETIS (N = 33.511 unit, D = 60 domain)               |
    |    - Kovariat unit: x_1di ~ N(1.6, 0.6^2), x_2di ~ Unif(0, 2.2)               |
    |    - Ketergantungan spasial: efek acak Besag ICAR (matriks ketetanggaan W)    |
    |    - Respons unit: y_di ~ Bernoulli(pi_di) dalam {0, 1}                       |
    |    - Nilai sebenarnya: P_d = sum(y_di) / N_d (rata-rata P_d = 27,4%)          |
    +-------------------------------------------------------------------------------+
                                          |
                                          v
    +-------------------------------------------------------------------------------+
    | 2. PENARIKAN SAMPEL SURVEI (R = 100 ulangan Monte Carlo)                      |
    |    - SRS tanpa pengembalian per domain (n_d dalam [22, 55], rata-rata 40)     |
    |    - Fraksi sampling f_d = n_d / N_d ~ 7,1%                                   |
    |    - Estimator langsung: p_hat_d = mean(y_sample_d)                           |
    |    - Varians sampling desain: V_hat_d = (1 - f_d) * p(1 - p) / (n_d - 1)      |
    +-------------------------------------------------------------------------------+
                                          |
                       +------------------+------------------+
                       |                                     |
                       v                                     v
    +------------------------------------+ +------------------------------------+
    | 3. fastsae::hb_area()              | | 4. tipsae::fit_sae()               |
    |    - INLA Bayesian (Laplace penuh) | |    - Hamiltonian Monte Carlo (Stan)|
    |    - Keluarga: "beta"              | |    - Likelihood: "beta"            |
    |    - Spasial: "besag"              | |    - Galat spasial: TRUE           |
    |    - Strategi: "laplace"           | |    - Waktu: ~12-40 detik / estimasi|
    |    - Prior bawaan (hb-priors)      | |                                    |
    |    - Waktu: ~6-9 detik / estimasi  | |                                    |
    +------------------------------------+ +------------------------------------+
                       |                                     |
                       +------------------+------------------+
                                          |
                                          v
    +-------------------------------------------------------------------------------+
    | 5. EVALUASI EMPIRIS TERHADAP NILAI SEBENARNYA P_d                             |
    |    - RMSE empiris, RRMSE empiris (%), dan bias relatif empiris (%)            |
    |    - Tingkat cakupan selang kredibel 95% pada 100 ulangan                     |
    |    - Peningkatan efisiensi vs estimator langsung (rasio MSE)                  |
    +-------------------------------------------------------------------------------+

------------------------------------------------------------------------

## Formulasi Matematis Sifat-Sifat Empiris

Untuk $`R = 100`$ ulangan sampel survei yang saling bebas, sifat-sifat
empiris untuk setiap domain $`d = 1, \dots, D`$ didefinisikan sebagai
berikut:

1.  **Rata-rata dan bias empiris**:
    ``` math
    \overline{\hat{\theta}}_d = \frac{1}{R}\sum_{r=1}^R \hat{\theta}_d^{(r)}, \quad \text{Bias}_d = \overline{\hat{\theta}}_d - P_d, \quad \text{RB}_d = \frac{\text{Bias}_d}{P_d} \times 100\%
    ```

2.  ***Root mean squared error* (RMSE) empiris**:
    ``` math
    \text{RMSE}_d = \sqrt{\frac{1}{R}\sum_{r=1}^R \left(\hat{\theta}_d^{(r)} - P_d\right)^2}, \quad \text{RRMSE}_d = \frac{\text{RMSE}_d}{P_d} \times 100\%
    ```

3.  **Tingkat cakupan (CR) empiris selang kredibel 95%**:
    ``` math
    \text{CR}_d = \frac{1}{R}\sum_{r=1}^R \mathbb{I}\left(P_d \in \left[\text{CI}_{\text{lower}, d}^{(r)}, \text{CI}_{\text{upper}, d}^{(r)}\right]\right) \times 100\%
    ```

4.  **Peningkatan efisiensi empiris dibandingkan estimator langsung**:
    ``` math
    \text{Eff}_d = \frac{\text{MSE}_d(\text{Direct})}{\text{MSE}_d(\text{Model})}
    ```

------------------------------------------------------------------------

## Hasil Tolok Ukur Empiris ($`R = 100`$ Ulangan, $`D = 60`$)

Tabel di bawah merangkum rata-rata makro metrik empiris di seluruh 60
domain dari simulasi Monte Carlo lengkap dengan $`R = 100`$. Yang
dinilai adalah
[`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
dengan *engine* dan prior bawaannya
([`inla_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md);
lihat [artikel prior
bawaan](https://ridsonap.is-a.dev/fastsae/articles/id/hb-priors.md))
serta [`tipsae::fit_sae()`](https://rdrr.io/pkg/tipsae/man/fit_sae.html)
dengan `type_disp = "var"` (satu rantai 600 iterasi). Waktu eksekusi
dihitung per estimasi, dengan lima ulangan dijalankan paralel.

*(Sel yang disorot hijau menandai kinerja terbaik untuk setiap kriteria
penilaian.)*

| Estimator / Model | RMSE empiris | RRMSE empiris (%) | MARB empiris (%) | RB empiris (%) | Cakupan empiris selang kredibel 95% | Rata-rata lebar selang kredibel | Peningkatan efisiensi (vs langsung) | Waktu (detik/estimasi) |
|:---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| **Estimator langsung (sampel $`\hat{p}_d`$)** | 0,0681 | 26,27% | 20,71% | +2,27% | — | — | 1,000x | — |
| **[`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md) (Besag spasial)** | 0,04913 | 19,31% | 16,09% | +2,36% | 92,13% | 0,1855 | 1,783x | 6,18 detik |
| **[`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html) (Besag spasial)** | 0,04911 | 19,16% | 15,89% | +1,06% | 93,18% | 0,1979 | 1,806x | 11,65 detik |

### Temuan Empiris Utama

1.  **Penurunan varians dan galat yang berarti**: Baik `fastsae` maupun
    `tipsae` menurunkan RMSE empiris dari **0,0681** menjadi **0,0491**,
    yaitu peningkatan efisiensi empiris sekitar **1,8x** dibandingkan
    estimasi survei langsung. Akurasi kedua estimasi Bayesian ini
    praktis sama.
2.  **Kuantifikasi ketidakpastian**: `fastsae` memberikan selang
    kredibel yang sedikit lebih sempit ($`0{,}186`$ vs $`0{,}198`$)
    dengan cakupan 92,1%, dibandingkan 93,2% untuk `tipsae`; keduanya
    sedikit di bawah nilai nominal 95%.
3.  **Skalabilitas komputasi**: Di sini
    [`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
    hanya memerlukan sekitar separuh waktu MCMC Stan; model Beta-nya
    dengan `vardir` diestimasi dua kali (lihat
    [`?hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)).

Pengujian sebelumnya untuk simulasi ini memakai prior PC dengan
`param = c(1, 0.05)` untuk simpangan baku area dan menetapkan presisi
likelihood Beta pada estimasi langsung; selangnya mencakup 88,3%
proporsi sebenarnya (RMSE 0,0485, rata-rata lebar 0,172). Nilai bawaan
yang sekarang memperlebar selang dan membawa cakupannya mendekati
cakupan `tipsae`.

------------------------------------------------------------------------

## Tolok Ukur Skala Nasional ($`D = 500`$ Domain, $`N \approx 250.000`$ Individu, $`R = 25`$)

Untuk menilai skalabilitas pada tingkat nasional (setara dengan rincian
$`\approx 514`$ kabupaten/kota di Indonesia atau ratusan *county* di
negara lain), dibangkitkan populasi sintetis nasional berisi
$`N = 250.824`$ individu pada kisi spasial $`20 \times 25`$ yang terdiri
atas $`D = 500`$ domain yang bersinggungan. Sampel survei berukuran
$`n_d \in [20, 50]`$ per domain ($`n \approx 17.500`$ unit survei)
ditarik dalam $`R = 25`$ ulangan Monte Carlo.

*(Sel yang disorot hijau menandai kinerja terbaik.)*

| Estimator / Model | RMSE empiris | RRMSE empiris (%) | MARB empiris (%) | RB empiris (%) | Cakupan selang kredibel 95% (%) | Rata-rata lebar selang kredibel | Peningkatan efisiensi (vs langsung) | Waktu (detik/estimasi) |
|:---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| **Estimator langsung (sampel $`\hat{p}_d`$)** | 0,0719 | 27,14% | 21,75% | +2,96% | — | — | 1,000x | — |
| **[`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md) (Besag spasial)** | 0,0489 | 18,85% | 15,82% | +2,81% | 93,38% | 0,1935 | 1,991x | 8,53 detik |
| **[`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html) (Besag spasial)** | 0,0487 | 18,68% | 15,66% | +2,00% | 93,69% | 0,2017 | 2,014x | 40,00 detik |

### Hasil Utama pada Skala Nasional

- **Percepatan 4,7x**: `fastsae` mengestimasi model spasial lengkap
  dengan 500 domain dalam **8,5 detik**, dibandingkan **40,0 detik**
  untuk MCMC Stan (satu rantai 500 iterasi).
- **Akurasi setara**: RMSE empiris kedua estimasi Bayesian praktis sama
  ($`0{,}0489`$ vs $`0{,}0487`$, dibandingkan $`0{,}0719`$ untuk
  estimator langsung).
- **Cakupan**: 93,4% untuk `fastsae` dan 93,7% untuk `tipsae`, dengan
  selang yang sedikit lebih sempit untuk `fastsae` ($`0{,}194`$ vs
  $`0{,}202`$). Pengujian sebelumnya dengan prior PC dan presisi yang
  ditetapkan pada estimasi langsung mencapai cakupan 89,0%.

------------------------------------------------------------------------

## Visualisasi

### 1. RMSE Empiris per Domain

Di seluruh 60 area kecil, baik `fastsae` maupun `tipsae` secara
konsisten memangkas varians estimasi menjadi sekitar separuh varians
estimasi survei langsung:

``` r

# Kode untuk mereproduksi visualisasi RMSE per domain
ggplot(df_rmse, aes(x = DomainIdx, y = RMSE, color = Method)) +
  geom_line() + geom_point() +
  labs(title = "RMSE empiris per domain (R = 100 ulangan)",
       x = "Indeks domain", y = "RMSE empiris terhadap P_d yang sebenarnya")
```

### 2. Cakupan Empiris Selang Kredibel 95%

Tingkat cakupan empiris setiap domain tersebar di sekitar garis rujukan
putus-putus pada nilai nominal 95%, dengan rata-rata 92,1% untuk
`fastsae` dan 93,2% untuk `tipsae`.

------------------------------------------------------------------------

## Cara Mereproduksi

Skrip lengkap tolok ukur populasi terhingga tersedia di repositori
paket:

``` bash
# Jalankan 100 ulangan Monte Carlo di 5 inti CPU
Rscript benchmarks/sim_mc_finite_population.R --R=100 --cores=5
```

------------------------------------------------------------------------

## Daftar Pustaka

- **Janicki, R.** (2020). Properties of the beta regression model for
  small area estimation of proportions and application to estimation of
  poverty rates. *Communications in Statistics - Theory and Methods*,
  49(9), 2264–2284.
- **De Nicolò, S., & Gardini, A.** (2024). The R Package tipsae: Tools
  for Mapping Proportions and Indicators on the Unit Interval. *Journal
  of Statistical Software*, 108(1), 1–36.
- **Rao, J. N. K., & Molina, I.** (2015). *Small Area Estimation* (2nd
  ed.). John Wiley & Sons.
- **Tzavidis, N., et al.** (2018). From start to finish: a framework for
  the production of small area official statistics. *Journal of the
  Royal Statistical Society: Series A*, 181(4), 927–979.
- **Rue, H., Martino, S., & Chopin, N.** (2009). Approximate Bayesian
  inference for latent Gaussian models by using integrated nested
  Laplace approximations. *Journal of the Royal Statistical Society:
  Series B*, 71(2), 319–392.
