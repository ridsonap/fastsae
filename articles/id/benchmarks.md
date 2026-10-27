# Performance & Scalability Benchmarks

*Terjemahan dari [versi bahasa
Inggris](https://ridsonap.is-a.dev/fastsae/articles/benchmarks.md).
Penjelasan istilah ada di [Daftar
Istilah](https://ridsonap.is-a.dev/fastsae/articles/id/glosarium.md).*

## Ringkasan Eksekutif

Small Area Estimation sering melibatkan registri administrasi yang besar
atau jaringan spasial yang luas dengan ratusan hingga ribuan domain.
Implementasi tradisional di R kerap bergantung pada loop yang
diinterpretasikan, rutin Fisher scoring dalam R murni, atau alokasi
memori besar untuk matriks kovarians padat $`D \times D`$.

**fastsae** mengestimasi model frequentist dalam C++ terkompilasi
(`RcppArmadillo`; faktorisasi Cholesky sparse untuk model SAR) dan model
Bayesian dengan INLA. Pada $`n = 1.000`$ domain:

- **Fay-Herriot**: 1.250x lebih cepat daripada `sae` dan 21.071x lebih
  cepat daripada `emdi`.
- **Spatial Fay-Herriot**: 46x lebih cepat daripada `sae` dan 49x lebih
  cepat daripada `emdi`.
- **Spatio-Temporal Fay-Herriot** ($`T = 5`$, point estimate): 29x lebih
  cepat daripada `sae`.

------------------------------------------------------------------------

## Tabel Hasil Benchmark

| Model ($`n = 1.000`$ domain, 5 covariate) | `fastsae` | `sae` | `emdi` |
|:---|:--:|:--:|:--:|
| **Fay-Herriot** | **0,000639 s** (1,91 MB) | 0,799 s (39,0 MB; 1.250x lebih lambat) | 13,5 s (222 MB; 21.071x lebih lambat) |
| **Spatial Fay-Herriot** | **0,518 s** (78,7 MB) | 24,1 s (239 MB; 46x lebih lambat) | 25,2 s (310 MB; 49x lebih lambat) |
| **Spatio-Temporal FH ($`T = 5`$), point estimate** | **39,8 s** (495 MB) | 1160 s (3,60 GB; 29x lebih lambat) | \- |

Waktu adalah median dari `bench::mark()`. Baris frequentist diukur pada
26-10-2026 dengan fastsae 1.0.0, sae 1.3, dan emdi 2.2.3 pada data
simulasi yang sama (Apple M4, R versi 4.6.1 (2026-06-24);
`benchmarks/bench_frequentist.R`), masing-masing dalam proses R
tersendiri. Memori adalah kenaikan memori residen puncak proses tersebut
selama estimasi, termasuk memori yang dialokasikan oleh kode
terkompilasi (kenaikan beberapa puluh MB bervariasi antarpenjalanan dan
tidak boleh dibandingkan). ‘sae’ tidak memiliki analytical MSE untuk
model spatio-temporal, sehingga baris itu membandingkan point estimate;
dengan analytical MSE-nya, fastsae membutuhkan 39,3 s. Baris `fastsae`
Bayesian diukur ulang pada 14-10-2026
(`benchmarks/bench_bayes_fastsae.R`), sedangkan baris `tipsae` (4 chain,
masing-masing 2.000 iterasi) berasal dari penjalanan sebelumnya pada
data yang sama.

------------------------------------------------------------------------

## Penjelajah Benchmark Interaktif

Gunakan kontrol interaktif di bawah untuk membandingkan waktu eksekusi,
memori yang dialokasikan, dan iterasi per detik pada berbagai ukuran
sampel:

Fay-Herriot (FH)

Spatial FH (SFH)

Spatio-Temporal (STFH)

Execution Time

Memory Allocated

Throughput (iter/s)

Execution time (median, seconds) — FH · log scale

Speedup Factor (fastsae vs competitors) — FH

------------------------------------------------------------------------

## Wawasan Arsitektur: Mengapa fastsae Begitu Cepat?

1.  **Linear Solver C++ Terkompilasi**: Alih-alih menginterpretasikan
    loop bersarang di R, `fastsae` mengimplementasikan Fisher scoring
    dan analytical MSE dalam Armadillo C++. Model Spatial Fay-Herriot
    tidak pernah membentuk invers $`D \times D`$: dengan
    $`A = (I - \rho W')(I - \rho W)`$ dan
    $`M = A + \sigma^2 E'\Psi^{-1}E`$, $`V^{-1}`$, $`g_1, \dots, g_4`$,
    dan suku kriging diperoleh dari penyelesaian Cholesky sparse dengan
    $`A`$ dan $`M`$ (Eigen, AMD ordering).

2.  **Operasi Matriks Zero-Copy**: Alokasi memori untuk struktur antara
    yang besar ($`V`$, $`V^{-1}`$, $`P`$) dihindari atau dipakai ulang
    di seluruh iterasi Fisher, bukan dialokasikan ulang di heap.

3.  **Bootstrap Multi-Thread dengan OpenMP**: Resampling bootstrap
    berjalan paralel secara native di seluruh inti CPU memakai
    `#pragma omp parallel for`, sehingga menghindari overhead
    serialisasi proses worker R
    ([`parallel::makeCluster`](https://rdrr.io/r/parallel/makeCluster.html)
    / `foreach`).

4.  **Woodbury Identity untuk Data Panel**: Pada `eblup_stfh`, inversi
    matriks kovarians blok $`(DT \times DT)`$$`V`$ direduksi menjadi
    operasi pada blok $`D \times D`$ dan $`T \times T`$ secara terpisah
    melalui dekomposisi Kronecker dan Woodbury, sehingga bottleneck
    $`O((DT)^3)`$ berubah menjadi langkah yang skalabel
    $`O(D^3 + T^3)`$.

5.  **Integrated Nested Laplace Approximations (INLA)**: Untuk model
    hierarkis tergeneralisasi dan spatio-temporal yang kompleks
    (`hb_area`), `fastsae` memakai INLA untuk menghitung marginal
    posterior analitik langsung dari Gaussian Markov Random Field yang
    sparse, sehingga sepenuhnya melewati sampling Monte Carlo.

------------------------------------------------------------------------

## Benchmark Bayesian Spatio-Temporal: `fastsae` (INLA) vs `tipsae` (Stan MCMC)

Untuk menilai kinerja small area estimation Bayesian,
[`fastsae::hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
dibandingkan dengan
[`tipsae::fit_sae()`](https://rdrr.io/pkg/tipsae/man/fit_sae.html),
paket Stan MCMC terkini untuk model small area Beta spatio-temporal.

Pengujian memakai data panel resmi `tipsae`: **`emilia`** (38 distrik
kesehatan di Emilia-Romagna selama 5 tahun, $`N = 190`$ domain
$`\times`$ tahun) dengan matriks kontiguitas poligon spasial $`W`$ dari
**`emilia_shp`**.

| Metrik | [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html) (Stan MCMC) | [`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md) (INLA) |
|:---|:---|:---|
| **Model** | Beta, Besag ICAR + RW(1) spesifik domain | `family = "beta"`, `spatial = "besag"`, `temporal = "rw1"`, `st_interaction = "domain-specific"`; struktur sama, prior berbeda |
| **Komputasi** | Hamiltonian Monte Carlo (NUTS) | Integrated nested Laplace approximation (mode `"compact"` bawaan INLA untuk Beta) |
| **Waktu eksekusi** | 22 menit (4 chain x 10.000 iterasi dengan 2.000 warm-up, paralel; $`\hat R`$ terbesar 1.007). Setiap transisi mencapai kedalaman pohon maksimum; 4 chain x 2.000 iterasi (120 s) belum konvergen ($`\hat R`$ 1.07, dan 1.12 pada runner Linux) | 1.1 s |
| **Posterior mean vs `tipsae` yang konvergen** | Acuan | $`r = 0.999`$, selisih absolut rata-rata 0.0017 (terbesar 0.006) |
| **Posterior SD vs `tipsae`** | Acuan | rasio 0.81 sampai 1.22 (prior berbeda) |

Diukur pada 16-10-2026 (acuan `tipsae` diestimasi ulang pada 17-10-2026)
pada satu mesin (Apple silicon, R 4.6.1, INLA 26.8.7, tipsae 1.0.3);
waktu bergantung pada perangkat keras. INLA tidak memerlukan chain atau
warm-up, tetapi hasil estimasinya tetap perlu diperiksa (`convergence`,
dari pencarian modus INLA).

### Penjelajah Benchmark SAE Beta Interaktif ($`n = 30`$ sampai $`n = 1.000`$)

Gunakan kontrol interaktif di bawah untuk membandingkan runtime,
penggunaan memori, dan faktor speedup antara
[`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
dan [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html)
pada model **Standard Beta**, **Spatial Beta (Besag ICAR)**, dan
**Spatio-Temporal Beta (Besag + RW1)**:

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

### Finite Population Simulation dan Validasi Ground Truth

Untuk studi simulasi design-based yang komprehensif tentang pemulihan
parameter dan properti empiris pada pengambilan sampel survei berulang
dari finite population ($`D = 60`$ dan skala nasional $`D = 500`$),
lihat artikel khusus: [Finite Population Simulation dan
Validasi](https://ridsonap.is-a.dev/fastsae/articles/id/finite-population-validation.md).
