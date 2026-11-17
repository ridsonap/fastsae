# Tolok Ukur Kinerja dan Skalabilitas

## Ringkasan

Estimasi area kecil sering melibatkan registri administratif yang besar
atau jaringan spasial yang luas dengan ratusan hingga ribuan domain.
Implementasi tradisional di R kerap bergantung pada perulangan yang
diinterpretasikan, rutin *Fisher scoring* dalam R murni, atau alokasi
memori besar untuk matriks kovarians padat berukuran $`D \times D`$.

**fastsae** mengestimasi model frekuentis dalam C++ terkompilasi
(`RcppArmadillo`; faktorisasi Cholesky matriks jarang untuk model SAR)
dan model Bayesian dengan INLA. Pada $`n = 1.000`$ domain:

- **Fay-Herriot**: 1.250x lebih cepat daripada `sae` dan 21.071x lebih
  cepat daripada `emdi`.
- **Fay-Herriot Spasial**: 46x lebih cepat daripada `sae` dan 49x lebih
  cepat daripada `emdi`.
- **Fay-Herriot Spatio-Temporal** ($`T = 5`$, estimasi titik): 29x lebih
  cepat daripada `sae`.

------------------------------------------------------------------------

## Tabel Hasil Tolok Ukur

| Model ($`n = 1.000`$ domain, 5 kovariat) | `fastsae` | `sae` | `emdi` |
|:---|:--:|:--:|:--:|
| **Fay-Herriot** | **0,000639 detik** (1,91 MB) | 0,799 detik (39,0 MB; 1.250x lebih lambat) | 13,5 detik (222 MB; 21.071x lebih lambat) |
| **Fay-Herriot Spasial** | **0,518 detik** (78,7 MB) | 24,1 detik (239 MB; 46x lebih lambat) | 25,2 detik (310 MB; 49x lebih lambat) |
| **FH Spatio-Temporal ($`T = 5`$), estimasi titik** | **39,8 detik** (495 MB) | 1160 detik (3,60 GB; 29x lebih lambat) | \- |

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
sebelumnya pada data yang sama.

------------------------------------------------------------------------

## Penjelajah Tolok Ukur Interaktif

Gunakan tombol-tombol di bawah untuk membandingkan waktu eksekusi,
memori yang dialokasikan, dan banyaknya iterasi per detik untuk berbagai
ukuran sampel:

Fay-Herriot (FH)

FH Spasial (SFH)

Spatio-Temporal (STFH)

Waktu Eksekusi

Memori yang Dialokasikan

Laju (iterasi/detik)

Waktu eksekusi (median, detik) — FH · skala log

Faktor Percepatan (fastsae vs paket lain) — FH

------------------------------------------------------------------------

## Di Balik Kecepatannya: Mengapa fastsae Begitu Cepat?

1.  **Pemecah sistem linear dalam C++ terkompilasi**: Alih-alih
    menginterpretasikan perulangan bersarang di R, `fastsae`
    mengimplementasikan *Fisher scoring* dan MSE analitik dengan
    Armadillo di C++. Model Fay-Herriot spasial tidak pernah membentuk
    invers berukuran $`D \times D`$: dengan
    $`A = (I - \rho W')(I - \rho W)`$ dan
    $`M = A + \sigma^2 E'\Psi^{-1}E`$, nilai $`V^{-1}`$,
    $`g_1, \dots, g_4`$, dan suku-suku *kriging* diperoleh dari
    penyelesaian Cholesky matriks jarang dengan $`A`$ dan $`M`$ (Eigen,
    pengurutan AMD).

2.  **Operasi matriks tanpa penyalinan**: Alokasi memori untuk struktur
    antara yang besar ($`V`$, $`V^{-1}`$, $`P`$) dihindari atau dipakai
    ulang di setiap iterasi Fisher, alih-alih dialokasikan ulang di
    *heap*.

3.  **Bootstrap *multi-thread* dengan OpenMP**: Pengambilan sampel ulang
    bootstrap dijalankan secara paralel langsung di seluruh inti CPU
    dengan `#pragma omp parallel for`, sehingga terhindar dari beban
    serialisasi proses pekerja R
    ([`parallel::makeCluster`](https://rdrr.io/r/parallel/makeCluster.html)
    / `foreach`).

4.  **Identitas Woodbury untuk data panel**: Pada `eblup_stfh`, invers
    matriks kovarians blok $`V`$ berukuran $`(DT \times DT)`$ direduksi
    menjadi operasi pada blok-blok $`D \times D`$ dan $`T \times T`$
    melalui dekomposisi Kronecker dan Woodbury, sehingga hambatan
    $`O((DT)^3)`$ berubah menjadi langkah-langkah $`O(D^3 + T^3)`$ yang
    dapat diskalakan.

5.  ***Integrated Nested Laplace Approximations* (INLA)**: Untuk model
    hierarkis terampat dan model spatio-temporal yang kompleks
    (`hb_area`), `fastsae` memakai INLA untuk menghitung distribusi
    marginal posterior secara analitik langsung dari medan acak Markov
    Gaussian yang jarang, tanpa penarikan sampel Monte Carlo sama
    sekali.

------------------------------------------------------------------------

## Tolok Ukur Spatio-Temporal Bayesian: `fastsae` (INLA) vs `tipsae` (MCMC Stan)

Untuk menilai kinerja estimasi area kecil Bayesian,
[`fastsae::hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
dibandingkan dengan
[`tipsae::fit_sae()`](https://rdrr.io/pkg/tipsae/man/fit_sae.html),
paket MCMC Stan terkini untuk model area kecil Beta spatio-temporal.

Pengujian memakai data panel resmi dari `tipsae`: **`emilia`** (38
distrik kesehatan di Emilia-Romagna selama 5 tahun, $`N = 190`$ domain
$`\times`$ tahun) dengan matriks persinggungan poligon spasial $`W`$
dari **`emilia_shp`**.

| Metrik | [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html) (MCMC Stan) | [`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md) (INLA) |
|:---|:---|:---|
| **Model** | Beta, Besag ICAR + RW(1) khusus per domain | `family = "beta"`, `spatial = "besag"`, `temporal = "rw1"`, `st_interaction = "domain-specific"`; struktur sama, prior berbeda |
| **Komputasi** | *Hamiltonian Monte Carlo* (NUTS) | *Integrated nested Laplace approximation* (mode `"compact"` bawaan INLA untuk Beta) |
| **Waktu eksekusi** | 22 menit (4 rantai x 10.000 iterasi dengan 2.000 iterasi *warm-up*, dijalankan paralel; $`\hat R`$ terbesar 1,007). Setiap transisi mencapai kedalaman pohon maksimum; dengan 4 rantai x 2.000 iterasi (120 detik), rantainya belum konvergen ($`\hat R`$ 1,07, dan 1,12 di mesin Linux) | 1,1 detik |
| **Rata-rata posterior vs `tipsae` yang sudah konvergen** | Rujukan | $`r = 0{,}999`$, rata-rata selisih absolut 0,0017 (terbesar 0,006) |
| **Simpangan baku posterior vs `tipsae`** | Rujukan | rasio 0,81 hingga 1,22 (priornya berbeda) |

Diukur pada 2026-10-16 (rujukan `tipsae` diestimasi ulang pada
2026-10-17) di satu mesin (Apple silicon, R 4.6.1, INLA 26.8.7, tipsae
1.0.3); waktu eksekusi bergantung pada perangkat kerasnya. INLA tidak
memerlukan rantai maupun *warm-up*, tetapi hasilnya tetap perlu
diperiksa (`convergence`, dari pencarian modus INLA).

### Penjelajah Tolok Ukur SAE Beta Interaktif ($`n = 30`$ hingga $`n = 1.000`$)

Gunakan tombol-tombol di bawah untuk membandingkan waktu eksekusi,
kebutuhan memori, dan faktor percepatan antara
[`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
dan [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html)
untuk model **Beta Standar**, **Beta Spasial (Besag ICAR)**, dan **Beta
Spatio-Temporal (Besag + RW1)**:

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

### Simulasi Populasi Terhingga dan Validasi terhadap Nilai Sebenarnya

Studi simulasi berbasis desain yang menilai pemulihan parameter dan
sifat-sifat empiris di bawah penarikan sampel survei berulang dari
populasi terhingga ($`D = 60`$ dan skala nasional $`D = 500`$) dibahas
di artikel tersendiri: [Simulasi dan Validasi Populasi
Terhingga](https://ridsonap.is-a.dev/fastsae/articles/id/finite-population-validation.md).
