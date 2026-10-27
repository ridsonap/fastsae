# Daftar Istilah

Dokumentasi berbahasa Indonesia ini **mempertahankan istilah teknis
dalam bahasa Inggris** (misalnya *small area estimation*, *direct
estimator*, *random effect*), karena itulah yang dipakai di literatur,
di dokumentasi fungsi (halaman referensi), dan di keluaran R. Kalimat
penjelasnya ditulis dalam bahasa Indonesia. Nama fungsi, argumen, dan
kode tidak diterjemahkan dan sama persis dengan dokumentasi berbahasa
Inggris.

Halaman ini memberi penjelasan singkat istilah yang paling sering
muncul.

## Konsep dasar

| Istilah | Penjelasan |
|:---|:---|
| small area estimation (SAE) | Pendugaan parameter untuk domain yang sampelnya terlalu kecil bagi estimasi langsung yang andal. |
| area / domain | Kelompok populasi yang diestimasi, misalnya kabupaten/kota. Dalam dokumentasi ini keduanya dipakai bergantian. |
| direct estimator | Estimator yang hanya memakai data sampel dari domain yang bersangkutan (mis. dari Susenas). |
| model-based estimator | Estimator yang meminjam informasi dari domain lain melalui model. |
| auxiliary variable / covariate | Variabel penyerta yang tersedia untuk seluruh domain, termasuk yang tidak tersampel. |
| sampled / unsampled area | Domain yang memiliki / tidak memiliki pengamatan survei (`y = NA`). Domain unsampled diprediksi dari model. |
| sampling variance (`vardir`) | Varians direct estimator akibat pengambilan sampel. |
| sampling error | Selisih direct estimator dengan nilai sebenarnya karena sampel. |
| relative standard error (RSE) | Simpangan baku dibagi estimasi, dalam persen. Dipakai BPS untuk menilai keandalan estimasi. |
| design-based / model-based | Inferensi berdasarkan desain sampling / berdasarkan model statistik. |

## Model

| Istilah | Penjelasan |
|:---|:---|
| area-level model | Model dengan domain sebagai unit analisis (direct estimator per domain), mis. Fay-Herriot. |
| unit-level model | Model dengan responden/elemen sampel sebagai unit analisis, mis. Battese-Harter-Fuller (BHF). |
| two-fold subarea model | Model dengan subarea yang bersarang dalam area. |
| random effect / fixed effect | Efek acak (spesifik domain) / efek tetap (koefisien regresi). |
| variance component | Komponen varians, mis. varians random effect $`\sigma_u^2`$. |
| spatial weights matrix | Matriks pembobot spasial yang menyatakan ketetanggaan antardomain. |
| SAR, CAR (Besag, BYM, BYM2, Leroux) | Struktur autokorelasi spasial untuk random effect. |
| RW1, RW2, AR(1) | Struktur temporal untuk random effect (random walk, autoregressive). |
| spatio-temporal | Gabungan efek spasial dan temporal pada data panel. |

## Estimator dan error

| Istilah | Penjelasan |
|:---|:---|
| EBLUP | *Empirical best linear unbiased predictor*: prediktor tak bias linear terbaik empiris. |
| MSE (analytical / bootstrap) | *Mean squared error* (kuadrat tengah galat) dari estimator. Dihitung dengan rumus (analytical) atau simulasi ulang (bootstrap). |
| REML / ML | *Restricted maximum likelihood* / *maximum likelihood*: metode estimasi variance component. |
| Fisher scoring | Algoritma iteratif untuk memaksimumkan likelihood. |
| shrinkage | Estimator model “menarik” direct estimator ke arah prediksi regresi (synthetic estimator); semakin besar sampling variance, semakin kuat tarikannya. |
| synthetic estimator | Prediksi regresi $`x_d^\top \hat{\beta}`$ tanpa random effect. |
| kriging | Prediksi untuk domain unsampled dari domain tetangganya. |
| benchmarking | Menyesuaikan estimasi domain agar konsisten dengan total/rata-rata pada tingkat yang lebih tinggi. |
| residual | Selisih pengamatan dan nilai terprediksi. |

## Bayesian

| Istilah | Penjelasan |
|:---|:---|
| hierarchical Bayes (HB) | Pendekatan Bayesian untuk model bertingkat. |
| prior / posterior | Distribusi parameter sebelum / sesudah melihat data. |
| credible interval | Selang kredibel: padanan Bayesian dari selang kepercayaan. |
| likelihood | Fungsi kemungkinan data. |
| INLA | *Integrated Nested Laplace Approximations*: inferensi Bayesian berbasis Laplace approximation, tanpa MCMC. |
| MCMC, chain | *Markov chain Monte Carlo*; rantai sampel yang dibangkitkan oleh sampler seperti NUTS (Stan, NumPyro). |
| engine | Mesin estimasi: [`inla_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md) atau [`jax_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md). |
| precision | Kebalikan varians. |
