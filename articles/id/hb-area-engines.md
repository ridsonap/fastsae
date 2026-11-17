# Engine hb_area(): INLA vs JAX (CPU, GPU Apple, CUDA) vs tipsae

[`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
dapat diestimasi dengan dua *engine*:

- **[`inla_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)**
  (bawaan): *integrated nested Laplace approximation* dengan
  [INLA](https://www.r-inla.org).
- **[`jax_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)**:
  MCMC penuh dengan *No-U-Turn sampler* dari NumPyro di atas
  [JAX](https://jax.readthedocs.io). *Engine* ini berjalan di CPU, di
  GPU NVIDIA (`device = "cuda"`), atau di GPU Apple silicon
  (`device = "metal"`).

Waktu INLA di halaman ini diukur dalam mode `"compact"`, mode bawaan
INLA. Sejak fastsae 1.0.0, model Gaussian diestimasi dalam mode
`"classic"` INLA, yang lebih akurat untuk model tersebut (lihat bagian
*Engines* di
[`?hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md))
dan dalam pemeriksaan kami memerlukan waktu yang kurang lebih sama;
keluarga lainnya tetap memakai mode *compact*.

Halaman ini membandingkan kecepatan kedua *engine*, di CPU maupun GPU,
serta seberapa dekat estimasi keduanya. Perbandingannya memakai model
tingkat area spasial. Untuk model Beta, perbandingan ini juga
menyertakan paket [tipsae](https://CRAN.R-project.org/package=tipsae),
yang mengestimasi model yang sama dengan MCMC di Stan.

> **Hasil GPU NVIDIA (CUDA) masih menyusul.** Hasilnya akan ditambahkan
> di sini dan dapat direproduksi dengan *notebook* yang tautannya ada di
> akhir halaman ini.

## Rancangan

**Model.**

- **Data:** kisi teratur berisi $`D`$ sel dengan ketetanggaan *queen*,
  disimulasikan dengan
  `sim_area_data(D = n, W = W, n_unsampled = 0, seed = 42 + n)`.
- **Gaussian:** `hb_area(y_gaussian ~ x1 + x2, spatial = "besag")`,
  yaitu model Fay-Herriot dengan efek area CAR intrinsik (Besag).
- **Beta:**
  `hb_area(y_beta ~ x1 + x2, family = "beta", spatial = "besag")`.
  tipsae mengestimasinya dengan
  `fit_sae(likelihood = "beta", spatial_error = TRUE)`, yang juga
  memakai efek ICAR.

tipsae hanya dijalankan untuk model Beta, sampai $`D = 2.500`$.

| Mesin | CPU | GPU | Perangkat lunak | Tanggal |
|:---|:---|:---|:---|:---|
| Apple M4 | Apple M4 (10 inti, RAM 16 GB) | Apple M4 | R 4.6.1; fastsae 1.0.0; INLA 26.8.7; tipsae 1.0.3; jax 0.5.0, jaxlib 0.5.0, jax-metal 0.1.1, numpyro 0.19.0 | 2026-10-24 |

Setiap *engine* dijalankan dalam proses R tersendiri, karena JAX
menetapkan perangkatnya pada estimasi pertama dalam sebuah sesi. Setiap
estimasi dijalankan 2 kali pada data yang sama. Putaran pertama mencakup
kompilasi JAX; tabel dan gambar menampilkan putaran yang lebih cepat,
yaitu putaran hangat setelah kompilasi (penjelajah interaktif juga dapat
menampilkan putaran pertama). JAX memakai pengaturan bawaan
[`jax_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)
(2 rantai × 1.000 tarikan setelah 500 *warm-up*), dan tipsae memakai
pengaturan bawaannya sendiri. Kodenya ada di
[`benchmarks/bench_hb_area_engines.R`](https://github.com/ridsonap/fastsae/blob/main/benchmarks/bench_hb_area_engines.R).

**Membandingkan antarmesin.** Waktu yang diukur di mesin yang berbeda
tidak dapat dibandingkan satu sama lain. Yang dapat dibandingkan adalah
rasio di dalam satu mesin: GPU vs CPU, dan INLA vs tipsae.

## Waktu komputasi

Pilih model, mesin, dan jenis waktu yang ingin ditampilkan. Arahkan
kursor ke sebuah titik untuk melihat waktunya, perangkatnya, dan R-hat
estimasinya.

Waktu komputasi

Seberapa cepat INLA dibandingkan engine lain

Gambar statis (semua mesin)

![Waktu komputasi putaran hangat (skala log). Makin rendah, makin cepat.
Satu baris untuk setiap
mesin.](hb-area-engines_files/figure-html/time-plot-1.png)

Waktu komputasi putaran hangat (skala log). Makin rendah, makin cepat.
Satu baris untuk setiap mesin.

| Mesin | Model | D | INLA | JAX CPU | JAX GPU | GPU vs CPU | tipsae | INLA vs tipsae |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|
| Apple M4 | Gaussian | 100 | 0,7 detik | 4,5 detik | 136 detik (metal) | 30,3x lebih lambat | \- | \- |
| Apple M4 | Gaussian | 500 | 3,1 detik | 32 detik | 245 detik (metal) | 7,5x lebih lambat | \- | \- |
| Apple M4 | Gaussian | 1.000 | 0,8 detik | 93 detik | 377 detik (metal) | 4,1x lebih lambat | \- | \- |
| Apple M4 | Gaussian | 2.500 | 1,4 detik | 841 detik | 725 detik (metal) | **1,2x lebih cepat** | \- | \- |
| Apple M4 | Beta | 100 | 2,6 detik | 24 detik | 221 detik (metal) | 9,2x lebih lambat | 67 detik | 26x lebih cepat |
| Apple M4 | Beta | 500 | 3,1 detik | 146 detik | 281 detik (metal) | 1,9x lebih lambat | 152 detik | 49x lebih cepat |
| Apple M4 | Beta | 1.000 | 3,7 detik | 374 detik | 376 detik (metal) | 1,0x lebih lambat | 235 detik | 63x lebih cepat |
| Apple M4 | Beta | 2.500 | 6,3 detik | 1.269 detik | 819 detik (metal) | **1,5x lebih cepat** | 620 detik | 98x lebih cepat |

Hasilnya:

- **INLA adalah *engine* tercepat pada semua ukuran.** INLA tidak
  menarik sampel. Sebagai contoh, di Apple M4 model Gaussian spasial
  dengan $`D = 2.500`$ memerlukan 1,4 detik.
- **tipsae.** Pada $`D = 2.500`$ di Apple M4, tipsae memerlukan 620
  detik untuk model Beta spasial, INLA 6,3 detik (98x lebih cepat), dan
  JAX di CPU 1.269 detik. Jadi, *engine* JAX di CPU lebih cepat daripada
  tipsae sampai $`D = 500`$ (hingga 2,8x), tetapi lebih lambat pada
  $`D`$ yang lebih besar (pada $`D = 2.500`$, waktunya 2,0x waktu
  tipsae).
- **Pada model spasial, waktu *engine* JAX melonjak seiring bertambahnya
  $`D`$.** *Engine* ini merepresentasikan efek ICAR melalui vektor eigen
  Laplacian graf berukuran $`D \times D`$. Menghitungnya memerlukan
  $`O(D^3)`$ satu kali, lalu setiap evaluasi gradien memerlukan
  $`O(D^2)`$. Komputasi matriks padat semacam inilah yang dipercepat
  oleh GPU.
- **CPU atau GPU?** Di Apple M4 (Metal): GPU lebih cepat mulai
  $`D = 2.500`$ untuk model Gaussian, sedangkan GPU lebih cepat mulai
  $`D = 2.500`$ untuk model Beta. GPU memiliki biaya tetap untuk setiap
  estimasi (pada Metal, setiap putaran dalam perulangan pembangunan
  pohon NUTS harus menunggu GPU dan CPU bersinkronisasi), sehingga GPU
  baru menguntungkan bila $`D`$ besar. `jax_engine(device = "auto")`
  memakai CUDA bila ada GPU NVIDIA dan CPU bila tidak ada, termasuk di
  macOS.

## Kesesuaian estimasi

Setiap *engine* dibandingkan dengan INLA pada mesin dan data yang sama.
`cor` adalah korelasi rata-rata posterior, `maks |selisih|` adalah
selisih absolut terbesar di antara keduanya, dan `rasio sd` adalah
median rasio simpangan baku posterior (*engine* / INLA).

| Mesin    | Model    | Engine    | D     | cor    | maks \|selisih\| | rasio sd |
|:---------|:---------|:----------|:------|:-------|:-----------------|:---------|
| Apple M4 | Gaussian | JAX CPU   | 100   | 1,0000 | 0,0088           | 1,000    |
| Apple M4 | Gaussian | JAX CPU   | 500   | 1,0000 | 0,02             | 1,000    |
| Apple M4 | Gaussian | JAX CPU   | 1.000 | 1,0000 | 0,015            | 1,000    |
| Apple M4 | Gaussian | JAX CPU   | 2.500 | 1,0000 | 0,02             | 0,999    |
| Apple M4 | Gaussian | JAX Metal | 100   | 1,0000 | 0,0098           | 0,998    |
| Apple M4 | Gaussian | JAX Metal | 500   | 1,0000 | 0,014            | 0,998    |
| Apple M4 | Gaussian | JAX Metal | 1.000 | 1,0000 | 0,019            | 1,000    |
| Apple M4 | Gaussian | JAX Metal | 2.500 | 1,0000 | 0,018            | 0,999    |
| Apple M4 | Beta     | JAX CPU   | 100   | 0,9999 | 0,014            | 1,014    |
| Apple M4 | Beta     | JAX CPU   | 500   | 0,9999 | 0,016            | 1,016    |
| Apple M4 | Beta     | JAX CPU   | 1.000 | 0,9999 | 0,014            | 1,013    |
| Apple M4 | Beta     | JAX CPU   | 2.500 | 0,9999 | 0,015            | 1,014    |
| Apple M4 | Beta     | JAX Metal | 100   | 0,9999 | 0,014            | 1,022    |
| Apple M4 | Beta     | JAX Metal | 500   | 0,9999 | 0,015            | 1,017    |
| Apple M4 | Beta     | JAX Metal | 1.000 | 0,9999 | 0,014            | 1,014    |
| Apple M4 | Beta     | JAX Metal | 2.500 | 0,9999 | 0,016            | 1,013    |
| Apple M4 | Beta     | tipsae    | 100   | 0,9865 | 0,32             | 1,068    |
| Apple M4 | Beta     | tipsae    | 500   | 0,9983 | 0,064            | 1,039    |
| Apple M4 | Beta     | tipsae    | 1.000 | 0,9982 | 0,058            | 1,038    |
| Apple M4 | Beta     | tipsae    | 2.500 | 0,9976 | 0,091            | 1,035    |

- ***Engine* JAX.** Di setiap perangkat, rata-rata posteriornya sesuai
  dengan INLA dengan korelasi paling sedikit 0,9999, dan hal yang sama
  berlaku antarperangkat. GPU mengevaluasi operasi bilangan
  *floating-point* dalam urutan yang berbeda dari CPU, sehingga
  rantainya sedikit berbeda, tetapi estimasinya tetap sama.
- **tipsae.** tipsae sesuai dengan INLA dengan korelasi paling sedikit
  0,986. Simpangan baku posteriornya berbeda (median rasio 1,03 hingga
  1,07), kemungkinan besar karena kedua paket memakai prior yang
  berbeda.
- **Konvergensi.** Di antara semua estimasi JAX, *split* R-hat terbesar
  adalah 1,016. Untuk tipsae nilainya 3,710 (pada $`D = 100`$ di Apple
  M4), di atas ambang lazim 1,05, sehingga panjang rantai bawaannya
  belum cukup pada ukuran itu.

## Memakai GPU

``` r

# Penyiapan lingkungan Python sekali saja (CUDA di Linux/Windows, Metal di Apple silicon)
fastsae::setup_numpyro_env(device = "cuda")   # atau device = "metal"

# Di sesi R yang baru: minta GPU sebelum estimasi JAX lainnya
fit <- hb_area(y ~ x1 + x2, vardir = "vardir", data = df, spatial = "besag", W = W,
               engine = jax_engine(device = "cuda"))
fit$device
```

**Metal.** jax-metal belum mengimplementasikan dua operasi yang dipakai
oleh *sampler* NUTS dari NumPyro: penghitung bit `popcnt`, dan sebagian
pernyataan `lax.cond` dengan keluaran konstan. Dengan
`device = "metal"`, fastsae menggantikan keduanya dengan operasi yang
setara dan menarik sampel rantai-rantainya secara berurutan. Bila sebuah
model tetap tidak dapat dikompilasi di Metal, estimasi beralih ke CPU
disertai peringatan.

**Mereproduksi hasil CUDA.** Dua *notebook* memasang fastsae, INLA,
tipsae, dan JAX dengan CUDA, lalu menjalankan tolok ukur yang sama:

- **Kaggle** (disarankan):
  [`benchmarks/kaggle/hb_area_engines_kaggle.ipynb`](https://github.com/ridsonap/fastsae/blob/main/benchmarks/kaggle/hb_area_engines_kaggle.ipynb).
  Impor ke *notebook* Kaggle. *Notebook* ini berbasis Python karena
  Kaggle menyediakan GPU untuk *notebook* Python, dan ia memasang R
  sendiri. Pilih akselerator GPU dan aktifkan akses internet. Setelah
  itu, gunakan *Save & Run All*, yang menjalankan seluruh tolok ukur di
  server Kaggle tanpa perlu membiarkan peramban tetap terbuka.
- **Google Colab**:
  [`benchmarks/colab/hb_area_engines_colab.ipynb`](https://colab.research.google.com/github/ridsonap/fastsae/blob/main/benchmarks/colab/hb_area_engines_colab.ipynb).
  Pilih *runtime* R dengan GPU T4 dan jalankan semua sel. Biarkan tab
  tetap terbuka sampai prosesnya selesai.

**Singkatnya:**

- Untuk sebagian besar survei (ratusan hingga beberapa ribu domain),
  [`inla_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)
  adalah pilihan tercepat.
- Untuk MCMC penuh, *engine* JAX di CPU saja sudah lebih cepat daripada
  Stan melalui tipsae untuk $`D`$ kecil hingga menengah, tetapi
  keunggulan itu hilang seiring bertambahnya $`D`$ (lihat tabel di
  atas).
- GPU baru sepadan untuk *engine* JAX pada model spasial yang besar.
