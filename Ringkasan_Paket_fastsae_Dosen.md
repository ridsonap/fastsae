# RINGKASAN EKSEKUTIF PENGEMBANGAN PAKET R: `fastsae`
**Fast & Scalable Small Area Estimation via Integrated Nested Laplace Approximations**

---

**Diajukan Oleh:** Mahasiswa Bimbingan / Tim Peneliti  
**Tujuan Dokumen:** Laporan kemajuan penelitian dan pengajuan ringkasan fungsionalitas paket R `fastsae` untuk persiapan publikasi ilmiah (*Journal of Statistical Software* / *The R Journal*) serta implementasi produksi data statistik resmi (*official statistics*).  
**Repositori Kode:** `https://github.com/ridsonap/fastsae`  
**Status Kesiapan:** 529 unit tests lulus (100%), dokumentasi manual lengkap (CRAN-compliant), naskah publikasi 17 halaman siap telaah.

---

## 1. Latar Belakang & Urgensi Masalah

Dalam penyusunan kebijakan publik berbasis bukti (*evidence-based policy*), pemerintah membutuhkan indikator sosial-ekonomi (seperti persentase kemiskinan, prevalensi stunting, tingkat pengangguran) hingga tingkat wilayah administrasi yang sangat kecil (kabupaten/kota, kecamatan, bahkan desa/kelurahan).

Namun, Badan Pusat Statistik (BPS) dan lembaga statistik resmi internasional menghadapi **tiga kendala fundamental**:
1. **Keterbatasan Ukuran Sampel Survei:** Penduga langsung (*direct estimator*) seperti survei Susenas/Sakernas memiliki ukuran sampel yang sangat kecil atau bahkan nol di level mikro, menghasilkan *Relative Standard Error* (RSE / CV) yang sangat tinggi (>30%) sehingga tidak reliabel untuk dipublikasikan.
2. **Hambatan Komputasi Metode Bayesian (MCMC):** Pendekatan Bayesian hierarkis selama ini mengandalkan algoritma *Markov Chain Monte Carlo* (MCMC melalui Stan/JAGS). Pada skala nasional dengan ratusan kabupaten/kota atau ribuan desa, dan terutama pada pemodelan spasio-temporal, MCMC membutuhkan waktu komputasi yang sangat lambat (berkisar antara puluhan menit hingga berjam-jam) serta sering menghadapi masalah konvergensi rantai (*chain convergence*).
3. **Ketiadaan Solusi Terpadu untuk Model Non-Gaussian & Kalibrasi Konsistensi:** 
   - Paket R SAE konvensional (misalnya paket `sae` dan `emdi`) sebagian besar mengasumsikan distribusi normal/Gaussian, yang sering menghasilkan estimasi tidak realistis (misalnya proporsi bernilai negatif atau >100%).
   - Angka estimasi area kecil sering kali tidak cocok jika dijumlahkan kembali ke angka agregat resmi nasional atau provinsi (*inconsistency issue*).

---

## 2. Solusi & Kebaruan yang Ditawarkan Paket `fastsae`

Paket **`fastsae`** dirancang sebagai ekosistem komputasi SAE modern yang menggabungkan kecepatan komputasi deterministik tinggi, fleksibilitas distribusi probabilitas, dan kepatuhan terhadap kaidah produksi statistik resmi:

1. **Akselerasi Berbasis INLA (*Integrated Nested Laplace Approximations*):**
   Mengganti simulasi MCMC dengan aproksimasi analitik Laplace bersarang dan aljabar matriks jarang (*sparse matrix* C++). Memberikan hasil posterior mean dan interval kredibel yang identik dengan MCMC, namun dengan kecepatan **10x hingga 100x lebih cepat** (konvergensi dalam hitungan detik).
2. **Cakupan Model Area-Level & Unit-Level Terluas:**
   - **Area-Level Gaussian:** Model standar Fay-Herriot (FH), Spasial FH (SAR, CAR, BYM2), dan Spasio-Temporal FH (AR(1) temporal autoregressive dengan struktur interaksi ruang-waktu Knorr-Held).
   - **Non-Gaussian Hierarchical Bayes (HB):** Model Beta (untuk indikator terbatas $0 < \theta < 1$ seperti proporsi kemiskinan), Binomial/Logit (prevalensi dengan denominator ukuran sampel), dan Poisson (frekuensi kejadian langka dengan offset populasi berisiko).
   - **Unit-Level Models:** Model Battese-Harter-Fuller (BHF) dengan metode analitik cepat Henderson III dan REML.
3. **Inovasi Kalibrasi: *Two-Stage Hierarchical Benchmarking* (`benchmark_sae`):**
   Memastikan estimasi level kabupaten/kota secara serentak konsisten dengan target provinsi, dan target provinsi konsisten dengan target nasional:
   $$\sum_{d \in \text{Prov}_g} \tilde{w}_{dg} \hat{\theta}_{dg}^{\text{BM}} = T_g^* \quad \text{dan} \quad \sum_g \tilde{W}_g T_g^* = T_{\text{nat}}$$
   Mendukung metode *Ratio*, *Difference*, *Optimal Variance-Weighted Quadratic Loss* (Datta et al. 2011; Steorts et al. 2014), dan *Constrained Logit* (Berg & Fuller 2014).
4. **Perangkat Diagnostik & Standar Mutu Statistik Otomatis:**
   Dilengkapi uji kalibrasi Brown et al. (2001) (uji regresi bias dan uji Wald $W$), otomatisasi klasifikasi reliabilitas indikator (Reliable / Caution / Unreliable), komparasi model otomatis (`compare_models`), visualisasi peta tematik (`map_sae`), serta ekspor data multi-format ber-metadata (`export_sae`).

---

## 3. Matriks Perbandingan `fastsae` terhadap Paket R Sejenis

Tabel berikut merangkum posisi keunggulan `fastsae` dibandingkan dengan paket-paket SAE mapan di CRAN:

| Fitur / Dimensi | `fastsae` | `sae` | `emdi` | `hbsae` | `tipsae` | `SUMMER` |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Engine Inferensi** | **INLA + C++** | REML / FH | REML / MCMC | MCMC | Stan MCMC | INLA |
| **Model Fay-Herriot Standar** | **Ya** | Ya | Ya | Ya | Tidak | Tidak |
| **Model Spasial (SAR/CAR/BYM2)** | **Ya** | Ya | Tidak | Tidak | Ya (CAR) | Ya (BYM2) |
| **Spatio-Temporal (AR1, Knorr-Held)** | **Ya** | Ya | Tidak | Tidak | Ya | Ya |
| **Distribusi Non-Gaussian (Beta/Bin/Pois)** | **Ya (Lengkap)**| Tidak | Tidak | Tidak | Hanya Beta | Bin/Pois |
| **Unit-Level Model (BHF)** | **Ya** | Ya | Ya | Ya | Tidak | Tidak |
| **Two-Stage Hierarchical Benchmarking** | **Ya** | Tidak | Hanya 1-Tahap | Tidak | Tidak | Tidak |
| **Uji Diagnostik Brown et al. (2001)** | **Ya (Built-in)**| Manual | Manual | Tidak | Tidak | Tidak |
| **Kecepatan pada Data Skala Besar** | **Sangat Cepat (<5s)** | Cepat | Cepat-Sedang | Lambat | Sangat Lambat | Cepat |
| **Pemetaan & Ekspor Produksi Terpadu** | **Ya (`map_sae`)**| Tidak | Ya | Tidak | Tidak | Parsial |

---

## 4. Alur Kerja Praktis Penggunaan Paket (*Workflow*)

Implementasi paket dirancang sangat intuitif menggunakan sintaks formula R standar:

```r
library(fastsae)

# 1. Fitting Model Spasio-Temporal Bayesian Area-Level via INLA
fit_model <- hb_area(
  formula = y_rate ~ poverty_proxy + illiteracy,
  family  = "beta",          # Menjamin estimasi tetap di rentang (0, 1)
  spatial = "bym2",          # Struktur korelasi spasial tetangga
  W       = adjacency_matrix,
  data    = data_kabupaten
)

# 2. Diagnostik Ilmiah & Evaluasi Kualitas Model
diag <- diagnose(fit_model)
summary(diag)                # Menampilkan reduksi RSE & uji Brown et al. (2001)

# 3. Two-Stage Hierarchical Benchmarking (Kabupaten -> Provinsi -> Nasional)
fit_calibrated <- benchmark_sae(
  fit_model,
  group           = "id_provinsi",
  weight          = "populasi",
  national_target = 9.36,     # Target resmi nasional (misal angka kemiskinan nasional)
  method          = "ratio"
)

# 4. Visualisasi Spasial & Ekspor Hasil
map_sae(fit_calibrated, sf_geom = peta_kabupaten, variable = "benchmarked")
export_sae(fit_calibrated, file = "Estimasi_Resmi_Kabupaten.xlsx")
```

---

## 5. Ringkasan Bukti Empiris & Validasi Numerik

1. **Efisiensi Komputasi (Waktu Eksekusi):**
   - Pada simulasi model Spatio-Temporal Beta SAE ($D = 100$ area, $T = 3$ periode waktu):
     - `tipsae` (Stan MCMC): **~420 detik (7 menit)** dengan 4 chains x 2000 iterasi.
     - `fastsae` (INLA): **~3.4 detik** (akselerasi lebih dari **120x lebih cepat**).
2. **Validasi Populasi Hingga (*Finite Population Monte Carlo Validation*):**
   - Berdasarkan eksperimen Monte Carlo skala nasional ($D = 500$ domain, populasi $N \approx 250.000$ individu, $R = 25$ replikasi):
     - **Relative Bias (RB):** Model `fastsae` menghasilkan bias mendekati 0% (rata-rata $|RB| < 0.8\%$), membuktikan sifat tak bias secara empiris.
     - **Efisiensi Estimasi:** Menghasilkan penurunan Relative Root Mean Squared Error (RRMSE) hingga **35% - 50%** dibandingkan penduga langsung (*Direct Estimator*), terutama pada area-area dengan ukuran sampel kecil.

---

## 6. Poin-Poin yang Dimohonkan Arahan dari Dosen Pembimbing

Untuk penyempurnaan naskah artikel ilmiah dan paket ini, beberapa poin diskusi yang dimohonkan masukan:
1. **Pemilihan Kasus Data Riil:** Apakah sebaiknya artikel memfokuskan studi kasus pada data riil kemiskinan kabupaten/kota di Indonesia (data BPS: Susenas & Podes), stunting balita (SKI/SSGI), atau data internasional terstandar?
2. **Spesifikasi Prior Bayesian:** Masukan terhadap pemilihan prior *default* (misalnya penalized complexity priors / PC-priors pada parameter dispersi dan korelasi spasial) agar sesuai dengan konsensus metodologis di literatur statistika terapan.
3. **Penyelarasan Format Publikasi:** Kesiapan draf naskah (*manuscript*) yang saat ini disusun menggunakan format *Journal of Statistical Software (JSS)* atau apakah ada target jurnal alternatif yang lebih diutamakan oleh Bapak/Ibu Dosen.

---

**Lampiran yang Tersedia:**
1. Naskah lengkap draf artikel jurnal: `paper/article.pdf` (17 halaman).
2. Dokumentasi fungsi paket dan vignettes: Direktori `man/` dan `vignettes/`.
3. Skrip reproduktibilitas benchmark simulasi: Direktori `benchmarks/`.
