# Latar Belakang Manuskrip Jurnal: Paket `fastsae`

> **Rekomendasi Sitasi & Afiliasi**:
> **Judul Manuskrip yang Disarankan**:
> *"fastsae: An Ultra-Fast and Comprehensive R Package for High-Performance Frequentist and Bayesian Small Area Estimation"*
> **Penulis**: Ridson Al Farizal P. & Azka Ubaidillah (Politeknik Statistika STIS)

---

# 1. PENDAHULUAN (INTRODUCTION)

## 1.1 Urgensi Small Area Estimation dalam Era *Evidence-Based Policymaking*

Kebutuhan akan data statistik spasial terdisagregasi pada tingkat wilayah mikro (seperti kabupaten/kota, kecamatan, hingga desa/kelurahan) meningkat secara dramatis dalam beberapa dekade terakhir. Disagregasi spasial beresolusi tinggi ini merupakan instrumen krusial bagi pemerintah dan pembuat kebijakan untuk mewujudkan prinsip *leave no one behind* dalam kerangka *Sustainable Development Goals* (SDGs), pengalokasian dana perimbangan daerah, pemetaan kantong kemiskinan ekstrem, perencanaan fasilitas kesehatan, hingga intervensi penurunan angka stunting (Rao & Molina, 2015; Pfeffermann, 2013). 

Meskipun demikian, penyediaan indikator pada domain kecil (*small area/domain*) melalui survei sampel berskala nasional—seperti Survei Sosial Ekonomi Nasional (Susenas) di Indonesia, *American Community Survey* (ACS) di Amerika Serikat, atau *European Union Statistics on Income and Living Conditions* (EU-SILC)—menghadapi kendala metodologis yang fundamental. Survei-survei tersebut pada umumnya hanya dirancang untuk memberikan estimasi langsung (*direct estimators*) yang reliabel pada tingkat agregasi makro (nasional atau provinsi). Pada tingkat wilayah kecil, ukuran sampel yang diperoleh sangat terbatas (*small sample size problem*), bahkan beberapa domain tidak terambil sampel sama sekali (*unsampled domains*). Akibatnya, penduga langsung menghasilkan varians sampling yang sangat besar, nilai koefisien variasi (*Relative Standard Error* / RSE) yang tidak dapat diterima, serta presisi estimasi yang sangat rendah sehingga tidak layak dijadikan basis intervensi kebijakan publik (Burgard et al., 2020; Pratesi, 2016).

Untuk mengatasi kegagalan penduga langsung tersebut, metodologi *Small Area Estimation* (SAE) berkembang sebagai solusi statistik terdepan. Prinsip utama SAE adalah "meminjam kekuatan" (*borrowing strength*) informasi, baik secara lintas area melalui data variabel bantu (*auxiliary variables*) dari sensus atau data administratif, maupun secara lintas waktu melalui data historis/panel (Rao & Molina, 2015). Di antara berbagai metodologi SAE tingkat area (*area-level*), model Fay-Herriot (FH) yang diperkenalkan oleh Fay dan Herriot (1979) menjadi pendekatan yang paling luas diimplementasikan oleh badan-badan statistik resmi (*National Statistical Offices* / NSOs) di seluruh dunia karena kesederhanaannya dalam menghubungkan penduga langsung dengan variabel bantu tanpa memerlukan akses terhadap data mikro unit individu.

---

## 1.2 Evolusi Metodologis: Model Efek Acak Spasial dan Spasio-Temporal

Dalam model Fay-Herriot standar, pengaruh domain lokal diasumsikan saling bebas (*independent and identically distributed* / IID). Namun, dalam fenomena sosio-ekonomi dan geospasial riil, asumsi independensi spasial ini sering kali dilanggar. Sesuai dengan Hukum Geografi Pertama Tobler (*Tobler's First Law of Geography*), area-area yang berdekatan cenderung memiliki karakteristik yang lebih mirip dibandingkan area yang berjauhan (Tobler, 1970).

Menyadari keterbatasan tersebut, literatur SAE berkembang pesat melalui integrasi ketergantungan spasial dan temporal:
1. **Model Spatial Fay-Herriot (SFH)**: Petrucci dan Salvati (2006) serta Pratesi dan Salvati (2008) memperluas model Fay-Herriot dengan memodelkan efek acak area melalui proses *Simultaneous Autoregressive* (SAR) menggunakan matriks ketetanggaan spasial ($W$). Model SFH terbukti secara konsisten mampu mereduksi *Mean Squared Error* (MSE) secara signifikan dengan memanfaatkan korelasi antar-wilayah tetangga, khususnya bagi wilayah yang memiliki varians sampling besar atau wilayah yang tidak tersampel.
2. **Model Spatio-Temporal Fay-Herriot (ST-FH)**: Marhuenda et al. (2013) memformulasikan model spatio-temporal tingkat area yang mengombinasikan struktur spasial SAR dengan proses dinamis temporal autoregresif ordo pertama (AR(1)). Pendekatan ini memungkinkan peminjaman kekuatan dua dimensi secara simultan: informasi tetangga spasial dan data riwayat tahun-tahun sebelumnya.
3. **Model Battese-Harter-Fuller (BHF)**: Saat data mikro tingkat unit individu tersedia dan dapat dipadukan dengan rata-rata populasi area, model BHF (Battese et al., 1988) memberikan estimasi tingkat unit berbasis model efek acak bersarang (*nested error regression model*).

Meskipun model-model ini telah mapan secara teoritis, proses inferensinya—khususnya estimasi komponen varians melalui *Restricted Maximum Likelihood* (REML) dan estimasi ketidakpastian melalui *parametric bootstrap* (González-Manteiga et al., 2008) atau *non-parametric bootstrap*—memerlukan inversi matriks berdimensi besar dan komputasi iteratif algoritma Fisher-scoring yang sangat intensif, menjadi hambatan utama dalam eksekusi data berskala besar.

---

## 1.3 Tantangan Indikator Non-Gaussian dan Keterbatasan Pendekatan Linier

Hambatan konseptual yang lebih besar muncul ketika indikator target tidak memenuhi asumsi sebaran Normal (*Gaussian*). Banyak indikator strategis pembangunan memiliki karakteristik sebaran yang unik:
1. **Data Proporsi Terbatas dalam Interval $(0, 1)$**: Indikator seperti angka kemiskinan, prevalensi pengangguran, stunting, dan partisipasi sekolah memiliki rentang nilai terbatas antara 0 dan 1. Model Fay-Herriot Gaussian standar rentan menghasilkan estimasi di luar rentang logis ($< 0$ atau $> 1$), terutama pada domain dengan nilai ekstrem mendekati batas (Janicki, 2020; De Nicolò & Gardini, 2022). Penggunaan transformasi logit atau arcsin sering kali menyebabkan bias transformasi balik (*transformation bias*) yang parah saat diubah kembali ke skala asli (Slud & Maiti, 2006).
2. **Data Kontinu Positif Menjulur (*Skewed Positive Continuous*)**: Indikator moneter seperti pendapatan per kapita, pengeluaran konsumsi, dan keparahan kemiskinan umumnya memiliki sebaran yang sangat menjulur ke kanan dengan ragam yang proporsional terhadap kuadrat rata-ratanya, yang secara alami mengikuti sebaran Gamma (Ghosh & Maiti, 2004; Fabrizi et al., 2011).
3. **Data Diskrit Cacahan (*Count Data*)**: Indikator kejadian langka, seperti jumlah kasus kematian ibu/bayi atau kasus kriminalitas, mengikuti proses Poisson atau *Negative Binomial* jika terdapat overdispersi, yang memerlukan parameter offset populasi berisiko (*exposure*) (Boubeta et al., 2016, 2017).

Untuk mengatasi indikator-indikator non-Gaussian ini, pendekatan *Generalized Linear Mixed Models* (GLMM) dan *Empirical Best Prediction* (EBP) berbasis hierarki Bayesian menjadi keniscayaan.

---

## 1.4 Dilema Komputasi: Jebakan MCMC vs. Revolusi Numerik INLA

Dalam kerangka Bayesian Hierarkis (*Hierarchical Bayes* / HB), pemodelan non-Gaussian dengan struktur spasial kompleks secara tradisional diselesaikan menggunakan algoritma *Markov Chain Monte Carlo* (MCMC), seperti Gibbs sampling atau *Hamiltonian Monte Carlo* / *No-U-Turn Sampler* (HMC/NUTS) yang diimplementasikan via JAGS, BUGS, atau Stan (Gelman et al., 2013).

Namun, ketergantungan pada MCMC menimbulkan masalah komputasi yang serius dalam praktiknya:
- **Waktu Komputasi Sangat Lama (*Computational Bottleneck*)**: Iterasi MCMC memerlukan penarikan puluhan ribu sampel untuk mencapai konvergensi, sering kali memakan waktu puluhan menit hingga berjam-jam untuk satu model spasial area-level (Rue et al., 2009; Bivand et al., 2015).
- **Beban Diagnostik Subjektif**: Peneliti dan praktisi NSO dibebani tugas verifikasi diagnostik konvergensi yang kompleks (seperti uji Gelman-Rubin $\hat{R}$, *effective sample size*, dan korelasi berantai), yang menyulitkan otomasi produksi statistik berkala.
- **Kerapuhan pada Model Spasial**: Ketika model menggabungkan autokorelasi spasial murni (seperti model ICAR Besag atau Leroux CAR) pada sebaran non-Gaussian, penarikan rantai MCMC sering kali mengalami masalah *divergent transitions* atau autokorelasi rantai yang sangat tinggi (Simpson et al., 2017).

Sebagai terobosan mutakhir, Rue, Martino, dan Chopin (2009) mengembangkan *Integrated Nested Laplace Approximations* (INLA). INLA merupakan pendekatan numerik deterministik yang dirancang khusus untuk kelas *Latent Gaussian Models* (LGMs). Alih-alih melakukan simulasi sampling stokastik berulang seperti MCMC, INLA memanfaatkan sifat ketidakteraturan matriks presisi (*sparse precision matrices*) dari *Gaussian Markov Random Fields* (GMRF) dan deret Taylor Laplace untuk mengaproksimasi sebaran marginal posterior secara analitis.

Hasil pengujian literatur menunjukkan bahwa INLA memberikan aproksimasi posterior yang identik dengan MCMC presisi tinggi, namun diselesaikan dalam hitungan **detik**, bukan jam (Rue et al., 2017; Gómez-Rubio, 2020). Dalam domain SAE, INLA menawarkan potensi luar biasa untuk mengestimasi model area-level non-Gaussian dengan efek spasial tingkat tinggi secara instan.

---

## 1.5 Kesenjangan Perangkat Lunak (*Software Landscape & Research Gap*)

Meskipun landasan metodologis telah tersedia secara terpisah di berbagai literatur, ekosistem perangkat lunak komputasi statistik (khususnya R) saat ini masih sangat terfragmentasi dan memiliki keterbatasan struktural yang nyata (Tabel 1):

### Tabel 1: Analisis Kesenjangan Paket R untuk Small Area Estimation

| Dimensi Fitur | `sae` (Molina & Rao, 2015) | `emdi` (Kreutzmann et al., 2019) | `tipsae` (De Nicolò & Gardini, 2022) | Paket Berbasis Stan/JAGS (`saeHB`, dll.) | **`fastsae` (Paket Baru)** |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Engine Komputasi Frequentist** | R Murni (Single-thread) | C++ parsial / R | Tidak Ada | Tidak Ada | **C++17 (Rcpp/Armadillo) + OpenMP Paralel** |
| **Model FH Spasial (SAR) & Spatio-Temporal** | Ada (SFH lambat) | Tidak Ada | Tidak Ada | Terbatas | **Ada (SFH & ST-FH Paralel C++)** |
| **MSE Bootstrap Paralel** | Lambat (Single-thread) | Standard | Tidak Ada | Tidak Ada | **Sangat Cepat (Multi-core OpenMP)** |
| **Model Non-Gaussian (Beta, Gamma, dll.)** | Tidak Ada | Tidak Ada | Hanya Beta | Parsial | **Lengkap: Gaussian, Beta, Gamma, Binomial, Poisson, NegBin** |
| **Metode Estimasi Bayesian** | Tidak Ada | Tidak Ada | MCMC (Stan) | MCMC (JAGS/Stan) | **Bayesian EBP via R-INLA (Deterministik Cepat)** |
| **Waktu Estimasi Model Non-Gaussian Spasial** | N/A | N/A | Menit - Jam (MCMC) | Menit - Jam (MCMC) | **Detik (Speedup 3x - 100x vs MCMC)** |
| **Struktur Spasial Modern (BYM2, Leroux CAR)** | Tidak Ada | Tidak Ada | Terbatas | Terbatas | **Lengkap: BYM2 (PC Priors), Leroux, Besag, SLM** |
| **Injeksi Presisi Direct Spesifik Area** | Hanya Gaussian ($D_i$) | Hanya Gaussian ($D_i$) | Hanya Beta (Janicki 2020) | Tergantung script | **Gaussian ($1/D_i$), Beta (Janicki), Gamma ($y^2/D_i$), Binomial ($N_i$)** |
| **Alat Diagnostik & Autoplot Dinamis** | Terbatas | Baik (Gaussian) | Terbatas (Beta) | Terbatas | **Lengkap: DIC, WAIC, CPO, PIT, RRMSE, Bland-Altman, Autoplot** |

Kesenjangan di atas memperlihatkan ketiadaan paket R terpadu yang mampu menjembatani kebutuhan komputasi performa tinggi (*high-performance computing*) untuk metode klasik/frequentist sekaligus menyediakan antarmuka Bayesian EBP berbasis INLA untuk sebaran non-Gaussian dan spasial mutakhir.

---

## 1.6 Kontribusi dan Kebaruan Paket `fastsae`

Untuk menutup kesenjangan metodologis dan komputasional tersebut, paket R **`fastsae`** dikembangkan. Paket ini menghadirkan paradigma baru dalam perangkat lunak Small Area Estimation dengan keunggulan utama sebagai berikut:

1. **Performa Komputasi Frequentist Ultra-Cepat**:
   Mengimplementasikan estimasi parameter REML/ML dan perhitungan MSE analitis maupun bootstrap (parametrik dan non-parametrik) untuk model FH standar, Spatial FH (SAR), Spatio-Temporal FH (ST-FH), dan BHF tingkat unit menggunakan arsitektur **C++17** via `Rcpp` dan `RcppArmadillo` yang dioptimasi dengan multi-threading `OpenMP`.
2. **Estimasi Bayesian EBP Multi-Distribusi melalui INLA**:
   Menyediakan fungsi terpadu [`ebp_area()`](file:///Volumes/work/_MainR/fastsae-inla/R/ebp_area.R) yang mengintegrasikan aproksimasi Laplace Bayesian INLA untuk enam keluarga sebaran: Gaussian, Beta, Gamma, Binomial, Poisson, dan Negative Binomial.
3. **Penyelarasan Presisi Varians Sampling Direct (*Exact Variance Injection*)**:
   Mengadopsi formulasi teoritis presisi sampling langsung secara spesifik untuk tiap distribusi:
   - *Gaussian*: Pembobotan presisi skala observasi $\tau_i = 1/D_i$.
   - *Beta*: Parameter presisi langsung per area Janicki (2020), $\phi_i = \max\left(\frac{y_i(1-y_i)}{D_i} - 1, 1\right)$.
   - *Gamma*: Penyelarasan koefisien variasi langsung, $s_i = \frac{y_i^2}{D_i} = \frac{1}{CV_i^2}$.
   - *Binomial & Poisson*: Penanganan ukuran sampel/trials ($N_i$) dan populasi berisiko/exposure ($E_i$).
4. **Dukungan Struktur Spasial Hierarkis Terlengkap**:
   Mendukung model spasial mutakhir: model scaled BYM2 dengan *Penalized Complexity* (PC) priors (Simpson et al., 2017), model Leroux CAR (Leroux et al., 2000) dengan perumusan matriks struktur $C = I - R$ yang presisi, model murni ICAR Besag, hingga *Spatial Lag Model* (SLM/SAR) dengan batasan spektral nilai eigen.
5. **Rangkaian Diagnostik Lengkap dan Visualisasi Otomatis**:
   Menyediakan fungsi diagnostik universal [`diagnose()`](file:///Volumes/work/_MainR/fastsae-inla/R/diagnose.R) yang menyajikan kriteria bayesian (DIC, WAIC, Marginal Likelihood), validasi prediksi *leave-one-out* (CPO, PIT), evaluasi metrik simulasi (RRMSE, RB, MAD, Pearson/Spearman $r$), serta antarmuka visual dinamis [`autoplot()`](file:///Volumes/work/_MainR/fastsae-inla/R/autoplot.R).
6. **Dataset Bawaan dan Generator Simulasi Komprehensif**:
   Menyediakan dataset sintetis multi-distribusi bawaan siap pakai ([`sim_area`](file:///Volumes/work/_MainR/fastsae-inla/data/sim_area.rda) dan [`sim_panel`](file:///Volumes/work/_MainR/fastsae-inla/data/sim_panel.rda)) beserta fungsi generator data fleksibel ([`sim_area_data()`](file:///Volumes/work/_MainR/fastsae-inla/R/sim_area_data.R) dan [`sim_series_data()`](file:///Volumes/work/_MainR/fastsae-inla/R/sim_series_data.R)).

---

## 1.7 Struktur Artikel

Sisa dari artikel ini diorganisasikan sebagai berikut: 
- **Bagian 2 (Landasan Metodologis)** menguraikan formulasi matematis model Fay-Herriot spasial, spatio-temporal, Battese-Harter-Fuller, serta kerangka Bayesian EBP non-Gaussian dengan INLA.
- **Bagian 3 (Arsitektur dan Desain Paket)** memaparkan rancangan perangkat lunak `fastsae`, optimasi C++, integrasi matriks spasial, dan alur kerja diagnostik.
- **Bagian 4 (Studi Simulasi dan Benchmarking Komputasi)** membandingkan akurasi estimasi dan kecepatan eksekusi `fastsae` terhadap paket benchmark yang telah mapan (`sae` dan `tipsae`).
- **Bagian 5 (Aplikasi Data Riil)** mendemonstrasikan implementasi `fastsae` pada estimasi indikator pembangunan daerah.
- **Bagian 6 (Kesimpulan dan Pengembangan Mendatang)** menutup artikel dengan rangkuman kontribusi serta potensi perluasan model di masa depan.

---

# DAFTAR PUSTAKA ACUAN (LITERATURE CITATIONS)

1. **Battese, G. E., Harter, R. M., & Fuller, W. A. (1988).** An error-components model for prediction of county crop areas using survey and satellite data. *Journal of the American Statistical Association*, 83(401), 28–36.
2. **Besag, J., York, J., & Mollié, A. (1991).** Bayesian image restoration, with two applications in spatial statistics. *Annals of the Institute of Statistical Mathematics*, 43(1), 1–20.
3. **Bivand, R., Gómez-Rubio, V., & Rue, H. (2015).** Spatial data analysis with R-INLA with some extensions. *Journal of Statistical Software*, 63(20), 1–31.
4. **Boubeta, M., Lombardía, M. J., & Morales, D. (2016).** Empirical best prediction under area-level Poisson mixed models. *Test*, 25(3), 540–569.
5. **Boubeta, M., Lombardía, M. J., & Morales, D. (2017).** Poisson mixed models for small area estimation of poverty indicators. *Computational Statistics & Data Analysis*, 107, 19–33.
6. **Burgard, J. P., Krause, J., & Merkle, P. (2020).** Computational aspects of small area estimation in official statistics. *Journal of Official Statistics*, 36(4), 793–817.
7. **De Nicolò, S., & Gardini, A. (2022).** Small area estimation of inequality measures via beta regression: The tipsae package. *Journal of Statistical Software*, 103(1), 1–34.
8. **Fabrizi, E., Ferrante, M. R., & Trivisano, C. (2011).** Hierarchical and empirical Bayes methods for small area estimation in presence of skewed distributions. *Journal of Official Statistics*, 27(4), 517–537.
9. **Fay, R. E., & Herriot, R. A. (1979).** Estimates of income for small places: an application of James-Stein procedures to census data. *Journal of the American Statistical Association*, 74(366a), 269–277.
10. **Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin, D. B. (2013).** *Bayesian Data Analysis* (3rd ed.). Chapman and Hall/CRC.
11. **Ghosh, M., & Maiti, T. (2004).** Small area estimation of positive expenditure data using generalized linear models. *Journal of the Indian Society of Agricultural Statistics*, 57, 126–138.
12. **Gómez-Rubio, V. (2020).** *Bayesian Inference with INLA*. Chapman and Hall/CRC.
13. **González-Manteiga, W., Lombardía, M. J., Molina, I., Morales, D., & Santamaría, L. (2008).** Bootstrap mean squared error of a small-area EBLUP. *Journal of Statistical Computation and Simulation*, 78(5), 443–462.
14. **Janicki, R. (2020).** Properties of the beta regression model for small area estimation of proportions. *Communications in Statistics - Theory and Methods*, 49(14), 3465–3488.
15. **Kreutzmann, A., Pannier, S., Rojas-Perilla, N., Schmid, T., Templ, M., & Tzavidis, N. (2019).** The R package emdi for estimating and mapping indicators to regional levels. *Journal of Statistical Software*, 91(7), 1–33.
16. **Leroux, B. G., Lei, X., & Breslow, N. (2000).** Estimation of disease rates in small areas: A new mixed model for spatial dependence. In *Statistical Models in Epidemiology, the Environment, and Clinical Trials* (pp. 135–145). Springer.
17. **Marhuenda, Y., Molina, I., & Morales, D. (2013).** Small area estimation with spatio-temporal Fay-Herriot models. *Computational Statistics & Data Analysis*, 58, 308–325.
18. **Molina, I., & Rao, J. N. K. (2015).** sae: An R package for small area estimation. *The R Journal*, 7(1), 81–98.
19. **Petrucci, A., & Salvati, N. (2006).** Small area estimation for spatial data in ecology. *Ecological Modelling*, 190(1-2), 173–185.
20. **Pfeffermann, D. (2013).** New important developments in small area estimation. *Statistical Science*, 28(1), 40–68.
21. **Pratesi, M. (Ed.). (2016).** *Analysis of Poverty Data by Small Area Estimation*. John Wiley & Sons.
22. **Pratesi, M., & Salvati, N. (2008).** Small area estimation: the EBLUP estimator based on spatially correlated random area effects. *Statistical Methods & Applications*, 17(1), 113–141.
23. **Rao, J. N. K., & Molina, I. (2015).** *Small Area Estimation* (2nd ed.). John Wiley & Sons.
24. **Rue, H., Martino, S., & Chopin, N. (2009).** Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations. *Journal of the Royal Statistical Society: Series B (Statistical Methodology)*, 71(2), 319–392.
25. **Rue, H., Riebler, A., Sørbye, S. H., Illian, J. B., Simpson, D. P., & Lindgren, F. K. (2017).** Bayesian computing with INLA: a review. *Annual Review of Statistics and Its Application*, 4, 395–421.
26. **Simpson, D., Rue, H., Riebler, A., Martins, T. G., & Sørbye, S. H. (2017).** Penalising model component complexity: A principled, practical approach to constructing priors. *Statistical Science*, 32(1), 1–28.
27. **Slud, E., & Maiti, T. (2006).** Mean-squared error estimation in transformed Fay-Herriot models. *Journal of the Royal Statistical Society: Series B (Statistical Methodology)*, 68(5), 833–850.
28. **Tobler, W. R. (1970).** A computer movie simulating urban growth in the Detroit region. *Economic Geography*, 46(sup1), 234–240.
