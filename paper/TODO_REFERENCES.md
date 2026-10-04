# TODO_REFERENCES — entri yang TIDAK dapat diverifikasi

Aturan: tidak ada DOI/entri yang dikarang. Apa pun yang tidak dapat dikonfirmasi ke
sumber penerbit/Crossref/CRAN dicatat di sini dan **tidak dikutip** di laporan.

Semua entri yang **berhasil diverifikasi** ada di `references.bib` (34 entri, DOI
dikonfirmasi lewat `api.crossref.org/works/<doi>` dan halaman penerbit).

## 1. Tidak diverifikasi → tidak dikutip

| Entri yang diminta | Status | Alasan |
|---|---|---|
| Di Ranzo & Marhuenda — makalah paket `tipsae` | **TIDAK DITEMUKAN** | Tidak ada penulis "Di Ranzo" pada item SAE mana pun di Crossref; kolaborasi Di Ranzo–Marhuenda tidak muncul di pencarian web. Sitasi kanonik `tipsae` ternyata **De Nicolò & Gardini (2024)**, JSS 108(1):1–36, DOI `10.18637/jss.v108.i01` (sudah dimasukkan). |
| "Baraldi & Petrucci" untuk beta Fay–Herriot | **TIDAK DITEMUKAN** | Crossref/web hanya mengembalikan karya Petrucci tentang FH spasial bersama Pratesi/Salvati. Tidak dikutip. |
| Makalah EBP untuk *Human Development Index* | **TIDAK DIVERIFIKASI** | Pekerjaan SAE untuk HDI ada, tetapi tidak ditemukan makalah EBP/beta yang dapat dikonfirmasi. Tidak dikutip. |
| Tzavidis dkk. — makalah beta-FH | **TIDAK DIVERIFIKASI** | Tidak ditemukan makalah beta-FH Tzavidis yang dapat dikonfirmasi. Sebagai gantinya dipakai entri terverifikasi: Janicki (2020), De Nicolò dkk. (2024), Esteban dkk. (2019). |
| Datta & Lahiri (2000), *Statistica Sinica* 10:613–627 | **Terverifikasi tanpa DOI** | Makalah dikonfirmasi di situs resmi *Statistica Sinica* (10(2), 613–627), tetapi tidak ada DOI yang ditetapkan (era sebelum DOI jurnal tersebut) → kolom `doi` dihilangkan, bukan ditebak. |
| Wang, Fuller & Qu (2008), *Survey Methodology* 34:29–36 | **Terverifikasi tanpa DOI** | Dikonfirmasi di situs Statistics Canada (*Survey Methodology* 34(1):29–36, Cat. 12-001-X); jurnal ini tidak terdaftar di Crossref → tanpa `doi`. |
| Item non-arsip: working paper AMSActa, buletin *The Survey Statistician* (2024) tentang `tipsae` | **TIDAK DIVERIFIKASI** | Metadata tidak konsisten antar sumber (nomor/bulan buletin bertentangan). Tidak dikutip. |

## 2. Koreksi terhadap daftar awal (bukan ketidakverifikasian)

- **emdi (Kreutzmann dkk. 2019)** — judul pada daftar awal ("Estimation of Indicators for
  Development Domains and Spatial Transformation of Data") **salah**. Judul sebenarnya:
  *"The R Package emdi for Estimating and Mapping Regionally Disaggregated Indicators"*,
  JSS 91(7):1–33, DOI `10.18637/jss.v091.i07`. Daftar penulis (6 orang, sesuai urutan
  terbit): Kreutzmann, Pannier, Rojas-Perilla, Schmid, Templ, Tzavidis — **bukan** daftar
  panjang (Germann, Rao, Marhuenda, Brumback, Nandrau, Salvati, Pratesi, …) yang muncul di
  daftar awal; daftar itu adalah penulis *paket* CRAN, bukan makalah JSS.
- **Riebler, Sørbye, Rue (2016)** — makalah punya **empat** penulis: Riebler, Sørbye,
  **Simpson**, Rue. Disesuaikan pada `references.bib`.
- **Besag (1974)** — pagination: badan artikel di metadata Crossref/OUP = 192–225, tetapi
  sitasi kanonik dalam literatur (termasuk diskusi) = **192–236**; dipertahankan 192–236
  dan dicatat di sini.
- **Banerjee, Carlin & Gelfand (2014)** — Crossref mencatat buku terbit 2014-09-12; beberapa
  katalog mencetak 2015. Dipakai 2014 (edisi ke-2).
- **`tipsae`** — deskripsi: model hierarki Bayes berbasis Stan/HMC untuk respons di interval
  unit (Beta klasik, Flexible Beta, zero/one-inflated Beta), bukan "empirical best predictor
  untuk indikator beta". Rumusan lama dibuang.

## 3. Referensi yang disebut kode tetapi tanpa DOI di kode

Blok `@references` di `R/*.R` tidak memuat satu pun DOI/URL (hanya URL repo INLA di
`R/inla_utils.R:12`). Semua DOI di `references.bib` berasal dari verifikasi eksternal,
bukan dari kode. Sumber roxygen per fungsi dicatat di `paper/_notes/notes_*.md`.

## 4. Entri terverifikasi dengan metadata yang tidak lengkap

| Entri | DOI terverifikasi | Yang tidak tersedia |
|---|---|---|
| Bates, Mächler, Bolker & Walker (2015), *Fitting Linear Mixed-Effects Models Using lme4* | `10.18637/jss.v067.i01` (Crossref + content negotiation DOI) | Rentang halaman tidak didaftarkan penerbit (JSS tidak mengisi kolom `page`) → kolom `pages` dihilangkan, bukan ditebak. |
| Lindgren & Rue (2015), *Bayesian Spatial Modelling with R-INLA* | `10.18637/jss.v063.i19` (sama) | Sama: tanpa rentang halaman di metadata resmi. |

## 5. Rujukan di dalam blok @references kode (belum terverifikasi)

Blok `@references` `R/benchmark.R:79-93` dan `R/diagnose.R:27-33` memuat lima
dan tiga entri; yang statusnya:

| Entri (sesuai kode) | Status | Tindakan |
|---|---|---|
| Rao & Molina (2015), *Small Area Estimation*, Bab 10 "Benchmarking and Other Practical Issues", pp. 297-315 | Terverifikasi sebagai `rao2015sae` (Wiley, DOI terdaftar); rentang halaman bab dikutip dari kode | Dikutip memakai `rao2015sae`; halaman bab dicatat apa adanya |
| You & Rao (2002), *Canad. J. Statist.* 30(3), 431-439 | Terverifikasi sebagai `you2002pseudo` (tanpa DOI terdaftar) | Dikutip |
| Datta, Ghosh, Steorts & Maples, *TEST* 20(3), 574-588, **2011** (kode) | **Terverifikasi**: DOI `10.1007/s11749-010-0218-y`, judul/vol/2/halaman cocok; tetapi tahun terdaftar di Crossref & kutipan resmi DOI = **2010**, bukan 2011 | Entri `datta2010bayesian` memakai 2010; selisih tahun dicatat di sini |
| Steorts, Hall & Ghosh (2014), *J. Survey Statist. Methodol.* 2(2), 173-193 | **TIDAK DITEMUKAN** di Crossref (judul maupun pencarian penulis) | Tidak dikutip; hanya dicatat di sini |
| Berg & Fuller (2014), *J. Survey Statist. Methodol.* 2(3), 256-283, judul "…with a constrained multinomial logit model" | **Terverifikasi dengan judul berbeda**: DOI `10.1093/jssam/smu011` = *Small Area Prediction of Proportions with Applications to the Canadian Labour Force Survey*, 2(3), **227-256** | Entri `berg2014small` memakai judul & halaman resmi; judul & halaman versi kode dicatat di sini |
| Brown, Chambers, Heady & Heasman (2001), *ONS Internal Report* | **TIDAK DITEMUKAN** (laporan internal, tanpa DOI/ISBN) | Tidak dikutip sebagai referensi bibliografis; rumus tes Brown dilaporkan apa adanya dari kode |
| Cliff & Ord (1981), *Spatial Processes: Models & Applications*, Pion London | Buku, tidak terdaftar DOI di Crossref | Tidak dikutip; hanya dicatat |
