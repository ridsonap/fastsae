# Eksplorasi & Rekomendasi Nama Fungsi Uji Diagnostik Estimasi SAE

Dokumen ini menyajikan eksplorasi nama-nama alternatif yang menarik, elegan, dan modern untuk fungsi pengujian/evaluasi hasil estimasi Small Area Estimation pada paket [`fastsae`](file:///Volumes/work/_MainR/fastsae-inla).

---

## 1. Perbandingan Alternatif Nama Fungsi

Berikut adalah kurasi nama-nama terbaik yang dikelompokkan berdasarkan gaya dan filosofi desain:

| No | Nama Fungsi | Gaya & Karakter | Keunggulan & Nuansa | Contoh Pemanggilan |
| :---: | :--- | :--- | :--- | :--- |
| **1** | **`check_sae()`** <br>*(Paling Direkomendasikan)* | **Modern & Action-Oriented** <br>(Tren *easystats* / *tidyverse*) | Sangat ramah pengguna (*user-friendly*), intuitif, mudah diketik, dan langsung terbaca sebagai perintah: *"periksa hasil SAE ini"*. | `check_sae(fit)` <br>`autoplot(check_sae(fit))` |
| **2** | **`diagnose()`** | **Sleek, Clean & Direct** | Sangat elegan, satu kata tanpa prefiks panjang. Langsung to-the-point dan terasa seperti fungsi kelas satu (*first-class API*). | `diagnose(fit)` <br>`autoplot(diagnose(fit))` |
| **3** | **`sae_eval()`** / **`evaluate()`** | **Concise & Methodological** | Lebih ringkas dari `diagnostic`. Menekankan evaluasi menyeluruh (presisi, efisiensi gain, dan uji kalibrasi bias). | `sae_eval(fit)` <br>`autoplot(sae_eval(fit))` |
| **4** | **`sae_audit()`** / **`audit_sae()`** | **Official Statistics & Quality Assurance** | Nuansa profesional yang sangat kuat, relevan dengan standar audit kelayakan publikasi data resmi di BPS / Eurostat. | `sae_audit(fit)` <br>`autoplot(sae_audit(fit))` |
| **5** | **`sae_validate()`** | **Rigorous & Scientific** | Menekankan pembuktian validitas statistik (apakah model valid dan tidak bias terhadap data survei). | `sae_validate(fit)` |
| **6** | **`sae_quality()`** | **Output-Oriented** | Berorientasi pada indikator mutu/kualitas statistik (RSE < 25%, MSE reduction, fit). | `sae_quality(fit)` |

---

## 2. Bedah 3 Kandidat Terkuat

### Kandidat A: `check_sae()` — *Pilihan Paling Populer & Modern*
* **Mengapa Menarik?**
  Mengikuti konvensi R modern (seperti paket terkemuka `performance::check_model()` atau `check_normality()`). Pengguna R era sekarang sangat menyukai fungsi dengan awalan kata kerja aktif `check_*`.
* **Kelebihan**:
  - Sangat alami diucapkan dan diketik.
  - Memiliki ruang untuk fungsi turunan jika kelak dibutuhkan (misal `check_residuals()`, `check_spatial()`).
* **Sintaks**:
  ```r
  chk <- check_sae(fit)
  print(chk)
  autoplot(chk)
  ```

---

### Kandidat B: `diagnose()` — *Pilihan Paling Bersih & Minimalis*
* **Mengapa Menarik?**
  Tidak ada redundansi kata `sae` (karena nama paketnya sendiri sudah `fastsae`). Menghindari penamaan yang panjang (*verbosity*).
* **Kelebihan**:
  - Super ringkas (hanya 8 huruf).
  - Taraf elegan setara dengan fungsi-fungsi generik R seperti `summary()`, `plot()`, `fitted()`.
* **Sintaks**:
  ```r
  diag <- diagnose(fit)
  print(diag)
  autoplot(diag)
  ```

---

### Kandidat C: `sae_audit()` — *Pilihan Paling Eksklusif untuk Lembaga Statistik*
* **Mengapa Menarik?**
  Dalam survei BPS dan statistik resmi, proses penentuan apakah estimasi level kabupaten/kota/kecamatan layak dirilis sering disebut **"Quality Audit"**.
* **Kelebihan**:
  - Memiliki diferensiasi unik dibanding paket SAE lain (seperti `sae` atau `emdi`).
  - Outputnya menyajikan checklist kelayakan rilis publikasi (*Publication Readiness Audit*).
* **Sintaks**:
  ```r
  audit <- sae_audit(fit)
  print(audit)
  ```

---

## 3. Strategi Desain API: Fungsi Utama + Alias Fleksibel

Kita dapat menetapkan **satu nama utama**, dan menyediakan **alias ramah** sehingga pengguna bebas menggunakan gaya yang mereka sukai tanpa kebingungan:

* **Nama Utama**: **`check_sae()`** (atau **`diagnose()`**)
* **Alias**:
  ```r
  # Fungsi utama
  check_sae <- function(object, W = NULL, truth = NULL, ...) { ... }

  # Alias kenyamanan
  diagnose <- check_sae
  sae_diagnostic <- check_sae
  ```
Dengan demikian, pengguna yang terbiasa mengetik `diagnose(fit)` maupun `check_sae(fit)` akan mendapatkan hasil yang sama persis!

---

## 4. User Review Required

> [!IMPORTANT]
> Dari pilihan di atas:
> 1. **`check_sae()`** (Modern & ramah)
> 2. **`diagnose()`** (Minimalis & ringkas)
> 3. **`sae_eval()`** (Padat & metodologis)
> 4. **`sae_audit()`** (Gaya jaminan mutu / audit data resmi)
>
> Manakah nama yang paling Anda sukai untuk dijadikan fungsi utama di paket `fastsae`?
