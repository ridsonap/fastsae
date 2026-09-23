# Analisis & Penjelasan Arsitektur: Perbedaan `ebp`, `ebp_area`, dan `ebp_spatial`

Dokumen ini menjelaskan latar belakang, perbandingan teknis, dan rancangan arsitektur antarmuka fungsi estimasi Bayesian / EBP pada paket [`fastsae`](file:///Volumes/work/_MainR/fastsae-inla).

---

## 1. Latar Belakang & Asal-Usul

Saat perancangan awal model Bayesian INLA di `fastsae`, terdapat dua kebutuhan desain:
1. **Pola Eksplisit Berbasis Level Data** (Konsisten dengan konvensi paket `fastsae` yang sudah ada: `eblup_fh`, `eblup_sfh`, `eblup_stfh`, `eblup_bhf`):
   - Di `fastsae`, nama fungsi membedakan tipe model dan level data secara jelas:
     - Non-spasial area: `eblup_fh()`
     - Spasial area: `eblup_sfh()`
   - Mengikuti analogi tersebut, dibentuklah:
     - `ebp_area()`: untuk model area-level umum.
     - `ebp_spatial_area()`: khusus untuk model spasial area-level dengan matriks bobot $W$.
2. **Ergonomi & Kemudahan Pengguna (Shorthand)**:
   - Pengguna sering kali menginginkan nama fungsi yang ringkas dan cepat diketik, mirip paket lain (seperti `emdi::ebp` atau `sae::eblupFH`).
   - Oleh karena itu, dibuatlah alias:
     - `ebp <- ebp_area`
     - `ebp_spatial <- ebp_spatial_area`

---

## 2. Diagram Hubungan & Hirarki Fungsi

Secara internal di dalam kode R saat ini, struktur hubungannya adalah sebagai berikut:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                                 ebp_area()                                  │
│                          (Core Engine Estimasi Area)                        │
│                                                                             │
│  - Parameter: spatial = c("none", "bym2", "bym", "besag", "generic1")      │
│  - Default: spatial = "none" (Model IID / Non-Spasial)                      │
│  - Matriks W: Opsional jika "none", Wajib jika spatial != "none"            │
│  - Method: "inla" (Bayesian) atau "laplace" (Frequentist GLMM via lme4)     │
└──────────────────────────────────────┬──────────────────────────────────────┘
                                       │
            ┌──────────────────────────┴──────────────────────────┐
            │                                                     │
     [Alias Langsung]                                     [Dedicated Wrapper]
            │                                                     │
            ▼                                                     ▼
    ┌───────────────┐                                 ┌───────────────────────┐
    │     ebp()     │                                 │   ebp_spatial_area()  │
    │  (= ebp_area) │                                 │                       │
    └───────────────┘                                 │ - W WAJIB diisi       │
                                                      │ - Default: BYM2       │
                                                      │ - Memanggil ebp_area()│
                                                      └───────────┬───────────┘
                                                                  │
                                                           [Alias Langsung]
                                                                  │
                                                                  ▼
                                                      ┌───────────────────────┐
                                                      │     ebp_spatial()     │
                                                      │ (= ebp_spatial_area)  │
                                                      └───────────────────────┘
```

---

## 3. Matriks Perbandingan Teknis

| Fitur / Karakteristik | `ebp_area()` | `ebp()` | `ebp_spatial_area()` | `ebp_spatial()` |
| :--- | :--- | :--- | :--- | :--- |
| **Status Kode** | **Core Function** | Alias dari `ebp_area` | **Spatial Wrapper** | Alias dari `ebp_spatial_area` |
| **Tujuan Utama** | Estimator area terpadu (non-spasial & spasial) | Alias singkat `ebp_area` | Estimator khusus model spasial | Alias singkat `ebp_spatial_area` |
| **Model Spasial Default** | `spatial = "none"` | `spatial = "none"` | `spatial = "bym2"` | `spatial = "bym2"` |
| **Kewajiban Argumen `W`** | Opsional (hanya wajib jika `spatial != "none"`) | Opsional (hanya wajib jika `spatial != "none"`) | **Wajib** (`missing(W)` $\rightarrow$ error) | **Wajib** (`missing(W)` $\rightarrow$ error) |
| **Pilihan Model Spasial** | `"none"`, `"bym2"`, `"bym"`, `"besag"`, `"generic1"` | `"none"`, `"bym2"`, `"bym"`, `"besag"`, `"generic1"` | `"bym2"`, `"bym"`, `"besag"`, `"generic1"`, `"slm"` | `"bym2"`, `"bym"`, `"besag"`, `"generic1"`, `"slm"` |
| **Dukungan Metode** | `"inla"` dan `"laplace"` | `"inla"` dan `"laplace"` | `"inla"` | `"inla"` |
| **Objek Output** | S3 `fastsae_ebp_area` | S3 `fastsae_ebp_area` | S3 `fastsae_ebp_area` | S3 `fastsae_ebp_area` |

---

## 4. Apakah Terjadi Redundansi?

Secara fungsional:
* `ebp` dan `ebp_area` adalah **100% objek yang sama** (`identical(ebp, ebp_area) == TRUE`).
* `ebp_spatial` dan `ebp_spatial_area` adalah **100% objek yang sama** (`identical(ebp_spatial, ebp_spatial_area) == TRUE`).
* `ebp_area()` sebenarnya **sudah mampu menangani model spasial** jika pengguna memberikan `spatial = "bym2"` dan `W = proxmat`.

Lalu mengapa `ebp_spatial()` tetap dibuat?
1. **Mencegah Kesalahan Pengguna**: Pada `ebp_spatial()`, jika pengguna lupa memasukkan `W`, fungsi langsung memberi pesan peringatan yang ramah dan jelas sebelum memanggil INLA.
2. **Ergonomi**: Pengguna tidak perlu mengetik `spatial = "bym2"` secara manual karena sudah menjadi nilai *default*.
3. **Keselarasan dengan EBLUP**: Pengguna yang biasa menggunakan pasangan `eblup_fh()` (non-spasial) dan `eblup_sfh()` (spasial) merasa alami jika memiliki pasangan `ebp()` (non-spasial) dan `ebp_spatial()` (spasial).

---

## 5. Opsi Desain untuk Pengguna (User Choice)

Ada 3 pilihan arah perapihan API paket yang dapat kita pilih:

### Opsi A: Pertahankan Status Quo (Rekomendasi)
* Tetap sediakan `ebp_area` & `ebp_spatial_area` sebagai nama fungsi formal/deskriptif, dan `ebp` & `ebp_spatial` sebagai *shortcut*.
* **Kelebihan**: Kompatibel dengan semua skrip yang sudah ditulis; pengguna bebas memilih gaya penamaan (panjang/formal vs pendek/ergonomis).
* **Kekurangan**: Ada 4 nama fungsi di halaman dokumentasi / auto-complete RStudio.

### Opsi B: Sederhanakan Hanya Menjadi 2 Fungsi Utama (`ebp` dan `ebp_spatial`)
* Jadikan `ebp()` dan `ebp_spatial()` sebagai fungsi utama.
* Hapus atau *deprecate* sufiks `_area` (`ebp_area` dan `ebp_spatial_area`).
* **Kelebihan**: Sangat bersih, intuitif, dan tidak ada duplikasi di auto-complete.
* **Catatan**: Jika kelak ditambahkan `ebp_unit()`, fungsi non-spasial area tetap bernama `ebp()`.

### Opsi C: Satukan Semuanya Menjadi 1 Fungsi Tunggal (`ebp` atau `ebp_area`)
* Cukup satu fungsi saja: `ebp(formula, data, spatial = "none", W = NULL, ...)`
* Hapus `ebp_spatial`. Jika ingin model spasial, cukup panggil `ebp(..., spatial = "bym2", W = proxmat)`.
* **Kelebihan**: Desain paling minimalis (*unified API*).
* **Kekurangan**: Berbeda pola dengan `eblup_fh` vs `eblup_sfh` yang terpisah di fungsi EBLUP.
