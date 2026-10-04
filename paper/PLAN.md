# PLAN — Kerangka laporan teknis `fastsae`

Bahasa: **Indonesia**, akademis tetapi mudah dibaca.
Mesin: `pdflatex` + `bibtex` (natbib, author-year), dikendalikan `latexmk` (`.latexmkrc`).
Kerangka kerja: `algorithm2e` (pseudocode), `booktabs` (tabel), `tikz` (grafik tetangga & DAG),
`listings` (cuplikan kode R).

## Keputusan struktur

- `main.tex` berisi preambul + `\input` tiap berkas di `sections/`.
- Bab model (bagian 5) dipecah menjadi satu berkas per keluarga model agar tiap berkas
  tetap kecil dan mudah diverifikasi terhadap kode:
  `05_fh`, `05_sfh`, `05_stfh`, `05_bhf`, `05_twofold`, `05_hb_area`, `05_hb_unit`,
  `05_hb_twofold`.
- Utilitas (benchmarking, diagnostik, perbandingan, ekspor, peta) = bagian 7 (`07_utils.tex`).
- Seluruh angka benchmark dihasilkan oleh skrip di `scripts/` yang membaca
  `inst/extdata/*.rds` → menghasilkan `figures/*.pdf` dan `tables/*.tex`.
  **Tidak ada angka yang diketik manual di teks.**
- Setiap tabel/figur diberi keterangan sumber (`\caption{... Sumber: ...}`).

## Outline akhir

| Bagian | Berkas | Isi |
|---|---|---|
| 1. Pendahuluan | `sections/01_intro.tex` | Masalah SAE, estimasi langsung vs model-based, kontribusi `fastsae`, arsitektur (R / Rcpp-Armadillo / OpenMP / INLA / S3), instalasi termasuk INLA, struktur repo |
| 2. Desain perangkat lunak | `sections/02_design.tex` | Inti C++: Fisher scoring, identitas Woodbury/blok (di mana & mana yang tidak), strategi memori $O(m^2)$, bootstrap paralel OpenMP, pembagian kerja R–C++, alur data, kelas `fastsae` |
| 3. Benchmark | `sections/03_benchmark.tex` | Metodologi (hardware/R/iterasi/seed/n, dari `benchmarks/` + `inst/extdata`), hasil waktu & memori FH/spatial/ST/Bayesian, skala dengan $n$ (skala log), ekuivalensi numerik vs pembanding, keadilan & kaveat, kontradiksi README |
| 4. Latar statistik | `sections/04_background.tex` | Tabel notasi; LMM; GLS & BLUP; EBLUP; ML & REML; Fisher scoring; dekomposisi MSE $g_1{+}g_2{+}g_3$; bootstrap parametrik/nonparametrik; domain unsampled & prediksi sintetis; benchmarking internal/external |
| 5. Per model | `sections/05_*.tex` | Tiap fungsi: (a) setting data; (b) persamaan model; (c) prediktor; (d) estimasi komponen varians; (e) MSE; (f) pseudocode `algorithm2e` yang dicocokkan baris-per-baris; (g) signature, argumen, objek return, contoh; (h) referensi kunci |
| 5.1 `eblup_fh` | `05_fh.tex` | FH area-level (Fisher scoring ML/REML, $g_1g_2g_3$, cabang unsampled) |
| 5.2 `eblup_sfh` | `05_sfh.tex` | FH spasial error-components/SAR |
| 5.3 `eblup_stfh` | `05_stfh.tex` | FH spasio-temporal (spasial konstan-waktu + AR(1)) |
| 5.4 `eblup_bhf` | `05_bhf.tex` | BHF unit-level (lme4 + prediktor C++, bootstrap) |
| 5.5 `eblup_twofold` | `05_twofold.tex` | Model sub-area two-fold (Torabi–Rao), derivasi matriks $g_3$ |
| 5.6 `hb_area` | `05_hb_area.tex` | Bayes area-level: 6 likelihood, 6 struktur spasial, 5 temporal, 7 interaksi |
| 5.7 `hb_unit` | `05_hb_unit.tex` | Bayes unit-level |
| 5.8 `hb_twofold` | `05_hb_twofold.tex` | Bayes two-fold |
| 6. Bab Bayes | `sections/06_bayesian.tex` | Latent Gaussian & INLA; GMRF; CAR/ICAR/proper CAR/SAR; BYM & BYM2; PC prior & log-gamma; prior default paket; likelihood & link; skala rata-rata area; unsampled; **TikZ grafik tetangga (W, D, Q)** & **DAG model hierarki** |
| 7. Utilitas | `sections/07_utils.tex` | Diagnostik (Brown, GOF, Moran's I, flag), perbandingan, benchmarking (rumus 4 metode + hierarki), ekspor, pemetaan |
| 8. Lampiran | `sections/08_appendix.tex` | Derivasi: Woodbury, skor REML & matriks informasi, $g_1/g_2/g_3$, $g_3$ turunan untuk two-fold; glosarium; reproduktibilitas (`sessionInfo`, perintah regenerasi) |

## Berkas pendukung

- `references.bib` — hanya referensi terverifikasi (batch 1–3, DOI dikonfirmasi via Crossref).
- `TODO_REFERENCES.md` — entri yang tidak dapat diverifikasi + alasannya.
- `DISCREPANCIES.md` — semua perbedaan implementasi vs rumus/dokumentasi publik.
- `CODE_MAP.md` — peta fungsi → file → model → estimasi → MSE → C++/OpenMP.
- `scripts/*.R` — regenerasi tabel & gambar; `figures/`, `tables/`.

## Urutan pengerjaan

1. ✅ Inventaris (CODE_MAP, dump data, catatan per file kode).
2. ✅ Verifikasi referensi batch 1–2; batch 3 (paket R) berjalan.
3. Susun DISCREPANCIES.md & TODO_REFERENCES.md.
4. Skrip `scripts/` → jalankan → `figures/` + `tables/`.
5. Tulis bagian 1–8.
6. `latexmk` sampai bersih (tanpa undefined citation/reference/citation undefined).
7. Checklist kualitas + ringkasan akhir.
