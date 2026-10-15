# Konteks untuk agen riset laporan teknis fastsae

Repo: `/Volumes/work/_MainR/fastsae-cov` (paket R `fastsae` v1.0.0).

Aturan keras untuk SEMUA agen:
- **Jangan mengubah apa pun di luar `paper/_notes/`.** Jangan sentuh `R/`, `src/`, `tests/`, `man/`, `DESCRIPTION`, `NAMESPACE`, `README*`, `vignettes/`, `benchmarks/`, `inst/`.
- Semua klaim HARUS berasal dari kode/dokumentasi yang benar-benar Anda baca. Sertakan nomor baris (`R/foo.R:123`).
- Jangan mengarang angka, rumus, atau referensi. Jika tidak yakin, tulis "TIDAK PASTI".
- Format output: Markdown, maksimal ~400 baris per file.
