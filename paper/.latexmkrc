# latexmk konfigurasi untuk laporan teknis fastsae
# Pakai:  latexmk          (dari direktori paper/)
#         latexmk -C       (bersihkan)

# 1. Bangun PDF dengan pdflatex
$pdf_mode = 1;

# 2. Jalankan bibtex bila ada \bibliography (nilai 2 = paksa bibtex)
$bibtex_use = 2;

# 3. Berkas utama
@default_files = ('main.tex');

# 4. pdflatex: berhenti pada error, tampilkan nomor baris berkas
$pdflatex = 'pdflatex -interaction=nonstopmode -halt-on-error -file-line-error %O %S';

# 5. bibtex: jangan berhenti diam pada error
$bibtex = 'bibtex %O %B';

# 6. Berkas tambahan yang ikut dibersihkan oleh `latexmk -c`
$clean_ext = 'bbl blg run.xml fdb_latexmk synctex.gz';

# 7. Tidak memakai direktori keluaran: semua \input{} relatif terhadap
#    paper/ (sections/, tables/, figures/) sehingga harus tetap satu tingkat.
$out_dir = '';

# 8. Nama keluaran sesuai spesifikasi laporan: setelah build sukses,
#    main.pdf disalin ke fastsae_technical_report.pdf (identik byte per byte).
$success_cmd = 'cp -f %R.pdf fastsae_technical_report.pdf';

# 9. Gerbang mutu: sitasi/referensi tak terdefinisi atau label ganda
#    dianggap gagal build (sesuai checklist kualitas laporan).
$warnings_as_errors = 1;
