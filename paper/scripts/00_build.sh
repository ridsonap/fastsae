#!/bin/sh
# scripts/00_build.sh -- membangun laporan teknis fastsae
#
# Pemakaian (dari folder paper/):
#   sh scripts/00_build.sh          # hanya membangun PDF dari sumber yang ada
#   sh scripts/00_build.sh --all    # regenerasi tabel & figur dulu, lalu membangun
#
# Keluaran: main.pdf dan fastsae_technical_report.pdf (salinan identik).
# Konfigurasi latexmk ada di .latexmkrc (pdf_mode=1, bibtex_use=2, nonstopmode
# + halt-on-error + file-line-error).

set -e
cd "$(dirname "$0")/.."

if [ "$1" = "--all" ]; then
  Rscript scripts/01_tables.R
  Rscript scripts/02_figures.R
fi

latexmk -pdf main.tex
cp -f main.pdf fastsae_technical_report.pdf
ls -l main.pdf fastsae_technical_report.pdf
