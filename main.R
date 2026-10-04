library(devtools)

# chmod +x configure
Rcpp::compileAttributes()
devtools::document()
devtools::check()
devtools::test()
load_all()

usethis::use_citation()


# Paper -------------------------------------------------------------------
# Pada hari Anda memutuskan untuk mengirimkan (submit) artikel:
#
# 1. Pastikan tanggal submisi pada berkas fastsae.Rmd:3 sesuai dengan tanggal pengiriman (date: 'YYYY-MM-DD').
# 2. Masuk ke direktori artikel di R console dan buat arsip ZIP submisi:
#   setwd("paper/rjournal")
# rjtools::zip_paper(".")
#
# 3. Lakukan submisi melalui fungsi pembantu resmi:
#   rjtools::submit_rjournal(".")
#   (Fungsi ini akan membuka Google Form pengajuan resmi The R Journal dan memandu pengunggahan arsip .zip yang telah
#   divalidasi).

