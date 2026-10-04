# scripts/01_tables.R
# Menghasilkan seluruh tabel benchmark laporan teknis fastsae
# dari file sumber di inst/extdata/*.rds.
#
# Jalankan dari folder paper/:   Rscript scripts/01_tables.R
# Output: tables/*.tex
#
# ATURAN: tidak ada angka yang diketik manual; semua dihitung dari .rds.
# Angka ditampilkan dengan formatC(format="g", digits=6) dan dicatat di caption.

root <- if (dir.exists("../inst/extdata")) "../inst/extdata" else "inst/extdata"
outdir <- if (dir.exists("tables") || dir.create("tables", showWarnings = FALSE)) "tables" else "tables"
dir.create(outdir, showWarnings = FALSE)

fmt <- function(x) formatC(as.numeric(x), format = "g", digits = 6)
SRC <- "Sumber: file \\texttt{%s} di \\texttt{inst/extdata} (repo \\texttt{fastsae}); angka dihitung ulang oleh \\texttt{scripts/01\\_tables.R}, ditampilkan hingga 6 digit signifikan."
# escape karakter khusus LaTeX untuk teks yang disisipkan ke dalam \texttt{...}
texesc <- function(s) gsub("_", "\\\\_", s, fixed = FALSE)

read_b <- function(f) as.data.frame(readRDS(file.path(root, f)))

# ---------------------------------------------------------------------------
# helper: tabel waktu/memori (median detik / MB) per n untuk tiap file
# ---------------------------------------------------------------------------
bench_table <- function(file, label, caption_extra = NULL, what = c("time", "mem"),
                        width = NULL) {
  what <- match.arg(what)
  d <- read_b(file)
  col <- if (what == "time") "Median (s)" else "Memory (MB)"
  meth <- unique(d$Method)
  ns <- sort(unique(d$n))
  M <- sapply(meth, function(m) {
    v <- d[d$Method == m, col][match(ns, d$n[d$Method == m])]
    v
  })
  if (is.null(dim(M))) M <- matrix(M, ncol = length(meth))
  colnames(M) <- meth
  hdr <- paste(meth, collapse = " & ")
  rows <- paste(ns, "&", apply(M, 1, function(r) paste(fmt(r), collapse = " & ")), "\\\\")
  body <- paste(rows, collapse = "\n")
  numc <- paste(rep("r", length(meth) + 1), collapse = "")
  cap <- if (what == "time")
    sprintf("Waktu median per iterasi (detik) menurut jumlah area $n$.")
  else sprintf("Alokasi memori total per iterasi (MB) menurut jumlah area $n$.")
  if (!is.null(caption_extra)) cap <- paste(cap, caption_extra)
  cap <- paste(cap, sprintf(SRC, texesc(file)))
  txt <- sprintf(
"\\begin{table}[htbp]
\\centering
\\caption{%s}
\\label{%s}
\\small
\\begin{tabular}{l%s}
\\toprule
$n$ & %s \\\\
\\midrule
%s
\\bottomrule
\\end{tabular}
\\end{table}
", cap, label, numc, hdr, body)
  writeLines(txt, file.path(outdir, paste0(gsub(":", "-", label), ".tex")))
  invisible(M)
}

# --- waktu: FH, spatial FH, spatio-temporal FH -------------------------------
bench_table("eblup_benchmark.rds", "tab:time-fh", what = "time")
bench_table("seblup_benchmark.rds", "tab:time-sfh", what = "time")
bench_table("steblup_benchmark.rds", "tab:time-stfh", what = "time",
            caption_extra = "Model spasio-temporal dengan $T$ periode waktu; waktu termasuk estimasi komponen varian.")

# --- memori: FH, spatial FH, spatio-temporal FH ------------------------------
bench_table("eblup_benchmark.rds", "tab:mem-fh", what = "mem")
bench_table("seblup_benchmark.rds", "tab:mem-sfh", what = "mem")
bench_table("steblup_benchmark.rds", "tab:mem-stfh", what = "mem")

# --- model beta -------------------------------------------------------------
bench_table("beta_benchmark.rds", "tab:time-beta", what = "time",
            caption_extra = "Model area-level beta (\\texttt{fastsae::hb\\_area} vs \\texttt{tipsae::fit\\_sae}).")
bench_table("spatial_beta_benchmark.rds", "tab:time-sbeta", what = "time",
            caption_extra = "Model area-level beta dengan efek spasial.")
bench_table("spatio_temporal_beta_benchmark.rds", "tab:time-stbeta", what = "time",
            caption_extra = "Model beta spasio-temporal ($T=5$); data tersedia hanya untuk $n=30$ dan $n=50$.")

# --- Bayes normal ------------------------------------------------------------
bench_table("hb_benchmark.rds", "tab:time-hb", what = "time",
            caption_extra = "Model hierarki Bayes area-level Gaussian (\\texttt{hb\\_area}) vs \\texttt{hbsae::fSAE.Area}.")
bench_table("hb_benchmark.rds", "tab:mem-hb", what = "mem",
            caption_extra = "Model hierarki Bayes area-level Gaussian.")

# ---------------------------------------------------------------------------
# Tabel rasio kecepatan (pembanding / fastsae)
# ---------------------------------------------------------------------------
ratio_table <- function() {
  specs <- list(
    list(f = "eblup_benchmark.rds",  nm = "FH area-level (\\texttt{eblup\\_fh})"),
    list(f = "seblup_benchmark.rds", nm = "FH spasial (\\texttt{eblup\\_sfh})"),
    list(f = "steblup_benchmark.rds", nm = "FH spasio-temporal (\\texttt{eblup\\_stfh})"),
    list(f = "beta_benchmark.rds",   nm = "Beta area-level (\\texttt{hb\\_area})"),
    list(f = "spatial_beta_benchmark.rds", nm = "Beta spasial"),
    list(f = "spatio_temporal_beta_benchmark.rds", nm = "Beta spasio-temporal"),
    list(f = "hb_benchmark.rds",     nm = "Hierarki Bayes Gaussian")
  )
  rows <- character(0)
  for (s in specs) {
    d <- read_b(s$f)
    ns <- sort(unique(d$n))
    for (m in setdiff(unique(d$Method), "fastsae")) {
      r <- sapply(ns, function(n) {
        a <- d$`Median (s)`[d$Method == "fastsae" & d$n == n]
        b <- d$`Median (s)`[d$Method == m & d$n == n]
        b / a
      })
      rows <- c(rows, sprintf("%s & %s & %s & %s & %s & %s \\\\",
                              s$nm, m, fmt(r[1]), fmt(r[length(r)]),
                              fmt(min(r)), fmt(max(r))))
      rows <- c(rows, sprintf(" & \\multicolumn{5}{r}{\\emph{rata-rata lintas $n$: %s$\\times$}} \\\\[-2pt]",
                              fmt(mean(r))))
    }
  }
  txt <- sprintf(
"\\begin{table}[htbp]
\\centering
\\caption{Rasio kecepatan $=$ median waktu pembanding $/$ median waktu \\texttt{fastsae} (semakin besar semakin cepat \\texttt{fastsae}; nilai $<1$ berarti pembanding lebih cepat). Rasio dihitung per $n$ lalu diringkas.}
\\label{tab:ratio}
\\small
\\begin{tabular}{llrrrr}
\\toprule
Model & Pembanding & $n$ pertama & $n$ terbesar & min & max \\\\
\\midrule
%s
\\bottomrule
\\end{tabular}
\\end{table}
", paste(rows, collapse = "\n"))
  writeLines(txt, file.path(outdir, "tab-ratio.tex"))
}
ratio_table()

# ---------------------------------------------------------------------------
# Tabel ekuivalensi numerik (fastsae vs tipsae, n = 1000)
# ---------------------------------------------------------------------------
equiv_table <- function() {
  d <- read_b("beta_estimates_n1000.rds")
  ef <- d[["estimasi fastsae"]]; et <- d[["estimasi tipsae"]]
  mf <- d[["mse fastsae"]];      mt <- d[["mse tipsae"]]
  val <- c(
    "Pearson $r$ (estimasi)"          = cor(ef, et, method = "pearson"),
    "Spearman $\\rho$ (estimasi)"     = cor(ef, et, method = "spearman"),
    "MAE (estimasi)"                  = mean(abs(ef - et)),
    "RMSD (estimasi)"                 = sqrt(mean((ef - et)^2)),
    "Galat absolut maksimum"          = max(abs(ef - et)),
    "Pearson $r$ (MSE)"               = cor(mf, mt, method = "pearson"),
    "Median selisih mutlak MSE"       = median(abs(mf - mt)),
    "Median rasio MSE fastsae/tipsae" = median(mf / mt)
  )
  rows <- sprintf("%s & %s \\\\", names(val), fmt(val))
  txt <- sprintf(
"\\begin{table}[htbp]
\\centering
\\caption{Ekuivalensi numerik antara \\texttt{fastsae} dan \\texttt{tipsae} pada $n=1000$ area (indikator beta). Semua metrik dihitung dari pasangan kolom pada file sumber.}
\\label{tab:equiv}
\\small
\\begin{tabular}{lr}
\\toprule
Metrik & Nilai \\\\
\\midrule
%s
\\bottomrule
\\end{tabular}
\\end{table}
", paste(rows, collapse = "\n"))
  writeLines(txt, file.path(outdir, "tab-equiv.tex"))
}
equiv_table()

# ---------------------------------------------------------------------------
# Tabel ringkasan memori lintas pembanding (mean per file)
# ---------------------------------------------------------------------------
mem_summary <- function() {
  files <- c("eblup_benchmark.rds", "seblup_benchmark.rds", "steblup_benchmark.rds",
             "beta_benchmark.rds", "spatial_beta_benchmark.rds",
             "spatio_temporal_beta_benchmark.rds", "hb_benchmark.rds")
  rows <- character(0)
  for (f in files) {
    d <- read_b(f)
    for (m in unique(d$Method)) {
      v <- d$`Memory (MB)`[d$Method == m]
      rows <- c(rows, sprintf("\\texttt{%s} & %s & %s & %s & %s \\\\",
                              texesc(sub("_benchmark\\.rds$", "", f)), m,
                              fmt(mean(v)), fmt(min(v)), fmt(max(v))))
    }
  }
  txt <- sprintf(
"\\begin{table}[htbp]
\\centering
\\caption{Ringkasan alokasi memori total per iterasi (MB): rerata, minimum dan maksimum lintas seluruh $n$ pada tiap file sumber. Kolom ``maks'' adalah nilai pada $n$ terbesar pada file tersebut.}
\\label{tab:mem-summary}
\\small
\\begin{tabular}{llrrr}
\\toprule
File & Metode & Rerata & Min & Maks \\\\
\\midrule
%s
\\bottomrule
\\end{tabular}
\\end{table}
", paste(rows, collapse = "\n"))
  writeLines(txt, file.path(outdir, "tab-mem-summary.tex"))
}
mem_summary()

cat("Selesai. Berkas tabel di", normalizePath(outdir), ":\n")
print(list.files(outdir, pattern = "\\.tex$"))
