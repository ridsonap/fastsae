# scripts/01_tables.R
# Menghasilkan seluruh tabel benchmark laporan teknis fastsae
# dari file sumber di inst/extdata/*.rds.
#
# Jalankan dari folder technical-report/:   Rscript scripts/01_tables.R
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

# ---------------------------------------------------------------------------
# Tabel perbandingan engine hb_area() pada model spasial (ICAR)
# sumber: hb_engine_benchmark*.rds (benchmarks/bench_hb_area_engines.R)
# ---------------------------------------------------------------------------
engine_tables <- function() {
  files <- list.files(root, pattern = "^hb_engine_benchmark.*\\.rds$")
  if (length(files) == 0) return(invisible(NULL))
  B <- lapply(files, function(f) readRDS(file.path(root, f)))
  tm <- do.call(rbind, lapply(B, `[[`, "timings"))
  ag <- do.call(rbind, lapply(B, `[[`, "agreement"))
  model_lab <- c(gaussian = "Gaussian", beta = "Beta")
  src <- sprintf(SRC, texesc(paste(files, collapse = ", ")))
  get <- function(m, d, e, n) {
    v <- tm$time_warm[tm$machine == m & tm$design == d & tm$engine == e & tm$n == n]
    if (length(v) == 0) "--" else fmt(v[1])
  }
  keys <- unique(tm[, c("machine", "design", "n")])
  keys <- keys[order(keys$machine, match(keys$design, names(model_lab)), keys$n), ]
  rows <- vapply(seq_len(nrow(keys)), function(i) {
    k <- keys[i, ]
    gpu <- intersect(c("metal", "cuda"), tm$engine[tm$machine == k$machine])
    gpu_txt <- if (length(gpu) == 0) "--" else sprintf("%s (%s)", get(k$machine, k$design, gpu[1], k$n), gpu[1])
    sprintf("%s & %s & %s & %s & %s & %s & %s \\\\", k$machine, model_lab[[k$design]], k$n,
            get(k$machine, k$design, "inla", k$n), get(k$machine, k$design, "cpu", k$n),
            gpu_txt, get(k$machine, k$design, "tipsae", k$n))
  }, character(1))
  txt <- sprintf(
"\\begin{table}[htbp]
\\centering
\\caption{Waktu komputasi (detik, run \\emph{warm}) \\code{hb\\_area()} dengan efek ICAR (Besag) pada grid beraturan: INLA, JAX/NumPyro NUTS di CPU dan di GPU, serta \\pkg{tipsae} (Stan, hanya model Beta). JAX: 2 rantai $\\times$ 1000 draw setelah 500 warmup. %s}
\\label{tab:engine-time}
\\small
\\begin{tabular}{llrrrrr}
\\toprule
Mesin & Model & $D$ & INLA & JAX CPU & JAX GPU & tipsae \\\\
\\midrule
%s
\\bottomrule
\\end{tabular}
\\end{table}
", src, paste(rows, collapse = "\n"))
  writeLines(txt, file.path(outdir, "tab-engine-time.tex"))

  eng_lab <- c(cpu = "JAX CPU", metal = "JAX Metal", cuda = "JAX CUDA", tipsae = "tipsae")
  g <- unique(ag[, c("machine", "design", "engine")])
  g <- g[order(g$machine, match(g$design, names(model_lab)), match(g$engine, names(eng_lab))), ]
  rows2 <- vapply(seq_len(nrow(g)), function(i) {
    k <- g[i, ]
    a <- ag[ag$machine == k$machine & ag$design == k$design & ag$engine == k$engine, ]
    r <- tm$rhat_max[tm$machine == k$machine & tm$design == k$design & tm$engine == k$engine]
    sprintf("%s & %s & %s & %s & %s--%s & %s \\\\", k$machine, model_lab[[k$design]], eng_lab[[k$engine]],
            fmt(min(a$cor_hb)), fmt(min(a$median_ratio_sd)), fmt(max(a$median_ratio_sd)), fmt(max(r, na.rm = TRUE)))
  }, character(1))
  txt2 <- sprintf(
"\\begin{table}[htbp]
\\centering
\\caption{Kesesuaian tiap engine terhadap INLA pada mesin dan data yang sama, lintas seluruh $D$: korelasi minimum rerata posterior, rentang median rasio sd posterior (engine/INLA), dan $\\hat R$ terpisah (\\emph{split}) terbesar. %s}
\\label{tab:engine-equiv}
\\small
\\begin{tabular}{lllrrr}
\\toprule
Mesin & Model & Engine & $\\min r$ & Rasio sd & $\\max \\hat R$ \\\\
\\midrule
%s
\\bottomrule
\\end{tabular}
\\end{table}
", src, paste(rows2, collapse = "\n"))
  writeLines(txt2, file.path(outdir, "tab-engine-equiv.tex"))
}
engine_tables()

# ---------------------------------------------------------------------------
# Tabel simulasi Monte Carlo eblup_stfh dengan sel tak tersampel
# sumber: benchmarks/output/sim_stfh_unsampled.rds (benchmarks/sim_stfh_unsampled.R)
# ---------------------------------------------------------------------------
stfh_unsampled_table <- function() {
  f <- file.path(root, "..", "..", "benchmarks", "output", "sim_stfh_unsampled.rds")
  if (!file.exists(f)) return(invisible(NULL))
  sim <- readRDS(f)
  sm <- sim$summary
  st <- sim$settings
  lab <- c(sampled = "tersampel", gap = "periode hilang",
           "unsampled domain" = "domain tanpa sampel")
  rows <- sprintf("%s & %d & %s & %s & %s & %s & %s \\\\", lab[sm$type], sm$cells,
                  fmt(sm$mse_eblup), fmt(sm$mse_synthetic), fmt(sm$rb_analytical_pct),
                  fmt(sm$rb_bootstrap_pct), fmt(sm$coverage_pct))
  txt <- sprintf(
"\\begin{table}[htbp]
\\centering
\\caption{Simulasi Monte Carlo berbasis model untuk \\code{eblup\\_stfh} dengan sel tanpa estimasi langsung: $D=%d$, $T=%d$, %d replikasi (bootstrap $B=%d$ pada %d replikasi pertama), $(\\sigma_1^2,\\rho_1,\\sigma_2^2,\\rho_2)=(%s)$. MSE empiris EBLUP dan estimator sintetis $x'\\hat\\beta$; bias relatif (\\%%) estimator MSE analitik dan bootstrap; cakupan (\\%%) interval $\\hat\\theta\\pm1{,}96\\sqrt{\\widehat{\\mathrm{MSE}}}$ analitik, dirata-ratakan per jenis sel. Sumber: \\texttt{benchmarks/output/sim\\_stfh\\_unsampled.rds} (\\texttt{benchmarks/sim\\_stfh\\_unsampled.R}); angka dihitung ulang oleh \\texttt{scripts/01\\_tables.R}.}
\\label{tab:stfh-unsampled-sim}
\\small
\\begin{tabular}{lrrrrrr}
\\toprule
Jenis sel & Sel & \\multicolumn{2}{c}{MSE empiris} & \\multicolumn{2}{c}{RB (\\%%)} & Cakupan \\\\
\\cmidrule(lr){3-4} \\cmidrule(lr){5-6}
 & & EBLUP & sintetis & analitik & bootstrap & (\\%%) \\\\
\\midrule
%s
\\bottomrule
\\end{tabular}
\\end{table}
", st$D, st$T, st$R, st$B, st$R_boot, paste(format(st$pars), collapse = ", "), paste(rows, collapse = "\n"))
  writeLines(txt, file.path(outdir, "tab-stfh-unsampled-sim.tex"))
}
stfh_unsampled_table()

cat("Selesai. Berkas tabel di", normalizePath(outdir), ":\n")
print(list.files(outdir, pattern = "\\.tex$"))
