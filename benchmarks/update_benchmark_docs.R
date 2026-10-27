#!/usr/bin/env Rscript

# ==============================================================================
# Write the benchmark numbers of inst/extdata/*_benchmark.rds into
#   - vignettes/benchmarks.Rmd: the RAW_DATA / BETA_DATA arrays of the
#     interactive charts and the summary between <!-- BENCH-SUMMARY:START/END -->
#   - vignettes/hb-area-inla.Rmd: the BETA_DATA array of its scalability chart
#   - README.Rmd and pkgdown/index.md (the site home page): the tables between
#     <!-- BENCH-README:START/END -->
# (inst/extdata is excluded from the package build, so the vignette cannot read
# the files itself.) Run from the package root after updating the .rds files:
#   Rscript benchmarks/update_benchmark_docs.R
# ==============================================================================

rd <- function(f) as.data.frame(readRDS(file.path("inst/extdata", f)))
freq <- list(FH = rd("eblup_benchmark.rds"), SFH = rd("seblup_benchmark.rds"), STFH = rd("steblup_benchmark.rds"))
beta <- list(STD = rd("beta_benchmark.rds"), SPATIAL = rd("spatial_beta_benchmark.rds"),
             ST = rd("spatio_temporal_beta_benchmark.rds"))

sig <- function(x, d = 4) formatC(signif(x, d), format = "fg", digits = d, flag = "#") |> trimws() |> sub(pattern = "\\.$", replacement = "")
num <- function(x, d = 4) sig(x, d)

js_rows <- function(df, algo = NULL) {
  # the charts compare packages; the extra fastsae row with the STFH analytical MSE is left out
  df <- df[df$Method != "fastsae_mse", ]
  df <- df[order(df$n, df$Method), ]
  sprintf("    { Method: \"%s\", Median: %s, Mem: %s, IPS: %s, n: %d%s }",
          df$Method, num(df$`Median (s)`), num(df$`Memory (MB)`), num(df$`Iterations/sec`), as.integer(df$n),
          if (is.null(algo)) "" else sprintf(", Algo: \"%s\"", algo))
}

replace_block <- function(lines, start_pat, end_pat, new) {
  s <- grep(start_pat, lines)[1]
  if (is.na(s)) stop("start marker not found: ", start_pat)
  e <- s + grep(end_pat, lines[(s + 1):length(lines)])[1]
  if (is.na(e)) stop("end marker not found: ", end_pat)
  c(lines[seq_len(s)], new, lines[e:length(lines)])
}

med <- function(df, m, n) { v <- df$`Median (s)`[df$Method == m & df$n == n]; if (length(v)) v[1] else NA_real_ }
mem <- function(df, m, n) { v <- df$`Memory (MB)`[df$Method == m & df$n == n]; if (length(v)) v[1] else NA_real_ }
# Increase of the peak resident memory of the R process (benchmarks/bench_frequentist.R)
peak <- function(df, m, n) { v <- df$`Peak RSS (MB)`[df$Method == m & df$n == n]; if (length(v)) v[1] else NA_real_ }
fmt_s <- function(x) ifelse(is.na(x), "\\-", paste0(sig(x, 3), " s"))
fmt_mb <- function(x) ifelse(is.na(x), "\\-", ifelse(x >= 1024, paste0(sig(x / 1024, 3), " GB"), paste0(sig(x, 3), " MB")))
fmt_x <- function(x) ifelse(is.na(x), "\\-", paste0(format(round(x), big.mark = ","), "x"))
N <- 1000

# ---- frequentist summary table (n = 1000) ------------------------------------
freq_rows <- unlist(lapply(names(freq), function(a) {
  df <- freq[[a]]
  lab <- c(FH = "Fay-Herriot", SFH = "Spatial Fay-Herriot",
           STFH = "Spatio-Temporal FH ($T = 5$), point estimates")[[a]]
  f <- med(df, "fastsae", N); s <- med(df, "sae", N); e <- med(df, "emdi", N)
  cell <- function(m, v) if (is.na(v)) "\\-" else
    sprintf("%s (%s; %s slower)", fmt_s(v), fmt_mb(peak(df, m, N)), fmt_x(v / f))
  sprintf("| **%s** | **%s** (%s) | %s | %s |",
          lab, fmt_s(f), fmt_mb(peak(df, "fastsae", N)), cell("sae", s), cell("emdi", e))
}))
freq_table <- c(
  "| Model ($n = 1,000$ domains, 5 covariates) | `fastsae` | `sae` | `emdi` |",
  "|:---|:--:|:--:|:--:|",
  freq_rows
)
# Times only: bench's mem_alloc counts neither the memory of the external inla program
# nor that of the Stan sampler, so it cannot compare the two packages
beta_rows <- unlist(lapply(names(beta), function(a) {
  df <- beta[[a]]
  n_use <- max(df$n)
  lab <- c(STD = "Beta", SPATIAL = "Spatial Beta (Besag ICAR)", ST = "Spatio-Temporal Beta (Besag + RW1, $T = 5$)")[[a]]
  f <- med(df, "fastsae", n_use); t <- med(df, "tipsae", n_use)
  sprintf("| **%s** | %s | **%s** | %s | %s |", lab, format(n_use, big.mark = ","),
          fmt_s(f), fmt_s(t), fmt_x(t / f))
}))
beta_table <- c(
  "| Model | $n$ | `fastsae::hb_area` (INLA) | `tipsae::fit_sae` (Stan) | Ratio of times |",
  "|:---|:--:|:--:|:--:|:--:|",
  beta_rows
)
fmeta <- readRDS("inst/extdata/frequentist_benchmark_meta.rds")
date_bayes <- "2026-10-14"  # benchmarks/bench_bayes_fastsae.R (update when re-running it)
mem_note <- paste0(
  "Times are medians of `bench::mark()`. The frequentist rows were measured on ", fmeta$date,
  " with fastsae ", fmeta$versions[["fastsae"]], ", sae ", fmeta$versions[["sae"]], " and emdi ",
  fmeta$versions[["emdi"]], " on the same simulated data sets (", fmeta$cpu, ", ", fmeta$r,
  "; `benchmarks/bench_frequentist.R`), each in its own R process; memory is the increase of the peak ",
  "resident memory of that process during the fits, which includes memory allocated by compiled code ",
  "(increases of a few tens of MB vary from run to run and should not be compared). ",
  "'sae' has no analytical MSE for the spatio-temporal model, so that row compares point estimates; with ",
  "its analytical MSE, fastsae took ", fmt_s(med(freq$STFH, "fastsae_mse", N)), "."
)
freq_note <- paste0(
  "The Bayesian `fastsae` rows were re-measured on ", date_bayes, " (`benchmarks/bench_bayes_fastsae.R`), ",
  "the `tipsae` rows (4 chains of 2,000 iterations) come from an earlier run on the same data."
)

# ---- vignette ------------------------------------------------------------------
vig <- readLines("vignettes/benchmarks.Rmd")
raw <- c(
  unlist(lapply(names(freq), function(a) c(sprintf("    // %s", a), paste0(js_rows(freq[[a]], a), ",")))),
  NULL
)
raw[length(raw)] <- sub(",$", "", raw[length(raw)])
vig <- replace_block(vig, "const RAW_DATA = \\[", "^  \\];", raw)
bd <- unlist(lapply(names(beta), function(a) {
  r <- js_rows(beta[[a]])
  r[-length(r)] <- paste0(r[-length(r)], ",")
  c(sprintf("    %s: [", a), r, "    ],")
}))
bd[length(bd)] <- "    ]"
vig <- replace_block(vig, "const BETA_DATA = \\{", "^  \\};", bd)

f1000 <- function(a, m) med(freq[[a]], m, N) / med(freq[[a]], "fastsae", N)
summary_block <- c(
  "**fastsae** fits the frequentist models in compiled C++ (`RcppArmadillo`; sparse Cholesky factorisations for",
  "the SAR models) and the Bayesian models with INLA. At $n = 1,000$ domains:",
  "",
  sprintf("- **Fay-Herriot**: %s faster than `sae` and %s faster than `emdi`.", fmt_x(f1000("FH", "sae")), fmt_x(f1000("FH", "emdi"))),
  sprintf("- **Spatial Fay-Herriot**: %s faster than `sae` and %s faster than `emdi`.", fmt_x(f1000("SFH", "sae")), fmt_x(f1000("SFH", "emdi"))),
  sprintf("- **Spatio-Temporal Fay-Herriot** ($T = 5$, point estimates): %s faster than `sae`.", fmt_x(f1000("STFH", "sae"))),
  "",
  "---",
  "",
  "## Benchmark Results Table",
  "",
  freq_table,
  "",
  paste(mem_note, freq_note),
  ""
)
vig <- replace_block(vig, "<!-- BENCH-SUMMARY:START -->", "<!-- BENCH-SUMMARY:END -->", summary_block)
writeLines(vig, "vignettes/benchmarks.Rmd")
hb_vig <- readLines("vignettes/hb-area-inla.Rmd")
hb_vig <- replace_block(hb_vig, "const BETA_DATA = \\{", "^  \\};", bd)
writeLines(hb_vig, "vignettes/hb-area-inla.Rmd")

# ---- README --------------------------------------------------------------------
rme <- readLines("README.Rmd")
readme_block <- c(
  "#### 1. Frequentist EBLUP with analytical MSE (time and peak memory)",
  "",
  freq_table,
  "",
  "#### 2. Bayesian Beta models: INLA vs Stan MCMC (largest $n$ benchmarked)",
  "",
  beta_table,
  "",
  paste(mem_note, freq_note, "See the [benchmark article](https://ridsonap.is-a.dev/fastsae/articles/benchmarks.html) for all sizes."),
  ""
)
rme <- replace_block(rme, "<!-- BENCH-README:START -->", "<!-- BENCH-README:END -->", readme_block)
writeLines(rme, "README.Rmd")
home <- readLines("pkgdown/index.md")
home <- replace_block(home, "<!-- BENCH-README:START -->", "<!-- BENCH-README:END -->", readme_block)
writeLines(home, "pkgdown/index.md")
cat("Updated vignettes/benchmarks.Rmd, vignettes/hb-area-inla.Rmd, README.Rmd and pkgdown/index.md",
    "(run devtools::build_readme() next)\n")
