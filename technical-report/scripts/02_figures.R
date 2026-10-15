# scripts/02_figures.R
# Menghasilkan seluruh figur benchmark laporan teknis fastsae
# dari file sumber di inst/extdata/*.rds.
#
# Jalankan dari folder technical-report/:   Rscript scripts/02_figures.R
# Output: figures/*.pdf
# Versi bahasa Inggris untuk artikel R Journal:
#   Rscript scripts/02_figures.R --lang=en --out=rjournal/figures
#
# ATURAN: tidak ada angka yang diketik manual; sumbu dan kurva dibangun
# langsung dari .rds. Skala log pada kedua sumbu bila diminta.

root <- if (dir.exists("../inst/extdata")) "../inst/extdata" else "inst/extdata"

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(name, default) {
  a <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(a)) sub(paste0("^--", name, "="), "", a[1]) else default
}
lang <- arg_val("lang", "id")
out_dir <- arg_val("out", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fig <- function(f) file.path(out_dir, f)

# label sumbu dan judul panel (id = laporan teknis, en = artikel R Journal)
labs <- list(
  id = c(n = "jumlah area n", time = "waktu median (detik, log)", mem = "memori (MB, log)",
         ratio = "rasio waktu (pembanding / fastsae)", fh = "(a) Fay-Herriot area-level",
         sfh = "(b) Fay-Herriot spasial", stfh = "(c) Fay-Herriot spasio-temporal",
         hb = "(d) Hierarki Bayes Gaussian", beta_a = "(a) Beta area-level",
         beta_b = "(b) Beta spasial", beta_c = "(c) Beta spasio-temporal",
         beta_d = "(d) Beta area-level", est_t = "estimasi tipsae", est_f = "estimasi fastsae",
         mse_t = "MSE tipsae (log)", mse_f = "MSE fastsae (log)", pred = "(a) Prediktor area",
         mse = "(b) MSE", done = "Selesai. Berkas figur di",
         D = "jumlah area D", eng_time = "waktu warm (detik, log)",
         eng_g = "(a) Fay-Herriot spasial (ICAR)", eng_b = "(b) Beta spasial (ICAR)"),
  en = c(n = "number of areas n", time = "median time (seconds, log)", mem = "memory allocated (MB, log)",
         ratio = "time ratio (comparator / fastsae)", fh = "(a) Fay-Herriot",
         sfh = "(b) Spatial Fay-Herriot", stfh = "(c) Spatio-temporal Fay-Herriot",
         hb = "(d) Hierarchical Bayes Gaussian", beta_a = "(a) Beta",
         beta_b = "(b) Spatial Beta", beta_c = "(c) Spatio-temporal Beta",
         beta_d = "(d) Beta", est_t = "tipsae estimate", est_f = "fastsae estimate",
         mse_t = "tipsae MSE (log)", mse_f = "fastsae MSE (log)", pred = "(a) Area predictions",
         mse = "(b) MSE", done = "Done. Figures in",
         D = "number of domains D", eng_time = "warm time (seconds, log)",
         eng_g = "(a) Spatial Fay-Herriot (ICAR)", eng_b = "(b) Spatial Beta (ICAR)")
)[[lang]]
L <- function(id) unname(labs[[id]])

read_b <- function(f) as.data.frame(readRDS(file.path(root, f)))

# warna metode yang konsisten di seluruh laporan
pal <- c(fastsae = "#1b6ca8", sae = "#c0392b", emdi = "#7f8c8d",
         tipsae = "#27ae60", hbsae = "#e67e22")
meth_col <- function(m) {
  o <- pal[as.character(m)]
  o[is.na(o)] <- "black"
  unname(o)
}

plot_series <- function(d, logxy = TRUE, ylab, main, leg = FALSE, ylim = NULL) {
  ns <- sort(unique(d$n))
  meths <- unique(d$Method)
  xs <- if (logxy) log10(ns) else ns
  ys <- lapply(meths, function(m) d$`Median (s)`[d$Method == m][match(ns, d$n[d$Method == m])])
  if (is.null(ylim)) ylim <- range(unlist(ys), finite = TRUE)
  plot(NA, xlim = range(xs), ylim = ylim,
       xaxt = if (logxy) "n" else "s",
       xlab = if (logxy) L("n") else "n",
       ylab = ylab, main = main, log = if (logxy) "y" else "")
  if (logxy) axis(1, at = log10(ns), labels = ns)
  grid(col = "grey85", lty = 1)
  for (i in seq_along(meths)) {
    lines(xs, ys[[i]], col = meth_col(meths[i]), lwd = 2, pch = 16, type = "b")
    points(xs, ys[[i]], col = meth_col(meths[i]), pch = 16)
  }
  if (leg) legend("topleft", legend = meths, col = meth_col(meths),
                  lwd = 2, pch = 16, bty = "n", cex = 0.85)
}

# ---------------------------------------------------------------------------
# Figur 1: skala waktu FH (area, spasial, spasio-temporal) -- log-log
# ---------------------------------------------------------------------------
pdf(fig("fig-time-fh.pdf"), width = 9, height = 6.5, useDingbats = FALSE)
op <- par(mfrow = c(2, 2), mar = c(4.2, 4.4, 2.6, 1.1), cex.lab = 1.05)
plot_series(read_b("eblup_benchmark.rds"), ylab = L("time"),
            main = L("fh"), leg = TRUE)
plot_series(read_b("seblup_benchmark.rds"), ylab = L("time"),
            main = L("sfh"), leg = TRUE)
plot_series(read_b("steblup_benchmark.rds"), ylab = L("time"),
            main = L("stfh"), leg = TRUE)
plot_series(read_b("hb_benchmark.rds"), ylab = L("time"),
            main = L("hb"), leg = TRUE)
par(op); dev.off()

# ---------------------------------------------------------------------------
# Figur 2: skala waktu model beta -- log-log
# ---------------------------------------------------------------------------
pdf(fig("fig-time-beta.pdf"), width = 9, height = 4.4, useDingbats = FALSE)
op <- par(mfrow = c(1, 3), mar = c(4.2, 4.4, 2.6, 1.1), cex.lab = 1.05)
plot_series(read_b("beta_benchmark.rds"), ylab = L("time"),
            main = L("beta_a"), leg = TRUE)
plot_series(read_b("spatial_beta_benchmark.rds"), ylab = L("time"),
            main = L("beta_b"), leg = TRUE)
plot_series(read_b("spatio_temporal_beta_benchmark.rds"),
            ylab = L("time"),
            main = L("beta_c"), leg = TRUE)
par(op); dev.off()

# ---------------------------------------------------------------------------
# Figur 3: memori vs n -- log-log
# ---------------------------------------------------------------------------
pdf(fig("fig-mem.pdf"), width = 9, height = 6.5, useDingbats = FALSE)
op <- par(mfrow = c(2, 2), mar = c(4.2, 4.4, 2.6, 1.1), cex.lab = 1.05)
mem_series <- function(file, main) {
  d <- read_b(file)
  ns <- sort(unique(d$n)); meths <- unique(d$Method)
  xs <- log10(ns)
  ys <- lapply(meths, function(m) d$`Memory (MB)`[d$Method == m][match(ns, d$n[d$Method == m])])
  plot(NA, xlim = range(xs), ylim = range(unlist(ys), finite = TRUE),
       xaxt = "n", xlab = L("n"),
       ylab = L("mem"), main = main, log = "y")
  axis(1, at = log10(ns), labels = ns)
  grid(col = "grey85", lty = 1)
  for (i in seq_along(meths)) {
    lines(xs, ys[[i]], col = meth_col(meths[i]), lwd = 2, type = "b", pch = 16)
    points(xs, ys[[i]], col = meth_col(meths[i]), pch = 16)
  }
  legend("topleft", legend = meths, col = meth_col(meths), lwd = 2, pch = 16,
         bty = "n", cex = 0.85)
}
mem_series("eblup_benchmark.rds", L("fh"))
mem_series("seblup_benchmark.rds", L("sfh"))
mem_series("steblup_benchmark.rds", L("stfh"))
mem_series("beta_benchmark.rds", L("beta_d"))
par(op); dev.off()

# ---------------------------------------------------------------------------
# Figur 4: rasio kecepatan vs n (pembanding / fastsae)
# ---------------------------------------------------------------------------
pdf(fig("fig-ratio.pdf"), width = 9, height = 4.4, useDingbats = FALSE)
op <- par(mfrow = c(1, 3), mar = c(4.4, 4.6, 2.6, 1.1), cex.lab = 1.05)
ratio_panel <- function(file, main, cols = NULL) {
  d <- read_b(file)
  ns <- sort(unique(d$n))
  meths <- setdiff(unique(d$Method), "fastsae")
  ft <- d$`Median (s)`[d$Method == "fastsae"][match(ns, d$n[d$Method == "fastsae"])]
  rs <- lapply(meths, function(m)
    d$`Median (s)`[d$Method == m][match(ns, d$n[d$Method == m])] / ft)
  plot(NA, xlim = range(log10(ns)), ylim = range(unlist(rs), finite = TRUE),
       xaxt = "n", xlab = L("n"),
       ylab = L("ratio"),
       main = main, log = "y")
  axis(1, at = log10(ns), labels = ns)
  abline(h = 1, lty = 2, col = "grey40")
  grid(col = "grey85", lty = 1)
  for (i in seq_along(meths)) {
    lines(log10(ns), rs[[i]], col = meth_col(meths[i]), lwd = 2, type = "b", pch = 16)
    points(log10(ns), rs[[i]], col = meth_col(meths[i]), pch = 16)
  }
  legend("topleft", legend = meths, col = meth_col(meths), lwd = 2, pch = 16,
         bty = "n", cex = 0.85)
}
ratio_panel("eblup_benchmark.rds", L("fh"))
ratio_panel("seblup_benchmark.rds", L("sfh"))
ratio_panel("beta_benchmark.rds", sub("^\\(a\\)", "(c)", L("beta_a")))
par(op); dev.off()

# ---------------------------------------------------------------------------
# Figur 5: ekuivalensi numerik fastsae vs tipsae (n = 1000)
# ---------------------------------------------------------------------------
d <- read_b("beta_estimates_n1000.rds")
ef <- d[["estimasi fastsae"]]; et <- d[["estimasi tipsae"]]
mf <- d[["mse fastsae"]];      mt <- d[["mse tipsae"]]
pdf(fig("fig-equiv.pdf"), width = 9, height = 4.4, useDingbats = FALSE)
op <- par(mfrow = c(1, 2), mar = c(4.4, 4.6, 2.6, 1.1), cex.lab = 1.05)
plot(et, ef, pch = 16, col = "#1b6ca8aa", cex = 0.7,
     xlab = L("est_t"), ylab = L("est_f"),
     main = L("pred"))
abline(0, 1, col = "red", lwd = 1.5, lty = 2)
plot(mt, mf, pch = 16, col = "#1b6ca8aa", cex = 0.7, log = "xy",
     xlab = L("mse_t"), ylab = L("mse_f"),
     main = L("mse"))
abline(0, 1, col = "red", lwd = 1.5, lty = 2)
par(op); dev.off()

# ---------------------------------------------------------------------------
# Figur 6: engine hb_area() pada model spasial (INLA, JAX CPU/GPU, tipsae)
# sumber: hb_engine_benchmark*.rds; garis putus-putus = mesin kedua dst.
# ---------------------------------------------------------------------------
eng_files <- list.files(root, pattern = "^hb_engine_benchmark.*\\.rds$", full.names = TRUE)
if (length(eng_files) > 0) {
  tm <- do.call(rbind, lapply(eng_files, function(f) readRDS(f)$timings))
  eng_col <- c(inla = "#1D6A5C", cpu = "#E69F00", metal = "#7B5EA7", cuda = "#76B900", tipsae = "#C0392B")
  eng_lab <- c(inla = "INLA", cpu = "JAX CPU", metal = "JAX Metal", cuda = "JAX CUDA", tipsae = "tipsae")
  machines <- unique(tm$machine)
  pdf(fig("fig-engines.pdf"), width = 9, height = 4.4, useDingbats = FALSE)
  op <- par(mfrow = c(1, 2), mar = c(4.4, 4.6, 2.6, 1.1), cex.lab = 1.05)
  for (d in c("gaussian", "beta")) {
    sub <- tm[tm$design == d, ]
    ns <- sort(unique(sub$n))
    plot(NA, xlim = range(log10(ns)), ylim = range(sub$time_warm), log = "y", xaxt = "n",
         xlab = L("D"), ylab = L("eng_time"), main = L(if (d == "gaussian") "eng_g" else "eng_b"))
    axis(1, at = log10(ns), labels = ns)
    grid(col = "grey85", lty = 1)
    for (m in seq_along(machines)) for (e in names(eng_col)) {
      r <- sub[sub$engine == e & sub$machine == machines[m], ]
      if (nrow(r) == 0) next
      r <- r[order(r$n), ]
      lines(log10(r$n), r$time_warm, col = eng_col[[e]], lwd = 2, lty = m, type = "b", pch = 16)
    }
    present <- intersect(names(eng_col), unique(sub$engine))
    legend("topleft", legend = eng_lab[present], col = eng_col[present], lwd = 2, pch = 16,
           bty = "n", cex = 0.85)
    if (length(machines) > 1) legend("bottomright", legend = machines, lty = seq_along(machines),
                                     lwd = 2, bty = "n", cex = 0.8)
  }
  par(op); dev.off()
}

cat(L("done"), normalizePath(out_dir), ":\n")
print(list.files(out_dir, pattern = "\\.pdf$"))
