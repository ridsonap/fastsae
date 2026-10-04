# scripts/02_figures.R
# Menghasilkan seluruh figur benchmark laporan teknis fastsae
# dari file sumber di inst/extdata/*.rds.
#
# Jalankan dari folder paper/:   Rscript scripts/02_figures.R
# Output: figures/*.pdf
#
# ATURAN: tidak ada angka yang diketik manual; sumbu dan kurva dibangun
# langsung dari .rds. Skala log pada kedua sumbu bila diminta.

root <- if (dir.exists("../inst/extdata")) "../inst/extdata" else "inst/extdata"
dir.create("figures", showWarnings = FALSE)

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
       xlab = if (logxy) expression("jumlah area "*n) else "n",
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
pdf("figures/fig-time-fh.pdf", width = 9, height = 6.5, useDingbats = FALSE)
op <- par(mfrow = c(2, 2), mar = c(4.2, 4.4, 2.6, 1.1), cex.lab = 1.05)
plot_series(read_b("eblup_benchmark.rds"), ylab = "waktu median (detik, log)",
            main = "(a) Fay-Herriot area-level", leg = TRUE)
plot_series(read_b("seblup_benchmark.rds"), ylab = "waktu median (detik, log)",
            main = "(b) Fay-Herriot spasial", leg = TRUE)
plot_series(read_b("steblup_benchmark.rds"), ylab = "waktu median (detik, log)",
            main = "(c) Fay-Herriot spasio-temporal", leg = TRUE)
plot_series(read_b("hb_benchmark.rds"), ylab = "waktu median (detik, log)",
            main = "(d) Hierarki Bayes Gaussian", leg = TRUE)
par(op); dev.off()

# ---------------------------------------------------------------------------
# Figur 2: skala waktu model beta -- log-log
# ---------------------------------------------------------------------------
pdf("figures/fig-time-beta.pdf", width = 9, height = 4.4, useDingbats = FALSE)
op <- par(mfrow = c(1, 3), mar = c(4.2, 4.4, 2.6, 1.1), cex.lab = 1.05)
plot_series(read_b("beta_benchmark.rds"), ylab = "waktu median (detik, log)",
            main = "(a) Beta area-level", leg = TRUE)
plot_series(read_b("spatial_beta_benchmark.rds"), ylab = "waktu median (detik, log)",
            main = "(b) Beta spasial", leg = TRUE)
plot_series(read_b("spatio_temporal_beta_benchmark.rds"),
            ylab = "waktu median (detik, log)",
            main = "(c) Beta spasio-temporal", leg = TRUE)
par(op); dev.off()

# ---------------------------------------------------------------------------
# Figur 3: memori vs n -- log-log
# ---------------------------------------------------------------------------
pdf("figures/fig-mem.pdf", width = 9, height = 6.5, useDingbats = FALSE)
op <- par(mfrow = c(2, 2), mar = c(4.2, 4.4, 2.6, 1.1), cex.lab = 1.05)
mem_series <- function(file, main) {
  d <- read_b(file)
  ns <- sort(unique(d$n)); meths <- unique(d$Method)
  xs <- log10(ns)
  ys <- lapply(meths, function(m) d$`Memory (MB)`[d$Method == m][match(ns, d$n[d$Method == m])])
  plot(NA, xlim = range(xs), ylim = range(unlist(ys), finite = TRUE),
       xaxt = "n", xlab = expression("jumlah area "*n),
       ylab = expression("memori (MB, log)"), main = main, log = "y")
  axis(1, at = log10(ns), labels = ns)
  grid(col = "grey85", lty = 1)
  for (i in seq_along(meths)) {
    lines(xs, ys[[i]], col = meth_col(meths[i]), lwd = 2, type = "b", pch = 16)
    points(xs, ys[[i]], col = meth_col(meths[i]), pch = 16)
  }
  legend("topleft", legend = meths, col = meth_col(meths), lwd = 2, pch = 16,
         bty = "n", cex = 0.85)
}
mem_series("eblup_benchmark.rds", "(a) Fay-Herriot area-level")
mem_series("seblup_benchmark.rds", "(b) Fay-Herriot spasial")
mem_series("steblup_benchmark.rds", "(c) Fay-Herriot spasio-temporal")
mem_series("beta_benchmark.rds", "(d) Beta area-level")
par(op); dev.off()

# ---------------------------------------------------------------------------
# Figur 4: rasio kecepatan vs n (pembanding / fastsae)
# ---------------------------------------------------------------------------
pdf("figures/fig-ratio.pdf", width = 9, height = 4.4, useDingbats = FALSE)
op <- par(mfrow = c(1, 3), mar = c(4.4, 4.6, 2.6, 1.1), cex.lab = 1.05)
ratio_panel <- function(file, main, cols = NULL) {
  d <- read_b(file)
  ns <- sort(unique(d$n))
  meths <- setdiff(unique(d$Method), "fastsae")
  ft <- d$`Median (s)`[d$Method == "fastsae"][match(ns, d$n[d$Method == "fastsae"])]
  rs <- lapply(meths, function(m)
    d$`Median (s)`[d$Method == m][match(ns, d$n[d$Method == m])] / ft)
  plot(NA, xlim = range(log10(ns)), ylim = range(unlist(rs), finite = TRUE),
       xaxt = "n", xlab = expression("jumlah area "*n),
       ylab = expression("rasio waktu (pembanding / fastsae)"),
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
ratio_panel("eblup_benchmark.rds", "(a) Fay-Herriot area-level")
ratio_panel("seblup_benchmark.rds", "(b) Fay-Herriot spasial")
ratio_panel("beta_benchmark.rds", "(c) Beta area-level")
par(op); dev.off()

# ---------------------------------------------------------------------------
# Figur 5: ekuivalensi numerik fastsae vs tipsae (n = 1000)
# ---------------------------------------------------------------------------
d <- read_b("beta_estimates_n1000.rds")
ef <- d[["estimasi fastsae"]]; et <- d[["estimasi tipsae"]]
mf <- d[["mse fastsae"]];      mt <- d[["mse tipsae"]]
pdf("figures/fig-equiv.pdf", width = 9, height = 4.4, useDingbats = FALSE)
op <- par(mfrow = c(1, 2), mar = c(4.4, 4.6, 2.6, 1.1), cex.lab = 1.05)
plot(et, ef, pch = 16, col = "#1b6ca8aa", cex = 0.7,
     xlab = "estimasi tipsae", ylab = "estimasi fastsae",
     main = "(a) Prediktor area")
abline(0, 1, col = "red", lwd = 1.5, lty = 2)
plot(mt, mf, pch = 16, col = "#1b6ca8aa", cex = 0.7, log = "xy",
     xlab = "MSE tipsae (log)", ylab = "MSE fastsae (log)",
     main = "(b) MSE")
abline(0, 1, col = "red", lwd = 1.5, lty = 2)
par(op); dev.off()

cat("Selesai. Berkas figur di", normalizePath("figures"), ":\n")
print(list.files("figures", pattern = "\\.pdf$"))
