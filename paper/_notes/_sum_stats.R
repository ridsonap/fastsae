# cross-check: medians, per-n ratios, extremes
d <- "/Volumes/work/_MainR/fastsae-cov/inst/extdata"
f <- function(x) format(x, digits = 12)
cat("--- seblup: median of Median(s) and median of Memory (MB) ---\n")
x <- readRDS(file.path(d, "seblup_benchmark.rds"))
print(f(tapply(x[["Median (s)"]], x$Method, median)))
print(f(tapply(x[["Memory (MB)"]], x$Method, median)))
cat("--- eblup: median of Median(s) ---\n")
y <- readRDS(file.path(d, "eblup_benchmark.rds"))
print(f(tapply(y[["Median (s)"]], y$Method, median)))
print(f(tapply(y[["Memory (MB)"]], y$Method, median)))
cat("--- steblup: median ---\n")
z <- readRDS(file.path(d, "steblup_benchmark.rds"))
print(f(tapply(z[["Median (s)"]], z$Method, median)))
print(f(tapply(z[["Memory (MB)"]], z$Method, median)))
cat("\n--- per-n ratio (comparator / fastsae), all files ---\n")
for (fn in list.files(d, pattern = "rds$")) {
  x <- readRDS(file.path(d, fn))
  if (!is.data.frame(x) || !all(c("Method", "Median (s)", "n") %in% names(x))) next
  if (!"fastsae" %in% x$Method) next
  cat(fn, "\n")
  for (nn in sort(unique(x$n))) {
    s <- x[x$n == nn, ]
    v <- s[["Median (s)"]]; names(v) <- s$Method
    rs <- setdiff(names(v), "fastsae")
    cat("  n=", nn, ": ", paste(sprintf("%s=%.4g", rs, v[rs] / v[["fastsae"]]), collapse = "  "), "\n", sep = "")
  }
}
cat("\n--- fastsae Memory (MB) above 10 and above 25, all files ---\n")
for (fn in list.files(d, pattern = "rds$")) {
  x <- readRDS(file.path(d, fn))
  if (!is.data.frame(x) || !"Memory (MB)" %in% names(x)) next
  s <- x[x$Method == "fastsae" & x[["Memory (MB)"]] > 10, ]
  if (nrow(s)) cat(fn, ":", paste(sprintf("n=%d:%.6gMB", s$n, s[["Memory (MB)"]]), collapse = ", "), "\n")
}
cat("\n--- n_iter column values per file ---\n")
for (fn in list.files(d, pattern = "rds$")) {
  x <- readRDS(file.path(d, fn))
  if (!is.data.frame(x) || !"n_iter" %in% names(x)) next
  cat(fn, ":", paste(sort(unique(x$n_iter)), collapse = ", "), "\n")
}
