# Ringkasan struktur benchmarks/output/*.rds (read-only)
options(digits = 17, width = 220, scipen = 999)
files <- sort(list.files("/Volumes/work/_MainR/fastsae-cov/benchmarks/output", pattern = "\\.rds$", full.names = TRUE))
for (f in files) {
  cat("\n==============================\nFILE:", basename(f), " size:", file.info(f)$size, "bytes\n")
  o <- readRDS(f)
  cat("class:", paste(class(o), collapse = ", "), "\n")
  if (is.list(o) && !is.data.frame(o)) cat("nama elemen:", paste(names(o), collapse = ", "), "\n")
  walk <- function(x, path, depth = 0) {
    if (depth > 3) return(invisible())
    if (is.data.frame(x)) {
      cat("--", path, ": data.frame", nrow(x), "x", ncol(x), "\n")
      cls <- vapply(x, function(c) paste(class(c), collapse = "/"), character(1))
      cat("   kolom:", paste(paste0(names(x), " [", cls, "]"), collapse = "; "), "\n")
      # print head/tail small tables fully if small
      if (nrow(x) <= 12) print(x, digits = 10) else { cat("   head:\n"); print(utils::head(x, 4), digits = 10); cat("   tail:\n"); print(utils::tail(x, 3), digits = 10) }
    } else if (is.list(x)) {
      cat("--", path, ": list (", length(x), ") nama:", paste(names(x), collapse = ", "), "\n")
      for (nm in names(x)) if (!is.null(nm) && nzchar(nm)) walk(x[[nm]], paste0(path, "$", nm), depth + 1)
      if (is.null(names(x))) for (i in seq_along(x)) walk(x[[i]], paste0(path, "[[", i, "]]"), depth + 1)
    } else if (is.atomic(x)) {
      cat("--", path, ":", paste(class(x), collapse = "/"), "len", length(x), "\n")
      if (length(x) <= 40) cat("   nilai:", paste(format(x, digits = 10), collapse = ", "), "\n")
      else cat("   head:", paste(format(utils::head(x, 5), digits = 10), collapse = ", "), "\n")
    }
  }
  if (is.data.frame(o)) walk(o, basename(o)) else for (nm in names(o)) walk(o[[nm]], paste0(basename(f), "$", nm))
}
