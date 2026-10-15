# Dump inst/extdata/*.rds to stdout for bench_data_dump.md
# READ-ONLY: only reads files, prints to stdout.
options(digits = 17, width = 200, scipen = 999)

files <- sort(list.files("/Volumes/work/_MainR/fastsae-cov/inst/extdata", pattern = "\\.rds$", full.names = TRUE))

fmt_num <- function(x) {
  if (is.numeric(x)) {
    out <- format(x, digits = 17, trim = TRUE, scientific = FALSE)
    # fall back for extremely large/small
    out[!is.finite(x)] <- as.character(x[!is.finite(x)])
    out
  } else {
    as.character(x)
  }
}

dump_df <- function(df, name) {
  cat("### Object:", name, "\n\n")
  cat("- class:", paste(class(df), collapse = ", "), "\n")
  if (is.data.frame(df)) {
    cat("- jumlah baris:", nrow(df), "\n")
    cat("- jumlah kolom:", ncol(df), "\n")
    cls <- vapply(df, function(c) paste(class(c), collapse = "/"), character(1))
    cat("- kolom + tipe:\n")
    for (i in seq_along(df)) cat("  ", i, ".", names(df)[i], ":", cls[i], "\n")
    cat("\n")
    # markdown table
    hdr <- paste0("| ", paste(names(df), collapse = " | "), " |")
    sep <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
    rows <- apply(df, 1, function(r) paste0("| ", paste(vapply(seq_along(r), function(j) fmt_num(r[[j]]), character(1)), collapse = " | "), " |"))
    cat(hdr, "\n", sep, "\n", sep = "\n")
    cat(rows, sep = "\n")
    cat("\n\n")
  } else {
    cat("- (bukan data.frame)\n\n")
    print(df)
    cat("\n")
  }
}

walk <- function(obj, path) {
  if (is.data.frame(obj)) {
    dump_df(obj, path)
  } else if (is.list(obj)) {
    cat("### Object:", path, "\n\n")
    cat("- class:", paste(class(obj), collapse = ", "), "\n")
    cat("- elemen list:", length(obj), "\n")
    if (!is.null(names(obj))) cat("- nama elemen:", paste(names(obj), collapse = ", "), "\n")
    if (!is.null(dim(obj))) cat("- dim:", paste(dim(obj), collapse = " x "), "\n")
    if (length(obj) == 0) { cat("(kosong)\n\n"); return(invisible()) }
    # atomic vector leaf?
    if (!is.null(names(obj)) && all(vapply(obj, function(x) is.atomic(x) && !is.list(x), logical(1)))) {
      cat("\n")
      for (nm in names(obj)) cat("  - ", nm, ": ", fmt_num(obj[[nm]]), "\n", sep = "")
      cat("\n")
      return(invisible())
    }
    cat("\n")
    if (is.null(names(obj))) {
      for (i in seq_along(obj)) walk(obj[[i]], paste0(path, "[[", i, "]]"))
    } else {
      for (nm in names(obj)) {
        if (is.data.frame(obj[[nm]]) || is.list(obj[[nm]])) walk(obj[[nm]], paste0(path, "$", nm))
        else {
          cat("- ", nm, " (", paste(class(obj[[nm]]), collapse = "/"), "): ", sep = "")
          v <- obj[[nm]]
          if (is.atomic(v) && length(v) <= 50) cat(paste(fmt_num(v), collapse = ", ")) else cat("len=", length(v))
          cat("\n")
        }
      }
      cat("\n")
    }
  } else if (is.matrix(obj) || is.array(obj)) {
    cat("### Object:", path, "\n\n")
    cat("- class:", paste(class(obj), collapse = ", "), "\n")
    cat("- dim:", paste(dim(obj), collapse = " x "), "\n\n")
    print(obj, digits = 17)
    cat("\n")
  } else if (is.atomic(obj)) {
    cat("### Object:", path, "\n\n")
    cat("- class:", paste(class(obj), collapse = ", "), "\n")
    cat("- length:", length(obj), "\n\n")
    print(obj, digits = 17)
    cat("\n")
  } else {
    cat("### Object:", path, "\n\n- class:", paste(class(obj), collapse = ", "), "(str)\n\n")
    str(obj)
    cat("\n")
  }
}

for (f in files) {
  cat("\n\n======================================================================\n")
  cat("## FILE:", basename(f), "\n")
  cat("path:", normalizePath(f), "\n")
  cat("size bytes:", file.info(f)$size, "\n")
  obj <- readRDS(f)
  cat("class:", paste(class(obj), collapse = ", "), "\n")
  if (!is.null(dim(obj))) cat("dim:", paste(dim(obj), collapse = " x "), "\n")
  cat("======================================================================\n\n")
  walk(obj, basename(f))
}
