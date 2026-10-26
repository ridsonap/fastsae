#!/usr/bin/env Rscript

# ==============================================================================
# Benchmark: frequentist EBLUP + analytical MSE (FH, SFH, STFH) of fastsae, sae and
# emdi on the same simulated data sets. sae has no analytical MSE for the
# spatio-temporal model, so STFH compares point estimates (fastsae vs
# sae::eblupSTFH()); the row "fastsae_mse" adds the analytical MSE of fastsae.
#
# Each (model, method, D) cell runs in its own R process (bench_frequentist_child.R)
# under /usr/bin/time, which gives
#   - the median time of bench::mark(),
#   - the memory allocated by R during one call (bench's mem_alloc; it does not see
#     memory allocated by compiled code), and
#   - the peak resident set size of the process minus that of a process that loads
#     the same packages and data without fitting ("Peak RSS"), which also counts
#     memory allocated by compiled code.
# Cells are cached in benchmarks/output/freq_cells/, so an interrupted run resumes.
# Writes inst/extdata/{eblup,seblup,steblup}_benchmark.rds and
# inst/extdata/frequentist_benchmark_meta.rds.
#
# Usage (from the package root, on an otherwise idle machine):
#   Rscript benchmarks/bench_frequentist.R [D ...]
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
n_args <- suppressWarnings(as.integer(args))
D_values <- if (length(n_args[!is.na(n_args)]) > 0) n_args[!is.na(n_args)] else c(30, 50, 100, 250, 500, 1000)

child <- "benchmarks/bench_frequentist_child.R"
cell_dir <- "benchmarks/output/freq_cells"
dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)

methods <- list(FH = c("fastsae", "sae", "emdi"), SFH = c("fastsae", "sae", "emdi"), STFH = c("fastsae", "fastsae_mse", "sae"))
files <- c(FH = "eblup_benchmark.rds", SFH = "seblup_benchmark.rds", STFH = "steblup_benchmark.rds")

# Iterations from the expected time of one call (earlier runs), so that slow cells
# are not repeated needlessly
expected <- function(model, method, D) {
  old <- tryCatch(as.data.frame(readRDS(file.path("inst/extdata", files[[model]]))), error = function(e) NULL)
  if (method == "fastsae_mse") method <- "fastsae"
  v <- if (is.null(old)) numeric() else old$`Median (s)`[old$Method == method & old$n == D]
  if (length(v)) v[1] else 1
}
n_iter_for <- function(t) if (t < 0.5) 20L else if (t < 5) 10L else if (t < 60) 5L else if (t < 600) 2L else 1L

run_cell <- function(model, method, D, iter) {
  out <- file.path(cell_dir, sprintf("%s_%s_%d.rds", model, method, D))
  if (file.exists(out)) return(readRDS(out))
  tmp <- tempfile(fileext = ".rds")
  log <- system2("/usr/bin/time", c("-l", "Rscript", child, model, method, D, iter, tmp),
                 stdout = TRUE, stderr = TRUE)
  if (!file.exists(tmp)) stop("cell ", model, "/", method, "/", D, " failed:\n", paste(log, collapse = "\n"))
  rss <- grep("maximum resident set size", log, value = TRUE)
  res <- readRDS(tmp)
  res$peak_rss_mb <- as.numeric(sub("^\\s*([0-9]+).*$", "\\1", rss[1])) / 1024^2
  saveRDS(res, out)
  res
}

rows <- list()
for (model in names(methods)) {
  for (D in D_values) {
    base <- run_cell(model, "baseline", D, 1L)
    for (m in methods[[model]]) {
      iter <- n_iter_for(expected(model, m, D))
      cat(sprintf(">>> %-4s D = %4d %-8s (%d iterations) ... ", model, D, m, iter))
      r <- run_cell(model, m, D, iter)
      cat(sprintf("%.4g s, peak RSS +%.1f MB, converged: %s\n", r$median,
                  r$peak_rss_mb - base$peak_rss_mb, r$converged))
      rows[[length(rows) + 1]] <- data.frame(
        model = model, Method = m, `Min (s)` = r$min, `Median (s)` = r$median,
        `Iterations/sec` = r$itr_sec, `Memory (MB)` = r$mem_alloc_mb,
        `Peak RSS (MB)` = max(r$peak_rss_mb - base$peak_rss_mb, 0),
        n_iter = r$n_iter, converged = r$converged, n = D, check.names = FALSE
      )
    }
  }
}
all_rows <- do.call(rbind, rows)

for (model in names(files)) {
  out <- all_rows[all_rows$model == model, setdiff(names(all_rows), "model")]
  out <- out[order(out$n, match(out$Method, c("fastsae", "fastsae_mse", "emdi", "sae"))), ]
  out <- tibble::as_tibble(out)
  class(out) <- c("bench_mark", "tbl_df", "tbl", "data.frame")
  saveRDS(out, file.path("inst/extdata", files[[model]]))
  cat("\n==", files[[model]], "\n")
  print(as.data.frame(out)[, c("Method", "Median (s)", "Memory (MB)", "Peak RSS (MB)", "converged", "n_iter", "n")], digits = 4)
}

meta <- list(
  date = Sys.Date(),
  r = R.version.string,
  versions = vapply(c("fastsae", "sae", "emdi", "bench"), function(p) as.character(packageVersion(p)), ""),
  cpu = tryCatch(system("sysctl -n machdep.cpu.brand_string", intern = TRUE), error = function(e) NA_character_),
  cores = parallel::detectCores(),
  blas = sessionInfo()$BLAS,
  omp_threads = Sys.getenv("OMP_NUM_THREADS", "unset")
)
saveRDS(meta, "inst/extdata/frequentist_benchmark_meta.rds")
