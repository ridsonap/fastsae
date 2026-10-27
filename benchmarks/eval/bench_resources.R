#!/usr/bin/env Rscript

# ==============================================================================
# Speed and memory of the evaluation article, one cell (model, method, D) per R
# process, run one after another on an otherwise idle machine.
#
#   sfh   fastsae / sae / emdi, REML + analytical MSE, D = 100, 500, 1000, 2500
#   beta  inla / jax (NumPyro, CPU) / metal (NumPyro, Apple GPU) / tipsae,
#         D = 50, 100, 500, 1000 (base scenario of sim_accuracy.R on a grid)
#
# Time: each cell fits the same data `reps` times in one process; the first time
# includes start-up costs (JAX compilation, Stan's first call), the median of the
# later ones is the warm time. For MCMC, the smallest bulk/tail ESS of each fit is
# kept so that time per effective draw can be reported.
#
# Memory: the parent polls the resident set size of the child process and all its
# descendants (INLA runs an external program) every 0.1 s with the ps package and
# keeps the peak; the peak of a child that loads the same packages and data
# without fitting ("baseline") is subtracted. Starting Python/JAX is part of the
# NumPyro cells. bench::mark()'s mem_alloc is not used
# because it counts only R's own allocations. Peaks shorter than the polling
# interval can be missed, and on Apple silicon GPU (Metal) memory is unified memory
# that RSS does not fully show.
#
# Cells are cached in benchmarks/output/eval/resources/, so a run resumes.
# Writes inst/extdata/eval_resources.rds.
#
# Usage: Rscript benchmarks/eval/bench_resources.R [--models=sfh,beta] [--reps=3]
#          [--D=50,100] [--methods=inla,tipsae]   (subsets, e.g. for a quick check)
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit[1])
}

# ---- child: fit one cell ---------------------------------------------------------
if ("--child" %in% args) {
  suppressPackageStartupMessages(library(fastsae))
  source("benchmarks/eval/eval_common.R")
  model <- arg_value("model"); method <- arg_value("method")
  D <- as.integer(arg_value("D")); reps <- as.integer(arg_value("reps"))
  if (model == "sfh") {
    sc <- list(map = "grid", D = D, rho = 0.5, sigma2 = 1, ratio = 1)
    des <- sfh_design(sc, 1000)
    dat <- sfh_replicate(des, sc, 1)
    # every cell (and the baseline) loads the same packages, so the peak difference is the fit
    loadNamespace("sae"); loadNamespace("emdi")
  } else {
    sc <- list(map = "grid", D = D, sigma_s = 0.5, sigma_v = 0.2, kappa = 25)
    des <- beta_design(sc, 1000)
    dat <- beta_replicate(des, sc, 1)
    loadNamespace("tipsae"); loadNamespace("INLA")
  }
  times <- ess <- rep(NA_real_, reps)
  if (method != "baseline") {
    for (r in seq_len(reps)) {
      times[r] <- system.time(fit <- suppressWarnings(suppressMessages(
        if (model == "sfh") fit_sfh(method, dat, des$W_std) else fit_beta(method, dat, des, seed = r)
      )))[["elapsed"]]
      if (model == "beta") ess[r] <- fit$ess_min
    }
  }
  saveRDS(list(times = times, ess = ess), arg_value("out"))
  quit(save = "no")
}

# ---- parent ---------------------------------------------------------------------
models <- strsplit(arg_value("models", "sfh,beta"), ",", fixed = TRUE)[[1]]
reps_default <- as.integer(arg_value("reps", "3"))
cell_dir <- file.path(Sys.getenv("EVAL_OUT", "benchmarks/output/eval"), "resources")
dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)
grid <- list(
  sfh = list(methods = c("fastsae", "sae", "emdi"), D = c(100, 500, 1000, 2500)),
  beta = list(methods = c("inla", "jax", "metal", "tipsae"), D = c(50, 100, 500, 1000))
)

peak_tree_rss <- function(proc) {
  peak <- 0
  while (proc$is_alive()) {
    h <- tryCatch(proc$as_ps_handle(), error = function(e) NULL)
    if (!is.null(h)) {
      rss <- tryCatch({
        kids <- ps::ps_children(h, recursive = TRUE)
        sum(vapply(c(list(h), kids), function(k)
          tryCatch(ps::ps_memory_info(k)[["rss"]], error = function(e) 0), 0))
      }, error = function(e) 0)
      peak <- max(peak, rss)
    }
    Sys.sleep(0.1)
  }
  peak / 1024^2
}

run_cell <- function(model, method, D, reps) {
  out <- file.path(cell_dir, sprintf("%s_%s_%d.rds", model, method, D))
  if (file.exists(out)) return(readRDS(out))
  tmp <- tempfile(fileext = ".rds")
  log <- tempfile(fileext = ".log")
  proc <- processx::process$new(
    "Rscript", c("benchmarks/eval/bench_resources.R", "--child", paste0("--model=", model),
                 paste0("--method=", method), paste0("--D=", D), paste0("--reps=", reps), paste0("--out=", tmp)),
    stdout = log, stderr = "2>&1"
  )
  peak <- peak_tree_rss(proc)
  if (!file.exists(tmp)) stop("cell ", model, "/", method, "/", D, " failed:\n", paste(readLines(log), collapse = "\n"))
  res <- readRDS(tmp)
  res$peak_mb <- peak
  saveRDS(res, out)
  res
}

`%||%` <- function(a, b) if (length(a) == 0) b else a
only_D <- arg_value("D"); only_methods <- arg_value("methods")
rows <- list()
for (model in models) {
  Ds <- if (is.null(only_D)) grid[[model]]$D else as.numeric(strsplit(only_D, ",", fixed = TRUE)[[1]])
  ms <- if (is.null(only_methods)) grid[[model]]$methods else intersect(strsplit(only_methods, ",", fixed = TRUE)[[1]], grid[[model]]$methods)
  for (D in Ds) {
    base <- run_cell(model, "baseline", D, 1L)
    for (m in ms) {
      # the slowest cells (several minutes per fit) are fitted twice
      reps <- if ((model == "sfh" && m != "fastsae" && D >= 2500) || (model == "beta" && m != "inla" && D >= 1000)) 2L else reps_default
      cat(sprintf(">>> %-4s %-8s D = %4d (%d fits) ... ", model, m, D, reps))
      r <- tryCatch(run_cell(model, m, D, reps), error = function(e) { cat("FAILED\n", conditionMessage(e), "\n"); NULL })
      if (is.null(r)) next
      warm <- if (reps > 1) r$times[-1] else r$times
      cat(sprintf("first %.2f s, warm %.2f s, peak +%.0f MB\n", r$times[1], stats::median(warm), r$peak_mb - base$peak_mb))
      rows[[length(rows) + 1]] <- data.frame(
        model = model, method = m, D = D, reps = reps,
        time_first = r$times[1], time_warm = stats::median(warm),
        time_warm_min = min(warm), time_warm_max = max(warm),
        ess_min = if (all(is.na(r$ess))) NA_real_ else stats::median(r$ess[-1] %||% r$ess, na.rm = TRUE),
        peak_mb = r$peak_mb, baseline_mb = base$peak_mb, peak_increase_mb = max(r$peak_mb - base$peak_mb, 0)
      )
    }
  }
}

py <- tryCatch({
  reticulate::use_virtualenv("r-numpyro-gpu", required = FALSE)
  md <- reticulate::import("importlib.metadata")
  vapply(c("jax", "jaxlib", "jax-metal", "numpyro"), function(p) tryCatch(md$version(p), error = function(e) NA_character_), "")
}, error = function(e) NA_character_)
meta <- list(
  date = Sys.Date(), r = R.version.string,
  cpu = tryCatch(system("sysctl -n machdep.cpu.brand_string", intern = TRUE), error = function(e) NA_character_),
  cores = parallel::detectCores(),
  ram_gb = tryCatch(as.numeric(system("sysctl -n hw.memsize", intern = TRUE)) / 1024^3, error = function(e) NA_real_),
  os = utils::sessionInfo()$running, blas = utils::sessionInfo()$BLAS,
  omp_threads = Sys.getenv("OMP_NUM_THREADS", "unset"),
  versions = vapply(c("fastsae", "sae", "emdi", "tipsae", "INLA", "rstan"), function(p)
    tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_), ""),
  python = py,
  settings = c(jax = "jax_engine() defaults: 2 chains x 1000 draws, 500 warmup",
               tipsae = "fit_sae() defaults: 4 chains x 2000 iterations (1000 warmup), chains run one after another",
               inla = "inla_engine() defaults; hb_area() runs INLA single-threaded",
               sfh = "REML, analytical MSE, fastsae n_threads = 1")
)
out <- list(cells = do.call(rbind, rows), meta = meta)
saveRDS(out, file.path(Sys.getenv("EVAL_EXTDATA", "inst/extdata"), "eval_resources.rds"))
print(out$cells, digits = 3)
