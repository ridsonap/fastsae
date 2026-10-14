#!/usr/bin/env Rscript

# ==============================================================================
# Benchmark: hb_area() engines
#
#   inla   inla_engine()                     INLA (Laplace approximation)
#   cpu    jax_engine(device = "cpu")        NumPyro NUTS on the CPU (JAX/XLA)
#   metal  jax_engine(device = "metal")      NumPyro NUTS on the Apple GPU (jax-metal)
#
# JAX fixes its device when it is first used in a process, so every engine is
# run in its own Rscript child process. Each fit is timed `--reps` times on the
# same data: the first time includes JAX's compilation of the sampler, the
# fastest later one is the warm time. JAX settings are the jax_engine()
# defaults (2 chains x 1000 draws, 500 warmup).
#
# Designs (data from sim_area_data(D = n, n_unsampled = 0, seed = 42 + n)):
#   gaussian  y_gaussian ~ x1 + x2, Fay-Herriot (family = "gaussian")
#   beta      y_beta ~ x1 + x2, Beta model for proportions (family = "beta")
#
# Writes inst/extdata/hb_engine_benchmark.rds with
#   $timings    design, engine, n, time_first, time_warm, device, rhat_max
#   $agreement  design, engine, n, cor_hb, max_abs_diff_hb, median_ratio_sd
#               (JAX engines against INLA on the same data)
#   $meta       date, machine, versions and settings
#
# Usage: Rscript benchmarks/bench_hb_area_engines.R [--n=100,1000,5000,10000]
#          [--designs=gaussian,beta] [--engines=inla,cpu,metal] [--reps=2]
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit[1])
}
split_arg <- function(x) strsplit(x, ",", fixed = TRUE)[[1]]

fit_design <- function(design, engine, df) {
  eng <- switch(engine,
    inla = fastsae::inla_engine(),
    cpu = fastsae::jax_engine(device = "cpu", progress = FALSE),
    metal = fastsae::jax_engine(device = "metal", progress = FALSE)
  )
  switch(design,
    gaussian = fastsae::hb_area(y_gaussian ~ x1 + x2, vardir = "vardir", data = df,
                                family = "gaussian", engine = eng, print_result = FALSE),
    beta = fastsae::hb_area(y_beta ~ x1 + x2, vardir = "vardir", data = df,
                            family = "beta", engine = eng, print_result = FALSE)
  )
}

# ---- Child: one engine, one design, all n --------------------------------------
if ("--child" %in% args) {
  suppressPackageStartupMessages(library(fastsae))
  engine <- arg_value("engine", "inla")
  design <- arg_value("design", "gaussian")
  ns <- as.integer(split_arg(arg_value("n", "100")))
  reps <- as.integer(arg_value("reps", "2"))
  out <- lapply(ns, function(n) {
    df <- sim_area_data(D = n, n_unsampled = 0, seed = 42 + n)$data
    times <- numeric(reps)
    for (r in seq_len(reps)) {
      times[r] <- system.time(fit <- fit_design(design, engine, df))[["elapsed"]]
    }
    cat(sprintf("%-8s %-6s n = %6d  first %8.2fs  warm %8.2fs  %s\n", design, engine, n,
                times[1], min(times[-1], times[1]), fit$device %||% "INLA"))
    list(
      timing = data.frame(
        design = design, engine = engine, n = n,
        time_first = times[1],
        time_warm = if (reps > 1) min(times[-1]) else times[1],
        device = fit$device %||% "cpu (INLA)",
        rhat_max = fit$rhat_max %||% NA_real_
      ),
      est = data.frame(hb = fit$df_hb$hb, sd = fit$df_hb$sd)
    )
  })
  names(out) <- ns
  if (engine != "inla") {
    md <- reticulate::import("importlib.metadata")
    pkgs <- c("jax", "jaxlib", "jax-metal", "numpyro")
    attr(out, "python") <- vapply(pkgs, function(p) {
      tryCatch(md$version(p), error = function(e) NA_character_)
    }, character(1))
  }
  saveRDS(out, arg_value("out", tempfile(fileext = ".rds")))
  quit(save = "no")
}

# ---- Driver -------------------------------------------------------------------
ns <- arg_value("n", "100,1000,5000,10000")
designs <- split_arg(arg_value("designs", "gaussian,beta"))
engines <- split_arg(arg_value("engines", "inla,cpu,metal"))
reps <- arg_value("reps", "2")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])

runs <- list()
for (design in designs) {
  for (engine in engines) {
    out <- tempfile(fileext = ".rds")
    status <- system2(file.path(R.home("bin"), "Rscript"),
      c(shQuote(script), "--child", paste0("--engine=", engine), paste0("--design=", design),
        paste0("--n=", ns), paste0("--reps=", reps), paste0("--out=", out)))
    if (status != 0 || !file.exists(out)) stop("Benchmark failed: ", design, " / ", engine)
    runs[[paste(design, engine)]] <- readRDS(out)
  }
}

timings <- do.call(rbind, lapply(runs, function(r) do.call(rbind, lapply(r, `[[`, "timing"))))
rownames(timings) <- NULL

agreement <- do.call(rbind, lapply(designs, function(design) {
  ref <- runs[[paste(design, "inla")]]
  do.call(rbind, lapply(setdiff(engines, "inla"), function(engine) {
    run <- runs[[paste(design, engine)]]
    do.call(rbind, lapply(names(run), function(n) {
      a <- run[[n]]$est
      b <- ref[[n]]$est
      data.frame(design = design, engine = engine, n = as.integer(n),
                 cor_hb = stats::cor(a$hb, b$hb),
                 max_abs_diff_hb = max(abs(a$hb - b$hb)),
                 median_ratio_sd = stats::median(a$sd / b$sd))
    }))
  }))
}))

py_versions <- Filter(Negate(is.null), lapply(runs, attr, "python"))
py_versions <- if (length(py_versions) > 0) py_versions[[length(py_versions)]] else NULL

meta <- list(
  date = "2026-10-14",
  machine = paste(Sys.info()[["sysname"]], Sys.info()[["machine"]],
                  tryCatch(system("sysctl -n machdep.cpu.brand_string", intern = TRUE),
                           error = function(e) "")),
  r_version = R.version.string,
  fastsae = as.character(utils::packageVersion("fastsae")),
  inla = tryCatch(as.character(utils::packageVersion("INLA")), error = function(e) NA_character_),
  python = py_versions,
  jax_settings = "jax_engine() defaults: 2 chains x 1000 draws, 500 warmup",
  reps = as.integer(reps)
)

res <- list(timings = timings, agreement = agreement, meta = meta)
dest <- file.path(dirname(script), "..", "inst", "extdata", "hb_engine_benchmark.rds")
saveRDS(res, dest)
print(timings)
print(agreement)
cat("Saved", normalizePath(dest), "\n")
