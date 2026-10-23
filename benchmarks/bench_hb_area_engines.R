#!/usr/bin/env Rscript

# ==============================================================================
# Benchmark: hb_area() engines and tipsae on spatial (ICAR) area-level models
#
#   inla    inla_engine()                     INLA (Laplace approximation)
#   cpu     jax_engine(device = "cpu")        NumPyro NUTS on the CPU (JAX/XLA)
#   metal   jax_engine(device = "metal")      NumPyro NUTS on the Apple GPU (jax-metal)
#   cuda    jax_engine(device = "cuda")       NumPyro NUTS on an NVIDIA GPU
#   tipsae  tipsae::fit_sae()                 Stan NUTS (Beta model only)
#
# JAX fixes its device when it is first used in a process, so every engine is
# run in its own Rscript child process. Each fit is timed `--reps` times on the
# same data: the first time includes compilation (JAX), the fastest later one
# is the warm time. JAX settings are the jax_engine() defaults (2 chains x 1000
# draws, 500 warmup); tipsae uses its defaults.
#
# Designs: a regular grid of n cells with queen adjacency (Besag ICAR random
# effect), data from sim_area_data(D = n, W = W, n_unsampled = 0, seed = 42 + n):
#   gaussian  y_gaussian ~ x1 + x2, spatial = "besag"
#   beta      y_beta ~ x1 + x2, family = "beta", spatial = "besag", vardir = y_beta (1 - y_beta) / 26
#             (tipsae: likelihood = "beta", spatial_error = TRUE)
#
# Writes `--out` (default inst/extdata/hb_engine_benchmark.rds) with
#   $timings    machine, design, engine, n, time_first, time_warm, device, rhat_max
#   $agreement  machine, design, engine, n, cor_hb, max_abs_diff_hb, median_ratio_sd
#               (every engine against INLA on the same machine and data)
#   $meta       machine, hardware, versions and settings
# The file is rewritten after every engine/design, and combinations already in
# it are skipped, so an interrupted run resumes where it stopped.
#
# Usage: Rscript benchmarks/bench_hb_area_engines.R [--n=100,500,1000,2500]
#          [--designs=gaussian,beta] [--engines=inla,cpu,metal,tipsae] [--reps=2]
#          [--tipsae-max-n=2500] [--machine="Apple M4"] [--out=path.rds]
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit[1])
}
split_arg <- function(x) strsplit(x, ",", fixed = TRUE)[[1]]

# Regular grid of n cells (as in bench_bayes_fastsae.R)
make_grid_sf <- function(n_val) {
  factors <- which(n_val %% 1:n_val == 0)
  nx <- factors[which.min(abs(factors - sqrt(n_val)))]
  ny <- n_val / nx
  if (nx > ny) {
    tmp <- nx; nx <- ny; ny <- tmp
  }
  poly <- sf::st_polygon(list(matrix(c(0, 0, 1, 0, 1, 1, 0, 1, 0, 0), ncol = 2, byrow = TRUE)))
  grid <- sf::st_make_grid(poly, n = c(nx, ny))
  sf::st_sf(domain = as.character(1:n_val), geometry = grid[1:n_val])
}

make_data <- function(n) {
  grid_sf <- make_grid_sf(n)
  nb <- spdep::poly2nb(grid_sf)
  W <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
  dimnames(W) <- list(grid_sf$domain, grid_sf$domain)
  df <- fastsae::sim_area_data(D = n, W = W, n_unsampled = 0, seed = 42 + n)$data
  df$domain <- as.character(df$domain)
  # Sampling variance of the Beta direct estimate, y (1 - y) / (kappa + 1) with the precision
  # kappa = 25 of the simulator, computed from the estimate as a survey would. The Gaussian
  # `vardir` exceeds mu (1 - mu) / 2 in about half of the areas, which floors the precision at 1.
  df$vardir_beta <- df$y_beta * (1 - df$y_beta) / 26
  list(df = df, nb = nb, grid_sf = grid_sf)
}

# Returns data.frame(domain, hb, sd) plus device / rhat_max attributes
fit_engine <- function(design, engine, d) {
  if (engine == "tipsae") {
    fit <- tipsae::fit_sae(
      formula_fixed = y_beta ~ x1 + x2, data = d$df, domains = "domain",
      disp_direct = "vardir_beta", type_disp = "var", likelihood = "beta",
      spatial_error = TRUE, spatial_df = d$grid_sf, domains_spatial_df = "domain",
      temporal_error = FALSE, refresh = 0, seed = 42
    )
    est <- summary(fit)$model_estimates
    out <- data.frame(domain = as.character(est$Domains), hb = est$mean, sd = est$sd)
    attr(out, "device") <- "cpu (Stan)"
    attr(out, "rhat_max") <- tryCatch(
      max(rstan::summary(fit$stanfit)$summary[, "Rhat"], na.rm = TRUE),
      error = function(e) NA_real_
    )
    return(out)
  }
  eng <- switch(engine,
    inla = fastsae::inla_engine(),
    fastsae::jax_engine(device = engine, progress = FALSE)
  )
  fml <- if (design == "gaussian") y_gaussian ~ x1 + x2 else y_beta ~ x1 + x2
  vardir <- if (design == "gaussian") "vardir" else "vardir_beta"
  fit <- fastsae::hb_area(fml, vardir = vardir, data = d$df, domain = "domain",
                          family = design, spatial = "besag", W = d$nb,
                          engine = eng, print_result = FALSE)
  out <- data.frame(domain = as.character(fit$df_hb$domain), hb = fit$df_hb$hb, sd = fit$df_hb$sd)
  attr(out, "device") <- if (engine == "inla") "cpu (INLA)" else fit$device
  attr(out, "rhat_max") <- fit$rhat_max %||% NA_real_
  out
}

# ---- Child: one engine, one design, all n --------------------------------------
if ("--child" %in% args) {
  suppressPackageStartupMessages(library(fastsae))
  engine <- arg_value("engine", "inla")
  design <- arg_value("design", "gaussian")
  ns <- as.integer(split_arg(arg_value("n", "100")))
  reps <- as.integer(arg_value("reps", "2"))
  out_file <- arg_value("out", tempfile(fileext = ".rds"))
  res <- list()
  for (n in ns) {
    d <- make_data(n)
    times <- numeric(reps)
    for (r in seq_len(reps)) {
      times[r] <- system.time(est <- fit_engine(design, engine, d))[["elapsed"]]
    }
    device <- attr(est, "device")
    if (engine %in% c("metal", "cuda") && !grepl(engine, device, ignore.case = TRUE)) {
      stop("Requested device '", engine, "' but sampling ran on '", device, "'.")
    }
    warm <- if (reps > 1) min(times[-1]) else times[1]
    cat(sprintf("%-8s %-6s n = %6d  first %8.2fs  warm %8.2fs  %s\n",
                design, engine, n, times[1], warm, device))
    res[[as.character(n)]] <- list(
      timing = data.frame(design = design, engine = engine, n = n,
                          time_first = times[1], time_warm = warm, device = device,
                          rhat_max = attr(est, "rhat_max")),
      est = data.frame(domain = est$domain, hb = est$hb, sd = est$sd)
    )
    saveRDS(res, out_file)  # partial results survive an interruption
  }
  if (engine %in% c("cpu", "metal", "cuda")) {
    md <- reticulate::import("importlib.metadata")
    pkgs <- c("jax", "jaxlib", "jax-metal", "jax-cuda12-plugin", "numpyro")
    attr(res, "python") <- vapply(pkgs, function(p) {
      tryCatch(md$version(p), error = function(e) NA_character_)
    }, character(1))
  }
  saveRDS(res, out_file)
  quit(save = "no")
}

# ---- Driver -------------------------------------------------------------------
ns <- as.integer(split_arg(arg_value("n", "100,500,1000,2500")))
designs <- split_arg(arg_value("designs", "gaussian,beta"))
engines <- split_arg(arg_value("engines", "inla,cpu,metal,tipsae"))
reps <- arg_value("reps", "2")
tipsae_max_n <- as.integer(arg_value("tipsae-max-n", "2500"))
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
dest <- arg_value("out", file.path(dirname(script), "..", "inst", "extdata", "hb_engine_benchmark.rds"))

sys_out <- function(cmd, args) {
  out <- tryCatch(suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = FALSE)),
                  error = function(e) character(0))
  if (length(out) == 0 || !is.null(attr(out, "status"))) NA_character_ else trimws(out[1])
}
cpu_name <- if (Sys.info()[["sysname"]] == "Darwin") {
  sys_out("sysctl", c("-n", "machdep.cpu.brand_string"))
} else {
  sub(".*:\\s*", "", grep("model name", readLines("/proc/cpuinfo", warn = FALSE), value = TRUE)[1])
}
gpu_name <- if (nzchar(Sys.which("nvidia-smi"))) {
  sys_out("nvidia-smi", c("--query-gpu=name", "--format=csv,noheader"))
} else if (Sys.info()[["sysname"]] == "Darwin") {
  sub(".*Chipset Model:\\s*", "", grep("Chipset Model", system2("system_profiler",
    "SPDisplaysDataType", stdout = TRUE), value = TRUE)[1])
} else {
  NA_character_
}
ram_gb <- if (Sys.info()[["sysname"]] == "Darwin") {
  as.numeric(sys_out("sysctl", c("-n", "hw.memsize"))) / 1024^3
} else {
  as.numeric(sub("\\D*(\\d+).*", "\\1", grep("MemTotal", readLines("/proc/meminfo"), value = TRUE))) / 1024^2
}
has_nvidia <- !is.na(gpu_name) && nzchar(Sys.which("nvidia-smi"))
machine <- arg_value("machine", if (dir.exists("/kaggle") && has_nvidia) {
  paste0("Kaggle, ", gpu_name)
} else if (dir.exists("/content") && has_nvidia) {
  paste0("Google Colab, ", gpu_name)
} else if (!is.na(cpu_name) && grepl("Apple", cpu_name)) {
  cpu_name
} else {
  Sys.info()[["nodename"]]
})

# Resume from an earlier, interrupted run of the same output file
runs <- list()
if (file.exists(dest)) {
  prev <- readRDS(dest)
  if (!is.null(prev$runs)) runs <- prev$runs
}

save_results <- function(runs) {
  all_runs <- runs
  # Report only the requested grid (an earlier run may have used other sizes)
  runs <- lapply(runs, function(r) {
    keep <- r[names(r) %in% as.character(ns)]
    attr(keep, "python") <- attr(r, "python")
    keep
  })
  runs <- runs[lengths(runs) > 0]
  timings <- do.call(rbind, lapply(runs, function(r) do.call(rbind, lapply(r, `[[`, "timing"))))
  rownames(timings) <- NULL
  timings <- cbind(machine = machine, timings)
  agreement <- do.call(rbind, lapply(names(runs), function(key) {
    parts <- strsplit(key, " ")[[1]]
    ref <- runs[[paste(parts[1], "inla")]]
    if (parts[2] == "inla" || is.null(ref)) return(NULL)
    run <- runs[[key]]
    do.call(rbind, lapply(intersect(names(run), names(ref)), function(n) {
      a <- run[[n]]$est
      b <- ref[[n]]$est
      i <- match(b$domain, a$domain)
      data.frame(machine = machine, design = parts[1], engine = parts[2], n = as.integer(n),
                 cor_hb = stats::cor(a$hb[i], b$hb),
                 max_abs_diff_hb = max(abs(a$hb[i] - b$hb)),
                 median_ratio_sd = stats::median(a$sd[i] / b$sd))
    }))
  }))
  py <- Filter(Negate(is.null), lapply(runs, attr, "python"))
  meta <- list(
    date = format(Sys.Date()),
    machine = machine,
    cpu = cpu_name,
    cpu_cores = parallel::detectCores(),
    gpu = gpu_name,
    ram_gb = round(ram_gb, 1),
    r_version = R.version.string,
    fastsae = as.character(utils::packageVersion("fastsae")),
    inla = tryCatch(as.character(utils::packageVersion("INLA")), error = function(e) NA_character_),
    tipsae = tryCatch(as.character(utils::packageVersion("tipsae")), error = function(e) NA_character_),
    python = if (length(py) > 0) py[[length(py)]] else NULL,
    jax_settings = "jax_engine() defaults: 2 chains x 1000 draws, 500 warmup",
    tipsae_settings = "tipsae::fit_sae() defaults",
    model = "Besag ICAR on a regular grid with queen adjacency",
    reps = as.integer(reps)
  )
  saveRDS(list(timings = timings, agreement = agreement, meta = meta, runs = all_runs), dest)
  invisible(list(timings = timings, agreement = agreement))
}

for (design in designs) {
  for (engine in engines) {
    if (engine == "tipsae" && design != "beta") next
    key <- paste(design, engine)
    n_todo <- if (engine == "tipsae") ns[ns <= tipsae_max_n] else ns
    n_todo <- setdiff(n_todo, as.integer(names(runs[[key]])))
    if (length(n_todo) == 0) next
    out <- tempfile(fileext = ".rds")
    status <- system2(file.path(R.home("bin"), "Rscript"),
      c(shQuote(script), "--child", paste0("--engine=", engine), paste0("--design=", design),
        paste0("--n=", paste(n_todo, collapse = ",")), paste0("--reps=", reps),
        paste0("--out=", out)))
    if (file.exists(out)) {
      new <- readRDS(out)
      runs[[key]] <- c(runs[[key]], new)
      attr(runs[[key]], "python") <- attr(new, "python")
      save_results(runs)
    }
    if (status != 0) stop("Benchmark failed: ", key, " (partial results kept in ", dest, ")")
  }
}

res <- save_results(runs)
print(res$timings)
print(res$agreement)
cat("Saved", normalizePath(dest), "\n")
