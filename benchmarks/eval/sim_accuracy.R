#!/usr/bin/env Rscript

# ==============================================================================
# Monte Carlo accuracy study of the evaluation article (ADEMP design, see
# vignettes/articles/fastsae-evaluation.Rmd):
#   sfh   eblup_sfh() vs sae::mseSFH() vs emdi::fh(correlation = "spatial"), REML,
#         analytical MSE; truth theta_d = x_d'beta + u_d
#   beta  hb_area(family = "beta", spatial = "bym") with INLA and NumPyro (CPU) vs
#         tipsae::fit_sae(spatial_error = TRUE); truth mu_d
#
# Every (model, scenario, method, replicate) is saved to
# benchmarks/output/eval/<model>/<scenario>_<method>_<rep>.rds and skipped when it
# exists, so an interrupted run resumes. Data depend only on the scenario and the
# replicate, so all methods see the same data sets. A failed fit is saved as a
# failure (with its error message), not dropped.
#
# Usage (from the package root, after R CMD INSTALL . --preclean):
#   Rscript benchmarks/eval/sim_accuracy.R --list=sfh|beta     # one task per line
#   Rscript benchmarks/eval/sim_accuracy.R --task=sfh:base:fastsae:1-100
#   Rscript benchmarks/eval/sim_accuracy.R --collect           # -> inst/extdata/eval_*_accuracy.rds
# The tasks are run in parallel by benchmarks/eval/run_overnight.sh.
# ==============================================================================

suppressPackageStartupMessages(library(fastsae))
source("benchmarks/eval/eval_common.R")

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit[1])
}
out_root <- Sys.getenv("EVAL_OUT", "benchmarks/output/eval")

# ---- scenarios (one factor changed at a time from the base scenario) --------------
sfh_base <- list(map = "grid", D = 100, rho = 0.5, sigma2 = 1, ratio = 1)
sfh_scen <- list(
  base = list(),
  rho_low = list(rho = 0.2), rho_high = list(rho = 0.8), rho_edge = list(rho = 0.95),
  psi_low = list(ratio = 0.5), psi_high = list(ratio = 2),
  D50 = list(D = 50), D500 = list(D = 500),
  mys = list(map = "mys")
)
beta_base <- list(map = "grid", D = 50, sigma_s = 0.5, sigma_v = 0.2, kappa = 25)
beta_scen <- list(
  base = list(),
  spat_low = list(sigma_s = 0.2), spat_high = list(sigma_s = 1),
  kappa_low = list(kappa = 10), kappa_high = list(kappa = 100),
  voronoi = list(map = "voronoi")
)
scenario <- function(model, name) {
  base <- if (model == "sfh") sfh_base else beta_base
  modifyList(base, (if (model == "sfh") sfh_scen else beta_scen)[[name]])
}
methods <- list(sfh = c("fastsae", "sae", "emdi"), beta = c("inla", "jax", "tipsae"))
R_reps <- c(sfh = 500L, beta = 100L)
block <- c(sfh = 100L, beta = 25L)
# fixed seeds: design 1000 + scenario index, data 1e5 * scenario index + replicate
scen_index <- function(model, name) match(name, names(if (model == "sfh") sfh_scen else beta_scen))

# ---- list tasks -----------------------------------------------------------------
if (!is.null(arg_value("list"))) {
  model <- arg_value("list")
  sc_names <- names(if (model == "sfh") sfh_scen else beta_scen)
  starts <- seq(1L, R_reps[[model]], by = block[[model]])
  tasks <- expand.grid(start = starts, method = methods[[model]], scen = sc_names, stringsAsFactors = FALSE)
  # slowest method first, so the long tasks do not start last
  tasks <- tasks[order(-match(tasks$method, methods[[model]])), ]
  cat(sprintf("%s:%s:%s:%d-%d\n", model, tasks$scen, tasks$method, tasks$start,
              pmin(tasks$start + block[[model]] - 1L, R_reps[[model]])), sep = "")
  quit(save = "no")
}

# ---- run one task -----------------------------------------------------------------
task <- arg_value("task")
if (!is.null(task)) {
  p <- strsplit(task, ":", fixed = TRUE)[[1]]
  model <- p[1]; scen <- p[2]; method <- p[3]
  reps <- as.integer(strsplit(p[4], "-", fixed = TRUE)[[1]])
  reps <- seq(reps[1], reps[2])
  sc <- scenario(model, scen)
  k <- scen_index(model, scen)
  dir <- file.path(out_root, model)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  des <- if (model == "sfh") sfh_design(sc, 1000 + k) else beta_design(sc, 1000 + k)
  for (r in reps) {
    f <- file.path(dir, sprintf("%s_%s_%03d.rds", scen, method, r))
    if (file.exists(f)) next
    dat <- if (model == "sfh") sfh_replicate(des, sc, 1e5 * k + r) else beta_replicate(des, sc, 1e5 * k + r)
    warns <- character()
    t <- system.time(fit <- withCallingHandlers(
      tryCatch(
        if (model == "sfh") fit_sfh(method, dat, des$W_std) else fit_beta(method, dat, des, seed = r),
        error = function(e) list(error = conditionMessage(e))
      ),
      warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }
    ))[["elapsed"]]
    fit$time <- t
    fit$warnings <- unique(warns)
    fit$truth <- if (model == "sfh") dat$theta else dat$mu
    fit$vardir <- dat$vardir
    saveRDS(fit, f)
  }
  cat(sprintf("done %s\n", task))
  quit(save = "no")
}

# ---- collect --------------------------------------------------------------------
if ("--collect" %in% args) {
  for (model in c("sfh", "beta")) {
    dir <- file.path(out_root, model)
    files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
    if (length(files) == 0) next
    rows <- list(); areas <- list()
    for (f in files) {
      x <- readRDS(f)
      key <- regmatches(basename(f), regexec("^(.*)_([a-z]+)_([0-9]+)\\.rds$", basename(f)))[[1]]
      scen <- key[2]; method <- key[3]; rep <- as.integer(key[4])
      sc <- scenario(model, scen)
      failed <- !is.null(x$error)
      rows[[length(rows) + 1]] <- data.frame(
        scenario = scen, method = method, rep = rep, failed = failed,
        error = if (failed) x$error else NA_character_,
        n_warnings = length(x$warnings), time = x$time,
        beta0 = if (failed) NA else x$beta[1], beta1 = if (failed) NA else x$beta[2],
        beta2 = if (failed) NA else x$beta[3],
        sigma2 = if (model == "sfh" && !failed) x$sigma2 else NA_real_,
        rho = if (model == "sfh" && !failed) x$rho else NA_real_,
        converged = if (model == "sfh" && !failed) x$converged else NA,
        sigma_v = if (model == "beta" && !failed) x$sigma_v else NA_real_,
        sigma_s = if (model == "beta" && !failed) x$sigma_s else NA_real_,
        rhat_max = if (model == "beta" && !failed) x$rhat_max else NA_real_,
        ess_min = if (model == "beta" && !failed) x$ess_min else NA_real_,
        divergences = if (model == "beta" && !failed) x$divergences else NA_real_
      )
      if (!failed) {
        half <- if (model == "sfh") 1.96 * sqrt(pmax(x$mse, 0)) else NULL
        areas[[length(areas) + 1]] <- data.frame(
          scenario = scen, method = method, rep = rep, area = seq_along(x$truth),
          truth = x$truth, est = x$est,
          lower = if (model == "sfh") x$est - half else x$lower,
          upper = if (model == "sfh") x$est + half else x$upper,
          mse = if (model == "sfh") x$mse else x$sd^2
        )
      }
    }
    scen_list <- if (model == "sfh") sfh_scen else beta_scen
    out <- list(
      fits = do.call(rbind, rows),
      areas = do.call(rbind, areas),
      scenarios = lapply(stats::setNames(names(scen_list), names(scen_list)), function(n) scenario(model, n)),
      truth_par = if (model == "sfh") c(beta0 = 1, beta1 = 0.8, beta2 = -0.5) else c(beta0 = -1, beta1 = 0.5, beta2 = -0.4),
      meta = list(
        date = Sys.Date(), r = R.version.string,
        cpu = tryCatch(system("sysctl -n machdep.cpu.brand_string", intern = TRUE), error = function(e) NA_character_),
        versions = vapply(c("fastsae", "sae", "emdi", "tipsae", "INLA", "rstan"), function(p)
          tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_), ""),
        R_reps = R_reps[[model]], design = "model-based; design (x, sampling variances) fixed per scenario"
      )
    )
    saveRDS(out, file.path(Sys.getenv("EVAL_EXTDATA", "inst/extdata"), sprintf("eval_%s_accuracy.rds", model)))
    cat(sprintf("%s: %d fits (%d failed) -> inst/extdata/eval_%s_accuracy.rds\n",
                model, nrow(out$fits), sum(out$fits$failed), model))
  }
}
