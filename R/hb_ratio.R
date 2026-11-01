#' Hierarchical Bayes Ratio Model for a Sub-Proportion of a Parent Indicator
#'
#' @description Estimates an indicator that is a subset of a parent indicator (for instance
#' extreme poverty within poverty, severe stunting within stunting) as the product of the parent
#' estimate and a modelled sub-proportion: \eqn{\theta_d = s_d P_d}, where \eqn{P_d} is a given
#' (e.g. officially released) estimate of the parent indicator and \eqn{s_d} is the share of the
#' parent population that has the sub-characteristic. The sub-proportion follows a binomial
#' area-level model (\code{\link{hb_area}} with \code{family = "binomial"}) whose trials are the
#' (effective) numbers of sampled parent cases. Since \eqn{0 \le s_d \le 1}, the estimates never
#' exceed the parent (\eqn{\theta_d \le P_d}), and an optional benchmark on the logit scale of
#' \eqn{s_d} keeps that property.
#'
#' @param formula Formula for the direct sub-proportion, \code{s ~ x1 + x2}: the response is the
#'   direct estimate of \eqn{s_d} (for instance the direct ratio of the sub-indicator to the
#'   parent), or the number of parent cases with the sub-characteristic.
#' @param data Data frame with one row per area.
#' @param trials Vector, column name, or formula: (effective) numbers of sampled parent cases,
#'   e.g. from \code{\link{effective_trials}} with the parent cases as \code{n}.
#' @param parent Vector, column name, or formula: the parent estimate \eqn{P_d} (non-negative).
#' @param parent_se Optional vector, column name, or formula: standard errors of \code{parent}.
#'   With them, \code{sd}, \code{mse}, \code{rse} and the intervals include the uncertainty of
#'   the parent (treated as independent of \eqn{s_d}); without them they are conditional on
#'   \code{parent}.
#' @param domain Optional column name or formula of area identifiers.
#' @param benchmark_target Optional benchmark target for the weighted mean of
#'   \eqn{\theta_d} (scalar, or named by \code{benchmark_group}).
#' @param benchmark_weight Vector, column name, or formula of area weights for the benchmark
#'   (e.g. population sizes); needed with \code{benchmark_target}.
#' @param benchmark_group Optional vector, column name, or formula: groups with their own target.
#' @param n_draws Number of posterior draws (default 2000).
#' @param seed Integer seed of the posterior draws and of the parent draws (default 1);
#'   \code{NULL} uses the session's random number stream.
#' @param print_result Logical. If \code{TRUE} (default), prints a summary.
#' @param ... Further arguments to \code{\link{hb_area}} (e.g. \code{spatial}, \code{W},
#'   \code{engine = inla_engine(...)}); \code{family} is \code{"binomial"}.
#'
#' @details The benchmark solves, for each posterior draw, for one shift \eqn{a_g} per group on
#' the logit scale such that \eqn{\sum_{d \in g} w_d P_d\, \mathrm{logit}^{-1}(\mathrm{logit}(s_d) +
#' a_g) / \sum_{d \in g} w_d = T_g}, so the benchmarked \eqn{s_d} stay in (0, 1) and the estimates
#' below their parent. The parent enters at its given value; the parent draws for the total
#' uncertainty are \eqn{\max(P_d + \mathrm{se}_d z, 0)}, \eqn{z \sim N(0, 1)}.
#'
#' Modelling the sub-proportion borrows strength for the part of the indicator that the parent
#' does not already explain, and it keeps the two indicators coherent. In an application to
#' extreme poverty within poverty for 514 Indonesian districts, the ratio model gave a relative
#' standard error of at most 50\% in 124 districts against 77 for a binomial HB model of the
#' extreme poverty rate itself, and no district estimate above its poverty rate (against 7).
#'
#' Only the INLA engine is supported.
#'
#' @return An object of class \code{c("fastsae_hb_ratio", "fastsae_hb", "fastsae")} with
#'   \describe{
#'     \item{df_hb}{One row per area: \code{domain}, \code{y} (direct sub-proportion times
#'       \code{parent}), \code{hb} (posterior mean of \eqn{\theta_d}), \code{sd}, \code{mse},
#'       \code{rse}, \code{ci_lower}, \code{ci_upper}, \code{rse_cond} (RSE given the parent),
#'       \code{ratio}, \code{ratio_sd} (posterior mean and SD of the possibly benchmarked
#'       \eqn{s_d}), \code{ratio_direct}, \code{parent}, \code{parent_se} and \code{trials}.}
#'     \item{draws}{Matrix of posterior draws of \eqn{\theta_d} (areas x draws) behind \code{df_hb}
#'       (for instance for \code{add_reliability_flags(draws = )}).}
#'     \item{ratio_fit}{The \code{\link{hb_area}} fit of the sub-proportion (use it for
#'       \code{diagnose()} of the model).}
#'     \item{benchmark}{\code{NULL}, or a list with the targets, groups and the weighted means
#'       of the benchmarked estimates.}
#'   }
#'   and \code{estcoef}, \code{random_effect_var}, \code{goodness}, \code{convergence},
#'   \code{inla_mode_status}, \code{formula}, \code{model} and \code{call}.
#'
#' @seealso \code{\link{hb_area}}, \code{\link{effective_trials}}, \code{\link{benchmark_sae}}.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   set.seed(1)
#'   D <- 60
#'   d <- data.frame(area = 1:D, x = rnorm(D), pop = rpois(D, 5e4))
#'   d$p0 <- plogis(-2 + 0.3 * rnorm(D))                 # released parent rates
#'   d$n_poor <- rpois(D, 40)                             # sampled parent cases
#'   s_true <- plogis(-1.5 + 0.5 * d$x + rnorm(D, 0, 0.3))
#'   d$s_dir <- rbinom(D, d$n_poor, s_true) / d$n_poor
#'   fit <- hb_ratio(s_dir ~ x,
#'     data = d, trials = "n_poor", parent = "p0", domain = "area",
#'     benchmark_target = 0.025, benchmark_weight = "pop", print_result = FALSE
#'   )
#'   head(fit$df_hb)
#'   all(fit$df_hb$hb <= d$p0)
#' }
#' }
#' @export
hb_ratio <- function(formula,
                     data,
                     trials,
                     parent,
                     parent_se = NULL,
                     domain = NULL,
                     benchmark_target = NULL,
                     benchmark_weight = NULL,
                     benchmark_group = NULL,
                     n_draws = 2000,
                     seed = 1L,
                     print_result = TRUE,
                     ...) {
  call_matched <- match.call()
  dots <- list(...)
  if (!is.null(dots$family) && !identical(dots$family, "binomial")) {
    cli::cli_abort("{.fn hb_ratio} models the sub-proportion with {.code family = \"binomial\"}.")
  }
  dots$family <- NULL
  if (!is.null(dots$self_benchmark) && isTRUE(dots$self_benchmark)) {
    cli::cli_abort("Use {.arg benchmark_target} of {.fn hb_ratio}, not {.arg self_benchmark}: it benchmarks the product with the parent.")
  }
  engine <- dots$engine %||% "inla"
  if (identical(engine, "jax") || inherits(engine, "fastsae_engine") && !identical(engine$name, "inla")) {
    cli::cli_abort("{.fn hb_ratio} supports the INLA engine only.")
  }
  if (!is.numeric(n_draws) || length(n_draws) != 1 || is.na(n_draws) || n_draws < 10) {
    cli::cli_abort("{.arg n_draws} must be a number of at least 10.")
  }
  n_draws <- as.integer(n_draws)
  data <- as.data.frame(data)
  D <- nrow(data)
  parent_vec <- as.numeric(.get_variable(data, parent))
  if (anyNA(parent_vec) || any(parent_vec < 0)) {
    cli::cli_abort("{.arg parent} must be non-negative and known in every area.")
  }
  parent_se_vec <- if (!is.null(parent_se)) as.numeric(.get_variable(data, parent_se))
  if (!is.null(parent_se_vec) && (anyNA(parent_se_vec) || any(parent_se_vec < 0))) {
    cli::cli_abort("{.arg parent_se} must be non-negative and known in every area.")
  }
  sb_info <- NULL
  if (!is.null(benchmark_target)) {
    if (is.null(benchmark_weight)) {
      cli::cli_abort("{.arg benchmark_target} needs {.arg benchmark_weight} (e.g. population sizes).")
    }
    w <- as.numeric(.get_variable(data, benchmark_weight))
    if (anyNA(w) || any(w < 0) || sum(w) <= 0) {
      cli::cli_abort("{.arg benchmark_weight} must be non-negative, known, and not all zero.")
    }
    grp <- if (!is.null(benchmark_group)) as.character(.get_variable(data, benchmark_group)) else rep("all", D)
    ug <- unique(grp)
    tgt <- if (length(benchmark_target) == 1 && is.null(names(benchmark_target))) {
      stats::setNames(rep(benchmark_target, length(ug)), ug)
    } else {
      benchmark_target[ug]
    }
    if (anyNA(tgt)) {
      cli::cli_abort("{.arg benchmark_target} must be one number or be named by every group of {.arg benchmark_group}.")
    }
    # w_d P_d / sum_g w: the weighted mean of s_d P_d is sum(w_norm * s_d)
    w_norm <- w * parent_vec / stats::ave(w, grp, FUN = sum)
    reach <- tapply(w_norm, grp, sum)[ug]
    if (any(tgt <= 0 | tgt >= reach)) {
      bad <- ug[tgt <= 0 | tgt >= reach]
      cli::cli_abort(c(
        "The benchmark target must be positive and below the weighted mean of {.arg parent} (the value at s = 1).",
        "x" = "Not so for group{?s} {.val {bad}}."
      ))
    }
    sb_info <- list(unique_groups = ug, group_vec = grp, w_norm = w_norm, target_vec = tgt)
  }

  fit_s <- do.call(hb_area, c(
    list(formula = formula, data = data, domain = domain, family = "binomial", trials = trials, print_result = FALSE),
    dots
  ))
  if (nrow(fit_s$df_hb) != D) {
    cli::cli_abort("{.fn hb_ratio} needs one row of {.arg data} per area.")
  }

  draw_seed <- if (is.null(seed)) 0L else as.integer(seed)
  s_draws <- .with_seed(seed, {
    smp <- .inla_draws(fit_s$fit, n_draws, draw_seed)
    rn <- rownames(smp[[1]]$latent)
    pred_rows <- grep("^Predictor:", rn)[seq_len(D)]
    eta <- vapply(smp, function(s) s$latent[pred_rows, 1], numeric(D))
    .inla_inv_link(eta, "binomial", fit_s$link)
  })
  if (!is.null(sb_info)) s_draws <- .benchmark_adjust(s_draws, sb_info, "logit")
  level_cond <- s_draws * parent_vec
  level <- if (!is.null(parent_se_vec)) {
    z <- .with_seed(if (is.null(seed)) NULL else seed + 1L, matrix(stats::rnorm(D * n_draws), D))
    s_draws * pmax(parent_vec + parent_se_vec * z, 0)
  } else {
    level_cond
  }

  # df_hb$y of a binomial hb_area() fit holds the numbers of cases
  ratio_direct <- fit_s$df_hb$y / fit_s$df_hb$trials
  est <- rowMeans(level_cond)
  sd_lv <- apply(level, 1, stats::sd)
  rse_of <- function(sd) ifelse(est > 0, 100 * sd / est, NA_real_)
  df <- data.frame(
    domain = fit_s$df_hb$domain,
    y = ratio_direct * parent_vec,
    hb = est,
    sd = sd_lv,
    mse = sd_lv^2,
    rse = rse_of(sd_lv),
    ci_lower = apply(level, 1, stats::quantile, probs = 0.025, names = FALSE),
    ci_upper = apply(level, 1, stats::quantile, probs = 0.975, names = FALSE),
    rse_cond = rse_of(apply(level_cond, 1, stats::sd)),
    ratio = rowMeans(s_draws),
    ratio_sd = apply(s_draws, 1, stats::sd),
    ratio_direct = ratio_direct,
    parent = parent_vec,
    parent_se = if (is.null(parent_se_vec)) NA_real_ else parent_se_vec,
    trials = fit_s$df_hb$trials,
    stringsAsFactors = FALSE
  )

  bm <- NULL
  if (!is.null(sb_info)) {
    w <- as.numeric(.get_variable(data, benchmark_weight))
    bm <- list(
      target = sb_info$target_vec,
      group = sb_info$group_vec,
      achieved = tapply(w * est, sb_info$group_vec, sum)[sb_info$unique_groups] /
        tapply(w, sb_info$group_vec, sum)[sb_info$unique_groups]
    )
  }

  out <- list(
    df_hb = df,
    hb = df,
    df_eblup = df,
    draws = level,
    ratio_fit = fit_s,
    benchmark = bm,
    estcoef = fit_s$estcoef,
    random_effect_var = fit_s$random_effect_var,
    goodness = fit_s$goodness,
    convergence = fit_s$convergence,
    inla_mode_status = fit_s$inla_mode_status,
    engine = fit_s$engine,
    family = "binomial",
    formula = formula,
    data = data,
    model = paste0("Ratio model: parent x ", fit_s$model, " sub-proportion"),
    call = call_matched
  )
  class(out) <- c("fastsae_hb_ratio", "fastsae_hb", "fastsae")
  if (print_result) print(out)
  out
}
