#' Estimation engines for the Hierarchical Bayes models
#'
#' Constructors for the \code{engine} argument of \code{\link{hb_area}},
#' \code{\link{hb_unit}} and \code{\link{hb_twofold}}. Each collects the
#' settings of one estimation engine, so that the model functions themselves
#' only describe the model. Passing the string \code{"inla"} or \code{"jax"}
#' is the same as passing \code{inla_engine()} or \code{jax_engine()} with
#' their defaults.
#'
#' @section INLA:
#' \code{inla_engine()} fits the model by integrated nested Laplace
#' approximation (\pkg{INLA}); it is available for all HB models. Priors are
#' given in \pkg{INLA} syntax. A prior or setting left at \code{NULL} takes the
#' default of the model function; a setting the model does not use (such as
#' \code{prior_rho_time} in \code{hb_unit()}) is ignored with a warning.
#'
#' @section JAX:
#' \code{jax_engine()} fits the model by MCMC (NumPyro's No-U-Turn sampler on
#' JAX, Phan et al. 2019) on CPU, CUDA GPU or Apple Metal, through
#' \pkg{reticulate}. It needs Python with \code{jax} and \code{numpyro}; see
#' \code{\link{setup_numpyro_env}} and \code{\link{check_numpyro_available}}.
#' The fixed effects get a normal or a regularized horseshoe prior (Piironen
#' and Vehtari 2017), and the variance components half-normal priors.
#' \code{spatial = "slm"} and \code{temporal = "rw2"} are not available.
#' \code{jax_engine()} is available for \code{hb_area()}, \code{hb_unit()} and \code{hb_twofold()}.
#'
#' @param strategy INLA approximation strategy: \code{"laplace"} (default, full
#'   Laplace approximation, highest accuracy), \code{"simplified.laplace"}
#'   (faster) or \code{"gaussian"} (fastest).
#' @param link Link function for the likelihood (\code{hb_area()} only).
#'   \code{NULL} (default) uses the canonical link of the family: identity for
#'   gaussian, logit for binomial and beta, log for poisson, nbinomial and gamma.
#'   Binomial and beta also accept \code{"probit"}, \code{"cloglog"} and
#'   \code{"loglog"}. The estimates (\code{hb}) are on the response scale for
#'   any link.
#' @param scale_model Logical; scale the intrinsic GMRF models (besag, bym,
#'   rw1, rw2) to a generalized variance of 1. Default \code{TRUE}.
#' @param prior_prec Prior for the precision of the area (domain) random
#'   effect. \code{NULL} (default): the PC prior
#'   \code{list(prior = "pc.prec", param = c(1, 0.01))} in \code{hb_area()} and
#'   \code{hb_twofold()}, \code{list(prior = "loggamma", param = c(0.01, 0.01))}
#'   in \code{hb_unit()}. The PC prior puts probability 0.01 on a standard
#'   deviation above 1, which suits area-level responses such as proportions,
#'   rates and log scales; unit-level responses are often on a raw scale (e.g.
#'   hectares), where it would over-shrink the area effects, hence the vague
#'   log-gamma prior. For an area-level response with a large scale, set
#'   \code{prior_prec = list(prior = "pc.prec", param = c(U, 0.01))} with
#'   \code{U} a plausible upper value of the area standard deviation.
#' @param prior_phi Prior for the mixing parameter of \code{spatial = "bym2"}.
#'   \code{NULL} (default): the PC prior \code{list(prior = "pc", param = c(0.5, 0.5))}.
#' @param prior_rho Prior for the spatial correlation of
#'   \code{spatial = "generic1"} or \code{"slm"} (\code{hb_area()} only);
#'   \code{NULL} uses the INLA default.
#' @param prior_prec_time Prior for the precision of the temporal random
#'   effect (\code{hb_area()} only); \code{NULL} uses the INLA default.
#' @param prior_rho_time Prior for the autocorrelation of
#'   \code{temporal = "ar1"} (\code{hb_area()} only); \code{NULL} uses the INLA
#'   default.
#' @param prior_prec_subarea Prior for the precision of the subarea random
#'   effect (\code{hb_twofold()} only). \code{NULL} (default):
#'   \code{list(prior = "pc.prec", param = c(1, 0.01))}.
#' @param draws Number of posterior draws from
#'   \code{INLA::inla.posterior.sample()}: for posterior self-benchmarking in
#'   \code{hb_area()} (\code{NULL}: 1000) and for the area aggregates of
#'   \code{hb_twofold()} (\code{NULL}: 200; \code{0} skips the sampling and
#'   approximates the aggregate uncertainty from independent sub-area effects).
#' @param seed Seed of \code{INLA::inla.posterior.sample()} for those draws, also used
#'   by \code{benchmark_sae()}. The default \code{1L} makes them reproducible (INLA
#'   then samples on one thread); \code{0L} is INLA's own default, a random seed.
#' @param ... Further arguments passed to \code{INLA::inla()}, such as
#'   \code{control.compute} or \code{num.threads}.
#'
#' @return An object of class \code{"fastsae_engine"}.
#'
#' @references
#' Phan, D., Pradhan, N. and Jankowiak, M. (2019). Composable effects for
#' flexible and accelerated probabilistic programming in NumPyro.
#' \emph{arXiv:1912.11554}.
#'
#' Piironen, J. and Vehtari, A. (2017). Sparsity information and
#' regularization in the horseshoe and other shrinkage priors.
#' \emph{Electronic Journal of Statistics}, 11(2), 5018--5051.
#'
#' @seealso \code{\link{hb_area}}, \code{\link{hb_unit}}, \code{\link{hb_twofold}}
#' @examples
#' inla_engine(strategy = "simplified.laplace")
#' jax_engine(samples = 2000, chains = 4, device = "cuda")
#' @name hb_engine
NULL

#' @rdname hb_engine
#' @export
inla_engine <- function(
  strategy = c("laplace", "simplified.laplace", "gaussian"),
  link = NULL,
  scale_model = TRUE,
  prior_prec = NULL,
  prior_phi = NULL,
  prior_rho = NULL,
  prior_prec_time = NULL,
  prior_rho_time = NULL,
  prior_prec_subarea = NULL,
  draws = NULL,
  seed = 1L,
  ...
) {
  strategy <- match.arg(tolower(strategy), choices = c("laplace", "simplified.laplace", "gaussian"))
  if (!is.null(link) && !(is.character(link) && length(link) == 1 && !is.na(link))) {
    cli::cli_abort("{.arg link} must be {.code NULL} or a single string.")
  }
  .check_flag(scale_model, "scale_model")
  for (nm in c("prior_prec", "prior_phi", "prior_rho", "prior_prec_time", "prior_rho_time", "prior_prec_subarea")) {
    if (!is.null(get(nm)) && !is.list(get(nm))) {
      cli::cli_abort("{.arg {nm}} must be {.code NULL} or a list in INLA prior syntax.")
    }
  }
  inla_args <- list(...)
  if (length(inla_args) > 0 && (is.null(names(inla_args)) || any(!nzchar(names(inla_args))))) {
    cli::cli_abort("Arguments in {.arg ...} must be named (they are passed to {.fn INLA::inla}).")
  }
  structure(
    list(
      name = "inla",
      strategy = strategy,
      link = link,
      scale_model = scale_model,
      prior_prec = prior_prec,
      prior_phi = prior_phi,
      prior_rho = prior_rho,
      prior_prec_time = prior_prec_time,
      prior_rho_time = prior_rho_time,
      prior_prec_subarea = prior_prec_subarea,
      draws = if (is.null(draws)) NULL else .check_count(draws, "draws", min = 0),
      seed = .check_count(seed, "seed", min = 0),
      inla_args = inla_args
    ),
    class = c("fastsae_engine_inla", "fastsae_engine")
  )
}

#' @rdname hb_engine
#' @param warmup Number of warmup (adaptation) iterations per chain. Default 500.
#' @param samples Number of posterior draws kept per chain. Default 1000.
#' @param chains Number of MCMC chains. Default 2. With more than one chain,
#'   the fit reports the largest split R-hat and \code{convergence} is
#'   \code{rhat_max < 1.05}.
#' @param device \code{"auto"} (default: CUDA when available, otherwise CPU),
#'   \code{"cpu"}, \code{"cuda"} or \code{"metal"} (Apple silicon GPU, through
#'   jax-metal). JAX fixes its device at the first fit of an R session.
#'   \code{"metal"} has a fixed cost of about two minutes per fit, so on Apple
#'   silicon it pays off only for large numbers of domains (thousands); for
#'   that reason \code{"auto"} uses the CPU on macOS.
#' @param seed Random seed of the sampler. Default 42.
#' @param prior_beta Prior for the regression coefficients: \code{"normal"}
#'   (default) or \code{"horseshoe"} (regularized horseshoe, for many
#'   candidate covariates).
#' @param robust Logical; Student-t instead of normal random effects (the
#'   area effects, and the sub-area effects of \code{hb_twofold()}), with the
#'   degrees of freedom estimated. Default \code{FALSE}.
#' @param keep_draws Logical; keep the posterior draws of the area parameters
#'   in \code{$draws} of the fit (needed by \code{benchmark_sae()} to propagate
#'   the posterior covariance). Default \code{TRUE}.
#' @param progress Logical; show NumPyro's progress bar. Default
#'   \code{interactive()}.
#' @export
jax_engine <- function(
  warmup = 500,
  samples = 1000,
  chains = 2,
  device = c("auto", "cpu", "cuda", "metal"),
  seed = 42,
  prior_beta = c("normal", "horseshoe"),
  robust = FALSE,
  keep_draws = TRUE,
  progress = interactive()
) {
  .check_flag(robust, "robust")
  .check_flag(keep_draws, "keep_draws")
  .check_flag(progress, "progress")
  structure(
    list(
      name = "jax",
      warmup = .check_count(warmup, "warmup", min = 0),
      samples = .check_count(samples, "samples", min = 1),
      chains = .check_count(chains, "chains", min = 1),
      device = match.arg(tolower(device), choices = c("auto", "cpu", "cuda", "metal")),
      seed = .check_count(seed, "seed", min = 0),
      prior_beta = match.arg(tolower(prior_beta), choices = c("normal", "horseshoe")),
      robust = robust,
      keep_draws = keep_draws,
      progress = progress
    ),
    class = c("fastsae_engine_jax", "fastsae_engine")
  )
}

#' @rdname hb_engine
#' @param x A \code{"fastsae_engine"} object.
#' @export
print.fastsae_engine <- function(x, ...) {
  if (x$name == "jax") {
    cli::cli_text("{.strong JAX engine} (NumPyro NUTS)")
    cli::cli_bullets(c(
      "*" = "{x$chains} chain{?s} x {x$samples} draws, {x$warmup} warmup, seed {x$seed}",
      "*" = "Device: {x$device}",
      "*" = "Prior on coefficients: {x$prior_beta}{if (x$robust) ', Student-t random effects' else ''}"
    ))
  } else {
    cli::cli_text("{.strong INLA engine}")
    bullets <- c(
      "*" = "Strategy: {x$strategy}",
      "*" = "Link: {if (is.null(x$link)) 'canonical' else x$link}"
    )
    if (length(x$inla_args) > 0) {
      inla_names <- names(x$inla_args)
      bullets <- c(bullets, "*" = "Passed to {.fn INLA::inla}: {.arg {inla_names}}")
    }
    cli::cli_bullets(bullets)
  }
  invisible(x)
}

#' Coerce the `engine` argument of hb_area() to an engine object
#' @noRd
.as_hb_engine <- function(engine) {
  if (inherits(engine, "fastsae_engine")) {
    return(engine)
  }
  if (is.character(engine) && length(engine) == 1 && !is.na(engine)) {
    name <- tolower(engine)
    if (name == "laplace") {
      cli::cli_abort(c(
        "{.code engine = \"laplace\"} has been removed.",
        "i" = "Use {.code engine = inla_engine(strategy = \"laplace\")}."
      ))
    }
    if (name == "inla") return(inla_engine())
    if (name == "jax") return(jax_engine())
  }
  cli::cli_abort(c(
    "{.arg engine} must be {.val inla}, {.val jax}, or an object from {.fn inla_engine} or {.fn jax_engine}.",
    "x" = "You supplied an object of class {.cls {class(engine)}}."
  ))
}

#' Map the deprecated `method` argument and INLA arguments given to a model function
#'
#' fastsae 0.1.0 took `method`, `strategy`, `link`, `scale_model`, `prior_*`,
#' `benchmark_draws`, `benchmark_seed`, `n_samples` and INLA arguments directly
#' in `hb_area()`, `hb_unit()` and `hb_twofold()`. They now belong to
#' `inla_engine()`; old code keeps working with a deprecation warning.
#' @param renamed Named character vector: old argument name -> `inla_engine()` field.
#' @param fun Name of the calling function, for the messages.
#' @return An engine object.
#' @noRd
.legacy_engine_args <- function(engine, engine_missing, method = NULL, dots = list(),
                                renamed = .legacy_inla_args$hb_area, fun = "hb_area") {
  if (!is.null(method)) {
    method <- match.arg(tolower(method), choices = c("inla", "laplace"))
    if (!engine_missing) {
      cli::cli_abort("Supply {.arg engine} only; {.arg method} is deprecated.")
    }
    new_call <- if (method == "laplace") "engine = inla_engine(strategy = \"laplace\")" else "engine = \"inla\""
    cli::cli_warn(c(
      "{.arg method} is deprecated; use {.arg engine} instead.",
      "i" = "{.code method = \"{method}\"} is now {.code {new_call}}.",
      if (method == "laplace") c("i" = "The lme4 GLMM route was removed; the model is fitted by INLA.")
    ))
    engine <- if (method == "laplace") inla_engine(strategy = "laplace") else inla_engine()
  }
  spec <- .as_hb_engine(engine)
  if (length(dots) == 0) {
    return(spec)
  }

  dot_names <- names(dots)
  if (is.null(dot_names) || any(!nzchar(dot_names))) {
    cli::cli_abort("Arguments in {.arg ...} must be named.")
  }
  jax_args <- intersect(dot_names, c(setdiff(names(formals(jax_engine)), "seed"), "control"))
  if (length(jax_args) > 0) {
    cli::cli_abort(c(
      "{.arg {jax_args}} {?is/are} not {?an argument/arguments} of {.fn {fun}}.",
      "i" = "Set MCMC options in {.code engine = jax_engine(...)}."
    ))
  }
  if (spec$name != "inla") {
    cli::cli_abort(c(
      "{.arg {dot_names}} {?is an INLA argument/are INLA arguments} and cannot be used with {.fn jax_engine}.",
      "i" = "INLA settings go in {.fn inla_engine}."
    ))
  }

  cli::cli_warn(c(
    "Passing {.arg {dot_names}} to {.fn {fun}} is deprecated.",
    "i" = "Use {.code engine = inla_engine(...)} instead."
  ))
  fields <- dots[intersect(dot_names, names(renamed))]
  names(fields) <- renamed[names(fields)]
  extra <- dots[setdiff(dot_names, names(renamed))]

  # Settings given both ways are ambiguous unless the engine is a plain default
  if (!engine_missing && inherits(engine, "fastsae_engine")) {
    defaults <- inla_engine()
    clash <- c(
      names(fields)[!vapply(names(fields), function(f) identical(engine[[f]], defaults[[f]]), logical(1))],
      intersect(names(extra), names(engine$inla_args))
    )
    if (length(clash) > 0) {
      cli::cli_abort("{.arg {clash}} {?is/are} set both in {.fn inla_engine} and in {.fn {fun}}.")
    }
  }
  args <- spec[setdiff(names(spec), c("name", "inla_args"))]
  args[names(fields)] <- fields
  do.call(inla_engine, c(args, utils::modifyList(spec$inla_args, extra)))
}

#' Old argument names of the model functions -> inla_engine() fields
#' @noRd
.legacy_inla_args <- list(
  hb_area = c(
    strategy = "strategy", link = "link", scale_model = "scale_model",
    prior_prec = "prior_prec", prior_phi = "prior_phi", prior_rho = "prior_rho",
    prior_prec_time = "prior_prec_time", prior_rho_time = "prior_rho_time",
    benchmark_draws = "draws", benchmark_seed = "seed"
  ),
  hb_unit = c(
    strategy = "strategy", scale_model = "scale_model",
    prior_prec = "prior_prec", prior_phi = "prior_phi"
  ),
  hb_twofold = c(
    strategy = "strategy", scale_model = "scale_model",
    prior_prec_area = "prior_prec", prior_prec_subarea = "prior_prec_subarea",
    prior_phi = "prior_phi", n_samples = "draws"
  )
)

#' Resolve the INLA engine of a model function
#'
#' Coerces `engine` (with the deprecated arguments in `dots`), rejects other
#' engines when `inla_only`, fills the settings left at `NULL` with the
#' function's `defaults`, and warns about settings the function does not use.
#' @param used `inla_engine()` fields the function reads (besides `strategy`,
#'   `seed` and `inla_args`).
#' @noRd
.resolve_inla_engine <- function(engine, engine_missing, dots, fun, defaults, used, method = NULL) {
  engine <- .legacy_engine_args(engine, engine_missing, method = method, dots = dots,
                                renamed = .legacy_inla_args[[fun]], fun = fun)
  if (engine$name != "inla") {
    return(engine)
  }
  base <- inla_engine()
  fields <- setdiff(names(base), c("name", "strategy", "seed", "inla_args", used))
  unused <- fields[!vapply(fields, function(f) identical(engine[[f]], base[[f]]), logical(1))]
  if (length(unused) > 0) {
    cli::cli_warn("{.arg {unused}} of {.fn inla_engine} {?is/are} not used by {.fn {fun}} and ignored.")
  }
  for (f in names(defaults)) {
    if (is.null(engine[[f]])) engine[f] <- list(defaults[[f]])
  }
  engine
}

#' Abort when a model function gets an engine it does not support
#' @noRd
.require_inla_engine <- function(engine, fun) {
  if (engine$name != "inla") {
    cli::cli_abort(c(
      "{.fn {fun}} is fitted with INLA only.",
      "i" = "Use {.code engine = \"inla\"} or {.fn inla_engine}."
    ))
  }
  invisible(engine)
}

#' Validate a single logical flag
#' @noRd
.check_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    cli::cli_abort("{.arg {arg}} must be {.code TRUE} or {.code FALSE}.")
  }
  invisible(x)
}

#' Validate a single whole number >= min and return it as integer
#' @noRd
.check_count <- function(x, arg, min = 0) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < min || x %% 1 != 0) {
    cli::cli_abort("{.arg {arg}} must be a whole number >= {min}.")
  }
  as.integer(x)
}
