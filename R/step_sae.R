#' Stepwise Selection of the Covariates of a Small Area Model
#'
#' Selects the covariates of a fitted \pkg{fastsae} model by stepwise search (backward,
#' forward or both) on an information criterion, refitting the model with the same
#' settings (family, spatial or temporal structure, engine, ...) at every step.
#'
#' @param object A fitted model: \code{eblup_fh()}, \code{eblup_sfh()}, \code{eblup_stfh()},
#'   \code{eblup_bhf()}, \code{eblup_twofold()}, \code{hb_area()}, \code{hb_unit()} or
#'   \code{hb_twofold()}.
#' @param scope The range of models searched: a formula with the largest set of covariates
#'   (the upper model), or a list with formulas \code{lower} and \code{upper}. Default: the
#'   covariates of \code{object} as upper model and the intercept-only model as lower model.
#'   Forward selection needs an upper model larger than \code{object}.
#' @param direction \code{"both"} (default), \code{"backward"} or \code{"forward"}.
#' @param criterion Information criterion to minimise. EBLUP models: \code{"AIC"},
#'   \code{"BIC"} (default) or \code{"KIC"}. HB models: \code{"LCPO"} (default with
#'   \code{inla_engine()}), \code{"WAIC"} (default with \code{jax_engine()}) or \code{"DIC"}.
#' @param steps Maximum number of steps. Default 1000.
#' @param trace Logical; print the criterion of each candidate model. Default \code{TRUE}.
#'
#' @details
#' \strong{EBLUP models.} The criteria are computed from the full Gaussian log-likelihood
#' \eqn{\ell} at the estimates, with \eqn{k} the number of coefficients plus variance
#' parameters: \eqn{\mathrm{AIC} = -2\ell + 2k}, \eqn{\mathrm{BIC} = -2\ell + k \log m} and
#' \eqn{\mathrm{KIC} = -2\ell + 3k} (Marhuenda, Morales and Pardo 2014; the analytic criteria of
#' \code{emdi::step()}). The restricted (REML) likelihood depends on the covariates through
#' \eqn{X} and cannot compare models with different covariates, so models with a
#' \code{method} argument are fitted by ML during the search; the selected model is then
#' refitted with the method of \code{object}. \code{eblup_stfh()} has no ML fit and is
#' compared through the likelihood at its REML estimates. Only sampled domains enter the
#' criteria.
#'
#' \strong{HB models} (every family). \code{"LCPO"}, the negative mean log conditional
#' predictive ordinate \eqn{-\frac{1}{m}\sum_d \log \mathrm{CPO}_d} (leave-one-out
#' cross-validation; Held, Schrödle and Rue 2010), and \code{"WAIC"} (Watanabe 2010) estimate
#' the out-of-sample predictive accuracy; \code{"DIC"} (Spiegelhalter et al. 2002) is less
#' reliable for hierarchical models. With one observation per area and an area random effect,
#' WAIC is unstable (Vehtari, Gelman and Gabry 2017): in a simulation of a binomial area-level
#' model (40 areas, two true and two noise covariates) WAIC selected exactly the true covariates
#' in 27\% of the replicates and dropped a true one in 13\%, LCPO in 70\%
#' (\code{benchmarks/validation/step_sae.R}). LCPO is therefore the default with
#' \code{inla_engine()}; \code{jax_engine()} has no CPO and uses WAIC. With \code{jax_engine()}
#' every candidate model is an MCMC fit, so the search can take long; the criteria then also
#' carry Monte Carlo error.
#'
#' \strong{Spatial and temporal models.} Every candidate is fitted with the spatial or temporal
#' structure of \code{object}. Selecting the covariates without the spatial effect and adding
#' it afterwards is not equivalent: a spatially smooth covariate and the spatial random effect
#' explain the same variation (spatial confounding), so a covariate kept by a non-spatial model
#' may be redundant, and its coefficient can change markedly, once the spatial effect is in the
#' model (Reich, Hodges and Zadnik 2006; Hodges and Reich 2010). In the simulation above, with a
#' spatial area effect and a spatially smooth noise covariate, the Fay-Herriot model kept that
#' covariate in 20\% of the replicates and the spatial Fay-Herriot model in 8\%. Check the residual spatial
#' autocorrelation of the selected model, e.g. with \code{\link{spatial_test}()}.
#'
#' The model is refitted by evaluating its call with a new formula in the environment where
#' \code{step_sae()} is called, so the data and other arguments of the call must be
#' available there.
#'
#' @return An object of class \code{"fastsae_step"}: a list with \code{model} (the selected
#'   model, fitted with the settings of \code{object}), \code{path} (a data frame with the
#'   term added or dropped at each step and the criterion), \code{criterion},
#'   \code{direction}, and the \code{initial} and \code{final} formulas.
#'
#' @references
#' Held, L., Schrödle, B. and Rue, H. (2010). Posterior and cross-validatory predictive checks:
#' a comparison of MCMC and INLA. In T. Kneib and G. Tutz (eds.), \emph{Statistical Modelling
#' and Regression Structures}, 91-110. Physica-Verlag.
#'
#' Hodges, J. S. and Reich, B. J. (2010). Adding spatially-correlated errors can mess up the
#' fixed effect you love. \emph{The American Statistician}, 64(4), 325-334.
#'
#' Marhuenda, Y., Morales, D. and Pardo, M. C. (2014). Information criteria for Fay-Herriot
#' model selection. \emph{Computational Statistics & Data Analysis}, 70, 268-280.
#'
#' Reich, B. J., Hodges, J. S. and Zadnik, V. (2006). Effects of residual smoothing on the
#' posterior of the fixed effects in disease-mapping models. \emph{Biometrics}, 62(4),
#' 1197-1206.
#'
#' Spiegelhalter, D. J., Best, N. G., Carlin, B. P. and van der Linde, A. (2002). Bayesian
#' measures of model complexity and fit. \emph{Journal of the Royal Statistical Society:
#' Series B}, 64(4), 583-639.
#'
#' Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model evaluation using
#' leave-one-out cross-validation and WAIC. \emph{Statistics and Computing}, 27(5), 1413-1432.
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable
#' information criterion in singular learning theory. \emph{Journal of Machine Learning
#' Research}, 11, 3571-3594.
#'
#' @seealso \code{\link{compare_sae}}, \code{\link{spatial_test}}
#' @export
#' @examples
#' data(mys)
#' fit <- eblup_fh(y ~ x1 + x2 + x3, vardir = "vardir", data = mys, print_result = FALSE)
#' sel <- step_sae(fit, criterion = "BIC")
#' sel
#' sel$model$formula
step_sae <- function(object, scope = NULL, direction = c("both", "backward", "forward"),
                     criterion = NULL, steps = 1000, trace = TRUE) {
  if (!inherits(object, "fastsae") || is.null(object$call)) {
    cli::cli_abort("{.arg object} must be a model fitted by {.pkg fastsae}.")
  }
  direction <- match.arg(direction)
  is_hb <- inherits(object, "fastsae_hb")
  allowed <- if (is_hb) c("LCPO", "WAIC", "DIC") else c("AIC", "BIC", "KIC")
  criterion <- criterion %||% if (!is_hb) "BIC" else if (identical(object$engine, "inla")) "LCPO" else "WAIC"
  if (!is.character(criterion) || length(criterion) != 1L || !criterion %in% allowed) {
    cli::cli_abort(c(
      "{.arg criterion} must be one of {.or {.val {allowed}}} for this model.",
      "i" = "EBLUP models use AIC, BIC or KIC; HB models WAIC, DIC or LCPO."
    ))
  }
  if (criterion == "LCPO" && !identical(object$engine, "inla")) {
    cli::cli_abort("{.code criterion = \"LCPO\"} needs a model fitted with {.fn inla_engine}.")
  }
  if (is_hb && identical(object$engine, "jax")) {
    cli::cli_inform(c(
      "i" = "Every candidate model is an MCMC fit with {.fn jax_engine}; the search can take long.",
      "i" = "WAIC is unstable with one observation per area; {.fn inla_engine} allows the leave-one-out criterion LCPO."
    ))
  }
  env <- parent.frame()

  # The search fits EBLUP models by ML (REML likelihoods are not comparable across X)
  call <- object$call
  fun_formals <- names(formals(eval(call[[1L]], env)))
  search_call <- call
  if (!is_hb && "method" %in% fun_formals) search_call$method <- "ML"
  search_call$print_result <- FALSE

  fit_formula <- function(fml) {
    cl <- search_call
    cl$formula <- fml
    suppressWarnings(suppressMessages(eval(cl, env)))
  }
  score <- function(fit) .sae_criterion(fit, criterion)

  start_formula <- stats::formula(object$formula %||% call$formula)
  response <- deparse(start_formula[[2L]])
  terms_of <- function(f) attr(stats::terms(f), "term.labels")
  upper <- terms_of(start_formula)
  lower <- character(0)
  if (!is.null(scope)) {
    if (is.list(scope)) {
      if (!is.null(scope$upper)) upper <- terms_of(scope$upper)
      if (!is.null(scope$lower)) lower <- terms_of(scope$lower)
    } else {
      upper <- terms_of(scope)
    }
  }
  current <- terms_of(start_formula)
  if (!all(lower %in% current)) {
    cli::cli_abort("The lower model must be contained in the model of {.arg object}.")
  }
  build <- function(trms) stats::reformulate(if (length(trms)) trms else "1", response = response)

  current_fit <- if (identical(search_call, call)) object else fit_formula(start_formula)
  current_score <- score(current_fit)
  path <- data.frame(Step = "", Criterion = current_score, stringsAsFactors = FALSE)
  if (trace) cli::cli_text("Start: {criterion} = {round(current_score, 4)} ({deparse(build(current))})")

  for (i in seq_len(steps)) {
    cands <- list()
    if (direction %in% c("both", "backward")) {
      for (t in setdiff(current, lower)) cands[[paste("-", t)]] <- setdiff(current, t)
    }
    if (direction %in% c("both", "forward")) {
      for (t in setdiff(upper, current)) cands[[paste("+", t)]] <- c(current, t)
    }
    if (length(cands) == 0L) break
    res <- lapply(cands, function(trms) {
      fit <- tryCatch(fit_formula(build(trms)), error = function(e) NULL)
      list(fit = fit, score = if (is.null(fit)) NA_real_ else score(fit))
    })
    scores <- vapply(res, `[[`, numeric(1), "score")
    if (trace) {
      tab <- data.frame(Step = c("<none>", names(scores)), Criterion = c(current_score, scores))
      tab <- tab[order(tab$Criterion), ]
      names(tab)[2] <- criterion
      print(tab, row.names = FALSE)
    }
    if (all(is.na(scores)) || min(scores, na.rm = TRUE) >= current_score) break
    best <- which.min(scores)
    current <- cands[[best]]
    current_fit <- res[[best]]$fit
    current_score <- scores[[best]]
    path <- rbind(path, data.frame(Step = names(cands)[best], Criterion = current_score))
    if (trace) cli::cli_text("Step: {names(cands)[best]}: {criterion} = {round(current_score, 4)}")
  }

  final_formula <- build(current)
  final_call <- call
  final_call$formula <- final_formula
  final_call$print_result <- FALSE
  final <- if (identical(search_call, call)) current_fit else suppressMessages(eval(final_call, env))
  names(path)[2] <- criterion

  structure(
    list(
      model = final, path = path, criterion = criterion, direction = direction,
      initial = build(terms_of(start_formula)), final = final_formula
    ),
    class = "fastsae_step"
  )
}

#' Information criterion of a fitted fastsae model, for step_sae()
#' @noRd
.sae_criterion <- function(fit, criterion) {
  g <- fit$goodness
  if (criterion == "LCPO") {
    cpo <- fit$fit$cpo$cpo
    cpo <- cpo[is.finite(cpo) & cpo > 0]
    if (length(cpo) == 0L) cli::cli_abort("The INLA fit has no conditional predictive ordinates.")
    return(-mean(log(cpo)))
  }
  # hb_area(): DIC, WAIC; hb_unit() and hb_twofold(): dic, waic
  if (criterion %in% c("WAIC", "DIC") && !criterion %in% names(g)) criterion <- tolower(criterion)
  if (is.null(g) || (!criterion %in% names(g) && criterion != "KIC")) {
    cli::cli_abort("The model does not report {.val {toupper(criterion)}}.")
  }
  if (criterion == "KIC") {
    # KIC = -2 loglik + 3k, with k from AIC = -2 loglik + 2k
    k <- (g[["AIC"]] + 2 * g[["loglikelihood"]]) / 2
    return(-2 * g[["loglikelihood"]] + 3 * k)
  }
  as.numeric(g[[criterion]])
}

#' @export
print.fastsae_step <- function(x, ...) {
  cli::cli_h1("Stepwise covariate selection ({x$criterion}, {x$direction})")
  cli::cli_text("{.strong Initial model}: {deparse(x$initial)}")
  cli::cli_text("{.strong Final model}: {deparse(x$final)}")
  cat("\n")
  print(x$path, row.names = FALSE)
  invisible(x)
}
