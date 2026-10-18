#' Hierarchical Bayes for Area-Level Small Area Estimation
#'
#' @description
#' Estimates small area parameters using area-level Hierarchical Bayesian models
#' under various distributions (Gaussian, Binomial, Poisson, Negative Binomial, Beta, Gamma)
#' and spatial random effect structures (non-spatial, BYM2, BYM, Besag/ICAR, Leroux, SAR).
#' Two estimation engines are available through \code{engine}:
#' \itemize{
#'   \item \code{\link{inla_engine}()} (default, or \code{"inla"}): Integrated Nested
#'     Laplace Approximations (\pkg{INLA}; Rue et al., 2009).
#'   \item \code{\link{jax_engine}()} (or \code{"jax"}): full MCMC with the No-U-Turn sampler
#'     of 'NumPyro' (Phan et al., 2019) on 'JAX', on the CPU or an NVIDIA GPU. Needs Python
#'     with \code{jax} and \code{numpyro} through \pkg{reticulate}; see
#'     \code{\link{setup_numpyro_env}}.
#' }
#' The engine constructors hold all engine settings (priors, approximation strategy and
#' further \code{INLA::inla()} arguments for INLA; MCMC length, device and coefficient prior
#' for JAX), so \code{hb_area()} itself only describes the model.
#'
#' @section Engines:
#' Both engines return the same structure, so \code{print}, \code{summary},
#' \code{\link{benchmark_sae}}, \code{\link{compare_sae}} and the other methods work
#' with either. Feature support:
#' \tabular{lcc}{
#'   \strong{Feature} \tab \strong{inla} \tab \strong{jax} \cr
#'   Families \tab all six \tab all six \cr
#'   \code{spatial} \tab all \tab all but \code{"slm"} \cr
#'   \code{temporal} \tab all \tab all but \code{"rw2"} \cr
#'   \code{st_interaction} \tab all \tab all \cr
#'   Self-benchmarking \tab posterior, wfq \tab posterior; wfq (gaussian) \cr
#'   Horseshoe prior, Student-t random effects \tab -- \tab yes \cr
#' }
#' The \code{"jax"} engine uses weakly informative priors (half-normal standard deviations,
#' uniform spatial mixing / autocorrelation parameters) instead of the PC priors of INLA, so its
#' estimates differ slightly. For \code{family = "gamma"} it estimates a common shape and does
#' not use \code{vardir}.
#'
#' With the \code{"inla"} engine, Gaussian models and the Leroux and SLM effects of any family are
#' fitted in INLA's \code{"classic"} mode (Rue, Martino and Chopin, 2009); the other models use
#' INLA's default \code{"compact"} mode (van Niekerk et al., 2023), which leaves the linear
#' predictors out of the latent field and is designed for data-rich models. Gaussian models with
#' a \code{"type1"} to \code{"type4"} interaction and three hyperparameters (a Besag or iid area
#' effect and a random walk or iid in time) also use a grid over the hyperparameters
#' (\code{int.strategy = "grid"}); with \code{temporal = "ar1"} or a BYM, BYM2, Leroux or SLM area
#' effect the grid grows exponentially (13 minutes for 12 areas with AR(1)) and was not validated,
#' so INLA's default integration is kept. These choices rest on comparisons with long Stan runs with the
#' same model and priors (Gaussian, binomial, Poisson, beta, gamma and negative binomial; iid,
#' Besag, BYM2, Leroux and SLM effects; RW1 in time and the four Knorr-Held interactions):
#' \itemize{
#'   \item Gaussian: compact overstated the posterior SDs by up to 18\% when the area variance was
#'     small relative to \code{vardir}; classic matched Stan (means within 0.07 posterior SD, SDs
#'     within 0.91-1.03 times) at the same speed. With an interaction both modes overstated the SDs
#'     by 6-13\%; classic with the grid matched Stan (means within 0.04 SD, SDs within 2\%) and
#'     took 5-8 times longer.
#'   \item Binomial, Poisson, beta and gamma with iid, Besag or BYM2 effects or an interaction:
#'     compact matched Stan (means within 0.16 SD, SDs within 0.95-1.08 times), and classic was
#'     much slower for many areas (beta, 2000 areas: 32 s against 1.3 s).
#'   \item Leroux and SLM: compact moved the Poisson posterior means by 0.12-0.15 SD; classic was
#'     within 0.03 SD.
#'   \item Negative binomial: with one count per area the overdispersion and the area effect are
#'     only weakly identified, and both modes were off (means up to 0.64 SD); the model gives a
#'     warning.
#' }
#' \code{inla_engine(inla.mode = , control.inla = )} overrides these settings. For Gaussian models
#' with a small area variance, the classic summary of that variance (\code{random_effect_var}) was
#' biased downwards, while the area estimates were accurate. The validation scripts are in the
#' package repository (\code{benchmarks/validation}).
#'
#' @param formula An object of class \code{formula} specifying the fixed-effects model
#'   (e.g., \code{y ~ x1 + x2}).
#' @param data A \code{data.frame} containing the area-level data (one row per area, or per
#'   area and period). Areas without a direct estimate (unsampled areas) are kept as rows
#'   with \code{NA} response and their covariates; they are predicted from the model. An
#'   area with no row in \code{data} has no covariates and is not predicted.
#' @param domain Vector, column name, or one-sided formula referencing the area/domain
#'   identifier in \code{data}. If \code{NULL}, domains are numbered consecutively.
#' @param family Character string specifying the response distribution / likelihood.
#'   Options:
#'   \itemize{
#'     \item \code{"gaussian"}: Continuous response (Fay-Herriot model). If \code{vardir}
#'       is provided, sampling variances are treated as known. If \code{vardir = NULL},
#'       residual variance is estimated.
#'     \item \code{"binomial"}: Binary / proportion response. Requires \code{trials}.
#'     \item \code{"poisson"}: Count response / disease rate. Supports \code{exposure}.
#'     \item \code{"nbinomial"}: Negative Binomial for overdispersed counts.
#'     \item \code{"beta"}: Continuous proportions strictly in (0, 1). If \code{vardir}
#'       is provided, area-specific precision is calculated via Janicki (2020) formula
#'       \eqn{\phi_i = y_i(1 - y_i)/V_i - 1} and injected via INLA's \code{scale} argument.
#'       If \code{trials} is provided, precision is scaled as \eqn{n_i - 1}.
#'     \item \code{"gamma"}: Skewed positive continuous response.
#'   }
#' @param spatial Character string specifying the spatial random effect structure:
#'   \itemize{
#'     \item \code{"none"}: Non-spatial IID area random intercept (default).
#'     \item \code{"bym2"}: Scaled Besag-York-Mollié 2 model (Riebler et al., 2016),
#'       decomposing area variance into spatial (\eqn{\phi}) and unstructured components.
#'     \item \code{"bym"}: Classic Besag-York-Mollié (1991) model.
#'     \item \code{"besag"}: Intrinsic Conditional Autoregressive (ICAR) model.
#'     \item \code{"generic1"}: Leroux spatial autoregressive model.
#'     \item \code{"sar"}: Simultaneous autoregressive (SAR) area effect,
#'       \eqn{\eta = X\beta + u}, \eqn{u = (I - \rho W)^{-1}\epsilon},
#'       \eqn{\epsilon \sim N(0, \sigma^2 I)}: the model of \code{\link{eblup_sfh}}
#'       (Petrucci and Salvati, 2006), with \eqn{\rho} uniform on \eqn{|\rho| \le 0.999}.
#'       With INLA it is the \code{"slm"} latent model with a zero covariate matrix, so that
#'       only \eqn{\epsilon} is filtered. Not combinable with the \code{"type3"},
#'       \code{"type4"} and \code{"separable"} interactions, which are CAR fields.
#'     \item \code{"slm"}: Spatial lag model, \eqn{\eta = (I - \rho W)^{-1}(X\beta + \epsilon)}
#'       (INLA \code{"slm"}). Note that the covariates are spatially filtered too, so this
#'       is \emph{not} the Bayesian counterpart of \code{\link{eblup_sfh}}, whose SAR
#'       structure applies to the random effect only: \eqn{X\beta + (I - \rho W)^{-1}\epsilon}.
#'       Use \code{"sar"} for that model.
#'   }
#' @param W Proximity or spatial adjacency matrix. Can be a square \code{matrix},
#'   \code{Matrix}, or \code{spdep} \code{nb} or \code{listw} object. Dimensions must
#'   match the total number of domains in \code{data}. Required when \code{spatial != "none"}.
#'   The CAR models (\code{"besag"}, \code{"bym"}, \code{"bym2"}, \code{"leroux"}) use only the
#'   neighbourhood pattern \code{W > 0}, not the weights: give a neighbourhood graph (for
#'   example \code{create_weights(method = "queen")} or \code{"knn"}). A dense matrix such as
#'   distance-decay weights becomes a nearly complete graph and gives a warning. \code{"sar"}
#'   and \code{"slm"} use the weights, row-standardised (as \code{\link{eblup_sfh}} expects).
#' @param vardir Vector, column name, or formula specifying the sampling variances of the
#'   direct estimator (for \code{family = "gaussian"} or \code{family = "beta"}).
#' @param trials Vector, column name, or formula specifying sample sizes / total trials
#'   per area (for \code{family = "binomial"}).
#' @param exposure Vector, column name, or formula specifying expected counts or population
#'   exposure offsets (for \code{family = "poisson"} or \code{"nbinomial"}).
#' @param engine Estimation engine: \code{\link{inla_engine}()} or \code{\link{jax_engine}()},
#'   or the strings \code{"inla"} (default) and \code{"jax"} for their default settings.
#'   See the Engines section.
#' @param time Vector, column name, or one-sided formula referencing the time period
#'   identifier in \code{data} (e.g. \code{time = "year"} or \code{~year}). Required
#'   when \code{temporal != "none"}.
#' @param temporal Character string specifying the temporal random effect structure:
#'   \itemize{
#'     \item \code{"none"}: No temporal random effect (default cross-sectional model).
#'     \item \code{"rw1"}: First-order random walk across time periods.
#'     \item \code{"rw2"}: Second-order random walk across time periods.
#'     \item \code{"ar1"}: First-order autoregressive process across time periods.
#'     \item \code{"iid"}: Unstructured independent time effects.
#'   }
#'   With \code{engine = jax_engine()}, the \code{"rw1"} and \code{"ar1"} effects sum to zero over
#'   the periods (\code{"ar1"} only when the model has an intercept), so the intercept is the
#'   mean level: over a few periods the level of an AR(1) effect otherwise trades off against
#'   the intercept and NUTS diverges. INLA fits \code{"ar1"} without this constraint, so the
#'   intercepts of the two engines differ by the mean of the time effects.
#' @param st_interaction Character string specifying the spatio-temporal structure:
#'   \itemize{
#'     \item \code{"none"}: Additive main spatial and temporal effects (default).
#'     \item \code{"domain-specific"}: Domain-specific temporal random walk / AR(1) dynamics with
#'       shared variance parameter (directly corresponds to the \pkg{tipsae} spatio-temporal model).
#'       Not available with \code{temporal = "rw2"}: without a common trend, the linear trend of
#'       a domain without direct estimates has a flat prior; use \code{"type2"}.
#'     \item \code{"separable"}: Grouped dynamic spatial field evolving over time via AR(1)
#'       (corresponds to classical Spatio-Temporal Fay-Herriot models, Marhuenda et al. 2013).
#'       Only \code{temporal = "ar1"} or \code{"iid"}: with a random walk in time every domain's
#'       level in the field has a flat prior, as the field has no spatial main effect; use
#'       \code{"type4"}.
#'     \item \code{"type1"}: Knorr-Held (2000) Type I interaction (unstructured space \eqn{\times} unstructured time).
#'     \item \code{"type2"}: Knorr-Held Type II interaction (unstructured space \eqn{\times} structured time):
#'       an independent \code{temporal} process for each domain.
#'     \item \code{"type3"}: Knorr-Held Type III interaction (structured space \eqn{\times} unstructured time):
#'       an independent ICAR field for each period.
#'     \item \code{"type4"}: Knorr-Held Type IV interaction (structured space \eqn{\times} structured time):
#'       an ICAR field evolving over time with the \code{temporal} structure.
#'   }
#'   The \code{"type*"} options add the interaction to the main spatial and temporal
#'   effects. \code{"type3"} and \code{"type4"} need a structured \code{spatial} effect
#'   (\code{"bym2"}, \code{"bym"}, \code{"besag"} or \code{"leroux"}); with INLA the interaction
#'   is an ICAR (\code{"besag"}) field on the graph of \code{W}. With \code{temporal = "rw1"} or
#'   \code{"rw2"}, the type IV interaction (and type II with \code{"rw2"}) uses the scaled
#'   Kronecker structure with all the identifiability constraints of Knorr-Held (2000): sums over
#'   the domains are 0 in every period, and sums over the periods (and, for \code{"rw2"}, linear
#'   trends) are 0 in every domain. INLA's default constraints leave the latter out; with them
#'   alone, a domain without direct estimates had a posterior SD of 36-58 (other domains about
#'   0.4).
#' @param self_benchmark Logical. If \code{TRUE}, the weighted aggregate of the area
#'   estimates is made to match the direct survey aggregate or \code{benchmark_target}; see
#'   \code{benchmark_method}.
#' @param benchmark_method How self-benchmarking is done when \code{self_benchmark = TRUE}:
#'   \code{"posterior"} (default) benchmarks every posterior draw of the area parameters
#'   (an additive shift for \code{"gaussian"}, a common shift on the logit scale for
#'   \code{"binomial"} and \code{"beta"}, a common ratio for \code{"poisson"},
#'   \code{"nbinomial"} and \code{"gamma"}), so the constraint holds exactly for the posterior
#'   mean and the reported SD / credible intervals refer to the benchmarked estimator (Bayesian
#'   benchmarking, cf. Datta et al., 2011). \code{"wfq"} uses the augmented-covariate approach of
#'   Wang, Fuller and Qu (2008), which is exact only for the linear Gaussian EBLUP; for the HB
#'   models it holds only approximately (non-flat prior on the coefficients, integration over
#'   the hyperparameters) and not at all for non-Gaussian likelihoods. With INLA,
#'   posterior draws need the package \pkg{sn} (used by
#'   \code{INLA::inla.posterior.sample()}); without it \code{"wfq"} is used with a warning.
#' @param benchmark_weight Optional vector, column name, or one-sided formula referencing the domain
#'   benchmarking weights in \code{data} (e.g. population sizes or shares). Defaults to equal weights.
#' @param benchmark_target Optional numeric scalar or named vector specifying fixed aggregate
#'   benchmark target(s). If \code{NULL} (default), calibrates to the weighted aggregate of the direct estimates.
#' @param benchmark_group Optional vector, column name, or one-sided formula referencing the
#'   grouping / stratum column in \code{data} (e.g. province) for group-specific self-benchmarking.
#' @param print_result Logical. If \code{TRUE} (default), prints a summary of results.
#' @param method Deprecated; use \code{engine}. \code{method = "laplace"} becomes
#'   \code{engine = inla_engine(strategy = "laplace")}; the former \pkg{lme4} GLMM route for
#'   binomial and Poisson models has been removed.
#' @param ... Deprecated. \code{strategy}, \code{link}, \code{scale_model}, \code{prior_*},
#'   \code{benchmark_draws}, \code{benchmark_seed} and arguments of \code{INLA::inla()} are
#'   now given to \code{\link{inla_engine}()}; passed here they are forwarded to it with a
#'   warning.
#'
#' @returns An object of class \code{c("fastsae_hb_area", "fastsae_hb", "fastsae")} containing:
#' \itemize{
#'   \item \code{df_hb}: Data frame with domain estimates, including \code{domain},
#'     \code{time} (if specified), observed \code{y}, predicted \code{hb}, linear predictor \code{linear_pred},
#'     posterior standard error \code{sd}, \code{mse}, relative error \code{rse} (\%),
#'     95\% credible interval (\code{ci_lower}, \code{ci_upper}), and \code{random_effect}.
#'   \item \code{estcoef}: Data frame of estimated regression coefficients.
#'   \item \code{hyperpar}: Data frame of hyperparameter posterior estimates.
#'   \item \code{random_effect_var}: Posterior mean of the area random effect variance (see
#'     the Engines section on its accuracy when the variance is small).
#'   \item \code{random_effect_var_time}: Estimated temporal random effect variance (if temporal).
#'   \item \code{phi}: Estimated spatial variance proportion (for BYM2) or spatial autocorrelation.
#'   \item \code{rho_time}: Estimated temporal autocorrelation parameter (for AR1).
#'   \item \code{goodness}: Model fit metrics (DIC, WAIC, Marginal Log-Likelihood).
#'   \item \code{family}: Response family used.
#'   \item \code{spatial}: Spatial model type used.
#'   \item \code{temporal}: Temporal model type used.
#'   \item \code{st_interaction}: Spatio-temporal interaction type used.
#'   \item \code{self_benchmark}: Logical indicating whether self-benchmarking was active.
#'   \item \code{benchmark_summary}: Summary table of benchmark calibration (if self-benchmarking).
#'   \item \code{fit}: Raw fitted model object (\code{inla}, or the list returned by the
#'     Python backend for the JAX engine).
#'   \item \code{engine}: Name of the estimation engine (\code{"inla"} or \code{"jax"});
#'     \code{engine_spec} holds its settings (\code{\link{inla_engine}()} or
#'     \code{\link{jax_engine}()} object).
#'   \item \code{draws}, \code{rhat_max}, \code{device}, \code{elapsed_seconds}: With the
#'     JAX engine, the posterior draws of the area parameters (rows of \code{df_hb}
#'     by draws), the largest split R-hat (\code{convergence} is \code{rhat_max < 1.05}), the
#'     JAX device and the sampling time.
#'   \item \code{call}: Matched function call.
#' }
#'
#' @references
#' \enumerate{
#'   \item Rao, J. N. K., and Molina, I. (2015). \emph{Small Area Estimation}. John Wiley & Sons.
#'   \item Petrucci, A., and Salvati, N. (2006). Small area estimation for spatial correlation in
#'     watershed erosion assessment. \emph{Journal of Agricultural, Biological, and Environmental
#'     Statistics}, 11(2), 169-182.
#'   \item Wang, J., Fuller, W. A., and Qu, Y. (2008). Small area estimation under a restriction.
#'     \emph{Survey Methodology}, 34(1), 29-36.
#'   \item Riebler, A., Sørbye, S. H., Simpson, D., and Rue, H. (2016). An intuitive Bayesian spatial model
#'     for disease mapping that accounts for scaling. \emph{Statistical Methods in Medical Research}, 25(4), 1145-1165.
#'   \item Marhuenda, Y., Molina, I., and Morales, D. (2013). Small area estimation with spatio-temporal Fay-Herriot models.
#'     \emph{Computational Statistics & Data Analysis}, 58, 308-325.
#'   \item De Nicolò, S., and Gardini, A. (2024). The R Package tipsae: Tools for Mapping Proportions and Indicators on the Unit Interval.
#'     \emph{Journal of Statistical Software}, 108(1), 1-36.
#'   \item Phan, D., Pradhan, N., and Jankowiak, M. (2019). Composable effects for flexible and
#'     accelerated probabilistic programming in NumPyro. \emph{arXiv:1912.11554}.
#'   \item Piironen, J., and Vehtari, A. (2017). Sparsity information and regularization in the
#'     horseshoe and other shrinkage priors. \emph{Electronic Journal of Statistics}, 11(2), 5018-5051.
#'   \item Rue, H., Martino, S., and Chopin, N. (2009). Approximate Bayesian inference for latent Gaussian
#'     models by using integrated nested Laplace approximations. \emph{Journal of the Royal Statistical Society: Series B}, 71(2), 319-392.
#'   \item Knorr-Held, L. (2000). Bayesian modelling of inseparable space-time variation in disease
#'     risk. \emph{Statistics in Medicine}, 19(17-18), 2555-2567.
#'   \item van Niekerk, J., Krainski, E., Rustand, D., and Rue, H. (2023). A new avenue for Bayesian
#'     inference with INLA. \emph{Computational Statistics & Data Analysis}, 181, 107692.
#' }
#'
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   library(fastsae)
#'
#'   # 1. Non-spatial Gaussian Fay-Herriot with INLA
#'   m_norm <- hb_area(
#'     y ~ x1 + x2 + x3,
#'     data = mys,
#'     vardir = "vardir",
#'     family = "gaussian"
#'   )
#'
#'   # 2. Self-benchmarking Fay-Herriot with INLA (benchmarked posterior draws)
#'   m_bench <- hb_area(
#'     y ~ x1 + x2 + x3,
#'     data = mys,
#'     vardir = "vardir",
#'     self_benchmark = TRUE,
#'     benchmark_weight = "n"
#'   )
#'
#'   # 3. Spatial BYM2 Gaussian Fay-Herriot with INLA, faster approximation;
#'   #    the CAR model needs a neighbourhood graph (4 nearest centroids)
#'   W_knn <- create_weights(cbind(sim_area$x_coord, sim_area$y_coord), method = "knn", k = 4)
#'   m_bym2 <- hb_area(
#'     y ~ x1 + x2 + x3,
#'     data = mys,
#'     vardir = "vardir",
#'     family = "gaussian",
#'     spatial = "bym2",
#'     W = W_knn,
#'     engine = inla_engine(strategy = "simplified.laplace")
#'   )
#' }
#'
#' # 4. The same Fay-Herriot model by MCMC (NumPyro NUTS on JAX)
#' if (check_numpyro_available()) {
#'   m_jax <- hb_area(
#'     y ~ x1 + x2 + x3,
#'     data = mys,
#'     vardir = "vardir",
#'     engine = jax_engine(samples = 500, chains = 2)
#'   )
#' }
#' }
hb_area <- function(
  formula,
  data,
  domain = NULL,
  time = NULL,
  family = c("gaussian", "binomial", "poisson", "nbinomial", "beta", "gamma"),
  spatial = c("none", "bym2", "bym", "besag", "generic1", "leroux", "sar", "slm"),
  temporal = c("none", "rw1", "rw2", "ar1", "iid"),
  st_interaction = c("none", "domain-specific", "separable", "type1", "type2", "type3", "type4"),
  W = NULL,
  vardir = NULL,
  trials = NULL,
  exposure = NULL,
  engine = "inla",
  self_benchmark = FALSE,
  benchmark_weight = NULL,
  benchmark_target = NULL,
  benchmark_group = NULL,
  benchmark_method = c("posterior", "wfq"),
  print_result = TRUE,
  method = NULL,
  ...
) {
  call_matched <- match.call()
  benchmark_method <- match.arg(benchmark_method)
  family <- match.arg(tolower(family), choices = c("gaussian", "binomial", "poisson", "nbinomial", "beta", "gamma"))
  spatial <- match.arg(tolower(spatial), choices = c("none", "bym2", "bym", "besag", "generic1", "leroux", "sar", "slm"))
  if (spatial == "leroux") spatial <- "generic1" # same model; INLA calls it generic1
  temporal <- match.arg(tolower(temporal), choices = c("none", "rw1", "rw2", "ar1", "iid"))
  st_interaction <- match.arg(tolower(st_interaction), choices = c("none", "domain-specific", "separable", "type1", "type2", "type3", "type4"))
  # A random walk in time leaves one level (rw1) or one linear trend (rw2) per domain
  # with a flat prior when nothing else carries it: the separable field has no
  # spatial main effect, and the domain-specific rw2 no common trend. A domain
  # without direct estimates then got a posterior SD of 40 or more.
  if (st_interaction == "separable" && temporal %in% c("rw1", "rw2")) {
    cli::cli_abort(c(
      "{.code st_interaction = \"separable\"} cannot be combined with {.code temporal = \"{temporal}\"}.",
      "x" = "A random walk in time gives every domain an unpenalised level in the separable field, which has no spatial main effect: domains without direct estimates cannot be predicted.",
      "i" = "Use {.code temporal = \"ar1\"} (a stationary field) or {.code st_interaction = \"type4\"} (spatial and temporal main effects plus their interaction)."
    ))
  }
  if (st_interaction == "domain-specific" && temporal == "rw2") {
    cli::cli_abort(c(
      "{.code st_interaction = \"domain-specific\"} cannot be combined with {.code temporal = \"rw2\"}.",
      "x" = "An rw2 per domain has an unpenalised linear trend in every domain and there is no common trend: domains without direct estimates cannot be predicted.",
      "i" = "Use {.code st_interaction = \"type2\"} (a common trend plus domain deviations) or {.code temporal = \"rw1\"}."
    ))
  }
  if (spatial == "slm" && st_interaction == "separable") {
    cli::cli_abort(c(
      "{.code spatial = \"slm\"} cannot be combined with {.code st_interaction = \"separable\"}.",
      "i" = "Use {.code spatial = \"bym2\"}, {.code \"besag\"} or {.code \"leroux\"} for a spatial field evolving over time."
    ))
  }
  if (spatial == "sar" && st_interaction %in% c("type3", "type4", "separable")) {
    cli::cli_abort(c(
      "{.code spatial = \"sar\"} cannot be combined with {.code st_interaction = \"{st_interaction}\"}.",
      "i" = "This interaction is a CAR field on the graph of {.arg W}; use {.code spatial = \"bym2\"}, {.code \"besag\"} or {.code \"leroux\"}."
    ))
  }
  if (st_interaction %in% c("type3", "type4") && !spatial %in% c("bym2", "bym", "besag", "generic1")) {
    cli::cli_abort(c(
      "{.code st_interaction = \"{st_interaction}\"} needs a structured spatial effect.",
      "i" = "Use {.code spatial = \"bym2\"}, {.code \"bym\"}, {.code \"besag\"} or {.code \"leroux\"} with {.arg W}."
    ))
  }
  engine <- .resolve_inla_engine(
    engine, missing(engine), list(...), "hb_area",
    defaults = list(
      prior_prec = list(prior = "pc.prec", param = c(1, 0.01)),
      prior_phi = list(prior = "pc", param = c(0.5, 0.5)),
      draws = 1000L
    ),
    used = c("link", "scale_model", "prior_prec", "prior_phi", "prior_rho", "prior_prec_time", "prior_rho_time", "draws"),
    method = method
  )
  if (engine$name == "jax") {
    .check_jax_model(
      spatial = spatial, temporal = temporal, family = family,
      self_benchmark = self_benchmark, benchmark_method = benchmark_method
    )
  } else {
    if (isTRUE(self_benchmark) && benchmark_method == "posterior" && engine$draws < 1) {
      cli::cli_abort("Posterior self-benchmarking needs {.code inla_engine(draws = )} of at least 1.")
    }
    .check_inla_link(engine$link, family)
    if (family == "nbinomial") {
      # With one count per area the negative binomial overdispersion and the area
      # effect both describe extra variation between areas; against Stan, INLA's
      # posterior means were off by 0.33-0.64 SD and the SDs by 0.57-1.24 times,
      # in both modes and with grid integration.
      cli::cli_warn(c(
        "With one count per area, the negative binomial overdispersion and the area random effect are only weakly identified, and INLA's approximation of this model can be inaccurate.",
        "i" = "Against long MCMC runs the posterior means were off by up to 0.6 posterior SD. Consider {.code family = \"poisson\"} (the area effect already captures overdispersion between areas) or {.code engine = \"jax\"} (MCMC)."
      ))
    }
  }

  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame or tibble.")
  }

  n_obs <- nrow(data)

  # 1. Extract domain & time identifiers
  if (is.null(domain)) {
    if (is.null(time)) {
      domain_vec <- seq_len(n_obs)
    } else {
      cli::cli_abort("When {.arg time} is specified, {.arg domain} must also be specified to identify domains across time periods.")
    }
  } else {
    domain_vec <- .get_variable(data, domain)
    if (is.null(time)) .check_unique_domain(domain_vec, shared_effect = TRUE)
  }

  unique_domains <- unique(domain_vec)
  n_unique_domains <- length(unique_domains)

  time_vec <- if (!is.null(time)) .get_variable(data, time) else NULL
  if (temporal != "none" && is.null(time_vec)) {
    cli::cli_abort("When {.code temporal != 'none'}, {.arg time} (time/period column) must be specified.")
  }

  unique_times <- if (!is.null(time_vec)) sort(unique(time_vec)) else NULL
  n_unique_times <- if (!is.null(unique_times)) length(unique_times) else 0

  if (temporal != "none" && n_unique_times < 2) {
    cli::cli_abort("Temporal modeling requires at least 2 distinct time periods in {.arg time}.")
  }

  # 2. Model frame & Response validation
  mf <- stats::model.frame(formula, data, na.action = stats::na.pass)
  y <- stats::model.response(mf)

  # Validate optional inputs
  vardir_vec <- if (!is.null(vardir)) .get_variable(data, vardir) else NULL
  trials_vec <- if (!is.null(trials)) .get_variable(data, trials) else NULL
  exposure_vec <- if (!is.null(exposure)) .get_variable(data, exposure) else NULL

  response_var <- all.vars(formula)[1]

  if (family == "binomial") {
    if (is.null(trials_vec)) {
      cli::cli_abort("For {.code family = 'binomial'}, {.arg trials} (sample sizes / total trials per area) must be specified.")
    }
    # Check if response contains proportions (values in [0, 1] with non-integers)
    y_non_na <- y[!is.na(y)]
    if (length(y_non_na) > 0 && all(y_non_na >= 0 & y_non_na <= 1) && any(y_non_na %% 1 != 0)) {
      cli::cli_alert_info("Response {.arg y} appears to be proportions; converting to integer counts: {.code round(y * trials)}.")
      y_counts <- as.integer(round(y * trials_vec))
      data[[response_var]] <- y_counts
      y <- y_counts
    } else {
      y_counts <- ifelse(is.na(y), NA_integer_, as.integer(round(y)))
      data[[response_var]] <- y_counts
      y <- y_counts
    }
  }

  if (family == "gamma") {
    y_non_na <- y[!is.na(y)]
    if (any(y_non_na <= 0, na.rm = TRUE)) {
      cli::cli_abort("For {.code family = 'gamma'}, response {.arg y} must be strictly positive (> 0).")
    }
  }

  # INLA treats a missing covariate as 0, so every row needs complete
  # covariates; the design of the sampled rows must identify beta.
  X_fixed <- stats::model.matrix(attr(mf, "terms"), mf)
  .check_covariates_na(X_fixed, "row")
  if (any(!is.na(y))) {
    .check_design(X_fixed[!is.na(y), , drop = FALSE], "sampled observations")
  }
  if (!is.null(vardir_vec) && family %in% c("gaussian", "beta", "gamma")) {
    .check_vardir(y, vardir_vec, if (is.null(time_vec)) "area" else "observation")
  }

  # 3. Spatial weights handling
  W_obj <- NULL
  if (spatial != "none") {
    if (is.null(W) && inherits(data, "sf")) {
      cli::cli_alert_info("Argument {.arg W} not provided; automatically constructing binary spatial adjacency matrix from {.cls sf} polygon geometry.")
      W <- create_weights(data, method = "queen", style = "B")
    }
    W_obj <- .convert_spatial_weights(W, n_domains = n_unique_domains, spatial = spatial, domain_names = unique_domains)
  }

  # 4. Self-benchmarking setup (Wang, Fuller, and Qu, 2008)
  sb_info <- NULL
  # Posterior benchmarking draws with INLA::inla.posterior.sample(), which needs
  # the package 'sn'; without it fall back to WFQ
  if (isTRUE(self_benchmark) && benchmark_method == "posterior" && engine$name == "inla" && !.has_sn()) {
    cli::cli_warn(c(
      "Package {.pkg sn} (needed by {.fn INLA::inla.posterior.sample}) is not installed; using {.code benchmark_method = \"wfq\"} instead.",
      "i" = "Install it with {.code install.packages(\"sn\")} for posterior benchmarking."
    ))
    benchmark_method <- "wfq"
  }
  if (isTRUE(self_benchmark)) {
    # Direct estimates on the scale of the area parameter (the scale of `hb`):
    # proportions for binomial, rates per unit exposure for counts. The
    # benchmark targets are computed from these.
    y_direct <- y
    if (benchmark_method == "posterior") {
      if (family == "binomial") {
        y_direct <- y / trials_vec
      } else if (family %in% c("poisson", "nbinomial") && !is.null(exposure_vec)) {
        y_direct <- y / exposure_vec
      }
    }
    sb_info <- .setup_self_benchmark(
      data = data,
      formula = formula,
      y = y_direct,
      vardir = vardir_vec,
      weight = benchmark_weight,
      target = benchmark_target,
      group = benchmark_group
    )
    if (benchmark_method == "wfq") {
      data <- sb_info$data_aug
      formula <- sb_info$formula_aug
    }
    sb_info$method <- benchmark_method
  }

  # --------------------------------------------------------------------------
  # Route estimation: NumPyro/JAX MCMC or INLA
  # --------------------------------------------------------------------------
  if (engine$name == "jax") {
    out <- .fit_jax_area(
      formula = formula,
      data = data,
      y = y,
      domain = domain_vec,
      time = time_vec,
      family = family,
      spatial = spatial,
      temporal = temporal,
      st_interaction = st_interaction,
      W_obj = W_obj,
      vardir = vardir_vec,
      trials = trials_vec,
      exposure = exposure_vec,
      engine = engine,
      sb_info = sb_info,
      benchmark_weight = benchmark_weight,
      benchmark_target = benchmark_target,
      benchmark_group = benchmark_group,
      call = call_matched
    )
  } else {
    .check_inla_installed()

    out <- do.call(.fit_inla_area, quote = TRUE, c(
      list(
        formula = formula,
        data = data,
        y = y,
        domain = domain_vec,
        time = time_vec,
        family = family,
        spatial = spatial,
        temporal = temporal,
        st_interaction = st_interaction,
        W_obj = W_obj,
        vardir = vardir_vec,
        trials = trials_vec,
        exposure = exposure_vec,
        strategy = engine$strategy,
        link = engine$link,
        scale_model = engine$scale_model,
        prior_prec = engine$prior_prec,
        prior_phi = engine$prior_phi,
        prior_rho = engine$prior_rho,
        prior_prec_time = engine$prior_prec_time,
        prior_rho_time = engine$prior_rho_time,
        sb_info = sb_info,
        benchmark_weight = benchmark_weight,
        benchmark_target = benchmark_target,
        benchmark_group = benchmark_group,
        benchmark_draws = engine$draws,
        benchmark_seed = engine$seed,
        call = call_matched
      ),
      engine$inla_args
    ))
  }

  out$engine <- engine$name
  out$engine_spec <- engine
  out$W <- W

  if (print_result) {
    print(out)
  }

  return(out)
}

#' Internal worker to fit area-level model via INLA
#' @noRd
.fit_inla_area <- function(
  formula,
  data,
  y,
  domain,
  time = NULL,
  family,
  spatial = "none",
  temporal = "none",
  st_interaction = "none",
  W_obj = NULL,
  vardir = NULL,
  trials = NULL,
  exposure = NULL,
  strategy = "laplace",
  link = NULL,
  scale_model = TRUE,
  prior_prec = list(prior = "pc.prec", param = c(1, 0.01)),
  prior_phi = list(prior = "pc", param = c(0.5, 0.5)),
  prior_rho = NULL,
  prior_prec_time = NULL,
  prior_rho_time = NULL,
  sb_info = NULL,
  benchmark_weight = NULL,
  benchmark_target = NULL,
  benchmark_group = NULL,
  benchmark_draws = 1000,
  benchmark_seed = 1L,
  call = NULL,
  ...
) {
  # Respect CRAN limits and prevent inla.mclapply / OpenMP from exceeding core limits
  chk_cores <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  is_cran_check <- (nzchar(chk_cores) && (chk_cores != "false")) || nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_", ""))

  if (is_cran_check) {
    old_mc <- getOption("mc.cores")
    options(mc.cores = 1L)
    on.exit(options(mc.cores = old_mc), add = TRUE)

    old_inla_threads <- tryCatch(INLA::inla.getOption("num.threads"), error = function(e) NULL)
    if (!is.null(old_inla_threads)) {
      INLA::inla.setOption("num.threads", "1:1")
      on.exit(INLA::inla.setOption("num.threads", old_inla_threads), add = TRUE)
    }
  }

  n_obs <- nrow(data)

  unique_domains <- unique(domain)
  n_domains <- length(unique_domains)
  domain_id <- as.integer(factor(domain, levels = unique_domains))

  inla_data <- data
  inla_data$..domain_id.. <- domain_id
  inla_data$..obs_id.. <- seq_len(n_obs)

  unique_times <- NULL
  n_times <- 1
  time_id <- NULL
  if (!is.null(time)) {
    unique_times <- sort(unique(time))
    n_times <- length(unique_times)
    time_id <- as.integer(factor(time, levels = unique_times))
    inla_data$..time_id.. <- time_id
    inla_data$..time_id2.. <- time_id
  }
  inla_data$..domain_id2.. <- domain_id

  # Map family name to INLA family
  inla_family <- switch(family,
    "gaussian"  = "gaussian",
    "binomial"  = "binomial",
    "poisson"   = "poisson",
    "nbinomial" = "nbinomial",
    "beta"      = "beta",
    "gamma"     = "gamma",
    "gaussian"
  )

  # Construct INLA formula with random effect
  response_var <- all.vars(formula)[1]
  terms_fixed <- attr(stats::terms(formula), "term.labels")
  fixed_str <- if (length(terms_fixed) > 0) paste(terms_fixed, collapse = " + ") else "1"

  X_mat <- NULL
  rho_range <- NULL

  # Setup temporal prior hyperparameters
  if (is.null(prior_prec_time)) {
    prior_prec_time <- list(prior = "pc.prec", param = c(1, 0.01))
  }
  hyper_time <- list(prec = prior_prec_time)
  if (temporal == "ar1" && !is.null(prior_rho_time)) {
    hyper_time$rho <- prior_rho_time
  }

  rand_terms <- c()

  # 1. Spatial random effect specification
  if (spatial == "none") {
    if (temporal == "none" || st_interaction %in% c("none", "domain-specific", "type1", "type2")) {
      rand_terms <- c(rand_terms, "f(..domain_id.., model = 'iid', hyper = list(prec = prior_prec))")
    }
  } else if (spatial == "bym2") {
    if (st_interaction != "separable") {
      rand_terms <- c(rand_terms, paste0(
        "f(..domain_id.., model = 'bym2', graph = W_obj$graph, scale.model = ",
        scale_model,
        ", hyper = list(prec = prior_prec, phi = prior_phi))"
      ))
    }
  } else if (spatial == "bym") {
    if (st_interaction != "separable") {
      rand_terms <- c(rand_terms, paste0(
        "f(..domain_id.., model = 'bym', graph = W_obj$graph, scale.model = ",
        scale_model,
        ", hyper = list(prec.unstruct = prior_prec, prec.spatial = prior_prec))"
      ))
    }
  } else if (spatial == "besag") {
    if (st_interaction != "separable") {
      rand_terms <- c(rand_terms, paste0(
        "f(..domain_id.., model = 'besag', graph = W_obj$graph, scale.model = ",
        scale_model,
        ", hyper = list(prec = prior_prec))"
      ))
    }
  } else if (spatial == "generic1") {
    adj <- W_obj$graph
    R_mat <- diag(rowSums(adj)) - adj
    C_mat <- diag(nrow(adj)) - R_mat
    hyper_generic1 <- list(prec = prior_prec)
    if (!is.null(prior_rho)) {
      hyper_generic1$beta <- prior_rho
    } else if (!is.null(prior_phi) && !identical(prior_phi$prior, "pc")) {
      hyper_generic1$beta <- prior_phi
    }
    if (st_interaction != "separable") {
      rand_terms <- c(rand_terms, "f(..domain_id.., model = 'generic1', Cmatrix = C_mat, hyper = hyper_generic1)")
    }
  } else if (spatial == "sar") {
    # SAR random effect of eblup_sfh(): u = (I - rho W)^{-1} eps. INLA's "slm" is
    # (I - rho W)^{-1} (X beta + eps); with a zero X only eps is filtered and the
    # covariates stay fixed effects. rho is uniform on |rho| <= 0.999, the range of
    # eblup_sfh() (W is row-standardised).
    rho_range <- c(-0.999, 0.999)
    n_sar <- nrow(W_obj$W_mat)
    args_sar <- list(
      rho.min = rho_range[1],
      rho.max = rho_range[2],
      W = Matrix::Matrix(W_obj$W_mat, sparse = TRUE),
      X = matrix(0, n_sar, 1),
      Q.beta = Matrix::Diagonal(1, 1)
    )
    hyper_sar <- list(
      prec = prior_prec,
      rho = if (!is.null(prior_rho)) prior_rho else list(prior = "logitbeta", param = c(1, 1))
    )
    rand_terms <- c(rand_terms, "f(..domain_id.., model = 'slm', args.slm = args_sar, hyper = hyper_sar)")
  } else if (spatial == "slm") {
    W_mat <- W_obj$W_mat
    rs <- rowSums(W_mat)
    rs_inv <- ifelse(rs > 0, 1 / rs, 0)
    W_std <- W_mat * rs_inv
    W_sparse <- Matrix::Matrix(W_std, sparse = TRUE)
    e <- eigen(W_std, only.values = TRUE)$values
    re_e <- Re(e[abs(Im(e)) < 1e-5])
    rho_min <- 1 / min(re_e)
    rho_max <- 1 / max(re_e)
    if (is.infinite(rho_min) || rho_min < -1) rho_min <- -0.999
    if (is.infinite(rho_max) || rho_max > 1) rho_max <- 0.999
    rho_range <- c(rho_min, rho_max)

    terms_no_y <- stats::delete.response(stats::terms(formula))
    mf_X <- stats::model.frame(terms_no_y, data = data, na.action = stats::na.pass)
    X_mat <- stats::model.matrix(terms_no_y, data = mf_X)
    p_cov <- ncol(X_mat)
    Q_beta <- Matrix::Diagonal(p_cov, 1e-4)

    args_slm <- list(
      rho.min = rho_min,
      rho.max = rho_max,
      W = W_sparse,
      X = X_mat,
      Q.beta = Q_beta
    )
    hyper_slm <- list(prec = prior_prec)
    if (!is.null(prior_rho)) {
      hyper_slm$rho <- prior_rho
    } else if (!is.null(prior_phi) && !identical(prior_phi$prior, "pc")) {
      hyper_slm$rho <- prior_phi
    } else {
      hyper_slm$rho <- list(prior = "logitbeta", param = c(1, 1))
    }

    rand_terms <- c(rand_terms, "f(..domain_id.., model = 'slm', args.slm = args_slm, hyper = hyper_slm)")
  }

  # 2. Temporal & Space-Time random effect specification
  if (temporal != "none") {
    scale_model_time <- if (temporal %in% c("rw1", "rw2")) paste0(", scale.model = ", scale_model) else ""

    if (st_interaction == "none") {
      # Additive temporal main effect
      rand_terms <- c(rand_terms, paste0(
        "f(..time_id.., model = '", temporal, "'", scale_model_time, ", hyper = hyper_time)"
      ))
    } else if (st_interaction == "domain-specific") {
      # Domain-specific temporal dynamics (directly corresponds to tipsae)
      rand_terms <- c(rand_terms, paste0(
        "f(..time_id.., model = '", temporal, "', group = ..domain_id.., control.group = list(model = 'iid')",
        scale_model_time, ", hyper = hyper_time)"
      ))
    } else if (st_interaction == "separable") {
      # Dynamic spatial field evolving over time
      if (spatial == "bym2") {
        rand_terms <- c(rand_terms, paste0(
          "f(..domain_id.., model = 'bym2', graph = W_obj$graph, group = ..time_id.., control.group = list(model = '",
          temporal, "'), scale.model = ", scale_model, ", hyper = list(prec = prior_prec, phi = prior_phi))"
        ))
      } else if (spatial == "besag") {
        rand_terms <- c(rand_terms, paste0(
          "f(..domain_id.., model = 'besag', graph = W_obj$graph, group = ..time_id.., control.group = list(model = '",
          temporal, "'), scale.model = ", scale_model, ", hyper = list(prec = prior_prec))"
        ))
      } else if (spatial == "bym") {
        rand_terms <- c(rand_terms, paste0(
          "f(..domain_id.., model = 'bym', graph = W_obj$graph, group = ..time_id.., control.group = list(model = '",
          temporal, "'), scale.model = ", scale_model, ", hyper = list(prec.unstruct = prior_prec, prec.spatial = prior_prec))"
        ))
      } else if (spatial == "generic1") {
        # Grouped generic1 needs a non-negative diagonal of Cmatrix (I - R has not);
        # besagproper2 is the same Leroux model, Q = tau((1 - lambda) I + lambda R)
        hyper_leroux <- list(prec = prior_prec)
        if (!is.null(hyper_generic1$beta)) hyper_leroux$lambda <- hyper_generic1$beta
        rand_terms <- c(rand_terms, paste0(
          "f(..domain_id.., model = 'besagproper2', graph = W_obj$graph, group = ..time_id.., control.group = list(model = '",
          temporal, "'), hyper = hyper_leroux)"
        ))
      } else {
        rand_terms <- c(rand_terms, paste0(
          "f(..domain_id.., model = 'iid', group = ..time_id.., control.group = list(model = '",
          temporal, "'), hyper = list(prec = prior_prec))"
        ))
      }
    } else if (st_interaction == "type1") {
      rand_terms <- c(rand_terms, paste0(
        "f(..time_id.., model = '", temporal, "'", scale_model_time, ", hyper = hyper_time)"
      ))
      rand_terms <- c(rand_terms, "f(..obs_id.., model = 'iid', hyper = list(prec = prior_prec))")
    } else if (st_interaction %in% c("type2", "type3", "type4")) {
      # Knorr-Held (2000) interactions: main effects plus a space-time term on
      # copies of the indices (INLA allows each covariate in one f() only)
      rand_terms <- c(rand_terms, paste0(
        "f(..time_id.., model = '", temporal, "'", scale_model_time, ", hyper = hyper_time)"
      ))
      # With a random walk in time, type IV (and type II with rw2) is a generic0
      # term with all Knorr-Held (2000) constraints (see .st_generic0())
      st_g0 <- .st_generic0(st_interaction, temporal,
        domain_id = domain_id, time_id = time_id,
        n_domains = n_domains, n_times = n_times, graph = W_obj$graph,
        scale_model = scale_model
      )
      if (!is.null(st_g0)) {
        inla_data$..st_id2.. <- st_g0$index
        rand_terms <- c(rand_terms, paste0(
          "f(..st_id2.., model = 'generic0', Cmatrix = st_g0$Cmatrix, rankdef = st_g0$rankdef, ",
          "constr = FALSE, extraconstr = st_g0$extraconstr, hyper = ",
          if (st_interaction == "type2") "hyper_time)" else "list(prec = prior_prec))"
        ))
      } else if (st_interaction == "type2") {
        # unstructured space x structured time: one temporal process per domain
        rand_terms <- c(rand_terms, paste0(
          "f(..time_id2.., model = '", temporal, "', group = ..domain_id.., control.group = list(model = 'iid')",
          scale_model_time, ", hyper = hyper_time)"
        ))
      } else {
        # structured space (ICAR) x unstructured (type3) or structured (type4) time
        group_model <- if (st_interaction == "type3") "iid" else temporal
        group_extra <- if (group_model %in% c("rw1", "rw2")) {
          paste0(", scale.model = ", scale_model)
        } else if (group_model == "ar1" && !is.null(prior_rho_time)) {
          ", hyper = list(rho = prior_rho_time)"
        } else {
          ""
        }
        rand_terms <- c(rand_terms, paste0(
          "f(..domain_id2.., model = 'besag', graph = W_obj$graph, scale.model = ", scale_model,
          ", group = ..time_id.., control.group = list(model = '", group_model, "'", group_extra,
          "), hyper = list(prec = prior_prec))"
        ))
      }
    }
  }

  if (spatial == "slm") {
    inla_formula_str <- paste(response_var, "~ -1 +", paste(rand_terms, collapse = " + "))
  } else {
    inla_formula_str <- paste(response_var, "~", fixed_str, "+", paste(rand_terms, collapse = " + "))
  }

  inla_formula <- stats::as.formula(inla_formula_str)

  # INLA controls
  control_compute <- list(
    dic = TRUE,
    waic = TRUE,
    cpo = TRUE,
    mlik = TRUE,
    config = TRUE
  )

  # Predictor control: compute fitted values on original scale
  control_predictor <- list(compute = TRUE, link = 1)

  # INLA strategy control
  control_inla <- list(strategy = strategy)
  # Gaussian models with a Knorr-Held interaction: INLA's default integration over
  # the three precisions overstated posterior SDs by 6-13% against Stan (posterior
  # means off by up to 0.33 SD) in both modes; a grid matched Stan (<= 0.04 SD,
  # SDs within 2%), about 5-8 times slower. The grid grows exponentially with the
  # number of hyperparameters (ar1 or a two-parameter spatial effect add one or
  # two: 30 s and 13 min for 12 areas against 1 s), and only the three-precision
  # models were validated, so larger models keep INLA's default (CCD).
  # inla_engine(control.inla = ) overrides.
  if (family == "gaussian" && st_interaction %in% c("type1", "type2", "type3", "type4") &&
    .n_hyper_st(spatial, temporal, st_interaction, vardir_known = !is.null(vardir)) <= 3) {
    control_inla <- c(control_inla, list(int.strategy = "grid", dz = 0.5, diff.logdens = 8))
  }

  # Family specific settings
  inla_args <- list(
    formula = inla_formula,
    family = inla_family,
    data = inla_data,
    control.compute = control_compute,
    control.predictor = control_predictor,
    control.inla = control_inla
  )
  # Gaussian models, and the Leroux and SLM effects of any family, fit in INLA's
  # "classic" mode. Against Stan, the default "compact" mode (van Niekerk et al.,
  # 2023) overstated the posterior SDs of Gaussian models by up to 18% when the
  # area variance is small relative to vardir, and moved the posterior means of
  # Poisson Leroux and SLM models by 0.12-0.15 SD; classic matched Stan. For the
  # other non-Gaussian models both modes matched Stan, and classic was much
  # slower for many areas (beta, n = 2000: 32 s against 1.3 s), so they keep
  # INLA's default. inla_engine(inla.mode = ) overrides.
  use_classic <- family == "gaussian" || spatial %in% c("generic1", "sar", "slm")
  if (use_classic && "inla.mode" %in% names(formals(INLA::inla))) inla_args$inla.mode <- "classic"
  # User settings (inla_engine(...)) override; control.* lists are merged
  inla_args <- utils::modifyList(inla_args, list(...))

  # WFQ self-benchmarking: diffuse fixed-effects prior so the augmented covariate is not shrunk
  if (!is.null(sb_info) && identical(sb_info$method, "wfq")) {
    if (is.null(inla_args$control.fixed)) {
      inla_args$control.fixed <- list(
        mean.intercept = 0,
        prec.intercept = 1e-4,
        mean = 0,
        prec = 1e-6
      )
    }
  }

  # Respect CRAN core limits and testing environments
  if (is.null(inla_args$num.threads)) {
    chk_cores <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
    if ((nzchar(chk_cores) && (chk_cores != "false")) || nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_", ""))) {
      inla_args$num.threads <- "1:1"
    }
  }

  # Gaussian with known sampling variances (Fay-Herriot)
  if (family == "gaussian" && !is.null(vardir)) {
    if (any(vardir[!is.na(y)] <= 0, na.rm = TRUE)) {
      cli::cli_abort("{.arg vardir} must be strictly positive for sampled domains.")
    }
    inla_args$scale <- 1 / vardir
    inla_args$control.family <- list(
      hyper = list(prec = list(initial = 0, fixed = TRUE))
    )
  }

  # Beta with area-specific sampling variances or sample sizes (Janicki 2020)
  if (family == "beta") {
    if (!is.null(vardir)) {
      if (any(vardir[!is.na(y)] <= 0, na.rm = TRUE)) {
        cli::cli_abort("{.arg vardir} must be strictly positive for sampled domains.")
      }
      phi_dir <- (y * (1 - y) / vardir) - 1
      phi_dir[phi_dir < 1] <- 1
      inla_args$scale <- phi_dir
      inla_args$control.family <- list(
        hyper = list(theta = list(initial = 0, fixed = TRUE))
      )
    } else if (!is.null(trials)) {
      phi_trials <- trials - 1
      phi_trials[phi_trials < 1] <- 1
      inla_args$scale <- phi_trials
      inla_args$control.family <- list(
        hyper = list(theta = list(initial = 0, fixed = TRUE))
      )
    }
  }

  # Gamma with area-specific sampling variances or CV (Gamma SAE)
  if (family == "gamma" && !is.null(vardir)) {
    if (any(vardir[!is.na(y)] <= 0, na.rm = TRUE)) {
      cli::cli_abort("{.arg vardir} must be strictly positive for sampled domains.")
    }
    s_gamma <- rep(1, n_obs)
    idx_valid <- which(!is.na(y) & !is.na(vardir) & vardir > 0)
    s_gamma[idx_valid] <- (y[idx_valid]^2) / vardir[idx_valid]
    s_gamma[s_gamma <= 0 | is.na(s_gamma)] <- 1
    inla_args$scale <- s_gamma
    inla_args$control.family <- list(
      hyper = list(theta = list(initial = 0, fixed = TRUE))
    )
  }

  # Binomial with trials
  if (family == "binomial" && !is.null(trials)) {
    inla_args$Ntrials <- trials
  }

  # Poisson / Negative Binomial with exposure/offset
  if (family %in% c("poisson", "nbinomial") && !is.null(exposure)) {
    inla_args$E <- exposure
  }

  # Link function of the likelihood (canonical when link = NULL)
  link <- .check_inla_link(link, family)
  inla_args$control.family$control.link <- list(model = link)

  # Run INLA fitting
  fit <- tryCatch(
    do.call(INLA::inla, inla_args),
    error = function(e) {
      cli::cli_abort(c(
        "Fitting model with INLA failed.",
        "x" = e$message
      ))
    }
  )

  # Extract formatted results
  res <- .extract_inla_results(
    fit = fit,
    data = data,
    y = y,
    domain = domain,
    time = time,
    family = family,
    spatial = spatial,
    temporal = temporal,
    st_interaction = st_interaction,
    vardir = vardir,
    trials = trials,
    exposure = exposure,
    X_mat = X_mat,
    rho_range = rho_range,
    domain_id = domain_id,
    time_id = time_id,
    unique_domains = unique_domains,
    unique_times = unique_times,
    call = call
  )

  if (!is.null(sb_info)) {
    if (identical(sb_info$method, "posterior")) {
      res <- .posterior_benchmark(res, sb_info, family, n_draws = benchmark_draws, seed = benchmark_seed, link = link)
    }
    res$self_benchmark <- TRUE
    res$benchmark_method <- sb_info$method
    res$benchmark_weight <- benchmark_weight
    res$benchmark_target <- benchmark_target
    res$benchmark_group <- benchmark_group
    res$benchmark_summary <- .compute_benchmark_summary(
      unique_groups = sb_info$unique_groups,
      group_vec = sb_info$group_vec,
      w_norm = sb_info$w_norm,
      y = sb_info$original_y,
      target_vec = sb_info$target_vec,
      estimates = res$df_hb$hb
    )
    if (identical(sb_info$method, "wfq")) {
      res$df_hb$y <- sb_info$original_y
      res$hb$y <- sb_info$original_y
      res$df_eblup$y <- sb_info$original_y
    }
  }

  res$link <- link

  if (inherits(data, "sf") && requireNamespace("sf", quietly = TRUE)) {
    if (!is.null(res$df_hb) && nrow(res$df_hb) == nrow(data)) {
      res$df_hb <- sf::st_set_geometry(res$df_hb, sf::st_geometry(data))
    }
  }

  return(res)
}
