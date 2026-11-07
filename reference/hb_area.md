# Hierarchical Bayes for Area-Level Small Area Estimation

Estimates small area parameters using area-level Hierarchical Bayesian
models under various distributions (Gaussian, Binomial, Poisson,
Negative Binomial, Beta, Gamma) and spatial random effect structures
(non-spatial, BYM2, BYM, Besag/ICAR, Leroux, SAR). Two estimation
engines are available through `engine`:

- [`inla_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)
  (default, or `"inla"`): Integrated Nested Laplace Approximations
  (INLA; Rue et al., 2009).

- [`jax_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)
  (or `"jax"`): full MCMC with the No-U-Turn sampler of 'NumPyro' (Phan
  et al., 2019) on 'JAX', on the CPU, an NVIDIA GPU or the Apple silicon
  GPU. Needs Python with `jax` and `numpyro` through reticulate; see
  [`setup_numpyro_env`](https://ridsonap.is-a.dev/fastsae/reference/setup_numpyro_env.md).

The engine constructors hold all engine settings (priors, approximation
strategy and further
[`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html) arguments for
INLA; MCMC length, device and coefficient prior for JAX), so `hb_area()`
itself only describes the model.

## Usage

``` r
hb_area(
  formula,
  data,
  domain = NULL,
  time = NULL,
  family = c("gaussian", "binomial", "poisson", "nbinomial", "beta", "gamma"),
  spatial = c("none", "bym2", "bym", "besag", "generic1", "leroux", "sar", "slm"),
  temporal = c("none", "rw1", "rw2", "ar1", "iid"),
  st_interaction = c("none", "domain-specific", "separable", "type1", "type2", "type3",
    "type4"),
  W = NULL,
  vardir = NULL,
  trials = NULL,
  exposure = NULL,
  var_x = NULL,
  me_model = c("structural", "functional"),
  engine = "inla",
  self_benchmark = FALSE,
  benchmark_weight = NULL,
  benchmark_target = NULL,
  benchmark_group = NULL,
  benchmark_method = c("posterior", "wfq"),
  print_result = TRUE,
  method = NULL,
  ...
)
```

## Arguments

- formula:

  An object of class `formula` specifying the fixed-effects model (e.g.,
  `y ~ x1 + x2`).

- data:

  A `data.frame` containing the area-level data (one row per area, or
  per area and period). Areas without a direct estimate (unsampled
  areas) are kept as rows with `NA` response and their covariates; they
  are predicted from the model. An area with no row in `data` has no
  covariates and is not predicted.

- domain:

  Vector, column name, or one-sided formula referencing the area/domain
  identifier in `data`. If `NULL`, domains are numbered consecutively.

- time:

  Vector, column name, or one-sided formula referencing the time period
  identifier in `data` (e.g. `time = "year"` or `~year`). Required when
  `temporal != "none"`.

- family:

  Character string specifying the response distribution / likelihood.
  Options:

  - `"gaussian"`: Continuous response (Fay-Herriot model). If `vardir`
    is provided, sampling variances are treated as known. If
    `vardir = NULL`, residual variance is estimated, with a
    half-normal(sd y) prior on its SD; it is then only weakly identified
    against the area effect.

  - `"binomial"`: Counts or proportions. Requires `trials`. Proportions
    (all in \[0, 1\], at least one not a whole number) are turned into
    numbers of cases \\y_i n_i\\, without rounding. For survey-weighted
    proportions give effective sample sizes \\n_i = n_i^{obs} /
    \mathrm{deff}\_i\\ in `trials`: then \\y_i n_i\\ is the effective
    number of cases and the binomial variance \\p_i(1 - p_i)/n_i\\
    matches the design variance of the direct estimate (Chen, Wakefield
    and Lumley 2014; Franco and Bell 2015). Direct estimates of 0 or 1
    are kept as such. Non-integer cases or trials use INLA's
    `"xbinomial"` likelihood, which equals the binomial one on whole
    numbers. The design effects are not estimated by the model: the
    intervals are only as good as the effective sample sizes.

  - `"poisson"`: Count response / disease rate. Supports `exposure`.

  - `"nbinomial"`: Negative Binomial for overdispersed counts. If
    `vardir` is provided, the size of each area is fixed by \\V_i =
    \mu_i + \mu_i^2/n_i\\ at \\n_i = \hat\mu_i^2/(V_i - \hat\mu_i)\\ (a
    near-Poisson size when \\V_i \le \hat\mu_i\\), with \\\hat\mu_i\\
    the posterior mean of a first fit that estimates one size (so the
    model is fitted twice); without `vardir` that common size is weakly
    identified against the area effect (a warning).

  - `"beta"`: Continuous proportions strictly in (0, 1). If `vardir` is
    provided, the precision of each area is fixed by \\V_i = \mu_i(1 -
    \mu_i)/(\phi_i + 1)\\ at \\\phi_i = \hat\mu_i(1 - \hat\mu_i)/V_i -
    1\\ (at least 1), with \\\hat\mu_i\\ the posterior mean of a first
    fit with the precision at the direct estimate \\y_i\\ (Janicki
    2020), so the model is fitted twice. If `trials` is provided
    instead, the precision is \\n_i - 1\\ (simple random sampling; give
    effective sample sizes for complex designs, which gives the
    beta-logistic model of Liu, Lahiri and Kalton 2014); with both,
    `trials` is ignored (a warning). Without either, one precision is
    estimated, which is weakly identified against the area effect.

  - `"gamma"`: Skewed positive continuous response. If `vardir` is
    provided, the precision of each area is scaled as \\y_i^2 / V_i\\;
    without it the precision is estimated from one observation per area,
    which is weakly identified against the area effect (a warning; with
    INLA the 95% intervals covered 83-88% of the true means in a
    simulation, 91-95% with `vardir`).

- spatial:

  Character string specifying the spatial random effect structure:

  - `"none"`: Non-spatial IID area random intercept (default).

  - `"bym2"`: Scaled Besag-York-Mollié 2 model (Riebler et al., 2016),
    decomposing area variance into spatial (\\\phi\\) and unstructured
    components.

  - `"bym"`: Classic Besag-York-Mollié (1991) model.

  - `"besag"`: Intrinsic Conditional Autoregressive (ICAR) model.

  - `"generic1"`: Leroux spatial autoregressive model.

  - `"sar"`: Simultaneous autoregressive (SAR) area effect, \\\eta =
    X\beta + u\\, \\u = (I - \rho W)^{-1}\epsilon\\, \\\epsilon \sim
    N(0, \sigma^2 I)\\: the model of
    [`eblup_sfh`](https://ridsonap.is-a.dev/fastsae/reference/eblup_sfh.md)
    (Petrucci and Salvati, 2006), with \\\rho\\ uniform on \\\|\rho\|
    \le 0.999\\. With INLA it is the `"slm"` latent model with a zero
    covariate matrix, so that only \\\epsilon\\ is filtered. Not
    combinable with the `"type3"`, `"type4"` and `"separable"`
    interactions, which are CAR fields.

  - `"slm"`: Spatial lag model, \\\eta = (I - \rho W)^{-1}(X\beta +
    \epsilon)\\ (INLA `"slm"`). Note that the covariates are spatially
    filtered too, so this is *not* the Bayesian counterpart of
    [`eblup_sfh`](https://ridsonap.is-a.dev/fastsae/reference/eblup_sfh.md),
    whose SAR structure applies to the random effect only: \\X\beta +
    (I - \rho W)^{-1}\epsilon\\. Use `"sar"` for that model.

- temporal:

  Character string specifying the temporal random effect structure:

  - `"none"`: No temporal random effect (default cross-sectional model).

  - `"rw1"`: First-order random walk across time periods.

  - `"rw2"`: Second-order random walk across time periods.

  - `"ar1"`: First-order autoregressive process across time periods.

  - `"iid"`: Unstructured independent time effects.

  With `engine = jax_engine()`, the `"rw1"` and `"ar1"` effects sum to
  zero over the periods (`"ar1"` only when the model has an intercept),
  so the intercept is the mean level: over a few periods the level of an
  AR(1) effect otherwise trades off against the intercept and NUTS
  diverges. INLA fits `"ar1"` without this constraint, so the intercepts
  of the two engines differ by the mean of the time effects.

- st_interaction:

  Character string specifying the spatio-temporal structure:

  - `"none"`: Additive main spatial and temporal effects (default).

  - `"domain-specific"`: Domain-specific temporal random walk / AR(1)
    dynamics with shared variance parameter (directly corresponds to the
    tipsae spatio-temporal model). Not available with
    `temporal = "rw2"`: without a common trend, the linear trend of a
    domain without direct estimates has a flat prior; use `"type2"`.

  - `"separable"`: Grouped dynamic spatial field evolving over time via
    AR(1) (corresponds to classical Spatio-Temporal Fay-Herriot models,
    Marhuenda et al. 2013). Only `temporal = "ar1"` or `"iid"`: with a
    random walk in time every domain's level in the field has a flat
    prior, as the field has no spatial main effect; use `"type4"`.

  - `"type1"`: Knorr-Held (2000) Type I interaction (unstructured space
    \\\times\\ unstructured time).

  - `"type2"`: Knorr-Held Type II interaction (unstructured space
    \\\times\\ structured time): an independent `temporal` process for
    each domain.

  - `"type3"`: Knorr-Held Type III interaction (structured space
    \\\times\\ unstructured time): an independent ICAR field for each
    period.

  - `"type4"`: Knorr-Held Type IV interaction (structured space
    \\\times\\ structured time): an ICAR field evolving over time with
    the `temporal` structure.

  The `"type*"` options add the interaction to the main spatial and
  temporal effects. `"type3"` and `"type4"` need a structured `spatial`
  effect (`"bym2"`, `"bym"`, `"besag"` or `"leroux"`); with INLA the
  interaction is an ICAR (`"besag"`) field on the graph of `W`. With
  `temporal = "rw1"` or `"rw2"`, the type IV interaction (and type II
  with `"rw2"`) uses the scaled Kronecker structure with all the
  identifiability constraints of Knorr-Held (2000): sums over the
  domains are 0 in every period, and sums over the periods (and, for
  `"rw2"`, linear trends) are 0 in every domain. INLA's default
  constraints leave the latter out; with them alone, a domain without
  direct estimates had a posterior SD of 36-58 (other domains about
  0.4).

- W:

  Proximity or spatial adjacency matrix. Can be a square `matrix`,
  `Matrix`, or `spdep` `nb` or `listw` object. Dimensions must match the
  total number of domains in `data`. Required when `spatial != "none"`.
  The CAR models (`"besag"`, `"bym"`, `"bym2"`, `"leroux"`) use only the
  neighbourhood pattern `W > 0`, not the weights: give a neighbourhood
  graph (for example `create_weights(method = "queen")` or `"knn"`). A
  dense matrix such as distance-decay weights becomes a nearly complete
  graph and gives a warning. `"sar"` and `"slm"` use the weights,
  row-standardised (as
  [`eblup_sfh`](https://ridsonap.is-a.dev/fastsae/reference/eblup_sfh.md)
  expects).

- vardir:

  Vector, column name, or formula specifying the sampling variances of
  the direct estimator (for `family = "gaussian"`, `"beta"`, `"gamma"`
  or `"nbinomial"`).

- trials:

  Vector, column name, or formula specifying sample sizes / total trials
  per area (for `family = "binomial"`); effective sample sizes (sample
  size divided by the design effect) need not be whole numbers.

- exposure:

  Vector, column name, or formula specifying expected counts or
  population exposure offsets (for `family = "poisson"` or
  `"nbinomial"`).

- var_x:

  Optional variances of the errors of covariates measured with error,
  for instance the squared standard errors of covariates that are
  themselves survey estimates: a vector or list named by the covariates,
  each element a column name, one-sided formula or numeric vector (one
  value per row of `data`), e.g. `var_x = c(x1 = "var_x1")`. The
  covariates must enter `formula` as main effects, not transformed or in
  interactions, and must not be constant. A variance of 0 marks a value
  without error. See the section Covariates measured with error.

- me_model:

  Model of the true values of the covariates in `var_x`: `"structural"`
  (default) or `"functional"`; see the section Covariates measured with
  error.

- engine:

  Estimation engine:
  [`inla_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)
  or
  [`jax_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md),
  or the strings `"inla"` (default) and `"jax"` for their default
  settings. See the Engines section.

- self_benchmark:

  Logical. If `TRUE`, the weighted aggregate of the area estimates is
  made to match the direct survey aggregate or `benchmark_target`; see
  `benchmark_method`.

- benchmark_weight:

  Optional vector, column name, or one-sided formula referencing the
  domain benchmarking weights in `data` (e.g. population sizes or
  shares). Defaults to equal weights.

- benchmark_target:

  Optional numeric scalar or named vector specifying fixed aggregate
  benchmark target(s). If `NULL` (default), calibrates to the weighted
  aggregate of the direct estimates.

- benchmark_group:

  Optional vector, column name, or one-sided formula referencing the
  grouping / stratum column in `data` (e.g. province) for group-specific
  self-benchmarking.

- benchmark_method:

  How self-benchmarking is done when `self_benchmark = TRUE`:
  `"posterior"` (default) benchmarks every posterior draw of the area
  parameters (an additive shift for `"gaussian"`, a common shift on the
  logit scale for `"binomial"` and `"beta"`, a common ratio for
  `"poisson"`, `"nbinomial"` and `"gamma"`), so the constraint holds
  exactly for the posterior mean and the reported SD / credible
  intervals refer to the benchmarked estimator (Bayesian benchmarking,
  cf. Datta et al., 2011). `"wfq"` uses the augmented-covariate approach
  of Wang, Fuller and Qu (2008), which is exact only for the linear
  Gaussian EBLUP; for the HB models it holds only approximately
  (non-flat prior on the coefficients, integration over the
  hyperparameters) and not at all for non-Gaussian likelihoods. With
  INLA, posterior draws need the package sn (used by
  [`INLA::inla.posterior.sample()`](https://rdrr.io/pkg/INLA/man/posterior.sample.html));
  without it `"wfq"` is used with a warning.

- print_result:

  Logical. If `TRUE` (default), prints a summary of results.

- method:

  Deprecated; use `engine`. `method = "laplace"` becomes
  `engine = inla_engine(strategy = "laplace")`; the former lme4 GLMM
  route for binomial and Poisson models has been removed.

- ...:

  Deprecated. `strategy`, `link`, `scale_model`, `prior_*`,
  `benchmark_draws`, `benchmark_seed` and arguments of
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html) are now given
  to
  [`inla_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md);
  passed here they are forwarded to it with a warning.

## Value

An object of class `c("fastsae_hb_area", "fastsae_hb", "fastsae")`
containing:

- `df_hb`: Data frame with domain estimates, including `domain`, `time`
  (if specified), the direct estimate `y` on the scale of `hb` (a
  proportion for the binomial family, a rate per unit exposure for
  counts with `exposure`; the fitted numbers of cases or counts are in
  `cases` or `count`), predicted `hb`, linear predictor `linear_pred`,
  posterior standard deviation `sd`, `mse` (the posterior variance
  \\sd^2\\, named as in the EBLUP results; it is not a frequentist MSE),
  relative error `rse` (%), 95% credible interval (`ci_lower`,
  `ci_upper`), and `random_effect`. With `vardir`, also `vardir` and,
  for the beta and gamma families, the `precision` of each area the
  model was fitted with.

- `estcoef`: Data frame of estimated regression coefficients.

- `hyperpar`: Data frame of hyperparameter posterior estimates.

- `random_effect_var`: Posterior mean of the area random effect variance
  (see the Engines section on its accuracy when the variance is small).

- `random_effect_var_time`: Estimated temporal random effect variance
  (if temporal).

- `phi`: Estimated spatial variance proportion (for BYM2) or spatial
  autocorrelation.

- `rho_time`: Estimated temporal autocorrelation parameter (for AR1).

- `measurement_error`: With `var_x`, a list with the `model`, the
  `covariates` measured with error and, for the structural model, the
  mean `mu_x` and variance `sigma2_x` of their true values; `NULL`
  otherwise. The coefficients of these covariates are in `estcoef` with
  the others.

- `goodness`: Model fit metrics (DIC, WAIC, Marginal Log-Likelihood).

- `family`: Response family used.

- `spatial`: Spatial model type used.

- `temporal`: Temporal model type used.

- `st_interaction`: Spatio-temporal interaction type used.

- `self_benchmark`: Logical indicating whether self-benchmarking was
  active.

- `benchmark_summary`: Summary table of benchmark calibration (if
  self-benchmarking).

- `fit`: Raw fitted model object (`inla`, or the list returned by the
  Python backend for the JAX engine).

- `engine`: Name of the estimation engine (`"inla"` or `"jax"`);
  `engine_spec` holds its settings
  ([`inla_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)
  or
  [`jax_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)
  object).

- `draws`, `rhat_max`, `device`, `elapsed_seconds`: With the JAX engine,
  the posterior draws of the area parameters (rows of `df_hb` by draws),
  the largest split R-hat (`convergence` is `rhat_max < 1.05`), the JAX
  device and the sampling time.

- `call`: Matched function call.

## Engines

Both engines return the same structure, so `print`, `summary`,
[`benchmark_sae`](https://ridsonap.is-a.dev/fastsae/reference/benchmark_sae.md),
[`compare_sae`](https://ridsonap.is-a.dev/fastsae/reference/compare_sae.md)
and the other methods work with either. Feature support:

|  |  |  |
|----|----|----|
| **Feature** | **inla** | **jax** |
| Families | all six | all six |
| `spatial` | all | all but `"slm"` |
| `temporal` | all | all but `"rw2"` |
| `st_interaction` | all | all |
| Self-benchmarking | posterior, wfq | posterior; wfq (gaussian) |
| Covariates measured with error (`var_x`) | yes | yes |
| Horseshoe prior, Student-t random effects | – | yes |

By default both engines use the same priors (see `prior_prec` of
[`inla_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)
and `priors` of
[`jax_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)),
so their estimates differ only by the INLA approximation and the Monte
Carlo error.

With the `"inla"` engine, Gaussian models and the Leroux and SLM effects
of any family are fitted in INLA's `"classic"` mode (Rue, Martino and
Chopin, 2009); the other models use INLA's default `"compact"` mode (van
Niekerk et al., 2023), which leaves the linear predictors out of the
latent field and is designed for data-rich models. Gaussian models with
a `"type1"` to `"type4"` interaction and three hyperparameters (a Besag
or iid area effect and a random walk or iid in time) also use a grid
over the hyperparameters (`int.strategy = "grid"`); with
`temporal = "ar1"` or a BYM, BYM2, Leroux or SLM area effect the grid
grows exponentially (13 minutes for 12 areas with AR(1)) and was not
validated, so INLA's default integration is kept. These choices rest on
comparisons with long Stan runs with the same model and priors
(Gaussian, binomial, Poisson, beta, gamma and negative binomial; iid,
Besag, BYM2, Leroux and SLM effects; RW1 in time and the four Knorr-Held
interactions):

- Gaussian: compact overstated the posterior SDs by up to 18% when the
  area variance was small relative to `vardir`; classic matched Stan
  (means within 0.07 posterior SD, SDs within 0.91-1.03 times) at the
  same speed. With an interaction both modes overstated the SDs by
  6-13%; classic with the grid matched Stan (means within 0.04 SD, SDs
  within 2%) and took 5-8 times longer.

- Binomial, Poisson, beta and gamma with iid, Besag or BYM2 effects or
  an interaction: compact matched Stan (means within 0.16 SD, SDs within
  0.95-1.08 times), and classic was much slower for many areas (beta,
  2000 areas: 32 s against 1.3 s).

- Leroux and SLM: compact moved the Poisson posterior means by 0.12-0.15
  SD; classic was within 0.03 SD.

- Negative binomial: with one count per area the overdispersion and the
  area effect are only weakly identified, and both modes were off (means
  up to 0.64 SD); the model gives a warning.

`inla_engine(inla.mode = , control.inla = )` overrides these settings.
For Gaussian models with a small area variance, the classic summary of
that variance (`random_effect_var`) was biased downwards, while the area
estimates were accurate. The validation scripts are in the package
repository (`benchmarks/validation`).

## Covariates measured with error

A covariate that is itself an estimate, for instance from another
survey, is observed as \\w\_{dk} = x\_{dk} + \eta\_{dk}\\, where
\\x\_{dk}\\ is its true value and \\\eta\_{dk} \sim N(0, C\_{dk})\\ an
error with known variance \\C\_{dk}\\ (`var_x`), independent of the
sampling error of \\y_d\\, of the area effect and of the errors of the
other covariates. Used as if it were the true value, \\w\_{dk}\\ gets an
attenuated coefficient and every area is shrunk alike, however large its
\\C\_{dk}\\, and the predictor can be worse than the direct estimate
(Ybarra and Lohr, 2008). With `var_x` the true values replace the
observed ones in the linear predictor:

- `"structural"` (default): \\x\_{dk} \sim N(\mu_k, \sigma_k^2)\\
  independently, the structural model of Bell et al. (2019) with a
  common mean, with \\\mu_k \sim N(\bar w_k, (100 s_k)^2)\\ and a
  half-normal(\\s_k\\) prior on \\\sigma_k\\, where \\\bar w_k\\ and
  \\s_k\\ are the mean and SD of the observed values. The coefficient is
  that of the true covariate.

- `"functional"`: a flat prior on the true values (Arima, Datta and
  Liseo, 2015), the Bayesian counterpart of the model of Ybarra and
  Lohr. Integrating the true values out gives \\y_d \sim N(w_d'\beta,
  \psi_d + \sigma_u^2 + \beta' C_d \beta)\\, the regression on the
  observed values with a larger variance: the area estimates agree with
  the Ybarra-Lohr EBLUP (`emdi::fh(method = "me")`), but the coefficient
  keeps the attenuation that the moment estimator of Ybarra and Lohr
  removes.

The coefficients of these covariates get the prior of the other
coefficients (with INLA `control.fixed`, \\N(0, 1000)\\ by default;
`prior_beta` of
[`jax_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)).
A variance of 0 marks a value without error. Rows of the same domain
with the same observed value and variance share one true value (an
estimate used in several periods); every other row has its own. Not
available with `spatial = "slm"`, whose latent effect contains the
covariates.

In a simulation with 60 areas (10 without a direct estimate) and 100
replicates per scenario (moderate and heavy errors, skewed true values,
and a mix of exact and very noisy covariates;
`benchmarks/validation/me_area_coverage.R`), the 95% intervals of the
structural model covered 94-97% of the true area means in every scenario
and its RMSE was the smallest, also against the Ybarra-Lohr EBLUP; the
naive model covered 88.7% of the areas with noisy covariates in the
mixed scenario (96.5% of those without error). For a true coefficient of
0.8, the structural estimates averaged 0.80-0.86 (their intervals
covered 85-97%), the naive and functional ones 0.36-0.67 and those of
Ybarra-Lohr 0.79-0.95.

With INLA, each covariate is a classical measurement error term
(`"mec"`; Muff et al., 2015) whose coefficient is a hyperparameter, and
the model is fitted in INLA's classic mode (compact moved the posterior
means of binomial models by up to 0.08 posterior SD against Stan,
classic 0.02). For the structural model, \\\mu_k\\ and \\\sigma_k\\ are
fixed at their posterior mode given the observed values, which \\y\\
hardly informs: against Stan, which integrates them out, the posterior
means of the areas were within 0.03 posterior SD and the SDs within
0.97-1.02 times with one such covariate (Gaussian and binomial models,
30-100 areas); with two, the means were within 0.11 SD and the posterior
SDs of the coefficients 0.87-0.98 times those of Stan. As two more
hyperparameters, \\\mu_k\\ and \\\sigma_k\\ moved the means by up to
0.22 SD. The JAX engine samples them.

## References

1.  Rao, J. N. K., and Molina, I. (2015). *Small Area Estimation*. John
    Wiley & Sons.

2.  Besag, J., York, J., and Mollie, A. (1991). Bayesian image
    restoration, with two applications in spatial statistics. *Annals of
    the Institute of Statistical Mathematics*, 43(1), 1-20.

3.  Janicki, R. (2020). Properties of the beta regression model for
    small area estimation of proportions and application to estimation
    of poverty rates. *Communications in Statistics - Theory and
    Methods*, 49(9), 2264-2284.

4.  Leroux, B. G., Lei, X., and Breslow, N. (2000). Estimation of
    disease rates in small areas: a new mixed model for spatial
    dependence. In *Statistical Models in Epidemiology, the Environment,
    and Clinical Trials*, 135-178. Springer.

5.  Petrucci, A., and Salvati, N. (2006). Small area estimation for
    spatial correlation in watershed erosion assessment. *Journal of
    Agricultural, Biological, and Environmental Statistics*, 11(2),
    169-182.

6.  Wang, J., Fuller, W. A., and Qu, Y. (2008). Small area estimation
    under a restriction. *Survey Methodology*, 34(1), 29-36.

7.  Chen, C., Wakefield, J., and Lumley, T. (2014). The use of sampling
    weights in Bayesian hierarchical models for small area estimation.
    *Spatial and Spatio-temporal Epidemiology*, 11, 33-43.

8.  Liu, B., Lahiri, P., and Kalton, G. (2014). Hierarchical Bayes
    modeling of survey-weighted small area proportions. *Survey
    Methodology*, 40(1), 1-13.

9.  Franco, C., and Bell, W. R. (2015). Borrowing information over time
    in binomial/logit normal models for small area estimation.
    *Statistics in Transition new series*, 16(4), 563-584.

10. Riebler, A., Sørbye, S. H., Simpson, D., and Rue, H. (2016). An
    intuitive Bayesian spatial model for disease mapping that accounts
    for scaling. *Statistical Methods in Medical Research*, 25(4),
    1145-1165.

11. Marhuenda, Y., Molina, I., and Morales, D. (2013). Small area
    estimation with spatio-temporal Fay-Herriot models. *Computational
    Statistics & Data Analysis*, 58, 308-325.

12. De Nicolò, S., and Gardini, A. (2024). The R Package tipsae: Tools
    for Mapping Proportions and Indicators on the Unit Interval.
    *Journal of Statistical Software*, 108(1), 1-36.

13. Phan, D., Pradhan, N., and Jankowiak, M. (2019). Composable effects
    for flexible and accelerated probabilistic programming in NumPyro.
    *arXiv:1912.11554*.

14. Piironen, J., and Vehtari, A. (2017). Sparsity information and
    regularization in the horseshoe and other shrinkage priors.
    *Electronic Journal of Statistics*, 11(2), 5018-5051.

15. Rue, H., Martino, S., and Chopin, N. (2009). Approximate Bayesian
    inference for latent Gaussian models by using integrated nested
    Laplace approximations. *Journal of the Royal Statistical Society:
    Series B*, 71(2), 319-392.

16. Knorr-Held, L. (2000). Bayesian modelling of inseparable space-time
    variation in disease risk. *Statistics in Medicine*, 19(17-18),
    2555-2567.

17. van Niekerk, J., Krainski, E., Rustand, D., and Rue, H. (2023). A
    new avenue for Bayesian inference with INLA. *Computational
    Statistics & Data Analysis*, 181, 107692.

18. Ybarra, L. M. R., and Lohr, S. L. (2008). Small area estimation when
    auxiliary information is measured with error. *Biometrika*, 95(4),
    919-931.

19. Arima, S., Datta, G. S., and Liseo, B. (2015). Bayesian estimators
    for small area models when auxiliary information is measured with
    error. *Scandinavian Journal of Statistics*, 42(2), 518-529.

20. Bell, W. R., Chung, H. C., Datta, G. S., and Franco, C. (2019).
    Measurement error in small area estimation: Functional versus
    structural versus naive models. *Survey Methodology*, 45(1), 61-80.

21. Muff, S., Riebler, A., Held, L., Rue, H., and Saner, P. (2015).
    Bayesian analysis of measurement error models using integrated
    nested Laplace approximations. *Journal of the Royal Statistical
    Society: Series C (Applied Statistics)*, 64(2), 231-252.

## Examples

``` r
# \donttest{
if (requireNamespace("INLA", quietly = TRUE)) {
  library(fastsae)

  # 1. Non-spatial Gaussian Fay-Herriot with INLA
  m_norm <- hb_area(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = "vardir",
    family = "gaussian"
  )

  # 2. Self-benchmarking Fay-Herriot with INLA (benchmarked posterior draws); the
  #    sample sizes n are used only for illustration (normally population shares)
  m_bench <- hb_area(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = "vardir",
    self_benchmark = TRUE,
    benchmark_weight = "n"
  )

  # 3. Spatial BYM2 Gaussian Fay-Herriot with INLA, faster approximation;
  #    the CAR model needs a neighbourhood graph (4 nearest centroids)
  W_knn <- create_weights(cbind(sim_area$x_coord, sim_area$y_coord), method = "knn", k = 4)
  m_bym2 <- hb_area(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = "vardir",
    family = "gaussian",
    spatial = "bym2",
    W = W_knn,
    engine = inla_engine(strategy = "simplified.laplace")
  )
}
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> hb_area(formula = y ~ x1 + x2 + x3, data = mys, family = "gaussian", vardir =
#> "vardir")
#> 
#> ✔ Convergence: Yes (INLA mode status 0)
#> Model: HB-GAUSSIAN (Non-spatial)
#> Engine: inla
#> Random effect variance (sigma2_u): 2.880427 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  3.1102e+00  7.9794e-01  3.8978e+00  9.7084e-05  1.5557e+00
#> x1          -1.9836e-03  1.0289e-02 -1.9279e-01  8.4712e-01 -2.2596e-02
#> x2           5.5695e-02  6.3488e-02  8.7724e-01  3.8035e-01 -6.8341e-02
#> x3           3.3458e-02  5.9934e-02  5.5825e-01  5.7667e-01 -8.5903e-02
#>             ci_upper
#> (Intercept)   4.7053
#> x1            0.0181
#> x2            0.1821
#> x3            0.1504
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred        sd       mse      rse ci_lower
#> 1      1 8.359527 7.608038    7.608038 0.7646542 0.5846960 10.05061 6.123799
#> 2      2 7.599650 6.779536    6.779536 0.8360702 0.6990134 12.33226 5.161203
#> 3      3 5.514137 5.183061    5.183061 0.8343359 0.6961165 16.09736 3.557406
#> 4      4 3.869326 4.201224    4.201224 0.7385938 0.5455208 17.58044 2.744928
#> 5      5 6.305063 6.322898    6.322898 0.9421055 0.8875627 14.89990 4.473097
#> 6      6 3.926807 4.048398    4.048398 0.5872604 0.3448748 14.50600 2.894942
#>   ci_upper random_effect    vardir
#> 1 9.121606     2.9389384 0.6618838
#> 2 8.439310     2.5424526 0.8374691
#> 3 6.830949     0.9605958 0.8822257
#> 4 5.642492    -1.3194786 0.6581716
#> 5 8.172684    -0.0398370 1.2788021
#> 6 5.198417    -0.8211113 0.3878004
#> ... and 36 more rows.
#> 
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> hb_area(formula = y ~ x1 + x2 + x3, data = mys, vardir = "vardir",
#> self_benchmark = TRUE, benchmark_weight = "n")
#> 
#> ✔ Convergence: Yes (INLA mode status 0)
#> Model: HB-GAUSSIAN (Non-spatial)
#> Engine: inla
#> Self-Benchmarking: Active (benchmarked posterior draws)
#> Random effect variance (sigma2_u): 2.88043 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  3.1102e+00  7.9794e-01  3.8978e+00  9.7084e-05  1.5557e+00
#> x1          -1.9836e-03  1.0289e-02 -1.9279e-01  8.4712e-01 -2.2596e-02
#> x2           5.5695e-02  6.3488e-02  8.7724e-01  3.8035e-01 -6.8341e-02
#> x3           3.3458e-02  5.9934e-02  5.5825e-01  5.7667e-01 -8.5903e-02
#>             ci_upper
#> (Intercept)   4.7053
#> x1            0.0181
#> x2            0.1821
#> x3            0.1504
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred        sd       mse      rse ci_lower
#> 1      1 8.359527 7.549761    7.608040 0.8082664 0.6532946 10.70585 5.951396
#> 2      2 7.599650 6.732247    6.779538 0.8590779 0.7380148 12.76064 5.021879
#> 3      3 5.514137 5.137559    5.183062 0.8642473 0.7469234 16.82214 3.486832
#> 4      4 3.869326 4.132869    4.201224 0.7690848 0.5914914 18.60898 2.654658
#> 5      5 6.305063 6.332665    6.322898 0.9177923 0.8423427 14.49299 4.515273
#> 6      6 3.926807 4.020022    4.048398 0.6805457 0.4631424 16.92890 2.718330
#>   ci_upper random_effect    vardir hb_unbenchmarked
#> 1 9.153156    2.93893997 0.6618838         7.608040
#> 2 8.338095    2.54245403 0.8374691         6.779538
#> 3 6.890517    0.96059616 0.8822257         5.183062
#> 4 5.626265   -1.31947946 0.6581716         4.201224
#> 5 8.100661   -0.03983705 1.2788021         6.322898
#> 6 5.352789   -0.82111186 0.3878004         4.048398
#> ... and 36 more rows.
#> 
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> hb_area(formula = y ~ x1 + x2 + x3, data = mys, family = "gaussian", spatial =
#> "bym2", W = W_knn, vardir = "vardir", engine = inla_engine(strategy =
#> "simplified.laplace"))
#> 
#> ✔ Convergence: Yes (INLA mode status 0)
#> Model: HB-GAUSSIAN (BYM2)
#> Engine: inla
#> Random effect variance (sigma2_u): 3.010629 
#> Spatial autocorrelation (rho): 0.1575 
#> Spatial mixing fraction (phi): 0.1575 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  3.08439701  0.79828816  3.86376396  0.00011165  1.52767840
#> x1          -0.00111530  0.01042069 -0.10702792  0.91476683 -0.02191738
#> x2           0.05038359  0.06424577  0.78423200  0.43290399 -0.07548008
#> x3           0.03824033  0.06058285  0.63120722  0.52790504 -0.08216187
#>             ci_upper
#> (Intercept)   4.6783
#> x1            0.0193
#> x2            0.1780
#> x3            0.1568
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred        sd       mse      rse ci_lower
#> 1      1 8.359527 7.577917    7.577917 0.7693246 0.5918603 10.15219 6.083137
#> 2      2 7.599650 6.776312    6.776312 0.8380426 0.7023153 12.36724 5.152850
#> 3      3 5.514137 5.157808    5.157808 0.8350899 0.6973751 16.19079 3.531204
#> 4      4 3.869326 4.202365    4.202365 0.7374892 0.5438904 17.54939 2.748301
#> 5      5 6.305063 6.286759    6.286759 0.9468196 0.8964673 15.06054 4.428799
#> 6      6 3.926807 4.070923    4.070923 0.5896402 0.3476756 14.48419 2.912835
#>   ci_upper random_effect    vardir
#> 1 9.099238    2.88745981 0.6618838
#> 2 8.438686    2.54028718 0.8374691
#> 3 6.807561    0.98861922 0.8822257
#> 4 5.641528   -1.27395010 0.6581716
#> 5 8.146576   -0.06698827 1.2788021
#> 6 5.225598   -0.74997802 0.3878004
#> ... and 36 more rows.
#> 

# 4. The same Fay-Herriot model by MCMC (NumPyro NUTS on JAX)
if (setup_numpyro_env(check_only = TRUE)) {
  m_jax <- hb_area(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = "vardir",
    engine = jax_engine(samples = 500, chains = 2)
  )
}
# }
```
