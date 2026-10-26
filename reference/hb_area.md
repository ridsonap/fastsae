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

  - `"binomial"`: Binary / proportion response. Requires `trials`.

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
    effective sample sizes for complex designs); with both, `trials` is
    ignored (a warning). Without either, one precision is estimated,
    which is weakly identified against the area effect.

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
  per area (for `family = "binomial"`).

- exposure:

  Vector, column name, or formula specifying expected counts or
  population exposure offsets (for `family = "poisson"` or
  `"nbinomial"`).

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
  (if specified), observed `y`, predicted `hb`, linear predictor
  `linear_pred`, posterior standard deviation `sd`, `mse` (the posterior
  variance \\sd^2\\, named as in the EBLUP results; it is not a
  frequentist MSE), relative error `rse` (%), 95% credible interval
  (`ci_lower`, `ci_upper`), and `random_effect`. With `vardir`, also
  `vardir` and, for the beta and gamma families, the `precision` of each
  area the model was fitted with.

- `estcoef`: Data frame of estimated regression coefficients.

- `hyperpar`: Data frame of hyperparameter posterior estimates.

- `random_effect_var`: Posterior mean of the area random effect variance
  (see the Engines section on its accuracy when the variance is small).

- `random_effect_var_time`: Estimated temporal random effect variance
  (if temporal).

- `phi`: Estimated spatial variance proportion (for BYM2) or spatial
  autocorrelation.

- `rho_time`: Estimated temporal autocorrelation parameter (for AR1).

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

7.  Riebler, A., Sørbye, S. H., Simpson, D., and Rue, H. (2016). An
    intuitive Bayesian spatial model for disease mapping that accounts
    for scaling. *Statistical Methods in Medical Research*, 25(4),
    1145-1165.

8.  Marhuenda, Y., Molina, I., and Morales, D. (2013). Small area
    estimation with spatio-temporal Fay-Herriot models. *Computational
    Statistics & Data Analysis*, 58, 308-325.

9.  De Nicolò, S., and Gardini, A. (2024). The R Package tipsae: Tools
    for Mapping Proportions and Indicators on the Unit Interval.
    *Journal of Statistical Software*, 108(1), 1-36.

10. Phan, D., Pradhan, N., and Jankowiak, M. (2019). Composable effects
    for flexible and accelerated probabilistic programming in NumPyro.
    *arXiv:1912.11554*.

11. Piironen, J., and Vehtari, A. (2017). Sparsity information and
    regularization in the horseshoe and other shrinkage priors.
    *Electronic Journal of Statistics*, 11(2), 5018-5051.

12. Rue, H., Martino, S., and Chopin, N. (2009). Approximate Bayesian
    inference for latent Gaussian models by using integrated nested
    Laplace approximations. *Journal of the Royal Statistical Society:
    Series B*, 71(2), 319-392.

13. Knorr-Held, L. (2000). Bayesian modelling of inseparable space-time
    variation in disease risk. *Statistics in Medicine*, 19(17-18),
    2555-2567.

14. van Niekerk, J., Krainski, E., Rustand, D., and Rue, H. (2023). A
    new avenue for Bayesian inference with INLA. *Computational
    Statistics & Data Analysis*, 181, 107692.

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
#> Random effect variance (sigma2_u): 2.880417 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  3.11016649  0.79793357  3.89777622  0.00009708  1.55566885
#> x1          -0.00198357  0.01028857 -0.19279357  0.84712064 -0.02259571
#> x2           0.05569461  0.06348802  0.87724594  0.38035307 -0.06834074
#> x3           0.03345853  0.05993420  0.55825429  0.57667075 -0.08590249
#>             ci_upper
#> (Intercept)   4.7053
#> x1            0.0181
#> x2            0.1821
#> x3            0.1504
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred        sd       mse      rse ci_lower
#> 1      1 8.359527 7.608034    7.608034 0.7646552 0.5846976 10.05063 6.123792
#> 2      2 7.599650 6.779532    6.779532 0.8360710 0.6990147 12.33228 5.161197
#> 3      3 5.514137 5.183059    5.183059 0.8343357 0.6961161 16.09736 3.557405
#> 4      4 3.869326 4.201226    4.201226 0.7385935 0.5455204 17.58043 2.744929
#> 5      5 6.305063 6.322898    6.322898 0.9421046 0.8875611 14.89989 4.473098
#> 6      6 3.926807 4.048399    4.048399 0.5872603 0.3448746 14.50599 2.894943
#>   ci_upper random_effect    vardir
#> 1 9.121603    2.93893449 0.6618838
#> 2 8.439307    2.54244887 0.8374691
#> 3 6.830947    0.96059492 0.8822257
#> 4 5.642493   -1.31947626 0.6581716
#> 5 8.172683   -0.03983681 1.2788021
#> 6 5.198417   -0.82110983 0.3878004
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
#> Random effect variance (sigma2_u): 2.880409 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  3.1102e+00  7.9793e-01  3.8978e+00  9.7077e-05  1.5557e+00
#> x1          -1.9836e-03  1.0289e-02 -1.9279e-01  8.4712e-01 -2.2596e-02
#> x2           5.5695e-02  6.3488e-02  8.7725e-01  3.8035e-01 -6.8341e-02
#> x3           3.3459e-02  5.9934e-02  5.5826e-01  5.7667e-01 -8.5902e-02
#>             ci_upper
#> (Intercept)   4.7053
#> x1            0.0181
#> x2            0.1821
#> x3            0.1504
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred        sd       mse      rse ci_lower
#> 1      1 8.359527 7.548517    7.608032 0.8084043 0.6535176 10.70945 5.989925
#> 2      2 7.599650 6.731801    6.779530 0.8638367 0.7462138 12.83218 5.004285
#> 3      3 5.514137 5.136292    5.183058 0.8620172 0.7430736 16.78287 3.459030
#> 4      4 3.869326 4.132899    4.201227 0.7687101 0.5909153 18.59978 2.686835
#> 5      5 6.305063 6.336474    6.322898 0.9149853 0.8371981 14.43998 4.555482
#> 6      6 3.926807 4.019001    4.048399 0.6802954 0.4628018 16.92698 2.745402
#>   ci_upper random_effect    vardir hb_unbenchmarked
#> 1 9.145367    2.93893253 0.6618838         7.608032
#> 2 8.421588    2.54244704 0.8374691         6.779530
#> 3 6.892511    0.96059462 0.8822257         5.183058
#> 4 5.634806   -1.31947499 0.6581716         4.201227
#> 5 8.204404   -0.03983666 1.2788021         6.322898
#> 6 5.352775   -0.82110901 0.3878004         4.048399
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
#> Random effect variance (sigma2_u): 3.010341 
#> Spatial autocorrelation (rho): 0.1571 
#> Spatial mixing fraction (phi): 0.1571 
#> 
#> Fixed Effects Coefficients:
#>                   beta  std.error     zvalue     pvalue   ci_lower ci_upper
#> (Intercept)  3.0843539  0.7981641  3.8643107  0.0001114  1.5278934   4.6781
#> x1          -0.0011155  0.0104187 -0.1070656  0.9147369 -0.0219140   0.0193
#> x2           0.0503857  0.0642349  0.7843970  0.4328072 -0.0754561   0.1780
#> x3           0.0382391  0.0605731  0.6312888  0.5278517 -0.0821454   0.1567
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred        sd       mse      rse ci_lower
#> 1      1 8.359527 7.577621    7.577621 0.7693028 0.5918268 10.15230 6.082931
#> 2      2 7.599650 6.775971    6.775971 0.8380187 0.7022753 12.36751 5.152607
#> 3      3 5.514137 5.157678    5.157678 0.8350475 0.6973044 16.19038 3.531165
#> 4      4 3.869326 4.202490    4.202490 0.7374601 0.5438474 17.54817 2.748474
#> 5      5 6.305063 6.286797    6.286797 0.9467412 0.8963189 15.05920 4.428984
#> 6      6 3.926807 4.070945    4.070945 0.5896216 0.3476536 14.48365 2.912893
#>   ci_upper random_effect    vardir
#> 1 9.098953    2.88719620 0.6618838
#> 2 8.438358    2.53997834 0.8374691
#> 3 6.807361    0.98850805 0.8822257
#> 4 5.641590   -1.27381628 0.6581716
#> 5 8.146456   -0.06694308 1.2788021
#> 6 5.225583   -0.74994440 0.3878004
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
