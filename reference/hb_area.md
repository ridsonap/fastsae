# Hierarchical Bayes for Area-Level Small Area Estimation

Estimates small area parameters using area-level Hierarchical Bayesian
models under various distributions (Gaussian, Binomial, Poisson,
Negative Binomial, Beta, Gamma) and spatial random effect structures
(non-spatial, BYM2, BYM, Besag/ICAR, Leroux). Two estimation engines are
available through `engine`:

- [`inla_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)
  (default, or `"inla"`): Integrated Nested Laplace Approximations
  (INLA; Rue et al., 2009).

- [`jax_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)
  (or `"jax"`): full MCMC with the No-U-Turn sampler of 'NumPyro' (Phan
  et al., 2019) on 'JAX', on the CPU or an NVIDIA GPU. Needs Python with
  `jax` and `numpyro` through reticulate; see
  [`setup_numpyro_env`](https://ridsonap.github.io/fastsae/reference/setup_numpyro_env.md).

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
  spatial = c("none", "bym2", "bym", "besag", "generic1", "leroux", "slm"),
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
    `vardir = NULL`, residual variance is estimated.

  - `"binomial"`: Binary / proportion response. Requires `trials`.

  - `"poisson"`: Count response / disease rate. Supports `exposure`.

  - `"nbinomial"`: Negative Binomial for overdispersed counts.

  - `"beta"`: Continuous proportions strictly in (0, 1). If `vardir` is
    provided, area-specific precision is calculated via Janicki (2020)
    formula \\\phi_i = y_i(1 - y_i)/V_i - 1\\ and injected via INLA's
    `scale` argument. If `trials` is provided, precision is scaled as
    \\n_i - 1\\.

  - `"gamma"`: Skewed positive continuous response.

- spatial:

  Character string specifying the spatial random effect structure:

  - `"none"`: Non-spatial IID area random intercept (default).

  - `"bym2"`: Scaled Besag-York-Mollié 2 model (Riebler et al., 2016),
    decomposing area variance into spatial (\\\phi\\) and unstructured
    components.

  - `"bym"`: Classic Besag-York-Mollié (1991) model.

  - `"besag"`: Intrinsic Conditional Autoregressive (ICAR) model.

  - `"generic1"`: Leroux spatial autoregressive model.

  - `"slm"`: Spatial lag model, \\\eta = (I - \rho W)^{-1}(X\beta +
    \epsilon)\\ (INLA `"slm"`). Note that the covariates are spatially
    filtered too, so this is *not* the Bayesian counterpart of
    [`eblup_sfh`](https://ridsonap.github.io/fastsae/reference/eblup_sfh.md),
    whose SAR structure applies to the random effect only: \\X\beta +
    (I - \rho W)^{-1}\epsilon\\.

- temporal:

  Character string specifying the temporal random effect structure:

  - `"none"`: No temporal random effect (default cross-sectional model).

  - `"rw1"`: First-order random walk across time periods.

  - `"rw2"`: Second-order random walk across time periods.

  - `"ar1"`: First-order autoregressive process across time periods.

  - `"iid"`: Unstructured independent time effects.

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
  graph and gives a warning. `"slm"` uses the row-standardised weights.

- vardir:

  Vector, column name, or formula specifying the sampling variances of
  the direct estimator (for `family = "gaussian"` or `family = "beta"`).

- trials:

  Vector, column name, or formula specifying sample sizes / total trials
  per area (for `family = "binomial"`).

- exposure:

  Vector, column name, or formula specifying expected counts or
  population exposure offsets (for `family = "poisson"` or
  `"nbinomial"`).

- engine:

  Estimation engine:
  [`inla_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)
  or
  [`jax_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md),
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
  [`inla_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md);
  passed here they are forwarded to it with a warning.

## Value

An object of class `c("fastsae_hb_area", "fastsae_hb", "fastsae")`
containing:

- `df_hb`: Data frame with domain estimates, including `domain`, `time`
  (if specified), observed `y`, predicted `hb`, linear predictor
  `linear_pred`, posterior standard error `sd`, `mse`, relative error
  `rse` (%), 95% credible interval (`ci_lower`, `ci_upper`), and
  `random_effect`.

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
  ([`inla_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)
  or
  [`jax_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)
  object).

- `draws`, `rhat_max`, `device`, `elapsed_seconds`: With the JAX engine,
  the posterior draws of the area parameters (rows of `df_hb` by draws),
  the largest split R-hat (`convergence` is `rhat_max < 1.05`), the JAX
  device and the sampling time.

- `call`: Matched function call.

## Engines

Both engines return the same structure, so `print`, `summary`,
[`benchmark_sae`](https://ridsonap.github.io/fastsae/reference/benchmark_sae.md),
[`compare_sae`](https://ridsonap.github.io/fastsae/reference/compare_sae.md)
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

The `"jax"` engine uses weakly informative priors (half-normal standard
deviations, uniform spatial mixing / autocorrelation parameters) instead
of the PC priors of INLA, so its estimates differ slightly. For
`family = "gamma"` it estimates a common shape and does not use
`vardir`.

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

2.  Wang, J., Fuller, W. A., and Qu, Y. (2008). Small area estimation
    under a restriction. *Survey Methodology*, 34(1), 29-36.

3.  Riebler, A., Sørbye, S. H., Simpson, D., and Rue, H. (2016). An
    intuitive Bayesian spatial model for disease mapping that accounts
    for scaling. *Statistical Methods in Medical Research*, 25(4),
    1145-1165.

4.  Marhuenda, Y., Molina, I., and Morales, D. (2013). Small area
    estimation with spatio-temporal Fay-Herriot models. *Computational
    Statistics & Data Analysis*, 58, 308-325.

5.  De Nicolò, S., and Gardini, A. (2024). The R Package tipsae: Tools
    for Mapping Proportions and Indicators on the Unit Interval.
    *Journal of Statistical Software*, 108(1), 1-36.

6.  Phan, D., Pradhan, N., and Jankowiak, M. (2019). Composable effects
    for flexible and accelerated probabilistic programming in NumPyro.
    *arXiv:1912.11554*.

7.  Piironen, J., and Vehtari, A. (2017). Sparsity information and
    regularization in the horseshoe and other shrinkage priors.
    *Electronic Journal of Statistics*, 11(2), 5018-5051.

8.  Rue, H., Martino, S., and Chopin, N. (2009). Approximate Bayesian
    inference for latent Gaussian models by using integrated nested
    Laplace approximations. *Journal of the Royal Statistical Society:
    Series B*, 71(2), 319-392.

9.  Knorr-Held, L. (2000). Bayesian modelling of inseparable space-time
    variation in disease risk. *Statistics in Medicine*, 19(17-18),
    2555-2567.

10. van Niekerk, J., Krainski, E., Rustand, D., and Rue, H. (2023). A
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

  # 2. Self-benchmarking Fay-Herriot with INLA (benchmarked posterior draws)
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
#> Random effect variance (sigma2_u): 1.884566 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  3.0427e+00  6.9228e-01  4.3952e+00  1.1066e-05  1.6969e+00
#> x1          -1.1730e-03  8.7917e-03 -1.3342e-01  8.9386e-01 -1.8740e-02
#> x2           5.1638e-02  5.5506e-02  9.3033e-01  3.5220e-01 -5.6669e-02
#> x3           3.7778e-02  5.2655e-02  7.1747e-01  4.7309e-01 -6.6888e-02
#>             ci_upper
#> (Intercept)   4.4243
#> x1            0.0159
#> x2            0.1619
#> x3            0.1404
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred        sd       mse      rse ci_lower
#> 1      1 8.359527 7.339901    7.339901 0.7461476 0.5567363 10.16564 5.898815
#> 2      2 7.599650 6.498467    6.498467 0.8065379 0.6505034 12.41120 4.946476
#> 3      3 5.514137 5.050324    5.050324 0.7981022 0.6369671 15.80299 3.500404
#> 4      4 3.869326 4.304521    4.304521 0.7105254 0.5048464 16.50649 2.900284
#> 5      5 6.305063 6.322252    6.322252 0.8809212 0.7760222 13.93366 4.591943
#> 6      6 3.926807 4.085714    4.085714 0.5728372 0.3281424 14.02049 2.959795
#>   ci_upper random_effect    vardir
#> 1 8.823250    2.67457157 0.6618838
#> 2 8.108175    2.29078118 0.8374691
#> 3 6.632282    0.90454028 0.8822257
#> 4 5.688281   -1.16017927 0.6581716
#> 5 8.053107   -0.02602763 1.2788021
#> 6 5.206918   -0.71974447 0.3878004
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
#> Random effect variance (sigma2_u): 1.884647 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  3.0428e+00  6.9232e-01  4.3951e+00  1.1074e-05  1.6968e+00
#> x1          -1.1734e-03  8.7922e-03 -1.3346e-01  8.9383e-01 -1.8742e-02
#> x2           5.1640e-02  5.5508e-02  9.3032e-01  3.5221e-01 -5.6673e-02
#> x3           3.7776e-02  5.2657e-02  7.1739e-01  4.7313e-01 -6.6894e-02
#>             ci_upper
#> (Intercept)   4.4244
#> x1            0.0159
#> x2            0.1619
#> x3            0.1404
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred        sd       mse      rse ci_lower
#> 1      1 8.359527 7.256192    7.340056 0.7836072 0.6140402 10.79915 5.775898
#> 2      2 7.599650 6.425748    6.498622 0.8350583 0.6973223 12.99550 4.823442
#> 3      3 5.514137 4.983210    5.050404 0.8253429 0.6811909 16.56248 3.372134
#> 4      4 3.869326 4.229123    4.304468 0.7325886 0.5366861 17.32247 2.857941
#> 5      5 6.305063 6.316739    6.322253 0.8528188 0.7272999 13.50093 4.654609
#> 6      6 3.926807 4.041630    4.085695 0.6628529 0.4393739 16.40063 2.808622
#>   ci_upper random_effect    vardir hb_unbenchmarked
#> 1 8.797418    2.67472105 0.6618838         7.340056
#> 2 8.110314    2.29091891 0.8374691         6.498622
#> 3 6.673834    0.90458097 0.8822257         5.050404
#> 4 5.636487   -1.16025892 0.6581716         4.304468
#> 5 7.993049   -0.02603191 1.2788021         6.322253
#> 6 5.315588   -0.71979537 0.3878004         4.085695
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
#> Random effect variance (sigma2_u): 1.937337 
#> Spatial autocorrelation (rho): 0.1222 
#> Spatial mixing fraction (phi): 0.1222 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  3.0275e+00  6.9084e-01  4.3824e+00  1.1739e-05  1.6838e+00
#> x1          -5.9700e-04  8.8271e-03 -6.7632e-02  9.4608e-01 -1.8201e-02
#> x2           4.7972e-02  5.5706e-02  8.6117e-01  3.8914e-01 -6.0873e-02
#> x3           4.1073e-02  5.2813e-02  7.7770e-01  4.3674e-01 -6.3814e-02
#>             ci_upper
#> (Intercept)   4.4054
#> x1            0.0166
#> x2            0.1585
#> x3            0.1441
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred        sd       mse      rse ci_lower
#> 1      1 8.359527 7.297146    7.297146 0.7516880 0.5650348 10.30112 5.843470
#> 2      2 7.599650 6.485992    6.485992 0.8084611 0.6536093 12.46473 4.929167
#> 3      3 5.514137 5.022828    5.022828 0.7981279 0.6370082 15.89001 3.473567
#> 4      4 3.869326 4.307604    4.307604 0.7084316 0.5018753 16.44607 2.907338
#> 5      5 6.305063 6.284932    6.284932 0.8840853 0.7816069 14.06674 4.549571
#> 6      6 3.926807 4.109193    4.109193 0.5750883 0.3307266 13.99516 2.978849
#>   ci_upper random_effect    vardir
#> 1 8.789626    2.61748921 0.6618838
#> 2 8.098339    2.27861174 0.8374691
#> 3 6.605361    0.91262521 0.8822257
#> 4 5.687144   -1.12684962 0.6581716
#> 5 8.022934   -0.05624276 1.2788021
#> 6 5.234751   -0.66297248 0.3878004
#> ... and 36 more rows.
#> 

# 4. The same Fay-Herriot model by MCMC (NumPyro NUTS on JAX)
if (check_numpyro_available()) {
  m_jax <- hb_area(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = "vardir",
    engine = jax_engine(samples = 500, chains = 2)
  )
}
# }
```
