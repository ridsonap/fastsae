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

  A `data.frame` containing the area-level data (one row per area).

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
    tipsae spatio-temporal model).

  - `"separable"`: Grouped dynamic spatial field evolving over time via
    AR(1) / RW (corresponds to classical Spatio-Temporal Fay-Herriot
    models, Marhuenda et al. 2013).

  - `"type1"`: Knorr-Held (2000) Type I interaction (unstructured space
    \\\times\\ unstructured time).

  - `"type2"`: Knorr-Held Type II interaction (unstructured space
    \\\times\\ structured time).

  - `"type3"`: Knorr-Held Type III interaction (structured space
    \\\times\\ unstructured time).

  - `"type4"`: Knorr-Held Type IV interaction (structured space
    \\\times\\ structured time).

- W:

  Proximity or spatial adjacency matrix. Can be a square `matrix`,
  `Matrix`, or `spdep` `nb` or `listw` object. Dimensions must match the
  total number of domains in `data`. Required when `spatial != "none"`.

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

An object of class `c("fastsae_hb_area", "fastsae")` containing:

- `df_hb`: Data frame with domain estimates, including `domain`, `time`
  (if specified), observed `y`, predicted `hb`, linear predictor
  `linear_pred`, posterior standard error `sd`, `mse`, relative error
  `rse` (%), 95% credible interval (`ci_lower`, `ci_upper`), and
  `random_effect`.

- `estcoef`: Data frame of estimated regression coefficients.

- `hyperpar`: Data frame of hyperparameter posterior estimates.

- `random_effect_var`: Estimated area random effect variance.

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

|                                   |                |                           |
|-----------------------------------|----------------|---------------------------|
| **Feature**                       | **inla**       | **jax**                   |
| Families                          | all six        | all six                   |
| `spatial`                         | all            | all but `"slm"`           |
| `temporal`                        | all            | all but `"rw2"`           |
| `st_interaction`                  | all            | all                       |
| Self-benchmarking                 | posterior, wfq | posterior; wfq (gaussian) |
| Horseshoe prior, Student-t errors | –              | yes                       |

The `"jax"` engine uses weakly informative priors (half-normal standard
deviations, uniform spatial mixing / autocorrelation parameters) instead
of the PC priors of INLA, so its estimates differ slightly. For
`family = "gamma"` it estimates a common shape and does not use
`vardir`.

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

  # 3. Spatial BYM2 Gaussian Fay-Herriot with INLA, faster approximation
  m_bym2 <- hb_area(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = "vardir",
    family = "gaussian",
    spatial = "bym2",
    W = mys_proxmat,
    engine = inla_engine(strategy = "simplified.laplace")
  )
}
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> hb_area(formula = y ~ x1 + x2 + x3, data = mys, family = "gaussian", vardir =
#> "vardir")
#> 
#> ✔ Convergence: Yes (in - iterations)
#> Model: HB-GAUSSIAN (Non-spatial)
#> Engine: inla
#> Random effect variance (sigma2_u): 1.886439 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  3.0429e+00  6.9339e-01  4.3885e+00  1.1414e-05  1.6955e+00
#> x1          -1.1779e-03  8.8063e-03 -1.3375e-01  8.9360e-01 -1.8784e-02
#> x2           5.1663e-02  5.5581e-02  9.2952e-01  3.5262e-01 -5.6772e-02
#> x3           3.7759e-02  5.2727e-02  7.1613e-01  4.7391e-01 -6.7095e-02
#>             ci_upper
#> (Intercept)   4.4275
#> x1            0.0160
#> x2            0.1621
#> x3            0.1405
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred        sd       mse      rse ci_lower
#> 1      1 8.359527 7.346485    7.346485 0.7431548 0.5522791 10.11579 5.909880
#> 2      2 7.599650 6.505052    6.505052 0.8038271 0.6461380 12.35697 4.956248
#> 3      3 5.514137 5.053768    5.053768 0.7981023 0.6369672 15.79222 3.502543
#> 4      4 3.869326 4.302273    4.302273 0.7107212 0.5051247 16.51967 2.898577
#> 5      5 6.305063 6.322332    6.322332 0.8821568 0.7782007 13.95303 4.589863
#> 6      6 3.926807 4.084922    4.084922 0.5731467 0.3284971 14.03079 2.958629
#>   ci_upper random_effect    vardir
#> 1 8.822834    2.67527841 0.6618838
#> 2 8.107624    2.29158429 0.8374691
#> 3 6.634360    0.90447903 0.8822257
#> 4 5.687244   -1.16086830 0.6581716
#> 5 8.055299   -0.02614769 1.2788021
#> 6 5.206913   -0.72017671 0.3878004
#> ... and 36 more rows.
#> 
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> hb_area(formula = y ~ x1 + x2 + x3, data = mys, vardir = "vardir",
#> self_benchmark = TRUE, benchmark_weight = "n")
#> 
#> ✔ Convergence: Yes (in - iterations)
#> Model: HB-GAUSSIAN (Non-spatial)
#> Engine: inla
#> Self-Benchmarking: Active (benchmarked posterior draws)
#> Random effect variance (sigma2_u): 1.886439 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  3.0429e+00  6.9338e-01  4.3885e+00  1.1412e-05  1.6955e+00
#> x1          -1.1778e-03  8.8062e-03 -1.3375e-01  8.9360e-01 -1.8784e-02
#> x2           5.1663e-02  5.5580e-02  9.2952e-01  3.5262e-01 -5.6772e-02
#> x3           3.7759e-02  5.2726e-02  7.1614e-01  4.7391e-01 -6.7094e-02
#>             ci_upper
#> (Intercept)   4.4275
#> x1            0.0160
#> x2            0.1621
#> x3            0.1405
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred        sd       mse      rse ci_lower
#> 1      1 8.359527 7.271818    7.346465 0.7848122 0.6159302 10.79252 5.731621
#> 2      2 7.599650 6.475403    6.505031 0.8239815 0.6789455 12.72479 4.947109
#> 3      3 5.514137 5.024318    5.053758 0.8336537 0.6949784 16.59238 3.457790
#> 4      4 3.869326 4.257796    4.302280 0.7358881 0.5415314 17.28331 2.820859
#> 5      5 6.305063 6.287153    6.322332 0.8485922 0.7201087 13.49724 4.596055
#> 6      6 3.926807 4.006735    4.084924 0.6253046 0.3910058 15.60634 2.757263
#>   ci_upper random_effect    vardir hb_unbenchmarked
#> 1 8.846776    2.67525510 0.6618838         7.346465
#> 2 8.083912    2.29156222 0.8374691         6.505031
#> 3 6.692077    0.90447345 0.8822257         5.053758
#> 4 5.691186   -1.16085488 0.6581716         4.302280
#> 5 7.978053   -0.02614671 1.2788021         6.322332
#> 6 5.246663   -0.72016817 0.3878004         4.084924
#> ... and 36 more rows.
#> 
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> hb_area(formula = y ~ x1 + x2 + x3, data = mys, family = "gaussian", spatial =
#> "bym2", W = mys_proxmat, vardir = "vardir", engine = inla_engine(strategy =
#> "simplified.laplace"))
#> 
#> ✔ Convergence: Yes (in - iterations)
#> Model: HB-GAUSSIAN (BYM2)
#> Engine: inla
#> Random effect variance (sigma2_u): 1.876016 
#> Spatial autocorrelation (rho): 0.4816 
#> Spatial mixing fraction (phi): 0.4816 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  3.0436e+00  6.7542e-01  4.5062e+00  6.5998e-06  1.7297e+00
#> x1          -1.1805e-03  8.7842e-03 -1.3439e-01  8.9309e-01 -1.8691e-02
#> x2           5.1683e-02  5.5471e-02  9.3171e-01  3.5149e-01 -5.6610e-02
#> x3           3.7718e-02  5.2621e-02  7.1679e-01  4.7351e-01 -6.6710e-02
#>             ci_upper
#> (Intercept)   4.3884
#> x1            0.0159
#> x2            0.1617
#> x3            0.1403
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred        sd       mse      rse ci_lower
#> 1      1 8.359527 7.348365    7.348365 0.7396355 0.5470607 10.06531 5.914503
#> 2      2 7.599650 6.506743    6.506743 0.8002489 0.6403983 12.29876 4.960544
#> 3      3 5.514137 5.054871    5.054871 0.7974464 0.6359208 15.77580 3.503512
#> 4      4 3.869326 4.301792    4.301792 0.7103434 0.5045878 16.51273 2.899965
#> 5      5 6.305063 6.322403    6.322403 0.8823151 0.7784800 13.95538 4.589856
#> 6      6 3.926807 4.084758    4.084758 0.5731818 0.3285374 14.03221 2.958635
#>   ci_upper random_effect    vardir
#> 1 8.813957    2.67914346 0.6618838
#> 2 8.097758    2.29504075 0.8374691
#> 3 6.632256    0.90609242 0.8822257
#> 4 5.686897   -1.16239487 0.6581716
#> 5 8.055402   -0.02604528 1.2788021
#> 6 5.207004   -0.72114836 0.3878004
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
