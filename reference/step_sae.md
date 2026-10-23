# Stepwise Selection of the Covariates of a Small Area Model

Selects the covariates of a fitted fastsae model by stepwise search
(backward, forward or both) on an information criterion, refitting the
model with the same settings (family, spatial or temporal structure,
engine, ...) at every step.

## Usage

``` r
step_sae(
  object,
  scope = NULL,
  direction = c("both", "backward", "forward"),
  criterion = NULL,
  steps = 1000,
  trace = TRUE
)
```

## Arguments

- object:

  A fitted model:
  [`eblup_fh()`](https://ridsonap.is-a.dev/fastsae/reference/eblup_fh.md),
  [`eblup_sfh()`](https://ridsonap.is-a.dev/fastsae/reference/eblup_sfh.md),
  [`eblup_stfh()`](https://ridsonap.is-a.dev/fastsae/reference/eblup_stfh.md),
  [`eblup_bhf()`](https://ridsonap.is-a.dev/fastsae/reference/eblup_bhf.md),
  [`eblup_twofold()`](https://ridsonap.is-a.dev/fastsae/reference/eblup_twofold.md),
  [`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md),
  [`hb_unit()`](https://ridsonap.is-a.dev/fastsae/reference/hb_unit.md)
  or
  [`hb_twofold()`](https://ridsonap.is-a.dev/fastsae/reference/hb_twofold.md).

- scope:

  The range of models searched: a formula with the largest set of
  covariates (the upper model), or a list with formulas `lower` and
  `upper`. Default: the covariates of `object` as upper model and the
  intercept-only model as lower model. Forward selection needs an upper
  model larger than `object`.

- direction:

  `"both"` (default), `"backward"` or `"forward"`.

- criterion:

  Information criterion to minimise. EBLUP models: `"AIC"`, `"BIC"`
  (default) or `"KIC"`. HB models: `"LCPO"` (default with
  [`inla_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)),
  `"WAIC"` (default with
  [`jax_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md))
  or `"DIC"`.

- steps:

  Maximum number of steps. Default 1000.

- trace:

  Logical; print the criterion of each candidate model. Default `TRUE`.

## Value

An object of class `"fastsae_step"`: a list with `model` (the selected
model, fitted with the settings of `object`), `path` (a data frame with
the term added or dropped at each step and the criterion), `criterion`,
`direction`, and the `initial` and `final` formulas.

## Details

**EBLUP models.** The criteria are computed from the full Gaussian
log-likelihood \\\ell\\ at the estimates, with \\k\\ the number of
coefficients plus variance parameters: \\\mathrm{AIC} = -2\ell + 2k\\,
\\\mathrm{BIC} = -2\ell + k \log m\\ and \\\mathrm{KIC} = -2\ell + 3k\\
(Marhuenda, Morales and Pardo 2014; the analytic criteria of
`emdi::step()`). The restricted (REML) likelihood depends on the
covariates through \\X\\ and cannot compare models with different
covariates, so models with a `method` argument are fitted by ML during
the search; the selected model is then refitted with the method of
`object`.
[`eblup_stfh()`](https://ridsonap.is-a.dev/fastsae/reference/eblup_stfh.md)
has no ML fit and is compared through the likelihood at its REML
estimates. Only sampled domains enter the criteria.

**HB models** (every family). `"LCPO"`, the negative mean log
conditional predictive ordinate \\-\frac{1}{m}\sum_d \log
\mathrm{CPO}\_d\\ (leave-one-out cross-validation; Held, Schrödle and
Rue 2010), and `"WAIC"` (Watanabe 2010) estimate the out-of-sample
predictive accuracy; `"DIC"` (Spiegelhalter et al. 2002) is less
reliable for hierarchical models. With one observation per area and an
area random effect, WAIC is unstable (Vehtari, Gelman and Gabry 2017):
in a simulation of a binomial area-level model (40 areas, two true and
two noise covariates) WAIC selected exactly the true covariates in 27%
of the replicates and dropped a true one in 13%, LCPO in 70%
(`benchmarks/validation/step_sae.R`). LCPO is therefore the default with
[`inla_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md);
[`jax_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)
has no CPO and uses WAIC. With
[`jax_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)
every candidate model is an MCMC fit, so the search can take long; the
criteria then also carry Monte Carlo error.

**Spatial and temporal models.** Every candidate is fitted with the
spatial or temporal structure of `object`. Selecting the covariates
without the spatial effect and adding it afterwards is not equivalent: a
spatially smooth covariate and the spatial random effect explain the
same variation (spatial confounding), so a covariate kept by a
non-spatial model may be redundant, and its coefficient can change
markedly, once the spatial effect is in the model (Reich, Hodges and
Zadnik 2006; Hodges and Reich 2010). In the simulation above, with a
spatial area effect and a spatially smooth noise covariate, the
Fay-Herriot model kept that covariate in 20% of the replicates and the
spatial Fay-Herriot model in 8%. Check the residual spatial
autocorrelation of the selected model, e.g. with
[`spatial_test()`](https://ridsonap.is-a.dev/fastsae/reference/spatial_test.md).

The model is refitted by evaluating its call with a new formula in the
environment where `step_sae()` is called, so the data and other
arguments of the call must be available there.

## References

Held, L., Schrödle, B. and Rue, H. (2010). Posterior and
cross-validatory predictive checks: a comparison of MCMC and INLA. In T.
Kneib and G. Tutz (eds.), *Statistical Modelling and Regression
Structures*, 91-110. Physica-Verlag.

Hodges, J. S. and Reich, B. J. (2010). Adding spatially-correlated
errors can mess up the fixed effect you love. *The American
Statistician*, 64(4), 325-334.

Marhuenda, Y., Morales, D. and Pardo, M. C. (2014). Information criteria
for Fay-Herriot model selection. *Computational Statistics & Data
Analysis*, 70, 268-280.

Reich, B. J., Hodges, J. S. and Zadnik, V. (2006). Effects of residual
smoothing on the posterior of the fixed effects in disease-mapping
models. *Biometrics*, 62(4), 1197-1206.

Spiegelhalter, D. J., Best, N. G., Carlin, B. P. and van der Linde, A.
(2002). Bayesian measures of model complexity and fit. *Journal of the
Royal Statistical Society: Series B*, 64(4), 583-639.

Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*, 27(5), 1413-1432.

Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
and widely applicable information criterion in singular learning theory.
*Journal of Machine Learning Research*, 11, 3571-3594.

## See also

[`compare_sae`](https://ridsonap.is-a.dev/fastsae/reference/compare_sae.md),
[`spatial_test`](https://ridsonap.is-a.dev/fastsae/reference/spatial_test.md)

## Examples

``` r
data(mys)
fit <- eblup_fh(y ~ x1 + x2 + x3, vardir = "vardir", data = mys, print_result = FALSE)
sel <- step_sae(fit, criterion = "BIC")
#> Start: BIC = 147.3513 (y ~ x1 + x2 + x3)
#>    Step      BIC
#>    - x1 143.9140
#>    - x3 144.2980
#>    - x2 144.7444
#>  <none> 147.3513
#> Step: - x1: BIC = 143.914
#>    Step      BIC
#>    - x3 141.0402
#>    - x2 141.6679
#>  <none> 143.9140
#>    + x1 147.3513
#> Step: - x3: BIC = 141.0402
#>    Step      BIC
#>  <none> 141.0402
#>    + x3 143.9140
#>    + x1 144.2980
#>    - x2 149.7081
sel
#> 
#> ── Stepwise covariate selection (BIC, both) ────────────────────────────────────
#> Initial model: y ~ x1 + x2 + x3
#> Final model: y ~ x2
#> 
#>  Step      BIC
#>       147.3513
#>  - x1 143.9140
#>  - x3 141.0402
sel$model$formula
#> y ~ x2
#> <environment: 0x55daff05a4b8>
```
