# Estimation engines for the Hierarchical Bayes models

Constructors for the `engine` argument of
[`hb_area`](https://ridsonap.github.io/fastsae/reference/hb_area.md),
[`hb_unit`](https://ridsonap.github.io/fastsae/reference/hb_unit.md) and
[`hb_twofold`](https://ridsonap.github.io/fastsae/reference/hb_twofold.md).
Each collects the settings of one estimation engine, so that the model
functions themselves only describe the model. Passing the string
`"inla"` or `"jax"` is the same as passing `inla_engine()` or
`jax_engine()` with their defaults.

## Usage

``` r
inla_engine(
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
  seed = 0L,
  ...
)

jax_engine(
  warmup = 500,
  samples = 1000,
  chains = 2,
  device = c("auto", "cpu", "cuda", "metal"),
  seed = 42,
  prior_beta = c("normal", "horseshoe"),
  robust = FALSE,
  keep_draws = TRUE,
  progress = interactive()
)

# S3 method for class 'fastsae_engine'
print(x, ...)
```

## Arguments

- strategy:

  INLA approximation strategy: `"laplace"` (default, full Laplace
  approximation, highest accuracy), `"simplified.laplace"` (faster) or
  `"gaussian"` (fastest).

- link:

  Link function for the likelihood
  ([`hb_area()`](https://ridsonap.github.io/fastsae/reference/hb_area.md)
  only). `NULL` (default) uses the canonical link of the family:
  identity for gaussian, logit for binomial and beta, log for poisson,
  nbinomial and gamma.

- scale_model:

  Logical; scale the intrinsic GMRF models (besag, bym, rw1, rw2) to a
  generalized variance of 1. Default `TRUE`.

- prior_prec:

  Prior for the precision of the area (domain) random effect. `NULL`
  (default): the PC prior `list(prior = "pc.prec", param = c(1, 0.01))`
  in
  [`hb_area()`](https://ridsonap.github.io/fastsae/reference/hb_area.md)
  and
  [`hb_twofold()`](https://ridsonap.github.io/fastsae/reference/hb_twofold.md),
  `list(prior = "loggamma", param = c(0.01, 0.01))` in
  [`hb_unit()`](https://ridsonap.github.io/fastsae/reference/hb_unit.md).

- prior_phi:

  Prior for the mixing parameter of `spatial = "bym2"`. `NULL`
  (default): the PC prior `list(prior = "pc", param = c(0.5, 0.5))`.

- prior_rho:

  Prior for the spatial correlation of `spatial = "generic1"` or `"slm"`
  ([`hb_area()`](https://ridsonap.github.io/fastsae/reference/hb_area.md)
  only); `NULL` uses the INLA default.

- prior_prec_time:

  Prior for the precision of the temporal random effect
  ([`hb_area()`](https://ridsonap.github.io/fastsae/reference/hb_area.md)
  only); `NULL` uses the INLA default.

- prior_rho_time:

  Prior for the autocorrelation of `temporal = "ar1"`
  ([`hb_area()`](https://ridsonap.github.io/fastsae/reference/hb_area.md)
  only); `NULL` uses the INLA default.

- prior_prec_subarea:

  Prior for the precision of the subarea random effect
  ([`hb_twofold()`](https://ridsonap.github.io/fastsae/reference/hb_twofold.md)
  only). `NULL` (default):
  `list(prior = "pc.prec", param = c(1, 0.01))`.

- draws:

  Number of posterior draws from
  [`INLA::inla.posterior.sample()`](https://rdrr.io/pkg/INLA/man/posterior.sample.html):
  for posterior self-benchmarking in
  [`hb_area()`](https://ridsonap.github.io/fastsae/reference/hb_area.md)
  (`NULL`: 1000) and for the area aggregates of
  [`hb_twofold()`](https://ridsonap.github.io/fastsae/reference/hb_twofold.md)
  (`NULL`: 200; `0` skips the sampling and approximates the aggregate
  uncertainty from independent sub-area effects).

- seed:

  Random seed of the sampler. Default 42.

- ...:

  Further arguments passed to
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html), such as
  `control.compute` or `num.threads`.

- warmup:

  Number of warmup (adaptation) iterations per chain. Default 500.

- samples:

  Number of posterior draws kept per chain. Default 1000.

- chains:

  Number of MCMC chains. Default 2. With more than one chain, the fit
  reports the largest split R-hat and `convergence` is
  `rhat_max < 1.05`.

- device:

  `"auto"` (default: CUDA when available, otherwise CPU), `"cpu"`,
  `"cuda"` or `"metal"` (Apple silicon GPU, through jax-metal). JAX
  fixes its device at the first fit of an R session. `"metal"` has a
  fixed cost of about two minutes per fit, so on Apple silicon it pays
  off only for large numbers of domains (thousands); for that reason
  `"auto"` uses the CPU on macOS.

- prior_beta:

  Prior for the regression coefficients: `"normal"` (default) or
  `"horseshoe"` (regularized horseshoe, for many candidate covariates).

- robust:

  Logical; Student-t sampling errors instead of normal (Gaussian family
  only). Default `FALSE`.

- keep_draws:

  Logical; keep the posterior draws of the area parameters in `$draws`
  of the fit (needed by
  [`benchmark_sae()`](https://ridsonap.github.io/fastsae/reference/benchmark_sae.md)
  to propagate the posterior covariance). Default `TRUE`.

- progress:

  Logical; show NumPyro's progress bar. Default
  [`interactive()`](https://rdrr.io/r/base/interactive.html).

- x:

  A `"fastsae_engine"` object.

## Value

An object of class `"fastsae_engine"`.

## INLA

`inla_engine()` fits the model by integrated nested Laplace
approximation (INLA); it is available for all HB models. Priors are
given in INLA syntax. A prior or setting left at `NULL` takes the
default of the model function; a setting the model does not use (such as
`prior_rho_time` in
[`hb_unit()`](https://ridsonap.github.io/fastsae/reference/hb_unit.md))
is ignored with a warning.

## JAX

`jax_engine()` fits the model by MCMC (NumPyro's No-U-Turn sampler on
JAX, Phan et al. 2019) on CPU, CUDA GPU or Apple Metal, through
reticulate. It needs Python with `jax` and `numpyro`; see
[`setup_numpyro_env`](https://ridsonap.github.io/fastsae/reference/setup_numpyro_env.md)
and
[`check_numpyro_available`](https://ridsonap.github.io/fastsae/reference/check_numpyro_available.md).
The fixed effects get a normal or a regularized horseshoe prior
(Piironen and Vehtari 2017), and the variance components half-normal
priors. `spatial = "slm"` and `temporal = "rw2"` are not available.
`jax_engine()` is available for
[`hb_area()`](https://ridsonap.github.io/fastsae/reference/hb_area.md)
only.

## References

Phan, D., Pradhan, N. and Jankowiak, M. (2019). Composable effects for
flexible and accelerated probabilistic programming in NumPyro.
*arXiv:1912.11554*.

Piironen, J. and Vehtari, A. (2017). Sparsity information and
regularization in the horseshoe and other shrinkage priors. *Electronic
Journal of Statistics*, 11(2), 5018–5051.

## See also

[`hb_area`](https://ridsonap.github.io/fastsae/reference/hb_area.md),
[`hb_unit`](https://ridsonap.github.io/fastsae/reference/hb_unit.md),
[`hb_twofold`](https://ridsonap.github.io/fastsae/reference/hb_twofold.md)

## Examples

``` r
inla_engine(strategy = "simplified.laplace")
#> INLA engine
#> • Strategy: simplified.laplace
#> • Link: canonical
jax_engine(samples = 2000, chains = 4, device = "cuda")
#> JAX engine (NumPyro NUTS)
#> • 4 chains x 2000 draws, 500 warmup, seed 42
#> • Device: cuda
#> • Prior on coefficients: normal
```
