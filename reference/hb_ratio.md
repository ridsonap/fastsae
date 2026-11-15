# Hierarchical Bayes Ratio Model for a Sub-Proportion of a Parent Indicator

Estimates an indicator that is a subset of a parent indicator (for
instance extreme poverty within poverty, severe stunting within
stunting) as the product of the parent estimate and a modelled
sub-proportion: \\\theta_d = s_d P_d\\, where \\P_d\\ is a given (e.g.
officially released) estimate of the parent indicator and \\s_d\\ is the
share of the parent population that has the sub-characteristic. The
sub-proportion follows a binomial area-level model
([`hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
with `family = "binomial"`) whose trials are the (effective) numbers of
sampled parent cases. Since \\0 \le s_d \le 1\\, the estimates never
exceed the parent (\\\theta_d \le P_d\\), and an optional benchmark on
the logit scale of \\s_d\\ keeps that property.

## Usage

``` r
hb_ratio(
  formula,
  data,
  trials,
  parent,
  parent_se = NULL,
  domain = NULL,
  engine = "inla",
  benchmark_target = NULL,
  benchmark_weight = NULL,
  benchmark_group = NULL,
  n_draws = 2000,
  seed = 1L,
  print_result = TRUE,
  ...
)
```

## Arguments

- formula:

  Formula for the direct sub-proportion, `s ~ x1 + x2`: the response is
  the direct estimate of \\s_d\\ (for instance the direct ratio of the
  sub-indicator to the parent), or the number of parent cases with the
  sub-characteristic.

- data:

  Data frame with one row per area.

- trials:

  Vector, column name, or formula: (effective) numbers of sampled parent
  cases, e.g. from
  [`effective_trials`](https://ridsonap.is-a.dev/fastsae/reference/effective_trials.md)
  with the parent cases as `n`.

- parent:

  Vector, column name, or formula: the parent estimate \\P_d\\
  (non-negative).

- parent_se:

  Optional vector, column name, or formula: standard errors of `parent`.
  With them, `sd`, `mse`, `rse` and the intervals include the
  uncertainty of the parent (treated as independent of \\s_d\\); without
  them they are conditional on `parent`.

- domain:

  Optional column name or formula of area identifiers.

- engine:

  Estimation engine of the sub-proportion model, as in
  [`hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md):
  [`inla_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)
  (default, or `"inla"`) or
  [`jax_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)
  (or `"jax"`, MCMC with NumPyro). See Details for the posterior draws
  of each engine.

- benchmark_target:

  Optional benchmark target for the weighted mean of \\\theta_d\\
  (scalar, or named by `benchmark_group`).

- benchmark_weight:

  Vector, column name, or formula of area weights for the benchmark
  (e.g. population sizes); needed with `benchmark_target`.

- benchmark_group:

  Optional vector, column name, or formula: groups with their own
  target.

- n_draws:

  Number of posterior draws from the INLA fit (default 2000). Not used
  with
  [`jax_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md),
  whose MCMC draws (`chains` x `samples`) are taken as they are.

- seed:

  Integer seed of the posterior draws from the INLA fit and of the
  parent draws (default 1); `NULL` uses the session's random number
  stream. The seed of the MCMC sampler is `jax_engine(seed = )`.

- print_result:

  Logical. If `TRUE` (default), prints a summary.

- ...:

  Further arguments to
  [`hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
  (e.g. `spatial`, `W`, `var_x`); `family` is `"binomial"`.

## Value

An object of class `c("fastsae_hb_ratio", "fastsae_hb", "fastsae")` with

- df_hb:

  One row per area: `domain`, `y` (direct sub-proportion times
  `parent`), `hb` (posterior mean of \\\theta_d\\), `sd`, `mse`, `rse`,
  `ci_lower`, `ci_upper`, `rse_cond` (RSE given the parent), `ratio`,
  `ratio_sd` (posterior mean and SD of the possibly benchmarked
  \\s_d\\), `ratio_direct`, `parent`, `parent_se` and `trials`.

- draws:

  Matrix of posterior draws of \\\theta_d\\ (areas x draws) behind
  `df_hb` (for instance for `add_reliability_flags(draws = )`).

- ratio_fit:

  The
  [`hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
  fit of the sub-proportion (use it for
  [`diagnose()`](https://ridsonap.is-a.dev/fastsae/reference/diagnose.md)
  of the model).

- benchmark:

  `NULL`, or a list with the targets, groups and the weighted means of
  the benchmarked estimates.

and `estcoef`, `random_effect_var`, `goodness`, `convergence`,
`inla_mode_status` (INLA) or the MCMC diagnostics (JAX, see Details),
`engine`, `formula`, `model` and `call`.

## Details

The benchmark solves, for each posterior draw, for one shift \\a_g\\ per
group on the logit scale such that \\\sum\_{d \in g} w_d P_d\\
\mathrm{logit}^{-1}(\mathrm{logit}(s_d) + a_g) / \sum\_{d \in g} w_d =
T_g\\, so the benchmarked \\s_d\\ stay in (0, 1) and the estimates below
their parent. The parent enters at its given value; the parent draws for
the total uncertainty are \\\max(P_d + \mathrm{se}\_d z, 0)\\, \\z \sim
N(0, 1)\\.

Modelling the sub-proportion borrows strength for the part of the
indicator that the parent does not already explain, and it keeps the two
indicators coherent. In an application to extreme poverty within poverty
for 514 Indonesian districts, the ratio model gave a relative standard
error of at most 50% in 124 districts against 77 for a binomial HB model
of the extreme poverty rate itself, and no district estimate above its
poverty rate (against 7).

Both engines of
[`hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md) fit
the sub-proportion model, with the same default priors. They differ in
where the posterior draws of \\s_d\\ come from. With
[`inla_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md),
`n_draws` draws are sampled from the fitted approximation
([`INLA::inla.posterior.sample()`](https://rdrr.io/pkg/INLA/man/posterior.sample.html),
which needs the package sn), so `ratio` agrees with `ratio_fit$df_hb$hb`
up to Monte Carlo error. With
[`jax_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md),
the draws are those of the sampler, so the two are equal without a
benchmark, and the object carries the convergence diagnostics of the
sub-proportion model (`rhat_max`, `ess_min`, `divergences`,
`treedepth_hits`, `ebfmi`, `mcmc_diagnostics`; see
[`jax_engine`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)).
`autoplot(fit, type = "trace")` and `type = "rank"` plot the chains of
that model, where `hb[j]` is \\s_d\\ of row `j`. `hb_ratio()` always
returns `draws`; `jax_engine(keep_draws = FALSE)` only drops the draws
of \\s_d\\ from `ratio_fit`. The logit link is the only one of the JAX
engine.

## See also

[`hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md),
[`inla_engine`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md),
[`jax_engine`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md),
[`effective_trials`](https://ridsonap.is-a.dev/fastsae/reference/effective_trials.md),
[`benchmark_sae`](https://ridsonap.is-a.dev/fastsae/reference/benchmark_sae.md).

## Examples

``` r
# \donttest{
if (requireNamespace("INLA", quietly = TRUE)) {
  set.seed(1)
  D <- 60
  d <- data.frame(area = 1:D, x = rnorm(D), pop = rpois(D, 5e4))
  d$p0 <- plogis(-2 + 0.3 * rnorm(D))                 # released parent rates
  d$n_poor <- rpois(D, 40)                             # sampled parent cases
  s_true <- plogis(-1.5 + 0.5 * d$x + rnorm(D, 0, 0.3))
  d$s_dir <- rbinom(D, d$n_poor, s_true) / d$n_poor
  fit <- hb_ratio(s_dir ~ x,
    data = d, trials = "n_poor", parent = "p0", domain = "area",
    benchmark_target = 0.025, benchmark_weight = "pop", print_result = FALSE
  )
  head(fit$df_hb)
  all(fit$df_hb$hb <= d$p0)
}
#> ℹ Response `y` appears to be proportions; using `y * trials` as the number of cases (not rounded).
#> [1] TRUE
# }
if (FALSE) { # \dontrun{
# MCMC (needs Python with jax and numpyro, see setup_numpyro_env())
fit_jax <- hb_ratio(s_dir ~ x,
  data = d, trials = "n_poor", parent = "p0", domain = "area",
  engine = jax_engine(chains = 4)
)
fit_jax$rhat_max
} # }
```
