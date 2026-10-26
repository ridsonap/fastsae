# Two-Fold Hierarchical Bayes for Sub-Area Level Small Area Estimation

Estimates small area parameters under the two-fold sub-area level model
of Torabi and Rao (2014), formulated as a hierarchical Bayes model,
using Integrated Nested Laplace Approximations (INLA) or, with
`engine = jax_engine()`, MCMC (NumPyro's No-U-Turn sampler on JAX). The
model accounts for nested random effects at both the primary area level
(\\v_i\\) and the nested sub-area level (\\u\_{ij}\\), supporting
Gaussian, Binomial (logistic), and Poisson responses, with optional
spatial correlation structures (BYM2, Besag) across areas and aggregate
area-level estimation.

## Usage

``` r
hb_twofold(
  formula,
  vardir = NULL,
  domain,
  subarea = NULL,
  data,
  family = c("gaussian", "binomial", "poisson"),
  spatial = c("none", "bym2", "besag"),
  W = NULL,
  weight = NULL,
  trials = NULL,
  exposure = NULL,
  engine = "inla",
  compute_area = TRUE,
  print_result = TRUE,
  ...
)
```

## Arguments

- formula:

  An object of class `formula` specifying the fixed-effects model (e.g.,
  `y ~ x1 + x2`).

- vardir:

  Character string or vector specifying direct sampling variances
  (\\\psi\_{ij}\\) for continuous (Gaussian) responses. Required when
  `family = "gaussian"`.

- domain:

  Character string, column name, or vector referencing the primary area
  identifier in `data`. Must be supplied.

- subarea:

  Optional character string, column name, or vector referencing the
  sub-area identifier in `data`. If `NULL`, sub-areas are numbered
  consecutively per domain.

- data:

  A `data.frame` or tibble containing one row per sub-area.

- family:

  Character string specifying the response likelihood. Options:

  - `"gaussian"`: Continuous response with known sampling variance
    (default).

  - `"binomial"`: Binary / count proportion with sample sizes specified
    by `trials`.

  - `"poisson"`: Count response with expected exposures specified by
    `exposure`.

- spatial:

  Character string specifying the spatial structure across primary
  areas:

  - `"none"`: Non-spatial independent and identically distributed (IID)
    area effects (default).

  - `"bym2"`: Scaled Besag-York-Mollié 2 spatial model (Riebler et al.,
    2016).

  - `"besag"`: Intrinsic Conditional Autoregressive (ICAR) spatial
    model.

- W:

  Proximity or spatial adjacency matrix across unique primary areas. Can
  be a square `matrix`, `Matrix`, or `spdep` `nb` or `listw` object.
  Dimensions must match the total number of unique primary areas in
  `domain`. Required when `spatial != "none"`.

- weight:

  Optional character string or numeric vector specifying sub-area
  population weights (\\w\_{ij}\\) used to aggregate sub-area
  predictions to primary area estimates \\\hat{\theta}\_i = \sum_j
  w\_{ij} \hat{\theta}\_{ij}\\. If `NULL`, equal weights per area are
  used.

- trials:

  Optional character string or vector specifying total number of trials
  per sub-area. Required when `family = "binomial"`.

- exposure:

  Optional character string or vector specifying expected baseline
  exposures per sub-area. Required when `family = "poisson"`.

- engine:

  Estimation engine: `"inla"` (default) or an
  [`inla_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)
  object, or `"jax"` or a
  [`jax_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)
  object.
  [`inla_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)
  holds the approximation strategy, the priors (`prior_prec` for the
  area effect \\\tau_v\\ and `prior_prec_subarea` for the sub-area
  effect \\\tau_u\\, both by default the prior of
  [`inla_engine`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)`(prior_prec = )`;
  `prior_phi` for BYM2), `scale_model`, the number of posterior `draws`
  for the area aggregates of the binomial and Poisson models
  (default 2000) and their `seed`, and further
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html) arguments.
  The Gaussian model is fitted in INLA's `"classic"` mode, the binomial
  and Poisson models in INLA's default `"compact"` mode, as in
  [`hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md);
  `inla_engine(inla.mode = )` selects a mode.
  [`jax_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)
  samples the same model by MCMC (see the Engines section).

- compute_area:

  Logical. If `TRUE` (default), computes aggregate primary area-level
  estimates and posterior standard errors: for the Gaussian model with
  INLA from the exact posterior of the weighted sums of the sub-area
  predictors (INLA linear combinations), otherwise from joint posterior
  draws (sampled from the INLA fit, or the MCMC draws with
  [`jax_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)).

- print_result:

  Logical. If `TRUE` (default), prints a summary of results.

- ...:

  Deprecated. `strategy`, `scale_model`, `prior_prec_area` (now
  `prior_prec`), `prior_prec_subarea`, `prior_phi`, `n_samples` (now
  `draws`) and arguments of
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html) are now given
  to
  [`inla_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md);
  passed here they are forwarded to it with a warning.

## Value

An object of class `c("fastsae_hb_twofold", "fastsae_hb", "fastsae")`
containing:

- `df_hb`: Data frame with sub-area level estimates:

  - `domain`: Primary area identifier.

  - `subarea`: Sub-area identifier.

  - `y`: Direct estimate / observed response (NA for unsampled).

  - `hb`: Posterior mean estimate of sub-area mean.

  - `linear_pred`: Posterior mean of linear predictor.

  - `vardir`: Direct sampling variance (for Gaussian).

  - `sd`: Posterior standard deviation (standard error).

  - `mse`: Posterior variance (\\\text{sd}^2\\), named as in the EBLUP
    results.

  - `rse`: Relative Standard Error (%).

  - `ci_lower`: 2.5% quantile of posterior credible interval.

  - `ci_upper`: 97.5% quantile of posterior credible interval.

  - `random_effect_area`: Posterior mean of area random effect
    (\\\hat{v}\_i\\).

  - `random_effect_subarea`: Posterior mean of sub-area random effect
    (\\\hat{u}\_{ij}\\).

- `df_area`: Data frame with aggregate primary area-level estimates:

  - `domain`: Primary area identifier.

  - `hb_area`: Aggregate area mean estimate.

  - `sd_area`: Posterior standard error of area estimate.

  - `mse_area`: Posterior variance of the area estimate.

  - `rse_area`: Relative Standard Error of area estimate (%).

  - `ci_lower_area`: 2.5% quantile of area credible interval.

  - `ci_upper_area`: 97.5% quantile of area credible interval.

  - `n_subareas`: Number of sub-areas in the area.

- `df_subarea`: Alias pointing to `df_hb`.

- `estcoef`: Data frame of estimated regression coefficients (\\\beta\\,
  standard error, z-value, p-value, credible intervals).

- `hyperpar`: Data frame of hyperparameter posterior estimates.

- `random_effect_var`: Named vector with the posterior means of the area
  (`sigma2_v`) and sub-area (`sigma2_u`) variances. When a variance is
  close to 0 its posterior is skewed and INLA's summary can be biased
  (on simulated data `sigma2_u` was 0.047 against an exact posterior
  mean of 0.072); the sub-area and area estimates are the reliable
  output.

- `phi`: Estimated spatial variance proportion (for BYM2).

- `goodness`: Model fit metrics (DIC, WAIC, Marginal Log-Likelihood).

- `family`: Response family used.

- `spatial`: Spatial model type used.

- `model`: Label of model (`"HB-TWOFOLD"`).

- `method`: Description of method.

- `convergence`: Logical indicating convergence.

- `level`: `"subarea"`.

- `fit`: Raw fitted INLA model object, or the list returned by the
  Python backend for the JAX engine.

- `engine`, `engine_spec`: Engine name (`"inla"` or `"jax"`) and its
  settings.

- `draws`, `rhat_max`, `device`, `elapsed_seconds`: With the JAX engine,
  the posterior draws of the sub-area parameters (rows of `df_hb` by
  draws, kept when `jax_engine(keep_draws = TRUE)`), the largest split
  R-hat (`convergence` is `rhat_max < 1.05`), the JAX device and the
  sampling time.

- `call`: Matched function call.

## Engines

Both engines fit the same model with the same default priors and return
the same structure. The standard deviations of the area effect (iid,
Besag or BYM2) and of the iid sub-area effect have half-normal priors
with scale \\s_y\\, the standard deviation of the sub-area direct
estimates, for `family = "gaussian"`, and scale 1 on the link scale
otherwise; the BYM2 mixing parameter has INLA's PC prior and the
coefficients \\N(0, 1000)\\ priors.
`jax_engine(priors = "weakly_informative")` uses half-normal(1) standard
deviations on the standardised scale instead. See
[`vignette("hb-priors", package = "fastsae")`](https://ridsonap.is-a.dev/fastsae/articles/hb-priors.md).
Under the earlier PC priors of
[`inla_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md),
on simulated data (20 areas, 56 sub-areas; Gaussian, binomial, Poisson
and Gaussian with BYM2) the area means of the two engines correlated
above 0.999 and differed by a median of 0.02-0.07 posterior SD, with
area SDs 0.85-1.14 times those of INLA. The posterior means of the
variances were larger with MCMC (area variance 0.93 against 0.68 for the
Gaussian model): their posteriors are skewed, and INLA's summary is
biased downwards for small variances. `jax_engine(robust = TRUE)` gives
Student-t area and sub-area effects and `prior_beta = "horseshoe"` a
regularized horseshoe prior on the coefficients.

## References

1.  Torabi, M., & Rao, J. N. K. (2014). On small area estimation under a
    sub-area level model. *Journal of Multivariate Analysis*, 127,
    36-55.

2.  Rao, J. N. K., & Molina, I. (2015). *Small Area Estimation* (2nd
    ed.). John Wiley & Sons.

3.  Riebler, A., Sørbye, S. H., Simpson, D., & Rue, H. (2016). An
    intuitive Bayesian spatial model for disease mapping that accounts
    for scaling. *Statistical Methods in Medical Research*, 25(4),
    1145-1165.

## Examples

``` r
# \donttest{
if (requireNamespace("INLA", quietly = TRUE)) {
  library(fastsae)
  set.seed(42)
  m <- 10
  dat <- do.call(rbind, lapply(1:m, function(d) {
    nd <- sample(2:4, 1)
    data.frame(
      area = d,
      subarea = paste0(d, "-", seq_len(nd)),
      x1 = rnorm(nd),
      vardir = runif(nd, 0.2, 0.8)
    )
  }))
  dat$y <- 1 + dat$x1 + rnorm(m)[dat$area] + rnorm(nrow(dat), sd = 0.5) +
    rnorm(nrow(dat), sd = sqrt(dat$vardir))

  fit <- hb_twofold(y ~ x1,
    vardir = "vardir", domain = "area",
    subarea = "subarea", data = dat
  )
  summary(fit)
}
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> hb_twofold(formula = y ~ x1, vardir = "vardir", domain = "area", subarea =
#> "subarea", data = dat)
#> 
#> ✔ Convergence: Yes (INLA mode status 0)
#> Model: HB-TWOFOLD-GAUSSIAN (Non-spatial)
#> Engine: inla
#> Method: INLA (laplace)
#> Area effect variance (sigma2_v): 1.263078 
#> Subarea effect variance (sigma2_u): 0.352764 
#> Spatial mixing fraction (phi): NA 
#> 
#> Fixed Effects Coefficients:
#>                   beta  std.error     zvalue     pvalue   ci_lower ci_upper
#> (Intercept) 9.1768e-01 3.8578e-01 2.3788e+00 1.7369e-02 1.4671e-01   1.6901
#> x1          9.4873e-01 1.6723e-01 5.6731e+00 1.4025e-08 6.0891e-01   1.2704
#> 
#> HB Estimates (First 6 domains):
#>   domain subarea         y        hb linear_pred    vardir        sd       mse
#> 1      1     1-1 2.5249687 2.5697364   2.5697364 0.5114576 0.5562661 0.3094320
#> 2      1     1-2 2.3333192 2.1251792   2.1251792 0.6419530 0.5815634 0.3382160
#> 3      2     2-1 1.2440865 1.2421247   1.2421247 0.4773757 0.5157415 0.2659893
#> 4      2     2-2 0.9335319 0.8007818   0.8007817 0.7640087 0.5769019 0.3328158
#> 5      2     2-3 2.0476243 2.2346991   2.2346992 0.7869359 0.5897514 0.3478067
#> 6      3     3-1 1.6429514 1.4666033   1.4666033 0.5085271 0.5139411 0.2641354
#>        rse   ci_lower ci_upper random_effect_area random_effect_subarea
#> 1 21.64681  1.4711497 3.659016         0.21027680          -0.010542152
#> 2 27.36538  0.9925328 3.284346         0.21027680           0.089974094
#> 3 41.52091  0.2270468 2.259164        -0.05803548          -0.001137452
#> 4 72.04234 -0.3246619 1.958577        -0.05803548           0.042063658
#> 5 26.39064  1.0474818 3.380134        -0.05803548          -0.059243463
#> 6 35.04295  0.4640277 2.494332         0.30833357           0.096278956
#> ... and 22 more rows.
#> 
#> Area-Level Aggregates (First 6 areas):
#>   domain    hb_area   sd_area  mse_area rse_area ci_lower_area ci_upper_area
#> 1      1  2.3474578 0.4832042 0.2334863 20.58415    1.40292015     3.2995334
#> 2      2  1.4258685 0.4277236 0.1829475 29.99741    0.58607025     2.2648080
#> 3      3  1.4775234 0.3967959 0.1574470 26.85547    0.70201868     2.2589307
#> 4      4  3.0259776 0.4045413 0.1636536 13.36894    2.24472718     3.8313999
#> 5      5  0.6768729 0.3570974 0.1275186 52.75694   -0.02280287     1.3784494
#> 6      6 -0.8734795 0.3886167 0.1510230 44.49065   -1.63855488    -0.1153071
#>   n_subareas
#> 1          2
#> 2          3
#> 3          3
#> 4          2
#> 5          3
#> 6          3
#> ... and 4 more areas.
#> 
#> 
#> ── Summary of fastsae Fit ──────────────────────────────────────────────────────
#> Call :
#> hb_twofold(formula = y ~ x1, vardir = "vardir", domain = "area", subarea =
#> "subarea", data = dat)
#> 
#> ✔ Convergence: Yes (in - iterations)
#> Model: HB-TWOFOLD-GAUSSIAN (Non-spatial)
#> Engine: inla
#> Method: INLA (laplace)
#> 
#> Variance Components:
#> sigma2_v (area): 1.263078 
#> sigma2_u (subarea): 0.352764 
#> phi (spatial fraction): NA 
#> 
#> Coefficients:
#>                   beta  std.error     zvalue     pvalue   ci_lower ci_upper
#> (Intercept) 9.1768e-01 3.8578e-01 2.3788e+00 1.7369e-02 1.4671e-01   1.6901
#> x1          9.4873e-01 1.6723e-01 5.6731e+00 1.4025e-08 6.0891e-01   1.2704
#> 
#> Hyperparameters:
#>                                  mean        sd 0.025quant 0.5quant 0.975quant
#> Precision for ..area_id..    1.175629 0.8048668  0.2753772 0.972060   3.279658
#> Precision for ..subarea_id.. 5.328773 5.7246477  0.9067990 3.643309  20.280658
#>                                   mode
#> Precision for ..area_id..    0.6551048
#> Precision for ..subarea_id.. 1.9881957
#> 
#> Goodness of Fit:
#> $dic
#> [1] 76.50284
#> 
#> $p_eff_dic
#> [1] 16.65468
#> 
#> $waic
#> [1] 75.47097
#> 
#> $p_eff_waic
#> [1] 11.83889
#> 
#> $mlik
#> log marginal-likelihood (integration) 
#>                             -49.23015 
#> 
#> 
#> HB Summary Statistics:
#>        hb           linear_pred            sd              mse        
#>  Min.   :-1.7358   Min.   :-1.7358   Min.   :0.4117   Min.   :0.1695  
#>  1st Qu.: 0.2104   1st Qu.: 0.2104   1st Qu.:0.5083   1st Qu.:0.2584  
#>  Median : 1.2595   Median : 1.2595   Median :0.5463   Median :0.2984  
#>  Mean   : 1.0570   Mean   : 1.0570   Mean   :0.5318   Mean   :0.2864  
#>  3rd Qu.: 2.0838   3rd Qu.: 2.0838   3rd Qu.:0.5774   3rd Qu.:0.3334  
#>  Max.   : 3.3671   Max.   : 3.3671   Max.   :0.6495   Max.   :0.4218  
#>       rse         
#>  Min.   :  15.33  
#>  1st Qu.:  27.12  
#>  Median :  35.43  
#>  Mean   : 173.44  
#>  3rd Qu.:  71.63  
#>  Max.   :2226.32  
#> 
# }
```
