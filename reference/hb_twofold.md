# Two-Fold Hierarchical Bayes for Sub-Area Level Small Area Estimation

Estimates small area parameters under a two-fold sub-area level model
(Torabi & Rao, 2014 Bayesian hierarchical framework) using Integrated
Nested Laplace Approximations (INLA). The model accounts for nested
random effects at both the primary area level (\\v_i\\) and the nested
sub-area level (\\u\_{ij}\\), supporting Gaussian, Binomial (logistic),
and Poisson responses, with optional spatial correlation structures
(BYM2, Besag) across areas and aggregate area-level estimation.

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
  [`inla_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)
  object holding the approximation strategy, the priors (`prior_prec`
  for the area effect \\\tau_v\\ and `prior_prec_subarea` for the
  sub-area effect \\\tau_u\\, both by default
  `list(prior = "pc.prec", param = c(1, 0.01))`; `prior_phi` for BYM2),
  `scale_model`, the number of posterior `draws` for the area aggregates
  (default 200) and their `seed`, and further
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html) arguments.

- compute_area:

  Logical. If `TRUE` (default), computes aggregate primary area-level
  estimates and posterior standard errors via posterior sampling.

- print_result:

  Logical. If `TRUE` (default), prints a summary of results.

- ...:

  Deprecated. `strategy`, `scale_model`, `prior_prec_area` (now
  `prior_prec`), `prior_prec_subarea`, `prior_phi`, `n_samples` (now
  `draws`) and arguments of
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html) are now given
  to
  [`inla_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md);
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

  - `mse`: Posterior Mean Squared Error (\\\text{sd}^2\\).

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

  - `mse_area`: Posterior Mean Squared Error of area estimate.

  - `rse_area`: Relative Standard Error of area estimate (%).

  - `ci_lower_area`: 2.5% quantile of area credible interval.

  - `ci_upper_area`: 97.5% quantile of area credible interval.

  - `n_subareas`: Number of sub-areas in the area.

- `df_subarea`: Alias pointing to `df_hb`.

- `estcoef`: Data frame of estimated regression coefficients (\\\beta\\,
  standard error, z-value, p-value, credible intervals).

- `hyperpar`: Data frame of hyperparameter posterior estimates.

- `random_effect_var`: Named vector with area (`sigma2_v`) and sub-area
  (`sigma2_u`) variances.

- `phi`: Estimated spatial variance proportion (for BYM2).

- `goodness`: Model fit metrics (DIC, WAIC, Marginal Log-Likelihood).

- `family`: Response family used.

- `spatial`: Spatial model type used.

- `model`: Label of model (`"HB-TWOFOLD"`).

- `method`: Description of method.

- `convergence`: Logical indicating convergence.

- `level`: `"subarea"`.

- `fit`: Raw fitted INLA model object.

- `engine`, `engine_spec`: Engine name (`"inla"`) and its settings.

- `call`: Matched function call.

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

  fit <- hb_twofold(y ~ x1, vardir = "vardir", domain = "area",
                    subarea = "subarea", data = dat)
  summary(fit)
}
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> hb_twofold(formula = y ~ x1, vardir = "vardir", domain = "area", subarea =
#> "subarea", data = dat)
#> 
#> ✔ Convergence: Yes (in - iterations)
#> Model: HB-TWOFOLD-GAUSSIAN (Non-spatial)
#> Engine: inla
#> Method: INLA (laplace)
#> Area effect variance (sigma2_v): 0.678317 
#> Subarea effect variance (sigma2_u): 0.170334 
#> Spatial mixing fraction (phi): NA 
#> 
#> Fixed Effects Coefficients:
#>                   beta  std.error     zvalue     pvalue   ci_lower ci_upper
#> (Intercept) 9.1340e-01 2.9996e-01 3.0450e+00 2.3265e-03 3.1811e-01   1.5085
#> x1          9.4362e-01 1.5006e-01 6.2883e+00 3.2088e-10 6.4220e-01   1.2337
#> 
#> HB Estimates (First 6 domains):
#>   domain subarea         y       hb linear_pred    vardir        sd       mse
#> 1      1     1-1 2.5249687 2.548213   2.5482132 0.5114576 0.5170028 0.2672919
#> 2      1     1-2 2.3333192 2.076201   2.0762009 0.6419530 0.5311853 0.2821578
#> 3      2     2-1 1.2440865 1.245955   1.2459545 0.4773757 0.4771358 0.2276586
#> 4      2     2-2 0.9335319 0.794548   0.7945479 0.7640087 0.5212408 0.2716920
#> 5      2     2-3 2.0476243 2.251509   2.2515091 0.7869359 0.5354532 0.2867102
#> 6      3     3-1 1.6429514 1.419123   1.4191227 0.5085271 0.4719517 0.2227384
#>        rse   ci_lower ci_upper random_effect_area random_effect_subarea
#> 1 20.28884  1.5285125 3.562187         0.18728915           0.002827562
#> 2 25.58448  1.0409240 3.133963         0.18728915           0.072594580
#> 3 38.29480  0.3066717 2.186587        -0.04683704          -0.002253527
#> 4 65.60218 -0.2234888 1.838030        -0.04683704           0.028497291
#> 5 23.78197  1.1770522 3.293326        -0.04683704          -0.041870040
#> 6 33.25658  0.4982664 2.362042         0.28213336           0.079711439
#> ... and 22 more rows.
#> 
#> Area-Level Aggregates (First 6 areas):
#>   domain    hb_area   sd_area  mse_area rse_area ci_lower_area ci_upper_area
#> 1      1  2.3122071 0.4599466 0.2115508 19.89210    1.43187251    3.12452178
#> 2      2  1.4306705 0.4058428 0.1647084 28.36731    0.62529126    2.25033485
#> 3      3  1.4455033 0.3944195 0.1555667 27.28596    0.72273624    2.26323607
#> 4      4  2.8920432 0.3846355 0.1479445 13.29979    2.10664358    3.61491131
#> 5      5  0.6615624 0.3410308 0.1163020 51.54929    0.03919234    1.34196770
#> 6      6 -0.7020753 0.3742269 0.1400458 53.30296   -1.43796274   -0.03893978
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
#> sigma2_v (area): 0.678317 
#> sigma2_u (subarea): 0.170334 
#> phi (spatial fraction): NA 
#> 
#> Coefficients:
#>                   beta  std.error     zvalue     pvalue   ci_lower ci_upper
#> (Intercept) 9.1340e-01 2.9996e-01 3.0450e+00 2.3265e-03 3.1811e-01   1.5085
#> x1          9.4362e-01 1.5006e-01 6.2883e+00 3.2088e-10 6.4220e-01   1.2337
#> 
#> Hyperparameters:
#>                                   mean        sd 0.025quant 0.5quant 0.975quant
#> Precision for ..area_id..     1.931676  1.077486  0.6063528 1.686734   4.701193
#> Precision for ..subarea_id.. 18.689438 36.954030  1.5468367 8.858982  99.353688
#>                                  mode
#> Precision for ..area_id..    1.286689
#> Precision for ..subarea_id.. 3.380244
#> 
#> Goodness of Fit:
#> $dic
#> [1] 76.67422
#> 
#> $p_eff_dic
#> [1] 14.29145
#> 
#> $waic
#> [1] 77.9494
#> 
#> $p_eff_waic
#> [1] 12.00838
#> 
#> $mlik
#> log marginal-likelihood (integration) 
#>                             -51.36875 
#> 
#> 
#> HB Summary Statistics:
#>        hb           linear_pred            sd              mse        
#>  Min.   :-1.4828   Min.   :-1.4828   Min.   :0.3912   Min.   :0.1530  
#>  1st Qu.: 0.1835   1st Qu.: 0.1835   1st Qu.:0.4690   1st Qu.:0.2200  
#>  Median : 1.1956   Median : 1.1956   Median :0.4950   Median :0.2450  
#>  Mean   : 1.0545   Mean   : 1.0545   Mean   :0.4897   Mean   :0.2421  
#>  3rd Qu.: 2.0911   3rd Qu.: 2.0911   3rd Qu.:0.5301   3rd Qu.:0.2810  
#>  Max.   : 3.3102   Max.   : 3.3102   Max.   :0.5838   Max.   :0.3408  
#>       rse          
#>  Min.   :   15.00  
#>  1st Qu.:   26.05  
#>  Median :   33.64  
#>  Mean   : 3663.99  
#>  3rd Qu.:   61.03  
#>  Max.   :99970.90  
#> 
# }
```
