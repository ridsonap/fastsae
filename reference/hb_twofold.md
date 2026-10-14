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
#> Area effect variance (sigma2_v): 0.678307 
#> Subarea effect variance (sigma2_u): 0.170329 
#> Spatial mixing fraction (phi): NA 
#> 
#> Fixed Effects Coefficients:
#>                   beta  std.error     zvalue     pvalue   ci_lower ci_upper
#> (Intercept) 9.1340e-01 2.9996e-01 3.0451e+00 2.3264e-03 3.1811e-01   1.5085
#> x1          9.4362e-01 1.5006e-01 6.2883e+00 3.2093e-10 6.4220e-01   1.2337
#> 
#> HB Estimates (First 6 domains):
#>   domain subarea         y        hb linear_pred    vardir        sd       mse
#> 1      1     1-1 2.5249687 2.5482114   2.5482114 0.5114576 0.5170036 0.2672927
#> 2      1     1-2 2.3333192 2.0762013   2.0762012 0.6419530 0.5311865 0.2821591
#> 3      2     2-1 1.2440865 1.2459548   1.2459548 0.4773757 0.4771370 0.2276597
#> 4      2     2-2 0.9335319 0.7945495   0.7945494 0.7640087 0.5212428 0.2716940
#> 5      2     2-3 2.0476243 2.2515073   2.2515074 0.7869359 0.5354552 0.2867123
#> 6      3     3-1 1.6429514 1.4191230   1.4191229 0.5085271 0.4719533 0.2227399
#>        rse   ci_lower ci_upper random_effect_area random_effect_subarea
#> 1 20.28888  1.5285090 3.562186         0.18728806           0.002828283
#> 2 25.58454  1.0409221 3.133966         0.18728806           0.072596142
#> 3 38.29489  0.3066695 2.186590        -0.04683652          -0.002253695
#> 4 65.60231 -0.2234912 1.838036        -0.04683652           0.028498144
#> 5 23.78208  1.1770459 3.293328        -0.04683652          -0.041871224
#> 6 33.25669  0.4982635 2.362046         0.28213124           0.079713080
#> ... and 22 more rows.
#> 
#> Area-Level Aggregates (First 6 areas):
#>   domain    hb_area   sd_area  mse_area rse_area ci_lower_area ci_upper_area
#> 1      1  2.3122064 0.4394125 0.1930834 19.00404   1.368339872    3.05173140
#> 2      2  1.4306705 0.3873340 0.1500276 27.07360   0.642547138    2.11116682
#> 3      3  1.4455023 0.3985903 0.1588743 27.57452   0.768847926    2.27042395
#> 4      4  2.8920432 0.4031087 0.1624967 13.93854   2.037701586    3.60785382
#> 5      5  0.6615616 0.3377652 0.1140853 51.05575  -0.005930511    1.30464801
#> 6      6 -0.7020719 0.3747730 0.1404548 53.38101  -1.410804595   -0.08154649
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
#> sigma2_v (area): 0.678307 
#> sigma2_u (subarea): 0.170329 
#> phi (spatial fraction): NA 
#> 
#> Coefficients:
#>                   beta  std.error     zvalue     pvalue   ci_lower ci_upper
#> (Intercept) 9.1340e-01 2.9996e-01 3.0451e+00 2.3264e-03 3.1811e-01   1.5085
#> x1          9.4362e-01 1.5006e-01 6.2883e+00 3.2093e-10 6.4220e-01   1.2337
#> 
#> Hyperparameters:
#>                                   mean        sd 0.025quant 0.5quant 0.975quant
#> Precision for ..area_id..     1.931688  1.077468  0.6063656 1.686755   4.701146
#> Precision for ..subarea_id.. 18.692396 36.965123  1.5468611 8.859484  99.375561
#>                                  mode
#> Precision for ..area_id..    1.286718
#> Precision for ..subarea_id.. 3.380214
#> 
#> Goodness of Fit:
#> $dic
#> [1] 76.67427
#> 
#> $p_eff_dic
#> [1] 14.29156
#> 
#> $waic
#> [1] 77.94939
#> 
#> $p_eff_waic
#> [1] 12.00841
#> 
#> $mlik
#> log marginal-likelihood (integration) 
#>                             -51.36874 
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
#>  Mean   : 3656.04  
#>  3rd Qu.:   61.03  
#>  Max.   :99748.25  
#> 
# }
```
