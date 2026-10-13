# Empirical Best Linear Unbiased Prediction under a Two-fold Fay-Herriot Model.

This function gives the Empirical Best Linear Unbiased Prediction
(EBLUP) under a two-fold Fay-Herriot model (Torabi & Rao, 2014):
subareas nested within areas with an area random effect \\v_i\\ and a
subarea random effect \\u\_{ij}\\.

## Usage

``` r
eblup_twofold(
  formula,
  vardir,
  domain = NULL,
  subarea = NULL,
  data,
  method = c("REML", "ML"),
  mse = c("analytical", "bootstrap"),
  B = 200,
  seed = NULL,
  maxiter = 100,
  precision = 1e-04,
  self_benchmark = FALSE,
  benchmark_weight = NULL,
  benchmark_target = NULL,
  benchmark_group = NULL,
  print_result = TRUE
)
```

## Arguments

- formula:

  an object of class formula that contains a description of the model to
  be fitted.

- vardir:

  vector or column names from data that contain variance sampling from
  the direct estimator.

- domain:

  vector, column name or one-sided formula referencing the area
  identifier column in `data`. Must be supplied.

- subarea:

  vector, column name or one-sided formula referencing the subarea
  identifier column in `data`. If NULL, rows are numbered consecutively.

- data:

  a data frame or a data frame extension (e.g. a tibble), one row per
  subarea.

- method:

  Fitting method can be chosen between 'REML' and 'ML'.

- mse:

  MSE estimation method: 'analytical' (Prasad-Rao g1+g2+2\*g3, fast; ML
  adds the bias correction) or 'bootstrap' (parametric bootstrap, slower
  but robust).

- B:

  number of parametric bootstrap replicates (only used if mse =
  'bootstrap').

- seed:

  integer seed for bootstrap resampling (NULL = no seed set).

- maxiter:

  maximum number of iterations allowed in the Fisher-scoring algorithm.

- precision:

  convergence tolerance limit for the Fisher-scoring algorithm.

- self_benchmark:

  Logical. If `TRUE`, fits a self-benchmarking model via the augmented
  model approach (Wang, Fuller, and Qu, 2008; Rao and Molina, 2015, Sec
  10.5.1), ensuring that the weighted aggregate of small area estimates
  automatically matches the direct survey aggregate or specified
  benchmark targets without post-hoc adjustment.

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
  grouping / stratum column in `data` (e.g. area identifier for
  subarea-to-area benchmarking) for group-specific self-benchmarking.

- print_result:

  print coefficient or not, default value is TRUE.

## Value

The function returns a list with the following objects: `estcoef` a data
frame with the estimated model coefficients, `random_effect_var` named
vector with the estimated area (`sigma2_v`) and subarea (`sigma2_u`)
variances (paper notation), `goodness` vector containing several
goodness-of-fit measures, `df_eblup` a data frame that contains domain,
subarea, y, eblup, vardir, random_effect_area, random_effect_subarea,
mse, and rse.  

## Details

The model has a form that is response ~ auxiliary variables, with one
row per subarea. When the response contains NA the subarea is treated as
non-sampled and predicted with the synthetic estimator. MSE is estimated
either analytically (Prasad-Rao g1+g2+2\*g3 for REML, minus the ML bias
term b'grad(g1) for ML, with an exact matrix derivation for g3) or by
parametric bootstrap.

## References

1.  Torabi, M., & Rao, J. N. K. (2014). On small area estimation under a
    sub-area level model. *Journal of Multivariate Analysis*, 127,
    36–55.

2.  Rao, J. N. K., & Molina, I. (2015). *Small Area Estimation*. John
    Wiley & Sons.

## Examples

``` r
library(fastsae)

# Two-fold Fay-Herriot model (simulated illustration)
set.seed(1)
m <- 20
dat <- do.call(rbind, lapply(1:m, function(d) {
  nd <- sample(2:5, 1)
  data.frame(
    area = d, subarea = paste0(d, "-", seq_len(nd)),
    x1 = rnorm(nd), y = NA_real_, vardir = runif(nd, 0.2, 1)
  )
}))
dat$y <- 1 + dat$x1 + rnorm(m)[dat$area] + rnorm(nrow(dat), sd = 0.5) +
  rnorm(nrow(dat), sd = sqrt(dat$vardir))
m1 <- eblup_twofold(y ~ x1, vardir = "vardir", domain = "area",
                    subarea = "subarea", data = dat)
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> eblup_twofold(formula = y ~ x1, vardir = "vardir", domain = "area", subarea =
#> "subarea", data = dat)
#> 
#> ✔ Convergence: Yes (in 5 iterations)
#> Model: Two-fold Fay-Herriot (Area/Subarea-level)
#> Method: REML
#> Area effect variance (sigma2_v): 0.923152 
#> Subarea effect variance (sigma2_u): 0.557106 
#> 
#> Fixed Effects Coefficients:
#>                beta std.error  zvalue pvalue
#> (Intercept) 0.91925   0.25146 3.65565  3e-04
#> x1          1.09990   0.19812 5.55167  0e+00
#> 
#> EBLUP Estimates (First 6 domains):
#>   domain subarea          y      eblup    vardir random_effect_area
#> 1      1     1-1  0.6970537  0.8609432 0.9187117        0.399903297
#> 2      1     1-2  3.7070229  3.1225076 0.9557402        0.399903297
#> 3      2     2-1  2.3966129  1.7498606 0.7740948        0.002731102
#> 4      2     2-2 -1.6190158 -0.5691516 0.9935249        0.002731102
#> 5      2     2-3  1.2446340  1.3460273 0.5040281        0.002731102
#> 6      2     2-4  2.3206222  1.9710147 0.8219562        0.002731102
#>   random_effect_subarea       mse       rse
#> 1           -0.09938256 0.5451035  85.75609
#> 2            0.34071720 0.5680651  24.13768
#> 3            0.46545954 0.4319259  37.55789
#> 4           -0.58869784 0.5027167 124.57580
#> 5           -0.11207085 0.3373565  43.15098
#> 6            0.23695732 0.4494131  34.01206
#> ... and 67 more rows.
#> 
```
