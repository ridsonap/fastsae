# Generalized Variance Function (GVF) Smoothing for Sampling Variances

Fits Generalized Variance Function (GVF) models to smooth unstable
direct sampling variance estimates (\\D_i\\) in Small Area Estimation.
Supports log-linear regression (Wolter, 2007), Gamma generalized linear
models (GLM with log link), and CV power models.

## Usage

``` r
gvf_smooth(
  vardir,
  y = NULL,
  n = NULL,
  method = c("log_linear", "gamma_glm", "cv_power"),
  data = NULL,
  formula = NULL,
  bias_correction = TRUE
)
```

## Arguments

- vardir:

  vector, column name, or one-sided formula referencing the direct
  sampling variances in `data`.

- y:

  optional vector, column name, or one-sided formula referencing the
  direct survey estimates.

- n:

  optional vector, column name, or one-sided formula referencing domain
  sample sizes.

- method:

  character string specifying the GVF smoothing model:

  `"log_linear"`

  :   Log-linear model: \\\log(D_i) = \alpha_0 + \alpha_1 \log(y_i) +
      \alpha_2 \log(n_i) + \varepsilon_i\\. Includes second-order bias
      correction upon exponentiation.

  `"gamma_glm"`

  :   Gamma GLM with log link: \\\mathrm{E}\[D_i\] = \mu_i, \log(\mu_i)
      = \alpha_0 + \dots\\. Naturally suited for strictly positive
      variances without back-transformation bias.

  `"cv_power"`

  :   Power model on direct coefficients of variation: \\\log(CV_i^2) =
      \alpha_0 + \dots\\.

  Default is `"log_linear"`.

- data:

  optional data frame containing the variables.

- formula:

  optional user-supplied formula overriding the default model
  specification.

- bias_correction:

  logical. For `"log_linear"`, whether to apply log-normal bias
  correction \\\tilde{D}\_i = \exp(\hat{\mu}\_i + \hat{\sigma}^2 / 2)\\.
  Default is `TRUE`.

## Value

An object of class `"gvf_smooth"` and `"fastsae_gvf"`, containing:

- `smooth_vardir`:

  numeric vector of smoothed sampling variances \\\tilde{D}\_i\\.

- `vardir`:

  original raw sampling variances.

- `fitted_model`:

  the fitted `lm` or `glm` model object.

- `method`:

  method used.

- `r_squared`:

  coefficient of determination \\R^2\\ (or pseudo-\\R^2\\ for GLM).

- `df`:

  data frame with raw `vardir`, `smooth_vardir`, and covariates.

- `call`:

  the matched function call.

## References

1.  Wolter, K. M. (2007). *Introduction to Variance Estimation* (2nd
    ed.). Springer.

2.  Rao, J. N. K., & Molina, I. (2015). *Small Area Estimation* (2nd
    ed.). John Wiley & Sons.

3.  Rivest, L. P., & Belmonte, E. (2000). A conditional mean squared
    error of small area estimators. *Survey Methodology*, 26(1), 67-78.

## Examples

``` r
library(fastsae)

# Smooth direct variances using direct estimates y and sample size n
gvf_res <- gvf_smooth(
  vardir = "vardir",
  y = "y",
  n = "n",
  data = mys,
  method = "log_linear"
)

# Inspect summary
summary(gvf_res)
#> -- Generalized Variance Function (GVF) Smoothing --
#> Method: "log_linear"
#> Model formula: log_vardir ~ log_y + log_n
#> R-squared / Pseudo R-squared: 0.2321
#> Number of domains: 42
#> Raw vardir range: [0, 10.12571]
#> Smoothed vardir range: [0, 3.24426]
#> 
#> -- Model Coefficients --
#>               Estimate Std. Error   t value   Pr(>|t|)
#> (Intercept)  2.2968261  1.1981633  1.916956 0.06514378
#> log_y        0.4782909  0.2327915  2.054589 0.04902744
#> log_n       -0.3356704  0.1418147 -2.366966 0.02482499

# Use smoothed variances directly in Fay-Herriot model
fit_fh <- eblup_fh(
  y ~ x1 + x2 + x3,
  data = mys,
  vardir = gvf_res$smooth_vardir
)
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> eblup_fh(formula = y ~ x1 + x2 + x3, vardir = gvf_res$smooth_vardir, data =
#> mys)
#> 
#> ✔ Convergence: Yes (in 5 iterations)
#> Model: Fay-Herriot (Area-level)
#> Method: eblup
#> Random effect variance (sigma2_u): 2.042318 
#> 
#> Fixed Effects Coefficients:
#>                   beta  std.error     zvalue pvalue
#> (Intercept)  2.9422920  0.7376432  3.9887742 0.0001
#> x1          -0.0047696  0.0100273 -0.4756575 0.6343
#> x2           0.0782809  0.0577486  1.3555478 0.1752
#> x3           0.0172653  0.0527029  0.3275958 0.7432
#> 
#> EBLUP Estimates (First 6 domains):
#>   domain        y    eblup   vardir random_effect       mse      rse
#> 1      1 8.359527 6.507267 1.858058    2.03594397 1.1309699 16.34282
#> 2      2 7.599650 5.676565 2.463466    1.59431901 1.2938894 20.03840
#> 3      3 5.514137 4.817897 2.478087    0.57380608 1.3523650 24.13734
#> 4      4 3.869326 4.655423 1.717030   -0.93502085 1.0972074 22.50014
#> 5      5 6.305063 6.331141 1.322897   -0.04025957 0.9197609 15.14801
#> 6      6 3.926807 4.421276 1.996022   -0.50593814 1.2099341 24.87902
#> ... and 36 more rows.
#> 
```
