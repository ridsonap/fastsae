# Spatial and Spatio-Temporal Models

## Overview

When small areas are geographic units (such as counties, regencies, or
districts), spatial correlation frequently occurs: neighboring areas
tend to have more similar outcomes than distant ones. Furthermore, when
surveys are repeated over multiple years, temporal correlation between
time periods is also present.

**fastsae** provides specialized, high-performance implementations for
both: 1. **Spatial Fay-Herriot (`eblup_sfh`)**: Simultaneous
Autoregressive SAR(1) random effects. 2. **Spatio-Temporal Fay-Herriot
(`eblup_stfh`)**: Combined spatial SAR(1) and temporal AR(1) random
effects.

------------------------------------------------------------------------

## 1. Spatial Fay-Herriot Model (`eblup_sfh`)

### Model Formulation

The Spatial Fay-Herriot model (Pratesi & Salvati, 2008) models area
random effects using a Simultaneous Autoregressive (SAR(1)) process:

``` math
y = X\beta + u + e
```

``` math
u = \rho_1 W u + \epsilon_1 \implies u = (I - \rho_1 W)^{-1} \epsilon_1
```

where: - $`W`$ is a known row-standardized spatial proximity matrix
($`D \times D`$). - $`\rho_1 \in (-1, 1)`$ is the spatial autoregressive
parameter. - $`\epsilon_1 \sim N(0, \sigma_{u1}^2 I_D)`$ is independent
innovation noise. - $`e \sim N(0, V_e)`$ is sampling error with diagonal
matrix $`V_e = \text{diag}(D_1, \dots, D_D)`$.

### Fitting the Model

``` r

library(fastsae)
data("mys")
data("mys_proxmat")

# Fit Spatial Fay-Herriot model with REML
fit_sfh <- eblup_sfh(
  formula = y ~ x1 + x2 + x3,
  vardir = ~vardir,
  data = mys,
  W = mys_proxmat,
  method = "REML",
  print_result = FALSE
)

summary(fit_sfh)
#> 
#> ── Summary of fastsae Fit ──────────────────────────────────────────────────────
#> Call :
#> eblup_sfh(formula = y ~ x1 + x2 + x3, vardir = ~vardir, data = mys, method =
#> "REML", W = mys_proxmat, print_result = FALSE)
#> 
#> ✔ Convergence: Yes (in 9 iterations)
#> Model: Spatial Fay-Herriot (Area-level SAR)
#> Method: REML
#> 
#> Variance & Correlation Components:
#>  parameter  estimate std.error
#>   sigma2_u  1.396229 0.6748491
#>        rho -1.000000        NA
#> rho (spatial): -1 
#> 
#> Coefficients:
#>                    beta   std.error      zvalue pvalue
#> (Intercept)  3.02111730  0.61861801  4.88365555 0.0000
#> x1          -0.00072863  0.00805609 -0.09044507 0.9279
#> x2           0.04916379  0.05175135  0.95000019 0.3421
#> x3           0.04007654  0.04917325  0.81500691 0.4151
#> 
#> Goodness of Fit:
#> loglikelihood           AIC           BIC 
#>     -65.25028     142.50056     151.29498 
#> 
#> EBLUP Summary Statistics:
#>       mse              rse        
#>  Min.   :0.2006   Min.   : 9.955  
#>  1st Qu.:0.6195   1st Qu.:14.123  
#>  Median :0.8687   Median :20.411  
#>  Mean   :1.0842   Mean   :22.426  
#>  3rd Qu.:1.6293   3rd Qu.:29.634  
#>  Max.   :2.2268   Max.   :45.090
```

### Unsampled Domains and Spatial Kriging

A major advantage of `eblup_sfh` in **fastsae** is its automatic
handling of unsampled domains (domains with `y = NA`). Instead of
throwing an error or requiring manual subsetting, `eblup_sfh`
automatically performs **full-spatial kriging**:

``` r

# Count sampled vs unsampled domains
table(is.na(mys$y))
#> 
#> FALSE  TRUE 
#>    32    10

# Inspect estimates for unsampled domains
head(fit_sfh$df_eblup[is.na(mys$y), ])
#>    domain  y vardir    eblup random_effect      mse      rse
#> 21     21 NA      0 3.649858   -0.01707331 1.798100 36.73931
#> 25     25 NA      0 4.148486   -0.03953465 1.754488 31.92902
#> 26     26 NA      0 5.293005   -0.01455237 2.086796 27.29214
#> 27     27 NA      0 4.283739   -0.01062303 1.766047 31.02260
#> 28     28 NA      0 5.024852   -0.01952204 1.718356 26.08755
#> 34     34 NA      0 3.479814   -0.01480941 1.748554 37.99999
```

For unsampled areas, prediction borrows strength from both the
regression synthetic component $`x_d^\top \hat{\beta}`$ and spatial
proximity to neighboring sampled areas through
$`(I - \hat{\rho}_1 W)^{-1}`$.

### Multi-Threaded Bootstrap MSE

In addition to analytical Prasad-Rao style MSE, `eblup_sfh` supports
multi-threaded OpenMP bootstrap MSE estimation:

- **Parametric Bootstrap (`mse_method = "pbmse"`)**
- **Non-Parametric Bootstrap (`mse_method = "npbmse"`)**

``` r

fit_sfh_pb <- eblup_sfh(
  formula = y ~ x1 + x2 + x3,
  vardir = ~vardir,
  data = mys,
  W = mys_proxmat,
  mse_method = "pbmse",
  B = 50,
  n_threads = 2,
  seed = 123,
  print_result = FALSE
)
#> ℹ 10 unsampled domain(s) detected. Bootstrap MSE (pbmse) is computed for the 32 sampled domain(s); unsampled domain(s) get a full-spatial synthetic (kriging) prediction and analytical MSE instead.

head(fit_sfh_pb$df_eblup[, c("domain", "y", "eblup", "mse", "mse_pb", "mse_pbbc")])
#>   domain        y    eblup       mse    mse_pb  mse_pbbc
#> 1      1 8.359527 7.229608 0.5180135 0.5500837 0.5912158
#> 2      2 7.599650 6.369903 0.6092215 0.6103917 0.6867114
#> 3      3 5.514137 4.996497 0.6504225 0.8841162 0.7758816
#> 4      4 3.869326 4.343987 0.5213006 0.5271967 0.5896981
#> 5      5 6.305063 6.344684 0.7994159 0.7407658 0.9094689
#> 6      6 3.926807 4.094482 0.3422232 0.4123618 0.4430582
```

------------------------------------------------------------------------

## 2. Spatio-Temporal Fay-Herriot Model (`eblup_stfh`)

### Model Formulation

The Spatio-Temporal model (Marhuenda, Molina, & Morales, 2013) combines
spatial correlation and temporal autoregression across a balanced panel
of $`D`$ areas observed over $`T`$ time periods:

``` math
y_{dt} = x_{dt}^\top \beta + u_{1d} + u_{2dt} + e_{dt}
```

where: - $`u_1 = (u_{11}, \dots, u_{1D})^\top`$ captures static spatial
effects: $`u_1 = \rho_1 W u_1 + \epsilon_1`$. -
$`u_{2d} = (u_{2d1}, \dots, u_{2dT})^\top`$ captures dynamic temporal
effects for area $`d`$:
$`u_{2dt} = \rho_2 u_{2d,t-1} + \epsilon_{2dt}`$. -
$`e_{dt} \sim \text{ind. } N(0, D_{dt})`$ are sampling errors.

### Fitting the Spatio-Temporal Model

``` r

data("mys_panel")

# Prepare balanced panel without missing values
panel_data <- mys_panel[!is.na(mys_panel$y) & mys_panel$year >= 2024, ]
W_sub <- mys_proxmat[-c(21, 25), -c(21, 25)]

# Fit Spatio-Temporal model with the analytical MSE
fit_stfh <- eblup_stfh(
  formula = y ~ x1 + x2 + x3,
  data = panel_data,
  vardir = ~vardir,
  domain = ~area,
  time = ~year,
  W = W_sub,
  model = "ST",
  compute_mse = TRUE,
  print_result = FALSE
)

summary(fit_stfh)
#> 
#> ── Summary of fastsae Fit ──────────────────────────────────────────────────────
#> Call :
#> eblup_stfh(formula = y ~ x1 + x2 + x3, vardir = ~vardir, data = panel_data,
#> domain = ~area, time = ~year, W = W_sub, model = "ST", compute_mse = TRUE,
#> print_result = FALSE)
#> 
#> ✔ Convergence: Yes (in 7 iterations)
#> Model: Spatio-Temporal Fay-Herriot (SAR + AR(1) time effects)
#> Method: REML
#> 
#> Variance & Correlation Components:
#>    estimate std.error
#>   2.6699676 0.7419111
#>  -0.9990000        NA
#>   0.0000000 0.1473516
#>   0.3477329        NA
#> 
#> Coefficients:
#>                   beta  std.error     zvalue pvalue
#> (Intercept)  3.7761285  0.5612917  6.7275685 0.0000
#> x1          -0.0021031  0.0069943 -0.3006916 0.7636
#> x2           0.0394660  0.0362503  1.0887079 0.2763
#> x3           0.0344259  0.0349429  0.9852041 0.3245
#> 
#> Goodness of Fit:
#>   loglike       AIC       BIC 
#> -198.2219  412.4438  434.7437 
#> 
#> EBLUP Summary Statistics:
#>       mse              rse         
#>  Min.   :0.1286   Min.   :  5.791  
#>  1st Qu.:0.2735   1st Qu.:  9.516  
#>  Median :0.3681   Median : 10.835  
#>  Mean   :0.4665   Mean   : 15.801  
#>  3rd Qu.:0.4572   3rd Qu.: 14.707  
#>  Max.   :1.7254   Max.   :119.247
```

The estimated variance and correlation components include: - `sigma21`:
Spatial random effect variance ($`\sigma_{u1}^2`$). - `rho1`: Spatial
autocorrelation ($`\rho_1`$). - `sigma22`: Temporal random effect
variance ($`\sigma_{u2}^2`$). - `rho2`: Temporal AR(1) autoregression
coefficient ($`\rho_2`$).

The analytical MSE is $`g_1 + g_2 + 2 g_3 - g_4`$ (Prasad-Rao type, with
the exact derivative of the EBLUP in $`g_3`$ and the curvature term
$`g_4`$ for the nonlinear $`\rho_1`$, $`\rho_2`$); the components are
returned as `fit_stfh$g1`, …, `fit_stfh$g4`. A parametric bootstrap MSE
is available with `mse_type = "bootstrap"` (and `B`, `seed`).

------------------------------------------------------------------------

## 3. Comparing Models with `autoplot`

You can directly compare EBLUP estimates and uncertainty across
different models using
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html):

``` r

# Compare Fay-Herriot vs Spatial Fay-Herriot
fit_fh <- eblup_fh(y ~ x1 + x2 + x3, vardir = ~vardir, data = na.omit(mys), print_result = FALSE)
W_clean <- mys_proxmat[!is.na(mys$y), !is.na(mys$y)]
fit_sfh_clean <- eblup_sfh(y ~ x1 + x2 + x3, vardir = ~vardir, data = na.omit(mys), W = W_clean, print_result = FALSE)

autoplot(list("Standard FH" = fit_fh, "Spatial FH" = fit_sfh_clean), type = "comparison")
```

![](spatial-temporal_files/figure-html/compare_models-1.png)

------------------------------------------------------------------------

## 4. Bayesian Spatio-Temporal Modeling (`hb_area`)

In addition to frequentist C++ EBLUP implementations (`eblup_sfh` and
`eblup_stfh`), **fastsae** provides
[`hb_area()`](https://ridsonap.github.io/fastsae/reference/hb_area.md)
for Bayesian spatio-temporal modeling powered by **INLA**. This
supports: - Bounded proportions via Beta regression (matching
`tipsae`) - Discrete counts via Poisson and Negative Binomial
regression - Binary aggregates via Binomial regression - Flexible
space-time structures including domain-specific random walks, separable
dynamic fields, and Knorr-Held Types I–IV.

For full details and a case study on Italian poverty data, see the
vignette [Bayesian Area-Level Small Area Estimation with
INLA](https://ridsonap.github.io/fastsae/articles/hb-area-inla.md).

------------------------------------------------------------------------

## References

- Marhuenda, Y., Molina, I., & Morales, D. (2013). Small area estimation
  with spatio-temporal Fay-Herriot models. *Computational Statistics &
  Data Analysis*, 58, 308–325.
- Pratesi, M., & Salvati, N. (2008). Small area estimation for spatially
  correlated data: A Fay-Herriot with the spatial linear spline model.
  *Journal of Applied Statistics*, 35(7), 781–794.
