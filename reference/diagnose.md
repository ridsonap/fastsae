# Diagnostic Evaluation of Small Area Estimation Models

Conducts a comprehensive diagnostic evaluation of small area estimation
models (HB or EBLUP). Evaluates estimation precision, efficiency gains
over direct estimators, external calibration and bias diagnostics (Brown
et al., 2001), goodness-of-fit tests, and residual spatial
autocorrelation (Moran's I).

## Usage

``` r
diagnose(
  object,
  W = NULL,
  truth = NULL,
  rse_threshold = 25,
  alpha_level = 0.05
)
```

## Arguments

- object:

  A fitted model object of class `"fastsae"` (e.g., from
  [`hb_area`](https://ridsonap.github.io/fastsae/reference/hb_area.md),
  [`eblup_fh`](https://ridsonap.github.io/fastsae/reference/eblup_fh.md),
  [`eblup_sfh`](https://ridsonap.github.io/fastsae/reference/eblup_sfh.md),
  or
  [`eblup_stfh`](https://ridsonap.github.io/fastsae/reference/eblup_stfh.md)).

- W:

  Optional spatial proximity or adjacency matrix of dimension \\D \times
  D\\. If `NULL` and the fitted model contains a spatial matrix (e.g.
  `eblup_sfh` or spatial `hb_area`), it is automatically extracted from
  `object`.

- truth:

  Optional numeric vector containing true parameter values \\\theta_d\\
  (useful for simulation studies).

- rse_threshold:

  Numeric. Threshold for reliable Relative Standard Error (default is
  25%).

- alpha_level:

  Numeric. Significance level for hypothesis tests (default is 0.05).

## Value

An object of class `c("fastsae_diagnose", "list")` containing:

- model_info:

  Basic model characteristics (model name, sample size, unsampled
  count).

- precision:

  Summary of direct vs SAE RSE, proportion of areas with RSE \<
  threshold, and MSE reduction ratios.

- brown_test:

  Brown et al. (2001) diagnostics: `bias_regression` (descriptive OLS
  fit of the direct estimates on the model estimates, no test; see
  Details) and `goodness_of_fit` (the \\W\\ statistic with its
  upper-tail \\\chi^2\\ p-value).

- spatial_test:

  Moran's I test for spatial autocorrelation in model residuals (if `W`
  is available).

- bayesian_metrics:

  For Bayesian fits: DIC, WAIC, marginal log-likelihood and, when they
  are available per domain, CPO summaries and a Kolmogorov-Smirnov test
  of the PIT values.

- simulation_metrics:

  Mean relative bias, mean absolute relative error, relative RMSE across
  domains, and coverage rate (if `truth` is supplied).

- df_diag:

  Data frame combining direct estimates, SAE predictions, MSE, RSE, and
  residuals for diagnostics.

- status:

  Overall assessment summary string.

## Details

**Goodness of fit.** \\W = \sum_d (y_d - \hat\theta_d)^2 / (\psi_d +
\mathrm{MSE}\_d)\\ over the sampled domains is compared with
\\\chi^2_m\\ (Brown et al., 2001); only large values indicate that the
model estimates differ from the direct estimates more than their
uncertainty allows, so the p-value is upper-tailed (as in `emdi`). The
statistic ignores the positive correlation between \\y_d\\ and
\\\hat\theta_d\\, so it is conservative: under a correct model \\W\\ is
typically well below its degrees of freedom.

**Bias regression.** The OLS fit \\y_d = \alpha + \beta\hat\theta_d\\ is
reported without a formal test of \\\alpha = 0, \beta = 1\\: for
shrinkage estimators such as the EBLUP the slope exceeds 1 even when the
model is correct (about \\1/\gamma_d\\ for the part of the variation not
explained by the covariates), so that test would reject correct models.
Inspect the regression (and `autoplot(x, type = "calibration")`) for
systematic deviations instead.

**Moran's I** on the residuals \\y_d - \hat\theta_d\\ of the sampled
domains uses the randomisation variance of Cliff and Ord (1981),
two-sided.

## References

1.  Brown, G., Chambers, R., Heady, P., & Heasman, D. (2001). Evaluation
    of small area estimation methods: An application to the British
    Labour Force Survey. *ONS Internal Report*.

2.  Rao, J. N. K., & Molina, I. (2015). *Small Area Estimation* (2nd
    ed.). John Wiley & Sons.

3.  Cliff, A. D., & Ord, J. K. (1981). *Spatial Processes: Models &
    Applications*. Pion London.

## Examples

``` r
library(fastsae)
data(mys)
data(mys_proxmat)

# 1. Fit Fay-Herriot model
fit_fh <- eblup_fh(y ~ x1 + x2, data = mys, vardir = ~vardir, domain = ~area)
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> eblup_fh(formula = y ~ x1 + x2, vardir = ~vardir, domain = ~area, data = mys)
#> 
#> ✔ Convergence: Yes (in 7 iterations)
#> Model: Fay-Herriot (Area-level)
#> Method: REML
#> Random effect variance (sigma2_u): 2.569182 
#> 
#> Fixed Effects Coefficients:
#>                   beta  std.error     zvalue pvalue
#> (Intercept)  3.0689476  0.7631917  4.0212018 0.0001
#> x1          -0.0041073  0.0090768 -0.4525074 0.6509
#> x2           0.0861173  0.0309576  2.7817847 0.0054
#> 
#> EBLUP Estimates (First 6 domains):
#>   domain        y    eblup    vardir random_effect       mse       rse
#> 1      1 8.359527 7.594807 0.6618838    2.96835204 0.5602338  9.855256
#> 2      2 7.599650 6.777056 0.8374691    2.52354833 0.6766501 12.137828
#> 3      3 5.514137 5.247512 0.8822257    0.77645304 0.7048881 15.999508
#> 4      4 3.869326 4.247072 0.6581716   -1.47453645 0.5556741 17.551750
#> 5      5 6.305063 6.353632 1.2788021   -0.09757891 0.9308379 15.185005
#> 6      6 3.926807 4.090117 0.3878004   -1.08192989 0.3497100 14.458337
#> ... and 36 more rows.
#> 
diag_fh <- diagnose(fit_fh)
print(diag_fh)
#> ── fastsae Small Area Estimation Diagnostic Report ─────────────────────────────
#> Model: "FH" | Domains: 42 (Sampled: 32, Unsampled: 10)
#> 
#> ── 1. Precision & Efficiency Gain ──
#> 
#> ! Domains with RSE < 25%: 52.4% (Caution: low precision)
#> ✔ Average RSE reduction: Direct "29.49%" -> SAE "25.59%" (Gain: 3.9%)
#> ✔ Variance reduction in 100% of areas (MSE ratio median: 1.33, max: 4.44)
#> 
#> ── 2. Brown et al. (2001) Diagnostics ──
#> 
#> ✔ Goodness-of-Fit Statistic W: W = 4.57 (df = 32), p-value = 1 [Consistent with direct estimates]
#> ℹ Bias regression (direct ~ model, descriptive): alpha = -0.723 (SE 0.3573), beta = 1.1587 (SE 0.0643), R^2 = 0.915
#> A slope above 1 is expected for shrinkage estimators; see ?diagnose.
#> ────────────────────────────────────────────────────────────────────────────────
#> ✖ Final Assessment: WARNING: High proportions of domains with RSE >= threshold.
#> ────────────────────────────────────────────────────────────────────────────────

# 2. Fit Spatial HB model (requires INLA)
if (requireNamespace("INLA", quietly = TRUE)) {
  fit_hb <- hb_area(y ~ x1 + x2, data = mys, vardir = "vardir")
  diag_hb <- diagnose(fit_hb)
  print(diag_hb)
}
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> hb_area(formula = y ~ x1 + x2, data = mys, vardir = "vardir")
#> 
#> ✔ Convergence: Yes (INLA mode status 0)
#> Model: HB-GAUSSIAN (Non-spatial)
#> Engine: inla
#> Random effect variance (sigma2_u): 2.829982 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  3.07155706  0.79098782  3.88319134  0.00010309  1.53306278
#> x1          -0.00415124  0.00943430 -0.44001600  0.65992552 -0.02298330
#> x2           0.08620897  0.03204007  2.69066112  0.00713106  0.02302954
#>             ci_upper
#> (Intercept)   4.6544
#> x1            0.0143
#> x2            0.1495
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred        sd       mse      rse ci_lower
#> 1      1 8.359527 7.590846    7.590846 0.7640224 0.5837303 10.06505 6.107741
#> 2      2 7.599650 6.775392    6.775392 0.8337948 0.6952138 12.30622 5.161056
#> 3      3 5.514137 5.245482    5.245482 0.8229080 0.6771775 15.68794 3.640202
#> 4      4 3.869326 4.247402    4.247402 0.7360698 0.5417987 17.32988 2.794586
#> 5      5 6.305063 6.352990    6.352990 0.9388749 0.8814860 14.77847 4.508182
#> 6      6 3.926807 4.091352    4.091352 0.5837763 0.3407948 14.26854 2.943485
#>   ci_upper random_effect    vardir
#> 1 9.103078    2.96522702 0.6618838
#> 2 8.430372    2.52170388 0.8374691
#> 3 6.869507    0.77265159 0.8822257
#> 4 5.682178   -1.47620938 0.6581716
#> 5 8.195014   -0.09886313 1.2788021
#> 6 5.233334   -1.08175985 0.3878004
#> ... and 36 more rows.
#> 
#> ── fastsae Small Area Estimation Diagnostic Report ─────────────────────────────
#> Model: "HB-GAUSSIAN (Non-spatial)" | Domains: 42 (Sampled: 32, Unsampled: 10)
#> 
#> ── 1. Precision & Efficiency Gain ──
#> 
#> ! Domains with RSE < 25%: 52.4% (Caution: low precision)
#> ✔ Average RSE reduction: Direct "29.49%" -> SAE "25.84%" (Gain: 3.65%)
#> ✔ Variance reduction in 100% of areas (MSE ratio median: 1.38, max: 4.4)
#> 
#> ── 2. Brown et al. (2001) Diagnostics ──
#> 
#> ✔ Goodness-of-Fit Statistic W: W = 4.51 (df = 32), p-value = 1 [Consistent with direct estimates]
#> ℹ Bias regression (direct ~ model, descriptive): alpha = -0.7264 (SE 0.3513), beta = 1.1591 (SE 0.0632), R^2 = 0.918
#> A slope above 1 is expected for shrinkage estimators; see ?diagnose.
#> 
#> ── 5. Bayesian Information Criteria & Predictive Diagnostics ──
#> 
#> ℹ WAIC: 114.2 (p_eff: 13.09) | DIC: 118.64 (p_eff: 22.2)
#> ✔ PIT Calibration Test vs Uniform(0,1): D = 0.079, p-value = 0.9802 [Well-calibrated predictive distribution]
#> ! Leave-One-Out CPO: 2 observation(s) flagged with numerical caution.
#> ────────────────────────────────────────────────────────────────────────────────
#> ✖ Final Assessment: WARNING: High proportions of domains with RSE >= threshold.
#> ────────────────────────────────────────────────────────────────────────────────
```
