# Benchmark Small Area Estimation Predictions to Aggregate Targets

Calibrates small area model predictions (from `hb_area`, `eblup_fh`,
`eblup_sfh`, `eblup_stfh`, or `eblup_bhf`) so that their weighted
aggregate matches known benchmark targets at higher administrative
levels (e.g. provincial or national totals/means), as required in
official statistics production.

Supports:

- **Single-Stage Benchmarking**: Calibrates small areas directly to an
  aggregate target (either globally or within sub-regions/provinces via
  `group`).

- **Two-Stage Hierarchical Benchmarking** (Rao and Molina, 2015, Sec
  10.4; Steorts et al., 2014): Calibrates small areas (e.g. 500
  Kabupaten/Kota) within overarching groups (e.g. 34 Provinces) while
  *simultaneously* harmonizing the provincial aggregates to match an
  overarching National Target (`national_target`). Guarantees
  simultaneous mathematical consistency across all administrative tiers:
  \\\sum\_{d \in \text{Prov}\_g} \tilde{w}\_{dg}
  \hat{\theta}\_{dg}^{\text{BM}} = T_g^\*\\ and \\\sum_g \tilde{W}\_g
  T_g^\* = T\_{\text{nat}}\\.

Methods implemented:

- `"ratio"`: Proportional / multiplicative adjustment. Preserves area
  proportions and non-negativity (You and Rao, 2002; Rao and Molina,
  2015).

- `"difference"`: Uniform additive adjustment (Rao and Molina, 2015,
  Section 10.3).

- `"optimal"`: Variance-weighted quadratic loss benchmarking (Datta et
  al., 2011; Steorts et al., 2014). Domains with higher MSE absorb
  larger adjustments while precise domains remain minimally perturbed.

- `"logit"`: Logit-scale additive shift for bounded indicators
  \\\theta_d \in (0, 1)\\ (Berg and Fuller, 2014).

## Usage

``` r
benchmark_sae(object, ...)

# S3 method for class 'fastsae'
benchmark_sae(
  object,
  target = NULL,
  weight = NULL,
  method = c("ratio", "difference", "optimal", "logit"),
  group = NULL,
  type = c("mean", "total"),
  national_target = NULL,
  outer_target = NULL,
  ...
)

# Default S3 method
benchmark_sae(
  object,
  target = NULL,
  weight = NULL,
  method = c("ratio", "difference", "optimal", "logit"),
  group = NULL,
  type = c("mean", "total"),
  national_target = NULL,
  outer_target = NULL,
  ...
)

benchmark(object, ...)

# S3 method for class 'fastsae_benchmark'
print(x, all_groups = FALSE, ...)

# S3 method for class 'fastsae_benchmark'
summary(object, ...)

# S3 method for class 'fastsae_benchmark'
autoplot(object, ...)

# S3 method for class 'fastsae_benchmark'
plot(x, y = NULL, ...)
```

## Arguments

- object:

  A fitted `fastsae` model object (e.g., from `hb_area`, `eblup_fh`,
  etc.) or a numeric vector of model predictions.

- ...:

  Additional arguments passed to methods.

- target:

  Numeric value, named vector, or data frame specifying the benchmark
  target(s) \\T_g\\. If `group` is specified, `target` can be:

  - A single scalar numeric (constant target for all groups).

  - A named numeric vector matching group identifiers.

  - A data frame with columns `group` and `target`.

  - `NULL` if `national_target` is provided in two-stage hierarchical
    mode (initial provincial targets are automatically derived from
    model aggregates).

- weight:

  Optional numeric vector or character string naming the domain
  benchmark weight column in `object$data` (e.g., population sizes
  \\N_d\\ or population shares \\N_d / \sum N_j\\). Defaults to equal
  weights across domains if `NULL`.

- method:

  Character string specifying the benchmarking calibration method:
  `"ratio"` (default), `"difference"`, `"optimal"`, or `"logit"`.

- group:

  Optional vector or character string naming the grouping/stratum column
  (e.g., province or region) for multi-level hierarchical calibration.

- type:

  Character string: `"mean"` (default, benchmark target represents
  weighted average/rate, \\\sum\_{d \in g} w_d \hat{\theta}\_d = T_g\\
  with normalized weights) or `"total"` (benchmark target represents
  population total, \\\sum\_{d \in g} w_d \hat{\theta}\_d = T_g\\).

- national_target:

  Optional numeric scalar defining the overarching national target
  (Level 0) for Two-Stage Hierarchical Benchmarking. When provided
  alongside `group`, first calibrates/harmonizes group (provincial)
  targets to match `national_target`, then calibrates small areas to the
  harmonized provincial targets.

- outer_target:

  Alias for `national_target`.

- x:

  An object of class `fastsae_benchmark` (for `print` and `plot`
  methods).

- all_groups:

  Logical; if `TRUE`, prints verification status for all groups in
  hierarchical benchmarking. Default is `FALSE`.

- y:

  Ignored argument for compatibility with the generic `plot` method.

## Value

An object of class `c("fastsae_benchmark", "data.frame")` containing:

- `domain`: Domain identifier.

- `group`: Stratum / overarching group identifier (if specified).

- `weight`: Benchmark weight \\w_d\\.

- `original`: Model-based prediction before benchmarking
  (\\\hat{\theta}\_d\\).

- `benchmarked`: Calibrated estimate satisfying the benchmark constraint
  (\\\hat{\theta}\_d^{\text{BM}}\\).

- `adjustment`: Absolute adjustment (\\\hat{\theta}\_d^{\text{BM}} -
  \hat{\theta}\_d\\).

- `rel_adjustment_pct`: Relative adjustment percentage (%).

- `target`: Corresponding group benchmark target \\T_g^\*\\.

- `national_target`: Overarching national target (if two-stage
  hierarchical).

## References

1.  Rao, J. N. K., and Molina, I. (2015). *Small Area Estimation* (2nd
    ed.). John Wiley & Sons. Chapter 10: "Benchmarking and Other
    Practical Issues", pp. 297-315.

2.  You, Y., and Rao, J. N. K. (2002). A pseudo-empirical best linear
    unbiased prediction approach to small area estimation using survey
    weights. *The Canadian Journal of Statistics*, 30(3), 431-439.

3.  Steorts, R. C., Hall, P., and Ghosh, M. (2014). General benchmarking
    under quadratic loss with applications to small area estimation.
    *Journal of Survey Statistics and Methodology*, 2(2), 173-193.

4.  Datta, G. S., Ghosh, M., Steorts, R., and Maples, J. (2011).
    Bayesian benchmarking with applications to small area estimation.
    *Test*, 20(3), 574-588.

5.  Berg, E., and Fuller, W. A. (2014). Small area prediction of
    proportions with a constrained multinomial logit model. *Journal of
    Survey Statistics and Methodology*, 2(3), 256-283.

## Examples

``` r
library(fastsae)
data(mys)

# Assign province groups for hierarchical calibration example
mys$province <- rep(c("Prov_A", "Prov_B", "Prov_C"), length.out = nrow(mys))

# 1. Fit Fay-Herriot model
fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> eblup_fh(formula = y ~ x1 + x2, vardir = "vardir", data = mys)
#> 
#> ✔ Convergence: Yes (in 7 iterations)
#> Model: Fay-Herriot (Area-level)
#> Method: eblup
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

# 2. Ratio benchmarking to overall target (e.g. target = 6.5)
bm_ratio <- benchmark_sae(fit_fh, target = 6.5, method = "ratio")
head(bm_ratio)
#> === fastsae Small Area Benchmark Calibration ===
#> Method: Ratio (Multiplicative Adjustment) 
#> Target Type: Weighted Mean / Rate 
#> Total Domains: 6 
#> 
#> Level 1 (Group Consistency): 1 group
#>   Target: Target = 6.50000 | Original = 5.09590 -> Calibrated = 6.50000 [CONSISTENT]
#> 
#> First 6 benchmarked domains (Level 2):
#>   domain weight original benchmarked adjustment rel_adjustment_pct target
#> 1      1      1 7.594807    9.687449   2.092641           27.55358    6.5
#> 2      2      1 6.777056    8.644377   1.867322           27.55358    6.5
#> 3      3      1 5.247512    6.693390   1.445878           27.55358    6.5
#> 4      4      1 4.247072    5.417293   1.170221           27.55358    6.5
#> 5      5      1 6.353632    8.104286   1.750653           27.55358    6.5
#> 6      6      1 4.090117    5.217090   1.126974           27.55358    6.5
#> 

# 3. Two-Stage Hierarchical Calibration (e.g. Regency -> Province -> National)
# Calibrate domains within provinces while matching overarching national target:
bm_hier <- benchmark_sae(
  fit_fh,
  group = "province",
  weight = mys$n,
  national_target = 6.2,
  method = "ratio"
)
print(bm_hier)
#> === fastsae Two-Stage Hierarchical Benchmark Calibration ===
#> Method: Ratio (Multiplicative Adjustment) 
#> Target Type: Weighted Mean / Rate 
#> Total Domains: 42 
#> 
#> Level 0 (National Target): Target = 6.20000 | Calibrated Aggregate = 6.20000 [CONSISTENT]
#> 
#> Level 1 (Group Consistency): 3 groups
#>   Group 'Prov_A': Target = 6.83687 | Original = 6.30670 -> Calibrated = 6.83687 [CONSISTENT]
#>   Group 'Prov_B': Target = 6.32663 | Original = 5.83603 -> Calibrated = 6.32663 [CONSISTENT]
#>   Group 'Prov_C': Target = 5.47339 | Original = 5.04896 -> Calibrated = 5.47339 [CONSISTENT]
#> 
#> First 6 benchmarked domains (Level 2):
#>   domain  group    weight original benchmarked adjustment rel_adjustment_pct
#> 1      1 Prov_A  7279.992 7.594807    8.233252  0.6384450           8.406336
#> 2      2 Prov_B  2743.181 6.777056    7.346758  0.5697021           8.406336
#> 3      3 Prov_C  1706.433 5.247512    5.688636  0.4411235           8.406336
#> 4      4 Prov_A  3073.056 4.247072    4.604096  0.3570232           8.406336
#> 5      5 Prov_B 13400.228 6.353632    6.887740  0.5341077           8.406336
#> 6      6 Prov_C  2004.006 4.090117    4.433946  0.3438290           8.406336
#>     target national_target
#> 1 6.836868             6.2
#> 2 6.326627             6.2
#> 3 5.473390             6.2
#> 4 6.836868             6.2
#> 5 6.326627             6.2
#> 6 5.473390             6.2
#> ... and 36 more domains.
#> 
```
