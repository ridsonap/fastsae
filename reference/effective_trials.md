# Effective Sample Sizes of Direct Proportions

Effective sample sizes \\n_d / \mathrm{deff}\_d\\ for binomial
area-level models (`trials` of
[`hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md) and
[`hb_twofold`](https://ridsonap.is-a.dev/fastsae/reference/hb_twofold.md)),
from the direct proportions, their design variances and the sample
sizes. The design effect of an area is undefined when its direct
proportion is 0 or 1 (no design variance); such areas take a design
effect pooled over their group, or a supplied fallback.

## Usage

``` r
effective_trials(
  y,
  vardir,
  n,
  group = NULL,
  data = NULL,
  method = c("pooled", "area"),
  deff_fallback = NULL,
  deff_min = 1
)
```

## Arguments

- y:

  Vector, column name, or one-sided formula: direct proportions in
  \\\[0, 1\]\\.

- vardir:

  Vector, column name, or one-sided formula: design variances of `y`.

- n:

  Vector, column name, or one-sided formula: sample sizes (number of
  sampled units in the denominator of `y`).

- group:

  Optional vector, column name, or one-sided formula: groups (e.g.
  provinces) within which design effects are pooled.

- data:

  Optional data frame containing the variables.

- method:

  `"pooled"` (default): every area takes the design effect of its group,
  \\\sum_d v_d / \sum_d \\y_d (1 - y_d) / n_d\\\\ over the areas of the
  group with \\0 \< y_d \< 1\\ (the ratio of the summed design variances
  to the summed simple random sampling variances). `"area"`: each area
  takes its own design effect \\v_d / \\y_d (1 - y_d) / n_d\\\\ where it
  is defined, and the pooled one otherwise.

- deff_fallback:

  Optional vector, column name, or one-sided formula: design effects for
  the areas whose own design effect is undefined, used before the pooled
  ones with `method = "area"` (for instance the design effect of a
  related indicator of the same survey, which is defined where `y` is
  0).

- deff_min:

  Lower bound of the design effects (default 1, no gain over simple
  random sampling).

## Value

A data frame with one row per area and columns `n`, `deff` (the design
effect used, at least `deff_min`), `deff_source` (`"area"`,
`"fallback"`, `"group"` or `"overall"`) and `trials` (\\n_d /
\mathrm{deff}\_d\\). Areas with missing `y` or `n` get `NA`.

## Details

Single-area design effects are noisy and undefined at \\y_d \in \\0,
1\\\\, which is common for rare indicators. In a simulation of a
two-stage cluster sample (12 provinces x 6 districts), effective sample
sizes with design effects pooled per province gave 95% interval coverage
of 0.925 (districts) and 0.937 (provinces), close to 0.929 and 0.946
with the true design effects, against 0.911 and 0.897 when the design
was ignored (`benchmarks/validation/binomial_neff.R` in the source
repository). A group without any area with \\0 \< y_d \< 1\\ takes the
design effect pooled over all areas.

## References

1.  Kish, L. (1965). *Survey Sampling*. John Wiley & Sons.

2.  Chen, C., Wakefield, J., and Lumley, T. (2014). The use of sampling
    weights in Bayesian hierarchical models for small area estimation.
    *Spatial and Spatio-temporal Epidemiology*, 11, 33–43.

## See also

[`hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
(`family = "binomial"`),
[`gvf_smooth`](https://ridsonap.is-a.dev/fastsae/reference/gvf_smooth.md).

## Examples

``` r
set.seed(1)
d <- data.frame(
  prov = rep(1:3, each = 5), n = rpois(15, 80),
  y = c(0, rbeta(14, 2, 30))
)
d$vardir <- 1.6 * d$y * (1 - d$y) / d$n
effective_trials(~y, vardir = ~vardir, n = ~n, group = ~prov, data = d)
#>     n deff deff_source trials
#> 1  74  1.6       group 46.250
#> 2  91  1.6       group 56.875
#> 3  91  1.6       group 56.875
#> 4  83  1.6       group 51.875
#> 5  66  1.6       group 41.250
#> 6  84  1.6       group 52.500
#> 7  86  1.6       group 53.750
#> 8  85  1.6       group 53.125
#> 9  77  1.6       group 48.125
#> 10 72  1.6       group 45.000
#> 11 74  1.6       group 46.250
#> 12 77  1.6       group 48.125
#> 13 79  1.6       group 49.375
#> 14 79  1.6       group 49.375
#> 15 88  1.6       group 55.000

# own design effects where defined, pooled ones for y = 0
effective_trials("y", "vardir", "n", group = "prov", data = d, method = "area")
#>     n deff deff_source trials
#> 1  74  1.6       group 46.250
#> 2  91  1.6        area 56.875
#> 3  91  1.6        area 56.875
#> 4  83  1.6        area 51.875
#> 5  66  1.6        area 41.250
#> 6  84  1.6        area 52.500
#> 7  86  1.6        area 53.750
#> 8  85  1.6        area 53.125
#> 9  77  1.6        area 48.125
#> 10 72  1.6        area 45.000
#> 11 74  1.6        area 46.250
#> 12 77  1.6        area 48.125
#> 13 79  1.6        area 49.375
#> 14 79  1.6        area 49.375
#> 15 88  1.6        area 55.000
```
