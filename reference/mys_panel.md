# mys_panel: mean years of schooling of people with disabilities, 2016-2026

A synthetic panel dataset with direct estimates of the mean years of
schooling of people with disabilities in the 42 areas of
[`mys`](https://ridsonap.is-a.dev/fastsae/reference/mys.md) over 11
years (2016-2026). Areas 21 and 25 are never sampled (`y = NA` in every
year). The rows are in time-major order (all areas of 2016, then 2017,
...);
[`eblup_stfh`](https://ridsonap.is-a.dev/fastsae/reference/eblup_stfh.md)
sorts them internally, but
[`sae::eblupSTFH()`](https://rdrr.io/pkg/sae/man/eblupSTFH.html) needs
domain-major order.

## Usage

``` r
mys_panel
```

## Format

A data frame (tibble) with 462 rows and 10 variables:

- area:

  area identifier, 1 to 42.

- year:

  year of observation (2016 to 2026).

- y:

  direct estimate of the mean years of schooling (`NA` for areas 21 and
  25).

- vardir:

  sampling variance of the direct estimate (`NA` where `y` is `NA`).

- rse:

  relative standard error of the direct estimate (%).

- x1:

  number of elementary schools.

- x2:

  number of junior high schools.

- x3:

  number of senior high schools.

- n:

  number of people with disabilities in the Susenas sample of the area
  in that year (eligible sample size; synthetic, hence not an integer,
  and positive also for areas 21 and 25). `vardir` was not derived from
  it.

- weight:

  sample share, `n / sum(n)` over all rows of the panel (the weights sum
  to 1 over the panel, not within a year).

## Source

Simulated panel data based on empirical characteristics of the National
Socio-Economic Survey (Susenas), BPS-Statistics Indonesia.

## Examples

``` r
data(mys_panel)
head(mys_panel)
#> # A tibble: 6 × 10
#>    area  year     y vardir   rse    x1    x2    x3      n   weight
#>   <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>    <dbl>
#> 1     1  2016  8.26  0.662  9.85 123.   24.8 12.9   7229. 0.00287 
#> 2     2  2016  7.48  0.837 12.2   93.2  17.8  9.14  2767. 0.00110 
#> 3     3  2016  5.55  0.882 16.9   54.0  18.0  5.10  1693. 0.000671
#> 4     4  2016  3.43  0.658 23.7   87.4  35.8 17.8   3140. 0.00125 
#> 5     5  2016  6.46  1.28  17.5  138.   45.7 27.1  13421. 0.00532 
#> 6     6  2016  3.61  0.388 17.3   95.7  29.2 10.6   1984. 0.000787
table(sampled = !is.na(mys_panel$y), mys_panel$year)[, 1:3]
#>        
#> sampled 2016 2017 2018
#>   FALSE    2    2    2
#>   TRUE    40   40   40
```
