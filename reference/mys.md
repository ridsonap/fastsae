# mys: mean years of schooling of people with disabilities

A synthetic area-level dataset with direct estimates of the mean years
of schooling of people with disabilities in 42 regencies/municipalities,
10 of which are not sampled (`y = NA`). Used in the examples of the
area-level models; its spatial weights are
[`mys_proxmat`](https://ridsonap.is-a.dev/fastsae/reference/mys_proxmat.md).

## Usage

``` r
mys
```

## Format

A data frame (tibble) with 42 rows and 9 variables:

- area:

  area (regency/municipality) identifier, 1 to 42.

- y:

  direct estimate of the mean years of schooling (`NA` for the 10
  unsampled areas).

- vardir:

  sampling variance of the direct estimate (0 for the unsampled areas,
  where it is not used).

- rse:

  relative standard error of the direct estimate (%),
  \\100\sqrt{vardir}/y\\ (0 for the unsampled areas).

- x1:

  number of elementary schools.

- x2:

  number of junior high schools.

- x3:

  number of senior high schools.

- n:

  number of people with disabilities in the Susenas sample of the area
  (eligible sample size; synthetic, hence not an integer, and positive
  also for the unsampled areas). `vardir` was not derived from it.

- weight:

  sample share of the area, `n / sum(n)`.

## Source

Simulated based on empirical characteristics of the National
Socio-Economic Survey (Susenas), BPS-Statistics Indonesia.

## Examples

``` r
data(mys)
head(mys)
#> # A tibble: 6 × 9
#>    area     y vardir   rse    x1    x2    x3      n  weight
#>   <int> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>   <dbl>
#> 1     1  8.36  0.662  9.73   124    24    14  7280. 0.0326 
#> 2     2  7.60  0.837 12.0     89    18     9  2743. 0.0123 
#> 3     3  5.51  0.882 17.0     57    19     5  1706. 0.00764
#> 4     4  3.87  0.658 21.0     88    35    19  3073. 0.0138 
#> 5     5  6.31  1.28  17.9    141    46    29 13400. 0.0600 
#> 6     6  3.93  0.388 15.9     96    29    10  2004. 0.00897
# unsampled areas
mys$area[is.na(mys$y)]
#>  [1] 21 25 26 27 28 34 35 36 37 41
```
