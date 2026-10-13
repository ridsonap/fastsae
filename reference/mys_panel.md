# mys_panel: mean years of schooling people with disabilities 2016 - 2026.

A synthetic panel dataset containing the mean years of schooling for
people with disabilities across 42 regencies/municipalities over 11 time
periods (2016 - 2026).

## Usage

``` r
mys_panel
```

## Format

A data frame with 462 rows and 10 variables.

- area:

  regency municipality identifier

- year:

  year of observation (2016 to 2026)

- y:

  mean years of schooling people with disabilities (NA for unsampled
  areas)

- vardir:

  variance sampling from the direct estimator for each area

- rse:

  relative standard error (%)

- x1:

  Number of Elementary Schools

- x2:

  Number of Junior High Schools

- x3:

  Number of Senior High Schools

- n:

  Number of eligible samples

- weight:

  Weight

## Source

Simulated panel data based on empirical characteristics of the National
Socio-Economic Survey (Susenas), BPS-Statistics Indonesia.
