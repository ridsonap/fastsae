# Plot method for fastsae objects

Provides standard R
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) dispatch for
`fastsae` models. Automatically routes to
[`map_sae`](https://ridsonap.is-a.dev/fastsae/reference/map_sae.md) if
spatial geometry is supplied, to two-model comparison if a second
`fastsae` model is provided, or to
[`autoplot`](https://ggplot2.tidyverse.org/reference/autoplot.html) for
diagnostic and estimate plots.

## Usage

``` r
# S3 method for class 'fastsae'
plot(x, y = NULL, ...)
```

## Arguments

- x:

  An object of class `fastsae`.

- y:

  Optional second object (e.g. another `fastsae` model for comparison).

- ...:

  Additional arguments passed to
  [`map_sae`](https://ridsonap.is-a.dev/fastsae/reference/map_sae.md) or
  [`autoplot`](https://ggplot2.tidyverse.org/reference/autoplot.html).

## Value

A `ggplot` object.

## Examples

``` r
m <- eblup_fh(y ~ x1 + x2 + x3, vardir = "vardir", data = mys, print_result = FALSE)
plot(m)
```
