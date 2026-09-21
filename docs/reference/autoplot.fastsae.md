# Autoplot Method for fastsae Objects

Creates diagnostic and comparison plots for Small Area Estimation (SAE)
models fitted with fastsae. This extends the generic
[`autoplot`](https://ggplot2.tidyverse.org/reference/autoplot.html) from
ggplot2.

## Usage

``` r
# S3 method for class 'fastsae'
autoplot(object, type = c("comparison", "mse", "estimates", "scatter"), ...)

# S3 method for class 'list'
autoplot(object, type = c("comparison", "mse", "scatter"), ...)
```

## Arguments

- object:

  An object of class `fastsae`, or a (named) list of `fastsae` objects.

- type:

  Type of plot to create.

  - For a single `fastsae` object: `"comparison"`, `"mse"`,
    `"estimates"`, or `"scatter"`.

  - For a list of `fastsae` objects: `"comparison"`, `"mse"`, or
    `"scatter"`.

- ...:

  Additional arguments passed to internal plotting helpers (e.g.
  `title`) or to ggplot2 layers.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) { # \dontrun{
library(fastsae)

# Single model plot
fit_fh <- eblup_fh(y ~ x1 + x2 + x3, data = mys, vardir = "vardir")
autoplot(fit_fh, type = "estimates")

# Compare two models
fit_sfh <- eblup_sfh(y ~ x1 + x2 + x3, data = mys, vardir = "vardir", W = mys_proxmat)
autoplot(list("Fay-Herriot" = fit_fh, "Spatial FH" = fit_sfh), type = "comparison")
} # }
```
