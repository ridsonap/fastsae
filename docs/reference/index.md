# Package index

## Area-Level Models

Fay-Herriot family of models for area-level small area estimation

- [`eblup_fh()`](https://ridsonap.github.io/fastsae/reference/eblup_fh.md)
  : EBLUPs based on a Fay-Herriot Model.
- [`eblup_sfh()`](https://ridsonap.github.io/fastsae/reference/eblup_sfh.md)
  : EBLUPs based on a Spatial Fay-Herriot Model.
- [`eblup_stfh()`](https://ridsonap.github.io/fastsae/reference/eblup_stfh.md)
  : EBLUPs based on a Spatio-Temporal Fay-Herriot Model.

## Unit-Level Models

Battese-Harter-Fuller model for unit-level survey data

- [`eblup_bhf()`](https://ridsonap.github.io/fastsae/reference/eblup_bhf.md)
  [`pbmse_unit()`](https://ridsonap.github.io/fastsae/reference/eblup_bhf.md)
  : Empirical Best Linear Unbiased Prediction (EBLUP) for the
  Battese-Harter-Fuller Model

## Model Diagnostics & Visualization

S3 methods and ggplot2-based diagnostic plots

- [`autoplot(`*`<fastsae>`*`)`](https://ridsonap.github.io/fastsae/reference/autoplot.fastsae.md)
  [`autoplot(`*`<list>`*`)`](https://ridsonap.github.io/fastsae/reference/autoplot.fastsae.md)
  : Autoplot Method for fastsae Objects
- [`summary(`*`<fastsae>`*`)`](https://ridsonap.github.io/fastsae/reference/summary.fastsae.md)
  : Summarize a fastsae object
- [`coef(`*`<fastsae>`*`)`](https://ridsonap.github.io/fastsae/reference/coef.fastsae.md)
  : Extract coefficients from a fastsae object
- [`fitted(`*`<fastsae>`*`)`](https://ridsonap.github.io/fastsae/reference/fitted.fastsae.md)
  : Extract fitted values (EBLUP) from a fastsae object
- [`residuals(`*`<fastsae>`*`)`](https://ridsonap.github.io/fastsae/reference/residuals.fastsae.md)
  : Extract residuals from a fastsae object
- [`print(`*`<fastsae>`*`)`](https://ridsonap.github.io/fastsae/reference/print.fastsae.md)
  : Print a fastsae object
- [`print(`*`<summary.fastsae>`*`)`](https://ridsonap.github.io/fastsae/reference/print.summary.fastsae.md)
  : Print summary of a fastsae object

## Example Datasets

Built-in survey datasets and spatial proximity matrices

- [`mys`](https://ridsonap.github.io/fastsae/reference/mys.md) : mys:
  mean years of schooling people with disabilities.
- [`mys_panel`](https://ridsonap.github.io/fastsae/reference/mys_panel.md)
  : mys: mean years of schooling people with disabilities disabilities
  2016 - 2026.
- [`mys_proxmat`](https://ridsonap.github.io/fastsae/reference/mys_proxmat.md)
  : Example proximity matrix
- [`cornsoybean`](https://ridsonap.github.io/fastsae/reference/cornsoybean.md)
  : Corn and Soybean Survey and Satellite Data in 12 Iowa Counties
- [`cornsoybeanmeans`](https://ridsonap.github.io/fastsae/reference/cornsoybeanmeans.md)
  : Corn and Soybean Mean Number of Pixels per Segment for 12 Iowa
  Counties
