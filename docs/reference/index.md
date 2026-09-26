# Package index

## Classical Area-Level Models (EBLUP)

High-performance C++ implementation of Fay-Herriot models

- [`eblup_fh()`](https://ridsonap.github.io/fastsae/reference/eblup_fh.md)
  : Empirical Best Linear Unbiased Prediction based on a Fay-Herriot
  Model.
- [`eblup_sfh()`](https://ridsonap.github.io/fastsae/reference/eblup_sfh.md)
  : Empirical Best Linear Unbiased Prediction based on a Spatial
  Fay-Herriot Model.
- [`eblup_stfh()`](https://ridsonap.github.io/fastsae/reference/eblup_stfh.md)
  : Empirical Best Linear Unbiased Prediction based on a Spatio-Temporal
  Fay-Herriot Model.

## Bayesian Area-Level Models (INLA / EBP)

Hierarchical Bayesian Small Area Estimation using INLA across 6
distributions with spatial and spatio-temporal effects

- [`ebp_area()`](https://ridsonap.github.io/fastsae/reference/ebp_area.md)
  : Empirical Best Prediction for Area-Level Small Area Estimation

## Unit-Level Models

Battese-Harter-Fuller model for unit-level survey data

- [`eblup_bhf()`](https://ridsonap.github.io/fastsae/reference/eblup_bhf.md)
  : Empirical Best Linear Unbiased Prediction (EBLUP) for the
  Battese-Harter-Fuller Model

## Model Diagnostics & Visualization

Universal diagnostic framework and ggplot2 visualization methods

- [`diagnose()`](https://ridsonap.github.io/fastsae/reference/diagnose.md)
  : Diagnostic Evaluation of Small Area Estimation Models
- [`autoplot(`*`<fastsae_diagnose>`*`)`](https://ridsonap.github.io/fastsae/reference/autoplot.fastsae_diagnose.md)
  : Autoplot Method for fastsae_diagnose Objects
- [`autoplot(`*`<fastsae>`*`)`](https://ridsonap.github.io/fastsae/reference/autoplot.fastsae.md)
  [`autoplot(`*`<list>`*`)`](https://ridsonap.github.io/fastsae/reference/autoplot.fastsae.md)
  : Autoplot Method for fastsae Objects
- [`summary(`*`<fastsae>`*`)`](https://ridsonap.github.io/fastsae/reference/summary.fastsae.md)
  : Summarize a fastsae object
- [`coef(`*`<fastsae>`*`)`](https://ridsonap.github.io/fastsae/reference/coef.fastsae.md)
  : Extract coefficients from a fastsae object
- [`fitted(`*`<fastsae>`*`)`](https://ridsonap.github.io/fastsae/reference/fitted.fastsae.md)
  : Extract fitted values (EBLUP or EBP) from a fastsae object
- [`residuals(`*`<fastsae>`*`)`](https://ridsonap.github.io/fastsae/reference/residuals.fastsae.md)
  : Extract residuals from a fastsae object
- [`print(`*`<fastsae>`*`)`](https://ridsonap.github.io/fastsae/reference/print.fastsae.md)
  : Print a fastsae object
- [`print(`*`<summary.fastsae>`*`)`](https://ridsonap.github.io/fastsae/reference/print.summary.fastsae.md)
  : Print summary of a fastsae object

## Simulation & Synthetic Data Generation

Utilities to simulate area-level cross-sectional, panel, and spatial
weight matrices

- [`sim_area_data()`](https://ridsonap.github.io/fastsae/reference/sim_area_data.md)
  : Simulate Multi-Distribution Small Area Data with Spatial
  Autocorrelation
- [`sim_series_data()`](https://ridsonap.github.io/fastsae/reference/sim_series_data.md)
  : Simulate Spatio-Temporal Multi-Distribution Small Area Data
- [`sim_spatial_weights()`](https://ridsonap.github.io/fastsae/reference/sim_spatial_weights.md)
  : Simulate Spatial Proximity and Adjacency Matrices

## Example & Synthetic Datasets

Built-in empirical and synthetic survey datasets

- [`mys`](https://ridsonap.github.io/fastsae/reference/mys.md) : mys:
  mean years of schooling people with disabilities.
- [`mys_panel`](https://ridsonap.github.io/fastsae/reference/mys_panel.md)
  : mys_panel: mean years of schooling people with disabilities 2016 -
  2026.
- [`mys_proxmat`](https://ridsonap.github.io/fastsae/reference/mys_proxmat.md)
  : Example proximity matrix
- [`cornsoybean`](https://ridsonap.github.io/fastsae/reference/cornsoybean.md)
  : Corn and Soybean Survey and Satellite Data in 12 Iowa Counties
- [`cornsoybeanmeans`](https://ridsonap.github.io/fastsae/reference/cornsoybeanmeans.md)
  : Corn and Soybean Mean Number of Pixels per Segment for 12 Iowa
  Counties
- [`sim_area`](https://ridsonap.github.io/fastsae/reference/sim_area.md)
  : Multi-Distribution Synthetic Area-Level Dataset
- [`sim_panel`](https://ridsonap.github.io/fastsae/reference/sim_panel.md)
  : Multi-Distribution Spatio-Temporal Panel Dataset
