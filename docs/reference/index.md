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

## Bayesian Area-Level Models (INLA)

Hierarchical Bayesian SAE using Integrated Nested Laplace Approximations
(INLA) across 6 distributions with spatial and spatio-temporal effects

- [`hb_area()`](https://ridsonap.github.io/fastsae/reference/hb_area.md)
  : Hierarchical Bayes for Area-Level Small Area Estimation

## Two-Fold Subarea Models

Nested area-subarea level models (Torabi & Rao, 2014) with both area and
subarea random effects via REML EBLUP and Bayesian INLA

- [`eblup_twofold()`](https://ridsonap.github.io/fastsae/reference/eblup_twofold.md)
  : Empirical Best Linear Unbiased Prediction under a Two-fold
  Fay-Herriot Model.
- [`hb_twofold()`](https://ridsonap.github.io/fastsae/reference/hb_twofold.md)
  : Two-Fold Hierarchical Bayes for Sub-Area Level Small Area Estimation

## Unit-Level Models

Battese-Harter-Fuller models for unit-level survey data via REML EBLUP
and Bayesian INLA

- [`eblup_bhf()`](https://ridsonap.github.io/fastsae/reference/eblup_bhf.md)
  : Empirical Best Linear Unbiased Prediction (EBLUP) for the
  Battese-Harter-Fuller Model
- [`hb_unit()`](https://ridsonap.github.io/fastsae/reference/hb_unit.md)
  : Hierarchical Bayes for Unit-Level Small Area Estimation

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
  : Extract fitted values (EBLUP or HB) from a fastsae object
- [`residuals(`*`<fastsae>`*`)`](https://ridsonap.github.io/fastsae/reference/residuals.fastsae.md)
  : Extract residuals from a fastsae object
- [`print(`*`<fastsae>`*`)`](https://ridsonap.github.io/fastsae/reference/print.fastsae.md)
  : Print a fastsae object
- [`print(`*`<summary.fastsae>`*`)`](https://ridsonap.github.io/fastsae/reference/print.summary.fastsae.md)
  : Print summary of a fastsae object

## Benchmarking & Calibration

Calibrate small area estimates to aggregate targets via Ratio,
Difference, Optimal, and Logit methods (Rao & Molina, 2015)

- [`benchmark_sae()`](https://ridsonap.github.io/fastsae/reference/benchmark_sae.md)
  [`benchmark()`](https://ridsonap.github.io/fastsae/reference/benchmark_sae.md)
  [`print(`*`<fastsae_benchmark>`*`)`](https://ridsonap.github.io/fastsae/reference/benchmark_sae.md)
  [`summary(`*`<fastsae_benchmark>`*`)`](https://ridsonap.github.io/fastsae/reference/benchmark_sae.md)
  [`autoplot(`*`<fastsae_benchmark>`*`)`](https://ridsonap.github.io/fastsae/reference/benchmark_sae.md)
  [`plot(`*`<fastsae_benchmark>`*`)`](https://ridsonap.github.io/fastsae/reference/benchmark_sae.md)
  : Benchmark Small Area Estimation Predictions to Aggregate Targets

## Model Comparison & Spatial Mapping

Unified spatial choropleth mapping (‘satu pintu’), model concordance
evaluation, and official publication exports

- [`map_sae()`](https://ridsonap.github.io/fastsae/reference/map_sae.md)
  : Unified Spatial Choropleth Mapping for Small Area Estimates
- [`compare_sae()`](https://ridsonap.github.io/fastsae/reference/compare_sae.md)
  [`print(`*`<fastsae_comparison>`*`)`](https://ridsonap.github.io/fastsae/reference/compare_sae.md)
  [`summary(`*`<fastsae_comparison>`*`)`](https://ridsonap.github.io/fastsae/reference/compare_sae.md)
  [`plot(`*`<fastsae_comparison>`*`)`](https://ridsonap.github.io/fastsae/reference/compare_sae.md)
  [`autoplot(`*`<fastsae_comparison>`*`)`](https://ridsonap.github.io/fastsae/reference/compare_sae.md)
  : Compare Two Small Area Estimation Models
- [`export_sae()`](https://ridsonap.github.io/fastsae/reference/export_sae.md)
  : Export Small Area Estimation Results to Excel or CSV
- [`add_reliability_flags()`](https://ridsonap.github.io/fastsae/reference/add_reliability_flags.md)
  : Add Official Statistical Reliability Flags
- [`plot(`*`<fastsae>`*`)`](https://ridsonap.github.io/fastsae/reference/plot.fastsae.md)
  : Plot method for fastsae objects

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
