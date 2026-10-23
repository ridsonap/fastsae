
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/logo.png" alt="fastsae logo" align="right" height="139"/>

# fastsae: High-Performance Frequentist and Bayesian Small Area Estimation in R

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/fastsae)](https://CRAN.R-project.org/package=fastsae)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/fastsae)](https://CRAN.R-project.org/package=fastsae)
[![Monthly
Downloads](https://cranlogs.r-pkg.org/badges/fastsae)](https://CRAN.R-project.org/package=fastsae)
[![R-CMD-check](https://github.com/ridsonap/fastsae/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ridsonap/fastsae/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/ridsonap/fastsae/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ridsonap/fastsae)
[![License:
GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

<!-- badges: end -->

**fastsae** provides a high-performance, unified framework for Small
Area Estimation (SAE) in R. It combines compiled **C++**
(`RcppArmadillo`, OpenMP, sparse Cholesky factorisations for the spatial
models) for frequentist EBLUP models with **INLA** (Integrated Nested
Laplace Approximations) for Bayesian hierarchical models. It is orders
of magnitude faster than classical packages (`sae`, `emdi`) and MCMC
(`tipsae`) at realistic numbers of domains (see the benchmarks below),
and reproduces the `sae` estimates where both use the same method.

------------------------------------------------------------------------

### Highlights

- ⚡ **Blazing Fast**: Optimized C++ Fisher-scoring and deterministic
  INLA Laplace approximations turn minutes of runtime into milliseconds.
- 🎯 **Validated Against Reference Packages**: REML estimates and
  analytical MSEs agree with `sae` for interior estimates; deliberate,
  documented differences where `sae` is not valid (e.g. $\hat\rho$ on
  the boundary, ML information for ML). Bayesian posteriors agree
  closely with Stan HMC ($r = 0.989$).
- 🔀 **Unified Interface**: Common S3 workflow (`summary`, `coef`,
  `fitted`, `residuals`, `autoplot`, `compare_sae`) for both paradigms.
- 🧵 **Multi-Threaded Bootstrap**: OpenMP-accelerated parametric and
  non-parametric bootstrap MSE.
- 🗺️ **Comprehensive Spatial & Temporal**: SAR, Besag ICAR, BYM2,
  Spatial Lag, AR(1), and RW(1/2).
- 💾 **Low Memory**: Structured and sparse linear algebra instead of
  dense $D \times D$ inverses (spatial FH with a few thousand domains
  runs in seconds).

------------------------------------------------------------------------

### Supported Models

#### Frequentist EBLUP (C++ / OpenMP)

| Function | Model Level | Structure | MSE Estimation |
|:---|:---|:---|:---|
| `eblup_fh()` | Area-level | Standard Fay-Herriot with IID domain effects | Analytical (Prasad-Rao) |
| `eblup_sfh()` | Area-level | Spatial SAR(1) autoregressive effects | Analytical, Bootstrap |
| `eblup_stfh()` | Area-level | Spatio-temporal SAR(1) + AR(1) dynamics | Analytical, Bootstrap |
| `eblup_bhf()` | Unit-level | Battese-Harter-Fuller nested error regression | Bootstrap |

#### Bayesian Hierarchical Models (INLA)

| Function | Model Level | Distributions | Random Effects |
|:---|:---|:---|:---|
| `hb_area()` | Area-level | Gaussian, Beta, Binomial, Poisson, NegBin, Gamma | IID, BYM2, Besag ICAR, RW(1/2), AR(1) |
| `hb_unit()` | Unit-level | Gaussian | Unit-level nested error |

------------------------------------------------------------------------

### Performance at a Glance

<!-- BENCH-README:START -->

#### 1. Frequentist EBLUP with analytical MSE

| Model ($n = 1,000$ domains, 5 covariates) | `fastsae` | `sae` | `emdi` |
|:---|:--:|:--:|:--:|
| **Fay-Herriot** | **0.00119 s** (0.503 MB) | 1.50 s (72.0 MB; 1,263x slower) | 51.0 s (3.69 GB; 42,874x slower) |
| **Spatial Fay-Herriot** | **1.00 s** (0.203 MB) | 67.2 s (1.88 GB; 67x slower) | 51.1 s (3.68 GB; 51x slower) |
| **Spatio-Temporal FH ($T = 5$)** | **75.3 s** (3.29 MB) | 1790 s (31.7 GB; 24x slower) | \- |

#### 2. Bayesian Beta models: INLA vs Stan MCMC (largest $n$ benchmarked)

| Model | $n$ | `fastsae::hb_area` (INLA) | `tipsae::fit_sae` (Stan) | Speedup |
|:---|:--:|:--:|:--:|:--:|
| **Beta** | 1,000 | **1.61 s** (25.9 MB) | 187 s (17.9 GB) | 117x |
| **Spatial Beta (Besag ICAR)** | 1,000 | **1.87 s** (161 MB) | 404 s (22.5 GB) | 217x |
| **Spatio-Temporal Beta (Besag + RW1, $T = 5$)** | 50 | **1.71 s** (17.1 MB) | 908 s (4.77 GB) | 530x |

Times are medians of `bench::mark()`; memory is the total memory
allocated by R during one call (`bench::mark()$mem_alloc`), not the peak
resident memory. The `fastsae` rows were re-measured with fastsae 1.0.0
(frequentist on 2026-10-13, `benchmarks/bench_frequentist_fastsae.R`;
Bayesian on 2026-10-14, `benchmarks/bench_bayes_fastsae.R`); the
`sae`/`emdi` rows come from an earlier run of the same design (simulated
data, 5 covariates) on different random draws, the `tipsae` rows from an
earlier run on the same data. See the [benchmark
article](https://ridsonap.is-a.dev/fastsae/articles/benchmarks.html)
for all sizes.

<!-- BENCH-README:END -->

------------------------------------------------------------------------

## Installation

Install from **CRAN**:

``` r
install.packages("fastsae")
```

Or install the development version from **GitHub**:

``` r
# install.packages("remotes")
remotes::install_github("ridsonap/fastsae")
```

To enable Bayesian models, install **INLA**:

``` r
install.packages(
  "INLA",
  repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
  dep = TRUE
)
```

------------------------------------------------------------------------

## Quick Start

``` r
library(fastsae)

# 1. Frequentist Fay-Herriot (C++)
fit_fh <- eblup_fh(y ~ x1 + x2 + x3, vardir = ~vardir, data = na.omit(mys))
summary(fit_fh)

# 2. Diagnostic & Visualization
autoplot(fit_fh, type = "estimates")
```

------------------------------------------------------------------------

## 📖 Documentation & Tutorials

Explore comprehensive guides, mathematical formulations, and worked
examples:

👉
[**https://ridsonap.is-a.dev/fastsae/**](https://ridsonap.is-a.dev/fastsae/)

- [Getting Started with
  fastsae](https://ridsonap.is-a.dev/fastsae/articles/getting-started.html)
- [Spatial & Spatio-Temporal
  Models](https://ridsonap.is-a.dev/fastsae/articles/spatial-temporal.html)
- [Bayesian Area-Level Models via
  INLA](https://ridsonap.is-a.dev/fastsae/articles/hb-area-inla.html)
- [Model Diagnostics &
  Calibration](https://ridsonap.is-a.dev/fastsae/articles/model-diagnostics.html)
- [Interactive
  Benchmarks](https://ridsonap.is-a.dev/fastsae/articles/benchmarks.html)

------------------------------------------------------------------------

## Citation

``` r
citation("fastsae")
```

Bug reports and feature requests are welcome on the [issue
tracker](https://github.com/ridsonap/fastsae/issues).
