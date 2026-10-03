
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
Area Estimation (SAE) in R. It combines ultra-fast compiled **C++**
(`RcppArmadillo` + OpenMP) for frequentist EBLUP models with **INLA**
(Integrated Nested Laplace Approximations) for fast Bayesian
hierarchical models—achieving **50× to 6,400× speedups** over classical
packages (`sae`, `emdi`) and **\>230× speedups** over MCMC (`tipsae`),
with near-zero memory footprint and exact numerical equivalence.

------------------------------------------------------------------------

### Highlights

- ⚡ **Blazing Fast**: Optimized C++ Fisher-scoring and deterministic
  INLA Laplace approximations turn minutes of runtime into milliseconds.
- 🎯 **Exact Concordance**: Frequentist estimates match `sae` to machine
  tolerance; Bayesian posteriors match Stan HMC ($r = 0.989$).
- 🔀 **Unified Interface**: Common S3 workflow (`summary`, `coef`,
  `fitted`, `residuals`, `autoplot`, `compare_sae`) for both paradigms.
- 🧵 **Multi-Threaded Bootstrap**: OpenMP-accelerated parametric and
  non-parametric bootstrap MSE.
- 🗺️ **Comprehensive Spatial & Temporal**: SAR, Besag ICAR, BYM2,
  Spatial Lag, AR(1), and RW(1/2).
- 💾 **Ultra-Low Memory**: Avoids large dense matrix allocations,
  scaling effortlessly to $n > 10,000$ domains.

------------------------------------------------------------------------

### Supported Models

#### Frequentist EBLUP (C++ / OpenMP)

| Function | Model Level | Structure | MSE Estimation |
|:---|:---|:---|:---|
| `eblup_fh()` | Area-level | Standard Fay-Herriot with IID domain effects | Analytical (Prasad-Rao) |
| `eblup_sfh()` | Area-level | Spatial SAR(1) autoregressive effects | Analytical, Bootstrap |
| `eblup_stfh()` | Area-level | Spatio-temporal SAR(1) + AR(1) dynamics | Bootstrap |
| `eblup_bhf()` | Unit-level | Battese-Harter-Fuller nested error regression | Analytical, Bootstrap |

#### Bayesian Hierarchical Models (INLA)

| Function | Model Level | Distributions | Random Effects |
|:---|:---|:---|:---|
| `hb_area()` | Area-level | Gaussian, Beta, Binomial, Poisson, NegBin, Gamma | IID, BYM2, Besag ICAR, RW(1/2), AR(1) |
| `hb_unit()` | Unit-level | Gaussian | Unit-level nested error |

------------------------------------------------------------------------

### Performance at a Glance

#### 1. Frequentist Fay-Herriot Benchmark ($n = 1,000$ areas, 5 covariates)

| Task | `fastsae` | `sae` (Molina & Marhuenda) | `emdi` |
|:---|:--:|:--:|:--:|
| **Standard Fay-Herriot** | $\color{#2ea44f}{\mathbf{0.0015\text{ s}}}$ | 0.291 s (~190x slower) | 9.64 s (~6,400x slower) |
| **Spatial Fay-Herriot** | $\color{#2ea44f}{\mathbf{0.165\text{ s}}}$ | 12.60 s (~76x slower) | 8.69 s (~53x slower) |
| **Peak RAM Usage** | $\color{#2ea44f}{\mathbf{< 10\text{ MB}}}$ | ~400 MB | ~800 MB |

#### 2. Bayesian Beta SAE Benchmark ($n = 1,000$ areas, compared against Stan MCMC)

| Model Specification | `fastsae::hb_area` (INLA) | `tipsae::fit_sae` (Stan HMC) | Speedup & Efficiency |
|:---|:--:|:--:|:--:|
| **Standard Beta** | $\color{#2ea44f}{\mathbf{1.64\text{ s (134 MB)}}}$ | 187.16 s (18.3 GB) | $\color{#2ea44f}{\mathbf{\sim 114\times\text{ faster, } >135\times\text{ less RAM}}}$ |
| **Spatial Beta (Besag ICAR)** | $\color{#2ea44f}{\mathbf{1.73\text{ s (161 MB)}}}$ | 412.39 s (22.9 GB) | $\color{#2ea44f}{\mathbf{\sim 238\times\text{ faster, } >140\times\text{ less RAM}}}$ |
| **Concordance vs MCMC** | $\color{#2ea44f}{\mathbf{r = 0.9893}}$ | Baseline | $\color{#2ea44f}{\mathbf{\text{MAE} = 0.003\text{ (Exact match)}}}$ |

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
[**https://ridsonap.github.io/fastsae/**](https://ridsonap.github.io/fastsae/)

- [Getting Started with
  fastsae](https://ridsonap.github.io/fastsae/articles/getting-started.html)
- [Spatial & Spatio-Temporal
  Models](https://ridsonap.github.io/fastsae/articles/spatial-temporal.html)
- [Bayesian Area-Level Models via
  INLA](https://ridsonap.github.io/fastsae/articles/hb-area-inla.html)
- [Model Diagnostics &
  Calibration](https://ridsonap.github.io/fastsae/articles/model-diagnostics.html)
- [Interactive
  Benchmarks](https://ridsonap.github.io/fastsae/articles/benchmarks.html)

------------------------------------------------------------------------

## Citation

``` r
citation("fastsae")
```

Bug reports and feature requests are welcome on the [issue
tracker](https://github.com/ridsonap/fastsae/issues).
