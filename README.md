
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

**fastsae** provides a unified framework for Small Area Estimation (SAE)
in R. It combines compiled **C++** (`RcppArmadillo`, OpenMP, sparse
Cholesky factorisations for the spatial models) for frequentist EBLUP
models with **INLA** (Integrated Nested Laplace Approximations) for
Bayesian hierarchical models, and can also fit the Bayesian models by
MCMC (NumPyro on JAX). With hundreds to thousands of domains it is much
faster than `sae`, `emdi` and Stan-based `tipsae` (see the benchmarks
below; with a few dozen domains `sae` is about as fast), and it
reproduces the `sae` estimates where both use the same method.

------------------------------------------------------------------------

### Highlights

- ⚡ **Fast**: compiled Fisher scoring for the EBLUP models and INLA for
  the Bayesian models; the spatial Fay-Herriot model with analytical MSE
  for 1,000 domains takes well under a second.
- 🎯 **Validated against reference packages**: on the same data, the
  Fay-Herriot, spatial, spatio-temporal and unit-level EBLUPs, variance
  components and analytical MSEs agree with `sae` to about $10^{-14}$
  (unit level: $10^{-9}$). The two deliberate differences are
  documented: a spatial correlation $\hat\rho$ on the boundary is
  treated as fixed in the MSE, and ML fits of the spatial model use the
  ML information.
- 🗺️ **Unsampled domains**: areas, sub-areas and periods without a
  direct estimate (`y = NA`) stay in the model and are predicted, with
  their MSE.
- 🔀 **Unified interface**: one result class for all models, with
  `summary`, `coef`, `fitted`, `residuals`, `autoplot`, `compare_sae`,
  `benchmark_sae`, `map_sae` and `export_sae`.
- 🧵 **Multi-threaded bootstrap**: OpenMP parametric and non-parametric
  bootstrap MSE, with results that do not depend on the number of
  threads.
- 🧮 **Spatial and temporal effects**: SAR (EBLUP and HB), Besag ICAR,
  BYM, BYM2, Leroux, spatial lag, RW1, RW2, AR(1) and space-time
  interactions.
- 💾 **Sparse and structured linear algebra**: no dense inverse of the
  $D \times D$ (or $DT \times DT$) covariance matrix.

------------------------------------------------------------------------

### Supported Models

#### Frequentist EBLUP (C++ / OpenMP)

| Function | Model level | Structure | MSE estimation |
|:---|:---|:---|:---|
| `eblup_fh()` | Area | Fay-Herriot, iid area effects (optional log transformation) | Analytical (Prasad-Rao / Datta-Lahiri) |
| `eblup_sfh()` | Area | Spatial Fay-Herriot, SAR area effects | Analytical, parametric and non-parametric bootstrap |
| `eblup_stfh()` | Area | Spatio-temporal Fay-Herriot, SAR + AR(1) effects | Analytical, parametric bootstrap |
| `eblup_twofold()` | Sub-area | Two-fold sub-area model (area + sub-area effects), with area means | Analytical, parametric bootstrap |
| `eblup_bhf()` | Unit | Battese-Harter-Fuller nested error regression | Parametric bootstrap |

#### Hierarchical Bayes models (INLA, or MCMC with `jax_engine()`)

| Function | Model level | Likelihoods | Random effects |
|:---|:---|:---|:---|
| `hb_area()` | Area | Gaussian, beta, binomial, Poisson, negative binomial, gamma | iid, Besag ICAR, BYM, BYM2, Leroux, SAR, spatial lag; RW1, RW2, AR(1); space-time interactions |
| `hb_unit()` | Unit | Gaussian, binomial, Poisson | iid, Besag ICAR or BYM2 area effects |
| `hb_twofold()` | Sub-area | Gaussian, binomial, Poisson | iid, Besag ICAR or BYM2 area effects + iid sub-area effects |

The default priors of the HB models, and the simulations behind them,
are documented in the [priors
article](https://ridsonap.is-a.dev/fastsae/articles/hb-priors.html).

------------------------------------------------------------------------

### Performance at a Glance

<!-- BENCH-README:START -->

#### 1. Frequentist EBLUP with analytical MSE (time and peak memory)

| Model ($n = 1,000$ domains, 5 covariates) | `fastsae` | `sae` | `emdi` |
|:---|:--:|:--:|:--:|
| **Fay-Herriot** | **0.000639 s** (1.91 MB) | 0.799 s (39.0 MB; 1,250x slower) | 13.5 s (222 MB; 21,071x slower) |
| **Spatial Fay-Herriot** | **0.518 s** (78.7 MB) | 24.1 s (239 MB; 46x slower) | 25.2 s (310 MB; 49x slower) |
| **Spatio-Temporal FH ($T = 5$), point estimates** | **39.8 s** (495 MB) | 1160 s (3.60 GB; 29x slower) | \- |

#### 2. Bayesian Beta models: INLA vs Stan MCMC (largest $n$ benchmarked)

| Model | $n$ | `fastsae::hb_area` (INLA) | `tipsae::fit_sae` (Stan) | Ratio of times |
|:---|:--:|:--:|:--:|:--:|
| **Beta** | 1,000 | **3.73 s** | 180 s | 48x |
| **Spatial Beta (Besag ICAR)** | 1,000 | **4.28 s** | 362 s | 85x |
| **Spatio-Temporal Beta (Besag + RW1, $T = 5$)** | 50 | **3.80 s** | 716 s | 189x |

Times are medians of `bench::mark()`. The frequentist rows were measured
on 2026-10-26 with fastsae 1.0.0, sae 1.3 and emdi 2.2.3 on the same
simulated data sets (Apple M4, R version 4.6.1 (2026-06-24);
`benchmarks/bench_frequentist.R`), each in its own R process; memory is
the increase of the peak resident memory of that process during the
fits, which includes memory allocated by compiled code (increases of a
few tens of MB vary from run to run and should not be compared). ‘sae’
has no analytical MSE for the spatio-temporal model, so that row
compares point estimates; with its analytical MSE, fastsae took 39.3 s.
The Bayesian `fastsae` rows were re-measured on 2026-10-14
(`benchmarks/bench_bayes_fastsae.R`), the `tipsae` rows (4 chains of
2,000 iterations) come from an earlier run on the same data. See the
[benchmark
article](https://ridsonap.is-a.dev/fastsae/articles/benchmarks.html) for
all sizes.

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

To enable Bayesian models, install **INLA** (not on CRAN):

``` r
install.packages(
  "INLA",
  repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
  dep = TRUE
)
```

Optionally, the HB models can also be fitted by MCMC (NumPyro on JAX,
CPU or GPU) through `reticulate`; `setup_numpyro_env()` installs the
Python side and `setup_numpyro_env(check_only = TRUE)` checks it.

------------------------------------------------------------------------

## Quick Start

`mys` holds direct estimates of the mean years of schooling of 42 areas,
10 of them unsampled (`y = NA`). The unsampled areas stay in the data
and are predicted:

``` r
library(fastsae)

fit_fh <- eblup_fh(y ~ x1 + x2 + x3, vardir = ~vardir, domain = ~area, data = mys,
                   print_result = FALSE)
summary(fit_fh)
```

    ## 

    ## ── Summary of fastsae Fit ──────────────────────────────────────────────────────

    ## Call :

    ## eblup_fh(formula = y ~ x1 + x2 + x3, vardir = ~vardir, domain = ~area, data =
    ## mys, print_result = FALSE)

    ## 

    ## ✔ Convergence: Yes (in 6 iterations)

    ## Model: Fay-Herriot (Area-level)

    ## Method: REML

    ## 
    ## Variance Components:
    ## sigma2_u: 2.608103 
    ## 
    ## Coefficients:
    ##                   beta  std.error     zvalue pvalue
    ## (Intercept)  3.1077510  0.7697687  4.0372527 0.0001
    ## x1          -0.0019323  0.0098886 -0.1954019 0.8451
    ## x2           0.0555184  0.0614129  0.9040187 0.3660
    ## x3           0.0335344  0.0580013  0.5781663 0.5632
    ## 
    ## Goodness of Fit:
    ## loglikelihood           AIC           BIC 
    ##     -65.14251     140.28502     147.61370 
    ## 
    ## EBLUP Summary Statistics:
    ##       mse              rse        
    ##  Min.   :0.2037   Min.   : 9.844  
    ##  1st Qu.:0.6886   1st Qu.:14.889  
    ##  Median :1.0156   Median :22.157  
    ##  Mean   :1.5434   Mean   :26.275  
    ##  3rd Qu.:2.4886   3rd Qu.:35.020  
    ##  Max.   :3.7074   Max.   :53.065

``` r
head(fit_fh$df_eblup)
```

    ##   domain        y    eblup    vardir random_effect       mse       rse
    ## 1      1 8.359527 7.612738 0.6618838    2.94266461 0.5616180  9.844182
    ## 2      2 7.599650 6.782316 0.8374691    2.54539588 0.6784873 12.144869
    ## 3      3 5.514137 5.187060 0.8822257    0.96692722 0.7190378 16.347620
    ## 4      4 3.869326 4.201545 0.6581716   -1.31646412 0.5619314 17.841554
    ## 5      5 6.305063 6.323679 1.2788021   -0.03796677 0.9370273 15.307573
    ## 6      6 3.926807 4.048590 0.3878004   -0.81904048 0.3548310 14.713194

A spatial Fay-Herriot model with parametric bootstrap MSE, and the
hierarchical Bayes counterpart with a BYM2 area effect (needs INLA; the
CAR models use the neighbourhood graph `mys_proxmat > 0`):

``` r
fit_sfh <- eblup_sfh(y ~ x1 + x2 + x3, vardir = ~vardir, domain = ~area, data = mys,
                     W = mys_proxmat, mse_type = "pbmse", B = 200, seed = 1)

fit_hb <- hb_area(y ~ x1 + x2 + x3, vardir = "vardir", domain = "area", data = mys,
                  spatial = "bym2", W = mys_proxmat)

compare_sae(fit_fh, fit_hb)
autoplot(fit_fh, type = "estimates")
```

------------------------------------------------------------------------

## 📖 Documentation & Tutorials

Explore comprehensive guides, mathematical formulations, and worked
examples:

👉
[**https://ridsonap.is-a.dev/fastsae/**](https://ridsonap.is-a.dev/fastsae/)

🇮🇩 Dokumentasi dalam Bahasa Indonesia (Indonesian translation):
[**Beranda**](https://ridsonap.is-a.dev/fastsae/articles/id/beranda.html)

- [Getting Started with
  fastsae](https://ridsonap.is-a.dev/fastsae/articles/getting-started.html)
- [Spatial & Spatio-Temporal Models
  (EBLUP)](https://ridsonap.is-a.dev/fastsae/articles/spatial-temporal.html)
- [Bayesian Area-Level Models
  (INLA)](https://ridsonap.is-a.dev/fastsae/articles/hb-area-inla.html)
- [Default Priors of the HB
  Models](https://ridsonap.is-a.dev/fastsae/articles/hb-priors.html)
- [Model Diagnostics &
  Calibration](https://ridsonap.is-a.dev/fastsae/articles/model-diagnostics.html)
- [Unit-Level Estimation (BHF and
  HB)](https://ridsonap.is-a.dev/fastsae/articles/unit-level-bhf.html)
- [Sub-Area Level Two-Fold
  Models](https://ridsonap.is-a.dev/fastsae/articles/twofold-subarea.html)
- [Performance
  Benchmarks](https://ridsonap.is-a.dev/fastsae/articles/benchmarks.html)
- [hb_area() Engines: INLA vs JAX vs
  tipsae](https://ridsonap.is-a.dev/fastsae/articles/hb-area-engines.html)
- [Finite Population
  Validation](https://ridsonap.is-a.dev/fastsae/articles/finite-population-validation.html)

------------------------------------------------------------------------

## Citation

``` r
citation("fastsae")
```

Bug reports and feature requests are welcome on the [issue
tracker](https://github.com/ridsonap/fastsae/issues).
