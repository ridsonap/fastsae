# fastsae: High-Performance Frequentist and Bayesian Small Area Estimation in R

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/fastsae)](https://CRAN.R-project.org/package=fastsae)
[![R-CMD-check](https://github.com/ridsonap/fastsae/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ridsonap/fastsae/actions/workflows/R-CMD-check.yaml)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**fastsae** implements widely used Small Area Estimation (SAE) models for survey data. Frequentist estimation runs in compiled **C++** (`RcppArmadillo`) with **OpenMP** multi-threading, and Bayesian estimation uses **INLA** (Integrated Nested Laplace Approximations). With hundreds to thousands of domains it is much faster than existing R implementations (see [Performance at a Glance](#performance-at-a-glance); with a few dozen domains `sae` is about as fast), while agreeing numerically with reference packages such as `sae`.

---

## Why fastsae?

- ⚡ **Speed**: with 1,000 domains, the C++ Fisher-scoring fits with analytical MSE are 29 to 1,250 times faster than `sae` (spatio-temporal model: 40 s instead of 19 minutes); with a few dozen domains the differences are a fraction of a second.
- 🎯 **Numerical agreement**: on the same data, coefficients, variance components, EBLUPs and analytical MSEs match `sae` (Molina & Marhuenda) to about $10^{-14}$; the documented exceptions are a spatial correlation on the boundary and the ML MSE of the spatial model.
- 🔀 **Frequentist and Bayesian in one package**: EBLUP models in C++ and hierarchical Bayesian models through INLA.
- 🧵 **Parallel bootstrap**: multi-core OpenMP parametric and non-parametric bootstrap MSE estimation.
- 🗺️ **Unsampled areas**: domains with missing observations (`y = NA`) are predicted automatically.
- 💾 **Low memory**: avoids unnecessary $O(m^2)$ / $O(m^3)$ allocations.
- 📊 **Standard R interface**: S3 methods (`summary`, `coef`, `fitted`, `residuals`) and publication-ready `autoplot()`.

---

## Understanding INLA: Fast Bayesian Inference without MCMC

**INLA** (Integrated Nested Laplace Approximations; Rue et al., 2009) provides a deterministic alternative to Markov Chain Monte Carlo (MCMC) for Bayesian inference in latent Gaussian models. Instead of drawing thousands of samples, INLA:

- **Approximates** posterior marginals analytically using Laplace approximations
- **Integrates** out non-Gaussian observations analytically
- **Delivers** results in seconds rather than hours

`fastsae` uses INLA through the [`INLA`](https://www.r-inla-download.org) R package to fit hierarchical Bayesian small area models with:

- **6 likelihood families**: Gaussian, Beta, Binomial, Poisson, Negative Binomial, Gamma
- **Spatial random effects**: Besag ICAR, BYM, BYM2, Leroux, SAR, spatial lag
- **Temporal random effects**: RW1, RW2, AR(1)
- **Spatio-temporal interactions**: separable and non-separable structures

**Compared with MCMC**: no chains to run and monitor, deterministic results, and fits in seconds. INLA is an approximation: fastsae chose its settings by comparing them with exact posteriors and long Stan runs, and `engine = jax_engine()` fits the same models by MCMC when a check is needed.

---

## Two-Fold Subarea Models: Nested Hierarchical Structure

The two-fold subarea model (Torabi & Rao, 2014) extends the standard Fay-Herriot model to hierarchical data where **subareas are nested within areas**:

$$y_{ij} = x_{ij}^\top \beta + v_i + u_{ij} + e_{ij}$$

where:
- $v_i \sim N(0, \sigma_v^2)$: **Area-level random effect**
- $u_{ij} \sim N(0, \sigma_u^2)$: **Subarea-level random effect**
- $e_{ij} \sim N(0, D_{ij})$: Sampling error with known variance

**Available implementations:**
- `eblup_twofold()`: Frequentist EBLUP via REML Fisher-scoring (C++)
- `hb_twofold()`: Bayesian hierarchical model via INLA with spatial support

---

## Comparison with Other Packages

| Feature | `fastsae` | `sae` (Molina & Rao) | `emdi` (Kreutzmann et al.) |
|:---|:---:|:---:|:---:|
| **Core Computation** | **C++ (RcppArmadillo)** | Pure R | R / lme4 / nlme |
| **Multi-Threading** | **Native OpenMP** (`n_threads`) | Single-threaded | Optional foreach/parallel |
| **Agreement with `sae` on the same data** | **About $10^{-14}$** | Reference | About $10^{-6}$ (convergence tolerance) |
| **Unsampled Area Support** | **Automatic (synthetic, kriging for spatial models)** | Manual subsetting required | Synthetic (out-of-sample domains) |
| **Bootstrap MSE** | **Parallel C++ (OpenMP)** | R | R |
| **S3 Methods Support** | `print`, `summary`, `coef`, `fitted`, `residuals`, `autoplot` | Custom lists | Standard S3 |

---

## Supported Models

### Frequentist Models (EBLUP, C++)

| Model | Function | Random Effect Structure | MSE Estimation |
|:---|:---|:---|:---|
| **Fay-Herriot** (Area-level) | `eblup_fh()` | Independent area effects | Analytical (Prasad-Rao) |
| **Spatial Fay-Herriot** | `eblup_sfh()` | SAR(1) spatial autocorrelation | Analytical, Bootstrap |
| **Spatio-Temporal Fay-Herriot** | `eblup_stfh()` | SAR(1) + Temporal AR(1) | Analytical, Bootstrap |
| **Battese-Harter-Fuller** (Unit-level) | `eblup_bhf()` | Random intercept nested in domains | Bootstrap |
| **Two-Fold Subarea** | `eblup_twofold()` | Area + Subarea nested effects | Analytical, Bootstrap |

### Bayesian Models (INLA)

| Model | Function | Distributions | Spatial/Temporal |
|:---|:---|:---|:---|
| **Area-level HB** | `hb_area()` | Gaussian, Beta, Binomial, Poisson, NB, Gamma | ✓ Full support |
| **Unit-level HB** | `hb_unit()` | Gaussian, Binomial, Poisson | ✓ BYM2/Besag |
| **Two-Fold HB** | `hb_twofold()` | Gaussian, Binomial, Poisson | ✓ BYM2/Besag |

The Bayesian models are fitted with `engine = inla_engine()` (default). `hb_area()`, `hb_unit()` and `hb_twofold()` can also be fitted by MCMC in Python (JAX/NumPyro) with `engine = jax_engine()`.

---

## Performance at a Glance

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

Times are medians of `bench::mark()`. The frequentist rows were measured on 2026-10-26 with fastsae 1.0.0, sae 1.3 and emdi 2.2.3 on the same simulated data sets (Apple M4, R version 4.6.1 (2026-06-24); `benchmarks/bench_frequentist.R`), each in its own R process; memory is the increase of the peak resident memory of that process during the fits, which includes memory allocated by compiled code (increases of a few tens of MB vary from run to run and should not be compared). 'sae' has no analytical MSE for the spatio-temporal model, so that row compares point estimates; with its analytical MSE, fastsae took 39.3 s. The Bayesian `fastsae` rows were re-measured on 2026-10-14 (`benchmarks/bench_bayes_fastsae.R`), the `tipsae` rows (4 chains of 2,000 iterations) come from an earlier run on the same data. See the [benchmark article](https://ridsonap.is-a.dev/fastsae/articles/benchmarks.html) for all sizes.

<!-- BENCH-README:END -->

---

## Installation

```r
# Install from CRAN (once published)
install.packages("fastsae")

# Or install the development version from GitHub
# install.packages("remotes")
remotes::install_github("ridsonap/fastsae")
```

---

## Documentation and User Guides

Explore detailed guides, worked examples, and performance analyses:

- [**Getting Started with fastsae**](articles/getting-started.html): Step-by-step walkthrough of model fitting, diagnostics, and S3 extraction.
- [**An Optimal SAE Workflow**](articles/optimal-workflow.html): The order of a complete analysis, from reliability of the direct estimates to benchmarking and reporting.
- [**Spatial & Spatio-Temporal Models (EBLUP)**](articles/spatial-temporal.html): Working with spatial weights matrices, spatial kriging, and spatio-temporal panels.
- [**Bayesian Area-Level Models (INLA)**](articles/hb-area-inla.html): Hierarchical Bayesian SAE using INLA with 6 distributions and spatial-temporal effects.
- [**Model Diagnostics & Calibration**](articles/model-diagnostics.html): Diagnostics and benchmarking to aggregate targets.
- [**Unit-Level Estimation (BHF)**](articles/unit-level-bhf.html): Fitting the Battese-Harter-Fuller model with individual survey data.
- [**Performance Benchmarks**](articles/benchmarks.html): Interactive charts and speedup comparisons across sample sizes.
- [**Function Reference**](reference/index.html): Comprehensive API reference and parameter documentation.

---

## References

- Fay, R. E., & Herriot, R. A. (1979). Estimates of income for small places: An application of James-Stein procedures to Census data. *Journal of the American Statistical Association*, 74(366), 269–277.
- Battese, G. E., Harter, R. M., & Fuller, W. A. (1988). An error-components model for prediction of county crop areas using survey and satellite data. *Journal of the American Statistical Association*, 83(401), 28–36.
- Marhuenda, Y., Molina, I., & Morales, D. (2013). Small area estimation with spatio-temporal Fay-Herriot models. *Computational Statistics & Data Analysis*, 58, 308–325.
- Pratesi, M., & Salvati, N. (2008). Small area estimation: The EBLUP estimator based on spatially correlated random area effects. *Statistical Methods & Applications*, 17(1), 113–141.
- Torabi, M., & Rao, J. N. K. (2014). On small area estimation under a sub-area level model. *Journal of Multivariate Analysis*, 127, 36–55.
- Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian inference for latent Gaussian models using integrated nested Laplace approximations (with discussion). *Journal of the Royal Statistical Society: Series B*, 71(2), 319–392.
- Rao, J. N. K., & Molina, I. (2015). *Small Area Estimation* (2nd ed.). John Wiley & Sons.
