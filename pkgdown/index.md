# fastsae: High-Performance Frequentist and Bayesian Small Area Estimation in R

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/fastsae)](https://CRAN.R-project.org/package=fastsae)
[![R-CMD-check](https://github.com/ridsonap/fastsae/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ridsonap/fastsae/actions/workflows/R-CMD-check.yaml)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**fastsae** implements widely used Small Area Estimation (SAE) models for survey data. Frequentist estimation runs in compiled **C++** (`RcppArmadillo`) with **OpenMP** multi-threading, and Bayesian estimation uses **INLA** (Integrated Nested Laplace Approximations). It is orders of magnitude faster than existing R implementations (see [Performance at a Glance](#performance-at-a-glance)) while agreeing numerically with reference packages such as `sae`.

---

## Why fastsae?

- ⚡ **Speed**: optimized C++ Fisher-scoring algorithms turn fits that take seconds in other R packages into milliseconds.
- 🎯 **Numerical agreement**: point estimates ($\hat{\beta}$, $\hat{\theta}$) and variance components ($\hat{\sigma}_u^2$, $\hat{\rho}$) match `sae` (Molina & Marhuenda) to numerical tolerance.
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
- **Spatial random effects**: BYM2, Besag ICAR, Spatial Lag
- **Temporal random effects**: RW1, RW2, AR(1)
- **Spatio-temporal interactions**: separable and non-separable structures

**Key advantage over MCMC**: No convergence diagnostics, no burn-in, reproducible results, and 10–100× speedups.

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
| **Numerical Consistency** | **Reference baseline** | Baseline | Approximations |
| **Unsampled Area Support** | **Automatic (Spatial Kriging)** | Manual subsetting required | Limited |
| **Bootstrap Speed** | **Parallel C++ (Ultra-Fast)** | Slow (Interpreted R loops) | Moderate |
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
| **Unit-level HB** | `hb_unit()` | Gaussian | ✓ Via nested error |
| **Two-Fold HB** | `hb_twofold()` | Gaussian, Binomial, Poisson | ✓ BYM2/Besag |

The Bayesian models are fitted with `engine = inla_engine()` (default). `hb_area()` can also be fitted by MCMC in Python (JAX/NumPyro) with `engine = jax_engine()`.

---

## Performance at a Glance

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

Times are medians of `bench::mark()`; memory is the total memory allocated by R during one call (`bench::mark()$mem_alloc`), not the peak resident memory. The `fastsae` rows were re-measured with fastsae 1.0.0 (frequentist on 2026-10-13, `benchmarks/bench_frequentist_fastsae.R`; Bayesian on 2026-10-14, `benchmarks/bench_bayes_fastsae.R`); the `sae`/`emdi` rows come from an earlier run of the same design (simulated data, 5 covariates) on different random draws, the `tipsae` rows from an earlier run on the same data. See the [benchmark article](https://ridsonap.is-a.dev/fastsae/articles/benchmarks.html) for all sizes.

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
- Pratesi, M., & Salvati, N. (2008). Small area estimation for spatially correlated data: A Fay-Herriot with the spatial linear spline model. *Journal of Applied Statistics*, 35(7), 781–794.
- Torabi, M., & Rao, J. N. K. (2014). On small area estimation under a semi-parametric mixed model. *Journal of Multivariate Analysis*, 124, 226–236.
- Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian inference for latent Gaussian models using integrated nested Laplace approximations (with discussion). *Journal of the Royal Statistical Society: Series B*, 71(2), 319–392.
- Rao, J. N. K., & Molina, I. (2015). *Small Area Estimation* (2nd ed.). John Wiley & Sons.
