# fastsae: High-Performance Frequentist and Bayesian Small Area Estimation in R

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/fastsae)](https://CRAN.R-project.org/package=fastsae)
[![R-CMD-check](https://github.com/ridsonap/fastsae/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ridsonap/fastsae/actions/workflows/R-CMD-check.yaml)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**fastsae** implements widely used Small Area Estimation (SAE) models for survey data. Frequentist estimation runs in compiled **C++** (`RcppArmadillo`) with **OpenMP** multi-threading, and Bayesian estimation uses **INLA** (Integrated Nested Laplace Approximations). It provides massive speedups (up to **6,400x faster**) while maintaining **exact numerical equivalence** with gold-standard implementations.

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
| **RAM Consumption** | **Minimal (< 25 MB)** | Moderate (~100 MB) | High (~800+ MB to 4 GB) |
| **S3 Methods Support** | `print`, `summary`, `coef`, `fitted`, `residuals`, `autoplot` | Custom lists | Standard S3 |

---

## Supported Models

### Frequentist Models (EBLUP, C++)

| Model | Function | Random Effect Structure | MSE Estimation |
|:---|:---|:---|:---|
| **Fay-Herriot** (Area-level) | `eblup_fh()` | Independent area effects | Analytical (Prasad-Rao) |
| **Spatial Fay-Herriot** | `eblup_sfh()` | SAR(1) spatial autocorrelation | Analytical, Bootstrap |
| **Spatio-Temporal Fay-Herriot** | `eblup_stfh()` | SAR(1) + Temporal AR(1) | Bootstrap |
| **Battese-Harter-Fuller** (Unit-level) | `eblup_bhf()` | Random intercept nested in domains | Bootstrap |
| **Two-Fold Subarea** | `eblup_twofold()` | Area + Subarea nested effects | Analytical, Bootstrap |

### Bayesian Models (INLA)

| Model | Function | Distributions | Spatial/Temporal |
|:---|:---|:---|:---|
| **Area-level HB** | `hb_area()` | Gaussian, Beta, Binomial, Poisson, NB, Gamma | ✓ Full support |
| **Unit-level HB** | `hb_unit()` | Gaussian | ✓ Via nested error |
| **Two-Fold HB** | `hb_twofold()` | Gaussian, Binomial, Poisson | ✓ BYM2/Besag |

---

## Performance at a Glance

Benchmark across $n = 1,000$ areas (5 covariates):

| Task | `fastsae` | `sae` | `emdi` |
|:---|:--:|:--:|:--:|
| **Standard Fay-Herriot** | **0.0015 s** | 0.291 s (~190x slower) | 9.64 s (~6,400x slower) |
| **Spatial Fay-Herriot** | **0.165 s** | 12.60 s (~76x slower) | 8.69 s (~53x slower) |
| **Peak RAM usage** | **< 10 MB** | ~400 MB | ~800 MB |

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
