# fastsae 0.1.1

### Major New Features & Enhancements
* **Hierarchical Bayes Area-Level Models (`hb_area`)**:
  * Fast Bayesian inference powered by Integrated Nested Laplace Approximations (`INLA`).
  * Supports 6 probability distributions for non-Gaussian area indicators: Gaussian, Binomial (logit link with trials), Poisson (log rate with exposure/offsets), Negative Binomial, Beta, and Gamma.
  * Spatial random effects: BYM2 (Besag-York-Mollié 2 with PC-priors) and Besag/CAR models.
  * Automatic handling of unsampled domains and domain-specific variance predictions.
* **Hierarchical Bayes Unit-Level Models (`hb_unit`)**:
  * Bayesian unit-level Battese-Harter-Fuller (BHF) and generalized linear mixed models via `INLA`.
  * Supports Gaussian, Binomial (binary unit outcomes), and Poisson (count outcomes).
  * Domain random effects with optional spatial structures (BYM2, Besag).
  * Finite population adjustments when domain population sizes ($N_d$) are supplied.
  * Seamless handling of both domain-level population means and micro-level census data.
* **Two-Stage Hierarchical Benchmarking (`benchmark_sae`)**:
  * Unified calibration framework for multi-level administrative structures (e.g., Domain -> Province -> National).
  * Implements Ratio, Difference, Optimal Variance-Weighted Quadratic Loss (Datta et al. 2011; Steorts et al. 2014), and Constrained Logit benchmarking (Berg & Fuller 2014).
* **Comprehensive Diagnostics & Visualization**:
  * `diagnose()`: External calibration tests (Brown et al. 2001), residual analysis, spatial autocorrelation tests (Moran's I), and official statistical reliability classification (BPS standards).
  * `compare_sae()`: Tabular model comparison across information criteria and error metrics.
  * `map_sae()`: Spatial thematic choropleth maps integrated with `sf`.
  * `export_sae()`: Export estimation results with full statistical metadata to Excel (`openxlsx`/`writexl`) and CSV.
* **Simulation Utilities**:
  * `sim_area_data()`, `sim_spatial_weights()`, and `sim_series_data()`.
* **Reliability & Performance**:
  * High-performance C++ core (`RcppArmadillo` and OpenMP) with strict validation checks.

# fastsae 0.1.0

* Initial development version.
* High-performance Area-level Fay-Herriot (FH), Spatial FH (SFH), and Spatio-Temporal FH (STFH) models.
* Unit-level Battese-Harter-Fuller (BHF) models via Henderson III and REML.
