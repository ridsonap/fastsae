# Performance & Scalability Benchmarks

## Executive Summary

Small Area Estimation often involves large administrative registries or
extensive spatial networks with hundreds or thousands of domains.
Traditional implementations in R often rely on interpreted loops, pure R
Fisher-scoring routines, or large memory allocations for $`D \times D`$
dense covariance matrices.

**fastsae** fits the frequentist models in compiled C++
(`RcppArmadillo`; sparse Cholesky factorisations for the SAR models) and
the Bayesian models with INLA. At $`n = 1,000`$ domains:

- **Fay-Herriot**: 1,250x faster than `sae` and 21,071x faster than
  `emdi`.
- **Spatial Fay-Herriot**: 46x faster than `sae` and 49x faster than
  `emdi`.
- **Spatio-Temporal Fay-Herriot** ($`T = 5`$, point estimates): 29x
  faster than `sae`.

------------------------------------------------------------------------

## Benchmark Results Table

| Model ($`n = 1,000`$ domains, 5 covariates) | `fastsae` | `sae` | `emdi` |
|:---|:--:|:--:|:--:|
| **Fay-Herriot** | **0.000639 s** (1.91 MB) | 0.799 s (39.0 MB; 1,250x slower) | 13.5 s (222 MB; 21,071x slower) |
| **Spatial Fay-Herriot** | **0.518 s** (78.7 MB) | 24.1 s (239 MB; 46x slower) | 25.2 s (310 MB; 49x slower) |
| **Spatio-Temporal FH ($`T = 5`$), point estimates** | **39.8 s** (495 MB) | 1160 s (3.60 GB; 29x slower) | \- |

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
2,000 iterations) come from an earlier run on the same data.

------------------------------------------------------------------------

## Interactive Benchmark Explorer

Use the interactive controls below to compare execution time, memory
allocated, and iterations per second across sample sizes:

Fay-Herriot (FH)

Spatial FH (SFH)

Spatio-Temporal (STFH)

Execution Time

Memory Allocated

Throughput (iter/s)

Execution time (median, seconds) — FH · log scale

Speedup Factor (fastsae vs competitors) — FH

------------------------------------------------------------------------

## Architectural Insights: Why is fastsae so Fast?

1.  **Compiled C++ Linear Solvers**: Instead of interpreting nested
    loops in R, `fastsae` implements Fisher scoring and the analytical
    MSE in Armadillo C++. The spatial Fay-Herriot model never forms a
    $`D \times D`$ inverse: with $`A = (I - \rho W')(I - \rho W)`$ and
    $`M = A + \sigma^2 E'\Psi^{-1}E`$, $`V^{-1}`$, $`g_1, \dots, g_4`$
    and the kriging terms follow from sparse Cholesky solves with $`A`$
    and $`M`$ (Eigen, AMD ordering).

2.  **Zero-Copy Matrix Operations**: Memory allocations for large
    intermediate structures ($`V`$, $`V^{-1}`$, $`P`$) are avoided or
    reused across Fisher iterations rather than reallocated on the heap.

3.  **OpenMP Multi-Threaded Bootstrap**: Bootstrap resampling runs
    natively in parallel across CPU cores using
    `#pragma omp parallel for`, avoiding the serialization overhead of R
    worker processes
    ([`parallel::makeCluster`](https://rdrr.io/r/parallel/makeCluster.html)
    / `foreach`).

4.  **Woodbury Identity for Panel Data**: In `eblup_stfh`, inversion of
    the block $`(DT \times DT)`$ covariance matrix $`V`$ is reduced to
    operations on individual $`D \times D`$ and $`T \times T`$ blocks
    via Kronecker and Woodbury decomposition, transforming an
    $`O((DT)^3)`$ bottleneck into scalable $`O(D^3 + T^3)`$ steps.

5.  **Integrated Nested Laplace Approximations (INLA)**: For complex
    hierarchical generalized and spatio-temporal models (`hb_area`),
    `fastsae` utilizes INLA to compute analytical posterior margins
    directly from sparse Gaussian Markov Random Fields, bypassing Monte
    Carlo sampling entirely.

------------------------------------------------------------------------

## Bayesian Spatio-Temporal Benchmark: `fastsae` (INLA) vs `tipsae` (Stan MCMC)

To assess Bayesian small area estimation performance,
[`fastsae::hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
was benchmarked against
[`tipsae::fit_sae()`](https://rdrr.io/pkg/tipsae/man/fit_sae.html), the
state-of-the-art Stan MCMC package for spatio-temporal Beta small area
models.

The test used the official `tipsae` panel dataset: **`emilia`** (38
health districts in Emilia-Romagna over 5 years, $`N = 190`$ domains
$`\times`$ years) with spatial polygon contiguity matrix $`W`$ from
**`emilia_shp`**.

| Metric | [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html) (Stan MCMC) | [`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md) (INLA) |
|:---|:---|:---|
| **Model** | Beta, Besag ICAR + domain-specific RW(1) | `family = "beta"`, `spatial = "besag"`, `temporal = "rw1"`, `st_interaction = "domain-specific"`; same structure, different priors |
| **Computation** | Hamiltonian Monte Carlo (NUTS) | Integrated nested Laplace approximation (INLA’s default `"compact"` mode for Beta) |
| **Execution time** | 22 min (4 chains x 10,000 iterations with 2,000 warm-up, in parallel; largest $`\hat R`$ 1.007). Every transition reaches the maximum tree depth; 4 chains x 2,000 iterations (120 s) have not converged ($`\hat R`$ 1.07, and 1.12 on a Linux runner) | 1.1 s |
| **Posterior means vs converged `tipsae`** | Reference | $`r = 0.999`$, mean absolute difference 0.0017 (largest 0.006) |
| **Posterior SD vs `tipsae`** | Reference | ratio 0.81 to 1.22 (the priors differ) |

Measured on 2026-10-16 (`tipsae` reference refitted on 2026-10-17) on
one machine (Apple silicon, R 4.6.1, INLA 26.8.7, tipsae 1.0.3); times
depend on the hardware. INLA needs no chains or warm-up, but its fit
should still be checked (`convergence`, from INLA’s mode search).

### Interactive Beta SAE Benchmark Explorer ($`n = 30`$ to $`n = 1,000`$)

Use the interactive controls below to compare runtime, memory footprint,
and speedup factors between
[`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
and [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html)
across **Standard Beta**, **Spatial Beta (Besag ICAR)**, and
**Spatio-Temporal Beta (Besag + RW1)** models:

Standard Beta

Spatial Beta (ICAR)

Spatio-Temporal Beta

Execution Time

Memory Allocated

Throughput (iter/s)

Execution time (median, seconds) — Standard Beta · log scale

Speedup Factor & Efficiency Multiplier (fastsae vs tipsae) — Standard
Beta

#### Agreement with Stan MCMC (Beta model, $`n = 1,000`$)

To check that the speed of
[`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
does not come at the cost of its estimates, both packages were fitted
with their defaults to the $`n = 1,000`$ data set of the Beta benchmark
(independent area effects, sampling variances from the direct estimates;
`benchmarks/bench_beta_concordance.R`). The two packages do not fit
identical models (their priors differ, and fastsae fixes the precision
at a first-stage estimate while tipsae makes it a function of the model
mean), so the comparison is between what a user obtains from each
package:

| Metric | [`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md) (INLA) | [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html) (Stan) |
|:---|:--:|:--:|
| RMSE against the true area means | 0.0678 | 0.0684 |
| Coverage of the 95% intervals | 93.3% | 91.1% |
| Mean width of the 95% intervals | 0.252 | 0.259 |
| Time | 4.3 s | 184 s |

The posterior means of the two packages differed by 0.0096 on average
(absolute difference), and the posterior SDs of fastsae were 0.99 times
those of tipsae in the median area (90% of the areas between 0.88 and
1.17). The largest $`\hat{R}`$ of the tipsae fit was 1.008.

------------------------------------------------------------------------

### Finite Population Simulation & Ground Truth Validation

For comprehensive design-based simulation studies evaluating parameter
recovery and empirical properties under repeated survey sampling from
finite populations ($`D = 60`$ and national-scale $`D = 500`$), see the
dedicated article: [Finite Population Simulation &
Validation](https://ridsonap.is-a.dev/fastsae/articles/finite-population-validation.md).
