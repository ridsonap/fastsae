# Changelog

## fastsae 1.0.0

#### Estimation engines for `hb_area()`, `hb_unit()` and `hb_twofold()`

- New argument `engine`, taking an engine object from the new
  constructors
  [`inla_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)
  and
  [`jax_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)
  (or the strings `"inla"`, the default, and `"jax"` for their default
  settings). Each constructor holds all settings of its engine, so
  [`hb_area()`](https://ridsonap.github.io/fastsae/reference/hb_area.md)
  itself only describes the model:
  - `inla_engine(strategy, link, scale_model, prior_prec, prior_phi, prior_rho, prior_prec_time, prior_rho_time, prior_prec_subarea, draws, seed, ...)`,
    where `...` goes to
    [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html).
    `control.*` lists given there are merged with the package defaults
    (before, they clashed with them in
    [`hb_area()`](https://ridsonap.github.io/fastsae/reference/hb_area.md)
    and replaced them in
    [`hb_unit()`](https://ridsonap.github.io/fastsae/reference/hb_unit.md)/[`hb_twofold()`](https://ridsonap.github.io/fastsae/reference/hb_twofold.md)).
    Priors and `draws` left at `NULL` take the default of the model
    function, so the defaults of each function are unchanged; settings a
    function does not use are ignored with a warning.
  - `jax_engine(warmup, samples, chains, device, seed, prior_beta, robust, keep_draws, progress)`
    fits the same area-level models by full MCMC (NumPyro NUTS on JAX:
    CPU, NVIDIA GPU or the Apple silicon GPU via `device = "metal"`)
    through ‘reticulate’ (Suggests); the engine formerly shipped in the
    separate ‘fastsaegpu’ package now lives in fastsae. It supports all
    six families, `bym2`/`bym`/`besag`/`generic1` (alias `"leroux"`),
    `rw1`/`ar1`/`iid` and every `st_interaction`, plus
    `prior_beta = "horseshoe"` and `robust = TRUE` (Student-t area
    effects). Results have the same structure as the INLA engine, with
    the posterior draws in `draws`, the largest split R-hat in
    `rhat_max` (`convergence` is `rhat_max < 1.05`), and DIC/WAIC in
    `goodness`. Posterior self-benchmarking (also by `benchmark_group`)
    adjusts the MCMC draws exactly as for INLA, and
    [`benchmark_sae()`](https://ridsonap.github.io/fastsae/reference/benchmark_sae.md)
    uses the draws for the MSE matrix.
  - The engine settings are returned in `engine_spec`; `engine` holds
    the engine name.
- [`hb_unit()`](https://ridsonap.github.io/fastsae/reference/hb_unit.md)
  and
  [`hb_twofold()`](https://ridsonap.github.io/fastsae/reference/hb_twofold.md)
  take the same `engine = inla_engine(...)` (INLA only). In
  [`hb_twofold()`](https://ridsonap.github.io/fastsae/reference/hb_twofold.md),
  `prior_prec_area` is now `inla_engine(prior_prec = )`, `n_samples` is
  `inla_engine(draws = )`, and the posterior draws for the area
  aggregates honour `seed`. Their results carry `engine` and
  `engine_spec` too.
- `jax_engine(device = "metal")` now samples on the Apple GPU. jax-metal
  cannot compile two operations of NumPyro’s NUTS sampler (the `popcnt`
  bit count, and `lax.cond` statements whose outputs are constants), so
  fits used to fall back to the CPU. fastsae now replaces both with
  equivalent operations that Metal can compile, and steps the chains
  from Python; the draws match the CPU sampler. On an Apple M4, Metal is
  faster than the CPU for large models (about 2x at 10,000 domains) and
  slower for small ones, where a fixed cost of about two minutes per fit
  dominates (see the “hb_area() engines” benchmark article).
  `device = "auto"` keeps using the CPU on macOS.
- New
  [`setup_numpyro_env()`](https://ridsonap.github.io/fastsae/reference/setup_numpyro_env.md)
  and
  [`check_numpyro_available()`](https://ridsonap.github.io/fastsae/reference/check_numpyro_available.md)
  to install and check the Python side. Nothing is installed or
  downloaded unless
  [`setup_numpyro_env()`](https://ridsonap.github.io/fastsae/reference/setup_numpyro_env.md)
  is called.
- The frequentist lme4 GLMM route (`method = "laplace"` for non-spatial
  binomial and Poisson models) has been removed; such models are now
  fitted by INLA. ‘lme4’ moved to Suggests.
- Deprecated: `method` (`method = "laplace"` becomes
  `engine = inla_engine(strategy = "laplace")`), and passing `strategy`,
  `link`, `scale_model`, `prior_*`, `benchmark_draws`, `benchmark_seed`,
  `n_samples` or
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html) arguments
  directly to
  [`hb_area()`](https://ridsonap.github.io/fastsae/reference/hb_area.md),
  [`hb_unit()`](https://ridsonap.github.io/fastsae/reference/hb_unit.md)
  or
  [`hb_twofold()`](https://ridsonap.github.io/fastsae/reference/hb_twofold.md).
  They still work, with a warning, and are forwarded to
  [`inla_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md).
- [`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)
  show the engine, and the self-benchmarking line now names the method
  used (posterior draws or Wang-Fuller-Qu) instead of always reporting
  Wang-Fuller-Qu.

#### Statistical corrections (mathematical audit)

- [`eblup_stfh()`](https://ridsonap.github.io/fastsae/reference/eblup_stfh.md):
  the analytical $`g_3`$ was computed with the wrong matrix orientation
  ($`(K_a \Psi V \Psi K_b)_{ii}`$ instead of
  $`\psi_i^2 (K_a V K_b)_{ii}`$), which is wrong whenever `vardir` is
  heterogeneous. $`g_3`$ now uses the exact derivative of the EBLUP,
  $`\partial\tilde\theta/\partial\theta_a = \Psi P V_a P y`$, and the
  REML information at the final estimates.

- [`eblup_stfh()`](https://ridsonap.github.io/fastsae/reference/eblup_stfh.md)
  bootstrap MSE: the variance components are now re-estimated in every
  replicate (before, they were held fixed, so the bootstrap missed the
  $`g_3`$ contribution), and the SAR effects are generated with the
  lower Cholesky factor (before,
  $`\mathrm{Var}(u_1) = RR' \ne \Omega_1`$).

- [`eblup_stfh()`](https://ridsonap.github.io/fastsae/reference/eblup_stfh.md)
  now sorts the panel internally into the domain-major order the model
  requires and returns results in the caller’s row order. Before,
  unsorted input was silently mis-specified – including the time-major
  `mys_panel` used in the examples, tests and vignette (where
  [`sae::eblupSTFH()`](https://rdrr.io/pkg/sae/man/eblupSTFH.html),
  which assumes domain-major data, gave the same wrong fit).

- [`eblup_stfh()`](https://ridsonap.github.io/fastsae/reference/eblup_stfh.md):
  when the classical Fisher scoring fails, or converges to a negative
  variance that would only be truncated to 0, the fit switches to
  constrained REML (variances \>= 0, \|rho\| \<= 0.999; projected
  scoring with step halving), which returns the actual constrained
  maximiser. Convergence in a simulation went from 63% to 99% of the
  replicates; fits on which the classical iteration converges to an
  interior point are unchanged (identical to
  [`sae::eblupSTFH()`](https://rdrr.io/pkg/sae/man/eblupSTFH.html)). A
  correlation whose variance is 0 is unidentified and gets an `NA`
  standard error.

- [`eblup_stfh()`](https://ridsonap.github.io/fastsae/reference/eblup_stfh.md)
  bootstrap: non-converged replicates are replaced by new draws until
  `B` valid replicates are obtained (as in the SFH bootstrap).

- [`eblup_twofold()`](https://ridsonap.github.io/fastsae/reference/eblup_twofold.md):
  the analytical MSE estimator is now $`g_1 + g_2 + 2 g_3`$ (REML) and
  $`g_1 + g_2 + 2 g_3 - b_{ML}'\nabla g_1`$ (ML); before it was
  $`g_1 + g_2 + g_3`$, which ignores the downward bias of the plug-in
  $`g_1`$. The REML Fisher information now includes the cross-area
  terms.

- [`eblup_sfh()`](https://ridsonap.github.io/fastsae/reference/eblup_sfh.md)
  with unsampled domains: the SAR process is defined on all domains of
  `W`, so the sampled domains are now fitted with the marginal
  covariance $`\sigma^2 \Sigma_{full}[s,s]`$ (before: the SAR built from
  the sub-matrix `W[s,s]`), consistent with the kriging of the unsampled
  domains and the bootstrap.

- [`eblup_sfh()`](https://ridsonap.github.io/fastsae/reference/eblup_sfh.md):
  an upper-boundary $`\hat\rho = 0.999`$ is no longer mapped to $`1`$
  (singular for a row-standardised `W`); bootstrap replicates are no
  longer rejected when estimates exceed `1e6` (scale-dependent check);
  `npbmse` no longer fails when $`\hat\sigma^2_u = 0`$.

- [`hb_unit()`](https://ridsonap.github.io/fastsae/reference/hb_unit.md)
  with `popsize_var`: the non-sampled part is now predicted at the
  non-sampled covariate mean
  $`\bar x_{d,r} = (N_d \bar X_d - n_d \bar x_{d,s})/(N_d - n_d)`$
  (before: at $`\bar X_d`$), and its posterior SD includes
  $`\sigma_e^2/(N_d - n_d)`$.

- [`hb_unit()`](https://ridsonap.github.io/fastsae/reference/hb_unit.md)
  (Gaussian, non-spatial): INLA’s posterior mode search now starts at
  the REML estimates of the nested error model. Under the vague
  log-gamma prior it could stall near $`\sigma^2_u = 0`$ depending on
  the number of threads (e.g. $`\hat\sigma^2_u \approx 1`$ instead of
  $`\approx 68`$ on `cornsoybean`). A `control.mode` passed through
  [`inla_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)
  takes precedence.

- `benchmark_sae(method = "optimal")`: the adjustment is now
  proportional to $`w_i\,\mathrm{MSE}_i`$, the minimiser of
  $`\sum_i (\theta_i' - \hat\theta_i)^2/\mathrm{MSE}_i`$ subject to the
  benchmark (before: $`\mathrm{MSE}_i / w_i`$). Results are unchanged
  for equal weights.

- [`eblup_twofold()`](https://ridsonap.github.io/fastsae/reference/eblup_twofold.md)
  with one subarea per area (only `s2v + s2u` identified): Fisher
  scoring uses the pseudo-inverse of the singular information.

- [`eblup_bhf()`](https://ridsonap.github.io/fastsae/reference/eblup_bhf.md):
  the bootstrap MSE now uses the same population means (`popnmean_xpop`)
  as the estimate.

- `eblup_fh(transform = "log")`: documentation now describes the
  back-transformation accurately (conditional-mean correction
  $`\exp(\tilde\mu + g_1/2)`$ with delta-method MSE; not the full
  Slud-Maiti second-order correction), and self-benchmarking warns that
  the constraint holds on the log scale only.

- `eblup_fh(transform = "log")`: unsampled domains now always use the
  synthetic correction $`\exp(x'\hat\beta + \hat\sigma^2_u/2)`$, even
  when their `vardir` is supplied.

- HB models: reported variance components are posterior means of the
  variance, $`E[1/\tau \mid y]`$, instead of $`1/E[\tau \mid y]`$.

- [`eblup_stfh()`](https://ridsonap.github.io/fastsae/reference/eblup_stfh.md)
  analytical MSE now includes the second-derivative term $`g_4`$ (V is
  nonlinear in $`\rho_1, \rho_2`$), as the SFH MSE does:
  `mse = g1 + g2 + 2 g3 - g4`, with `g4` returned. Simulation: relative
  bias -1.7% (was +2.4% without $`g_4`$, -6.1% without $`g_3`$).

- [`eblup_sfh()`](https://ridsonap.github.io/fastsae/reference/eblup_sfh.md):
  when the classical Fisher scoring fails or ends with
  $`\hat\sigma^2 < 0`$, the (RE)ML likelihood is maximised under the
  constraints (as for
  [`eblup_stfh()`](https://ridsonap.github.io/fastsae/reference/eblup_stfh.md)),
  also inside the bootstrap. An estimate at $`\rho = -0.999`$ is
  reported as $`-1`$ only when $`I + W`$ is non-singular (not for
  bipartite graphs such as a rook grid, where it made the bootstrap MSE
  explode).

- [`eblup_sfh()`](https://ridsonap.github.io/fastsae/reference/eblup_sfh.md)
  unsampled (kriged) domains: the analytical MSE now includes $`g_3`$
  (exact derivative of the kriging predictor) and the second-order bias
  of the plug-in $`g_1`$ (and the ML bias term):
  `g1 + g2 + g3 - tr(I^{-1} Hess g1)/2`. Simulation (REML): relative
  bias -3.8% (was -7.5%).

- [`eblup_sfh()`](https://ridsonap.github.io/fastsae/reference/eblup_sfh.md)
  /
  [`eblup_stfh()`](https://ridsonap.github.io/fastsae/reference/eblup_stfh.md)
  analytical MSE: a correlation estimated on the boundary
  ($`|\hat\rho| \ge 0.999`$, or reported as $`-1`$) is treated as fixed
  – only the information block of the other parameters is inverted –
  because the Prasad-Rao expansion is not valid there and could explode.
  The SE of such a $`\rho`$ is `NA`. For SFH this departs from
  [`sae::mseSFH`](https://rdrr.io/pkg/sae/man/mseSFH.html) in that case
  only (interior estimates still agree exactly).

- `eblup_sfh(method = "ML")` analytical MSE: $`g_3`$ and the ML bias
  term $`b_{ML} = \frac12 I^{-1} h`$ now use the ML information instead
  of the REML one (as
  [`sae::mseSFH`](https://rdrr.io/pkg/sae/man/mseSFH.html) does). With
  $`\hat\sigma^2 \approx 0`$ and $`\hat\rho \approx 1`$ the REML choice
  made the bias term explode (MSE up to $`10^6`$ in simulations); with
  few domains the ML information also reduced the median MSE bias
  (e.g. +3.6% vs +6.5%, $`m = 36`$). The ML MSE therefore no longer
  equals [`sae::mseSFH`](https://rdrr.io/pkg/sae/man/mseSFH.html); REML
  is unchanged.

- [`eblup_stfh()`](https://ridsonap.github.io/fastsae/reference/eblup_stfh.md)
  is much faster: REML scores and Fisher information are computed from
  the block structure of $`V`$ (domain blocks plus a rank-$`D`$ spatial
  part) without forming $`M \times M`$ matrices ($`O(T D^3 + M(D+p)^2)`$
  instead of $`O(M^3)`$ per iteration; identical results). Fitting the
  $`40 \times 11`$`mys_panel` went from about 7.7 s to 0.4 s, and the
  parametric bootstrap (which re-estimates the model in every replicate)
  from about 2.4 s to 0.4 s per replicate.

- `hb_area(self_benchmark = TRUE)` now benchmarks every posterior draw
  by default (`benchmark_method = "posterior"`; difference / logit /
  ratio adjustment by family), so the constraint is exact for every
  likelihood and SD / credible intervals refer to the benchmarked
  estimator. Benchmark targets are computed on the parameter scale
  (proportions for binomial, rates for counts with exposure); before,
  the WFQ targets used counts while the estimates were proportions /
  rates. `benchmark_method = "wfq"` keeps the previous augmented-model
  approach.

- `eblup_fh(transform = "log", self_benchmark = TRUE)`: ratio
  benchmarking of the back-transformed estimates on the original scale
  (exact), instead of a WFQ constraint that only held on the log scale.

- [`benchmark_sae()`](https://ridsonap.github.io/fastsae/reference/benchmark_sae.md):
  the group MSE used by the hierarchical “optimal” method now includes
  the covariances between area prediction errors (full MSE matrix for
  FH, SFH, STFH and INLA area-level HB fits), and the output has a new
  `mse_benchmarked` column, the first-order MSE of the benchmarked
  estimates ($`mse_i - 2 k_i (Mw)_i + adj_i^2`$; within about 1% of the
  empirical MSE in a simulation).

- [`diagnose()`](https://ridsonap.github.io/fastsae/reference/diagnose.md):
  the Brown et al. (2001) goodness-of-fit statistic $`W`$ now uses an
  upper-tail $`\chi^2`$ decision (as
  [`emdi::compare.fh`](https://rdrr.io/pkg/emdi/man/compare.html)); the
  former two-sided interval flagged every correct FH fit in a
  simulation, because $`W`$ is conservative. The bias regression of the
  direct on the model estimates is now descriptive
  (`brown_test$bias_regression`: `alpha`, `beta`, standard errors,
  $`R^2`$) instead of an F test of $`\alpha = 0, \beta = 1`$, which
  rejected 93% of correct FH fits (the slope exceeds 1 for shrinkage
  estimators); it no longer affects `status`. `simulation_metrics$rrmse`
  is the relative RMSE across domains (`mean_rrmse` duplicated
  `mean_abs_rb`). DIC/WAIC/marginal likelihood are now also read for
  [`hb_unit()`](https://ridsonap.github.io/fastsae/reference/hb_unit.md)
  and
  [`hb_twofold()`](https://ridsonap.github.io/fastsae/reference/hb_twofold.md)
  fits, and CPO/PIT are only used when they are per domain.

- [`spatial_test()`](https://ridsonap.github.io/fastsae/reference/spatial_test.md):
  the FH vs spatial FH comparison uses ML fits, so the likelihood ratio
  and AIC/BIC are computed from ML log-likelihoods at ML estimates
  (before: at REML estimates). The recommendation now rests on a single
  test at level `alpha` (`decision_test`: the likelihood ratio test when
  `vardir` is given, otherwise the LM-error test); requiring any of
  three tests recommended a spatial model for 9% of non-spatial data
  sets at `alpha = 0.05`. The Moran’s I, LM-error and LM-lag statistics
  are unchanged and agree with `spdep`.

- [`diagnose()`](https://ridsonap.github.io/fastsae/reference/diagnose.md):
  standardised residuals are $`(y - \hat\theta)/\sqrt{vardir - MSE}`$
  (the variance of $`y - \hat\theta`$) instead of $`/\sqrt{vardir}`$,
  which shrank them (sd 0.76 under a correct FH model) and hid outliers
  in the residual and Q-Q plots.

- Reliability tiers
  ([`add_reliability_flags()`](https://ridsonap.github.io/fastsae/reference/add_reliability_flags.md),
  [`export_sae()`](https://ridsonap.github.io/fastsae/reference/export_sae.md),
  [`map_sae()`](https://ridsonap.github.io/fastsae/reference/map_sae.md),
  [`compare_sae()`](https://ridsonap.github.io/fastsae/reference/compare_sae.md),
  `autoplot(type = "rse")`) default to the BPS-Statistics Indonesia
  criteria: RSE \<= 25% reliable, 25% \< RSE \<= 50% use with caution,
  RSE \> 50% unreliable (was 20%/30%, documented as the BPS standard).
  `map_sae(type = "reliability")` now uses its `thresholds` argument
  (they were ignored).

- `benchmark_sae(type = "total")` in a hierarchical benchmark without
  group targets: the group aggregates are weighted totals (they were
  weighted means, on the wrong scale for stage 1).

#### Breaking changes and clean-up

- [`benchmark()`](https://rdrr.io/pkg/tipsae/man/benchmark.html), a
  plain alias of
  [`benchmark_sae()`](https://ridsonap.github.io/fastsae/reference/benchmark_sae.md),
  is removed; use
  [`benchmark_sae()`](https://ridsonap.github.io/fastsae/reference/benchmark_sae.md).
- The leftover handling of `df_ebp`/`ebp` result tables (from the
  removed `ebp_area()`) is removed from the methods and utilities;
  results are in `df_hb` (HB) or `df_eblup` (EBLUP).
- [`diagnose()`](https://ridsonap.github.io/fastsae/reference/diagnose.md)
  returns `brown_test$bias_regression` instead of `brown_test$wald_test`
  and `simulation_metrics$rrmse` instead of `mean_rrmse` (see above).

#### Other fixes in the utilities

- [`map_sae()`](https://ridsonap.github.io/fastsae/reference/map_sae.md)
  works with any name of the `sf` geometry column (it assumed
  `"geometry"`), labels HB fits correctly (they were labelled “EBP-”),
  and stops with a clear message for panel or sub-area tables with
  repeated domains instead of mis-merging them.
- [`compare_sae()`](https://ridsonap.github.io/fastsae/reference/compare_sae.md)
  and
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  plots label panel and sub-area rows uniquely (repeated domains made
  some plots fail or overlap);
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  of an unnamed list of models names them.
- [`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)
  show the estimation method of
  [`eblup_fh()`](https://ridsonap.github.io/fastsae/reference/eblup_fh.md)
  and
  [`eblup_twofold()`](https://ridsonap.github.io/fastsae/reference/eblup_twofold.md)
  (REML or ML; was “eblup”) and describe `eblup_stfh(model = "S")`
  correctly.
- [`export_sae()`](https://ridsonap.github.io/fastsae/reference/export_sae.md)
  lists both variance components of two-fold models in the model
  summary.

#### Performance

- [`benchmark_sae()`](https://ridsonap.github.io/fastsae/reference/benchmark_sae.md):
  the MSE covariance between areas is applied as a matrix-free operator
  (FH: $`O(n p^2)`$, no $`n \times n`$ matrices). FH with $`m = 2000`$:
  9.4 s / 490 MB -\> 0.002 s / 1.7 MB; $`m = 20000`$ now runs in about
  0.01 s.
- [`eblup_bhf()`](https://ridsonap.github.io/fastsae/reference/eblup_bhf.md):
  the REML/ML fit and the bootstrap MSE now run in C++
  (`src/eblup_bhf_core.cpp`). The likelihood is profiled to
  $`\gamma = \sigma_u^2/\sigma_e^2 \ge 0`$ and evaluated from
  within-domain cross-products and domain means ($`O(Dp^2)`$ per
  evaluation, independent of the number of units); the maximum is found
  from a log-grid scan plus Brent refinement, so $`\sigma_u^2 = 0`$ is
  handled exactly. Estimates agree with `lme4`/`sae` up to their
  optimiser tolerance (deviance never worse), and the bootstrap uses the
  same random draws as before, so MSEs agree with the previous
  `lme4`-based implementation to about 1e-5 for a given seed. Bootstrap:
  cornsoybean $`B = 200`$ 0.64 s -\> 0.005 s; $`D = 300`$, $`n = 6000`$,
  $`B = 100`$ 0.96 s -\> 0.035 s. `n_threads` now runs replicates in
  parallel (results do not depend on it); for small models R’s
  sequential RNG dominates the run time.
- [`eblup_bhf()`](https://ridsonap.github.io/fastsae/reference/eblup_bhf.md)
  output: `fit$lme` (the `lmerMod` object) is replaced by
  `fit$residuals` (unit-level conditional residuals, used by
  [`residuals()`](https://rdrr.io/r/stats/residuals.html)); `estcoef`
  standard errors are $`(X'\hat V^{-1}X)^{-1}`$ as in `lme4`.
  Rank-deficient designs now give an explicit error.
- [`eblup_sfh()`](https://ridsonap.github.io/fastsae/reference/eblup_sfh.md)
  (fit, analytical MSE, kriging and the bootstrap refits) no longer
  forms any $`D \times D`$ inverse or $`D^3`$ product. With
  $`A = (I - \rho W')(I - \rho W)`$ and
  $`M = A + \sigma^2 E'\Psi^{-1}E`$ (both sparse), $`V^{-1}`$,
  $`GV^{-1}`$, $`g_1 = \sigma^2 \mathrm{diag}(M^{-1})`$, the
  scores/information, $`g_3`$, $`g_4`$, the ML bias and the kriging
  terms follow from Cholesky solves with $`A`$ and $`M`$
  (`src/sar_engine.h`; Eigen sparse $`LDL'`$ with AMD ordering, dense
  factorisations when $`W`$ itself is dense, e.g. distance weights).
  Results are unchanged (relative differences below 2e-10). Analytical:
  $`D = 400`$ 2.2 s -\> 0.10 s, $`D = 1600`$ 80 s -\> 2.5 s,
  $`D = 2500`$ about 6 s; bootstrap ($`D = 400`$, $`B = 20`$): 16 s -\>
  1.3 s. `RcppEigen` is added to `LinkingTo`.
- [`eblup_stfh()`](https://ridsonap.github.io/fastsae/reference/eblup_stfh.md):
  the final step (Fisher information at the estimates, SEs,
  $`g_1`$-$`g_4`$ of the analytical MSE) no longer forms $`M \times M`$
  matrices ($`M = DT`$); it uses the same structured representation of
  $`V`$ as the iterations ($`A^{-1}`$ block diagonal,
  $`P = A^{-1} - U W U'`$ with $`U = [A^{-1}Z, V^{-1}X]`$ whose
  $`A^{-1}Z`$ part has one non-zero per row). The scoring step also
  exploits the diagonal blocks of $`U'Z`$ and $`U'EU`$ and the sparsity
  of $`C = 2\rho W'W - W - W'`$. Results are unchanged (relative
  differences below 3e-11). $`T = 5`$: $`D = 250`$ 28 s -\> 1.6 s,
  $`D = 500`$ about 4 min -\> 10 s, $`D = 1000`$ 75 s.
- [`create_weights()`](https://ridsonap.github.io/fastsae/reference/create_weights.md):
  distances for `knn`, `distance` and `idw` are computed in blocks of
  rows (at most about 1M distances at a time) instead of a full
  $`D \times D`$ distance matrix, and the weight matrix is filled once
  from the neighbour pairs (no intermediate 0/1 matrix and no
  normalisation copies; `idw` styles are applied per block). The result
  is the same matrix as before (identical, or within 1e-16 for `idw`).
  Peak memory at $`D = 10000`$: queen 2.3 GB -\> 0.8 GB, knn 3.1 GB -\>
  0.8 GB, distance 3.8 GB -\> 0.8 GB, i.e. essentially the returned
  matrix itself (764 MB).

#### Bug fixes

- [`create_weights()`](https://ridsonap.github.io/fastsae/reference/create_weights.md)
  with a coordinate matrix whose columns are named `X`/`Y` (e.g. from
  [`sf::st_coordinates()`](https://r-spatial.github.io/sf/reference/st_coordinates.html))
  failed with “undefined columns selected”; coordinate column names are
  now matched case-insensitively.

#### Major New Features & Enhancements

- **Two-Fold Subarea EBLUP (`eblup_twofold`)**:
  - Implements nested two-fold sub-area level models (Torabi &
    Rao, 2014) with area effects $`v_i`$ and sub-area effects
    $`u_{ij}`$.
  - Ultra-fast Fisher-scoring algorithms in C++ for REML and ML
    estimation.
  - Analytical Prasad-Rao MSE ($`g_1 + g_2 + 2 g_3`$ for REML, plus the
    bias correction for ML) with exact matrix derivation and
    multi-threaded parametric bootstrap MSE.
  - Automatic prediction for non-sampled subareas using area-level
    random effects and regression synthetic estimation.
- **Hierarchical Bayes Two-Fold Subarea Models (`hb_twofold`)**:
  - Bayesian nested two-fold sub-area models via `INLA`.
  - Supports Gaussian, Binomial (logit link), and Poisson (log rate)
    likelihoods.
  - Area-level random effects support IID or spatial structures (BYM2,
    Besag).
  - Simultaneous posterior estimation for subareas (`df_hb`) and
    weighted area aggregates (`df_area`) with full Monte Carlo
    uncertainty propagation.
- **Hierarchical Bayes Area-Level Models (`hb_area`)**:
  - Fast Bayesian inference powered by Integrated Nested Laplace
    Approximations (`INLA`).
  - Supports 6 probability distributions for non-Gaussian area
    indicators: Gaussian, Binomial (logit link with trials), Poisson
    (log rate with exposure/offsets), Negative Binomial, Beta, and
    Gamma.
  - Spatial random effects: BYM2 (Besag-York-Mollié 2 with PC-priors)
    and Besag/CAR models.
  - Automatic handling of unsampled domains and domain-specific variance
    predictions.
- **Hierarchical Bayes Unit-Level Models (`hb_unit`)**:
  - Bayesian unit-level Battese-Harter-Fuller (BHF) and generalized
    linear mixed models via `INLA`.
  - Supports Gaussian, Binomial (binary unit outcomes), and Poisson
    (count outcomes).
  - Domain random effects with optional spatial structures (BYM2,
    Besag).
  - Finite population adjustments when domain population sizes ($`N_d`$)
    are supplied.
  - Seamless handling of both domain-level population means and
    micro-level census data.
- **Two-Stage Hierarchical Benchmarking (`benchmark_sae`)**:
  - Unified calibration framework for multi-level administrative
    structures (e.g., Domain -\> Province -\> National).
  - Implements Ratio, Difference, Optimal Variance-Weighted Quadratic
    Loss (Datta et al. 2011; Steorts et al. 2014), and Constrained Logit
    benchmarking (Berg & Fuller 2014).
- **Comprehensive Diagnostics & Visualization**:
  - [`diagnose()`](https://ridsonap.github.io/fastsae/reference/diagnose.md):
    Brown et al. (2001) goodness-of-fit and bias diagnostics, residual
    analysis, spatial autocorrelation tests (Moran’s I), and official
    statistical reliability classification (BPS standards).
  - [`compare_sae()`](https://ridsonap.github.io/fastsae/reference/compare_sae.md):
    Tabular model comparison across information criteria and error
    metrics.
  - [`map_sae()`](https://ridsonap.github.io/fastsae/reference/map_sae.md):
    Spatial thematic choropleth maps integrated with `sf`.
  - [`export_sae()`](https://ridsonap.github.io/fastsae/reference/export_sae.md):
    Export estimation results with full statistical metadata to Excel
    (`openxlsx`/`writexl`) and CSV.
- **Simulation Utilities**:
  - [`sim_area_data()`](https://ridsonap.github.io/fastsae/reference/sim_area_data.md),
    [`sim_spatial_weights()`](https://ridsonap.github.io/fastsae/reference/sim_spatial_weights.md),
    and
    [`sim_series_data()`](https://ridsonap.github.io/fastsae/reference/sim_series_data.md).
- **Reliability & Performance**:
  - High-performance C++ core (`RcppArmadillo` and OpenMP) with strict
    validation checks.

## fastsae 0.1.0

CRAN release: 2026-09-30

- Initial development version.
- High-performance Area-level Fay-Herriot (FH), Spatial FH (SFH), and
  Spatio-Temporal FH (STFH) models.
- Unit-level Battese-Harter-Fuller (BHF) models via Henderson III and
  REML.
