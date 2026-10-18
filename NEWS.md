# fastsae 1.0.0

### Estimation engines for `hb_area()`, `hb_unit()` and `hb_twofold()`
* New argument `engine`, taking an engine object from the new constructors `inla_engine()` and
  `jax_engine()` (or the strings `"inla"`, the default, and `"jax"` for their default settings).
  Each constructor holds all settings of its engine, so `hb_area()` itself only describes the
  model:
  * `inla_engine(strategy, link, scale_model, prior_prec, prior_phi, prior_rho, prior_prec_time,
    prior_rho_time, prior_prec_subarea, draws, seed, ...)`, where `...` goes to `INLA::inla()`.
    `control.*` lists given there are merged with the package defaults (before, they clashed with
    them in `hb_area()` and replaced them in `hb_unit()`/`hb_twofold()`). Priors and `draws` left at `NULL`
    take the default of the model function, so the defaults of each function are unchanged;
    settings a function does not use are ignored with a warning.
  * `jax_engine(warmup, samples, chains, device, seed, prior_beta, robust, keep_draws, progress)`
    fits the same area-level models by full MCMC (NumPyro NUTS on JAX: CPU, NVIDIA GPU or the
    Apple silicon GPU via `device = "metal"`) through
    'reticulate' (Suggests); the engine formerly shipped in the separate 'fastsaegpu' package now
    lives in fastsae. It supports all six families, `bym2`/`bym`/`besag`/`generic1` (alias
    `"leroux"`), `rw1`/`ar1`/`iid` and every `st_interaction`, plus `prior_beta = "horseshoe"` and
    `robust = TRUE` (Student-t area effects). Results have the same structure as the INLA engine,
    with the posterior draws in `draws`, the largest split R-hat in `rhat_max` (`convergence` is
    `rhat_max < 1.05`), and DIC/WAIC in `goodness`. Posterior self-benchmarking (also by
    `benchmark_group`) adjusts the MCMC draws exactly as for INLA, and `benchmark_sae()` uses the
    draws for the MSE matrix.
  * The engine settings are returned in `engine_spec`; `engine` holds the engine name.
* `hb_twofold()` also takes `engine = jax_engine()` (or `"jax"`): the two-fold model with an
  iid, Besag or BYM2 area effect and an iid sub-area effect is sampled by NumPyro NUTS, for the
  Gaussian, binomial and Poisson families, and the area aggregates in `df_area` come from the
  MCMC draws. Against INLA on simulated data (20 areas) the area means correlated above 0.999
  and differed by a median of 0.02-0.07 posterior SD (the priors differ: half-normal standard
  deviations instead of PC priors); a fit took 7-11 s against 11-31 s for INLA with the area
  aggregates. The result has the structure of the INLA engine plus `draws`, `rhat_max`, `device`
  and `elapsed_seconds`.
* `hb_unit()` also takes `engine = jax_engine()` (or `"jax"`): the nested error model
  (Gaussian, Bernoulli or Poisson units; iid, Besag or BYM2 domain effects) is sampled by
  NumPyro NUTS from the sampled units, and the domain means, including the finite population
  mean with `popsize_var`, are summarised from the MCMC draws. Checked against a Stan fit of the
  same model and priors (30 domains; domain means within 0.04 posterior SD, SDs within 5%), and,
  with their own priors, against `hbsae::fSAE(method = "HB")` (Gaussian; within 0.1 SD) and the
  INLA engine (within 0.06 SD). On `cornsoybean` (12 counties) the priors matter: the county SDs
  are 1.1-1.4 times INLA's and 0.8-1.0 times hbsae's.
* The NumPyro backend returns the draws of the residual variance (`sigma2_e`) and can return the
  draws of selected rows only (`draw_rows`).
* Convergence diagnostics for `jax_engine()` fits (Vehtari et al., 2021), identical to
  `posterior::rhat()`, `ess_bulk()` and `ess_tail()` (checked to 1e-10): the rank-normalized split
  R-hat and bulk/tail ESS of every hyperparameter, coefficient and area parameter in
  `mcmc_diagnostics`, with `ess_min`, `divergences` (divergent transitions) and `treedepth_hits`.
  `rhat_max` used to cover only the coefficients and the area parameters, with the classic split
  R-hat, so a variance whose chains had not mixed went unnoticed: on `mys` without `vardir`
  (2 x 300 draws) the residual variance had R-hat 1.13 and a bulk ESS of 14, with 13 divergent
  transitions. Fits warn for R-hat >= 1.05, ESS below 100 per chain and divergent transitions
  (class `fastsae_mcmc_warning`), inform for R-hat between 1.01 and 1.05, and `print()` shows
  the ESS and divergences.
* `autoplot(fit, type = "trace")` and `type = "rank"` plot the chains of a `jax_engine()` fit:
  the hyperparameters, coefficients and the area with the largest R-hat by default, or any
  `pars`, with R-hat and ESS in the panel labels.
* `jax_engine()`: the WAIC of a Gaussian model with an unknown residual variance (`hb_area()`
  without `vardir`, `hb_unit()`) used the posterior mean of the residual SD in every draw; it
  now uses the SD of each draw (on `mys`: 129.3 instead of 135.9).
* `jax_engine()`: the `temporal = "ar1"` effect sums to zero over the periods (as `"rw1"` already
  did) when the model has an intercept, which then carries the mean level. Over a few periods
  the level of a stationary AR(1) effect traded off against the intercept along a ridge whose
  width depends on rho, and NUTS reported divergent transitions: on `sim_panel` (5 years,
  4 x 1000 draws, 4 seeds) 30-201 before and 1-6 after, with the area estimates 0.11 instead of
  0.21 posterior SD from INLA; with `st_interaction = "type2"`, 14-57 before and 0 after. With
  `spatial = "bym2"` one of 4 seeds still gave 38 (3-4 for the others; 30-154 before), which a
  longer `warmup` reduces. The INLA engine keeps its unconstrained AR(1) effect, so the
  intercepts of the two engines differ by the mean of the time effects (the area estimates are
  not affected by this).
* `hb_unit()` checks that the unit responses are 0 or 1 for `family = "binomial"` and non-negative
  integers for `"poisson"`.
* `print()` failed on a fit whose call holds braces, e.g. from `do.call()` with a function object:
  the call was read as cli markup.
* `jax_engine(robust = TRUE)` is documented as what it does: Student-t random effects (it was
  described as Student-t sampling errors for the Gaussian family).
* `hb_unit()` and `hb_twofold()` take the same `engine = inla_engine(...)`. In
  `hb_twofold()`, `prior_prec_area` is now `inla_engine(prior_prec = )`, `n_samples` is
  `inla_engine(draws = )`, and the posterior draws for the area aggregates honour `seed`. Their
  results carry `engine` and `engine_spec` too.
* `jax_engine(device = "metal")` now samples on the Apple GPU. jax-metal cannot compile two
  operations of NumPyro's NUTS sampler (the `popcnt` bit count, and `lax.cond` statements whose
  outputs are constants), so fits used to fall back to the CPU. fastsae now replaces both with
  equivalent operations that Metal can compile, and steps the chains from Python; the draws match
  the CPU sampler. On an Apple M4, Metal is faster than the CPU for large models (about 2x at
  10,000 domains) and slower for small ones, where a fixed cost of about two minutes per fit
  dominates (see the "hb_area() engines" benchmark article). `device = "auto"` keeps using the CPU
  on macOS.
* New `setup_numpyro_env()` and `check_numpyro_available()` to install and check the Python side.
  Nothing is installed or downloaded unless `setup_numpyro_env()` is called. The check looks only at
  an initialised Python, `RETICULATE_PYTHON` / `RETICULATE_PYTHON_ENV` and the `r-numpyro-gpu`
  environment, so it answers in about a second instead of searching the whole system for Python
  installations (which took over 20 s on Windows).
* The frequentist lme4 GLMM route (`method = "laplace"` for non-spatial binomial and Poisson
  models) has been removed; such models are now fitted by INLA. 'lme4' moved to Suggests.
* Deprecated: `method` (`method = "laplace"` becomes `engine = inla_engine(strategy = "laplace")`),
  and passing `strategy`, `link`, `scale_model`, `prior_*`, `benchmark_draws`, `benchmark_seed`,
  `n_samples` or `INLA::inla()` arguments directly to `hb_area()`, `hb_unit()` or `hb_twofold()`.
  They still work, with a warning, and are forwarded to `inla_engine()`.
* `print()`/`summary()` show the engine, and the self-benchmarking line now names the method used
  (posterior draws or Wang-Fuller-Qu) instead of always reporting Wang-Fuller-Qu.

### Unsampled domains and missing periods in `eblup_stfh()`
* `eblup_stfh()` no longer requires a complete panel. Cells whose response is `NA` -- missing
  periods of a sampled domain, or domains that are never sampled -- stay in the data and are
  predicted. The model is defined on the full domain x time grid (the SAR effect on all domains
  of `W`, the AR(1) effect on all periods); the variance components are estimated by REML from
  the sampled cells (the marginal likelihood with the missing cells integrated out), and every
  cell gets the best linear unbiased predictor given the sampled data (Saei and Chambers, 2005;
  Benedetti et al., 2024). A missing period borrows from the other periods of its domain and
  from the neighbours; a domain never sampled is predicted from its neighbours (spatial kriging).
* Analytical MSE for those cells: $g_1 + g_2 + g_3$, with $g_3$ from the exact derivative of
  the predictor in the variance components. The further correction of the plug-in $g_1$ is not
  applied: it is numerically unstable when $\rho_1$ approaches $\pm 1$ or a variance component is
  weakly identified, and in Monte Carlo $g_1 + g_2 + g_3$ has a relative bias within a few
  percent and close to 95% coverage for these cells. The parametric
  bootstrap now draws the random effects on the full grid, refits on the sampled cells and
  covers every cell. `df_eblup` gains a `sampled` column; `print()` reports the number of
  predicted cells. The auxiliary variables must be known in every cell, and `W` must include
  the unsampled domains.
* Results on complete panels are unchanged (bit-for-bit). The predictions and every MSE term
  are tested against a dense implementation, and `benchmarks/sim_stfh_unsampled.R` checks the
  MSE estimators by Monte Carlo.
* The examples now use the full `mys_panel` and `mys_proxmat` (domains 21 and 25 are unsampled).

### Area means in `eblup_twofold()`
* `eblup_twofold()` now also estimates the area means $\theta_i = \sum_j w_{ij}\theta_{ij}$,
  returned in `df_area` (`eblup_area`, `mse_area`, `rse_area`, `n_subareas`, `n_sampled`). The
  new argument `weight` gives the sub-area weights (normalised within each area, equal by
  default, as in `hb_twofold()`); `compute_area = FALSE` skips them.
  * The analytical MSE is the second-order estimator of the sub-areas applied to the linear
    target ($g_1 + g_2 + 2g_3$, minus the ML bias term). It includes the covariance of the
    sub-area errors through the area effect and $\hat\beta$, so it is not
    $\sum_j w_{ij}^2 \mathrm{mse}_{ij}$. Areas without a sampled sub-area get the synthetic
    estimator. `mse_type = "bootstrap"` computes the area MSE from the same replicates.
  * Tested against a dense implementation (finite-difference $g_3$, REML information) and, with
    a weight on one sub-area, against the sub-area MSE. `benchmarks/validation/twofold_area_mse.R`
    checks the MSE estimators by Monte Carlo.

### MSE argument, convergence and input checks
* All EBLUP functions choose the MSE with `mse_type`, which includes `"none"` (no MSE; `mse` and
  `rse` are `NA`): `eblup_fh()` (`"analytical"`, `"none"`), `eblup_sfh()` (`"analytical"`,
  `"pbmse"`, `"npbmse"`, `"none"`), `eblup_stfh()` and `eblup_twofold()` (`"analytical"`,
  `"bootstrap"`, `"none"`) and `eblup_bhf()` (`"bootstrap"`, `"none"`). The chosen value is
  returned in `mse_type`.
  * `eblup_sfh(mse_method = )` and `compute_mse = ` of `eblup_stfh()` and `eblup_bhf()` are
    deprecated and still work with a warning (`compute_mse = FALSE` is `mse_type = "none"`).
  * The MSE is now computed by default: `eblup_stfh()` gives the analytical MSE and `eblup_bhf()`
    the bootstrap MSE (both gave none before). Use `mse_type = "none"` for the previous output.
* Convergence is reported as a warning that can be caught. Before, `eblup_fh()`, `eblup_sfh()`,
  `eblup_stfh()` and `eblup_twofold()` only printed a message and `eblup_bhf()` said nothing.
  INLA fits of `hb_area()`, `hb_unit()` and `hb_twofold()` now read `convergence` from the INLA
  mode search (`ok`, `mode$mode.status`; it was always `TRUE`), and `hb_area(engine = "jax")`
  warns when the largest split R-hat is 1.05 or more.
* Inputs are checked before fitting, with a clear message:
  * `vardir` must be finite and positive in every sampled area (a missing value gave an error
    from the compiled code);
  * linearly dependent covariates stop the fit and name the aliased column (the compiled code
    used an approximate solution of a singular system);
  * there must be more sampled areas than coefficients (all MSEs were `NA`), and at least two
    sampled areas in `eblup_twofold()` and `eblup_bhf()`;
  * `W` of `eblup_sfh()` and `eblup_stfh()` must be non-negative; a zero row (a domain without
    neighbours) or rows that do not sum to 1 give a warning, since the SAR model with
    `rho` in (-1, 1) assumes a row-standardised `W`, as 'sae' requires;
  * the same checks now cover the HB models. A missing covariate is an error in `hb_area()`,
    `hb_twofold()` and in the `Xpop` means of `hb_unit()`: INLA treats it as 0, so the estimate
    of that area was wrong without any message. `hb_twofold()` replaced a missing, zero or
    negative `vardir` of a sampled subarea by 1, and `hb_area()` stopped with an INLA error on a
    missing one; both now stop with the message of the EBLUP functions. `hb_area()` also rejects
    linearly dependent covariates;
  * a missing population mean in `eblup_bhf()` names the domain (it gave an internal error on
    the number of rows);
  * duplicated `domain` values give a warning in `eblup_fh()`, `eblup_sfh()` and `hb_area()`
    without `time` (in `hb_area()` such rows share one area effect).
* A fixed `seed` of `eblup_sfh()`, `eblup_stfh()`, `eblup_bhf()` and `eblup_twofold()` no longer
  changes the random number stream of the session; the bootstrap MSEs are unchanged.
* `eblup_twofold()` uses `B = 100` and `seed = -1` (no seed) by default, as the other EBLUP
  functions (it used `B = 200` and `seed = NULL`; `NULL` still means no seed).
* Gaussian models of `hb_area()`, `hb_unit()` and `hb_twofold()`, and the Leroux and SLM effects
  of `hb_area()` for any family, are fitted in INLA's `"classic"` mode (Rue, Martino and Chopin,
  2009); the other models keep INLA's default `"compact"` mode (van Niekerk et al., 2023), which
  leaves the linear predictors out of the latent field and is designed for data-rich models.
  Gaussian models of `hb_area()` with a `"type1"` to `"type4"` interaction and three
  hyperparameters (a Besag or iid area effect, a random walk or iid in time) also integrate over a
  grid of hyperparameters. Larger models (`temporal = "ar1"`, or a BYM, BYM2, Leroux or SLM area
  effect) keep INLA's default integration: the grid grows exponentially with the number of
  hyperparameters (with AR(1), 30 s and 13 min for 12 areas against about 1 s) and was validated
  only for three. The choices rest on comparisons with the exact posterior (numerical
  integration over the variances) and with long Stan runs with the same model and priors (scripts
  in `benchmarks/validation`):
  * Gaussian: the compact mode overstated the posterior SDs by up to 18% (`hb_area()`, when the
    area variance is small relative to `vardir`; it was accurate on `mys`), 8-15% (`hb_unit()`,
    `cornsoybean`) and 10% (`hb_twofold()`), and `hb_unit()` gave different county means on
    identical calls (up to 3 hectares and 13% in the SDs). In the classic mode the posterior means
    were within 0.07 posterior SD and the SDs within 0.91-1.03 times the reference, at the same
    speed. With a space-time interaction both modes overstated the SDs by 6-13%; classic with the
    grid matched Stan (means within 0.04 SD, SDs within 2%) and takes 5-8 times longer.
  * Binomial, Poisson, beta and gamma: both modes matched Stan (compact: means within 0.16
    posterior SD, SDs within 0.95-1.08 times), and the classic mode was much slower for many
    areas or units (beta with 2000 areas: 32 s against 1.3 s; `hb_unit()` binomial with 4000
    units: 72 s against 2.6 s). The exception is the Leroux and SLM effects, where compact moved
    Poisson means by 0.12-0.15 SD (classic: within 0.03 SD).
  * With a small area variance, the classic summary of that variance (`random_effect_var`) is
    biased downwards while the area estimates are accurate; this is documented.
    `inla_engine(inla.mode = , control.inla = )` overrides the settings.
* `hb_area()` space-time interactions with a random walk in time are identified as in Knorr-Held
  (2000). INLA's default constraints left one level (rw1) or one linear trend (rw2) per domain
  with a flat prior: a domain without direct estimates had a posterior SD of 36-58 (other domains
  about 0.4) in the type IV interaction (rw1, rw2) and the type II interaction (rw2). These
  interactions are now a `generic0` term with the scaled Kronecker structure and all the
  constraints (sums over domains per period and over periods per domain, and linear trends per
  domain for rw2); they match Stan with the same constraints. The `"separable"` field with a
  random walk in time, which has no spatial main effect to carry the domain levels, and the
  `"domain-specific"` rw2, which has no common trend, now stop with an error suggesting
  `temporal = "ar1"`, `"type4"` or `"type2"`.
* `hb_area(family = "nbinomial")` warns that with one count per area the overdispersion and the
  area effect are only weakly identified: INLA's posterior means were off by up to 0.6 posterior
  SD against Stan in every mode. `"poisson"` or `engine = "jax"` (MCMC) is suggested.
* CAR models (`spatial = "besag"`, `"bym"`, `"bym2"`, `"leroux"`) of `hb_area()`, `hb_unit()` and
  `hb_twofold()` use only the pattern `W > 0`. A dense `W` such as distance-decay weights becomes
  a (nearly) complete graph and an almost unstructured spatial effect; it now gives a warning.
  `mys_proxmat` is such a matrix (all 861 pairs linked): its documentation said it was based on
  contiguous boundaries, and the examples and vignettes now fit CAR models on a 4-nearest-
  neighbour graph of the `sim_area` centroids. `mys_proxmat` remains suited to the SAR models.
* `hb_unit(family = "binomial" / "poisson")` warns that the domain means are approximate. The
  non-sampled part is predicted at the mean covariates, `inv_link(xbar' beta + u_d)`, which is
  biased for the mean of `inv_link(x' beta + u_d)` over the units: in a simulation -8%
  (binomial) and -16% (Poisson) for the sampled domains, with 90% and 78% coverage of the 95%
  intervals. Unbiased estimation needs unit-level population covariates (Molina and Rao, 2010);
  see `?hb_unit`. The Gaussian model is not affected.
* `hb_twofold(weight = )` warns about missing or negative weights (they get weight 0 in
  `df_area`, as before, but silently).
* `compare_sae()` between a Hierarchical Bayes and an EBLUP model no longer calls the model with
  the smaller posterior variance "more efficient": it reports which one has the lower reported
  uncertainty and notes that a posterior variance and an MSE are not the same quantity.
* Posterior draws of the INLA engine are reproducible: `inla_engine(seed = 1L)` is the new
  default (it was `0L`, INLA's random seed), and the seed is also set on the R side, since INLA
  draws the hyperparameter configurations with R's RNG; the session's random number stream is
  restored. This affects the area aggregates of `hb_twofold()` (their SD changed by up to 9%
  between identical calls), posterior self-benchmarking in `hb_area()` and `benchmark_sae()`.
* New vignette "Sub-Area Level Estimation with Two-Fold Models" (`eblup_twofold()`,
  `hb_twofold()`); the unit-level vignette covers `hb_unit()`.

### Fixes before the 1.0.0 release
* `benchmark_sae()` stops when a `target` or `national_target` is missing or not finite, or a
  `weight` is missing. Before, a target of `NA` (for example a `weighted.mean()` of direct
  estimates that include unsampled domains) silently made every benchmarked estimate `NA`.
* `hb_area(st_interaction = "type2" / "type3" / "type4")` now fit with INLA. They always failed
  before (the interaction reused the index of a main effect, which INLA rejects), and with
  `spatial = "none"` the type III/IV interaction was silently dropped. The interaction is now a
  separate term on copied indices: one `temporal` process per domain (type II), or an ICAR field
  per period (type III) or evolving with `temporal` (type IV). Types III and IV need a
  structured `spatial` effect and stop with a clear message otherwise.
* `hb_area(spatial = "leroux", st_interaction = "separable")` fits instead of crashing inside
  INLA: the grouped Leroux field uses INLA's `besagproper2`, the same model as `generic1`.
* `inla_engine(link = )` is used. It was accepted but ignored; it now sets the link of the
  likelihood (`"probit"`, `"cloglog"` and `"loglog"` for binomial and beta), checked against the
  family. The posterior draws for self-benchmarking and `benchmark_sae()` use the matching inverse
  link, and the link is returned in `link`.
* `hb_area(spatial = "bym")`: `random_effect_var` is the variance of the spatial component (it was
  the iid one). Interaction terms no longer feed `random_effect_var`, `rho` or `rho_time`, and with
  `st_interaction = "separable"` the AR(1) group correlation is returned in `rho_time`.
* `hb_area()` objects have class `c("fastsae_hb_area", "fastsae_hb", "fastsae")`, like
  `hb_unit()` and `hb_twofold()`.
* The documentation of `inla_engine(prior_prec = )` explains the different defaults: the PC prior
  `pc.prec(1, 0.01)` of `hb_area()` and `hb_twofold()` assumes an area standard deviation of
  order 1 (proportions, rates, log scales), while `hb_unit()` keeps the vague
  `loggamma(0.01, 0.01)` because unit-level responses are often on a raw scale.
* `goodness` of `eblup_fh()`, `eblup_sfh()`, `eblup_stfh()` and `eblup_twofold()` follows one
  convention. `loglikelihood` is the full Gaussian log-likelihood at the estimates, and `AIC`/`BIC`
  are computed from it for ML and REML fits, as in 'sae' and Marhuenda, Morales and Pardo (2014):
  the restricted likelihood depends on the covariates, so it cannot compare models with different
  covariates and is not reported. `goodness` always has these three elements.
  * `eblup_twofold()`: `loglikelihood`, `AIC` and `BIC` of REML fits were based on the restricted
    log-likelihood and now use the full one; `BIC` uses the number of sampled subareas (it used
    the number of areas).
  * `eblup_stfh()`: the element `loglike` is renamed `loglikelihood`, as in the other models.
* `eblup_twofold()`: a negative second-order MSE estimate (possible with the ML bias correction)
  is set to `NA` with a warning instead of giving a `NaN` RSE.
* `eblup_stfh()` reports `convergence = FALSE` when the fit ends at invalid variance components
  (it reported `TRUE` with no estimates). The documentation explains why the model is fitted by
  REML only (as in Marhuenda, Molina and Morales 2013 and `sae::eblupSTFH()`) and that `tvalue`
  is a z-statistic.

### Statistical corrections (mathematical audit)
* `eblup_stfh()`: the analytical $g_3$ was computed with the wrong matrix orientation
  ($(K_a \Psi V \Psi K_b)_{ii}$ instead of $\psi_i^2 (K_a V K_b)_{ii}$), which is wrong whenever
  `vardir` is heterogeneous. $g_3$ now uses the exact derivative of the EBLUP,
  $\partial\tilde\theta/\partial\theta_a = \Psi P V_a P y$, and the REML information at the final estimates.
* `eblup_stfh()` bootstrap MSE: the variance components are now re-estimated in every replicate
  (before, they were held fixed, so the bootstrap missed the $g_3$ contribution), and the SAR effects
  are generated with the lower Cholesky factor (before, $\mathrm{Var}(u_1) = RR' \ne \Omega_1$).
* `eblup_stfh()` now sorts the panel internally into the domain-major order the model requires and returns
  results in the caller's row order. Before, unsorted input was silently mis-specified -- including the
  time-major `mys_panel` used in the examples, tests and vignette (where `sae::eblupSTFH()`, which assumes
  domain-major data, gave the same wrong fit).
* `eblup_stfh()`: when the classical Fisher scoring fails, or converges to a negative variance that would
  only be truncated to 0, the fit switches to constrained REML (variances >= 0, |rho| <= 0.999; projected
  scoring with step halving), which returns the actual constrained maximiser. Convergence in a simulation
  went from 63% to 99% of the replicates; fits on which the classical iteration converges to an interior
  point are unchanged (identical to `sae::eblupSTFH()`). A correlation whose variance is 0 is unidentified
  and gets an `NA` standard error.
* `eblup_stfh()` bootstrap: non-converged replicates are replaced by new draws until `B` valid replicates
  are obtained (as in the SFH bootstrap).
* `eblup_twofold()`: the analytical MSE estimator is now $g_1 + g_2 + 2 g_3$ (REML) and
  $g_1 + g_2 + 2 g_3 - b_{ML}'\nabla g_1$ (ML); before it was $g_1 + g_2 + g_3$, which ignores the downward bias
  of the plug-in $g_1$. The REML Fisher information now includes the cross-area terms.
* `eblup_sfh()` with unsampled domains: the SAR process is defined on all domains of `W`, so the sampled
  domains are now fitted with the marginal covariance $\sigma^2 \Sigma_{full}[s,s]$ (before: the SAR built
  from the sub-matrix `W[s,s]`), consistent with the kriging of the unsampled domains and the bootstrap.
* `eblup_sfh()`: an upper-boundary $\hat\rho = 0.999$ is no longer mapped to $1$ (singular for a
  row-standardised `W`); bootstrap replicates are no longer rejected when estimates exceed `1e6`
  (scale-dependent check); `npbmse` no longer fails when $\hat\sigma^2_u = 0$.
* `hb_unit()` with `popsize_var`: the non-sampled part is now predicted at the non-sampled covariate mean
  $\bar x_{d,r} = (N_d \bar X_d - n_d \bar x_{d,s})/(N_d - n_d)$ (before: at $\bar X_d$), and its posterior
  SD includes $\sigma_e^2/(N_d - n_d)$.
* `hb_unit()` (Gaussian, non-spatial): INLA's posterior mode search now starts at the REML estimates
  of the nested error model. Under the vague log-gamma prior it could stall near $\sigma^2_u = 0$
  depending on the number of threads (e.g. $\hat\sigma^2_u \approx 1$ instead of $\approx 68$ on
  `cornsoybean`). A `control.mode` passed through `inla_engine()` takes precedence.
* `benchmark_sae(method = "optimal")`: the adjustment is now proportional to $w_i\,\mathrm{MSE}_i$, the
  minimiser of $\sum_i (\theta_i' - \hat\theta_i)^2/\mathrm{MSE}_i$ subject to the benchmark (before:
  $\mathrm{MSE}_i / w_i$). Results are unchanged for equal weights.
* `eblup_twofold()` with one subarea per area (only `s2v + s2u` identified): Fisher scoring uses the
  pseudo-inverse of the singular information.
* `eblup_bhf()`: the bootstrap MSE now uses the same population means (`popnmean_xpop`) as the estimate.
* `eblup_fh(transform = "log")`: documentation now describes the back-transformation accurately (conditional-mean
  correction $\exp(\tilde\mu + g_1/2)$ with delta-method MSE; not the full Slud-Maiti second-order correction),
  and self-benchmarking warns that the constraint holds on the log scale only.
* `eblup_fh(transform = "log")`: unsampled domains now always use the synthetic correction
  $\exp(x'\hat\beta + \hat\sigma^2_u/2)$, even when their `vardir` is supplied.
* HB models: reported variance components are posterior means of the variance, $E[1/\tau \mid y]$,
  instead of $1/E[\tau \mid y]$.

* `eblup_stfh()` analytical MSE now includes the second-derivative term $g_4$ (V is nonlinear in
  $\rho_1, \rho_2$), as the SFH MSE does: `mse = g1 + g2 + 2 g3 - g4`, with `g4` returned. Simulation:
  relative bias -1.7% (was +2.4% without $g_4$, -6.1% without $g_3$).
* `eblup_sfh()`: when the classical Fisher scoring fails or ends with $\hat\sigma^2 < 0$, the (RE)ML
  likelihood is maximised under the constraints (as for `eblup_stfh()`), also inside the bootstrap. An
  estimate at $\rho = -0.999$ is reported as $-1$ only when $I + W$ is non-singular (not for bipartite
  graphs such as a rook grid, where it made the bootstrap MSE explode).
* `eblup_sfh()` unsampled (kriged) domains: the analytical MSE now includes $g_3$ (exact derivative of
  the kriging predictor) and the second-order bias of the plug-in $g_1$ (and the ML bias term):
  `g1 + g2 + g3 - tr(I^{-1} Hess g1)/2`. Simulation (REML): relative bias -3.8% (was -7.5%).
* `eblup_sfh()` / `eblup_stfh()` analytical MSE: a correlation estimated on the boundary
  ($|\hat\rho| \ge 0.999$, or reported as $-1$) is treated as fixed -- only the information block of the
  other parameters is inverted -- because the Prasad-Rao expansion is not valid there and could explode.
  The SE of such a $\rho$ is `NA`. For SFH this departs from `sae::mseSFH` in that case only (interior
  estimates still agree exactly).
* `eblup_sfh(method = "ML")` analytical MSE: $g_3$ and the ML bias term $b_{ML} = \frac12 I^{-1} h$ now use
  the ML information instead of the REML one (as `sae::mseSFH` does). With $\hat\sigma^2 \approx 0$ and
  $\hat\rho \approx 1$ the REML choice made the bias term explode (MSE up to $10^6$ in simulations); with few
  domains the ML information also reduced the median MSE bias (e.g. +3.6% vs +6.5%, $m = 36$). The ML
  MSE therefore no longer equals `sae::mseSFH`; REML is unchanged.
* `eblup_stfh()` is much faster: REML scores and Fisher information are computed from the block
  structure of $V$ (domain blocks plus a rank-$D$ spatial part) without forming $M \times M$ matrices
  ($O(T D^3 + M(D+p)^2)$ instead of $O(M^3)$ per iteration; identical results). Fitting the
  $40 \times 11$ `mys_panel` went from about 7.7 s to 0.4 s, and the parametric bootstrap (which
  re-estimates the model in every replicate) from about 2.4 s to 0.4 s per replicate.
* `hb_area(self_benchmark = TRUE)` now benchmarks every posterior draw by default
  (`benchmark_method = "posterior"`; difference / logit / ratio adjustment by family), so the constraint
  is exact for every likelihood and SD / credible intervals refer to the benchmarked estimator. Benchmark
  targets are computed on the parameter scale (proportions for binomial, rates for counts with exposure);
  before, the WFQ targets used counts while the estimates were proportions / rates. `benchmark_method =
  "wfq"` keeps the previous augmented-model approach.
* `eblup_fh(transform = "log", self_benchmark = TRUE)`: ratio benchmarking of the back-transformed
  estimates on the original scale (exact), instead of a WFQ constraint that only held on the log scale.
* `benchmark_sae()`: the group MSE used by the hierarchical "optimal" method now includes the covariances
  between area prediction errors (full MSE matrix for FH, SFH, STFH and INLA area-level HB fits), and the
  output has a new `mse_benchmarked` column, the first-order MSE of the benchmarked estimates
  ($mse_i - 2 k_i (Mw)_i + adj_i^2$; within about 1% of the empirical MSE in a simulation).
* `diagnose()`: the Brown et al. (2001) goodness-of-fit statistic $W$ now uses an upper-tail
  $\chi^2$ decision (as `emdi::compare.fh`); the former two-sided interval flagged every correct
  FH fit in a simulation, because $W$ is conservative. The bias regression of the direct on the
  model estimates is now descriptive (`brown_test$bias_regression`: `alpha`, `beta`, standard
  errors, $R^2$) instead of an F test of $\alpha = 0, \beta = 1$, which rejected 93% of correct
  FH fits (the slope exceeds 1 for shrinkage estimators); it no longer affects `status`.
  `simulation_metrics$rrmse` is the relative RMSE across domains (`mean_rrmse` duplicated
  `mean_abs_rb`). DIC/WAIC/marginal likelihood are now also read for `hb_unit()` and
  `hb_twofold()` fits, and CPO/PIT are only used when they are per domain.
* `spatial_test()`: the FH vs spatial FH comparison uses ML fits, so the likelihood ratio and
  AIC/BIC are computed from ML log-likelihoods at ML estimates (before: at REML estimates). The
  recommendation now rests on a single test at level `alpha` (`decision_test`: the likelihood
  ratio test when `vardir` is given, otherwise the LM-error test); requiring any of three tests
  recommended a spatial model for 9% of non-spatial data sets at `alpha = 0.05`. The Moran's I,
  LM-error and LM-lag statistics are unchanged and agree with `spdep`.
* `diagnose()`: standardised residuals are $(y - \hat\theta)/\sqrt{vardir - MSE}$ (the variance of
  $y - \hat\theta$) instead of $/\sqrt{vardir}$, which shrank them (sd 0.76 under a correct FH
  model) and hid outliers in the residual and Q-Q plots.
* Reliability tiers (`add_reliability_flags()`, `export_sae()`, `map_sae()`, `compare_sae()`,
  `autoplot(type = "rse")`) default to the BPS-Statistics Indonesia criteria: RSE <= 25% reliable,
  25% < RSE <= 50% use with caution, RSE > 50% unreliable (was 20%/30%, documented as the BPS
  standard). `map_sae(type = "reliability")` now uses its `thresholds` argument (they were ignored).
* `benchmark_sae(type = "total")` in a hierarchical benchmark without group targets: the group
  aggregates are weighted totals (they were weighted means, on the wrong scale for stage 1).

### Breaking changes and clean-up
* `mys_proxmat` is now the row-standardised 4-nearest-neighbour graph of the 42 areas
  (symmetric, 4 to 8 neighbours per area, built from the `sim_area` centroids by
  `data-raw/mys_proxmat.R`). It was a dense distance-decay matrix with almost equal weights
  (largest weight per row 0.03-0.04, against 1/41 = 0.024 for equal weights), which the CAR models
  of `hb_area()` read as a complete graph. Results that use `mys_proxmat` change. On `mys` the
  SAR correlation of `eblup_sfh()` stays at its bound -1 with either matrix: the residuals of
  `mys` have no positive spatial autocorrelation. `sim_area` and `sim_panel` are unchanged
  (simulated with the old matrix).
* `benchmark()`, a plain alias of `benchmark_sae()`, is removed; use `benchmark_sae()`.
* The leftover handling of `df_ebp`/`ebp` result tables (from the removed `ebp_area()`) is removed
  from the methods and utilities; results are in `df_hb` (HB) or `df_eblup` (EBLUP).
* `diagnose()` returns `brown_test$bias_regression` instead of `brown_test$wald_test` and
  `simulation_metrics$rrmse` instead of `mean_rrmse` (see above).

### Other fixes in the utilities
* `map_sae()` works with any name of the `sf` geometry column (it assumed `"geometry"`), labels
  HB fits correctly (they were labelled "EBP-<family>"), and stops with a clear message for panel or
  sub-area tables with repeated domains instead of mis-merging them.
* `compare_sae()` and `autoplot()` plots label panel and sub-area rows uniquely (repeated domains
  made some plots fail or overlap); `autoplot()` of an unnamed list of models names them.
* `print()`/`summary()` show the estimation method of `eblup_fh()` and `eblup_twofold()` (REML or
  ML; was "eblup") and describe `eblup_stfh(model = "S")` correctly.
* `export_sae()` lists both variance components of two-fold models in the model summary.

### Performance
* `benchmark_sae()`: the MSE covariance between areas is applied as a matrix-free operator
  (FH: $O(n p^2)$, no $n \times n$ matrices). FH with $m = 2000$: 9.4 s / 490 MB -> 0.002 s / 1.7 MB;
  $m = 20000$ now runs in about 0.01 s.
* `eblup_bhf()`: the REML/ML fit and the bootstrap MSE now run in C++ (`src/eblup_bhf_core.cpp`).
  The likelihood is profiled to $\gamma = \sigma_u^2/\sigma_e^2 \ge 0$ and evaluated from within-domain
  cross-products and domain means ($O(Dp^2)$ per evaluation, independent of the number of units); the
  maximum is found from a log-grid scan plus Brent refinement, so $\sigma_u^2 = 0$ is handled exactly.
  Estimates agree with `lme4`/`sae` up to their optimiser tolerance (deviance never worse), and the
  bootstrap uses the same random draws as before, so MSEs agree with the previous `lme4`-based
  implementation to about 1e-5 for a given seed. Bootstrap: cornsoybean $B = 200$ 0.64 s -> 0.005 s;
  $D = 300$, $n = 6000$, $B = 100$ 0.96 s -> 0.035 s. `n_threads` now runs replicates in parallel
  (results do not depend on it); for small models R's sequential RNG dominates the run time.
* `eblup_bhf()` output: `fit$lme` (the `lmerMod` object) is replaced by `fit$residuals` (unit-level
  conditional residuals, used by `residuals()`); `estcoef` standard errors are $(X'\hat V^{-1}X)^{-1}$
  as in `lme4`. Rank-deficient designs now give an explicit error.
* `eblup_sfh()` (fit, analytical MSE, kriging and the bootstrap refits) no longer forms any
  $D \times D$ inverse or $D^3$ product. With $A = (I - \rho W')(I - \rho W)$ and
  $M = A + \sigma^2 E'\Psi^{-1}E$ (both sparse), $V^{-1}$, $GV^{-1}$, $g_1 = \sigma^2 \mathrm{diag}(M^{-1})$,
  the scores/information, $g_3$, $g_4$, the ML bias and the kriging terms follow from Cholesky
  solves with $A$ and $M$ (`src/sar_engine.h`; Eigen sparse $LDL'$ with AMD ordering, dense
  factorisations when $W$ itself is dense, e.g. distance weights). Results are unchanged (relative
  differences below 2e-10). Analytical: $D = 400$ 2.2 s -> 0.10 s, $D = 1600$ 80 s -> 2.5 s,
  $D = 2500$ about 6 s; bootstrap ($D = 400$, $B = 20$): 16 s -> 1.3 s. `RcppEigen` is added to
  `LinkingTo`.
* `eblup_stfh()`: the final step (Fisher information at the estimates, SEs, $g_1$-$g_4$ of the analytical
  MSE) no longer forms $M \times M$ matrices ($M = DT$); it uses the same structured representation of $V$
  as the iterations ($A^{-1}$ block diagonal, $P = A^{-1} - U W U'$ with $U = [A^{-1}Z, V^{-1}X]$ whose
  $A^{-1}Z$ part has one non-zero per row). The scoring step also exploits the diagonal blocks of $U'Z$ and
  $U'EU$ and the sparsity of $C = 2\rho W'W - W - W'$. Results are unchanged (relative differences below
  3e-11). $T = 5$: $D = 250$ 28 s -> 1.6 s, $D = 500$ about 4 min -> 10 s, $D = 1000$ 75 s.
* `create_weights()`: distances for `knn`, `distance` and `idw` are computed in blocks of rows (at most
  about 1M distances at a time) instead of a full $D \times D$ distance matrix, and the weight matrix is
  filled once from the neighbour pairs (no intermediate 0/1 matrix and no normalisation copies; `idw`
  styles are applied per block). The result is the same matrix as before (identical, or within 1e-16
  for `idw`). Peak memory at $D = 10000$: queen 2.3 GB -> 0.8 GB, knn 3.1 GB -> 0.8 GB, distance
  3.8 GB -> 0.8 GB, i.e. essentially the returned matrix itself (764 MB).

### Bug fixes
* `create_weights()` with a coordinate matrix whose columns are named `X`/`Y` (e.g. from
  `sf::st_coordinates()`) failed with "undefined columns selected"; coordinate column names are now
  matched case-insensitively.

### Major New Features & Enhancements
* **Two-Fold Subarea EBLUP (`eblup_twofold`)**:
  * Implements nested two-fold sub-area level models (Torabi & Rao, 2014) with area effects $v_i$ and sub-area effects $u_{ij}$.
  * Ultra-fast Fisher-scoring algorithms in C++ for REML and ML estimation.
  * Analytical Prasad-Rao MSE ($g_1 + g_2 + 2 g_3$ for REML, plus the bias correction for ML) with exact matrix derivation and multi-threaded parametric bootstrap MSE.
  * Automatic prediction for non-sampled subareas using area-level random effects and regression synthetic estimation.
* **Hierarchical Bayes Two-Fold Subarea Models (`hb_twofold`)**:
  * Bayesian nested two-fold sub-area models via `INLA`.
  * Supports Gaussian, Binomial (logit link), and Poisson (log rate) likelihoods.
  * Area-level random effects support IID or spatial structures (BYM2, Besag).
  * Simultaneous posterior estimation for subareas (`df_hb`) and weighted area aggregates (`df_area`) with full Monte Carlo uncertainty propagation.
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
  * `diagnose()`: Brown et al. (2001) goodness-of-fit and bias diagnostics, residual analysis, spatial autocorrelation tests (Moran's I), and official statistical reliability classification (BPS standards).
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
