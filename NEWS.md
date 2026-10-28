# fastsae 1.0.0

fastsae 1.0.0 adds hierarchical Bayes (HB) models, the two-fold sub-area model and
post-estimation tools to the frequentist EBLUP functions of 0.1.0, and corrects several
estimators of 0.1.0. The changes below are relative to fastsae 0.1.0 on CRAN.

## Breaking changes for users of 0.1.0

* All EBLUP functions choose the MSE with `mse_type`, which includes `"none"` (no MSE;
  `mse` and `rse` are `NA`): `eblup_fh()` (`"analytical"`, `"none"`), `eblup_sfh()`
  (`"analytical"`, `"pbmse"`, `"npbmse"`, `"none"`), `eblup_stfh()` and `eblup_twofold()`
  (`"analytical"`, `"bootstrap"`, `"none"`) and `eblup_bhf()` (`"bootstrap"`, `"none"`). The
  chosen value is returned in `mse_type`.
  * `eblup_sfh(mse_method = )` and `compute_mse = ` of `eblup_stfh()` and `eblup_bhf()` are
    deprecated and still work with a warning (`compute_mse = FALSE` is `mse_type = "none"`).
  * The MSE is now computed by default: `eblup_stfh()` gives the analytical MSE and
    `eblup_bhf()` the bootstrap MSE (both gave none before). Use `mse_type = "none"` for the
    previous output.
  * New arguments were inserted before existing ones (`transform` and `mse_type` in
    `eblup_fh()`, `mse_type` in `eblup_sfh()`), so calls that pass the later arguments by
    position change meaning; calls by name are unaffected.
* `eblup_bhf()`: `fit$lme` (the `lmerMod` object) is replaced by `fit$residuals` (unit-level
  conditional residuals, used by `residuals()`), since the model is now fitted in C++ (see
  below); the standard errors in `estcoef` are $(X'\hat V^{-1}X)^{-1}$ as in 'lme4'.
* `goodness` of `eblup_fh()`, `eblup_sfh()` and `eblup_stfh()` follows one convention:
  `loglikelihood` is the full Gaussian log-likelihood at the estimates, and `AIC`/`BIC` are
  computed from it for ML and REML fits, as in 'sae' and Marhuenda, Morales and Pardo (2014).
  In `eblup_stfh()` the element `loglike` is renamed `loglikelihood`.
* Convergence problems are reported as a warning that can be caught. Before, `eblup_fh()`,
  `eblup_sfh()` and `eblup_stfh()` only printed a message and `eblup_bhf()` said nothing.
* Inputs are checked before fitting, so some calls that ran before now stop with a clear
  message: `vardir` must be finite and positive in every sampled area; linearly dependent
  covariates stop the fit and name the aliased column (the compiled code used an approximate
  solution of a singular system); there must be more sampled areas than coefficients (all
  MSEs were `NA`); `W` of `eblup_sfh()` and `eblup_stfh()` must be non-negative (a zero row or
  rows that do not sum to 1 give a warning, since the SAR model with `rho` in (-1, 1) assumes
  a row-standardised `W`, as 'sae' requires); a missing population mean in `eblup_bhf()`
  names the domain. Duplicated `domain` values give a warning.
* A fixed `seed` of `eblup_sfh()`, `eblup_stfh()` and `eblup_bhf()` no longer changes the
  random number stream of the session; the bootstrap MSEs are unchanged.
* `mys_proxmat` is now the row-standardised 4-nearest-neighbour graph of the 42 areas
  (symmetric, 4 to 8 neighbours per area, built by `data-raw/mys_proxmat.R`). It was a dense
  distance-decay matrix with almost equal weights, which a CAR model reads as a complete
  graph, and its documentation wrongly said it was based on contiguous boundaries. Results
  that use `mys_proxmat` change.
* Several estimators of 0.1.0 were corrected, which changes their results (see
  "Changes to the models of 0.1.0"); the most visible are the spatio-temporal model on
  unsorted panels such as `mys_panel`, its analytical and bootstrap MSE, and the spatial model
  with unsampled domains.

## New models

### Hierarchical Bayes models

* `hb_area()`: area-level HB models with Gaussian, beta, binomial, Poisson, negative binomial
  and gamma likelihoods; iid, Besag (ICAR), BYM, BYM2, Leroux (`"leroux"`/`"generic1"`), SAR
  (`"sar"`, the Bayesian counterpart of `eblup_sfh()`) and spatial lag (`"slm"`) area
  effects; RW1, RW2, AR(1) and iid time effects; and the four space-time interactions of
  Knorr-Held (2000), area-specific trends and separable fields. Areas without a direct
  estimate are predicted. With `vardir`, the dispersion of each area follows from its sampling
  variance:
  * beta: the precision $\phi_d = \hat\mu_d(1 - \hat\mu_d)/V_d - 1$ at the posterior mean of a
    first fit (two stages). At the direct estimate, the model variance was
    $V_d \mu_d(1 - \mu_d)/(y_d(1 - y_d))$, not $V_d$, and the 95% intervals covered 81-84% of
    true proportions near 0.08, against 90-94% (`benchmarks/validation/beta_vardir_phi.R`);
  * negative binomial: the size $n = \hat\mu^2/(V - \hat\mu)$ at the posterior mean of a first
    fit (a near-Poisson size when $V \le \hat\mu$); the 95% intervals covered 90-95%
    (`benchmarks/validation/dispersion_vardir.R`, `nb_size_prior.R`);
  * gamma: the shape $y^2/V$ of each area.

  Without `vardir`, the gamma, negative binomial and Gaussian dispersions are only weakly
  identified with one observation per area, and `hb_area()` warns. `inla_engine(link = )`
  sets the link (`"probit"`, `"cloglog"` and `"loglog"` for binomial and beta).
* `hb_unit()`: unit-level HB models (Gaussian, Bernoulli or Poisson units; iid, Besag or BYM2
  domain effects), with finite population means when `popsize_var` is given: the non-sampled
  part is predicted at the non-sampled covariate mean
  $\bar x_{d,r} = (N_d \bar X_d - n_d \bar x_{d,s})/(N_d - n_d)$, and its posterior SD includes
  $\sigma_e^2/(N_d - n_d)$. For binomial and Poisson units the domain means are approximate
  (the non-sampled part is predicted at the mean covariates; -8% and -16% bias in a
  simulation) and `hb_unit()` warns; unbiased estimation needs unit-level population
  covariates (Molina and Rao, 2010). The responses are checked (0/1 for binomial,
  non-negative integers for Poisson).
* `hb_twofold()`: the two-fold sub-area model with an iid, Besag or BYM2 area effect and an iid
  sub-area effect (Gaussian, binomial, Poisson), with weighted area aggregates in `df_area`
  (Gaussian: exact INLA linear combinations; binomial and Poisson: draws from the Gaussian
  approximation of the latent field at INLA's configurations, `inla_engine(draws = 2000)` by
  default). Missing or negative weights give a warning.
* Binomial models of `hb_area()` and `hb_twofold()` take survey-weighted proportions with
  effective sample sizes $n^* = n/\mathrm{deff}$ in `trials` (the binomial/logit-normal model
  of Chen, Wakefield and Lumley, 2014, and Franco and Bell, 2015): proportions become the
  effective numbers of cases $n^* \hat p$, which are not rounded, and direct estimates of 0
  stay sampling zeros. Non-integer cases or trials use INLA's `"xbinomial"` likelihood (equal
  to the binomial one on whole numbers, which keep `"binomial"`); `jax_engine()` accepts them
  too. A response outside $[0, n]$ is an error.
* Estimation engines: `inla_engine()` (default) and `jax_engine()` hold all settings of an
  engine, so the model functions only describe the model; `engine = "inla"` or `"jax"` take the
  defaults. `inla_engine(strategy, link, scale_model, prior_*, draws, seed, ...)` passes `...`
  to `INLA::inla()` and merges `control.*` lists with the package defaults. `jax_engine()` fits
  the same models by NumPyro NUTS on JAX (CPU, NVIDIA GPU, or the Apple GPU with
  `device = "metal"`) through 'reticulate', with `prior_beta = "horseshoe"` and
  `robust = TRUE` (Student-t random effects) as options. MCMC fits report the rank-normalized
  split R-hat and bulk/tail ESS of every parameter (Vehtari et al., 2021; identical to
  'posterior'), divergences and tree-depth hits, warn for R-hat >= 1.05, ESS below 100 per
  chain and divergent transitions, and `autoplot(type = "trace")` / `"rank"` plot the chains.
  The engine settings are returned in `engine_spec`.
* `setup_numpyro_env()` installs the Python side, and `setup_numpyro_env(check_only = TRUE)`
  checks it without installing or searching the system for Python (about one second).

### Default priors and INLA settings of the HB models

* The random-effect standard deviations get a half-normal prior with the scale of the
  response, sd(y), for Gaussian models and half-normal(1) on the link scale otherwise, instead
  of INLA's PC prior $P(\sigma > 1) = 0.01$, which depends on the unit of the response (on
  `mys` with the response in tenths it moved the estimates by up to 4.1 posterior SD). In
  simulations the 95% intervals covered 93-97% with the half-normal prior (EBLUP: 93-96%) and
  85-95% with PC priors scaled to sd(y) (`benchmarks/validation/default_sd_prior*.R`).
* The Gaussian residual SD without `vardir` gets the half-normal(sd y) prior (INLA's default
  put it near zero, with 0-4% coverage); the negative binomial size and the gamma precision a
  half-normal(10) prior (INLA's defaults put them at 110-180 and 76-84 whatever the truth);
  the beta precision without `vardir` or `trials` INLA's log-gamma(1, 0.1); Leroux's $\lambda$
  and the SAR/SLM $\rho$ a uniform prior. `jax_engine(priors = "inla")`, the default, uses the
  same priors, so the two engines agree within Monte Carlo error and the INLA approximation;
  `priors = "weakly_informative"` keeps the earlier NumPyro priors.
* Gaussian models, and the Leroux and SLM effects for any family, are fitted in INLA's
  `"classic"` mode, the other models in its `"compact"` mode; Gaussian models with a type I-IV
  interaction and three hyperparameters integrate over a grid. In comparisons with exact
  posteriors and long Stan runs (`benchmarks/validation`), the compact mode overstated the
  posterior SDs of Gaussian models by up to 18%, while the classic mode was within 0.07
  posterior SD (means) and 0.91-1.03 times (SDs); for the other families both modes matched
  Stan and compact was much faster. `inla_engine(inla.mode = , control.inla = )` overrides this.
* Space-time interactions with a random walk in time carry all the constraints of Knorr-Held
  (2000), so domains without direct estimates are identified; the `"separable"` field with a
  random walk and the `"domain-specific"` RW2 stop with an error.
* CAR models (`"besag"`, `"bym"`, `"bym2"`, `"leroux"`) use only the pattern `W > 0`; a dense
  `W` such as distance-decay weights gives a warning.
* Variance components are reported as posterior means of the variance, $E[1/\tau \mid y]$.
* Posterior draws of the INLA engine are reproducible: `inla_engine(seed = 1L)` is the default
  and the seed is also set on the R side; the session's random number stream is restored.
* The new vignette "Default Priors of the Hierarchical Bayes Models" documents the priors of
  both engines and the simulations behind them.

### Frequentist models and options

* `eblup_twofold()`: the two-fold sub-area model of Torabi and Rao (2014) by REML or ML, with
  the analytical MSE $g_1 + g_2 + 2g_3$ (minus the ML bias term) or a parametric bootstrap,
  and the area means $\theta_i = \sum_j w_{ij}\theta_{ij}$ in `df_area` with an MSE that
  includes the covariances of the sub-area errors (`weight`, `compute_area`). With one
  sub-area per area, Fisher scoring uses the pseudo-inverse of the singular information; a
  negative second-order MSE is set to `NA` with a warning.
* `eblup_fh(transform = "log")`: Fay-Herriot model on the log scale, back-transformed with the
  conditional-mean correction $\exp(\tilde\mu + g_1/2)$ and a delta-method MSE (not the full
  Slud-Maiti second-order correction); unsampled domains get $\exp(x'\hat\beta + \hat\sigma^2_u/2)$.
* Self-benchmarking (`self_benchmark`, `benchmark_weight`, `benchmark_target`,
  `benchmark_group`) in the EBLUP functions (augmented model of Wang, Fuller and Qu, 2008) and
  in `hb_area()` (every posterior draw is benchmarked by default, `benchmark_method =
  "posterior"`, so the constraint is exact and the intervals refer to the benchmarked
  estimator; `"wfq"` uses the augmented model). With `transform = "log"`, `eblup_fh()`
  ratio-benchmarks on the original scale.

## New post-estimation tools and utilities

* `benchmark_sae()`: ratio, difference, MSE-weighted ("optimal") and logit benchmarking to
  one or two levels of targets (Datta et al., 2011; Bell, Datta and Ghosh, 2013; two-stage
  benchmarking in the spirit of Ghosh and Steorts, 2013). The optimal adjustment is
  proportional to $w_i\,\mathrm{MSE}_i$ and uses the full MSE matrix of FH, SFH, STFH and
  area-level HB fits; `mse_benchmarked` is a first-order MSE of the benchmarked estimates.
  Missing or non-finite targets and missing weights stop with an error.
* `diagnose()`: the goodness-of-fit statistic of Brown et al. (2001) with an upper-tail
  $\chi^2$ decision (as `emdi::compare.fh`), a descriptive bias regression, standardised
  residuals $(y - \hat\theta)/\sqrt{vardir - MSE}$, Moran's I of the residuals, DIC, WAIC and
  CPO/PIT for HB fits, and simulation metrics when the truth is known.
* `compare_sae()`: agreement of two models by the mean difference and 95% limits of agreement
  (Bland and Altman, 1986), Lin's (1989) concordance correlation, differences relative to the
  RMSE, rank changes and the overlap of the highest and lowest domains (`top`). Pearson's
  correlation is not reported, since it measures association, not agreement. Between an HB and
  an EBLUP fit it reports which has the lower reported uncertainty and notes that a posterior
  variance and an MSE are not the same quantity. `autoplot(type = "difference")` is a
  Bland-Altman plot.
* `step_sae()`: stepwise covariate selection for every model, refitting with the same
  settings: AIC, BIC (default) or KIC from the full log-likelihood for EBLUP models
  (Marhuenda, Morales and Pardo, 2014; models with a `method` are searched by ML and refitted
  with the original method; on `mys` the path equals `emdi::step()`), and LCPO, WAIC or DIC
  for HB models (`benchmarks/validation/step_sae.R`).
* `add_reliability_flags()`, `export_sae()`, `map_sae()`: reliability tiers by RSE (default
  cut-offs of BPS-Statistics Indonesia: RSE <= 25% reliable, 25-50% use with caution, > 50%
  unreliable), export of estimates and model summaries to Excel or CSV, and choropleth maps of
  estimates, RSE, reliability, differences and benchmarking (any `sf` geometry column).
* `create_weights()` builds spatial weights (queen, rook, k-nearest neighbours, distance bands,
  inverse distance) from `sf` objects, coordinates or neighbour lists, computing distances in
  blocks of rows (peak memory at 10,000 domains about the size of the returned matrix);
  `spatial_test()` gives Moran's I, LM-error and LM-lag tests (as 'spdep') and a recommendation
  from a single test at level `alpha` (likelihood ratio test of FH against spatial FH, both
  fitted by ML, when `vardir` is given); `gvf_smooth()` smooths sampling variances by a
  generalized variance function.
* `sim_area_data()` and `sim_series_data()` simulate area-level and panel data for all six
  families with spatial (and AR(1)) effects, returning the spatial weights; the random number
  stream of the session is restored after `seed`. New datasets `sim_area` and `sim_panel`.
* `eblup_bhf()` returns `goodness` (full log-likelihood, AIC and BIC); with `method = "ML"`
  they equal those of `lme4::lmer(REML = FALSE)`.
* `print()`/`summary()` show the estimation method and engine, and `autoplot()` labels panel
  and sub-area rows uniquely.

## Changes to the models of 0.1.0

### `eblup_stfh()`

* The panel is sorted internally into the domain-major order the model requires and results
  are returned in the caller's row order. Before, unsorted input was silently mis-specified,
  including the time-major `mys_panel` of the examples (where `sae::eblupSTFH()`, which
  assumes domain-major data, gives the same wrong fit).
* Unsampled domains and missing periods: cells whose response is `NA` stay in the data and are
  predicted. The model is defined on the full domain x time grid, the variance components are
  estimated by REML from the sampled cells, and every cell gets the best linear unbiased
  predictor given the sampled data (Saei and Chambers, 2005; Benedetti et al., 2024), with the
  analytical MSE $g_1 + g_2 + g_3$ for the predicted cells and a parametric bootstrap that
  covers every cell. `df_eblup` gains a `sampled` column; `W` must include the unsampled
  domains. Results on complete panels are unchanged.
* Analytical MSE: $g_3$ used the wrong matrix orientation (wrong whenever `vardir` is
  heterogeneous); it now uses the exact derivative of the EBLUP and the REML information at the
  final estimates, and the second-derivative term $g_4$ is included:
  `mse = g1 + g2 + 2 g3 - g4` (relative bias -1.7% in a simulation; +2.4% without $g_4$).
* Bootstrap MSE: the variance components are re-estimated in every replicate (they were held
  fixed, so $g_3$ was missed), the SAR effects are generated with the lower Cholesky factor,
  and non-converged replicates are replaced.
* When the classical Fisher scoring fails or converges to a negative variance, the fit switches
  to constrained REML (convergence in a simulation went from 63% to 99%; interior fits are
  identical to `sae::eblupSTFH()`); `convergence = FALSE` is reported when the fit ends at
  invalid variance components. Correlations on the boundary are treated as fixed in the MSE.
* The model is fitted by REML only, as in Marhuenda, Molina and Morales (2013); `tvalue` is a
  z-statistic.

### `eblup_sfh()`

* With unsampled domains, the SAR process is defined on all domains of `W`, so the sampled
  domains are fitted with the marginal covariance (before: the SAR built from `W[s, s]`),
  consistent with the kriging of the unsampled domains and with the bootstrap. Their analytical
  MSE includes $g_3$ (exact derivative of the kriging predictor) and the second-order bias of
  the plug-in $g_1$ (relative bias -3.8%, was -7.5%).
* When the classical Fisher scoring fails or ends with $\hat\sigma^2 < 0$, the (RE)ML
  likelihood is maximised under the constraints, also inside the bootstrap. A boundary
  $\hat\rho = 0.999$ is no longer mapped to 1, and $-0.999$ is reported as $-1$ only when
  $I + W$ is non-singular.
* A correlation on the boundary ($|\hat\rho| \ge 0.999$) is treated as fixed in the analytical
  MSE, where the Prasad-Rao expansion is not valid; for ML fits, $g_3$ and the bias term use the
  ML information (the REML information made the MSE explode for $\hat\sigma^2 \approx 0$ and
  $\hat\rho \approx 1$). These are the only cases where the MSE differs from `sae::mseSFH()`.
* Bootstrap replicates are no longer rejected when estimates exceed `1e6`, and `npbmse` no
  longer fails when $\hat\sigma^2_u = 0$.

### `eblup_fh()` and `eblup_sfh()`

* A random effect variance estimated as 0 (truncated, as in 'sae') is reported with a message:
  the EBLUPs are then synthetic estimates and the MSE estimates can understate their error
  (Li and Lahiri, 2010).

### `eblup_bhf()`

* The REML/ML fit and the bootstrap MSE run in C++, with the likelihood profiled to
  $\gamma = \sigma_u^2/\sigma_e^2 \ge 0$ and evaluated from within-domain cross-products, so
  $\sigma_u^2 = 0$ is handled exactly. Estimates agree with 'lme4' and 'sae' up to their
  optimiser tolerance; the bootstrap uses the same random draws as before, `n_threads` runs
  replicates in parallel (results do not depend on it), and it uses the same population means
  (`popnmean_xpop`) as the estimate. Rank-deficient designs give an explicit error.

## Performance

* `eblup_sfh()` (fit, analytical MSE, kriging and the bootstrap refits) forms no $D \times D$
  inverse: all terms follow from sparse Cholesky solves with $(I - \rho W')(I - \rho W)$ and
  $M = A + \sigma^2 E'\Psi^{-1}E$ (`RcppEigen` added to `LinkingTo`). Analytical MSE:
  $D = 1600$ 80 s -> 2.5 s; results unchanged (relative differences below 2e-10).
* `eblup_stfh()` computes the REML scores, the Fisher information and the analytical MSE from
  the block structure of $V$, without $M \times M$ matrices ($M = DT$): `mys_panel` 7.7 s ->
  0.4 s; $T = 5$, $D = 500$ about 4 min -> 10 s; results unchanged.
* `eblup_bhf()` bootstrap: `cornsoybean` with $B = 200$ 0.64 s -> 0.005 s.

## Documentation

* New vignettes on the HB area-level models, default priors, unit-level and two-fold models,
  model diagnostics, the `hb_area()` engines, finite population validation and benchmarks;
  every exported function and S3 method has examples.
* The website moved to <https://ridsonap.is-a.dev/fastsae/>.

# fastsae 0.1.0

* Initial CRAN release.
* High-performance area-level Fay-Herriot (FH), spatial FH (SFH) and spatio-temporal FH
  (STFH) models.
* Unit-level Battese-Harter-Fuller (BHF) models via Henderson III and REML.
