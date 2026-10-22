# Default Priors of the Hierarchical Bayes Models

[`hb_area()`](https://ridsonap.github.io/fastsae/reference/hb_area.md),
[`hb_unit()`](https://ridsonap.github.io/fastsae/reference/hb_unit.md)
and
[`hb_twofold()`](https://ridsonap.github.io/fastsae/reference/hb_twofold.md)
fit their models with INLA
([`inla_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md),
the default) or with NUTS in NumPyro
([`jax_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)).
Both engines use the same default priors, so a model gives the same
estimates whichever engine fits it. This article lists those priors,
explains why they were chosen and summarises the simulations behind
them.

With few areas, or with one observation per area, the prior of a
variance or dispersion parameter is not a technicality: the data hardly
inform it, and its prior decides how much the estimates are shrunk and
how wide the intervals are.

## Priors of the model parameters

| Parameter | Default ([`inla_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md), `jax_engine(priors = "inla")`) | `jax_engine(priors = "weakly_informative")` |
|----|----|----|
| Random-effect SD $`\sigma`$ (area, sub-area, time, area-by-time) | Half-normal$`(s)`$: $`s = \mathrm{sd}(y)`$ for Gaussian models, $`s = 1`$ on the link scale otherwise | Half-normal$`(1)`$ on the standardised response or the link scale |
| BYM2 mixing parameter $`\phi`$ | PC prior with $`P(\phi < 0.5) = 0.5`$ | Uniform$`(0, 1)`$ |
| AR(1) correlation $`\rho`$ | $`N(0, 1/0.15)`$ on $`\log\{(1 + \rho)/(1 - \rho)\}`$ | Uniform$`(-1, 1)`$ |
| Leroux $`\lambda`$, SAR $`\rho`$ | Uniform | Uniform |
| Coefficients $`\beta`$ | $`N(0, 1000)`$; intercept flat (INLA) or nearly flat (JAX) | $`N(0, 10^2)`$ on the standardised scale ($`N(0, 2.5^2)`$ for logit links) |

Here $`\mathrm{sd}(y)`$ is the standard deviation of the observed direct
estimates (of the unit responses for
[`hb_unit()`](https://ridsonap.github.io/fastsae/reference/hb_unit.md)).
[`inla_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)
writes the half-normal prior as an `expression:` prior on the log
precision $`\theta = -2 \log \sigma`$, so that INLA and NumPyro use the
same density.

[`jax_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)
has two priors that INLA does not have. With `robust = TRUE` the random
effects are Student-$`t`$ with degrees of freedom $`\nu \sim`$
Uniform$`(2.5, 30)`$. With `prior_beta = "horseshoe"` the coefficients
get a regularised horseshoe prior (Piironen and Vehtari 2017): local
scales $`\lambda_j \sim C^+(1)`$, global scale $`\tau \sim C^+(1)`$ and
slab $`c^2 \sim`$ Inv-Gamma$`(2, 8)`$.

## Priors of the likelihood

| Parameter | When | Prior |
|----|----|----|
| Sampling variance (Gaussian), gamma precision | `vardir` given | Fixed per area: $`1/V_d`$; $`y_d^2 / V_d`$ |
| Beta precision | `vardir` given | Fixed per area at $`\hat\mu_d(1 - \hat\mu_d)/V_d - 1`$, with $`\hat\mu_d`$ from a first fit |
| Beta precision | `trials` given (no `vardir`) | Fixed per area at $`n_d - 1`$ |
| Negative binomial size | `vardir` given | Fixed per area at $`\hat\mu_d^2 / (V_d - \hat\mu_d)`$, with $`\hat\mu_d`$ from a first fit |
| Gamma precision | no `vardir` | Half-normal$`(10)`$, both engines |
| Negative binomial size | no `vardir` | Half-normal$`(10)`$, both engines |
| Gaussian residual SD | no `vardir` (and the unit errors of [`hb_unit()`](https://ridsonap.github.io/fastsae/reference/hb_unit.md)) | Half-normal$`(\mathrm{sd}(y))`$, both engines |
| Beta precision $`\phi`$ | no `vardir`, no `trials` | Log-gamma$`(1, 0.1)`$ on $`\log\phi`$, i.e. $`\phi \sim`$ Exponential$`(0.1)`$, both engines |

## Why half-normal standard deviations

INLA’s penalised complexity (PC) prior on an SD,
$`P(\sigma > U) = \alpha`$, is a good default when $`U`$ is chosen for
the problem (Simpson et al. 2017). The common choice $`U = 1`$,
$`\alpha = 0.01`$ is not: it is tied to the unit of the response. On the
`mys` data, the REML estimate of the area SD is 1.62 and
$`\mathrm{sd}(y) = 2.21`$, but PC$`(1, 0.01)`$ gives
$`P(\sigma > 1.62) = 0.0006`$, and the prior shrinks the area effects
far more than the data suggest. In a Fay-Herriot simulation (30 areas)
the 95% intervals under PC$`(1, 0.01)`$ covered 94%, 85% and 38% of the
true means with the same data in three units (`y` times 0.1, 1 and 10).

Scaling the prior to the response removes the dependence on the unit.
Among the scaled priors, the half-normal covered as the EBLUP intervals,
while the PC priors, which put much mass near $`\sigma = 0`$, still
shrank the area SD (`benchmarks/validation/default_sd_prior.R`; 150
replicates, known `vardir`):

| Areas | Area SD | EBLUP | Half-normal$`(\mathrm{sd}\,y)`$ | PC, $`\alpha = 0.01`$ | PC, $`\alpha = 0.05`$ |
|----|----|----|----|----|----|
| 15 | 0.4 | 0.960 | 0.966 | 0.917 | 0.942 |
| 30 | 0.4 | 0.940 | 0.958 | 0.905 | 0.929 |
| 60 | 0.4 | 0.927 | 0.943 | 0.898 | 0.918 |
| 15 | 0.8 | 0.930 | 0.932 | 0.853 | 0.890 |
| 30 | 0.8 | 0.934 | 0.929 | 0.888 | 0.909 |
| 60 | 0.8 | 0.942 | 0.939 | 0.919 | 0.928 |

For models with a link function the random effects are on the logit or
log scale, where an SD of 1 is already large, so the scale is 1. Against
PC$`(1, 0.01)`$ and half-normal$`(2.5)`$ (the default of tipsae),
half-normal$`(1)`$ covered at least as well, and it estimated the area
SD with less bias (`default_sd_prior_glm.R`,
`default_sd_prior_glm_more.R`; mean coverage over the scenarios, worst
scenario in brackets):

| Response | Half-normal$`(1)`$ | Half-normal$`(2.5)`$ | PC$`(1, 0.01)`$ |
|----|----|----|----|
| Binomial, $`n_d`$ 20-200 | 0.927-0.950 | 0.927-0.952 | 0.918-0.944 |
| Binomial, $`n_d`$ 5-20 | 0.950 (0.932) | 0.953 (0.938) | 0.919 (0.864) |
| Beta, $`n_d`$ 5-20 | 0.938 (0.925) | 0.940 (0.929) | 0.908 (0.873) |
| Poisson | 0.947 (0.944) | 0.948 (0.944) | 0.943 (0.934) |
| Binomial BYM2 | 0.947 (0.941) | 0.947 (0.942) | 0.942 (0.934) |

The half-normal prior is also the weakly informative choice of Gelman
(2006) for hierarchical SDs. It suits NUTS better than INLA’s vague
log-gamma prior on the precision, whose ridge towards $`\sigma = 0`$
slows the sampler down.

The BYM2 and AR(1) priors are INLA’s defaults (Riebler et al. 2016),
reproduced in JAX:
[`jax_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)
evaluates the PC density of $`\phi`$ from the eigenvalues of the scaled
ICAR precision, as `INLA:::inla.pc.bym.phi()` does. With these priors
the two engines agree on `mys` within 0.03 (iid) and 0.05 (BYM2)
posterior SD.

## Dispersion with one observation per area

In an area-level model each area contributes one direct estimate. For
the gamma and negative binomial likelihoods the variance of that
estimate depends on a dispersion parameter as well as on the area
effect:

- gamma: $`\mathrm{CV}^2 \approx 1/\text{shape} + \sigma^2`$;
- negative binomial:
  $`\mathrm{Var} \approx \mu + \mu^2/\text{size} + \mu^2(e^{\sigma^2} - 1)`$.

With one observation per area the data identify only the sum, so the
dispersion follows its prior. INLA’s own defaults favour little sampling
noise: log-gamma$`(1, 0.01)`$ on the gamma precision and
pc.mgamma$`(7)`$ on the size, which shrinks towards Poisson. They put
the estimate at 76-84 (gamma) and 110-180 (size) whatever the truth, and
the intervals then missed the true means when the sampling noise was
real. Both engines therefore use a half-normal$`(10)`$ prior on the
gamma precision and on the size. Coverage of 95% intervals, sampled
areas, area SD 0.3 / 0.8 (`gamma_shape_prior.R`, `nb_size_prior.R`; 30
areas, 100 replicates):

| Likelihood        | True dispersion | Half-normal$`(10)`$ | INLA’s default |
|-------------------|-----------------|---------------------|----------------|
| Gamma             | shape 5         | 0.930 / 0.925       | 0.593 / 0.522  |
| Gamma             | shape 20        | 0.959 / 0.986       | 0.862 / 0.826  |
| Gamma             | shape 100       | 0.950 / 0.999       | 0.979 / 0.993  |
| Negative binomial | size 3          | 0.928 / 0.896       | 0.816 / 0.756  |
| Negative binomial | size 10         | 0.949 / 0.965       | 0.913 / 0.915  |
| Negative binomial | size 100        | 0.939 / 0.990       | 0.965 / 0.980  |

The half-normal$`(10)`$ prior has a cost: when the sampling noise is
small (shape or size 100) its intervals are too wide and its errors
larger.
[`hb_area()`](https://ridsonap.github.io/fastsae/reference/hb_area.md)
therefore warns when a gamma or negative binomial model is fitted
without `vardir`. The sampling variances fix the dispersion of each area
and give 92-94% (gamma) and 91-95% (negative binomial) coverage
(`dispersion_vardir.R`). A Poisson model is not an alternative for
overdispersed counts: it ignores the overdispersion and covered 68-75%.

The Gaussian and beta likelihoods have the same problem without
`vardir`. For the Gaussian residual SD, INLA’s own default,
log-gamma$`(1, 5 \times 10^{-5})`$ on the precision, is disastrous here:
it put the residual SD near zero (precision about $`1.5 \times 10^4`$),
the area means collapsed onto the direct estimates and their intervals
covered almost none of the true means. The half-normal prior with the
scale of the response, the prior of the area SD, covered nominally in
every unit of the response (`likelihood_priors.R`; 30 areas, 60
replicates, coverage of 95% intervals for sampled areas, $`y`$ in its
unit / times 10):

| Area SD | Residual SD | Half-normal$`(\mathrm{sd}\,y)`$ | INLA’s default |
|---------|-------------|---------------------------------|----------------|
| 0.5     | 0.5         | 0.943 / 0.943                   | 0.038 / 0.006  |
| 0.5     | 1           | 0.949 / 0.947                   | 0.024 / 0.001  |
| 1       | 0.5         | 0.955 / 0.968                   | 0.033 / 0.003  |
| 1       | 1           | 0.953 / 0.938                   | 0.019 / 0.002  |

For the beta precision, INLA’s default log-gamma$`(1, 0.1)`$ covered as
well as half-normal$`(10)`$ or better (area SD 0.3 / 0.8):

| True $`\phi`$ | Log-gamma$`(1, 0.1)`$ | Half-normal$`(10)`$ |
|---------------|-----------------------|---------------------|
| 10            | 0.960 / 0.937         | 0.965 / 0.938       |
| 50            | 0.931 / 0.903         | 0.913 / 0.869       |
| 200           | 0.875 / 0.925         | 0.841 / 0.871       |

so both engines use it. Neither prior covers well when the sampling
noise is small ($`\phi = 200`$), and the areas without a direct estimate
then covered only 73-80%; give `vardir` or `trials` when they are known.

### Beta precision from `vardir`

The beta variance is $`\mu_d(1 - \mu_d)/(\phi_d + 1)`$, so a known
sampling variance fixes $`\phi_d = \mu_d(1 - \mu_d)/V_d - 1`$ at the
unknown true proportion. Putting the direct estimate $`y_d`$ in place of
$`\mu_d`$ (Janicki 2020) gives the model a sampling variance of
$`V_d\,\mu_d(1 - \mu_d)/\{y_d(1 - y_d)\}`$ instead of $`V_d`$, too small
wherever $`y_d`$ happens to fall far from $`\mu_d`$.
[`hb_area()`](https://ridsonap.github.io/fastsae/reference/hb_area.md)
therefore fits the model twice and takes $`\phi_d`$ at the posterior
mean of the first fit. tipsae (`type_disp = "var"`) instead lets
$`\phi_d`$ follow the model mean inside the MCMC; that keeps the
variance at $`V_d`$ exactly, but a low mean then also lowers the
precision, which gives the likelihood a second mode near 0. Coverage of
95% intervals for sampled areas, area SD 0.3 / 0.8 (`beta_vardir_phi.R`;
30 areas, 60 replicates, $`V_d`$ the true sampling variance):

| Proportions near | Direct estimate | Two stages    | tipsae        | True $`\phi_d`$ |
|------------------|-----------------|---------------|---------------|-----------------|
| 0.08             | 0.844 / 0.812   | 0.938 / 0.896 | 0.938 / 0.931 | 0.953 / 0.918   |
| 0.23             | 0.922 / 0.931   | 0.948 / 0.933 | 0.951 / 0.954 | 0.950 / 0.944   |

tipsae covered as well, but 11 of its 240 fits did not converge (R-hat
up to 8.7).

## Changing the priors

[`inla_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)
takes INLA prior specifications for every hyperparameter. A prior given
there replaces the default.

``` r

library(fastsae)
data("mys")

# PC prior on the area SD instead of the half-normal default
fit <- hb_area(y ~ x1 + x2 + x3,
  data = mys, vardir = "vardir",
  engine = inla_engine(prior_prec = list(prior = "pc.prec", param = c(1, 0.01)))
)

# INLA's own prior on the negative binomial size (d: counts y, covariate x, exposure E)
fit_nb <- hb_area(y ~ x,
  data = d, family = "nbinomial", exposure = "E",
  engine = inla_engine(control.family = list(hyper = list(theta = list(prior = "pc.mgamma", param = 7))))
)
```

The other arguments are `prior_phi` (BYM2), `prior_rho` (Leroux, SAR),
`prior_prec_time` and `prior_rho_time` (temporal effect) and
`prior_prec_subarea`
([`hb_twofold()`](https://ridsonap.github.io/fastsae/reference/hb_twofold.md)).
[`jax_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)
offers the two prior sets of the first table:

``` r

fit_jax <- hb_area(y ~ x1 + x2 + x3,
  data = mys, vardir = "vardir",
  engine = jax_engine(priors = "weakly_informative")
)
```

## References

- Gelman, A. (2006). Prior distributions for variance parameters in
  hierarchical models. *Bayesian Analysis*, 1(3), 515-534.
- Janicki, R. (2020). Properties of the beta regression model for small
  area estimation of proportions and application to estimation of
  poverty rates. *Communications in Statistics - Theory and Methods*,
  49(9), 2264-2284.
- Piironen, J. and Vehtari, A. (2017). Sparsity information and
  regularization in the horseshoe and other shrinkage priors.
  *Electronic Journal of Statistics*, 11(2), 5018-5051.
- Riebler, A., Sørbye, S. H., Simpson, D. and Rue, H. (2016). An
  intuitive Bayesian spatial model for disease mapping that accounts for
  scaling. *Statistical Methods in Medical Research*, 25(4), 1145-1165.
- Simpson, D., Rue, H., Riebler, A., Martins, T. G. and Sørbye, S. H.
  (2017). Penalising model component complexity: a principled approach
  to constructing priors. *Statistical Science*, 32(1), 1-28.
