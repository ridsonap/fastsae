"""
GPU-Accelerated Hierarchical Bayesian Small Area Estimation (SAE) via NumPyro (JAX).
Supports Apple Silicon (Metal) and NVIDIA (CUDA) with graceful CPU fallback.
Implements non-centered parameterization, stable AR(1)/RW(1) operators,
and fast Kronecker tensor contraction for spatio-temporal interactions.
"""

import os
import sys
import numpy as np
from functools import partial
from scipy.special import logsumexp
from scipy.sparse import csr_matrix

def _backend_devices(jax, name):
    """Devices of a JAX backend; jax-metal registers its platform as 'METAL'."""
    for n in (name, name.upper()):
        try:
            devs = jax.devices(n)
            if len(devs) > 0:
                return devs
        except Exception:
            pass
    return []

def _popcount32(x):
    """Bit count of non-negative int32 values with shifts and masks only."""
    import jax.numpy as jnp
    x = x.astype(jnp.int32)
    x = x - ((x >> 1) & 0x55555555)
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333)
    x = (x + (x >> 4)) & 0x0F0F0F0F
    x = x + (x >> 8)
    x = x + (x >> 16)
    return x & 0x3F

def _leaf_idx_to_ckpt_idxs(n):
    """numpyro.infer.hmc_util._leaf_idx_to_ckpt_idxs without jnp.bitwise_count."""
    idx_max = _popcount32(n >> 1)
    num_subtrees = _popcount32((~n & (n + 1)) - 1)
    return idx_max - num_subtrees + 1, idx_max

def _patch_numpyro_for_metal():
    """Work around the two jax-metal gaps that stop NumPyro's NUTS from compiling.

    1. jax-metal cannot legalize ``popcnt`` (``jnp.bitwise_count``), which NUTS uses
       to index its U-turn checkpoints; it is replaced by an equivalent shift/mask count.
    2. jax-metal cannot serialize a ``lax.cond`` in which an output is a literal in
       every branch ("Unable to serialize MPS module"); each branch output is routed
       through ``lax.optimization_barrier``, which leaves the values unchanged.
    """
    from jax import lax
    from numpyro.infer import hmc_util
    if getattr(lax.cond, "_fastsae_metal", False):
        return
    import jax
    import jax.numpy as jnp
    hmc_util._leaf_idx_to_ckpt_idxs = _leaf_idx_to_ckpt_idxs
    cond = lax.cond

    def _barrier(fun):
        return lambda *a: jax.tree_util.tree_map(
            lambda v: lax.optimization_barrier(jnp.asarray(v)), fun(*a)
        )

    def metal_cond(pred, *args, **kwargs):
        if len(args) == 4 and not callable(args[0]):
            # NumPyro's form: cond(pred, true_operand, true_fun, false_operand, false_fun)
            t_op, t_fun, f_op, f_fun = args
            return cond(pred, t_op, _barrier(t_fun), f_op, _barrier(f_fun))
        t_fun, f_fun, *operands = args
        return cond(pred, _barrier(t_fun), _barrier(f_fun), *operands, **kwargs)

    metal_cond._fastsae_metal = True
    lax.cond = metal_cond

def init_jax_environment(device="auto"):
    """Initializes JAX on the requested accelerator with graceful CPU fallback."""
    preferred = device
    if preferred == "auto":
        # On Linux/Windows, try CUDA. On macOS, default to optimized CPU unless metal requested
        if sys.platform != "darwin":
            preferred = "cuda"
        else:
            preferred = "cpu"

    # JAX_PLATFORMS is only read when JAX initializes its backends
    if "jax" not in sys.modules:
        if preferred == "cpu":
            os.environ["JAX_PLATFORMS"] = "cpu"
        elif preferred == "metal":
            os.environ["JAX_PLATFORMS"] = "METAL,cpu"
        elif preferred in ("cuda", "gpu"):
            os.environ["JAX_PLATFORMS"] = "cuda,cpu"

    import jax
    target = "metal" if preferred == "metal" else ("cuda" if preferred in ("cuda", "gpu") else "cpu")
    devs = _backend_devices(jax, target) if target != "cpu" else []
    if devs:
        active_dev_str = f"{target}:{devs[0]}"
    else:
        devs = _backend_devices(jax, "cpu") or jax.devices()
        active_dev_str = "cpu" if devs[0].platform.lower() == "cpu" else f"{devs[0].platform}:{devs[0]}"
    # Pin computation to the selected device; JAX's default backend may differ
    # (e.g. METAL becomes the default whenever jax-metal is installed)
    jax.config.update("jax_default_device", devs[0])
    if active_dev_str.startswith("metal"):
        _patch_numpyro_for_metal()

    import jax.numpy as jnp
    import numpyro
    import numpyro.distributions as dist
    from numpyro.infer import MCMC, NUTS, init_to_median
    
    return jax, jnp, numpyro, dist, MCMC, NUTS, init_to_median, active_dev_str

def horseshoe_shrinkage(tau, lam, c2, xp=np):
    """Regularized horseshoe (Piironen & Vehtari 2017): local scale lambda_tilde and
    shrinkage factor kappa = 1 / (1 + tau^2 lambda_tilde^2), from the scale actually used by the prior."""
    lam_tilde = xp.sqrt(c2) * lam / xp.sqrt(c2 + tau ** 2 * lam ** 2)
    return lam_tilde, 1.0 / (1.0 + (tau * lam_tilde) ** 2)

def init_coefficients_at_zero(site=None, init_to_median=None):
    """NumPyro init strategy: coefficients at 0, every other site at its prior median.
    The median of 15 prior draws of an N(0, 1000) coefficient can be far out (|beta| ~ 10
    on a log scale), where a sharp likelihood (a gamma shape known from vardir) left
    chains stuck."""
    if site is None:
        return partial(init_coefficients_at_zero, init_to_median=init_to_median)
    if site["type"] == "sample" and not site["is_observed"] and site["name"] in ("beta", "beta_0"):
        import jax.numpy as jnp_
        return jnp_.zeros(site["fn"].shape())
    return init_to_median(site)


def pc_bym2_phi_setup(laplacian_eigenvalues, scale_factor, u=0.5, alpha=0.5):
    """PC prior of the BYM2 mixing parameter phi (Riebler et al., 2016), as INLA's
    inla.pc.bym.phi(): an exponential prior with rate lam on the distance
    d(phi) = sqrt(2 KLD), 2 KLD = sum_i [phi (g_i - 1) - log(1 + phi (g_i - 1))], where
    g_i are the eigenvalues of the scaled ICAR covariance (0 in its null space), and
    lam such that P(phi < u) = alpha. Returns (g - 1, lam)."""
    ev = np.asarray(laplacian_eigenvalues, dtype=np.float64)
    g = np.where(ev > 1e-6, 1.0 / np.maximum(ev, 1e-6) / scale_factor, 0.0)
    gm1 = g - 1.0
    x = u * gm1
    d_u = np.sqrt(max(np.sum(x - np.log1p(x)), 0.0))
    return gm1, float(-np.log(1.0 - alpha) / d_u)


def pc_bym2_phi_logdens(phi, gm1, lam, xp=np):
    """Log density of the PC prior of phi (see pc_bym2_phi_setup)."""
    x = phi * gm1
    d = xp.sqrt(xp.maximum(xp.sum(x - xp.log1p(x)), 1e-30))
    dd = 0.5 * xp.sum(gm1 * x / (1.0 + x)) / d  # d'(phi)
    return xp.log(lam) - lam * d + xp.log(xp.abs(dd))


def build_temporal_operator(T, temporal, rho, jnp):
    """
    Constructs a numerically stable non-centered 1D temporal operator (T x T).
    Avoids 0**0 gradient NaN singularity in JAX reverse-mode AD.
    """
    if temporal in ("none", "iid") or T <= 1:
        return jnp.eye(T)
    elif temporal == "rw1":
        # Cumulative sum with sum-to-zero projection across time periods
        # Ensures strict identifiability with fixed intercept
        L_cumsum = jnp.tril(jnp.ones((T, T)))
        P_zero = jnp.eye(T) - (1.0 / T) * jnp.ones((T, T))
        return jnp.dot(P_zero, L_cumsum)
    elif temporal == "ar1":
        # Lower triangular Toeplitz AR(1) operator
        idx = jnp.arange(T)
        diff = idx[:, None] - idx[None, :]
        mask = diff > 0
        # Clip powers to avoid rho**0 singularity at rho = 0
        powers_clipped = jnp.where(mask, diff, 1)
        # Explicit 1.0 on diagonal
        L = jnp.where(mask, rho ** powers_clipped, jnp.eye(T))
        # Scale column 0 for stationary distribution
        scale_col0 = 1.0 / jnp.sqrt(jnp.maximum(1.0 - rho**2, 1e-4))
        L = L.at[:, 0].multiply(scale_col0)
        return L
    return jnp.eye(T)

# ---------------------------------------------------------------------------
# Convergence diagnostics of Vehtari, Gelman, Simpson, Carpenter and Buerkner (2021),
# "Rank-normalization, folding, and localization: An improved R-hat for assessing
# convergence of MCMC", Bayesian Analysis 16(2). Same computations as rhat(),
# ess_bulk() and ess_tail() of the R package posterior. Draws are arrays of
# shape (chains, iterations, K), one column per parameter.
# ---------------------------------------------------------------------------

def _split_chains(x):
    """Each chain split in halves (the middle iteration dropped when odd)."""
    n = x.shape[1]
    if n == 1:
        return x
    return np.concatenate([x[:, :n // 2], x[:, int(np.ceil(n / 2)):]], axis=0)

def _is_constant(x):
    return ~np.all(np.isfinite(x), axis=(0, 1)) | (np.ptp(x.reshape(-1, x.shape[2]), axis=0) <= np.finfo(float).eps)

def _rhat_basic(x):
    m, n = x.shape[0], x.shape[1]
    chain_mean = x.mean(axis=1)
    var_between = n * np.var(chain_mean, axis=0, ddof=1)
    var_within = np.mean(np.var(x, axis=1, ddof=1), axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.sqrt((var_between / var_within + n - 1) / n)

def _z_scale(x):
    from scipy.stats import rankdata
    from scipy.special import ndtri
    m, n, k = x.shape
    r = rankdata(x.reshape(m * n, k), method="average", axis=0)
    return ndtri((r - 0.375) / (m * n - 0.75 + 1.0)).reshape(m, n, k)

def _autocovariance(x):
    """Biased autocovariances along axis 1 (mean-centred, divided by N), via the FFT."""
    n = x.shape[1]
    size = 1 << int(np.ceil(np.log2(2 * n)))
    yc = x - x.mean(axis=1, keepdims=True)
    f = np.fft.rfft(yc, n=size, axis=1)
    return np.fft.irfft(f * np.conj(f), n=size, axis=1)[:, :n] / n

def _ess_one(acov_means, chain_means, n, m):
    """ESS of one parameter from its chain-averaged autocovariances (Geyer's initial
    monotone sequence), as posterior:::.ess."""
    mean_var = acov_means[0] * n / (n - 1)
    var_plus = mean_var * (n - 1) / n
    if m > 1:
        var_plus += np.var(chain_means, ddof=1)
    rho = np.zeros(n)
    t = 0
    rho_even = 1.0
    rho[0] = rho_even
    rho_odd = 1.0 - (mean_var - acov_means[1]) / var_plus
    rho[1] = rho_odd
    while t < n - 5 and not np.isnan(rho_even + rho_odd) and rho_even + rho_odd > 0:
        t += 2
        rho_even = 1.0 - (mean_var - acov_means[t]) / var_plus
        rho_odd = 1.0 - (mean_var - acov_means[t + 1]) / var_plus
        if rho_even + rho_odd >= 0:
            rho[t] = rho_even
            rho[t + 1] = rho_odd
    max_t = t
    if rho_even > 0:
        rho[max_t] = rho_even
    t = 0
    while t <= max_t - 4:
        t += 2
        if rho[t] + rho[t + 1] > rho[t - 2] + rho[t - 1]:
            rho[t] = (rho[t - 2] + rho[t - 1]) / 2.0
            rho[t + 1] = rho[t]
    ess = m * n
    tau = -1.0 + 2.0 * np.sum(rho[:max_t]) + rho[max_t]
    return ess / max(tau, 1.0 / np.log10(ess))

def _ess(x):
    m, n, k = x.shape
    out = np.full(k, np.nan)
    if n < 3:
        return out
    acov = _autocovariance(x).mean(axis=0)
    chain_means = x.mean(axis=1)
    const = _is_constant(x)
    for j in range(k):
        if not const[j]:
            out[j] = _ess_one(acov[:, j], chain_means[:, j], n, m)
    return out

def rank_rhat(x):
    """Rank-normalized split R-hat: the larger of the bulk and the folded (tail) R-hat."""
    folded = np.abs(x - np.median(x.reshape(-1, x.shape[2]), axis=0))
    out = np.maximum(_rhat_basic(_z_scale(_split_chains(x))), _rhat_basic(_z_scale(_split_chains(folded))))
    out[_is_constant(x)] = np.nan
    return out

def ess_bulk(x):
    out = _ess(_z_scale(_split_chains(x)))
    out[_is_constant(x)] = np.nan
    return out

def ess_tail(x):
    flat = x.reshape(-1, x.shape[2])
    res = []
    for prob in (0.05, 0.95):
        q = np.quantile(flat, prob, axis=0)
        res.append(_ess(_split_chains((x <= q).astype(np.float64))))
    out = np.minimum(res[0], res[1])
    out[_is_constant(x)] = np.nan
    return out

# Scalar parameters reported with their diagnostics (backend name -> reported name)
_DIAG_SCALARS = (
    ("sigma2_u", "sigma2_u"), ("sigma2_spatial", "sigma2_spatial"), ("sigma2_iid", "sigma2_iid"),
    ("sigma2_t", "sigma2_t"), ("sigma2_subarea", "sigma2_subarea"), ("sigma2_e", "sigma2_e"),
    ("phi_est", "phi"), ("rho_spatial", "rho_spatial"), ("rho_t", "rho_t"), ("nu_u", "nu"),
    ("tau_horseshoe", "tau_horseshoe"), ("alpha_dispersion", "alpha_dispersion"),
    ("shape_param", "shape_param"), ("phi_beta", "phi_beta"),
)

def mcmc_diagnostics(samples, num_chains, hb_rows=None, max_trace_beta=8):
    """R-hat, bulk and tail ESS of the hyperparameters, coefficients and area parameters
    (rows hb_rows of hb_est, or all), and the per-chain draws of the hyperparameters,
    the first coefficients and the area parameter with the largest R-hat (trace)."""
    blocks = []
    for key, name in _DIAG_SCALARS:
        if key in samples:
            blocks.append(([name], np.asarray(samples[key], dtype=np.float64).reshape(-1, 1)))
    n_hyper = sum(len(b[0]) for b in blocks)
    beta = np.asarray(samples["beta"], dtype=np.float64)
    blocks.append((["beta[%d]" % (j + 1) for j in range(beta.shape[1])], beta))
    hb = np.asarray(samples["hb_est"], dtype=np.float64)
    if hb_rows is not None:
        hb = hb[:, np.asarray(hb_rows, dtype=np.int64)]
    blocks.append((["hb[%d]" % (j + 1) for j in range(hb.shape[1])], hb))

    names = [nm for b in blocks for nm in b[0]]
    flat = np.concatenate([b[1] for b in blocks], axis=1)
    x = flat.reshape(num_chains, -1, flat.shape[1])
    rhat = rank_rhat(x)
    table = {
        "parameter": names,
        "mean": flat.mean(axis=0),
        "sd": flat.std(axis=0, ddof=1),
        "rhat": rhat,
        "ess_bulk": ess_bulk(x),
        "ess_tail": ess_tail(x),
    }

    n_beta = beta.shape[1]
    trace_idx = list(range(n_hyper + min(n_beta, max_trace_beta)))
    hb_rhat = rhat[n_hyper + n_beta:]
    if hb_rhat.size > 0:
        worst = int(np.nanargmax(hb_rhat)) if np.any(np.isfinite(hb_rhat)) else 0
        trace_idx.append(n_hyper + n_beta + worst)
    trace = {"names": [names[i] for i in trace_idx], "draws": x[:, :, trace_idx]}
    return table, trace

def fit_numpyro_hb(
    y,
    X,
    domain_idx,
    time_idx=None,
    subarea_idx=None,
    num_subareas=1,
    vardir=None,
    trials=None,
    exposure=None,
    D=1,
    T=1,
    family="gaussian",
    spatial="none",
    temporal="none",
    st_interaction="none",
    W_adj=None,
    scale_factor=1.0,
    eig_values=None,
    eig_vectors=None,
    W_sar=None,
    nb_size=None,
    num_warmup=500,
    num_samples=1000,
    num_chains=2,
    device="auto",
    seed=42,
    benchmark=False,
    benchmark_weights=None,
    benchmark_target=None,
    benchmark_method="logit",
    prior_beta="normal",
    priors="inla",
    pc_phi=(0.5, 0.5),
    robust=False,
    twofold_weights=None,
    return_draws=False,
    draw_rows=None,
    progress_bar=True,
    max_tree_depth=10
):
    jax, jnp, numpyro, dist, MCMC, NUTS, init_to_median, dev_str = init_jax_environment(device)

    N, P = X.shape
    y_np = np.asarray(y, dtype=np.float32)
    X_jnp = jnp.asarray(X, dtype=jnp.float32)
    domain_jnp = jnp.asarray(domain_idx, dtype=jnp.int32)
    time_jnp = jnp.asarray(time_idx, dtype=jnp.int32) if time_idx is not None else None
    subarea_jnp = jnp.asarray(subarea_idx, dtype=jnp.int32) if subarea_idx is not None else None
    
    # Missing / unsampled domains mask
    observed_mask = ~np.isnan(y_np)
    obs_idx = np.where(observed_mask)[0]
    y_obs = jnp.asarray(y_np[obs_idx], dtype=jnp.float32)

    # Convert optional vectors safely
    vardir_np = np.asarray(vardir, dtype=np.float32) if vardir is not None else None
    trials_np = np.asarray(trials, dtype=np.float32) if trials is not None else None
    exposure_np = np.asarray(exposure, dtype=np.float32) if exposure is not None else None

    vardir_obs = jnp.asarray(vardir_np[obs_idx], dtype=jnp.float32) if vardir_np is not None else None

    # Gaussian: sample on a standardized scale, y_std = (y - y_c) / y_s, so the
    # unit-scale priors and NUTS adaptation work for any data scale (e.g. incomes ~1e6).
    # Draws are mapped back to the original scale right after sampling.
    X_np = np.asarray(X, dtype=np.float64)
    int_cols = np.where(np.all(X_np == 1.0, axis=0))[0]
    int_mask = np.zeros(P)
    if len(int_cols) > 0:
        int_mask[int_cols[0]] = 1.0
    y_c, y_s = 0.0, 1.0
    if family == "gaussian" and len(obs_idx) > 1:
        if len(int_cols) > 0:
            y_c = float(np.mean(y_np[obs_idx]))
        y_s = float(np.std(y_np[obs_idx]))
        if not np.isfinite(y_s) or y_s <= 0:
            y_s = 1.0
        y_obs = (y_obs - y_c) / y_s
        if vardir_obs is not None:
            vardir_obs = vardir_obs / (y_s ** 2)
    trials_obs = jnp.asarray(trials_np[obs_idx], dtype=jnp.float32) if trials_np is not None else None
    exposure_obs = jnp.asarray(exposure_np[obs_idx], dtype=jnp.float32) if exposure_np is not None else None

    # Spatial Spectral Decomposition (ICAR & Leroux basis)
    # Reuse eigendecomposition from R if available to avoid duplicate O(D^3) eigh
    icar_basis = None
    icar_rank = 0
    evals_all = None
    evecs_all = None
    if spatial in ("besag", "bym2", "bym", "leroux") and W_adj is not None:
        if eig_values is not None and eig_vectors is not None:
            evals = np.asarray(eig_values, dtype=np.float64)
            evecs = np.asarray(eig_vectors, dtype=np.float64)
        else:
            deg = np.sum(W_adj, axis=1)
            L_spatial = np.diag(deg) - W_adj
            evals, evecs = np.linalg.eigh(L_spatial)
        pos = evals > 1e-6
        evals_pos = evals[pos]
        evecs_pos = evecs[:, pos]
        # Unit-variance scaled ICAR basis (Riebler et al., 2016)
        icar_basis = jnp.asarray(evecs_pos / np.sqrt(evals_pos * scale_factor), dtype=jnp.float32)
        icar_rank = icar_basis.shape[1]
        evals_all = jnp.asarray(evals, dtype=jnp.float32)
        evecs_all = jnp.asarray(evecs, dtype=jnp.float32)

    # SAR: u = (I - rho W)^{-1} eps. When W has a real spectrum, W = V diag(lam) V^{-1}
    # (e.g. a row-standardised symmetric graph) and the solve is two matrix-vector
    # products for any rho; otherwise it is a linear solve per evaluation.
    sar_V = sar_Vinv = sar_lam = sar_W = None
    if spatial == "sar" and W_sar is not None:
        W_np = np.asarray(W_sar, dtype=np.float64)
        lam, V = np.linalg.eig(W_np)
        if (np.max(np.abs(lam.imag)) < 1e-8 and np.max(np.abs(V.imag)) < 1e-8
                and np.linalg.cond(V.real) < 1e6):
            sar_V = jnp.asarray(V.real, dtype=jnp.float32)
            sar_Vinv = jnp.asarray(np.linalg.inv(V.real), dtype=jnp.float32)
            sar_lam = jnp.asarray(lam.real, dtype=jnp.float32)
        else:
            sar_W = jnp.asarray(W_np, dtype=jnp.float32)

    pen_idx = np.where(int_mask == 0.0)[0]

    # Negative binomial with vardir: known size per area, computed in R
    nb_size_np = np.asarray(nb_size, dtype=np.float64)[obs_idx] if (family == "nbinomial" and nb_size is not None) else None
    nb_size_obs = jnp.asarray(nb_size_np, dtype=jnp.float32) if nb_size_np is not None else None

    # Gamma with vardir: known shape y^2 / vardir per area (INLA's scale; 1 when invalid)
    gamma_shape_np = None
    gamma_shape_obs = None
    if family == "gamma" and vardir_np is not None:
        y_o = y_np[obs_idx].astype(np.float64)
        v_o = vardir_np[obs_idx].astype(np.float64)
        gamma_shape_np = np.where(np.isfinite(v_o) & (v_o > 0), y_o ** 2 / np.where(v_o > 0, v_o, 1.0), 1.0)
        gamma_shape_np[~np.isfinite(gamma_shape_np) | (gamma_shape_np <= 0)] = 1.0
        gamma_shape_obs = jnp.asarray(gamma_shape_np, dtype=jnp.float32)

    # Default priors. "inla": those of inla_engine(), so both engines fit the same model
    # (.default_prior_prec() in R). Random-effect SDs: half-normal with scale sd(y) for
    # Gaussian models (sigma_std = sigma / y_s ~ HalfNormal(sd(y) / y_s) on the standardised
    # scale) and 1 on the link scale otherwise (sd_scale = y_s = 1).
    # "weakly_informative": half-normal with scale 1 on the (standardised) model scale.
    inla_priors = priors == "inla"
    if family == "gaussian" and len(obs_idx) > 1:
        sd_scale = float(np.std(y_np[obs_idx].astype(np.float64), ddof=1))
        if not np.isfinite(sd_scale) or sd_scale <= 0:
            sd_scale = 1.0
    else:
        sd_scale = 1.0
    def sample_sd(name):
        return numpyro.sample(name, dist.HalfNormal(sd_scale / y_s if inla_priors else 1.0))

    pc_phi_gm1 = None
    pc_phi_lam = 1.0
    if inla_priors and spatial == "bym2" and evals_all is not None:
        gm1_np, pc_phi_lam = pc_bym2_phi_setup(np.asarray(evals_all), scale_factor, pc_phi[0], pc_phi[1])
        pc_phi_gm1 = jnp.asarray(gm1_np, dtype=jnp.float32)

    def model():
        # Regression coefficients
        if inla_priors:
            # INLA: N(0, 1000) on the original scale (its flat intercept also gets it:
            # a wider one let chains drift with the area effects on the centred scale)
            beta_scale = np.sqrt(1000.0) / y_s
        else:
            beta_scale = 2.5 if family in ("beta", "binomial") else 10.0
        if prior_beta == "horseshoe" and len(pen_idx) > 0:
            # Regularized Horseshoe (Finnish Horseshoe) prior (Carvalho et al. 2010; Piironen & Vehtari 2017)
            # Penalized coefficients: every column except the intercept (if any)
            P_pen = len(pen_idx)
            z_beta = numpyro.sample("z_beta", dist.Normal(0.0, 1.0).expand([P_pen]))
            lambda_beta = numpyro.sample("lambda_beta", dist.HalfCauchy(1.0).expand([P_pen]))
            tau_hs = numpyro.sample("tau_hs", dist.HalfCauchy(1.0))
            c2_hs = numpyro.sample("c2_hs", dist.InverseGamma(2.0, 8.0))

            lambda_tilde, kappa = horseshoe_shrinkage(tau_hs, lambda_beta, c2_hs, xp=jnp)
            beta = jnp.zeros(P).at[pen_idx].set(z_beta * tau_hs * lambda_tilde)
            if len(int_cols) > 0:
                # Unpenalized weakly informative normal intercept
                beta_0 = numpyro.sample("beta_0", dist.Normal(0.0, beta_scale))
                beta = beta.at[int_cols[0]].set(beta_0)
            numpyro.deterministic("beta", beta)
            numpyro.deterministic("tau_horseshoe", tau_hs)
            # Full length P; 0 for the unpenalized intercept
            numpyro.deterministic("kappa_shrinkage", jnp.zeros(P).at[pen_idx].set(kappa))
        else:
            beta = numpyro.sample("beta", dist.Normal(0.0, beta_scale).expand([P]))
        linpred_fixed = jnp.dot(X_jnp, beta)

        # 1. Spatial / Domain Random Effect
        u_spatial = jnp.zeros(D)
        
        # Base random effect distribution (Gaussian vs Robust Student-t)
        if robust:
            nu_u = numpyro.sample("nu_u", dist.Uniform(2.5, 30.0))
            dist_rand = dist.StudentT(df=nu_u, loc=0.0, scale=1.0)
        else:
            dist_rand = dist.Normal(0.0, 1.0)

        # Avoid double-counting spatial field if separable spatio-temporal structure is active
        if st_interaction != "separable":
            if spatial == "none":
                sigma_s = sample_sd("sigma_s")
                z_s = numpyro.sample("z_s", dist_rand.expand([D]))
                u_spatial = sigma_s * z_s
                numpyro.deterministic("sigma2_u", sigma_s ** 2)
            elif spatial == "besag":
                sigma_s = sample_sd("sigma_s")
                z_icar = numpyro.sample("z_icar", dist_rand.expand([icar_rank]))
                u_spatial = sigma_s * jnp.dot(icar_basis, z_icar)
                numpyro.deterministic("sigma2_u", sigma_s ** 2)
            elif spatial == "bym":
                sigma_s = sample_sd("sigma_s")
                sigma_iid = sample_sd("sigma_iid")
                z_icar = numpyro.sample("z_icar", dist_rand.expand([icar_rank]))
                z_iid = numpyro.sample("z_iid", dist_rand.expand([D]))
                u_spatial = sigma_s * jnp.dot(icar_basis, z_icar) + sigma_iid * z_iid
                numpyro.deterministic("sigma2_spatial", sigma_s ** 2)
                numpyro.deterministic("sigma2_iid", sigma_iid ** 2)
                numpyro.deterministic("sigma2_u", sigma_s ** 2 + sigma_iid ** 2)
            elif spatial == "bym2":
                sigma_s = sample_sd("sigma_s")
                if pc_phi_gm1 is not None:
                    # On logit(phi) in (-15, 12), the support of INLA's tabulated prior: the
                    # untruncated PC prior puts about 10% of its mass at phi > 1 - 1e-8.
                    phi_logit = numpyro.sample("phi_logit", dist.Uniform(-15.0, 12.0))
                    phi = jax.nn.sigmoid(phi_logit)
                    numpyro.factor(
                        "phi_pc_prior",
                        pc_bym2_phi_logdens(phi, pc_phi_gm1, pc_phi_lam, xp=jnp) + jnp.log(phi) + jnp.log1p(-phi)
                    )
                else:
                    phi = numpyro.sample("phi", dist.Beta(1.0, 1.0))
                z_iid = numpyro.sample("z_iid", dist_rand.expand([D]))
                z_icar = numpyro.sample("z_icar", dist_rand.expand([icar_rank]))
                u_icar = jnp.dot(icar_basis, z_icar)
                u_spatial = sigma_s * (jnp.sqrt(1.0 - phi) * z_iid + jnp.sqrt(phi) * u_icar)
                numpyro.deterministic("sigma2_u", sigma_s ** 2)
                numpyro.deterministic("phi_est", phi)
            elif spatial == "leroux":
                sigma_s = sample_sd("sigma_s")
                rho_leroux = numpyro.sample("rho_leroux", dist.Beta(1.0, 1.0))
                prec_evals = rho_leroux * evals_all + (1.0 - rho_leroux)
                scale_leroux = 1.0 / jnp.sqrt(jnp.maximum(prec_evals, 1e-6))
                z_leroux = numpyro.sample("z_leroux", dist_rand.expand([D]))
                u_spatial = sigma_s * jnp.dot(evecs_all, scale_leroux * z_leroux)
                numpyro.deterministic("sigma2_u", sigma_s ** 2)
                numpyro.deterministic("rho_spatial", rho_leroux)
            elif spatial == "sar":
                sigma_s = sample_sd("sigma_s")
                rho_sar = numpyro.sample("rho_sar", dist.Uniform(-0.999, 0.999))
                z_sar = numpyro.sample("z_sar", dist_rand.expand([D]))
                if sar_V is not None:
                    u_sar = jnp.dot(sar_V, jnp.dot(sar_Vinv, z_sar) / (1.0 - rho_sar * sar_lam))
                else:
                    u_sar = jnp.linalg.solve(jnp.eye(D) - rho_sar * sar_W, z_sar)
                u_spatial = sigma_s * u_sar
                numpyro.deterministic("sigma2_u", sigma_s ** 2)
                numpyro.deterministic("rho_spatial", rho_sar)

        # 2. Temporal Random Effect
        u_temporal = jnp.zeros(T)
        rho_t = 0.0
        if temporal != "none" and T > 1:
            sigma_t = sample_sd("sigma_t")
            numpyro.deterministic("sigma2_t", sigma_t ** 2)
            if temporal == "ar1":
                if inla_priors:
                    # INLA: N(0, precision 0.15) on log((1 + rho) / (1 - rho))
                    x_rho = numpyro.sample("rho_t_internal", dist.Normal(0.0, 1.0 / np.sqrt(0.15)))
                    rho_t = numpyro.deterministic("rho_t", jnp.tanh(x_rho / 2.0))
                else:
                    rho_t = numpyro.sample("rho_t", dist.Uniform(-0.99, 0.99))
            C_T = build_temporal_operator(T, temporal, rho_t, jnp)
            z_t = numpyro.sample("z_t", dist.Normal(0.0, 1.0).expand([T]))
            u_temporal = sigma_t * jnp.dot(C_T, z_t)
            if temporal == "ar1" and len(int_cols) > 0:
                # Sum-to-zero, as for rw1: the level of a stationary AR(1) over a few
                # periods trades off against the intercept along a ridge whose width
                # depends on rho, which NUTS cannot adapt to (divergent transitions).
                # The intercept then carries the mean level over the periods.
                u_temporal = u_temporal - jnp.mean(u_temporal)

        # 3. Spatio-Temporal Interaction
        u_st = jnp.zeros(N)
        if st_interaction != "none" and temporal != "none" and T > 1:
            sigma_st = sample_sd("sigma_st")
            C_T = build_temporal_operator(T, temporal, rho_t, jnp)

            if st_interaction == "type1":
                # Unstructured space x Unstructured time (I_D (x) I_T)
                z_st = numpyro.sample("z_st", dist.Normal(0.0, 1.0).expand([D, T]))
                u_st_mat = sigma_st * z_st
                u_st = u_st_mat[domain_jnp, time_jnp]
            elif st_interaction in ("type2", "domain-specific"):
                # Unstructured space x Structured time (I_D (x) Q_t)
                z_st = numpyro.sample("z_st", dist.Normal(0.0, 1.0).expand([D, T]))
                u_st_mat = sigma_st * jnp.matmul(z_st, C_T.T)
                u_st = u_st_mat[domain_jnp, time_jnp]
            elif st_interaction == "type3":
                # Structured space x Unstructured time (Q_s (x) I_T)
                if icar_basis is not None:
                    z_st = numpyro.sample("z_st", dist.Normal(0.0, 1.0).expand([icar_rank, T]))
                    u_st_mat = sigma_st * jnp.matmul(icar_basis, z_st)
                    u_st = u_st_mat[domain_jnp, time_jnp]
            elif st_interaction in ("type4", "separable"):
                # Structured space x Structured time (Q_s (x) Q_t via fast Kronecker contraction)
                if icar_basis is not None:
                    z_st = numpyro.sample("z_st", dist.Normal(0.0, 1.0).expand([icar_rank, T]))
                    temp_filt = jnp.matmul(z_st, C_T.T)
                    u_st_mat = sigma_st * jnp.matmul(icar_basis, temp_filt)
                    u_st = u_st_mat[domain_jnp, time_jnp]

        # Total latent linear predictor
        u_area_rows = u_spatial[domain_jnp]
        rand_total = u_area_rows

        # Nested Sub-Area Random Effect (Torabi & Rao, 2014)
        if subarea_jnp is not None and num_subareas > 1:
            sigma_sub = sample_sd("sigma_subarea")
            z_sub = numpyro.sample("z_subarea", dist_rand.expand([num_subareas]))
            u_subarea = sigma_sub * z_sub
            numpyro.deterministic("sigma2_subarea", sigma_sub ** 2)
            rand_total = rand_total + u_subarea[subarea_jnp]

        if time_jnp is not None:
            rand_total = rand_total + u_temporal[time_jnp] + u_st
        else:
            rand_total = rand_total + u_st

        eta_raw = linpred_fixed + rand_total
        eta = jnp.clip(eta_raw, -8.0, 8.0) if family in ("beta", "binomial") else eta_raw
        # Only S x N sites that are not derivable from others are recorded (each is kept on
        # the device for every draw): linear_pred = link(hb_est) is computed after sampling.
        numpyro.deterministic("rand_eff", rand_total)
        # Split components for hb_twofold df_hb (random_effect_area/subarea); nested models only
        if subarea_jnp is not None and num_subareas > 1:
            numpyro.deterministic("rand_eff_area", u_area_rows)
            numpyro.deterministic("rand_eff_subarea", u_subarea[subarea_jnp])

        # Likelihood and observation sampling
        if family == "gaussian":
            mu = eta
            numpyro.deterministic("hb_est", mu)
            if vardir_obs is not None:
                scale_obs = jnp.sqrt(vardir_obs)
            else:
                sigma_e = numpyro.sample("sigma_e", dist.HalfNormal(1.0))
                scale_obs = sigma_e
                numpyro.deterministic("sigma2_e", sigma_e ** 2)
            numpyro.sample("y_obs", dist.Normal(mu[obs_idx], scale_obs), obs=y_obs)

        elif family == "binomial":
            p = jax.nn.sigmoid(eta)
            numpyro.deterministic("hb_est", p)
            numpyro.sample("y_obs", dist.Binomial(total_count=trials_obs, probs=p[obs_idx]), obs=y_obs)

        elif family == "poisson":
            rate = jnp.exp(eta)
            numpyro.deterministic("hb_est", rate)
            lam = rate[obs_idx] * (exposure_obs if exposure_obs is not None else 1.0)
            numpyro.sample("y_obs", dist.Poisson(rate=lam), obs=y_obs)

        elif family == "beta":
            p = jax.nn.sigmoid(eta)
            numpyro.deterministic("hb_est", p)
            if vardir_obs is not None:
                phi_beta = jnp.maximum((y_obs * (1.0 - y_obs) / vardir_obs) - 1.0, 1.0)
            elif trials_obs is not None:
                phi_beta = jnp.maximum(trials_obs - 1.0, 1.0)
            else:
                phi_beta = numpyro.sample("phi_beta", dist.HalfNormal(10.0))
            a = jnp.maximum(p[obs_idx] * phi_beta, 1e-4)
            b = jnp.maximum((1.0 - p[obs_idx]) * phi_beta, 1e-4)
            y_obs_clipped = jnp.clip(y_obs, 1e-5, 1.0 - 1e-5)
            numpyro.sample("y_obs", dist.Beta(a, b), obs=y_obs_clipped)

        elif family == "nbinomial":
            if nb_size_obs is not None:
                # As INLA: the size of each area from vardir (.nb_size_from_vardir in R)
                alpha_nb = nb_size_obs
            else:
                alpha_nb = numpyro.sample("alpha_nb", dist.HalfNormal(10.0))
                numpyro.deterministic("alpha_dispersion", alpha_nb)
            rate = jnp.exp(eta)
            numpyro.deterministic("hb_est", rate)
            lam = rate[obs_idx] * (exposure_obs if exposure_obs is not None else 1.0)
            numpyro.sample("y_obs", dist.NegativeBinomial2(mean=lam, concentration=alpha_nb), obs=y_obs)

        elif family == "gamma":
            if vardir_obs is not None:
                # As INLA: the shape (precision) of each area is y^2 / vardir
                shape_gamma = gamma_shape_obs
            else:
                shape_gamma = numpyro.sample("shape_gamma", dist.HalfNormal(10.0))
                numpyro.deterministic("shape_param", shape_gamma)
            mu_gamma = jnp.exp(eta)
            numpyro.deterministic("hb_est", mu_gamma)
            rate_gamma = shape_gamma / jnp.maximum(mu_gamma[obs_idx], 1e-6)
            numpyro.sample("y_obs", dist.Gamma(concentration=shape_gamma, rate=rate_gamma), obs=y_obs)

    # Run MCMC on GPU/Accelerator
    # On Apple Metal, sequential avoids vmap control-flow shader limits (the
    # vectorized initialisation fails there)
    on_metal = dev_str.lower().startswith("metal")
    pref_chain_method = "sequential" if on_metal else ("vectorized" if num_chains > 1 else "sequential")
    rng_key = jax.random.PRNGKey(seed)
    kernel = NUTS(model, init_strategy=partial(init_coefficients_at_zero, init_to_median=init_to_median),
                  target_accept_prob=0.85, max_tree_depth=max_tree_depth)
    # divergent transitions and the number of leapfrog steps of each iteration (tree depth)
    extra_fields = ("diverging", "num_steps")

    # On Metal, NumPyro must step the chain from Python (its progress-bar loop):
    # the alternative fuses the whole chain into one kernel, which jax-metal
    # never finishes. The bar itself is silenced when no progress is wanted.
    import numpyro.util as numpyro_util
    tqdm_module = numpyro_util.tqdm
    if on_metal and not progress_bar:
        import types
        numpyro_util.tqdm = types.SimpleNamespace(
            trange=lambda *a, **k: tqdm_module.trange(*a, **{**k, "disable": True})
        )
    try:
        mcmc = MCMC(
            kernel,
            num_warmup=num_warmup,
            num_samples=num_samples,
            num_chains=num_chains,
            chain_method=pref_chain_method,
            progress_bar=progress_bar or on_metal
        )
        mcmc.run(rng_key, extra_fields=extra_fields)
    except Exception as e:
        err_msg = str(e)
        if "legalize" in err_msg or "popcnt" in err_msg or "bytecode" in err_msg or "MPS" in err_msg:
            dev_str = "cpu (fallback)"
            cpu_devs = jax.devices("cpu")
            if len(cpu_devs) > 0:
                with jax.default_device(cpu_devs[0]):
                    mcmc = MCMC(
                        kernel,
                        num_warmup=num_warmup,
                        num_samples=num_samples,
                        num_chains=num_chains,
                        chain_method="sequential",
                        progress_bar=progress_bar
                    )
                    mcmc.run(rng_key, extra_fields=extra_fields)
            else:
                raise e
        elif pref_chain_method == "vectorized":
            mcmc = MCMC(
                kernel,
                num_warmup=num_warmup,
                num_samples=num_samples,
                num_chains=num_chains,
                chain_method="sequential",
                progress_bar=progress_bar
            )
            mcmc.run(rng_key, extra_fields=extra_fields)
        else:
            raise e
    finally:
        numpyro_util.tqdm = tqdm_module

    samples = mcmc.get_samples()

    # Divergent transitions and iterations that reached the maximum tree depth
    extra = mcmc.get_extra_fields()
    num_divergent = int(np.sum(np.asarray(extra["diverging"]))) if "diverging" in extra else None
    treedepth_hits = (int(np.sum(np.asarray(extra["num_steps"]) >= 2 ** max_tree_depth - 1))
                      if "num_steps" in extra else None)

    if y_c != 0.0 or y_s != 1.0:
        samples = {k: np.asarray(v) for k, v in samples.items()}
        samples["hb_est"] = y_c + y_s * samples["hb_est"]
        for k in ("rand_eff", "rand_eff_area", "rand_eff_subarea", "sigma_e"):
            if k in samples:
                samples[k] = y_s * samples[k]
        samples["beta"] = y_s * samples["beta"] + y_c * int_mask
        for k in list(samples):
            if k.startswith("sigma2_"):
                samples[k] = (y_s ** 2) * samples[k]

    # Convergence: rank-normalized split R-hat and bulk/tail ESS of every hyperparameter,
    # coefficient and area parameter (samples are chain-major, as get_samples() flattens them)
    diag_table, diag_trace = mcmc_diagnostics(samples, num_chains, hb_rows=draw_rows)
    rhats = np.asarray(diag_table["rhat"])
    rhats = rhats[np.isfinite(rhats)]
    rhat_max = float(np.max(rhats)) if num_chains > 1 and rhats.size > 0 else None

    # Extract posterior statistics — keep as numpy arrays for zero-copy reticulate transfer
    beta_samples = np.asarray(samples["beta"])
    beta_mean = np.mean(beta_samples, axis=0)
    beta_sd = np.std(beta_samples, axis=0)
    beta_ci_lower, beta_ci_upper = np.percentile(beta_samples, [2.5, 97.5], axis=0)

    hb_samples = np.asarray(samples["hb_est"])
    hb_mean = np.mean(hb_samples, axis=0)
    hb_sd = np.std(hb_samples, axis=0)
    hb_ci_lower, hb_ci_upper = np.percentile(hb_samples, [2.5, 97.5], axis=0)

    # In-Model Benchmark / Calibration across MCMC sample draws
    benchmarked = False
    hb_bench_mean = None
    hb_bench_sd = None
    hb_bench_ci_lower = None
    hb_bench_ci_upper = None
    actual_target = None
    used_method = str(benchmark_method).lower() if benchmark_method is not None else "logit"

    if benchmark and benchmark_weights is not None:
        benchmarked = True
        w = np.asarray(benchmark_weights, dtype=np.float64)
        sum_w = np.sum(w)
        w_norm = (w / sum_w) if sum_w > 0 else (np.ones(N) / N)

        if benchmark_target is not None:
            actual_target = float(benchmark_target)
        else:
            # Self-benchmarking: direct survey weighted aggregate
            actual_target = float(np.sum(w_norm * y_np))

        if used_method == "logit":
            p_clip = np.clip(hb_samples, 1e-7, 1.0 - 1e-7)
            logit_s = np.log(p_clip / (1.0 - p_clip))
            delta = np.zeros(hb_samples.shape[0], dtype=np.float64)
            for _ in range(15):
                cur = 1.0 / (1.0 + np.exp(-(logit_s + delta[:, None])))
                f = np.sum(cur * w_norm[None, :], axis=1) - actual_target
                df = np.sum(cur * (1.0 - cur) * w_norm[None, :], axis=1)
                step = f / np.maximum(df, 1e-9)
                delta -= step
                if np.max(np.abs(f)) < 1e-8:
                    break
            hb_bench_samples = 1.0 / (1.0 + np.exp(-(logit_s + delta[:, None])))
        elif used_method == "optimal":
            # Posterior variance where the sampling variance is unknown (unsampled domains)
            psi = np.var(hb_samples, axis=0)
            if vardir_np is not None:
                psi = np.where(np.isfinite(vardir_np), vardir_np, psi)
            denom = np.sum((w_norm ** 2) * psi)
            if denom > 1e-9:
                agg_s = np.sum(hb_samples * w_norm[None, :], axis=1, keepdims=True)
                lambda_s = (actual_target - agg_s) / denom
                hb_bench_samples = hb_samples + lambda_s * (w_norm[None, :] * psi[None, :])
            else:
                agg_s = np.sum(hb_samples * w_norm[None, :], axis=1, keepdims=True)
                hb_bench_samples = hb_samples * (actual_target / np.maximum(agg_s, 1e-8))
        elif used_method == "difference":
            agg_s = np.sum(hb_samples * w_norm[None, :], axis=1, keepdims=True)
            hb_bench_samples = hb_samples + (actual_target - agg_s)
        else:  # ratio
            agg_s = np.sum(hb_samples * w_norm[None, :], axis=1, keepdims=True)
            hb_bench_samples = hb_samples * (actual_target / np.maximum(agg_s, 1e-8))

        hb_bench_mean = np.mean(hb_bench_samples, axis=0)
        hb_bench_sd = np.std(hb_bench_samples, axis=0)
        hb_bench_ci_lower, hb_bench_ci_upper = np.percentile(hb_bench_samples, [2.5, 97.5], axis=0)

    # Twofold subarea -> area aggregation (Torabi & Rao 2014; Rao & Molina 2015 Ch.8)
    # Area mean: theta_j. = sum_k W_jk * theta_jk, weights normalized per area.
    # ponytail: gaussian-only aggregation; extend when non-gaussian twofold needed.
    area_mean = None
    area_sd = None
    area_ci_lower = None
    area_ci_upper = None
    if twofold_weights is not None:
        w_tf = np.asarray(twofold_weights, dtype=np.float64)
        w_tf = np.where(np.isnan(w_tf), 0.0, w_tf)
        domain_np = np.asarray(domain_idx, dtype=np.int64)
        # Normalized weights per area (equal weights if an area's weights sum to 0), applied to
        # all draws at once through a sparse N x D aggregation matrix
        w_sum = np.bincount(domain_np, weights=w_tf, minlength=int(D))[domain_np]
        n_area = np.bincount(domain_np, minlength=int(D))[domain_np]
        w_norm_tf = np.where(w_sum > 0, w_tf / np.where(w_sum > 0, w_sum, 1.0), 1.0 / n_area)
        agg = csr_matrix((w_norm_tf, (np.arange(len(domain_np)), domain_np)), shape=(len(domain_np), int(D)))
        area_draws = np.asarray(agg.T.dot(hb_samples.T).T, dtype=np.float64)
        area_mean = np.mean(area_draws, axis=0)
        area_sd = np.std(area_draws, axis=0)
        area_ci_lower, area_ci_upper = np.percentile(area_draws, [2.5, 97.5], axis=0)

    # Linear predictor on the link scale, recovered from hb_est = inverse_link(eta)
    if family == "gaussian":
        linpred_mean = hb_mean
    elif family in ("beta", "binomial"):
        p64 = np.clip(hb_samples.astype(np.float64), 1e-12, 1.0 - 1e-12)
        linpred_mean = np.mean(np.log(p64) - np.log1p(-p64), axis=0)
    else:
        linpred_mean = np.mean(np.log(np.maximum(hb_samples.astype(np.float64), 1e-300)), axis=0)

    rand_eff_samples = np.asarray(samples["rand_eff"])
    rand_eff_mean = np.mean(rand_eff_samples, axis=0)
    # ponytail: means only; full draws for split effects skipped (recompute from hb draws when needed).
    if "rand_eff_area" in samples:
        rand_eff_area_mean = np.mean(np.asarray(samples["rand_eff_area"]), axis=0)
    else:
        rand_eff_area_mean = None
    if "rand_eff_subarea" in samples:
        rand_eff_subarea_mean = np.mean(np.asarray(samples["rand_eff_subarea"]), axis=0)
    else:
        rand_eff_subarea_mean = None

    # Hyperparameters
    hyperparams = {}
    if "sigma2_u" in samples:
        hyperparams["sigma2_u"] = float(np.mean(samples["sigma2_u"]))
    if "sigma2_spatial" in samples:
        hyperparams["sigma2_spatial"] = float(np.mean(samples["sigma2_spatial"]))
    if "sigma2_iid" in samples:
        hyperparams["sigma2_iid"] = float(np.mean(samples["sigma2_iid"]))
    if "rho_spatial" in samples:
        hyperparams["rho_spatial"] = float(np.mean(samples["rho_spatial"]))
    if "sigma2_t" in samples:
        hyperparams["sigma2_t"] = float(np.mean(samples["sigma2_t"]))
    if "phi_est" in samples:
        hyperparams["phi"] = float(np.mean(samples["phi_est"]))
    if "rho_t" in samples:
        hyperparams["rho_t"] = float(np.mean(samples["rho_t"]))
    if "alpha_dispersion" in samples:
        hyperparams["alpha_dispersion"] = float(np.mean(samples["alpha_dispersion"]))
    if "shape_param" in samples:
        hyperparams["shape_param"] = float(np.mean(samples["shape_param"]))
    if "phi_beta" in samples:
        hyperparams["phi_beta"] = float(np.mean(samples["phi_beta"]))
    if "nu_u" in samples:
        hyperparams["nu_degrees_of_freedom"] = float(np.mean(samples["nu_u"]))
    if "tau_horseshoe" in samples:
        hyperparams["tau_horseshoe"] = float(np.mean(samples["tau_horseshoe"]))
    if "sigma2_e" in samples:
        hyperparams["sigma2_e"] = float(np.mean(samples["sigma2_e"]))
    if "sigma2_subarea" in samples:
        hyperparams["sigma2_subarea"] = float(np.mean(samples["sigma2_subarea"]))
        if "sigma2_u" in hyperparams:
            s2_u = hyperparams["sigma2_u"]
            s2_sub = hyperparams["sigma2_subarea"]
            hyperparams["icc_nested"] = float(s2_u / (s2_u + s2_sub + 1e-8))

    shrinkage_weights = None
    if "kappa_shrinkage" in samples:
        shrinkage_weights = np.mean(np.asarray(samples["kappa_shrinkage"]), axis=0)

    # Accurate Pointwise Log-Likelihood, WAIC, and DIC
    S_total = hb_samples.shape[0]
    y_obs_val = y_np[obs_idx]
    
    if family == "gaussian":
        if vardir_np is not None:
            sigma_val = np.sqrt(vardir_np[obs_idx])
            sigma_hat = sigma_val
        else:
            # Unknown residual SD (unit-level models): per draw for WAIC and the mean deviance,
            # its posterior mean for the plug-in deviance of DIC
            sigma_val = np.asarray(samples["sigma_e"], dtype=np.float64)[:, None]
            sigma_hat = float(np.mean(sigma_val))
        mu_post = hb_samples[:, obs_idx]
        ll = -0.5 * np.log(2.0 * np.pi * sigma_val**2) - 0.5 * ((y_obs_val - mu_post) / sigma_val)**2
        mu_mean = np.mean(mu_post, axis=0)
        ll_mean = -0.5 * np.log(2.0 * np.pi * sigma_hat**2) - 0.5 * ((y_obs_val - mu_mean) / sigma_hat)**2
    elif family == "binomial":
        p_post = np.clip(hb_samples[:, obs_idx], 1e-6, 1.0 - 1e-6)
        n_val = trials_np[obs_idx]
        from scipy.special import gammaln
        log_comb = gammaln(n_val + 1) - gammaln(y_obs_val + 1) - gammaln(n_val - y_obs_val + 1)
        ll = log_comb + y_obs_val * np.log(p_post) + (n_val - y_obs_val) * np.log(1.0 - p_post)
        p_mean = np.mean(p_post, axis=0)
        ll_mean = log_comb + y_obs_val * np.log(p_mean) + (n_val - y_obs_val) * np.log(1.0 - p_mean)
    elif family == "poisson":
        rate_post = np.clip(hb_samples[:, obs_idx], 1e-6, 1e8)
        e_val = exposure_np[obs_idx] if exposure_np is not None else 1.0
        lam_post = rate_post * e_val
        from scipy.special import gammaln
        ll = y_obs_val * np.log(lam_post) - lam_post - gammaln(y_obs_val + 1)
        lam_mean = np.mean(lam_post, axis=0)
        ll_mean = y_obs_val * np.log(lam_mean) - lam_mean - gammaln(y_obs_val + 1)
    elif family == "beta":
        y_val_clip = np.clip(y_obs_val, 1e-5, 1.0 - 1e-5)
        p_post = np.clip(hb_samples[:, obs_idx], 1e-5, 1.0 - 1e-5)
        if vardir_np is not None:
            phi_val = np.maximum((y_obs_val * (1.0 - y_obs_val) / vardir_np[obs_idx]) - 1.0, 1.0)
        else:
            phi_val = float(np.mean(samples.get("phi_beta", 10.0)))
        a_post = np.maximum(p_post * phi_val, 1e-4)
        b_post = np.maximum((1.0 - p_post) * phi_val, 1e-4)
        from scipy.special import betaln
        ll = (a_post - 1.0) * np.log(y_val_clip) + (b_post - 1.0) * np.log(1.0 - y_val_clip) - betaln(a_post, b_post)
        p_mean = np.mean(p_post, axis=0)
        a_mean = np.maximum(p_mean * phi_val, 1e-4)
        b_mean = np.maximum((1.0 - p_mean) * phi_val, 1e-4)
        ll_mean = (a_mean - 1.0) * np.log(y_val_clip) + (b_mean - 1.0) * np.log(1.0 - y_val_clip) - betaln(a_mean, b_mean)
    elif family == "nbinomial":
        from scipy.special import gammaln
        rate_post = np.clip(hb_samples[:, obs_idx], 1e-6, 1e8)
        e_val = exposure_np[obs_idx] if exposure_np is not None else 1.0
        lam_post = rate_post * e_val
        if nb_size_np is not None:
            r_post = nb_size_np[None, :]  # known size per area
        else:
            r_post = np.asarray(samples["alpha_dispersion"])[:, None]  # shape (S_total, 1)
        ll = (gammaln(y_obs_val + r_post) - gammaln(y_obs_val + 1.0) - gammaln(r_post)
              + r_post * (np.log(r_post) - np.log(r_post + lam_post))
              + y_obs_val * (np.log(lam_post) - np.log(r_post + lam_post)))
        lam_mean = np.mean(lam_post, axis=0)
        r_mean = np.mean(r_post, axis=0)
        ll_mean = (gammaln(y_obs_val + r_mean) - gammaln(y_obs_val + 1.0) - gammaln(r_mean)
                   + r_mean * (np.log(r_mean) - np.log(r_mean + lam_mean))
                   + y_obs_val * (np.log(lam_mean) - np.log(r_mean + lam_mean)))
    elif family == "gamma":
        from scipy.special import gammaln
        mu_post = np.clip(hb_samples[:, obs_idx], 1e-6, 1e8)
        if gamma_shape_np is not None:
            alpha_post = gamma_shape_np[None, :]  # known shape per area
        else:
            alpha_post = np.asarray(samples["shape_param"])[:, None]  # shape (S_total, 1)
        rate_post = alpha_post / mu_post
        ll = (alpha_post * np.log(rate_post) - gammaln(alpha_post)
              + (alpha_post - 1.0) * np.log(y_obs_val) - rate_post * y_obs_val)
        mu_mean = np.mean(mu_post, axis=0)
        alpha_mean = np.mean(alpha_post, axis=0)
        rate_mean = alpha_mean / mu_mean
        ll_mean = (alpha_mean * np.log(rate_mean) - gammaln(alpha_mean)
                   + (alpha_mean - 1.0) * np.log(y_obs_val) - rate_mean * y_obs_val)
    else:
        ll = np.zeros((S_total, len(obs_idx)))
        ll_mean = np.zeros(len(obs_idx))

    lppd = float(np.sum(logsumexp(ll, axis=0) - np.log(S_total)))
    p_waic = float(np.sum(np.var(ll, axis=0, ddof=1)))
    waic = float(-2.0 * (lppd - p_waic))

    D_bar = float(-2.0 * np.mean(np.sum(ll, axis=1)))
    D_hat = float(-2.0 * np.sum(ll_mean))
    p_d = float(D_bar - D_hat)
    dic = float(D_bar + p_d)

    return {
        "beta_mean": beta_mean,
        "beta_sd": beta_sd,
        "beta_ci_lower": beta_ci_lower,
        "beta_ci_upper": beta_ci_upper,
        "hb_mean": hb_mean,
        "hb_sd": hb_sd,
        "hb_ci_lower": hb_ci_lower,
        "hb_ci_upper": hb_ci_upper,
        "linpred_mean": linpred_mean,
        "rand_eff_mean": rand_eff_mean,
        "rand_eff_area_mean": rand_eff_area_mean,
        "rand_eff_subarea_mean": rand_eff_subarea_mean,
        "hyperparameters": hyperparams,
        "waic": waic,
        "p_waic": p_waic,
        "dic": dic,
        "p_dic": p_d,
        "device_used": dev_str,
        "benchmarked": benchmarked,
        "benchmark_target": actual_target,
        "benchmark_method": used_method if benchmarked else None,
        "hb_bench_mean": hb_bench_mean,
        "hb_bench_sd": hb_bench_sd,
        "hb_ci_lower_bench": hb_bench_ci_lower,
        "hb_ci_upper_bench": hb_bench_ci_upper,
        "area_mean": area_mean,
        "area_sd": area_sd,
        "area_ci_lower": area_ci_lower,
        "area_ci_upper": area_ci_upper,
        "shrinkage_weights": shrinkage_weights,
        "rhat_max": rhat_max,
        "diagnostics": diag_table,
        "trace": diag_trace,
        "num_divergent": num_divergent,
        "treedepth_hits": treedepth_hits,
        "max_tree_depth": max_tree_depth,
        "hb_draws": _draw_columns(hb_samples, draw_rows) if return_draws else None,
        "sigma2_e_draws": np.asarray(samples["sigma2_e"], dtype=np.float64) if return_draws and "sigma2_e" in samples else None
    }


def _draw_columns(hb_samples, draw_rows):
    """Draws of all rows, or of the rows in draw_rows only (e.g. the D prediction rows of a
    unit-level model, which leaves the many sampled units on the Python side)."""
    if draw_rows is not None:
        hb_samples = hb_samples[:, np.asarray(draw_rows, dtype=np.int64)]
    return np.asarray(hb_samples, dtype=np.float64)
