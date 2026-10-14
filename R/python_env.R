#' Check whether the NumPyro / JAX engine is available
#'
#' Checks that \code{hb_area(engine = "jax")} can run: the package \pkg{reticulate}
#' and a Python with the modules \code{jax} and \code{numpyro}. It looks for the
#' \code{r-numpyro-gpu} virtualenv (see \code{\link{setup_numpyro_env}}) or any
#' Python already configured for \pkg{reticulate}. The check never installs or
#' downloads anything: the managed ('uv') Python of \pkg{reticulate} is not
#' provisioned unless the \code{RETICULATE_USE_MANAGED_VENV} environment variable is set.
#'
#' @param device Character: \code{"auto"}, \code{"metal"}, \code{"cuda"}, or \code{"cpu"}.
#' @return Logical indicating whether the environment is properly configured.
#' @seealso \code{\link{setup_numpyro_env}}, \code{\link{hb_area}}
#' @examples
#' check_numpyro_available()
#' @export
check_numpyro_available <- function(device = "auto") {
  device <- match.arg(device, c("auto", "metal", "cuda", "cpu"))
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    cli::cli_alert_warning("Package {.pkg reticulate} is not installed.")
    return(FALSE)
  }

  if (device == "cpu" || (device == "auto" && Sys.info()["sysname"] == "Darwin")) {
    if (!nzchar(Sys.getenv("JAX_PLATFORMS"))) {
      old_plat <- Sys.getenv("JAX_PLATFORMS", unset = NA)
      on.exit({
        if (is.na(old_plat) || !nzchar(old_plat)) {
          Sys.unsetenv("JAX_PLATFORMS")
        } else {
          Sys.setenv(JAX_PLATFORMS = old_plat)
        }
      }, add = TRUE)
      Sys.setenv(JAX_PLATFORMS = "cpu")
    }
  }

  if (!reticulate::py_available()) {
    try(reticulate::use_virtualenv("r-numpyro-gpu", required = FALSE), silent = TRUE)
    # A check must not download anything: without this, reticulate (>= 1.41)
    # provisions a 'uv'-managed Python into the user cache when none is configured.
    if (is.na(Sys.getenv("RETICULATE_USE_MANAGED_VENV", unset = NA))) {
      Sys.setenv(RETICULATE_USE_MANAGED_VENV = "no")
      on.exit(Sys.unsetenv("RETICULATE_USE_MANAGED_VENV"), add = TRUE)
    }
  }

  has_np <- tryCatch(reticulate::py_module_available("numpyro"), error = function(e) FALSE)
  has_jax <- tryCatch(reticulate::py_module_available("jax"), error = function(e) FALSE)
  if (!has_np || !has_jax) {
    return(FALSE)
  }

  tryCatch({
    jax <- reticulate::import("jax")
    devs <- jax$devices()
    cli::cli_alert_success("JAX devices found: {paste(vapply(devs, as.character, character(1)), collapse = ', ')}")
    TRUE
  }, error = function(e) {
    FALSE
  })
}

#' Set up the Python environment for the NumPyro / JAX engine
#'
#' Creates (or updates) a Python environment with \code{jax} and \code{numpyro} for
#' \code{hb_area(engine = "jax")} and activates it for \pkg{reticulate}.
#'
#' @param envname Name of the conda or virtualenv environment. Default \code{"r-numpyro-gpu"}.
#' @param method \code{"auto"} (virtualenv), \code{"virtualenv"}, or \code{"conda"}.
#' @param device \code{"auto"} (CUDA if \code{nvidia-smi} is found, otherwise CPU), \code{"metal"}
#'   (Apple Silicon), \code{"cuda"} (NVIDIA), or \code{"cpu"}. Versions are pinned to tested
#'   combinations: \code{jax==0.5.0} and \code{numpyro==0.19.0} for CPU/CUDA; \code{jax==0.4.34},
#'   \code{jaxlib==0.4.34}, \code{jax-metal==0.1.1} for Metal.
#' @return Invisible logical \code{TRUE} upon successful setup.
#' @seealso \code{\link{check_numpyro_available}}, \code{\link{hb_area}}
#' @examples
#' \dontrun{
#' setup_numpyro_env(device = "auto")
#' }
#' @export
setup_numpyro_env <- function(envname = "r-numpyro-gpu", method = "auto", device = "auto") {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg reticulate} is required. Install it with {.code install.packages(\"reticulate\")}.")
  }
  method <- match.arg(method, c("auto", "virtualenv", "conda"))
  device <- match.arg(device, c("auto", "metal", "cuda", "cpu"))

  if (device == "auto") {
    # Not "metal" on Apple Silicon: NUTS on Metal is slower than on the CPU
    device <- if (nzchar(Sys.which("nvidia-smi"))) "cuda" else "cpu"
  }

  cli::cli_h2("Configuring NumPyro environment for target device: {.val {device}}")

  # Pinned, tested combinations. jax-metal 0.1.1 (its last release) is built for
  # jax/jaxlib 0.4.34; unpinned installs pull a jax it cannot load.
  pkgs <- c("numpy", "scipy", "numpyro==0.19.0")
  if (device == "metal") {
    pkgs <- c(pkgs, "jax==0.4.34", "jaxlib==0.4.34", "jax-metal==0.1.1")
  } else if (device == "cuda") {
    pkgs <- c(pkgs, "jax[cuda12]==0.5.0")
  } else {
    pkgs <- c(pkgs, "jax==0.5.0", "jaxlib==0.5.0")
  }

  cli::cli_alert_info("Creating/updating {method} environment {.val {envname}} with: {paste(pkgs, collapse = ', ')}")
  if (method == "conda") {
    reticulate::conda_create(envname)
    reticulate::conda_install(envname, packages = pkgs, pip = TRUE)
    reticulate::use_condaenv(envname, required = TRUE)
  } else {
    reticulate::virtualenv_create(envname)
    reticulate::virtualenv_install(envname, packages = pkgs)
    reticulate::use_virtualenv(envname, required = TRUE)
  }

  cli::cli_alert_success("NumPyro environment successfully configured!")
  invisible(TRUE)
}

#' Load the Python NumPyro backend shipped in inst/python
#' @noRd
.get_numpyro_backend <- function(device = "auto") {
  if (!check_numpyro_available(device = device)) {
    cli::cli_abort(c(
      "{.code engine = \"jax\"} needs Python with {.pkg jax} and {.pkg numpyro}, which were not found.",
      "i" = "Run {.code setup_numpyro_env(device = \"{device}\")} or configure a Python environment with jax and numpyro for {.pkg reticulate}."
    ))
  }

  candidates <- c(
    system.file("python", "hb_area_numpyro.py", package = "fastsae"),
    Sys.getenv("FASTSAE_PYTHON_PATH", "")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (length(candidates) == 0) {
    cli::cli_abort("Backend script {.file hb_area_numpyro.py} not found. Reinstall {.pkg fastsae}.")
  }

  reticulate::import_from_path("hb_area_numpyro", path = dirname(candidates[1]))
}

#' Warn when an explicitly requested accelerator was not used
#' @noRd
.warn_device_fallback <- function(device, device_used) {
  if (device == "metal" && grepl("fallback", as.character(device_used))) {
    cli::cli_warn(c(
      "The Metal backend could not compile the model; sampling ran on CPU.",
      "i" = "jax-metal does not implement every JAX operation. Check the installation with {.code setup_numpyro_env(device = 'metal')}."
    ))
  } else if (device %in% c("metal", "cuda") && !grepl(device, as.character(device_used), ignore.case = TRUE)) {
    cli::cli_warn(c(
      "Requested {.code device = '{device}'} is not available to JAX in this session; sampling ran on {.val {device_used}}.",
      "i" = "JAX fixes its platforms when first used in an R session: restart R and request {.val {device}} first, or check the installation with {.code setup_numpyro_env(device = '{device}')}."
    ))
  }
  invisible(NULL)
}
