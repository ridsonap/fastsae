# Fits with few MCMC draws (fast tests that check agreement, not convergence) may
# warn about low ESS or divergent transitions: drop those warnings and messages only.
quiet_mcmc <- function(expr) {
  suppressWarnings(suppressMessages(expr), classes = "fastsae_mcmc_warning")
}
