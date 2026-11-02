#' Effective Sample Sizes of Direct Proportions
#'
#' @description Effective sample sizes \eqn{n_d / \mathrm{deff}_d} for binomial area-level models
#' (\code{trials} of \code{\link{hb_area}} and \code{\link{hb_twofold}}), from the direct
#' proportions, their design variances and the sample sizes. The design effect of an area is
#' undefined when its direct proportion is 0 or 1 (no design variance); such areas take a
#' design effect pooled over their group, or a supplied fallback.
#'
#' @param y Vector, column name, or one-sided formula: direct proportions in \eqn{[0, 1]}.
#' @param vardir Vector, column name, or one-sided formula: design variances of \code{y}.
#' @param n Vector, column name, or one-sided formula: sample sizes (number of sampled units
#'   in the denominator of \code{y}).
#' @param group Optional vector, column name, or one-sided formula: groups (e.g. provinces)
#'   within which design effects are pooled.
#' @param data Optional data frame containing the variables.
#' @param method \code{"pooled"} (default): every area takes the design effect of its group,
#'   \eqn{\sum_d v_d / \sum_d \{y_d (1 - y_d) / n_d\}} over the areas of the group with
#'   \eqn{0 < y_d < 1} (the ratio of the summed design variances to the summed simple random
#'   sampling variances). \code{"area"}: each area takes its own design effect
#'   \eqn{v_d / \{y_d (1 - y_d) / n_d\}} where it is defined, and the pooled one otherwise.
#' @param deff_fallback Optional vector, column name, or one-sided formula: design effects for
#'   the areas whose own design effect is undefined, used before the pooled ones with
#'   \code{method = "area"} (for instance the design effect of a related indicator of the same
#'   survey, which is defined where \code{y} is 0).
#' @param deff_min Lower bound of the design effects (default 1, no gain over simple random
#'   sampling).
#'
#' @details Single-area design effects are noisy and undefined at \eqn{y_d \in \{0, 1\}}, which is
#' common for rare indicators. In a simulation of a two-stage cluster sample (12 provinces x 6
#' districts), effective sample sizes with design effects pooled per province gave 95\% interval
#' coverage of 0.925 (districts) and 0.937 (provinces), close to 0.929 and 0.946 with the true
#' design effects, against 0.911 and 0.897 when the design was ignored
#' (\code{benchmarks/validation/binomial_neff.R} in the source repository). A group without any
#' area with \eqn{0 < y_d < 1} takes the design effect pooled over all areas.
#'
#' @return A data frame with one row per area and columns \code{n}, \code{deff} (the design
#'   effect used, at least \code{deff_min}), \code{deff_source} (\code{"area"},
#'   \code{"fallback"}, \code{"group"} or \code{"overall"}) and \code{trials}
#'   (\eqn{n_d / \mathrm{deff}_d}). Areas with missing \code{y} or \code{n} get \code{NA}.
#'
#' @references
#' \enumerate{
#'   \item Kish, L. (1965). \emph{Survey Sampling}. John Wiley & Sons.
#'   \item Chen, C., Wakefield, J., and Lumley, T. (2014). The use of sampling weights in Bayesian
#'     hierarchical models for small area estimation. \emph{Spatial and Spatio-temporal
#'     Epidemiology}, 11, 33--43.
#' }
#'
#' @seealso \code{\link{hb_area}} (\code{family = "binomial"}), \code{\link{gvf_smooth}}.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(
#'   prov = rep(1:3, each = 5), n = rpois(15, 80),
#'   y = c(0, rbeta(14, 2, 30))
#' )
#' d$vardir <- 1.6 * d$y * (1 - d$y) / d$n
#' effective_trials(~y, vardir = ~vardir, n = ~n, group = ~prov, data = d)
#'
#' # own design effects where defined, pooled ones for y = 0
#' effective_trials("y", "vardir", "n", group = "prov", data = d, method = "area")
#' @export
effective_trials <- function(y,
                             vardir,
                             n,
                             group = NULL,
                             data = NULL,
                             method = c("pooled", "area"),
                             deff_fallback = NULL,
                             deff_min = 1) {
  method <- match.arg(method)
  get <- function(v) {
    if (is.null(v) || is.null(data)) {
      if (inherits(v, "formula") || (is.character(v) && length(v) == 1)) {
        cli::cli_abort("Column names and formulas need {.arg data}.")
      }
      return(v)
    }
    .get_variable(data, v)
  }
  y <- as.numeric(get(y))
  vardir <- as.numeric(get(vardir))
  n <- as.numeric(get(n))
  group <- get(group)
  deff_fallback <- if (!is.null(deff_fallback)) as.numeric(get(deff_fallback))
  D <- length(y)
  lens <- c(
    vardir = length(vardir), n = length(n),
    group = if (is.null(group)) D else length(group),
    deff_fallback = if (is.null(deff_fallback)) D else length(deff_fallback)
  )
  if (any(lens != D)) {
    bad <- names(lens)[lens != D]
    cli::cli_abort("{.arg {bad}} must have the length of {.arg y} ({D}).")
  }
  if (!is.numeric(deff_min) || length(deff_min) != 1 || is.na(deff_min) || deff_min <= 0) {
    cli::cli_abort("{.arg deff_min} must be a single positive number.")
  }
  obs <- !is.na(y) & !is.na(n)
  if (any(y[obs] < 0 | y[obs] > 1)) {
    cli::cli_abort("{.arg y} must hold proportions in [0, 1].")
  }
  if (any(n[obs] <= 0)) {
    cli::cli_abort("{.arg n} must be positive in every area with a direct estimate.")
  }
  if (any(vardir[obs] < 0, na.rm = TRUE)) {
    cli::cli_abort("{.arg vardir} must not be negative.")
  }
  if (is.null(group)) group <- rep(1L, D)
  group <- as.character(group)

  srs <- y * (1 - y) / n
  ok <- obs & !is.na(vardir) & vardir > 0 & srs > 0
  if (!any(ok)) {
    cli::cli_abort("No area has a direct proportion strictly between 0 and 1 with a positive {.arg vardir}; design effects cannot be estimated.")
  }
  # Pooled design effects: ratio of summed design variances to summed SRS variances
  overall <- sum(vardir[ok]) / sum(srs[ok])
  pooled <- tapply(vardir[ok], group[ok], sum) / tapply(srs[ok], group[ok], sum)
  deff <- unname(pooled[group])
  source <- ifelse(is.na(deff), "overall", "group")
  deff[is.na(deff)] <- overall
  if (method == "area") {
    if (!is.null(deff_fallback)) {
      fb <- !ok & !is.na(deff_fallback) & deff_fallback > 0
      deff[fb] <- deff_fallback[fb]
      source[fb] <- "fallback"
    }
    deff[ok] <- vardir[ok] / srs[ok]
    source[ok] <- "area"
  }
  deff <- pmax(deff, deff_min)
  trials <- n / deff
  trials[!obs] <- NA_real_
  deff[!obs] <- NA_real_
  source[!obs] <- NA_character_
  data.frame(n = n, deff = deff, deff_source = source, trials = trials, stringsAsFactors = FALSE)
}
