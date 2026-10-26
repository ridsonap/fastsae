#' @title mys: mean years of schooling of people with disabilities
#' @description A synthetic area-level dataset with direct estimates of the mean years of
#'   schooling of people with disabilities in 42 regencies/municipalities, 10 of which are not
#'   sampled (\code{y = NA}). Used in the examples of the area-level models; its spatial
#'   weights are \code{\link{mys_proxmat}}.
#' @format
#' A data frame (tibble) with 42 rows and 9 variables:
#'
#' \describe{
#'   \item{area}{area (regency/municipality) identifier, 1 to 42.}
#'   \item{y}{direct estimate of the mean years of schooling (\code{NA} for the 10 unsampled
#'     areas).}
#'   \item{vardir}{sampling variance of the direct estimate (0 for the unsampled areas, where
#'     it is not used).}
#'   \item{rse}{relative standard error of the direct estimate (\%),
#'     \eqn{100\sqrt{vardir}/y} (0 for the unsampled areas).}
#'   \item{x1}{number of elementary schools.}
#'   \item{x2}{number of junior high schools.}
#'   \item{x3}{number of senior high schools.}
#'   \item{n}{number of people with disabilities in the Susenas sample of the area (eligible
#'     sample size; synthetic, hence not an integer, and positive also for the unsampled
#'     areas). \code{vardir} was not derived from it.}
#'   \item{weight}{sample share of the area, \code{n / sum(n)}.}
#' }
#' @source Simulated based on empirical characteristics of the National Socio-Economic Survey (Susenas), BPS-Statistics Indonesia.
#' @keywords datasets
#' @examples
#' data(mys)
#' head(mys)
#' # unsampled areas
#' mys$area[is.na(mys$y)]
"mys"

#' @title mys_panel: mean years of schooling of people with disabilities, 2016-2026
#' @description A synthetic panel dataset with direct estimates of the mean years of schooling
#'   of people with disabilities in the 42 areas of \code{\link{mys}} over 11 years
#'   (2016-2026). Areas 21 and 25 are never sampled (\code{y = NA} in every year). The rows
#'   are in time-major order (all areas of 2016, then 2017, ...); \code{\link{eblup_stfh}}
#'   sorts them internally, but \code{sae::eblupSTFH()} needs domain-major order.
#' @format
#' A data frame (tibble) with 462 rows and 10 variables:
#'
#' \describe{
#'   \item{area}{area identifier, 1 to 42.}
#'   \item{year}{year of observation (2016 to 2026).}
#'   \item{y}{direct estimate of the mean years of schooling (\code{NA} for areas 21 and 25).}
#'   \item{vardir}{sampling variance of the direct estimate (\code{NA} where \code{y} is
#'     \code{NA}).}
#'   \item{rse}{relative standard error of the direct estimate (\%).}
#'   \item{x1}{number of elementary schools.}
#'   \item{x2}{number of junior high schools.}
#'   \item{x3}{number of senior high schools.}
#'   \item{n}{number of people with disabilities in the Susenas sample of the area in that
#'     year (eligible sample size; synthetic, hence not an integer, and positive also for
#'     areas 21 and 25). \code{vardir} was not derived from it.}
#'   \item{weight}{sample share, \code{n / sum(n)} over all rows of the panel (the weights sum
#'     to 1 over the panel, not within a year).}
#' }
#' @source Simulated panel data based on empirical characteristics of the National Socio-Economic Survey (Susenas), BPS-Statistics Indonesia.
#' @keywords datasets
#' @examples
#' data(mys_panel)
#' head(mys_panel)
#' table(sampled = !is.na(mys_panel$y), mys_panel$year)[, 1:3]
"mys_panel"

#' @title Example proximity matrix
#' @description A 42 by 42 row-standardized spatial weight matrix for the 42 areas in \code{mys}:
#'   each area is linked to its 4 nearest areas by centroid distance
#'   (\code{sim_area$x_coord}, \code{sim_area$y_coord}), a link in either direction counts,
#'   and each row is divided by its number of neighbours (4 to 8). The weights suit the SAR
#'   models (\code{\link{eblup_sfh}}, \code{\link{eblup_stfh}}, \code{hb_area(spatial = "slm")}),
#'   and the pattern \code{W > 0} is the neighbourhood graph the CAR models of
#'   \code{\link{hb_area}} use (\code{"besag"}, \code{"bym"}, \code{"bym2"}, \code{"leroux"}).
#'   Until fastsae 1.0.0 it was a dense distance-decay matrix that linked every pair of areas,
#'   which a CAR model reads as a complete graph.
#' @format A square numeric matrix with 42 rows and 42 columns with values between 0 and 1;
#'   every row sums to 1.
#' @source Built from the \code{sim_area} centroids by \code{data-raw/mys_proxmat.R} in the
#'   package sources.
#' @keywords datasets
#' @examples
#' data(mys_proxmat)
#' range(rowSums(mys_proxmat))
#' # neighbours per area (the graph used by the CAR models)
#' table(rowSums(mys_proxmat > 0))
"mys_proxmat"


#' Corn and Soybean Survey and Satellite Data in 12 Iowa Counties
#'
#' Survey and satellite data for corn and soy beans in 12 Iowa counties,
#' originally obtained from the 1978 June Enumerative Survey of the U.S.
#' Department of Agriculture and from LANDSAT satellite observations during
#' the 1978 growing season.
#'
#' This dataset is included for demonstration purposes and is originally
#' provided in the \pkg{sae} package.
#'
#' @usage data(cornsoybean)
#'
#' @format A data frame with 37 observations on the following 5 variables:
#' \describe{
#'   \item{County}{numeric county code.}
#'   \item{CornHec}{reported hectares of corn from the survey.}
#'   \item{SoyBeansHec}{reported hectares of soy beans from the survey.}
#'   \item{CornPix}{number of pixels of corn in the sample segment within county, from satellite data.}
#'   \item{SoyBeansPix}{number of pixels of soy beans in the sample segment within county, from satellite data.}
#' }
#'
#' @source
#' Battese, G.E., Harter, R.M., and Fuller, W.A. (1988).
#' *An Error-Components Model for Prediction of County Crop Areas Using Survey and Satellite Data.*
#' *Journal of the American Statistical Association*, 83, 28–36.
#'
#'
#' @keywords datasets
#' @examples
#' data(cornsoybean)
#' head(cornsoybean)
"cornsoybean"


#' Corn and Soybean Mean Number of Pixels per Segment for 12 Iowa Counties
#'
#' County means of number of pixels per segment of corn and soy beans,
#' from satellite data, for 12 counties in Iowa. The dataset includes
#' population size, sample size, and means of auxiliary variables used in
#' the dataset \code{\link{cornsoybean}}.
#'
#' This dataset is provided for demonstration purposes and is originally
#' distributed with the \pkg{sae} package.
#'
#' @usage data(cornsoybeanmeans)
#'
#' @format A data frame with 12 observations on the following 6 variables:
#' \describe{
#'   \item{CountyIndex}{numeric county code.}
#'   \item{CountyName}{name of the county.}
#'   \item{SampSegments}{number of sample segments in the county (sample size).}
#'   \item{PopnSegments}{number of population segments in the county (population size).}
#'   \item{MeanCornPixPerSeg}{mean number of corn pixels per segment in the county.}
#'   \item{MeanSoyBeansPixPerSeg}{mean number of soy beans pixels per segment in the county.}
#' }
#'
#' @source
#' Battese, G.E., Harter, R.M., and Fuller, W.A. (1988).
#' *An Error-Components Model for Prediction of County Crop Areas Using Survey and Satellite Data.*
#' *Journal of the American Statistical Association*, 83, 28–36.
#'
#' @keywords datasets
#' @examples
#' data(cornsoybeanmeans)
#' head(cornsoybeanmeans)
"cornsoybeanmeans"


#' Multi-Distribution Synthetic Area-Level Dataset
#'
#' A synthetic dataset containing 42 domains with multiple response variables
#' representing various probability distributions (Gaussian, Poisson, Binomial,
#' Beta, Negative Binomial, and Gamma). Designed for testing and benchmarking
#' generalized linear and spatial small area estimation models (HB / INLA).
#' Its centroids give the neighbourhood graph of \code{\link{mys_proxmat}}.
#'
#' @usage data(sim_area)
#' @aliases sae_area_multi
#'
#' @format A data frame with 42 rows and 14 variables:
#' \describe{
#'   \item{area}{Integer domain identifier (1 to 42).}
#'   \item{x1}{Continuous explanatory covariate.}
#'   \item{x2}{Uniform explanatory covariate.}
#'   \item{y_gaussian}{Gaussian direct estimator response with unsampled domains as NA.}
#'   \item{vardir}{Direct sampling variance for Gaussian Fay-Herriot model.}
#'   \item{y_poisson}{Count response (Poisson) with unsampled domains as NA.}
#'   \item{exposure}{Expected population count / exposure offset for Poisson model.}
#'   \item{y_binomial}{Number of successes (Binomial) with unsampled domains as NA.}
#'   \item{trials}{Sample size / number of trials for Binomial model.}
#'   \item{y_beta}{Continuous proportion response in (0, 1) for Beta regression.}
#'   \item{y_nbinomial}{Overdispersed count response for Negative Binomial model.}
#'   \item{y_gamma}{Skewed positive continuous response for Gamma regression.}
#'   \item{x_coord}{Centroid x coordinate.}
#'   \item{y_coord}{Centroid y coordinate.}
#' }
#'
#' @source Simulated using \code{\link{sim_area_data}} with the dense distance-decay
#'   weights that \code{\link{mys_proxmat}} held until fastsae 1.0.0; the centroids are a
#'   two-dimensional scaling of those weights.
#'
#' @keywords datasets
#' @examples
#' data(sim_area)
#' head(sim_area)
"sim_area"


#' Multi-Distribution Spatio-Temporal Panel Dataset
#'
#' A synthetic panel dataset containing 42 domains observed over 5 time periods
#' (2022 to 2026, total 210 observations). Features multiple response variables
#' representing various probability distributions (Gaussian, Poisson, Binomial,
#' Beta, Negative Binomial, and Gamma) driven by domain-level spatial autocorrelation
#' and first-order autoregressive AR(1) temporal dynamics. Its 42 areas and centroids are
#' those of \code{\link{sim_area}} (see \code{\link{mys_proxmat}}), in domain-major order
#' for direct use in \code{\link{eblup_stfh}}.
#'
#' @usage data(sim_panel)
#' @aliases sae_panel_multi
#'
#' @format A data frame with 210 rows and 15 variables in domain-major order:
#' \describe{
#'   \item{area}{Integer domain identifier (1 to 42).}
#'   \item{year}{Year of observation (2022 to 2026).}
#'   \item{x1}{Dynamic explanatory covariate with spatial & temporal variation.}
#'   \item{x2}{Dynamic uniform explanatory covariate.}
#'   \item{y_gaussian}{Gaussian direct estimator response with unsampled domains as NA.}
#'   \item{vardir}{Direct sampling variance for Gaussian Fay-Herriot model.}
#'   \item{y_poisson}{Count response (Poisson) with unsampled domains as NA.}
#'   \item{exposure}{Expected population count / exposure offset for Poisson model.}
#'   \item{y_binomial}{Number of successes (Binomial) with unsampled domains as NA.}
#'   \item{trials}{Sample size / number of trials for Binomial model.}
#'   \item{y_beta}{Continuous proportion response in (0, 1) for Beta regression.}
#'   \item{y_nbinomial}{Overdispersed count response for Negative Binomial model.}
#'   \item{y_gamma}{Skewed positive continuous response for Gamma regression.}
#'   \item{x_coord}{Centroid x coordinate.}
#'   \item{y_coord}{Centroid y coordinate.}
#' }
#'
#' @source Simulated using \code{\link{sim_series_data}} with the dense distance-decay
#'   weights that \code{\link{mys_proxmat}} held until fastsae 1.0.0.
#'
#' @keywords datasets
#' @examples
#' data(sim_panel)
#' head(sim_panel)
"sim_panel"
