#' @title mys: mean years of schooling people with disabilities.
#' @description A synthetic dataset containing the mean years of schooling for people with disabilities across 42 regencies/municipalities.
#' @format
#' A data frame with 42 rows and 9 variables with 10 domains as non-sampled areas.
#'
#' \describe{
#'   \item{area}{regency municipality identifier}
#'   \item{y}{mean years of schooling people with disabilities (NA for unsampled areas)}
#'   \item{vardir}{variance sampling from the direct estimator for each area}
#'   \item{rse}{relative standard error (\%)}
#'   \item{x1}{Number of Elementary Schools}
#'   \item{x2}{Number of Junior High Schools}
#'   \item{x3}{Number of Senior High Schools}
#'   \item{n}{Number of eligible samples}
#'   \item{weight}{Weight}
#' }
#' @source Simulated based on empirical characteristics of the National Socio-Economic Survey (Susenas), BPS-Statistics Indonesia.
"mys"

#' @title mys_panel: mean years of schooling people with disabilities 2016 - 2026.
#' @description A synthetic panel dataset containing the mean years of schooling for people with disabilities across 42 regencies/municipalities over 11 time periods (2016 - 2026).
#' @format
#' A data frame with 462 rows and 10 variables.
#'
#' \describe{
#'   \item{area}{regency municipality identifier}
#'   \item{year}{year of observation (2016 to 2026)}
#'   \item{y}{mean years of schooling people with disabilities (NA for unsampled areas)}
#'   \item{vardir}{variance sampling from the direct estimator for each area}
#'   \item{rse}{relative standard error (\%)}
#'   \item{x1}{Number of Elementary Schools}
#'   \item{x2}{Number of Junior High Schools}
#'   \item{x3}{Number of Senior High Schools}
#'   \item{n}{Number of eligible samples}
#'   \item{weight}{Weight}
#' }
#' @source Simulated panel data based on empirical characteristics of the National Socio-Economic Survey (Susenas), BPS-Statistics Indonesia.
"mys_panel"

#' @title Example proximity matrix
#' @description A 42 by 42 row-standardized spatial proximity matrix for the 42 areas in \code{mys}.
#' @format A square numeric matrix with 42 rows and 42 columns with values between 0 and 1.
#' @source Simulated based on contiguous administrative boundaries.
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
"cornsoybean"


#' Corn and Soybean Mean Number of Pixels per Segment for 12 Iowa Counties
#'
#' County means of number of pixels per segment of corn and soy beans,
#' from satellite data, for 12 counties in Iowa. The dataset includes
#' population size, sample size, and means of auxiliary variables used in
#' the dataset \code{\link[sae]{cornsoybean}}.
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
"cornsoybeanmeans"


#' Multi-Distribution Synthetic Area-Level Dataset
#'
#' A synthetic dataset containing 42 domains with multiple response variables
#' representing various probability distributions (Gaussian, Poisson, Binomial,
#' Beta, Negative Binomial, and Gamma). Designed for testing and benchmarking
#' generalized linear and spatial small area estimation models (EBP / INLA).
#' It shares the 42-domain spatial structure of \code{\link{mys_proxmat}}.
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
#' @source Simulated using \code{\link{sim_area_data}} based on the spatial
#'   proximity structure of \code{\link{mys_proxmat}}.
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
#' and first-order autoregressive AR(1) temporal dynamics. It shares the 42-domain
#' spatial structure of \code{\link{mys_proxmat}} and is formatted in domain-major order
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
#' @source Simulated using \code{\link{sim_series_data}} based on the spatial
#'   proximity structure of \code{\link{mys_proxmat}}.
#'
#' @keywords datasets
#' @examples
#' data(sim_panel)
#' head(sim_panel)
"sim_panel"
