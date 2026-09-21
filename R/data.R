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
