# Corn and Soybean Survey and Satellite Data in 12 Iowa Counties

Survey and satellite data for corn and soy beans in 12 Iowa counties,
originally obtained from the 1978 June Enumerative Survey of the U.S.
Department of Agriculture and from LANDSAT satellite observations during
the 1978 growing season.

## Usage

``` r
data(cornsoybean)
```

## Format

A data frame with 37 observations on the following 5 variables:

- County:

  numeric county code.

- CornHec:

  reported hectares of corn from the survey.

- SoyBeansHec:

  reported hectares of soy beans from the survey.

- CornPix:

  number of pixels of corn in the sample segment within county, from
  satellite data.

- SoyBeansPix:

  number of pixels of soy beans in the sample segment within county,
  from satellite data.

## Source

Battese, G.E., Harter, R.M., and Fuller, W.A. (1988). \*An
Error-Components Model for Prediction of County Crop Areas Using Survey
and Satellite Data.\* \*Journal of the American Statistical
Association\*, 83, 28–36.

## Details

This dataset is included for demonstration purposes and is originally
provided in the sae package.

## Examples

``` r
data(cornsoybean)
head(cornsoybean)
#>   County CornHec SoyBeansHec CornPix SoyBeansPix
#> 1      1  165.76        8.09     374          55
#> 2      2   96.32      106.03     209         218
#> 3      3   76.08      103.60     253         250
#> 4      4  185.35        6.47     432          96
#> 5      4  116.43       63.82     367         178
#> 6      5  162.08       43.50     361         137
```
