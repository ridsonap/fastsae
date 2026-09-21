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
