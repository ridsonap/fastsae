# Corn and Soybean Mean Number of Pixels per Segment for 12 Iowa Counties

County means of number of pixels per segment of corn and soy beans, from
satellite data, for 12 counties in Iowa. The dataset includes population
size, sample size, and means of auxiliary variables used in the dataset
[`cornsoybean`](https://ridsonap.is-a.dev/fastsae/reference/cornsoybean.md).

## Usage

``` r
data(cornsoybeanmeans)
```

## Format

A data frame with 12 observations on the following 6 variables:

- CountyIndex:

  numeric county code.

- CountyName:

  name of the county.

- SampSegments:

  number of sample segments in the county (sample size).

- PopnSegments:

  number of population segments in the county (population size).

- MeanCornPixPerSeg:

  mean number of corn pixels per segment in the county.

- MeanSoyBeansPixPerSeg:

  mean number of soy beans pixels per segment in the county.

## Source

Battese, G.E., Harter, R.M., and Fuller, W.A. (1988). \*An
Error-Components Model for Prediction of County Crop Areas Using Survey
and Satellite Data.\* \*Journal of the American Statistical
Association\*, 83, 28–36.

## Details

This dataset is provided for demonstration purposes and is originally
distributed with the sae package.

## Examples

``` r
data(cornsoybeanmeans)
head(cornsoybeanmeans)
#>   CountyIndex CountyName SampSegments PopnSegments MeanCornPixPerSeg
#> 1           1 CerroGordo            1          545            295.29
#> 2           2   Hamilton            1          566            300.40
#> 3           3      Worth            1          394            289.60
#> 4           4   Humboldt            2          424            290.74
#> 5           5   Franklin            3          564            318.21
#> 6           6 Pocahontas            3          570            257.17
#>   MeanSoyBeansPixPerSeg
#> 1                189.70
#> 2                196.65
#> 3                205.28
#> 4                220.22
#> 5                188.06
#> 6                247.13
```
