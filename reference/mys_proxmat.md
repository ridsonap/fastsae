# Example proximity matrix

A 42 by 42 row-standardized spatial proximity matrix for the 42 areas in
`mys`. Every pair of areas has a positive weight, larger for areas whose
centroids (`sim_area$x_coord`, `sim_area$y_coord`) are close, so it
suits the SAR models
([`eblup_sfh`](https://ridsonap.github.io/fastsae/reference/eblup_sfh.md),
[`eblup_stfh`](https://ridsonap.github.io/fastsae/reference/eblup_stfh.md),
`hb_area(spatial = "slm")`). The CAR models of
[`hb_area`](https://ridsonap.github.io/fastsae/reference/hb_area.md) use
only the pattern `W > 0`, which for this matrix is a complete graph;
build a neighbourhood graph for them, for example
`create_weights(cbind(sim_area$x_coord, sim_area$y_coord), method = "knn", k = 4)`.

## Usage

``` r
mys_proxmat
```

## Format

A square numeric matrix with 42 rows and 42 columns with values between
0 and 1.

## Source

Simulated proximity weights between the 42 areas.
