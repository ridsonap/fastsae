# Example proximity matrix

A 42 by 42 row-standardized spatial weight matrix for the 42 areas in
`mys`: each area is linked to its 4 nearest areas by centroid distance
(`sim_area$x_coord`, `sim_area$y_coord`), a link in either direction
counts, and each row is divided by its number of neighbours (4 to 8).
The weights suit the SAR models
([`eblup_sfh`](https://ridsonap.is-a.dev/fastsae/reference/eblup_sfh.md),
[`eblup_stfh`](https://ridsonap.is-a.dev/fastsae/reference/eblup_stfh.md),
`hb_area(spatial = "slm")`), and the pattern `W > 0` is the
neighbourhood graph the CAR models of
[`hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md) use
(`"besag"`, `"bym"`, `"bym2"`, `"leroux"`). Until fastsae 1.0.0 it was a
dense distance-decay matrix that linked every pair of areas, which a CAR
model reads as a complete graph.

## Usage

``` r
mys_proxmat
```

## Format

A square numeric matrix with 42 rows and 42 columns with values between
0 and 1; every row sums to 1.

## Source

Built from the `sim_area` centroids by `data-raw/mys_proxmat.R` in the
package sources.
