# Retrieve stored cell fraction estimates

Retrieve stored cell fraction estimates

## Usage

``` r
get_cell_fractions(x)
```

## Arguments

- x:

  A VISTA object.

## Value

A data.frame of cell fraction estimates with samples in rows.

## Examples

``` r
v <- example_vista()
if (requireNamespace('xCell2', quietly = TRUE)) {
  vx <- try(run_cell_deconvolution(v, method = 'xCell2'), silent = TRUE)
  if (!inherits(vx, 'try-error')) head(get_cell_fractions(vx))
}
#> Starting xCell2 Analysis...
```
