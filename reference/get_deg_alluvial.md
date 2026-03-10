# Plot alluvial diagram showing gene regulation transitions across comparisons

Plot alluvial diagram showing gene regulation transitions across
comparisons

## Usage

``` r
get_deg_alluvial(x, sample_comparisons, show_other = FALSE)
```

## Arguments

- x:

  A VISTA object (DE results in `comparisons(x)`).

- sample_comparisons:

  Character vector of comparison names to include (\>= 2).

- show_other:

  Logical; include "Other" genes. Default FALSE.

## Value

A ggplot object.

## Examples

``` r
v <- example_vista()
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
if (requireNamespace('ggalluvial', quietly = TRUE)) {
  p <- get_deg_alluvial(v, sample_comparisons = names(comparisons(v)))
  print(p)
}
#> Error in get_deg_alluvial(v, sample_comparisons = names(comparisons(v))): At least two comparisons required.
```
