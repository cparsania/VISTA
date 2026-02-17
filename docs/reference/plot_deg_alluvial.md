# Plot alluvial diagram showing gene regulation transitions across comparisons

Plot alluvial diagram showing gene regulation transitions across
comparisons

## Usage

``` r
plot_deg_alluvial(vista_obj, sample_comparisons, show_other = FALSE)
```

## Arguments

- vista_obj:

  A VISTA object (DE results in metadata(vista_obj)\$de_results).

- sample_comparisons:

  Character vector of comparison names to include (\>= 2).

- show_other:

  Logical; include "Other" genes. Default FALSE.

## Value

A ggplot object.
