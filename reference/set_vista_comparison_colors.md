# Set manual comparison colors in a VISTA object

Updates `metadata(x)$comparison$colors` using a user-supplied named
color vector. Comparison colors are used in fold-change and
comparison-level plots. The color map must include all currently active
comparisons.

## Usage

``` r
set_vista_comparison_colors(object, color_map)
```

## Arguments

- object:

  A `VISTA` object.

- color_map:

  Named character vector of colors, with names equal to comparison names
  (for example `"A_VS_B"`).

## Value

A modified `VISTA` object with updated comparison colors.
