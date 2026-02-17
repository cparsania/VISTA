# Set manual group colors in a VISTA object

Updates `metadata(x)$group$colors` using a user-supplied named color
vector. This controls group-level coloring across VISTA plots that use
group mapping (for example PCA, MDS, and expression plots). The color
map must include all groups currently present in the object.

## Usage

``` r
set_vista_group_colors(object, color_map)
```

## Arguments

- object:

  A `VISTA` object.

- color_map:

  Named character vector of colors, with names equal to group labels.

## Value

A modified `VISTA` object with updated group colors.
