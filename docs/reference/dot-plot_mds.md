# Plot MDS Coordinates

Internal plotting function for rendering MDS results using `ggplot2`.

## Usage

``` r
.plot_mds(
  mds_df,
  group_col,
  circle_size,
  label_replicates,
  sample_colors,
  color_vals
)
```

## Arguments

- mds_df:

  A data frame with MDS coordinates and sample metadata.

- group_col:

  The grouping column used for coloring points.

- circle_size:

  Size of the points.

- label_replicates:

  Logical; whether to label each point with the sample name.

- sample_colors:

  Logical; whether to color points by group.

- color_vals:

  Named color vector for groups.

## Value

A `ggplot` object.
