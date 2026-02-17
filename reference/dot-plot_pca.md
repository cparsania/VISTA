# Plot PCA results using ggplot2

Internal function that creates a PCA plot using ggplot2 with optional
coloring, labeling, clustering ellipses, and customizable color
palettes.

## Usage

``` r
.plot_pca(
  pca_df,
  group_col,
  circle_size,
  label_replicates,
  sample_colors,
  show_clusters,
  color_vals,
  pca
)
```

## Arguments

- pca_df:

  A data frame of PCA coordinates and metadata (from
  [`.prepare_pca_dataframe()`](dot-prepare_pca_dataframe.md)).

- group_col:

  Column name used for grouping/coloring samples.

- circle_size:

  Point size for the scatter plot.

- label_replicates:

  Logical; whether to display sample labels.

- sample_colors:

  Logical; whether to color by group.

- show_clusters:

  Logical; whether to draw ellipses for group clusters.

- color_vals:

  Named color vector for groups.

- pca:

  Original PCA object used to compute variance explained for axis
  labels.

## Value

A ggplot2 object.
