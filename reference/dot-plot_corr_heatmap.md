# Plot a correlation heatmap from a matrix

Plot a correlation heatmap from a matrix

## Usage

``` r
.plot_corr_heatmap(
  cor_mat,
  vis_method = "square",
  plot_type = "full",
  show_diagonal = TRUE,
  show_corr_values = TRUE,
  col_corr_values = "black",
  size_corr_values = 4,
  cluster_samples = TRUE,
  scale_range = NULL
)
```

## Arguments

- cor_mat:

  A symmetric correlation matrix.

- vis_method:

  Type of shape for visualization: "square" or "circle".

- plot_type:

  Type of plot: "full", "lower", or "upper".

- show_diagonal:

  Logical; whether to show diagonal values.

- show_corr_values:

  Logical; whether to label correlation values.

- col_corr_values:

  Color for text labels.

- size_corr_values:

  Numeric size of text labels.

- cluster_samples:

  Logical; whether to cluster samples hierarchically.

- scale_range:

  Optional numeric vector of length 2 to fix the color scale.

## Value

A ggplot2 heatmap.
