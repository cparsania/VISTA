# Draw a sample correlation heatmap

Plots sample-sample correlation matrix derived from normalized counts
with optional clustering and annotations.

## Usage

``` r
get_corr_heatmap(
  x,
  sample_group = NULL,
  group_column = NULL,
  genes = NULL,
  corr_method = "pearson",
  vis_method = "square",
  plot_type = "full",
  cluster_samples = TRUE,
  show_diagonal = TRUE,
  show_corr_values = TRUE,
  col_corr_values = "black",
  size_corr_values = 4,
  scale_range = NULL,
  viridis_option = "viridis",
  viridis_direction = 1,
  viridis_begin = 0,
  viridis_end = 1
)
```

## Arguments

- x:

  A `VISTA` object.

- sample_group:

  Optional character vector of groups (referencing `group_column`) to
  include.

- group_column:

  Optional column name in `sample_info` defining the grouping used for
  filtering.

- genes:

  Optional character vector of gene IDs to limit the matrix.

- corr_method:

  Correlation method passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) (e.g.,
  `"pearson"`).

- vis_method:

  Currently unused; retained for backward compatibility.

- plot_type:

  Either `"full"`, `"lower"`, or `"upper"` to control which triangle is
  drawn.

- cluster_samples:

  Logical; hierarchically cluster samples before plotting when `TRUE`.

- show_diagonal:

  Logical; include the correlation diagonal when `TRUE`.

- show_corr_values:

  Logical; overlay correlation coefficients as text.

- col_corr_values:

  Color for the text labels.

- size_corr_values:

  Numeric text size multiplier.

- scale_range:

  Optional numeric vector of length two giving limits for the color
  scale.

- viridis_option:

  Character viridis palette name.

- viridis_direction:

  Integer (1 or -1) controlling palette direction.

- viridis_begin, viridis_end:

  Palette endpoints between 0 and 1.

## Value

An object returned by this function.

## Examples

``` r
NULL
#> NULL
```
