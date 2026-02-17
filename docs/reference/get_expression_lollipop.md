# Plot mean expression per group as a lollipop chart

Summarizes expression per group for a handful of genes using a
stem-and-dot (lollipop) plot. Values are averaged across replicates in
each group; this function does not show individual replicates. For
per-sample display or pairwise tests, use
[`get_expression_barplot()`](get_expression_barplot.md).

## Usage

``` r
get_expression_lollipop(
  x,
  genes,
  sample_group = NULL,
  group_column = NULL,
  log_transform = TRUE,
  facet_scale = "free_y",
  point_size = 6,
  line_size = 1.2,
  label = TRUE,
  label_digits = 1,
  display_id = NULL,
  display_from = NULL,
  display_orgdb = NULL
)
```

## Arguments

- x:

  A `VISTA` object.

- genes:

  Character vector (≤15 genes) to plot.

- sample_group:

  Optional character vector of groups (from `group_column`) to include.

- group_column:

  Optional column name in `sample_info` to use for grouping samples.

- log_transform:

  Logical; log2-transform expression before plotting.

- facet_scale:

  Scaling option passed to `facet_wrap()` when plotting multiple genes.

- point_size:

  Numeric size of the dots.

- line_size:

  Numeric size of the stems.

- label:

  Logical; draw numeric labels above the dots.

- label_digits:

  Integer; digits to show in labels when `label = TRUE`.

- display_id:

  Optional ID/column name to use for labels/facets. If supplied and
  present in `rowData(x)`, those values are used; otherwise falls back
  to ID mapping.

- display_from:

  Optional source ID type for mapping (used when `display_id` is not
  found in `rowData`).

- display_orgdb:

  Optional `OrgDb` object used for ID mapping when `display_id` is set
  but not found in `rowData`.
