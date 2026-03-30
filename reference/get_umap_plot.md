# Generate a UMAP plot for samples in a VISTA object

Runs UMAP on normalized counts, optionally restricting to selected
groups or genes. UMAP is intended for exploratory sample-level
structure.

## Usage

``` r
get_umap_plot(
  x,
  sample_group = NULL,
  group_column = NULL,
  color_by = NULL,
  genes = NULL,
  top_n_genes = NULL,
  label = FALSE,
  label_size = 3,
  point_size = 10,
  shape_by = NULL,
  shape_values = NULL,
  n_neighbors = 15,
  min_dist = 0.1,
  metric = "euclidean",
  seed = 123,
  use_vista_colors = NULL,
  palette = NULL,
  colors = NULL,
  use_group_colors = TRUE
)
```

## Arguments

- x:

  A `VISTA` object.

- sample_group:

  Optional character vector of groups to include (based on
  `group_column`).

- group_column:

  Optional column name in `sample_info` used for filtering/grouping.
  Defaults to the stored grouping column.

- color_by:

  Optional column name in `sample_info` used for point color. Defaults
  to `group_column`.

- genes:

  Optional character vector of gene identifiers to restrict the matrix.

- top_n_genes:

  Optional integer selecting top variable genes to include.

- label:

  Logical; draw sample labels when `TRUE`.

- label_size:

  Numeric label size when `label = TRUE`.

- point_size:

  Numeric point size.

- shape_by:

  Optional column name in `sample_info` used to map point shape.

- shape_values:

  Optional vector passed to `scale_shape_manual()` when `shape_by` is
  set.

- n_neighbors:

  UMAP `n_neighbors` parameter.

- min_dist:

  UMAP `min_dist` parameter.

- metric:

  UMAP distance metric.

- seed:

  Integer random seed passed to UMAP.

- use_vista_colors:

  Deprecated alias for `use_group_colors`. When supplied, it overrides
  `use_group_colors`.

- palette:

  Optional qualitative palette name used when generating colours for
  non-group metadata levels.

- colors:

  Optional named character vector of manual colours overriding both
  `palette` and stored VISTA colours.

- use_group_colors:

  Logical; when `TRUE`, prefer the stored VISTA group colours when
  colouring by the grouping column.

## Value

A ggplot object with UMAP1/UMAP2 coordinates.

## Examples

``` r
if (requireNamespace("uwot", quietly = TRUE)) {
  vista <- example_vista()
  get_umap_plot(vista, top_n_genes = 50)
}
#> Warning: `n_neighbors` (15) must be smaller than sample size (6); using 5.
```
