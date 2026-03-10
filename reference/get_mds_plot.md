# Generate an MDS plot for samples in a VISTA object

Runs classical multidimensional scaling on normalized counts, optionally
restricting to groups or genes.

## Usage

``` r
get_mds_plot(
  x,
  sample_group = NULL,
  group_column = NULL,
  genes = NULL,
  top_n_genes = NULL,
  label = FALSE,
  label_size = 3,
  point_size = 10,
  shape_by = NULL,
  shape_values = NULL,
  color_by = NULL,
  use_vista_colors = TRUE,
  palette = NULL,
  colors = NULL
)
```

## Arguments

- x:

  A `VISTA` object.

- sample_group:

  Optional character vector of groups to include (based on the column
  specified by `group_column`).

- group_column:

  Optional column name in `sample_info` to use for grouping/filtering.

- genes:

  Optional character vector of gene identifiers to restrict the matrix.

- top_n_genes:

  Optional integer selecting the top variable genes to include.

- label:

  Logical; draw sample labels when `TRUE`.

- label_size:

  Numeric size of sample labels when `label = TRUE`.

- point_size:

  Numeric size for points.

- shape_by:

  Optional column name in `sample_info` used to map point shape. When
  `NULL`, shapes are not mapped.

- shape_values:

  Optional vector of shapes passed to `scale_shape_manual()` when
  `shape_by` is set. Use a named vector to map shapes to specific
  levels.

- color_by:

  Optional column name in `sample_info` used for point colour. Defaults
  to the active grouping column.

- use_vista_colors:

  Logical; when `TRUE`, prefer the stored VISTA group colours when
  colouring by the grouping column.

- palette:

  Optional qualitative palette name used when generating colours for
  non-group metadata levels.

- colors:

  Optional named character vector of manual colours overriding both
  `palette` and stored VISTA colours.

## Value

An object returned by this function.

## Examples

``` r
v <- example_vista()
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
p <- get_mds_plot(v)
print(p)
```
