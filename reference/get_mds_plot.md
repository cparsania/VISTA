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
  label_replicates = FALSE,
  label_size = 3,
  circle_size = 10,
  sample_colors = TRUE,
  shape_by = NULL,
  shape_values = NULL
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

- label_replicates:

  Logical; draw sample labels when `TRUE`.

- label_size:

  Numeric size of replicate labels when `label_replicates = TRUE`.

- circle_size:

  Numeric size for points.

- sample_colors:

  Logical; apply stored group colors when `TRUE`.

- shape_by:

  Optional column name in `sample_info` used to map point shape. When
  `NULL`, shapes are not mapped.

- shape_values:

  Optional vector of shapes passed to `scale_shape_manual()` when
  `shape_by` is set. Use a named vector to map shapes to specific
  levels.

## Value

An object returned by this function.

## Examples

``` r
NULL
#> NULL
```
