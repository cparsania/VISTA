# Plot expression as a lollipop chart

Displays selected genes as stem-and-dot (lollipop) plots. By default,
expression is summarized per group (mean across replicates). With
`by = "sample"`, the plot shows individual samples.

## Usage

``` r
get_expression_lollipop(
  x,
  genes,
  sample_group = NULL,
  group_column = NULL,
  by = c("group", "sample"),
  sample_order = c("input", "group", "expression"),
  facet_by = c("auto", "gene", "none"),
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

- by:

  One of `"group"` (default; summarize replicates by group) or
  `"sample"` (show individual samples).

- sample_order:

  Ordering used when `by = "sample"`: `"input"` preserves the current
  sample order, `"group"` groups samples by `group_column`, and
  `"expression"` ranks samples by mean expression across the selected
  genes.

- facet_by:

  Faceting mode: `"auto"` (default; facet by gene when more than one
  gene is requested), `"gene"`, or `"none"`. For multiple genes,
  `"none"` falls back to `"gene"` to avoid unreadable combined panels.

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

## Value

An object returned by this function.

A `ggplot2` object.

## Examples

``` r
NULL
#> NULL
data("count_data", package = "VISTA")
data("sample_metadata", package = "VISTA")

vista <- create_vista(
  counts = count_data[1:200, ],
  sample_info = sample_metadata[1:6, ],
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control"
)
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing

genes <- rownames(vista)[1:3]
get_expression_lollipop(vista, genes = genes)

get_expression_lollipop(vista, genes = genes[1:2], by = "sample")

```
