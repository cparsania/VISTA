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
  facet_nrow = NULL,
  facet_ncol = NULL,
  point_size = 6,
  linewidth = 1.2,
  label = TRUE,
  label_digits = 1,
  display_id = NULL,
  display_from = NULL,
  display_orgdb = NULL,
  facet_scales = facet_scale,
  max_genes = 15,
  line_size = NULL
)
```

## Arguments

- x:

  A `VISTA` object.

- genes:

  Character vector (\<=15 genes) to plot.

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

  Deprecated; use `facet_scales`.

- facet_nrow, facet_ncol:

  Optional layout passed to `facet_wrap()` when faceting.

- point_size:

  Numeric size of the dots.

- linewidth:

  Numeric width of the stems.

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

- facet_scales:

  Scaling option passed to `facet_wrap()` when plotting multiple genes.

- max_genes:

  Maximum number of genes accepted in one call (default 15). Previously
  an undocumented hard cap.

- line_size:

  Deprecated; use `linewidth`.

## Value

A `ggplot2` object drawing expression as lollipops.

A `ggplot2` object.

## Examples

``` r
v <- example_vista()
genes <- head(rownames(v), 5)
p <- get_expression_lollipop(v, genes = genes)
print(p)

data("count_data", package = "VISTA")
data("sample_metadata", package = "VISTA")

vista <- create_vista(
  counts = count_data[seq_len(200), ],
  sample_info = sample_metadata[seq_len(6), ],
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

genes <- rownames(vista)[seq_len(3)]
get_expression_lollipop(vista, genes = genes)

get_expression_lollipop(vista, genes = genes[seq_len(2)], by = "sample")

```
