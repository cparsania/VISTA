# Plot gene expression means with optional statistics

Summarizes expression per group for a handful of genes via barplots with
optional ggpubr comparisons.

## Usage

``` r
get_expression_barplot(
  x,
  genes,
  sample_group = NULL,
  group_column = NULL,
  log_transform = TRUE,
  stats_group = FALSE,
  facet_scale = "free_y",
  facet_scales = facet_scale,
  p.label = "p.signif",
  comparisons = NULL,
  display_id = NULL,
  display_from = NULL,
  display_orgdb = NULL
)
```

## Arguments

- x:

  A `VISTA` object.

- genes:

  Character vector (≤10 genes) to plot.

- sample_group:

  Optional character vector of groups (from `group_column`) to include.

- group_column:

  Optional column name in `sample_info` to use for grouping samples.

- log_transform:

  Logical; log2-transform expression before plotting.

- stats_group:

  Logical; add statistical comparisons between groups when `TRUE`.

- facet_scale:

  Scaling option passed to `facet_wrap()` (deprecated; use
  `facet_scales`).

- facet_scales:

  Facet scales argument passed to `facet_wrap()` when faceting by gene.

- p.label:

  Label format for
  [`ggpubr::stat_compare_means()`](https://rpkgs.datanovia.com/ggpubr/reference/stat_compare_means.html).

- comparisons:

  Optional list of specific group comparisons for
  `stat_compare_means()`.

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

## Examples

``` r
# Create VISTA object
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

# Plot expression for select genes
genes <- rownames(vista)[1:3]
get_expression_barplot(vista, genes = genes)


# With statistics
get_expression_barplot(vista, genes = genes, stats_group = TRUE)


# Without log transformation
get_expression_barplot(vista, genes = genes, log_transform = FALSE)

```
