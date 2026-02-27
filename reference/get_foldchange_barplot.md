# Plot fold-change barplots across comparisons for selected genes

Plot fold-change barplots across comparisons for selected genes

## Usage

``` r
get_foldchange_barplot(
  x,
  genes,
  sample_comparisons = NULL,
  facet = TRUE,
  coord_flip = FALSE,
  display_id = NULL,
  sort_by = c("input", "log2fc", "abs_log2fc"),
  facet_comparison = FALSE,
  facet_scales = "free_y"
)
```

## Arguments

- x:

  A `VISTA` object containing differential expression results.

- genes:

  Character vector of gene IDs to plot.

- sample_comparisons:

  Optional character vector of comparison names to include; defaults to
  all available.

- facet:

  Logical; facet by gene when `TRUE`.

- coord_flip:

  Logical; flip axes when `TRUE`.

- display_id:

  Optional column in `rowData(x)` to use for gene labels. Input gene
  matching still uses `gene_id`.

- sort_by:

  How to order genes when faceting: `"input"` (use supplied order),
  `"log2fc"` (descending log2FC of the first comparison), or
  `"abs_log2fc"` (descending max absolute log2FC across comparisons).

- facet_comparison:

  Logical; facet by comparison (x = gene) instead of faceting by gene (x
  = comparison).

- facet_scales:

  Facet scales argument passed to `facet_wrap()` when faceting (default
  `"free_y"`).

## Value

An object returned by this function.

## Examples

``` r
NULL
#> NULL
```
