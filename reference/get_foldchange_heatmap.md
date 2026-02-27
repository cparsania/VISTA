# Fold-change heatmap

Visualizes log2 fold-change matrices across comparisons with
ComplexHeatmap, supporting clustering and annotations.

## Usage

``` r
get_foldchange_heatmap(
  vista_obj,
  sample_comparisons,
  genes,
  repair_genes = FALSE,
  color_default = TRUE,
  col = NULL,
  show_row_names = FALSE,
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  row_names_font_size = 10,
  label_specific_rows = NULL,
  label_specific_rows_gp = grid::gpar(fontsize = 5),
  show_column_names = TRUE,
  cluster_columns = TRUE,
  show_heatmap_legend = TRUE,
  kmeans_k = NULL,
  return_type = c("heatmap", "clusters", "both"),
  annotate_columns = FALSE,
  column_anno_palette = "Set2",
  heatmap_name = NULL,
  display_id = NULL,
  display_from = NULL,
  display_orgdb = NULL,
  ...
)
```

## Arguments

- vista_obj:

  A `VISTA` object with stored differential expression results.

- sample_comparisons:

  Character vector of comparison names to include.

- genes:

  Character vector of gene identifiers to display.

- repair_genes:

  Logical; attempt to simplify `gene_id` strings by removing prefixes.

- color_default:

  Logical; use the default diverging palette when `TRUE`. Set to `FALSE`
  to supply `col`.

- col:

  Optional
  [`circlize::colorRamp2`](https://rdrr.io/pkg/circlize/man/colorRamp2.html)
  color function used when `color_default = FALSE`.

- show_row_names:

  Logical; draw row (gene) names.

- cluster_rows:

  Logical; cluster rows.

- show_row_dend:

  Logical; display the row dendrogram.

- row_names_font_size:

  Numeric font size for row names.

- label_specific_rows:

  Optional character vector of genes to highlight with `anno_mark()`.

- label_specific_rows_gp:

  [`grid::gpar()`](https://rdrr.io/r/grid/gpar.html) object controlling
  highlighted labels.

- show_column_names:

  Logical; draw column labels.

- cluster_columns:

  Logical; cluster columns.

- show_heatmap_legend:

  Logical; display the heatmap legend.

- kmeans_k:

  Optional integer specifying the number of k-means clusters for rows.

- return_type:

  `"heatmap"`, `"clusters"`, or `"both"` selecting the returned value.

- annotate_columns:

  Logical; add annotation bars keyed to the sample grouping column.

- column_anno_palette:

  Qualitative palette name used for column annotations.

- heatmap_name:

  Optional legend title.

- display_id:

  Optional ID/column name to use for plot labels. If supplied

- display_from:

  Optional source ID type for mapping (used when `display_id`

- display_orgdb:

  Optional `OrgDb` object used for ID mapping when

- ...:

  Additional arguments forwarded to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).

## Value

An object returned by this function.

## Examples

``` r
NULL
#> NULL
```
