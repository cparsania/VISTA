# Expression heatmap

Summarizes expression for selected genes/groups via ComplexHeatmap with
optional transformations and annotations.

## Usage

``` r
get_expression_heatmap(
  vista_obj,
  samples,
  genes,
  value_transform = c("zscore", "log2", "raw"),
  repair_genes = FALSE,
  color_default = TRUE,
  col = NULL,
  summarise_replicates = TRUE,
  summarise_method = c("mean", "median"),
  convert_rowmeans = FALSE,
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
  cluster_by = NULL,
  column_anno_palette = "Dark 3",
  heatmap_name = NULL,
  display_id = NULL,
  display_from = NULL,
  display_orgdb = NULL,
  group_column = NULL,
  ...
)
```

## Arguments

- vista_obj:

  A `VISTA` object.

- samples:

  Character vector of group labels specifying which samples to include
  (based on the selected grouping column).

- genes:

  Character vector of gene identifiers to display.

- value_transform:

  One of `"zscore"`, `"log2"`, or `"raw"` controlling how expression
  values are transformed.

- repair_genes:

  Logical; if `TRUE`, split `gene_id` strings such as `ID:SYMBOL` to
  display the symbol.

- color_default:

  Logical; use the default blue-white-red palette when `TRUE`. Set to
  `FALSE` to supply `col`.

- col:

  Optional
  [`circlize::colorRamp2`](https://rdrr.io/pkg/circlize/man/colorRamp2.html)
  function used when `color_default = FALSE`.

- summarise_replicates:

  Logical; average replicates per group before plotting when `TRUE`.

- summarise_method:

  `"mean"` or `"median"` summarization used when
  `summarise_replicates = TRUE`.

- convert_rowmeans:

  Logical; subtract row means prior to plotting.

- show_row_names:

  Logical; display row names (genes) beside the heatmap.

- cluster_rows:

  Logical; cluster rows when drawing the heatmap.

- show_row_dend:

  Logical; display the row dendrogram.

- row_names_font_size:

  Numeric font size for row names.

- label_specific_rows:

  Optional character vector of row names to highlight via `anno_mark()`.

- label_specific_rows_gp:

  [`grid::gpar()`](https://rdrr.io/r/grid/gpar.html) object controlling
  the highlighted row labels.

- show_column_names:

  Logical; draw column names when `TRUE`.

- cluster_columns:

  Logical; cluster columns.

- show_heatmap_legend:

  Logical; display the heatmap legend.

- kmeans_k:

  Optional integer specifying the number of k-means clusters to compute
  for rows.

- return_type:

  `"heatmap"`, `"clusters"`, or `"both"` selecting the returned object.

- annotate_columns:

  Logical or character vector. `TRUE` adds one column annotation using
  `group_column`; a character vector adds multiple annotations from
  `sample_info`.

- cluster_by:

  Optional annotation column used to split/cluster columns. Defaults to
  the first annotation column when `annotate_columns` is enabled.

- column_anno_palette:

  Qualitative palette name used for the column annotation.

- heatmap_name:

  Optional legend title.

- display_id:

  Optional ID/column name to use for row labels. If supplied

- display_from:

  Optional source ID type for mapping (used when `display_id`

- display_orgdb:

  Optional `OrgDb` object used for ID mapping when

- group_column:

  Optional column name in `sample_info` used to interpret `samples`.

- ...:

  Additional arguments passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).

## Value

An object returned by this function.

## Examples

``` r
NULL
#> NULL
```
