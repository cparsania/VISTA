# Raincloud plot of fold-change distributions across comparisons

Uses
[`ggrain::geom_rain()`](https://rdrr.io/pkg/ggrain/man/geom_rain.html)
to display log2 fold-change distributions for selected comparisons, with
optional statistical testing across comparisons.

## Usage

``` r
get_foldchange_raincloud(
  x,
  genes = NULL,
  sample_comparisons = NULL,
  facet_by = c("auto", "comparison", "none"),
  rain_side = c("r", "l", "f", "f1x1", "f2x2"),
  id.long.var = NULL,
  alpha = 0.5,
  point_size = 1.5,
  p.label = "p.signif",
  stats_group = FALSE,
  stats_method = "t.test",
  label = FALSE,
  label_column = "gene_id",
  label_size = 3,
  label_max_overlaps = 50,
  display_id = NULL,
  display_from = NULL,
  display_orgdb = NULL
)
```

## Arguments

- x:

  A `VISTA` object containing differential expression results.

- genes:

  Optional character vector of gene IDs to include.

- sample_comparisons:

  Optional character vector of comparison names to plot.

- facet_by:

  Faceting mode: `"auto"` (default), `"comparison"`, or `"none"`.

- rain_side:

  Side specification passed to
  [`ggrain::geom_rain()`](https://rdrr.io/pkg/ggrain/man/geom_rain.html);
  one of `"r"`, `"l"`, `"f"`, `"f1x1"`, or `"f2x2"`.

- id.long.var:

  Optional column name passed to
  [`ggrain::geom_rain()`](https://rdrr.io/pkg/ggrain/man/geom_rain.html)
  as `id.long.var` to identify repeated measurements.

- alpha:

  Alpha for jittered points.

- point_size:

  Point size for jittered points.

- p.label:

  Label type passed to
  [`ggpubr::stat_compare_means()`](https://rdrr.io/pkg/ggpubr/man/stat_compare_means.html).

- stats_group:

  Logical; add pairwise statistical tests when `TRUE`.

- stats_method:

  Statistical method passed to
  [`ggpubr::stat_compare_means()`](https://rdrr.io/pkg/ggpubr/man/stat_compare_means.html).

- label:

  Logical; add text labels to points using `ggrepel`.

- label_column:

  Column name in the plotting data used for labels. Defaults to
  `"gene_id"` for fold-change raincloud plots.

- label_size:

  Text size for point labels.

- label_max_overlaps:

  Maximum overlaps passed to
  [`ggrepel::geom_text_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html).

- display_id:

  Optional ID/column name to use for labels. If supplied and present in
  `rowData(x)`, those values are used; otherwise falls back to ID
  mapping.

- display_from:

  Optional source ID type for mapping (used when `display_id` is not
  found in `rowData`).

- display_orgdb:

  Optional `OrgDb` object used for ID mapping when `display_id` is set
  but not found in `rowData`.

## Value

A `ggplot2` object.

## Details

`id.long.var` controls which repeated unit is connected by lines in
[`ggrain::geom_rain()`](https://rdrr.io/pkg/ggrain/man/geom_rain.html).

Recommended usage for fold-change raincloud plots:

- `id.long.var = NULL` (default): best for clean distribution summaries.

- `id.long.var = "gene_id"`: most useful option; connects each gene
  across comparisons.

- `id.long.var = "comparison"` is generally not useful because
  comparison is already on the x-axis.

- Continuous value columns (e.g. `log2FoldChange`) are not suitable
  identifiers for line connections.

- Point labels (`label = TRUE`) work best with `facet_by = "none"`
  unless only a small set of genes is shown.

For identifier display consistency with other VISTA plotting functions,
set `display_id` (for example, `"SYMBOL"`). When provided, `genes` can
be given in that ID space, and default point labels use the mapped
display IDs.

## Examples

``` r
v <- example_vista()
comp <- names(comparisons(v))[1]
genes <- head(as.character(comparisons(v)[[comp]]$gene_id), 20)
p <- get_foldchange_raincloud(v, sample_comparison = comp, genes = genes)
print(p)
```
