# Raincloud plot of expression values

Uses
[`ggrain::geom_rain()`](https://rdrr.io/pkg/ggrain/man/geom_rain.html)
to combine a half-violin, boxplot, and jittered points per sample/group
to show distribution, summary, and individual values.

## Usage

``` r
get_expression_raincloud(
  x,
  genes = NULL,
  sample_group = NULL,
  group_column = NULL,
  by = "group",
  value_transform = c("log2", "zscore", "none"),
  summarise = FALSE,
  facet_by = c("auto", "gene", "none"),
  sample_order = c("input", "group", "expression"),
  rain_side = c("r", "l", "f", "f1x1", "f2x2"),
  id.long.var = NULL,
  alpha = 0.5,
  point_size = 1.5,
  p.label = "p.signif",
  stats_group = FALSE,
  stats_method = "t.test",
  label = FALSE,
  label_column = "gene",
  label_size = 3,
  label_max_overlaps = 50,
  display_id = NULL,
  display_from = NULL,
  display_orgdb = NULL
)
```

## Arguments

- x:

  A `VISTA` object.

- genes:

  Optional character vector of gene IDs to include; defaults to all.

- sample_group:

  Optional subset of groups (values of `group_column`) to keep.

- group_column:

  Grouping column in `sample_info`; defaults to the stored grouping.

- by:

  Plot unit. Violin plots currently support only `"group"`, because a
  violin needs replicate-level distributions within groups rather than
  single values per sample.

- value_transform:

  One of `"log2"`, `"zscore"`, or `"none"`.

- summarise:

  Logical; when `TRUE`, collapse replicates within each group for each
  gene (one value per gene per group). This is useful for pooled
  multi-gene raincloud plots where each dot represents one gene-level
  summary in a group.

- facet_by:

  Faceting mode: `"auto"` (default), `"gene"`, or `"none"`.

- sample_order:

  Ordering for sample-level x-axis display: `"input"`, `"group"`, or
  `"expression"` before values are grouped into violins.

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
  [`ggpubr::stat_compare_means()`](https://rpkgs.datanovia.com/ggpubr/reference/stat_compare_means.html).

- stats_group:

  Logical; add pairwise statistical tests when `TRUE`.

- stats_method:

  Statistical method passed to
  [`ggpubr::stat_compare_means()`](https://rpkgs.datanovia.com/ggpubr/reference/stat_compare_means.html).

- label:

  Logical; add text labels to points using `ggrepel`.

- label_column:

  Column name in the plotting data used for labels. Defaults to `"gene"`
  for expression raincloud plots.

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

Recommended usage for expression raincloud plots:

- `id.long.var = NULL` (default): best for clean distribution summaries.

- `id.long.var = "gene"`: best when plotting a small number of genes and
  showing gene-level trajectories across x levels.

- `id.long.var = "<subject_id_column>"`: best for
  paired/repeated-measure designs when a subject ID exists in
  `sample_info`.

- `id.long.var = "sample"` or the grouping variable is usually less
  informative and can over-connect points.

- Point labels (`label = TRUE`) work best with `facet_by = "none"` or a
  small number of genes.

For identifier display consistency with other VISTA plotting functions,
set `display_id` (for example, `"SYMBOL"`). When provided, `genes` can
be given in that ID space, and default point labels use the mapped
display IDs.

## Examples

``` r
v <- example_vista()
genes <- head(rownames(v), 5)
p <- get_expression_raincloud(v, genes = genes, summarise = TRUE)
#> Warning: `summarise = TRUE` with `facet_by = 'gene'` gives one summarized value per
#> group in each facet, so raincloud distributions are not informative. Consider
#> `facet_by = 'none'` with `id.long.var = 'gene'`.
print(p)
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: Computation failed in `stat_half_ydensity()`.
#> Caused by error in `$<-.data.frame`:
#> ! replacement has 1 row, data has 0
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: Computation failed in `stat_half_ydensity()`.
#> Caused by error in `$<-.data.frame`:
#> ! replacement has 1 row, data has 0
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: Computation failed in `stat_half_ydensity()`.
#> Caused by error in `$<-.data.frame`:
#> ! replacement has 1 row, data has 0
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: Computation failed in `stat_half_ydensity()`.
#> Caused by error in `$<-.data.frame`:
#> ! replacement has 1 row, data has 0
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: Computation failed in `stat_half_ydensity()`.
#> Caused by error in `$<-.data.frame`:
#> ! replacement has 1 row, data has 0
```
