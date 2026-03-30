# Chromosome plot for expression

Convenience wrapper around
[`get_chromosome_plot()`](get_chromosome_plot.md) for expression-based
colouring (optional group mean, rowData columns, or assay columns).

## Usage

``` r
get_expression_chromosome_plot(
  x,
  txdb,
  keytype = "GENEID",
  id_column = NULL,
  genes = NULL,
  sample_group = NULL,
  group_column = NULL,
  value_column = NULL,
  value_from = c("rowdata", "assay"),
  value_assay = "norm_counts",
  facet_value_columns = FALSE,
  side_by_side_groups = FALSE,
  paginate = TRUE,
  panels_per_page = 24,
  summarise_replicates = FALSE,
  summarise_method = c("mean", "median"),
  group_value = NULL,
  label_n = 20,
  label_genes = c("top", "all", "none"),
  display_id = NULL,
  line_length = 0.02,
  line_width = 0.6,
  filter_chrom = NULL,
  log_transform = TRUE,
  value_label = "log2(mean expr)"
)
```

## Arguments

- x:

  A `VISTA` object.

- txdb:

  A TxDb object (e.g., from GenomicFeatures).

- keytype:

  Key type in the TxDb matching `id_column` (default `"GENEID"`).

- id_column:

  Optional column in `rowData(x)` used to match to TxDb keys. When
  `NULL`, rownames(x) are used as keys.

- genes:

  Optional character vector of gene IDs to label (alternative to
  `label_top_n`). When provided, all genes are plotted but only these
  are labeled. Defaults to `NULL` (no explicit label set).

- sample_group:

  Optional character vector of groups (from `group_column`) used to
  subset assay columns when `value_from = "assay"`.

- group_column:

  Optional column name in `sample_info` used for `sample_group`
  selection and replicate summarisation.

- value_column:

  Optional column in `rowData(x)` used for colouring.

- value_from:

  Source for `value_column` data: `"rowdata"` (default) or `"assay"`.
  When `"assay"`, the selected assay column is copied into `rowData`
  temporarily for colouring.

- value_assay:

  Assay name to pull values from when `value_from = "assay"`. Default
  `"norm_counts"`.

- facet_value_columns:

  Ignored (kept for compatibility); multiple `value_column`s are always
  arranged in a chromosome-by-column grid.

- side_by_side_groups:

  Logical; when plotting multiple `value_column`s, arrange one row per
  chromosome with one panel per selected sample/group. This is useful
  for direct side-by-side group comparison.

- paginate:

  Logical; when `TRUE`, split large multi-panel chromosome plots into
  pages to avoid patchwork viewport-size errors.

- panels_per_page:

  Maximum number of panels per page when `paginate = TRUE`. Set to
  `NULL` to disable automatic paging.

- summarise_replicates:

  Logical; when `TRUE` and `value_from = "assay"`, aggregate replicates
  by `group_column` before plotting.

- summarise_method:

  `"mean"` or `"median"` summarisation used when
  `summarise_replicates = TRUE`.

- group_value:

  Optional group label (from `group_column`); when supplied, uses mean
  expression for that group for colouring (assay `norm_counts`).

- label_n:

  Integer; number of genes with the largest absolute values to label
  when `genes` is not supplied. Set to 0 to disable automatic labels.

- label_genes:

  One of `"top"` (default), `"all"`, or `"none"` controlling chromosome
  text labels when `genes` is not supplied.

- display_id:

  Optional column in `rowData(x)` to use for point labels (fallback to
  gene_id/rownames).

- line_length:

  Horizontal half-length (in megabases) of the tick used to mark each
  gene position. Default `0.02`. Increase for longer ticks.

- line_width:

  Line width of the tick marks. Default `0.6`.

- filter_chrom:

  Optional character vector of chromosomes to keep (e.g.,
  `c("chr1","chr2")`). When `NULL`, all chromosomes returned by the TxDb
  are shown.

- log_transform:

  Logical; when `group_value` is used (and `value_column` is `NULL`) or
  when `value_from = "assay"`, apply log2(x + 1) before coloring.

- value_label:

  Optional legend title override for the colour scale.

## Value

A `ggplot2` object or a list of `ggplot2` objects when multiple
`value_column`s are provided and patchwork is unavailable.

## Details

- If multiple `value_column`s are supplied, one plot is produced per
  column with its own colour legend titled by that column. Chromosomes
  are laid out top-down in a single column: for each chromosome, plots
  for all `value_column`s appear sequentially. Legends for each column
  are shown only on the first occurrence (at the bottom) to avoid
  overlap. Requires patchwork; otherwise a list of plots is returned.

- Set `side_by_side_groups = TRUE` to place selected groups/samples in
  separate columns for each chromosome (same-row comparison).

- `paginate = TRUE` automatically splits large panel collections into a
  list of patchwork pages, which prevents "viewport too small" rendering
  errors in IDE plot panes.

- For `value_from = "assay"`, colours follow the VISTA ecosystem: group
  colours are used when `summarise_replicates = TRUE`, while sample
  colours are derived from group colours when plotting individual
  replicates.

- For `value_from = "assay"` with `value_column = NULL`, all selected
  assay columns are used (samples after `sample_group` filtering, or
  groups after replicate summarisation).

- Labels are kept consistent across `value_column`s: if `genes` is
  provided it is used for all panels; otherwise the top `label_n` (by
  absolute value in the first `value_column`) are used for all panels.

- When `value_from = "assay"`, the specified assay column is copied into
  `rowData` on the fly so it can be used for colouring.

- When `group_value` is provided (and no `value_column`), colouring is
  based on log2 mean expression for that group (assay `norm_counts`),
  using a data-driven colour scale. Ignored when `value_column` is
  supplied.

## Examples

``` r
v <- example_vista()
p <- try(get_expression_chromosome_plot(v), silent = TRUE)
if (!inherits(p, 'try-error')) print(p)
```
