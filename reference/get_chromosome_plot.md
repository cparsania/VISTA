# Plot gene positions along chromosomes using a TxDb

Retrieves gene coordinates on the fly from a user-supplied TxDb and
plots selected genes along chromosomes, optionally colouring by a
numeric value and labeling the most variable genes.

## Usage

``` r
get_chromosome_plot(
  x,
  txdb,
  keytype = "GENEID",
  id_column = NULL,
  genes = NULL,
  value_column = NULL,
  comparison = NULL,
  group_value = NULL,
  label_top_n = 20,
  display_id = NULL,
  display_from = NULL,
  display_orgdb = NULL,
  line_length = 0.02,
  line_width = 0.6,
  filter_chrom = NULL,
  value_label = NULL,
  use_data_range = FALSE,
  force_fc_limits = FALSE,
  scale_mode = c("diverging", "sequential"),
  scale_low = NULL,
  scale_mid = NULL,
  scale_high = NULL,
  scale_midpoint = 0
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

- value_column:

  Optional column in `rowData(x)` used for colouring.

- comparison:

  Optional comparison name; when supplied, uses `log2fc` from
  `metadata(x)$de_results[[comparison]]` for colouring. If multiple
  comparisons are provided, one panel per comparison is shown (log2FC
  clipped to +/-2).

- group_value:

  Optional group label (from `group_column`); when supplied, uses mean
  expression for that group for colouring (assay `norm_counts`).

- label_top_n:

  Integer; number of genes with largest \|value\| (or random if no
  value) to label. Ignored when `genes` is provided. Set to 0 to disable
  labels.

- display_id:

  Optional column in `rowData(x)` to use for point labels (fallback to
  gene_id/rownames).

- display_from:

  Optional source ID type for mapping when `display_id` is not present
  in `rowData(x)`.

- display_orgdb:

  Optional `OrgDb` object used for identifier mapping when `display_id`
  is not present in `rowData(x)`.

- line_length:

  Horizontal half-length (in megabases) of the tick used to mark each
  gene position. Default `0.02`. Increase for longer ticks.

- line_width:

  Line width of the tick marks. Default `0.6`.

- filter_chrom:

  Optional character vector of chromosomes to keep (e.g.,
  `c("chr1","chr2")`). When `NULL`, all chromosomes returned by the TxDb
  are shown.

- value_label:

  Optional legend title override for the colour scale.

- scale_mode:

  Colour scale mode: `"diverging"` (default) or `"sequential"`.

- scale_low:

  Optional low-end colour override.

- scale_mid:

  Optional midpoint colour override (used with diverging scale).

- scale_high:

  Optional high-end colour override.

- scale_midpoint:

  Numeric midpoint passed to
  [`ggplot2::scale_color_gradient2()`](https://ggplot2.tidyverse.org/reference/scale_gradient.html)
  when `scale_mode = "diverging"`.

## Value

A `ggplot2` object.

## Details

The `genes` argument is only used as an explicit label set (all genes
are still plotted). Values in `genes` must match either the rownames of
`x` or the values in `id_column` when that is supplied. For example, if
`id_column = "ENTREZID"`, then `genes` should contain Entrez IDs to be
labeled. When multiple comparisons are supplied, `value_column` and
`group_value` are ignored and the plot is facetted by comparison with a
fixed +/-2 log2FC scale.
