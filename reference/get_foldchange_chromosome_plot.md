# Chromosome plot for fold change

Convenience wrapper around
[`get_chromosome_plot()`](get_chromosome_plot.md) for fold-change
colouring. When multiple comparisons are supplied, panels are facetted
by comparison with log2FC clipped to ±2.

## Usage

``` r
get_foldchange_chromosome_plot(
  x,
  txdb,
  keytype = "GENEID",
  id_column = NULL,
  genes = NULL,
  comparison = NULL,
  value_column = NULL,
  label_top_n = 20,
  display_id = NULL,
  line_length = 0.02,
  line_width = 0.6,
  filter_chrom = NULL
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

- comparison:

  Optional comparison name; when supplied, uses `log2fc` from
  `metadata(x)$de_results[[comparison]]` for colouring. If multiple
  comparisons are provided, one panel per comparison is shown (log2FC
  clipped to ±2).

- value_column:

  Optional column in `rowData(x)` used for colouring.

- label_top_n:

  Integer; number of genes with largest \|value\| (or random if no
  value) to label. Ignored when `genes` is provided. Set to 0 to disable
  labels.

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

## Value

A `ggplot2` object.
