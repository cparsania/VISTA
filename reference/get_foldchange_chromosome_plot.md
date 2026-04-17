# Chromosome plot for fold change

Convenience wrapper around
[`get_chromosome_plot()`](get_chromosome_plot.md) for fold-change
colouring. When multiple comparisons are supplied, panels are facetted
by comparison with log2FC clipped to +/-2.

## Usage

``` r
get_foldchange_chromosome_plot(
  x,
  txdb,
  keytype = "GENEID",
  id_column = NULL,
  genes = NULL,
  sample_comparison = NULL,
  value_column = NULL,
  label_n = 20,
  display_id = NULL,
  display_from = NULL,
  display_orgdb = NULL,
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

- sample_comparison:

  Optional comparison name (or vector of names) used for fold-change
  colouring.

- value_column:

  Optional column in `rowData(x)` used for colouring.

- label_n:

  Integer; number of genes with the largest absolute fold-change to
  label when `genes` is not supplied. Set to 0 to disable labels.

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

## Examples

``` r
v <- example_vista()
p <- try(get_foldchange_chromosome_plot(v, sample_comparison = names(comparisons(v))[1]), silent = TRUE)
if (!inherits(p, 'try-error')) print(p)
```
