# VISTA S4 Class Definition

The `VISTA` class extends `SummarizedExperiment`. All analysis metadata
(DE results, summaries, cutoffs, grouping info, etc.) is stored in the
`metadata()` slot of the object.

## Details

Core elements stored in `metadata(v)`:

- `$de_results`: named `SimpleList` of DE result tables.

- `$de_summary`: named `SimpleList` of DEG summary tables.

- `$de_cutoffs`: named list of thresholds (log2FC, p-value type/cutoff,
  method settings).

- `$group`: list with `column`, `palette`, `colors` describing the
  grouping/fill scheme.

- `$provenance`: list with constructor version, timestamp, session info.

Feature-level annotations live in `rowData(v)`, sample metadata in
`colData(v)`, and normalized counts (and any additional assays) in
`assay(v, "norm_counts")` by default. Use
[`create_vista()`](create_vista.md) for standard construction from
counts + sample_info; use [`vista()`](vista.md) for low-level
instantiation when you already have normalized assays/metadata.
