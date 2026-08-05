# VISTA 1.1.2

Targeted at Bioconductor 3.24 (VISTA 1.2.0). This release fixes several
defects that produced silently wrong results, makes VISTA behave like a
proper Bioconductor object, and begins a deprecation cycle for a number of
inconsistent argument names. **No existing script stops working**: every
renamed argument still functions and warns.

## Data-integrity fixes

Check whether you were affected:

- **`get_expression_matrix(genes = ..., summarise = TRUE)` mislabelled rows.**
  Rows were subset in assay order but labelled in the order `genes` was
  supplied, so whenever those differed every row carried another gene's
  identifier and values.

  *Scope:* no other VISTA function called `get_expression_matrix()`, so no
  plot, heatmap, DEG table, export or generated report was affected. Internal
  group averaging goes through `norm_counts(summarise = TRUE)`, which labels
  rows from the matrix itself and was never affected by this. Only code that
  called `get_expression_matrix()` directly could be wrong.

  *Were you affected?* Only if you passed `genes` **and** `summarise = TRUE`
  **and** your gene vector was in some order other than the object's own row
  order. Passing genes already in `rownames(x)` order (or omitting `genes`)
  always produced correct output. Note that ranked inputs — for example
  `get_genes_by_regulation(top_n = ...)`, which sorts by absolute fold change —
  are normally *not* in row order, so those calls were affected. Re-run them.

  The bug dates back to the first release, so it is present in 1.0.0 and every
  0.99.x.
- **GSEA ranked vectors slid scores onto the wrong genes.** `get_gsea()`
  reassigned `names(rank_vec)` from a mapping helper that drops unmapped and
  duplicate identifiers, so the shorter vector was padded with `NA` and every
  score after the first gap moved to a different gene. The resulting NES
  values and leading-edge gene lists were meaningless.

  *Scope:* only `set_type = "kegg"` was affected in practice. That path always
  converts to ENTREZID, and any real dataset contains identifiers with no
  ENTREZ mapping, so the vector always shortened. `set_type = "msigdb"`
  converted an identifier type to *itself*, which short-circuits before
  querying the annotation database and cannot shorten. `set_type = "go"`
  never reassigned the names at all. Verified on a 22-gene input containing
  two identifiers absent from the OrgDb: `"kegg"` returned 20 names for 22
  scores, `"msigdb"` returned 22.

  No other VISTA function called `get_gsea()`, so nothing else in the package
  consumed a corrupted ranking. Re-run any `get_gsea(set_type = "kegg")`
  results; MSigDB and GO results are unaffected.
- **`get_corr_heatmap()` drew a staircase, not a triangle.** The `triangle`
  mask was computed against input order while the axes were re-levelled to
  clustered order, under the default `cluster_by = "correlation"`.
- **Multi-file count importers bound columns positionally.** The STAR, HTSeq
  and RSEM importers took gene identifiers from the first file only and never
  checked that row *i* described the same gene in the others, so files with
  equal row counts but different ordering were silently misaligned and every
  count for those samples was attributed to the wrong gene.

  *Scope:* only `read_vista_counts()` given a **vector of two or more file
  paths** with `format` of `"star"`, `"htseq"` or `"rsem"`, and only when
  those files did not already list the same genes in the same order. A single
  file, a data frame, a matrix, `"featurecounts"` and `"tximport"` all take a
  name-keyed path and were never affected. Files written by one pipeline run
  normally do share an order, so a uniform set of inputs was fine; the risk is
  mixing files produced at different times or by different pipeline versions,
  where HTSeq's trailing `__no_feature`/`__ambiguous` rows or STAR's four
  leading `N_*` rows may differ.

  Imports now match by gene identifier and **error**, naming the file and up
  to five missing identifiers, when a file does not cover the reference gene
  set. This can surface an error where a script previously appeared to work;
  it was producing wrong numbers. Anything built from a multi-file STAR,
  HTSeq or RSEM import should be re-imported and re-analysed.
- **`run_vista_report()` reported its own defaults as your parameters.** When a
  prebuilt object was supplied via `vista_rds`, the report printed
  deseq2 / 1 / 0.05 / padj regardless of how the object was built.
- **`run_cell_deconvolution()`'s SYMBOL retry told xCell2 the wrong identifier
  type.** When the first scoring attempt failed, VISTA retried with a
  SYMBOL-keyed matrix built from `rowData`, but carried the original
  `gene_id_type` (for example `"ensembl"`) into the retry -- handing xCell2
  symbols while declaring them Ensembl IDs. The retry therefore either failed
  with the real error masked behind a generic message, or scored against a
  near-empty gene overlap.
- **Deconvolution was broken for `ENSG...:SYMBOL` rownames** — the common case
  reached automatically by `gene_id_type = "auto"`.
- Fixed `get_ma_plot()` reading a non-existent metadata key for its threshold
  fallback, `get_expression_lollipop()` failing to resolve symbols through an
  OrgDb, comparison colours collapsing to `character(0)` and crashing
  `get_foldchange_lineplot()`, `NaN` columns when summarising over a factor
  with unused levels, and YAML/filename escaping in report and asset output.
- **Hardened the DE-table alignment helper**, which reindexes every comparison
  onto the counts matrix and so underpins the whole package. It now rejects a
  duplicated reference gene list instead of emitting rows whose names R had
  silently de-duplicated (`g2` -> `g2.1`) away from the `gene_id` column they
  came from, coerces tibble input to a plain data frame — a tibble ignores
  `rownames<-`, so every gene read as absent and was replaced by an `NA` row
  with a perfectly correct label — and asserts as a post-condition that the
  returned rownames and `gene_id` column both match the reference exactly.
  No shipped code path reached the tibble case; both are closed against future
  callers, because a wrong answer here is invisible in the output's shape.

## Reproducibility

- **Pathway gene capping is now deterministic.** `get_pathway_heatmap(max_genes)`
  and `get_enrichment_chord(max_genes)` truncated the pathway gene set
  positionally. Pathway membership is a *set*, and clusterProfiler's ordering of
  it is arbitrary and, verified here, differs between R sessions -- so the same
  call could plot a different subset of genes on different runs. The heatmap now
  ranks candidates by expression variance (most informative rows first) and the
  chord by pathway participation, both breaking ties on the gene identifier.
- **Deconvolution scores are no longer aligned to samples by position when
  their own labels disagree.** `run_cell_deconvolution()` accepted scores whose
  row count matched the sample count and stamped the object's sample names onto
  them, even when the scores carried different labels of their own -- which
  would attach every sample's cell fractions to the wrong sample. Alignment is
  now by name whenever names are available in either orientation; contradictory
  labels are an error; and position is used only when the scores carry no
  labels at all, with a message saying so.
- The airway vignette sets a seed. GSEA estimates p-values by permutation and
  ComplexHeatmap's row k-means is unseeded, so the rendered document previously
  changed on every build for reasons unrelated to the code.

## Analysis changes

- The edgeR and limma backends now filter with the same predicate as DESeq2
  (a gene is kept when at least `min_replicates` samples reach `min_counts`);
  they previously ignored `min_counts` and used a hardcoded `cpm > 1`. They
  also filter *before* `calcNormFactors()`, per the edgeR user's guide. edgeR
  and limma results shift slightly relative to 1.1.1.
- `min_replicates` is documented accurately: it counts samples across the
  whole experiment, never within each group.

## New features

- **Raw counts are retained** as a `counts` assay, reachable with `counts(v)`.
  `as_deseq_dataset()` converts a VISTA object back into a `DESeqDataSet`.
  Roughly doubles the assay footprint; opt out with
  `create_vista(keep_raw_counts = FALSE)`.
- **`[` method.** Subsetting now keeps metadata consistent: row subsetting
  reindexes every DE table and recounts DEG summaries, column subsetting
  prunes orphaned colour maps and warns when a comparison loses its samples.
  Previously `v[1:10, ]` left the DE tables describing all original genes and
  `validObject()` still returned `TRUE`.
- **`show()` method.** Objects now display as `class: VISTA` with their
  grouping, comparisons, active DE source, cutoffs and schema version.
- **`updateObject()`** migrates objects from older metadata schemas.
- **`?VISTA`** resolves, and there is a package overview page.
- Unknown arguments passed through `...` are now rejected with a did-you-mean
  suggestion. A typo such as `get_expression_heatmap(v, gene = my_genes)`
  previously plotted the default gene set instead of erroring.
- `get_foldchange_lineplot()` gains `return_type`. Its default stays `"both"`
  (the list it has always returned); pass `"plot"` for a bare ggplot.
- `get_enrichment_chord()` gains `return_type`. Its default stays `"data"`
  because it draws to the active device and has always returned its table
  invisibly; `"plot"` returns a recorded plot that `save_vista_plot()` can now
  write to a file. VISTA enables the device display list before drawing, since
  `png()` and `pdf()` leave it off and would otherwise yield an empty
  recording that silently saved nothing.

## Deprecations

All still work and warn once per session. See `?VISTA-deprecated` for the
full table and timelines; each becomes defunct in 1.4.0.

| Function | Deprecated | Use instead |
|---|---|---|
| `get_deg_count_pieplot()`, `get_deg_count_donutplot()` | `label` | `label_type` |
| `get_corr_heatmap()` | `cluster_by` | `order_by` |
| `get_corr_heatmap()` | `show_corr_values`, `col_corr_values` | `label`, `label_color` |
| `get_celltype_group_dotplot()` | `error` | `errorbar` |
| four `get_expression_*()` plots | `comparisons` | `stat_comparisons` |
| `get_pca_plot()`, `get_mds_plot()`, `get_umap_plot()` | `top_n_genes` | `top_n` |
| `get_expression_lollipop()`, `get_foldchange_lollipop()` | `line_size` | `linewidth` |
| `get_volcano_plot()` | `col_up`, `col_down`, `col_other(s)`, `lab_size` | `colors`, `label_size` |
| `get_pca_plot()` | `sample.seed` | removed; it never had any effect |
| `get_expression_barplot()`, `get_expression_lollipop()` | `facet_scale` | `facet_scales` |
| `get_expression_violinplot()`, `get_expression_lineplot()` | `value_transform` | `log_transform` |

`return_type` values are unified on `c("plot", "data", "both")`. The legacy
spellings still work and warn: `"heatmap"`/`"clusters"` (both heatmaps),
`"matrix"` (`get_celltype_heatmap()`), `"genes"` (`get_pathway_heatmap()`).
Arguments that select the *shape* of a purely tabular result --
`get_pathway_genes()`, `get_genes_by_regulation()`, `read_vista_counts()`,
`derive_vista_metadata()` -- keep their own vocabularies, because they have no
plot component.

Undocumented gene caps (20 for embeddings, 15 for lollipop, 25 for barplot)
are now the `max_genes` argument, with the existing values as defaults.

## Testing

- Added label/value integrity tests covering every exported function that
  accepts a gene set. They assert two properties the 1.0.0 defect violated:
  results do not depend on the order genes were supplied in, and the value
  reported for a gene is that gene's own value, checked against an independent
  recomputation from the object. The previous suite asserted shape --
  `is.matrix()`, `nrow()`, class -- all of which the broken output satisfied.
- Verified by mutation testing: reintroducing the exact 1.0.0 defect now fails
  five tests across three files, and replacing a name-keyed lookup with a
  positional one elsewhere fails nine.
- The DE-table alignment helper is now covered end to end: identifier resolution
  and its fallbacks, DESeq2/edgeR/limma column canonicalisation, type-preserving
  `NA` padding across numeric, integer, character, logical, factor, `Date` and
  `POSIXct` columns, and the reordering guarantee itself. Both faults fixed above
  were found by these tests.

## Internal

- Accessor generics dispatch only on the object rather than on every formal.
- `Depends: R (>= 4.6.0)`; `viridis` dropped from Imports (unused);
  `BiocGenerics` added. `grDevices`, `graphics`, `stats`, `tools` and `utils`
  are now declared — all five were already used via `::` without being listed.
  The maintainer's ORCID iD is recorded in `Authors@R`.
- `print.vista` removed — no object could carry that class.
- Plot layer data is pinned by snapshot tests so API renames cannot silently
  change output.

# VISTA 1.1.1

## Bug fixes

- `get_volcano_plot()` now inherits `log2fc_cutoff`, `pval_cutoff`, and the new
  `p_value_type` argument from `cutoffs(x)` instead of hard-coding the raw
  `pvalue` column at 0.05. Previously the volcano could colour a different set
  of genes than `deg_summary()`, `get_deg_count_barplot()`, and
  `get_genes_by_regulation()` reported for the same object. Pass the arguments
  explicitly to restore the old behaviour
  (`p_value_type = "pvalue"`, `pval_cutoff = 0.05`). The y-axis label now
  reflects which p-value column is plotted.

- Consensus DE tables (`method = "both"`) no longer overwrite `pvalue` and
  `padj` with `1` for every gene that is not called by both backends. Each gene
  now carries the less-significant of the two backend values, so the columns
  stay continuous and usable for volcano/MA plots and ranking. Per-backend
  values remain available in `pvalue_deseq2`/`pvalue_edger`/`padj_deseq2`/
  `padj_edger`, and DEG calls are unchanged.

- `get_ma_plot()` read its threshold fallback from the non-existent
  `metadata(x)$cutoffs` key, silently using package defaults instead of the
  object's stored cutoffs.

## Reproducibility

- **Pathway gene capping is now deterministic.** `get_pathway_heatmap(max_genes)`
  and `get_enrichment_chord(max_genes)` truncated the pathway gene set
  positionally. Pathway membership is a *set*, and clusterProfiler's ordering of
  it is arbitrary and, verified here, differs between R sessions -- so the same
  call could plot a different subset of genes on different runs. The heatmap now
  ranks candidates by expression variance (most informative rows first) and the
  chord by pathway participation, both breaking ties on the gene identifier.
- **Deconvolution scores are no longer aligned to samples by position when
  their own labels disagree.** `run_cell_deconvolution()` accepted scores whose
  row count matched the sample count and stamped the object's sample names onto
  them, even when the scores carried different labels of their own -- which
  would attach every sample's cell fractions to the wrong sample. Alignment is
  now by name whenever names are available in either orientation; contradictory
  labels are an error; and position is used only when the scores carry no
  labels at all, with a message saying so.
- The airway vignette sets a seed. GSEA estimates p-values by permutation and
  ComplexHeatmap's row k-means is unseeded, so the rendered document previously
  changed on every build for reasons unrelated to the code.

## Analysis changes

- The edgeR and limma-voom backends now filter genes with the same predicate as
  the DESeq2 backend — a gene is kept when at least `min_replicates` samples
  reach `min_counts`. They previously ignored `min_counts` at this step and used
  a hard-coded `cpm > 1` threshold, so the same arguments produced different
  feature sets across backends.

- The edgeR and limma-voom backends now filter *before* calling
  `edgeR::calcNormFactors()`, per the edgeR user's guide; normalization factors
  were previously estimated from genes that were about to be discarded. Both
  changes shift edgeR/limma results slightly relative to 1.1.0.

## Documentation

- `min_replicates` is documented accurately: it counts samples across the whole
  experiment, not within each group, which is what all backends have always
  implemented.
- `create_vista()` now documents the consensus table's column semantics,
  including the `support` column.
- Removed a reference to `get_enrichment_network()` from the README; no such
  function exists. Corrected the citation block and listed the missing
  expression-plot entries.

## Testing

- Added label/value integrity tests covering every exported function that
  accepts a gene set. They assert two properties the 1.0.0 defect violated:
  results do not depend on the order genes were supplied in, and the value
  reported for a gene is that gene's own value, checked against an independent
  recomputation from the object. The previous suite asserted shape --
  `is.matrix()`, `nrow()`, class -- all of which the broken output satisfied.
- Verified by mutation testing: reintroducing the exact 1.0.0 defect now fails
  five tests across three files, and replacing a name-keyed lookup with a
  positional one elsewhere fails nine.

## Internal

- Removed duplicate definitions of `.prepare_sample_metadata()` and
  `.filter_genes()`, which existed in both `utils-internal.R` and
  `viz_related.R`, with collation order silently deciding which one ran.
- Removed unreachable helpers `.plot_pca()`, `.plot_mds()`,
  `.prepare_mds_dataframe()`, `.prepare_corr_matrix()`, `.plot_corr_heatmap()`,
  and `.cluster_log2fc_matrix()` along with their generated man pages.
- pkgdown-only articles under `vignettes/guides/` and `vignettes/workflows/`
  are no longer shipped in the source tarball.

# VISTA 1.0.0

- First Bioconductor release (Bioconductor 3.23).

# VISTA 0.99.8

# VISTA 0.99.7

# VISTA 0.99.6

# VISTA 0.99.5

- Added `read_vista_counts()`, `read_vista_metadata()`, and
  `match_vista_inputs()` to standardize common RNA-seq input formats without
  changing the existing `create_vista()` API.
- Added `derive_vista_metadata()` to bootstrap starter sample metadata from
  count-derived sample names using split- or regex-based parsing.
- Added lightweight import support for plain count tables, featureCounts,
  STAR gene counts, HTSeq-count, tximport-like inputs, and RSEM gene result
  files.

# VISTA 0.99.4

# VISTA 0.99.3

# VISTA 0.99.2

# VISTA 0.99.1

- `example_vista()` now uses a precomputed default object to reduce example,
  test, and package-check runtime while preserving the existing API.

# VISTA 0.99.0

*Submitted to Bioconductor 2026-02-11*

## Overview

VISTA (Visualization Toolkit for Transcriptomic Analysis) provides a unified S4-based framework for differential expression analysis of RNA-seq data, wrapping DESeq2 and edgeR workflows with consistent metadata management and rich visualization capabilities.

## Key Features

### Core Infrastructure
- S4 `VISTA` class extending `SummarizedExperiment` for standardized data management
- Unified differential expression workflow supporting DESeq2 and edgeR backends
- Consistent metadata structure for comparisons, cutoffs, and group information
- Flexible color palette system for visualizations

### Visualization Suite (28+ functions)
- **Dimension reduction**: PCA, MDS plots with customizable aesthetics
- **DE results**: Volcano plots, MA plots, DEG count barplots
- **Expression**: Barplots, boxplots, violin plots, density plots, joyplots, heatmaps
- **Comparisons**: Venn diagrams, alluvial plots, correlation heatmaps, pairwise plots
- **Fold-change**: Scatter plots, barplots, matrix visualizations, chromosome plots

### Functional Analysis
- MSigDB enrichment with flexible ID mapping (SYMBOL, ENSEMBL, ENTREZID)
- GO enrichment analysis (BP, MF, CC ontologies)
- KEGG pathway enrichment
- GSEA support with customizable gene sets
- Integrated visualization functions for enrichment results

### Optional Features
- Cell-type deconvolution via xCell2 integration
- Automated report generation with Quarto
- Accessor functions for all metadata components

## Implementation Details

- Comprehensive input validation and edge case handling
- Extensive test suite (>70% coverage)
- Complete roxygen2 documentation with runnable examples
- BiocStyle vignettes demonstrating complete workflows
- Proper namespace management and import declarations

## Bug Fixes

- Fixed contradictory roxygen documentation markers in internal utilities
- Added missing `@importFrom` declarations across all modules
- Improved error messages for invalid inputs
- Enhanced edge case handling in visualization functions
- Heatmap utilities now validate non-character `genes` input explicitly, support
  minimal-call defaults, and allow custom colours for multi-column annotations
- DEG count pie/donut plots now optionally include non-DE genes as an `"Other"`
  slice and support configurable label text colour
- `get_genes_by_regulation()` now supports top-gene ranking by `abs(log2fc)`
  and optional annotated table output
- PCA/MDS/UMAP plots now accept the standardized `use_group_colors` argument
  while keeping `use_vista_colors` as a deprecated compatibility alias

---
