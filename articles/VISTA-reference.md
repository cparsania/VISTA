# VISTA Function Reference

## Overview

This reference vignette summarizes the VISTA API with argument details.
It is intended as a high‑level map of functions and parameters. For
authoritative signatures, run `?function_name` in R (or open the
generated Rd page) for the current version.

## Core construction

### `create_vista()`

Creates a VISTA object and runs DE analysis.

Arguments:

- `counts`: Raw counts (matrix/data.frame) with a gene ID column +
  sample columns.
- `sample_info`: Sample metadata with `sample_names` or matching
  rownames.
- `column_geneid`: Column name in `counts` that holds gene identifiers.
- `group_column`: Column in `sample_info` used to group samples.
- `group_numerator`: Numerator group(s) for comparisons.
- `group_denominator`: Denominator group(s), same length/order as
  numerator.
- `method`: `"deseq2"` (default), `"edger"`, or `"both"`.
- `min_counts`: Minimum total counts per gene to retain.
- `min_replicates`: Minimum samples per group meeting `min_counts`.
- `log2fc_cutoff`: Absolute LFC threshold for DEG calling.
- `pval_cutoff`: P‑value cutoff.
- `p_value_type`: `"padj"` or `"pvalue"` (for
  [`create_vista()`](../reference/create_vista.md) this controls DE
  calling across all methods).
- `covariates`: Optional additive covariate columns from `sample_info`.
- `design_formula`: Optional explicit design formula including
  `group_column`.
- `consensus_mode`: For `method = "both"`, `"intersection"` or
  `"union"`.
- `consensus_log2fc`: For consensus tables, how `log2fc` is populated.
- `result_source`: Active DE source (`"consensus"`, `"deseq2"`,
  `"edger"`).
- `group_palette`: Colorspace palette name for group colors.
- `comparison_palette`: Palette name for comparison colors.

### `vista()`

Low‑level constructor for a pre‑normalized `SummarizedExperiment`.

Arguments:

- `norm_counts`: Normalized expression matrix.
- `sample_info`: Sample metadata.
- `row_data`: Feature metadata.
- `comparisons`: Named list of DE tables.
- `deg_summary`: Named list of summary tables.
- `cutoffs`: Named list of thresholds.
- `group_column`: Grouping column name.
- `group_palette`: Group palette name.

## Accessors

- `comparisons(x, source = "active")`: DE results list (active or
  method-specific source).
- `deg_summary(x, source = "active")`: DEG summary list (active or
  method-specific source).
- `cutoffs(x)`: DE threshold list.
- `norm_counts(x, summarise = FALSE)`: Normalized counts (optionally by
  group).
- `sample_info(x)`: Sample metadata.
- `row_data(x)`: Feature metadata.
- `group_colors(x)`: Named group colors.
- `group_palette(x)`: Palette name.
- `set_de_source(x, source)`: Switch active DE source used by plotting
  functions.
- `set_vista_group_colors(x, color_map)`: Set manual named group colors.
- `set_vista_comparison_colors(x, color_map)`: Set manual named
  comparison colors.

## DE wrappers

- `run_deseq_analysis(...)`: DESeq2 pipeline; arguments mirror
  [`create_vista()`](../reference/create_vista.md).
- `run_edger_analysis(...)`: edgeR pipeline; arguments mirror
  [`create_vista()`](../reference/create_vista.md).

## Expression matrices

### `get_expression_matrix()`

Arguments:

- `x`: VISTA object.
- `assay_name`: Assay to extract (default `"norm_counts"`).
- `genes`: Optional gene IDs to keep.
- `samples`: Optional sample IDs (or group labels when
  `summarise = TRUE`).
- `group_column`: Grouping column for summarization.
- `summarise`: Average replicates per group when `TRUE`.
- `transform`: `"none"`, `"log2"`, or `"zscore"`.

## Expression plots

### `get_expression_boxplot()`

Arguments:

- `x`, `genes`, `sample_group`, `group_column`, `facet`,
  `log_transform`.
- `display_id`, `display_from`, `display_orgdb`.
- `facet_scales`, `stats_group`, `p.label`, `comparisons`.
- `pool_genes`, `x_by`, `facet_by`, `fill_by`.

Notes: `genes` capped at 20 when provided. If `facet_by = "gene"` and
`genes` is NULL, top 20 variable genes are selected automatically.

### `get_expression_barplot()`

Arguments:

- `x`, `genes` (≤10), `sample_group`, `group_column`, `log_transform`.
- `stats_group`, `facet_scale`, `p.label`, `comparisons`.
- `display_id`, `display_from`, `display_orgdb`.

### `get_expression_lollipop()`

Arguments:

- `x`, `genes` (≤15), `sample_group`, `group_column`, `log_transform`.
- `facet_scale`, `point_size`, `line_size`, `label`, `label_digits`.
- `display_id`, `display_from`, `display_orgdb`.

### `get_expression_lineplot()`

Arguments:

- `x`, `genes`, `sample_group`, `group_column`.
- `value_transform` (`"log2"`, `"zscore"`, `"none"`).
- `summarise_groups`, `facet_by_group`, `color_palette`.

### Distribution helpers

- `get_expression_violinplot(...)`
- `get_expression_densityplot(...)`
- `get_expression_raincloud(...)`

Common arguments: `genes`, `sample_group`, `group_column`,
`value_transform`, `summarise`, `facet` and plotting controls.

## Dimensionality reduction

### `get_pca_plot()`

Arguments:

- `x`, `sample_group`, `group_column`, `genes`, `top_n_genes`.
- `label_replicates`, `label_size`, `circle_size`, `sample_colors`.
- `sample.seed`, `show_clusters`.

### `get_mds_plot()`

Arguments:

- `x`, `sample_group`, `group_column`, `genes`, `top_n_genes`.
- `label_replicates`, `label_size`, `circle_size`, `sample_colors`.
- `sample.seed`, `show_clusters`.

## Correlation plots

- `get_pairwise_corr_plot(...)`: GGally ggpairs for selected
  samples/genes.
- `get_corr_heatmap(...)`: correlation heatmap (with optional
  clustering).

## DEG summaries and comparison plots

- `get_deg_count_barplot(...)`: DEG counts by comparison/regulation.
- `get_deg_count_pieplot(...)`: circular DEG composition summary (pie).
- `get_deg_count_donutplot(...)`: circular DEG composition summary
  (donut).
- `get_deg_venn_diagram(...)`: overlap of DEG sets (2–4 comparisons).
- `plot_deg_alluvial(...)`: transitions of regulation across
  comparisons.

## Volcano/MA

- `get_volcano_plot(...)`: volcano plot with optional label mapping.
- `get_ma_plot(...)`: MA plot with optional label mapping.

Common label arguments: `display_id`, `display_from`, `display_orgdb`.

## Fold change utilities

- `get_foldchange_matrix(...)`: log2FC matrix for comparisons.
- `get_foldchange_heatmap(...)`: heatmap of log2FC.
- `get_foldchange_lineplot(...)`: line plot across comparisons.
- `get_foldchange_boxplot(...)`: distribution across comparisons.
- `get_foldchange_barplot(...)`: barplot for selected genes.
- `get_foldchange_scatter(...)`: comparison‑comparison scatter with
  concordance colors.

## Enrichment

- `enrichMsigDB(...)`: base ORA using MSigDB.
- `get_msigdb_enrichment(...)`: convenience wrapper for VISTA
  comparisons.
- `get_go_enrichment(...)`: GO ORA for a comparison.
- `get_kegg_enrichment(...)`: KEGG ORA for a comparison.
- `get_gsea(...)`: GSEA for MSigDB/GO/KEGG.
- `get_enrichment_plot(...)`: dot plot of enrichment results.
- `get_pathway_genes(...)`: extract per-pathway genes from enrichment
  outputs.
- `get_pathway_heatmap(...)`: visualize pathway-derived genes with
  `get_expression_heatmap(...)`.
- `get_enrichment_chord(...)`: chord diagram of gene–pathway
  relationships (requires circlize); supports
  `gene_order_by = "foldchange"` or `"abs_foldchange"`.

## Deconvolution

- `run_cell_deconvolution(...)`: xCell2 workflow.
- `get_cell_fractions(x)`: retrieve stored fractions.
- `plot_celltype_barplot(...)`: barplot of cell fractions.

## Reporting

- `run_vista_report(config, output_file)`: render Quarto report.

## Row metadata helper

- `set_rowdata(...)`: attach/append annotations using `orgdb` or
  provided table.

## Utility

- `get_genes_by_regulation(...)`: extract regulated gene sets for a
  comparison.

## Notes

This reference is intentionally exhaustive but compact. For full
parameter validation rules and return types, consult each function’s
help page.
