# Package index

## VISTA Class & Constructors

Create and manage VISTA objects for differential expression analysis.

- [`VISTA-class`](VISTA-class.md) : VISTA S4 Class Definition
- [`create_vista()`](create_vista.md) : Create a VISTA Object with
  Internal DE Analysis
- [`as_vista()`](as_vista.md) : Coerce SummarizedExperiment to VISTA
- [`run_deseq_analysis()`](run_deseq_analysis.md)
  [`run_edger_analysis()`](run_deseq_analysis.md)
  [`run_limma_analysis()`](run_deseq_analysis.md) : Run Differential
  Expression Analysis with DESeq2, edgeR, or limma-voom
- [`set_rowdata()`](set_rowdata.md) : Set or append rowData annotations
  on a VISTA object
- [`print(`*`<VISTA>`*`)`](print.VISTA.md) : Print a VISTA object like a
  SummarizedExperiment

## Accessor Functions

Extract components from a VISTA object.

- [`comparisons(`*`<VISTA>`*`)`](VISTA-accessors.md)
  [`deg_summary(`*`<VISTA>`*`)`](VISTA-accessors.md)
  [`cutoffs(`*`<VISTA>`*`)`](VISTA-accessors.md)
  [`norm_counts(`*`<VISTA>`*`)`](VISTA-accessors.md)
  [`sample_info(`*`<VISTA>`*`)`](VISTA-accessors.md)
  [`row_data(`*`<VISTA>`*`)`](VISTA-accessors.md)
  [`group_colors(`*`<VISTA>`*`)`](VISTA-accessors.md)
  [`group_palette(`*`<VISTA>`*`)`](VISTA-accessors.md) : Accessor
  Methods for VISTA Object
- [`set_de_source()`](set_de_source.md) : Set active DE source in a
  VISTA object
- [`set_vista_group_colors()`](set_vista_group_colors.md) : Set manual
  group colors in a VISTA object
- [`set_vista_comparison_colors()`](set_vista_comparison_colors.md) :
  Set manual comparison colors in a VISTA object

## Input Preparation

Standardize count matrices and sample metadata before constructing a
VISTA object.

- [`read_vista_counts()`](read_vista_counts.md) : Read and standardize
  count inputs for VISTA
- [`derive_vista_metadata()`](derive_vista_metadata.md) : Derive starter
  sample metadata from count sample names
- [`read_vista_metadata()`](read_vista_metadata.md) : Read and
  standardize sample metadata for VISTA
- [`match_vista_inputs()`](match_vista_inputs.md) : Match count and
  metadata inputs for VISTA

## Quality Control Plots

Sample-level quality control and exploratory visualizations.

- [`get_pca_plot()`](get_pca_plot.md) : PCA plot
- [`get_mds_plot()`](get_mds_plot.md) : Generate an MDS plot for samples
  in a VISTA object
- [`get_umap_plot()`](get_umap_plot.md) : Generate a UMAP plot for
  samples in a VISTA object
- [`get_corr_heatmap()`](get_corr_heatmap.md) : Draw a sample
  correlation heatmap
- [`get_pairwise_corr_plot()`](get_pairwise_corr_plot.md) : Plot
  pairwise correlations between samples

## Differential Expression Plots

Visualize differential expression results across comparisons.

- [`get_volcano_plot()`](get_volcano_plot.md) : Generate a volcano plot
  for a comparison in a VISTA object
- [`get_ma_plot()`](get_ma_plot.md) : Generate MA plot from a VISTA
  object
- [`get_deg_count_barplot()`](get_deg_count_barplot.md) : Barplot of DEG
  counts (Up/Down) across comparisons
- [`get_deg_count_pieplot()`](get_deg_count_pieplot.md) : Pie chart of
  DEG counts (Up/Down) across comparisons
- [`get_deg_count_donutplot()`](get_deg_count_donutplot.md) : Donut
  chart of DEG counts (Up/Down) across comparisons
- [`get_deg_venn_diagram()`](get_deg_venn_diagram.md) : DEG Venn diagram
- [`get_deg_alluvial()`](get_deg_alluvial.md) : Plot alluvial diagram
  showing gene regulation transitions across comparisons
- [`get_genes_by_regulation()`](get_genes_by_regulation.md) : Get Genes
  by Regulation

## Expression Plots

Visualize gene expression patterns across samples and groups.

- [`get_expression_heatmap()`](get_expression_heatmap.md) : Expression
  heatmap
- [`get_expression_barplot()`](get_expression_barplot.md) : Plot gene
  expression as barplots
- [`get_expression_boxplot()`](get_expression_boxplot.md) : Plot gene
  expression distributions as boxplots
- [`get_expression_violinplot()`](get_expression_violinplot.md) : Violin
  plot of expression values
- [`get_expression_density()`](get_expression_density.md) : Plot
  expression distributions as density curves
- [`get_expression_scatter()`](get_expression_scatter.md) : Compare
  normalized expression between two samples or groups
- [`get_expression_lineplot()`](get_expression_lineplot.md) : Gene
  expression line plot
- [`get_expression_lollipop()`](get_expression_lollipop.md) : Plot
  expression as a lollipop chart
- [`get_expression_joyplot()`](get_expression_joyplot.md) : Plot
  expression distributions as ridgelines
- [`get_expression_raincloud()`](get_expression_raincloud.md) :
  Raincloud plot of expression values
- [`get_expression_matrix()`](get_expression_matrix.md) : Retrieve an
  expression matrix from a VISTA object
- [`get_expression_chromosome_plot()`](get_expression_chromosome_plot.md)
  : Chromosome plot for expression

## Fold-Change Plots

Visualize and compare log2 fold-changes across comparisons.

- [`get_foldchange_heatmap()`](get_foldchange_heatmap.md) : Fold-change
  heatmap
- [`get_foldchange_scatter()`](get_foldchange_scatter.md) : Fold-change
  scatterplot between two comparisons
- [`get_foldchange_barplot()`](get_foldchange_barplot.md) : Plot
  fold-change barplots across comparisons for selected genes
- [`get_foldchange_boxplot()`](get_foldchange_boxplot.md) : Plot
  fold-change distributions across comparisons
- [`get_foldchange_raincloud()`](get_foldchange_raincloud.md) :
  Raincloud plot of fold-change distributions across comparisons
- [`get_foldchange_lineplot()`](get_foldchange_lineplot.md) :
  Fold-change line plot across comparisons
- [`get_foldchange_lollipop()`](get_foldchange_lollipop.md) :
  Fold-change plotting helpers (overview)
- [`get_foldchange_matrix()`](get_foldchange_matrix.md) : Extract a log2
  fold-change matrix
- [`get_foldchange_chromosome_plot()`](get_foldchange_chromosome_plot.md)
  : Chromosome plot for fold change

## Functional Enrichment

Gene set enrichment and pathway analysis.

- [`get_msigdb_enrichment()`](get_msigdb_enrichment.md) : Run MSigDB
  enrichment directly from a VISTA comparison
- [`get_go_enrichment()`](get_go_enrichment.md) : Run GO enrichment
  directly from a VISTA comparison
- [`get_kegg_enrichment()`](get_kegg_enrichment.md) : Run KEGG
  enrichment directly from a VISTA comparison
- [`get_gsea()`](get_gsea.md) : Gene set enrichment analysis (GSEA) from
  a VISTA comparison
- [`enrichMsigDB()`](enrichMsigDB.md) : Perform MSigDB
  over-representation analysis on a VISTA object
- [`get_enrichment_plot()`](get_enrichment_plot.md) : Plot enrichment
  results using -log10(FDR)
- [`get_pathway_genes()`](get_pathway_genes.md) : Extract genes from
  enriched pathways
- [`get_pathway_heatmap()`](get_pathway_heatmap.md) : Plot
  pathway-specific expression heatmaps from enrichment output
- [`get_enrichment_chord()`](get_enrichment_chord.md) : Chord diagram of
  enrichment gene–pathway relationships

## Cell-Type Deconvolution

Estimate and visualize cell-type composition from bulk RNA-seq.

- [`run_cell_deconvolution()`](run_cell_deconvolution.md) : Run Cell
  Deconvolution on Bulk RNA-seq from VISTA Object
- [`get_cell_fractions()`](get_cell_fractions.md) : Retrieve stored cell
  fraction estimates
- [`get_celltype_barplot()`](get_celltype_barplot.md) : Plot cell-type
  composition as stacked bars
- [`get_celltype_group_dotplot()`](get_celltype_group_dotplot.md) : Plot
  group-level deconvolution scores as dot plot
- [`get_celltype_heatmap()`](get_celltype_heatmap.md) : Plot cell-type
  deconvolution heatmap

## Reporting

Generate automated analysis reports.

- [`run_vista_report()`](run_vista_report.md) : Generate a
  publication-ready VISTA workflow report

## Validation & Benchmarking

Verify VISTA objects and benchmark backend fidelity against standalone
workflows.

- [`validate_vista()`](validate_vista.md) : Validate a VISTA object
- [`validate_vista_deep()`](validate_vista_deep.md) : Deep validation of
  VISTA differential-expression fidelity
- [`benchmark_vista_equivalence()`](benchmark_vista_equivalence.md) :
  Benchmark VISTA against standalone differential-expression backends

## Exporting

Save publication-ready plots, tables, and reproducible asset bundles.

- [`save_vista_plot()`](save_vista_plot.md) : Save a VISTA plot object
  to disk
- [`save_vista_data()`](save_vista_data.md) : Save VISTA tabular outputs
  to disk
- [`export_vista_assets()`](export_vista_assets.md) : Export a complete
  VISTA asset bundle

## Example Data

Built-in datasets for examples and testing.

- [`example_vista()`](example_vista.md) : Build a small example VISTA
  object
- [`count_data`](count_data.md) : Example RNA-seq count matrix shipped
  with VISTA
- [`sample_metadata`](sample_metadata.md) : Sample metadata accompanying
  the VISTA airway example counts
