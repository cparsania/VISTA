# Package index

## VISTA Class & Constructors

Create and manage VISTA objects for differential expression analysis.

- [`VISTA-class`](https://cparsania.github.io/VISTA/reference/VISTA-class.md)
  [`VISTA`](https://cparsania.github.io/VISTA/reference/VISTA-class.md)
  : VISTA S4 Class Definition
- [`create_vista()`](https://cparsania.github.io/VISTA/reference/create_vista.md)
  : Create a VISTA Object with Internal DE Analysis
- [`as_vista()`](https://cparsania.github.io/VISTA/reference/as_vista.md)
  : Coerce SummarizedExperiment to VISTA
- [`run_deseq_analysis()`](https://cparsania.github.io/VISTA/reference/run_deseq_analysis.md)
  [`run_edger_analysis()`](https://cparsania.github.io/VISTA/reference/run_deseq_analysis.md)
  [`run_limma_analysis()`](https://cparsania.github.io/VISTA/reference/run_deseq_analysis.md)
  : Run Differential Expression Analysis with DESeq2, edgeR, or
  limma-voom
- [`set_rowdata()`](https://cparsania.github.io/VISTA/reference/set_rowdata.md)
  : Set or append rowData annotations on a VISTA object
- [`print(`*`<VISTA>`*`)`](https://cparsania.github.io/VISTA/reference/print.VISTA.md)
  : Print a VISTA object

## Object Methods

Bioconductor methods a VISTA object supports: display, subsetting,
schema migration, and hand-off to DESeq2.

- [`show(`*`<VISTA>`*`)`](https://cparsania.github.io/VISTA/reference/VISTA-show.md)
  : Display a VISTA object
- [`` `[`( ``*`<VISTA>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://cparsania.github.io/VISTA/reference/VISTA-subset.md)
  : Subset a VISTA object
- [`updateObject(`*`<VISTA>`*`)`](https://cparsania.github.io/VISTA/reference/updateObject-VISTA-method.md)
  : Update a VISTA object to the current metadata schema
- [`as_deseq_dataset()`](https://cparsania.github.io/VISTA/reference/as_deseq_dataset.md)
  : Convert a VISTA object to a DESeqDataSet

## Accessor Functions

Extract components from a VISTA object.

- [`comparisons(`*`<VISTA>`*`)`](https://cparsania.github.io/VISTA/reference/VISTA-accessors.md)
  [`deg_summary(`*`<VISTA>`*`)`](https://cparsania.github.io/VISTA/reference/VISTA-accessors.md)
  [`cutoffs(`*`<VISTA>`*`)`](https://cparsania.github.io/VISTA/reference/VISTA-accessors.md)
  [`norm_counts(`*`<VISTA>`*`)`](https://cparsania.github.io/VISTA/reference/VISTA-accessors.md)
  [`sample_info(`*`<VISTA>`*`)`](https://cparsania.github.io/VISTA/reference/VISTA-accessors.md)
  [`row_data(`*`<VISTA>`*`)`](https://cparsania.github.io/VISTA/reference/VISTA-accessors.md)
  [`group_colors(`*`<VISTA>`*`)`](https://cparsania.github.io/VISTA/reference/VISTA-accessors.md)
  [`group_palette(`*`<VISTA>`*`)`](https://cparsania.github.io/VISTA/reference/VISTA-accessors.md)
  : Accessor Methods for VISTA Object
- [`counts(`*`<VISTA>`*`)`](https://cparsania.github.io/VISTA/reference/counts.md)
  : Raw counts from a VISTA object
- [`set_de_source()`](https://cparsania.github.io/VISTA/reference/set_de_source.md)
  : Set active DE source in a VISTA object
- [`set_vista_group_colors()`](https://cparsania.github.io/VISTA/reference/set_vista_group_colors.md)
  : Set manual group colors in a VISTA object
- [`set_vista_comparison_colors()`](https://cparsania.github.io/VISTA/reference/set_vista_comparison_colors.md)
  : Set manual comparison colors in a VISTA object

## Input Preparation

Standardize count matrices and sample metadata before constructing a
VISTA object.

- [`read_vista_counts()`](https://cparsania.github.io/VISTA/reference/read_vista_counts.md)
  : Read and standardize count inputs for VISTA
- [`derive_vista_metadata()`](https://cparsania.github.io/VISTA/reference/derive_vista_metadata.md)
  : Derive starter sample metadata from count sample names
- [`read_vista_metadata()`](https://cparsania.github.io/VISTA/reference/read_vista_metadata.md)
  : Read and standardize sample metadata for VISTA
- [`match_vista_inputs()`](https://cparsania.github.io/VISTA/reference/match_vista_inputs.md)
  : Match count and metadata inputs for VISTA

## Quality Control Plots

Sample-level quality control and exploratory visualizations.

- [`get_pca_plot()`](https://cparsania.github.io/VISTA/reference/get_pca_plot.md)
  : PCA plot
- [`get_mds_plot()`](https://cparsania.github.io/VISTA/reference/get_mds_plot.md)
  : Generate an MDS plot for samples in a VISTA object
- [`get_umap_plot()`](https://cparsania.github.io/VISTA/reference/get_umap_plot.md)
  : Generate a UMAP plot for samples in a VISTA object
- [`get_corr_heatmap()`](https://cparsania.github.io/VISTA/reference/get_corr_heatmap.md)
  : Draw a sample correlation heatmap
- [`get_pairwise_corr_plot()`](https://cparsania.github.io/VISTA/reference/get_pairwise_corr_plot.md)
  : Plot pairwise correlations between samples

## Differential Expression Plots

Visualize differential expression results across comparisons.

- [`get_volcano_plot()`](https://cparsania.github.io/VISTA/reference/get_volcano_plot.md)
  : Generate a volcano plot for a comparison in a VISTA object
- [`get_ma_plot()`](https://cparsania.github.io/VISTA/reference/get_ma_plot.md)
  : Generate MA plot from a VISTA object
- [`get_deg_count_barplot()`](https://cparsania.github.io/VISTA/reference/get_deg_count_barplot.md)
  : Barplot of DEG counts (Up/Down) across comparisons
- [`get_deg_count_pieplot()`](https://cparsania.github.io/VISTA/reference/get_deg_count_pieplot.md)
  : Pie chart of DEG counts (Up/Down) across comparisons
- [`get_deg_count_donutplot()`](https://cparsania.github.io/VISTA/reference/get_deg_count_donutplot.md)
  : Donut chart of DEG counts (Up/Down) across comparisons
- [`get_deg_venn_diagram()`](https://cparsania.github.io/VISTA/reference/get_deg_venn_diagram.md)
  : DEG Venn diagram
- [`get_deg_alluvial()`](https://cparsania.github.io/VISTA/reference/get_deg_alluvial.md)
  : Plot alluvial diagram showing gene regulation transitions across
  comparisons
- [`get_genes_by_regulation()`](https://cparsania.github.io/VISTA/reference/get_genes_by_regulation.md)
  : Get Genes by Regulation

## Expression Plots

Visualize gene expression patterns across samples and groups.

- [`get_expression_heatmap()`](https://cparsania.github.io/VISTA/reference/get_expression_heatmap.md)
  : Expression heatmap
- [`get_expression_barplot()`](https://cparsania.github.io/VISTA/reference/get_expression_barplot.md)
  : Plot gene expression as barplots
- [`get_expression_boxplot()`](https://cparsania.github.io/VISTA/reference/get_expression_boxplot.md)
  : Plot gene expression distributions as boxplots
- [`get_expression_violinplot()`](https://cparsania.github.io/VISTA/reference/get_expression_violinplot.md)
  : Violin plot of expression values
- [`get_expression_density()`](https://cparsania.github.io/VISTA/reference/get_expression_density.md)
  : Plot expression distributions as density curves
- [`get_expression_scatter()`](https://cparsania.github.io/VISTA/reference/get_expression_scatter.md)
  : Compare normalized expression between two samples or groups
- [`get_expression_lineplot()`](https://cparsania.github.io/VISTA/reference/get_expression_lineplot.md)
  : Gene expression line plot
- [`get_expression_lollipop()`](https://cparsania.github.io/VISTA/reference/get_expression_lollipop.md)
  : Plot expression as a lollipop chart
- [`get_expression_joyplot()`](https://cparsania.github.io/VISTA/reference/get_expression_joyplot.md)
  : Plot expression distributions as ridgelines
- [`get_expression_raincloud()`](https://cparsania.github.io/VISTA/reference/get_expression_raincloud.md)
  : Raincloud plot of expression values
- [`get_expression_matrix()`](https://cparsania.github.io/VISTA/reference/get_expression_matrix.md)
  : Retrieve an expression matrix from a VISTA object
- [`get_expression_chromosome_plot()`](https://cparsania.github.io/VISTA/reference/get_expression_chromosome_plot.md)
  : Chromosome plot for expression

## Fold-Change Plots

Visualize and compare log2 fold-changes across comparisons.

- [`get_foldchange_heatmap()`](https://cparsania.github.io/VISTA/reference/get_foldchange_heatmap.md)
  : Fold-change heatmap
- [`get_foldchange_scatter()`](https://cparsania.github.io/VISTA/reference/get_foldchange_scatter.md)
  : Fold-change scatterplot between two comparisons
- [`get_foldchange_barplot()`](https://cparsania.github.io/VISTA/reference/get_foldchange_barplot.md)
  : Plot fold-change barplots across comparisons for selected genes
- [`get_foldchange_boxplot()`](https://cparsania.github.io/VISTA/reference/get_foldchange_boxplot.md)
  : Plot fold-change distributions across comparisons
- [`get_foldchange_raincloud()`](https://cparsania.github.io/VISTA/reference/get_foldchange_raincloud.md)
  : Raincloud plot of fold-change distributions across comparisons
- [`get_foldchange_lineplot()`](https://cparsania.github.io/VISTA/reference/get_foldchange_lineplot.md)
  : Fold-change line plot across comparisons
- [`get_foldchange_lollipop()`](https://cparsania.github.io/VISTA/reference/get_foldchange_lollipop.md)
  : Fold-change plotting helpers (overview)
- [`get_foldchange_matrix()`](https://cparsania.github.io/VISTA/reference/get_foldchange_matrix.md)
  : Extract a log2 fold-change matrix
- [`get_foldchange_chromosome_plot()`](https://cparsania.github.io/VISTA/reference/get_foldchange_chromosome_plot.md)
  : Chromosome plot for fold change

## Functional Enrichment

Gene set enrichment and pathway analysis.

- [`get_msigdb_enrichment()`](https://cparsania.github.io/VISTA/reference/get_msigdb_enrichment.md)
  : Run MSigDB enrichment directly from a VISTA comparison
- [`get_go_enrichment()`](https://cparsania.github.io/VISTA/reference/get_go_enrichment.md)
  : Run GO enrichment directly from a VISTA comparison
- [`get_kegg_enrichment()`](https://cparsania.github.io/VISTA/reference/get_kegg_enrichment.md)
  : Run KEGG enrichment directly from a VISTA comparison
- [`get_gsea()`](https://cparsania.github.io/VISTA/reference/get_gsea.md)
  : Gene set enrichment analysis (GSEA) from a VISTA comparison
- [`enrichMsigDB()`](https://cparsania.github.io/VISTA/reference/enrichMsigDB.md)
  : Perform MSigDB over-representation analysis on a VISTA object
- [`get_enrichment_plot()`](https://cparsania.github.io/VISTA/reference/get_enrichment_plot.md)
  : Plot enrichment results using -log10(FDR)
- [`get_pathway_genes()`](https://cparsania.github.io/VISTA/reference/get_pathway_genes.md)
  : Extract genes from enriched pathways
- [`get_pathway_heatmap()`](https://cparsania.github.io/VISTA/reference/get_pathway_heatmap.md)
  : Plot pathway-specific expression heatmaps from enrichment output
- [`get_enrichment_chord()`](https://cparsania.github.io/VISTA/reference/get_enrichment_chord.md)
  : Chord diagram of enrichment gene–pathway relationships

## Cell-Type Deconvolution

Estimate and visualize cell-type composition from bulk RNA-seq.

- [`run_cell_deconvolution()`](https://cparsania.github.io/VISTA/reference/run_cell_deconvolution.md)
  : Run Cell Deconvolution on Bulk RNA-seq from VISTA Object
- [`get_cell_fractions()`](https://cparsania.github.io/VISTA/reference/get_cell_fractions.md)
  : Retrieve stored cell fraction estimates
- [`get_celltype_barplot()`](https://cparsania.github.io/VISTA/reference/get_celltype_barplot.md)
  : Plot cell-type composition as stacked bars
- [`get_celltype_group_dotplot()`](https://cparsania.github.io/VISTA/reference/get_celltype_group_dotplot.md)
  : Plot group-level deconvolution scores as dot plot
- [`get_celltype_heatmap()`](https://cparsania.github.io/VISTA/reference/get_celltype_heatmap.md)
  : Plot cell-type deconvolution heatmap

## Reporting

Generate automated analysis reports.

- [`run_vista_report()`](https://cparsania.github.io/VISTA/reference/run_vista_report.md)
  : Generate a publication-ready VISTA workflow report

## Validation & Benchmarking

Verify VISTA objects and benchmark backend fidelity against standalone
workflows.

- [`validate_vista()`](https://cparsania.github.io/VISTA/reference/validate_vista.md)
  : Validate a VISTA object
- [`validate_vista_deep()`](https://cparsania.github.io/VISTA/reference/validate_vista_deep.md)
  : Deep validation of VISTA differential-expression fidelity
- [`benchmark_vista_equivalence()`](https://cparsania.github.io/VISTA/reference/benchmark_vista_equivalence.md)
  : Benchmark VISTA against standalone differential-expression backends

## Exporting

Save publication-ready plots, tables, and reproducible asset bundles.

- [`save_vista_plot()`](https://cparsania.github.io/VISTA/reference/save_vista_plot.md)
  : Save a VISTA plot object to disk
- [`save_vista_data()`](https://cparsania.github.io/VISTA/reference/save_vista_data.md)
  : Save VISTA tabular outputs to disk
- [`export_vista_assets()`](https://cparsania.github.io/VISTA/reference/export_vista_assets.md)
  : Export a complete VISTA asset bundle

## Example Data

Built-in datasets for examples and testing.

- [`example_vista()`](https://cparsania.github.io/VISTA/reference/example_vista.md)
  : Build a small example VISTA object
- [`count_data`](https://cparsania.github.io/VISTA/reference/count_data.md)
  : Example RNA-seq count matrix shipped with VISTA
- [`sample_metadata`](https://cparsania.github.io/VISTA/reference/sample_metadata.md)
  : Sample metadata accompanying the VISTA airway example counts

## Deprecated

Renamed arguments that still work, and the release each becomes defunct
in.

- [`VISTA-deprecated`](https://cparsania.github.io/VISTA/reference/VISTA-deprecated.md)
  : Deprecated and defunct arguments in VISTA
