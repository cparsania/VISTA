# VISTA <img src="man/figures/logo.png" align="right" height="240"/>

> **V**isualization and **I**ntegrated **S**ystem for **T**ranscriptomic **A**nalysis

<!-- badges: start -->

[![Lifecycle: maturing](https://lifecycle.r-lib.org/articles/figures/lifecycle-maturing.svg)](https://lifecycle.r-lib.org/articles/stages.html#maturing)
[![R-CMD-check](https://github.com/cparsania/VISTA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cparsania/VISTA/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/cparsania/VISTA/graph/badge.svg)](https://app.codecov.io/gh/cparsania/VISTA)
[![pkgdown site](https://img.shields.io/badge/pkgdown-site-blue)](https://cparsania.github.io/VISTA/)

<!-- badges: end -->

## Overview

VISTA is a Bioconductor-oriented framework for RNA-seq differential expression analysis that keeps data, statistics, annotations, and visualization in a single `SummarizedExperiment`-based object. It is designed for analysts who want a reproducible workflow without rebuilding the same plotting and result-handling code for every project.

VISTA provides:

- A single entry point for `DESeq2`, `edgeR`, `limma-voom`, and a DESeq2/edgeR consensus workflow.
- Design-aware modeling through covariates or an explicit design formula.
- Publication-ready visualizations for QC, differential expression, expression patterns, fold-change structure, pathway enrichment, and deconvolution.
- Consistent color management across comparisons and groups.
- Export helpers and a YAML-driven Quarto reporting workflow for shareable analyses.

## Why VISTA

Most RNA-seq projects repeat the same sequence of steps: normalize counts, fit differential models, extract contrasts, make QC plots, label genes, summarize pathways, and assemble figures for collaborators. The friction is usually not the statistics alone; it is the repeated glue code that connects them.

VISTA addresses that by organizing the workflow around one validated object:

- `assay(x)` stores normalized expression data.
- `rowData(x)` stores feature annotations.
- `colData(x)` stores sample metadata.
- `comparisons(x)` stores differential expression tables.
- VISTA accessors and plotting functions operate directly on that object.

That design makes it easier to move from raw counts to a consistent analysis narrative without switching data structures between each step.

## Installation

### Bioconductor

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("VISTA")
```

### Development version

```r
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}

pak::pak("cparsania/VISTA")
```

## Quick Start

```r
library(VISTA)

# Example data shipped with the package
# count_data: gene-by-sample count matrix
# sample_metadata: sample annotations

data("count_data", package = "VISTA")
data("sample_metadata", package = "VISTA")

vista <- create_vista(
  counts = count_data,
  sample_info = sample_metadata,
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control",
  method = "deseq2",
  log2fc_cutoff = 1,
  pval_cutoff = 0.05
)

comp <- names(comparisons(vista))[1]

# Inspect the object and the first differential comparison
vista
head(comparisons(vista)[[comp]][, c("gene_id", "log2fc", "padj")])

# QC
get_pca_plot(vista, label_replicates = TRUE)
get_mds_plot(vista)
get_corr_heatmap(vista)
# Optional nonlinear view (requires uwot)
# get_umap_plot(vista, color_by = "cell")

# Differential expression
get_volcano_plot(vista, sample_comparison = comp)
get_ma_plot(vista, sample_comparison = comp)
get_deg_count_barplot(vista)
get_deg_count_donutplot(vista)

# Expression-focused views
up_genes <- get_genes_by_regulation(
  vista,
  sample_comparisons = comp,
  regulation = "Up"
)[[comp]]

get_expression_heatmap(vista, genes = up_genes[1:40], kmeans_k = 3)
get_expression_barplot(vista, genes = up_genes[1:3], by = "sample", facet_by = "gene")
get_expression_lollipop(vista, genes = up_genes[1:3], by = "sample", facet_by = "gene")

# Fold-change views
get_foldchange_barplot(vista, genes = up_genes[1:6], sample_comparisons = comp, facet_by = "gene")
get_foldchange_lollipop(vista, sample_comparison = comp, genes = up_genes[1:6], facet_by = "gene")
```

## Core Workflow

### 1. Build the analysis object

```r
vista <- create_vista(
  counts = count_data,
  sample_info = sample_metadata,
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control",
  method = "limma"  # or "deseq2", "edger", "both"
)
```

Use `method = "both"` when you want a DESeq2/edgeR consensus workflow:

```r
vista_consensus <- create_vista(
  counts = count_data,
  sample_info = sample_metadata,
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control",
  method = "both",
  result_source = "consensus"
)

# Switch the active DE table if needed
vista_consensus <- set_de_source(vista_consensus, "edger")
```

### 2. Add covariates or an explicit design

```r
# Add sample-level covariates to adjust the DE model
vista_cov <- create_vista(
  counts = count_data,
  sample_info = sample_metadata,
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control",
  covariates = "cell"
)

# Or provide a full design formula
vista_design <- create_vista(
  counts = count_data,
  sample_info = sample_metadata,
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control",
  design_formula = ~ cell + cond_long
)
```

### 3. Add feature annotations

```r
# Example for human data
vista <- set_rowdata(
  vista,
  orgdb = org.Hs.eg.db,
  columns = c("SYMBOL", "GENENAME", "ENTREZID")
)
```

### 4. Explore, export, and report

```r
# Export selected tables / matrices / plots
export_vista_assets(
  vista,
  out_dir = "vista_assets",
  include_data = c("comparison", "norm_counts", "sample_info")
)

# Render a parameterized HTML report from YAML
file.copy(
  system.file("reports", "vista-report-template.yml", package = "VISTA"),
  "vista-report.yml"
)
run_vista_report("vista-report.yml")
```

## Visualization Coverage

### Quality control and sample structure

- `get_pca_plot()`
- `get_mds_plot()`
- `get_umap_plot()`
- `get_corr_heatmap()`
- `get_pairwise_corr_plot()`

### Differential expression summaries

- `get_volcano_plot()`
- `get_ma_plot()`
- `get_deg_count_barplot()`
- `get_deg_count_pieplot()`
- `get_deg_count_donutplot()`
- `get_deg_venn_diagram()`
- `plot_deg_alluvial()`

### Expression patterns

- `get_expression_heatmap()`
- `get_expression_boxplot()`
- `get_expression_violinplot()`
- `get_expression_barplot()`
- `get_expression_lollipop()`
- `get_expression_scatter()`
- `get_expression_lineplot()`
- `get_expression_density()`
- `get_expression_joyplot()`
- `get_expression_raincloud()`

### Fold-change structure

- `get_foldchange_scatter()`
- `get_foldchange_barplot()`
- `get_foldchange_lollipop()`
- `get_foldchange_boxplot()`
- `get_foldchange_lineplot()`
- `get_foldchange_heatmap()`
- `get_foldchange_matrix()`
- `get_foldchange_chromosome_plot()`

### Pathway and enrichment views

- `get_msigdb_enrichment()`
- `get_go_enrichment()`
- `get_kegg_enrichment()`
- `get_gsea()`
- `get_enrichment_plot()`
- `get_enrichment_chord()`
- `get_enrichment_network()`
- `get_pathway_genes()`
- `get_pathway_heatmap()`

### Deconvolution (optional workflow)

- `run_cell_deconvolution()`
- `plot_celltype_barplot()`
- `get_celltype_group_dotplot()`
- `get_celltype_heatmap()`

## Example Analyses

### Enrichment from a VISTA comparison

```r
# Human example: requires org.Hs.eg.db
msig <- get_msigdb_enrichment(
  vista,
  sample_comparison = comp,
  regulation = "Up",
  orgdb = org.Hs.eg.db,
  species = "Homo sapiens",
  msigdb_category = "H"
)

go_bp <- get_go_enrichment(
  vista,
  sample_comparison = comp,
  regulation = "Up",
  ont = "BP",
  orgdb = org.Hs.eg.db,
  species = "Homo sapiens"
)

kegg <- get_kegg_enrichment(
  vista,
  sample_comparison = comp,
  regulation = "Up",
  orgdb = org.Hs.eg.db,
  species = "Homo sapiens"
)

get_enrichment_plot(msig$enrich)
```

### Multi-comparison overlap

```r
# Example only when your VISTA object contains multiple contrasts
get_deg_venn_diagram(
  vista_consensus,
  sample_comparisons = c("treatment1_VS_control", "control_VS_treatment1"),
  regulation = "Up"
)
```

### Consistent color control

```r
group_colors(vista)
vista <- set_vista_group_colors(
  vista,
  c(control = "#264653", treatment1 = "#E76F51")
)

# For multi-comparison objects
vista_consensus <- set_vista_comparison_colors(
  vista_consensus,
  c(treatment1_VS_control = "#6C5CE7")
)
```

## Bioconductor-Compatible Object Design

VISTA extends `SummarizedExperiment`, so standard Bioconductor workflows remain available.

```r
# Standard access
assay(vista)[1:5, 1:5]
rowData(vista)
colData(vista)
metadata(vista)

# VISTA accessors
comparisons(vista)
deg_summary(vista)
cutoffs(vista)
norm_counts(vista, summarise = TRUE)

# Validation helpers
validate_vista(vista, level = "full")
```

Advanced users can also coerce an existing `SummarizedExperiment`:

```r
# se <- your existing SummarizedExperiment
# vista_from_se <- as_vista(se, group_column = "cond_long")
# validate_vista(vista_from_se)
```

## Documentation

- Package website: <https://cparsania.github.io/VISTA/>
- Introduction article: <https://cparsania.github.io/VISTA/articles/VISTA-airway.html>
- Comparison workflow: <https://cparsania.github.io/VISTA/articles/VISTA-comparison.html>
- Color management: <https://cparsania.github.io/VISTA/articles/VISTA-colors.html>
- Function reference article: <https://cparsania.github.io/VISTA/articles/VISTA-reference.html>
- Deconvolution article: <https://cparsania.github.io/VISTA/articles/VISTA-deconvolution.html>

Local help remains available through standard R documentation, for example:

```r
?create_vista
?get_expression_heatmap
?get_go_enrichment
?run_vista_report
```

## Citation

If you use VISTA in published work, cite the package release used in your analysis.

```text
Parsania C (2026). VISTA: Visualization and Integrated System for Transcriptomic Analysis.
R package version 0.99.0.
```

## Support

- Issues: <https://github.com/cparsania/VISTA/issues>
- Bioconductor Support: <https://support.bioconductor.org>

## Contributing

Contributions should preserve reproducibility and backward compatibility.

Before opening a pull request:

- Add or update tests for functional changes.
- Run `devtools::document()` if roxygen comments changed.
- Run `devtools::test()` and `R CMD check`.
- Update `NEWS.md` for user-visible changes.

## License

GPL-3
