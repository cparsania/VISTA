# VISTA <img src="man/figures/logo.png" align="right" height="240" />

> **V**isualization and **I**ntegrated **S**ystem for **T**ranscriptomic **A**nalysis

<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check](https://img.shields.io/badge/R%20CMD%20check-passing-brightgreen)](https://github.com/chiragparsania/VISTA/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

VISTA streamlines differential expression (DE) analysis workflows for RNA-seq data by providing a unified, **SummarizedExperiment**-based framework that wraps **DESeq2** and **edgeR**. The package offers:

- 🧬 **Unified DE workflow** - Single interface for both DESeq2 and edgeR methods
- ⚖️ **Consensus DE mode** - Optional intersection/union consensus across DESeq2 and edgeR
- 🧪 **Design-aware modeling** - Add covariates or explicit model formula in DE analysis
- 📊 **30+ visualization functions** - PCA, volcano plots, heatmaps, expression plots, and more
- 🔬 **Functional enrichment** - Integrated MSigDB, GO, and KEGG pathway analysis
- 🎨 **Consistent aesthetics** - Harmonized color palettes and publication-ready themes
- 📝 **Automated reporting** - Generate comprehensive analysis reports with Quarto
- ✅ **Extensive validation** - Input checking and clear error messages
- 🧪 **Comprehensive tests** - Broad test coverage for core workflows

## Installation

### From Bioconductor (Recommended - when available)

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("VISTA")
```

### Development Version

```r
# Install from GitHub
if (!requireNamespace("pak", quietly = TRUE))
    install.packages("pak")

pak::pak("chiragparsania/VISTA")
```

## Quick Start

```r
library(VISTA)

# Load your count matrix and sample metadata
data("count_data", package = "VISTA")
data("sample_metadata", package = "VISTA")
# Example datasets are shipped as tibbles
class(count_data)
class(sample_metadata)

# Create VISTA object with differential expression analysis
vista <- create_vista(
  counts = count_data,
  sample_info = sample_metadata,
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control",
  log2fc_cutoff = 1.0,
  pval_cutoff = 0.05,
  method = "deseq2"  # or "edger"
)

# Explore results
vista
names(comparisons(vista))
# compact summary preview
lapply(deg_summary(vista), head, 3)

# Quality control visualizations
get_pca_plot(vista, label_replicates = TRUE)
get_corr_heatmap(vista)

# Differential expression visualizations
get_volcano_plot(vista, sample_comparison = "treatment1_VS_control")
get_ma_plot(vista, sample_comparison = "treatment1_VS_control")
get_deg_count_barplot(vista)
get_deg_count_pieplot(vista)

# Optional: adjust for covariates in DE model
vista_cov <- create_vista(
  counts = count_data,
  sample_info = sample_metadata,
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control",
  method = "deseq2",
  covariates = "cell"
)

# Optional: DESeq2 + edgeR consensus workflow
vista_both <- create_vista(
  counts = count_data,
  sample_info = sample_metadata,
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control",
  method = "both",
  result_source = "consensus"
)

# Switch active DE source post-hoc
vista_both <- set_de_source(vista_both, "edger")

# Expression pattern analysis
genes_of_interest <- get_genes_by_regulation(
  vista,
  sample_comparisons = "treatment1_VS_control",
  regulation = "Up"
)$treatment1_VS_control

get_expression_heatmap(vista, genes = genes_of_interest[1:50], kmeans_k = 3)
get_expression_barplot(vista, genes = genes_of_interest[1:5], stats_group = TRUE)

# Functional enrichment
enrichment_results <- get_msigdb_enrichment(
  vista,
  sample_comparison = "treatment1_VS_control",
  regulation = "Up",
  msigdb_category = "H"  # Hallmark gene sets
)

get_enrichment_plot(enrichment_results$enrich)
```

## Key Features

### 🎯 Complete Analysis Workflow

```r
# 1. Set up analysis
vista <- create_vista(counts, sample_info, ...)

# 2. Add gene annotations
vista <- set_rowdata(
  vista,
  orgdb = org.Hs.eg.db,
  columns = c("SYMBOL", "GENENAME", "ENTREZID")
)

# 3. Quality control
get_pca_plot(vista)
get_mds_plot(vista)
get_corr_heatmap(vista)

# 4. Differential expression
get_volcano_plot(vista, sample_comparison = "A_VS_B")
get_deg_venn_diagram(vista, sample_comparison = c("A_VS_B", "C_VS_B"))

# 5. Expression patterns
get_expression_heatmap(vista, genes = my_genes, kmeans_k = 4)
get_expression_boxplot(vista, genes = my_genes)

# 6. Functional analysis
go_results <- get_go_enrichment(vista, regulation = "Up", ont = "BP")
kegg_results <- get_kegg_enrichment(vista, regulation = "Down")
gsea_results <- get_gsea(vista, sample_comparison = "A_VS_B")
```

### 📊 Visualization Suite

**Dimension Reduction:**
- `get_pca_plot()` - Principal component analysis
- `get_mds_plot()` - Multidimensional scaling

**DE Results:**
- `get_volcano_plot()` - Volcano plots with EnhancedVolcano integration
- `get_ma_plot()` - MA plots with customizable labeling
- `get_deg_count_barplot()` - DEG summary statistics
- `get_deg_count_pieplot()`, `get_deg_count_donutplot()` - Circular DEG composition summaries

**Expression Patterns:**
- `get_expression_barplot()`, `get_expression_boxplot()`, `get_expression_violinplot()`
- `get_expression_density()`, `get_expression_joyplot()`, `get_expression_raincloud()`
- `get_expression_heatmap()` - Complex heatmaps with k-means clustering
- `get_expression_scatter()`, `get_expression_lineplot()`, `get_expression_lollipop()`

**Comparisons:**
- `get_corr_heatmap()` - Sample correlation heatmaps
- `get_deg_venn_diagram()` - Overlap between comparisons
- `plot_deg_alluvial()` - Alluvial plots for DEG transitions

**Fold-change Analysis:**
- `get_foldchange_scatter()`, `get_foldchange_barplot()`, `get_foldchange_matrix()`
- `get_foldchange_heatmap()`, `get_foldchange_chromosome_plot()`

### 🔬 Functional Enrichment

```r
# MSigDB enrichment (Hallmark, C2, C5, etc.)
msig <- get_msigdb_enrichment(
  vista,
  sample_comparison = "Treatment_VS_Control",
  regulation = "Up",
  msigdb_category = "H",
  from_type = "ENSEMBL"
)

# GO enrichment
go_bp <- get_go_enrichment(vista, regulation = "Up", ont = "BP")
go_mf <- get_go_enrichment(vista, regulation = "Up", ont = "MF")

# KEGG pathway enrichment
kegg <- get_kegg_enrichment(vista, regulation = "Down")

# Gene Set Enrichment Analysis (GSEA)
gsea <- get_gsea(
  vista,
  sample_comparison = "Treatment_VS_Control",
  set_type = "msigdb"
)

# Visualize results
get_enrichment_plot(msig$enrich)
enrichplot::dotplot(go_bp$enrich)
enrichplot::cnetplot(kegg$enrich)
```

### 🎨 Customizable Aesthetics

```r
# Use different color palettes
data("count_data", package = "VISTA")
data("sample_metadata", package = "VISTA")

vista <- create_vista(
  counts = count_data[1:500, ],
  sample_info = sample_metadata,
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control",
  group_palette = "Harmonic"  # Dark 2, Set 2, Warm, Cold, etc.
)

# Access and modify colors
group_colors(vista)
group_palette(vista)
vista <- set_vista_group_colors(
  vista,
  c(control = "#264653", treatment1 = "#E76F51")
)

# Optional: set manual colors for multi-comparison objects
vista_multi <- create_vista(
  counts = count_data[1:500, ],
  sample_info = sample_metadata,
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = c("treatment1", "control"),
  group_denominator = c("control", "treatment1")
)
vista_multi <- set_vista_comparison_colors(
  vista_multi,
  c(treatment1_VS_control = "#6C5CE7", control_VS_treatment1 = "#00A896")
)

# All plots respect the stored color scheme
get_pca_plot(vista, sample_colors = TRUE)
```

### 📦 S4 Class Design

VISTA extends `SummarizedExperiment`, making it compatible with Bioconductor workflows:

```r
# Standard SummarizedExperiment methods work
assay(vista, "norm_counts")[1:5, 1:5]
colData(vista)
rowData(vista)
metadata(vista)

# Plus VISTA-specific accessors
names(comparisons(vista))
deg_summary(vista)
cutoffs(vista)
norm_counts(vista, summarise = TRUE)  # Group-averaged expression
```

## Documentation

- **Vignettes:**
  - [Introduction to VISTA](vignettes/VISTA-airway.Rmd) - Complete workflow tutorial
  - [Color and Palette Design](vignettes/VISTA-colors.Rmd) - Consistent and customizable color workflows
  - [Function Reference](vignettes/VISTA-reference.Rmd) - Comprehensive function guide

- **Help Pages:** `?create_vista`, `?get_pca_plot`, `?get_msigdb_enrichment`, etc.

- **Example Data:** `data("count_data")`, `data("sample_metadata")`

## Automated Report Template

VISTA ships with a YAML-driven report workflow template for publication-ready
single-comparison analysis reports.

```r
# 1) Copy and edit the template
file.copy(
  system.file("reports", "vista-report-template.yml", package = "VISTA"),
  "vista-report.yml"
)

# 2) Render the report
run_vista_report("vista-report.yml")
```

The generated HTML report includes:
- QC plots (PCA, MDS, correlation heatmap)
- Differential expression views (volcano, MA, DEG summaries)
- Expression heatmap and fold-change views
- MSigDB/GO/KEGG enrichment (when available)
- Downloadable tables and plot assets
- Interactive tables with export buttons when `DT` is installed

## Citation

If you use VISTA in published research, please cite:

```
Parsania C (2025). VISTA: Visualization Toolkit for Transcriptomic Analysis.
R package version 0.99.0.
```

## Getting Help

- **Issues:** [GitHub Issues](https://github.com/chiragparsania/VISTA/issues)
- **Questions:** Open a discussion or issue on GitHub
- **Bioconductor Support:** [support.bioconductor.org](https://support.bioconductor.org)

## Development

VISTA is actively developed with:
- ✅ Comprehensive input validation
- ✅ Broad test coverage
- ✅ Complete roxygen2 documentation
- ✅ BiocCheck compliance
- ✅ CI/CD integration ready

### Package Structure

```
VISTA/
├── R/                      # R source code
│   ├── VISTA-class.R      # S4 class definition
│   ├── accessors.R        # Accessor methods
│   ├── rnaseq_related.R   # DE analysis wrappers
│   ├── viz_related.R      # Visualization functions
│   ├── functional_analysis.R  # Enrichment analysis
│   └── utils-internal.R   # Internal utilities
├── tests/testthat/        # Test suite
├── vignettes/             # Package vignettes
├── data/                  # Example datasets
└── man/                   # Documentation
```

## License

GPL-3

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure `R CMD check` passes
5. Submit a pull request

## Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

---

**Built with** ❤️ **for the Bioconductor community**
