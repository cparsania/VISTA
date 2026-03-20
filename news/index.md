# Changelog

## VISTA 0.99.0

*Submitted to Bioconductor 2026-02-11*

### Overview

VISTA (Visualization Toolkit for Transcriptomic Analysis) provides a
unified S4-based framework for differential expression analysis of
RNA-seq data, wrapping DESeq2 and edgeR workflows with consistent
metadata management and rich visualization capabilities.

### Key Features

#### Core Infrastructure

- S4 `VISTA` class extending `SummarizedExperiment` for standardized
  data management
- Unified differential expression workflow supporting DESeq2 and edgeR
  backends
- Consistent metadata structure for comparisons, cutoffs, and group
  information
- Flexible color palette system for visualizations

#### Visualization Suite (28+ functions)

- **Dimension reduction**: PCA, MDS plots with customizable aesthetics
- **DE results**: Volcano plots, MA plots, DEG count barplots
- **Expression**: Barplots, boxplots, violin plots, density plots,
  joyplots, heatmaps
- **Comparisons**: Venn diagrams, alluvial plots, correlation heatmaps,
  pairwise plots
- **Fold-change**: Scatter plots, barplots, matrix visualizations,
  chromosome plots

#### Functional Analysis

- MSigDB enrichment with flexible ID mapping (SYMBOL, ENSEMBL, ENTREZID)
- GO enrichment analysis (BP, MF, CC ontologies)
- KEGG pathway enrichment
- GSEA support with customizable gene sets
- Integrated visualization functions for enrichment results

#### Optional Features

- Cell-type deconvolution via xCell2 integration
- Automated report generation with Quarto
- Accessor functions for all metadata components

### Implementation Details

- Comprehensive input validation and edge case handling
- Extensive test suite (\>70% coverage)
- Complete roxygen2 documentation with runnable examples
- BiocStyle vignettes demonstrating complete workflows
- Proper namespace management and import declarations

### Bug Fixes

- Fixed contradictory roxygen documentation markers in internal
  utilities
- Added missing `@importFrom` declarations across all modules
- Improved error messages for invalid inputs
- Enhanced edge case handling in visualization functions
- Heatmap utilities now validate non-character `genes` input explicitly,
  support minimal-call defaults, and allow custom colours for
  multi-column annotations
- DEG count pie/donut plots now optionally include non-DE genes as an
  `"Other"` slice and support configurable label text colour
- [`get_genes_by_regulation()`](../reference/get_genes_by_regulation.md)
  now supports top-gene ranking by `abs(log2fc)` and optional annotated
  table output
- PCA/MDS/UMAP plots now accept the standardized `use_group_colors`
  argument while keeping `use_vista_colors` as a deprecated
  compatibility alias

------------------------------------------------------------------------
