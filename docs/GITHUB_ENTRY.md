# VISTA: A Reproducible, End-to-End Framework for RNA-seq Differential Expression and Visualization

**VISTA** (Visualization Toolkit for Transcriptomic Analysis) is an
R/Bioconductor-style package that unifies statistical differential
expression workflows and publication-ready visualization for bulk
RNA-seq data.

It is designed for three audiences: - **Experimental biologists** who
need clear, interpretable outputs without manually stitching multiple
tools. - **Bioinformaticians** who need robust defaults with room for
customization. - **Core facilities and collaborative teams** who need
reproducible analysis objects and consistent reporting.

------------------------------------------------------------------------

## Why VISTA?

RNA-seq analysis often fragments across separate scripts for
normalization, differential testing, plotting, and enrichment. This
increases the risk of inconsistent sample handling, undocumented
parameter drift, and difficult-to-reproduce figures.

VISTA addresses this by: - Wrapping established differential expression
engines (`DESeq2`, `edgeR`) behind one coherent interface. - Storing
counts, sample metadata, gene metadata, results, and cutoffs in a single
`SummarizedExperiment`-based object. - Providing a broad visualization
layer (PCA, MDS, volcano, MA, fold-change, correlation, heatmaps,
expression distributions). - Integrating enrichment workflows (MSigDB,
GO, KEGG, GSEA) through `clusterProfiler` and `msigdbr`.

------------------------------------------------------------------------

## Scientific Foundation

VISTA is built on methods that are widely accepted in transcriptomics: -
**Negative binomial modeling for RNA-seq counts** via `DESeq2` and
`edgeR` for differential expression inference. - **False discovery rate
control** (Benjamini-Hochberg) for multiple testing. - **Gene-set and
pathway enrichment paradigms** (`clusterProfiler`, GSEA-style workflows,
MSigDB gene sets). - **Bioconductor data structures**
(`SummarizedExperiment`) for interoperability and transparent metadata
tracking.

This design choice prioritizes methodological rigor while reducing
user-side engineering overhead.

------------------------------------------------------------------------

## What You Can Do with VISTA

### 1. Build an analysis object once

- Validate count and sample inputs.
- Run DE analysis across one or multiple contrasts.
- Persist DE outputs and summary statistics in the same object.

### 2. Explore quality and structure

- Principal component analysis (`get_pca_plot`)
- Multidimensional scaling (`get_mds_plot`)
- Pairwise and matrix-level sample correlations
  (`get_pairwise_corr_plot`, `get_corr_heatmap`)

### 3. Interpret differential expression

- Volcano and MA plots (`get_volcano_plot`, `get_ma_plot`)
- DEG counts and overlap (`get_deg_count_barplot`,
  `get_deg_venn_diagram`, `plot_deg_alluvial`)
- Fold-change focused visualizations (`get_foldchange_*`)

### 4. Profile expression behavior

- Expression heatmaps (`get_expression_heatmap`)
- Distribution plots (`get_expression_boxplot`,
  `get_expression_violinplot`, `get_expression_density`,
  `get_expression_raincloud`)
- Gene-wise trend and scatter views (`get_expression_lineplot`,
  `get_expression_scatter`, `get_expression_lollipop`)

### 5. Extend to biological interpretation

- MSigDB enrichment (`get_msigdb_enrichment`)
- GO enrichment (`get_go_enrichment`)
- KEGG enrichment (`get_kegg_enrichment`)
- Rank-based gene set enrichment (`get_gsea`)

------------------------------------------------------------------------

## Minimal Reproducible Example

``` r

library(VISTA)

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

# Core outputs
names(comparisons(vista))
deg_summary(vista)

# QC and DE interpretation
get_pca_plot(vista, label_replicates = TRUE)
get_volcano_plot(vista, sample_comparison = names(comparisons(vista))[1])

# Enrichment
msig <- get_msigdb_enrichment(
  vista,
  sample_comparison = names(comparisons(vista))[1],
  regulation = "Up",
  msigdb_category = "H"
)
get_enrichment_plot(msig$enrich)
```

------------------------------------------------------------------------

## Reproducibility and Reporting

VISTA supports reproducible practice by: - Keeping analysis artifacts in
one object (`assays`, `colData`, `rowData`, `metadata`). - Making
parameter choices explicit at object creation and plotting time. -
Providing report generation utilities for consistent communication in
collaborative settings.

For manuscript-grade workflows, users should still pair VISTA with
versioned scripts, locked package environments, and explicit session
information.

------------------------------------------------------------------------

## Installation

``` r

if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
pak::pak("chiragparsania/VISTA")
```

------------------------------------------------------------------------

## Citation and Method References

If VISTA contributes to published work, cite both VISTA and its core
statistical/methodological dependencies.

- Love MI, Huber W, Anders S. **Moderated estimation of fold change and
  dispersion for RNA-seq data with DESeq2.** Genome Biology (2014).\
  <https://doi.org/10.1186/s13059-014-0550-8>
- Robinson MD, McCarthy DJ, Smyth GK. **edgeR: a Bioconductor package
  for differential expression analysis of digital gene expression
  data.** Bioinformatics (2010).\
  <https://doi.org/10.1093/bioinformatics/btp616>
- McCarthy DJ, Chen Y, Smyth GK. **Differential expression analysis of
  multifactor RNA-Seq experiments with respect to biological
  variation.** Nucleic Acids Research (2012).\
  <https://doi.org/10.1093/nar/gks042>
- Benjamini Y, Hochberg Y. **Controlling the false discovery rate: a
  practical and powerful approach to multiple testing.** JRSS-B (1995).\
  <https://doi.org/10.1111/j.2517-6161.1995.tb02031.x>
- Huber W et al. **Orchestrating high-throughput genomic analysis with
  Bioconductor.** Nature Methods (2015).\
  <https://doi.org/10.1038/nmeth.3252>
- Morgan M et al. **SummarizedExperiment package (Bioconductor).**\
  <https://bioconductor.org/packages/SummarizedExperiment>
- Yu G et al. **clusterProfiler: an R package for comparing biological
  themes among gene clusters.** OMICS (2012).\
  <https://doi.org/10.1089/omi.2011.0118>
- Wu T et al. **clusterProfiler 4.0: A universal enrichment tool for
  interpreting omics data.** The Innovation (2021).\
  <https://doi.org/10.1016/j.xinn.2021.100141>
- Liberzon A et al. **The Molecular Signatures Database Hallmark Gene
  Set Collection.** Cell Systems (2015).\
  <https://doi.org/10.1016/j.cels.2015.12.004>
- Subramanian A et al. **Gene set enrichment analysis: a knowledge-based
  approach for interpreting genome-wide expression profiles.** PNAS
  (2005).\
  <https://doi.org/10.1073/pnas.0506580102>

------------------------------------------------------------------------

## Contact

- **Maintainer:** Chirag Parsania\
- **Issues / feature requests:** GitHub Issues in this repository
