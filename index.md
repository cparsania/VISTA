# VISTA

> **V**isualization and **I**ntegrated **S**ystem for **T**ranscriptomic
> **A**nalysis

VISTA is a Bioconductor framework for RNA-seq differential expression
that keeps counts, statistics, annotations, and figures in a single
validated `SummarizedExperiment`-based object — so you can go from raw
counts to a publication-ready narrative without rebuilding the same glue
code each project.

``` r

BiocManager::install("VISTA")
```

- 📦 **One object, one grammar** — `DESeq2`, `edgeR`, `limma-voom`, and
  a DESeq2/edgeR consensus behind a single entry point.
- 🎨 **40+ publication-ready plots** — QC, differential expression,
  expression patterns, fold-change structure, enrichment, and
  deconvolution.
- 🧬 **Bioconductor-native** — extends `SummarizedExperiment`, so
  `assay()`, `rowData()`, `colData()`, and `[` all work as expected.
- 📄 **Reproducible reporting** — export helpers plus a YAML-driven
  Quarto workflow.

------------------------------------------------------------------------

## Contents

[Installation](#installation) · [Quick start](#quick-start) · [Why
VISTA](#why-vista) · [Core workflow](#core-workflow) · [Plot
catalogue](#plot-catalogue) · [Object
design](#bioconductor-compatible-object-design) ·
[Documentation](#documentation) · [Citation](#citation)

------------------------------------------------------------------------

## Installation

**Release** — from Bioconductor:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("VISTA")
```

**Development version** (GitHub):

``` r

if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}

pak::pak("cparsania/VISTA")
```

## Quick start

``` r

library(VISTA)

data("count_data", package = "VISTA")        # gene-by-sample counts
data("sample_metadata", package = "VISTA")   # sample annotations

prepared_counts  <- read_vista_counts(count_data, format = "matrix", gene_id_column = "gene_id")
prepared_samples <- read_vista_metadata(sample_metadata)
matched_inputs   <- match_vista_inputs(prepared_counts, prepared_samples)

vista <- create_vista(
  counts          = matched_inputs$counts,
  sample_info     = matched_inputs$sample_info,
  column_geneid   = matched_inputs$column_geneid,
  group_column    = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control",
  method          = "deseq2",
  log2fc_cutoff   = 1,
  pval_cutoff     = 0.05
)

comp <- names(comparisons(vista))[1]

vista                                                    # object summary
head(comparisons(vista)[[comp]][, c("gene_id", "log2fc", "padj")])
```

Already have a `SummarizedExperiment`? Skip the import step:

``` r

vista <- as_vista(se, group_column = "cond_long")
```

### Explore it

``` r

# Quality control
get_pca_plot(vista, label = TRUE)
get_mds_plot(vista, use_group_colors = TRUE)
get_corr_heatmap(vista)

# Differential expression
get_volcano_plot(vista, sample_comparison = comp)
get_ma_plot(vista, sample_comparison = comp)
get_deg_count_barplot(vista)

# Expression and fold-change views
up_genes <- get_genes_by_regulation(
  vista, sample_comparisons = comp, regulation = "Up", top_n = 40
)[[comp]]

get_expression_heatmap(vista, genes = up_genes, kmeans_k = 3)
get_expression_barplot(vista, genes = up_genes[1:3], by = "sample", facet_by = "gene")
get_foldchange_heatmap(vista)
get_foldchange_lollipop(vista, sample_comparison = comp, genes = up_genes[1:6], facet_by = "gene")
```

## Why VISTA

Most RNA-seq projects repeat the same sequence: normalize counts, fit
models, extract contrasts, plot QC, label genes, summarize pathways,
assemble figures. The friction is rarely the statistics — it is the glue
code between them.

VISTA organizes that work around one validated object:

| Component                      | Accessor                        |
|--------------------------------|---------------------------------|
| Normalized expression          | `assay(x)` / `norm_counts(x)`   |
| Raw filtered counts            | `counts(x)`                     |
| Feature annotations            | `rowData(x)`                    |
| Sample metadata                | `colData(x)` / `sample_info(x)` |
| Differential expression tables | `comparisons(x)`                |
| Analysis parameters            | `cutoffs(x)`                    |

Every accessor and plotting function reads that same object, so you move
from raw counts to a consistent analysis narrative without switching
data structures between steps.

## Core workflow

### 1. Prepare counts and metadata

``` r

prepared_counts  <- read_vista_counts(count_data, format = "matrix", gene_id_column = "gene_id")
prepared_samples <- read_vista_metadata(sample_metadata)
matched_inputs   <- match_vista_inputs(prepared_counts, prepared_samples)
```

These helpers import plain matrices and data frames, featureCounts, STAR
gene counts, HTSeq-count, tximport-like inputs, and RSEM gene results.
No metadata sheet yet? Bootstrap one from the count sample names:

``` r

starter_metadata <- derive_vista_metadata(
  matched_inputs$counts,
  column_geneid = matched_inputs$column_geneid,
  parser = "auto"
)
```

### 2. Build the analysis object

``` r

vista <- create_vista(
  counts = matched_inputs$counts,
  sample_info = matched_inputs$sample_info,
  column_geneid = matched_inputs$column_geneid,
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control",
  method = "limma"                       # or "deseq2", "edger", "both"
)
```

Use `method = "both"` for a DESeq2/edgeR consensus:

``` r

vista_consensus <- create_vista(
  ...,
  method = "both",
  result_source = "consensus"
)

vista_consensus <- set_de_source(vista_consensus, "edger")   # switch the active table
```

### 3. Adjust the model

``` r

# Sample-level covariates
vista_cov <- create_vista(..., covariates = "cell")

# Or a full design formula
vista_design <- create_vista(..., design_formula = ~ cell + cond_long)
```

### 4. Add feature annotations

``` r

vista <- set_rowdata(
  vista,
  orgdb = org.Hs.eg.db,
  columns = c("SYMBOL", "GENENAME", "ENTREZID")
)
```

### 5. Export and report

``` r

export_vista_assets(
  vista,
  out_dir = "vista_assets",
  include_data = c("comparison", "norm_counts", "sample_info")
)

file.copy(
  system.file("reports", "vista-report-template.yml", package = "VISTA"),
  "vista-report.yml"
)
run_vista_report("vista-report.yml")
```

## Plot catalogue

**Quality control and sample structure**

[`get_pca_plot()`](https://cparsania.github.io/VISTA/reference/get_pca_plot.md)
·
[`get_mds_plot()`](https://cparsania.github.io/VISTA/reference/get_mds_plot.md)
·
[`get_umap_plot()`](https://cparsania.github.io/VISTA/reference/get_umap_plot.md)
·
[`get_corr_heatmap()`](https://cparsania.github.io/VISTA/reference/get_corr_heatmap.md)
·
[`get_pairwise_corr_plot()`](https://cparsania.github.io/VISTA/reference/get_pairwise_corr_plot.md)

**Differential expression summaries**

[`get_volcano_plot()`](https://cparsania.github.io/VISTA/reference/get_volcano_plot.md)
·
[`get_ma_plot()`](https://cparsania.github.io/VISTA/reference/get_ma_plot.md)
·
[`get_deg_count_barplot()`](https://cparsania.github.io/VISTA/reference/get_deg_count_barplot.md)
·
[`get_deg_count_pieplot()`](https://cparsania.github.io/VISTA/reference/get_deg_count_pieplot.md)
·
[`get_deg_count_donutplot()`](https://cparsania.github.io/VISTA/reference/get_deg_count_donutplot.md)
·
[`get_deg_venn_diagram()`](https://cparsania.github.io/VISTA/reference/get_deg_venn_diagram.md)
·
[`get_deg_alluvial()`](https://cparsania.github.io/VISTA/reference/get_deg_alluvial.md)

**Expression patterns**

[`get_expression_heatmap()`](https://cparsania.github.io/VISTA/reference/get_expression_heatmap.md)
·
[`get_expression_boxplot()`](https://cparsania.github.io/VISTA/reference/get_expression_boxplot.md)
·
[`get_expression_violinplot()`](https://cparsania.github.io/VISTA/reference/get_expression_violinplot.md)
·
[`get_expression_barplot()`](https://cparsania.github.io/VISTA/reference/get_expression_barplot.md)
·
[`get_expression_lollipop()`](https://cparsania.github.io/VISTA/reference/get_expression_lollipop.md)
·
[`get_expression_scatter()`](https://cparsania.github.io/VISTA/reference/get_expression_scatter.md)
·
[`get_expression_lineplot()`](https://cparsania.github.io/VISTA/reference/get_expression_lineplot.md)
·
[`get_expression_density()`](https://cparsania.github.io/VISTA/reference/get_expression_density.md)
·
[`get_expression_joyplot()`](https://cparsania.github.io/VISTA/reference/get_expression_joyplot.md)
·
[`get_expression_raincloud()`](https://cparsania.github.io/VISTA/reference/get_expression_raincloud.md)
·
[`get_expression_chromosome_plot()`](https://cparsania.github.io/VISTA/reference/get_expression_chromosome_plot.md)
·
[`get_expression_matrix()`](https://cparsania.github.io/VISTA/reference/get_expression_matrix.md)

**Fold-change structure**

[`get_foldchange_scatter()`](https://cparsania.github.io/VISTA/reference/get_foldchange_scatter.md)
·
[`get_foldchange_barplot()`](https://cparsania.github.io/VISTA/reference/get_foldchange_barplot.md)
·
[`get_foldchange_lollipop()`](https://cparsania.github.io/VISTA/reference/get_foldchange_lollipop.md)
·
[`get_foldchange_boxplot()`](https://cparsania.github.io/VISTA/reference/get_foldchange_boxplot.md)
·
[`get_foldchange_lineplot()`](https://cparsania.github.io/VISTA/reference/get_foldchange_lineplot.md)
·
[`get_foldchange_heatmap()`](https://cparsania.github.io/VISTA/reference/get_foldchange_heatmap.md)
·
[`get_foldchange_matrix()`](https://cparsania.github.io/VISTA/reference/get_foldchange_matrix.md)
·
[`get_foldchange_chromosome_plot()`](https://cparsania.github.io/VISTA/reference/get_foldchange_chromosome_plot.md)

**Pathway and enrichment**

[`get_msigdb_enrichment()`](https://cparsania.github.io/VISTA/reference/get_msigdb_enrichment.md)
·
[`get_go_enrichment()`](https://cparsania.github.io/VISTA/reference/get_go_enrichment.md)
·
[`get_kegg_enrichment()`](https://cparsania.github.io/VISTA/reference/get_kegg_enrichment.md)
·
[`get_gsea()`](https://cparsania.github.io/VISTA/reference/get_gsea.md)
·
[`get_enrichment_plot()`](https://cparsania.github.io/VISTA/reference/get_enrichment_plot.md)
·
[`get_enrichment_chord()`](https://cparsania.github.io/VISTA/reference/get_enrichment_chord.md)
·
[`get_pathway_genes()`](https://cparsania.github.io/VISTA/reference/get_pathway_genes.md)
·
[`get_pathway_heatmap()`](https://cparsania.github.io/VISTA/reference/get_pathway_heatmap.md)

**Deconvolution (optional)**

[`run_cell_deconvolution()`](https://cparsania.github.io/VISTA/reference/run_cell_deconvolution.md)
·
[`get_celltype_barplot()`](https://cparsania.github.io/VISTA/reference/get_celltype_barplot.md)
·
[`get_celltype_group_dotplot()`](https://cparsania.github.io/VISTA/reference/get_celltype_group_dotplot.md)
·
[`get_celltype_heatmap()`](https://cparsania.github.io/VISTA/reference/get_celltype_heatmap.md)

## Harmonized plot API

VISTA plotting functions share one argument grammar, so the same concept
always uses the same name across plot families:

| Argument | Controls |
|----|----|
| `sample_group`, `group_column` | Sample filtering and grouping |
| `sample_comparison` / `sample_comparisons` | One contrast / several contrasts |
| `genes`, `top_n` | Which features to show, and how many |
| `by` | Group-level vs sample-level view |
| `facet_by` | Layout by gene, group, comparison, or none |
| `sample_order` | Sample sequencing for per-sample plots |
| `display_id` | User-facing gene labels |
| `summarise` | Collapse replicates to group means |
| `color_by`, `palette`, `colors` | Colour control |
| `return_type` | `"plot"`, `"data"`, or `"both"` |

Older argument names keep working and warn with the release in which
they become defunct — see
[`?"VISTA-deprecated"`](https://cparsania.github.io/VISTA/reference/VISTA-deprecated.md).

## Bioconductor-compatible object design

VISTA extends `SummarizedExperiment`, so standard Bioconductor workflows
apply.

``` r

# Standard Bioconductor access
assay(vista)[1:5, 1:5]
rowData(vista)
colData(vista)
metadata(vista)
vista[1:100, ]                      # subsetting keeps DE tables aligned

# VISTA accessors
comparisons(vista)
deg_summary(vista)
cutoffs(vista)
norm_counts(vista, summarise = TRUE)

# Validation
validate_vista(vista, level = "full")
```

Raw filtered counts are retained alongside the normalized assay, so an
object can go straight back into DESeq2:

``` r

counts(vista)                                     # integer counts
dds <- as_deseq_dataset(vista, design = ~ cond_long)
```

> **Note**
> [`counts()`](https://cparsania.github.io/VISTA/reference/counts.md),
> [`as_deseq_dataset()`](https://cparsania.github.io/VISTA/reference/as_deseq_dataset.md),
> and `[` are available in the development version and are scheduled for
> the next Bioconductor release.

## Example analyses

### Enrichment from a comparison

``` r

msig <- get_msigdb_enrichment(
  vista, sample_comparison = comp, regulation = "Up",
  orgdb = org.Hs.eg.db, species = "Homo sapiens", msigdb_category = "H"
)

go_bp <- get_go_enrichment(
  vista, sample_comparison = comp, regulation = "Up",
  ont = "BP", orgdb = org.Hs.eg.db, species = "Homo sapiens"
)

get_enrichment_plot(msig$enrich)
```

### Consistent colour control

``` r

group_colors(vista)

vista <- set_vista_group_colors(
  vista,
  c(control = "#264653", treatment1 = "#E76F51")
)

vista_consensus <- set_vista_comparison_colors(
  vista_consensus,
  c(treatment1_VS_control = "#6C5CE7")
)
```

## Documentation

📖 **[Package website](https://cparsania.github.io/VISTA/)** ·
**[Function
reference](https://cparsania.github.io/VISTA/reference/index.html)**

**Workflows**

| Article | Description |
|----|----|
| [Complete RNA-seq workflow](https://cparsania.github.io/VISTA/articles/VISTA-airway.html) | End-to-end analysis of the `airway` dataset |
| [DESeq2 vs edgeR](https://cparsania.github.io/VISTA/articles/workflows/VISTA-comparison.html) | Comparing backends and building a consensus |
| [Code economy](https://cparsania.github.io/VISTA/articles/workflows/VISTA-code-economy.html) | VISTA against a standard R workflow |
| [Cell-type deconvolution](https://cparsania.github.io/VISTA/articles/workflows/VISTA-deconvolution.html) | Optional deconvolution workflow |

**Visualization guides**

| Guide | Description |
|----|----|
| [Preparing counts and metadata](https://cparsania.github.io/VISTA/articles/guides/VISTA-input-preparation.html) | Importing common count formats |
| [Colour and palette design](https://cparsania.github.io/VISTA/articles/guides/VISTA-colors.html) | Consistent colours across comparisons |
| [Enrichment chord diagrams](https://cparsania.github.io/VISTA/articles/guides/VISTA-chord.html) | Pathway–gene chord plots |
| [Raincloud plots](https://cparsania.github.io/VISTA/articles/guides/VISTA-raincloud.html) | Distribution views |

Local help works as usual:

``` r

?create_vista
?get_expression_heatmap
?get_go_enrichment
?run_vista_report
```

## Citation

If you use VISTA in published work, cite the release used in your
analysis:

``` r

citation("VISTA")
```

``` text
Parsania C (2026). VISTA: Visualization and Integrated System for Transcriptomic Analysis.
doi:10.18129/B9.bioc.VISTA, https://bioconductor.org/packages/VISTA/
```

## Support

- 🐛 [Report an issue](https://github.com/cparsania/VISTA/issues)
- 💬 [Bioconductor Support](https://support.bioconductor.org/tag/VISTA/)

## Contributing

Contributions should preserve reproducibility and backward
compatibility. Before opening a pull request:

- Add or update tests for functional changes.
- Run `devtools::document()` if roxygen comments changed.
- Run `devtools::test()` and `R CMD check`.
- Update `NEWS.md` for user-visible changes.

## License

[GPL-3](https://www.gnu.org/licenses/gpl-3.0)
