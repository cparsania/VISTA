# Complete RNA-seq Analysis Workflow with VISTA

## Introduction

This vignette demonstrates a complete RNA-seq differential expression
workflow using **VISTA** (Visualization and Integrated System for
Transcriptomic Analysis). We’ll use the well-known **airway** dataset
from Bioconductor, which contains RNA-seq data from airway smooth muscle
cells treated with dexamethasone.

### What is VISTA?

VISTA streamlines RNA-seq analysis by:

- **Unifying DE workflows**: Wraps DESeq2 and edgeR with consistent
  output
- **Simplifying visualization**: 28+ publication-ready plotting
  functions
- **Integrating enrichment**: Built-in MSigDB, GO, and KEGG analysis
- **Ensuring reproducibility**: S4 class structure with comprehensive
  metadata

### Dataset Overview

The **airway** dataset (Himes et al. 2014) includes:

- **Samples**: 8 human airway smooth muscle cell lines
- **Treatment**: 4 treated with dexamethasone, 4 untreated
- **Sequencing**: Illumina HiSeq 2000
- **Features**: ~64,000 genes

**Reference**: Himes BE et al. (2014). “RNA-Seq transcriptome profiling
identifies CRISPLD2 as a glucocorticoid responsive gene that modulates
cytokine function in airway smooth muscle cells.” *PLoS One* 9(6):
e99625.

## Installation and Setup

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("VISTA", "airway", "org.Hs.eg.db"))
install.packages("ggplot2")
```

``` r

# Load required packages
library(VISTA)
library(ggplot2)         # For plotting functions
library(airway)          # Dataset
library(org.Hs.eg.db)    # Human gene annotations
library(magrittr)    # %>% 
```

If required packages are missing, install them before rendering this
vignette.

## Data Preparation

### Load the airway dataset

``` r

# Load the SummarizedExperiment object
data("airway", package = "airway")

# Examine the structure
airway
#> class: RangedSummarizedExperiment 
#> dim: 63677 8 
#> metadata(1): ''
#> assays(1): counts
#> rownames(63677): ENSG00000000003 ENSG00000000005 ... ENSG00000273492
#>   ENSG00000273493
#> rowData names(10): gene_id gene_name ... seq_coord_system symbol
#> colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
#> colData names(9): SampleName cell ... Sample BioSample
```

### Extract counts and metadata

``` r

# Extract count matrix
counts_matrix <- assay(airway, "counts")

# Preview counts (first 5 genes, first 4 samples)
counts_matrix[1:5, 1:4]
#>                 SRR1039508 SRR1039509 SRR1039512 SRR1039513
#> ENSG00000000003        679        448        873        408
#> ENSG00000000005          0          0          0          0
#> ENSG00000000419        467        515        621        365
#> ENSG00000000457        260        211        263        164
#> ENSG00000000460         60         55         40         35

# Extract sample metadata
sample_metadata <- as.data.frame(colData(airway))
sample_metadata
#>            SampleName    cell   dex albut        Run avgLength Experiment
#> SRR1039508 GSM1275862  N61311 untrt untrt SRR1039508       126  SRX384345
#> SRR1039509 GSM1275863  N61311   trt untrt SRR1039509       126  SRX384346
#> SRR1039512 GSM1275866 N052611 untrt untrt SRR1039512       126  SRX384349
#> SRR1039513 GSM1275867 N052611   trt untrt SRR1039513        87  SRX384350
#> SRR1039516 GSM1275870 N080611 untrt untrt SRR1039516       120  SRX384353
#> SRR1039517 GSM1275871 N080611   trt untrt SRR1039517       126  SRX384354
#> SRR1039520 GSM1275874 N061011 untrt untrt SRR1039520       101  SRX384357
#> SRR1039521 GSM1275875 N061011   trt untrt SRR1039521        98  SRX384358
#>               Sample    BioSample
#> SRR1039508 SRS508568 SAMN02422669
#> SRR1039509 SRS508567 SAMN02422675
#> SRR1039512 SRS508571 SAMN02422678
#> SRR1039513 SRS508572 SAMN02422670
#> SRR1039516 SRS508575 SAMN02422682
#> SRR1039517 SRS508576 SAMN02422673
#> SRR1039520 SRS508579 SAMN02422683
#> SRR1039521 SRS508580 SAMN02422677

# Simplify treatment labels for clarity
sample_metadata$treatment <- ifelse(
  sample_metadata$dex == "trt",
  "Dexamethasone",
  "Untreated"
)
```

### Prepare data for VISTA

VISTA includes helper functions that standardize counts and sample
metadata before calling
[`create_vista()`](https://cparsania.github.io/VISTA/reference/create_vista.md):

``` r

prepared_counts <- read_vista_counts(
  counts_matrix,
  format = "matrix"
)

prepared_samples <- read_vista_metadata(
  sample_metadata %>%
    tibble::as_tibble() %>%
    dplyr::rename("sample_names" = "Run")
)

matched_inputs <- match_vista_inputs(prepared_counts, prepared_samples)

# Preview the create_vista-ready count table
matched_inputs$counts[1:5, 1:5]
#>                         gene_id SRR1039508 SRR1039509 SRR1039512 SRR1039513
#> ENSG00000000003 ENSG00000000003        679        448        873        408
#> ENSG00000000005 ENSG00000000005          0          0          0          0
#> ENSG00000000419 ENSG00000000419        467        515        621        365
#> ENSG00000000457 ENSG00000000457        260        211        263        164
#> ENSG00000000460 ENSG00000000460         60         55         40         35
```

For a more complete guide covering file-derived sample-name repair and
starter metadata generation with
[`derive_vista_metadata()`](https://cparsania.github.io/VISTA/reference/derive_vista_metadata.md),
see the pkgdown article `Preparing Counts and Metadata for VISTA`.

The matched sample sheet now has stable `sample_names` aligned to the
count columns:

``` r

matched_inputs$sample_info[, c("sample_names", "cell", "treatment", "dex")]
#>            sample_names    cell     treatment   dex
#> SRR1039508   SRR1039508  N61311     Untreated untrt
#> SRR1039509   SRR1039509  N61311 Dexamethasone   trt
#> SRR1039512   SRR1039512 N052611     Untreated untrt
#> SRR1039513   SRR1039513 N052611 Dexamethasone   trt
#> SRR1039516   SRR1039516 N080611     Untreated untrt
#> SRR1039517   SRR1039517 N080611 Dexamethasone   trt
#> SRR1039520   SRR1039520 N061011     Untreated untrt
#> SRR1039521   SRR1039521 N061011 Dexamethasone   trt
```

## Create VISTA Object

### Using DESeq2 backend

The primary method for creating a VISTA object:

``` r

# Create VISTA object with DESeq2 backend
vista <- create_vista(
  counts = matched_inputs$counts,
  sample_info = matched_inputs$sample_info,
  column_geneid = matched_inputs$column_geneid,
  group_column = "treatment",
  group_numerator = "Dexamethasone",
  group_denominator = "Untreated",
  method = "deseq2",
  min_counts = 10,
  min_replicates = 2,
  log2fc_cutoff = 1.0,
  pval_cutoff = 0.05,
  p_value_type = "padj"
)

# Examine the VISTA object
vista
#> class: VISTA 
#> dim: 17199 8 
#> metadata(12): de_results de_summary ... design comparison
#> assays(2): norm_counts counts
#> rownames(17199): ENSG00000000003 ENSG00000000419 ... ENSG00000273487
#>   ENSG00000273488
#> rowData names(1): baseMean
#> colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
#> colData names(11): SampleName cell ... sizeFactor sample_names
#> -------- VISTA --------
#> group column: treatment (Untreated, Dexamethasone)
#> comparisons: Dexamethasone_VS_Untreated
#> DE source: deseq2
#> cutoffs: |log2FC| >= 1, padj <= 0.05
#> raw counts: available via counts()
#> schema: 1.1.0
```

The VISTA object stores:

- **Normalized counts** in `assays`
- **Sample metadata** in `colData`
- **Gene annotations** in `rowData`
- **DE results** in `metadata`

### Validate object integrity

[`create_vista()`](https://cparsania.github.io/VISTA/reference/create_vista.md)
runs validation by default (`validate = TRUE`). You can also run it
explicitly:

``` r

validate_vista(vista, level = "full")
```

For advanced users importing a pre-built `SummarizedExperiment`, use
[`as_vista()`](https://cparsania.github.io/VISTA/reference/as_vista.md)
and then validate:

``` r

se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(norm_counts = norm_counts(vista)),
  colData = S4Vectors::DataFrame(sample_info(vista), row.names = sample_info(vista)$sample_names),
  rowData = S4Vectors::DataFrame(row_data(vista), row.names = rownames(norm_counts(vista)))
)
vista2 <- as_vista(se, group_column = "treatment")
validate_vista(vista2, level = "full")
```

### Alternative: Using edgeR backend

``` r

# Create VISTA object with edgeR backend
vista_edger <- create_vista(
  counts = matched_inputs$counts,
  sample_info = matched_inputs$sample_info,
  column_geneid = matched_inputs$column_geneid,
  group_column = "treatment",
  group_numerator = "Dexamethasone",
  group_denominator = "Untreated",
  method = "edger",  # Use edgeR instead of DESeq2
  min_counts = 10,
  min_replicates = 2,
  log2fc_cutoff = 1.0,
  pval_cutoff = 0.05,
  p_value_type = "padj"
)

vista_edger
#> class: VISTA 
#> dim: 17199 8 
#> metadata(12): de_results de_summary ... design comparison
#> assays(2): norm_counts counts
#> rownames(17199): ENSG00000000003 ENSG00000000419 ... ENSG00000273487
#>   ENSG00000273488
#> rowData names(1): baseMean
#> colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
#> colData names(10): SampleName cell ... treatment sample_names
#> -------- VISTA --------
#> group column: treatment (Untreated, Dexamethasone)
#> comparisons: Dexamethasone_VS_Untreated
#> DE source: edger
#> cutoffs: |log2FC| >= 1, padj <= 0.05
#> raw counts: available via counts()
#> schema: 1.1.0
```

### Alternative: Using limma-voom backend

``` r

# Create VISTA object with limma-voom backend
vista_limma <- create_vista(
  counts = matched_inputs$counts,
  sample_info = matched_inputs$sample_info,
  column_geneid = matched_inputs$column_geneid,
  group_column = "treatment",
  group_numerator = "Dexamethasone",
  group_denominator = "Untreated",
  method = "limma",
  min_counts = 10,
  min_replicates = 2,
  log2fc_cutoff = 1.0,
  pval_cutoff = 0.05,
  p_value_type = "padj"
)

vista_limma
#> class: VISTA 
#> dim: 17199 8 
#> metadata(12): de_results de_summary ... design comparison
#> assays(2): norm_counts counts
#> rownames(17199): ENSG00000000003 ENSG00000000419 ... ENSG00000273487
#>   ENSG00000273488
#> rowData names(1): baseMean
#> colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
#> colData names(10): SampleName cell ... treatment sample_names
#> -------- VISTA --------
#> group column: treatment (Untreated, Dexamethasone)
#> comparisons: Dexamethasone_VS_Untreated
#> DE source: limma
#> cutoffs: |log2FC| >= 1, padj <= 0.05
#> raw counts: available via counts()
#> schema: 1.1.0
```

All three backends share the same filtering rule: a gene is kept when at
least `min_replicates` samples reach `min_counts`.

``` r

# The three backends model the same feature set
data.frame(
  backend = c("deseq2", "edger", "limma"),
  genes   = c(nrow(vista), nrow(vista_edger), nrow(vista_limma)),
  DEGs    = c(
    sum(comparisons(vista)[[1]]$regulation %in% c("Up", "Down")),
    sum(comparisons(vista_edger)[[1]]$regulation %in% c("Up", "Down")),
    sum(comparisons(vista_limma)[[1]]$regulation %in% c("Up", "Down"))
  )
)
#>   backend genes DEGs
#> 1  deseq2 17199  853
#> 2   edger 17199  892
#> 3   limma 17199  660
```

### Advanced: covariates, design formula, and consensus mode

``` r

# Covariate-adjusted model (additive design)
vista_cov <- create_vista(
  counts = matched_inputs$counts,
  sample_info = matched_inputs$sample_info,
  column_geneid = matched_inputs$column_geneid,
  group_column = "treatment",
  group_numerator = "Dexamethasone",
  group_denominator = "Untreated",
  method = "deseq2",
  covariates = c("cell")
)

# Equivalent explicit model formula
vista_formula <- create_vista(
  counts = matched_inputs$counts,
  sample_info = matched_inputs$sample_info,
  column_geneid = matched_inputs$column_geneid,
  group_column = "treatment",
  group_numerator = "Dexamethasone",
  group_denominator = "Untreated",
  method = "deseq2",
  design_formula = "~ cell + treatment"
)

# Run both DESeq2 and edgeR and keep consensus as active source
vista_both <- create_vista(
  counts = matched_inputs$counts,
  sample_info = matched_inputs$sample_info,
  column_geneid = matched_inputs$column_geneid,
  group_column = "treatment",
  group_numerator = "Dexamethasone",
  group_denominator = "Untreated",
  method = "both",
  consensus_mode = "intersection",  # or "union"
  result_source = "consensus"       # or "deseq2"/"edger"
)

# Access source-specific outputs. Each returns one table per contrast; count
# the DEG calls rather than printing 17,000 rows.
vapply(
  c("consensus", "deseq2", "edger"),
  function(src) sum(comparisons(vista_both, source = src)[[1]]$regulation %in% c("Up", "Down")),
  integer(1)
)
#> consensus    deseq2     edger 
#>       845       867       942

# The consensus table keeps each backend's own estimates alongside the
# consensus call, plus a `support` column recording which backends agreed.
table(comparisons(vista_both, source = "consensus")[[1]]$support)
#> 
#>        both deseq2_only  edger_only        none 
#>         845          22          97       17122

# Switch the active source used by plotting functions
vista_both <- set_de_source(vista_both, "edger")
vista_both
#> class: VISTA 
#> dim: 18086 8 
#> metadata(12): de_results de_summary ... design comparison
#> assays(2): norm_counts counts
#> rownames(18086): ENSG00000000003 ENSG00000000419 ... ENSG00000273487
#>   ENSG00000273488
#> rowData names(1): baseMean
#> colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
#> colData names(11): SampleName cell ... sizeFactor sample_names
#> -------- VISTA --------
#> group column: treatment (Untreated, Dexamethasone)
#> comparisons: Dexamethasone_VS_Untreated
#> DE source: edger (of deseq2, edger, consensus)
#> cutoffs: |log2FC| >= 1, padj <= 0.05
#> raw counts: available via counts()
#> schema: 1.1.0
```

### Add gene annotations

Enhance the object with gene symbols and descriptions:

``` r

vista <- set_rowdata(
  vista,
  orgdb = org.Hs.eg.db,
  columns = c("SYMBOL", "GENENAME", "ENTREZID"),
  keytype = "ENSEMBL"
)

# View updated gene annotations
head(rowData(vista))
#> DataFrame with 6 rows and 4 columns
#>                  baseMean      SYMBOL               GENENAME    ENTREZID
#>                 <numeric> <character>            <character> <character>
#> ENSG00000000003  709.7752      TSPAN6          tetraspanin 6        7105
#> ENSG00000000419  521.1362        DPM1 dolichyl-phosphate m..        8813
#> ENSG00000000457  237.5619       SCYL3 SCY1 like pseudokina..       57147
#> ENSG00000000460   58.0338       FIRRM FIGNL1 interacting r..       55732
#> ENSG00000000971 5826.6231         CFH    complement factor H        3075
#> ENSG00000001036 1284.4108       FUCA2   alpha-L-fucosidase 2        2519
```

## Explore the Results

### Access differential expression results

``` r

# Get comparison names
comp_names <- names(comparisons(vista))
comp_names
#> [1] "Dexamethasone_VS_Untreated"

# Get DE results for the comparison
de_results <- comparisons(vista)[[1]]
head(de_results)
#>                         gene_id   baseMean      log2fc      lfcSE       stat
#> ENSG00000000003 ENSG00000000003  709.77518 -0.38027500 0.17108885 -2.2226756
#> ENSG00000000419 ENSG00000000419  521.13616  0.20227695 0.09533472  2.1217555
#> ENSG00000000457 ENSG00000000457  237.56187  0.03272066 0.12278673  0.2664837
#> ENSG00000000460 ENSG00000000460   58.03383 -0.11851609 0.30815777 -0.3845955
#> ENSG00000000971 ENSG00000000971 5826.62312  0.43942357 0.25944519  1.6937048
#> ENSG00000001036 ENSG00000001036 1284.41081 -0.24322719 0.11682555 -2.0819691
#>                     pvalue      padj regulation
#> ENSG00000000003 0.02623769 0.1206592      Other
#> ENSG00000000419 0.03385828 0.1448222      Other
#> ENSG00000000457 0.78986672 0.9229371      Other
#> ENSG00000000460 0.70053714 0.8844114      Other
#> ENSG00000000971 0.09032139 0.2843077      Other
#> ENSG00000001036 0.03734529 0.1541257      Other

# Get DEG summary
deg_summary(vista)
#> $Dexamethasone_VS_Untreated
#>   regulation     n
#> 1       Down   388
#> 2      Other 16346
#> 3         Up   465

# Get analysis cutoffs
cutoffs(vista)
#> $log2fc
#> [1] 1
#> 
#> $pval
#> [1] 0.05
#> 
#> $p_value_type
#> [1] "padj"
#> 
#> $method
#> [1] "deseq2"
#> 
#> $min_counts
#> [1] 10
#> 
#> $min_replicates
#> [1] 2
#> 
#> $covariates
#> character(0)
#> 
#> $design_formula
#> NULL
#> 
#> $consensus_mode
#> NULL
#> 
#> $consensus_log2fc
#> NULL
#> 
#> $active_source
#> [1] "deseq2"
```

### Count significant genes

``` r

# Extract upregulated genes
up_genes <- get_genes_by_regulation(
  vista,
  sample_comparisons = comp_names[1],
  regulation = "Up",
  #top_n = 50
)

# Extract downregulated genes
down_genes <- get_genes_by_regulation(
  vista,
  sample_comparisons = comp_names[1],
  regulation = "Down",
  #top_n = 50
)

# Summary
cat("Upregulated genes:", length(up_genes[[1]]), "\n")
#> Upregulated genes: 465
cat("Downregulated genes:", length(down_genes[[1]]), "\n")
#> Downregulated genes: 388
```

## Quality Control Visualizations

### Sample Correlation Heatmap

Check sample relationships and potential batch effects.

#### Basic correlation heatmap

``` r

  get_corr_heatmap(vista, label_size = 18,base_size = 18, viridis_direction = -1)
```

![](VISTA-airway_files/figure-html/corr-heatmap-basic-1.png)

#### Customize color scheme

``` r

# Reverse viridis color direction
get_corr_heatmap(
  vista,
  viridis_direction = -1,
  viridis_option = "plasma",
  label_size = 18,
  base_size = 18
)
```

![](VISTA-airway_files/figure-html/corr-heatmap-colors-1.png)

#### Show correlation values

``` r

# Display correlation coefficients
get_corr_heatmap(
  vista,
  viridis_direction = -1,
  label = TRUE,
  label_color = 'white',
  viridis_option = "mako",
  label_size = 14,
  base_size = 14
) 
```

![](VISTA-airway_files/figure-html/corr-heatmap-values-1.png)

### Principal Component Analysis (PCA)

Visualize sample clustering and variation.

#### Basic PCA with labels

``` r

get_pca_plot(
  vista,
  label = TRUE,label_size = 5
)
```

![](VISTA-airway_files/figure-html/pca-basic-1.png)

#### PCA colored by different metadata

``` r

# Shape points by cell line
get_pca_plot(
  vista,
  label = TRUE,
  label_size = 5,
  shape_by = "cell"
)
```

![](VISTA-airway_files/figure-html/pca-shape-1.png)

#### PCA with top variable genes

``` r

# Use top 500 most variable genes
get_pca_plot(
  vista,
  top_n = 500,
  show_clusters = TRUE,
  shape_by = "cell"
)
```

![](VISTA-airway_files/figure-html/pca-top-genes-1.png)

#### PCA with custom circle size

``` r

# Larger points for better visibility
get_pca_plot(
  vista,
  label = TRUE,
  point_size = 15,label_size = 5
)
```

![](VISTA-airway_files/figure-html/pca-circle-size-1.png)

#### PCA without labels

``` r

# Clean plot without sample labels
get_pca_plot(
  vista,
  label = FALSE,
  point_size = 12
)
```

![](VISTA-airway_files/figure-html/pca-no-labels-1.png)

### Multidimensional Scaling (MDS)

Alternative dimensionality reduction method.

#### Basic MDS plot

``` r

get_mds_plot(
  vista,
  label = TRUE
)
```

![](VISTA-airway_files/figure-html/mds-basic-1.png)

#### MDS with top variable genes

``` r

get_mds_plot(
  vista,
  top_n = 500,
  label = TRUE
)
```

![](VISTA-airway_files/figure-html/mds-top-genes-1.png)

#### MDS with custom shapes

``` r

# Shape points by cell line
get_mds_plot(
  vista,
  label = TRUE,
  shape_by = "cell"
)
```

![](VISTA-airway_files/figure-html/mds-shapes-1.png)

### Uniform Manifold Approximation and Projection (UMAP)

Non-linear sample embedding for exploratory structure.

#### Basic UMAP plot

``` r

get_umap_plot(
  vista,
  label = TRUE
)
```

![](VISTA-airway_files/figure-html/umap-basic-1.png)

#### UMAP colored by a user-defined metadata column

``` r

get_umap_plot(
  vista,
  color_by = "cell",
  shape_by = "treatment",
  label = TRUE
)
```

![](VISTA-airway_files/figure-html/umap-color-by-cell-1.png)

## Differential Expression Visualizations

### DEG Count Summary

#### Basic count barplot

``` r

get_deg_count_barplot(vista)
```

![](VISTA-airway_files/figure-html/deg-count-basic-1.png)

#### Faceted by regulation

``` r

get_deg_count_barplot(
  vista,
  facet_by = "regulation"
)
```

![](VISTA-airway_files/figure-html/deg-count-facet-1.png)

### Volcano Plot

Classic volcano plot showing log2FC against significance. By default the
axis and thresholds are inherited from the cutoffs the object was built
with, so the genes coloured here are exactly the ones
[`deg_summary()`](https://cparsania.github.io/VISTA/reference/VISTA-accessors.md)
reports. Because this object used `p_value_type = "padj"`, the y-axis is
-log10(adjusted p-value); pass `p_value_type = "pvalue"` to plot raw
p-values instead.

#### Basic volcano plot

``` r

get_volcano_plot(
  vista,
  sample_comparison = comp_names[1]
)
```

![](VISTA-airway_files/figure-html/volcano-basic-1.png)

#### Customize cutoffs and labels

``` r

get_volcano_plot(
  vista,
  sample_comparison = comp_names[1],
  log2fc_cutoff = 1.5,
  pval_cutoff = 0.01,
  label_size = 3,
  display_id = "SYMBOL"
)
```

![](VISTA-airway_files/figure-html/volcano-custom-1.png)

#### Custom colors

``` r

# Custom colors for up/down regulated genes
get_volcano_plot(
  vista,
  sample_comparison = comp_names[1],
  colors = c(Up = "red", Down = "blue", Other = "grey80"),
  display_id = "SYMBOL"
)
```

![](VISTA-airway_files/figure-html/volcano-colors-1.png)

### MA Plot

Mean expression vs log2 fold-change.

#### Basic MA plot

``` r

get_ma_plot(
  vista,
  sample_comparison = comp_names[1]
)
```

![](VISTA-airway_files/figure-html/ma-basic-1.png)

#### Label top genes

``` r

get_ma_plot(
  vista,
  sample_comparison = comp_names[1],
  label_n = 10,
  point_size = 2,
  display_id = "SYMBOL"
)
```

![](VISTA-airway_files/figure-html/ma-labeled-1.png)

#### Custom cutoffs

``` r

get_ma_plot(
  vista,
  sample_comparison = comp_names[1],
  label_n = 5,
  display_id = "SYMBOL",
  point_size = 1.5,
  alpha = 0.8
)
```

![](VISTA-airway_files/figure-html/ma-custom-1.png)

## Expression Pattern Analysis

### Prepare gene sets

``` r

# Get top 50 DEGs by adjusted p-value
de_table <- comparisons(vista)[[1]]

top_degs <- get_genes_by_regulation(vista, names(comparisons(vista))[1],regulation = "Both",top_n = 50)[[1]]  # top 50 by abs fold change 

top_up <- get_genes_by_regulation(vista, names(comparisons(vista))[1],regulation = "Up",top_n = 50)[[1]]  

top_down <- get_genes_by_regulation(vista, names(comparisons(vista))[1],regulation = "Down",top_n = 50)[[1]]  


cat("Top upregulated genes:\n")
#> Top upregulated genes:
print(top_up)
#>  [1] "ENSG00000179593" "ENSG00000109906" "ENSG00000250978" "ENSG00000132518"
#>  [5] "ENSG00000171819" "ENSG00000127954" "ENSG00000249364" "ENSG00000137673"
#>  [9] "ENSG00000100033" "ENSG00000168481" "ENSG00000168309" "ENSG00000264868"
#> [13] "ENSG00000152583" "ENSG00000163884" "ENSG00000177575" "ENSG00000127324"
#> [17] "ENSG00000101342" "ENSG00000270689" "ENSG00000268894" "ENSG00000152779"
#> [21] "ENSG00000128045" "ENSG00000096060" "ENSG00000152463" "ENSG00000173838"
#> [25] "ENSG00000273259" "ENSG00000101347" "ENSG00000187288" "ENSG00000219565"
#> [29] "ENSG00000211445" "ENSG00000143127" "ENSG00000128917" "ENSG00000170214"
#> [33] "ENSG00000163083" "ENSG00000178723" "ENSG00000248187" "ENSG00000157152"
#> [37] "ENSG00000170323" "ENSG00000231246" "ENSG00000233117" "ENSG00000157514"
#> [41] "ENSG00000189221" "ENSG00000165995" "ENSG00000182836" "ENSG00000112936"
#> [45] "ENSG00000269289" "ENSG00000174697" "ENSG00000179094" "ENSG00000187193"
#> [49] "ENSG00000006788" "ENSG00000102760"
cat("\nTop downregulated genes:\n")
#> 
#> Top downregulated genes:
print(top_down)
#>  [1] "ENSG00000128285" "ENSG00000267339" "ENSG00000019186" "ENSG00000183454"
#>  [5] "ENSG00000146006" "ENSG00000122679" "ENSG00000155897" "ENSG00000143494"
#>  [9] "ENSG00000141469" "ENSG00000108700" "ENSG00000162692" "ENSG00000175489"
#> [13] "ENSG00000183092" "ENSG00000250657" "ENSG00000136267" "ENSG00000214814"
#> [17] "ENSG00000261121" "ENSG00000105989" "ENSG00000122877" "ENSG00000188176"
#> [21] "ENSG00000131771" "ENSG00000165272" "ENSG00000184564" "ENSG00000079101"
#> [25] "ENSG00000188501" "ENSG00000119714" "ENSG00000223811" "ENSG00000130487"
#> [29] "ENSG00000166670" "ENSG00000165388" "ENSG00000013293" "ENSG00000123405"
#> [33] "ENSG00000145777" "ENSG00000140600" "ENSG00000124134" "ENSG00000146250"
#> [37] "ENSG00000116991" "ENSG00000126878" "ENSG00000197046" "ENSG00000128165"
#> [41] "ENSG00000084710" "ENSG00000173110" "ENSG00000123689" "ENSG00000106003"
#> [45] "ENSG00000181634" "ENSG00000154864" "ENSG00000182732" "ENSG00000136999"
#> [49] "ENSG00000015520" "ENSG00000095585"
```

### Expression Heatmaps

#### Basic heatmap

``` r

get_expression_heatmap(vista)
```

![](VISTA-airway_files/figure-html/heatmap-basic-1.png)

#### Heatmap with explicit gene set

``` r

get_expression_heatmap(
  vista,
  sample_group = levels(colData(vista)$treatment),
  genes = top_degs,
  display_id = "SYMBOL",
  show_row_names = TRUE
)
```

![](VISTA-airway_files/figure-html/heatmap-explicit-1.png)

#### Heatmap with k-means clustering

``` r

get_expression_heatmap(
  vista,
  sample_group = levels(colData(vista)$treatment),
  genes = top_degs,
  display_id = "SYMBOL",
  kmeans_k = 3,
  show_row_names = TRUE
)
```

![](VISTA-airway_files/figure-html/heatmap-kmeans-1.png)

#### Heatmap with column annotations

``` r

get_expression_heatmap(
  vista,
  sample_group = levels(colData(vista)$treatment),
  genes = top_degs,
  show_row_names = TRUE,
  display_id = "SYMBOL",
  kmeans_k = 3,
  cluster_row_slices = FALSE,
  summarise_replicates = FALSE,
  annotate_columns = TRUE
)
```

![](VISTA-airway_files/figure-html/heatmap-annotated-1.png)

#### Heatmap with multiple column annotations and `cluster_by`

``` r

# Use multiple sample-level columns in top annotation.
# By default, columns are split by the first annotation column.


get_expression_heatmap(
  vista,
  sample_group = levels(colData(vista)$treatment),
  genes = top_degs,
  show_row_names = FALSE,
  display_id = "SYMBOL",
  summarise_replicates = FALSE,
  annotate_columns = c("treatment", "cell"),
  cluster_by = "cell"
)
```

![](VISTA-airway_files/figure-html/heatmap-annotated-multi-1.png)

#### Heatmap showing each replicate

``` r

get_expression_heatmap(
  vista,
  sample_group = levels(colData(vista)$treatment),
  genes = top_degs,
  display_id = "SYMBOL",
  summarise_replicates = FALSE,
  show_row_names = FALSE,
  annotate_columns = TRUE,
  kmeans_k = 2
)
```

![](VISTA-airway_files/figure-html/heatmap-summarized-1.png)

### Expression Barplots

#### Basic barplot

``` r

get_expression_barplot(
  vista,
  genes = top_up[1:4],
  display_id = "SYMBOL",
)+theme_minimal(base_size = 16)
```

![](VISTA-airway_files/figure-html/barplot-basic-1.png)

#### Log-transformed with statistics

``` r

# Add statistical comparisons between groups
get_expression_barplot(
  vista,
  genes = top_up[1:4],
  display_id = "SYMBOL",
  log_transform = TRUE,
  stats_group = TRUE  # Enable statistical annotations
) + theme_minimal(base_size = 16)
```

![](VISTA-airway_files/figure-html/barplot-log-stats-1.png)

#### Per-sample barplot for selected genes

``` r

get_expression_barplot(
  vista,
  genes = top_up[1:2],
  display_id = "SYMBOL",
  by = "sample",
  sample_order = "group"
)+ theme(text = element_text(size = 16))
```

![](VISTA-airway_files/figure-html/barplot-per-sample-1.png)

#### Compare up and down regulated genes

``` r

# Compare expression of both up- and down-regulated genes
selected_genes <- c(top_up[1:3], top_down[1:3])
get_expression_barplot(
  vista,
  genes = selected_genes,
  display_id = "SYMBOL",
  log_transform = TRUE,
  stats_group = TRUE
)+ theme(text = element_text(size = 16))
```

![](VISTA-airway_files/figure-html/barplot-comparison-1.png)

### Expression Boxplots

#### Basic boxplot

``` r

get_expression_boxplot(
  vista,
  genes = top_up[1:4],
  display_id = "SYMBOL"
)+ theme(text = element_text(size = 16))
```

![](VISTA-airway_files/figure-html/boxplot-basic-1.png)

#### Boxplot without faceting

``` r

# All genes overlaid on same plot
get_expression_boxplot(
  vista,
  genes = top_up[1:3],
  display_id = "SYMBOL",
  facet_by = "none"
) + theme(text = element_text(size = 16))
```

![](VISTA-airway_files/figure-html/boxplot-no-facet-1.png)

#### Boxplot with faceting by gene

``` r

# Each gene in separate panel - must specify facet_by = "gene"
get_expression_boxplot(
  vista,
  genes = top_up[1:3],
  display_id = "SYMBOL",
  facet_by = "gene",
  facet_scales = "free_y"
) + theme(text = element_text(size = 16))
```

![](VISTA-airway_files/figure-html/boxplot-facets-1.png)

#### Boxplot with gene facets AND statistics

``` r

# Each gene in separate panel WITH statistical comparisons
get_expression_boxplot(
  vista,
  genes = top_up[1:4],
  display_id = "SYMBOL",
  log_transform = TRUE,
  facet_by = "gene",
  facet_scales = "free_y",
  stats_group = TRUE,  # Add statistics to each gene panel
  p.label = "p.signif"
)+ theme(text = element_text(size = 16))
```

![](VISTA-airway_files/figure-html/boxplot-facets-stats-1.png)

#### Pooled genes with statistics

``` r

# Pool all genes together for group comparison with statistical test
get_expression_boxplot(
  vista,
  genes = top_up[1:5],
  display_id = "SYMBOL",
  log_transform = TRUE,
  pool_genes = TRUE,
  facet_by = "none",
  stats_group = TRUE,  # Required for statistical annotations
  p.label = "p.signif"
)+ theme(text = element_text(size = 16))
```

![](VISTA-airway_files/figure-html/boxplot-pooled-1.png)

#### Log-transformed with p-values

``` r

# Show statistical comparisons between treatment groups
get_expression_boxplot(
  vista,
  genes = top_up[1:4],
  display_id = "SYMBOL",
  log_transform = TRUE,
  stats_group = TRUE,  # Enable statistical annotations
  p.label = "p.signif"
)+ theme(text = element_text(size = 16))
```

![](VISTA-airway_files/figure-html/boxplot-pvalues-1.png)

### Expression Violin Plots

#### Basic violin plot

``` r

get_expression_violinplot(
  vista,
  genes = top_up[1:4],
  display_id = "SYMBOL"
)+ theme(text = element_text(size = 16))
```

![](VISTA-airway_files/figure-html/violin-basic-1.png)

#### Violin with log2 transformation

``` r

get_expression_violinplot(
  vista,
  genes = top_up[1:4],,
  display_id = "SYMBOL",
  value_transform = "none"
)+ theme(text = element_text(size = 16))
```

![](VISTA-airway_files/figure-html/violin-log-1.png)

#### Violin with z-score transformation

``` r

get_expression_violinplot(
  vista,
  genes = top_up[1:4],
  value_transform = "zscore"
)+ theme(text = element_text(size = 16))
```

![](VISTA-airway_files/figure-html/violin-zscore-1.png)

### Additional Expression Plots

#### Density plot

``` r

get_expression_density(
  vista,
  genes = top_up[1:50],
  log_transform = TRUE
)+ theme(text = element_text(size = 16))
```

![](VISTA-airway_files/figure-html/density-plot-1.png)

#### Scatter plot (sample vs sample)

``` r

# Compare two samples
samples <- colnames(vista)
get_expression_scatter(
  vista,
  sample_x = samples[1],
  sample_y = samples[2],
  log_transform = TRUE,
  label_n = 50,display_id = "SYMBOL",label_size = 4
)+ theme(text = element_text(size = 16))
```

![](VISTA-airway_files/figure-html/scatter-plot-1.png)

#### Line plot (expression across samples)

``` r

get_expression_lineplot(
  vista,
  genes = top_up[1:3],
  log_transform = TRUE,display_id = "SYMBOL",
  by = "sample",facet_by = "none",
  group_column = "treatment",sample_group = c("Untreated","Dexamethasone")
)
```

![](VISTA-airway_files/figure-html/line-plot-1.png)

#### Lollipop plot

``` r

get_expression_lollipop(
  vista,
  genes = top_up[1:4],
  display_id = "SYMBOL",
  log_transform = TRUE
)
```

![](VISTA-airway_files/figure-html/lollipop-plot-1.png)

#### Per-sample lollipop plot

``` r

get_expression_lollipop(
  vista,
  genes = top_up[1:2],
  display_id = "SYMBOL",
  by = "sample",
  sample_order = "expression"
)
```

![](VISTA-airway_files/figure-html/lollipop-per-sample-1.png)

#### Joyplot by treatment group

``` r

# Ridges by treatment group - shows distribution for each group
get_expression_joyplot(
  vista,
  genes = top_up[1:5],
  log_transform = TRUE,
  y_by = "group",      # Each treatment group gets a ridge
  color_by = "group"   # Color by treatment group
)
```

![](VISTA-airway_files/figure-html/joyplot-1.png)

#### Joyplot by sample

``` r

# Ridges by individual sample - shows distribution for each sample
get_expression_joyplot(
  vista,
  genes = top_up[1:3],
  log_transform = TRUE,
  y_by = "sample",     # Each sample gets a ridge
  color_by = "group"   # Color by treatment group
)
```

![](VISTA-airway_files/figure-html/joyplot-sample-1.png)

#### Raincloud plot (expression)

``` r

get_expression_raincloud(
  vista,
  genes = top_up,
  value_transform = "log2",
  summarise = TRUE,
  facet_by = "none",
  id.long.var = "gene",
  stats_group = TRUE
)
```

![](VISTA-airway_files/figure-html/raincloud-expression-1.png)

## Functional Enrichment Analysis

### MSigDB Enrichment

#### Hallmark gene sets - Upregulated

``` r

msig_up <- get_msigdb_enrichment(
  vista,
  sample_comparison = comp_names[1],
  regulation = "Up",
  msigdb_category = "H",  # Hallmark gene sets
  species = "Homo sapiens",
  from_type = "ENSEMBL"
)

# View top enriched pathways
if (!is.null(msig_up$enrich) && nrow(msig_up$enrich@result) > 0) {
  head(msig_up$enrich@result[, c("Description", "pvalue", "p.adjust", "Count")])
}
#>                                                                           Description
#> HALLMARK_TNFA_SIGNALING_VIA_NFKB                     HALLMARK_TNFA_SIGNALING_VIA_NFKB
#> HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
#> HALLMARK_ADIPOGENESIS                                           HALLMARK_ADIPOGENESIS
#> HALLMARK_UV_RESPONSE_DN                                       HALLMARK_UV_RESPONSE_DN
#> HALLMARK_HYPOXIA                                                     HALLMARK_HYPOXIA
#> HALLMARK_P53_PATHWAY                                             HALLMARK_P53_PATHWAY
#>                                                  pvalue     p.adjust Count
#> HALLMARK_TNFA_SIGNALING_VIA_NFKB           3.153757e-13 1.545341e-11    28
#> HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION 1.136118e-06 2.783488e-05    19
#> HALLMARK_ADIPOGENESIS                      2.272727e-04 3.712120e-03    15
#> HALLMARK_UV_RESPONSE_DN                    3.664291e-04 4.488757e-03    12
#> HALLMARK_HYPOXIA                           7.260526e-04 5.929430e-03    14
#> HALLMARK_P53_PATHWAY                       7.260526e-04 5.929430e-03    14
```

#### Hallmark gene sets - Downregulated

``` r

msig_down <- get_msigdb_enrichment(
  vista,
  sample_comparison = comp_names[1],
  regulation = "Down",
  msigdb_category = "H",
  species = "Homo sapiens",
  from_type = "ENSEMBL"
)

if (!is.null(msig_down$enrich) && nrow(msig_down$enrich@result) > 0) {
  head(msig_down$enrich@result[, c("Description", "pvalue", "p.adjust", "Count")])
}
#>                                                                           Description
#> HALLMARK_P53_PATHWAY                                             HALLMARK_P53_PATHWAY
#> HALLMARK_TNFA_SIGNALING_VIA_NFKB                     HALLMARK_TNFA_SIGNALING_VIA_NFKB
#> HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
#> HALLMARK_MTORC1_SIGNALING                                   HALLMARK_MTORC1_SIGNALING
#> HALLMARK_MYOGENESIS                                               HALLMARK_MYOGENESIS
#> HALLMARK_APOPTOSIS                                                 HALLMARK_APOPTOSIS
#>                                                  pvalue     p.adjust Count
#> HALLMARK_P53_PATHWAY                       1.669778e-06 7.180044e-05    17
#> HALLMARK_TNFA_SIGNALING_VIA_NFKB           4.166736e-04 8.958483e-03    13
#> HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION 1.377933e-03 1.975037e-02    12
#> HALLMARK_MTORC1_SIGNALING                  4.197724e-03 3.610043e-02    11
#> HALLMARK_MYOGENESIS                        4.197724e-03 3.610043e-02    11
#> HALLMARK_APOPTOSIS                         8.312106e-03 5.673959e-02     9
```

### Enrichment Visualizations

#### VISTA dotplot (default)

``` r

# VISTA's wrapper function
if (!is.null(msig_up$enrich) && nrow(msig_up$enrich@result) > 0) {
  get_enrichment_plot(msig_up$enrich, top_n = 10)
}
```

![](VISTA-airway_files/figure-html/enrichment-dotplot-1.png)

#### Barplot (clusterProfiler native)

``` r

# Use generic barplot with enrichResult method
if (!is.null(msig_up$enrich) && nrow(msig_up$enrich@result) > 0) {
  barplot(msig_up$enrich, showCategory = 10)
}
```

![](VISTA-airway_files/figure-html/enrichment-barplot-1.png)

#### Dotplot with customization

``` r

# Customized dotplot with more categories
if (!is.null(msig_up$enrich) && nrow(msig_up$enrich@result) > 0) {
  enrichplot::dotplot(msig_up$enrich, showCategory = 20, font.size = 12)
}
```

![](VISTA-airway_files/figure-html/enrichment-dotplot-native-1.png)

#### Network plot (clusterProfiler native)

``` r

# Gene-concept network showing gene-pathway relationships
if (!is.null(msig_up$enrich) && nrow(msig_up$enrich@result) > 0) {
  enrichplot::cnetplot(msig_up$enrich, showCategory = 5)
}
```

![](VISTA-airway_files/figure-html/enrichment-network-1.png)

#### Chord diagram (gene–pathway relationships)

The chord diagram reveals which **hub genes** drive multiple enriched
pathways and how much redundancy exists across terms. Chords can be
coloured by fold-change when a VISTA object is supplied.

``` r

# Pathway-coloured chord diagram (no VISTA object needed)
if (!is.null(msig_up$enrich) && nrow(msig_up$enrich@result) > 0) {
  get_enrichment_chord(
    msig_up,
    top_n    = 6,
    color_by = "pathway",
    title    = "Gene-Pathway Membership (Up-regulated, Hallmark)"
  )
}
```

![](VISTA-airway_files/figure-html/enrichment-chord-pathway-1.png)

``` r

# Fold-change coloured chords
if (!is.null(msig_up$enrich) && nrow(msig_up$enrich@result) > 0) {
  get_enrichment_chord(
    msig_up,
    vista  = vista,
    sample_comparison = names(comparisons(vista))[1],
    top_n      = 6,
    color_by   = "foldchange",
    display_id = "SYMBOL",
    title      = "Gene-Pathway Chord (coloured by log2FC)"
  )
}
```

![](VISTA-airway_files/figure-html/enrichment-chord-fc-1.png)

``` r

# Show only hub genes shared across 2+ pathways
if (!is.null(msig_up$enrich) && nrow(msig_up$enrich@result) > 0) {
  chord_result <- get_enrichment_chord(
    msig_up,
    vista    = vista,
    sample_comparison   = names(comparisons(vista))[1],
    top_n        = 8,
    min_pathways = 2,
    color_by     = "regulation",
    display_id   = "SYMBOL",
    title        = "Hub Genes Bridging Multiple Pathways"
  )

  # Inspect hub genes returned invisibly
  if (length(chord_result$hub_genes)) {
    cat("Hub genes:", paste(head(chord_result$hub_genes, 10), collapse = ", "), "\n")
  }
}
```

![](VISTA-airway_files/figure-html/enrichment-chord-hub-1.png)

    #> Hub genes: ENSG00000003402, ENSG00000060718, ENSG00000067082, ENSG00000099860, ENSG00000102804, ENSG00000118689, ENSG00000120129, ENSG00000122641, ENSG00000130066, ENSG00000131979

### Pathway-Specific Expression Heatmaps

Use enrichment output to extract pathway genes and visualize their
expression directly.

#### Extract genes from top pathways

``` r

if (!is.null(msig_up$enrich) && nrow(msig_up$enrich@result) > 0) {
  pathway_gene_list <- get_pathway_genes(
    msig_up$enrich,
    top_n = 3,
    return_type = "list"
  )

  # Preview a few genes per pathway. Pathway membership is a SET -- the order
  # clusterProfiler returns genes in is not meaningful and is not stable across
  # R sessions -- so sort before taking a head().
  lapply(pathway_gene_list, function(g) head(sort(g), 5))
}
#> $HALLMARK_TNFA_SIGNALING_VIA_NFKB
#> [1] "ENSG00000003402" "ENSG00000067082" "ENSG00000099860" "ENSG00000102804"
#> [5] "ENSG00000105835"
#> 
#> $HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
#> [1] "ENSG00000060718" "ENSG00000070404" "ENSG00000099860" "ENSG00000105664"
#> [5] "ENSG00000107796"
#> 
#> $HALLMARK_ADIPOGENESIS
#> [1] "ENSG00000079435" "ENSG00000095637" "ENSG00000127083" "ENSG00000128311"
#> [5] "ENSG00000132170"
```

#### Heatmap of genes from top enriched pathways

``` r

if (!is.null(msig_up$enrich) && nrow(msig_up$enrich@result) > 0) {
  get_pathway_heatmap(
    x = vista,
    enrichment = msig_up$enrich,
    sample_group = c("Untreated", "Dexamethasone"),
    top_n = 2,
    gene_mode = "union",
    max_genes = 60,
    value_transform = "zscore",
    display_id = "SYMBOL",
    annotate_columns = TRUE,
    summarise_replicates = FALSE,
    show_row_names = FALSE
  )
}
```

![](VISTA-airway_files/figure-html/pathway-heatmap-1.png)

### GO Enrichment

#### Biological Process

``` r

go_bp <- get_go_enrichment(
  vista,
  sample_comparison = comp_names[1],
  regulation = "Up",
  ont = "BP",  # Biological Process
  species = "Homo sapiens",
  from_type = "ENSEMBL"
)

if (!is.null(go_bp$enrich) && nrow(go_bp$enrich@result) > 0) {
  head(go_bp$enrich@result[, c("Description", "pvalue", "p.adjust", "Count")], n = 10)
}
#>                                                                                    Description
#> GO:0030198                                                   extracellular matrix organization
#> GO:0043062                                                extracellular structure organization
#> GO:0045229                                       external encapsulating structure organization
#> GO:0060986                                                         endocrine hormone secretion
#> GO:0032970                                          regulation of actin filament-based process
#> GO:0071375                                       cellular response to peptide hormone stimulus
#> GO:0050886                                                                   endocrine process
#> GO:0003012                                                               muscle system process
#> GO:0035360 positive regulation of peroxisome proliferator activated receptor signaling pathway
#> GO:0043500                                                                   muscle adaptation
#>                  pvalue     p.adjust Count
#> GO:0030198 3.993655e-07 0.0005015893    20
#> GO:0043062 4.228448e-07 0.0005015893    20
#> GO:0045229 4.475811e-07 0.0005015893    20
#> GO:0060986 3.589831e-06 0.0025081955     8
#> GO:0032970 4.017769e-06 0.0025081955    22
#> GO:0071375 4.476256e-06 0.0025081955    20
#> GO:0050886 7.259704e-06 0.0034867322    10
#> GO:0003012 1.573717e-05 0.0066135461    24
#> GO:0035360 2.458483e-05 0.0091838006     4
#> GO:0043500 3.708999e-05 0.0121426416    10
```

#### GO Visualization

``` r

if (!is.null(go_bp$enrich) && nrow(go_bp$enrich@result) > 0) {
  get_enrichment_plot(go_bp$enrich, top_n = 20)
}
```

![](VISTA-airway_files/figure-html/go-plot-1.png)

### Gene Set Enrichment Analysis (GSEA)

GSEA uses ranked gene lists to identify pathways enriched at the top or
bottom of the ranking. VISTA automatically prepares the ranked list from
your differential expression results.

#### GSEA with MSigDB Hallmark gene sets

``` r

# Run GSEA using VISTA's native function
set.seed(20260101)
gsea_results <- get_gsea(
  vista,
  sample_comparison = comp_names[1],
  set_type = "msigdb",
  from_type = "ENSEMBL",
  species = "Homo sapiens",
  msigdb_category = "H",  # Hallmark gene sets
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# Show results
if (!is.null(gsea_results$enrich) && nrow(gsea_results$enrich@result) > 0) {
  head(gsea_results$enrich@result[, c("Description", "NES", "pvalue", "p.adjust")], n = 10)
}
#>                                                                       Description
#> HALLMARK_ADIPOGENESIS                                       HALLMARK_ADIPOGENESIS
#> HALLMARK_ANDROGEN_RESPONSE                             HALLMARK_ANDROGEN_RESPONSE
#> HALLMARK_TNFA_SIGNALING_VIA_NFKB                 HALLMARK_TNFA_SIGNALING_VIA_NFKB
#> HALLMARK_XENOBIOTIC_METABOLISM                     HALLMARK_XENOBIOTIC_METABOLISM
#> HALLMARK_OXIDATIVE_PHOSPHORYLATION             HALLMARK_OXIDATIVE_PHOSPHORYLATION
#> HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY
#> HALLMARK_APICAL_JUNCTION                                 HALLMARK_APICAL_JUNCTION
#> HALLMARK_BILE_ACID_METABOLISM                       HALLMARK_BILE_ACID_METABOLISM
#> HALLMARK_IL2_STAT5_SIGNALING                         HALLMARK_IL2_STAT5_SIGNALING
#> HALLMARK_HYPOXIA                                                 HALLMARK_HYPOXIA
#>                                               NES       pvalue     p.adjust
#> HALLMARK_ADIPOGENESIS                    2.011813 1.210205e-07 6.051025e-06
#> HALLMARK_ANDROGEN_RESPONSE               1.654001 1.601434e-03 2.669057e-02
#> HALLMARK_TNFA_SIGNALING_VIA_NFKB         1.615979 6.403582e-04 1.600895e-02
#> HALLMARK_XENOBIOTIC_METABOLISM           1.507815 3.901630e-03 3.901630e-02
#> HALLMARK_OXIDATIVE_PHOSPHORYLATION       1.491037 3.071152e-03 3.838940e-02
#> HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY 1.483494 4.844167e-02 1.513802e-01
#> HALLMARK_APICAL_JUNCTION                 1.473142 5.503714e-03 3.931224e-02
#> HALLMARK_BILE_ACID_METABOLISM            1.457941 2.406557e-02 1.002732e-01
#> HALLMARK_IL2_STAT5_SIGNALING             1.453587 4.734781e-03 3.931224e-02
#> HALLMARK_HYPOXIA                         1.450031 1.590027e-02 7.227397e-02
```

#### GSEA with GO Biological Process

``` r

# Run GSEA with GO terms
set.seed(20260101)
gsea_go <- get_gsea(
  vista,
  sample_comparison = comp_names[1],
  set_type = "go",
  from_type = "ENSEMBL",
  orgdb = org.Hs.eg.db,
  ont = "BP",  # Biological Process
  pvalueCutoff = 0.05
)

# Show results
if (!is.null(gsea_go$enrich) && nrow(gsea_go$enrich@result) > 0) {
  head(gsea_go$enrich@result[, c("Description", "NES", "pvalue", "p.adjust")], n = 10)
}
#>                                                                           Description
#> GO:0014888                                                 striated muscle adaptation
#> GO:0035357               peroxisome proliferator activated receptor signaling pathway
#> GO:0071276                                           cellular response to cadmium ion
#> GO:0071294                                              cellular response to zinc ion
#> GO:0035358 regulation of peroxisome proliferator activated receptor signaling pathway
#> GO:0046688                                                     response to copper ion
#> GO:0044060                                            regulation of endocrine process
#> GO:0097501                                               stress response to metal ion
#> GO:0071280                                            cellular response to copper ion
#> GO:0060986                                                endocrine hormone secretion
#>                 NES       pvalue    p.adjust
#> GO:0014888 2.173415 2.226967e-06 0.005564077
#> GO:0035357 2.105064 6.608012e-05 0.022798908
#> GO:0071276 2.101350 1.369619e-05 0.009777125
#> GO:0071294 2.094459 8.048658e-06 0.009727138
#> GO:0035358 2.039947 4.082184e-04 0.059684360
#> GO:0046688 2.034456 2.290070e-04 0.049038254
#> GO:0044060 2.030552 1.839123e-04 0.046606174
#> GO:0097501 2.012979 7.586127e-06 0.009727138
#> GO:0071280 2.012772 2.337896e-04 0.049038254
#> GO:0060986 2.002814 5.223680e-04 0.059684360
```

#### GSEA enrichment overview

``` r

# Show all significant pathways using VISTA's visualization
if (!is.null(gsea_results$enrich) && nrow(gsea_results$enrich@result) > 0) {
  get_enrichment_plot(gsea_results$enrich, top_n = 15)
}
```

![](VISTA-airway_files/figure-html/gsea-dotplot-1.png)

#### GSEA plot for top pathway

``` r

# Show enrichment plot for the top pathway
if (!is.null(gsea_results$enrich) && nrow(gsea_results$enrich@result) > 0) {
  # Create GSEA enrichment plot with running enrichment score
  enrichplot::gseaplot2(
    gsea_results$enrich,
    geneSetID = 1,  # Top pathway
    title = gsea_results$enrich@result$Description[1],
    pvalue_table = TRUE
  )
}
```

![](VISTA-airway_files/figure-html/gsea-pathway-1.png)

#### GSEA plot for multiple pathways

``` r

# Show top 3 pathways together
if (!is.null(gsea_results$enrich) && nrow(gsea_results$enrich@result) > 0) {
  enrichplot::gseaplot2(
    gsea_results$enrich,
    geneSetID = 1:3,  # Top 3 pathways
    pvalue_table = TRUE,
    ES_geom = "line"
  )
}
```

![](VISTA-airway_files/figure-html/gsea-multi-1.png)

#### GSEA with GO visualization

``` r

# Visualize GO GSEA results
if (!is.null(gsea_go$enrich) && nrow(gsea_go$enrich@result) > 0) {
  get_enrichment_plot(gsea_go$enrich, top_n = 15)
}
```

![](VISTA-airway_files/figure-html/gsea-go-plot-1.png)

### KEGG Pathway Enrichment

#### KEGG upregulated genes

``` r

kegg_up <- get_kegg_enrichment(
  vista,
  sample_comparison = comp_names[1],
  regulation = "Up",
  species = "Homo sapiens",
  from_type = "ENSEMBL"
)

if (!is.null(kegg_up$enrich) && nrow(kegg_up$enrich@result) > 0) {
  head(kegg_up$enrich@result[, c("Description", "pvalue", "p.adjust", "Count")], n = 10)
}
```

#### KEGG downregulated genes

``` r

kegg_down <- get_kegg_enrichment(
  vista,
  sample_comparison = comp_names[1],
  regulation = "Down",
  species = "Homo sapiens",
  from_type = "ENSEMBL"
)

if (!is.null(kegg_down$enrich) && nrow(kegg_down$enrich@result) > 0) {
  head(kegg_down$enrich@result[, c("Description", "pvalue", "p.adjust", "Count")])
}
```

#### KEGG Visualization

``` r

if (!is.null(kegg_up$enrich) && nrow(kegg_up$enrich@result) > 0) {
  get_enrichment_plot(kegg_up$enrich, top_n = 15)
}
```

## Fold-Change Analysis

### Fold-change Matrix

Useful for comparing multiple comparisons:

``` r

fc_matrix <- get_foldchange_matrix(vista)
head(fc_matrix, n = 10)
#>                 Dexamethasone_VS_Untreated
#> ENSG00000000003                -0.38027500
#> ENSG00000000419                 0.20227695
#> ENSG00000000457                 0.03272066
#> ENSG00000000460                -0.11851609
#> ENSG00000000971                 0.43942357
#> ENSG00000001036                -0.24322719
#> ENSG00000001084                -0.03060707
#> ENSG00000001167                -0.49118392
#> ENSG00000001460                -0.13476909
#> ENSG00000001461                -0.04363732
```

### Fold-change Barplot and Lollipop

#### Per-gene fold-change barplot

``` r

get_foldchange_barplot(
  vista,
  genes = top_up[1:3],
  sample_comparisons = comp_names,
  display_id = "SYMBOL",
  facet_by = "gene"
)
```

![](VISTA-airway_files/figure-html/fc-barplot-gene-1.png)

#### Per-gene fold-change lollipop

``` r

get_foldchange_lollipop(
  vista,
  sample_comparison = comp_names[1],
  genes = top_up[1:3],
  display_id = "SYMBOL",
  facet_by = "gene"
)
```

![](VISTA-airway_files/figure-html/fc-lollipop-gene-1.png)

### Fold-change Raincloud

``` r

get_foldchange_raincloud(
  vista,
  sample_comparisons = comp_names,
  facet_by = "none",
  id.long.var = "gene_id",
  stats_group = TRUE
)
```

![](VISTA-airway_files/figure-html/raincloud-foldchange-1.png)

### Fold-change Heatmap

#### Basic FC heatmap

``` r

get_foldchange_heatmap(vista)
```

![](VISTA-airway_files/figure-html/fc-heatmap-basic-1.png)

#### FC heatmap for selected genes

``` r

# Select genes with large fold-changes
fc_genes <- rownames(fc_matrix)[abs(fc_matrix[, 1]) > 2][1:30]

get_foldchange_heatmap(
  vista,
  sample_comparisons = comp_names,
  genes = fc_genes,
  display_id = "SYMBOL"
)
```

![](VISTA-airway_files/figure-html/fc-heatmap-explicit-1.png)

#### FC heatmap with gene names

``` r

get_foldchange_heatmap(
  vista,
  sample_comparisons = comp_names,
  genes = fc_genes[1:25],
  show_row_names = TRUE,
  display_id = "SYMBOL"
)
```

![](VISTA-airway_files/figure-html/fc-heatmap-names-1.png)

#### FC heatmap for specific gene set

``` r

# Use top upregulated genes
get_foldchange_heatmap(
  vista,
  sample_comparisons = comp_names,
  genes = top_up,
  show_row_names = TRUE,
  display_id = "SYMBOL"
)
```

![](VISTA-airway_files/figure-html/fc-heatmap-custom-1.png)

## Export Results

### Export DE results to file

``` r

# Export complete DE table with annotations
de_annotated <- merge(
  comparisons(vista)[[1]],
  as.data.frame(rowData(vista)),
  by.x = "gene_id",
  by.y = "row.names"
)

# Write to CSV
write.csv(
  de_annotated,
  file = "airway_dexamethasone_vs_untreated_DE_results.csv",
  row.names = FALSE
)

# Export significant genes only
sig_genes <- de_annotated[de_annotated$regulation %in% c("Up", "Down"), ]
write.csv(
  sig_genes,
  file = "airway_significant_DEGs.csv",
  row.names = FALSE
)
```

### Save VISTA object

``` r

# Save the complete VISTA object for later use
saveRDS(vista, file = "airway_vistaect.rds")

# Load it back
# vista <- readRDS("airway_vistaect.rds")
```

## Summary

### Workflow Completed

In this workflow, we:

1.  ✅ Loaded the airway RNA-seq dataset
2.  ✅ Created a VISTA object with DESeq2 analysis
3.  ✅ Added gene annotations from org.Hs.eg.db
4.  ✅ Performed quality control (PCA, MDS, UMAP, correlation)
5.  ✅ Visualized differential expression (volcano, MA plots)
6.  ✅ Analyzed expression patterns (heatmaps, barplots, boxplots,
    violin/raincloud plots, and more)
7.  ✅ Performed functional enrichment (MSigDB, GO, KEGG)
8.  ✅ Explored fold-change patterns
9.  ✅ Generated publication-ready visualizations

### Key Features Demonstrated

- **Single-function workflow**:
  [`create_vista()`](https://cparsania.github.io/VISTA/reference/create_vista.md)
  handles DE analysis
- **Consistent interface**: All plot functions follow the same pattern
- **Flexible visualizations**: Easy to customize colors, labels,
  thresholds
- **Multiple plot types**: 29+ plotting functions for every analysis
  need
- **Integrated enrichment**: No need to wrangle gene IDs manually
- **Publication-ready**: All plots are ggplot2/ComplexHeatmap objects

### Plotting Functions Used

#### QC Plots

- [`get_corr_heatmap()`](https://cparsania.github.io/VISTA/reference/get_corr_heatmap.md) -
  Sample correlation
- [`get_pca_plot()`](https://cparsania.github.io/VISTA/reference/get_pca_plot.md) -
  Principal component analysis
- [`get_mds_plot()`](https://cparsania.github.io/VISTA/reference/get_mds_plot.md) -
  Multidimensional scaling
- [`get_umap_plot()`](https://cparsania.github.io/VISTA/reference/get_umap_plot.md) -
  Nonlinear sample embedding

#### DE Visualization

- [`get_deg_count_barplot()`](https://cparsania.github.io/VISTA/reference/get_deg_count_barplot.md) -
  DEG summary counts
- [`get_volcano_plot()`](https://cparsania.github.io/VISTA/reference/get_volcano_plot.md) -
  Volcano plots
- [`get_ma_plot()`](https://cparsania.github.io/VISTA/reference/get_ma_plot.md) -
  MA plots

#### Expression Plots

- [`get_expression_heatmap()`](https://cparsania.github.io/VISTA/reference/get_expression_heatmap.md) -
  Expression heatmaps
- [`get_expression_barplot()`](https://cparsania.github.io/VISTA/reference/get_expression_barplot.md) -
  Expression barplots
- [`get_expression_boxplot()`](https://cparsania.github.io/VISTA/reference/get_expression_boxplot.md) -
  Expression boxplots
- [`get_expression_violinplot()`](https://cparsania.github.io/VISTA/reference/get_expression_violinplot.md) -
  Violin plots
- [`get_expression_density()`](https://cparsania.github.io/VISTA/reference/get_expression_density.md) -
  Density plots
- [`get_expression_scatter()`](https://cparsania.github.io/VISTA/reference/get_expression_scatter.md) -
  Sample-vs-sample scatter
- [`get_expression_lineplot()`](https://cparsania.github.io/VISTA/reference/get_expression_lineplot.md) -
  Expression across samples
- [`get_expression_lollipop()`](https://cparsania.github.io/VISTA/reference/get_expression_lollipop.md) -
  Lollipop plots
- [`get_expression_joyplot()`](https://cparsania.github.io/VISTA/reference/get_expression_joyplot.md) -
  Ridgeline plots

#### Enrichment Plots

- [`get_enrichment_plot()`](https://cparsania.github.io/VISTA/reference/get_enrichment_plot.md) -
  Generic enrichment visualization
- [`get_msigdb_enrichment()`](https://cparsania.github.io/VISTA/reference/get_msigdb_enrichment.md) -
  MSigDB enrichment
- [`get_go_enrichment()`](https://cparsania.github.io/VISTA/reference/get_go_enrichment.md) -
  GO enrichment
- [`get_kegg_enrichment()`](https://cparsania.github.io/VISTA/reference/get_kegg_enrichment.md) -
  KEGG pathway enrichment
- [`get_pathway_genes()`](https://cparsania.github.io/VISTA/reference/get_pathway_genes.md) -
  Extract genes driving enriched pathways
- [`get_pathway_heatmap()`](https://cparsania.github.io/VISTA/reference/get_pathway_heatmap.md) -
  Plot pathway-derived expression heatmaps
- [`get_enrichment_chord()`](https://cparsania.github.io/VISTA/reference/get_enrichment_chord.md) -
  Chord diagram of gene-pathway relationships

#### Fold-Change

- [`get_foldchange_matrix()`](https://cparsania.github.io/VISTA/reference/get_foldchange_matrix.md) -
  Extract FC matrix
- [`get_foldchange_heatmap()`](https://cparsania.github.io/VISTA/reference/get_foldchange_heatmap.md) -
  Visualize FC patterns

### Next Steps

- Try with your own data
- Explore edgeR backend: `method = "edger"`
- Explore limma-voom backend: `method = "limma"`
- Test multiple comparisons simultaneously
- Customize plots with ggplot2 themes
- Generate automated reports with
  [`run_vista_report()`](https://cparsania.github.io/VISTA/reference/run_vista_report.md)
- Integrate with downstream tools

## Session Information

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] magrittr_2.0.5              org.Hs.eg.db_3.23.1        
#>  [3] AnnotationDbi_1.74.0        airway_1.32.0              
#>  [5] SummarizedExperiment_1.42.0 Biobase_2.72.0             
#>  [7] GenomicRanges_1.64.0        Seqinfo_1.2.0              
#>  [9] IRanges_2.46.0              S4Vectors_0.50.1           
#> [11] BiocGenerics_0.58.1         generics_0.1.4             
#> [13] MatrixGenerics_1.24.0       matrixStats_1.5.0          
#> [15] ggplot2_4.0.3               VISTA_1.1.5                
#> [17] BiocStyle_2.40.0           
#> 
#> loaded via a namespace (and not attached):
#>   [1] splines_4.6.1           ggplotify_0.1.3         tibble_3.3.1           
#>   [4] ggpp_0.6.1              polyclip_1.10-7         enrichit_0.2.1         
#>   [7] lifecycle_1.0.5         httr2_1.3.0             rstatix_1.1.0          
#>  [10] edgeR_4.10.1            doParallel_1.0.17       processx_3.9.0         
#>  [13] lattice_0.22-9          MASS_7.3-65             backports_1.5.1        
#>  [16] limma_3.68.4            sass_0.4.10             rmarkdown_2.31         
#>  [19] jquerylib_0.1.4         yaml_2.3.12             otel_0.2.0             
#>  [22] ggtangle_0.1.2          EnhancedVolcano_1.30.0  DBI_1.3.0              
#>  [25] RColorBrewer_1.1-3      abind_1.4-8             purrr_1.2.2            
#>  [28] msigdbr_26.1.0          yulab.utils_0.2.4       tweenr_2.0.3           
#>  [31] rappdirs_0.3.4          aisdk_1.4.12            gdtools_0.5.1          
#>  [34] circlize_0.4.18         enrichplot_1.32.0       ggrepel_0.9.8          
#>  [37] tidytree_0.4.8          RSpectra_0.16-2         pkgdown_2.2.1          
#>  [40] codetools_0.2-20        DelayedArray_0.38.2     DOSE_4.6.0             
#>  [43] ggforce_0.5.0           tidyselect_1.2.1        shape_1.4.6.1          
#>  [46] aplot_0.3.1             farver_2.1.2            jsonlite_2.0.0         
#>  [49] GetoptLong_1.1.1        Formula_1.2-6           ggridges_0.5.7         
#>  [52] iterators_1.0.14        systemfonts_1.3.2       foreach_1.5.2          
#>  [55] tools_4.6.1             ggnewscale_0.5.2        treeio_1.36.1          
#>  [58] ragg_1.5.2              Rcpp_1.1.2              glue_1.8.1             
#>  [61] gridExtra_2.3.1         SparseArray_1.12.2      xfun_0.60              
#>  [64] DESeq2_1.52.0           qvalue_2.44.0           dplyr_1.2.1            
#>  [67] withr_3.0.3             BiocManager_1.30.27     fastmap_1.2.0          
#>  [70] GGally_2.4.0            ggpointdensity_0.2.1    callr_3.8.0            
#>  [73] digest_0.6.39           R6_2.6.1                gridGraphics_0.5-1     
#>  [76] textshaping_1.0.5       colorspace_2.1-3        GO.db_3.23.1           
#>  [79] RSQLite_3.53.3          ggrain_0.1.2            tidyr_1.3.2            
#>  [82] fontLiberation_0.1.0    FNN_1.1.4.1             httr_1.4.8             
#>  [85] htmlwidgets_1.6.4       S4Arrays_1.12.0         scatterpie_0.2.6       
#>  [88] ggstats_0.13.0          uwot_0.2.4              pkgconfig_2.0.3        
#>  [91] gtable_0.3.6            blob_1.3.0              ComplexHeatmap_2.28.0  
#>  [94] S7_0.2.2                XVector_0.52.0          clusterProfiler_4.20.0 
#>  [97] htmltools_0.5.9         carData_3.0-6           fontBitstreamVera_0.1.1
#> [100] bookdown_0.47           clue_0.3-68             scales_1.4.0           
#> [103] png_0.1-9               ggfun_0.2.1             knitr_1.51             
#> [106] reshape2_1.4.5          rjson_0.2.23            nlme_3.1-169           
#> [109] curl_7.1.0              cachem_1.1.0            GlobalOptions_0.1.4    
#> [112] stringr_1.6.0           parallel_4.6.1          desc_1.4.3             
#> [115] pillar_1.11.1           grid_4.6.1              vctrs_0.7.3            
#> [118] ggpubr_1.0.0            car_3.1-5               tidydr_0.0.6           
#> [121] cluster_2.1.8.2         evaluate_1.0.5          cli_3.6.6              
#> [124] locfit_1.5-9.12         compiler_4.6.1          rlang_1.3.0            
#> [127] crayon_1.5.3            ggsignif_0.6.4          labeling_0.4.3         
#> [130] ps_1.9.3                forcats_1.0.1           plyr_1.8.9             
#> [133] fs_2.1.0                ggiraph_0.9.6           stringi_1.8.9          
#> [136] viridisLite_0.4.3       BiocParallel_1.46.0     assertthat_0.2.1       
#> [139] babelgene_22.9          Biostrings_2.80.1       lazyeval_0.2.3         
#> [142] GOSemSim_2.38.3         fontquiver_0.2.1        Matrix_1.7-5           
#> [145] patchwork_1.3.2         bit64_4.8.2             KEGGREST_1.52.2        
#> [148] statmod_1.5.2           igraph_2.3.3            broom_1.0.13           
#> [151] memoise_2.0.1           bslib_0.12.0            ggtree_4.2.0           
#> [154] bit_4.6.0               ape_5.8-1               gson_0.2.1             
#> [157] polynom_1.4-1
```

## References

- Himes BE et al. (2014). RNA-Seq Transcriptome Profiling Identifies
  CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine
  Function in Airway Smooth Muscle Cells. *PLoS One*, 9(6), e99625.
- Love MI, Huber W, Anders S (2014). Moderated estimation of fold change
  and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15,
  550.
- Robinson MD, McCarthy DJ, Smyth GK (2010). edgeR: a Bioconductor
  package for differential expression analysis of digital gene
  expression data. *Bioinformatics*, 26(1), 139-140.
- Subramanian A et al. (2005). Gene set enrichment analysis: a
  knowledge-based approach for interpreting genome-wide expression
  profiles. *PNAS*, 102(43), 15545-15550.
- Yu G et al. (2012). clusterProfiler: an R package for comparing
  biological themes among gene clusters. *OMICS*, 16(5), 284-287.
