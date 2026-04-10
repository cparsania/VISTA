# Preparing Counts and Metadata for VISTA

## Overview

VISTA expects two main inputs for
[`create_vista()`](../../reference/create_vista.md):

1.  A count table with one gene identifier column plus one column per
    sample.
2.  A sample metadata table (`sample_info`) containing `sample_names`
    and the grouping variables used in the analysis.

In practice, users often start from count tables whose column names are
file paths, alignment outputs, or sequencing-derived sample identifiers.
Metadata can also be incomplete or missing entirely. VISTA’s input
preparation helpers are designed to standardize these objects before
differential analysis starts:

- [`read_vista_counts()`](../../reference/read_vista_counts.md)
- [`derive_vista_metadata()`](../../reference/derive_vista_metadata.md)
- [`read_vista_metadata()`](../../reference/read_vista_metadata.md)
- [`match_vista_inputs()`](../../reference/match_vista_inputs.md)

## Load Example Data

``` r
library(VISTA)
library(dplyr)
library(tibble)

data("count_data", package = "VISTA")
data("sample_metadata", package = "VISTA")

dim(count_data)
#> [1] 63677     9
head(sample_metadata[, c("sample_names", "cond_long")])
#> # A tibble: 6 × 2
#>   sample_names cond_long 
#>   <chr>        <chr>     
#> 1 SRR1039508   control   
#> 2 SRR1039509   treatment1
#> 3 SRR1039512   control   
#> 4 SRR1039513   treatment1
#> 5 SRR1039516   control   
#> 6 SRR1039517   treatment1
```

## 1. Standardize Counts

For an ordinary in-memory count table,
[`read_vista_counts()`](../../reference/read_vista_counts.md) keeps the
existing API contract and returns a
[`create_vista()`](../../reference/create_vista.md)-ready structure.

``` r
prepared_counts <- read_vista_counts(
  count_data[seq_len(50), ],
  format = "matrix",
  gene_id_column = "gene_id"
)

names(prepared_counts)
#> [1] "counts"          "row_data"        "column_geneid"   "sample_names"   
#> [5] "sample_name_map" "input_format"    "report"
head(prepared_counts$counts[, seq_len(4)])
#>           gene_id SRR1039508 SRR1039509 SRR1039512
#> 1 ENSG00000000003        679        448        873
#> 2 ENSG00000000005          0          0          0
#> 3 ENSG00000000419        467        515        621
#> 4 ENSG00000000457        260        211        263
#> 5 ENSG00000000460         60         55         40
#> 6 ENSG00000000938          0          0          2
head(prepared_counts$row_data)
#>                         gene_id
#> ENSG00000000003 ENSG00000000003
#> ENSG00000000005 ENSG00000000005
#> ENSG00000000419 ENSG00000000419
#> ENSG00000000457 ENSG00000000457
#> ENSG00000000460 ENSG00000000460
#> ENSG00000000938 ENSG00000000938
```

The returned object includes:

- `counts`: standardized count table with `gene_id`
- `row_data`: feature metadata aligned to the counts
- `sample_names`: standardized sample columns
- `sample_name_map`: mapping from original to repaired sample names

## 2. Repair File-Derived Sample Names

RNA-seq count tables often inherit sample columns from alignment or
quantification files. VISTA can repair these names conservatively when
they are clearly file-derived and the repaired names remain unique.

``` r
counts_paths <- count_data[seq_len(20), c("gene_id", sample_metadata$sample_names[seq_len(4)]), drop = FALSE]

colnames(counts_paths)[2:5] <- c(
  "/proj/run/03_alignment/HET_1_U/HET_1_U_star_alignAligned.sortedByCoord.out.bam",
  "/proj/run/03_alignment/HET_1_ovary/HET_1_ovary_star_alignAligned.sortedByCoord.out.bam",
  "/proj/run/03_alignment/WT_1_U/WT_1_U_star_alignAligned.sortedByCoord.out.bam",
  "/proj/run/03_alignment/WT_1_ovary/WT_1_ovary_star_alignAligned.sortedByCoord.out.bam"
)

prepared_paths <- read_vista_counts(
  counts_paths,
  format = "matrix",
  gene_id_column = "gene_id",
  repair_sample_names = "auto"
)

prepared_paths$sample_name_map
#>                                                                                 original
#> 1         /proj/run/03_alignment/HET_1_U/HET_1_U_star_alignAligned.sortedByCoord.out.bam
#> 2 /proj/run/03_alignment/HET_1_ovary/HET_1_ovary_star_alignAligned.sortedByCoord.out.bam
#> 3           /proj/run/03_alignment/WT_1_U/WT_1_U_star_alignAligned.sortedByCoord.out.bam
#> 4   /proj/run/03_alignment/WT_1_ovary/WT_1_ovary_star_alignAligned.sortedByCoord.out.bam
#>      repaired
#> 1     HET_1_U
#> 2 HET_1_ovary
#> 3      WT_1_U
#> 4  WT_1_ovary
```

Automatic repair currently handles patterns such as:

- full file paths
- STAR outputs (`Aligned.sortedByCoord.out.bam`, `ReadsPerGene.out.tab`)
- RSEM outputs (`.genes.results`, `.isoforms.results`)
- generic quantification files such as `quant.sf`
- lane/read suffixes like `_S1_L001_R1_001`

## 3. Derive Starter Metadata from Sample Names

If you do not yet have a metadata sheet,
[`derive_vista_metadata()`](../../reference/derive_vista_metadata.md)
can create a starter `sample_info` table from the count sample names.

### Split-based parsing

``` r
derived_split <- derive_vista_metadata(
  prepared_paths,
  parser = "split",
  split = "_",
  fields = c("genotype", "replicate", "tissue")
)

derived_split
#>   sample_names genotype replicate tissue
#> 1      HET_1_U      HET         1      U
#> 2  HET_1_ovary      HET         1  ovary
#> 3       WT_1_U       WT         1      U
#> 4   WT_1_ovary       WT         1  ovary
```

### Regex-based parsing

For sequencing accessions or other structured labels, regex parsing
gives more control:

``` r
derived_regex <- derive_vista_metadata(
  counts = NULL,
  sample_names = c("SRR1039508", "SRR1039509", "SRR1039512"),
  parser = "regex",
  pattern = "SRR(\\d+)",
  fields = "run_id"
)

derived_regex
#>   sample_names  run_id
#> 1   SRR1039508 1039508
#> 2   SRR1039509 1039509
#> 3   SRR1039512 1039512
```

### Template mode

Template mode adds empty placeholders for columns users commonly fill in
next, such as `group` and `batch`.

``` r
derive_vista_metadata(
  counts = NULL,
  sample_names = c("sampleA", "sampleB"),
  parser = "none",
  return_type = "template"
)
#>   sample_names group batch
#> 1      sampleA  <NA>  <NA>
#> 2      sampleB  <NA>  <NA>
```

## 4. Standardize Existing Metadata

If a metadata table already exists,
[`read_vista_metadata()`](../../reference/read_vista_metadata.md)
standardizes it into the form expected by
[`create_vista()`](../../reference/create_vista.md).

``` r
prepared_samples <- read_vista_metadata(sample_metadata)

head(prepared_samples[, c("sample_names", "cond_long")])
#>            sample_names  cond_long
#> SRR1039508   SRR1039508    control
#> SRR1039509   SRR1039509 treatment1
#> SRR1039512   SRR1039512    control
#> SRR1039513   SRR1039513 treatment1
#> SRR1039516   SRR1039516    control
#> SRR1039517   SRR1039517 treatment1
```

VISTA will infer `sample_names` from:

- an existing `sample_names` column
- non-default rownames
- common aliases such as `sample`, `sample_id`, or `Run`

## 5. Match Counts and Metadata

[`match_vista_inputs()`](../../reference/match_vista_inputs.md) aligns
count columns and metadata rows so they can be passed directly into
[`create_vista()`](../../reference/create_vista.md).

``` r
matched_inputs <- match_vista_inputs(prepared_counts, prepared_samples)

matched_inputs$column_geneid
#> [1] "gene_id"
identical(matched_inputs$sample_info$sample_names, colnames(matched_inputs$counts)[-1])
#> [1] TRUE
```

When appropriate, you can drop unmatched samples instead of erroring:

``` r
sample_extra <- prepared_samples[seq_len(4), , drop = FALSE]
extra_row <- sample_extra[1, , drop = FALSE]
extra_row$sample_names <- "extra_sample"
if ("cond_long" %in% colnames(extra_row)) extra_row$cond_long <- "control"
if ("dex" %in% colnames(extra_row)) extra_row$dex <- "untrt"
if ("cell" %in% colnames(extra_row)) extra_row$cell <- "N999"
rownames(extra_row) <- "extra_sample"
sample_extra <- rbind(sample_extra, extra_row)

matched_drop <- match_vista_inputs(
  prepared_counts,
  sample_extra,
  drop_unmatched = TRUE
)

matched_drop$report
#> $n_genes
#> [1] 50
#> 
#> $n_samples
#> [1] 4
#> 
#> $dropped_from_counts
#> [1] "SRR1039516" "SRR1039517" "SRR1039520" "SRR1039521"
#> 
#> $dropped_from_sample_info
#> [1] "extra_sample"
```

## 6. Create a VISTA Object

Once inputs are standardized,
[`create_vista()`](../../reference/create_vista.md) is unchanged.

``` r
vista <- create_vista(
  counts = matched_inputs$counts,
  sample_info = matched_inputs$sample_info,
  column_geneid = matched_inputs$column_geneid,
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control",
  method = "deseq2",
  min_counts = 5,
  min_replicates = 1
)

vista
#> class: VISTA 
#> dim: 44 8 
#> metadata(12): de_results de_summary ... design comparison
#> assays(1): norm_counts
#> rownames(44): ENSG00000000003 ENSG00000000419 ... ENSG00000004399
#>   ENSG00000004455
#> rowData names(1): baseMean
#> colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
#> colData names(14): SampleName cell ... sizeFactor sample_names
```

## Practical Guidance

- Use [`read_vista_counts()`](../../reference/read_vista_counts.md)
  first whenever count columns come from external tools or file names.
- Use
  [`derive_vista_metadata()`](../../reference/derive_vista_metadata.md)
  when you only have sample names and need a reviewable starter metadata
  sheet.
- Use [`read_vista_metadata()`](../../reference/read_vista_metadata.md)
  when you already have a metadata table and want VISTA to standardize
  it.
- Use [`match_vista_inputs()`](../../reference/match_vista_inputs.md)
  before [`create_vista()`](../../reference/create_vista.md) to catch
  alignment issues early.

## Session Info

``` r
sessionInfo()
#> R version 4.5.3 (2026-03-11)
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
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] tibble_3.3.1     dplyr_1.2.1      VISTA_0.99.6     BiocStyle_2.38.0
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3          jsonlite_2.0.0             
#>   [3] tidydr_0.0.6                magrittr_2.0.5             
#>   [5] ggtangle_0.1.1              farver_2.1.2               
#>   [7] rmarkdown_2.31              fs_2.0.1                   
#>   [9] ragg_1.5.2                  vctrs_0.7.2                
#>  [11] memoise_2.0.1               ggtree_4.0.5               
#>  [13] rstatix_0.7.3               htmltools_0.5.9            
#>  [15] S4Arrays_1.10.1             curl_7.0.0                 
#>  [17] broom_1.0.12                Formula_1.2-5              
#>  [19] SparseArray_1.10.10         gridGraphics_0.5-1         
#>  [21] sass_0.4.10                 bslib_0.10.0               
#>  [23] htmlwidgets_1.6.4           desc_1.4.3                 
#>  [25] plyr_1.8.9                  cachem_1.1.0               
#>  [27] igraph_2.2.3                lifecycle_1.0.5            
#>  [29] pkgconfig_2.0.3             gson_0.1.0                 
#>  [31] Matrix_1.7-4                R6_2.6.1                   
#>  [33] fastmap_1.2.0               MatrixGenerics_1.22.0      
#>  [35] digest_0.6.39               aplot_0.2.9                
#>  [37] enrichplot_1.30.5           colorspace_2.1-2           
#>  [39] ggnewscale_0.5.2            GGally_2.4.0               
#>  [41] patchwork_1.3.2             AnnotationDbi_1.72.0       
#>  [43] S4Vectors_0.48.1            DESeq2_1.50.2              
#>  [45] textshaping_1.0.5           GenomicRanges_1.62.1       
#>  [47] RSQLite_2.4.6               ggpubr_0.6.3               
#>  [49] polyclip_1.10-7             httr_1.4.8                 
#>  [51] abind_1.4-8                 compiler_4.5.3             
#>  [53] withr_3.0.2                 bit64_4.6.0-1              
#>  [55] fontquiver_0.2.1            backports_1.5.1            
#>  [57] S7_0.2.1                    BiocParallel_1.44.0        
#>  [59] carData_3.0-6               DBI_1.3.0                  
#>  [61] ggstats_0.13.0              ggforce_0.5.0              
#>  [63] R.utils_2.13.0              ggsignif_0.6.4             
#>  [65] MASS_7.3-65                 rappdirs_0.3.4             
#>  [67] DelayedArray_0.36.1         tools_4.5.3                
#>  [69] otel_0.2.0                  scatterpie_0.2.6           
#>  [71] ape_5.8-1                   msigdbr_26.1.0             
#>  [73] R.oo_1.27.1                 glue_1.8.0                 
#>  [75] nlme_3.1-168                GOSemSim_2.36.0            
#>  [77] grid_4.5.3                  cluster_2.1.8.2            
#>  [79] reshape2_1.4.5              fgsea_1.36.2               
#>  [81] generics_0.1.4              gtable_0.3.6               
#>  [83] R.methodsS3_1.8.2           tidyr_1.3.2                
#>  [85] data.table_1.18.2.1         utf8_1.2.6                 
#>  [87] car_3.1-5                   XVector_0.50.0             
#>  [89] BiocGenerics_0.56.0         ggrepel_0.9.8              
#>  [91] pillar_1.11.1               stringr_1.6.0              
#>  [93] babelgene_22.9              limma_3.66.0               
#>  [95] yulab.utils_0.2.4           splines_4.5.3              
#>  [97] tweenr_2.0.3                treeio_1.34.0              
#>  [99] lattice_0.22-9              bit_4.6.0                  
#> [101] tidyselect_1.2.1            fontLiberation_0.1.0       
#> [103] GO.db_3.22.0                locfit_1.5-9.12            
#> [105] Biostrings_2.78.0           knitr_1.51                 
#> [107] fontBitstreamVera_0.1.1     bookdown_0.46              
#> [109] IRanges_2.44.0              Seqinfo_1.0.0              
#> [111] edgeR_4.8.2                 SummarizedExperiment_1.40.0
#> [113] stats4_4.5.3                xfun_0.57                  
#> [115] Biobase_2.70.0              statmod_1.5.1              
#> [117] matrixStats_1.5.0           stringi_1.8.7              
#> [119] lazyeval_0.2.3              ggfun_0.2.0                
#> [121] yaml_2.3.12                 evaluate_1.0.5             
#> [123] codetools_0.2-20            gdtools_0.5.0              
#> [125] qvalue_2.42.0               BiocManager_1.30.27        
#> [127] ggplotify_0.1.3             cli_3.6.5                  
#> [129] systemfonts_1.3.2           jquerylib_0.1.4            
#> [131] Rcpp_1.1.1                  png_0.1-9                  
#> [133] parallel_4.5.3              assertthat_0.2.1           
#> [135] pkgdown_2.2.0               ggplot2_4.0.2              
#> [137] blob_1.3.0                  clusterProfiler_4.18.4     
#> [139] DOSE_4.4.0                  tidytree_0.4.7             
#> [141] ggiraph_0.9.6               scales_1.4.0               
#> [143] purrr_1.2.1                 crayon_1.5.3               
#> [145] rlang_1.2.0                 cowplot_1.2.0              
#> [147] fastmatch_1.1-8             KEGGREST_1.50.0
```
