# Accessor Methods for VISTA Object

These accessor functions expose the analysis components stored in a
`VISTA` object. Core expression matrices and annotations live in the
underlying `SummarizedExperiment`, while differential expression
results, summaries, and configuration details are kept inside
`metadata(x)`.

## Usage

``` r
# S4 method for class 'VISTA'
comparisons(object)

# S4 method for class 'VISTA'
deg_summary(object)

# S4 method for class 'VISTA'
cutoffs(object)

# S4 method for class 'VISTA'
norm_counts(object, summarise = FALSE)

# S4 method for class 'VISTA'
sample_info(object)

# S4 method for class 'VISTA'
row_data(object)

# S4 method for class 'VISTA'
group_colors(object)

# S4 method for class 'VISTA'
group_palette(object)
```

## Arguments

- object:

  An object of class `VISTA`.

- summarise:

  Logical. If TRUE, returns mean-normalized counts grouped by the
  grouping column stored in the `VISTA` object (e.g., condition or
  treatment). Default is FALSE.

## Value

The content of the respective slot or processed data:

- comparisons:

  A named list of differential expression tables stored in
  `metadata(x)$de_results`.

- deg_summary:

  A named list of DEG summary tables stored in `metadata(x)$de_summary`.

- cutoffs:

  A list of analysis thresholds held in `metadata(x)$de_cutoffs` (empty
  list if absent).

- norm_counts:

  A matrix of normalized counts, optionally averaged by group.

- sample_info:

  A `DataFrame` of sample metadata.

- row_data:

  A `DataFrame` of gene-level annotation (e.g., baseMean, gene ID).

- group_colors:

  A named character vector of colours from `metadata(x)$group$colors`.

- group_palette:

  The qualitative palette name stored in `metadata(x)$group$palette`.

## See also

[`vista()`](vista.md), [`run_deseq_analysis()`](run_deseq_analysis.md)

## Examples

``` r
# Create example VISTA object
data("count_data", package = "VISTA")
data("sample_metadata", package = "VISTA")

vista <- create_vista(
  counts = count_data[1:100, ],
  sample_info = sample_metadata[1:6, ],
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control"
)
#> Warning: some variables in design formula are characters, converting to factors
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing

# Access differential expression comparisons
comps <- comparisons(vista)
names(comps)
#> [1] "treatment1_VS_control"

# View DEG summary statistics
deg_summary(vista)
#> $treatment1_VS_control
#> # A tibble: 2 × 2
#>   regulation     n
#>   <chr>      <int>
#> 1 Other         83
#> 2 Up             2
#> 

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
#> [1] 1
#> 

# Access normalized counts
nc <- norm_counts(vista)
head(nc)
#>                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000000003  683.40828  520.64761  760.37986  625.21466 1005.78678
#> ENSG00000000419  470.03191  598.51232  540.88876  559.32194  518.80214
#> ENSG00000000457  261.68800  245.21573  229.07205  251.31178  216.53582
#> ENSG00000000460   60.38954   63.91879   34.83986   53.63361   68.93793
#> ENSG00000000971 3272.10648 4275.58608 5380.14477 6515.71752 5940.15198
#> ENSG00000001036 1442.30347 1234.21376 1509.43676 1350.03460 1258.55921
#>                 SRR1039517
#> ENSG00000000003  765.83276
#> ENSG00000000419  584.43207
#> ENSG00000000457  242.11141
#> ENSG00000000460   46.08163
#> ENSG00000000971 8065.74772
#> ENSG00000001036 1052.56289

# Get group-summarized counts
nc_summary <- norm_counts(vista, summarise = TRUE)
head(nc_summary)
#>                    control treatment1
#> ENSG00000000003  816.52497  637.23168
#> ENSG00000000419  509.90761  580.75544
#> ENSG00000000457  235.76529  246.21297
#> ENSG00000000460   54.72244   54.54468
#> ENSG00000000971 4864.13441 6285.68377
#> ENSG00000001036 1403.43315 1212.27042

# Access sample metadata
sample_info(vista)
#> DataFrame with 6 rows and 13 columns
#>            SampleName     cell      dex    albut        Run avgLength
#>              <factor> <factor> <factor> <factor>   <factor> <integer>
#> SRR1039508 GSM1275862  N61311     untrt    untrt SRR1039508       126
#> SRR1039509 GSM1275863  N61311     trt      untrt SRR1039509       126
#> SRR1039512 GSM1275866  N052611    untrt    untrt SRR1039512       126
#> SRR1039513 GSM1275867  N052611    trt      untrt SRR1039513        87
#> SRR1039516 GSM1275870  N080611    untrt    untrt SRR1039516       120
#> SRR1039517 GSM1275871  N080611    trt      untrt SRR1039517       126
#>            Experiment    Sample    BioSample  cond_long  cond_short      groups
#>              <factor>  <factor>     <factor>   <factor> <character> <character>
#> SRR1039508  SRX384345 SRS508568 SAMN02422669 control           CTRL     control
#> SRR1039509  SRX384346 SRS508567 SAMN02422675 treatment1      TREAT1  treatment1
#> SRR1039512  SRX384349 SRS508571 SAMN02422678 control           CTRL     control
#> SRR1039513  SRX384350 SRS508572 SAMN02422670 treatment1      TREAT1  treatment1
#> SRR1039516  SRX384353 SRS508575 SAMN02422682 control           CTRL     control
#> SRR1039517  SRX384354 SRS508576 SAMN02422673 treatment1      TREAT1  treatment1
#>            sizeFactor
#>             <numeric>
#> SRR1039508   0.993550
#> SRR1039509   0.860467
#> SRR1039512   1.148110
#> SRR1039513   0.652576
#> SRR1039516   1.131453
#> SRR1039517   1.367139

# Access gene-level annotations
head(row_data(vista))
#> DataFrame with 6 rows and 1 column
#>                  baseMean
#>                 <numeric>
#> ENSG00000000003  726.8783
#> ENSG00000000419  545.3315
#> ENSG00000000457  240.9891
#> ENSG00000000460   54.6336
#> ENSG00000000971 5574.9091
#> ENSG00000001036 1307.8518

# Get group colors
group_colors(vista)
#>    control treatment1 
#>  "#C87A8A"  "#00A396" 

# Get palette name
group_palette(vista)
#> [1] "Dark 2"
```
