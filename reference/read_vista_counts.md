# Read and standardize count inputs for VISTA

`read_vista_counts()` helps standardize common RNA-seq count inputs into
a count table that can be passed directly to
[`create_vista()`](create_vista.md). It supports plain matrices/data
frames, featureCounts outputs, STAR gene counts, HTSeq-count outputs,
tximport-like lists, and RSEM gene result files.

## Usage

``` r
read_vista_counts(
  x,
  format = c("auto", "matrix", "featurecounts", "star", "htseq", "tximport", "rsem"),
  gene_id_column = NULL,
  sample_columns = NULL,
  sample_names = NULL,
  annotation_columns = NULL,
  count_column = NULL,
  tx2gene = NULL,
  counts_from = c("counts", "abundance", "length"),
  drop_technical = TRUE,
  remove_special_rows = TRUE,
  make_unique_ids = FALSE,
  repair_sample_names = c("auto", "none"),
  return_type = c("list", "data.frame", "matrix"),
  verbose = TRUE
)
```

## Arguments

- x:

  Count input. Supported values depend on `format` and include a matrix,
  data frame, single file path, vector of file paths, or a tximport-like
  list with `counts`, `abundance`, and/or `length`.

- format:

  Input format. One of `"auto"`, `"matrix"`, `"featurecounts"`,
  `"star"`, `"htseq"`, `"tximport"`, or `"rsem"`.

- gene_id_column:

  Optional gene identifier column in tabular inputs. When omitted, VISTA
  uses common names such as `gene_id`/`Geneid`, or falls back to
  rownames for matrices/data frames with unique rownames.

- sample_columns:

  Optional character vector of sample count columns to retain from
  tabular inputs.

- sample_names:

  Optional sample names to use when `x` is a vector of per-sample files.

- annotation_columns:

  Optional feature annotation columns to retain in the returned
  `row_data`.

- count_column:

  Optional count column selector for formats that expose multiple count
  choices. For STAR, use one of `"unstranded"`, `"stranded_first"`, or
  `"stranded_second"`. For RSEM, this defaults to `"expected_count"`.

- tx2gene:

  Optional two-column mapping used to summarize transcript-level
  tximport-like inputs to genes. The first column should contain
  transcript IDs and the second column gene IDs.

- counts_from:

  Which matrix to extract from a tximport-like input: `"counts"`,
  `"abundance"`, or `"length"`.

- drop_technical:

  Logical; when `TRUE`, drop known technical summary rows from
  STAR/HTSeq inputs.

- remove_special_rows:

  Logical; alias for `drop_technical`, retained for clarity in
  file-based imports.

- make_unique_ids:

  Logical; if `TRUE`, duplicate gene IDs are repaired with
  [`make.unique()`](https://rdrr.io/r/base/make.unique.html). Otherwise
  duplicated gene IDs raise an error.

- repair_sample_names:

  Strategy for repairing sample column names. `"auto"` (default) strips
  common file-path and alignment/count suffixes when the repaired names
  are unique, while `"none"` leaves sample columns unchanged. In
  automatic mode VISTA currently:

  - strips directory paths to the basename

  - uses the parent directory for generic quantification files such as
    `quant.sf` or `abundance.tsv`

  - removes common RNA-seq output suffixes such as
    `Aligned.sortedByCoord.out.bam`, `ReadsPerGene.out.tab`,
    `.genes.results`, `.isoforms.results`, `.bam`, and `.fastq.gz`

  - removes common lane/read suffixes such as `_S1_L001_R1_001`,
    `_L001_R2_001`, `_R1`, and `_R2`

  Repaired names are only applied when they remain non-empty and unique.
  Otherwise VISTA keeps the original count column names and records the
  unchanged mapping in `sample_name_map`.

- return_type:

  Return `"list"` (default), standardized `"data.frame"`, or numeric
  `"matrix"`.

- verbose:

  Logical; print an informational import summary.

## Value

If `return_type = "list"`, a list with:

- counts:

  A standardized count table with `gene_id` plus sample columns.

- row_data:

  Feature metadata aligned to the count table.

- column_geneid:

  Always `"gene_id"` for the standardized output.

- sample_names:

  Sample columns in the standardized count table.

- sample_name_map:

  A two-column mapping of original and repaired sample names.

- input_format:

  Resolved import format.

- report:

  Basic import summary.

If `return_type = "data.frame"`, returns the standardized count table.
If `return_type = "matrix"`, returns a numeric matrix with gene IDs as
rownames.

## Details

Internally, VISTA uses a format-specific importer for each supported
input type, then normalizes the result into a common structure with:

- a count table with a `gene_id` column plus sample columns

- optional feature metadata in `row_data`

- sample names inferred from columns or file names

- an auditable `sample_name_map` showing original and repaired names

## Examples

``` r
data("count_data", package = "VISTA")

cnt <- read_vista_counts(
  count_data[seq_len(25), ],
  format = "matrix",
  gene_id_column = "gene_id"
)
#> Imported 25 features and 8 samples from "matrix" input.

head(cnt$counts[, seq_len(4)])
#>           gene_id SRR1039508 SRR1039509 SRR1039512
#> 1 ENSG00000000003        679        448        873
#> 2 ENSG00000000005          0          0          0
#> 3 ENSG00000000419        467        515        621
#> 4 ENSG00000000457        260        211        263
#> 5 ENSG00000000460         60         55         40
#> 6 ENSG00000000938          0          0          2
cnt$sample_names
#> [1] "SRR1039508" "SRR1039509" "SRR1039512" "SRR1039513" "SRR1039516"
#> [6] "SRR1039517" "SRR1039520" "SRR1039521"
```
