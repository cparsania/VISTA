# Set or append rowData annotations on a VISTA object

Accepts a data.frame/tibble/DataFrame of gene-level annotations, aligns
it to the VISTA row order, and stores it in `rowData(x)`. Rows are
matched by a key column (default: tries `gene_id` or rownames); Ensembl
version suffixes can be stripped for matching.

## Usage

``` r
set_rowdata(
  x,
  annotations = NULL,
  orgdb = NULL,
  key_col = NULL,
  keytype = NULL,
  columns = c("SYMBOL", "GENENAME", "ENSEMBL", "ENTREZID", "TXCHROM", "TXSTART", "TXEND"),
  drop_version = TRUE,
  overwrite = FALSE
)
```

## Arguments

- x:

  A `VISTA` object.

- annotations:

  Optional data.frame/tibble/DataFrame with one row per gene and a
  column containing the gene IDs to match against `rownames(x)`. If
  omitted, annotations are pulled from `orgdb`.

- orgdb:

  Optional `OrgDb` object; when supplied (and `annotations` is `NULL`),
  annotations are retrieved via
  [`AnnotationDbi::select()`](https://rdrr.io/pkg/AnnotationDbi/man/AnnotationDb-class.html).

- key_col:

  Name of the column in `annotations` that holds the gene IDs. If
  `NULL`, the function will try `gene_id`, `gene`, `ENSEMBL`, `SYMBOL`,
  or use `rownames(annotations)`. Ignored when `annotations` is `NULL`
  and `orgdb` is used.

- keytype:

  Key type for `orgdb` lookups (e.g., `"ENSEMBL"`, `"SYMBOL"`). If
  `NULL`, inferred from `rownames(x)` (`ENSEMBL` if they start with
  "ENS", otherwise `SYMBOL`).

- columns:

  Character vector of OrgDb columns to retrieve when using `orgdb`.
  Default:
  `c("SYMBOL","GENENAME","ENSEMBL","ENTREZID","TXCHROM","TXSTART","TXEND")`.
  The `TXCHROM`/`TXSTART`/`TXEND` fields carry basic genomic coordinates
  when available in the OrgDb.

- drop_version:

  Logical; if `TRUE`, strips Ensembl version suffixes (e.g., `.1`) from
  both the VISTA rownames and the key column/keys before matching.

- overwrite:

  Logical; if `TRUE`, replaces existing `rowData`. If `FALSE`, new
  columns are appended (overwriting by name when names collide).

## Value

The updated `VISTA` object with `rowData` populated/appended.

## Details

OrgDb packages rarely include full genomic coordinates; the default
`TXCHROM/TXSTART/TXEND` columns may therefore be `NA` unless your OrgDb
provides them. For reliable coordinates, fetch them from an
`EnsDb`/`TxDb` (via `genes()` or `biomaRt`/AnnotationHub), build an
annotation table keyed on your gene IDs, and supply that via the
`annotations` argument. When fetching from an `OrgDb`, only columns
available in that database will be filled.

## Examples

``` r
vista <- example_vista()
custom_annot <- data.frame(
  gene_id = rownames(vista)[seq_len(10)],
  custom_info = paste0("Info_", seq_len(10))
)
vista2 <- set_rowdata(vista, annotations = custom_annot, key_col = "gene_id")
#> Warning: Missing annotations for 113 genes; filling those rows with NA.
head(SummarizedExperiment::rowData(vista2)$custom_info)
#> [1] "Info_1" "Info_2" "Info_3" "Info_4" "Info_5" "Info_6"

# \donttest{
# Load example VISTA object
data("count_data", package = "VISTA")
data("sample_metadata", package = "VISTA")

vista <- create_vista(
  counts = count_data[seq_len(100), ],
  sample_info = sample_metadata[seq_len(6), ],
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control"
)
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing

# Add annotations from OrgDb (human)
if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  vista <- set_rowdata(
    vista,
    orgdb = org.Hs.eg.db::org.Hs.eg.db,
    columns = c("SYMBOL", "GENENAME", "ENTREZID")
  )

  # View updated rowData
  head(SummarizedExperiment::rowData(vista))
}
#> 'select()' returned 1:many mapping between keys and columns
#> DataFrame with 6 rows and 4 columns
#>                  baseMean      SYMBOL               GENENAME    ENTREZID
#>                 <numeric> <character>            <character> <character>
#> ENSG00000000003  726.8783      TSPAN6          tetraspanin 6        7105
#> ENSG00000000419  545.3315        DPM1 dolichyl-phosphate m..        8813
#> ENSG00000000457  240.9891       SCYL3 SCY1 like pseudokina..       57147
#> ENSG00000000460   54.6336       FIRRM FIGNL1 interacting r..       55732
#> ENSG00000000971 5574.9091         CFH    complement factor H        3075
#> ENSG00000001036 1307.8518       FUCA2   alpha-L-fucosidase 2        2519

# Or provide custom annotations
custom_annot <- data.frame(
  gene_id = rownames(vista)[seq_len(10)],
  custom_info = paste0("Info_", seq_len(10))
)
vista <- set_rowdata(vista, annotations = custom_annot, key_col = "gene_id")
#> Warning: Missing annotations for 75 genes; filling those rows with NA.
# }
```
