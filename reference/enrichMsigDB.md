# Perform MSigDB over-representation analysis on a VISTA object

This function performs ORA on a set of genes from a VISTA object. Gene
identifiers stored as “ENSG00000187634:SAMD11” will be split at the
colon and the appropriate part extracted based on `from_type`. Genes are
then optionally converted via `orgdb` to match the identifier type
required by the MSigDB gene sets.

## Usage

``` r
enrichMsigDB(
  x,
  gene_list,
  from_type = "SYMBOL",
  orgdb,
  msigdb_category = "H",
  msigdb_subcategory = NULL,
  species = "Mus musculus",
  background = NULL,
  col_genetype = "GENETYPE",
  feature_type = "protein-coding",
  ...
)
```

## Arguments

- x:

  A VISTA object.

- gene_list:

  A non-empty character vector of gene identifiers. This argument is
  required and must not be `NULL`.

- from_type:

  Identifier type of the genes in `gene_list`. One of `"SYMBOL"`,
  `"ENSEMBL"`, or `"ENTREZID"` (default `"SYMBOL"`). Ensembl IDs have
  version suffixes stripped automatically.

- orgdb:

  An `OrgDb` object used for ID conversion. If omitted, a default is
  chosen based on `species`: `org.Mm.eg.db` for mouse, `org.Hs.eg.db`
  for human.

- msigdb_category:

  MSigDB category (e.g. `"H"`, `"C2"`, `"C5"`). Default `"H"`.

- msigdb_subcategory:

  Optional MSigDB sub-collection (e.g. `"BP"` for C5). Default `NULL`.

- species:

  Species name for MSigDB (default `"Mus musculus"`).

- background:

  Optional character vector of background gene IDs. If `NULL`, all
  features in the VISTA object (optionally filtered by `feature_type`)
  are used.

- col_genetype:

  Column name in `rowData(x)` used to filter the background by gene type
  (default `"GENETYPE"`).

- feature_type:

  Gene type to retain in the background when `background` is `NULL`
  (default `"protein-coding"`).

- ...:

  Additional arguments to pass to
  [`clusterProfiler::enricher()`](https://rdrr.io/pkg/clusterProfiler/man/enricher.html).

## Value

A list containing a single `enrichResult` object named `"enrich"`.
