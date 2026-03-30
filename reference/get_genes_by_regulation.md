# Get Genes by Regulation

Extract gene IDs by regulation class from selected comparisons in a
VISTA object.

## Usage

``` r
get_genes_by_regulation(
  x,
  sample_comparisons,
  regulation = "Both",
  top_n = NULL,
  display_id = NULL,
  return_type = c("list", "table")
)
```

## Arguments

- x:

  A VISTA object.

- sample_comparisons:

  Character vector of comparison names to include.

- regulation:

  One of "Up", "Down", "Both", or "All" (default: "Both").

- top_n:

  Optional integer limiting each comparison to the top DE genes ranked
  by absolute log2 fold change after regulation filtering.

- display_id:

  Optional column name in `rowData(x)` to append when
  `return_type = "table"`, for example `"SYMBOL"`.

- return_type:

  Either `"list"` (default) to return gene ID vectors or `"table"` to
  return one data frame per comparison with gene metadata.

## Value

A named list of character vectors (one per comparison) when
`return_type = "list"`, or a named list of data frames when
`return_type = "table"`.

## Examples

``` r
v <- example_vista()
comp <- names(comparisons(v))[1]
genes <- get_genes_by_regulation(v, sample_comparisons = comp, regulation = 'Up', top_n = 20)
str(genes, max.level = 1)
#> List of 1
#>  $ treatment1_VS_control: chr [1:2] "ENSG00000004799" "ENSG00000003402"
```
