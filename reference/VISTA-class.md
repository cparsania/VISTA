# VISTA S4 Class Definition

The `VISTA` class extends `SummarizedExperiment`. All analysis metadata
(DE results, summaries, cutoffs, grouping info, etc.) is stored in the
`metadata()` slot of the object.

## Value

A `VISTA` S4 object.

## Details

Core elements stored in `metadata(v)`:

- `$de_results`: named `SimpleList` of DE result tables.

- `$de_summary`: named `SimpleList` of DEG summary tables.

- `$de_cutoffs`: named list of thresholds (log2FC, p-value type/cutoff,
  method settings).

- `$group`: list with `column`, `palette`, `colors` describing the
  grouping/fill scheme.

- `$provenance`: list with constructor version, timestamp, session info.

- `$vista_schema_version`: metadata schema tag used for compatibility
  checks.

Feature-level annotations live in `rowData(v)`, sample metadata in
`colData(v)`, and normalized counts (and any additional assays) in
`assay(v, "norm_counts")` by default. Use
[`create_vista()`](create_vista.md) as the primary end-user constructor.
Advanced users can convert an existing `SummarizedExperiment` with
[`as_vista()`](as_vista.md).

## See also

[`create_vista`](create_vista.md), [`as_vista`](as_vista.md),
[`validate_vista`](validate_vista.md)

## Examples

``` r
methods::showClass('VISTA')
#> Class "VISTA" [package "VISTA"]
#> 
#> Slots:
#>                                                                               
#> Name:            colData            assays             NAMES   elementMetadata
#> Class:         DataFrame    Assays_OR_NULL character_OR_NULL         DataFrame
#>                         
#> Name:           metadata
#> Class:              list
#> 
#> Extends: 
#> Class "SummarizedExperiment", directly
#> Class "RectangularData", by class "SummarizedExperiment", distance 2
#> Class "Vector", by class "SummarizedExperiment", distance 2
#> Class "Annotated", by class "SummarizedExperiment", distance 3
#> Class "vector_OR_Vector", by class "SummarizedExperiment", distance 3
```
