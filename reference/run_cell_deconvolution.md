# Run Cell Deconvolution on Bulk RNA-seq from VISTA Object

Estimates cell-type proportions in bulk RNA-seq using single-cell
reference or xCell2.

## Usage

``` r
run_cell_deconvolution(
  x,
  method = c("xCell2"),
  single_cell_reference = NULL,
  reference_labels = NULL,
  gene_id_type = c("auto", "symbol", "ensembl", "ensembl_symbol"),
  xcell2_reference = NULL,
  transform = c("log2", "raw"),
  ...
)
```

## Arguments

- x:

  A VISTA object.

- method:

  Deconvolution method. Currently only `"xCell2"` is supported.

- single_cell_reference:

  Reserved for future reference-based methods (ignored).

- reference_labels:

  Reserved for future reference-based methods (ignored).

- gene_id_type:

  Type of gene identifiers: `"auto"`, `"symbol"`, `"ensembl"`, or
  `"ensembl_symbol"`.

- xcell2_reference:

  Optional xCell2 reference object or dataset name (e.g.,
  `"DICE_demo.xCell2Ref"`). Used when xCell2 exposes `xCell2Analysis()`.

- transform:

  Expression transformation: "log2" or "raw".

- ...:

  Additional arguments passed to the specific method.

## Value

VISTA object with cell_fractions added to metadata.
