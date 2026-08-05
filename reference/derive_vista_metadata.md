# Derive starter sample metadata from count sample names

`derive_vista_metadata()` creates a starter `sample_info` table from
count sample names. It is intended for projects where users have count
columns but do not yet have a separate metadata sheet. The derived table
can be edited, passed through
[`read_vista_metadata()`](https://cparsania.github.io/VISTA/reference/read_vista_metadata.md),
and then aligned with
[`match_vista_inputs()`](https://cparsania.github.io/VISTA/reference/match_vista_inputs.md).

## Usage

``` r
derive_vista_metadata(
  counts,
  column_geneid = NULL,
  sample_names = NULL,
  parser = c("auto", "split", "regex", "none"),
  split = "_",
  fields = NULL,
  pattern = NULL,
  sample_column = "sample_names",
  repair_sample_names = c("auto", "none"),
  return_type = c("data.frame", "template"),
  verbose = TRUE
)
```

## Arguments

- counts:

  Count input accepted by
  [`read_vista_counts()`](https://cparsania.github.io/VISTA/reference/read_vista_counts.md),
  or the list returned by
  [`read_vista_counts()`](https://cparsania.github.io/VISTA/reference/read_vista_counts.md).

- column_geneid:

  Optional gene identifier column for raw tabular count inputs. Ignored
  when `counts` is the list output of
  [`read_vista_counts()`](https://cparsania.github.io/VISTA/reference/read_vista_counts.md).

- sample_names:

  Optional explicit sample names to derive metadata from. When supplied,
  these override names extracted from `counts`.

- parser:

  Metadata parsing mode. `"auto"` tries a simple delimiter-based split
  when sample names have a consistent structure. `"split"` uses `split`
  explicitly. `"regex"` uses `pattern`. `"none"` returns only the
  `sample_names` column.

- split:

  Delimiter used when `parser = "split"` or when `"auto"` chooses
  split-based parsing.

- fields:

  Optional field names for parsed metadata columns. When omitted, VISTA
  uses `part_1`, `part_2`, etc.

- pattern:

  Regular expression used when `parser = "regex"`. Capture groups are
  mapped to `fields` in order.

- sample_column:

  Name of the sample identifier column in the returned metadata. Default
  is `"sample_names"`.

- repair_sample_names:

  Strategy passed to
  [`read_vista_counts()`](https://cparsania.github.io/VISTA/reference/read_vista_counts.md)
  when sample names are taken from `counts`. One of `"auto"` or
  `"none"`.

- return_type:

  Return `"data.frame"` (default) or `"template"`. Both return a data
  frame; `"template"` adds empty placeholder columns for `group` and
  `batch`.

- verbose:

  Logical; print an informational derivation summary.

## Value

A data frame containing `sample_names` plus any parsed metadata columns.

## Examples

``` r
data("count_data", package = "VISTA")
data("sample_metadata", package = "VISTA")

counts_in <- count_data[seq_len(8), c("gene_id", sample_metadata$sample_names[seq_len(6)]), drop = FALSE]
meta <- derive_vista_metadata(
  counts_in,
  column_geneid = "gene_id",
  parser = "regex",
  pattern = "SRR(\\d+)",
  fields = "run_id"
)
#> Derived metadata for 6 samples using parser "regex".
head(meta)
#>   sample_names  run_id
#> 1   SRR1039508 1039508
#> 2   SRR1039509 1039509
#> 3   SRR1039512 1039512
#> 4   SRR1039513 1039513
#> 5   SRR1039516 1039516
#> 6   SRR1039517 1039517
```
