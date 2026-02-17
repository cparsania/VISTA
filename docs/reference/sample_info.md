# Sample metadata accompanying the VISTA example counts

Metadata describing each column in `count_data`, including friendly
sample names and the experimental group assignment used throughout the
vignettes and tests.

## Usage

``` r
data(sample_info)
```

## Format

A tibble with 63 rows and 2 columns:

- sample_names:

  Character identifiers matching the column names in `count_data`.

- groups:

  High-level grouping labels (e.g., `"Scr"`, `"ZNF219_sh5KD"`).

## Source

Simulated metadata bundled with the VISTA package.

## See also

[count_data](count_data.md)
