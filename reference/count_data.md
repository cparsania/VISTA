# Example RNA-seq count matrix shipped with VISTA

A count matrix derived from the Bioconductor `airway` dataset and
formatted for VISTA examples. The first column stores Ensembl gene IDs;
the remaining columns are sample-level counts.

## Usage

``` r
data(count_data)
```

## Format

A tibble with 63,677 rows and 9 columns.

- gene_id:

  Character column storing Ensembl gene identifiers.

- SRR1039508:

  Numeric count values for airway sample SRR1039508.

- SRR1039509:

  Numeric count values for airway sample SRR1039509.

- SRR1039512:

  Numeric count values for airway sample SRR1039512.

- SRR1039513:

  Numeric count values for airway sample SRR1039513.

- SRR1039516:

  Numeric count values for airway sample SRR1039516.

- SRR1039517:

  Numeric count values for airway sample SRR1039517.

- SRR1039520:

  Numeric count values for airway sample SRR1039520.

- SRR1039521:

  Numeric count values for airway sample SRR1039521.

## Source

Derived from the `airway` Bioconductor dataset.
