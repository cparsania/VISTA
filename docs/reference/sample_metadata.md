# Sample metadata accompanying the VISTA airway example counts

Metadata describing each column in `count_data`, including friendly
sample names and the experimental group assignment used throughout the
vignettes and tests.

## Usage

``` r
data(sample_metadata)
```

## Format

A data frame with 8 rows and 13 columns:

- sample_names:

  Character identifiers matching the column names in `count_data`.

- groups:

  Group labels used in examples (`"control"`, `"treatment1"`).

- cond_long:

  Long-form condition labels (`"control"`, `"treatment1"`).

- cond_short:

  Short-form condition labels (`"CTRL"`, `"TREAT1"`).

- SampleName:

  Original sample label from the airway colData.

- dex:

  Original airway treatment labels (`"untrt"`/`"trt"`).

- cell:

  Donor/cell-line identifiers from airway.

- albut:

  Original airway albuterol indicator.

- Run:

  SRA run identifier (also copied into `sample_names`).

- avgLength:

  Average transcript length metadata from airway.

- Experiment:

  SRA experiment accession from airway.

- Sample:

  SRA sample accession from airway.

- BioSample:

  NCBI BioSample accession from airway.

## Source

Derived from the `airway` Bioconductor dataset.

## See also

[count_data](count_data.md)
