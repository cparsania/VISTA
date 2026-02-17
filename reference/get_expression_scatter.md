# Compare normalized expression between two samples or groups

Plots gene-level expression for two selected samples or group means,
colours points by local density (viridis), labels the most divergent
genes, and reports Pearson/Spearman correlation.

## Usage

``` r
get_expression_scatter(
  x,
  sample_x,
  sample_y,
  by = c("sample", "group"),
  group_column = NULL,
  genes = NULL,
  log_transform = TRUE,
  label_top_n = 20,
  method = c("pearson", "spearman"),
  display_id = NULL
)
```

## Arguments

- x:

  A `VISTA` object.

- sample_x:

  First sample or group to plot (character scalar).

- sample_y:

  Second sample or group to plot (character scalar).

- by:

  One of `"sample"` (use individual samples) or `"group"` (average
  replicates within `group_column` before plotting). Default `"sample"`.

- group_column:

  Column in `sample_info` used when `by = "group"` (defaults to stored
  grouping column).

- genes:

  Optional character vector of genes to include; defaults to all.

- log_transform:

  Logical; apply log2(x + 1) transform. Default `TRUE`.

- label_top_n:

  Integer; number of most divergent genes to label (ranked by \|x -
  y\|). Set to 0 to disable labels.

- method:

  Correlation method for the subtitle; `"pearson"` (default) or
  `"spearman"`.

- display_id:

  Optional column in `rowData(x)` to use for point labels (fallback to
  gene_id/rownames when not available).

## Value

A `ggplot2` object.
