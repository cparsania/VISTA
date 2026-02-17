# Enhanced volcano plot with smart column detection & coloring

Internal helper that wraps EnhancedVolcano and colors points by
regulation.

- If `regulation_col` is present (e.g., "Up"/"Down"/"Other"), it is
  used.

- Otherwise, regulation is derived from FC/P cutoffs.

## Usage

``` r
.EnhancedVolcano2(
  toptable,
  lab,
  x = NULL,
  y = NULL,
  pCutoff = 1e-04,
  FCcutoff = 1.5,
  col_by_regul = TRUE,
  regulation_col = NULL,
  col_up = "#b2182b",
  col_down = "#2166ac",
  col_others = "#e0e0e0",
  return_keyvals = FALSE,
  ...
)
```

## Arguments

- toptable:

  data.frame with DE results.

- lab:

  Character vector of labels (same length as nrow(toptable)).

- x:

  FC column name (auto-detected if NULL).

- y:

  P/P-adj column name (auto-detected if NULL).

- pCutoff:

  Numeric; P-value cutoff. Default 1e-4.

- FCcutoff:

  Numeric; symmetric scalar or c(down, up). Default 1.5.

- col_by_regul:

  Logical; color by regulation. Default TRUE.

- regulation_col:

  Optional column in `toptable` to use directly (expects levels
  Up/Down/Other).

- col_up, :

  col_down, col_others Colors.

- return_keyvals:

  Logical; invisibly return the computed `keyvals`. Default FALSE.

- ...:

  Passed to EnhancedVolcano::EnhancedVolcano()

## Value

ggplot object (and invisibly the `keyvals` if `return_keyvals=TRUE`)
