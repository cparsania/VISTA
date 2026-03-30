# DEG Venn diagram

Visualizes overlaps between DEG sets for two to four comparisons.

## Usage

``` r
get_deg_venn_diagram(
  x,
  sample_comparisons,
  regulation = "Up",
  palette = "Set 2",
  auto_scale = FALSE,
  show_percentage = TRUE,
  ...
)
```

## Arguments

- x:

  A `VISTA` object.

- sample_comparisons:

  Character vector of 2–4 comparison names to include in the Venn
  diagram.

- regulation:

  One of `"Up"`, `"Down"`, `"Both"`, or `"All"` selecting which genes to
  include.

- palette:

  Qualitative palette name passed to
  [`colorspace::qualitative_hcl()`](https://colorspace.R-Forge.R-project.org/reference/hcl_palettes.html)
  for fill colors.

- auto_scale:

  Logical; pass through to
  [`ggvenn::ggvenn()`](https://yanlinlin82.github.io/ggvenn/reference/ggvenn.html)
  to scale circles by size.

- show_percentage:

  Logical; request percentage labels from
  [`ggvenn::ggvenn()`](https://yanlinlin82.github.io/ggvenn/reference/ggvenn.html).

- ...:

  Additional arguments forwarded to
  [`ggvenn::ggvenn()`](https://yanlinlin82.github.io/ggvenn/reference/ggvenn.html).

## Value

An object returned by this function.

## Examples

``` r
v <- example_vista()
comps <- names(comparisons(v))
if (length(comps) >= 2) {
  p <- get_deg_venn_diagram(v, sample_comparisons = comps[1:2])
  print(p)
}
```
