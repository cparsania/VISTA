# Shared plotting defaults ---------------------------------------------------
#
# One place where every cross-function default lives. VISTA 1.2.0 does not
# change any value here relative to 1.1.x -- this file exists so that the
# planned 1.4.0 default harmonization is a single atomic edit rather than a
# scatter of magic numbers across ~60 plotting functions.
#
# Values are read through `.vista_default()` rather than inlined as formals, so
# `options(vista.legacy_defaults = TRUE)` can restore 1.x values for one
# deprecation cycle when the flip happens.

#' @keywords internal
#' @noRd
.vista_defaults <- function() {
  list(
    # Point / line geometry
    point_size_embedding = 10,
    point_size_scatter = 1.5,
    point_size_lollipop = 6,
    point_size_volcano = 1,
    point_size_ma = 1.2,
    linewidth_lollipop = 1.2,

    # Text
    label_size = 3,
    label_size_corr = 4,
    row_label_fontsize = 10,
    base_size = 12,

    # Selection
    top_n_heatmap = 50,
    top_n_foldchange_heatmap = 10,
    top_n_enrichment = 10,
    top_n_chord = 8,
    max_genes_embedding = 20,
    max_genes_lollipop = 15,
    max_genes_barplot = 25,
    max_genes_boxplot = 20,

    # Transparency
    alpha_embedding = 0.85,
    alpha_density = 0.4,
    alpha_distribution = 0.5,

    # Palettes
    group_palette = "Dark 2",
    comparison_palette = "Dark 3",
    anno_palette_expression = "Dark 3",
    anno_palette_foldchange = "Set2",

    # Regulation colours
    up_color = "#a40000",
    down_color = "#007e2f",
    other_color = "grey70"
  )
}

#' Look up a shared VISTA default
#'
#' @param name Character scalar naming an entry in [.vista_defaults()].
#' @return The default value.
#' @keywords internal
#' @noRd
.vista_default <- function(name) {
  defaults <- .vista_defaults()
  if (!name %in% names(defaults)) {
    cli::cli_abort("Unknown VISTA default {.val {name}}.")
  }
  defaults[[name]]
}
