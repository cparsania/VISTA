#' @title Accessor Methods for VISTA Object
#' @description These accessor functions expose the analysis components stored in a `VISTA` object.
#' Core expression matrices and annotations live in the underlying
#' `SummarizedExperiment`, while differential expression results, summaries, and
#' configuration details are kept inside `metadata(x)`.
#'
#' @param object An object of class `VISTA`.
#' @param summarise Logical. If TRUE, returns mean-normalized counts grouped by the grouping column
#' stored in the `VISTA` object (e.g., condition or treatment). Default is FALSE.
#' @param source Which DE result source to use for `comparisons()`/`deg_summary()`.
#'   One of `"active"`, `"deseq2"`, `"edger"`, `"limma"`, or `"consensus"`.
#'   `"active"` uses the currently selected source stored in metadata.
#' @aliases comparisons deg_summary cutoffs norm_counts sample_info row_data group_colors group_palette
#'
#' @return The content of the respective slot or processed data:
#' \describe{
#'   \item{comparisons}{A named list of differential expression tables stored in `metadata(x)$de_results`.}
#'   \item{deg_summary}{A named list of DEG summary tables stored in `metadata(x)$de_summary`.}
#'   \item{cutoffs}{A list of analysis thresholds held in `metadata(x)$de_cutoffs` (empty list if absent).}
#'   \item{norm_counts}{A matrix of normalized counts, optionally averaged by group.}
#'   \item{sample_info}{A `DataFrame` of sample metadata.}
#'   \item{row_data}{A `DataFrame` of gene-level annotation (e.g., baseMean, gene ID).}
#'   \item{group_colors}{A named character vector of colours from `metadata(x)$group$colors`.}
#'   \item{group_palette}{The qualitative palette name stored in `metadata(x)$group$palette`.}
#' }
#'
#' @examples
#' # Create example VISTA object
#' data("count_data", package = "VISTA")
#' data("sample_metadata", package = "VISTA")
#'
#' vista <- create_vista(
#'   counts = count_data[1:100, ],
#'   sample_info = sample_metadata[1:6, ],
#'   column_geneid = "gene_id",
#'   group_column = "cond_long",
#'   group_numerator = "treatment1",
#'   group_denominator = "control"
#' )
#'
#' # Access differential expression comparisons
#' comps <- comparisons(vista)
#' names(comps)
#' comparisons(vista, source = "active")
#'
#' # View DEG summary statistics
#' deg_summary(vista)
#'
#' # Get analysis cutoffs
#' cutoffs(vista)
#'
#' # Access normalized counts
#' nc <- norm_counts(vista)
#' head(nc)
#'
#' # Get group-summarized counts
#' nc_summary <- norm_counts(vista, summarise = TRUE)
#' head(nc_summary)
#'
#' # Access sample metadata
#' sample_info(vista)
#'
#' # Access gene-level annotations
#' head(row_data(vista))
#'
#' # Get group colors
#' group_colors(vista)
#'
#' # Get palette name
#' group_palette(vista)
#'
#' # For method = "both", inspect method-specific outputs
#' \donttest{
#' vista_both <- create_vista(
#'   counts = count_data,
#'   sample_info = sample_metadata,
#'   column_geneid = "gene_id",
#'   group_column = "cond_long",
#'   group_numerator = "treatment1",
#'   group_denominator = "control",
#'   method = "both",
#'   result_source = "consensus"
#' )
#' comparisons(vista_both, source = "deseq2")
#' comparisons(vista_both, source = "edger")
#' vista_both <- set_de_source(vista_both, "edger")
#' }
#'
#' @name VISTA-accessors
#' @seealso [create_vista()], [as_vista()], [run_deseq_analysis()]
#' @importFrom SummarizedExperiment assay colData rowData
#' @import methods
NULL

#' @keywords internal
.resolve_de_source <- function(object, source = c("active", "deseq2", "edger", "limma", "consensus"), which = c("results", "summary")) {
  source <- match.arg(source)
  which <- match.arg(which)
  md <- S4Vectors::metadata(object)
  by_method_key <- if (which == "results") "de_results_by_method" else "de_summary_by_method"
  legacy_key <- if (which == "results") "de_results" else "de_summary"

  if (source == "active") {
    active <- md$de_active_source %||% NULL
    if (!is.null(active) && !is.null(md[[by_method_key]]) && active %in% names(md[[by_method_key]])) {
      return(as.list(md[[by_method_key]][[active]]))
    }
    legacy <- md[[legacy_key]]
    if (is.null(legacy)) return(list())
    return(as.list(legacy))
  }

  by_method <- md[[by_method_key]]
  if (is.null(by_method) || !source %in% names(by_method)) {
    available <- if (is.null(by_method)) character(0) else names(by_method)
    cli::cli_abort("DE source {.val {source}} is not available in this object. Available: {.val {available}}")
  }
  as.list(by_method[[source]])
}

#' @export
setGeneric("comparisons", function(object, source = "active") standardGeneric("comparisons"))
#' @rdname VISTA-accessors
#' @export
setMethod("comparisons", "VISTA", function(object, source = "active") {
  .resolve_de_source(object, source = source, which = "results")
})


#' @export
setGeneric("deg_summary", function(object, source = "active") standardGeneric("deg_summary"))
#' @rdname VISTA-accessors
#' @export
setMethod("deg_summary", "VISTA", function(object, source = "active") {
  .resolve_de_source(object, source = source, which = "summary")
})


#' @export
setGeneric("cutoffs", function(object) standardGeneric("cutoffs"))
#' @rdname VISTA-accessors
#' @export
setMethod("cutoffs", "VISTA", function(object) {
  md <- S4Vectors::metadata(object)
  rlang::`%||%`(md$de_cutoffs, list())
})


#' @export
setGeneric("norm_counts", function(object, summarise = FALSE) standardGeneric("norm_counts"))
#' @rdname VISTA-accessors
#' @export
setMethod("norm_counts", "VISTA", function(object, summarise = FALSE) {
  mat <- assay(object, "norm_counts")
  if (!summarise) {
    return(mat)
  }

  si <- as.data.frame(colData(object))
  md <- S4Vectors::metadata(object)
  group_info <- md$group
  group_column <- if (!is.null(group_info)) group_info$column else NULL

  if (is.null(group_column) || !group_column %in% colnames(si)) {
    cli::cli_abort("Grouping column {.val {group_column}} not found in sample_info for summarization.")
  }

  si$sample <- rownames(si)
  group_map <- split(si$sample, si[[group_column]])
  summarized <- vapply(group_map, function(samples) {
    rowMeans(mat[, samples, drop = FALSE])
  }, FUN.VALUE = numeric(nrow(mat)))

  rownames(summarized) <- rownames(mat)
  return(summarized)
})


#' @export
setGeneric("sample_info", function(object) standardGeneric("sample_info"))
#' @rdname VISTA-accessors
#' @export
setMethod("sample_info", "VISTA", function(object) colData(object))


#' @export
setGeneric("row_data", function(object) standardGeneric("row_data"))
#' @rdname VISTA-accessors
#' @export
setMethod("row_data", "VISTA", function(object) rowData(object))


#' @export
setGeneric("group_colors", function(object) standardGeneric("group_colors"))
#' @rdname VISTA-accessors
#' @export
setMethod("group_colors", "VISTA", function(object) {
  md <- S4Vectors::metadata(object)
  group_info <- md$group
  colors <- if (!is.null(group_info)) group_info$colors else NULL
  rlang::`%||%`(colors, character())
})




#' @export
setGeneric("group_palette", function(object) standardGeneric("group_palette"))
#' @rdname VISTA-accessors
#' @export
setMethod("group_palette", "VISTA", function(object) {
  md <- S4Vectors::metadata(object)
  group_info <- md$group
  palette <- if (!is.null(group_info)) group_info$palette else NULL
  rlang::`%||%`(palette, character())
})

#' @keywords internal
.validate_named_color_map <- function(color_map, arg_name = "color_map") {
  if (!is.character(color_map) || !length(color_map)) {
    cli::cli_abort("{.arg {arg_name}} must be a non-empty named character vector of colors.")
  }
  nms <- names(color_map)
  if (is.null(nms) || any(!nzchar(nms)) || anyDuplicated(nms)) {
    cli::cli_abort("{.arg {arg_name}} must have unique, non-empty names.")
  }
  bad <- nms[!vapply(color_map, function(clr) {
    !inherits(try(grDevices::col2rgb(clr), silent = TRUE), "try-error")
  }, logical(1))]
  if (length(bad)) {
    cli::cli_abort("Invalid color value(s) in {.arg {arg_name}} for name(s): {.val {bad}}")
  }
  color_map
}

#' Set manual group colors in a VISTA object
#'
#' Updates `metadata(x)$group$colors` using a user-supplied named color vector.
#' This controls group-level coloring across VISTA plots that use group mapping
#' (for example PCA, MDS, and expression plots). The color map must include all
#' groups currently present in the object.
#'
#' @param object A `VISTA` object.
#' @param color_map Named character vector of colors, with names equal to group labels.
#'
#' @return A modified `VISTA` object with updated group colors.
#' @export
set_vista_group_colors <- function(object, color_map) {
  stopifnot(inherits(object, "VISTA"))
  color_map <- .validate_named_color_map(color_map, "color_map")

  md <- S4Vectors::metadata(object)
  group_info <- md$group %||% list()
  group_col <- group_info$column %||% colnames(SummarizedExperiment::colData(object))[1]

  if (is.null(group_col) || !group_col %in% colnames(SummarizedExperiment::colData(object))) {
    cli::cli_abort("Could not resolve grouping column in {.arg object}.")
  }

  groups <- unique(as.character(SummarizedExperiment::colData(object)[[group_col]]))
  missing_groups <- setdiff(groups, names(color_map))
  if (length(missing_groups)) {
    cli::cli_abort("Missing group color(s) for: {.val {missing_groups}}")
  }

  extra_groups <- setdiff(names(color_map), groups)
  if (length(extra_groups)) {
    cli::cli_warn("Ignoring color_map entries not present in object groups: {.val {extra_groups}}")
  }

  group_info$column <- group_col
  group_info$colors <- color_map[groups]
  group_info$palette <- "manual"
  md$group <- group_info
  S4Vectors::metadata(object) <- md

  validate_vista(object, level = "core", error = TRUE)
  object
}

#' Set manual comparison colors in a VISTA object
#'
#' Updates `metadata(x)$comparison$colors` using a user-supplied named color
#' vector. Comparison colors are used in fold-change and comparison-level plots.
#' The color map must include all currently active comparisons.
#'
#' @param object A `VISTA` object.
#' @param color_map Named character vector of colors, with names equal to
#'   comparison names (for example `"A_VS_B"`).
#'
#' @return A modified `VISTA` object with updated comparison colors.
#' @export
set_vista_comparison_colors <- function(object, color_map) {
  stopifnot(inherits(object, "VISTA"))
  color_map <- .validate_named_color_map(color_map, "color_map")

  comp_names <- names(comparisons(object))
  if (!length(comp_names)) {
    cli::cli_abort("No active comparisons found in {.arg object}.")
  }

  missing_comparisons <- setdiff(comp_names, names(color_map))
  if (length(missing_comparisons)) {
    cli::cli_abort("Missing comparison color(s) for: {.val {missing_comparisons}}")
  }

  extra_comparisons <- setdiff(names(color_map), comp_names)
  if (length(extra_comparisons)) {
    cli::cli_warn("Ignoring color_map entries not present in active comparisons: {.val {extra_comparisons}}")
  }

  md <- S4Vectors::metadata(object)
  comp_info <- md$comparison %||% list()
  comp_info$colors <- color_map[comp_names]
  comp_info$palette <- "manual"
  md$comparison <- comp_info
  S4Vectors::metadata(object) <- md

  validate_vista(object, level = "core", error = TRUE)
  object
}

#' Set active DE source in a VISTA object
#'
#' Switches the DE result source used by all downstream VISTA plotting and
#' accessor functions that read the active metadata slots.
#'
#' @param object A `VISTA` object.
#' @param source One of `"active"`, `"deseq2"`, `"edger"`, `"limma"`, or `"consensus"`.
#'   When `"active"`, the currently active DE source is kept.
#'
#' @return A modified `VISTA` object with updated active DE source.
#' @export
set_de_source <- function(object, source = c("deseq2", "edger", "limma", "consensus", "active")) {
  stopifnot(inherits(object, "VISTA"))
  source <- match.arg(source)
  md <- S4Vectors::metadata(object)
  if (identical(source, "active")) {
    source <- md$de_active_source %||% NULL
    if (is.null(source)) {
      # Backward compatibility for older VISTA objects without DE source metadata.
      if (!is.null(md$de_results) && !is.null(md$de_summary)) {
        return(object)
      }
      cli::cli_abort("{.arg source} = {.val active} is not available because this object has no active DE source.")
    }
  }
  if (is.null(md$de_results_by_method) || !source %in% names(md$de_results_by_method)) {
    available <- if (is.null(md$de_results_by_method)) character(0) else names(md$de_results_by_method)
    cli::cli_abort("DE source {.val {source}} is not available in this object. Available: {.val {available}}")
  }

  md$de_results <- md$de_results_by_method[[source]]
  md$de_summary <- md$de_summary_by_method[[source]]
  md$de_active_source <- source
  if (is.list(md$de_cutoffs)) {
    md$de_cutoffs$active_source <- source
  }
  S4Vectors::metadata(object) <- md
  validate_vista(object, level = "full", error = TRUE)
  object
}
