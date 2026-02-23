#' @title VISTA S4 Class Definition
#'
#' @description The `VISTA` class extends `SummarizedExperiment`. All analysis
#' metadata (DE results, summaries, cutoffs, grouping info, etc.) is stored in
#' the `metadata()` slot of the object.
#'
#' @details Core elements stored in `metadata(v)`:
#' \itemize{
#'   \item \code{$de_results}: named \code{SimpleList} of DE result tables.
#'   \item \code{$de_summary}: named \code{SimpleList} of DEG summary tables.
#'   \item \code{$de_cutoffs}: named list of thresholds (log2FC, p-value type/cutoff, method settings).
#'   \item \code{$group}: list with \code{column}, \code{palette}, \code{colors} describing the grouping/fill scheme.
#'   \item \code{$provenance}: list with constructor version, timestamp, session info.
#'   \item \code{$vista_schema_version}: metadata schema tag used for compatibility checks.
#' }
#' Feature-level annotations live in \code{rowData(v)}, sample metadata in
#' \code{colData(v)}, and normalized counts (and any additional assays) in
#' \code{assay(v, "norm_counts")} by default. Use \code{create_vista()} as the
#' primary end-user constructor. Advanced users can convert an existing
#' \code{SummarizedExperiment} with \code{as_vista()}.
#' @name VISTA-class
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @exportClass VISTA
#' @seealso \code{\link{create_vista}}, \code{\link{as_vista}}, \code{\link{validate_vista}}
NULL

#' @keywords internal
#' @noRd
.VISTA_SCHEMA_VERSION <- "1.0.0"

.vista_validity_messages <- function(object) {
  msgs <- character()

  an <- SummarizedExperiment::assayNames(object)
  if (!"norm_counts" %in% an) {
    msgs <- c(msgs, "Assay 'norm_counts' is required.")
    return(msgs)
  }

  mat <- SummarizedExperiment::assay(object, "norm_counts")
  if (!is.numeric(mat)) {
    msgs <- c(msgs, "Assay 'norm_counts' must be numeric.")
  }

  rn <- rownames(mat)
  cn <- colnames(mat)
  if (is.null(rn) || any(!nzchar(rn)) || anyDuplicated(rn)) {
    msgs <- c(msgs, "'norm_counts' rownames must be unique, non-empty gene IDs.")
  }
  if (is.null(cn) || any(!nzchar(cn)) || anyDuplicated(cn)) {
    msgs <- c(msgs, "'norm_counts' colnames must be unique, non-empty sample IDs.")
  }

  cd <- as.data.frame(SummarizedExperiment::colData(object), stringsAsFactors = FALSE)
  rd <- as.data.frame(SummarizedExperiment::rowData(object), stringsAsFactors = FALSE)
  cd_rn <- rownames(cd)
  rd_rn <- rownames(rd)

  if (is.null(cd_rn) || any(!nzchar(cd_rn)) || anyDuplicated(cd_rn)) {
    msgs <- c(msgs, "colData rownames must be unique, non-empty sample IDs.")
  }
  if (is.null(rd_rn) || any(!nzchar(rd_rn)) || anyDuplicated(rd_rn)) {
    msgs <- c(msgs, "rowData rownames must be unique, non-empty feature IDs.")
  }

  if (!is.null(cn) && !is.null(cd_rn) && !identical(cd_rn, cn)) {
    msgs <- c(msgs, "colData rownames must exactly match colnames(norm_counts).")
  }
  if (!is.null(rn) && !is.null(rd_rn) && !identical(rd_rn, rn)) {
    msgs <- c(msgs, "rowData rownames must exactly match rownames(norm_counts).")
  }

  md <- S4Vectors::metadata(object)
  required_md <- c("de_results", "de_summary", "de_cutoffs", "group", "provenance", "vista_schema_version")
  missing_md <- setdiff(required_md, names(md))
  if (length(missing_md) > 0) {
    msgs <- c(msgs, paste0("Missing metadata key(s): ", paste(missing_md, collapse = ", "), "."))
  }

  if ("group" %in% names(md) && is.list(md$group)) {
    gcol <- md$group$column
    if (!is.character(gcol) || length(gcol) != 1L || !nzchar(gcol)) {
      msgs <- c(msgs, "metadata(x)$group$column must be a non-empty character scalar.")
    } else if (!gcol %in% colnames(cd)) {
      msgs <- c(msgs, sprintf("Grouping column '%s' is not present in colData.", gcol))
    } else {
      grp <- cd[[gcol]]
      if (anyNA(grp)) {
        msgs <- c(msgs, sprintf("Grouping column '%s' contains NA values.", gcol))
      }
      gcols <- md$group$colors
      if (!is.null(gcols)) {
        if (!is.character(gcols) || is.null(names(gcols))) {
          msgs <- c(msgs, "metadata(x)$group$colors must be a named character vector.")
        } else {
          missing_colors <- setdiff(unique(as.character(grp)), names(gcols))
          if (length(missing_colors) > 0) {
            msgs <- c(msgs, paste0("Missing group color(s) for: ", paste(missing_colors, collapse = ", "), "."))
          }
        }
      }
    }
  } else if ("group" %in% names(md)) {
    msgs <- c(msgs, "metadata(x)$group must be a list.")
  }

  if (!is.null(md$vista_schema_version) &&
      (!is.character(md$vista_schema_version) || length(md$vista_schema_version) != 1L || !nzchar(md$vista_schema_version))) {
    msgs <- c(msgs, "metadata(x)$vista_schema_version must be a non-empty character scalar.")
  }

  if (length(msgs) == 0) TRUE else msgs
}

setClass(
  Class = "VISTA",
  contains = "SummarizedExperiment"
)

setValidity("VISTA", .vista_validity_messages)


#' @title Coerce SummarizedExperiment to VISTA
#' @description Convert a pre-existing `SummarizedExperiment` into a valid
#'   `VISTA` object.
#'
#' @param se A `SummarizedExperiment` object.
#' @param assay_name Assay in `se` to use as `norm_counts` (default: `"norm_counts"`).
#' @param group_column Grouping column in `colData(se)`.
#' @param comparisons Optional named list of differential expression tables.
#' @param deg_summary Optional named list of DEG summary tables.
#' @param cutoffs Optional named list of threshold metadata.
#' @param group_palette Palette name for group colors.
#' @param validate Logical; run object validation before returning.
#'
#' @return A `VISTA` object.
#' @examples
#' mat <- matrix(rnorm(60), nrow = 10)
#' rownames(mat) <- paste0("gene", seq_len(nrow(mat)))
#' colnames(mat) <- paste0("sample", seq_len(ncol(mat)))
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(norm_counts = mat),
#'   colData = S4Vectors::DataFrame(
#'     cond = rep(c("A", "B"), each = 3),
#'     row.names = colnames(mat)
#'   ),
#'   rowData = S4Vectors::DataFrame(
#'     gene_id = rownames(mat),
#'     row.names = rownames(mat)
#'   )
#' )
#' v <- as_vista(se, group_column = "cond")
#' v
#' @import methods
#' @export
as_vista <- function(se,
                     assay_name = "norm_counts",
                     group_column,
                     comparisons = list(),
                     deg_summary = list(),
                     cutoffs = list(),
                     group_palette = "Dark 3",
                     validate = TRUE) {
  if (!inherits(se, "SummarizedExperiment")) {
    cli::cli_abort("{.arg se} must inherit from {.cls SummarizedExperiment}.")
  }
  if (!assay_name %in% SummarizedExperiment::assayNames(se)) {
    cli::cli_abort(
      "Assay {.val {assay_name}} not found. Available assays: {.val {SummarizedExperiment::assayNames(se)}}"
    )
  }

  mat <- SummarizedExperiment::assay(se, assay_name)
  if (!is.numeric(mat)) {
    cli::cli_abort("Assay {.val {assay_name}} must be numeric.")
  }
  if (is.null(rownames(mat)) || any(!nzchar(rownames(mat))) || anyDuplicated(rownames(mat))) {
    cli::cli_abort("Assay rownames must be unique, non-empty gene identifiers.")
  }
  if (is.null(colnames(mat)) || any(!nzchar(colnames(mat))) || anyDuplicated(colnames(mat))) {
    cli::cli_abort("Assay colnames must be unique, non-empty sample identifiers.")
  }

  sample_info <- as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)
  if (is.null(rownames(sample_info)) || any(!nzchar(rownames(sample_info))) || anyDuplicated(rownames(sample_info))) {
    rownames(sample_info) <- colnames(mat)
  }
  if (!identical(rownames(sample_info), colnames(mat))) {
    if (all(colnames(mat) %in% rownames(sample_info))) {
      sample_info <- sample_info[colnames(mat), , drop = FALSE]
    } else {
      cli::cli_abort("`colData(se)` rownames must match assay column names.")
    }
  }
  sample_info$sample_names <- rownames(sample_info)

  row_data <- as.data.frame(SummarizedExperiment::rowData(se), stringsAsFactors = FALSE)
  if (is.null(rownames(row_data)) || any(!nzchar(rownames(row_data))) || anyDuplicated(rownames(row_data))) {
    rownames(row_data) <- rownames(mat)
  }
  if (!identical(rownames(row_data), rownames(mat))) {
    if (all(rownames(mat) %in% rownames(row_data))) {
      row_data <- row_data[rownames(mat), , drop = FALSE]
    } else {
      cli::cli_abort("`rowData(se)` rownames must match assay row names.")
    }
  }
  if (!"gene_id" %in% colnames(row_data)) {
    row_data$gene_id <- rownames(row_data)
  }

  .vista(
    norm_counts = mat,
    sample_info = sample_info,
    row_data = row_data,
    comparisons = comparisons,
    deg_summary = deg_summary,
    cutoffs = cutoffs,
    group_column = group_column,
    group_palette = group_palette,
    validate = validate
  )
}

#' @keywords internal
#' @noRd
.vista <- function(norm_counts,
                   sample_info,
                   row_data,
                   comparisons = list(),
                   deg_summary = list(),
                   cutoffs = list(),
                   group_column = character(1),
                   group_palette = "Dark 3",
                   validate = TRUE) {
  # ---- Basic structure checks ----
  if (is.null(rownames(norm_counts)))
    cli::cli_abort("norm_counts must have rownames.")
  if (is.null(colnames(norm_counts)))
    cli::cli_abort("norm_counts must have colnames (sample names).")
  if (is.null(rownames(sample_info)))
    cli::cli_abort("sample_info must have rownames (sample names).")
  if (is.null(rownames(row_data)))
    cli::cli_abort("row_data must have rownames (feature IDs).")

  # match dimensions / identities
  if (!identical(rownames(row_data), rownames(norm_counts))) {
    cli::cli_abort("row_data rownames must exactly match rownames(norm_counts).")
  }
  if (!identical(rownames(sample_info), colnames(norm_counts))) {
    cli::cli_abort("sample_info rownames must exactly match colnames(norm_counts).")
  }

  # ---- group_column (required) ----
  if (length(group_column) != 1L || is.na(group_column) || nchar(group_column) == 0)
    cli::cli_abort("'group_column' is required and cannot be empty.")
  if (!group_column %in% colnames(sample_info)) {
    cli::cli_abort("The specified group_column '{group_column}' is not found in sample_info.")
  }

  # ---- palette selection ----
  palette_choices <- c("Pastel 1", "Dark 2", "Dark 3", "Set 2", "Set 3",
                       "Warm", "Cold", "Harmonic", "Dynamic")
  group_palette <- match.arg(group_palette, palette_choices)

  # ---- compute group colors ----
  grp <- sample_info[[group_column]]
  if (anyNA(grp)) {
    cli::cli_abort("group_column contains NA values; please fix sample_info.")
  }
  groups <- unique(as.character(grp))
  group_colors <- colorspace::qualitative_hcl(length(groups), palette = group_palette)
  names(group_colors) <- groups

  # ---- build SE core ----
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(norm_counts = as.matrix(norm_counts)),
    colData = S4Vectors::DataFrame(sample_info),
    rowData = S4Vectors::DataFrame(row_data)
  )

  # ---- coerce lists for metadata ----
  de_results_sl <- if (length(comparisons)) S4Vectors::SimpleList(comparisons) else S4Vectors::SimpleList()
  de_summary_sl <- if (length(deg_summary)) S4Vectors::SimpleList(deg_summary) else S4Vectors::SimpleList()

  # validate each DE table's rownames align with features
  if (length(de_results_sl)) {
    rn <- rownames(se)
    for (nm in names(de_results_sl)) {
      dfi <- de_results_sl[[nm]]
      if (is.null(rownames(dfi)) || !identical(rownames(dfi), rn)) {
        cli::cli_abort("Each DE result ('{nm}') must have rownames identical to rownames(norm_counts).")
      }
    }
  }

  S4Vectors::metadata(se) <- list(
    de_results = de_results_sl,
    de_summary = de_summary_sl,
    de_cutoffs = cutoffs,
    group = list(
      column = group_column,
      palette = group_palette,
      colors = group_colors
    ),
    provenance = list(
      constructor = ".vista()",
      created_at = as.character(Sys.time())
    ),
    vista_schema_version = .VISTA_SCHEMA_VERSION
  )

  v <- methods::new("VISTA", se)
  if (isTRUE(validate)) {
    validate_vista(v, level = "core", error = TRUE)
  }
  v
}

# Backward-compatibility internal alias; not exported.
#' @keywords internal
#' @noRd
vista <- function(...) {
  .vista(...)
}
