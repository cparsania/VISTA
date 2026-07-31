# Raw counts -----------------------------------------------------------------
#
# Before 1.2.0 a VISTA object kept only normalized counts, so counts(v) failed
# and the object could not be handed back to DESeq2 or edgeR -- an odd gap for a
# package that wraps them. The raw matrix is now retained as a second assay
# named "counts" (registered after norm_counts, so an unqualified assay(x) is
# unchanged).

#' Raw counts from a VISTA object
#'
#' @description
#' Returns the filtered raw count matrix stored alongside the normalized
#' counts. Objects built before VISTA 1.2.0, or with
#' `create_vista(keep_raw_counts = FALSE)`, do not carry one.
#'
#' @details
#' Normalized counts are not invertible, so [BiocGenerics::updateObject()]
#' cannot back-fill this assay for an older object; it has to be rebuilt with
#' [create_vista()].
#'
#' @param object A `VISTA` object.
#' @param ... Ignored.
#'
#' @return A numeric matrix of raw counts with genes in rows and samples in
#'   columns.
#' @examples
#' v <- example_vista()
#' counts(v)[seq_len(3), seq_len(3)]
#'
#' @name counts
#' @aliases counts,VISTA-method
#' @importFrom BiocGenerics counts
#' @exportMethod counts
setMethod("counts", "VISTA", function(object, ...) {
  an <- SummarizedExperiment::assayNames(object)
  if (!"counts" %in% an) {
    cli::cli_abort(c(
      "This VISTA object does not carry a raw {.field counts} assay.",
      "i" = "Objects built before VISTA 1.2.0, or with {.code keep_raw_counts = FALSE}, store only normalized counts.",
      "x" = "Normalized counts cannot be inverted, so this cannot be recovered by {.fn updateObject}.",
      "i" = "Rebuild the object with {.fn create_vista} to obtain raw counts."
    ))
  }
  SummarizedExperiment::assay(object, "counts")
})

#' Convert a VISTA object to a DESeqDataSet
#'
#' @description
#' Builds a `DESeq2::DESeqDataSet` from the raw counts retained in a `VISTA`
#' object, so an analysis can be continued or re-run with DESeq2 directly.
#'
#' @details
#' This is deliberately explicit rather than making
#' `DESeq2::DESeqDataSet(v, design)` work implicitly: the inherited
#' `SummarizedExperiment` method would pick up whichever assay comes first,
#' which is `norm_counts`, and fail because those values are not integers.
#'
#' Counts are rounded to integers, as DESeq2 requires.
#'
#' @param x A `VISTA` object carrying a raw `counts` assay.
#' @param design A model formula. Defaults to `~ <group column>` using the
#'   grouping column stored in the object.
#'
#' @return A `DESeqDataSet`.
#' @examples
#' v <- example_vista()
#' if (requireNamespace("DESeq2", quietly = TRUE)) {
#'   dds <- as_deseq_dataset(v)
#'   class(dds)
#' }
#' @export
as_deseq_dataset <- function(x, design = NULL) {
  stopifnot(inherits(x, "VISTA"))
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg DESeq2} is required to build a {.cls DESeqDataSet}.")
  }

  raw <- counts(x)
  storage.mode(raw) <- "double"
  raw <- round(raw)

  group_col <- S4Vectors::metadata(x)$group$column
  if (is.null(design)) {
    if (is.null(group_col)) {
      cli::cli_abort("No grouping column stored in {.arg x}; supply {.arg design} explicitly.")
    }
    design <- stats::as.formula(paste("~", group_col))
  }
  design <- .coerce_design_formula(design)

  col_data <- SummarizedExperiment::colData(x)
  for (v in intersect(all.vars(design), colnames(col_data))) {
    if (is.character(col_data[[v]])) col_data[[v]] <- factor(col_data[[v]])
  }

  DESeq2::DESeqDataSetFromMatrix(
    countData = raw,
    colData = col_data,
    design = design
  )
}
