# Display --------------------------------------------------------------------
#
# S4 objects display through show(). VISTA defined no show() method, so it fell
# back to SummarizedExperiment's, and every VISTA object introduced itself as
# "class: SummarizedExperiment" with none of its analysis state visible.

#' @keywords internal
#' @noRd
.vista_show_summary_lines <- function(object) {
  md <- S4Vectors::metadata(object)
  lines <- character()

  add <- function(label, value) {
    if (length(value) && nzchar(paste(value, collapse = ""))) {
      lines <<- c(lines, sprintf("%s: %s", label, paste(value, collapse = ", ")))
    }
  }

  group_info <- md$group
  if (is.list(group_info) && !is.null(group_info$column)) {
    cd <- SummarizedExperiment::colData(object)
    levs <- if (group_info$column %in% colnames(cd)) {
      unique(as.character(cd[[group_info$column]]))
    } else {
      character()
    }
    add("group column", sprintf("%s (%s)", group_info$column, paste(levs, collapse = ", ")))
  }

  comps <- tryCatch(names(comparisons(object)), error = function(e) NULL)
  if (length(comps)) {
    shown <- utils::head(comps, 3)
    more <- length(comps) - length(shown)
    add(
      "comparisons",
      sprintf("%s%s", paste(shown, collapse = ", "), if (more > 0) sprintf(" (+%d more)", more) else "")
    )
  } else {
    add("comparisons", "none")
  }

  if (!is.null(md$de_active_source)) {
    available <- md$de_available_sources
    add(
      "DE source",
      sprintf(
        "%s%s", md$de_active_source,
        if (length(available) > 1L) sprintf(" (of %s)", paste(available, collapse = ", ")) else ""
      )
    )
  }

  cuts <- md$de_cutoffs
  if (is.list(cuts) && length(cuts)) {
    add("cutoffs", sprintf(
      "|log2FC| >= %s, %s <= %s",
      cuts$log2fc %||% "?", cuts$p_value_type %||% "p", cuts$pval %||% "?"
    ))
  }

  an <- SummarizedExperiment::assayNames(object)
  if ("counts" %in% an) add("raw counts", "available via counts()")

  add("schema", md$vista_schema_version %||% "unknown")
  lines
}

#' Display a VISTA object
#'
#' Prints the `SummarizedExperiment` header followed by a short summary of the
#' analysis state -- grouping, comparisons, active differential-expression
#' source, thresholds and schema version -- none of which was visible before.
#'
#' @param object A `VISTA` object.
#' @return `object`, invisibly. Called for its output.
#' @examples
#' v <- example_vista()
#' v
#' @name VISTA-show
#' @aliases show,VISTA-method
#' @exportMethod show
setMethod("show", "VISTA", function(object) {
  se_show <- methods::selectMethod("show", "SummarizedExperiment")
  out <- utils::capture.output(se_show(methods::as(object, "SummarizedExperiment")))
  # The inherited header names the parent class; this object is a VISTA.
  out[1] <- sub("^class: SummarizedExperiment", "class: VISTA", out[1])
  cat(out, sep = "\n")

  extra <- .vista_show_summary_lines(object)
  if (length(extra)) {
    cat("-------- VISTA --------", extra, sep = "\n")
  }
  invisible(object)
})
