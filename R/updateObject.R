# Schema versioning ----------------------------------------------------------
#
# metadata(x)$vista_schema_version records the metadata layout an object was
# built with. Before 1.2.0 it was written once and never read: nothing compared
# it to .VISTA_SCHEMA_VERSION and there was no migration path, so an object
# serialized under an older layout deserialized silently and could be missing
# keys the current code expects.

#' Compare an object's schema version to the running package
#'
#' @param x A `VISTA` object.
#' @return One of `"current"`, `"older"`, `"newer"`, or `"unknown"`.
#' @keywords internal
#' @noRd
.vista_schema_compare <- function(x) {
  stored <- S4Vectors::metadata(x)$vista_schema_version
  if (is.null(stored) || !is.character(stored) || length(stored) != 1L || !nzchar(stored)) {
    return("unknown")
  }
  stored_v <- tryCatch(package_version(stored), error = function(e) NULL)
  if (is.null(stored_v)) return("unknown")

  current_v <- package_version(.VISTA_SCHEMA_VERSION)
  if (stored_v == current_v) "current" else if (stored_v < current_v) "older" else "newer"
}

#' Update a VISTA object to the current metadata schema
#'
#' @description
#' Brings a `VISTA` object created by an older version of the package up to the
#' current metadata layout. Missing metadata keys are back-filled with their
#' documented defaults, the schema tag is stamped to the running version, and
#' the migration is recorded in `metadata(x)$provenance$updates`.
#'
#' @details
#' Some content cannot be recovered by migration. In particular, objects built
#' before VISTA 1.2.0 did not retain raw counts, and normalized counts are not
#' invertible -- [counts()] on such an object reports that a rebuild is
#' required rather than inventing values.
#'
#' Objects carrying a *newer* schema than the running package are left
#' untouched and reported by [validate_vista()] as an issue, because newer
#' metadata may carry semantics this version would misread.
#'
#' @param object A `VISTA` object.
#' @param ... Passed to other methods (unused).
#' @param verbose Logical; report what was migrated.
#'
#' @return A `VISTA` object stamped with the current schema version.
#' @examples
#' v <- example_vista()
#'
#' # Simulate an object written by an older release.
#' S4Vectors::metadata(v)$vista_schema_version <- "0.9.0"
#' v2 <- BiocGenerics::updateObject(v)
#' S4Vectors::metadata(v2)$vista_schema_version
#'
#' @importFrom BiocGenerics updateObject
#' @exportMethod updateObject
setMethod("updateObject", "VISTA", function(object, ..., verbose = FALSE) {
  state <- .vista_schema_compare(object)
  md <- S4Vectors::metadata(object)
  before <- md$vista_schema_version

  current_schema <- .VISTA_SCHEMA_VERSION

  if (identical(state, "newer")) {
    cli::cli_warn(c(
      "This object's schema ({.val {before}}) is newer than VISTA's ({.val {current_schema}}).",
      "i" = "Update the VISTA package rather than the object."
    ))
    return(object)
  }

  defaults <- list(
    de_results = S4Vectors::SimpleList(),
    de_summary = S4Vectors::SimpleList(),
    de_cutoffs = list(),
    provenance = list(),
    group = list(
      column = colnames(SummarizedExperiment::colData(object))[1],
      palette = .vista_default("group_palette"),
      colors = NULL
    )
  )
  filled <- names(defaults)[vapply(names(defaults), function(k) is.null(md[[k]]), logical(1))]
  for (key in filled) md[[key]] <- defaults[[key]]

  prov <- md$provenance
  if (!is.list(prov)) prov <- list()
  prov$updates <- c(
    prov$updates,
    list(list(
      from = before %||% NA_character_,
      to = .VISTA_SCHEMA_VERSION,
      filled = filled,
      at = as.character(Sys.time())
    ))
  )
  md$provenance <- prov
  md$vista_schema_version <- .VISTA_SCHEMA_VERSION

  S4Vectors::metadata(object) <- md

  if (isTRUE(verbose)) {
    cli::cli_inform(c(
      "Updated VISTA schema {.val {before %||% 'unknown'}} -> {.val {current_schema}}.",
      if (length(filled)) "i" = "Back-filled metadata key{?s}: {.field {filled}}."
    ))
  }

  object
})
