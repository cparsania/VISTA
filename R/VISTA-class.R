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
#' }
#' Feature-level annotations live in \code{rowData(v)}, sample metadata in
#' \code{colData(v)}, and normalized counts (and any additional assays) in
#' \code{assay(v, "norm_counts")} by default. Use \code{create_vista()} for
#' standard construction from counts + sample_info; use \code{vista()} for
#' low-level instantiation when you already have normalized assays/metadata.
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @exportClass VISTA
setClass(
  Class = "VISTA",
  contains = "SummarizedExperiment"
)




#' @title VISTA constructor
#' @description Create a VISTA object that extends `SummarizedExperiment`.
#'
#' @param norm_counts A normalized expression matrix (or \code{data.frame})
#'   with genes as rows and samples as columns. Must have rownames.
#' @param sample_info A \code{data.frame} with sample metadata. Row names must
#'   match \code{colnames(norm_counts)}.
#' @param row_data A \code{data.frame} with feature metadata. Row names must
#'   match \code{rownames(norm_counts)}.
#' @param comparisons A named list of differential expression result tables
#'   (e.g. \code{S4Vectors::DataFrame}s) keyed by comparison name. Default: empty list.
#' @param deg_summary A named list of per-comparison DEG summary tables. Default: empty list.
#' @param cutoffs A named list of thresholds used for DEG classification (e.g. \code{alpha}, \code{lfc}). Default: empty list.
#' @param group_column A **non-empty** character scalar naming the grouping column in \code{sample_info}.
#'   This is **required** and must exist in \code{sample_info}.
#' @param group_palette A colorspace qualitative palette name used to assign colors to groups in \code{group_column}.
#'   One of \code{c("Pastel 1","Dark 2","Dark 3","Set 2","Set 3","Warm","Cold","Harmonic","Dynamic")}.
#'   Default: \code{"Dark 3"}.
#'
#' @return A \code{VISTA} object. The following keys are stored in \code{metadata(vista)}:
#' \itemize{
#'   \item \code{$de_results}: named \code{SimpleList} of DE result tables.
#'   \item \code{$de_summary}: named \code{SimpleList} of DEG summary tables.
#'   \item \code{$de_cutoffs}: named list of thresholds.
#'   \item \code{$group}: list with \code{column}, \code{palette}, \code{colors}.
#'   \item \code{$provenance}: list with constructor version, timestamp, session info.
#' }
#'
#' @examples
#' # Low-level constructor (advanced usage)
#' # Most users should use create_vista() instead
#'
#' # Create sample data
#' norm_counts_mat <- matrix(
#'   rnorm(100), nrow = 10, ncol = 10,
#'   dimnames = list(
#'     paste0("gene", 1:10),
#'     paste0("sample", 1:10)
#'   )
#' )
#'
#' sample_info_df <- data.frame(
#'   groups = rep(c("A", "B"), each = 5),
#'   row.names = colnames(norm_counts_mat)
#' )
#'
#' row_data_df <- data.frame(
#'   gene_id = rownames(norm_counts_mat),
#'   row.names = rownames(norm_counts_mat)
#' )
#'
#' # Create VISTA object
#' v <- vista(
#'   norm_counts = norm_counts_mat,
#'   sample_info = sample_info_df,
#'   row_data = row_data_df,
#'   group_column = "groups"
#' )
#'
#' # Examine object
#' v
#' norm_counts(v)[1:5, 1:5]
#'
#' @import SummarizedExperiment S4Vectors methods
#' @export
vista <- function(norm_counts,
                  sample_info,
                  row_data,
                  comparisons  = list(),
                  deg_summary  = list(),
                  cutoffs      = list(),
                  group_column = character(1),
                  group_palette = "Dark 3") {



  # ---- Basic structure checks ----
  if (is.null(rownames(norm_counts)))
    stop("norm_counts must have rownames.")
  if (is.null(colnames(norm_counts)))
    stop("norm_counts must have colnames (sample names).")
  if (is.null(rownames(sample_info)))
    stop("sample_info must have rownames (sample names).")
  if (is.null(rownames(row_data)))
    stop("row_data must have rownames (feature IDs).")

  # match dimensions / identities
  if (!identical(rownames(row_data), rownames(norm_counts))) {
    stop("row_data rownames must exactly match rownames(norm_counts).")
  }
  if (!identical(rownames(sample_info), colnames(norm_counts))) {
    stop("sample_info rownames must exactly match colnames(norm_counts).")
  }

  # ---- group_column (required) ----
  if (length(group_column) != 1L || is.na(group_column) || nchar(group_column) == 0)
    stop("'group_column' is required and cannot be empty.")
  if (!group_column %in% colnames(sample_info)) {
    stop("The specified group_column '", group_column, "' is not found in sample_info.")
  }

  # ---- palette selection ----
  palette_choices <- c("Pastel 1", "Dark 2", "Dark 3", "Set 2", "Set 3",
                       "Warm", "Cold", "Harmonic", "Dynamic")
  group_palette <- match.arg(group_palette, palette_choices)

  # ---- compute group colors ----
  grp <- sample_info[[group_column]]
  if (anyNA(grp))
    stop("group_column contains NA values; please fix sample_info.")
  groups <- unique(as.character(grp))
  group_colors <- colorspace::qualitative_hcl(length(groups), palette = group_palette)
  names(group_colors) <- groups

  # ---- build SE core ----
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(norm_counts = as.matrix(norm_counts)),
    colData = S4Vectors::DataFrame(sample_info),
    rowData = S4Vectors::DataFrame(row_data)
  )

  # ---- coerce lists for metadata ----
  # Prefer SimpleList so it's light and plays nicely with Bioc serialization
  de_results_sl <- if (length(comparisons)) S4Vectors::SimpleList(comparisons) else S4Vectors::SimpleList()
  de_summary_sl <- if (length(deg_summary)) S4Vectors::SimpleList(deg_summary) else S4Vectors::SimpleList()

  # (optional) validate each DE table's rownames align with features
  if (length(de_results_sl)) {
    rn <- rownames(se)
    for (nm in names(de_results_sl)) {
      dfi <- de_results_sl[[nm]]
      if (is.null(rownames(dfi)) || !identical(rownames(dfi), rn)) {
        stop("Each DE result ('", nm, "') must have rownames identical to rownames(norm_counts).")
      }
    }
  }

  # ---- stuff everything into metadata() ----
  S4Vectors::metadata(se) <- list(
    de_results = de_results_sl,
    de_summary = de_summary_sl,
    de_cutoffs = cutoffs,
    group = list(
      column  = group_column,
      palette = group_palette,
      colors  = group_colors
    ),
    provenance = list(
      constructor = "vista()",
      created_at = as.character(Sys.time())
    )
  )

  # ---- return as VISTA ----
  methods::new("VISTA", se)
}
