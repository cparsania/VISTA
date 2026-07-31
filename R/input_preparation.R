#' Read and standardize count inputs for VISTA
#'
#' `read_vista_counts()` helps standardize common RNA-seq count inputs into a
#' count table that can be passed directly to `create_vista()`. It supports
#' plain matrices/data frames, featureCounts outputs, STAR gene counts,
#' HTSeq-count outputs, tximport-like lists, and RSEM gene result files.
#'
#' Internally, VISTA uses a format-specific importer for each supported input
#' type, then normalizes the result into a common structure with:
#' \itemize{
#'   \item a count table with a `gene_id` column plus sample columns
#'   \item optional feature metadata in `row_data`
#'   \item sample names inferred from columns or file names
#'   \item an auditable `sample_name_map` showing original and repaired names
#' }
#'
#' @param x Count input. Supported values depend on `format` and include a
#'   matrix, data frame, single file path, vector of file paths, or a
#'   tximport-like list with `counts`, `abundance`, and/or `length`.
#' @param format Input format. One of `"auto"`, `"matrix"`,
#'   `"featurecounts"`, `"star"`, `"htseq"`, `"tximport"`, or `"rsem"`.
#' @param gene_id_column Optional gene identifier column in tabular inputs.
#'   When omitted, VISTA uses common names such as `gene_id`/`Geneid`, or falls
#'   back to rownames for matrices/data frames with unique rownames.
#' @param sample_columns Optional character vector of sample count columns to
#'   retain from tabular inputs.
#' @param sample_names Optional sample names to use when `x` is a vector of
#'   per-sample files.
#' @param annotation_columns Optional feature annotation columns to retain in the
#'   returned `row_data`.
#' @param count_column Optional count column selector for formats that expose
#'   multiple count choices. For STAR, use one of `"unstranded"`,
#'   `"stranded_first"`, or `"stranded_second"`. For RSEM, this defaults to
#'   `"expected_count"`.
#' @param tx2gene Optional two-column mapping used to summarize transcript-level
#'   tximport-like inputs to genes. The first column should contain transcript
#'   IDs and the second column gene IDs.
#' @param counts_from Which matrix to extract from a tximport-like input:
#'   `"counts"`, `"abundance"`, or `"length"`.
#' @param drop_technical Logical; when `TRUE`, drop known technical summary rows
#'   from STAR/HTSeq inputs.
#' @param remove_special_rows Logical; alias for `drop_technical`, retained for
#'   clarity in file-based imports.
#' @param make_unique_ids Logical; if `TRUE`, duplicate gene IDs are repaired
#'   with [make.unique()]. Otherwise duplicated gene IDs raise an error.
#' @param repair_sample_names Strategy for repairing sample column names.
#'   `"auto"` (default) strips common file-path and alignment/count suffixes
#'   when the repaired names are unique, while `"none"` leaves sample columns
#'   unchanged. In automatic mode VISTA currently:
#'   \itemize{
#'     \item strips directory paths to the basename
#'     \item uses the parent directory for generic quantification files such as
#'       `quant.sf` or `abundance.tsv`
#'     \item removes common RNA-seq output suffixes such as
#'       `Aligned.sortedByCoord.out.bam`, `ReadsPerGene.out.tab`,
#'       `.genes.results`, `.isoforms.results`, `.bam`, and `.fastq.gz`
#'     \item removes common lane/read suffixes such as `_S1_L001_R1_001`,
#'       `_L001_R2_001`, `_R1`, and `_R2`
#'   }
#'   Repaired names are only applied when they remain non-empty and unique.
#'   Otherwise VISTA keeps the original count column names and records the
#'   unchanged mapping in `sample_name_map`.
#' @param return_type Return `"list"` (default), standardized `"data.frame"`,
#'   or numeric `"matrix"`.
#' @param verbose Logical; print an informational import summary.
#'
#' @importFrom cli cli_abort cli_inform cli_warn
#' @importFrom rlang `%||%`
#'
#' @return If `return_type = "list"`, a list with:
#' \describe{
#'   \item{counts}{A standardized count table with `gene_id` plus sample columns.}
#'   \item{row_data}{Feature metadata aligned to the count table.}
#'   \item{column_geneid}{Always `"gene_id"` for the standardized output.}
#'   \item{sample_names}{Sample columns in the standardized count table.}
#'   \item{sample_name_map}{A two-column mapping of original and repaired sample names.}
#'   \item{input_format}{Resolved import format.}
#'   \item{report}{Basic import summary.}
#' }
#'
#' If `return_type = "data.frame"`, returns the standardized count table. If
#' `return_type = "matrix"`, returns a numeric matrix with gene IDs as rownames.
#'
#' @examples
#' data("count_data", package = "VISTA")
#'
#' cnt <- read_vista_counts(
#'   count_data[seq_len(25), ],
#'   format = "matrix",
#'   gene_id_column = "gene_id"
#' )
#'
#' head(cnt$counts[, seq_len(4)])
#' cnt$sample_names
#' @export
read_vista_counts <- function(x,
                              format = c("auto", "matrix", "featurecounts", "star", "htseq", "tximport", "rsem"),
                              gene_id_column = NULL,
                              sample_columns = NULL,
                              sample_names = NULL,
                              annotation_columns = NULL,
                              count_column = NULL,
                              tx2gene = NULL,
                              counts_from = c("counts", "abundance", "length"),
                              drop_technical = TRUE,
                              remove_special_rows = TRUE,
                              make_unique_ids = FALSE,
                              repair_sample_names = c("auto", "none"),
                              return_type = c("list", "data.frame", "matrix"),
                              verbose = TRUE) {
  format <- match.arg(format)
  counts_from <- match.arg(counts_from)
  repair_sample_names <- match.arg(repair_sample_names)
  return_type <- match.arg(return_type)
  drop_technical <- isTRUE(drop_technical) && isTRUE(remove_special_rows)

  resolved_format <- .detect_vista_count_format(x, format = format)
  parsed <- switch(
    resolved_format,
    matrix = .import_counts_matrix(
      x = x,
      gene_id_column = gene_id_column,
      sample_columns = sample_columns,
      annotation_columns = annotation_columns
    ),
    featurecounts = .import_counts_featurecounts(
      x = x,
      gene_id_column = gene_id_column,
      sample_columns = sample_columns,
      annotation_columns = annotation_columns
    ),
    star = .import_counts_star(
      x = x,
      sample_names = sample_names,
      count_column = count_column,
      drop_technical = drop_technical
    ),
    htseq = .import_counts_htseq(
      x = x,
      sample_names = sample_names,
      drop_technical = drop_technical
    ),
    tximport = .import_counts_tximport(
      x = x,
      tx2gene = tx2gene,
      counts_from = counts_from
    ),
    rsem = .import_counts_rsem(
      x = x,
      sample_names = sample_names,
      count_column = count_column,
      annotation_columns = annotation_columns
    ),
    cli::cli_abort("Unsupported {.arg format}: {.val {resolved_format}}.")
  )

  out <- .standardize_counts_import(
    parsed = parsed,
    input_format = resolved_format,
    make_unique_ids = make_unique_ids,
    repair_sample_names = repair_sample_names
  )

  if (isTRUE(verbose)) {
    cli::cli_inform(
      "Imported {.field {nrow(out$counts)}} features and {.field {length(out$sample_names)}} samples from {.val {resolved_format}} input."
    )
  }

  switch(
    return_type,
    list = out,
    `data.frame` = out$counts,
    matrix = .counts_df_to_matrix(out$counts, column_geneid = out$column_geneid)
  )
}

#' Read and standardize sample metadata for VISTA
#'
#' `read_vista_metadata()` standardizes a sample sheet for use as
#' `sample_info` in `create_vista()`. It infers or creates the required
#' `sample_names` column using the same conventions VISTA already accepts in the
#' constructor.
#'
#' @param x Sample metadata as a data frame or file path.
#' @param sample_column Optional column to use as `sample_names`. If omitted,
#'   VISTA uses an existing `sample_names` column, non-default rownames, or
#'   common aliases such as `sample`, `sample_id`, or `Run`.
#' @param required_columns Optional character vector of columns that must be
#'   present after import.
#' @param drop_empty Logical; if `TRUE`, remove columns that are entirely `NA`
#'   or empty strings.
#' @param standardize_names Logical; if `TRUE`, coerce the final `sample_names`
#'   column to character and set rownames to match it.
#' @param verbose Logical; print an informational import summary.
#'
#' @return A data frame suitable for use as `sample_info` in `create_vista()`.
#'
#' @examples
#' data("sample_metadata", package = "VISTA")
#'
#' si <- read_vista_metadata(sample_metadata[seq_len(6), ])
#' head(si$sample_names)
#' @export
read_vista_metadata <- function(x,
                                sample_column = NULL,
                                required_columns = NULL,
                                drop_empty = TRUE,
                                standardize_names = TRUE,
                                verbose = TRUE) {
  sample_info <- .read_tabular_input(x)

  if (!is.null(sample_column)) {
    if (!is.character(sample_column) || length(sample_column) != 1L || !sample_column %in% colnames(sample_info)) {
      cli::cli_abort("{.arg sample_column} must name a column present in {.arg x}.")
    }
    sample_info$sample_names <- as.character(sample_info[[sample_column]])
  }

  if (isTRUE(drop_empty) && ncol(sample_info)) {
    keep_cols <- vapply(sample_info, function(col) {
      vals <- as.character(col)
      !all(is.na(vals) | !nzchar(trimws(vals)))
    }, logical(1))
    sample_info <- sample_info[, keep_cols, drop = FALSE]
  }

  sample_info <- .normalize_sample_info(sample_info)

  if (!is.null(required_columns)) {
    if (!is.character(required_columns)) {
      cli::cli_abort("{.arg required_columns} must be NULL or a character vector.")
    }
    missing_cols <- setdiff(required_columns, colnames(sample_info))
    if (length(missing_cols)) {
      cli::cli_abort("Required metadata column(s) missing: {.val {missing_cols}}")
    }
  }

  if (isTRUE(standardize_names)) {
    sample_info$sample_names <- as.character(sample_info$sample_names)
    rownames(sample_info) <- sample_info$sample_names
  }

  if (isTRUE(verbose)) {
    cli::cli_inform(
      "Imported {.field {nrow(sample_info)}} samples with {.field {ncol(sample_info)}} metadata columns."
    )
  }

  sample_info
}

#' @keywords internal
#' @noRd
.default_metadata_fields <- function(n) {
  paste0("part_", seq_len(n))
}

#' @keywords internal
#' @noRd
.derive_vista_split_fields <- function(sample_names, split = "_", fields = NULL) {
  if (!is.character(split) || length(split) != 1L || !nzchar(split)) {
    cli::cli_abort("{.arg split} must be a single non-empty delimiter.")
  }
  pieces <- strsplit(sample_names, split = split, fixed = TRUE)
  widths <- vapply(pieces, length, integer(1))
  if (length(unique(widths)) != 1L || widths[[1]] <= 1L) {
    cli::cli_abort(
      "Could not split sample names consistently with delimiter {.val {split}}."
    )
  }
  n_fields <- widths[[1]]
  if (is.null(fields)) {
    fields <- .default_metadata_fields(n_fields)
  }
  if (!is.character(fields) || length(fields) != n_fields) {
    cli::cli_abort(
      "{.arg fields} must have length {.val {n_fields}} for the supplied sample names."
    )
  }

  parsed <- do.call(rbind, pieces)
  parsed <- as.data.frame(parsed, stringsAsFactors = FALSE, check.names = FALSE)
  names(parsed) <- fields
  parsed
}

#' @keywords internal
#' @noRd
.derive_vista_regex_fields <- function(sample_names, pattern, fields = NULL) {
  if (!is.character(pattern) || length(pattern) != 1L || !nzchar(pattern)) {
    cli::cli_abort("{.arg pattern} must be a single non-empty regular expression.")
  }
  matches <- regexec(pattern, sample_names, perl = TRUE)
  groups <- regmatches(sample_names, matches)
  matched <- vapply(groups, length, integer(1)) > 0L
  if (!all(matched)) {
    cli::cli_abort(
      "Some sample names did not match {.arg pattern}: {.val {utils::head(sample_names[!matched], 10)}}"
    )
  }

  n_fields <- length(groups[[1]]) - 1L
  if (n_fields <= 0L) {
    cli::cli_abort("{.arg pattern} must define at least one capture group.")
  }
  if (is.null(fields)) {
    fields <- .default_metadata_fields(n_fields)
  }
  if (!is.character(fields) || length(fields) != n_fields) {
    cli::cli_abort(
      "{.arg fields} must have length {.val {n_fields}} for the supplied regex captures."
    )
  }

  parsed <- do.call(rbind, lapply(groups, function(x) x[-1]))
  parsed <- as.data.frame(parsed, stringsAsFactors = FALSE, check.names = FALSE)
  names(parsed) <- fields
  parsed
}

#' @keywords internal
#' @noRd
.guess_vista_split_delimiter <- function(sample_names) {
  candidates <- c("_", "-", ".")
  widths <- lapply(candidates, function(delim) {
    vapply(strsplit(sample_names, split = delim, fixed = TRUE), length, integer(1))
  })
  valid <- vapply(widths, function(x) length(unique(x)) == 1L && x[[1]] > 1L, logical(1))
  if (!any(valid)) {
    return(NULL)
  }
  delim <- candidates[which(valid)[1]]
  list(delimiter = delim, n_fields = widths[[which(valid)[1]]][[1]])
}

#' @keywords internal
#' @noRd
.derive_vista_metadata_fields <- function(sample_names,
                                          parser = c("auto", "split", "regex", "none"),
                                          split = "_",
                                          fields = NULL,
                                          pattern = NULL) {
  parser <- match.arg(parser)
  sample_names <- as.character(sample_names)

  if (identical(parser, "none")) {
    return(NULL)
  }
  if (identical(parser, "split")) {
    return(.derive_vista_split_fields(sample_names, split = split, fields = fields))
  }
  if (identical(parser, "regex")) {
    return(.derive_vista_regex_fields(sample_names, pattern = pattern, fields = fields))
  }

  guessed <- .guess_vista_split_delimiter(sample_names)
  if (is.null(guessed)) {
    return(NULL)
  }
  .derive_vista_split_fields(
    sample_names = sample_names,
    split = guessed$delimiter,
    fields = fields
  )
}

#' Derive starter sample metadata from count sample names
#'
#' `derive_vista_metadata()` creates a starter `sample_info` table from count
#' sample names. It is intended for projects where users have count columns but
#' do not yet have a separate metadata sheet. The derived table can be edited,
#' passed through `read_vista_metadata()`, and then aligned with
#' `match_vista_inputs()`.
#'
#' @param counts Count input accepted by `read_vista_counts()`, or the list
#'   returned by `read_vista_counts()`.
#' @param column_geneid Optional gene identifier column for raw tabular count
#'   inputs. Ignored when `counts` is the list output of `read_vista_counts()`.
#' @param sample_names Optional explicit sample names to derive metadata from.
#'   When supplied, these override names extracted from `counts`.
#' @param parser Metadata parsing mode. `"auto"` tries a simple delimiter-based
#'   split when sample names have a consistent structure. `"split"` uses
#'   `split` explicitly. `"regex"` uses `pattern`. `"none"` returns only the
#'   `sample_names` column.
#' @param split Delimiter used when `parser = "split"` or when `"auto"` chooses
#'   split-based parsing.
#' @param fields Optional field names for parsed metadata columns. When omitted,
#'   VISTA uses `part_1`, `part_2`, etc.
#' @param pattern Regular expression used when `parser = "regex"`. Capture
#'   groups are mapped to `fields` in order.
#' @param sample_column Name of the sample identifier column in the returned
#'   metadata. Default is `"sample_names"`.
#' @param repair_sample_names Strategy passed to `read_vista_counts()` when
#'   sample names are taken from `counts`. One of `"auto"` or `"none"`.
#' @param return_type Return `"data.frame"` (default) or `"template"`. Both
#'   return a data frame; `"template"` adds empty placeholder columns for
#'   `group` and `batch`.
#' @param verbose Logical; print an informational derivation summary.
#'
#' @return A data frame containing `sample_names` plus any parsed metadata
#'   columns.
#'
#' @examples
#' data("count_data", package = "VISTA")
#' data("sample_metadata", package = "VISTA")
#'
#' counts_in <- count_data[seq_len(8), c("gene_id", sample_metadata$sample_names[seq_len(6)]), drop = FALSE]
#' meta <- derive_vista_metadata(
#'   counts_in,
#'   column_geneid = "gene_id",
#'   parser = "regex",
#'   pattern = "SRR(\\d+)",
#'   fields = "run_id"
#' )
#' head(meta)
#' @export
derive_vista_metadata <- function(counts,
                                  column_geneid = NULL,
                                  sample_names = NULL,
                                  parser = c("auto", "split", "regex", "none"),
                                  split = "_",
                                  fields = NULL,
                                  pattern = NULL,
                                  sample_column = "sample_names",
                                  repair_sample_names = c("auto", "none"),
                                  return_type = c("data.frame", "template"),
                                  verbose = TRUE) {
  parser <- match.arg(parser)
  repair_sample_names <- match.arg(repair_sample_names)
  return_type <- match.arg(return_type)

  if (!is.null(sample_names)) {
    if (!is.character(sample_names) || !length(sample_names)) {
      cli::cli_abort("{.arg sample_names} must be NULL or a non-empty character vector.")
    }
    resolved_names <- as.character(sample_names)
    sample_name_map <- data.frame(
      original = resolved_names,
      repaired = resolved_names,
      stringsAsFactors = FALSE
    )
  } else {
    payload <- .coerce_vista_count_payload(counts, column_geneid = column_geneid)
    if (!identical(repair_sample_names, "auto")) {
      payload <- read_vista_counts(
        x = payload$counts,
        format = "matrix",
        gene_id_column = payload$column_geneid,
        repair_sample_names = repair_sample_names,
        verbose = FALSE
      )
    }
    sample_name_map <- payload$sample_name_map
    resolved_names <- if (!is.null(sample_name_map) && nrow(sample_name_map)) {
      as.character(sample_name_map$repaired)
    } else {
      setdiff(colnames(payload$counts), payload$column_geneid)
    }
  }

  df <- data.frame(sample_names = resolved_names, stringsAsFactors = FALSE, check.names = FALSE)
  names(df)[1] <- sample_column

  parsed <- .derive_vista_metadata_fields(
    sample_names = resolved_names,
    parser = parser,
    split = split,
    fields = fields,
    pattern = pattern
  )
  if (!is.null(parsed) && ncol(parsed)) {
    df <- cbind(df, parsed, stringsAsFactors = FALSE)
  }

  if (identical(return_type, "template")) {
    if (!"group" %in% colnames(df)) {
      df$group <- NA_character_
    }
    if (!"batch" %in% colnames(df)) {
      df$batch <- NA_character_
    }
  }

  attr(df, "sample_name_map") <- sample_name_map

  if (isTRUE(verbose)) {
    cli::cli_inform(
      "Derived metadata for {.field {nrow(df)}} samples using parser {.val {parser}}."
    )
  }

  df
}

#' Match count and metadata inputs for VISTA
#'
#' `match_vista_inputs()` aligns standardized counts and sample metadata so they
#' can be passed directly to `create_vista()`. It accepts the raw output from
#' `read_vista_counts()` or a count data frame/matrix plus sample metadata.
#'
#' @param counts Standardized counts from `read_vista_counts()`, or a compatible
#'   count matrix/data frame.
#' @param sample_info Sample metadata from `read_vista_metadata()` or a data
#'   frame coercible to that format.
#' @param column_geneid Optional gene identifier column for raw tabular counts.
#'   Ignored when `counts` is the list output of `read_vista_counts()`.
#' @param sample_column Optional sample identifier column in `sample_info`.
#' @param reorder Logical; if `TRUE` (default), reorder `sample_info` to match
#'   count columns.
#' @param drop_unmatched Logical; if `TRUE`, keep only the intersection of count
#'   samples and metadata samples. Otherwise mismatches raise an error.
#' @param verbose Logical; print an informational alignment summary.
#'
#' @return A list with standardized `counts`, aligned `sample_info`,
#'   `column_geneid`, `sample_names`, `sample_name_map`, `row_data`, and a small `report`.
#'
#' @examples
#' data("count_data", package = "VISTA")
#' data("sample_metadata", package = "VISTA")
#'
#' cnt <- read_vista_counts(
#'   count_data[seq_len(25), ],
#'   format = "matrix",
#'   gene_id_column = "gene_id",
#'   verbose = FALSE
#' )
#' si <- read_vista_metadata(
#'   sample_metadata[sample_metadata$sample_names %in% cnt$sample_names, ],
#'   verbose = FALSE
#' )
#' matched <- match_vista_inputs(cnt, si, verbose = FALSE)
#'
#' matched$column_geneid
#' identical(matched$sample_info$sample_names, colnames(matched$counts)[-1])
#' @export
match_vista_inputs <- function(counts,
                               sample_info,
                               column_geneid = NULL,
                               sample_column = NULL,
                               reorder = TRUE,
                               drop_unmatched = FALSE,
                               verbose = TRUE) {
  count_payload <- .coerce_vista_count_payload(counts, column_geneid = column_geneid)
  sample_info <- read_vista_metadata(
    x = sample_info,
    sample_column = sample_column,
    verbose = FALSE
  )
  sample_info <- .normalize_sample_info(sample_info = sample_info)

  count_samples <- setdiff(colnames(count_payload$counts), count_payload$column_geneid)
  meta_samples <- sample_info$sample_names
  common_samples <- intersect(count_samples, meta_samples)
  missing_in_counts <- setdiff(meta_samples, count_samples)
  missing_in_meta <- setdiff(count_samples, meta_samples)

  if (!isTRUE(drop_unmatched) && (length(missing_in_counts) || length(missing_in_meta))) {
    cli::cli_abort(c(
      "Counts and sample metadata do not align.",
      if (length(missing_in_counts)) "x" = "Present in {.arg sample_info} only: {.val {utils::head(missing_in_counts, 10)}}",
      if (length(missing_in_meta)) "i" = "Present in {.arg counts} only: {.val {utils::head(missing_in_meta, 10)}}"
    ))
  }

  if (isTRUE(drop_unmatched)) {
    keep_cols <- c(count_payload$column_geneid, common_samples)
    count_payload$counts <- count_payload$counts[, keep_cols, drop = FALSE]
    sample_info <- sample_info[sample_info$sample_names %in% common_samples, , drop = FALSE]
    count_samples <- setdiff(colnames(count_payload$counts), count_payload$column_geneid)
    if (!length(count_samples)) {
      cli::cli_abort("No overlapping samples remain after dropping unmatched inputs.")
    }
  } else {
    sample_info <- .normalize_sample_info(
      sample_info = sample_info,
      counts = count_payload$counts,
      column_geneid = count_payload$column_geneid
    )
  }

  if (isTRUE(reorder)) {
    sample_info <- sample_info[match(count_samples, sample_info$sample_names), , drop = FALSE]
    rownames(sample_info) <- sample_info$sample_names
  }

  out <- list(
    counts = count_payload$counts,
    sample_info = sample_info,
    column_geneid = count_payload$column_geneid,
    sample_names = count_samples,
    sample_name_map = count_payload$sample_name_map,
    row_data = count_payload$row_data,
    input_format = count_payload$input_format,
    report = list(
      n_genes = nrow(count_payload$counts),
      n_samples = length(count_samples),
      dropped_from_counts = missing_in_meta,
      dropped_from_sample_info = missing_in_counts
    )
  )

  if (isTRUE(verbose)) {
    cli::cli_inform(
      "Matched {.field {out$report$n_samples}} samples between counts and sample metadata."
    )
  }

  out
}

#' @keywords internal
#' @noRd
.read_tabular_input <- function(x, header = TRUE, comment.char = "", stringsAsFactors = FALSE) {
  if (is.data.frame(x)) {
    return(as.data.frame(x, stringsAsFactors = stringsAsFactors, check.names = FALSE))
  }
  if (is.matrix(x)) {
    return(as.data.frame(x, stringsAsFactors = stringsAsFactors, check.names = FALSE))
  }
  if (!is.character(x) || length(x) != 1L || !file.exists(x)) {
    cli::cli_abort("{.arg x} must be a data frame, matrix, or an existing file path.")
  }

  ext <- tolower(tools::file_ext(x))
  reader <- if (identical(ext, "csv")) utils::read.csv else utils::read.delim
  reader(
    file = x,
    header = header,
    stringsAsFactors = stringsAsFactors,
    check.names = FALSE,
    comment.char = comment.char
  )
}

#' @keywords internal
#' @noRd
.detect_vista_count_format <- function(x, format = c("auto", "matrix", "featurecounts", "star", "htseq", "tximport", "rsem")) {
  format <- match.arg(format)
  if (!identical(format, "auto")) {
    return(format)
  }

  if (is.list(x) && !is.data.frame(x) && !is.matrix(x) && "counts" %in% names(x)) {
    return("tximport")
  }
  if (is.matrix(x)) {
    return("matrix")
  }
  if (is.data.frame(x)) {
    cols <- colnames(x)
    if (all(c("Geneid", "Chr", "Start", "End", "Strand", "Length") %in% cols)) {
      return("featurecounts")
    }
    if ("expected_count" %in% cols && "gene_id" %in% cols) {
      return("rsem")
    }
    return("matrix")
  }
  if (is.character(x) && length(x) >= 1L && all(file.exists(x))) {
    if (length(x) > 1L) {
      first_lines <- readLines(x[[1]], n = 8L, warn = FALSE)
      if (any(grepl("^N_unmapped", first_lines)) || any(grepl("ReadsPerGene.out.tab$", basename(x)))) {
        return("star")
      }
      if (any(grepl("^__", first_lines))) {
        return("htseq")
      }
      if (grepl("genes\\.results$", basename(x[[1]]))) {
        return("rsem")
      }
      return("htseq")
    }

    first_lines <- readLines(x, n = 8L, warn = FALSE)
    if (any(grepl("^N_unmapped", first_lines)) || grepl("ReadsPerGene.out.tab$", basename(x))) {
      return("star")
    }
    if (any(grepl("^__", first_lines))) {
      return("htseq")
    }

    tab <- tryCatch(.read_tabular_input(x, comment.char = "#"), error = function(e) NULL)
    if (!is.null(tab)) {
      cols <- colnames(tab)
      if (all(c("Geneid", "Chr", "Start", "End", "Strand", "Length") %in% cols)) {
        return("featurecounts")
      }
      if ("expected_count" %in% cols && "gene_id" %in% cols) {
        return("rsem")
      }
    }
    return("matrix")
  }

  cli::cli_abort("Could not infer a supported count {.arg format} from {.arg x}.")
}

#' @keywords internal
#' @noRd
.resolve_gene_id_column <- function(df, gene_id_column = NULL) {
  if (!is.null(gene_id_column)) {
    if (!is.character(gene_id_column) || length(gene_id_column) != 1L || !gene_id_column %in% colnames(df)) {
      cli::cli_abort("{.arg gene_id_column} must name a column present in the count input.")
    }
    return(gene_id_column)
  }

  candidates <- intersect(c("gene_id", "Geneid", "gene", "GENEID", "GeneID", "ENSEMBL"), colnames(df))
  if (length(candidates)) {
    return(candidates[[1]])
  }

  rn <- rownames(df)
  if (!is.null(rn) && length(rn) == nrow(df) && all(nzchar(rn)) && !anyDuplicated(rn)) {
    return(NULL)
  }

  non_numeric <- colnames(df)[!vapply(df, is.numeric, logical(1))]
  if (length(non_numeric)) {
    return(non_numeric[[1]])
  }

  cli::cli_abort(
    "Could not determine the gene identifier column. Supply {.arg gene_id_column} or use rownames."
  )
}

#' @keywords internal
#' @noRd
.infer_sample_columns <- function(df, gene_id_column = NULL, sample_columns = NULL, known_annotation = NULL) {
  if (!is.null(sample_columns)) {
    if (!is.character(sample_columns) || !length(sample_columns)) {
      cli::cli_abort("{.arg sample_columns} must be NULL or a character vector.")
    }
    missing_cols <- setdiff(sample_columns, colnames(df))
    if (length(missing_cols)) {
      cli::cli_abort("Sample column(s) not found: {.val {missing_cols}}")
    }
    return(sample_columns)
  }

  drop_cols <- c(gene_id_column, known_annotation)
  candidates <- setdiff(colnames(df), drop_cols)
  if (!length(candidates)) {
    cli::cli_abort("No sample columns were detected in the count input.")
  }
  candidates
}

#' @keywords internal
#' @noRd
.coerce_numeric_counts <- function(df, sample_columns) {
  out <- df
  for (nm in sample_columns) {
    vals <- suppressWarnings(as.numeric(out[[nm]]))
    if (anyNA(vals) && any(!is.na(out[[nm]]))) {
      cli::cli_abort("Sample column {.field {nm}} contains non-numeric values.")
    }
    out[[nm]] <- vals
  }
  out
}

#' @keywords internal
#' @noRd
.counts_df_to_matrix <- function(counts_df, column_geneid = "gene_id") {
  mat <- as.matrix(counts_df[, setdiff(colnames(counts_df), column_geneid), drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- as.character(counts_df[[column_geneid]])
  mat
}

#' @keywords internal
#' @noRd
.import_counts_matrix <- function(x,
                                  gene_id_column = NULL,
                                  sample_columns = NULL,
                                  annotation_columns = NULL) {
  df <- .read_tabular_input(x)
  gene_col <- .resolve_gene_id_column(df, gene_id_column = gene_id_column)

  if (is.null(gene_col)) {
    df$gene_id <- rownames(df)
    gene_col <- "gene_id"
  }

  sample_cols <- .infer_sample_columns(
    df = df,
    gene_id_column = gene_col,
    sample_columns = sample_columns,
    known_annotation = annotation_columns
  )
  df <- .coerce_numeric_counts(df, sample_cols)

  row_cols <- unique(c(gene_col, annotation_columns))
  row_cols <- row_cols[row_cols %in% colnames(df)]

  list(
    counts = df[, c(gene_col, sample_cols), drop = FALSE],
    row_data = df[, row_cols, drop = FALSE]
  )
}

#' @keywords internal
#' @noRd
.import_counts_featurecounts <- function(x,
                                         gene_id_column = NULL,
                                         sample_columns = NULL,
                                         annotation_columns = NULL) {
  df <- if (is.character(x) && length(x) == 1L && file.exists(x)) {
    .read_tabular_input(x, comment.char = "#")
  } else {
    .read_tabular_input(x)
  }
  gene_col <- gene_id_column %||% "Geneid"
  known_annotation <- c("Chr", "Start", "End", "Strand", "Length")
  annotation_cols <- annotation_columns %||% known_annotation
  sample_cols <- .infer_sample_columns(
    df = df,
    gene_id_column = gene_col,
    sample_columns = sample_columns,
    known_annotation = annotation_cols
  )
  df <- .coerce_numeric_counts(df, sample_cols)

  keep_row <- unique(c(gene_col, annotation_cols))
  keep_row <- keep_row[keep_row %in% colnames(df)]

  list(
    counts = df[, c(gene_col, sample_cols), drop = FALSE],
    row_data = df[, keep_row, drop = FALSE]
  )
}

#' @keywords internal
#' @noRd
.parse_star_column <- function(count_column = NULL) {
  count_column <- count_column %||% "unstranded"
  if (is.numeric(count_column) && length(count_column) == 1L) {
    idx <- as.integer(count_column)
    if (!idx %in% c(2L, 3L, 4L)) {
      cli::cli_abort("STAR {.arg count_column} index must be 2, 3, or 4.")
    }
    return(idx)
  }
  idx <- match.arg(as.character(count_column), c("unstranded", "stranded_first", "stranded_second"))
  c(unstranded = 2L, stranded_first = 3L, stranded_second = 4L)[[idx]]
}

#' @keywords internal
#' @noRd
.infer_file_sample_names <- function(paths) {
  repaired <- .repair_vista_sample_names(paths)
  out <- repaired$repaired
  out[!nzchar(out)] <- vapply(paths[!nzchar(out)], function(path) {
    tools::file_path_sans_ext(basename(path))
  }, character(1))
  make.unique(out)
}

#' Bind per-sample count files by gene identifier
#'
#' Shared by the STAR, HTSeq and RSEM importers. Each of those previously used
#' `vapply(FUN.VALUE = numeric(length(genes)))`, which enforces only that every
#' file has the same *number* of rows -- never that row `i` describes the same
#' gene. Files with equal row counts but different ordering were silently pasted
#' onto the first file's gene labels.
#'
#' Rows are taken from the first file and every other file is indexed by
#' `match()` against it. A file that does not cover the reference gene set is an
#' error rather than a silent zero-fill, because a silent zero is the same class
#' of bug this function exists to prevent.
#'
#' @param tabs List of data frames, one per sample.
#' @param id_col Column name or index holding gene identifiers.
#' @param value_col Column name or index holding counts.
#' @param sample_names Character vector of sample names, parallel to `tabs`.
#' @param files Optional file paths used in diagnostics; falls back to `sample_names`.
#' @param format Short label for the input format, used in messages.
#'
#' @return A `data.frame` with `gene_id` plus one numeric column per sample.
#' @keywords internal
#' @noRd
.bind_count_files <- function(tabs, id_col, value_col, sample_names,
                              files = NULL, format = "count") {
  labels <- if (!is.null(files) && length(files) == length(tabs)) {
    basename(as.character(files))
  } else {
    as.character(sample_names)
  }

  has_col <- function(df, col) {
    if (is.numeric(col)) ncol(df) >= col else col %in% colnames(df)
  }
  for (i in seq_along(tabs)) {
    if (!has_col(tabs[[i]], id_col) || !has_col(tabs[[i]], value_col)) {
      cli::cli_abort(
        "{format} file {.file {labels[[i]]}} is missing the expected identifier or count column."
      )
    }
  }

  ref_ids <- as.character(tabs[[1]][[id_col]])
  if (!length(ref_ids)) {
    cli::cli_abort("{format} file {.file {labels[[1]]}} contains no rows.")
  }

  as_counts <- function(df, label) {
    vals <- suppressWarnings(as.numeric(df[[value_col]]))
    if (anyNA(vals) && any(!is.na(df[[value_col]]))) {
      cli::cli_abort("{format} file {.file {label}} contains non-numeric counts.")
    }
    vals
  }

  ids_list <- lapply(tabs, function(df) as.character(df[[id_col]]))
  identical_ids <- all(vapply(ids_list, identical, logical(1), ref_ids))

  if (identical_ids) {
    # Fast path: every file already lists the same genes in the same order, so
    # positional binding is correct even when identifiers repeat.
    counts <- vapply(
      seq_along(tabs),
      function(i) as_counts(tabs[[i]], labels[[i]]),
      numeric(length(ref_ids))
    )
  } else {
    if (anyDuplicated(ref_ids)) {
      cli::cli_abort(c(
        "{format} file {.file {labels[[1]]}} has duplicated gene identifiers and the files are not in a common order.",
        "i" = "Deduplicate the identifiers, or supply the files in a matching row order."
      ))
    }
    counts <- matrix(NA_real_, nrow = length(ref_ids), ncol = length(tabs))
    for (i in seq_along(tabs)) {
      ids <- ids_list[[i]]
      if (anyDuplicated(ids)) {
        cli::cli_abort(
          "{format} file {.file {labels[[i]]}} has duplicated gene identifiers, so it cannot be matched by identifier."
        )
      }
      idx <- match(ref_ids, ids)
      if (anyNA(idx)) {
        missing_ids <- ref_ids[is.na(idx)]
        cli::cli_abort(c(
          "Gene identifiers in {format} file {.file {labels[[i]]}} do not cover the reference gene set.",
          "x" = "Missing ({length(missing_ids)} of {length(ref_ids)}): {.val {utils::head(missing_ids, 5)}}",
          "i" = "The reference set is taken from {.file {labels[[1]]}}."
        ))
      }
      extra <- length(ids) - length(ref_ids)
      if (extra > 0L) {
        cli::cli_inform(
          "Dropped {extra} identifier{?s} present in {.file {labels[[i]]}} but absent from {.file {labels[[1]]}}."
        )
      }
      counts[, i] <- as_counts(tabs[[i]], labels[[i]])[idx]
    }
  }

  counts <- as.matrix(counts)
  colnames(counts) <- sample_names
  data.frame(
    gene_id = ref_ids,
    counts,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

#' @keywords internal
#' @noRd
.import_counts_star <- function(x,
                                sample_names = NULL,
                                count_column = NULL,
                                drop_technical = TRUE) {
  star_col <- .parse_star_column(count_column)

  if (is.character(x) && length(x) >= 1L && all(file.exists(x))) {
    files <- x
    sample_names <- sample_names %||% .infer_file_sample_names(files)
    if (length(sample_names) != length(files)) {
      cli::cli_abort("{.arg sample_names} must have the same length as the STAR file vector.")
    }

    mats <- lapply(files, function(path) {
      utils::read.delim(path, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
    })
    for (i in seq_along(mats)) {
      if (ncol(mats[[i]]) < 4L) {
        cli::cli_abort("STAR gene count files must have four columns.")
      }
    }
    out <- .bind_count_files(
      tabs = mats,
      id_col = 1L,
      value_col = star_col,
      sample_names = sample_names,
      files = files,
      format = "STAR"
    )
  } else {
    df <- .read_tabular_input(x)
    if (ncol(df) < 4L) {
      cli::cli_abort("STAR tabular input must contain gene identifiers plus three count columns.")
    }
    colnames(df)[1] <- "gene_id"
    if (ncol(df) == 4L) {
      colnames(df)[2:4] <- c("unstranded", "stranded_first", "stranded_second")
    }
    star_name <- names(c(unstranded = 2L, stranded_first = 3L, stranded_second = 4L))[match(star_col, c(2L, 3L, 4L))]
    out <- data.frame(
      gene_id = df$gene_id,
      sample1 = suppressWarnings(as.numeric(df[[star_name]])),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  if (isTRUE(drop_technical)) {
    keep <- !out$gene_id %in% c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
    out <- out[keep, , drop = FALSE]
  }

  list(
    counts = out,
    row_data = out["gene_id"]
  )
}

#' @keywords internal
#' @noRd
.import_counts_htseq <- function(x,
                                 sample_names = NULL,
                                 drop_technical = TRUE) {
  if (is.character(x) && length(x) >= 1L && all(file.exists(x))) {
    files <- x
    sample_names <- sample_names %||% .infer_file_sample_names(files)
    if (length(sample_names) != length(files)) {
      cli::cli_abort("{.arg sample_names} must have the same length as the HTSeq file vector.")
    }

    mats <- lapply(files, function(path) {
      utils::read.delim(path, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
    })
    out <- .bind_count_files(
      tabs = mats,
      id_col = 1L,
      value_col = 2L,
      sample_names = sample_names,
      files = files,
      format = "HTSeq"
    )
  } else {
    df <- .read_tabular_input(x)
    gene_col <- .resolve_gene_id_column(df)
    sample_cols <- .infer_sample_columns(df, gene_id_column = gene_col)
    df <- .coerce_numeric_counts(df, sample_cols)
    out <- df[, c(gene_col, sample_cols), drop = FALSE]
    names(out)[1] <- "gene_id"
  }

  if (isTRUE(drop_technical)) {
    out <- out[!startsWith(as.character(out$gene_id), "__"), , drop = FALSE]
  }

  list(
    counts = out,
    row_data = out["gene_id"]
  )
}

#' @keywords internal
#' @noRd
.import_counts_tximport <- function(x,
                                    tx2gene = NULL,
                                    counts_from = c("counts", "abundance", "length")) {
  counts_from <- match.arg(counts_from)
  if (!is.list(x) || is.data.frame(x) || is.matrix(x)) {
    cli::cli_abort("tximport input must be a list-like object.")
  }
  if (!counts_from %in% names(x)) {
    cli::cli_abort("Requested {.arg counts_from} matrix {.val {counts_from}} is not present in the tximport input.")
  }

  mat <- x[[counts_from]]
  if (!is.matrix(mat) && !is.data.frame(mat)) {
    cli::cli_abort("The selected tximport matrix must be a matrix or data frame.")
  }
  mat <- as.matrix(mat)

  if (!is.null(tx2gene)) {
    tx2gene <- as.data.frame(tx2gene, stringsAsFactors = FALSE)
    if (ncol(tx2gene) < 2L) {
      cli::cli_abort("{.arg tx2gene} must contain at least two columns: transcript ID and gene ID.")
    }
    tx_map <- stats::setNames(as.character(tx2gene[[2]]), as.character(tx2gene[[1]]))
    gene_ids <- tx_map[rownames(mat)]
    keep <- !is.na(gene_ids) & nzchar(gene_ids)
    mat <- mat[keep, , drop = FALSE]
    gene_ids <- gene_ids[keep]
    mat <- rowsum(mat, group = gene_ids, reorder = FALSE)
  }

  out <- data.frame(
    gene_id = rownames(mat),
    as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  list(
    counts = out,
    row_data = out["gene_id"]
  )
}

#' @keywords internal
#' @noRd
.import_counts_rsem <- function(x,
                                sample_names = NULL,
                                count_column = NULL,
                                annotation_columns = NULL) {
  count_column <- count_column %||% "expected_count"

  if (is.character(x) && length(x) >= 1L && all(file.exists(x))) {
    files <- x
    sample_names <- sample_names %||% .infer_file_sample_names(files)
    if (length(sample_names) != length(files)) {
      cli::cli_abort("{.arg sample_names} must have the same length as the RSEM file vector.")
    }

    tabs <- lapply(files, .read_tabular_input)
    first <- tabs[[1]]
    if (!all(c("gene_id", count_column) %in% colnames(first))) {
      cli::cli_abort("RSEM input must contain {.field gene_id} and {.field {count_column}} columns.")
    }
    for (df in tabs) {
      if (!all(c("gene_id", count_column) %in% colnames(df))) {
        cli::cli_abort("Every RSEM file must contain {.field gene_id} and {.field {count_column}}.")
      }
    }
    out_counts <- .bind_count_files(
      tabs = tabs,
      id_col = "gene_id",
      value_col = count_column,
      sample_names = sample_names,
      files = files,
      format = "RSEM"
    )
    row_cols <- unique(c("gene_id", annotation_columns %||% setdiff(colnames(first), count_column)))
    row_cols <- row_cols[row_cols %in% colnames(first)]
    row_data <- first[, row_cols, drop = FALSE]
  } else {
    df <- .read_tabular_input(x)
    if (!all(c("gene_id", count_column) %in% colnames(df))) {
      cli::cli_abort("RSEM tabular input must contain {.field gene_id} and {.field {count_column}}.")
    }
    sample_cols <- sample_names %||% NULL
    if (!is.null(sample_cols) && !all(sample_cols %in% colnames(df))) {
      cli::cli_abort("{.arg sample_names} must refer to columns present in the RSEM table when importing a merged table.")
    }
    if (is.null(sample_cols)) {
      sample_cols <- count_column
    }
    if (identical(sample_cols, count_column)) {
      out_counts <- df[, c("gene_id", count_column), drop = FALSE]
      names(out_counts)[2] <- "sample1"
    } else {
      out_counts <- df[, c("gene_id", sample_cols), drop = FALSE]
    }
    out_counts <- .coerce_numeric_counts(out_counts, setdiff(names(out_counts), "gene_id"))
    row_cols <- unique(c("gene_id", annotation_columns %||% setdiff(colnames(df), count_column)))
    row_cols <- row_cols[row_cols %in% colnames(df)]
    row_data <- df[, row_cols, drop = FALSE]
  }

  list(
    counts = out_counts,
    row_data = row_data
  )
}

#' @keywords internal
#' @noRd
.standardize_counts_import <- function(parsed,
                                       input_format,
                                       make_unique_ids = FALSE,
                                       repair_sample_names = c("auto", "none")) {
  repair_sample_names <- match.arg(repair_sample_names)
  counts <- as.data.frame(parsed$counts, stringsAsFactors = FALSE, check.names = FALSE)
  row_data <- as.data.frame(parsed$row_data, stringsAsFactors = FALSE, check.names = FALSE)

  gene_col <- names(counts)[1]
  gene_ids <- as.character(counts[[gene_col]])
  if (anyNA(gene_ids) || any(!nzchar(gene_ids))) {
    cli::cli_abort("Imported gene IDs must be non-empty and non-missing.")
  }
  if (anyDuplicated(gene_ids)) {
    if (isTRUE(make_unique_ids)) {
      gene_ids <- make.unique(gene_ids)
    } else {
      dupes <- unique(gene_ids[duplicated(gene_ids)])
      cli::cli_abort("Duplicated gene IDs detected: {.val {utils::head(dupes, 10)}}")
    }
  }

  names(counts)[1] <- "gene_id"
  counts$gene_id <- gene_ids
  sample_names <- setdiff(colnames(counts), "gene_id")
  if (!length(sample_names)) {
    cli::cli_abort("No sample columns were detected after standardizing the count input.")
  }
  sample_name_map <- data.frame(
    original = sample_names,
    repaired = sample_names,
    stringsAsFactors = FALSE
  )
  if (identical(repair_sample_names, "auto")) {
    repaired <- .repair_vista_sample_names(sample_names)
    sample_name_map$repaired <- repaired$repaired
    if (isTRUE(repaired$applied)) {
      names(counts)[match(sample_names, names(counts))] <- repaired$repaired
      sample_names <- unname(repaired$repaired)
    } else if (isTRUE(repaired$changed)) {
      cli::cli_warn(repaired$message)
    }
  }
  sample_names <- unname(sample_names)
  if (anyDuplicated(sample_names)) {
    cli::cli_abort("Sample column names must be unique.")
  }
  counts <- .coerce_numeric_counts(counts, sample_names)

  if (!"gene_id" %in% colnames(row_data)) {
    row_data$gene_id <- gene_ids
  } else {
    row_data$gene_id <- gene_ids
  }
  rownames(row_data) <- gene_ids

  list(
    counts = counts[, c("gene_id", sample_names), drop = FALSE],
    row_data = row_data,
    column_geneid = "gene_id",
    sample_names = sample_names,
    sample_name_map = sample_name_map,
    input_format = input_format,
    report = list(
      n_genes = length(gene_ids),
      n_samples = length(sample_names),
      repaired_sample_names = sum(sample_name_map$original != sample_name_map$repaired)
    )
  )
}

#' @keywords internal
#' @noRd
.coerce_vista_count_payload <- function(counts, column_geneid = NULL) {
  if (is.list(counts) && all(c("counts", "column_geneid") %in% names(counts))) {
    return(list(
      counts = as.data.frame(counts$counts, stringsAsFactors = FALSE, check.names = FALSE),
      column_geneid = counts$column_geneid,
      row_data = counts$row_data %||% data.frame(gene_id = counts$counts[[counts$column_geneid]], stringsAsFactors = FALSE),
      sample_name_map = counts$sample_name_map %||% data.frame(
        original = setdiff(colnames(counts$counts), counts$column_geneid),
        repaired = setdiff(colnames(counts$counts), counts$column_geneid),
        stringsAsFactors = FALSE
      ),
      input_format = counts$input_format %||% "matrix"
    ))
  }

  standardized <- read_vista_counts(
    x = counts,
    format = "matrix",
    gene_id_column = column_geneid,
    verbose = FALSE
  )
  list(
    counts = standardized$counts,
    column_geneid = standardized$column_geneid,
    row_data = standardized$row_data,
    sample_name_map = standardized$sample_name_map,
    input_format = standardized$input_format
  )
}

#' @title Repair file-derived sample names for VISTA imports
#'
#' @description
#' Internal helper used by `read_vista_counts()` to derive stable sample names
#' from count column labels that look like file paths or pipeline-generated
#' file names.
#'
#' Repair rules are applied conservatively and in order:
#' \itemize{
#'   \item strip directory paths to the basename
#'   \item for generic quantification files such as `quant.sf` or
#'   `abundance.tsv`, use the parent directory name
#'   \item strip common RNA-seq output suffixes such as
#'   `Aligned.sortedByCoord.out.bam`, `ReadsPerGene.out.tab`,
#'   `.genes.results`, `.isoforms.results`, `.bam`, `.fastq.gz`
#'   \item strip common lane/read suffixes such as `_S1_L001_R1_001`,
#'   `_L001_R2_001`, `_R1`, or `_R2`
#'   \item normalize whitespace and repeated punctuation
#' }
#'
#' Repaired names are only applied automatically when they are non-empty and
#' unique; otherwise the original names are retained and a warning is issued.
#'
#' @param sample_names Character vector of sample column names.
#'
#' @return A list with `original`, `repaired`, `changed`, `applied`, and
#'   `message`.
#' @keywords internal
#' @noRd
.repair_vista_sample_names <- function(sample_names) {
  if (!is.character(sample_names)) {
    cli::cli_abort("{.arg sample_names} must be a character vector.")
  }

  original <- sample_names
  repaired <- sample_names

  derive_one <- function(x) {
    if (is.na(x) || !nzchar(x)) {
      return(x)
    }

    is_path <- grepl("[/\\\\]", x)
    base <- if (is_path) basename(x) else x
    parent <- if (is_path) basename(dirname(x)) else ""
    candidate <- base

    if (tolower(candidate) %in% c("quant.sf", "abundance.tsv", "abundance.h5", "readspergene.out.tab")) {
      if (nzchar(parent) && !identical(parent, ".") && !identical(parent, "/")) {
        candidate <- parent
      }
    }

    suffix_patterns <- c(
      "(_star_align)?Aligned\\.sortedByCoord\\.out\\.(bam|sam|cram)$",
      "_ReadsPerGene\\.out\\.tab$",
      "ReadsPerGene\\.out\\.tab$",
      "\\.genes\\.results$",
      "\\.isoforms\\.results$",
      "\\.(bam|sam|cram)$",
      "\\.(fastq|fq)(\\.gz)?$",
      "\\.(csv|tsv|txt|tab)$"
    )
    for (pat in suffix_patterns) {
      candidate <- sub(pat, "", candidate, perl = TRUE)
    }

    lane_patterns <- c(
      "_S[0-9]+_L[0-9]{3}_R[12]_001$",
      "_L[0-9]{3}_R[12]_001$",
      "_R[12]_001$",
      "_R[12]$"
    )
    for (pat in lane_patterns) {
      candidate <- sub(pat, "", candidate, perl = TRUE)
    }

    candidate <- gsub("%20", "_", candidate, fixed = TRUE)
    candidate <- gsub("[[:space:]]+", "_", candidate, perl = TRUE)
    candidate <- gsub("[()\\[\\]{}]+", "", candidate, perl = TRUE)
    candidate <- gsub("[^A-Za-z0-9._-]+", "_", candidate, perl = TRUE)
    candidate <- gsub("_{2,}", "_", candidate, perl = TRUE)
    candidate <- gsub("^[._-]+|[._-]+$", "", candidate, perl = TRUE)

    if (!nzchar(candidate)) {
      candidate <- if (nzchar(parent) && !identical(parent, ".") && !identical(parent, "/")) {
        parent
      } else {
        tools::file_path_sans_ext(base)
      }
      candidate <- gsub("[^A-Za-z0-9._-]+", "_", candidate, perl = TRUE)
      candidate <- gsub("_{2,}", "_", candidate, perl = TRUE)
      candidate <- gsub("^[._-]+|[._-]+$", "", candidate, perl = TRUE)
    }

    candidate
  }

  repaired <- vapply(sample_names, derive_one, character(1))
  changed <- !is.na(original) & !is.na(repaired) & original != repaired
  can_apply <- all(!is.na(repaired)) && all(nzchar(repaired)) && !anyDuplicated(repaired)

  message <- NULL
  if (any(changed) && !can_apply) {
    message <- paste(
      "Automatic sample-name repair was skipped because repaired names were not unique or became empty.",
      "Provide explicit sample names if you want to override the original columns."
    )
    repaired <- original
  }

  list(
    original = original,
    repaired = repaired,
    changed = any(changed),
    applied = any(changed) && can_apply,
    message = message
  )
}
