#' Create a VISTA Object with Internal DE Analysis
#'
#' This wrapper performs differential expression (DE) analysis (DESeq2, edgeR, limma, or both)
#' and returns a fully initialized \code{VISTA} object. The object stores
#' expression matrices and annotations in the SummarizedExperiment core, while
#' all DE outputs and configuration live in \code{metadata(vista)}:
#' \itemize{
#'   \item \code{$de_results}: named \code{SimpleList} of per-contrast DE tables
#'   \item \code{$de_summary}: named \code{SimpleList} of summary tables
#'   \item \code{$de_cutoffs}: list of thresholds/method options
#'   \item \code{$group}: list with \code{column}, \code{palette}, \code{colors}
#' }
#'
#' @param counts Raw counts (matrix/data.frame) with a gene-id column and sample columns.
#' @param sample_info Data frame with sample metadata. Must contain \code{sample_names}
#'   (or have rownames equal to sample columns in \code{counts}) and the \code{group_column}.
#' @param column_geneid Column name in \code{counts} that contains gene identifiers.
#' @param group_column Column in \code{sample_info} used to group samples.
#' @param group_numerator Character vector of numerator groups for pairwise comparisons.
#' @param group_denominator Character vector of denominator groups (same length/order as numerator).
#' @param method \code{"deseq2"}, \code{"edger"}, \code{"limma"}, or \code{"both"}.
#' @param min_counts Minimum total counts per gene to retain (default: 10).
#' @param min_replicates Minimum samples per group meeting \code{min_counts} (default: 1).
#' @param log2fc_cutoff Absolute LFC threshold for DEG calling (default: 1).
#' @param pval_cutoff p-value (raw or adjusted) threshold (default: 0.05).
#' @param p_value_type Which p-value column to use (\code{"padj"} or \code{"pvalue"}). Default: \code{"padj"}.
#' @param covariates Optional character vector of additional sample_info columns to adjust for.
#'   These are included as additive terms in the DE design.
#' @param design_formula Optional model formula (or formula string) overriding automatic
#'   construction from \code{group_column} + \code{covariates}. Must include \code{group_column}.
#' @param consensus_mode When \code{method = "both"}, how to define consensus calls:
#'   \code{"intersection"} (both methods significant in same direction) or
#'   \code{"union"} (either method significant; discordant directions excluded).
#' @param consensus_log2fc When \code{method = "both"}, how to populate consensus
#'   \code{log2fc}: \code{"mean"}, \code{"deseq2"}, or \code{"edger"}.
#' @param result_source Active DE source used in \code{metadata(v)$de_results}. For
#'   \code{method = "both"}, one of \code{"consensus"}, \code{"deseq2"}, \code{"edger"}.
#'   For single-method runs, this must match \code{method}.
#' @param group_palette Qualitative palette name for \code{colorspace::qualitative_hcl()}.
#'   One of \code{c("Pastel 1","Dark 2","Dark 3","Set 2","Set 3","Warm","Cold","Harmonic","Dynamic")}.
#'   Default: \code{"Dark 2"}.
#' @param comparison_palette Qualitative palette name used to assign colors per
#'   comparison (stored in \code{metadata(v)$comparison$colors}). Defaults to \code{"Dark 3"}.
#' @param validate Logical; if \code{TRUE} (default), run full \code{validate_vista()}
#'   checks before returning the object.
#'
#' @return A \code{VISTA} object:
#'   \itemize{
#'     \item \strong{assays(v)}: \code{norm_counts} (matrix)
#'     \item \strong{colData(v)}: \code{sample_info} (DataFrame)
#'     \item \strong{rowData(v)}: \code{row_data} (DataFrame)
#'     \item \strong{metadata(v)}: \code{de_results}, \code{de_summary}, \code{de_cutoffs}, \code{group}, \code{comparison}, \code{provenance}
#'   }
#'
#' @details
#' Contrast names follow \code{"numerator_VS_denominator"}.
#' Each DE table must have rownames identical to the final \code{norm_counts} rownames.
#' When \code{method = "both"}, method-specific and consensus DE tables are stored in
#' \code{metadata(v)$de_results_by_method} and \code{metadata(v)$de_summary_by_method},
#' and the active source is tracked in \code{metadata(v)$de_active_source}.
#' @seealso \link{as_vista}, \link{VISTA-class}, \link[colorspace]{qualitative_hcl}
#'
#' @examples
#' # Load example data
#' data("count_data", package = "VISTA")
#' data("sample_metadata", package = "VISTA")
#'
#' # Create VISTA object with DESeq2 (default method)
#' vista <- create_vista(
#'   counts = count_data[1:100, ],
#'   sample_info = sample_metadata[1:6, ],
#'   column_geneid = "gene_id",
#'   group_column = "cond_long",
#'   group_numerator = "treatment1",
#'   group_denominator = "control",
#'   log2fc_cutoff = 0.6,
#'   pval_cutoff = 0.05
#' )
#'
#' # Examine the VISTA object
#' vista
#'
#' # Access comparisons
#' names(comparisons(vista))
#'
#' # View DEG summary
#' deg_summary(vista)
#'
#' # View cutoffs used
#' cutoffs(vista)
#'
#' # Multiple comparisons example
#' \donttest{
#' vista_multi <- create_vista(
#'   counts = count_data,
#'   sample_info = sample_metadata,
#'   column_geneid = "gene_id",
#'   group_column = "cell",
#'   group_numerator = c("N052611", "N080611"),
#'   group_denominator = c("N61311", "N61311"),
#'   method = "edger",
#'   log2fc_cutoff = 1.0,
#'   pval_cutoff = 0.01
#' )
#' }
#'
#' @importFrom dplyr %>% mutate
#' @importFrom tibble as_tibble column_to_rownames rownames_to_column tibble
#' @importFrom DESeq2 DESeq DESeqDataSetFromMatrix results counts
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment colData rowData
#' @export
create_vista <- function(counts,
                         sample_info,
                         column_geneid,
                         group_column,
                         group_numerator,
                         group_denominator,
                         method = c("deseq2", "edger", "limma", "both"),
                         min_counts = 10,
                         min_replicates = 1,
                         log2fc_cutoff = 1,
                         pval_cutoff = 0.05,
                         p_value_type = "padj",
                         covariates = NULL,
                         design_formula = NULL,
                         consensus_mode = c("intersection", "union"),
                         consensus_log2fc = c("mean", "deseq2", "edger"),
                         result_source = NULL,
                         group_palette = "Dark 2",
                         comparison_palette = "Dark 3",
                         validate = TRUE) {

  method <- match.arg(method)
  p_value_type <- match.arg(p_value_type, c("padj", "pvalue"))
  consensus_mode <- match.arg(consensus_mode)
  consensus_log2fc <- match.arg(consensus_log2fc)
  sample_info <- .normalize_sample_info(
    sample_info = sample_info,
    counts = counts,
    column_geneid = column_geneid
  )
  design_formula <- .coerce_design_formula(design_formula)
  covariates <- .validate_covariates(sample_info, group_column, covariates)
  if (is.factor(group_numerator)) group_numerator <- as.character(group_numerator)
  if (is.factor(group_denominator)) group_denominator <- as.character(group_denominator)

  # Validate group vectors
  if (!is.character(group_numerator) || length(group_numerator) == 0) {
    cli::cli_abort("{.arg group_numerator} must be a non-empty character vector.")
  }
  if (!is.character(group_denominator) || length(group_denominator) == 0) {
    cli::cli_abort("{.arg group_denominator} must be a non-empty character vector.")
  }
  if (length(group_numerator) != length(group_denominator)) {
    cli::cli_abort(
      "{.arg group_numerator} and {.arg group_denominator} must have equal length.
      Got {length(group_numerator)} and {length(group_denominator)}."
    )
  }

  # Validate groups exist in sample_info
  if (!group_column %in% colnames(sample_info)) {
    cli::cli_abort("{.arg group_column} '{group_column}' not found in {.arg sample_info}.")
  }

  unique_groups <- unique(sample_info[[group_column]])
  missing_num <- setdiff(group_numerator, unique_groups)
  missing_den <- setdiff(group_denominator, unique_groups)

  if (length(missing_num) > 0) {
    cli::cli_abort(
      "{.arg group_numerator} contains groups not in data: {.val {missing_num}}"
    )
  }
  if (length(missing_den) > 0) {
    cli::cli_abort(
      "{.arg group_denominator} contains groups not in data: {.val {missing_den}}"
    )
  }

  # Warn about low replication
  groups_in_use <- unique(c(group_numerator, group_denominator))
  groups_reps <- table(sample_info[[group_column]][sample_info[[group_column]] %in% groups_in_use])
  if (any(groups_reps < 2)) {
    low_rep_groups <- names(groups_reps)[groups_reps < 2]
    cli::cli_warn(
      "Group(s) {.val {low_rep_groups}} have fewer than 2 replicates.
      Statistical testing may be unreliable."
    )
  }

  if (method == "both") {
    allowed_sources <- c("consensus", "deseq2", "edger")
    if (is.null(result_source)) {
      result_source <- "consensus"
    } else {
      result_source <- match.arg(result_source, allowed_sources)
    }
  } else {
    if (is.null(result_source)) {
      result_source <- method
    }
    if (!identical(result_source, method)) {
      cli::cli_abort("{.arg result_source} must match {.arg method} when {.arg method} is not {.val both}.")
    }
  }

  run_backend <- function(backend) {
    switch(
      backend,
      deseq2 = run_deseq_analysis(
        counts = counts,
        sample_info = sample_info,
        column_geneid = column_geneid,
        group_column = group_column,
        group_numerator = group_numerator,
        group_denominator = group_denominator,
        covariates = covariates,
        design_formula = design_formula,
        min_counts = min_counts,
        min_replicates = min_replicates,
        log2fc_cutoff = log2fc_cutoff,
        pval_cutoff = pval_cutoff,
        p_value_type = p_value_type
      ),
      edger = run_edger_analysis(
        counts = counts,
        sample_info = sample_info,
        column_geneid = column_geneid,
        group_column = group_column,
        group_numerator = group_numerator,
        group_denominator = group_denominator,
        covariates = covariates,
        design_formula = design_formula,
        min_counts = min_counts,
        min_replicates = min_replicates,
        log2fc_cutoff = log2fc_cutoff,
        pval_cutoff = pval_cutoff,
        p_value_type = ifelse(p_value_type == "padj", "FDR", "PValue")
      ),
      limma = run_limma_analysis(
        counts = counts,
        sample_info = sample_info,
        column_geneid = column_geneid,
        group_column = group_column,
        group_numerator = group_numerator,
        group_denominator = group_denominator,
        covariates = covariates,
        design_formula = design_formula,
        min_counts = min_counts,
        min_replicates = min_replicates,
        log2fc_cutoff = log2fc_cutoff,
        pval_cutoff = pval_cutoff,
        p_value_type = ifelse(p_value_type == "padj", "FDR", "PValue")
      )
    )
  }

  if (method == "both") {
    de_results_deseq <- run_backend("deseq2")
    de_results_edger <- run_backend("edger")
    norm_source <- if (result_source == "edger") "edger" else "deseq2"
    de_results_base <- if (norm_source == "edger") de_results_edger else de_results_deseq
  } else {
    de_results_base <- run_backend(method)
    de_results_deseq <- NULL
    de_results_edger <- NULL
  }

  ref_rn <- rownames(de_results_base$norm_counts)
  .std_comparisons <- function(comp_list) {
    purrr::map(comp_list, ~{
      df <- .tidy_de_results(.x, rowname_col = "gene_id")
      .align_de_to_counts(df, ref_rn)
    })
  }

  if (method == "both") {
    comps_deseq <- .std_comparisons(de_results_deseq$comparisons)
    comps_edger <- .std_comparisons(de_results_edger$comparisons)

    p_type_consensus <- if (p_value_type == "padj") "padj" else "pvalue"
    consensus_results <- .build_consensus_results(
      deseq_comparisons = comps_deseq,
      edger_comparisons = comps_edger,
      ref_rn = ref_rn,
      log2fc_cutoff = log2fc_cutoff,
      pval_cutoff = pval_cutoff,
      p_value_type = p_type_consensus,
      consensus_mode = consensus_mode,
      consensus_log2fc = consensus_log2fc
    )

    comparisons_by_source <- list(
      deseq2 = comps_deseq,
      edger = comps_edger,
      consensus = consensus_results$comparisons
    )
    deg_summary_by_source <- list(
      deseq2 = purrr::map(comps_deseq, ~ dplyr::count(.x, regulation)),
      edger = purrr::map(comps_edger, ~ dplyr::count(.x, regulation)),
      consensus = consensus_results$deg_summary
    )
  } else {
    standardized <- .std_comparisons(de_results_base$comparisons)
    comparisons_by_source <- list()
    comparisons_by_source[[method]] <- standardized
    deg_summary_by_source <- list()
    deg_summary_by_source[[method]] <- purrr::map(standardized, ~ dplyr::count(.x, regulation))
  }

  active_comparisons <- comparisons_by_source[[result_source]]
  active_deg_summary <- deg_summary_by_source[[result_source]]

  # --- Align core pieces before calling .vista() ---

  # 1) sample_info: ensure rownames and order == colnames(norm_counts)
  si <- as.data.frame(de_results_base$sample_info, stringsAsFactors = FALSE)
  if (is.null(rownames(si)) && "sample_names" %in% colnames(si)) {
    rownames(si) <- as.character(si$sample_names)
  }
  if (!all(colnames(de_results_base$norm_counts) %in% rownames(si))) {
    cli::cli_abort("Some sample columns in normalized counts are not present in sample metadata rownames.")
  }
  si <- si[colnames(de_results_base$norm_counts), , drop = FALSE]
  # Keep explicit sample_names for downstream examples/tests and user ergonomics.
  si$sample_names <- rownames(si)

  # 2) row_data: promote gene_id to rownames, align to counts rows, pad missing
  rd <- as.data.frame(de_results_base$row_data, stringsAsFactors = FALSE)
  if (is.null(rownames(rd)) && "gene_id" %in% colnames(rd)) {
    rownames(rd) <- as.character(rd$gene_id)
  }
  if (is.null(rownames(rd))) {
    if (nrow(rd) == length(ref_rn)) {
      rownames(rd) <- ref_rn
    } else {
      cli::cli_abort("Could not resolve row names for {.arg row_data}.")
    }
  }

  missing_ids <- setdiff(ref_rn, rownames(rd))
  if (length(missing_ids)) {
    na_block <- matrix(NA_real_, nrow = length(missing_ids), ncol = ncol(rd),
                       dimnames = list(missing_ids, colnames(rd)))
    rd <- rbind(rd, na_block)
  }
  rd <- rd[ref_rn, , drop = FALSE]

  # --- Now it's safe to construct VISTA ---
  v <- .vista(
    norm_counts   = de_results_base$norm_counts,
    sample_info   = si,
    row_data      = rd,
    comparisons   = active_comparisons,
    deg_summary   = active_deg_summary,
    cutoffs       = list(
      log2fc       = log2fc_cutoff,
      pval         = pval_cutoff,
      p_value_type = p_value_type,
      method       = method,
      min_counts   = min_counts,
      min_replicates = min_replicates,
      covariates = covariates,
      design_formula = if (is.null(design_formula)) NULL else paste(deparse(design_formula), collapse = ""),
      consensus_mode = if (method == "both") consensus_mode else NULL,
      consensus_log2fc = if (method == "both") consensus_log2fc else NULL,
      active_source = result_source
    ),
    group_column  = group_column,
    group_palette = group_palette,
    validate      = FALSE
  )

  meta <- S4Vectors::metadata(v)
  meta$de_results_by_method <- lapply(comparisons_by_source, S4Vectors::SimpleList)
  meta$de_summary_by_method <- lapply(deg_summary_by_source, S4Vectors::SimpleList)
  meta$de_active_source <- result_source
  meta$de_available_sources <- names(comparisons_by_source)
  meta$design <- list(
    formula = if (is.null(design_formula)) NULL else paste(deparse(design_formula), collapse = ""),
    covariates = covariates
  )

  # store comparison colors
  comp_names <- names(active_comparisons)
  if (!is.null(comp_names) && length(comp_names)) {
    comp_cols <- colorspace::qualitative_hcl(length(comp_names), palette = comparison_palette)
    names(comp_cols) <- comp_names
    meta$comparison <- list(
      palette = comparison_palette,
      colors = comp_cols
    )
  }
  S4Vectors::metadata(v) <- meta

  # Optional, dev-only check
  if (getOption("vista.strict", FALSE)) {
    rn <- rownames(v)
    meta <- S4Vectors::metadata(v)
    for (nm in names(meta$de_results)) {
      dfi <- meta$de_results[[nm]]
      if (!identical(rownames(dfi), rn)) {
        warning(sprintf("Rownames mismatch in metadata(v)$de_results[['%s']].", nm))
      }
    }
  }

  if (isTRUE(validate)) {
    validate_vista(v, level = "full", error = TRUE)
  }

  v
}

#' Retrieve an expression matrix from a VISTA object
#'
#' Returns the specified assay (default `"norm_counts"`) with optional gene/sample
#' subsetting, replicate summarization by group, and simple transformations.
#'
#' @param x A `VISTA` object.
#' @param assay_name Assay to extract (default: `"norm_counts"`).
#' @param genes Optional character vector of gene IDs to keep.
#' @param samples Optional character vector of sample IDs (or group labels when
#'   `summarise = TRUE`) to keep.
#' @param group_column Optional column in `sample_info` used when summarising
#'   replicates. Defaults to the stored `group_column`.
#' @param summarise Logical; if `TRUE`, averages replicates per `group_column`.
#' @param transform One of `"none"`, `"log2"`, or `"zscore"` applied after
#'   subsetting and optional summarisation. Z-scores are computed within the
#'   selected samples/columns.
#'
#' @return A numeric matrix with genes in rows and samples (or groups) in columns.
#' @export
get_expression_matrix <- function(x,
                                  assay_name = "norm_counts",
                                  genes = NULL,
                                  samples = NULL,
                                  group_column = NULL,
                                  summarise = FALSE,
                                  transform = c("none", "log2", "zscore")) {
  stopifnot(inherits(x, "VISTA"))
  transform <- match.arg(transform)

  mat <- SummarizedExperiment::assay(x, assay_name)

  # subset genes
  if (!is.null(genes)) {
    mat <- mat[intersect(rownames(mat), genes), , drop = FALSE]
  }

  # subset samples or summarise by group
  meta <- as.data.frame(SummarizedExperiment::colData(x))
  meta$sample <- rownames(meta)
  group_info <- S4Vectors::metadata(x)$group %||% list()
  group_col <- group_column %||% group_info$column %||% "group"

  if (summarise) {
    stopifnot(group_col %in% names(meta))
    if (!is.null(samples)) {
      # treat samples as group labels when summarising
      meta <- meta[meta[[group_col]] %in% samples, , drop = FALSE]
    }
    group_map <- split(meta$sample, meta[[group_col]])
    mat <- vapply(group_map, function(samps) rowMeans(mat[, samps, drop = FALSE]),
                  FUN.VALUE = numeric(nrow(mat)))
    mat <- as.matrix(mat)
    rownames(mat) <- if (!is.null(genes)) genes[genes %in% rownames(mat)] else rownames(mat)
  } else {
    if (!is.null(samples)) {
      keep <- intersect(colnames(mat), samples)
      mat <- mat[, keep, drop = FALSE]
    }
  }

  # apply transform within the selected columns
  if (transform == "log2") {
    mat <- log2(mat + 1)
  } else if (transform == "zscore") {
    mat <- t(scale(t(mat)))
  }

  mat
}




#' @title Get Genes by Regulation
#' @description Extract gene IDs by regulation class from selected comparisons in a VISTA object.
#' @param x A VISTA object.
#' @param sample_comparisons Character vector of comparison names to include.
#' @param regulation One of "Up", "Down", "Both", or "All" (default: "Both").
#' @return A named list of character vectors (one per comparison).
#' @export
get_genes_by_regulation <- function(x, sample_comparisons, regulation = "Both") {
  stopifnot(inherits(x, "VISTA"))

  comps <- S4Vectors::metadata(x)$de_results
  if (is.null(comps) || !length(comps)) {
    cli::cli_abort("No DE results found in metadata(x)$de_results.")
  }
  bad <- setdiff(sample_comparisons, names(comps))
  if (length(bad)) cli::cli_abort("Unknown comparisons: {.val {bad}}")

  regulation <- match.arg(regulation, choices = c("Up", "Down", "Both", "All"))

  # cutoffs from metadata if available, else defaults
  cuts <- cutoffs(x)
  lfc_cut <- if (!is.null(cuts) && !is.null(cuts$log2fc)) cuts$log2fc else 1
  p_cut   <- if (!is.null(cuts) && !is.null(cuts$pval)) cuts$pval else 0.05

  .ensure_regul <- function(df) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    # gene_id
    if (!"gene_id" %in% names(df)) {
      id_candidates <- c("gene","geneID","GeneID","symbol","SYMBOL","ensembl","ensembl_id","ID","id")
      hit <- id_candidates[id_candidates %in% names(df)][1]
      if (!is.na(hit)) names(df)[names(df) == hit] <- "gene_id"
    }
    if (!"gene_id" %in% names(df)) {
      rn <- rownames(df)
      if (!is.null(rn) && length(rn) == nrow(df) && all(nzchar(rn)) && !anyDuplicated(rn))
        df$gene_id <- rn
    }
    if (!"gene_id" %in% names(df)) cli::cli_abort("A DE table lacks a 'gene_id' column and valid rownames.")

    # FC & p columns
    fc_col <- if ("log2fc" %in% names(df)) "log2fc" else
      if ("log2FoldChange" %in% names(df)) "log2FoldChange" else
        if ("logFC" %in% names(df)) "logFC" else NA_character_
    p_col  <- if ("pvalue" %in% names(df)) "pvalue" else
      if ("padj" %in% names(df)) "padj" else NA_character_
    if (is.na(fc_col) || is.na(p_col)) {
      if (!"regulation" %in% names(df)) {
        cli::cli_abort("Cannot compute regulation: missing log2FC/p-value columns and no 'regulation' present.")
      }
    } else if (!"regulation" %in% names(df)) {
      xv <- suppressWarnings(as.numeric(gsub(",", "", as.character(df[[fc_col]]), fixed = TRUE)))
      yv <- suppressWarnings(as.numeric(gsub(",", "", as.character(df[[p_col]]),  fixed = TRUE)))
      yv[!is.finite(yv) | is.na(yv)] <- 1
      df$regulation <- ifelse(xv >=  lfc_cut & yv <= p_cut, "Up",
                              ifelse(xv <= -lfc_cut & yv <= p_cut, "Down", "Other"))
    }
    df
  }

  out <- purrr::map(sample_comparisons, function(comp) {
    df <- .ensure_regul(comps[[comp]])
    if (regulation == "All") {
      genes <- df$gene_id
    } else if (regulation == "Both") {
      genes <- df$gene_id[df$regulation %in% c("Up", "Down")]
    } else {
      genes <- df$gene_id[df$regulation %in% regulation]
    }
    unique(stats::na.omit(as.character(genes)))
  })
  names(out) <- sample_comparisons
  out
}
