#' Run Differential Expression Analysis with DESeq2, edgeR, or limma-voom
#'
#' Perform differential expression (DE) analysis across multiple group comparisons
#' using DESeq2, edgeR, or limma-voom. These functions process raw count data,
#' normalize it, execute pairwise group-level tests, and return standardized DEG outputs
#' compatible with \code{VISTA}-based visualization and analysis.
#'
#' @description
#' These functions encapsulate the standard RNA-seq analysis workflow using
#' DESeq2 (\link{run_deseq_analysis}), edgeR (\link{run_edger_analysis}), or
#' limma-voom (\link{run_limma_analysis}), including:
#' gene filtering, design matrix setup, normalization, model fitting, differential testing,
#' DEG classification (\code{"Up"}, \code{"Down"}, \code{"Other"}), and result formatting.
#'
#' Both methods return output in a harmonized structure ready for downstream use
#' in \link{create_vista} or standalone DEG summaries.
#'
#' @param counts A data frame or matrix of raw counts with one gene per row.
#'               Must include a column defined by \code{column_geneid}, and column names
#'               must match entries in \code{sample_info$sample_names}.
#' @param column_geneid A string identifying the column name containing gene identifiers.
#' @param sample_info A data frame with sample metadata. Must contain \code{sample_names} and
#'        the specified grouping column.
#' @param group_column The name of the column in \code{sample_info} that defines experimental groups.
#' @param group_numerator A character vector of numerator group(s) for fold-change comparisons.
#' @param group_denominator A character vector of denominator group(s) for fold-change comparisons.
#' @param covariates Optional character vector of additional sample_info columns to adjust for.
#' @param design_formula Optional model formula (or formula string). When provided,
#'   it overrides automatic design construction from \code{group_column} +
#'   \code{covariates}. Must include \code{group_column}.
#' @param min_counts Minimum per-sample read count a gene must reach to count as
#'   detected in that sample. Also applied as a minimum row total across all
#'   samples in an initial pre-filter. Default: \code{10}.
#' @param min_replicates Minimum number of samples in which a gene must reach
#'   \code{min_counts} for that gene to enter the model. Counted across the whole
#'   experiment, not within each group. Default: \code{1}.
#' @param log2fc_cutoff Absolute log2 fold-change threshold to define DEGs. Default: \code{1}.
#' @param pval_cutoff P-value or adjusted p-value cutoff for significance. Default: \code{0.05}.
#' @param p_value_type For DESeq2: one of \code{"padj"} or \code{"pvalue"}.
#'   For edgeR/limma: one of \code{"FDR"} or \code{"PValue"}.
#'
#' @return A named list with components:
#' \itemize{
#'   \item \code{norm_counts}: Matrix of normalized expression values
#'     (CPM for edgeR/limma, DESeq2-normalized counts).
#'   \item \code{sample_info}: Updated sample metadata.
#'   \item \code{row_data}: Gene-level metadata, including mean expression.
#'   \item \code{comparisons}: Named list of DEG result tibbles (one per comparison), each containing
#'         standardized columns: \code{gene_id}, \code{log2fc}, \code{pvalue}, \code{p.adj}, and \code{regulation}.
#'   \item \code{deg_summary}: List of summary tables showing DEG regulation counts.
#' }
#'
#' @details
#' \itemize{
#'   \item For DESeq2, normalization is performed via \link[DESeq2]{DESeq}, and DE testing uses \link[DESeq2]{results}.
#'   \item For edgeR, normalization uses \link[edgeR]{calcNormFactors}, and testing uses \link[edgeR]{glmLRT}.
#'   \item For limma, normalization uses \link[edgeR]{calcNormFactors} + \link[limma]{voom},
#'     and testing uses \link[limma]{eBayes}.
#' }
#' Low-abundance filtering is applied before model fitting, and for the
#' edgeR/limma backends before \link[edgeR]{calcNormFactors} so that
#' normalization factors are not estimated from genes that will be discarded.
#' All three backends use the same filtering predicate: a gene is retained when
#' at least \code{min_replicates} samples reach \code{min_counts}.
#' Gene regulation status is determined via \code{.categorize_deg_results()} based on user thresholds.
#'
#' All output comparison results are internally standardized via \code{.tidy_de_results()} to ensure
#' a uniform column schema compatible with VISTA plotting tools.
#'
#' @seealso \link{create_vista}, \link[DESeq2]{DESeq}, \link[edgeR]{glmLRT}, \link[limma]{voom}
#'
#' @examples
#' v <- example_vista()
#' si <- as.data.frame(sample_info(v))
#' data("count_data", package = "VISTA")
#' counts_small <- count_data[seq_len(200), c("gene_id", si$sample_names), drop = FALSE]
#' limma_results <- run_limma_analysis(
#'   counts = counts_small,
#'   sample_info = si,
#'   column_geneid = "gene_id",
#'   group_column = "cond_long",
#'   group_numerator = "treatment1",
#'   group_denominator = "control",
#'   min_counts = 5,
#'   min_replicates = 1
#' )
#' names(limma_results$comparisons)
#'
#' @importFrom DESeq2 DESeq DESeqDataSetFromMatrix results counts estimateSizeFactors
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmQLFit glmQLFTest topTags
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment assay
#' @importFrom tibble as_tibble rownames_to_column column_to_rownames
#' @importFrom dplyr filter mutate arrange left_join select
#' @importFrom cli cli_abort cli_warn
#' @importFrom rlang `%||%` abort
#' @importFrom stats setNames
#'
#' @name run_deseq_analysis
NULL

#' @keywords internal
#' @noRd
.coerce_design_formula <- function(design_formula) {
  if (is.null(design_formula)) return(NULL)
  if (inherits(design_formula, "formula")) return(design_formula)
  if (is.character(design_formula) && length(design_formula) == 1L && nzchar(design_formula)) {
    return(stats::as.formula(design_formula))
  }
  cli::cli_abort("{.arg design_formula} must be NULL, a formula, or a single formula string.")
}

#' @keywords internal
.normalize_sample_info <- function(sample_info, counts = NULL, column_geneid = NULL) {
  if (!is.data.frame(sample_info)) {
    sample_info <- as.data.frame(sample_info, stringsAsFactors = FALSE)
  } else {
    sample_info <- as.data.frame(sample_info, stringsAsFactors = FALSE)
  }

  .is_default_row_names <- function(rn) {
    !is.null(rn) &&
      length(rn) == nrow(sample_info) &&
      identical(as.character(rn), as.character(seq_len(nrow(sample_info))))
  }

  .infer_sample_name_col <- function(sample_info, counts = NULL, column_geneid = NULL) {
    candidate_cols <- intersect(
      c("sample_names", "sample_name", "sample", "sample_id", "Run", "run", "Sample", "SampleName"),
      colnames(sample_info)
    )
    if (!length(candidate_cols)) {
      return(NULL)
    }

    valid_ids <- function(x) {
      x <- as.character(x)
      length(x) == nrow(sample_info) &&
        !anyNA(x) &&
        all(nzchar(x)) &&
        !anyDuplicated(x)
    }

    if (!is.null(counts)) {
      sample_cols <- colnames(counts)
      if (!is.null(column_geneid) && nzchar(column_geneid)) {
        sample_cols <- setdiff(sample_cols, column_geneid)
      }
      for (nm in candidate_cols) {
        vals <- as.character(sample_info[[nm]])
        if (valid_ids(vals) && all(vals %in% sample_cols)) {
          return(nm)
        }
      }
    }

    for (nm in candidate_cols) {
      vals <- as.character(sample_info[[nm]])
      if (valid_ids(vals)) {
        return(nm)
      }
    }
    NULL
  }

  if (!"sample_names" %in% colnames(sample_info)) {
    rn <- rownames(sample_info)
    has_valid_rn <- !is.null(rn) &&
      length(rn) == nrow(sample_info) &&
      !anyNA(rn) &&
      all(nzchar(rn)) &&
      !anyDuplicated(rn)

    if (has_valid_rn && !.is_default_row_names(rn)) {
      sample_info$sample_names <- as.character(rn)
    } else {
      inferred_col <- .infer_sample_name_col(
        sample_info = sample_info,
        counts = counts,
        column_geneid = column_geneid
      )
      if (!is.null(inferred_col)) {
        sample_info$sample_names <- as.character(sample_info[[inferred_col]])
      } else {
        cli::cli_abort(
          c(
            "{.arg sample_info} must contain a {.field sample_names} column, or have unique non-default rownames.",
            "i" = "Alternatively, provide a sample identifier column (e.g. {.field Run}) with values matching count columns."
          )
        )
      }
    }
  }

  sample_info$sample_names <- as.character(sample_info$sample_names)
  if (anyNA(sample_info$sample_names) || any(!nzchar(sample_info$sample_names)) || anyDuplicated(sample_info$sample_names)) {
    cli::cli_abort("{.field sample_names} in {.arg sample_info} must be unique, non-empty, and non-NA.")
  }

  if (!is.null(counts)) {
    sample_cols <- colnames(counts)
    if (!is.null(column_geneid) && nzchar(column_geneid)) {
      sample_cols <- setdiff(sample_cols, column_geneid)
    }
    missing_samples <- setdiff(sample_info$sample_names, sample_cols)
    if (length(missing_samples)) {
      cli::cli_abort(
        c(
          "{.field sample_names} in {.arg sample_info} must match count column names.",
          "x" = "Missing in {.arg counts}: {.val {utils::head(missing_samples, 10)}}",
          "i" = "First available count columns: {.val {utils::head(sample_cols, 10)}}"
        )
      )
    }
  }

  sample_info
}

#' @keywords internal
.validate_covariates <- function(sample_info, group_column, covariates = NULL) {
  if (is.null(covariates)) return(character(0))
  if (is.character(covariates) && !length(covariates)) return(character(0))
  if (!is.character(covariates)) {
    cli::cli_abort("{.arg covariates} must be NULL or a character vector.")
  }
  covariates <- unique(covariates)
  missing_cov <- setdiff(covariates, colnames(sample_info))
  if (length(missing_cov)) {
    cli::cli_abort("Covariate column(s) not found in {.arg sample_info}: {.val {missing_cov}}")
  }
  setdiff(covariates, group_column)
}

#' @keywords internal
.resolve_design_formula <- function(sample_info,
                                    group_column,
                                    covariates = NULL,
                                    design_formula = NULL,
                                    backend = c("deseq2", "edger", "limma")) {
  backend <- match.arg(backend)
  design_formula <- .coerce_design_formula(design_formula)
  covariates <- .validate_covariates(sample_info, group_column, covariates)

  if (is.null(design_formula)) {
    terms <- unique(c(covariates, group_column))
    f_txt <- if (backend %in% c("edger", "limma")) {
      paste("~ 0 +", paste(terms, collapse = " + "))
    } else {
      paste("~", paste(terms, collapse = " + "))
    }
    return(stats::as.formula(f_txt))
  }

  design_vars <- all.vars(design_formula)
  if (!group_column %in% design_vars) {
    cli::cli_abort("{.arg design_formula} must include {.field {group_column}}.")
  }
  missing_vars <- setdiff(design_vars, colnames(sample_info))
  if (length(missing_vars)) {
    cli::cli_abort("Variable(s) in {.arg design_formula} not found in {.arg sample_info}: {.val {missing_vars}}")
  }

  if (backend %in% c("edger", "limma")) {
    f_txt <- paste(deparse(design_formula), collapse = "")
    if (grepl("[:*]", f_txt)) {
      cli::cli_abort("Interactions in {.arg design_formula} are not currently supported for edgeR. Use additive terms only.")
    }
    rhs_txt <- trimws(strsplit(f_txt, "~", fixed = TRUE)[[1]][2])
    has_no_intercept <- grepl("(^|\\+)\\s*0\\s*\\+|\\-\\s*1", rhs_txt)
    if (!has_no_intercept) {
      design_formula <- stats::as.formula(paste("~ 0 +", rhs_txt))
    }
  }

  design_formula
}

#' @keywords internal
.edgeR_group_contrast <- function(design_colnames, group_column, numerator, denominator) {
  contrast <- numeric(length(design_colnames))
  std_names <- make.names(design_colnames)
  num_coef <- make.names(paste0(group_column, numerator))
  den_coef <- make.names(paste0(group_column, denominator))
  idx_num <- which(std_names == num_coef)
  idx_den <- which(std_names == den_coef)

  if (length(idx_num) > 1L || length(idx_den) > 1L) {
    cli::cli_abort(
      paste0(
        "Ambiguous edgeR contrast for ", numerator, " vs ", denominator, ". ",
        "Expected unique design columns for ", group_column, " levels. ",
        "Available design columns: ", paste(design_colnames, collapse = ", ")
      )
    )
  }

  # With treatment coding, one side can be the implicit reference level and
  # therefore absent from the design matrix.
  if (length(idx_num) == 1L && length(idx_den) == 1L) {
    contrast[idx_num] <- 1
    contrast[idx_den] <- -1
    return(contrast)
  }
  if (length(idx_num) == 1L && length(idx_den) == 0L) {
    contrast[idx_num] <- 1
    return(contrast)
  }
  if (length(idx_num) == 0L && length(idx_den) == 1L) {
    contrast[idx_den] <- -1
    return(contrast)
  }

  if (length(idx_num) == 0L && length(idx_den) == 0L) {
    cli::cli_abort(
      paste0(
        "Could not construct edgeR contrast for ",
        numerator, " vs ", denominator, ". ",
        "Expected design columns for ", group_column, " levels were not found. ",
        "Available design columns: ", paste(design_colnames, collapse = ", ")
      )
    )
  }
  contrast
}

#' @rdname run_deseq_analysis
#' @export
run_deseq_analysis <- function(
    counts,
    sample_info,
    column_geneid,
    group_column,
    group_numerator,
    group_denominator,
    covariates = NULL,
    design_formula = NULL,
    min_counts = 10,
    min_replicates = 1,
    log2fc_cutoff = 1,
    pval_cutoff = 0.05,
    p_value_type = "padj"
) {

  p_value_type <- match.arg(p_value_type, c("padj", "pvalue"))
  sample_info <- .normalize_sample_info(
    sample_info = sample_info,
    counts = counts,
    column_geneid = column_geneid
  )
  if (is.factor(group_numerator)) group_numerator <- as.character(group_numerator)
  if (is.factor(group_denominator)) group_denominator <- as.character(group_denominator)

  # validate arguments

  .validate_de_inputs(
    counts = counts,
    sample_info = sample_info,
    column_geneid = column_geneid,
    group_column = group_column,
    group_numerator = group_numerator,
    group_denominator = group_denominator,
    p_value_type = p_value_type
  )

  # Preprocess count data
  count_data <- counts %>%
    dplyr::select(tidyselect::all_of(c(column_geneid, sample_info$sample_names))) %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    tibble::column_to_rownames(var = column_geneid)
  count_data[is.na(count_data)] <- 0
  count_data <- count_data[rowSums(count_data) >= min_counts, ]

  # Prepare DESeq2 input
  col_data <- sample_info %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    tibble::column_to_rownames("sample_names")
  col_data[[group_column]] <- factor(col_data[[group_column]])

  de_design <- .resolve_design_formula(
    sample_info = col_data,
    group_column = group_column,
    covariates = covariates,
    design_formula = design_formula,
    backend = "deseq2"
  )

  # Avoid DESeq2 warnings by explicitly casting character design variables.
  design_vars <- intersect(all.vars(de_design), colnames(col_data))
  for (vv in design_vars) {
    if (is.character(col_data[[vv]])) {
      col_data[[vv]] <- factor(col_data[[vv]])
    }
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = as.matrix(count_data),
    colData = col_data,
    design = de_design
  )
  dds <- dds[rowSums(DESeq2::counts(dds) >= min_counts) >= min_replicates, ]
  dds <- DESeq2::DESeq(dds)

  # Normalized counts
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)

  sample_annot <- as.data.frame(SummarizedExperiment::colData(dds))

  # Create rowData
  row_data <- S4Vectors::DataFrame(
    baseMean = rowMeans(norm_counts),
    row.names = rownames(norm_counts)
  )


  # perform all DE comparison
  deseq_comps <- .run_deseq_comparisons(dds = dds,
                                        group_column = group_column,
                                        group_numerator = group_numerator,
                                        group_denominator = group_denominator,
                                        log2fc_cutoff = log2fc_cutoff,
                                        pval_cutoff = pval_cutoff,
                                        p_value_type = p_value_type)

  # put all deseq related stuff required for vista in a list

  deseq_results <- list(norm_counts = norm_counts,
                        raw_counts = DESeq2::counts(dds, normalized = FALSE),
                        sample_info = sample_annot,
                        row_data = row_data,
                        comparisons = deseq_comps$comparisons,
                        deg_summary = deseq_comps$deg_summary)

  deseq_results

}

#' @rdname run_deseq_analysis
#' @export
run_edger_analysis <- function(
    counts,
    sample_info,
    column_geneid,
    group_column,
    group_numerator,
    group_denominator,
    covariates = NULL,
    design_formula = NULL,
    min_counts = 10,
    min_replicates = 1,
    log2fc_cutoff = 1,
    pval_cutoff = 0.05,
    p_value_type = "FDR"
) {
  p_value_type <- match.arg(p_value_type, c("FDR", "PValue"))
  sample_info <- .normalize_sample_info(
    sample_info = sample_info,
    counts = counts,
    column_geneid = column_geneid
  )
  if (is.factor(group_numerator)) group_numerator <- as.character(group_numerator)
  if (is.factor(group_denominator)) group_denominator <- as.character(group_denominator)

  # Validate inputs
  .validate_de_inputs(
    counts = counts,
    sample_info = sample_info,
    column_geneid = column_geneid,
    group_column = group_column,
    group_numerator = group_numerator,
    group_denominator = group_denominator,
    p_value_type = p_value_type
  )

  if (!requireNamespace("edgeR", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg edgeR} is required for edgeR analysis. Please install it or choose method = 'deseq2'.")
  }

  # Preprocess counts
  count_data <- counts %>%
    dplyr::select(tidyselect::all_of(c(column_geneid, sample_info$sample_names))) %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    tibble::column_to_rownames(column_geneid)
  count_data[is.na(count_data)] <- 0
  count_data <- count_data[rowSums(count_data) >= min_counts, ]

  # Metadata
  col_data <- sample_info %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    tibble::column_to_rownames("sample_names")
  col_data[[group_column]] <- factor(col_data[[group_column]])

  de_design <- .resolve_design_formula(
    sample_info = col_data,
    group_column = group_column,
    covariates = covariates,
    design_formula = design_formula,
    backend = "edger"
  )

  # Keep edgeR design terms explicit and stable across runs.
  design_vars <- intersect(all.vars(de_design), colnames(col_data))
  for (vv in design_vars) {
    if (is.character(col_data[[vv]])) {
      col_data[[vv]] <- factor(col_data[[vv]])
    }
  }

  dge <- edgeR::DGEList(counts = count_data, group = col_data[[group_column]])

  # Filter before calcNormFactors(): normalization factors estimated on
  # unfiltered low-count genes are biased (see the edgeR user's guide).
  # The predicate matches the DESeq2 backend so `min_counts`/`min_replicates`
  # mean the same thing across all three engines.
  keep <- rowSums(dge$counts >= min_counts) >= min_replicates
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)

  design <- stats::model.matrix(de_design, data = col_data)

  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmFit(dge, design)

  norm_counts <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)

  sample_annot <- col_data


  row_data <- S4Vectors::DataFrame(
    baseMean = rowMeans(norm_counts),
    row.names = rownames(norm_counts)
  )


  # Run all comparisons
  comps <- list()
  deg_summary <- list()

  for (i in seq_along(group_numerator)) {
    num <- group_numerator[[i]]
    den <- group_denominator[[i]]
    if (identical(num, den)) next

    contrast_vector <- .edgeR_group_contrast(
      design_colnames = colnames(design),
      group_column = group_column,
      numerator = num,
      denominator = den
    )
    lrt <- edgeR::glmLRT(fit, contrast = contrast_vector)
    res <- edgeR::topTags(lrt, n = Inf)$table

    comp_name <- paste(num, den, sep = "_VS_")
    res$gene_id <- rownames(res)

    categorized <- .categorize_deg_results(
      de_results = res,
      log2fc_cutoff = log2fc_cutoff,
      pval_cutoff = pval_cutoff,
      p_value_type = p_value_type
    )

    comps[[comp_name]] <- categorized
    deg_summary[[comp_name]] <- categorized %>%
      dplyr::count(regulation)
  }

  list(
    norm_counts = norm_counts,
    raw_counts = as.matrix(dge$counts),
    sample_info = sample_annot,
    row_data = row_data,
    comparisons = comps,
    deg_summary = deg_summary
  )
}

#' @rdname run_deseq_analysis
#' @export
run_limma_analysis <- function(
    counts,
    sample_info,
    column_geneid,
    group_column,
    group_numerator,
    group_denominator,
    covariates = NULL,
    design_formula = NULL,
    min_counts = 10,
    min_replicates = 1,
    log2fc_cutoff = 1,
    pval_cutoff = 0.05,
    p_value_type = "FDR"
) {
  p_value_type <- match.arg(p_value_type, c("FDR", "PValue"))
  sample_info <- .normalize_sample_info(
    sample_info = sample_info,
    counts = counts,
    column_geneid = column_geneid
  )
  if (is.factor(group_numerator)) group_numerator <- as.character(group_numerator)
  if (is.factor(group_denominator)) group_denominator <- as.character(group_denominator)

  .validate_de_inputs(
    counts = counts,
    sample_info = sample_info,
    column_geneid = column_geneid,
    group_column = group_column,
    group_numerator = group_numerator,
    group_denominator = group_denominator,
    p_value_type = p_value_type
  )

  if (!requireNamespace("limma", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg limma} is required for limma analysis. Please install it or choose another method.")
  }
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg edgeR} is required for limma-voom analysis. Please install it or choose another method.")
  }

  count_data <- counts %>%
    dplyr::select(tidyselect::all_of(c(column_geneid, sample_info$sample_names))) %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    tibble::column_to_rownames(column_geneid)
  count_data[is.na(count_data)] <- 0
  count_data <- count_data[rowSums(count_data) >= min_counts, , drop = FALSE]

  col_data <- sample_info %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    tibble::column_to_rownames("sample_names")
  col_data[[group_column]] <- factor(col_data[[group_column]])

  de_design <- .resolve_design_formula(
    sample_info = col_data,
    group_column = group_column,
    covariates = covariates,
    design_formula = design_formula,
    backend = "limma"
  )

  design_vars <- intersect(all.vars(de_design), colnames(col_data))
  for (vv in design_vars) {
    if (is.character(col_data[[vv]])) {
      col_data[[vv]] <- factor(col_data[[vv]])
    }
  }

  dge <- edgeR::DGEList(counts = count_data)

  # Filter before calcNormFactors(); see run_edger_analysis() for rationale.
  keep <- rowSums(dge$counts >= min_counts) >= min_replicates
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)

  design <- stats::model.matrix(de_design, data = col_data)
  v <- limma::voom(dge, design = design, plot = FALSE)
  fit <- limma::lmFit(v, design)

  comps <- list()
  deg_summary <- list()

  for (i in seq_along(group_numerator)) {
    num <- group_numerator[[i]]
    den <- group_denominator[[i]]
    if (identical(num, den)) next

    contrast_vector <- .edgeR_group_contrast(
      design_colnames = colnames(design),
      group_column = group_column,
      numerator = num,
      denominator = den
    )
    contrast_matrix <- matrix(contrast_vector, ncol = 1, dimnames = list(colnames(design), "contrast"))
    fit2 <- limma::contrasts.fit(fit, contrasts = contrast_matrix)
    fit2 <- limma::eBayes(fit2)
    res <- limma::topTable(fit2, coef = 1, number = Inf, sort.by = "none")

    comp_name <- paste(num, den, sep = "_VS_")
    res$gene_id <- rownames(res)
    res$PValue <- res$P.Value
    res$FDR <- res$adj.P.Val

    categorized <- .categorize_deg_results(
      de_results = res,
      log2fc_cutoff = log2fc_cutoff,
      pval_cutoff = pval_cutoff,
      p_value_type = p_value_type
    )

    comps[[comp_name]] <- categorized
    deg_summary[[comp_name]] <- categorized %>%
      dplyr::count(regulation)
  }

  norm_counts <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)

  list(
    norm_counts = norm_counts,
    raw_counts = as.matrix(dge$counts),
    sample_info = col_data,
    row_data = S4Vectors::DataFrame(
      baseMean = rowMeans(norm_counts),
      row.names = rownames(dge$counts)
    ),
    comparisons = comps,
    deg_summary = deg_summary
  )
}




#' @keywords internal
.validate_de_inputs <- function(
    counts,
    sample_info,
    column_geneid,
    group_column,
    group_numerator,
    group_denominator,
    p_value_type
) {
  if (is.factor(group_numerator)) group_numerator <- as.character(group_numerator)
  if (is.factor(group_denominator)) group_denominator <- as.character(group_denominator)

  if (!is.character(group_numerator) || !length(group_numerator)) {
    cli::cli_abort("{.arg group_numerator} must be a non-empty character vector.")
  }
  if (!is.character(group_denominator) || !length(group_denominator)) {
    cli::cli_abort("{.arg group_denominator} must be a non-empty character vector.")
  }
  if (length(group_numerator) != length(group_denominator)) {
    cli::cli_abort(
      "{.arg group_numerator} and {.arg group_denominator} must have equal length. Got {length(group_numerator)} and {length(group_denominator)}."
    )
  }

  # Check column_geneid
  if (!is.character(column_geneid) || length(column_geneid) != 1 || !column_geneid %in% colnames(counts)) {
    cli::cli_abort(
      c(
        "!" = "{.arg column_geneid} must be a single character string and must exist in {.arg counts}.",
        "x" = "Column {.val {column_geneid}} not found in {.arg counts}.",
        "i" = "Available columns in {.arg counts}: {.val {colnames(counts)}}"
      )
    )
  }

  # Check sample_info structure
  if (!is.character(group_column) || length(group_column) != 1) {
    cli::cli_abort("{.arg group_column} must be a single character string.")
  }

  if (!group_column %in% colnames(sample_info)) {
    cli::cli_abort("Column {.val {group_column}} not found in {.arg sample_info}.")
  }

  if (!"sample_names" %in% colnames(sample_info)) {
    cli::cli_abort("{.arg sample_info} must contain a column named {.field sample_names}.")
  }

  # Check p-value type
  if (!p_value_type %in% c("padj", "pvalue","FDR","PValue")) {
    cli::cli_abort("{.arg p_value_type} must be one of {.val {'padj', 'pvalue', 'FDR', 'PValue'}}.")
  }

  # Check if all requested groups exist
  all_groups <- unique(sample_info[[group_column]])
  requested_groups <- unique(c(group_numerator, group_denominator))
  missing_groups <- setdiff(requested_groups, all_groups)

  if (length(missing_groups) > 0) {
    cli::cli_abort(
      c(
        "!" = "Some groups specified in {.field group_numerator} or {.field group_denominator} are not found in column {.field {group_column}} of {.arg sample_info}.",
        "x" = "Missing group(s): {.val {missing_groups}}",
        "i" = "Available group(s) in column {.field {group_column}}: {.val {all_groups}}"
      )
    )
  }

  invisible(TRUE)
}



# Internal helper for fold-change gene validation
.available_fc_genes <- function(x, sample_comparisons = NULL) {
  comps <- comparisons(x)
  if (!length(comps)) {
    return(character())
  }

  if (is.null(sample_comparisons)) {
    sample_comparisons <- names(comps)
  }

  available <- unique(unlist(lapply(sample_comparisons, function(comp_name) {
    tbl <- as.data.frame(comps[[comp_name]], stringsAsFactors = FALSE)
    if ("gene_id" %in% names(tbl)) {
      ids <- as.character(tbl$gene_id)
    } else {
      ids <- rownames(tbl)
    }
    ids[!is.na(ids) & nzchar(ids)]
  }), use.names = FALSE))

  if (!length(available)) {
    rd <- as.data.frame(row_data(x), stringsAsFactors = FALSE)
    if ("gene_id" %in% names(rd)) {
      available <- as.character(rd$gene_id)
    } else {
      available <- rownames(x)
    }
  }

  unique(available[!is.na(available) & nzchar(available)])
}

# Internal validator for get_foldchange_matrix
.validate_fc_inputs <- function(x, sample_comparisons, genes) {
  comps <- comparisons(x)
  all_comp_names <- names(comps)

  if (!is.null(sample_comparisons)) {
    missing_comps <- setdiff(sample_comparisons, all_comp_names)
    if (length(missing_comps) > 0) {
      cli::cli_abort(c(
        "!" = "Some {.arg sample_comparisons} not found in {.fun comparisons(x)}.",
        "x" = "Missing: {.val {missing_comps}}",
        "i" = "Available comparisons: {.val {all_comp_names}}"
      ))
    }
  }

  all_genes <- .available_fc_genes(x, sample_comparisons = sample_comparisons)
  if (!is.null(genes)) {
    missing_genes <- setdiff(genes, all_genes)
    if (length(missing_genes) > 0) {
      cli::cli_abort(c(
        "!" = "Some {.arg genes} not found in the VISTA object.",
        "x" = "Missing: {.val {missing_genes}}",
        "i" = "Available genes: e.g. {.val {utils::head(all_genes, 5)}} ..."
      ))
    }
  }
}


#' Prepare PCA result as a tidy data frame
#'
#' Converts PCA output and metadata into a tidy format suitable for ggplot2.
#' Includes PC1 and PC2 scores and merges with sample metadata.
#'
#' @param pca A PCA object returned by [stats::prcomp()].
#' @param meta A sample metadata data frame, typically from `.prepare_sample_metadata()`.
#'
#' @return A data frame containing PCA components and sample metadata.
#' @keywords internal
.prepare_pca_dataframe <- function(pca, meta) {
  as.data.frame(pca$x[, seq_len(2), drop = FALSE]) %>%
    tibble::rownames_to_column("sample") %>%
    dplyr::left_join(meta, by = "sample")
}

#' Enhanced volcano plot with smart column detection & coloring
#'
#' Internal helper that wraps EnhancedVolcano and colors points by regulation.
#' - If `regulation_col` is present (e.g., "Up"/"Down"/"Other"), it is used.
#' - Otherwise, regulation is derived from FC/P cutoffs.
#'
#' @param toptable data.frame with DE results.
#' @param lab Character vector of labels (same length as nrow(toptable)).
#' @param x FC column name (auto-detected if NULL).
#' @param y P/P-adj column name (auto-detected if NULL).
#' @param pCutoff Numeric; P-value cutoff. Default 1e-4.
#' @param FCcutoff Numeric; symmetric scalar or c(down, up). Default 1.5.
#' @param col_by_regul Logical; color by regulation. Default TRUE.
#' @param regulation_col Optional column in `toptable` to use directly (expects levels Up/Down/Other).
#' @param col_up, col_down, col_others Colors.
#' @param return_keyvals Logical; invisibly return the computed `keyvals`. Default FALSE.
#' @param ... Passed to EnhancedVolcano::EnhancedVolcano()
#' @return ggplot object (and invisibly the `keyvals` if `return_keyvals=TRUE`)
#' @keywords internal
# replace your helper with this version
.EnhancedVolcano2 <- function(toptable,
                              lab,
                              x = NULL,
                              y = NULL,
                              pCutoff = 1e-4,
                              FCcutoff = 1.5,
                              col_by_regul = TRUE,
                              regulation_col = NULL,
                              col_up = "#b2182b",
                              col_down = "#2166ac",
                              col_others = "#e0e0e0",
                              return_keyvals = FALSE,
                              ...) {

  if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg EnhancedVolcano} must be installed to draw volcano plots.")
  }
  toptable <- as.data.frame(toptable, stringsAsFactors = FALSE)
  dots <- list(...)

  # Backward-compatible argument aliases commonly used in older VISTA vignettes.
  if ("lab_size" %in% names(dots) && !"labSize" %in% names(dots)) {
    dots$labSize <- dots$lab_size
    dots$lab_size <- NULL
  }
  if ("point_size" %in% names(dots) && !"pointSize" %in% names(dots)) {
    dots$pointSize <- dots$point_size
    dots$point_size <- NULL
  }

  # --- auto-detect columns if not supplied
  if (is.null(x)) {
    x_candidates <- c("log2fc","log2FoldChange","LFC","lfc")
    x <- x_candidates[x_candidates %in% names(toptable)][1]
    if (is.na(x)) stop("FC column not found; tried: ", paste(x_candidates, collapse = ", "))
  }
  if (is.null(y)) {
    y_candidates <- c("padj","FDR","qvalue","svalue","pvalue","P.Value","p_val")
    y <- y_candidates[y_candidates %in% names(toptable)][1]
    if (is.na(y)) stop("P-value column not found; tried: ", paste(y_candidates, collapse = ", "))
  }

  # --- robust numeric coercion for FC & P columns
  to_num <- function(v) {
    if (is.numeric(v)) return(v)
    v_chr <- as.character(v)
    v_chr <- gsub(",", "", v_chr, fixed = TRUE)  # remove thousand separators if any
    suppressWarnings(as.numeric(v_chr))
  }

  xv <- to_num(toptable[[x]])
  yv <- to_num(toptable[[y]])

  # if any Inf/-Inf or NA after coercion, make EV-safe:
  if (any(!is.finite(xv), na.rm = TRUE)) {
    finite <- is.finite(xv)
    if (any(finite)) {
      cap <- max(abs(xv[finite]), na.rm = TRUE)
      cap <- if (is.finite(cap) && cap > 0) cap else 1
      xv[!finite] <- sign(xv[!finite]) * cap  # push to finite boundary
    } else {
      xv[!finite] <- 0
    }
  }
  # p-values: coerce NAs to 1 to avoid false positives
  yv[!is.finite(yv) | is.na(yv)] <- 1

  toptable[[x]] <- xv
  toptable[[y]] <- yv

  # --- normalize FCcutoff (allow c(down, up))
  if (length(FCcutoff) == 1) {
    fc_down <- -abs(FCcutoff); fc_up <- abs(FCcutoff)
  } else if (length(FCcutoff) == 2) {
    fc_down <- -abs(FCcutoff[1]); fc_up <- abs(FCcutoff[2])
  } else stop("FCcutoff must be length 1 or 2.")

  # --- build keyvals (use regulation column if present)
  keyvals <- rep(col_others, nrow(toptable)); names(keyvals) <- rep("Other", nrow(toptable))
  if (col_by_regul) {
    if (!is.null(regulation_col) && regulation_col %in% names(toptable)) {
      reg <- as.character(toptable[[regulation_col]])
      keyvals[reg == "Up"]   <- col_up
      keyvals[reg == "Down"] <- col_down
      keyvals[is.na(reg) | !(reg %in% c("Up","Down","Other"))] <- col_others
      names(keyvals)[keyvals == col_up]    <- "Up"
      names(keyvals)[keyvals == col_down]  <- "Down"
      names(keyvals)[keyvals == col_others]<- "Other"
    } else {
      up_idx   <- xv >= fc_up   & yv <= pCutoff
      down_idx <- xv <= fc_down & yv <= pCutoff
      keyvals[up_idx]   <- col_up
      keyvals[down_idx] <- col_down
      names(keyvals)[keyvals == col_up]    <- "Up"
      names(keyvals)[keyvals == col_down]  <- "Down"
      names(keyvals)[keyvals == col_others]<- "Other"
    }
  }

  # --- call EnhancedVolcano (it expects symmetric FCcutoff)
  ev_args <- c(
    list(
      toptable = toptable,
      lab = lab,
      x = x,
      y = y,
      pCutoff = pCutoff,
      FCcutoff = max(abs(fc_up), abs(fc_down)),
      colCustom = if (col_by_regul) keyvals else NULL
    ),
    dots
  )
  plt <- do.call(EnhancedVolcano::EnhancedVolcano, ev_args)

  if (isTRUE(return_keyvals)) {
    return(structure(plt, keyvals = keyvals))
  }
  plt
}



#' Perform pairwise DE comparisons using DESeq2 results
#'
#' @param dds A DESeq2 object after running DESeq().
#' @param group_column The column used for grouping in the design formula.
#' @param group_numerator Vector of numerator group names.
#' @param group_denominator Vector of denominator group names.
#' @param log2fc_cutoff Absolute log2 fold-change cutoff for DEG classification.
#' @param pval_cutoff P-value threshold for DEG classification.
#' @param p_value_type Column to use for p-value filtering ("padj" or "pvalue").
#'
#' @return A list with:
#' \itemize{
#'   \item \code{comparisons}: List of DEG result tibbles (with gene ID and fold changes)
#'   \item \code{deg_summary}: Summary table of DEG counts by regulation
#' }
#' @keywords internal
.run_deseq_comparisons <- function(dds,
                                   group_column,
                                   group_numerator,
                                   group_denominator,
                                   log2fc_cutoff,
                                   pval_cutoff,
                                   p_value_type) {

  comparisons <- list()
  deg_summary <- list()

  for (i in seq_along(group_numerator)) {
    num <- group_numerator[[i]]
    den <- group_denominator[[i]]
    if (identical(num, den)) next

    contrast <- c(group_column, num, den)
    res <- DESeq2::results(dds, contrast = contrast)

    df_deg <- .categorize_deg_results(
      de_results = res,
      log2fc_cutoff = log2fc_cutoff,
      pval_cutoff = pval_cutoff,
      p_value_type = p_value_type
    )

    comp_name <- paste(num, den, sep = "_VS_")
    comparisons[[comp_name]] <- df_deg
    deg_summary[[comp_name]] <- df_deg %>% dplyr::count(regulation)
  }

  list(comparisons = comparisons, deg_summary = deg_summary)
}


#' Categorize Differential Expression Results (DESeq2 or edgeR)
#'
#' This function classifies genes as "Up", "Down", or "Other" based on fold-change and p-value thresholds.
#' It supports input from DESeq2 (`DESeqResults`), edgeR (e.g., `topTags()`), or a general tibble/data frame.
#'
#' @param de_results A differential expression results object. Can be a DESeqResults object,
#'   an edgeR `topTags()` table, or a data frame/tibble containing fold changes and p-values.
#' @param log2fc_cutoff Numeric; minimum absolute log2 fold change to classify a DEG. Default is 1.
#' @param pval_cutoff Numeric; p-value or FDR threshold to define significance. Default is 0.05.
#' @param p_value_type String; which column to use for filtering significance. One of `"padj"`, `"pvalue"`, `"FDR"`, `"PValue"`. Default is `"padj"`.
#'
#' @return A tibble with gene IDs, fold changes, p-values, and a new `regulation` column,
#'   indicating the category of each gene:
#'   - `"Up"`: Genes with log2 fold change >= `log2fc_cutoff` and p-value <= `pval_cutoff`.
#'   - `"Down"`: Genes with log2 fold change <= -`log2fc_cutoff` and p-value <= `pval_cutoff`.
#'   - `"Other"`: Genes not meeting the above criteria.
#' @keywords internal
.categorize_deg_results <- function(de_results,
                                    log2fc_cutoff = 1,
                                    pval_cutoff = 0.05,
                                    p_value_type = c("padj", "pvalue", "FDR", "PValue")) {


  p_value_type <- match.arg(p_value_type)

  # Convert to tibble and preserve gene_id
  if (inherits(de_results, "DESeqResults")) {
    results_tbl <- as_tibble(de_results, rownames = "gene_id")
  } else if (is.data.frame(de_results)) {
    results_tbl <- as_tibble(de_results)
    if (!"gene_id" %in% colnames(results_tbl)) {
      if (!is.null(rownames(de_results))) {
        results_tbl <- tibble::rownames_to_column(results_tbl, "gene_id")
      } else {
        stop("No 'gene_id' column or rownames found to identify genes.")
      }
    }
  } else {
    stop("Unsupported input: must be a DESeqResults object or data.frame.")
  }

  # Determine fold change column
  if ("log2FoldChange" %in% colnames(results_tbl)) {
    fc_col <- "log2FoldChange"
  } else if ("logFC" %in% colnames(results_tbl)) {
    fc_col <- "logFC"
  } else if ("log2fc" %in% colnames(results_tbl)) {
    fc_col <- "log2fc"
  } else {
    stop("No fold change column ('log2FoldChange', 'logFC', or 'log2fc') found.")
  }

  if (!p_value_type %in% colnames(results_tbl)) {
    stop(paste("The selected p-value type", p_value_type, "is not present in the results."))
  }

  results_tbl <- results_tbl %>%
    mutate(
      regulation = dplyr::case_when(
        .data[[fc_col]] >= log2fc_cutoff & .data[[p_value_type]] <= pval_cutoff ~ "Up",
        .data[[fc_col]] <= -log2fc_cutoff & .data[[p_value_type]] <= pval_cutoff ~ "Down",
        TRUE ~ "Other"
      )
    )

  return(results_tbl)
}

#' @title Standardize Differential Expression Results
#' @description Harmonize DE results (DESeq2/edgeR) to a common schema and ensure
#' gene IDs are stored as rownames. Returns S4Vectors::DataFrame.
#' @param tbl data.frame/tibble or S4Vectors::DataFrame with DE results.
#' @param rowname_col optional column to promote to rownames if missing.
#' @return S4Vectors::DataFrame with rownames = gene IDs and columns like
#'   log2FC, pvalue, padj (if present), baseMean/AveExpr mapped to baseMean, etc.
#' @keywords internal
.tidy_de_results <- function(tbl, rowname_col = NULL) {
  # accept S4Vectors::DataFrame or base/tibble
  if (inherits(tbl, "DataFrame")) {
    df <- as.data.frame(tbl, stringsAsFactors = FALSE)
  } else if (is.data.frame(tbl)) {
    df <- as.data.frame(tbl, stringsAsFactors = FALSE)
  } else {
    cli::cli_abort("{.arg tbl} must be a data.frame/tibble or S4Vectors::DataFrame.")
  }

  # --- resolve gene_id (KEEP as an explicit column; do NOT delete) -------------
  # priority: user-specified -> existing 'gene_id' -> common id-like columns -> rownames
  id_candidates <- c("gene","geneID","GeneID","symbol","SYMBOL","ensembl","ensembl_id","ID","id")
  id_col <- NULL
  if (!is.null(rowname_col) && rowname_col %in% names(df)) {
    id_col <- rowname_col
  } else if ("gene_id" %in% names(df)) {
    id_col <- "gene_id"
  } else if (length(intersect(id_candidates, names(df))) > 0) {
    id_col <- intersect(id_candidates, names(df))[1]
  } else {
    rn <- rownames(df)
    if (!is.null(rn) && length(rn) == nrow(df) &&
        all(nzchar(rn)) && !anyDuplicated(rn) &&
        !all(rn == as.character(seq_len(nrow(df))))) {
      df$gene_id <- rn
      id_col <- "gene_id"
    }
  }
  if (is.null(id_col)) {
    cli::cli_abort(
      "Gene identifiers not found. Supply {.arg rowname_col} or include a {.field gene_id}/{.field symbol} column."
    )
  }
  # rename to canonical 'gene_id' if needed
  if (id_col != "gene_id") names(df)[names(df) == id_col] <- "gene_id"
  df$gene_id <- as.character(df$gene_id)

  # ensure valid IDs
  if (anyNA(df$gene_id) || any(df$gene_id == "") || anyDuplicated(df$gene_id)) {
    cli::cli_abort("Column {.field gene_id} must be unique, non-empty, and non-NA.")
  }

  # (optional) set rownames too, but KEEP the gene_id column
  rownames(df) <- df$gene_id

  # --- standardize column names ------------------------------------------------
  rename_one <- function(x, candidates, to) {
    hit <- intersect(candidates, colnames(x))
    if (length(hit) >= 1) colnames(x)[match(hit[1], colnames(x))] <- to
    x
  }
  # use lower-case 'log2fc' to match your plotting code
  df <- rename_one(df, c("log2fc","log2FoldChange","logFC","LFC","lfc"), "log2fc")
  df <- rename_one(df, c("pvalue","PValue","P.Value","p_val","p"),          "pvalue")
  df <- rename_one(df, c("padj","FDR","qvalue","svalue"),                    "padj")
  df <- rename_one(df, c("baseMean","AveExpr"),                              "baseMean")
  df <- rename_one(df, c("stat","LR","t"),                                   "stat")

  # --- robust numeric coercion -------------------------------------------------
  to_num <- function(v) {
    if (is.numeric(v)) return(v)
    v <- gsub(",", "", as.character(v), fixed = TRUE)
    suppressWarnings(as.numeric(v))
  }
  num_like <- intersect(c("baseMean","log2fc","lfcSE","stat","pvalue","padj","logCPM"), names(df))
  for (nm in num_like) df[[nm]] <- to_num(df[[nm]])
  # p-values: make non-finite NA -> 1
  if ("pvalue" %in% names(df)) df$pvalue[!is.finite(df$pvalue) | is.na(df$pvalue)] <- 1

  # --- canonical column order (gene_id first) ---------------------------------
  canon <- c("gene_id","baseMean","log2fc","lfcSE","stat","pvalue","padj","regulation")
  df <- df[, c(intersect(canon, names(df)), setdiff(names(df), canon)), drop = FALSE]

  # Return a plain data.frame (keeps EnhancedVolcano happy). Wrap with DataFrame() where needed.
  df
}





#' Align a DE table to reference gene order; add NA rows for missing genes.
#' Keeps a real `gene_id` column and returns a base data.frame.
#' Canonical column order: gene_id, baseMean, log2fc, lfcSE, stat, pvalue, padj, regulation, ...
#' @keywords internal
.align_de_to_counts <- function(df,
                                ref_rn,
                                id_col = NULL,
                                strict = FALSE,
                                warn_missing = TRUE) {
  # ---- coerce inputs ---------------------------------------------------------
  if (inherits(df, "DataFrame")) df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (!is.data.frame(df)) {
    cli::cli_abort("{.arg df} must be a data.frame/tibble or S4Vectors::DataFrame.")
  }
  # This function aligns by rownames, but a tibble silently ignores
  # `rownames<-`, so every gene would look absent and be replaced by an NA row:
  # correct labels, no data. Demote to a plain data.frame before touching them.
  if (!identical(class(df), "data.frame")) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)
  }
  ref_rn <- as.character(ref_rn)

  # The reference identifiers become the returned rownames. If they repeat,
  # `df[ref_rn, ]` silently de-duplicates them ("g2" -> "g2.1"), leaving the
  # rownames disagreeing with the gene_id column they were built from -- the
  # same label/value split that produced the 1.0.0 defect.
  if (!length(ref_rn)) {
    cli::cli_abort("{.arg ref_rn} must contain at least one gene identifier.")
  }
  if (anyNA(ref_rn) || any(!nzchar(ref_rn))) {
    cli::cli_abort("{.arg ref_rn} must not contain NA or empty gene identifiers.")
  }
  if (anyDuplicated(ref_rn)) {
    dups <- unique(ref_rn[duplicated(ref_rn)])
    cli::cli_abort(c(
      "{.arg ref_rn} must be unique.",
      "x" = "Duplicated: {.val {utils::head(dups, 5)}}",
      "i" = "Aligning to repeated identifiers would produce rownames that disagree with the {.field gene_id} column."
    ))
  }

  # ---- resolve / create canonical gene_id column -----------------------------
  pick <- NULL
  if (!is.null(id_col) && id_col %in% names(df)) {
    pick <- id_col
  } else if ("gene_id" %in% names(df)) {
    pick <- "gene_id"
  } else {
    # try rownames
    rn <- rownames(df)
    if (!is.null(rn) && length(rn) == nrow(df) &&
        all(nzchar(rn)) && !anyDuplicated(rn) &&
        !all(rn == as.character(seq_len(nrow(df))))) {
      df$gene_id <- rn
      pick <- "gene_id"
    } else {
      # as last resort: first column if clean character IDs
      fc <- names(df)[1]
      if (is.character(df[[fc]]) &&
          all(nzchar(df[[fc]])) && !anyDuplicated(df[[fc]]) && !any(is.na(df[[fc]]))) {
        pick <- fc
      }
    }
  }
  if (is.null(pick)) {
    cli::cli_abort(
      "Could not determine gene identifiers. Provide {.arg id_col} or include a {.field gene_id} column."
    )
  }
  if (pick != "gene_id") names(df)[names(df) == pick] <- "gene_id"
  df$gene_id <- as.character(df$gene_id)
  if (anyNA(df$gene_id) || any(df$gene_id == "") || anyDuplicated(df$gene_id)) {
    cli::cli_abort("{.field gene_id} must be unique, non-empty, and non-NA.")
  }
  # set rownames but KEEP gene_id column
  rownames(df) <- df$gene_id

  # ---- standardize column names (lowercase 'log2fc' to match plotting code) --
  rename_one <- function(x, candidates, to) {
    hit <- intersect(candidates, names(x))
    if (length(hit) >= 1) names(x)[match(hit[1], names(x))] <- to
    x
  }
  df <- rename_one(df, c("log2fc","log2FoldChange","logFC","LFC","lfc"), "log2fc")
  df <- rename_one(df, c("pvalue","PValue","P.Value","p_val","p"),       "pvalue")
  df <- rename_one(df, c("padj","FDR","qvalue","svalue"),                 "padj")
  df <- rename_one(df, c("baseMean","AveExpr"),                           "baseMean")
  df <- rename_one(df, c("lfcSE","SE","SE_lfc","lfc_se"),                 "lfcSE")
  df <- rename_one(df, c("stat","LR","t","W","F","z"),                    "stat")

  # ---- canonical column order ------------------------------------------------
  canon <- c("gene_id","baseMean","log2fc","lfcSE","stat","pvalue","padj","regulation")
  ordered_cols <- c(intersect(canon, names(df)), setdiff(names(df), canon))
  df <- df[, ordered_cols, drop = FALSE]

  # ---- add NA rows for genes missing from DE table ---------------------------
  missing_ids <- setdiff(ref_rn, rownames(df))
  if (length(missing_ids)) {
    if (strict) {
      cli::cli_abort("Missing {length(missing_ids)} genes from DE table; set {.code strict = FALSE} to fill with NA.")
    }
    if (warn_missing) {
      cli::cli_warn("Added {length(missing_ids)} NA rows to DE table to match reference genes.")
    }
    make_na_col <- function(x, n) {
      if (is.factor(x)) {
        factor(rep(NA, n), levels = levels(x))
      } else if (inherits(x, "Date")) {
        as.Date(rep(NA, n))
      } else if (inherits(x, "POSIXct")) {
        as.POSIXct(rep(NA, n), origin = "1970-01-01", tz = attr(x, "tzone") %||% "UTC")
      } else if (is.character(x)) {
        rep(NA_character_, n)
      } else if (is.integer(x)) {
        rep(NA_integer_, n)
      } else if (is.logical(x)) {
        rep(NA, n)
      } else {
        rep(NA_real_, n)  # numeric/double and fallback
      }
    }
    na_block <- as.data.frame(lapply(df, make_na_col, n = length(missing_ids)), stringsAsFactors = FALSE)
    rownames(na_block) <- missing_ids
    # ensure gene_id column is filled with rownames for NA rows
    if ("gene_id" %in% names(na_block)) na_block$gene_id <- rownames(na_block)
    df <- rbind(df, na_block)
  }

  # ---- final alignment to ref order (drops extra genes) ----------------------
  df <- df[ref_rn, , drop = FALSE]

  # Post-condition. This function exists to guarantee that a DE table lines up
  # with the counts matrix, so assert it rather than trusting the steps above:
  # every downstream plot and summary reads these two in parallel.
  if (!identical(rownames(df), ref_rn) ||
      !identical(as.character(df$gene_id), ref_rn)) {
    cli::cli_abort(c(
      "Internal error: aligned DE table does not match the reference gene order.",
      "i" = "Please report this at {.url https://github.com/cparsania/VISTA/issues}."
    ))
  }

  df
}

#' @keywords internal
.safe_numeric <- function(x, default = NA_real_) {
  out <- suppressWarnings(as.numeric(x))
  out[!is.finite(out)] <- default
  out
}

#' @keywords internal
.p_col_for_cutoff <- function(df, p_value_type = c("padj", "pvalue")) {
  p_value_type <- match.arg(p_value_type)
  if (p_value_type %in% names(df)) return(p_value_type)
  if (p_value_type == "padj" && "pvalue" %in% names(df)) return("pvalue")
  if (p_value_type == "pvalue" && "padj" %in% names(df)) return("padj")
  NA_character_
}

#' @keywords internal
.build_consensus_results <- function(deseq_comparisons,
                                     edger_comparisons,
                                     ref_rn,
                                     log2fc_cutoff = 1,
                                     pval_cutoff = 0.05,
                                     p_value_type = c("padj", "pvalue"),
                                     consensus_mode = c("intersection", "union"),
                                     consensus_log2fc = c("mean", "deseq2", "edger")) {
  p_value_type <- match.arg(p_value_type)
  consensus_mode <- match.arg(consensus_mode)
  consensus_log2fc <- match.arg(consensus_log2fc)

  common_comparisons <- intersect(names(deseq_comparisons), names(edger_comparisons))
  if (!length(common_comparisons)) {
    cli::cli_abort("No overlapping comparison names between DESeq2 and edgeR results for consensus building.")
  }

  consensus_comparisons <- list()
  consensus_summary <- list()

  for (comp_name in common_comparisons) {
    d <- .align_de_to_counts(
      df = .tidy_de_results(deseq_comparisons[[comp_name]], rowname_col = "gene_id"),
      ref_rn = ref_rn,
      warn_missing = FALSE
    )
    e <- .align_de_to_counts(
      df = .tidy_de_results(edger_comparisons[[comp_name]], rowname_col = "gene_id"),
      ref_rn = ref_rn,
      warn_missing = FALSE
    )

    d_log2fc <- .safe_numeric(d$log2fc, default = NA_real_)
    e_log2fc <- .safe_numeric(e$log2fc, default = NA_real_)
    d_pcol <- .p_col_for_cutoff(d, p_value_type = p_value_type)
    e_pcol <- .p_col_for_cutoff(e, p_value_type = p_value_type)

    d_sig <- if (!is.na(d_pcol)) {
      pvals <- .safe_numeric(d[[d_pcol]], default = 1)
      !is.na(d_log2fc) & abs(d_log2fc) >= log2fc_cutoff & pvals <= pval_cutoff
    } else {
      rep(FALSE, length(ref_rn))
    }
    e_sig <- if (!is.na(e_pcol)) {
      pvals <- .safe_numeric(e[[e_pcol]], default = 1)
      !is.na(e_log2fc) & abs(e_log2fc) >= log2fc_cutoff & pvals <= pval_cutoff
    } else {
      rep(FALSE, length(ref_rn))
    }

    same_direction <- sign(d_log2fc) == sign(e_log2fc) &
      sign(d_log2fc) != 0 &
      !is.na(d_log2fc) &
      !is.na(e_log2fc)

    support <- ifelse(
      d_sig & e_sig & same_direction, "both",
      ifelse(
        d_sig & !e_sig, "deseq2_only",
        ifelse(
          !d_sig & e_sig, "edger_only",
          ifelse(d_sig & e_sig & !same_direction, "discordant", "none")
        )
      )
    )

    consensus_significant <- if (consensus_mode == "intersection") {
      support == "both"
    } else {
      support %in% c("both", "deseq2_only", "edger_only")
    }

    log2fc_consensus <- switch(
      consensus_log2fc,
      mean = rowMeans(cbind(d_log2fc, e_log2fc), na.rm = TRUE),
      deseq2 = d_log2fc,
      edger = e_log2fc
    )
    # For method-specific calls in union mode, use the contributing method FC.
    log2fc_consensus[support == "deseq2_only"] <- d_log2fc[support == "deseq2_only"]
    log2fc_consensus[support == "edger_only"] <- e_log2fc[support == "edger_only"]
    log2fc_consensus[!is.finite(log2fc_consensus)] <- NA_real_

    d_pvalue <- if ("pvalue" %in% names(d)) .safe_numeric(d$pvalue, default = 1) else rep(1, length(ref_rn))
    e_pvalue <- if ("pvalue" %in% names(e)) .safe_numeric(e$pvalue, default = 1) else rep(1, length(ref_rn))
    d_padj <- if ("padj" %in% names(d)) .safe_numeric(d$padj, default = 1) else d_pvalue
    e_padj <- if ("padj" %in% names(e)) .safe_numeric(e$padj, default = 1) else e_pvalue

    # Every gene carries the less-significant (more conservative) of the two
    # backends, so the consensus p-value columns stay continuous and usable for
    # volcano/MA plots and ranking. Genes called by a single backend then take
    # that backend's value below.
    pvalue_consensus <- pmax(d_pvalue, e_pvalue, na.rm = TRUE)
    padj_consensus <- pmax(d_padj, e_padj, na.rm = TRUE)

    d_only_idx <- support == "deseq2_only"
    e_only_idx <- support == "edger_only"

    pvalue_consensus[d_only_idx] <- d_pvalue[d_only_idx]
    padj_consensus[d_only_idx] <- d_padj[d_only_idx]
    pvalue_consensus[e_only_idx] <- e_pvalue[e_only_idx]
    padj_consensus[e_only_idx] <- e_padj[e_only_idx]

    regulation <- rep("Other", length(ref_rn))
    up_idx <- consensus_significant & !is.na(log2fc_consensus) & log2fc_consensus >= log2fc_cutoff
    down_idx <- consensus_significant & !is.na(log2fc_consensus) & log2fc_consensus <= -log2fc_cutoff
    regulation[up_idx] <- "Up"
    regulation[down_idx] <- "Down"

    base_mean <- rowMeans(
      cbind(
        if ("baseMean" %in% names(d)) .safe_numeric(d$baseMean, default = NA_real_) else rep(NA_real_, length(ref_rn)),
        if ("baseMean" %in% names(e)) .safe_numeric(e$baseMean, default = NA_real_) else rep(NA_real_, length(ref_rn))
      ),
      na.rm = TRUE
    )
    base_mean[!is.finite(base_mean)] <- NA_real_

    out <- data.frame(
      gene_id = ref_rn,
      baseMean = base_mean,
      log2fc = log2fc_consensus,
      pvalue = pvalue_consensus,
      padj = padj_consensus,
      regulation = regulation,
      support = support,
      log2fc_deseq2 = d_log2fc,
      log2fc_edger = e_log2fc,
      pvalue_deseq2 = d_pvalue,
      pvalue_edger = e_pvalue,
      padj_deseq2 = d_padj,
      padj_edger = e_padj,
      stringsAsFactors = FALSE,
      row.names = ref_rn
    )

    consensus_comparisons[[comp_name]] <- out
    consensus_summary[[comp_name]] <- out %>% dplyr::count(regulation)
  }

  list(comparisons = consensus_comparisons, deg_summary = consensus_summary)
}






# Map gene IDs to display labels using an orgdb or supplied map
.map_gene_ids <- function(ids,
                          from_type = NULL,
                          to_type = "SYMBOL",
                          orgdb = NULL,
                          id_map = NULL) {
  ids <- as.character(ids)
  # `to_type == from_type` returns NA when from_type is NA and logical(0) when
  # it is NULL, and `||` can evaluate neither -- a display_id supplied without a
  # display_from failed with "missing value where TRUE/FALSE needed" instead of
  # saying what was missing. identical() is safe for both.
  if (is.null(to_type) || !length(to_type) || anyNA(to_type)) return(ids)
  if (any(to_type %in% c("", "none"))) return(ids)
  if (identical(as.character(to_type), as.character(from_type))) return(ids)

  if (!is.null(id_map)) {
    hit <- id_map[ids]
    return(ifelse(!is.na(hit) & nzchar(hit), hit, ids))
  }

  if (is.null(orgdb)) {
    cli::cli_abort("ID mapping requested but {.arg display_orgdb} was not supplied.")
  }
  if (is.null(from_type)) from_type <- to_type

  kts <- tryCatch(AnnotationDbi::keytypes(orgdb), error = function(e) NULL)
  if (is.null(kts)) cli::cli_abort("Could not retrieve key types from provided {.arg display_orgdb}.")
  if (!from_type %in% kts) cli::cli_abort("Unknown {.arg display_from} {.val {from_type}}. Available keytypes: {.val {kts}}")
  if (!to_type %in% kts) cli::cli_abort("Unknown {.arg display_id} {.val {to_type}}. Available keytypes: {.val {kts}}")

  res <- tryCatch(
    AnnotationDbi::select(orgdb, keys = ids, keytype = from_type, columns = to_type),
    error = function(e) NULL
  )
  if (is.null(res) || !nrow(res)) return(ids)
  # keep first mapping per input id
  res <- res[!duplicated(res[[from_type]]), , drop = FALSE]
  res <- res[!is.na(res[[to_type]]), , drop = FALSE]
  if (!nrow(res)) return(ids)
  map <- stats::setNames(res[[to_type]], res[[from_type]])
  out <- map[ids]
  ifelse(!is.na(out) & nzchar(out), out, ids)
}

#' Set or append rowData annotations on a VISTA object
#'
#' Accepts a data.frame/tibble/DataFrame of gene-level annotations, aligns it to
#' the VISTA row order, and stores it in `rowData(x)`. Rows are matched by a
#' key column (default: tries `gene_id` or rownames); Ensembl version suffixes
#' can be stripped for matching.
#'
#' @param x A `VISTA` object.
#' @param annotations Optional data.frame/tibble/DataFrame with one row per gene
#'   and a column containing the gene IDs to match against `rownames(x)`. If
#'   omitted, annotations are pulled from `orgdb`.
#' @param orgdb Optional `OrgDb` object; when supplied (and `annotations` is
#'   `NULL`), annotations are retrieved via `AnnotationDbi::select()`.
#' @param key_col Name of the column in `annotations` that holds the gene IDs.
#'   If `NULL`, the function will try `gene_id`, `gene`, `ENSEMBL`, `SYMBOL`, or
#'   use `rownames(annotations)`. Ignored when `annotations` is `NULL` and
#'   `orgdb` is used.
#' @param keytype Key type for `orgdb` lookups (e.g., `"ENSEMBL"`, `"SYMBOL"`).
#'   If `NULL`, inferred from `rownames(x)` (`ENSEMBL` if they start with "ENS",
#'   otherwise `SYMBOL`).
#' @param columns Character vector of OrgDb columns to retrieve when using
#'   `orgdb`. Default:
#'   `c("SYMBOL","GENENAME","ENSEMBL","ENTREZID","TXCHROM","TXSTART","TXEND")`.
#'   The `TXCHROM`/`TXSTART`/`TXEND` fields carry basic genomic coordinates when
#'   available in the OrgDb.
#' @param drop_version Logical; if `TRUE`, strips Ensembl version suffixes (e.g.,
#'   `.1`) from both the VISTA rownames and the key column/keys before matching.
#' @param overwrite Logical; if `TRUE`, replaces existing `rowData`. If `FALSE`,
#'   new columns are appended (overwriting by name when names collide).
#'
#' @details OrgDb packages rarely include full genomic coordinates; the default
#'   `TXCHROM/TXSTART/TXEND` columns may therefore be `NA` unless your OrgDb
#'   provides them. For reliable coordinates, fetch them from an `EnsDb`/`TxDb`
#'   (via `genes()` or `biomaRt`/AnnotationHub), build an annotation table keyed
#'   on your gene IDs, and supply that via the `annotations` argument. When
#'   fetching from an `OrgDb`, only columns available in that database will be
#'   filled.
#'
#' @return The updated `VISTA` object with `rowData` populated/appended.
#'
#' @examples
#' vista <- example_vista()
#' custom_annot <- data.frame(
#'   gene_id = rownames(vista)[seq_len(10)],
#'   custom_info = paste0("Info_", seq_len(10))
#' )
#' vista2 <- set_rowdata(vista, annotations = custom_annot, key_col = "gene_id")
#' head(SummarizedExperiment::rowData(vista2)$custom_info)
#'
#' \donttest{
#' # Load example VISTA object
#' data("count_data", package = "VISTA")
#' data("sample_metadata", package = "VISTA")
#'
#' vista <- create_vista(
#'   counts = count_data[seq_len(100), ],
#'   sample_info = sample_metadata[seq_len(6), ],
#'   column_geneid = "gene_id",
#'   group_column = "cond_long",
#'   group_numerator = "treatment1",
#'   group_denominator = "control"
#' )
#'
#' # Add annotations from OrgDb (human)
#' if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   vista <- set_rowdata(
#'     vista,
#'     orgdb = org.Hs.eg.db::org.Hs.eg.db,
#'     columns = c("SYMBOL", "GENENAME", "ENTREZID")
#'   )
#'
#'   # View updated rowData
#'   head(SummarizedExperiment::rowData(vista))
#' }
#'
#' # Or provide custom annotations
#' custom_annot <- data.frame(
#'   gene_id = rownames(vista)[seq_len(10)],
#'   custom_info = paste0("Info_", seq_len(10))
#' )
#' vista <- set_rowdata(vista, annotations = custom_annot, key_col = "gene_id")
#' }
#'
#' @export
set_rowdata <- function(x,
                        annotations = NULL,
                        orgdb = NULL,
                        key_col = NULL,
                        keytype = NULL,
                        columns = c("SYMBOL", "GENENAME", "ENSEMBL", "ENTREZID", "TXCHROM", "TXSTART", "TXEND"),
                        drop_version = TRUE,
                        overwrite = FALSE) {
  stopifnot(inherits(x, "VISTA"))
  if (is.null(annotations) && is.null(orgdb)) {
    cli::cli_abort("Provide either {.arg annotations} or {.arg orgdb} for lookup.")
  }

  rn <- rownames(x)
  strip_ver <- function(z) sub("\\..*$", "", as.character(z))
  rn_keys <- if (drop_version) strip_ver(rn) else rn

  if (is.null(annotations)) {
    if (is.null(keytype)) {
      keytype <- if (length(rn_keys) && grepl("^ENS", rn_keys[1])) "ENSEMBL" else "SYMBOL"
    }
    ann_df <- AnnotationDbi::select(
      orgdb,
      keys = unique(rn_keys),
      keytype = keytype,
      columns = columns
    )
    if (is.null(ann_df) || !nrow(ann_df)) cli::cli_abort("No annotations returned from {.arg orgdb}.")
    # keep first mapping per key
    ann_df <- ann_df[!duplicated(ann_df[[keytype]]), , drop = FALSE]
    key_col <- keytype
  } else {
    ann_df <- as.data.frame(annotations, stringsAsFactors = FALSE)
    if (is.null(key_col)) {
      candidates <- c("gene_id", "gene", "ENSEMBL", "SYMBOL")
      hit <- candidates[candidates %in% names(ann_df)][1]
      if (is.na(hit) && !is.null(rownames(ann_df)) && all(nzchar(rownames(ann_df)))) {
        ann_df$gene_id <- rownames(ann_df)
        hit <- "gene_id"
      }
      if (is.na(hit)) cli::cli_abort("Could not detect a gene ID column in {.arg annotations}; please supply {.arg key_col}.")
      key_col <- hit
    }
    if (!key_col %in% names(ann_df)) {
      cli::cli_abort("Column {.val {key_col}} not found in {.arg annotations}.")
    }
  }

  keys <- ann_df[[key_col]]
  if (drop_version) keys <- strip_ver(keys)
  ann_df[[key_col]] <- keys
  # keep first occurrence per key
  ann_df <- ann_df[!duplicated(ann_df[[key_col]]), , drop = FALSE]

  ord <- match(rn_keys, ann_df[[key_col]])
  if (anyNA(ord)) {
    missing <- rn[is.na(ord)]
    cli::cli_warn("Missing annotations for {length(missing)} genes; filling those rows with NA.")
  }

  ann_df <- ann_df[ord, , drop = FALSE]
  ann_df[[key_col]] <- NULL

  rd_new <- S4Vectors::DataFrame(ann_df)
  if (!overwrite && ncol(SummarizedExperiment::rowData(x)) > 0) {
    rd_existing <- SummarizedExperiment::rowData(x)
    dup_cols <- intersect(colnames(rd_existing), colnames(rd_new))
    if (length(dup_cols)) {
      rd_existing <- rd_existing[, setdiff(colnames(rd_existing), dup_cols), drop = FALSE]
    }
    SummarizedExperiment::rowData(x) <- cbind(rd_existing, rd_new)
  } else {
    SummarizedExperiment::rowData(x) <- rd_new
  }

  validate_vista(x, level = "core", error = TRUE)
  x
}
