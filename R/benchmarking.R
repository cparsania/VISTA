#' Benchmark VISTA against standalone differential-expression backends
#'
#' Run VISTA and direct DESeq2, edgeR, and limma pipelines with matched
#' preprocessing, contrast definitions, and DEG thresholds, then compare the
#' resulting tables, DEG calls, normalized matrices, and critical plot inputs.
#'
#' @param counts Raw counts (matrix/data.frame) with a gene-id column and sample columns.
#' @param sample_info Data frame with sample metadata.
#' @param column_geneid Column name in `counts` that contains gene identifiers.
#' @param group_column Column in `sample_info` used to group samples.
#' @param group_numerator Character vector of numerator groups for pairwise comparisons.
#' @param group_denominator Character vector of denominator groups.
#' @param methods Character vector of backends to benchmark. Any subset of
#'   `c("deseq2", "edger", "limma")`.
#' @param min_counts Minimum total counts per gene to retain.
#' @param min_replicates Minimum samples per group meeting filtering criteria.
#' @param log2fc_cutoff Absolute log2 fold-change threshold for DEG calling.
#' @param pval_cutoff P-value (or adjusted p-value) threshold.
#' @param p_value_type Either `"padj"` or `"pvalue"`.
#' @param covariates Optional character vector of additional sample_info columns.
#' @param design_formula Optional model formula (or formula string) including `group_column`.
#' @param tolerance Numeric tolerance used for floating-point comparisons.
#' @param return_plots Logical; if `TRUE`, return paired VISTA/reference plots
#'   for MA, volcano, DEG count, and PCA views.
#'
#' @return A list with fields `valid`, `comparison_summary`, `visual_summary`,
#'   and `methods`. Each element of `methods` contains the VISTA object, direct
#'   backend results, structural validation output, self-consistency checks, and
#'   optional plot objects.
#' @examples
#' \donttest{
#' data("count_data", package = "VISTA")
#' data("sample_metadata", package = "VISTA")
#'
#' target_groups <- c("control", "treatment1")
#' sample_subset <- sample_metadata[sample_metadata$cond_long %in% target_groups, ]
#' count_subset <- count_data[1:150, c("gene_id", sample_subset$sample_names)]
#'
#' bm <- benchmark_vista_equivalence(
#'   counts = count_subset,
#'   sample_info = sample_subset,
#'   column_geneid = "gene_id",
#'   group_column = "cond_long",
#'   group_numerator = "treatment1",
#'   group_denominator = "control",
#'   methods = c("deseq2", "edger"),
#'   min_counts = 5,
#'   min_replicates = 1
#' )
#'
#' bm$comparison_summary
#' bm$visual_summary
#' }
#' @export
benchmark_vista_equivalence <- function(counts,
                                        sample_info,
                                        column_geneid,
                                        group_column,
                                        group_numerator,
                                        group_denominator,
                                        methods = c("deseq2", "edger", "limma"),
                                        min_counts = 10,
                                        min_replicates = 1,
                                        log2fc_cutoff = 1,
                                        pval_cutoff = 0.05,
                                        p_value_type = "padj",
                                        covariates = NULL,
                                        design_formula = NULL,
                                        tolerance = 1e-8,
                                        return_plots = FALSE) {
  methods <- match.arg(methods, c("deseq2", "edger", "limma"), several.ok = TRUE)
  p_value_type <- match.arg(p_value_type, c("padj", "pvalue"))
  design_formula <- .coerce_design_formula(design_formula)
  sample_info <- .normalize_sample_info(
    sample_info = sample_info,
    counts = counts,
    column_geneid = column_geneid
  )
  covariates <- .validate_covariates(sample_info, group_column, covariates)

  method_results <- stats::setNames(vector("list", length(methods)), methods)
  comparison_summary <- vector("list", length(methods))
  visual_summary <- vector("list", length(methods))

  for (i in seq_along(methods)) {
    method <- methods[[i]]

    vista_obj <- create_vista(
      counts = counts,
      sample_info = sample_info,
      column_geneid = column_geneid,
      group_column = group_column,
      group_numerator = group_numerator,
      group_denominator = group_denominator,
      method = method,
      min_counts = min_counts,
      min_replicates = min_replicates,
      log2fc_cutoff = log2fc_cutoff,
      pval_cutoff = pval_cutoff,
      p_value_type = p_value_type,
      covariates = covariates,
      design_formula = design_formula,
      validate = TRUE
    )

    reference <- .run_reference_backend(
      method = method,
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
    )

    structural_validation <- validate_vista(vista_obj, level = "full", error = FALSE)
    self_consistency <- .validate_vista_self_consistency(
      x = vista_obj,
      p_value_type = p_value_type,
      log2fc_cutoff = log2fc_cutoff,
      pval_cutoff = pval_cutoff
    )

    method_check <- .compare_vista_and_reference(
      vista_obj = vista_obj,
      reference = reference,
      method = method,
      p_value_type = p_value_type,
      log2fc_cutoff = log2fc_cutoff,
      pval_cutoff = pval_cutoff,
      tolerance = tolerance,
      return_plots = return_plots
    )

    method_results[[method]] <- c(
      list(
        vista = vista_obj,
        reference = reference,
        structural_validation = structural_validation,
        self_consistency = self_consistency
      ),
      method_check
    )

    comparison_summary[[i]] <- method_check$comparison_summary |>
      dplyr::mutate(
        structural_valid = structural_validation$valid,
        self_consistency_valid = length(self_consistency$issues) == 0,
        all_checks_pass = .data$all_checks_pass &
          .data$structural_valid &
          .data$self_consistency_valid
      )

    visual_summary[[i]] <- method_check$visual_summary |>
      dplyr::mutate(
        structural_valid = structural_validation$valid,
        self_consistency_valid = length(self_consistency$issues) == 0,
        pass = .data$pass & .data$structural_valid & .data$self_consistency_valid
      )
  }

  comparison_summary <- dplyr::bind_rows(comparison_summary)
  visual_summary <- dplyr::bind_rows(visual_summary)
  valid <- all(comparison_summary$all_checks_pass) && all(visual_summary$pass)

  list(
    valid = valid,
    comparison_summary = comparison_summary,
    visual_summary = visual_summary,
    methods = method_results
  )
}

#' Deep validation of VISTA differential-expression fidelity
#'
#' This validator combines `validate_vista()` with backend-to-backend numerical
#' equivalence checks against standalone DESeq2, edgeR, and limma runs.
#'
#' @inheritParams benchmark_vista_equivalence
#' @param error Logical; if `TRUE`, abort when any discrepancy is detected.
#'
#' @return Invisibly returns the full benchmark report.
#' @examples
#' \donttest{
#' data("count_data", package = "VISTA")
#' data("sample_metadata", package = "VISTA")
#'
#' target_groups <- c("control", "treatment1")
#' sample_subset <- sample_metadata[sample_metadata$cond_long %in% target_groups, ]
#' count_subset <- count_data[1:150, c("gene_id", sample_subset$sample_names)]
#'
#' validate_vista_deep(
#'   counts = count_subset,
#'   sample_info = sample_subset,
#'   column_geneid = "gene_id",
#'   group_column = "cond_long",
#'   group_numerator = "treatment1",
#'   group_denominator = "control",
#'   methods = c("deseq2", "edger"),
#'   min_counts = 5,
#'   min_replicates = 1
#' )
#' }
#' @export
validate_vista_deep <- function(counts,
                                sample_info,
                                column_geneid,
                                group_column,
                                group_numerator,
                                group_denominator,
                                methods = c("deseq2", "edger", "limma"),
                                min_counts = 10,
                                min_replicates = 1,
                                log2fc_cutoff = 1,
                                pval_cutoff = 0.05,
                                p_value_type = "padj",
                                covariates = NULL,
                                design_formula = NULL,
                                tolerance = 1e-8,
                                return_plots = FALSE,
                                error = TRUE) {
  report <- benchmark_vista_equivalence(
    counts = counts,
    sample_info = sample_info,
    column_geneid = column_geneid,
    group_column = group_column,
    group_numerator = group_numerator,
    group_denominator = group_denominator,
    methods = methods,
    min_counts = min_counts,
    min_replicates = min_replicates,
    log2fc_cutoff = log2fc_cutoff,
    pval_cutoff = pval_cutoff,
    p_value_type = p_value_type,
    covariates = covariates,
    design_formula = design_formula,
    tolerance = tolerance,
    return_plots = return_plots
  )

  if (!report$valid) {
    issues <- c(
      .collect_benchmark_issues(report$comparison_summary, type = "comparison"),
      .collect_benchmark_issues(report$visual_summary, type = "visual"),
      unlist(lapply(report$methods, function(x) x$self_consistency$issues), use.names = FALSE),
      unlist(lapply(report$methods, function(x) x$structural_validation$issues), use.names = FALSE)
    )
    issues <- unique(stats::na.omit(issues))
    bullets <- stats::setNames(as.character(issues), rep("x", length(issues)))

    if (isTRUE(error)) {
      cli::cli_abort(c("Deep VISTA validation failed.", bullets))
    } else {
      cli::cli_warn(c("Deep VISTA validation reported issues.", bullets))
    }
  }

  invisible(report)
}

#' @keywords internal
.collect_benchmark_issues <- function(df, type = c("comparison", "visual")) {
  type <- match.arg(type)
  if (is.null(df) || !nrow(df)) {
    return(character())
  }

  if (type == "comparison") {
    failed <- df[!df$all_checks_pass, , drop = FALSE]
    if (!nrow(failed)) {
      return(character())
    }
    apply(failed, 1, function(row) {
      paste0(
        row[["method"]], " / ", row[["comparison"]], ": ",
        "regulation_identical=", row[["regulation_identical"]],
        ", deg_sets_identical=", row[["deg_sets_identical"]],
        ", norm_counts_identical=", row[["norm_counts_identical"]],
        ", baseMean_within_tolerance=", row[["baseMean_within_tolerance"]],
        ", log2fc_within_tolerance=", row[["log2fc_within_tolerance"]],
        ", pvalue_within_tolerance=", row[["pvalue_within_tolerance"]],
        ", padj_within_tolerance=", row[["padj_within_tolerance"]]
      )
    })
  } else {
    failed <- df[!df$pass, , drop = FALSE]
    if (!nrow(failed)) {
      return(character())
    }
    apply(failed, 1, function(row) {
      comp <- if (!is.na(row[["comparison"]]) && nzchar(row[["comparison"]])) {
        paste0(" / ", row[["comparison"]])
      } else {
        ""
      }
      paste0(
        row[["method"]], comp, " / ", row[["visual"]], ": ",
        row[["detail"]]
      )
    })
  }
}

#' @keywords internal
.run_reference_backend <- function(method,
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
                                   p_value_type = "padj") {
  switch(
    method,
    deseq2 = .run_reference_deseq2(
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
    edger = .run_reference_edger(
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
      p_value_type = if (identical(p_value_type, "padj")) "FDR" else "PValue"
    ),
    limma = .run_reference_limma(
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
      p_value_type = if (identical(p_value_type, "padj")) "FDR" else "PValue"
    )
  )
}

#' @keywords internal
.prepare_reference_inputs <- function(counts,
                                      sample_info,
                                      column_geneid,
                                      group_column,
                                      covariates = NULL,
                                      design_formula = NULL,
                                      min_counts = 10,
                                      backend = c("deseq2", "edger", "limma")) {
  backend <- match.arg(backend)

  count_data <- counts |>
    dplyr::select(tidyselect::all_of(c(column_geneid, sample_info$sample_names))) |>
    as.data.frame(stringsAsFactors = FALSE) |>
    tibble::as_tibble() |>
    tibble::column_to_rownames(column_geneid)
  count_data[is.na(count_data)] <- 0
  count_data <- count_data[rowSums(count_data) >= min_counts, , drop = FALSE]

  col_data <- sample_info |>
    as.data.frame(stringsAsFactors = FALSE) |>
    tibble::as_tibble() |>
    tibble::column_to_rownames("sample_names")
  col_data[[group_column]] <- factor(col_data[[group_column]])

  design_formula <- .resolve_design_formula(
    sample_info = col_data,
    group_column = group_column,
    covariates = covariates,
    design_formula = design_formula,
    backend = backend
  )

  design_vars <- intersect(all.vars(design_formula), colnames(col_data))
  for (vv in design_vars) {
    if (is.character(col_data[[vv]])) {
      col_data[[vv]] <- factor(col_data[[vv]])
    }
  }

  list(
    count_data = count_data,
    col_data = col_data,
    design_formula = design_formula
  )
}

#' @keywords internal
.run_reference_deseq2 <- function(counts,
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
                                  p_value_type = "padj") {
  inputs <- .prepare_reference_inputs(
    counts = counts,
    sample_info = sample_info,
    column_geneid = column_geneid,
    group_column = group_column,
    covariates = covariates,
    design_formula = design_formula,
    min_counts = min_counts,
    backend = "deseq2"
  )

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = as.matrix(inputs$count_data),
    colData = inputs$col_data,
    design = inputs$design_formula
  )
  dds <- dds[rowSums(DESeq2::counts(dds) >= min_counts) >= min_replicates, ]
  dds <- DESeq2::DESeq(dds)
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)

  comparisons <- list()
  deg_summary <- list()
  for (i in seq_along(group_numerator)) {
    num <- group_numerator[[i]]
    den <- group_denominator[[i]]
    if (identical(num, den)) {
      next
    }

    res <- DESeq2::results(dds, contrast = c(group_column, num, den))
    categorized <- .categorize_deg_results(
      de_results = res,
      log2fc_cutoff = log2fc_cutoff,
      pval_cutoff = pval_cutoff,
      p_value_type = p_value_type
    )
    comp_name <- paste(num, den, sep = "_VS_")
    comparisons[[comp_name]] <- categorized
    deg_summary[[comp_name]] <- categorized |>
      dplyr::count(.data$regulation)
  }

  list(
    norm_counts = norm_counts,
    sample_info = as.data.frame(SummarizedExperiment::colData(dds), stringsAsFactors = FALSE),
    row_data = S4Vectors::DataFrame(
      baseMean = rowMeans(norm_counts),
      row.names = rownames(norm_counts)
    ),
    comparisons = comparisons,
    deg_summary = deg_summary
  )
}

#' @keywords internal
.run_reference_edger <- function(counts,
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
                                 p_value_type = "FDR") {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg edgeR} is required for edgeR benchmarking.")
  }

  inputs <- .prepare_reference_inputs(
    counts = counts,
    sample_info = sample_info,
    column_geneid = column_geneid,
    group_column = group_column,
    covariates = covariates,
    design_formula = design_formula,
    min_counts = min_counts,
    backend = "edger"
  )

  dge <- edgeR::DGEList(counts = inputs$count_data, group = inputs$col_data[[group_column]])
  dge <- edgeR::calcNormFactors(dge)
  keep <- rowSums(edgeR::cpm(dge) > 1) >= min_replicates
  dge <- dge[keep, , keep.lib.sizes = FALSE]

  design <- stats::model.matrix(inputs$design_formula, data = inputs$col_data)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmFit(dge, design)
  norm_counts <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)

  comparisons <- list()
  deg_summary <- list()
  for (i in seq_along(group_numerator)) {
    num <- group_numerator[[i]]
    den <- group_denominator[[i]]
    if (identical(num, den)) {
      next
    }

    contrast_vector <- .edgeR_group_contrast(
      design_colnames = colnames(design),
      group_column = group_column,
      numerator = num,
      denominator = den
    )
    lrt <- edgeR::glmLRT(fit, contrast = contrast_vector)
    res <- edgeR::topTags(lrt, n = Inf)$table
    res$gene_id <- rownames(res)

    categorized <- .categorize_deg_results(
      de_results = res,
      log2fc_cutoff = log2fc_cutoff,
      pval_cutoff = pval_cutoff,
      p_value_type = p_value_type
    )
    comp_name <- paste(num, den, sep = "_VS_")
    comparisons[[comp_name]] <- categorized
    deg_summary[[comp_name]] <- categorized |>
      dplyr::count(.data$regulation)
  }

  list(
    norm_counts = norm_counts,
    sample_info = inputs$col_data,
    row_data = S4Vectors::DataFrame(
      baseMean = rowMeans(norm_counts),
      row.names = rownames(norm_counts)
    ),
    comparisons = comparisons,
    deg_summary = deg_summary
  )
}

#' @keywords internal
.run_reference_limma <- function(counts,
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
                                 p_value_type = "FDR") {
  if (!requireNamespace("limma", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg limma} is required for limma benchmarking.")
  }
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg edgeR} is required for limma benchmarking.")
  }

  inputs <- .prepare_reference_inputs(
    counts = counts,
    sample_info = sample_info,
    column_geneid = column_geneid,
    group_column = group_column,
    covariates = covariates,
    design_formula = design_formula,
    min_counts = min_counts,
    backend = "limma"
  )

  dge <- edgeR::DGEList(counts = inputs$count_data)
  dge <- edgeR::calcNormFactors(dge)
  keep <- rowSums(edgeR::cpm(dge) > 1) >= min_replicates
  dge <- dge[keep, , keep.lib.sizes = FALSE]

  design <- stats::model.matrix(inputs$design_formula, data = inputs$col_data)
  v <- limma::voom(dge, design = design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  norm_counts <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)

  comparisons <- list()
  deg_summary <- list()
  for (i in seq_along(group_numerator)) {
    num <- group_numerator[[i]]
    den <- group_denominator[[i]]
    if (identical(num, den)) {
      next
    }

    contrast_vector <- .edgeR_group_contrast(
      design_colnames = colnames(design),
      group_column = group_column,
      numerator = num,
      denominator = den
    )
    contrast_matrix <- matrix(
      contrast_vector,
      ncol = 1,
      dimnames = list(colnames(design), "contrast")
    )
    fit2 <- limma::contrasts.fit(fit, contrasts = contrast_matrix)
    fit2 <- limma::eBayes(fit2)
    res <- limma::topTable(fit2, coef = 1, number = Inf, sort.by = "none")
    res$gene_id <- rownames(res)
    res$PValue <- res$P.Value
    res$FDR <- res$adj.P.Val

    categorized <- .categorize_deg_results(
      de_results = res,
      log2fc_cutoff = log2fc_cutoff,
      pval_cutoff = pval_cutoff,
      p_value_type = p_value_type
    )
    comp_name <- paste(num, den, sep = "_VS_")
    comparisons[[comp_name]] <- categorized
    deg_summary[[comp_name]] <- categorized |>
      dplyr::count(.data$regulation)
  }

  list(
    norm_counts = norm_counts,
    sample_info = inputs$col_data,
    row_data = S4Vectors::DataFrame(
      baseMean = rowMeans(norm_counts),
      row.names = rownames(norm_counts)
    ),
    comparisons = comparisons,
    deg_summary = deg_summary
  )
}

#' @keywords internal
.compare_vista_and_reference <- function(vista_obj,
                                         reference,
                                         method,
                                         p_value_type = "padj",
                                         log2fc_cutoff = 1,
                                         pval_cutoff = 0.05,
                                         tolerance = 1e-8,
                                         return_plots = FALSE) {
  vista_norm <- as.matrix(norm_counts(vista_obj))
  ref_norm <- as.matrix(reference$norm_counts)
  ref_norm <- ref_norm[rownames(vista_norm), colnames(vista_norm), drop = FALSE]
  norm_check <- .compare_numeric_matrix(vista_norm, ref_norm, tolerance = tolerance)

  vista_comps <- comparisons(vista_obj)
  ref_comps <- reference$comparisons
  comp_names <- intersect(names(vista_comps), names(ref_comps))
  ref_rn <- rownames(vista_norm)

  comparison_rows <- vector("list", length(comp_names))
  visual_rows <- vector("list", length(comp_names) + 2L)
  plot_bundle <- if (isTRUE(return_plots)) list(ma = list(), volcano = list()) else NULL

  for (i in seq_along(comp_names)) {
    comp_name <- comp_names[[i]]
    vista_tbl <- .align_de_to_counts(
      df = .tidy_de_results(vista_comps[[comp_name]], rowname_col = "gene_id"),
      ref_rn = ref_rn,
      warn_missing = FALSE
    )
    ref_tbl <- .align_de_to_counts(
      df = .tidy_de_results(ref_comps[[comp_name]], rowname_col = "gene_id"),
      ref_rn = ref_rn,
      warn_missing = FALSE
    )

    comp_check <- .compare_de_tables(
      vista_tbl = vista_tbl,
      ref_tbl = ref_tbl,
      method = method,
      comparison = comp_name,
      norm_check = norm_check,
      tolerance = tolerance
    )
    comparison_rows[[i]] <- comp_check$summary

    visual_rows[[i]] <- .compare_visual_payload(
      vista_obj = vista_obj,
      reference = reference,
      method = method,
      comparison = comp_name,
      visual = "ma",
      vista_df = .ma_payload_from_vista(vista_obj, comp_name),
      ref_df = .ma_payload_from_reference(reference, comp_name),
      tolerance = tolerance
    )

    visual_rows[[length(comp_names) + i]] <- .compare_visual_payload(
      vista_obj = vista_obj,
      reference = reference,
      method = method,
      comparison = comp_name,
      visual = "volcano",
      vista_df = .volcano_payload_from_vista(
        vista_obj,
        comp_name,
        p_value_type = p_value_type,
        log2fc_cutoff = log2fc_cutoff,
        pval_cutoff = pval_cutoff
      ),
      ref_df = .volcano_payload_from_reference(
        reference,
        comp_name,
        p_value_type = p_value_type,
        log2fc_cutoff = log2fc_cutoff,
        pval_cutoff = pval_cutoff
      ),
      tolerance = tolerance
    )

    if (isTRUE(return_plots)) {
      plot_bundle$ma[[comp_name]] <- list(
        vista = get_ma_plot(vista_obj, sample_comparison = comp_name),
        reference = .plot_reference_ma(.ma_payload_from_reference(reference, comp_name))
      )
      plot_bundle$volcano[[comp_name]] <- list(
        vista = .plot_reference_volcano(
          .volcano_payload_from_vista(
            vista_obj,
            comp_name,
            p_value_type = p_value_type,
            log2fc_cutoff = log2fc_cutoff,
            pval_cutoff = pval_cutoff
          )
        ),
        reference = .plot_reference_volcano(
          .volcano_payload_from_reference(
            reference,
            comp_name,
            p_value_type = p_value_type,
            log2fc_cutoff = log2fc_cutoff,
            pval_cutoff = pval_cutoff
          )
        )
      )
    }
  }

  visual_rows[[length(comp_names) * 2 + 1L]] <- .compare_visual_payload(
    vista_obj = vista_obj,
    reference = reference,
    method = method,
    comparison = NA_character_,
    visual = "deg_count",
    vista_df = .deg_count_payload_from_vista(vista_obj),
    ref_df = .deg_count_payload_from_reference(reference),
    tolerance = tolerance
  )
  visual_rows[[length(comp_names) * 2 + 2L]] <- .compare_visual_payload(
    vista_obj = vista_obj,
    reference = reference,
    method = method,
    comparison = NA_character_,
    visual = "pca",
    vista_df = .pca_payload_from_vista(vista_obj),
    ref_df = .pca_payload_from_reference(reference, vista_obj),
    tolerance = tolerance,
    allow_sign_flip = TRUE
  )

  if (isTRUE(return_plots)) {
    plot_bundle$deg_count <- list(
      vista = get_deg_count_barplot(vista_obj, label = TRUE),
      reference = .plot_reference_deg_count(.deg_count_payload_from_reference(reference))
    )
    plot_bundle$pca <- list(
      vista = get_pca_plot(vista_obj),
      reference = .plot_reference_pca(.pca_payload_from_reference(reference, vista_obj), vista_obj)
    )
  }

  list(
    comparison_summary = dplyr::bind_rows(comparison_rows),
    visual_summary = dplyr::bind_rows(visual_rows),
    plots = plot_bundle
  )
}

#' @keywords internal
.compare_de_tables <- function(vista_tbl,
                               ref_tbl,
                               method,
                               comparison,
                               norm_check,
                               tolerance = 1e-8) {
  base_check <- .compare_numeric_vectors(vista_tbl$baseMean, ref_tbl$baseMean, tolerance = tolerance)
  log2fc_check <- .compare_numeric_vectors(vista_tbl$log2fc, ref_tbl$log2fc, tolerance = tolerance)
  pvalue_check <- .compare_numeric_vectors(vista_tbl$pvalue, ref_tbl$pvalue, tolerance = tolerance)
  padj_check <- .compare_numeric_vectors(vista_tbl$padj, ref_tbl$padj, tolerance = tolerance)

  vista_reg <- as.character(vista_tbl$regulation)
  ref_reg <- as.character(ref_tbl$regulation)
  regulation_identical <- identical(vista_reg, ref_reg)
  up_identical <- identical(
    vista_tbl$gene_id[vista_reg == "Up"],
    ref_tbl$gene_id[ref_reg == "Up"]
  )
  down_identical <- identical(
    vista_tbl$gene_id[vista_reg == "Down"],
    ref_tbl$gene_id[ref_reg == "Down"]
  )

  summary <- tibble::tibble(
    method = method,
    comparison = comparison,
    n_genes = nrow(vista_tbl),
    up_genes_identical = up_identical,
    down_genes_identical = down_identical,
    deg_sets_identical = up_identical && down_identical,
    regulation_identical = regulation_identical,
    norm_counts_identical = norm_check$identical,
    baseMean_within_tolerance = base_check$identical,
    log2fc_within_tolerance = log2fc_check$identical,
    pvalue_within_tolerance = pvalue_check$identical,
    padj_within_tolerance = padj_check$identical,
    max_abs_norm_counts_diff = norm_check$max_abs_diff,
    max_abs_baseMean_diff = base_check$max_abs_diff,
    max_abs_log2fc_diff = log2fc_check$max_abs_diff,
    max_abs_pvalue_diff = pvalue_check$max_abs_diff,
    max_abs_padj_diff = padj_check$max_abs_diff,
    all_checks_pass = regulation_identical &&
      up_identical &&
      down_identical &&
      norm_check$identical &&
      base_check$identical &&
      log2fc_check$identical &&
      pvalue_check$identical &&
      padj_check$identical
  )

  list(summary = summary)
}

#' @keywords internal
.compare_visual_payload <- function(vista_obj,
                                    reference,
                                    method,
                                    comparison,
                                    visual,
                                    vista_df,
                                    ref_df,
                                    tolerance = 1e-8,
                                    allow_sign_flip = FALSE) {
  if (allow_sign_flip) {
    ref_df <- .align_pca_signs(vista_df, ref_df)
  }

  shared_cols <- intersect(names(vista_df), names(ref_df))
  vista_df <- vista_df[, shared_cols, drop = FALSE]
  ref_df <- ref_df[, shared_cols, drop = FALSE]

  numeric_cols <- shared_cols[vapply(vista_df, is.numeric, logical(1))]
  other_cols <- setdiff(shared_cols, numeric_cols)
  numeric_checks <- lapply(numeric_cols, function(col) {
    .compare_numeric_vectors(vista_df[[col]], ref_df[[col]], tolerance = tolerance)
  })
  names(numeric_checks) <- numeric_cols

  non_numeric_identical <- TRUE
  if (length(other_cols)) {
    non_numeric_identical <- all(vapply(other_cols, function(col) {
      identical(as.character(vista_df[[col]]), as.character(ref_df[[col]]))
    }, logical(1)))
  }

  numeric_identical <- all(vapply(numeric_checks, `[[`, logical(1), "identical"))
  max_abs_diff <- 0
  if (length(numeric_checks)) {
    max_abs_diff <- max(vapply(numeric_checks, `[[`, numeric(1), "max_abs_diff"), na.rm = TRUE)
    if (!is.finite(max_abs_diff)) {
      max_abs_diff <- 0
    }
  }

  detail <- if (numeric_identical && non_numeric_identical) {
    "identical"
  } else {
    paste0(
      "non_numeric_identical=", non_numeric_identical,
      ", numeric_identical=", numeric_identical,
      ", max_abs_diff=", signif(max_abs_diff, 4)
    )
  }

  tibble::tibble(
    method = method,
    comparison = if (is.na(comparison)) NA_character_ else comparison,
    visual = visual,
    pass = numeric_identical && non_numeric_identical,
    detail = detail
  )
}

#' @keywords internal
.compare_numeric_vectors <- function(x, y, tolerance = 1e-8) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  both_na <- is.na(x) & is.na(y)
  equal <- both_na | (is.finite(x) & is.finite(y) & abs(x - y) <= tolerance)
  equal[is.na(equal)] <- FALSE

  diff <- abs(x - y)
  diff[both_na] <- 0
  diff[!is.finite(diff)] <- Inf
  max_abs_diff <- suppressWarnings(max(diff, na.rm = TRUE))
  if (!is.finite(max_abs_diff)) {
    max_abs_diff <- 0
  }

  list(
    identical = all(equal),
    max_abs_diff = max_abs_diff
  )
}

#' @keywords internal
.compare_numeric_matrix <- function(x, y, tolerance = 1e-8) {
  check <- .compare_numeric_vectors(as.vector(x), as.vector(y), tolerance = tolerance)
  list(
    identical = check$identical &&
      identical(rownames(x), rownames(y)) &&
      identical(colnames(x), colnames(y)),
    max_abs_diff = check$max_abs_diff
  )
}

#' @keywords internal
.validate_vista_self_consistency <- function(x,
                                             p_value_type = "padj",
                                             log2fc_cutoff = 1,
                                             pval_cutoff = 0.05) {
  comps <- comparisons(x)
  summaries <- deg_summary(x)
  ref_rn <- rownames(norm_counts(x))
  issues <- character()

  for (comp_name in names(comps)) {
    tbl <- .align_de_to_counts(
      df = .tidy_de_results(comps[[comp_name]], rowname_col = "gene_id"),
      ref_rn = ref_rn,
      warn_missing = FALSE
    )
    p_col <- .p_col_for_cutoff(tbl, p_value_type = p_value_type)
    p_vals <- if (!is.na(p_col)) .safe_numeric(tbl[[p_col]], default = 1) else rep(1, nrow(tbl))
    expected_reg <- ifelse(
      !is.na(tbl$log2fc) & tbl$log2fc >= log2fc_cutoff & p_vals <= pval_cutoff,
      "Up",
      ifelse(
        !is.na(tbl$log2fc) & tbl$log2fc <= -log2fc_cutoff & p_vals <= pval_cutoff,
        "Down",
        "Other"
      )
    )

    if (!identical(as.character(tbl$regulation), expected_reg)) {
      issues <- c(issues, sprintf("Comparison '%s' has a regulation column inconsistent with the configured thresholds.", comp_name))
    }

    expected_summary <- tbl |>
      dplyr::filter(.data$regulation %in% c("Up", "Down")) |>
      dplyr::count(.data$regulation, name = "n") |>
      dplyr::mutate(regulation = factor(.data$regulation, levels = c("Up", "Down"))) |>
      dplyr::arrange(.data$regulation)

    observed_summary <- summaries[[comp_name]] |>
      as.data.frame(stringsAsFactors = FALSE) |>
      dplyr::filter(.data$regulation %in% c("Up", "Down")) |>
      dplyr::mutate(regulation = factor(.data$regulation, levels = c("Up", "Down"))) |>
      dplyr::arrange(.data$regulation)

    expected_summary <- tidyr::complete(
      expected_summary,
      regulation = factor(c("Up", "Down"), levels = c("Up", "Down")),
      fill = list(n = 0)
    ) |>
      dplyr::arrange(.data$regulation)

    observed_summary <- tidyr::complete(
      observed_summary,
      regulation = factor(c("Up", "Down"), levels = c("Up", "Down")),
      fill = list(n = 0)
    ) |>
      dplyr::arrange(.data$regulation)

    if (!identical(as.integer(expected_summary$n), as.integer(observed_summary$n))) {
      issues <- c(issues, sprintf("Comparison '%s' has a DEG summary that does not match its regulation calls.", comp_name))
    }
  }

  list(valid = length(issues) == 0, issues = unique(issues))
}

#' @keywords internal
.ma_payload_from_vista <- function(x, comparison) {
  tbl <- .align_de_to_counts(
    df = .tidy_de_results(comparisons(x)[[comparison]], rowname_col = "gene_id"),
    ref_rn = rownames(norm_counts(x)),
    warn_missing = FALSE
  )
  tbl$baseMean <- rowMeans(norm_counts(x))[match(tbl$gene_id, rownames(norm_counts(x)))]
  tbl |>
    dplyr::transmute(
      gene_id = .data$gene_id,
      baseMean = .data$baseMean,
      log10_baseMean = log10(.data$baseMean + 1),
      log2fc = .data$log2fc,
      regulation = as.character(.data$regulation)
    )
}

#' @keywords internal
.ma_payload_from_reference <- function(reference, comparison) {
  tbl <- .align_de_to_counts(
    df = .tidy_de_results(reference$comparisons[[comparison]], rowname_col = "gene_id"),
    ref_rn = rownames(reference$norm_counts),
    warn_missing = FALSE
  )
  tbl$baseMean <- rowMeans(reference$norm_counts)[match(tbl$gene_id, rownames(reference$norm_counts))]
  tbl |>
    dplyr::transmute(
      gene_id = .data$gene_id,
      baseMean = .data$baseMean,
      log10_baseMean = log10(.data$baseMean + 1),
      log2fc = .data$log2fc,
      regulation = as.character(.data$regulation)
    )
}

#' @keywords internal
.volcano_payload_from_vista <- function(x,
                                        comparison,
                                        p_value_type = "padj",
                                        log2fc_cutoff = 1,
                                        pval_cutoff = 0.05) {
  tbl <- .align_de_to_counts(
    df = .tidy_de_results(comparisons(x)[[comparison]], rowname_col = "gene_id"),
    ref_rn = rownames(norm_counts(x)),
    warn_missing = FALSE
  )
  .volcano_payload_from_table(tbl, p_value_type, log2fc_cutoff, pval_cutoff)
}

#' @keywords internal
.volcano_payload_from_reference <- function(reference,
                                            comparison,
                                            p_value_type = "padj",
                                            log2fc_cutoff = 1,
                                            pval_cutoff = 0.05) {
  tbl <- .align_de_to_counts(
    df = .tidy_de_results(reference$comparisons[[comparison]], rowname_col = "gene_id"),
    ref_rn = rownames(reference$norm_counts),
    warn_missing = FALSE
  )
  .volcano_payload_from_table(tbl, p_value_type, log2fc_cutoff, pval_cutoff)
}

#' @keywords internal
.volcano_payload_from_table <- function(tbl,
                                        p_value_type = "padj",
                                        log2fc_cutoff = 1,
                                        pval_cutoff = 0.05) {
  p_col <- .p_col_for_cutoff(tbl, p_value_type = p_value_type)
  p_vals <- if (!is.na(p_col)) .safe_numeric(tbl[[p_col]], default = 1) else rep(1, nrow(tbl))
  neg_log10 <- -log10(p_vals + 1e-300)
  tbl |>
    dplyr::transmute(
      gene_id = .data$gene_id,
      log2fc = .data$log2fc,
      pvalue = p_vals,
      neg_log10_p = neg_log10,
      significant = abs(.data$log2fc) >= log2fc_cutoff & p_vals <= pval_cutoff,
      regulation = as.character(.data$regulation)
    )
}

#' @keywords internal
.deg_count_payload_from_vista <- function(x) {
  .collect_deg_count_data(x) |>
    as.data.frame(stringsAsFactors = FALSE) |>
    dplyr::arrange(.data$sample_comparisons, .data$regulation)
}

#' @keywords internal
.deg_count_payload_from_reference <- function(reference) {
  df <- purrr::imap_dfr(reference$deg_summary, function(tbl, comp_name) {
    as.data.frame(tbl, stringsAsFactors = FALSE) |>
      dplyr::filter(.data$regulation %in% c("Up", "Down")) |>
      dplyr::transmute(
        sample_comparisons = comp_name,
        regulation = as.character(.data$regulation),
        n = as.numeric(.data$n)
      )
  })

  if (!nrow(df)) {
    return(df)
  }

  df$sample_comparisons <- factor(df$sample_comparisons, levels = unique(df$sample_comparisons))
  df$regulation <- factor(df$regulation, levels = c("Up", "Down"))
  tidyr::complete(
    df,
    sample_comparisons,
    regulation,
    fill = list(n = 0)
  ) |>
    dplyr::arrange(.data$sample_comparisons, .data$regulation)
}

#' @keywords internal
.pca_payload_from_vista <- function(x) {
  mat <- norm_counts(x)
  meta <- as.data.frame(SummarizedExperiment::colData(x), stringsAsFactors = FALSE)
  meta$sample <- meta$sample_names
  pca <- stats::prcomp(t(mat), center = TRUE, scale. = TRUE)
  .prepare_pca_dataframe(pca, meta) |>
    dplyr::select("sample", "PC1", "PC2", tidyselect::everything()) |>
    dplyr::arrange(.data$sample)
}

#' @keywords internal
.pca_payload_from_reference <- function(reference, vista_obj) {
  mat <- as.matrix(reference$norm_counts)
  meta <- as.data.frame(reference$sample_info, stringsAsFactors = FALSE)
  if (!"sample_names" %in% names(meta)) {
    meta$sample_names <- rownames(meta)
  }
  meta$sample <- meta$sample_names
  pca <- stats::prcomp(t(mat), center = TRUE, scale. = TRUE)
  pca_df <- .prepare_pca_dataframe(pca, meta) |>
    dplyr::select("sample", "PC1", "PC2", tidyselect::everything()) |>
    dplyr::arrange(.data$sample)

  vista_group <- tryCatch(S4Vectors::metadata(vista_obj)$group$column, error = function(e) NULL)
  if (!is.null(vista_group) && vista_group %in% names(pca_df)) {
    pca_df[[vista_group]] <- as.character(pca_df[[vista_group]])
  }
  pca_df
}

#' @keywords internal
.align_pca_signs <- function(vista_df, ref_df) {
  out <- ref_df
  for (pc in c("PC1", "PC2")) {
    if (!pc %in% names(vista_df) || !pc %in% names(ref_df)) {
      next
    }
    corr <- suppressWarnings(stats::cor(vista_df[[pc]], ref_df[[pc]], use = "pairwise.complete.obs"))
    if (is.finite(corr) && corr < 0) {
      out[[pc]] <- -out[[pc]]
    }
  }
  out
}

#' @keywords internal
.plot_reference_ma <- function(df,
                               colors = c(Up = "#a40000", Down = "#16317d", Other = "gray70")) {
  ggplot2::ggplot(df, ggplot2::aes(x = .data$log10_baseMean, y = .data$log2fc, color = .data$regulation)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.2) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(x = "log10(baseMean + 1)", y = "log2 fold change", color = "Regulation")
}

#' @keywords internal
.plot_reference_volcano <- function(df,
                                    colors = c(Up = "#a40000", Down = "#007e2f", Other = "grey70")) {
  ggplot2::ggplot(df, ggplot2::aes(x = .data$log2fc, y = .data$neg_log10_p, color = .data$regulation)) +
    ggplot2::geom_point(size = 1, alpha = 0.8) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(x = "log2 fold change", y = "-log10(p-value)", color = "Regulation")
}

#' @keywords internal
.plot_reference_deg_count <- function(df,
                                      colors = c(Up = "red4", Down = "blue4")) {
  ggplot2::ggplot(df, ggplot2::aes(x = .data$sample_comparisons, y = .data$n, fill = .data$regulation)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = NULL, y = "Gene Count", fill = "Regulation")
}

#' @keywords internal
.plot_reference_pca <- function(df, vista_obj) {
  group_col <- tryCatch(S4Vectors::metadata(vista_obj)$group$column, error = function(e) NULL)
  if (is.null(group_col) || !group_col %in% names(df)) {
    group_col <- names(df)[!names(df) %in% c("sample", "PC1", "PC2")][1]
  }
  color_vals <- group_colors(vista_obj)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$PC1, y = .data$PC2, color = .data[[group_col]])) +
    ggplot2::geom_point(size = 10) +
    ggplot2::scale_color_manual(values = color_vals) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "PCA Plot", x = "PC1", y = "PC2", color = group_col)
}
