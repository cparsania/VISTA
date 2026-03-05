test_that("expression heatmap supports rowData display_id", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  vista <- make_small_vista()
  rowData(vista)$SYMBOL <- paste0("SYM", seq_len(nrow(vista)))
  genes <- rownames(vista)[1:10]
  groups <- unique(SummarizedExperiment::colData(vista)$cond_long)

  hm <- get_expression_heatmap(
    vista,
    sample_group = groups,
    genes = genes,
    display_id = "SYMBOL",
    show_row_names = TRUE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    summarise_replicates = FALSE,
    return_type = "heatmap"
  )
  expect_true(inherits(hm, "Heatmap") || inherits(hm, "HeatmapList"))
})

test_that("expression heatmap supports column annotation with summarised replicates", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  vista <- make_small_vista()
  genes <- rownames(vista)[1:12]
  groups <- unique(SummarizedExperiment::colData(vista)$cond_long)

  hm <- get_expression_heatmap(
    vista,
    sample_group = groups,
    genes = genes,
    annotate_columns = TRUE,
    summarise_replicates = TRUE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    return_type = "heatmap"
  )
  expect_true(inherits(hm, "Heatmap") || inherits(hm, "HeatmapList"))
})

test_that("expression heatmap supports multiple column annotations and custom cluster_by", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  vista <- make_small_vista()
  genes <- rownames(vista)[1:12]
  groups <- unique(SummarizedExperiment::colData(vista)$cond_long)
  SummarizedExperiment::colData(vista)$sample_batch <- rep(c("B1", "B2"), length.out = ncol(vista))

  hm <- get_expression_heatmap(
    vista,
    sample_group = groups,
    genes = genes,
    annotate_columns = c("cond_long", "sample_batch"),
    cluster_by = "sample_batch",
    summarise_replicates = FALSE,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    return_type = "heatmap"
  )
  expect_true(inherits(hm, "Heatmap") || inherits(hm, "HeatmapList"))
})

test_that("expression heatmap validates cluster_by against active annotation columns", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  vista <- make_small_vista()
  genes <- rownames(vista)[1:12]
  groups <- unique(SummarizedExperiment::colData(vista)$cond_long)

  expect_error(
    get_expression_heatmap(
      vista,
      sample_group = groups,
      genes = genes,
      annotate_columns = c("cond_long"),
      cluster_by = "nonexistent_col",
      summarise_replicates = FALSE,
      return_type = "heatmap"
    ),
    "cluster_by"
  )
})

test_that("expression heatmap matrix values and row order match source data", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  vista <- make_small_vista()
  genes <- rownames(vista)[1:12]
  si <- as.data.frame(SummarizedExperiment::colData(vista))
  si$sample <- rownames(si)
  group_col <- S4Vectors::metadata(vista)$group$column
  groups <- unique(si[[group_col]])
  sample_ids <- unlist(lapply(groups, function(g) si$sample[si[[group_col]] == g]))

  hm <- get_expression_heatmap(
    vista,
    genes = genes,
    sample_group = groups,
    value_transform = "raw",
    summarise_replicates = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    return_type = "heatmap"
  )

  expect_s4_class(hm, "Heatmap")
  expected <- as.matrix(SummarizedExperiment::assay(vista)[genes, sample_ids, drop = FALSE])
  expect_equal(rownames(hm@matrix), genes)
  expect_equal(colnames(hm@matrix), sample_ids)
  expect_equal(unname(hm@matrix), unname(expected))
})

test_that("expression heatmap kmeans clusters match matrix clustering partition", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  vista <- make_small_vista()
  genes <- rownames(vista)[1:20]
  si <- as.data.frame(SummarizedExperiment::colData(vista))
  si$sample <- rownames(si)
  group_col <- S4Vectors::metadata(vista)$group$column
  groups <- unique(si[[group_col]])
  sample_ids <- unlist(lapply(groups, function(g) si$sample[si[[group_col]] == g]))
  mat <- as.matrix(SummarizedExperiment::assay(vista)[genes, sample_ids, drop = FALSE])

  set.seed(2026)
  out <- get_expression_heatmap(
    vista,
    genes = genes,
    sample_group = groups,
    value_transform = "raw",
    summarise_replicates = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    kmeans_k = 3,
    return_type = "clusters"
  )

  expect_true(all(c("gene", "cluster") %in% colnames(out)))
  expect_equal(sort(out$gene), sort(genes))
  expect_false(anyNA(out$cluster))
  expect_true(all(out$cluster %in% 1:3))

  set.seed(2026)
  km <- stats::kmeans(mat, centers = 3, nstart = 1)

  obs <- out$cluster[match(genes, out$gene)]
  man <- km$cluster[match(genes, rownames(mat))]
  obs_same <- outer(obs, obs, `==`)
  man_same <- outer(man, man, `==`)
  expect_equal(unname(obs_same), unname(man_same))
})

.prepare_foldchange_fixture <- function(vista, min_genes = 3L) {
  comps <- names(comparisons(vista))
  stopifnot(length(comps) >= 1L)
  comp_tbl <- as.data.frame(comparisons(vista)[[comps[1]]], stringsAsFactors = FALSE)

  if (!"gene_id" %in% colnames(comp_tbl)) {
    comp_tbl$gene_id <- rownames(comp_tbl)
  }
  if (!"log2fc" %in% colnames(comp_tbl)) {
    if ("log2FoldChange" %in% colnames(comp_tbl)) {
      comp_tbl$log2fc <- comp_tbl$log2FoldChange
    } else if ("logFC" %in% colnames(comp_tbl)) {
      comp_tbl$log2fc <- comp_tbl$logFC
    } else {
      stop("No log2fc-like column available in comparison table.")
    }
  }

  # Align row_data gene_id to comparison gene IDs so internal validation
  # uses the same identifier space as fold-change tables.
  comp_gene_map <- stats::setNames(as.character(comp_tbl$gene_id), rownames(comp_tbl))
  aligned_gene_id <- unname(comp_gene_map[rownames(vista)])
  SummarizedExperiment::rowData(vista)$gene_id <- aligned_gene_id

  valid_genes <- unique(as.character(comp_tbl$gene_id))
  valid_genes <- valid_genes[!is.na(valid_genes) & nzchar(valid_genes)]
  valid_genes <- intersect(valid_genes, unique(as.character(SummarizedExperiment::rowData(vista)$gene_id)))

  fc_full <- get_foldchange_matrix(vista, sample_comparisons = comps, genes = valid_genes)
  complete_genes <- rownames(fc_full)[rowSums(is.na(fc_full)) == 0]
  complete_genes <- unique(complete_genes[!is.na(complete_genes) & nzchar(complete_genes)])

  if (length(complete_genes) < min_genes) {
    skip("Insufficient complete fold-change genes in treatment vs control comparison.")
  }

  list(
    vista = vista,
    comps = comps,
    genes = complete_genes
  )
}

test_that("fold-change heatmap matrix values and row order match foldchange matrix", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  fixture <- .prepare_foldchange_fixture(make_small_vista(), min_genes = 2L)
  vista <- fixture$vista
  comps <- fixture$comps
  genes <- head(fixture$genes, 15)

  hm <- get_foldchange_heatmap(
    vista,
    sample_comparisons = comps,
    genes = genes,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    return_type = "heatmap"
  )

  expect_s4_class(hm, "Heatmap")
  expected <- get_foldchange_matrix(vista, sample_comparisons = comps, genes = genes)
  expect_equal(rownames(hm@matrix), genes)
  expect_equal(colnames(hm@matrix), colnames(expected))
  expect_equal(unname(hm@matrix), unname(expected))
})

test_that("fold-change heatmap kmeans clusters match matrix clustering partition", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  fixture <- .prepare_foldchange_fixture(make_small_vista(), min_genes = 3L)
  vista <- fixture$vista
  comps <- fixture$comps
  genes <- head(fixture$genes, 20)
  mat <- get_foldchange_matrix(vista, sample_comparisons = comps, genes = genes)

  set.seed(2026)
  out <- get_foldchange_heatmap(
    vista,
    sample_comparisons = comps,
    genes = genes,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    kmeans_k = 3,
    return_type = "clusters"
  )

  expect_true(all(c("gene", "cluster") %in% colnames(out)))
  expect_equal(sort(out$gene), sort(genes))
  expect_false(anyNA(out$cluster))
  expect_true(all(out$cluster %in% 1:3))

  set.seed(2026)
  km <- stats::kmeans(mat, centers = 3, nstart = 1)

  obs <- out$cluster[match(genes, out$gene)]
  man <- km$cluster[match(genes, rownames(mat))]
  obs_same <- outer(obs, obs, `==`)
  man_same <- outer(man, man, `==`)
  expect_equal(unname(obs_same), unname(man_same))
})
