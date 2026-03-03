test_that("get_pca_plot produces valid plot", {
  vista <- make_small_vista()
  p <- get_pca_plot(vista)
  expect_s3_class(p, "ggplot")
})

test_that("get_pca_plot handles label_replicates parameter", {
  vista <- make_small_vista()

  p_no_labels <- get_pca_plot(vista, label_replicates = FALSE)
  expect_s3_class(p_no_labels, "ggplot")

  p_with_labels <- get_pca_plot(vista, label_replicates = TRUE)
  expect_s3_class(p_with_labels, "ggplot")
})

test_that("get_pca_plot handles top_n_genes parameter", {
  vista <- make_small_vista()

  p <- get_pca_plot(vista, top_n_genes = 50)
  expect_s3_class(p, "ggplot")
})

test_that("get_mds_plot produces valid plot", {
  vista <- make_small_vista()
  p <- get_mds_plot(vista)
  expect_s3_class(p, "ggplot")
})

test_that("get_mds_plot handles top_n_genes parameter", {
  vista <- make_small_vista()

  p <- get_mds_plot(vista, top_n_genes = 100)
  expect_s3_class(p, "ggplot")
})

test_that("get_umap_plot produces valid plot", {
  skip_if_not_installed("uwot")
  vista <- make_small_vista()
  nn <- max(2, ncol(vista) - 1)
  p <- get_umap_plot(vista, n_neighbors = nn)
  expect_s3_class(p, "ggplot")
})

test_that("get_umap_plot supports color_by metadata", {
  skip_if_not_installed("uwot")
  vista <- make_small_vista()
  nn <- max(2, ncol(vista) - 1)
  p <- get_umap_plot(vista, color_by = "cell", shape_by = "cond_long", n_neighbors = nn)
  expect_s3_class(p, "ggplot")
})

test_that("get_volcano_plot returns ggplot", {
  vista <- make_small_vista()
  comps <- comparisons(vista)

  # EnhancedVolcano currently emits ggplot2 deprecation warnings internally.
  p <- suppressWarnings(get_volcano_plot(vista, sample_comparison = names(comps)[1]))
  expect_s3_class(p, "ggplot")
})

test_that("get_volcano_plot handles display_id parameter", {
  vista <- make_small_vista()
  rowData(vista)$SYMBOL <- paste0("GENE", seq_len(nrow(vista)))
  comps <- comparisons(vista)

  p <- suppressWarnings(get_volcano_plot(
    vista,
    sample_comparison = names(comps)[1],
    display_id = "SYMBOL"
  ))
  expect_s3_class(p, "ggplot")
})

test_that("get_ma_plot returns ggplot", {
  vista <- make_small_vista()
  comps <- comparisons(vista)

  p <- get_ma_plot(vista, sample_comparison = names(comps)[1])
  expect_s3_class(p, "ggplot")
})

test_that("get_ma_plot handles topn parameter", {
  vista <- make_small_vista()
  comps <- comparisons(vista)

  p <- get_ma_plot(vista, sample_comparison = names(comps)[1], topn = 5)
  expect_s3_class(p, "ggplot")
})

test_that("get_corr_heatmap works with all samples", {
  vista <- make_small_vista()
  p <- get_corr_heatmap(vista)
  expect_s3_class(p, c("gg", "ggplot"))
})

test_that("get_corr_heatmap works with sample subsetting", {
  vista <- make_small_vista()
  groups <- unique(SummarizedExperiment::colData(vista)$cond_long)

  p <- get_corr_heatmap(vista, sample_group = groups[1])
  expect_s3_class(p, c("gg", "ggplot"))
})

test_that("get_deg_count_barplot summarizes correctly", {
  vista <- make_small_vista()
  p <- get_deg_count_barplot(vista)
  expect_s3_class(p, "ggplot")
})

test_that("get_deg_count_barplot handles facet options", {
  vista <- make_small_vista()
  comps <- comparisons(vista)

  p_reg <- get_deg_count_barplot(vista, facet_by = "regulation")
  expect_s3_class(p_reg, "ggplot")
})

test_that("get_deg_count_pieplot summarizes correctly", {
  vista <- make_small_vista()
  p <- get_deg_count_pieplot(vista)
  expect_s3_class(p, "ggplot")
})

test_that("get_deg_count_pieplot supports single pie mode", {
  vista <- make_small_vista()
  comp_name <- names(comparisons(vista))[1]
  p <- get_deg_count_pieplot(
    vista,
    sample_comparisons = comp_name,
    facet_by = "none"
  )
  expect_s3_class(p, "ggplot")
})

test_that("get_deg_count_donutplot summarizes correctly", {
  vista <- make_small_vista()
  p <- get_deg_count_donutplot(vista)
  expect_s3_class(p, "ggplot")
})

test_that("get_deg_count_donutplot handles label toggles", {
  vista <- make_small_vista()
  p <- get_deg_count_donutplot(vista, show_counts = FALSE, show_percent = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("get_expression_heatmap returns ComplexHeatmap object", {
  skip_if_not_installed("ComplexHeatmap")
  vista <- make_small_vista()
  genes <- rownames(vista)[1:20]
  groups <- unique(SummarizedExperiment::colData(vista)$cond_long)

  hm <- get_expression_heatmap(vista, genes = genes, samples = groups)
  expect_true(!is.null(hm))
})

test_that("get_expression_heatmap handles kmeans clustering", {
  skip_if_not_installed("ComplexHeatmap")
  vista <- make_small_vista()
  genes <- rownames(vista)[1:30]
  groups <- unique(SummarizedExperiment::colData(vista)$cond_long)

  hm <- get_expression_heatmap(vista, genes = genes, samples = groups, kmeans_k = 3)
  expect_true(!is.null(hm))
})
