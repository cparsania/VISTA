test_that("expression heatmap supports rowData display_id", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  vista <- make_small_vista()
  rowData(vista)$SYMBOL <- paste0("SYM", seq_len(nrow(vista)))
  genes <- rownames(vista)[1:10]
  groups <- unique(SummarizedExperiment::colData(vista)$cond_long)

  hm <- get_expression_heatmap(
    vista,
    samples = groups,
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
    samples = groups,
    genes = genes,
    annotate_columns = TRUE,
    summarise_replicates = TRUE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    return_type = "heatmap"
  )
  expect_true(inherits(hm, "Heatmap") || inherits(hm, "HeatmapList"))
})
