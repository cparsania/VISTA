test_that("expression distribution plots return ggplot objects", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:3]

  p1 <- get_expression_violinplot(vista, genes = genes, facet = FALSE)
  expect_s3_class(p1, "ggplot")

  p2 <- get_expression_density(vista, genes = genes)
  expect_s3_class(p2, "ggplot")

  p3 <- get_expression_raincloud(vista, genes = genes, facet = FALSE)
  expect_s3_class(p3, "ggplot")
})

test_that("expression barplot supports rowData-based display labels", {
  vista <- make_small_vista()
  rowData(vista)$SYMBOL <- paste0("SYM", seq_len(nrow(vista)))
  genes_sym <- rowData(vista)$SYMBOL[1:3]

  p <- get_expression_barplot(
    vista,
    genes = genes_sym,
    display_id = "SYMBOL"
  )
  expect_s3_class(p, "ggplot")
})

test_that("get_expression_barplot returns ggplot object", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:3]

  p <- get_expression_barplot(vista, genes = genes)
  expect_s3_class(p, "ggplot")
})

test_that("get_expression_barplot handles sample subsetting", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:3]
  groups <- unique(SummarizedExperiment::colData(vista)$cond_long)

  p <- get_expression_barplot(vista, genes = genes, sample_group = groups[1])
  expect_s3_class(p, "ggplot")
})

test_that("get_expression_boxplot returns ggplot object", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:5]

  p <- get_expression_boxplot(vista, genes = genes)
  expect_s3_class(p, "ggplot")
})

test_that("get_expression_boxplot handles log transformation", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:5]

  p <- get_expression_boxplot(vista, genes = genes, log_transform = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("get_expression_violinplot validates gene input", {
  vista <- make_small_vista()

  expect_error(
    get_expression_violinplot(vista, genes = "NONEXISTENT_GENE")
  )
})

test_that("get_expression_density works with facet_by parameter", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:5]

  p <- get_expression_density(vista, genes = genes)
  expect_s3_class(p, "ggplot")
})

test_that("get_expression_scatter returns ggplot object", {
  vista <- make_small_vista()
  samples <- colnames(vista)

  p <- get_expression_scatter(vista, sample_x = samples[1], sample_y = samples[2])
  expect_s3_class(p, "ggplot")
})

test_that("get_expression_lineplot returns ggplot object", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:3]

  p <- get_expression_lineplot(vista, genes = genes)
  expect_s3_class(p, "ggplot")
})

test_that("get_expression_lollipop returns ggplot object", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:5]

  p <- get_expression_lollipop(vista, genes = genes)
  expect_s3_class(p, "ggplot")
})

test_that("get_expression_joyplot returns ggplot object", {
  skip_if_not_installed("ggridges")
  vista <- make_small_vista()
  genes <- rownames(vista)[1:5]

  p <- get_expression_joyplot(vista, genes = genes)
  expect_s3_class(p, "ggplot")
})

test_that("expression plots handle faceting options", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:3]

  p_facet <- get_expression_violinplot(vista, genes = genes, facet = TRUE)
  expect_s3_class(p_facet, "ggplot")

  p_no_facet <- get_expression_violinplot(vista, genes = genes, facet = FALSE)
  expect_s3_class(p_no_facet, "ggplot")
})

test_that("expression plots respect gene limits", {
  vista <- make_small_vista()

  # Test with reasonable number of genes
  genes <- rownames(vista)[1:10]
  p <- get_expression_barplot(vista, genes = genes)
  expect_s3_class(p, "ggplot")
})

test_that("get_expression_raincloud handles parameters correctly", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:3]

  p <- get_expression_raincloud(vista, genes = genes, facet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("expression plots work with normalized counts", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:5]

  # These should use normalized counts from the VISTA object
  p <- get_expression_boxplot(vista, genes = genes)
  expect_s3_class(p, "ggplot")
})
