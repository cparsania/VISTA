test_that("expression distribution plots return ggplot objects", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:3]

  p1 <- get_expression_violinplot(vista, genes = genes, facet_by = "none")
  expect_s3_class(p1, "ggplot")

  p2 <- get_expression_density(vista, genes = genes)
  expect_s3_class(p2, "ggplot")

  skip_if_not_installed("ggrain")
  p3 <- get_expression_raincloud(vista, genes = genes, facet_by = "none")
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

test_that("get_expression_barplot bar heights match group means", {
  vista <- make_small_vista()
  gene <- rownames(vista)[1]
  p <- get_expression_barplot(
    vista,
    genes = gene,
    log_transform = FALSE,
    by = "group",
    facet_by = "none"
  )
  expect_s3_class(p, "ggplot")

  si <- as.data.frame(SummarizedExperiment::colData(vista))
  si$sample <- rownames(si)
  group_col <- S4Vectors::metadata(vista)$group$column
  expr <- as.numeric(SummarizedExperiment::assay(vista)[gene, si$sample])
  expected_means <- tapply(expr, si[[group_col]], mean, na.rm = TRUE)

  gb <- ggplot2::ggplot_build(p)
  bar_y <- gb$data[[1]]$y
  expect_equal(sort(as.numeric(bar_y)), sort(as.numeric(expected_means)))
})

test_that("get_expression_barplot handles sample subsetting", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:3]
  groups <- unique(SummarizedExperiment::colData(vista)$cond_long)

  p <- get_expression_barplot(vista, genes = genes, sample_group = groups[1])
  expect_s3_class(p, "ggplot")
})

test_that("get_expression_barplot supports per-sample mode", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:2]

  p <- get_expression_barplot(
    vista,
    genes = genes,
    by = "sample",
    sample_order = "group"
  )
  expect_s3_class(p, "ggplot")
})

test_that("get_expression_barplot rejects stats_group in per-sample mode", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1]

  expect_error(
    get_expression_barplot(
      vista,
      genes = genes,
      by = "sample",
      stats_group = TRUE
    ),
    "only supported when"
  )
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

test_that("get_expression_boxplot supports harmonized by alias and sample ordering", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:3]

  p <- get_expression_boxplot(
    vista,
    genes = genes,
    pool_genes = TRUE,
    by = "sample",
    sample_order = "group",
    facet_by = "auto"
  )
  expect_s3_class(p, "ggplot")
})

test_that("get_expression_violinplot includes the boxplot-style arguments", {
  box_args <- names(formals(get_expression_boxplot))
  violin_args <- names(formals(get_expression_violinplot))

  expect_true(all(box_args %in% violin_args))
})

test_that("get_expression_violinplot validates gene input", {
  vista <- make_small_vista()

  expect_error(
    get_expression_violinplot(vista, genes = "NONEXISTENT_GENE")
  )
})

test_that("get_expression_violinplot supports boxplot-style display and stats arguments", {
  vista <- make_small_vista()
  rowData(vista)$SYMBOL <- paste0("SYM", seq_len(nrow(vista)))
  genes_sym <- rowData(vista)$SYMBOL[1:3]

  p <- get_expression_violinplot(
    vista,
    genes = genes_sym,
    display_id = "SYMBOL",
    facet_by = "none"
  )
  expect_s3_class(p, "ggplot")

  skip_if_not_installed("ggpubr")
  p_stats <- get_expression_violinplot(
    vista,
    genes = rownames(vista)[1:2],
    stats_group = TRUE,
    facet_by = "none"
  )
  expect_s3_class(p_stats, "ggplot")
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

test_that("get_expression_lineplot includes the boxplot-style arguments", {
  box_args <- names(formals(get_expression_boxplot))
  line_args <- names(formals(get_expression_lineplot))

  expect_true(all(box_args %in% line_args))
})

test_that("get_expression_lineplot supports boxplot-style display arguments", {
  vista <- make_small_vista()
  rowData(vista)$SYMBOL <- paste0("SYM", seq_len(nrow(vista)))
  genes_sym <- rowData(vista)$SYMBOL[1:3]

  p <- get_expression_lineplot(
    vista,
    genes = genes_sym,
    display_id = "SYMBOL",
    facet_by = "none"
  )
  expect_s3_class(p, "ggplot")

  p_pool <- get_expression_lineplot(
    vista,
    genes = rownames(vista)[1:5],
    pool_genes = TRUE
  )
  expect_s3_class(p_pool, "ggplot")
})

test_that("get_expression_lollipop returns ggplot object", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:5]

  p <- get_expression_lollipop(vista, genes = genes)
  expect_s3_class(p, "ggplot")
})

test_that("get_expression_lollipop supports per-sample mode", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:3]

  p <- get_expression_lollipop(
    vista,
    genes = genes,
    by = "sample",
    sample_order = "expression"
  )
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

  p_facet <- get_expression_violinplot(vista, genes = genes, facet_by = "auto")
  expect_s3_class(p_facet, "ggplot")

  p_no_facet <- get_expression_violinplot(vista, genes = genes, facet_by = "none")
  expect_s3_class(p_no_facet, "ggplot")

  p_gene <- get_expression_violinplot(vista, genes = genes, facet_by = "gene")
  expect_s3_class(p_gene, "ggplot")
})

test_that("expression violin plot enforces group-based semantics", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:3]

  expect_error(
    get_expression_violinplot(vista, genes = genes, by = "sample"),
    "supports only"
  )
})

test_that("expression raincloud enforces group-based semantics", {
  skip_if_not_installed("ggrain")
  vista <- make_small_vista()
  genes <- rownames(vista)[1:3]

  expect_error(
    get_expression_raincloud(vista, genes = genes, by = "sample"),
    "supports only"
  )
})

test_that("expression plots respect gene limits", {
  vista <- make_small_vista()

  # Test with reasonable number of genes
  genes <- rownames(vista)[1:10]
  p <- get_expression_barplot(vista, genes = genes)
  expect_s3_class(p, "ggplot")
})

test_that("get_expression_raincloud handles parameters correctly", {
  skip_if_not_installed("ggrain")
  vista <- make_small_vista()
  genes <- rownames(vista)[1:3]

  p <- get_expression_raincloud(vista, genes = genes, facet_by = "gene", label = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("get_expression_raincloud supports gene-level group summaries", {
  skip_if_not_installed("ggrain")
  vista <- make_small_vista()
  genes <- rownames(vista)[1:4]

  p <- get_expression_raincloud(
    vista,
    genes = genes,
    summarise = TRUE,
    facet_by = "none",
    id.long.var = "gene"
  )
  expect_s3_class(p, "ggplot")
})

test_that("get_expression_raincloud warns when summarised values are faceted by gene", {
  skip_if_not_installed("ggrain")
  vista <- make_small_vista()
  genes <- rownames(vista)[1:4]

  expect_warning(
    get_expression_raincloud(
      vista,
      genes = genes,
      summarise = TRUE,
      facet_by = "gene"
    ),
    "not informative"
  )
})

test_that("distribution plots support sample ordering and palette overrides", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:4]

  p_density <- get_expression_density(
    vista,
    genes = genes,
    color_by = "sample",
    sample_order = "expression",
    palette = "Set 2"
  )
  expect_s3_class(p_density, "ggplot")

  skip_if_not_installed("ggridges")
  p_joy <- get_expression_joyplot(
    vista,
    genes = genes,
    color_by = "sample",
    sample_order = "group",
    palette = "Set 2"
  )
  expect_s3_class(p_joy, "ggplot")
})

test_that("get_expression_scatter supports harmonized label aliases", {
  vista <- make_small_vista()
  samples <- colnames(vista)

  p <- get_expression_scatter(
    vista,
    sample_x = samples[1],
    sample_y = samples[2],
    label_n = 5,
    label_size = 2.5
  )
  expect_s3_class(p, "ggplot")
})

test_that("expression plots work with normalized counts", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:5]

  # These should use normalized counts from the VISTA object
  p <- get_expression_boxplot(vista, genes = genes)
  expect_s3_class(p, "ggplot")
})
