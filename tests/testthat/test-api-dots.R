# `...` used to swallow unknown names silently, so a typo like
# get_expression_heatmap(v, gene = my_genes) plotted the default gene set
# instead of erroring.

test_that("unknown ... arguments are rejected with a suggestion", {
  v <- make_small_vista()
  comp <- names(comparisons(v))[[1]]
  genes <- rownames(v)[seq_len(5)]

  expect_error(
    suppressWarnings(get_volcano_plot(v, sample_comparison = comp, definitely_not_real = 1)),
    "unknown argument"
  )
  expect_error(
    get_expression_heatmap(v, genes = genes, definitely_not_real = 1),
    "unknown argument"
  )
  expect_error(
    get_foldchange_heatmap(v, genes = genes, definitely_not_real = 1),
    "unknown argument"
  )
  expect_error(
    save_vista_plot(get_pca_plot(v), file = tempfile(fileext = ".png"), definitely_not_real = 1),
    "unknown argument"
  )
})

test_that("a near-miss argument name suggests the intended one", {
  v <- make_small_vista()
  genes <- rownames(v)[seq_len(5)]

  msg <- tryCatch(
    get_expression_heatmap(v, genes = genes, kmeans_kk = 3),
    error = function(e) conditionMessage(e)
  )
  expect_match(msg, "did you mean", ignore.case = TRUE)
  expect_match(msg, "kmeans_k", fixed = TRUE)
})

test_that("legitimately forwarded arguments still reach the plotting engine", {
  skip_if_not_installed("ComplexHeatmap")
  v <- make_small_vista()
  genes <- rownames(v)[seq_len(10)]

  # column_title is a real ComplexHeatmap::Heatmap argument.
  expect_no_error(
    get_expression_heatmap(v, genes = genes, column_title = "kept working")
  )
  expect_no_error(
    get_foldchange_heatmap(v, genes = genes, column_title = "kept working")
  )
})

test_that("managed arguments cannot be smuggled through ...", {
  v <- make_small_vista()
  genes <- rownames(v)[seq_len(5)]

  expect_error(
    get_expression_heatmap(v, genes = genes, matrix = matrix(1)),
    "managed by"
  )
})

test_that("get_pathway_heatmap keeps its original explicit guard", {
  v <- make_small_vista()
  expect_error(
    get_pathway_heatmap(v, enrichment = list(), genes = "x"),
    "managed by|are managed"
  )
})
