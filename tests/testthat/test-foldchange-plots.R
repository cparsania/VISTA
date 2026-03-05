test_that("get_foldchange_barplot supports facet_by gene", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:3]

  p <- get_foldchange_barplot(
    vista,
    genes = genes,
    facet_by = "gene"
  )

  expect_s3_class(p, "ggplot")
})

test_that("get_foldchange_barplot values match DE table values", {
  vista <- make_small_vista()
  comp <- names(comparisons(vista))[1]
  genes <- rownames(vista)[1:4]

  p <- get_foldchange_barplot(
    vista,
    genes = genes,
    sample_comparisons = comp,
    facet_by = "none",
    sort_by = "input"
  )
  expect_s3_class(p, "ggplot")

  de_tbl <- as.data.frame(comparisons(vista)[[comp]])
  expected <- de_tbl$log2fc[match(genes, de_tbl$gene_id)]

  gb <- ggplot2::ggplot_build(p)
  shown <- gb$data[[1]]$y
  expect_equal(sort(as.numeric(shown)), sort(as.numeric(expected)))
})

test_that("get_foldchange_lollipop supports facet_by gene", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:3]
  comp <- names(comparisons(vista))[1]

  p <- get_foldchange_lollipop(
    vista,
    sample_comparison = comp,
    genes = genes,
    facet_by = "gene"
  )

  expect_s3_class(p, "ggplot")
})

test_that("get_foldchange_barplot supports facet_by comparison", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())
  cell_levels <- unique(sample_metadata$cell)

  vista_multi <- create_vista(
    counts = count_data[1:150, ],
    sample_info = sample_metadata,
    column_geneid = "gene_id",
    group_column = "cell",
    group_numerator = cell_levels[2:3],
    group_denominator = rep(cell_levels[1], 2),
    min_counts = 5,
    min_replicates = 1
  )

  p <- get_foldchange_barplot(
    vista_multi,
    genes = rownames(vista_multi)[1:3],
    sample_comparisons = names(comparisons(vista_multi))[1:2],
    facet_by = "comparison"
  )

  expect_s3_class(p, "ggplot")
})

test_that("get_foldchange_boxplot supports harmonized comparison and facet arguments", {
  vista <- make_small_vista()

  p <- get_foldchange_boxplot(
    vista,
    sample_comparisons = names(comparisons(vista))[1],
    facet_by = "auto"
  )

  expect_s3_class(p, "ggplot")
})

test_that("get_foldchange_raincloud supports harmonized label and facet arguments", {
  skip_if_not_installed("ggrain")
  vista <- make_small_vista()

  p <- get_foldchange_raincloud(
    vista,
    sample_comparisons = names(comparisons(vista))[1],
    facet_by = "none",
    label = FALSE
  )

  expect_s3_class(p, "ggplot")
})
