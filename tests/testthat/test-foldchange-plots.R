test_that("get_foldchange_barplot supports facet_by gene without breaking legacy defaults", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:3]

  p <- get_foldchange_barplot(
    vista,
    genes = genes,
    facet_by = "gene"
  )

  expect_s3_class(p, "ggplot")
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

test_that("legacy facet_comparison still works for foldchange barplot", {
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
    facet = TRUE,
    facet_comparison = TRUE
  )

  expect_s3_class(p, "ggplot")
})
