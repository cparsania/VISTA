test_that("foldchange scatter returns ggplot for two comparisons", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())
  cell_levels <- unique(sample_metadata$cell)

  # Need 2 comparisons for foldchange scatter
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

  comp_names <- names(comparisons(vista_multi))[1:2]
  p <- get_foldchange_scatter(vista_multi, sample_comparisons = comp_names)
  expect_s3_class(p, "ggplot")
})
