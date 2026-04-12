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

test_that("foldchange barplot, lollipop, and boxplot support display mapping and facet layout", {
  vista <- make_small_vista()
  rowData(vista)$SYMBOL <- paste0("SYM", seq_len(nrow(vista)))
  genes_sym <- rowData(vista)$SYMBOL[1:3]
  comp <- names(comparisons(vista))[1]

  p_bar <- get_foldchange_barplot(
    vista,
    genes = genes_sym,
    sample_comparisons = comp,
    display_id = "SYMBOL",
    facet_by = "none"
  )
  expect_s3_class(p_bar, "ggplot")

  p_lollipop <- get_foldchange_lollipop(
    vista,
    sample_comparison = comp,
    genes = genes_sym,
    display_id = "SYMBOL",
    facet_by = "gene",
    facet_ncol = 2
  )
  expect_s3_class(p_lollipop, "ggplot")
  expect_equal(p_lollipop$facet$params$ncol, 2)

  p_box <- get_foldchange_boxplot(
    vista,
    genes = genes_sym,
    sample_comparisons = comp,
    display_id = "SYMBOL",
    facet_by = "comparison",
    facet_ncol = 1
  )
  expect_s3_class(p_box, "ggplot")
  expect_equal(p_box$facet$params$ncol, 1)
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

test_that("foldchange plot APIs retain display and facet arguments where applicable", {
  expect_true(all(c("display_id", "display_from", "display_orgdb") %in% names(formals(get_foldchange_chromosome_plot))))
  expect_true(all(c("display_id", "display_from", "display_orgdb", "facet_nrow", "facet_ncol") %in% names(formals(get_foldchange_boxplot))))
  expect_true(all(c("display_id", "display_from", "display_orgdb", "facet_nrow", "facet_ncol") %in% names(formals(get_foldchange_raincloud))))
  expect_true(all(c("genes", "display_id", "display_from", "display_orgdb") %in% names(formals(get_foldchange_scatter))))
  expect_true(all(c("display_id", "display_from", "display_orgdb", "facet_nrow", "facet_ncol") %in% names(formals(get_foldchange_lollipop))))
  expect_true(all(c("display_id", "display_from", "display_orgdb", "facet_nrow", "facet_ncol") %in% names(formals(get_foldchange_lineplot))))
  expect_true(all(c("display_id", "display_from", "display_orgdb", "facet_nrow", "facet_ncol") %in% names(formals(get_foldchange_barplot))))
  expect_true(all(c("display_id", "display_from", "display_orgdb") %in% names(formals(get_foldchange_heatmap))))
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

test_that("get_foldchange_lineplot supports display mapping and cluster facet layout", {
  vista <- make_small_vista()
  rowData(vista)$SYMBOL <- paste0("SYM", seq_len(nrow(vista)))
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())
  cell_levels <- unique(sample_metadata$cell)

  vista_multi <- create_vista(
    counts = count_data[seq_len(150), ],
    sample_info = sample_metadata,
    column_geneid = "gene_id",
    group_column = "cell",
    group_numerator = cell_levels[2:3],
    group_denominator = rep(cell_levels[1], 2),
    min_counts = 5,
    min_replicates = 1
  )
  rowData(vista_multi)$SYMBOL <- paste0("SYM", seq_len(nrow(vista_multi)))
  genes_sym <- rowData(vista_multi)$SYMBOL[1:6]

  out <- get_foldchange_lineplot(
    vista_multi,
    sample_comparisons = names(comparisons(vista_multi))[1:2],
    genes = genes_sym,
    display_id = "SYMBOL",
    km = 2,
    facet_by = "cluster",
    facet_ncol = 1
  )

  expect_type(out, "list")
  expect_s3_class(out$plot, "ggplot")
  expect_true(all(c("gene_id", "display_gene", "cluster") %in% colnames(out$clustered_data)))
  expect_equal(out$plot$facet$params$ncol, 1)
})
