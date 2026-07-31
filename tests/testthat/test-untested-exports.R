# Coverage for exports that previously had no direct test.

test_that("get_expression_matrix subsets, summarises, and transforms", {
  v <- make_small_vista()
  genes <- rownames(v)[seq_len(5)]

  m <- get_expression_matrix(v, genes = genes)
  expect_true(is.matrix(m))
  expect_setequal(rownames(m), genes)
  expect_identical(colnames(m), colnames(norm_counts(v)))

  # Summarising collapses replicates onto group labels.
  grp <- as.character(sample_info(v)[[S4Vectors::metadata(v)$group$column]])
  ms <- get_expression_matrix(v, genes = genes, summarise = TRUE)
  expect_setequal(colnames(ms), unique(grp))
  expect_equal(nrow(ms), length(genes))

  # log2 transform is applied after subsetting.
  ml <- get_expression_matrix(v, genes = genes, transform = "log2")
  expect_equal(ml, log2(m + 1))

  # zscore rows are centred.
  mz <- get_expression_matrix(v, genes = genes, transform = "zscore")
  expect_true(all(abs(rowMeans(mz)) < 1e-8))

  # Sample subsetting keeps only the requested columns.
  keep <- colnames(m)[seq_len(2)]
  expect_identical(colnames(get_expression_matrix(v, genes = genes, sample_names = keep)), keep)
})

test_that("get_deg_venn_diagram returns a plot for two comparisons", {
  skip_if_not_installed("ggvenn")

  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  si <- sample_metadata[sample_metadata$cond_long %in% c("control", "treatment1"), ]
  cnt <- count_data[seq_len(300), c("gene_id", si$sample_names)]

  v <- suppressMessages(create_vista(
    counts = cnt,
    sample_info = si,
    column_geneid = "gene_id",
    group_column = "cond_long",
    group_numerator = c("treatment1", "control"),
    group_denominator = c("control", "treatment1"),
    min_counts = 5,
    min_replicates = 1
  ))

  comps <- names(comparisons(v))
  expect_length(comps, 2L)

  p <- get_deg_venn_diagram(v, sample_comparisons = comps, regulation = "Up")
  expect_s3_class(p, "ggplot")
})

test_that("get_deg_alluvial returns a plot across comparisons", {
  skip_if_not_installed("ggalluvial")

  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  si <- sample_metadata[sample_metadata$cond_long %in% c("control", "treatment1"), ]
  cnt <- count_data[seq_len(300), c("gene_id", si$sample_names)]

  v <- suppressMessages(create_vista(
    counts = cnt,
    sample_info = si,
    column_geneid = "gene_id",
    group_column = "cond_long",
    group_numerator = c("treatment1", "control"),
    group_denominator = c("control", "treatment1"),
    min_counts = 5,
    min_replicates = 1
  ))

  p <- get_deg_alluvial(v, sample_comparisons = names(comparisons(v)))
  expect_s3_class(p, "ggplot")
})

test_that("get_pairwise_corr_plot returns a GGally matrix", {
  skip_if_not_installed("GGally")

  v <- make_small_vista()
  p <- get_pairwise_corr_plot(v)
  expect_s3_class(p, "ggmatrix")
})

test_that("run_deseq_analysis returns the harmonized standalone structure", {
  v <- make_small_vista()
  si <- as.data.frame(sample_info(v))
  data("count_data", package = "VISTA", envir = environment())
  cnt <- count_data[seq_len(300), c("gene_id", si$sample_names), drop = FALSE]

  res <- suppressMessages(run_deseq_analysis(
    counts = cnt,
    sample_info = si,
    column_geneid = "gene_id",
    group_column = "cond_long",
    group_numerator = "treatment1",
    group_denominator = "control",
    min_counts = 5,
    min_replicates = 1
  ))

  expect_named(
    res,
    c("norm_counts", "sample_info", "row_data", "comparisons", "deg_summary")
  )
  expect_identical(names(res$comparisons), "treatment1_VS_control")
  expect_true(all(c("gene_id", "regulation") %in% colnames(res$comparisons[[1]])))
  expect_identical(colnames(res$norm_counts), si$sample_names)
})

test_that("run_edger_analysis and run_limma_analysis share the same schema", {
  v <- make_small_vista()
  si <- as.data.frame(sample_info(v))
  data("count_data", package = "VISTA", envir = environment())
  cnt <- count_data[seq_len(300), c("gene_id", si$sample_names), drop = FALSE]

  args <- list(
    counts = cnt,
    sample_info = si,
    column_geneid = "gene_id",
    group_column = "cond_long",
    group_numerator = "treatment1",
    group_denominator = "control",
    min_counts = 5,
    min_replicates = 1
  )

  res_edger <- suppressMessages(do.call(run_edger_analysis, args))
  res_limma <- suppressMessages(do.call(run_limma_analysis, args))

  for (res in list(res_edger, res_limma)) {
    expect_named(
      res,
      c("norm_counts", "sample_info", "row_data", "comparisons", "deg_summary")
    )
    expect_identical(names(res$comparisons), "treatment1_VS_control")
    expect_true(all(c("gene_id", "regulation") %in% colnames(res$comparisons[[1]])))
  }

  # Both backends filter identically, so they model the same feature set.
  expect_identical(rownames(res_edger$norm_counts), rownames(res_limma$norm_counts))
})
