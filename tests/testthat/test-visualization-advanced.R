test_that("get_pca_plot produces valid plot", {
  vista <- make_small_vista()
  p <- get_pca_plot(vista)
  expect_s3_class(p, "ggplot")
})

test_that("get_pca_plot handles label_replicates parameter", {
  vista <- make_small_vista()

  p_no_labels <- get_pca_plot(vista, label = FALSE)
  expect_s3_class(p_no_labels, "ggplot")

  p_with_labels <- get_pca_plot(vista, label = TRUE)
  expect_s3_class(p_with_labels, "ggplot")
})

test_that("get_pca_plot handles top_n_genes parameter", {
  vista <- make_small_vista()

  p <- get_pca_plot(vista, top_n_genes = 50)
  expect_s3_class(p, "ggplot")
})

test_that("embedding plots support harmonized color and label arguments", {
  vista <- make_small_vista()

  p_pca <- get_pca_plot(
    vista,
    color_by = "cell",
    label = TRUE,
    use_group_colors = FALSE,
    palette = "Set 2"
  )
  expect_s3_class(p_pca, "ggplot")

  p_mds <- get_mds_plot(
    vista,
    color_by = "cell",
    label = FALSE,
    use_group_colors = FALSE,
    colors = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
  )
  expect_s3_class(p_mds, "ggplot")
})

test_that("embedding plots still accept deprecated use_vista_colors alias", {
  vista <- make_small_vista()

  expect_warning(
    p <- get_mds_plot(vista, color_by = "cell", use_vista_colors = FALSE),
    "deprecated"
  )
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

  p <- get_ma_plot(vista, sample_comparison = names(comps)[1], label_n = 5)
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

test_that("get_deg_count_pieplot can include non-DE genes and set text color", {
  vista <- make_small_vista()
  p <- get_deg_count_pieplot(
    vista,
    show_other = TRUE,
    other_color = "grey80",
    text_color = "navy"
  )

  expect_s3_class(p, "ggplot")
  expect_true("Other" %in% as.character(p$data$regulation))
  expect_identical(p$layers[[2]]$aes_params$colour, "navy")
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
  p <- get_deg_count_donutplot(vista, label = "percent")
  expect_s3_class(p, "ggplot")
})

test_that("get_deg_count_donutplot can include non-DE genes", {
  vista <- make_small_vista()
  p <- get_deg_count_donutplot(vista, show_other = TRUE, facet_by = "comparison")

  expect_s3_class(p, "ggplot")
  expect_true("Other" %in% as.character(p$data$regulation))
})

test_that("get_expression_heatmap returns ComplexHeatmap object", {
  skip_if_not_installed("ComplexHeatmap")
  vista <- make_small_vista()
  genes <- rownames(vista)[1:20]
  groups <- unique(SummarizedExperiment::colData(vista)$cond_long)

  hm <- get_expression_heatmap(vista, sample_group = groups, genes = genes)
  expect_true(!is.null(hm))
})

test_that("get_expression_heatmap handles kmeans clustering", {
  skip_if_not_installed("ComplexHeatmap")
  vista <- make_small_vista()
  genes <- rownames(vista)[1:30]
  groups <- unique(SummarizedExperiment::colData(vista)$cond_long)

  hm <- get_expression_heatmap(vista, genes = genes, sample_group = groups, kmeans_k = 3)
  expect_true(!is.null(hm))
})

test_that(".vista_comparison_colors falls back to a palette when nothing matches (B5)", {
  v <- make_small_vista()

  # An empty intersection used to yield character(0), which walked straight past
  # the internal `is.null(cols)` fallback and past every caller's guard. It must
  # now produce a usable palette covering what was asked for -- the same
  # behaviour .vista_group_colors() has always had.
  cols <- VISTA:::.vista_comparison_colors(v, comparisons_present = "NOT_A_COMPARISON")
  expect_length(cols, 1L)
  expect_named(cols, "NOT_A_COMPARISON")
  expect_false(is.na(cols[["NOT_A_COMPARISON"]]))

  # The two sibling helpers must stay in step.
  gcols <- VISTA:::.vista_group_colors(v, groups_present = "NOT_A_GROUP")
  expect_length(gcols, 1L)
  expect_named(gcols, "NOT_A_GROUP")

  # A real comparison still resolves to the stored colour.
  comp <- names(comparisons(v))[[1]]
  stored <- VISTA:::.vista_comparison_colors(v, comparisons_present = comp)
  expect_true(comp %in% names(stored))
  expect_identical(
    unname(stored[[comp]]),
    unname(S4Vectors::metadata(v)$comparison$colors[[comp]])
  )
})

test_that("plots fall back to a default palette when stored comparison colors drift", {
  v <- make_small_vista()

  # Simulate an object whose stored colour map no longer names its comparisons
  # (hand-edited metadata, or comparisons added after construction).
  md <- S4Vectors::metadata(v)
  md$comparison$colors <- c(SOME_STALE_COMPARISON = "#123456")
  S4Vectors::metadata(v) <- md

  # This previously died with "subscript out of bounds" at pal[[1]].
  expect_no_error(
    get_foldchange_lineplot(v, sample_comparisons = names(comparisons(v)))
  )
})

test_that("corr heatmap triangle mask matches the rendered axis order (B3)", {
  v <- make_small_vista()

  # The mask used to be computed against input-order factor codes while the axes
  # were re-levelled to hclust order afterwards, so the retained cells formed a
  # staircase rather than a triangle.
  count_per_row <- function(p) {
    d <- p$data
    lv <- levels(d$Var1)
    as.integer(table(factor(as.character(d$Var1), levels = lv)))
  }

  for (cb in c("correlation", "group", "input", "none")) {
    lower <- get_corr_heatmap(v, triangle = "lower", cluster_by = cb)
    expect_identical(
      count_per_row(lower), seq_len(ncol(v)),
      info = paste("lower /", cb)
    )

    upper <- get_corr_heatmap(v, triangle = "upper", cluster_by = cb)
    expect_identical(
      count_per_row(upper), rev(seq_len(ncol(v))),
      info = paste("upper /", cb)
    )
  }
})

test_that("corr heatmap full triangle and show_diagonal stay consistent", {
  v <- make_small_vista()
  n <- ncol(v)

  full <- get_corr_heatmap(v, triangle = "full")
  expect_identical(nrow(full$data), as.integer(n * n))

  no_diag <- get_corr_heatmap(v, triangle = "full", show_diagonal = FALSE)
  expect_identical(nrow(no_diag$data), as.integer(n * n - n))

  lower_no_diag <- get_corr_heatmap(v, triangle = "lower", show_diagonal = FALSE)
  expect_identical(nrow(lower_no_diag$data), as.integer(n * (n - 1) / 2))

  # Values must still belong to their own sample pair after re-levelling.
  mat <- log2(SummarizedExperiment::assay(v, "norm_counts") + 1)
  cm <- stats::cor(mat, method = "pearson", use = "pairwise.complete.obs")
  d <- full$data
  for (i in sample(seq_len(nrow(d)), min(10L, nrow(d)))) {
    expect_equal(
      d$value[[i]],
      cm[as.character(d$Var1[[i]]), as.character(d$Var2[[i]])]
    )
  }
})
