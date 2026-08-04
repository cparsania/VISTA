# VISTA grew three incompatible return_type vocabularies for the same choice of
# "give me the picture, the numbers, or both". They are unified on
# c("plot", "data", "both"); the legacy spellings still work and warn.

test_that("the canonical vocabulary is accepted everywhere", {
  skip_if_not_installed("ComplexHeatmap")
  v <- make_small_vista()
  genes <- rownames(v)[seq_len(10)]

  hm <- get_expression_heatmap(v, genes = genes, return_type = "plot")
  expect_s4_class(hm, "Heatmap")

  fc <- get_foldchange_heatmap(v, genes = genes, return_type = "plot")
  expect_s4_class(fc, "Heatmap")

  both <- get_expression_heatmap(v, genes = genes, kmeans_k = 2, return_type = "both")
  expect_type(both, "list")
  expect_named(both, c("heatmap", "clusters"))

  clusters <- get_expression_heatmap(v, genes = genes, kmeans_k = 2, return_type = "data")
  expect_s3_class(clusters, "tbl_df")
  expect_true(all(c("gene", "cluster") %in% colnames(clusters)))
})

test_that("legacy return_type values still work and warn", {
  skip_if_not_installed("ComplexHeatmap")
  v <- make_small_vista()
  genes <- rownames(v)[seq_len(10)]

  expect_warning(
    old <- get_expression_heatmap(v, genes = genes, return_type = "heatmap"),
    class = "vista_deprecated_value"
  )
  new <- get_expression_heatmap(v, genes = genes, return_type = "plot")
  expect_identical(class(old), class(new))

  expect_warning(
    old_cl <- get_expression_heatmap(v, genes = genes, kmeans_k = 2, return_type = "clusters"),
    class = "vista_deprecated_value"
  )
  new_cl <- get_expression_heatmap(v, genes = genes, kmeans_k = 2, return_type = "data")
  expect_identical(colnames(old_cl), colnames(new_cl))
  expect_identical(nrow(old_cl), nrow(new_cl))

  # Legacy warnings carry the standard deprecation classes too.
  expect_warning(
    get_foldchange_heatmap(v, genes = genes, return_type = "heatmap"),
    class = "deprecatedWarning"
  )
})

test_that("get_celltype_heatmap and get_pathway_heatmap use the same vocabulary", {
  v <- make_small_vista()
  md <- S4Vectors::metadata(v)
  md$cell_fractions <- data.frame(
    A = stats::runif(ncol(v)), B = stats::runif(ncol(v)), C = stats::runif(ncol(v)),
    row.names = colnames(v)
  )
  S4Vectors::metadata(v) <- md

  mat <- get_celltype_heatmap(v, top_n = 2, return_type = "data")
  expect_true(is.matrix(mat))

  expect_warning(
    legacy <- get_celltype_heatmap(v, top_n = 2, return_type = "matrix"),
    class = "vista_deprecated_value"
  )
  expect_equal(legacy, mat)

  p <- get_celltype_heatmap(v, top_n = 2, return_type = "plot")
  expect_s3_class(p, "ggplot")
})

test_that("an unknown return_type is rejected with the canonical options", {
  v <- make_small_vista()
  genes <- rownames(v)[seq_len(5)]

  msg <- tryCatch(
    get_expression_heatmap(v, genes = genes, return_type = "nonsense"),
    error = function(e) conditionMessage(e)
  )
  expect_match(msg, "plot", fixed = TRUE)
  expect_match(msg, "data", fixed = TRUE)
  expect_match(msg, "both", fixed = TRUE)
})

test_that("get_foldchange_lineplot gains return_type without changing its default", {
  v <- make_small_vista()
  comps <- names(comparisons(v))

  # Default is "both" -- the list this function has always returned.
  default <- get_foldchange_lineplot(v, sample_comparisons = comps)
  expect_type(default, "list")
  expect_named(default, c("plot", "clustered_data"))

  explicit <- get_foldchange_lineplot(v, sample_comparisons = comps, return_type = "both")
  expect_named(explicit, c("plot", "clustered_data"))

  p <- get_foldchange_lineplot(v, sample_comparisons = comps, return_type = "plot")
  expect_s3_class(p, "ggplot")
  expect_equal(ggplot2::ggplot_build(p)$data, ggplot2::ggplot_build(default$plot)$data)

  d <- get_foldchange_lineplot(v, sample_comparisons = comps, return_type = "data")
  expect_s3_class(d, "data.frame")
  expect_equal(d, default$clustered_data)
})

test_that("shape-selecting return_type arguments are deliberately untouched", {
  # These pick the shape of a purely tabular result and have no plot component,
  # so they keep their own vocabularies.
  expect_identical(
    eval(formals(get_pathway_genes)$return_type),
    c("long", "list", "vector")
  )
  expect_identical(
    eval(formals(get_genes_by_regulation)$return_type),
    c("list", "table")
  )
  expect_identical(
    eval(formals(read_vista_counts)$return_type),
    c("list", "data.frame", "matrix")
  )
})

test_that("get_enrichment_chord returns data by default and can return a plot", {
  skip_if_not_installed("circlize")
  skip_if_not_installed("msigdbr")

  v <- make_small_vista()
  msig <- try(
    suppressMessages(get_msigdb_enrichment(
      v, sample_comparison = names(comparisons(v))[[1]],
      regulation = "Both", msigdb_category = "H", from_type = "ENSEMBL"
    )),
    silent = TRUE
  )
  skip_if(inherits(msig, "try-error") || is.null(msig$enrich), "enrichment unavailable")
  skip_if(nrow(as.data.frame(msig$enrich)) == 0, "no enriched terms in fixture")

  # A device must be open for base graphics to be recordable.
  f <- tempfile(fileext = ".png")
  on.exit(unlink(f), add = TRUE)
  grDevices::png(f)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  VISTA:::.vista_enable_display_list()

  # Default is "data": the invisible list this function has always returned.
  res <- get_enrichment_chord(msig, top_n = 3)
  expect_type(res, "list")
  expect_true(all(c("gene_data", "hub_genes") %in% names(res)))

  both <- get_enrichment_chord(msig, top_n = 3, return_type = "both")
  expect_true(all(c("gene_data", "hub_genes", "plot") %in% names(both)))
})

test_that("save_vista_plot writes a recorded base-graphics plot", {
  # circlize output is a base-graphics side effect, so save_vista_plot needs a
  # replay path rather than a ggsave path.
  f_dev <- tempfile(fileext = ".png")
  grDevices::png(f_dev)
  # png() keeps its display list off by default, which is exactly why
  # get_enrichment_chord() enables it before drawing.
  VISTA:::.vista_enable_display_list()
  plot(1:10, 1:10)
  rp <- VISTA:::.vista_record_plot()
  grDevices::dev.off()
  unlink(f_dev)

  skip_if(is.null(rp), "graphics device does not support display lists")
  expect_s3_class(rp, "recordedplot")

  out <- tempfile(fileext = ".png")
  on.exit(unlink(out), add = TRUE)
  path <- save_vista_plot(rp, file = out, width = 5, height = 4)

  expect_true(file.exists(path))
  expect_gt(file.size(path), 0)

  out_pdf <- tempfile(fileext = ".pdf")
  on.exit(unlink(out_pdf), add = TRUE)
  expect_true(file.exists(save_vista_plot(rp, file = out_pdf, width = 5, height = 4)))
})

test_that("an empty recording is reported rather than silently saved as nothing", {
  # Without dev.control("enable") a png device records an empty display list.
  # recordPlot() still returns a "recordedplot", but replaying it draws nothing
  # and writes no file -- so .vista_record_plot() must reject it.
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  plot(1:10, 1:10)
  empty <- VISTA:::.vista_record_plot()
  grDevices::dev.off()
  unlink(f)

  expect_null(empty)
})
