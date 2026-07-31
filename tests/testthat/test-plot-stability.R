# Regression net for the Phase 3 API harmonization.
#
# Phase 3 adds argument aliases and renames but must not move a single plotted
# value. These snapshots pin the layer data of representative plots across the
# main families, taken before any renaming, so an accidental behaviour change
# during the rename shows up as a snapshot diff rather than passing silently.

# Round to keep snapshots stable across BLAS/platform noise while still
# catching real changes.
layer_digest <- function(p, layer = 1L) {
  d <- ggplot2::ggplot_build(p)$data[[layer]]
  keep <- intersect(c("x", "y", "xmin", "xmax", "ymin", "ymax", "fill", "colour", "group"), names(d))
  d <- d[, keep, drop = FALSE]
  num <- vapply(d, is.numeric, logical(1))
  d[num] <- lapply(d[num], function(v) round(v, 6))
  d[order(do.call(paste, c(d, sep = "|"))), , drop = FALSE]
}

test_that("PCA layer data is stable", {
  v <- make_small_vista()
  expect_snapshot_value(layer_digest(get_pca_plot(v)), style = "json2", tolerance = 1e-6)
})

test_that("MDS layer data is stable", {
  v <- make_small_vista()
  expect_snapshot_value(layer_digest(get_mds_plot(v)), style = "json2", tolerance = 1e-6)
})

test_that("correlation heatmap layer data is stable", {
  v <- make_small_vista()
  expect_snapshot_value(
    layer_digest(get_corr_heatmap(v, triangle = "lower")),
    style = "json2", tolerance = 1e-6
  )
})

test_that("DEG count barplot layer data is stable", {
  v <- make_small_vista()
  expect_snapshot_value(layer_digest(get_deg_count_barplot(v)), style = "json2", tolerance = 1e-6)
})

test_that("MA plot layer data is stable", {
  v <- make_small_vista()
  comp <- names(comparisons(v))[[1]]
  expect_snapshot_value(
    layer_digest(get_ma_plot(v, sample_comparison = comp)),
    style = "json2", tolerance = 1e-6
  )
})

test_that("expression boxplot layer data is stable", {
  v <- make_small_vista()
  genes <- rownames(v)[seq_len(3)]
  expect_snapshot_value(
    layer_digest(get_expression_boxplot(v, genes = genes)),
    style = "json2", tolerance = 1e-6
  )
})

test_that("expression barplot layer data is stable", {
  v <- make_small_vista()
  genes <- rownames(v)[seq_len(3)]
  expect_snapshot_value(
    layer_digest(get_expression_barplot(v, genes = genes)),
    style = "json2", tolerance = 1e-6
  )
})

test_that("fold-change barplot layer data is stable", {
  v <- make_small_vista()
  genes <- rownames(v)[seq_len(5)]
  expect_snapshot_value(
    layer_digest(get_foldchange_barplot(v, genes = genes)),
    style = "json2", tolerance = 1e-6
  )
})
