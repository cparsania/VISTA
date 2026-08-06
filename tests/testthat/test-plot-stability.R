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

# PCA and MDS axes come from an eigendecomposition, and the sign of an
# eigenvector is arbitrary: prcomp() and cmdscale() may return v or -v for the
# same input, and which one you get depends on the LAPACK/BLAS build. A runner
# image update flipped both MDS axes and turned this suite red without anything
# in VISTA changing. Magnitudes are identical, so canonicalise each axis to a
# fixed orientation before comparing -- a real change still moves the numbers.
canon_axis <- function(v) {
  if (!is.numeric(v) || !length(v)) return(v)
  i <- which.max(abs(v))                       # largest coordinate is unambiguous
  if (length(i) == 1L && !is.na(v[i]) && v[i] < 0) -v else v
}

embedding_digest_df <- function(d) {
  for (ax in intersect(c("x", "y"), names(d))) d[[ax]] <- canon_axis(d[[ax]])
  keep <- intersect(c("x", "y", "fill", "colour", "group"), names(d))
  d <- d[, keep, drop = FALSE]
  num <- vapply(d, is.numeric, logical(1))
  d[num] <- lapply(d[num], function(v) round(v, 6))
  # sort after canonicalising, so row order does not inherit the flip either
  d <- d[order(do.call(paste, c(d, sep = "|"))), , drop = FALSE]
  rownames(d) <- NULL
  d
}

embedding_digest <- function(p, layer = 1L) {
  embedding_digest_df(ggplot2::ggplot_build(p)$data[[layer]])
}

test_that("PCA layer data is stable", {
  v <- make_small_vista()
  expect_snapshot_value(embedding_digest(get_pca_plot(v)), style = "json2", tolerance = 1e-6)
})

test_that("MDS layer data is stable", {
  v <- make_small_vista()
  expect_snapshot_value(embedding_digest(get_mds_plot(v)), style = "json2", tolerance = 1e-6)
})

test_that("embedding digests survive an axis sign flip", {
  # The property the snapshots depend on: negating an eigenvector -- exactly
  # what a different BLAS build does -- must not change the digest. This is
  # what turned the suite red with nothing in VISTA changed.
  v <- make_small_vista()

  for (nm in c("get_pca_plot", "get_mds_plot")) {
    d <- ggplot2::ggplot_build(getExportedValue("VISTA", nm)(v))$data[[1]]
    flipped <- d
    flipped$x <- -flipped$x
    flipped$y <- -flipped$y

    expect_equal(embedding_digest_df(d), embedding_digest_df(flipped), info = nm)

    # ...but a real change to a coordinate must still be caught, or the
    # canonicalisation would have made the snapshot meaningless.
    moved <- d
    moved$x[[1]] <- moved$x[[1]] + 1
    expect_false(
      isTRUE(all.equal(embedding_digest_df(d), embedding_digest_df(moved))),
      info = nm
    )
  }
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
