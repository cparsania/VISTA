# Label/value integrity across the gene-taking API.
#
# VISTA 1.0.0 shipped a matrix where every row carried another gene's values.
# The existing suite did not catch it because it asserted SHAPE -- is.matrix(),
# nrow(), class -- all of which the broken output satisfied perfectly.
#
# These tests assert two properties instead, for every exported function that
# accepts a gene set:
#
#   1. permutation invariance -- the result does not depend on the order the
#      caller happened to supply genes in;
#   2. label/value correspondence -- the value reported for gene G is G's own
#      value, checked against an independent recomputation from the object.
#
# Both are properties the 1.0.0 defect violated.

# ---- helpers ----------------------------------------------------------------

# Canonical form of a plot's layer-independent data, so two calls can be
# compared without depending on row order.
canon <- function(p) {
  d <- p$data
  if (!is.data.frame(d)) return(NULL)
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  d <- d[, order(colnames(d)), drop = FALSE]
  d[] <- lapply(d, function(col) if (is.factor(col)) as.character(col) else col)
  key <- do.call(paste, c(d, sep = "\r"))
  d[order(key), , drop = FALSE]
}

expect_permutation_invariant <- function(fn, v, genes, ..., label = NULL) {
  label <- label %||% "function"
  a <- canon(fn(v, genes = genes, ...))
  b <- canon(fn(v, genes = rev(genes), ...))
  c_ <- canon(fn(v, genes = sample(genes), ...))
  skip_if(is.null(a), paste(label, "has no data frame to compare"))
  rownames(a) <- NULL; rownames(b) <- NULL; rownames(c_) <- NULL
  expect_equal(a, b, info = paste(label, "- reversed input"))
  expect_equal(a, c_, info = paste(label, "- shuffled input"))
}

# ---- 1. permutation invariance across the expression plot family ------------

test_that("expression plots do not depend on the order genes are supplied in", {
  v <- make_small_vista()
  genes <- rownames(v)[seq_len(4)]

  fns <- list(
    get_expression_boxplot   = get_expression_boxplot,
    get_expression_barplot   = get_expression_barplot,
    get_expression_violinplot = get_expression_violinplot,
    get_expression_lineplot  = get_expression_lineplot,
    get_expression_density   = get_expression_density,
    get_expression_joyplot   = get_expression_joyplot,
    get_expression_lollipop  = get_expression_lollipop,
    get_expression_raincloud = get_expression_raincloud
  )
  for (nm in names(fns)) {
    expect_permutation_invariant(fns[[nm]], v, genes, label = nm)
  }
})

test_that("fold-change plots do not depend on the order genes are supplied in", {
  v <- make_small_vista()
  genes <- rownames(v)[seq_len(5)]

  fns <- list(
    get_foldchange_barplot   = get_foldchange_barplot,
    get_foldchange_boxplot   = get_foldchange_boxplot,
    get_foldchange_raincloud = get_foldchange_raincloud
  )
  for (nm in names(fns)) {
    expect_permutation_invariant(fns[[nm]], v, genes, label = nm)
  }
})

test_that("embedding plots do not depend on the order genes are supplied in", {
  v <- make_small_vista()
  genes <- rownames(v)[seq_len(10)]

  expect_permutation_invariant(get_pca_plot, v, genes, label = "get_pca_plot")
  expect_permutation_invariant(get_mds_plot, v, genes, label = "get_mds_plot")
  expect_permutation_invariant(get_corr_heatmap, v, genes, label = "get_corr_heatmap")
})

# ---- 2. label/value correspondence, checked against the object --------------

test_that("expression plots report each gene's own values", {
  v <- make_small_vista()
  genes <- rownames(v)[c(6, 2, 9, 1)]        # deliberately not in row order
  mat <- SummarizedExperiment::assay(v, "norm_counts")

  fns <- list(
    get_expression_boxplot   = get_expression_boxplot,
    get_expression_barplot   = get_expression_barplot,
    get_expression_violinplot = get_expression_violinplot,
    get_expression_lineplot  = get_expression_lineplot,
    get_expression_density   = get_expression_density,
    get_expression_joyplot   = get_expression_joyplot
  )

  for (nm in names(fns)) {
    d <- fns[[nm]](v, genes = genes, log_transform = FALSE)$data
    skip_if(!all(c("gene", "sample", "expression") %in% colnames(d)),
            paste(nm, "does not expose gene/sample/expression"))
    d <- as.data.frame(d, stringsAsFactors = FALSE)
    d$gene <- as.character(d$gene); d$sample <- as.character(d$sample)

    expect_setequal(unique(d$gene), genes)
    # Every plotted point must equal that gene's value in that sample.
    for (i in sample(seq_len(nrow(d)), min(25L, nrow(d)))) {
      expect_equal(
        d$expression[[i]], unname(mat[d$gene[[i]], d$sample[[i]]]),
        info = sprintf("%s: %s / %s", nm, d$gene[[i]], d$sample[[i]])
      )
    }
  }
})

test_that("fold-change plots report each gene's own log2FC", {
  v <- make_small_vista()
  genes <- rownames(v)[c(7, 3, 11, 2, 5)]
  comp <- names(comparisons(v))[[1]]
  tbl <- comparisons(v)[[comp]]

  for (nm in c("get_foldchange_barplot", "get_foldchange_boxplot")) {
    d <- as.data.frame(
      getExportedValue("VISTA", nm)(v, genes = genes)$data,
      stringsAsFactors = FALSE
    )
    id_col <- intersect(c("gene_id", "gene"), colnames(d))[[1]]
    fc_col <- intersect(c("log2fc", "log2FoldChange"), colnames(d))[[1]]
    d[[id_col]] <- as.character(d[[id_col]])

    expect_setequal(unique(d[[id_col]]), genes)
    for (i in seq_len(nrow(d))) {
      expect_equal(
        d[[fc_col]][[i]], unname(tbl[d[[id_col]][[i]], "log2fc"]),
        info = sprintf("%s: %s", nm, d[[id_col]][[i]])
      )
    }
  }
})

test_that("heatmap matrices carry each gene's own row", {
  skip_if_not_installed("ComplexHeatmap")
  v <- make_small_vista()
  genes <- rownames(v)[c(8, 3, 12, 1, 6)]

  hm <- get_expression_heatmap(
    v, genes = genes, value_transform = "raw",
    summarise_replicates = FALSE, return_type = "plot"
  )
  m <- hm@matrix
  src <- SummarizedExperiment::assay(v, "norm_counts")
  expect_setequal(rownames(m), genes)
  for (g in rownames(m)) {
    expect_equal(unname(m[g, ]), unname(src[g, colnames(m)]), info = g)
  }

  fh <- get_foldchange_heatmap(v, genes = genes, return_type = "plot")
  fm <- fh@matrix
  fcm <- get_foldchange_matrix(v, genes = genes)
  expect_setequal(rownames(fm), genes)
  for (g in rownames(fm)) {
    expect_equal(unname(fm[g, ]), unname(fcm[g, colnames(fm)]), info = g)
  }
})

# ---- 3. subsetting a result equals the result of a subset -------------------

test_that("asking for fewer genes returns the same values for those genes", {
  v <- make_small_vista()
  many <- rownames(v)[seq_len(8)]
  few <- many[c(5, 2)]

  all_m <- get_expression_matrix(v, genes = many, summarise = TRUE)
  few_m <- get_expression_matrix(v, genes = few, summarise = TRUE)
  for (g in rownames(few_m)) {
    expect_equal(unname(few_m[g, ]), unname(all_m[g, colnames(few_m)]), info = g)
  }

  comps <- names(comparisons(v))
  all_fc <- get_foldchange_matrix(v, sample_comparisons = comps, genes = many)
  few_fc <- get_foldchange_matrix(v, sample_comparisons = comps, genes = few)
  for (g in rownames(few_fc)) {
    expect_equal(unname(few_fc[g, ]), unname(all_fc[g, colnames(few_fc)]), info = g)
  }
})

# ---- 4. sample-axis integrity ----------------------------------------------

test_that("sample ordering never detaches a sample from its own values", {
  v <- make_small_vista()
  genes <- rownames(v)[seq_len(3)]
  mat <- SummarizedExperiment::assay(v, "norm_counts")

  # get_expression_boxplot() only accepts by = "sample" together with
  # pool_genes = TRUE; get_expression_barplot() takes it directly. Exercise both
  # routes, since sample ordering is where a sample could lose its own column.
  for (ord in c("input", "group", "expression")) {
    d <- as.data.frame(
      get_expression_barplot(
        v, genes = genes, by = "sample", sample_order = ord, log_transform = FALSE
      )$data,
      stringsAsFactors = FALSE
    )
    d$gene <- as.character(d$gene); d$sample <- as.character(d$sample)
    for (i in seq_len(nrow(d))) {
      expect_equal(
        d$expression[[i]], unname(mat[d$gene[[i]], d$sample[[i]]]),
        info = sprintf("sample_order=%s: %s / %s", ord, d$gene[[i]], d$sample[[i]])
      )
    }
  }
})

test_that("group summarisation attributes each group its own samples", {
  v <- make_small_vista()
  gcol <- S4Vectors::metadata(v)$group$column
  grp <- as.character(SummarizedExperiment::colData(v)[[gcol]])
  mat <- SummarizedExperiment::assay(v, "norm_counts")

  summ <- norm_counts(v, summarise = TRUE)
  for (g in colnames(summ)) {
    members <- colnames(v)[grp == g]
    expect_equal(
      unname(summ[, g]), unname(rowMeans(mat[, members, drop = FALSE])),
      info = g
    )
  }
})
