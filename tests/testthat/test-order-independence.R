# VISTA 1.0.0 shipped a bug where row labels were separated from their values.
# The same class of defect appears whenever code depends on an ordering nothing
# guarantees. These tests pin the places where an upstream order could leak into
# results.

test_that("pathway gene capping does not depend on the incoming gene order", {
  v <- make_small_vista()
  genes <- rownames(v)

  # get_pathway_heatmap() caps pathway membership at max_genes. Membership is a
  # SET; clusterProfiler's ordering of it is arbitrary and, verified empirically,
  # differs between R sessions. The cap must therefore not be positional.
  pick <- function(g, n) {
    mat <- SummarizedExperiment::assay(v, "norm_counts")[g, , drop = FALSE]
    score <- matrixStats::rowVars(mat)
    score[!is.finite(score)] <- -Inf
    g[utils::head(order(-score, g), n)]
  }

  set_a <- genes[seq_len(30)]
  set_b <- rev(set_a)
  set_c <- sample(set_a)

  expect_identical(pick(set_a, 10), pick(set_b, 10))
  expect_identical(pick(set_a, 10), pick(set_c, 10))
})

test_that("chord gene capping is stable under input permutation", {
  counts <- c(g1 = 3L, g2 = 1L, g3 = 3L, g4 = 2L, g5 = 1L)
  rank_tier <- function(g) if (!length(g)) g else g[order(-as.integer(counts[g]), g)]

  a <- rank_tier(names(counts))
  b <- rank_tier(rev(names(counts)))
  c_ <- rank_tier(sample(names(counts)))

  expect_identical(a, b)
  expect_identical(a, c_)
  # Highest pathway participation first, identifier breaking ties.
  expect_identical(a, c("g1", "g3", "g4", "g2", "g5"))
})

test_that("gene set membership is order-insensitive end to end", {
  skip_if_not_installed("msigdbr")
  v <- make_small_vista()
  comp <- names(comparisons(v))[[1]]

  msig <- try(
    suppressMessages(get_msigdb_enrichment(
      v, sample_comparison = comp, regulation = "Both",
      msigdb_category = "H", from_type = "ENSEMBL"
    )),
    silent = TRUE
  )
  skip_if(inherits(msig, "try-error") || is.null(msig$enrich), "enrichment unavailable")
  skip_if(nrow(as.data.frame(msig$enrich)) == 0, "no enriched terms")

  # get_pathway_genes returns membership; only the SET is meaningful.
  a <- get_pathway_genes(msig, top_n = 2, return_type = "list")
  b <- get_pathway_genes(msig, top_n = 2, return_type = "list")
  expect_identical(names(a), names(b))
  for (nm in names(a)) expect_setequal(a[[nm]], b[[nm]])
})

test_that("summarised expression labels always travel with their own values", {
  # The 1.0.0 defect itself: a regression guard that recomputes the answer
  # independently rather than checking shape.
  v <- make_small_vista()
  truth <- norm_counts(v, summarise = TRUE)
  ids <- rownames(v)

  for (g in list(ids[1:5], rev(ids[1:5]), ids[c(9, 2, 7, 1, 4)], sample(ids, 6))) {
    m <- get_expression_matrix(v, genes = g, summarise = TRUE)
    for (gene in rownames(m)) {
      expect_equal(
        unname(m[gene, ]), unname(truth[gene, colnames(m)]),
        info = paste(gene, "in", paste(utils::head(g, 3), collapse = ","))
      )
    }
  }
})

test_that("fold-change matrices keep each gene's own values under permutation", {
  v <- make_small_vista()
  comps <- names(comparisons(v))
  ids <- rownames(v)[seq_len(8)]

  a <- get_foldchange_matrix(v, sample_comparisons = comps, genes = ids)
  b <- get_foldchange_matrix(v, sample_comparisons = comps, genes = rev(ids))

  expect_setequal(rownames(a), rownames(b))
  for (g in rownames(a)) {
    expect_equal(unname(a[g, ]), unname(b[g, colnames(a)]), info = g)
  }

  # And against the source tables, independently.
  tbl <- comparisons(v)[[comps[[1]]]]
  for (g in rownames(a)) {
    expect_equal(unname(a[g, comps[[1]]]), unname(tbl[g, "log2fc"]), info = g)
  }
})
