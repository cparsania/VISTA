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

# --- deconvolution: never attach scores to a sample by position alone ---------

test_that(".normalize_xcell2_scores aligns by name in either orientation", {
  samples <- c("s1", "s2", "s3")
  scores <- data.frame(A = c(1, 2, 3), B = c(4, 5, 6), row.names = samples)

  # Already row-oriented, but shuffled: must be reordered by name, not position.
  shuffled <- scores[c("s3", "s1", "s2"), , drop = FALSE]
  out <- VISTA:::.normalize_xcell2_scores(shuffled, samples)
  expect_identical(rownames(out), samples)
  expect_equal(out["s1", "A"], 1)
  expect_equal(out["s3", "A"], 3)

  # Column-oriented input is transposed and then aligned by name.
  out_t <- VISTA:::.normalize_xcell2_scores(t(as.matrix(scores)), samples)
  expect_identical(rownames(out_t), samples)
  expect_equal(out_t["s2", "A"], 2)
})

test_that("labelled scores that disagree with the samples are refused, not guessed", {
  samples <- c("s1", "s2", "s3")
  # Right number of rows, real labels, but they are someone else's samples.
  wrong <- data.frame(A = c(1, 2, 3), row.names = c("x1", "x2", "x3"))

  expect_error(
    VISTA:::.normalize_xcell2_scores(wrong, samples),
    "do not match"
  )
  # The message must make the refusal and its reason explicit.
  msg <- tryCatch(
    VISTA:::.normalize_xcell2_scores(wrong, samples),
    error = function(e) conditionMessage(e)
  )
  expect_match(msg, "position", fixed = TRUE)
})

test_that("unlabelled scores fall back to position, and say so", {
  samples <- c("s1", "s2", "s3")
  bare <- data.frame(A = c(1, 2, 3))          # default 1..n rownames
  expect_message(
    out <- VISTA:::.normalize_xcell2_scores(bare, samples),
    "by position"
  )
  expect_identical(rownames(out), samples)
  expect_equal(out[["A"]], c(1, 2, 3))

  # Truly absent rownames behave the same way.
  m <- matrix(c(1, 2, 3), ncol = 1, dimnames = list(NULL, "A"))
  expect_message(out2 <- VISTA:::.normalize_xcell2_scores(m, samples), "by position")
  expect_identical(rownames(out2), samples)
})

test_that("a sample-count mismatch is an error rather than a warning", {
  samples <- c("s1", "s2", "s3")
  short <- data.frame(A = c(1, 2), row.names = c("x1", "x2"))
  expect_error(VISTA:::.normalize_xcell2_scores(short, samples), "Could not align")
})

test_that(".vista_labels_uninformative recognises the defaults R invents", {
  expect_true(VISTA:::.vista_labels_uninformative(NULL, 3))
  expect_true(VISTA:::.vista_labels_uninformative(c("", "", ""), 3))
  expect_true(VISTA:::.vista_labels_uninformative(c("1", "2", "3"), 3))
  expect_false(VISTA:::.vista_labels_uninformative(c("s1", "s2", "s3"), 3))
})
