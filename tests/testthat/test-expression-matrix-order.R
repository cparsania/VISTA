# get_expression_matrix() subsets rows in assay order but, before 1.2.0, relabelled
# them in the caller's argument order when summarise = TRUE. Every row then carried
# a different gene's identifier.

test_that("summarised matrix rows keep their own identity regardless of gene order", {
  v <- make_small_vista()
  g <- rownames(v)[seq_len(5)]
  g_rev <- rev(g)

  truth <- norm_counts(v, summarise = TRUE)

  # Requesting the genes in reverse must not change which values sit on which row.
  m <- get_expression_matrix(v, genes = g_rev, summarise = TRUE)
  expect_equal(m[rownames(m), colnames(m)], truth[rownames(m), colnames(m)])

  # And spot-check one gene against an independently computed group mean.
  grp <- as.character(sample_info(v)[[S4Vectors::metadata(v)$group$column]])
  ctrl <- colnames(v)[grp == "control"]
  expected <- mean(SummarizedExperiment::assay(v, "norm_counts")[g_rev[[1]], ctrl])
  expect_equal(unname(m[g_rev[[1]], "control"]), expected)
})

test_that("row order and labels agree for every gene ordering", {
  v <- make_small_vista()
  g <- rownames(v)[c(9, 2, 7, 1, 4)]

  m <- get_expression_matrix(v, genes = g, summarise = TRUE)
  truth <- norm_counts(v, summarise = TRUE)

  expect_setequal(rownames(m), g)
  for (gene in g) {
    expect_equal(
      unname(m[gene, ]), unname(truth[gene, colnames(m)]),
      info = gene
    )
  }
})

test_that("duplicated and missing genes are handled without a dimnames error", {
  v <- make_small_vista()
  g <- rownames(v)[seq_len(3)]

  # Duplicates previously made the replacement rowname vector longer than the
  # matrix, producing "length of 'dimnames' not equal to array extent".
  m_dup <- get_expression_matrix(v, genes = c(g, g), summarise = TRUE)
  expect_equal(nrow(m_dup), length(g))
  expect_setequal(rownames(m_dup), g)

  m_missing <- get_expression_matrix(
    v, genes = c(g, "NOT_A_REAL_GENE"), summarise = TRUE
  )
  expect_setequal(rownames(m_missing), g)
})

test_that("the unsummarised path is unaffected and still returns assay order", {
  v <- make_small_vista()
  g <- rev(rownames(v)[seq_len(5)])

  m <- get_expression_matrix(v, genes = g, summarise = FALSE)
  expect_identical(rownames(m), intersect(rownames(v), g))
  expect_identical(colnames(m), colnames(norm_counts(v)))
  expect_equal(
    m, SummarizedExperiment::assay(v, "norm_counts")[rownames(m), , drop = FALSE]
  )
})

test_that("transforms are applied to correctly-labelled rows", {
  v <- make_small_vista()
  g <- rev(rownames(v)[seq_len(4)])

  raw <- get_expression_matrix(v, genes = g, summarise = TRUE)
  logged <- get_expression_matrix(v, genes = g, summarise = TRUE, transform = "log2")
  expect_equal(logged, log2(raw + 1))

  z <- get_expression_matrix(v, genes = g, summarise = TRUE, transform = "zscore")
  expect_identical(rownames(z), rownames(raw))
})
