# .align_de_to_counts() is the machinery that pads and reorders an arbitrary DE
# table onto the counts matrix. Every downstream plot reads the aligned table and
# the counts matrix in parallel, so a fault here reproduces the 1.0.0 class of
# defect -- values under the wrong gene labels -- across the whole package.
#
# These tests cover the three things that can go wrong: picking the wrong column
# as the identifier, losing type or content while NA-padding, and emitting rows
# whose position no longer matches the reference.

align <- function(...) VISTA:::.align_de_to_counts(...)

de_tbl <- function(ids, ...) {
  data.frame(
    gene_id = ids, log2fc = seq_along(ids), ...,
    stringsAsFactors = FALSE
  )
}

# ---- identifier resolution --------------------------------------------------

test_that("gene identifiers are taken from id_col, gene_id, rownames, then column 1", {
  ref <- c("g1", "g2")

  # explicit id_col wins
  d <- data.frame(my_id = ref, gene_id = c("wrong1", "wrong2"), log2fc = c(1, 2),
                  stringsAsFactors = FALSE)
  out <- align(d, ref, id_col = "my_id")
  expect_identical(rownames(out), ref)
  expect_identical(out$gene_id, ref)

  # gene_id when no id_col
  expect_identical(rownames(align(de_tbl(ref), ref)), ref)

  # rownames when neither
  d <- data.frame(log2fc = c(1, 2), row.names = ref)
  expect_identical(rownames(align(d, ref)), ref)

  # first character column when there are no usable rownames
  d <- data.frame(id = ref, log2fc = c(1, 2), stringsAsFactors = FALSE)
  expect_identical(rownames(align(d, ref)), ref)

  # an unusable id_col falls through rather than erroring
  expect_identical(rownames(align(de_tbl(ref), ref, id_col = "absent")), ref)
})

test_that("default integer rownames are not mistaken for gene identifiers", {
  # data.frame() rownames are "1","2",... -- using them as gene IDs would align
  # the table to positions instead of genes.
  d <- data.frame(sym = c("g1", "g2"), log2fc = c(1, 2), stringsAsFactors = FALSE)
  expect_identical(rownames(d), c("1", "2"))
  out <- align(d, c("g1", "g2"))
  expect_identical(out$gene_id, c("g1", "g2"))
})

test_that("unusable identifiers are rejected rather than guessed at", {
  ref <- c("g1", "g2")

  # no character column and no meaningful rownames
  expect_error(align(data.frame(a = c(1, 2), b = c(3, 4)), ref),
               "Could not determine gene identifiers")
  # duplicated / empty / NA identifiers
  expect_error(align(de_tbl(c("g1", "g1")), ref), "must be unique")
  expect_error(align(de_tbl(c("g1", "")), ref), "non-empty")
  expect_error(align(de_tbl(c("g1", NA)), ref), "non-NA|unique")
  # first-column fallback must not accept a duplicated or NA column
  expect_error(align(data.frame(id = c("a", "a"), log2fc = c(1, 2),
                                stringsAsFactors = FALSE), ref),
               "Could not determine gene identifiers")
})

test_that("input must be a data.frame-like object", {
  expect_error(align(1:3, "g1"), "must be a data.frame")
  expect_error(align(matrix(1, 2, 2), c("g1", "g2")), "must be a data.frame")

  # S4Vectors::DataFrame is coerced
  d <- S4Vectors::DataFrame(gene_id = c("g1", "g2"), log2fc = c(1, 2))
  out <- align(d, c("g1", "g2"))
  expect_s3_class(out, "data.frame")
  expect_identical(rownames(out), c("g1", "g2"))

  # A tibble ignores `rownames<-`, so aligning one in place would read every
  # gene as absent and replace it with an NA row -- right labels, no data.
  skip_if_not_installed("tibble")
  tb <- tibble::tibble(gene_id = c("g1", "g2"), log2fc = c(10, 20))
  out <- align(tb, c("g1", "g2", "g3"), warn_missing = FALSE)
  expect_identical(rownames(out), c("g1", "g2", "g3"))
  expect_identical(out$gene_id, c("g1", "g2", "g3"))
  expect_equal(out$log2fc, c(10, 20, NA_real_))

  # ...including when the identifier has to come from the first column
  out <- align(tibble::tibble(sym = c("g1", "g2"), log2fc = c(5, 6)), c("g1", "g2"))
  expect_identical(out$gene_id, c("g1", "g2"))
  expect_equal(out$log2fc, c(5, 6))

  # a grouped tibble is the realistic dplyr hand-off
  g <- dplyr::group_by(
    tibble::tibble(gene_id = c("g1", "g2"), log2fc = c(1, 2), grp = c("a", "b")), grp
  )
  expect_equal(align(g, c("g1", "g2"))$log2fc, c(1, 2))

  # data.frame subclasses that do carry rownames still work
  d <- data.frame(gene_id = c("g1", "g2"), log2fc = c(1, 2), stringsAsFactors = FALSE)
  class(d) <- c("my_de_table", "data.frame")
  out <- align(d, c("g2", "g1"))
  expect_identical(rownames(out), c("g2", "g1"))
  expect_equal(out$log2fc, c(2, 1))
})

# ---- reference identifiers --------------------------------------------------

test_that("a duplicated reference is rejected instead of silently de-duplicated", {
  # R's df[c("g2","g2"), ] renames the second row "g2.1", leaving rownames
  # disagreeing with the gene_id column they were built from.
  expect_error(align(de_tbl(c("g1", "g2")), c("g1", "g2", "g2")), "must be unique")
  expect_error(align(de_tbl("g1"), character()), "at least one")
  expect_error(align(de_tbl("g1"), c("g1", NA)), "NA or empty")
  expect_error(align(de_tbl("g1"), c("g1", "")), "NA or empty")
})

# ---- column canonicalisation ------------------------------------------------

test_that("backend-specific column names are mapped to the canonical schema", {
  ref <- c("g1", "g2")

  edger <- data.frame(gene_id = ref, logFC = c(1, -1), PValue = c(.01, .2),
                      FDR = c(.02, .3), LR = c(5, 1), stringsAsFactors = FALSE)
  out <- align(edger, ref)
  expect_true(all(c("log2fc", "pvalue", "padj", "stat") %in% names(out)))
  expect_identical(out$log2fc, c(1, -1))
  expect_identical(out$padj, c(.02, .3))

  limma <- data.frame(gene_id = ref, logFC = c(2, -2), P.Value = c(.03, .4),
                      adj.P.Val = c(.05, .5), AveExpr = c(7, 8), t = c(3, -3),
                      stringsAsFactors = FALSE)
  out <- align(limma, ref)
  expect_identical(out$log2fc, c(2, -2))
  expect_identical(out$pvalue, c(.03, .4))
  expect_identical(out$baseMean, c(7, 8))
  expect_identical(out$stat, c(3, -3))

  deseq <- data.frame(gene_id = ref, log2FoldChange = c(1, 2), lfcSE = c(.1, .2),
                      stat = c(9, 9), pvalue = c(.01, .02), padj = c(.03, .04),
                      baseMean = c(50, 60), stringsAsFactors = FALSE)
  out <- align(deseq, ref)
  expect_identical(
    head(names(out), 7),
    c("gene_id", "baseMean", "log2fc", "lfcSE", "stat", "pvalue", "padj")
  )
})

test_that("renaming never collides with a column that already holds the canonical name", {
  # A table carrying both spellings must keep the canonical one and not end up
  # with two columns of the same name -- df$log2fc would then be ambiguous.
  ref <- c("g1", "g2")
  d <- data.frame(gene_id = ref, log2fc = c(1, 2), logFC = c(99, 99),
                  padj = c(.1, .2), FDR = c(.9, .9), stringsAsFactors = FALSE)
  out <- align(d, ref)
  expect_false(anyDuplicated(names(out)) > 0)
  expect_identical(out$log2fc, c(1, 2))   # not the logFC decoys
  expect_identical(out$padj, c(.1, .2))
})

test_that("non-canonical columns are retained after the canonical ones", {
  ref <- c("g1", "g2")
  d <- data.frame(gene_id = ref, log2fc = c(1, 2), custom = c("a", "b"),
                  stringsAsFactors = FALSE)
  out <- align(d, ref)
  expect_identical(names(out), c("gene_id", "log2fc", "custom"))
  expect_identical(out$custom, c("a", "b"))
})

# ---- NA padding -------------------------------------------------------------

test_that("padded rows keep every column's type", {
  # rbind() of mismatched types silently coerces whole columns -- an integer
  # count column becoming character would break every numeric consumer.
  d <- data.frame(
    gene_id = c("g1", "g2"), num = c(1.5, 2.5), int = c(1L, 2L),
    chr = c("a", "b"), lgl = c(TRUE, FALSE),
    fct = factor(c("x", "y"), levels = c("x", "y", "z")),
    stringsAsFactors = FALSE
  )
  d$date <- as.Date(c("2020-01-01", "2020-01-02"))
  d$time <- as.POSIXct(c("2020-01-01 10:00", "2020-01-02 10:00"), tz = "UTC")

  out <- align(d, c("g2", "g3", "g1"), warn_missing = FALSE)

  expect_identical(rownames(out), c("g2", "g3", "g1"))
  expect_type(out$num, "double")
  expect_type(out$int, "integer")
  expect_type(out$chr, "character")
  expect_type(out$lgl, "logical")
  expect_s3_class(out$fct, "factor")
  expect_identical(levels(out$fct), c("x", "y", "z"))
  expect_s3_class(out$date, "Date")
  expect_s3_class(out$time, "POSIXct")

  # the padded row is NA everywhere except its identifier
  pad <- out["g3", setdiff(names(out), "gene_id")]
  expect_true(all(vapply(pad, is.na, logical(1))))
  expect_identical(out["g3", "gene_id"], "g3")

  # and the real rows kept their own values
  expect_identical(out["g1", "chr"], "a")
  expect_identical(out["g2", "int"], 2L)
})

test_that("strict = TRUE refuses to pad", {
  expect_error(align(de_tbl("g1"), c("g1", "g2"), strict = TRUE),
               "Missing 1 gene")
  # nothing missing means nothing to complain about
  expect_no_error(align(de_tbl(c("g1", "g2")), c("g1", "g2"), strict = TRUE))
})

test_that("padding warns unless the caller opts out", {
  expect_warning(align(de_tbl("g1"), c("g1", "g2")), "Added 1 NA row")
  expect_no_warning(align(de_tbl("g1"), c("g1", "g2"), warn_missing = FALSE))
  expect_no_warning(align(de_tbl(c("g1", "g2")), c("g1", "g2")))
})

# ---- the alignment guarantee itself -----------------------------------------

test_that("output rows follow the reference exactly, whatever the input order", {
  ref <- paste0("g", 1:6)
  vals <- stats::setNames(c(10, 20, 30, 40, 50, 60), ref)

  for (perm in list(ref, rev(ref), sample(ref), ref[c(4, 1, 6, 2, 5, 3)])) {
    d <- data.frame(gene_id = perm, log2fc = unname(vals[perm]),
                    stringsAsFactors = FALSE)
    out <- align(d, ref)
    expect_identical(rownames(out), ref)
    expect_identical(out$gene_id, ref)
    # each gene kept its own value through the reordering
    expect_equal(out$log2fc, unname(vals[ref]))
  }
})

test_that("extra genes are dropped and missing genes padded, together", {
  d <- data.frame(gene_id = c("x1", "g2", "g1", "x2"), log2fc = c(-9, 2, 1, -9),
                  stringsAsFactors = FALSE)
  out <- align(d, c("g1", "g2", "g3"), warn_missing = FALSE)

  expect_identical(out$gene_id, c("g1", "g2", "g3"))
  expect_equal(out$log2fc, c(1, 2, NA_real_))
  expect_false(any(c("x1", "x2") %in% rownames(out)))
})

test_that("rownames and the gene_id column always agree", {
  # The invariant the 1.0.0 defect broke: two identifier sources that disagree.
  cases <- list(
    align(de_tbl(c("g2", "g1")), c("g1", "g2")),
    align(de_tbl(c("g1")), c("g1", "g2"), warn_missing = FALSE),
    align(de_tbl(c("g1", "g2", "g3")), c("g2", "g1"))
  )
  for (out in cases) {
    expect_identical(rownames(out), as.character(out$gene_id))
    expect_false(anyNA(rownames(out)))
  }
})

test_that("aligned DE tables satisfy the VISTA rownames contract", {
  # The consumer-level guarantee: what this function returns can be stored on a
  # VISTA object without tripping validation.
  v <- make_small_vista()
  ref <- rownames(v)
  tbl <- comparisons(v)[[1]]

  shuffled <- tbl[sample(nrow(tbl)), , drop = FALSE]
  out <- align(shuffled, ref)

  expect_identical(rownames(out), ref)
  expect_identical(rownames(out), rownames(norm_counts(v)))
  # values travelled with their genes
  expect_equal(out[ref, "log2fc"], tbl[ref, "log2fc"])
})
