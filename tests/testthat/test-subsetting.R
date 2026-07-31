# Before 1.2.0 VISTA inherited SummarizedExperiment's `[` unchanged, so
# subsetting produced an object whose metadata described genes and samples it no
# longer contained -- and validObject() reported it as valid.

test_that("row subsetting reindexes every stored DE table", {
  v <- make_small_vista()
  keep <- rownames(v)[seq_len(10)]
  vs <- v[seq_len(10), ]

  expect_s4_class(vs, "VISTA")
  expect_identical(nrow(vs), 10L)
  expect_identical(rownames(vs), keep)

  for (tbl in comparisons(vs)) {
    expect_identical(nrow(tbl), 10L)
    expect_identical(rownames(tbl), keep)
  }

  # Every by-method table follows too, not just the active one.
  md <- S4Vectors::metadata(vs)
  for (src in names(md$de_results_by_method)) {
    for (tbl in md$de_results_by_method[[src]]) {
      expect_identical(rownames(as.data.frame(tbl)), keep, info = src)
    }
  }

  expect_true(validate_vista(vs, level = "full", error = FALSE)$valid)
})

test_that("row subsetting preserves the requested order and values", {
  v <- make_small_vista()
  idx <- c(9L, 2L, 7L, 1L)
  keep <- rownames(v)[idx]
  vs <- v[idx, ]

  expect_identical(rownames(vs), keep)
  expect_identical(rownames(comparisons(vs)[[1]]), keep)

  # Values must travel with their genes.
  full <- comparisons(v)[[1]]
  sub <- comparisons(vs)[[1]]
  expect_equal(sub$log2fc, full[keep, "log2fc"])

  expect_equal(
    SummarizedExperiment::assay(vs, "norm_counts"),
    SummarizedExperiment::assay(v, "norm_counts")[keep, , drop = FALSE]
  )
})

test_that("logical and character row subscripts behave", {
  v <- make_small_vista()

  keep_lgl <- rowMeans(SummarizedExperiment::assay(v, "norm_counts")) > 500
  skip_if(sum(keep_lgl) < 2, "not enough high-expression genes in the fixture")
  vl <- v[keep_lgl, ]
  expect_identical(nrow(comparisons(vl)[[1]]), sum(keep_lgl))
  expect_true(validate_vista(vl, level = "full", error = FALSE)$valid)

  ids <- rownames(v)[c(3, 5)]
  vc <- v[ids, ]
  expect_identical(rownames(vc), ids)
  expect_identical(rownames(comparisons(vc)[[1]]), ids)
})

test_that("head() and single-gene subsets stay valid", {
  v <- make_small_vista()

  vh <- utils::head(v, 3)
  expect_identical(nrow(vh), 3L)
  expect_identical(nrow(comparisons(vh)[[1]]), 3L)
  expect_true(validate_vista(vh, level = "full", error = FALSE)$valid)

  v1 <- v[1, ]
  expect_identical(nrow(comparisons(v1)[[1]]), 1L)
  expect_true(validate_vista(v1, level = "full", error = FALSE)$valid)
})

test_that("DEG summaries are recounted from the retained rows", {
  v <- make_small_vista()
  vs <- v[seq_len(20), ]

  summ <- deg_summary(vs)[[1]]
  expect_identical(sum(summ$n), 20L)

  tbl <- comparisons(vs)[[1]]
  expected <- as.data.frame(table(tbl$regulation), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(expected))) {
    got <- summ$n[match(expected$Var1[[i]], as.character(summ$regulation))]
    expect_identical(as.integer(got), as.integer(expected$Freq[[i]]))
  }
})

test_that("column subsetting prunes colour maps and does not produce NaN", {
  v <- make_small_vista()
  grp <- as.character(sample_info(v)$cond_long)
  keep <- grp == "control"

  vs <- suppressWarnings(v[, keep])

  expect_identical(ncol(vs), sum(keep))
  expect_setequal(names(group_colors(vs)), "control")
  expect_false("treatment1" %in% names(group_colors(vs)))

  nc <- norm_counts(vs, summarise = TRUE)
  expect_false(anyNA(nc))
  expect_identical(colnames(nc), "control")

  # DE tables are results, not per-sample data, so they are retained intact.
  expect_identical(nrow(comparisons(vs)[[1]]), nrow(v))
})

test_that("dropping a comparison's samples warns", {
  v <- make_small_vista()
  keep <- as.character(sample_info(v)$cond_long) == "control"

  expect_warning(v[, keep], "no longer represented")
})

test_that("subsetting both dimensions works and is recorded in provenance", {
  v <- make_small_vista()
  vs <- suppressWarnings(v[seq_len(8), seq_len(3)])

  expect_identical(dim(vs), c(8L, 3L))
  expect_identical(nrow(comparisons(vs)[[1]]), 8L)

  hist <- S4Vectors::metadata(vs)$provenance$subset_history
  expect_true(length(hist) >= 1L)
  expect_identical(hist[[length(hist)]]$genes, 8L)
  expect_identical(hist[[length(hist)]]$samples, 3L)
})

test_that("a stale-metadata object is now caught by the core check", {
  v <- make_small_vista()

  # Hand-build the corruption `[` used to produce: shrink the assay but leave
  # the DE tables describing every original gene.
  broken <- v
  md <- S4Vectors::metadata(broken)
  broken <- methods::as(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(norm_counts = SummarizedExperiment::assay(v, "norm_counts")[seq_len(5), ]),
      colData = SummarizedExperiment::colData(v),
      rowData = SummarizedExperiment::rowData(v)[seq_len(5), , drop = FALSE]
    ),
    "VISTA"
  )
  S4Vectors::metadata(broken) <- md

  res <- suppressWarnings(validate_vista(broken, level = "core", error = FALSE))
  expect_false(res$valid)
  expect_true(any(grepl("rownames must match", res$issues)))
})

test_that("downstream plots reflect the subset rather than the original", {
  v <- make_small_vista()
  reg <- comparisons(v)[[1]]$regulation

  # Pick a subset that actually contains DEGs, otherwise the barplot has
  # nothing to draw and the assertion would be vacuous.
  de_idx <- which(reg %in% c("Up", "Down"))
  expect_gt(length(de_idx), 0L)
  # Keep one DEG and several non-DEGs, so the plotted total must shrink.
  idx <- sort(unique(c(de_idx[[1]], which(reg == "Other")[seq_len(5)])))

  vs <- v[idx, ]
  expect_identical(nrow(vs), length(idx))

  p <- get_deg_count_barplot(vs)
  expect_s3_class(p, "ggplot")

  # Counts must describe the subset, not the 123 genes of the original.
  sub_reg <- comparisons(vs)[[1]]$regulation
  expect_equal(
    sum(p$data$n),
    sum(sub_reg %in% c("Up", "Down"))
  )
  expect_lt(sum(p$data$n), sum(reg %in% c("Up", "Down")))
})
