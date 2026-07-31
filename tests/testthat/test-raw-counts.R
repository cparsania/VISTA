# Before 1.2.0 a VISTA object discarded its own input: counts(v) failed and the
# object could not be handed back to DESeq2 or edgeR.

test_that("VISTA objects retain raw counts as a second assay", {
  v <- make_small_vista()

  expect_true("counts" %in% SummarizedExperiment::assayNames(v))
  # norm_counts must stay FIRST so an unqualified assay(x) is unchanged.
  expect_identical(SummarizedExperiment::assayNames(v)[[1]], "norm_counts")
  expect_equal(
    SummarizedExperiment::assay(v),
    SummarizedExperiment::assay(v, "norm_counts")
  )
})

test_that("counts() returns integer-valued raw counts aligned with the object", {
  v <- make_small_vista()
  raw <- counts(v)

  expect_true(is.matrix(raw))
  expect_identical(dim(raw), dim(norm_counts(v)))
  expect_identical(rownames(raw), rownames(v))
  expect_identical(colnames(raw), colnames(v))

  expect_true(all(raw == round(raw)))
  expect_true(all(raw >= 0))

  # Raw and normalized are genuinely different matrices.
  expect_false(isTRUE(all.equal(as.vector(raw), as.vector(norm_counts(v)))))
})

test_that("all three backends retain raw counts", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())
  si <- sample_metadata[sample_metadata$cond_long %in% c("control", "treatment1"), ]
  cnt <- count_data[seq_len(300), c("gene_id", si$sample_names)]

  for (m in c("deseq2", "edger", "limma")) {
    v <- suppressMessages(create_vista(
      counts = cnt, sample_info = si, column_geneid = "gene_id",
      group_column = "cond_long", group_numerator = "treatment1",
      group_denominator = "control", method = m,
      min_counts = 5, min_replicates = 1
    ))
    expect_true("counts" %in% SummarizedExperiment::assayNames(v), info = m)
    expect_identical(dim(counts(v)), dim(norm_counts(v)), info = m)
    expect_true(all(counts(v) == round(counts(v))), info = m)
  }
})

test_that("keep_raw_counts = FALSE omits the assay and counts() explains why", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())
  si <- sample_metadata[sample_metadata$cond_long %in% c("control", "treatment1"), ]
  cnt <- count_data[seq_len(200), c("gene_id", si$sample_names)]

  v <- suppressMessages(create_vista(
    counts = cnt, sample_info = si, column_geneid = "gene_id",
    group_column = "cond_long", group_numerator = "treatment1",
    group_denominator = "control", min_counts = 5, min_replicates = 1,
    keep_raw_counts = FALSE
  ))

  expect_false("counts" %in% SummarizedExperiment::assayNames(v))
  expect_error(counts(v), "does not carry a raw")
  expect_error(counts(v), "keep_raw_counts")
  # The message must not promise updateObject() can fix it.
  expect_error(counts(v), "cannot be inverted")
})

test_that("raw counts survive subsetting in both dimensions", {
  v <- make_small_vista()

  vr <- v[seq_len(10), ]
  expect_identical(dim(counts(vr)), c(10L, ncol(v)))
  expect_identical(rownames(counts(vr)), rownames(vr))

  vc <- suppressWarnings(v[, seq_len(3)])
  expect_identical(dim(counts(vc)), c(nrow(v), 3L))
  expect_identical(colnames(counts(vc)), colnames(vc))
})

test_that("as_vista accepts raw counts and validates their dimensions", {
  mat <- matrix(as.numeric(rpois(60, 20)), nrow = 10)
  rownames(mat) <- paste0("gene", seq_len(10))
  colnames(mat) <- paste0("sample", seq_len(6))
  raw <- mat + 1

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(norm_counts = mat),
    colData = S4Vectors::DataFrame(
      cond = rep(c("A", "B"), each = 3), row.names = colnames(mat)
    ),
    rowData = S4Vectors::DataFrame(
      gene_id = rownames(mat), row.names = rownames(mat)
    )
  )

  v <- as_vista(se, group_column = "cond", raw_counts = raw)
  expect_true("counts" %in% SummarizedExperiment::assayNames(v))
  expect_equal(counts(v), raw)

  # Reordered rows are realigned rather than silently mismatched.
  v2 <- as_vista(se, group_column = "cond", raw_counts = raw[rev(rownames(raw)), ])
  expect_equal(counts(v2), raw)

  expect_error(
    as_vista(se, group_column = "cond", raw_counts = raw[seq_len(3), ]),
    "same genes and samples"
  )
})

test_that("as_deseq_dataset round-trips a VISTA object back into DESeq2", {
  skip_if_not_installed("DESeq2")
  v <- make_small_vista()

  dds <- as_deseq_dataset(v)
  expect_s4_class(dds, "DESeqDataSet")
  expect_identical(dim(dds), dim(v))
  expect_equal(DESeq2::counts(dds), round(counts(v)), ignore_attr = "storage.mode")

  # Default design uses the stored grouping column.
  expect_identical(
    as.character(DESeq2::design(dds)),
    as.character(stats::as.formula(paste("~", S4Vectors::metadata(v)$group$column)))
  )

  # An explicit design is honoured.
  dds2 <- as_deseq_dataset(v, design = ~ cell + cond_long)
  expect_s4_class(dds2, "DESeqDataSet")
})

test_that("the schema version records the new assay layout", {
  v <- make_small_vista()
  expect_identical(
    S4Vectors::metadata(v)$vista_schema_version,
    VISTA:::.VISTA_SCHEMA_VERSION
  )
  expect_identical(VISTA:::.vista_schema_compare(v), "current")
})
