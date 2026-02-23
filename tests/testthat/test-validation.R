test_that("as_vista builds a valid VISTA object", {
  mat <- matrix(rnorm(60), nrow = 10)
  rownames(mat) <- paste0("gene", seq_len(nrow(mat)))
  colnames(mat) <- paste0("sample", seq_len(ncol(mat)))

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(norm_counts = mat),
    colData = S4Vectors::DataFrame(
      cond = rep(c("A", "B"), each = 3),
      row.names = colnames(mat)
    ),
    rowData = S4Vectors::DataFrame(
      gene_id = rownames(mat),
      row.names = rownames(mat)
    )
  )

  v <- as_vista(se, group_column = "cond")
  expect_s4_class(v, "VISTA")

  report <- validate_vista(v, level = "core", error = FALSE)
  expect_true(report$valid)
  expect_true(is.character(S4Vectors::metadata(v)$vista_schema_version))
})

test_that("validate_vista reports non-VISTA objects", {
  expect_warning(
    report <- validate_vista(list(a = 1), error = FALSE),
    "validation reported issues"
  )
  expect_false(report$valid)
  expect_true(any(grepl("inherit from class 'VISTA'", report$issues)))
})

test_that("validate_vista(full) catches inconsistent method metadata", {
  v <- make_small_vista()
  md <- S4Vectors::metadata(v)
  md$de_results_by_method <- list(
    deseq2 = md$de_results
  )
  md$de_summary_by_method <- list(
    deseq2 = md$de_summary
  )
  md$de_active_source <- "consensus"
  S4Vectors::metadata(v) <- md

  expect_warning(
    report <- validate_vista(v, level = "full", error = FALSE),
    "validation reported issues"
  )
  expect_false(report$valid)
  expect_true(any(grepl("de_active_source", report$issues)))
})
