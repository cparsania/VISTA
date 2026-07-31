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

test_that("schema version is compared, not merely present", {
  v <- make_small_vista()
  expect_identical(VISTA:::.vista_schema_compare(v), "current")

  older <- v
  S4Vectors::metadata(older)$vista_schema_version <- "0.0.1"
  expect_identical(VISTA:::.vista_schema_compare(older), "older")

  newer <- v
  S4Vectors::metadata(newer)$vista_schema_version <- "99.0.0"
  expect_identical(VISTA:::.vista_schema_compare(newer), "newer")

  missing <- v
  S4Vectors::metadata(missing)$vista_schema_version <- NULL
  expect_identical(VISTA:::.vista_schema_compare(missing), "unknown")
})

test_that("an older schema informs and a newer schema is a validation issue", {
  v <- make_small_vista()

  older <- v
  S4Vectors::metadata(older)$vista_schema_version <- "0.0.1"
  expect_message(
    res <- validate_vista(older, level = "full", error = FALSE),
    "updateObject"
  )
  expect_true(res$valid)

  newer <- v
  S4Vectors::metadata(newer)$vista_schema_version <- "99.0.0"
  res_new <- suppressWarnings(validate_vista(newer, level = "full", error = FALSE))
  expect_false(res_new$valid)
  expect_true(any(grepl("newer than this version", res_new$issues)))
})

test_that("updateObject migrates an older object and records the migration", {
  v <- make_small_vista()
  older <- v
  S4Vectors::metadata(older)$vista_schema_version <- "0.9.0"
  S4Vectors::metadata(older)$de_cutoffs <- NULL

  expect_message(
    updated <- BiocGenerics::updateObject(older, verbose = TRUE),
    "Updated VISTA schema"
  )

  md <- S4Vectors::metadata(updated)
  expect_identical(md$vista_schema_version, VISTA:::.VISTA_SCHEMA_VERSION)
  expect_false(is.null(md$de_cutoffs))
  expect_true(length(md$provenance$updates) >= 1L)
  expect_identical(md$provenance$updates[[1]]$from, "0.9.0")
  expect_true("de_cutoffs" %in% md$provenance$updates[[1]]$filled)

  expect_true(suppressMessages(validate_vista(updated, level = "full", error = FALSE))$valid)
})

test_that("updateObject refuses to downgrade a newer object", {
  v <- make_small_vista()
  S4Vectors::metadata(v)$vista_schema_version <- "99.0.0"

  expect_warning(out <- BiocGenerics::updateObject(v), "newer than")
  expect_identical(S4Vectors::metadata(out)$vista_schema_version, "99.0.0")
})
