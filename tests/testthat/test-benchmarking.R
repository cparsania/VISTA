test_that("benchmark_vista_equivalence matches standalone DESeq2 and edgeR", {
  inputs <- make_benchmark_inputs()

  report <- do.call(
    benchmark_vista_equivalence,
    c(inputs, list(methods = c("deseq2", "edger"), tolerance = 1e-8, return_plots = FALSE))
  )

  expect_true(report$valid)
  expect_true(all(report$comparison_summary$deg_sets_identical))
  expect_true(all(report$comparison_summary$regulation_identical))
  expect_true(all(report$comparison_summary$norm_counts_identical))
  expect_true(all(report$comparison_summary$log2fc_within_tolerance))
  expect_true(all(report$comparison_summary$pvalue_within_tolerance))
  expect_true(all(report$comparison_summary$padj_within_tolerance))
  expect_true(all(report$visual_summary$pass))
})

test_that("benchmark_vista_equivalence matches standalone limma", {
  skip_if_not_installed("limma")
  inputs <- make_benchmark_inputs()

  report <- do.call(
    benchmark_vista_equivalence,
    c(inputs, list(methods = "limma", tolerance = 1e-8, return_plots = FALSE))
  )

  expect_true(report$valid)
  expect_true(all(report$comparison_summary$deg_sets_identical))
  expect_true(all(report$comparison_summary$regulation_identical))
  expect_true(all(report$comparison_summary$norm_counts_identical))
  expect_true(all(report$comparison_summary$log2fc_within_tolerance))
  expect_true(all(report$comparison_summary$pvalue_within_tolerance))
  expect_true(all(report$comparison_summary$padj_within_tolerance))
  expect_true(all(report$visual_summary$pass))
})

test_that("validate_vista_deep returns a passing report on benchmark inputs", {
  inputs <- make_benchmark_inputs()

  report <- do.call(
    validate_vista_deep,
    c(inputs, list(methods = "deseq2", tolerance = 1e-8, error = FALSE))
  )

  expect_true(report$valid)
})
