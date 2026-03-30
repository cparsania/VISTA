test_that("example_vista default uses the precomputed object", {
  vista <- example_vista()

  expect_s4_class(vista, "VISTA")
  expect_identical(vista, get("vista_example_default", envir = asNamespace("VISTA")))
})

test_that("example_vista rebuilds for non-default arguments", {
  vista_small <- example_vista(n_genes = 50)
  vista_both <- example_vista(method = "both")

  expect_s4_class(vista_small, "VISTA")
  expect_lte(nrow(vista_small), 50)
  expect_lt(nrow(vista_small), nrow(example_vista()))
  expect_s4_class(vista_both, "VISTA")
  expect_true("de_results_by_method" %in% names(S4Vectors::metadata(vista_both)))
})
