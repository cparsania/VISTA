test_that("run_cell_deconvolution requires xCell2", {
  vista <- make_small_vista()
  expect_error(run_cell_deconvolution(vista, method = "MuSiC"), '"xCell2"')
  skip_if_not_installed("xCell2")

  xcell2_exports <- getNamespaceExports("xCell2")
  has_supported_api <- "xCell2Score" %in% xcell2_exports || "xCell2Analysis" %in% xcell2_exports
  if (!has_supported_api) {
    skip("xCell2 does not expose xCell2Score() or xCell2Analysis() in this version.")
  }

  vista_cf <- tryCatch(
    run_cell_deconvolution(
      vista,
      method = "xCell2",
      xcell2_reference = "DICE_demo.xCell2Ref"
    ),
    error = function(e) e
  )
  if (inherits(vista_cf, "error")) {
    skip(paste0("xCell2 runtime failed in this environment: ", conditionMessage(vista_cf)))
  }

  expect_s4_class(vista_cf, "VISTA")
  fractions <- get_cell_fractions(vista_cf)
  expect_s3_class(fractions, "data.frame")
  expect_equal(nrow(fractions), ncol(norm_counts(vista_cf)))
})
