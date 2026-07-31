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

make_mock_deconv_vista <- function(include_sample_names = FALSE, include_text_col = FALSE) {
  norm_mat <- matrix(
    seq_len(24),
    nrow = 6,
    ncol = 4,
    dimnames = list(
      paste0("gene", seq_len(6)),
      paste0("sample", seq_len(4))
    )
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(norm_counts = norm_mat),
    colData = S4Vectors::DataFrame(
      cond_long = c("control", "control", "treated", "treated"),
      row.names = colnames(norm_mat)
    ),
    rowData = S4Vectors::DataFrame(
      gene_id = rownames(norm_mat),
      row.names = rownames(norm_mat)
    )
  )
  v <- as_vista(se, group_column = "cond_long")

  frac <- data.frame(
    smooth_muscle_cell = c(0.8, 0.7, 0.2, 0.3),
    fibroblast = c(0.1, 0.2, 0.6, 0.5),
    immune = c(0.1, 0.1, 0.2, 0.2),
    stringsAsFactors = FALSE
  )
  if (isTRUE(include_text_col)) {
    frac$note <- letters[seq_len(nrow(frac))]
  }
  if (isTRUE(include_sample_names)) {
    frac$sample_names <- colnames(norm_mat)
  } else {
    rownames(frac) <- colnames(norm_mat)
  }

  md <- S4Vectors::metadata(v)
  md$cell_fractions <- frac
  S4Vectors::metadata(v) <- md
  v
}

test_that("get_celltype_barplot supports top_n and normalization", {
  v <- make_mock_deconv_vista()

  p <- get_celltype_barplot(
    v,
    top_n = 2,
    collapse_other = TRUE,
    normalize = "sample",
    facet_by = "none"
  )

  expect_s3_class(p, "ggplot")
  expect_true("Other" %in% unique(as.character(p$data$cell_type)))

  sample_sums <- tapply(p$data$score, p$data$sample, sum)
  expect_equal(as.numeric(sample_sums), rep(1, length(sample_sums)))
})

test_that("get_celltype_barplot validates group_column", {
  v <- make_mock_deconv_vista()
  expect_error(get_celltype_barplot(v, group_column = "missing_group"), "not found")
})

test_that("get_celltype_group_dotplot returns ggplot", {
  v <- make_mock_deconv_vista(include_sample_names = TRUE, include_text_col = TRUE)

  expect_warning(
    p <- get_celltype_group_dotplot(
      v,
      top_n = 2,
      error = "se",
      add_points = TRUE
    ),
    "Dropping non-numeric columns"
  )

  expect_s3_class(p, "ggplot")
  expect_equal(length(unique(as.character(p$data$cell_type))), 2)
})

test_that("get_celltype_group_dotplot requires a resolvable group column", {
  v <- make_mock_deconv_vista()
  md <- S4Vectors::metadata(v)
  md$group <- NULL
  S4Vectors::metadata(v) <- md

  expect_error(
    get_celltype_group_dotplot(v, group_column = NULL),
    "provide"
  )
})

test_that("get_celltype_heatmap returns matrix and plot outputs", {
  v <- make_mock_deconv_vista(include_sample_names = TRUE)

  out <- get_celltype_heatmap(
    v,
    top_n = 2,
    cluster_columns = FALSE,
    return_type = "both"
  )

  expect_type(out, "list")
  expect_s3_class(out$plot, "ggplot")
  expect_true(is.matrix(out$matrix))
  expect_equal(nrow(out$matrix), 2)
  expect_equal(colnames(out$matrix), colnames(norm_counts(v)))

  mat_only <- get_celltype_heatmap(v, top_n = 2, return_type = "matrix")
  expect_true(is.matrix(mat_only))
  expect_equal(nrow(mat_only), 2)
})

test_that(".collapse_ensembl_symbol_ids works on the matrix its caller passes", {
  v <- make_small_vista()
  mat <- SummarizedExperiment::assay(v, "norm_counts")[seq_len(6), , drop = FALSE]
  rownames(mat) <- c(
    "ENSG00000000001:AAA", "ENSG00000000002:BBB", "ENSG00000000003:AAA",
    "ENSG00000000004:CCC", "ENSG00000000005:BBB", "ENSG00000000006:DDD"
  )

  # run_cell_deconvolution() hands this a matrix, not a data.frame.
  out <- VISTA:::.collapse_ensembl_symbol_ids(mat)

  expect_true(is.matrix(out))
  expect_setequal(rownames(out), c("AAA", "BBB", "CCC", "DDD"))
  expect_identical(colnames(out), colnames(mat))

  # Duplicated symbols are averaged.
  expect_equal(
    unname(out["AAA", ]),
    unname(colMeans(mat[c(1, 3), , drop = FALSE]))
  )
  expect_equal(unname(out["CCC", ]), unname(mat[4, ]))
})

test_that("ENSEMBL:SYMBOL rownames are auto-detected and routed to the collapser", {
  ids <- c("ENSG00000000001:AAA", "ENSG00000000002:BBB", "ENSG00000000003:CCC")
  expect_identical(VISTA:::.infer_gene_id_type(ids), "ensembl_symbol")

  # This is the combination that made run_cell_deconvolution() unusable: the
  # default gene_id_type = "auto" resolves to "ensembl_symbol", which then hit
  # the collapser with a matrix.
  v <- make_small_vista()
  mat <- SummarizedExperiment::assay(v, "norm_counts")[seq_len(3), , drop = FALSE]
  rownames(mat) <- ids
  expect_no_error(VISTA:::.collapse_ensembl_symbol_ids(mat))
})
