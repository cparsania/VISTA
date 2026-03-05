make_mock_vista <- function(with_comparisons = TRUE) {
  mat <- matrix(rpois(60, lambda = 20), nrow = 10)
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

  if (isTRUE(with_comparisons)) {
    de <- data.frame(
      gene_id = rownames(mat),
      log2fc = rnorm(nrow(mat)),
      pvalue = stats::runif(nrow(mat)),
      padj = stats::runif(nrow(mat)),
      stringsAsFactors = FALSE,
      row.names = rownames(mat)
    )
    de$regulation <- ifelse(de$log2fc > 0.5, "Up", ifelse(de$log2fc < -0.5, "Down", "Other"))

    md <- S4Vectors::metadata(v)
    md$de_results <- S4Vectors::SimpleList(A_vs_B = de)
    md$de_summary <- S4Vectors::SimpleList(A_vs_B = as.data.frame(table(de$regulation), stringsAsFactors = FALSE))
    md$de_cutoffs <- list(log2fc_cutoff = 1, pval_cutoff = 0.05, p_value_type = "padj")
    S4Vectors::metadata(v) <- md
  }

  v
}

test_that("as_vista adds required columns and preserves aligned dimensions", {
  mat <- matrix(rnorm(24), nrow = 6)
  rownames(mat) <- paste0("gene", seq_len(nrow(mat)))
  colnames(mat) <- paste0("sample", seq_len(ncol(mat)))

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(norm_counts = mat),
    colData = S4Vectors::DataFrame(
      cond = rep(c("A", "B"), each = 2),
      row.names = colnames(mat)
    ),
    rowData = S4Vectors::DataFrame(
      symbol = paste0("SYM", seq_len(nrow(mat))),
      row.names = rownames(mat)
    )
  )

  v <- as_vista(se, group_column = "cond")

  expect_s4_class(v, "VISTA")
  expect_identical(rownames(as.data.frame(sample_info(v))), colnames(norm_counts(v)))
  expect_identical(rownames(as.data.frame(row_data(v))), rownames(norm_counts(v)))
  expect_true("sample_names" %in% colnames(as.data.frame(sample_info(v))))
  expect_true("gene_id" %in% colnames(as.data.frame(row_data(v))))
})

test_that(".vista validates core constructor inputs", {
  mat <- matrix(rnorm(20), nrow = 5)
  rownames(mat) <- paste0("g", 1:5)
  colnames(mat) <- paste0("s", 1:4)
  si <- data.frame(cond = c("A", "A", "B", "B"), row.names = colnames(mat), stringsAsFactors = FALSE)
  rd <- data.frame(gene_id = rownames(mat), row.names = rownames(mat), stringsAsFactors = FALSE)

  expect_error(
    VISTA:::.vista(unname(mat), si, rd, group_column = "cond"),
    "rownames"
  )
  expect_error(
    VISTA:::.vista(mat, si[rev(rownames(si)), , drop = FALSE], rd, group_column = "cond"),
    "sample_info rownames must exactly match"
  )
  expect_error(
    VISTA:::.vista(mat, si, rd, group_column = "missing"),
    "not found"
  )

  bad_cmp <- data.frame(log2fc = rnorm(5), row.names = paste0("x", 1:5), stringsAsFactors = FALSE)
  expect_error(
    VISTA:::.vista(mat, si, rd, group_column = "cond", comparisons = list(c1 = bad_cmp)),
    "must have rownames identical"
  )
})

test_that("print methods return invisibly and produce SE-like output", {
  vista <- make_mock_vista(with_comparisons = TRUE)

  printed <- capture.output(out <- withVisible(print(vista)))
  expect_false(out$visible)
  expect_identical(out$value, vista)
  expect_true(length(printed) > 0)

  printed_legacy <- capture.output(out_legacy <- withVisible(print.vista(vista)))
  expect_false(out_legacy$visible)
  expect_identical(out_legacy$value, vista)
  expect_true(length(printed_legacy) > 0)
})

test_that("run_vista_report validates config inputs before rendering", {
  skip_if_not_installed("quarto")

  expect_error(
    run_vista_report(config = 1L, output_file = tempfile(fileext = ".html")),
    "must be a YAML path or a named list"
  )

  expect_error(
    run_vista_report(config = list(sidebar_position = "center"), output_file = tempfile(fileext = ".html")),
    "sidebar_position"
  )

  expect_error(
    run_vista_report(config = list(), output_file = ""),
    "output_file"
  )

  expect_error(
    run_vista_report(config = list(sidebar_position = "left"), output_file = tempfile(fileext = ".html")),
    "Provide either"
  )
})

test_that("run_vista_report validates required fields and comparison names", {
  skip_if_not_installed("quarto")

  counts <- data.frame(gene_id = c("g1", "g2"), s1 = c(10, 20), s2 = c(30, 40), check.names = FALSE)
  si <- data.frame(sample_names = c("s1", "s2"), cond = c("A", "B"), stringsAsFactors = FALSE)

  expect_error(
    run_vista_report(
      config = list(counts = counts, sample_info = si),
      output_file = tempfile(fileext = ".html")
    ),
    "Missing required configuration fields"
  )

  vista <- make_mock_vista(with_comparisons = TRUE)
  vista_rds <- tempfile(fileext = ".rds")
  saveRDS(vista, vista_rds)

  expect_error(
    run_vista_report(
      config = list(vista_rds = vista_rds, primary_comparison = "not_a_comp"),
      output_file = tempfile(fileext = ".html")
    ),
    "primary_comparison"
  )
})

test_that("run_vista_report errors when VISTA has no comparisons", {
  skip_if_not_installed("quarto")

  vista_empty <- make_mock_vista(with_comparisons = FALSE)
  vista_rds <- tempfile(fileext = ".rds")
  saveRDS(vista_empty, vista_rds)

  expect_error(
    run_vista_report(
      config = list(vista_rds = vista_rds),
      output_file = tempfile(fileext = ".html")
    ),
    "No differential comparisons found"
  )
})
