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

# --- B8/B9: the report must describe the analysis that actually ran ---

test_that("report parameters come from the object, not the config defaults (B8)", {
  skip_if_not_installed("quarto")
  skip_if_not_installed("yaml")

  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())
  si <- sample_metadata[sample_metadata$cond_long %in% c("control", "treatment1"), ]
  cnt <- count_data[seq_len(300), c("gene_id", si$sample_names)]

  # Built with settings that differ from every run_vista_report() default.
  v <- suppressMessages(create_vista(
    counts = cnt, sample_info = si, column_geneid = "gene_id",
    group_column = "cond_long", group_numerator = "treatment1",
    group_denominator = "control",
    method = "edger", log2fc_cutoff = 2, pval_cutoff = 0.01,
    p_value_type = "pvalue", min_counts = 5, min_replicates = 1
  ))

  cuts <- cutoffs(v)
  expect_identical(cuts$method, "edger")
  expect_equal(cuts$log2fc, 2)
  expect_equal(cuts$pval, 0.01)
  expect_identical(cuts$p_value_type, "pvalue")

  rds <- tempfile(fileext = ".rds")
  on.exit(unlink(rds), add = TRUE)
  saveRDS(v, rds)

  out_dir <- tempfile()
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  # Rendering needs the quarto CLI, which may be absent; the reconciliation we
  # care about happens before that and writes summary_*.csv either way.
  try(
    suppressWarnings(suppressMessages(run_vista_report(
      list(
        vista_rds = rds,
        assets_dir = "assets",
        output_file = file.path(out_dir, "r.html"),
        include_msigdb = FALSE, include_go = FALSE, include_kegg = FALSE,
        include_pathway_heatmap = FALSE
      )
    ))),
    silent = TRUE
  )

  summary_csv <- list.files(
    out_dir, pattern = "^summary_.*\\.csv$", recursive = TRUE, full.names = TRUE
  )
  skip_if(length(summary_csv) == 0, "report did not reach the summary table")

  tbl <- utils::read.csv(summary_csv[[1]], stringsAsFactors = FALSE)
  val <- function(m) tbl$value[match(m, tbl$metric)]

  # Previously these reported the config defaults deseq2 / 1 / 0.05 / padj.
  expect_identical(val("DE method"), "edger")
  expect_identical(val("LFC cutoff"), "2")
  expect_identical(val("P-value cutoff"), "0.01")
  expect_identical(val("P-value type"), "pvalue")

  expect_true("source" %in% colnames(tbl))
  expect_identical(tbl$source[match("DE method", tbl$metric)], "object")
})

test_that("enrichment identifier type is detected, not taken from display_id (B9)", {
  v <- make_small_vista()
  de_tbl <- comparisons(v)[[1]]
  expect_true(any(grepl("^ENS", de_tbl$gene_id)))

  # Reproduce the resolution rule now used in run_vista_report().
  resolve <- function(cfg) {
    cfg$from_type %||%
      (if (any(grepl("^ENS", de_tbl$gene_id))) "ENSEMBL" else "SYMBOL")
  }

  # The problem pairing: Ensembl IDs with SYMBOL labels and no from_type.
  expect_identical(resolve(list(display_id = "SYMBOL")), "ENSEMBL")
  # An explicit from_type still wins.
  expect_identical(resolve(list(from_type = "SYMBOL", display_id = "SYMBOL")), "SYMBOL")
  # display_id is irrelevant to the decision now.
  expect_identical(
    resolve(list(display_id = "SYMBOL")),
    resolve(list(display_id = "GENENAME"))
  )

  # Guard the implementation itself: display_id must not be in the chain.
  body_txt <- paste(deparse(body(run_vista_report)), collapse = " ")
  chain <- regmatches(
    body_txt, regexpr("enrich_from <- [^\n]*?(?=\\s{2,}if \\(is\\.null)", body_txt, perl = TRUE)
  )
  expect_length(chain, 1L)
  expect_false(grepl("display_id", chain, fixed = TRUE))
  expect_true(grepl("from_type", chain, fixed = TRUE))
})
