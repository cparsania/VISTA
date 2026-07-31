test_that("get_volcano_plot inherits the object's p_value_type and cutoffs", {
  v <- make_small_vista()
  comp <- names(comparisons(v))[1]
  cuts <- cutoffs(v)

  expect_identical(cuts$p_value_type, "padj")

  # EnhancedVolcano still uses the retired `size` line aesthetic upstream.
  p <- suppressWarnings(get_volcano_plot(v, sample_comparison = comp))
  expect_s3_class(p, "ggplot")

  # The y-axis must be the adjusted p-value, matching how the object was built.
  ylab <- p$labels$y
  expect_true(grepl("adj", paste(deparse(ylab), collapse = "")))
})

test_that("volcano colouring agrees with the stored regulation calls", {
  v <- make_small_vista()
  comp <- names(comparisons(v))[1]
  tbl <- comparisons(v)[[comp]]
  cuts <- cutoffs(v)

  # Replay the colour assignment the plot performs.
  keyvals <- attr(
    .EnhancedVolcano2(
      tbl,
      lab = tbl$gene_id,
      x = "log2fc",
      y = cuts$p_value_type,
      pCutoff = cuts$pval,
      FCcutoff = cuts$log2fc,
      col_by_regul = TRUE,
      return_keyvals = TRUE
    ),
    "keyvals"
  )

  plotted <- table(factor(names(keyvals), levels = c("Up", "Down", "Other")))
  stored <- table(factor(tbl$regulation, levels = c("Up", "Down", "Other")))

  expect_equal(as.vector(plotted), as.vector(stored))
})

test_that("get_volcano_plot arguments still override the stored cutoffs", {
  v <- make_small_vista()
  comp <- names(comparisons(v))[1]

  p <- suppressWarnings(
    get_volcano_plot(v, sample_comparison = comp, p_value_type = "pvalue")
  )
  expect_s3_class(p, "ggplot")
  expect_false(grepl("adj", paste(deparse(p$labels$y), collapse = "")))

  # Explicit thresholds must win over cutoffs(v).
  p2 <- suppressWarnings(get_volcano_plot(
    v,
    sample_comparison = comp,
    log2fc_cutoff = 3,
    pval_cutoff = 1e-6
  ))
  expect_s3_class(p2, "ggplot")
})

test_that("volcano falls back with a warning when the requested column is absent", {
  v <- make_small_vista()
  comp <- names(comparisons(v))[1]

  # comparisons() resolves through de_results_by_method for the active source,
  # so the active method's table is the one that must lose its padj column.
  md <- S4Vectors::metadata(v)
  active <- md$de_active_source
  stripped <- as.data.frame(md$de_results_by_method[[active]][[comp]], stringsAsFactors = FALSE)
  stripped$padj <- NULL
  md$de_results_by_method[[active]][[comp]] <- stripped
  md$de_results[[comp]] <- stripped
  S4Vectors::metadata(v) <- md

  expect_false("padj" %in% colnames(comparisons(v)[[comp]]))

  expect_warning(
    get_volcano_plot(v, sample_comparison = comp, p_value_type = "padj"),
    "not found"
  )
})

test_that("consensus p-values stay continuous and conservative", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  si <- sample_metadata[sample_metadata$cond_long %in% c("control", "treatment1"), ]
  cnt <- count_data[seq_len(400), c("gene_id", si$sample_names)]

  v <- suppressMessages(create_vista(
    counts = cnt,
    sample_info = si,
    column_geneid = "gene_id",
    group_column = "cond_long",
    group_numerator = "treatment1",
    group_denominator = "control",
    method = "both",
    result_source = "consensus",
    min_counts = 5,
    min_replicates = 1
  ))

  tbl <- comparisons(v)[[1]]

  # Previously ~95% of genes were forced to exactly 1.
  expect_lt(mean(tbl$padj == 1, na.rm = TRUE), 0.5)
  expect_lt(mean(tbl$pvalue == 1, na.rm = TRUE), 0.5)

  # Genes not uniquely called by one backend carry the less-significant value.
  shared <- tbl$support %in% c("both", "discordant", "none")
  expect_equal(
    tbl$padj[shared],
    pmax(tbl$padj_deseq2, tbl$padj_edger)[shared]
  )
  expect_equal(
    tbl$pvalue[shared],
    pmax(tbl$pvalue_deseq2, tbl$pvalue_edger)[shared]
  )

  # Single-backend calls keep the contributing backend's value.
  d_only <- tbl$support == "deseq2_only"
  if (any(d_only)) {
    expect_equal(tbl$padj[d_only], tbl$padj_deseq2[d_only])
  }
  e_only <- tbl$support == "edger_only"
  if (any(e_only)) {
    expect_equal(tbl$padj[e_only], tbl$padj_edger[e_only])
  }

  # Per-backend columns must survive intact.
  expect_true(all(
    c("log2fc_deseq2", "log2fc_edger", "pvalue_deseq2", "pvalue_edger",
      "padj_deseq2", "padj_edger", "support") %in% colnames(tbl)
  ))
})

test_that("min_counts/min_replicates filter identically across backends", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  si <- sample_metadata[sample_metadata$cond_long %in% c("control", "treatment1"), ]
  cnt <- count_data[seq_len(400), c("gene_id", si$sample_names)]

  build <- function(method) {
    suppressMessages(create_vista(
      counts = cnt,
      sample_info = si,
      column_geneid = "gene_id",
      group_column = "cond_long",
      group_numerator = "treatment1",
      group_denominator = "control",
      method = method,
      min_counts = 5,
      min_replicates = 3
    ))
  }

  v_deseq <- build("deseq2")
  v_edger <- build("edger")
  v_limma <- build("limma")

  expect_identical(rownames(v_edger), rownames(v_deseq))
  expect_identical(rownames(v_limma), rownames(v_deseq))

  # And the retained set is exactly the documented predicate.
  raw <- as.data.frame(cnt, stringsAsFactors = FALSE)
  rownames(raw) <- raw$gene_id
  raw$gene_id <- NULL
  raw <- as.matrix(raw)
  raw[is.na(raw)] <- 0
  prefiltered <- raw[rowSums(raw) >= 5, , drop = FALSE]
  expected <- rownames(prefiltered)[rowSums(prefiltered >= 5) >= 3]

  expect_setequal(rownames(v_deseq), expected)
})

test_that("stricter min_replicates removes genes for every backend", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  si <- sample_metadata[sample_metadata$cond_long %in% c("control", "treatment1"), ]
  cnt <- count_data[seq_len(400), c("gene_id", si$sample_names)]

  build <- function(method, min_replicates) {
    suppressMessages(create_vista(
      counts = cnt,
      sample_info = si,
      column_geneid = "gene_id",
      group_column = "cond_long",
      group_numerator = "treatment1",
      group_denominator = "control",
      method = method,
      min_counts = 5,
      min_replicates = min_replicates
    ))
  }

  for (m in c("deseq2", "edger", "limma")) {
    lenient <- build(m, 1)
    strict <- build(m, 4)
    expect_lt(nrow(strict), nrow(lenient))
    expect_true(all(rownames(strict) %in% rownames(lenient)))
  }
})
