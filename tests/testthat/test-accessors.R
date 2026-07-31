test_that("accessors expose metadata correctly", {
  vista <- make_small_vista()

  summarized <- norm_counts(vista, summarise = TRUE)
  expect_true(is.matrix(summarized))
  expect_equal(ncol(summarized), length(unique(SummarizedExperiment::colData(vista)$cond_long)))

  summary_list <- deg_summary(vista)
  expect_true(all(vapply(summary_list, function(df) "regulation" %in% names(df), logical(1))))

  cut <- cutoffs(vista)
  expect_true(all(c("log2fc", "pval", "p_value_type", "method", "min_counts", "min_replicates") %in% names(cut)))
})

test_that("get_genes_by_regulation returns genes", {
  vista <- make_small_vista()
  comp_name <- names(comparisons(vista))[1]

  res <- get_genes_by_regulation(vista, comp_name, regulation = "Up")
  expect_true(is.list(res))
  expect_named(res)
})

test_that("get_genes_by_regulation supports top_n ranking and table output", {
  vista <- make_small_vista()
  rowData(vista)$SYMBOL <- paste0("GENE", seq_len(nrow(vista)))
  comp_name <- names(comparisons(vista))[1]

  up_full <- get_genes_by_regulation(vista, comp_name, regulation = "Up")[[1]]
  up_top <- get_genes_by_regulation(vista, comp_name, regulation = "Up", top_n = 3)[[1]]
  expect_lte(length(up_top), 3)
  expect_true(all(up_top %in% up_full))

  comp_tbl <- comparisons(vista)[[comp_name]]
  up_tbl <- comp_tbl[comp_tbl$gene_id %in% up_full, , drop = FALSE]
  expected_top <- up_tbl$gene_id[order(abs(up_tbl$log2fc), decreasing = TRUE, na.last = NA)]
  expect_identical(up_top, unique(utils::head(expected_top, n = length(up_top))))

  res_tbl <- get_genes_by_regulation(
    vista,
    comp_name,
    regulation = "Up",
    top_n = 3,
    display_id = "SYMBOL",
    return_type = "table"
  )[[1]]
  expect_s3_class(res_tbl, "tbl_df")
  expect_true(all(c("gene_id", "SYMBOL") %in% colnames(res_tbl)))
  expect_lte(nrow(res_tbl), 3)
})

test_that("get_cell_fractions reports absence", {
  vista <- make_small_vista()
  expect_error(get_cell_fractions(vista), "No cell fraction estimates")
})

test_that("get_foldchange_matrix constructs matrix", {
  vista <- make_small_vista()
  mat <- get_foldchange_matrix(vista)
  expect_true(is.matrix(mat))
  expect_equal(colnames(mat), names(comparisons(vista)))
  expect_gt(nrow(mat), 0)
})

test_that("norm_counts returns normalized count matrix", {
  vista <- make_small_vista()

  nc <- norm_counts(vista)
  expect_true(is.matrix(nc))
  expect_gt(nrow(nc), 0)
  expect_gt(ncol(nc), 0)
})

test_that("norm_counts with summarise=TRUE computes group means", {
  vista <- make_small_vista()

  nc_summary <- norm_counts(vista, summarise = TRUE)
  expect_true(is.matrix(nc_summary))
  expect_equal(ncol(nc_summary), length(unique(SummarizedExperiment::colData(vista)$cond_long)))
})

test_that("comparisons returns proper list structure", {
  vista <- make_small_vista()

  comps <- comparisons(vista)
  expect_type(comps, "list")
  expect_true(all(vapply(comps, is.data.frame, logical(1))))
  expect_true(all(vapply(comps, nrow, integer(1)) > 0))
  expect_named(comps)
})

test_that("deg_summary contains expected structure", {
  vista <- make_small_vista()

  summary <- deg_summary(vista)
  expect_type(summary, "list")
  expect_length(summary, length(comparisons(vista)))
  expect_true(all(vapply(summary, function(df) "regulation" %in% names(df), logical(1))))
})

test_that("cutoffs returns all required elements", {
  vista <- make_small_vista()

  cut <- cutoffs(vista)
  expect_type(cut, "list")
  expect_true(all(c("log2fc", "pval", "p_value_type", "method", "min_counts", "min_replicates") %in% names(cut)))
  expect_true(is.numeric(cut$log2fc))
  expect_true(is.numeric(cut$pval))
})

test_that("sample_info returns DataFrame", {
  vista <- make_small_vista()

  si <- sample_info(vista)
  expect_s4_class(si, "DataFrame")
  expect_gt(nrow(si), 0)
})

test_that("row_data returns DataFrame", {
  vista <- make_small_vista()

  rd <- row_data(vista)
  expect_s4_class(rd, "DataFrame")
  expect_equal(nrow(rd), nrow(vista))
})

test_that("group_colors returns named color vector", {
  vista <- make_small_vista()

  colors <- group_colors(vista)
  expect_type(colors, "character")
  expect_named(colors)
  expect_length(colors, length(unique(SummarizedExperiment::colData(vista)$cond_long)))
})

test_that("group_palette returns palette name", {
  vista <- make_small_vista()

  palette <- group_palette(vista)
  expect_type(palette, "character")
  expect_length(palette, 1)
})

test_that("accessors work after object modification", {
  vista <- make_small_vista()

  # Add new column to rowData
  rowData(vista)$test_col <- "test_value"

  # Accessors should still work
  rd <- row_data(vista)
  expect_true("test_col" %in% colnames(rd))

  comps <- comparisons(vista)
  expect_type(comps, "list")
})

test_that("norm_counts handles different assay names", {
  vista <- make_small_vista()

  # Should have norm_counts assay
  expect_true("norm_counts" %in% SummarizedExperiment::assayNames(vista))

  nc <- norm_counts(vista)
  expect_true(is.matrix(nc))
})

test_that("comparisons and deg_summary support source argument", {
  vista <- make_small_vista()
  comps_active <- comparisons(vista, source = "active")
  summary_active <- deg_summary(vista, source = "active")

  expect_type(comps_active, "list")
  expect_type(summary_active, "list")
  expect_named(comps_active)
})

test_that("set_de_source switches active DE source when available", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  vista_both <- create_vista(
    counts = count_data[1:200, ],
    sample_info = sample_metadata,
    column_geneid = "gene_id",
    group_column = "cond_long",
    group_numerator = "treatment1",
    group_denominator = "control",
    method = "both",
    result_source = "consensus",
    min_counts = 5,
    min_replicates = 1
  )

  vista_switched <- set_de_source(vista_both, "edger")
  expect_identical(S4Vectors::metadata(vista_switched)$de_active_source, "edger")
  expect_named(comparisons(vista_switched))

  vista_same <- set_de_source(vista_switched, "active")
  expect_identical(S4Vectors::metadata(vista_same)$de_active_source, "edger")
})

test_that("set_vista_group_colors updates group colors", {
  vista <- make_small_vista()
  custom_cols <- c(
    control = "#264653",
    treatment1 = "#E76F51"
  )

  vista_custom <- set_vista_group_colors(vista, custom_cols)
  expect_identical(group_colors(vista_custom), custom_cols[names(group_colors(vista_custom))])
  expect_identical(group_palette(vista_custom), "manual")
})

test_that("set_vista_group_colors validates complete mapping", {
  vista <- make_small_vista()
  expect_error(
    set_vista_group_colors(vista, c(control = "#264653")),
    "Missing group color"
  )
})

test_that("set_vista_comparison_colors updates comparison colors", {
  vista <- make_small_vista()
  comp_name <- names(comparisons(vista))
  custom_cols <- stats::setNames("#6C5CE7", comp_name)

  vista_custom <- set_vista_comparison_colors(vista, custom_cols)
  expect_identical(S4Vectors::metadata(vista_custom)$comparison$colors[comp_name], custom_cols[comp_name])
  expect_identical(S4Vectors::metadata(vista_custom)$comparison$palette, "manual")
})

test_that("set_vista_comparison_colors validates complete mapping", {
  vista <- make_small_vista()
  expect_error(
    set_vista_comparison_colors(vista, c(wrong_name = "#6C5CE7")),
    "Missing comparison color"
  )
})

test_that("group summarisation drops empty factor levels instead of emitting NaN", {
  v <- make_small_vista()
  group_col <- S4Vectors::metadata(v)$group$column

  # Force the grouping column to a factor that retains a level with no samples,
  # which is what subsetting an object produces.
  cd <- SummarizedExperiment::colData(v)
  cd[[group_col]] <- factor(
    as.character(cd[[group_col]]),
    levels = c(unique(as.character(cd[[group_col]])), "ghost_group")
  )
  SummarizedExperiment::colData(v) <- cd

  nc <- norm_counts(v, summarise = TRUE)
  expect_false(anyNA(nc))
  expect_false("ghost_group" %in% colnames(nc))
  expect_setequal(colnames(nc), unique(as.character(cd[[group_col]])))

  em <- get_expression_matrix(v, summarise = TRUE)
  expect_false(anyNA(em))
  expect_false("ghost_group" %in% colnames(em))

  # The surviving groups still carry the right values.
  ctrl <- rownames(cd)[as.character(cd[[group_col]]) == "control"]
  expect_equal(
    unname(nc[, "control"]),
    unname(rowMeans(SummarizedExperiment::assay(v, "norm_counts")[, ctrl, drop = FALSE]))
  )
})

test_that("accessor generics dispatch only on the object (2e)", {
  for (g in c("comparisons", "deg_summary", "cutoffs", "norm_counts",
              "sample_info", "row_data", "group_colors", "group_palette")) {
    expect_identical(
      methods::getGeneric(g)@signature, "object",
      info = g
    )
  }
  expect_identical(methods::getGeneric("enrichMsigDB")@signature, "x")
})

test_that("non-class arguments are no longer forced at dispatch", {
  v <- make_small_vista()

  # Every formal used to join the dispatch signature, so `summarise` was
  # evaluated before the method body ran.
  expect_silent(nc <- norm_counts(v))
  expect_true(is.matrix(nc))

  # A lazy default that would error if forced at dispatch is fine now.
  f <- function(obj, s = stop("forced at dispatch")) norm_counts(obj)
  expect_true(is.matrix(f(v)))
})
