test_that("create_vista returns a valid VISTA object", {
  vista <- make_small_vista()

  expect_s4_class(vista, "VISTA")
  expect_true("norm_counts" %in% assayNames(vista))

  md <- S4Vectors::metadata(vista)
  expect_s4_class(md$de_results, "SimpleList")
  expect_s4_class(md$de_summary, "SimpleList")
  expect_type(md$de_cutoffs, "list")

  comps <- comparisons(vista)
  expect_length(comps, 1L)
  expect_true(all(vapply(comps, nrow, integer(1)) > 0))

  nc <- norm_counts(vista)
  expect_equal(ncol(nc), 6)
  expect_gt(nrow(nc), 0)

  colors <- group_colors(vista)
  expect_named(colors)
  expect_length(colors, length(unique(SummarizedExperiment::colData(vista)$cond_long)))
})

test_that("create_vista validates group vector lengths", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  expect_error(
    create_vista(
      counts = count_data[1:100, ],
      sample_info = sample_metadata,
      column_geneid = "gene_id",
      group_column = "cond_long",
      group_numerator = c("treatment1", "control"),
      group_denominator = "control",  # Length mismatch
      min_counts = 5
    ),
    "must have equal length"
  )
})

test_that("create_vista validates groups exist in data", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  expect_error(
    create_vista(
      counts = count_data[1:100, ],
      sample_info = sample_metadata,
      column_geneid = "gene_id",
      group_column = "cond_long",
      group_numerator = "NONEXISTENT_GROUP",
      group_denominator = "control",
      min_counts = 5
    ),
    "not in data"
  )
})

test_that("create_vista validates empty group vectors", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  expect_error(
    create_vista(
      counts = count_data[1:100, ],
      sample_info = sample_metadata,
      column_geneid = "gene_id",
      group_column = "cond_long",
      group_numerator = character(0),
      group_denominator = "control",
      min_counts = 5
    ),
    "non-empty"
  )
})

test_that("create_vista validates group_column exists", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  expect_error(
    create_vista(
      counts = count_data[1:100, ],
      sample_info = sample_metadata,
      column_geneid = "gene_id",
      group_column = "nonexistent_column",
      group_numerator = "treatment1",
      group_denominator = "control",
      min_counts = 5
    ),
    "not found"
  )
})

test_that("create_vista warns about low replication", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  # Create single-replicate subset (one sample per group)
  single_rep <- sample_metadata[c(1, 4), ]
  counts_sub <- count_data[1:100, c("gene_id", single_rep$sample_names)]

  # Expect the warning about low replication (DESeq2 may error after the warning)
  expect_warning(
    tryCatch(
      create_vista(
        counts = counts_sub,
        sample_info = single_rep,
        column_geneid = "gene_id",
        group_column = "cond_long",
        group_numerator = single_rep$cond_long[2],
        group_denominator = single_rep$cond_long[1],
        min_counts = 1,
        min_replicates = 1
      ),
      error = function(e) NULL
    ),
    "fewer than 2 replicates"
  )
})

test_that("create_vista works with edgeR method", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  target_groups <- c("treatment1", "control")
  sample_subset <- sample_metadata[sample_metadata$cond_long %in% target_groups, ]
  counts_subset <- count_data[1:150, c("gene_id", sample_subset$sample_names)]

  vista_edger <- create_vista(
    counts = counts_subset,
    sample_info = sample_subset,
    column_geneid = "gene_id",
    group_column = "cond_long",
    group_numerator = "treatment1",
    group_denominator = "control",
    method = "edger",
    min_counts = 5,
    min_replicates = 1
  )

  expect_s4_class(vista_edger, "VISTA")
  expect_true("norm_counts" %in% assayNames(vista_edger))
})

test_that("create_vista works with limma method", {
  skip_if_not_installed("limma")
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  target_groups <- c("treatment1", "control")
  sample_subset <- sample_metadata[sample_metadata$cond_long %in% target_groups, ]
  counts_subset <- count_data[1:150, c("gene_id", sample_subset$sample_names)]

  vista_limma <- create_vista(
    counts = counts_subset,
    sample_info = sample_subset,
    column_geneid = "gene_id",
    group_column = "cond_long",
    group_numerator = "treatment1",
    group_denominator = "control",
    method = "limma",
    min_counts = 5,
    min_replicates = 1
  )

  expect_s4_class(vista_limma, "VISTA")
  expect_true("norm_counts" %in% assayNames(vista_limma))
  expect_named(comparisons(vista_limma))
})

test_that("create_vista handles multiple comparisons", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())
  cell_levels <- unique(sample_metadata$cell)

  vista_multi <- create_vista(
    counts = count_data[1:150, ],
    sample_info = sample_metadata,
    column_geneid = "gene_id",
    group_column = "cell",
    group_numerator = cell_levels[2:3],
    group_denominator = rep(cell_levels[1], 2),
    min_counts = 5,
    min_replicates = 1
  )

  expect_s4_class(vista_multi, "VISTA")
  comps <- comparisons(vista_multi)
  expect_length(comps, 2)  # Should have 2 comparisons
})

test_that("create_vista supports covariate-adjusted design", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  vista_cov <- create_vista(
    counts = count_data[1:200, ],
    sample_info = sample_metadata,
    column_geneid = "gene_id",
    group_column = "cond_long",
    group_numerator = "treatment1",
    group_denominator = "control",
    method = "deseq2",
    covariates = "cell",
    min_counts = 5,
    min_replicates = 1
  )

  expect_s4_class(vista_cov, "VISTA")
  expect_true("cell" %in% cutoffs(vista_cov)$covariates)
})

test_that("create_vista supports explicit design formula", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  vista_formula <- create_vista(
    counts = count_data[1:200, ],
    sample_info = sample_metadata,
    column_geneid = "gene_id",
    group_column = "cond_long",
    group_numerator = "treatment1",
    group_denominator = "control",
    method = "deseq2",
    design_formula = "~ cell + cond_long",
    min_counts = 5,
    min_replicates = 1
  )

  expect_s4_class(vista_formula, "VISTA")
  expect_match(cutoffs(vista_formula)$design_formula, "cond_long")
})

test_that("create_vista supports combined method with consensus active source", {
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

  expect_s4_class(vista_both, "VISTA")
  md <- S4Vectors::metadata(vista_both)
  expect_true(all(c("deseq2", "edger", "consensus") %in% names(md$de_results_by_method)))
  expect_identical(md$de_active_source, "consensus")
  expect_true(all(vapply(comparisons(vista_both), is.data.frame, logical(1))))
})

test_that("create_vista supports covariates with edgeR backend", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  vista_edger_cov <- create_vista(
    counts = count_data[1:200, ],
    sample_info = sample_metadata,
    column_geneid = "gene_id",
    group_column = "cond_long",
    group_numerator = "treatment1",
    group_denominator = "control",
    method = "edger",
    covariates = "cell",
    min_counts = 5,
    min_replicates = 1
  )

  expect_s4_class(vista_edger_cov, "VISTA")
  expect_true("cell" %in% cutoffs(vista_edger_cov)$covariates)
})

test_that("create_vista supports combined method with covariates", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  vista_both_cov <- create_vista(
    counts = count_data[1:200, ],
    sample_info = sample_metadata,
    column_geneid = "gene_id",
    group_column = "cond_long",
    group_numerator = "treatment1",
    group_denominator = "control",
    method = "both",
    covariates = "cell",
    result_source = "consensus",
    min_counts = 5,
    min_replicates = 1
  )

  expect_s4_class(vista_both_cov, "VISTA")
  expect_named(comparisons(vista_both_cov, source = "deseq2"))
  expect_named(comparisons(vista_both_cov, source = "edger"))
  expect_named(comparisons(vista_both_cov, source = "consensus"))
})

test_that("create_vista accepts sample_info rownames when sample_names column is absent", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  sample_no_names <- as.data.frame(sample_metadata)
  rownames(sample_no_names) <- sample_no_names$sample_names
  sample_no_names[["sample_names"]] <- NULL

  vista_rn <- create_vista(
    counts = count_data[1:150, ],
    sample_info = sample_no_names,
    column_geneid = "gene_id",
    group_column = "cond_long",
    group_numerator = "treatment1",
    group_denominator = "control",
    min_counts = 5,
    min_replicates = 1
  )

  expect_s4_class(vista_rn, "VISTA")
  expect_true("sample_names" %in% colnames(as.data.frame(sample_info(vista_rn))))
})
