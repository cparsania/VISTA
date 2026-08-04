test_that("enrichment helpers error clearly without orgdb", {
  vista <- make_small_vista()
  comp <- names(comparisons(vista))[1]
  expect_error(get_go_enrichment(vista, comp, species = "Unknownus"), "orgdb")
  expect_error(get_kegg_enrichment(vista, comp, species = "Unknownus"), "orgdb")
})

test_that("get_gsea errors if no comparison", {
  vista <- make_small_vista()
  expect_error(get_gsea(vista, sample_comparison = "missing"), "Comparison")
})

test_that("get_msigdb_enrichment handles invalid comparisons", {
  vista <- make_small_vista()

  expect_error(
    get_msigdb_enrichment(vista, sample_comparison = "NONEXISTENT", regulation = "Up"),
    "Comparison|comparison|not found"
  )
})

test_that("get_msigdb_enrichment validates regulation parameter", {
  skip_if_not_installed("msigdbr")
  vista <- make_small_vista()
  comps <- names(comparisons(vista))

  # With small test data, there may be no significant DEGs
  # so we verify it either succeeds or gives a clear error about no genes
  tryCatch(
    {
      result <- get_msigdb_enrichment(
        vista,
        sample_comparison = comps[1],
        regulation = "Up",
        from_type = "ENSEMBL"
      )
      expect_type(result, "list")
    },
    error = function(e) {
      expect_match(conditionMessage(e), "No genes found|not found")
    }
  )
})

test_that("enrichMsigDB validates gene_list input", {
  vista <- make_small_vista()

  expect_error(
    enrichMsigDB(vista, gene_list = character(0)),
    "gene_list.*required|empty"
  )

  expect_error(
    enrichMsigDB(vista, gene_list = NULL),
    "gene_list.*required"
  )
})

test_that("enrichMsigDB returns proper structure", {
  skip_if_not_installed("msigdbr")
  vista <- make_small_vista()
  genes <- rownames(vista)[1:20]

  result <- enrichMsigDB(
    vista,
    gene_list = genes,
    from_type = "ENSEMBL",
    msigdb_category = "H"
  )

  expect_type(result, "list")
  expect_true("enrich" %in% names(result))
})

test_that("get_genes_by_regulation extracts correct gene sets", {
  vista <- make_small_vista()
  comps <- names(comparisons(vista))

  up_genes <- get_genes_by_regulation(
    vista,
    sample_comparisons = comps[1],
    regulation = "Up"
  )

  expect_type(up_genes, "list")
  expect_named(up_genes)

  down_genes <- get_genes_by_regulation(
    vista,
    sample_comparisons = comps[1],
    regulation = "Down"
  )

  expect_type(down_genes, "list")
  expect_named(down_genes)
})

test_that("get_enrichment_plot handles enrichResult objects", {
  skip_if_not_installed("msigdbr")
  skip_if_not_installed("enrichplot")

  vista <- make_small_vista()
  genes <- rownames(vista)[1:30]

  result <- enrichMsigDB(
    vista,
    gene_list = genes,
    from_type = "ENSEMBL",
    msigdb_category = "H"
  )

  # Should not error even if no significant results
  expect_no_error({
    p <- get_enrichment_plot(result$enrich)
  })
})

test_that("enrichment functions handle from_type parameter", {
  skip_if_not_installed("msigdbr")
  vista <- make_small_vista()

  # Add mock SYMBOL data
  rowData(vista)$SYMBOL <- paste0("GENE", seq_len(nrow(vista)))
  genes_symbol <- rowData(vista)$SYMBOL[1:20]

  result <- enrichMsigDB(
    vista,
    gene_list = genes_symbol,
    from_type = "SYMBOL",
    msigdb_category = "H"
  )

  expect_type(result, "list")
})

test_that("get_go_enrichment validates ontology parameter", {
  vista <- make_small_vista()
  comps <- names(comparisons(vista))

  # Valid ontologies should pass ontology argument validation.
  # Depending on environment this may continue or fail later (e.g., missing orgdb),
  # but it should not fail because ont is invalid.
  bp_err <- tryCatch(
    {
      get_go_enrichment(vista, sample_comparison = comps[1], ont = "BP")
      NULL
    },
    error = function(e) e
  )
  if (!is.null(bp_err)) {
    expect_false(grepl("one of.*BP.*MF.*CC|\\bont\\b", conditionMessage(bp_err), ignore.case = TRUE))
  }

  mf_err <- tryCatch(
    {
      get_go_enrichment(vista, sample_comparison = comps[1], ont = "MF")
      NULL
    },
    error = function(e) e
  )
  if (!is.null(mf_err)) {
    expect_false(grepl("one of.*BP.*MF.*CC|\\bont\\b", conditionMessage(mf_err), ignore.case = TRUE))
  }

  # Invalid ontology should fail fast.
  expect_error(
    get_go_enrichment(vista, sample_comparison = comps[1], ont = "INVALID"),
    "one of.*BP.*MF.*CC"
  )
})

test_that("get_pathway_genes extracts pathway members in multiple formats", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:20]

  term2gene <- data.frame(
    term = c(rep("Pathway_A", 10), rep("Pathway_B", 10)),
    gene = genes,
    stringsAsFactors = FALSE
  )
  enr <- clusterProfiler::enricher(gene = genes[1:12], TERM2GENE = term2gene)

  long_tbl <- get_pathway_genes(enr, top_n = 2, return_type = "long")
  expect_s3_class(long_tbl, "data.frame")
  expect_true(all(c("pathway_id", "pathway", "gene") %in% colnames(long_tbl)))

  gene_list <- get_pathway_genes(enr, pathways = "Pathway_A", return_type = "list")
  expect_type(gene_list, "list")
  expect_named(gene_list, "Pathway_A")
  expect_true(length(gene_list[[1]]) > 0)

  gene_vec <- get_pathway_genes(enr, top_n = 1, return_type = "vector")
  expect_type(gene_vec, "character")
  expect_true(length(gene_vec) > 0)
})

test_that("get_pathway_heatmap returns mapped genes and supports SYMBOL mapping", {
  vista <- make_small_vista()
  rowData(vista)$SYMBOL <- paste0("SYM", seq_len(nrow(vista)))

  symbols <- rowData(vista)$SYMBOL[1:20]
  term2gene <- data.frame(
    term = c(rep("Pathway_A", 10), rep("Pathway_B", 10)),
    gene = symbols,
    stringsAsFactors = FALSE
  )
  enr <- clusterProfiler::enricher(gene = symbols[1:12], TERM2GENE = term2gene)

  mapped_genes <- get_pathway_heatmap(
    x = vista,
    enrichment = enr,
    sample_group = unique(SummarizedExperiment::colData(vista)$cond_long),
    top_n = 1,
    gene_id_column = "SYMBOL",
    return_type = "data"
  )

  expect_type(mapped_genes, "character")
  expect_true(length(mapped_genes) > 0)
  expect_true(all(mapped_genes %in% rownames(vista)))
})

test_that("get_pathway_heatmap returns heatmap object when optional deps installed", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  vista <- make_small_vista()
  genes <- rownames(vista)[1:20]
  term2gene <- data.frame(
    term = c(rep("Pathway_A", 10), rep("Pathway_B", 10)),
    gene = genes,
    stringsAsFactors = FALSE
  )
  enr <- clusterProfiler::enricher(gene = genes[1:12], TERM2GENE = term2gene)

  out <- get_pathway_heatmap(
    x = vista,
    enrichment = enr,
    sample_group = unique(SummarizedExperiment::colData(vista)$cond_long),
    pathways = "Pathway_A",
    return_type = "both",
    value_transform = "zscore",
    summarise_replicates = TRUE,
    show_row_names = FALSE
  )

  expect_type(out, "list")
  expect_true(all(c("heatmap", "genes", "pathway_genes") %in% names(out)))
  expect_true(inherits(out$heatmap, "Heatmap") || inherits(out$heatmap, "HeatmapList"))
})

test_that("get_pathway_heatmap blocks genes and samples override via ...", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:20]
  term2gene <- data.frame(
    term = c(rep("Pathway_A", 10), rep("Pathway_B", 10)),
    gene = genes,
    stringsAsFactors = FALSE
  )
  enr <- clusterProfiler::enricher(gene = genes[1:12], TERM2GENE = term2gene)

  expect_error(
    get_pathway_heatmap(
      x = vista,
      enrichment = enr,
      sample_group = unique(SummarizedExperiment::colData(vista)$cond_long),
      pathways = "Pathway_A",
      genes = genes[1:5],
      return_type = "data"
    ),
    "managed by"
  )
})

test_that("get_pathway_genes validates pathway filters and top_n", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:20]

  term2gene <- data.frame(
    term = c(rep("Pathway_A", 10), rep("Pathway_B", 10)),
    gene = genes,
    stringsAsFactors = FALSE
  )
  enr <- clusterProfiler::enricher(gene = genes[1:12], TERM2GENE = term2gene)

  expect_error(
    get_pathway_genes(enr, pathways = "Not_A_Pathway"),
    "were not found"
  )

  expect_error(
    get_pathway_genes(enr, top_n = 0),
    "top_n"
  )

  id_tbl <- get_pathway_genes(enr, top_n = 1, pathway_column = "ID", return_type = "long")
  expect_true(all(id_tbl$pathway_id %in% as.data.frame(enr@result)$ID))
})

test_that("get_pathway_heatmap supports intersection and max_genes cap", {
  vista <- make_small_vista()
  genes <- rownames(vista)[1:30]

  term2gene <- data.frame(
    term = c(rep("Pathway_A", 15), rep("Pathway_B", 15)),
    gene = c(genes[1:15], genes[8:22]),
    stringsAsFactors = FALSE
  )
  enr <- clusterProfiler::enricher(gene = genes[1:20], TERM2GENE = term2gene)

  expected <- intersect(genes[1:15], genes[8:22])
  got <- get_pathway_heatmap(
    x = vista,
    enrichment = enr,
    sample_group = unique(SummarizedExperiment::colData(vista)$cond_long),
    top_n = 2,
    gene_mode = "intersection",
    max_genes = 5,
    return_type = "data"
  )

  expect_lte(length(got), 5)
  expect_true(all(got %in% expected))
})

test_that("get_gsea errors when no fold-change column is available", {
  vista <- make_small_vista()
  comp <- names(comparisons(vista))[1]

  md <- S4Vectors::metadata(vista)
  tbl <- as.data.frame(md$de_results[[comp]], stringsAsFactors = FALSE)
  tbl$log2fc <- NULL
  tbl$log2FoldChange <- NULL
  tbl$logFC <- NULL
  md$de_results_by_method <- NULL
  md$de_summary_by_method <- NULL
  md$de_active_source <- NULL
  md$de_results[[comp]] <- tbl
  S4Vectors::metadata(vista) <- md

  expect_error(
    get_gsea(vista, sample_comparison = comp, set_type = "msigdb"),
    "No log2FC"
  )
})

# --- B2: GSEA ranked vectors must keep every score attached to its own gene ---

test_that(".vista_map_ids_strict preserves length and position", {
  m <- c(TP53 = "7157", BRCA1 = "672", EGFR = "1956")
  ids <- c("TP53", "NOPE", "BRCA1", "ALSO_NOPE", "EGFR")

  out <- VISTA:::.vista_map_ids_strict(ids, "SYMBOL", "ENTREZID", id_map = m)

  expect_length(out, length(ids))
  expect_identical(out, c("7157", NA, "672", NA, "1956"))

  # Same key type is a pass-through.
  expect_identical(
    VISTA:::.vista_map_ids_strict(ids, "SYMBOL", "SYMBOL"), ids
  )
})

test_that(".vista_build_rank_vector drops unmapped genes without shifting scores", {
  m <- c(TP53 = "7157", BRCA1 = "672", EGFR = "1956")
  ids <- c("TP53", "NOPE", "BRCA1", "ALSO_NOPE", "EGFR")
  scores <- c(5, 4, 3, 2, 1)

  rv <- suppressMessages(VISTA:::.vista_build_rank_vector(
    scores = scores, ids = ids,
    from_type = "SYMBOL", to_type = "ENTREZID", id_map = m
  ))

  # Assigning a shorter name vector used to pad with NA and slide each score
  # onto the next gene: BRCA1 received 4 instead of 3, EGFR received 3 not 1.
  expect_length(rv, 3L)
  expect_false(anyNA(names(rv)))
  expect_equal(rv[["7157"]], 5)
  expect_equal(rv[["672"]], 3)
  expect_equal(rv[["1956"]], 1)

  # Sorted decreasing, as GSEA requires.
  expect_identical(rv, sort(rv, decreasing = TRUE))
})

test_that(".vista_build_rank_vector collapses many-to-one by largest absolute score", {
  m <- c(A = "1", B = "1", C = "2")
  rv <- suppressMessages(VISTA:::.vista_build_rank_vector(
    scores = c(0.5, -3, 2), ids = c("A", "B", "C"),
    from_type = "SYMBOL", to_type = "ENTREZID", id_map = m
  ))

  expect_length(rv, 2L)
  expect_equal(rv[["1"]], -3)   # |-3| > |0.5|
  expect_equal(rv[["2"]], 2)
})

test_that(".vista_build_rank_vector drops non-finite scores and errors when nothing maps", {
  m <- c(A = "1", B = "2")
  rv <- suppressMessages(VISTA:::.vista_build_rank_vector(
    scores = c(1, NA, 2), ids = c("A", "B", "A"),
    from_type = "SYMBOL", to_type = "ENTREZID", id_map = m
  ))
  expect_false(anyNA(rv))
  expect_length(rv, 1L)

  expect_error(
    suppressMessages(VISTA:::.vista_build_rank_vector(
      scores = c(1, 2), ids = c("X", "Y"),
      from_type = "SYMBOL", to_type = "ENTREZID", id_map = c(A = "1")
    )),
    "No gene identifiers could be mapped"
  )
})

test_that("get_gsea builds a correctly aligned rank vector end to end", {
  skip_if_not_installed("org.Hs.eg.db")
  v <- make_small_vista()
  comp <- names(comparisons(v))[[1]]
  tbl <- comparisons(v)[[comp]]

  rv <- suppressMessages(VISTA:::.vista_build_rank_vector(
    scores = tbl$log2fc, ids = tbl$gene_id,
    from_type = "ENSEMBL", to_type = "ENTREZID",
    orgdb = org.Hs.eg.db::org.Hs.eg.db
  ))

  expect_false(anyNA(names(rv)))
  expect_false(anyNA(rv))
  expect_false(anyDuplicated(names(rv)) > 0L)
  expect_identical(rv, sort(rv, decreasing = TRUE))

  # Every retained score must still equal the value for its own source gene.
  map <- VISTA:::.vista_map_ids_strict(
    tbl$gene_id, "ENSEMBL", "ENTREZID", orgdb = org.Hs.eg.db::org.Hs.eg.db
  )
  lookup <- stats::setNames(tbl$log2fc, map)
  for (nm in utils::head(names(rv), 10)) {
    expect_true(rv[[nm]] %in% tbl$log2fc[which(map == nm)], info = nm)
  }
})
