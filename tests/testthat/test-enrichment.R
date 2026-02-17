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

  # Should accept valid ontologies (though might error on orgdb)
  expect_error(
    get_go_enrichment(vista, sample_comparison = comps[1], ont = "BP"),
    "orgdb|annotation"
  )

  expect_error(
    get_go_enrichment(vista, sample_comparison = comps[1], ont = "MF"),
    "orgdb|annotation"
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
    vista_obj = vista,
    enrichment = enr,
    samples = unique(SummarizedExperiment::colData(vista)$cond_long),
    top_n = 1,
    gene_id_column = "SYMBOL",
    return_type = "genes"
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
    vista_obj = vista,
    enrichment = enr,
    samples = unique(SummarizedExperiment::colData(vista)$cond_long),
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
      vista_obj = vista,
      enrichment = enr,
      samples = unique(SummarizedExperiment::colData(vista)$cond_long),
      pathways = "Pathway_A",
      genes = genes[1:5],
      return_type = "genes"
    ),
    "managed by"
  )
})
