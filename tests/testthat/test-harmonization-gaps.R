# Phase 3g: the last plots that hand-rolled logic the rest of the family shares.
#
# Routing them through the shared resolvers surfaced three defects that the
# duplicated code had been hiding. These tests pin each one.

test_that("fold-change heatmap annotates columns by comparison, not by sample", {
  skip_if_not_installed("ComplexHeatmap")
  v <- make_small_vista()

  # Columns of a fold-change heatmap are comparisons. The previous code filtered
  # colData by `sample %in% colnames(fc_mat)`, which can never match, so the
  # annotation was built from an empty frame: a track with no levels and no
  # colours. It looked fine and told the reader nothing.
  fc_cols <- colnames(get_foldchange_matrix(v))
  expect_length(intersect(fc_cols, colnames(v)), 0)

  hm <- get_foldchange_heatmap(v, annotate_columns = TRUE, return_type = "plot")
  anno <- hm@top_annotation
  expect_false(is.null(anno))
  expect_identical(names(anno@anno_list), "comparison")

  cm <- anno@anno_list[[1]]@color_mapping
  expect_setequal(as.character(cm@levels), fc_cols)
  expect_length(cm@colors, length(fc_cols))
  expect_true(all(nzchar(cm@colors)))
})

test_that("fold-change heatmap column colours can be overridden", {
  skip_if_not_installed("ComplexHeatmap")
  v <- make_small_vista()
  comp <- names(comparisons(v))[[1]]

  hm <- get_foldchange_heatmap(
    v, annotate_columns = TRUE, return_type = "plot",
    column_anno_colors = stats::setNames(list(stats::setNames("#FF0000", comp)), "comparison")
  )
  cols <- hm@top_annotation@anno_list[[1]]@color_mapping@colors
  expect_true(any(grepl("^#FF0000", cols, ignore.case = TRUE)))

  expect_error(
    get_foldchange_heatmap(v, annotate_columns = TRUE, column_anno_colors = "red"),
    "named list"
  )

  # off by default, and unchanged when off
  expect_null(get_foldchange_heatmap(v, return_type = "plot")@top_annotation)
})

test_that("a display_id without a display_orgdb says what is missing", {
  # `to_type == from_type` is NA when from_type is NA and logical(0) when it is
  # NULL; `||` can evaluate neither, so this failed with the opaque
  # "missing value where TRUE/FALSE needed".
  v <- make_small_vista()
  expect_error(
    get_foldchange_raincloud(v, genes = rownames(v)[1:5], display_id = "NOPE"),
    "display_orgdb"
  )

  # the passthrough cases must keep short-circuiting
  ids <- c("a", "b")
  expect_identical(VISTA:::.map_gene_ids(ids, to_type = NULL), ids)
  expect_identical(VISTA:::.map_gene_ids(ids, to_type = "none"), ids)
  expect_identical(VISTA:::.map_gene_ids(ids, to_type = ""), ids)
  expect_identical(VISTA:::.map_gene_ids(ids, from_type = "X", to_type = "X"), ids)
  expect_identical(VISTA:::.map_gene_ids(ids, to_type = NA_character_), ids)
  expect_identical(
    unname(VISTA:::.map_gene_ids(ids, to_type = "SYMBOL", id_map = c(a = "A", b = "B"))),
    c("A", "B")
  )
})

test_that("boxplot faceting keeps its own auto rule after adopting the resolver", {
  # The shared resolver answers on the genes actually plotted; this plot draws
  # every gene when `genes = NULL`, where faceting all of them is never what
  # "auto" meant. Pooling collapses genes, so its auto answer is "group".
  v <- make_small_vista()
  g3 <- rownames(v)[1:3]

  facet_of <- function(p) {
    if (inherits(p$facet, "FacetNull")) return("none")
    paste(names(p$facet$params$facets), collapse = "+")
  }
  gcol <- S4Vectors::metadata(v)$group$column

  expect_identical(facet_of(get_expression_boxplot(v, genes = g3)), "gene")
  expect_identical(facet_of(get_expression_boxplot(v, genes = g3[1])), "none")
  expect_identical(facet_of(get_expression_boxplot(v)), "none")            # all genes
  expect_identical(facet_of(get_expression_boxplot(v, genes = g3, facet_by = "none")), "none")
  expect_identical(
    facet_of(get_expression_boxplot(v, genes = g3, pool_genes = TRUE, by = "group")), gcol
  )
  expect_identical(
    facet_of(get_expression_boxplot(v, genes = g3, pool_genes = TRUE, by = "group",
                                    facet_by = "none")), "none"
  )
})

test_that("every gene-taking plot resolves identifiers through one helper", {
  # Three byte-similar copies of the mapping had drifted apart; one silently
  # ignored display_from/display_orgdb. Guard against a fourth appearing.
  # Only meaningful against the sources, which are absent from an install.
  src_file <- testthat::test_path("..", "..", "R", "viz_related.R")
  skip_if_not(file.exists(src_file), "package sources not available")

  src <- readLines(src_file, warn = FALSE)
  expect_false(any(grepl("^\\s*genes_use <- genes\\s*$", src)))
  # and the resolver really is the single point of resolution
  expect_gt(sum(grepl("\\.resolve_foldchange_gene_ids\\(", src)), 5L)
})
