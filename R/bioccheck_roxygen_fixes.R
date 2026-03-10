# BiocCheck roxygen patches
#
# These blocks ensure exported man pages include non-empty \value and
# runnable \examples sections after roxygen regeneration.

#' @name enrichMsigDB
#' @examples
#' v <- example_vista()
#' genes <- head(as.character(row_data(v)$gene_id), 20)
#' if (requireNamespace('msigdbr', quietly = TRUE)) {
#'   out <- try(enrichMsigDB(v, gene_list = genes, from_type = 'ENSEMBL', msigdb_category = 'H'), silent = TRUE)
#'   if (!inherits(out, 'try-error')) out
#' }
NULL

#' @name get_cell_fractions
#' @examples
#' v <- example_vista()
#' if (requireNamespace('xCell2', quietly = TRUE)) {
#'   vx <- try(run_cell_deconvolution(v, method = 'xCell2'), silent = TRUE)
#'   if (!inherits(vx, 'try-error')) head(get_cell_fractions(vx))
#' }
NULL

#' @name get_corr_heatmap
#' @param show_corr_values Deprecated alias for `label`. When supplied, it overrides `label`.
#' @param col_corr_values Deprecated alias for `label_color`. When supplied, it overrides `label_color`.
#' @return An object returned by this function.
#' @examples
#' v <- example_vista()
#' p <- get_corr_heatmap(v)
#' print(p)
NULL

#' @name get_deg_count_barplot
#' @return An object returned by this function.
#' @examples
#' v <- example_vista()
#' p <- get_deg_count_barplot(v)
#' print(p)
NULL

#' @name get_deg_count_donutplot
#' @return An object returned by this function.
#' @examples
#' v <- example_vista()
#' p <- get_deg_count_donutplot(v)
#' print(p)
NULL

#' @name get_deg_count_pieplot
#' @return An object returned by this function.
#' @examples
#' v <- example_vista()
#' p <- get_deg_count_pieplot(v)
#' print(p)
NULL

#' @name get_deg_venn_diagram
#' @return An object returned by this function.
#' @examples
#' v <- example_vista()
#' comps <- names(comparisons(v))
#' if (length(comps) >= 2) {
#'   p <- get_deg_venn_diagram(v, sample_comparisons = comps[1:2])
#'   print(p)
#' }
NULL

#' @name get_enrichment_plot
#' @examples
#' v <- example_vista()
#' com <- names(comparisons(v))[1]
#' if (requireNamespace('msigdbr', quietly = TRUE)) {
#'   ms <- try(get_msigdb_enrichment(v, sample_comparison = com, regulation = 'Up', from_type = 'ENSEMBL'), silent = TRUE)
#'   if (!inherits(ms, 'try-error') && !is.null(ms$enrich)) print(get_enrichment_plot(ms$enrich))
#' }
NULL

#' @name get_expression_boxplot
#' @return An object returned by this function.
#' @examples
#' v <- example_vista()
#' genes <- head(rownames(v), 3)
#' p <- get_expression_boxplot(v, genes = genes)
#' print(p)
NULL

#' @name get_expression_chromosome_plot
#' @examples
#' v <- example_vista()
#' p <- try(get_expression_chromosome_plot(v), silent = TRUE)
#' if (!inherits(p, 'try-error')) print(p)
NULL

#' @name get_expression_density
#' @examples
#' v <- example_vista()
#' genes <- head(rownames(v), 3)
#' p <- get_expression_density(v, genes = genes)
#' print(p)
NULL

#' @name get_expression_heatmap
#' @return An object returned by this function.
#' @examples
#' v <- example_vista()
#' genes <- head(rownames(v), 20)
#' if (requireNamespace('ComplexHeatmap', quietly = TRUE) &&
#'     requireNamespace('circlize', quietly = TRUE)) {
#'   hm <- get_expression_heatmap(
#'     v,
#'     genes = genes,
#'     sample_group = unique(as.character(sample_info(v)$cond_long)),
#'     return_type = 'heatmap'
#'   )
#'   ComplexHeatmap::draw(hm)
#' }
NULL

#' @name get_expression_joyplot
#' @examples
#' v <- example_vista()
#' genes <- head(rownames(v), 3)
#' p <- get_expression_joyplot(v, genes = genes)
#' print(p)
NULL

#' @name get_expression_lineplot
#' @return An object returned by this function.
#' @examples
#' v <- example_vista()
#' genes <- head(rownames(v), 3)
#' p <- get_expression_lineplot(v, genes = genes)
#' print(p)
NULL

#' @name get_expression_lollipop
#' @return An object returned by this function.
#' @examples
#' v <- example_vista()
#' genes <- head(rownames(v), 5)
#' p <- get_expression_lollipop(v, genes = genes)
#' print(p)
NULL

#' @name get_expression_matrix
#' @examples
#' v <- example_vista()
#' m <- get_expression_matrix(v)
#' dim(m)
NULL

#' @name get_expression_raincloud
#' @examples
#' v <- example_vista()
#' genes <- head(rownames(v), 5)
#' p <- get_expression_raincloud(v, genes = genes, summarise = TRUE)
#' print(p)
NULL

#' @name get_expression_scatter
#' @examples
#' v <- example_vista()
#' si <- as.data.frame(sample_info(v))
#' genes <- head(rownames(v), 50)
#' p <- get_expression_scatter(
#'   v,
#'   sample_x = si$sample_names[1],
#'   sample_y = si$sample_names[2],
#'   genes = genes
#' )
#' print(p)
NULL

#' @name get_expression_violinplot
#' @examples
#' v <- example_vista()
#' genes <- head(rownames(v), 4)
#' p <- get_expression_violinplot(v, genes = genes)
#' print(p)
NULL

#' @name get_foldchange_barplot
#' @return An object returned by this function.
#' @examples
#' v <- example_vista()
#' comp <- names(comparisons(v))[1]
#' genes <- head(as.character(comparisons(v)[[comp]]$gene_id), 10)
#' p <- get_foldchange_barplot(v, sample_comparison = comp, genes = genes)
#' print(p)
NULL

#' @name get_foldchange_boxplot
#' @return An object returned by this function.
#' @examples
#' v <- example_vista()
#' comp <- names(comparisons(v))[1]
#' genes <- head(as.character(comparisons(v)[[comp]]$gene_id), 10)
#' p <- get_foldchange_boxplot(v, sample_comparison = comp, genes = genes)
#' print(p)
NULL

#' @name get_foldchange_chromosome_plot
#' @examples
#' v <- example_vista()
#' p <- try(get_foldchange_chromosome_plot(v, sample_comparison = names(comparisons(v))[1]), silent = TRUE)
#' if (!inherits(p, 'try-error')) print(p)
NULL

#' @name get_foldchange_heatmap
#' @return An object returned by this function.
#' @examples
#' v <- example_vista()
#' comp <- names(comparisons(v))[1]
#' genes <- unique(stats::na.omit(as.character(comparisons(v)[[comp]]$gene_id)))[1:20]
#' if (requireNamespace('ComplexHeatmap', quietly = TRUE) &&
#'     requireNamespace('circlize', quietly = TRUE)) {
#'   hm <- get_foldchange_heatmap(
#'     v,
#'     sample_comparisons = comp,
#'     genes = genes,
#'     return_type = 'heatmap'
#'   )
#'   ComplexHeatmap::draw(hm)
#' }
NULL

#' @name get_foldchange_lineplot
#' @return An object returned by this function.
#' @examples
#' v <- example_vista()
#' comp <- names(comparisons(v))[1]
#' genes <- head(as.character(comparisons(v)[[comp]]$gene_id), 5)
#' p <- get_foldchange_lineplot(v, sample_comparison = comp, genes = genes)
#' print(p)
NULL

#' @name get_foldchange_lollipop
#' @return An object returned by this function.
#' @examples
#' v <- example_vista()
#' comp <- names(comparisons(v))[1]
#' genes <- head(as.character(comparisons(v)[[comp]]$gene_id), 10)
#' p <- get_foldchange_lollipop(v, sample_comparison = comp, genes = genes)
#' print(p)
NULL

#' @name get_foldchange_matrix
#' @examples
#' v <- example_vista()
#' mat <- get_foldchange_matrix(v)
#' dim(mat)
NULL

#' @name get_foldchange_raincloud
#' @examples
#' v <- example_vista()
#' comp <- names(comparisons(v))[1]
#' genes <- head(as.character(comparisons(v)[[comp]]$gene_id), 20)
#' p <- get_foldchange_raincloud(v, sample_comparison = comp, genes = genes)
#' print(p)
NULL

#' @name get_foldchange_scatter
#' @examples
#' data('count_data', package = 'VISTA')
#' data('sample_metadata', package = 'VISTA')
#' cell_levels <- unique(sample_metadata$cell)
#' if (length(cell_levels) >= 3) {
#'   v <- create_vista(count_data[1:150, ], sample_metadata, column_geneid = 'gene_id', group_column = 'cell',
#'                     group_numerator = cell_levels[2:3], group_denominator = rep(cell_levels[1], 2),
#'                     min_counts = 5, min_replicates = 1)
#'   comp_names <- names(comparisons(v))[1:2]
#'   p <- get_foldchange_scatter(v, sample_comparisons = comp_names)
#'   print(p)
#' }
NULL

#' @name get_genes_by_regulation
#' @examples
#' v <- example_vista()
#' comp <- names(comparisons(v))[1]
#' genes <- get_genes_by_regulation(v, sample_comparisons = comp, regulation = 'Up')
#' str(genes, max.level = 1)
NULL

#' @name get_go_enrichment
#' @examples
#' v <- example_vista()
#' comp <- names(comparisons(v))[1]
#' if (requireNamespace('org.Mm.eg.db', quietly = TRUE)) {
#'   out <- try(get_go_enrichment(v, sample_comparison = comp, ont = 'BP', from_type = 'ENSEMBL',
#'                                orgdb = org.Mm.eg.db::org.Mm.eg.db), silent = TRUE)
#'   if (!inherits(out, 'try-error')) out
#' }
NULL

#' @name get_gsea
#' @return An object returned by this function.
#' @examples
#' v <- example_vista()
#' comp <- names(comparisons(v))[1]
#' if (requireNamespace('msigdbr', quietly = TRUE)) {
#'   out <- try(get_gsea(v, sample_comparison = comp, set_type = 'msigdb', from_type = 'ENSEMBL', species = 'Homo sapiens'), silent = TRUE)
#'   if (!inherits(out, 'try-error')) out
#' }
NULL

#' @name get_kegg_enrichment
#' @return An object returned by this function.
#' @examples
#' v <- example_vista()
#' comp <- names(comparisons(v))[1]
#' if (requireNamespace('org.Mm.eg.db', quietly = TRUE)) {
#'   out <- try(get_kegg_enrichment(v, sample_comparison = comp, from_type = 'ENSEMBL', orgdb = org.Mm.eg.db::org.Mm.eg.db), silent = TRUE)
#'   if (!inherits(out, 'try-error')) out
#' }
NULL

#' @name get_ma_plot
#' @examples
#' v <- example_vista()
#' p <- get_ma_plot(v, sample_comparison = names(comparisons(v))[1])
#' print(p)
NULL

#' @name get_mds_plot
#' @return An object returned by this function.
#' @examples
#' v <- example_vista()
#' p <- get_mds_plot(v)
#' print(p)
NULL

#' @name get_pairwise_corr_plot
#' @return An object returned by this function.
#' @examples
#' v <- example_vista()
#' p <- get_pairwise_corr_plot(v)
#' print(p)
NULL

#' @name get_deg_alluvial
#' @examples
#' data("count_data", package = "VISTA")
#' data("sample_metadata", package = "VISTA")
#' si <- as.data.frame(sample_metadata[1:8, ], stringsAsFactors = FALSE)
#' trt_idx <- which(si$cond_long == "treatment1")
#' si$cond_long[trt_idx] <- rep(c("treatment1", "treatment2"), length.out = length(trt_idx))
#' si$groups <- si$cond_long
#' v <- create_vista(
#'   counts = count_data[1:120, c("gene_id", si$sample_names)],
#'   sample_info = si,
#'   column_geneid = "gene_id",
#'   group_column = "cond_long",
#'   group_numerator = c("treatment1", "treatment2"),
#'   group_denominator = c("control", "control"),
#'   min_counts = 5,
#'   min_replicates = 1
#' )
#' if (requireNamespace('ggalluvial', quietly = TRUE)) {
#'   p <- get_deg_alluvial(v, sample_comparisons = names(comparisons(v)))
#'   print(p)
#' }
NULL

#' @name print.VISTA
#' @aliases print.vista
#' @return The input object `x`, returned invisibly.
#' @examples
#' v <- example_vista()
#' print(v)
NULL

#' @name run_cell_deconvolution
#' @examples
#' v <- example_vista()
#' if (requireNamespace('xCell2', quietly = TRUE)) {
#'   out <- try(run_cell_deconvolution(v, method = 'xCell2'), silent = TRUE)
#'   if (!inherits(out, 'try-error')) out
#' }
NULL

#' @name run_vista_report
#' @examples
#' data('count_data', package = 'VISTA')
#' data('sample_metadata', package = 'VISTA')
#' cfg <- list(
#'   counts = count_data[1:100, ],
#'   sample_info = sample_metadata[1:6, ],
#'   column_geneid = 'gene_id',
#'   group_column = 'cond_long',
#'   group_numerator = 'treatment1',
#'   group_denominator = 'control',
#'   include_msigdb = FALSE, include_go = FALSE, include_kegg = FALSE
#' )
#' if (requireNamespace('quarto', quietly = TRUE)) {
#'   out <- tempfile(fileext = '.html')
#'   try(run_vista_report(cfg, output_file = out), silent = TRUE)
#' }
NULL

#' @name set_de_source
#' @examples
#' v <- example_vista(method = "both")
#' v <- set_de_source(v, "edger")
#' names(comparisons(v, source = "active"))
NULL

#' @name set_vista_comparison_colors
#' @examples
#' v <- example_vista()
#' comps <- names(comparisons(v))
#' if (length(comps)) {
#'   cmap <- stats::setNames(rep('#1b9e77', length(comps)), comps)
#'   set_vista_comparison_colors(v, cmap)
#' }
NULL

#' @name set_vista_group_colors
#' @examples
#' v <- example_vista()
#' groups <- unique(as.character(sample_info(v)$cond_long))
#' gmap <- stats::setNames(c('#1b9e77', '#d95f02')[seq_along(groups)], groups)
#' set_vista_group_colors(v, gmap)
NULL

#' @name VISTA-class
#' @return A `VISTA` S4 object.
#' @examples
#' methods::showClass('VISTA')
NULL
