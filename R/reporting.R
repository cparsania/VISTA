#' Generate a publication-ready VISTA workflow report
#'
#' Builds a `VISTA` object (or uses a precomputed one), computes a comprehensive
#' single-comparison analysis panel, exports publication assets, and renders an
#' automated Quarto HTML/PDF report from YAML-driven parameters.
#'
#' The report focuses on one differential comparison (`primary_comparison`) and
#' includes:
#' \itemize{
#'   \item QC plots (PCA, MDS, correlation heatmap),
#'   \item DE plots (volcano, MA, DEG bar/pie/donut summaries),
#'   \item expression-focused views (boxplot, fold-change barplot, expression heatmap),
#'   \item enrichment outputs (MSigDB/GO/KEGG when available),
#'   \item downloadable artifacts (tables + plots + optional zip bundle),
#'   \item interactive HTML tables via \pkg{DT} when installed.
#' }
#'
#' @param config Path to a YAML config file (see template at
#'   `inst/reports/vista-report-template.yml`) or a named list of parameters.
#' @param output_file Optional output filename overriding `config$output_file`.
#'
#' @return Invisibly, the normalized output report path.
#' @export
run_vista_report <- function(config, output_file = "vista-report.html") {
  if (!requireNamespace("quarto", quietly = TRUE)) {
    cli::cli_abort("Quarto is required. Please install the {.pkg quarto} CLI.")
  }

  cfg <- if (is.character(config)) {
    if (!requireNamespace("yaml", quietly = TRUE)) {
      cli::cli_abort("Please install the {.pkg yaml} package to read YAML configs.")
    }
    yaml::read_yaml(config)
  } else {
    config
  }
  if (is.null(cfg) || !is.list(cfg)) {
    cli::cli_abort("{.arg config} must be a YAML path or a named list.")
  }

  cfg$report_title <- cfg$report_title %||% "VISTA RNA-seq Analysis Report"
  cfg$report_subtitle <- cfg$report_subtitle %||% NULL
  cfg$report_author <- cfg$report_author %||% NULL
  cfg$assets_dir <- cfg$assets_dir %||% "vista_report_assets"
  cfg$output_format <- tolower(cfg$output_format %||% "html")
  cfg$sidebar_position <- tolower(cfg$sidebar_position %||% "left")
  cfg$side_menu <- isTRUE(cfg$side_menu %||% TRUE)
  cfg$toc_depth <- as.integer(cfg$toc_depth %||% 3L)
  cfg$html_theme <- cfg$html_theme %||% "cosmo"
  cfg$embed_resources <- isTRUE(cfg$embed_resources %||% FALSE)

  cfg$de_method <- cfg$de_method %||% "deseq2"
  cfg$min_counts <- cfg$min_counts %||% 10
  cfg$min_replicates <- cfg$min_replicates %||% 1
  cfg$log2fc_cutoff <- cfg$log2fc_cutoff %||% 1
  cfg$pval_cutoff <- cfg$pval_cutoff %||% 0.05
  cfg$p_value_type <- cfg$p_value_type %||% "padj"
  cfg$covariates <- cfg$covariates %||% NULL
  if (!is.null(cfg$covariates)) {
    if (!length(cfg$covariates)) {
      cfg$covariates <- NULL
    } else {
      cfg$covariates <- as.character(unlist(cfg$covariates, use.names = FALSE))
    }
  }
  cfg$design_formula <- cfg$design_formula %||% NULL
  cfg$consensus_mode <- cfg$consensus_mode %||% "intersection"
  cfg$consensus_log2fc <- cfg$consensus_log2fc %||% "mean"
  cfg$de_source <- cfg$de_source %||% NULL
  cfg$group_palette <- cfg$group_palette %||% "Dark 2"

  cfg$heatmap_genes <- as.integer(cfg$heatmap_genes %||% 60L)
  cfg$heatmap_kmeans <- as.integer(cfg$heatmap_kmeans %||% 4L)
  cfg$annotate_columns <- isTRUE(cfg$annotate_columns %||% TRUE)
  cfg$summarise_replicates <- isTRUE(cfg$summarise_replicates %||% FALSE)
  cfg$top_table_n <- as.integer(cfg$top_table_n %||% 200L)
  cfg$top_genes_plot <- as.integer(cfg$top_genes_plot %||% 20L)
  cfg$include_deg_pie <- isTRUE(cfg$include_deg_pie %||% TRUE)
  cfg$include_deg_donut <- isTRUE(cfg$include_deg_donut %||% TRUE)

  cfg$include_msigdb <- isTRUE(cfg$include_msigdb %||% TRUE)
  cfg$include_go <- isTRUE(cfg$include_go %||% TRUE)
  cfg$include_kegg <- isTRUE(cfg$include_kegg %||% TRUE)
  cfg$msigdb_category <- cfg$msigdb_category %||% "H"
  cfg$msigdb_subcategory <- cfg$msigdb_subcategory %||% NULL
  cfg$species <- cfg$species %||% "Mus musculus"
  cfg$from_type <- cfg$from_type %||% NULL
  cfg$top_n_enrich <- as.integer(cfg$top_n_enrich %||% 15L)
  cfg$include_pathway_heatmap <- isTRUE(cfg$include_pathway_heatmap %||% TRUE)
  cfg$pathway_top_n <- as.integer(cfg$pathway_top_n %||% 3L)
  cfg$pathway_gene_mode <- cfg$pathway_gene_mode %||% "union"
  cfg$pathway_max_genes <- as.integer(cfg$pathway_max_genes %||% 80L)
  cfg$display_id <- cfg$display_id %||% NULL
  cfg$export_norm_counts <- isTRUE(cfg$export_norm_counts %||% FALSE)

  if (!cfg$sidebar_position %in% c("left", "right")) {
    cli::cli_abort("{.arg sidebar_position} must be one of {.val left} or {.val right}.")
  }

  if (missing(output_file) || is.null(output_file)) {
    output_file <- cfg$output_file %||% "vista-report.html"
  }
  if (!is.character(output_file) || length(output_file) != 1L || !nzchar(output_file)) {
    cli::cli_abort("{.arg output_file} must be a non-empty file path.")
  }

  out_dir <- dirname(output_file)
  out_file <- basename(output_file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  read_any <- function(path) {
    if (!file.exists(path)) cli::cli_abort("Input file not found: {.file {path}}")
    ext <- tolower(tools::file_ext(path))
    if (ext %in% c("rds", "rda", "rdata")) {
      readRDS(path)
    } else if (ext %in% c("tsv", "txt")) {
      utils::read.delim(path, check.names = FALSE)
    } else {
      utils::read.csv(path, check.names = FALSE)
    }
  }

  if (!is.null(cfg$vista_rds)) {
    vista_obj <- readRDS(cfg$vista_rds)
  } else {
    counts <- if (!is.null(cfg$counts_file)) read_any(cfg$counts_file) else cfg$counts
    sample_info_df <- if (!is.null(cfg$sample_info_file)) read_any(cfg$sample_info_file) else cfg$sample_info
    if (is.null(counts) || is.null(sample_info_df)) {
      cli::cli_abort(
        "Provide either {.arg vista_rds} or both {.arg counts_file}/{.arg sample_info_file} (or in-memory {.arg counts}/{.arg sample_info})."
      )
    }
    required <- c("column_geneid", "group_column", "group_numerator", "group_denominator")
    missing_required <- required[vapply(required, function(x) is.null(cfg[[x]]), logical(1))]
    if (length(missing_required) > 0) {
      cli::cli_abort("Missing required configuration fields: {.val {missing_required}}")
    }
    vista_obj <- create_vista(
      counts = counts,
      sample_info = sample_info_df,
      column_geneid = cfg$column_geneid,
      group_column = cfg$group_column,
      group_numerator = cfg$group_numerator,
      group_denominator = cfg$group_denominator,
      method = cfg$de_method,
      min_counts = cfg$min_counts,
      min_replicates = cfg$min_replicates,
      log2fc_cutoff = cfg$log2fc_cutoff,
      pval_cutoff = cfg$pval_cutoff,
      p_value_type = cfg$p_value_type,
      covariates = cfg$covariates,
      design_formula = cfg$design_formula,
      consensus_mode = cfg$consensus_mode,
      consensus_log2fc = cfg$consensus_log2fc,
      result_source = cfg$de_source,
      group_palette = cfg$group_palette
    )
  }

  comps <- names(comparisons(vista_obj))
  if (!length(comps)) cli::cli_abort("No differential comparisons found in the VISTA object.")

  primary_comp <- cfg$primary_comparison %||% comps[[1]]
  if (!primary_comp %in% comps) {
    cli::cli_abort(
      "{.arg primary_comparison} '{primary_comp}' not found. Available: {.val {comps}}"
    )
  }

  sanitize_name <- function(x) {
    out <- gsub("[^A-Za-z0-9_\\-]+", "_", as.character(x))
    out <- gsub("^_+|_+$", "", out)
    if (!nzchar(out)) "item" else out
  }
  comp_tag <- sanitize_name(primary_comp)

  assets_rel <- cfg$assets_dir
  plots_rel <- file.path(assets_rel, "plots")
  tables_rel <- file.path(assets_rel, "tables")
  plots_abs <- file.path(out_dir, plots_rel)
  tables_abs <- file.path(out_dir, tables_rel)
  dir.create(plots_abs, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_abs, recursive = TRUE, showWarnings = FALSE)

  save_plot <- function(obj, filename, width = 8, height = 6) {
    rel <- file.path(plots_rel, filename)
    abs <- file.path(out_dir, rel)
    if (inherits(obj, c("Heatmap", "HeatmapList"))) {
      grDevices::png(abs, width = width, height = height, units = "in", res = 300)
      ComplexHeatmap::draw(obj)
      grDevices::dev.off()
    } else {
      ggplot2::ggsave(abs, plot = obj, width = width, height = height, dpi = 300)
    }
    rel
  }

  write_table <- function(df, filename) {
    rel <- file.path(tables_rel, filename)
    abs <- file.path(out_dir, rel)
    utils::write.csv(df, abs, row.names = FALSE)
    rel
  }

  pick_first <- function(candidates, nms) {
    hit <- candidates[candidates %in% nms]
    if (length(hit)) hit[[1]] else NA_character_
  }

  de_tbl <- as.data.frame(comparisons(vista_obj)[[primary_comp]], stringsAsFactors = FALSE)
  if (!"gene_id" %in% colnames(de_tbl)) {
    de_tbl$gene_id <- rownames(de_tbl)
  }
  if (!"gene_id" %in% colnames(de_tbl)) {
    cli::cli_abort("DE results for {.val {primary_comp}} must contain a {.val gene_id} column.")
  }

  fc_col <- pick_first(c("log2fc", "log2FoldChange", "logFC"), colnames(de_tbl))
  p_col <- pick_first(c("p.adj", "padj", "FDR", "pvalue", "PValue"), colnames(de_tbl))

  if (!"regulation" %in% colnames(de_tbl) && !is.na(fc_col) && !is.na(p_col)) {
    fc_vals <- suppressWarnings(as.numeric(de_tbl[[fc_col]]))
    p_vals <- suppressWarnings(as.numeric(de_tbl[[p_col]]))
    p_vals[!is.finite(p_vals)] <- 1
    de_tbl$regulation <- ifelse(
      fc_vals >= cfg$log2fc_cutoff & p_vals <= cfg$pval_cutoff, "Up",
      ifelse(fc_vals <= -cfg$log2fc_cutoff & p_vals <= cfg$pval_cutoff, "Down", "Other")
    )
  }
  if (!"regulation" %in% colnames(de_tbl)) de_tbl$regulation <- "Other"

  order_de <- function(df) {
    if (!nrow(df)) return(df)
    pvals <- if (!is.na(p_col) && p_col %in% colnames(df)) suppressWarnings(as.numeric(df[[p_col]])) else rep(Inf, nrow(df))
    fcv <- if (!is.na(fc_col) && fc_col %in% colnames(df)) suppressWarnings(as.numeric(df[[fc_col]])) else rep(0, nrow(df))
    pvals[!is.finite(pvals)] <- Inf
    fcv[!is.finite(fcv)] <- 0
    df[order(pvals, -abs(fcv)), , drop = FALSE]
  }

  up_tbl <- order_de(de_tbl[de_tbl$regulation == "Up", , drop = FALSE])
  down_tbl <- order_de(de_tbl[de_tbl$regulation == "Down", , drop = FALSE])
  de_sorted <- order_de(de_tbl)

  top_up_genes <- head(unique(up_tbl$gene_id), cfg$top_genes_plot)
  top_down_genes <- head(unique(down_tbl$gene_id), cfg$top_genes_plot)
  heatmap_genes <- unique(c(top_up_genes, top_down_genes))
  if (length(heatmap_genes) < cfg$heatmap_genes && !is.na(fc_col) && fc_col %in% colnames(de_tbl)) {
    fc_abs <- suppressWarnings(abs(as.numeric(de_tbl[[fc_col]])))
    fc_abs[!is.finite(fc_abs)] <- 0
    fc_fill <- de_tbl$gene_id[order(fc_abs, decreasing = TRUE)]
    heatmap_genes <- unique(c(heatmap_genes, fc_fill))
  }
  heatmap_genes <- head(heatmap_genes, cfg$heatmap_genes)
  plot_genes <- unique(c(head(top_up_genes, 10), head(top_down_genes, 10)))
  if (!length(plot_genes)) plot_genes <- head(heatmap_genes, 10)
  if (!length(plot_genes)) plot_genes <- head(de_tbl$gene_id, 10)

  si <- as.data.frame(sample_info(vista_obj), stringsAsFactors = FALSE)
  group_col <- (S4Vectors::metadata(vista_obj)$group$column) %||% cfg$group_column %||% colnames(si)[1]
  if (!group_col %in% colnames(si)) {
    cli::cli_abort("Group column {.val {group_col}} not found in sample metadata.")
  }
  group_values <- unique(as.character(si[[group_col]]))
  group_values <- group_values[!is.na(group_values) & nzchar(group_values)]
  if (!length(group_values)) cli::cli_abort("No valid group values found for {.val {group_col}}.")

  summary_tbl <- data.frame(
    metric = c(
      "Generated",
      "Primary comparison",
      "Number of comparisons",
      "Samples",
      "Genes tested",
      "Upregulated genes",
      "Downregulated genes",
      "DE method",
      "LFC cutoff",
      "P-value cutoff",
      "P-value type"
    ),
    value = c(
      as.character(Sys.time()),
      primary_comp,
      as.character(length(comps)),
      as.character(ncol(norm_counts(vista_obj))),
      as.character(nrow(de_tbl)),
      as.character(nrow(up_tbl)),
      as.character(nrow(down_tbl)),
      as.character(cfg$de_method),
      as.character(cfg$log2fc_cutoff),
      as.character(cfg$pval_cutoff),
      as.character(cfg$p_value_type)
    ),
    stringsAsFactors = FALSE
  )

  table_paths <- list()
  table_paths$summary <- write_table(summary_tbl, paste0("summary_", comp_tag, ".csv"))
  table_paths$sample_info <- write_table(si, "sample_info.csv")
  table_paths$de_all <- write_table(de_tbl, paste0("de_results_", comp_tag, ".csv"))
  table_paths$de_top <- write_table(head(de_sorted, cfg$top_table_n), paste0("de_top_", comp_tag, ".csv"))
  table_paths$up <- write_table(head(up_tbl, cfg$top_table_n), paste0("de_up_", comp_tag, ".csv"))
  table_paths$down <- write_table(head(down_tbl, cfg$top_table_n), paste0("de_down_", comp_tag, ".csv"))
  if (cfg$export_norm_counts) {
    nc <- as.data.frame(norm_counts(vista_obj), stringsAsFactors = FALSE)
    nc$gene_id <- rownames(norm_counts(vista_obj))
    table_paths$norm_counts <- write_table(nc, "norm_counts.csv")
  }

  download_manifest <- data.frame(type = character(), label = character(), path = character(), stringsAsFactors = FALSE)
  register_download <- function(type, label, path) {
    if (!is.null(path) && is.character(path) && length(path) == 1L && nzchar(path)) {
      download_manifest <<- rbind(
        download_manifest,
        data.frame(type = type, label = label, path = path, stringsAsFactors = FALSE)
      )
    }
  }

  for (nm in names(table_paths)) {
    register_download("table", paste("Table:", nm), table_paths[[nm]])
  }

  plot_messages <- list()
  safe_plot <- function(slot, filename, plot_fun, width = 8, height = 6) {
    obj <- tryCatch(plot_fun(), error = function(e) {
      plot_messages[[slot]] <<- e$message
      NULL
    })
    if (is.null(obj)) return(NULL)
    tryCatch(save_plot(obj, filename, width = width, height = height), error = function(e) {
      plot_messages[[slot]] <<- e$message
      NULL
    })
  }

  plot_paths <- list()
  plot_paths$pca <- safe_plot("pca", "qc_pca.png", function() {
    get_pca_plot(vista_obj, label = TRUE)
  }, width = 9, height = 6)
  plot_paths$mds <- safe_plot("mds", "qc_mds.png", function() {
    get_mds_plot(vista_obj, label = TRUE)
  }, width = 9, height = 6)
  plot_paths$corr_heatmap <- safe_plot("corr_heatmap", "qc_corr_heatmap.png", function() {
    get_corr_heatmap(vista_obj)
  }, width = 8, height = 7)
  plot_paths$deg_bar <- safe_plot("deg_bar", "de_deg_count_barplot.png", function() {
    get_deg_count_barplot(vista_obj, facet_by = "comparison")
  }, width = 10, height = 6)
  if (cfg$include_deg_pie) {
    plot_paths$deg_pie <- safe_plot("deg_pie", "de_deg_count_pieplot.png", function() {
      get_deg_count_pieplot(vista_obj, facet_by = "comparison")
    }, width = 10, height = 5)
  }
  if (cfg$include_deg_donut) {
    plot_paths$deg_donut <- safe_plot("deg_donut", "de_deg_count_donutplot.png", function() {
      get_deg_count_donutplot(vista_obj, facet_by = "comparison")
    }, width = 10, height = 5)
  }
  plot_paths$volcano <- safe_plot("volcano", paste0("de_volcano_", comp_tag, ".png"), function() {
    get_volcano_plot(vista_obj, sample_comparison = primary_comp, display_id = cfg$display_id)
  }, width = 8, height = 6)
  plot_paths$ma <- safe_plot("ma", paste0("de_ma_", comp_tag, ".png"), function() {
    get_ma_plot(vista_obj, sample_comparison = primary_comp, display_id = cfg$display_id)
  }, width = 8, height = 6)
  plot_paths$expression_box <- safe_plot("expression_box", "expr_boxplot_top_genes.png", function() {
    get_expression_boxplot(vista_obj, genes = head(plot_genes, 8), display_id = cfg$display_id)
  }, width = 11, height = 6)
  plot_paths$foldchange_bar <- safe_plot("foldchange_bar", "de_foldchange_barplot_top_genes.png", function() {
    get_foldchange_barplot(
      vista_obj,
      genes = head(plot_genes, cfg$top_genes_plot),
      sample_comparisons = primary_comp,
      display_id = cfg$display_id
    )
  }, width = 12, height = 7)

  cluster_tbl <- NULL
  hm_return_type <- if (!is.na(cfg$heatmap_kmeans) && cfg$heatmap_kmeans >= 2) "both" else "heatmap"
  hm_annotate <- if (cfg$summarise_replicates) FALSE else cfg$annotate_columns
  if (length(heatmap_genes) >= 2) {
    hm_out <- tryCatch(
      get_expression_heatmap(
        x = vista_obj,
        sample_group = group_values,
        genes = heatmap_genes,
        value_transform = "zscore",
        annotate_columns = hm_annotate,
        summarise_replicates = cfg$summarise_replicates,
        show_row_names = FALSE,
        kmeans_k = if (hm_return_type == "both") cfg$heatmap_kmeans else NULL,
        return_type = hm_return_type,
        display_id = cfg$display_id
      ),
      error = function(e) {
        plot_messages$expression_heatmap <<- e$message
        NULL
      }
    )
    if (!is.null(hm_out)) {
      hm_obj <- hm_out
      if (is.list(hm_out)) {
        hm_obj <- hm_out$heatmap
        cluster_tbl <- hm_out$clusters
      }
      plot_paths$expression_heatmap <- tryCatch(
        save_plot(hm_obj, paste0("expr_heatmap_", comp_tag, ".png"), width = 8, height = 10),
        error = function(e) {
          plot_messages$expression_heatmap <<- e$message
          NULL
        }
      )
    }
  }
  if (!is.null(cluster_tbl) && nrow(cluster_tbl) > 0) {
    table_paths$clusters <- write_table(cluster_tbl, paste0("expr_heatmap_clusters_", comp_tag, ".csv"))
    register_download("table", "Table: clusters", table_paths$clusters)
  }

  for (nm in names(plot_paths)) {
    register_download("plot", paste("Plot:", nm), plot_paths[[nm]])
  }

  enrichment <- list(
    msigdb = list(plot = NULL, table = NULL, message = NULL),
    go = list(plot = NULL, table = NULL, message = NULL),
    kegg = list(plot = NULL, table = NULL, message = NULL)
  )
  enrichment_tables <- list(msigdb = NULL, go = NULL, kegg = NULL, pathway_genes = NULL)

  enrich_from <- cfg$from_type %||% cfg$display_id %||% if (any(grepl("^ENS", de_tbl$gene_id))) "ENSEMBL" else "SYMBOL"

  collect_enrichment <- function(kind, enrich_fun, file_stub) {
    out <- tryCatch(enrich_fun(), error = function(e) {
      enrichment[[kind]]$message <<- e$message
      NULL
    })
    if (is.null(out) || is.null(out$enrich)) {
      if (is.null(enrichment[[kind]]$message)) enrichment[[kind]]$message <<- "No enrichment output returned."
      return(NULL)
    }
    enr <- out$enrich
    if (!inherits(enr, c("enrichResult", "gseaResult")) || is.null(enr@result) || nrow(enr@result) == 0) {
      enrichment[[kind]]$message <<- "No significant enrichment terms."
      return(NULL)
    }
    tbl <- as.data.frame(enr@result, stringsAsFactors = FALSE)
    enrichment_tables[[kind]] <<- tbl
    enrichment[[kind]]$table <<- write_table(tbl, paste0(file_stub, "_", comp_tag, ".csv"))
    enrichment[[kind]]$plot <<- safe_plot(
      paste0("enrichment_", kind),
      paste0(file_stub, "_", comp_tag, ".png"),
      function() get_enrichment_plot(enr, top_n = cfg$top_n_enrich, title = paste(toupper(kind), "-", primary_comp)),
      width = 10,
      height = 7
    )
    register_download("table", paste("Table: enrichment", kind), enrichment[[kind]]$table)
    register_download("plot", paste("Plot: enrichment", kind), enrichment[[kind]]$plot)
    enr
  }

  msig_enr <- NULL
  if (cfg$include_msigdb) {
    msig_enr <- collect_enrichment("msigdb", function() {
      get_msigdb_enrichment(
        x = vista_obj,
        sample_comparison = primary_comp,
        regulation = "Both",
        from_type = enrich_from,
        orgdb = cfg$orgdb,
        msigdb_category = cfg$msigdb_category,
        msigdb_subcategory = cfg$msigdb_subcategory,
        species = cfg$species
      )
    }, "enrichment_msigdb")
  }
  if (cfg$include_go) {
    collect_enrichment("go", function() {
      get_go_enrichment(
        x = vista_obj,
        sample_comparison = primary_comp,
        regulation = "Both",
        from_type = enrich_from,
        orgdb = cfg$orgdb,
        species = cfg$species
      )
    }, "enrichment_go")
  }
  if (cfg$include_kegg) {
    collect_enrichment("kegg", function() {
      get_kegg_enrichment(
        x = vista_obj,
        sample_comparison = primary_comp,
        regulation = "Both",
        from_type = enrich_from,
        orgdb = cfg$orgdb,
        species = cfg$species
      )
    }, "enrichment_kegg")
  }

  if (!is.null(msig_enr) && "get_pathway_genes" %in% getNamespaceExports("VISTA")) {
    pathway_tbl <- tryCatch(
      get_pathway_genes(msig_enr, top_n = cfg$pathway_top_n, return_type = "long"),
      error = function(e) {
        enrichment$msigdb$message <- paste(c(enrichment$msigdb$message, paste("Pathway extraction:", e$message)), collapse = " | ")
        NULL
      }
    )
    if (!is.null(pathway_tbl) && nrow(pathway_tbl) > 0) {
      enrichment_tables$pathway_genes <- pathway_tbl
      table_paths$pathway_genes <- write_table(pathway_tbl, paste0("pathway_genes_", comp_tag, ".csv"))
      register_download("table", "Table: pathway genes", table_paths$pathway_genes)
    }
  }

  if (!is.null(msig_enr) &&
      cfg$include_pathway_heatmap &&
      "get_pathway_heatmap" %in% getNamespaceExports("VISTA")) {
    pathway_hm <- tryCatch(
      get_pathway_heatmap(
        x = vista_obj,
        enrichment = msig_enr,
        sample_group = group_values,
        top_n = cfg$pathway_top_n,
        gene_mode = cfg$pathway_gene_mode,
        max_genes = cfg$pathway_max_genes,
        value_transform = "zscore",
        display_id = cfg$display_id,
        annotate_columns = hm_annotate,
        summarise_replicates = cfg$summarise_replicates,
        show_row_names = FALSE
      ),
      error = function(e) {
        plot_messages$pathway_heatmap <- e$message
        NULL
      }
    )
    if (!is.null(pathway_hm)) {
      plot_paths$pathway_heatmap <- tryCatch(
        save_plot(pathway_hm, paste0("expr_pathway_heatmap_", comp_tag, ".png"), width = 9, height = 10),
        error = function(e) {
          plot_messages$pathway_heatmap <- e$message
          NULL
        }
      )
      register_download("plot", "Plot: pathway_heatmap", plot_paths$pathway_heatmap)
    }
  }

  zip_rel <- NULL
  zip_msg <- NULL
  rel_files <- unique(download_manifest$path)
  if (length(rel_files) > 0) {
    zip_rel <- file.path(assets_rel, paste0("download_bundle_", comp_tag, ".zip"))
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(out_dir)
    tryCatch(
      utils::zip(zipfile = zip_rel, files = rel_files),
      error = function(e) {
        zip_msg <<- e$message
      }
    )
    if (!file.exists(file.path(out_dir, zip_rel))) {
      zip_rel <- NULL
    }
  }
  if (!is.null(zip_rel)) {
    register_download("archive", "Download all artifacts (zip)", zip_rel)
  }

  bundle <- list(
    report = list(
      title = cfg$report_title,
      subtitle = cfg$report_subtitle,
      author = cfg$report_author,
      generated = as.character(Sys.time()),
      primary_comp = primary_comp,
      group_column = group_col,
      de_method = cfg$de_method
    ),
    summary_table = summary_tbl,
    plots = plot_paths,
    tables = list(
      sample_info = si,
      de_top = head(de_sorted, cfg$top_table_n),
      de_up = head(up_tbl, cfg$top_table_n),
      de_down = head(down_tbl, cfg$top_table_n),
      clusters = cluster_tbl,
      enrichment_msigdb = enrichment_tables$msigdb,
      enrichment_go = enrichment_tables$go,
      enrichment_kegg = enrichment_tables$kegg,
      pathway_genes = enrichment_tables$pathway_genes
    ),
    table_paths = table_paths,
    enrichment = enrichment,
    messages = list(plot = plot_messages, zip = zip_msg),
    downloads = download_manifest,
    zip_bundle = zip_rel
  )

  data_rds <- tempfile(fileext = ".rds")
  saveRDS(bundle, data_rds)
  data_rds <- normalizePath(data_rds, winslash = "/", mustWork = FALSE)

  esc_yaml <- function(x) gsub("\"", "\\\\\"", as.character(x), fixed = TRUE)
  fmt_lines <- if (cfg$output_format %in% c("html", "html-default", "html4")) {
    Filter(Negate(is.null), c(
      "format:",
      "  html:",
      paste0("    toc: ", tolower(as.character(cfg$side_menu))),
      if (cfg$side_menu) paste0("    toc-location: ", cfg$sidebar_position) else NULL,
      paste0("    toc-depth: ", cfg$toc_depth),
      "    number-sections: true",
      paste0("    theme: ", cfg$html_theme),
      "    code-fold: false",
      paste0("    embed-resources: ", tolower(as.character(cfg$embed_resources)))
    ))
  } else {
    paste0("format: ", cfg$output_format)
  }

  qmd <- file.path(out_dir, paste0("vista-report-", format(Sys.time(), "%Y%m%d%H%M%S"), ".qmd"))
  qmd_lines <- Filter(Negate(is.null), c(
    "---",
    paste0("title: \"", esc_yaml(cfg$report_title), "\""),
    if (!is.null(cfg$report_subtitle) && nzchar(cfg$report_subtitle)) paste0("subtitle: \"", esc_yaml(cfg$report_subtitle), "\"") else NULL,
    if (!is.null(cfg$report_author) && nzchar(cfg$report_author)) paste0("author: \"", esc_yaml(cfg$report_author), "\"") else NULL,
    fmt_lines,
    "resources:",
    paste0("  - ", assets_rel, "/**"),
    "params:",
    "  bundle_rds: null",
    "---",
    "",
    "```{r setup}",
    "bundle <- readRDS(params$bundle_rds)",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "render_table <- function(df, caption = NULL) {",
    "  if (is.null(df) || !nrow(df)) {",
    "    cat('No data available.\\n\\n')",
    "    return(invisible(NULL))",
    "  }",
    "  if (requireNamespace('DT', quietly = TRUE)) {",
    "    print(DT::datatable(",
    "      df,",
    "      rownames = FALSE,",
    "      caption = caption,",
    "      extensions = 'Buttons',",
    "      options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel'),",
    "                     pageLength = 15, scrollX = TRUE)",
    "    ))",
    "  } else {",
    "    print(knitr::kable(utils::head(df, 50), caption = caption))",
    "  }",
    "}",
    "include_plot <- function(path) {",
    "  if (is.null(path) || !nzchar(path) || !file.exists(path)) {",
    "    cat('Plot unavailable.\\n\\n')",
    "    return(invisible(NULL))",
    "  }",
    "  knitr::include_graphics(path)",
    "}",
    "print_links <- function(df) {",
    "  if (is.null(df) || !nrow(df)) {",
    "    cat('No downloadable artifacts available.\\n\\n')",
    "    return(invisible(NULL))",
    "  }",
    "  for (i in seq_len(nrow(df))) {",
    "    cat(sprintf('- [%s](%s)\\n', df$label[i], df$path[i]))",
    "  }",
    "  cat('\\n')",
    "}",
    "```",
    "",
    "## Executive Summary",
    "",
    "```{r}",
    "cat('**Primary comparison:** ', bundle$report$primary_comp, '\\n\\n', sep = '')",
    "cat('**DE method:** ', bundle$report$de_method, '\\n\\n', sep = '')",
    "cat('**Group column:** ', bundle$report$group_column, '\\n\\n', sep = '')",
    "render_table(bundle$summary_table, caption = 'Analysis summary')",
    "```",
    "",
    "## Quality Control",
    "",
    "### PCA",
    "```{r}",
    "include_plot(bundle$plots$pca)",
    "```",
    "",
    "### MDS",
    "```{r}",
    "include_plot(bundle$plots$mds)",
    "```",
    "",
    "### Correlation Heatmap",
    "```{r}",
    "include_plot(bundle$plots$corr_heatmap)",
    "```",
    "",
    "## Differential Expression",
    "",
    "### DEG Count Summary",
    "```{r}",
    "include_plot(bundle$plots$deg_bar)",
    "```",
    "",
    "### DEG Composition (Pie)",
    "```{r}",
    "include_plot(bundle$plots$deg_pie)",
    "```",
    "",
    "### DEG Composition (Donut)",
    "```{r}",
    "include_plot(bundle$plots$deg_donut)",
    "```",
    "",
    "### Volcano Plot",
    "```{r}",
    "include_plot(bundle$plots$volcano)",
    "```",
    "",
    "### MA Plot",
    "```{r}",
    "include_plot(bundle$plots$ma)",
    "```",
    "",
    "### Top Differentially Expressed Genes",
    "```{r}",
    "render_table(bundle$tables$de_top, caption = 'Top DE genes')",
    "```",
    "",
    "## Expression-Centric Views",
    "",
    "### Expression Heatmap",
    "```{r}",
    "include_plot(bundle$plots$expression_heatmap)",
    "```",
    "",
    "### Pathway Gene Heatmap",
    "```{r}",
    "include_plot(bundle$plots$pathway_heatmap)",
    "```",
    "",
    "### Fold-Change Barplot",
    "```{r}",
    "include_plot(bundle$plots$foldchange_bar)",
    "```",
    "",
    "### Expression Boxplot",
    "```{r}",
    "include_plot(bundle$plots$expression_box)",
    "```",
    "",
    "## Functional Enrichment",
    "",
    "### MSigDB",
    "```{r}",
    "include_plot(bundle$enrichment$msigdb$plot)",
    "if (!is.null(bundle$enrichment$msigdb$message) && nzchar(bundle$enrichment$msigdb$message)) {",
    "  cat('**Note:** ', bundle$enrichment$msigdb$message, '\\n\\n', sep = '')",
    "}",
    "render_table(bundle$tables$enrichment_msigdb, caption = 'MSigDB enrichment')",
    "```",
    "",
    "### GO",
    "```{r}",
    "include_plot(bundle$enrichment$go$plot)",
    "if (!is.null(bundle$enrichment$go$message) && nzchar(bundle$enrichment$go$message)) {",
    "  cat('**Note:** ', bundle$enrichment$go$message, '\\n\\n', sep = '')",
    "}",
    "render_table(bundle$tables$enrichment_go, caption = 'GO enrichment')",
    "```",
    "",
    "### KEGG",
    "```{r}",
    "include_plot(bundle$enrichment$kegg$plot)",
    "if (!is.null(bundle$enrichment$kegg$message) && nzchar(bundle$enrichment$kegg$message)) {",
    "  cat('**Note:** ', bundle$enrichment$kegg$message, '\\n\\n', sep = '')",
    "}",
    "render_table(bundle$tables$enrichment_kegg, caption = 'KEGG enrichment')",
    "```",
    "",
    "### Pathway Gene Membership",
    "```{r}",
    "render_table(bundle$tables$pathway_genes, caption = 'Genes driving top enriched pathways')",
    "```",
    "",
    "## Supplementary Tables",
    "",
    "### Upregulated Genes",
    "```{r}",
    "render_table(bundle$tables$de_up, caption = 'Top upregulated genes')",
    "```",
    "",
    "### Downregulated Genes",
    "```{r}",
    "render_table(bundle$tables$de_down, caption = 'Top downregulated genes')",
    "```",
    "",
    "### Heatmap Clusters",
    "```{r}",
    "render_table(bundle$tables$clusters, caption = 'Expression heatmap clusters')",
    "```",
    "",
    "## Downloads",
    "",
    "```{r}",
    "print_links(bundle$downloads)",
    "if (!is.null(bundle$zip_bundle) && nzchar(bundle$zip_bundle) && file.exists(bundle$zip_bundle)) {",
    "  cat(sprintf('- [Download all artifacts (zip)](%s)\\n\\n', bundle$zip_bundle))",
    "}",
    "if (!is.null(bundle$messages$zip) && nzchar(bundle$messages$zip)) {",
    "  cat('Zip bundle note: ', bundle$messages$zip, '\\n\\n', sep = '')",
    "}",
    "```",
    "",
    "## Session Info",
    "",
    "```{r}",
    "sessionInfo()",
    "```"
  ))
  writeLines(qmd_lines, qmd)

  quarto::quarto_render(
    input = qmd,
    output_file = out_file,
    execute_params = list(bundle_rds = data_rds)
  )

  generated <- file.path(out_dir, out_file)
  if (!file.exists(generated) && file.exists(out_file)) generated <- out_file
  invisible(normalizePath(generated, winslash = "/", mustWork = FALSE))
}
