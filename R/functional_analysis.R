#' @importFrom clusterProfiler enricher GSEA gseGO gseKEGG enrichGO enrichKEGG
#' @importFrom msigdbr msigdbr
#' @importFrom AnnotationDbi mapIds keys
#' @importFrom methods setGeneric setMethod
#' @importFrom SummarizedExperiment rowData
#' @importFrom cli cli_abort cli_warn
NULL

#' @keywords internal
.vista_msigdbr <- function(species, msigdb_category = "H", msigdb_subcategory = NULL) {
  msig_args <- list(species = species)
  msig_formals <- names(formals(msigdbr::msigdbr))

  if ("collection" %in% msig_formals) {
    msig_args$collection <- msigdb_category
  } else {
    msig_args$category <- msigdb_category
  }

  if (!is.null(msigdb_subcategory)) {
    if ("subcollection" %in% msig_formals) {
      msig_args$subcollection <- msigdb_subcategory
    } else {
      msig_args$subcategory <- msigdb_subcategory
    }
  }

  do.call(msigdbr::msigdbr, msig_args)
}

#' Perform MSigDB over-representation analysis on a VISTA object
#'
#' This function performs ORA on a set of genes from a VISTA object.  Gene
#' identifiers stored as “ENSG00000187634:SAMD11” will be split at the colon and
#' the appropriate part extracted based on `from_type`.  Genes are then
#' optionally converted via `orgdb` to match the identifier type required by
#' the MSigDB gene sets.
#'
#' @param x A VISTA object.
#' @param gene_list A non-empty character vector of gene identifiers.  This
#'   argument is required and must not be `NULL`.
#' @param from_type Identifier type of the genes in `gene_list`. One of
#'   `"SYMBOL"`, `"ENSEMBL"`, or `"ENTREZID"` (default `"SYMBOL"`). Ensembl IDs
#'   have version suffixes stripped automatically.
#' @param orgdb An `OrgDb` object used for ID conversion. If omitted, a default
#'   is chosen based on `species`: `org.Mm.eg.db` for mouse, `org.Hs.eg.db` for
#'   human.
#' @param msigdb_category MSigDB category (e.g. `"H"`, `"C2"`, `"C5"`). Default `"H"`.
#' @param msigdb_subcategory Optional MSigDB sub-collection (e.g. `"BP"` for C5). Default `NULL`.
#' @param species Species name for MSigDB (default `"Mus musculus"`).
#' @param background Optional character vector of background gene IDs. If `NULL`,
#'   all features in the VISTA object (optionally filtered by `feature_type`)
#'   are used.
#' @param col_genetype Column name in `rowData(x)` used to filter the background
#'   by gene type (default `"GENETYPE"`).
#' @param feature_type Gene type to retain in the background when `background`
#'   is `NULL` (default `"protein-coding"`).
#' @param ... Additional arguments to pass to `clusterProfiler::enricher()`.
#'
#' @return A list containing a single `enrichResult` object named `"enrich"`.
#' @aliases enrichMsigDB enrichMsigDB,VISTA-method
#' @export
setGeneric("enrichMsigDB", function(x,
                                    gene_list,
                                    from_type = "SYMBOL",
                                    orgdb,
                                    msigdb_category = "H",
                                    msigdb_subcategory = NULL,
                                    species = "Mus musculus",
                                    background = NULL,
                                    col_genetype = "GENETYPE",
                                    feature_type = "protein-coding",
                                    ...) {
  standardGeneric("enrichMsigDB")
}, signature = "x")

setMethod(
  "enrichMsigDB", "VISTA",
  function(x,
           gene_list,
           from_type = "SYMBOL",
           orgdb,
                                   msigdb_category  = "H",
                                   msigdb_subcategory = NULL,
                                   species   = "Mus musculus",
                                   background = NULL,
                                   col_genetype = "GENETYPE",
                                   feature_type = "protein-coding",
                                   ...) {

    if (!inherits(x, "VISTA")) stop("`x` must be a VISTA object.")
    if (missing(gene_list) || is.null(gene_list) || length(gene_list) == 0) {
      stop("`gene_list` is required and must not be NULL or empty.")
    }

    ## Choose a default OrgDb based on species if not supplied
    if (missing(orgdb) || is.null(orgdb)) {
      if (species %in% c("Mus musculus", "mouse", "mus musculus")) {
        orgdb <- tryCatch(getNamespace("org.Mm.eg.db")$org.Mm.eg.db,
                          error = function(e) NULL)
      } else if (species %in% c("Homo sapiens", "human", "homo sapiens")) {
        orgdb <- tryCatch(getNamespace("org.Hs.eg.db")$org.Hs.eg.db,
                          error = function(e) NULL)
      } else {
        stop("`orgdb` must be provided for species ", species)
      }
    }

    ## Retrieve MSigDB gene sets and pick the appropriate gene-field
    msig_df <- .vista_msigdbr(
      species = species,
      msigdb_category = msigdb_category,
      msigdb_subcategory = msigdb_subcategory
    )
    gene_field <- switch(toupper(from_type),
                         "SYMBOL"   = "gene_symbol",
                         "ENSEMBL"  = "ensembl_gene",
                         "ENTREZID" = "ncbi_gene",
                         "gene_symbol")
    msig_t2g <- unique(msig_df[c("gs_name", gene_field)])
    colnames(msig_t2g) <- c("term", "gene")

    ## Parse colon-separated IDs.  For SYMBOL, take the second element; otherwise, take the first.
    parse_ids <- function(ids, type) {
      if (!any(grepl(":", ids, fixed = TRUE))) return(ids)
      split <- strsplit(ids, ":", fixed = TRUE)
      if (toupper(type) == "SYMBOL") {
        vapply(split, function(p) if (length(p) > 1) p[2] else p[1], character(1), USE.NAMES = FALSE)
      } else {
        vapply(split, `[[`, character(1), 1L, USE.NAMES = FALSE)
      }
    }

    ## Convert IDs to match the msig_t2g$gene type
    convert_ids <- function(ids, from_type, to_type) {
      ids_base <- parse_ids(ids, from_type)
      # drop Ensembl version suffixes (e.g., ENSMUSG... .1)
      if (toupper(from_type) == "ENSEMBL") {
        ids_base <- sub("\\..*$", "", ids_base)
      }
      if (toupper(from_type) == "SYMBOL" && to_type == "gene_symbol") {
        return(unique(ids_base))
      }
      if (toupper(from_type) == "ENSEMBL" && to_type == "ensembl_gene") {
        return(unique(ids_base))
      }
      if (toupper(from_type) == "ENTREZID" && to_type == "ncbi_gene") {
        return(unique(ids_base))
      }
      ## Otherwise, convert via orgdb
      keytype  <- toupper(from_type)
      target   <- if (to_type == "gene_symbol") "SYMBOL"
      else if (to_type == "ensembl_gene") "ENSEMBL"
      else if (to_type == "ncbi_gene") "ENTREZID"
      else stop("Unsupported target type: ", to_type)
      res <- AnnotationDbi::select(orgdb, keys = ids_base,
                                   columns = target, keytype = keytype)
      return(unique(na.omit(res[[target]])))
    }

    ## Build the universe (background)
    if (is.null(background)) {
      ## Use all rownames as background
      if ("norm_counts" %in% assayNames(x)) {
        bg_raw <- rownames(assay(x, "norm_counts"))
      } else {
        anms   <- assayNames(x)
        bg_raw <- rownames(assay(x, anms[[1]]))
      }
      ## Optionally filter by gene type
      rd <- rowData(x)
      if (!is.null(rd) && col_genetype %in% colnames(rd)) {
        keep <- which(rd[[col_genetype]] == feature_type)
        if (length(keep) > 0) bg_raw <- bg_raw[keep]
      }
    } else {
      bg_raw <- unique(background)
    }
    bg_conv <- convert_ids(bg_raw, from_type, gene_field)
    bg_conv <- unique(bg_conv[!is.na(bg_conv) & nzchar(bg_conv)])
    # background handed to enricher; extend TERM2GENE with UN_ANNOTATED bucket so no genes are dropped
    term_genes <- unique(msig_t2g$gene)
    universe <- union(bg_conv, term_genes)
    msig_t2g_full <- rbind(
      msig_t2g,
      data.frame(term = "UN_ANNOTATED", gene = bg_conv, stringsAsFactors = FALSE)
    )

    ## Map the provided gene list to the appropriate ID type
    gl_ids <- convert_ids(gene_list, from_type, gene_field)
    if (length(gl_ids) == 0) {
      stop("`gene_list` could not be mapped to the requested identifier type.")
    }

    enr <- clusterProfiler::enricher(gene    = gl_ids,
                                     TERM2GENE = msig_t2g_full,
                                     universe  = NULL,
                                     ...)
    return(list(enrich = enr))
  }
)

#' Plot enrichment results using -log10(FDR)
#'
#' Generates a dot plot of enrichment results for `enrichResult`, `gseaResult`,
#' or `compareClusterResult` objects (including those returned by
#' [enrichMsigDB()]). Points are sized by gene/set count and coloured by
#' -log10(FDR). For compareCluster results, the plot is faceted by cluster with
#' top terms selected per cluster.
#'
#' @param x An object of class `enrichResult`, `gseaResult`, or `compareClusterResult`.
#' @param top_n Integer; number of top terms to plot (per cluster for compareCluster).
#' @param title Optional plot title.
#'
#' @return A `ggplot2` object.
#' @export
get_enrichment_plot <- function(x, top_n = 10, title = NULL) {
  if (inherits(x, "compareClusterResult")) {
    df <- x@compareClusterResult
    type <- "compare"
  } else if (inherits(x, "gseaResult")) {
    df <- x@result
    type <- "gsea"
  } else if (inherits(x, "enrichResult")) {
    df <- x@result
    type <- "enrich"
  } else {
    stop("Input must be enrichResult, gseaResult, or compareClusterResult.")
  }

  df <- df |>
    dplyr::mutate(logFDR = -log10(.data$p.adjust))

  if (type == "compare") {
    df <- df |>
      dplyr::group_by(.data$Cluster) |>
      dplyr::arrange(dplyr::desc(.data$logFDR), .by_group = TRUE) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::ungroup() |>
      dplyr::mutate(Description = forcats::fct_reorder(.data$Description, .data$logFDR))

    p <- ggplot2::ggplot(df, ggplot2::aes(x = logFDR, y = Description)) +
      ggplot2::geom_point(ggplot2::aes(size = .data$Count, color = .data$logFDR)) +
      ggplot2::facet_wrap(~Cluster, scales = "free_y") +
      ggplot2::scale_color_gradient(low = "blue", high = "red", name = expression(-log[10](FDR))) +
      ggplot2::scale_size(name = "Gene Count") +
      ggplot2::labs(
        x = expression(-log[10](FDR)),
        y = NULL,
        title = title %||% "Top Enriched Terms by Group"
      ) +
      ggplot2::theme_minimal(base_size = 14)

  } else {
    size_col <- if (type == "gsea") "setSize" else "Count"

    df <- df |>
      dplyr::arrange(dplyr::desc(.data$logFDR)) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::mutate(
        Description = forcats::fct_reorder(.data$Description, .data$logFDR),
        .plot_size = .data[[size_col]]
      )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = logFDR, y = Description)) +
      ggplot2::geom_point(ggplot2::aes(size = .data$.plot_size, color = .data$logFDR)) +
      ggplot2::scale_color_gradient(low = "blue", high = "red", name = expression(-log[10](FDR))) +
      ggplot2::scale_size(name = "Gene Count") +
      ggplot2::labs(
        x = expression(-log[10](FDR)),
        y = NULL,
        title = title %||% "Top Enriched Terms"
      ) +
      ggplot2::theme_minimal(base_size = 14)
  }

  p
}

#' Extract genes from enriched pathways
#'
#' Parses pathway-level gene members from an enrichment result and returns
#' either a long table, pathway-indexed list, or a unique vector of genes.
#'
#' @param x An `enrichResult`/`gseaResult`, or a list containing element
#'   `enrich` (e.g. output from `get_msigdb_enrichment()`).
#' @param pathways Optional character vector of pathway names to keep. Matches
#'   against `pathway_column`.
#' @param top_n Number of top pathways to use when `pathways` is `NULL`.
#'   Ranking uses `p.adjust` (then `pvalue`) when available. Default: `10`.
#' @param pathway_column Which enrichment column to match pathway names against:
#'   `"Description"` (default) or `"ID"`.
#' @param gene_column Which column stores pathway members. `"auto"` (default)
#'   uses `"geneID"` when present, otherwise `"core_enrichment"`.
#' @param gene_sep Delimiter used in pathway gene strings (default `"/"`).
#' @param return_type One of `"long"`, `"list"`, or `"vector"`.
#'
#' @return
#' Depending on `return_type`:
#' \itemize{
#'   \item `"long"`: data frame with `pathway_id`, `pathway`, and `gene`.
#'   \item `"list"`: named list of character vectors (genes per pathway).
#'   \item `"vector"`: unique character vector of genes.
#' }
#'
#' @examples
#' if (requireNamespace("msigdbr", quietly = TRUE)) {
#'   vista <- example_vista()
#'   msig <- get_msigdb_enrichment(
#'     vista,
#'     sample_comparison = names(comparisons(vista))[1],
#'     regulation = "Both",
#'     msigdb_category = "H",
#'     from_type = "ENSEMBL"
#'   )
#'   pathway_tbl <- get_pathway_genes(msig, top_n = 5, return_type = "long")
#'   head(pathway_tbl)
#' }
#'
#' \donttest{
#' data("count_data", package = "VISTA")
#' data("sample_metadata", package = "VISTA")
#'
#' vista <- create_vista(
#'   counts = count_data,
#'   sample_info = sample_metadata,
#'   column_geneid = "gene_id",
#'   group_column = "cond_long",
#'   group_numerator = "treatment1",
#'   group_denominator = "control"
#' )
#'
#' msig <- get_msigdb_enrichment(
#'   vista,
#'   sample_comparison = names(comparisons(vista))[1],
#'   regulation = "Up",
#'   species = "Homo sapiens",
#'   from_type = "ENSEMBL"
#' )
#'
#' pathway_tbl <- get_pathway_genes(msig, top_n = 5, return_type = "long")
#' head(pathway_tbl)
#' }
#'
#' @export
get_pathway_genes <- function(x,
                              pathways = NULL,
                              top_n = 10,
                              pathway_column = c("Description", "ID"),
                              gene_column = c("auto", "geneID", "core_enrichment"),
                              gene_sep = "/",
                              return_type = c("long", "list", "vector")) {
  pathway_column <- match.arg(pathway_column)
  gene_column <- match.arg(gene_column)
  # Not the plot-vs-data vocabulary: this selects the SHAPE of a purely tabular
  # result, so it is deliberately left alone.
  return_type <- match.arg(return_type)

  if (is.list(x) && !inherits(x, c("enrichResult", "gseaResult")) && "enrich" %in% names(x)) {
    x <- x$enrich
  }
  if (!inherits(x, c("enrichResult", "gseaResult"))) {
    cli::cli_abort(
      "{.arg x} must be an {.cls enrichResult}/{.cls gseaResult}, or a list containing {.val enrich}."
    )
  }

  res <- as.data.frame(x@result, stringsAsFactors = FALSE)
  if (!nrow(res)) {
    cli::cli_abort("No enrichment rows available in {.arg x}.")
  }
  if (!pathway_column %in% colnames(res)) {
    cli::cli_abort(
      "Column {.val {pathway_column}} not found in enrichment results. Available: {.val {colnames(res)}}"
    )
  }

  if (identical(gene_column, "auto")) {
    gene_column <- if ("geneID" %in% colnames(res)) {
      "geneID"
    } else if ("core_enrichment" %in% colnames(res)) {
      "core_enrichment"
    } else {
      NA_character_
    }
  }
  if (is.na(gene_column) || !gene_column %in% colnames(res)) {
    cli::cli_abort(
      "Could not determine a pathway gene column. Expected {.val geneID} or {.val core_enrichment}."
    )
  }

  if (!is.null(pathways)) {
    pathways <- unique(as.character(pathways))
    available <- unique(as.character(res[[pathway_column]]))
    missing_paths <- setdiff(pathways, available)
    if (length(missing_paths) > 0) {
      cli::cli_abort(
        "Some {.arg pathways} were not found: {.val {missing_paths}}"
      )
    }
    path_levels <- pathways
    res <- res[res[[pathway_column]] %in% pathways, , drop = FALSE]
    res <- res[order(match(res[[pathway_column]], path_levels)), , drop = FALSE]
  } else {
    if (!is.null(top_n)) {
      if (!is.numeric(top_n) || length(top_n) != 1 || is.na(top_n) || top_n <= 0) {
        cli::cli_abort("{.arg top_n} must be a single positive number.")
      }
      if ("p.adjust" %in% colnames(res)) {
        if ("pvalue" %in% colnames(res)) {
          res <- res[order(res$p.adjust, res$pvalue, na.last = TRUE), , drop = FALSE]
        } else {
          res <- res[order(res$p.adjust, na.last = TRUE), , drop = FALSE]
        }
      } else if ("pvalue" %in% colnames(res)) {
        res <- res[order(res$pvalue, na.last = TRUE), , drop = FALSE]
      }
      res <- utils::head(res, top_n)
    }
    path_levels <- unique(as.character(res[[pathway_column]]))
  }

  pathway_values <- as.character(res[[pathway_column]])
  pathway_ids <- if ("ID" %in% colnames(res)) as.character(res$ID) else pathway_values

  split_genes <- lapply(as.character(res[[gene_column]]), function(s) {
    if (length(s) == 0 || is.na(s) || !nzchar(s)) return(character(0))
    genes <- trimws(strsplit(s, split = gene_sep, fixed = TRUE)[[1]])
    genes <- genes[nzchar(genes)]
    unique(genes)
  })

  long_parts <- Map(
    function(path_id, path_label, genes) {
      if (!length(genes)) return(NULL)
      data.frame(
        pathway_id = path_id,
        pathway = path_label,
        gene = genes,
        stringsAsFactors = FALSE
      )
    },
    pathway_ids,
    pathway_values,
    split_genes
  )
  long_df <- dplyr::bind_rows(long_parts)

  if (!nrow(long_df)) {
    cli::cli_abort("No genes could be parsed from selected pathways.")
  }

  long_df <- long_df |>
    dplyr::distinct(.data$pathway_id, .data$pathway, .data$gene)

  # Keep pathway order stable for list/vector outputs.
  keep_levels <- unique(path_levels[path_levels %in% unique(long_df$pathway)])
  long_df$pathway <- factor(long_df$pathway, levels = keep_levels)
  long_df <- long_df[order(long_df$pathway), , drop = FALSE]
  long_df$pathway <- as.character(long_df$pathway)

  if (return_type == "long") {
    return(long_df)
  }

  pathway_list <- split(long_df$gene, factor(long_df$pathway, levels = keep_levels))
  pathway_list <- lapply(pathway_list, unique)

  if (return_type == "list") {
    return(pathway_list)
  }

  unique(unlist(pathway_list, use.names = FALSE))
}

#' Plot pathway-specific expression heatmaps from enrichment output
#'
#' This wrapper bridges enrichment results and expression heatmaps. It extracts
#' genes from selected pathways (via [get_pathway_genes()]), maps them to the
#' `VISTA` feature IDs, and forwards to [get_expression_heatmap()].
#'
#' @param x A `VISTA` object.
#' @param enrichment An `enrichResult`/`gseaResult`, or a list with element
#'   `enrich` as returned by `get_*_enrichment()` helpers.
#' @param sample_group Character vector of group labels to include (same semantics as
#'   [get_expression_heatmap()]).
#' @param pathways Optional pathway names to include. When `NULL`, top pathways
#'   are selected using `top_n`.
#' @param top_n Number of top pathways used when `pathways = NULL`. Default: `5`.
#' @param pathway_column Pathway matching column, `"Description"` (default) or `"ID"`.
#' @param gene_column Pathway gene-member column. `"auto"` uses `"geneID"` or
#'   `"core_enrichment"` based on availability.
#' @param gene_sep Delimiter used to parse pathway gene strings (default `"/"`).
#' @param gene_mode How to combine pathway genes for plotting: `"union"` (default)
#'   or `"intersection"`.
#' @param gene_id_column Optional column in `rowData(x)` used to map
#'   enrichment genes back to VISTA rownames (e.g., `"SYMBOL"` or `"ENTREZID"`).
#'   Leave `NULL` when enrichment genes already match VISTA rownames.
#' @param max_genes Optional cap on the number of genes passed to the heatmap.
#' @param return_type One of `"plot"` (default), `"data"`, or `"both"`. The
#'   legacy values `"heatmap"` and `"genes"` are still accepted and warn.
#' @param ... Additional arguments passed to [get_expression_heatmap()].
#'
#' @return
#' Depending on `return_type`:
#' \itemize{
#'   \item `"heatmap"`: a `ComplexHeatmap` object from [get_expression_heatmap()].
#'   \item `"both"`: list with `heatmap`, `genes`, and `pathway_genes`.
#'   \item `"genes"`: character vector of mapped genes selected for plotting.
#' }
#'
#' @examples
#' if (requireNamespace("msigdbr", quietly = TRUE)) {
#'   vista <- example_vista()
#'   msig <- get_msigdb_enrichment(
#'     vista,
#'     sample_comparison = names(comparisons(vista))[1],
#'     regulation = "Both",
#'     msigdb_category = "H",
#'     from_type = "ENSEMBL"
#'   )
#'   genes <- get_pathway_heatmap(
#'     vista,
#'     enrichment = msig,
#'     sample_group = c("control", "treatment1"),
#'     top_n = 3,
#'     return_type = "data"
#'   )
#'   head(genes)
#' }
#'
#' \donttest{
#' data("count_data", package = "VISTA")
#' data("sample_metadata", package = "VISTA")
#'
#' vista <- create_vista(
#'   counts = count_data,
#'   sample_info = sample_metadata,
#'   column_geneid = "gene_id",
#'   group_column = "cond_long",
#'   group_numerator = "treatment1",
#'   group_denominator = "control"
#' )
#'
#' msig <- get_msigdb_enrichment(
#'   vista,
#'   sample_comparison = names(comparisons(vista))[1],
#'   regulation = "Up",
#'   species = "Homo sapiens",
#'   from_type = "ENSEMBL"
#' )
#'
#' get_pathway_heatmap(
#'   vista,
#'   enrichment = msig,
#'   sample_group = c("control", "treatment1"),
#'   top_n = 3,
#'   value_transform = "zscore",
#'   annotate_columns = TRUE,
#'   summarise_replicates = FALSE
#' )
#' }
#'
#' @export
get_pathway_heatmap <- function(x,
                                enrichment,
                                sample_group = NULL,
                                pathways = NULL,
                                top_n = 5,
                                pathway_column = c("Description", "ID"),
                                gene_column = c("auto", "geneID", "core_enrichment"),
                                gene_sep = "/",
                                gene_mode = c("union", "intersection"),
                                gene_id_column = NULL,
                                max_genes = NULL,
                                return_type = c("plot", "data", "both"),
                                ...) {
  stopifnot(inherits(x, "VISTA"))
  gene_mode <- match.arg(gene_mode)
  return_type <- .vista_resolve_return_type(
    return_type, fun = "get_pathway_heatmap",
    legacy = c(heatmap = "plot", genes = "data")
  )

  dots <- list(...)
  blocked <- intersect(c("genes", "sample_group"), names(dots))
  if (length(blocked) > 0) {
    cli::cli_abort(
      "Argument(s) {.val {blocked}} are managed by {.fun get_pathway_heatmap}. Use {.arg pathways}/{.arg top_n} and {.arg sample_group}."
    )
  }

  pathway_genes <- get_pathway_genes(
    x = enrichment,
    pathways = pathways,
    top_n = top_n,
    pathway_column = pathway_column,
    gene_column = gene_column,
    gene_sep = gene_sep,
    return_type = "list"
  )

  selected_genes <- switch(
    gene_mode,
    union = unique(unlist(pathway_genes, use.names = FALSE)),
    intersection = {
      if (length(pathway_genes) == 1L) pathway_genes[[1]] else Reduce(intersect, pathway_genes)
    }
  )

  if (!length(selected_genes)) {
    cli::cli_abort(
      "No genes remained after combining selected pathways with {.arg gene_mode} = {.val {gene_mode}}."
    )
  }

  mapped_genes <- selected_genes
  if (!is.null(gene_id_column)) {
    rd <- as.data.frame(SummarizedExperiment::rowData(x), stringsAsFactors = FALSE)
    if (!gene_id_column %in% colnames(rd)) {
      cli::cli_abort(
        "{.arg gene_id_column} '{gene_id_column}' not found in rowData(x)."
      )
    }
    lookup <- as.character(rd[[gene_id_column]])
    names(lookup) <- rownames(x)
    mapped_genes <- names(lookup)[match(selected_genes, lookup)]
    mapped_genes <- unique(mapped_genes[!is.na(mapped_genes)])
  }

  genes_in_object <- intersect(mapped_genes, rownames(x))
  if (!length(genes_in_object)) {
    cli::cli_abort(
      "No selected pathway genes matched rownames(x). If IDs differ, supply {.arg gene_id_column}."
    )
  }

  if (!is.null(max_genes)) {
    if (!is.numeric(max_genes) || length(max_genes) != 1 || is.na(max_genes) || max_genes <= 0) {
      cli::cli_abort("{.arg max_genes} must be a single positive number.")
    }
    if (length(genes_in_object) > max_genes) {
      # Pathway membership is a SET. clusterProfiler's gene order is arbitrary
      # and is not stable across R sessions, so taking the positional head would
      # make the plotted genes differ between runs of the same analysis. Rank by
      # expression variance -- the most informative rows for a heatmap -- and
      # break ties on the identifier so the selection is fully deterministic.
      mat <- SummarizedExperiment::assay(x, "norm_counts")[genes_in_object, , drop = FALSE]
      score <- matrixStats::rowVars(mat)
      score[!is.finite(score)] <- -Inf
      ord <- order(-score, genes_in_object)
      genes_in_object <- genes_in_object[utils::head(ord, max_genes)]
    }
  }

  if (return_type == "data") {
    return(genes_in_object)
  }

  heatmap_obj <- do.call(
    get_expression_heatmap,
    c(
      list(
        x = x,
        sample_group = sample_group,
        genes = genes_in_object
      ),
      dots
    )
  )

  if (return_type == "plot") {
    return(heatmap_obj)
  }

  list(
    heatmap = heatmap_obj,
    genes = genes_in_object,
    pathway_genes = pathway_genes
  )
}


#' Run MSigDB enrichment directly from a VISTA comparison
#'
#' Convenience wrapper that pulls regulated genes from a stored differential
#' expression comparison in a `VISTA` object and runs [enrichMsigDB()] on them.
#'
#' @param x A `VISTA` object with DE results in `metadata(x)$de_results`.
#' @param sample_comparison Character scalar naming the comparison to use.
#' @param regulation One of `"Up"`, `"Down"`, `"Both"`, or `"All"`; selects which
#'   genes to send to enrichment.
#' @param from_type Identifier type of the genes in the DE table (passed to
#'   [enrichMsigDB()], default `"SYMBOL"`). Ensembl versions are stripped automatically.
#' @param orgdb An `OrgDb` object used for ID conversion (passed through). If omitted,
#'   the default is chosen from `species` (mouse/human).
#' @param msigdb_category MSigDB category (e.g., `"H"`, `"C2"`, `"C5"`). Default `"H"`.
#' @param msigdb_subcategory Optional MSigDB sub-collection. Default `NULL`.
#' @param species Species name for MSigDB (default `"Mus musculus"`).
#' @param background Optional background gene set (passed to [enrichMsigDB()]). Default `NULL`
#'   uses all features in `x` (optionally filtered by `feature_type`).
#' @param col_genetype Column in `rowData(x)` used to filter background by gene type. Default `"GENETYPE"`.
#' @param feature_type Gene type to retain in the background when filtering. Default `"protein-coding"`.
#' @param ... Additional arguments forwarded to [enrichMsigDB()].
#'
#' @return A list with `enrich` containing an `enrichResult`.
#'
#' @examples
#' if (requireNamespace("msigdbr", quietly = TRUE)) {
#'   vista <- example_vista()
#'   comp <- names(comparisons(vista))[1]
#'   msig <- get_msigdb_enrichment(
#'     vista,
#'     sample_comparison = comp,
#'     regulation = "Both",
#'     msigdb_category = "H",
#'     from_type = "ENSEMBL"
#'   )
#'   class(msig$enrich)
#' }
#'
#' \donttest{
#' # Create VISTA object
#' data("count_data", package = "VISTA")
#' data("sample_metadata", package = "VISTA")
#'
#' vista <- create_vista(
#'   counts = count_data[seq_len(200), ],
#'   sample_info = sample_metadata[seq_len(6), ],
#'   column_geneid = "gene_id",
#'   group_column = "cond_long",
#'   group_numerator = "treatment1",
#'   group_denominator = "control"
#' )
#' comp <- names(comparisons(vista))[1]
#'
#' # Run MSigDB enrichment on upregulated genes
#' msig_up <- get_msigdb_enrichment(
#'   vista,
#'   sample_comparison = comp,
#'   regulation = "Up",
#'   msigdb_category = "H",  # Hallmark gene sets
#'   from_type = "ENSEMBL"
#' )
#'
#' if (!is.null(msig_up$enrich)) {
#'   # View results
#'   head(msig_up$enrich)
#'   # Visualize enrichment
#'   get_enrichment_plot(msig_up$enrich)
#' }
#'
#' # Enrichment for downregulated genes
#' msig_down <- get_msigdb_enrichment(
#'   vista,
#'   sample_comparison = comp,
#'   regulation = "Down",
#'   msigdb_category = "C2"  # Curated gene sets
#' )
#' }
#'
#' @export
get_msigdb_enrichment <- function(x,
                                  sample_comparison,
                                  regulation = c("Up", "Down", "Both", "All"),
                                  from_type = "SYMBOL",
                                  orgdb,
                                  msigdb_category  = "H",
                                  msigdb_subcategory = NULL,
                                  species   = "Mus musculus",
                                  background = NULL,
                                  col_genetype = "GENETYPE",
                                  feature_type = "protein-coding",
                                  ...) {
  regulation <- match.arg(regulation)
  genes_list <- get_genes_by_regulation(x, sample_comparisons = sample_comparison, regulation = regulation)
  genes <- unique(unlist(genes_list, use.names = FALSE))
  if (!length(genes)) {
    cli::cli_abort("No genes found for {.arg sample_comparison} with regulation {.val {regulation}}.")
  }

  enrichMsigDB(
    x = x,
    gene_list = genes,
    from_type = from_type,
    orgdb = orgdb,
    msigdb_category = msigdb_category,
    msigdb_subcategory = msigdb_subcategory,
    species = species,
    background = background,
    col_genetype = col_genetype,
    feature_type = feature_type,
    ...
  )
}

# ──────────────────────────────────────────────────────────────────────────────
# GO / KEGG / GSEA helpers
# ──────────────────────────────────────────────────────────────────────────────

.vista_default_orgdb <- function(species) {
  if (missing(species) || is.null(species)) species <- "Mus musculus"
  if (species %in% c("Mus musculus", "mouse", "mus musculus")) {
    tryCatch(getNamespace("org.Mm.eg.db")$org.Mm.eg.db, error = function(e) NULL)
  } else if (species %in% c("Homo sapiens", "human", "homo sapiens")) {
    tryCatch(getNamespace("org.Hs.eg.db")$org.Hs.eg.db, error = function(e) NULL)
  } else {
    NULL
  }
}

#' Map identifiers preserving length and position
#'
#' Unlike `.vista_convert_ids()`, which de-duplicates and drops unmapped
#' identifiers because its callers want a gene *set*, this returns a vector the
#' same length as `ids` with `NA` where no mapping exists. That property is
#' essential anywhere identifiers are used as keys against a parallel vector of
#' values.
#'
#' It also differs from `.map_gene_ids()`, which falls back to the *input*
#' identifier when a mapping is missing. That is right for display labels and
#' wrong for keys -- handing an unmapped SYMBOL to a function expecting ENTREZ
#' IDs trades a visible gap for a silent mismatch.
#'
#' @param ids Character vector of identifiers.
#' @param from_type,to_type Key types understood by the supplied `orgdb`.
#' @param orgdb An `OrgDb` object.
#' @param id_map Optional named character vector used instead of `orgdb`,
#'   primarily so the behaviour can be unit-tested without an annotation package.
#'
#' @return A character vector the same length as `ids`, with `NA` for unmapped
#'   entries.
#' @keywords internal
#' @noRd
.vista_map_ids_strict <- function(ids, from_type, to_type, orgdb = NULL, id_map = NULL) {
  ids <- as.character(ids)
  if (identical(toupper(from_type), toupper(to_type))) {
    return(ids)
  }
  if (toupper(from_type) == "ENSEMBL") ids <- sub("\\..*$", "", ids)

  if (!is.null(id_map)) {
    return(unname(id_map[ids]))
  }
  if (is.null(orgdb)) {
    cli::cli_abort("An {.arg orgdb} is required to map {.val {from_type}} to {.val {to_type}}.")
  }

  keytype <- toupper(from_type)
  target <- toupper(to_type)
  res <- tryCatch(
    suppressMessages(
      AnnotationDbi::select(orgdb, keys = unique(ids), keytype = keytype, columns = target)
    ),
    error = function(e) NULL
  )
  if (is.null(res) || !nrow(res)) {
    return(rep(NA_character_, length(ids)))
  }
  res <- res[!duplicated(res[[keytype]]), , drop = FALSE]
  map <- stats::setNames(as.character(res[[target]]), res[[keytype]])
  out <- unname(map[ids])
  out[!nzchar(out) | is.na(out)] <- NA_character_
  out
}

#' Build a GSEA ranked vector with names that stay attached to their scores
#'
#' Maps `ids`, drops entries that fail to map (removing name and score
#' together), resolves many-to-one collisions by keeping the largest absolute
#' score, and returns a decreasing-sorted named numeric vector.
#'
#' @param scores Numeric vector of ranking statistics (typically log2FC).
#' @param ids Character vector of identifiers, parallel to `scores`.
#' @inheritParams .vista_map_ids_strict
#' @param verbose Logical; report dropped and collapsed counts.
#'
#' @return A named numeric vector sorted in decreasing order.
#' @keywords internal
#' @noRd
.vista_build_rank_vector <- function(scores, ids, from_type, to_type,
                                     orgdb = NULL, id_map = NULL, verbose = TRUE) {
  stopifnot(length(scores) == length(ids))

  keep <- is.finite(scores)
  scores <- scores[keep]
  ids <- as.character(ids)[keep]

  mapped <- .vista_map_ids_strict(
    ids, from_type = from_type, to_type = to_type, orgdb = orgdb, id_map = id_map
  )

  # Drop name and score together. Assigning a shorter vector to names() would
  # pad with NA and slide every remaining score onto the wrong identifier.
  ok <- !is.na(mapped) & nzchar(mapped)
  n_dropped <- sum(!ok)
  mapped <- mapped[ok]
  scores <- scores[ok]

  if (!length(mapped)) {
    cli::cli_abort(
      "No gene identifiers could be mapped from {.val {from_type}} to {.val {to_type}}."
    )
  }

  # Many-to-one mappings: keep the most extreme score for each target ID.
  n_collapsed <- 0L
  if (anyDuplicated(mapped)) {
    ord <- order(abs(scores), decreasing = TRUE)
    mapped_ord <- mapped[ord]
    scores_ord <- scores[ord]
    first <- !duplicated(mapped_ord)
    n_collapsed <- sum(!first)
    mapped <- mapped_ord[first]
    scores <- scores_ord[first]
  }

  if (isTRUE(verbose) && (n_dropped > 0L || n_collapsed > 0L)) {
    cli::cli_inform(c(
      "Prepared {length(scores)} ranked gene{?s} for GSEA.",
      if (n_dropped > 0L) "i" = "Dropped {n_dropped} gene{?s} with no {.val {to_type}} mapping.",
      if (n_collapsed > 0L) "i" = "Collapsed {n_collapsed} duplicate {.val {to_type}} mapping{?s}, keeping the largest absolute score."
    ))
  }

  sort(stats::setNames(scores, mapped), decreasing = TRUE)
}

.vista_convert_ids <- function(ids, from_type, to_type, orgdb) {
  ids <- as.character(ids)
  if (identical(from_type, to_type) || is.null(to_type) || to_type %in% c("", "none")) return(unique(ids))
  if (toupper(from_type) == "ENSEMBL") ids <- sub("\\..*$", "", ids)
  keytype <- toupper(from_type)
  target  <- toupper(to_type)
  res <- AnnotationDbi::select(orgdb, keys = ids, keytype = keytype, columns = target)
  if (is.null(res) || !nrow(res)) return(unique(ids))
  res <- res[!duplicated(res[[keytype]]), , drop = FALSE]
  res <- res[!is.na(res[[target]]), , drop = FALSE]
  if (!nrow(res)) return(unique(ids))
  map <- stats::setNames(res[[target]], res[[keytype]])
  out <- map[ids]
  unique(out[!is.na(out) & nzchar(out)])
}

#' Run GO enrichment directly from a VISTA comparison
#'
#' @param x A `VISTA` object with DE results.
#' @param sample_comparison Comparison name to use.
#' @param regulation One of `"Up"`, `"Down"`, `"Both"`, or `"All"`; selects genes.
#' @param ont GO ontology: `"BP"`, `"MF"`, or `"CC"`.
#' @param from_type Identifier type in the DE tables (default `"SYMBOL"`).
#' @param orgdb OrgDb object; defaults to mouse/human based on `species`.
#' @param species Species name to infer default OrgDb.
#' @param background Optional background gene set; default uses all features.
#' @param ... Passed to `clusterProfiler::enrichGO()`.
#'
#' @return A list with `enrich` containing an `enrichResult`.
#' @export
get_go_enrichment <- function(x,
                              sample_comparison,
                              regulation = c("Up", "Down", "Both", "All"),
                              ont = c("BP", "MF", "CC"),
                              from_type = "SYMBOL",
                              orgdb = NULL,
                              species = "Mus musculus",
                              background = NULL,
                              ...) {
  regulation <- match.arg(regulation)
  ont <- match.arg(ont)
  orgdb <- orgdb %||% .vista_default_orgdb(species)
  if (is.null(orgdb)) cli::cli_abort("An {.arg orgdb} is required for GO enrichment.")

  genes_list <- get_genes_by_regulation(x, sample_comparisons = sample_comparison, regulation = regulation)
  genes <- unique(unlist(genes_list, use.names = FALSE))
  if (!length(genes)) cli::cli_abort("No genes found for {.arg sample_comparison} with regulation {.val {regulation}}.")

  gene_ids <- .vista_convert_ids(genes, from_type = from_type, to_type = from_type, orgdb = orgdb)
  bg_ids <- NULL
  if (!is.null(background)) {
    bg_ids <- .vista_convert_ids(background, from_type = from_type, to_type = from_type, orgdb = orgdb)
  }

  .vista_check_dots(
    list(...), fun = "get_go_enrichment",
    allowed = names(formals(clusterProfiler::enrichGO)),
    blocked = c("gene", "OrgDb", "keyType", "ont", "universe")
  )
  res <- clusterProfiler::enrichGO(gene          = gene_ids,
                                   OrgDb         = orgdb,
                                   keyType       = toupper(from_type),
                                   ont           = ont,
                                   universe      = bg_ids,
                                   ...)
  list(enrich = res)
}

#' Run KEGG enrichment directly from a VISTA comparison
#'
#' @inheritParams get_go_enrichment
#' @param kegg_species KEGG organism code (e.g., `"mmu"` or `"hsa"`). If `NULL`,
#'   inferred from `species`.
#'
#' @export
get_kegg_enrichment <- function(x,
                                sample_comparison,
                                regulation = c("Up", "Down", "Both", "All"),
                                from_type = "SYMBOL",
                                orgdb = NULL,
                                species = "Mus musculus",
                                kegg_species = NULL,
                                background = NULL,
                                ...) {
  regulation <- match.arg(regulation)
  orgdb <- orgdb %||% .vista_default_orgdb(species)
  if (is.null(orgdb)) cli::cli_abort("An {.arg orgdb} is required for KEGG enrichment.")
  if (is.null(kegg_species)) {
    kegg_species <- if (species %in% c("Homo sapiens", "human", "homo sapiens")) "hsa" else "mmu"
  }

  genes_list <- get_genes_by_regulation(x, sample_comparisons = sample_comparison, regulation = regulation)
  genes <- unique(unlist(genes_list, use.names = FALSE))
  if (!length(genes)) cli::cli_abort("No genes found for {.arg sample_comparison} with regulation {.val {regulation}}.")

  gene_ids <- .vista_convert_ids(genes, from_type = from_type, to_type = "ENTREZID", orgdb = orgdb)
  bg_ids <- NULL
  if (!is.null(background)) {
    bg_ids <- .vista_convert_ids(background, from_type = from_type, to_type = "ENTREZID", orgdb = orgdb)
  }

  .vista_check_dots(
    list(...), fun = "get_kegg_enrichment",
    allowed = names(formals(clusterProfiler::enrichKEGG)),
    blocked = c("gene", "organism", "universe")
  )
  res <- clusterProfiler::enrichKEGG(gene = gene_ids,
                                     organism = kegg_species,
                                     universe = bg_ids,
                                     ...)
  list(enrich = res)
}

#' Gene set enrichment analysis (GSEA) from a VISTA comparison
#'
#' @inheritParams get_go_enrichment
#' @param set_type One of `"msigdb"`, `"go"`, or `"kegg"` selecting the gene set source.
#' @param msigdb_category,msigdb_subcategory Passed to `msigdbr::msigdbr()` when `set_type = "msigdb"`.
#' @param ... Additional arguments forwarded to the underlying GSEA function:
#'   `clusterProfiler::GSEA()` (msigdb TERM2GENE), `gseGO()`, or `gseKEGG()`
#'   depending on `set_type`.
#'
#' @export
get_gsea <- function(x,
                     sample_comparison,
                     set_type = c("msigdb", "go", "kegg"),
                     from_type = "SYMBOL",
                     orgdb = NULL,
                     species = "Mus musculus",
                     msigdb_category = "H",
                     msigdb_subcategory = NULL,
                     ...) {
  set_type <- match.arg(set_type)
  orgdb <- orgdb %||% .vista_default_orgdb(species)

  comps <- comparisons(x)
  if (!sample_comparison %in% names(comps)) cli::cli_abort("Comparison {.val {sample_comparison}} not found.")
  tbl <- comps[[sample_comparison]]
  if (!"log2fc" %in% names(tbl)) {
    if ("log2FoldChange" %in% names(tbl)) names(tbl)[names(tbl) == "log2FoldChange"] <- "log2fc"
    else if ("logFC" %in% names(tbl)) names(tbl)[names(tbl) == "logFC"] <- "log2fc"
  }
  fc <- tbl$log2fc
  if (is.null(fc)) cli::cli_abort("No log2FC column found for {.val {sample_comparison}}.")
  rn <- tbl$gene_id %||% rownames(tbl)
  if (is.null(rn)) cli::cli_abort("DE table needs gene IDs for GSEA.")

  fc <- suppressWarnings(as.numeric(fc))

  if (set_type == "msigdb") {
    msig_df <- .vista_msigdbr(
      species = species,
      msigdb_category = msigdb_category,
      msigdb_subcategory = msigdb_subcategory
    )
    gene_field <- switch(toupper(from_type),
                         "SYMBOL" = "gene_symbol",
                         "ENSEMBL" = "ensembl_gene",
                         "ENTREZID" = "ncbi_gene",
                         "gene_symbol")
    msig_t2g <- msig_df[, c("gs_name", gene_field)]
    colnames(msig_t2g) <- c("term", "gene")
    rank_vec <- .vista_build_rank_vector(
      scores = fc, ids = rn,
      from_type = from_type, to_type = toupper(from_type), orgdb = orgdb
    )
    res <- clusterProfiler::GSEA(geneList = rank_vec,
                                 TERM2GENE = msig_t2g,
                                 ...)
  } else if (set_type == "go") {
    keytype <- toupper(from_type)
    rank_vec <- .vista_build_rank_vector(
      scores = fc, ids = rn,
      from_type = from_type, to_type = keytype, orgdb = orgdb
    )
    res <- clusterProfiler::gseGO(geneList = rank_vec,
                                  OrgDb = orgdb,
                                  keyType = keytype,
                                  ...)
  } else {
    org_code <- if (species %in% c("Homo sapiens", "human", "homo sapiens")) "hsa" else "mmu"
    rank_vec <- .vista_build_rank_vector(
      scores = fc, ids = rn,
      from_type = from_type, to_type = "ENTREZID", orgdb = orgdb
    )
    res <- clusterProfiler::gseKEGG(geneList = rank_vec,
                                    organism = org_code,
                                    ...)
  }
  list(enrich = res)
}


# ──────────────────────────────────────────────────────────────────────────────
# Enrichment chord diagram
# ──────────────────────────────────────────────────────────────────────────────

#' Chord diagram of enrichment gene--pathway relationships
#'
#' Draws a chord diagram linking genes to the enriched pathways they belong to.
#' Chords can be coloured by fold-change, regulation direction, or source
#' pathway.
#'
#' The plot reveals **hub genes** (those driving multiple enriched terms) and
#' pathway redundancy (terms sharing many genes).
#' This complements [get_enrichment_plot()] (which shows significance) and
#' [get_pathway_heatmap()] (which shows expression patterns).
#'
#' @param x An `enrichResult`, `gseaResult`, or `compareClusterResult` from
#'   clusterProfiler, or a list containing an `enrich` element (e.g. output of
#'   [get_msigdb_enrichment()]).
#' @param vista Optional `VISTA` object.
#'   Required when `color_by` is `"foldchange"` or `"regulation"`.
#' @param sample_comparison Character scalar naming the DE comparison in `vista`
#'   to pull log2FC values from.
#'   Required when `vista` is supplied.
#' @param pathways Optional character vector of pathway names to include.
#'   Matches against `pathway_column`.
#' @param top_n Number of top pathways to display when `pathways` is `NULL`
#'   (default `8`).
#' @param pathway_column Column in the enrichment result to match pathway names:
#'   `"Description"` (default) or `"ID"`.
#' @param gene_column Column storing gene members: `"auto"` (default),
#'   `"geneID"`, or `"core_enrichment"`.
#' @param gene_sep Delimiter separating genes in `gene_column` (default `"/"`).
#' @param min_pathways Minimum number of pathways a gene must appear in to be
#'   shown. Set to `2` to display only hub genes shared across terms. Default
#'   `1` (show all genes).
#' @param max_genes Maximum number of genes to display (default `50`). A safety
#'   cap for readability.
#' @param gene_order_by Order of gene sectors in the chord plot:
#'   `"none"` (default), `"foldchange"` (descending log2FC),
#'   or `"abs_foldchange"` (descending absolute log2FC).
#'   Fold-change based ordering requires `vista` + `sample_comparison`.
#' @param gene_id_column Column in `rowData(vista)` used to map enrichment
#'   gene IDs to `vista` rownames (for FC lookup).
#' @param display_id Column in `rowData(vista)` providing display-friendly
#'   gene names.
#' @param color_by How to colour chords: `"foldchange"` (continuous gradient),
#'   `"regulation"` (Up / Down / Other), or `"pathway"` (source pathway).
#'   Falls back to `"pathway"` when `vista` is `NULL`.
#' @param up_color Colour for up-regulated genes (default `"#D73027"`).
#' @param down_color Colour for down-regulated genes (default `"#1A9850"`).
#' @param other_color Colour for non-significant genes (default `"grey70"`).
#' @param pathway_colors Optional named vector of colours for pathway sectors.
#'   When `NULL`, colours are generated from an HCL palette.
#' @param transparency Chord transparency, 0--1 (default `0.4`).
#' @param gap_degree Gap between sectors in degrees (default `2`).
#' @param label_cex Text size for sector labels (default `0.7`).
#' @param title Optional plot title.
#' @param return_type One of `"data"` (default), `"plot"` or `"both"`. The
#'   default is `"data"` because this function draws to the active device and
#'   has always returned its table invisibly. `"plot"` returns a recorded plot
#'   that [save_vista_plot()] can write to a file.
#'
#' @return With `return_type = "data"` (the default), invisibly a list with:
#'   \describe{
#'     \item{gene_data}{Tibble of genes with pathway membership and (optionally)
#'       fold-change values.}
#'     \item{hub_genes}{Character vector of genes appearing in two or more
#'       pathways.}
#'   }
#'   The chord diagram is drawn as a side effect.
#'
#' @examples
#' v <- example_vista()
#' msig <- get_msigdb_enrichment(
#'   v,
#'   sample_comparison = names(comparisons(v))[1],
#'   regulation = "Both",
#'   msigdb_category = "H",
#'   from_type = "ENSEMBL"
#' )
#' get_enrichment_chord(msig, top_n = 5)
#'
#' \donttest{
#' data("count_data", package = "VISTA")
#' data("sample_metadata", package = "VISTA")
#'
#' vista <- create_vista(
#'   counts = count_data[seq_len(200), ],
#'   sample_info = sample_metadata[seq_len(6), ],
#'   column_geneid = "gene_id",
#'   group_column = "cond_long",
#'   group_numerator = "treatment1",
#'   group_denominator = "control"
#' )
#'
#' msig <- get_msigdb_enrichment(
#'   vista,
#'   sample_comparison = names(comparisons(vista))[1],
#'   regulation = "Up", from_type = "ENSEMBL"
#' )
#'
#' # Simple: pathway-coloured chords
#' get_enrichment_chord(msig)
#'
#' # With fold-change colouring
#' get_enrichment_chord(
#'   msig, vista = vista,
#'   sample_comparison = names(comparisons(vista))[1],
#'   color_by = "foldchange"
#' )
#'
#' # Hub genes only (shared across 2+ pathways)
#' pw_long <- get_pathway_genes(msig, return_type = "long")
#' if (any(table(pw_long$gene) >= 2)) {
#'   get_enrichment_chord(msig, min_pathways = 2)
#' }
#' }
#'
#' @export
get_enrichment_chord <- function(x,
                                 vista           = NULL,
                                 sample_comparison = NULL,
                                 pathways        = NULL,
                                 top_n           = 8,
                                 pathway_column  = c("Description", "ID"),
                                 gene_column     = c("auto", "geneID", "core_enrichment"),
                                 gene_sep        = "/",
                                 min_pathways    = 1,
                                 max_genes       = 50,
                                 gene_order_by   = c("none", "foldchange", "abs_foldchange"),
                                 gene_id_column  = NULL,
                                 display_id      = NULL,
                                 color_by        = c("foldchange", "regulation", "pathway"),
                                 up_color        = "#D73027",
                                 down_color      = "#1A9850",
                                 other_color     = "grey70",
                                 pathway_colors  = NULL,
                                 transparency    = 0.4,
                                 gap_degree      = 2,
                                 label_cex       = 0.7,
                                 title           = NULL,
                                 return_type     = c("data", "plot", "both")) {

  # Default is "data" rather than "plot" because this function has always drawn
  # to the active device and returned its table invisibly; changing that would
  # break existing code.
  return_type <- .vista_resolve_return_type(
    return_type, fun = "get_enrichment_chord", default = "data"
  )
  want_plot <- return_type %in% c("plot", "both")
  if (want_plot) {
    # Off-screen devices (png, pdf) keep their display list off by default, so
    # recordPlot() would return an empty recording that silently saves nothing.
    # This has to happen before circlize draws.
    .vista_enable_display_list()
  }

  if (!requireNamespace("circlize", quietly = TRUE)) {
    cli::cli_abort(
      "Package {.pkg circlize} is required for chord diagrams. Install it with {.code install.packages(\"circlize\")}."
    )
  }

  color_by <- match.arg(color_by)
  gene_order_by <- match.arg(gene_order_by)

  # --- extract enrichment result from wrapper lists --------------------------
  if (is.list(x) && !inherits(x, c("enrichResult", "gseaResult", "compareClusterResult")) &&
      "enrich" %in% names(x)) {
    x <- x$enrich
  }
  if (!inherits(x, c("enrichResult", "gseaResult", "compareClusterResult"))) {
    cli::cli_abort(
      "{.arg x} must be an {.cls enrichResult}, {.cls gseaResult}, {.cls compareClusterResult}, or a list containing {.val enrich}."
    )
  }

  # --- gene-pathway membership table via get_pathway_genes -------------------
  long_tbl <- get_pathway_genes(
    x              = x,
    pathways       = pathways,
    top_n          = top_n,
    pathway_column = match.arg(pathway_column),
    gene_column    = match.arg(gene_column),
    gene_sep       = gene_sep,
    return_type    = "long"
  )

  if (!nrow(long_tbl)) {
    cli::cli_abort("No gene-pathway pairs found in the enrichment result.")
  }

  # --- filter to genes appearing in >= min_pathways pathways -----------------
  gene_counts <- table(long_tbl$gene)
  keep_genes  <- names(gene_counts[gene_counts >= min_pathways])

  if (!length(keep_genes)) {
    cli::cli_abort(
      "No genes appear in {.val {min_pathways}} or more pathways. Try reducing {.arg min_pathways}."
    )
  }

  long_tbl <- long_tbl[long_tbl$gene %in% keep_genes, , drop = FALSE]

  hub_genes <- names(gene_counts[gene_counts >= 2])

  # --- cap gene count for readability ----------------------------------------
  if (length(keep_genes) > max_genes) {
    # Within each tier, order by how many pathways a gene participates in and
    # then by identifier. Without this the cap keeps whichever genes happened to
    # come first in clusterProfiler's arbitrary ordering, which is not stable
    # across R sessions, so the same call could draw different genes.
    rank_tier <- function(g) {
      if (!length(g)) return(g)
      g[order(-as.integer(gene_counts[g]), g)]
    }
    hub_first  <- rank_tier(intersect(hub_genes, keep_genes))
    non_hub    <- rank_tier(setdiff(keep_genes, hub_first))
    keep_genes <- c(hub_first, non_hub)[seq_len(max_genes)]
    long_tbl   <- long_tbl[long_tbl$gene %in% keep_genes, , drop = FALSE]
    cli::cli_warn(
      "Gene count capped at {.val {max_genes}} for readability (hub genes prioritised)."
    )
  }

  # --- resolve fold-change if requested --------------------------------------
  fc_lookup <- NULL
  if (color_by %in% c("foldchange", "regulation")) {
    if (is.null(vista) || is.null(sample_comparison)) {
      cli::cli_warn(
        "{.arg vista} and {.arg sample_comparison} are required for {.val {color_by}} colouring. Falling back to {.val pathway}."
      )
      color_by <- "pathway"
    } else {
      stopifnot(inherits(vista, "VISTA"))
      comps <- comparisons(vista)
      if (!sample_comparison %in% names(comps)) {
        cli::cli_abort("Comparison {.val {sample_comparison}} not found in {.arg vista}.")
      }
      de_tbl <- comps[[sample_comparison]]

      fc_col <- intersect(c("log2fc", "log2FoldChange", "logFC"), names(de_tbl))[1]
      if (is.na(fc_col)) {
        cli::cli_warn("No fold-change column found in DE results. Falling back to {.val pathway} colouring.")
        color_by <- "pathway"
      } else {
        fc_ids  <- de_tbl$gene_id %||% rownames(de_tbl)
        fc_vals <- de_tbl[[fc_col]]
        fc_lookup <- stats::setNames(fc_vals, fc_ids)

        # map enrichment gene IDs -> vista rownames when gene_id_column supplied
        if (!is.null(gene_id_column)) {
          rd <- as.data.frame(SummarizedExperiment::rowData(vista), stringsAsFactors = FALSE)
          if (gene_id_column %in% colnames(rd)) {
            id_map <- stats::setNames(rownames(rd), rd[[gene_id_column]])
            enrichment_genes <- unique(long_tbl$gene)
            mapped <- id_map[enrichment_genes]
            mapped <- mapped[!is.na(mapped)]
            fc_rekeyed <- fc_lookup[mapped]
            names(fc_rekeyed) <- names(mapped)
            fc_lookup <- fc_rekeyed
          }
        }
      }
    }
  }

  # --- resolve display IDs --------------------------------------------------
  gene_labels <- stats::setNames(unique(long_tbl$gene), unique(long_tbl$gene))
  if (!is.null(display_id) && !is.null(vista)) {
    rd <- as.data.frame(SummarizedExperiment::rowData(vista), stringsAsFactors = FALSE)
    if (display_id %in% colnames(rd)) {
      lab_map <- stats::setNames(rd[[display_id]], rownames(rd))
      if (!is.null(gene_id_column) && gene_id_column %in% colnames(rd)) {
        id_map <- stats::setNames(rownames(rd), rd[[gene_id_column]])
        for (g in names(gene_labels)) {
          rn <- id_map[g]
          if (!is.na(rn) && rn %in% names(lab_map) && !is.na(lab_map[rn]) && nzchar(lab_map[rn])) {
            gene_labels[g] <- lab_map[rn]
          }
        }
      } else {
        for (g in names(gene_labels)) {
          if (g %in% names(lab_map) && !is.na(lab_map[g]) && nzchar(lab_map[g])) {
            gene_labels[g] <- lab_map[g]
          }
        }
      }
    }
  }

  # --- build data frame for chordDiagram ------------------------------------
  gene_order <- unique(long_tbl$gene)
  if (gene_order_by != "none") {
    if (is.null(fc_lookup)) {
      cli::cli_warn(
        "{.arg gene_order_by} = {.val {gene_order_by}} requires fold-change lookup ({.arg vista} + {.arg sample_comparison}). Using default gene order."
      )
    } else {
      gene_fc <- fc_lookup[gene_order]
      ord <- if (gene_order_by == "foldchange") {
        order(is.na(gene_fc), -gene_fc, gene_order)
      } else {
        order(is.na(gene_fc), -abs(gene_fc), -gene_fc, gene_order)
      }
      gene_order <- gene_order[ord]
    }
  }

  chord_df <- data.frame(
    from  = factor(long_tbl$pathway, levels = unique(long_tbl$pathway)),
    to    = factor(long_tbl$gene, levels = gene_order),
    value = 1,
    stringsAsFactors = FALSE
  )

  # --- colour setup ----------------------------------------------------------
  all_pathways <- levels(chord_df$from)
  all_genes    <- levels(chord_df$to)
  from_chr <- as.character(chord_df$from)
  to_chr <- as.character(chord_df$to)

  # pathway sector colours
  if (is.null(pathway_colors)) {
    pathway_colors <- colorspace::qualitative_hcl(length(all_pathways), palette = "Dark 3")
    names(pathway_colors) <- all_pathways
  } else {
    missing_pw <- setdiff(all_pathways, names(pathway_colors))
    if (length(missing_pw)) {
      extra <- colorspace::qualitative_hcl(length(missing_pw), palette = "Dark 3")
      names(extra) <- missing_pw
      pathway_colors <- c(pathway_colors, extra)
    }
  }

  # gene sectors: neutral grey (chords carry the colour information)
  gene_sector_colors <- rep("grey90", length(all_genes))
  names(gene_sector_colors) <- all_genes

  grid_col <- c(pathway_colors[all_pathways], gene_sector_colors)

  # chord (link) colours
  if (color_by == "pathway") {
    link_colors <- grDevices::adjustcolor(pathway_colors[from_chr], alpha.f = 1 - transparency)

  } else if (color_by == "regulation" && !is.null(fc_lookup)) {
    reg_colors <- vapply(to_chr, function(g) {
      fc <- fc_lookup[g]
      if (is.na(fc)) return(other_color)
      if (fc > 0) return(up_color)
      if (fc < 0) return(down_color)
      other_color
    }, character(1))
    link_colors <- grDevices::adjustcolor(reg_colors, alpha.f = 1 - transparency)

  } else if (color_by == "foldchange" && !is.null(fc_lookup)) {
    fc_for_genes <- fc_lookup[to_chr]
    fc_for_genes[is.na(fc_for_genes)] <- 0
    max_abs <- max(abs(fc_for_genes), na.rm = TRUE)
    if (max_abs == 0) max_abs <- 1
    col_fun <- circlize::colorRamp2(
      c(-max_abs, 0, max_abs),
      c(down_color, "white", up_color)
    )
    link_colors <- grDevices::adjustcolor(col_fun(fc_for_genes), alpha.f = 1 - transparency)

  } else {
    link_colors <- grDevices::adjustcolor(pathway_colors[from_chr], alpha.f = 1 - transparency)
  }

  # --- draw chord diagram ----------------------------------------------------
  sector_labels <- c(
    stats::setNames(
      ifelse(nchar(all_pathways) > 40,
             paste0(substr(all_pathways, 1, 37), "..."),
             all_pathways),
      all_pathways
    ),
    gene_labels[all_genes]
  )

  op <- graphics::par(mar = c(1, 1, 2, 1))
  on.exit(graphics::par(op), add = TRUE)

  circlize::circos.par(
    gap.degree   = gap_degree,
    start.degree = 90,
    clock.wise   = TRUE
  )

  circlize::chordDiagram(
    chord_df,
    order             = c(all_pathways, all_genes),
    grid.col          = grid_col,
    col               = link_colors,
    annotationTrack   = "grid",
    preAllocateTracks = list(track.height = 0.05)
  )

  # --- sector labels ---------------------------------------------------------
  circlize::circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      sector_index <- circlize::get.cell.meta.data("sector.index")
      xlim <- circlize::get.cell.meta.data("xlim")
      ylim <- circlize::get.cell.meta.data("ylim")
      lbl  <- sector_labels[sector_index]
      if (is.na(lbl)) lbl <- sector_index

      is_pathway <- sector_index %in% all_pathways
      circlize::circos.text(
        x      = mean(xlim),
        y      = ylim[1] + 0.5,
        labels = lbl,
        facing = "clockwise",
        niceFacing = TRUE,
        adj    = c(0, 0.5),
        cex    = if (is_pathway) label_cex else label_cex * 0.85,
        font   = if (is_pathway) 2L else 1L
      )
    },
    bg.border = NA
  )

  # --- title -----------------------------------------------------------------
  if (!is.null(title)) {
    graphics::title(main = title, cex.main = 1.1, line = 0.5)
  }

  # --- legend for FC / regulation colouring ----------------------------------
  if (color_by == "foldchange" && !is.null(fc_lookup) &&
      requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    fc_for_legend <- fc_lookup[all_genes]
    fc_for_legend <- fc_for_legend[!is.na(fc_for_legend)]
    if (length(fc_for_legend)) {
      max_abs <- max(abs(fc_for_legend), na.rm = TRUE)
      if (max_abs == 0) max_abs <- 1
      lgd <- ComplexHeatmap::Legend(
        col_fun   = circlize::colorRamp2(c(-max_abs, 0, max_abs), c(down_color, "white", up_color)),
        title     = "log2FC",
        direction = "vertical"
      )
      ComplexHeatmap::draw(lgd, x = grid::unit(2, "mm"), y = grid::unit(2, "mm"), just = c("left", "bottom"))
    }
  } else if (color_by == "regulation" &&
             requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    lgd <- ComplexHeatmap::Legend(
      labels    = c("Up", "Down", "Other"),
      legend_gp = grid::gpar(fill = c(up_color, down_color, other_color)),
      title     = "Regulation"
    )
    ComplexHeatmap::draw(lgd, x = grid::unit(2, "mm"), y = grid::unit(2, "mm"), just = c("left", "bottom"))
  }

  circlize::circos.clear()

  # --- return data invisibly -------------------------------------------------
  gene_data <- tibble::tibble(
    gene       = long_tbl$gene,
    pathway    = long_tbl$pathway,
    display    = gene_labels[long_tbl$gene],
    log2fc     = if (!is.null(fc_lookup)) fc_lookup[long_tbl$gene] else NA_real_,
    n_pathways = as.integer(gene_counts[long_tbl$gene])
  )
  if (gene_order_by != "none" && length(gene_order)) {
    gene_data$gene <- factor(gene_data$gene, levels = gene_order)
    gene_data <- gene_data[order(gene_data$gene, gene_data$pathway), , drop = FALSE]
    gene_data$gene <- as.character(gene_data$gene)
  }

  # circlize draws to the active device rather than returning an object, so
  # capture the result to let save_vista_plot() handle it like any other plot.
  recorded <- if (want_plot) .vista_record_plot() else NULL

  if (identical(return_type, "plot")) {
    if (is.null(recorded)) {
      cli::cli_abort(c(
        "Could not capture the chord diagram from the active graphics device.",
        "i" = "The device recorded an empty display list. Open a device before calling, e.g. {.code png(f); get_enrichment_chord(..., return_type = \"plot\"); dev.off()}."
      ))
    }
    return(recorded)
  }

  out <- list(gene_data = gene_data, hub_genes = hub_genes)
  if (identical(return_type, "both")) {
    if (is.null(recorded)) {
      cli::cli_warn(c(
        "Could not capture the chord diagram; {.field plot} will be NULL.",
        "i" = "Open a device before calling if you need the recorded plot."
      ))
    }
    out$plot <- recorded
    return(out)
  }
  invisible(out)
}
