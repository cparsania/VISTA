# ──────────────────────────────────────────────────────────────────────────────
# Helpers for new VISTA object layout
# ──────────────────────────────────────────────────────────────────────────────

#' @importFrom GGally ggpairs
#' @importFrom ggpubr ggbarplot stat_compare_means
#' @importFrom stringr str_replace str_detect str_c
#' @importFrom ggplot2 ggplot aes geom_point geom_density geom_violin geom_boxplot
#'   geom_jitter geom_line geom_segment geom_bar geom_tile scale_color_manual
#'   scale_fill_manual theme_minimal theme_bw facet_wrap facet_grid labs coord_flip
#'   element_text element_blank scale_x_discrete scale_y_continuous position_jitter
#' @importFrom tibble rownames_to_column as_tibble tibble
#' @importFrom dplyr filter left_join arrange mutate select group_by summarise pull
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment assay colData rowData assayNames
#' @importFrom colorspace qualitative_hcl
#' @importFrom stats prcomp cmdscale dist setNames
#' @importFrom cli cli_abort cli_warn
#' @importFrom rlang `%||%`
NULL

`%||%` <- function(a, b) if (!is.null(a)) a else b

.vista_group_info <- function(x) {
  meta <- S4Vectors::metadata(x)
  info <- meta$group
  if (is.null(info)) {
    return(list(
      column = colnames(SummarizedExperiment::colData(x))[1] %||% "group",
      colors = NULL,
      palette = "Dark 3"
    ))
  }
  info
}

.vista_group_col <- function(x) {
  info <- .vista_group_info(x)
  info$column %||% "group"
}

.vista_group_colors <- function(x, groups_present = NULL) {
  info <- .vista_group_info(x)
  cols <- info$colors
  if (!is.null(cols) && !is.null(groups_present)) {
    gp_levels <- unique(as.character(groups_present))
    cols <- cols[intersect(gp_levels, names(cols))]
    if (!length(cols)) cols <- NULL
  }
  if (is.null(cols) && !is.null(groups_present)) {
    pal_name <- (info$palette %||% "Dark 3")[1]
    pal <- colorspace::qualitative_hcl(length(unique(groups_present)), palette = pal_name)
    cols <- stats::setNames(pal, unique(groups_present))
  }
  cols
}

.vista_comparisons <- function(x) {
  comparisons(x)
}

.vista_comparison_colors <- function(x, comparisons_present = NULL) {
  meta <- S4Vectors::metadata(x)
  info <- meta$comparison %||% list()
  cols <- info$colors
  if (!is.null(cols) && !is.null(comparisons_present)) {
    comps <- unique(as.character(comparisons_present))
    cols <- cols[intersect(comps, names(cols))]
  }
  if (is.null(cols) && !is.null(comparisons_present)) {
    pal_name <- info$palette %||% "Dark 3"
    cols <- colorspace::qualitative_hcl(length(unique(comparisons_present)), palette = pal_name)
    cols <- stats::setNames(cols, unique(comparisons_present))
  }
  cols
}

.vista_deg_summary <- function(x) {
  deg_summary(x)
}

.resolve_plot_label_flag <- function(label = NULL,
                                     legacy = FALSE,
                                     legacy_arg = "label_replicates") {
  if (is.null(label)) {
    return(isTRUE(legacy))
  }
  if (!is.logical(label) || length(label) != 1L || is.na(label)) {
    cli::cli_abort("{.arg label} must be TRUE or FALSE.")
  }
  if (!identical(isTRUE(label), isTRUE(legacy))) {
    cli::cli_warn(
      "{.arg label} overrides the legacy {.arg {legacy_arg}} argument."
    )
  }
  isTRUE(label)
}

.resolve_plot_color_column <- function(meta,
                                       group_col,
                                       color_by = NULL,
                                       arg = "color_by") {
  color_col <- color_by %||% group_col
  if (!color_col %in% colnames(meta)) {
    cli::cli_abort(
      "Column {.val {color_col}} not found in {.arg sample_info}; cannot map {.arg {arg}}."
    )
  }
  color_col
}

.resolve_palette_values <- function(levels,
                                    colors = NULL,
                                    palette = NULL,
                                    default = "Dark 3") {
  levels <- unique(as.character(levels))
  if (!length(levels)) return(NULL)

  if (!is.null(colors)) {
    if (!is.character(colors) || !length(colors)) {
      cli::cli_abort("{.arg colors} must be a non-empty character vector when supplied.")
    }
    if (is.null(names(colors))) {
      if (length(colors) < length(levels)) {
        cli::cli_abort(
          "Not enough {.arg colors} supplied for levels: {.val {levels}}."
        )
      }
      colors <- colors[seq_along(levels)]
      names(colors) <- levels
    } else {
      missing_levels <- setdiff(levels, names(colors))
      if (length(missing_levels)) {
        cli::cli_abort(
          "{.arg colors} is missing values for levels: {.val {missing_levels}}."
        )
      }
      colors <- colors[levels]
    }
    return(colors)
  }

  pal <- colorspace::qualitative_hcl(length(levels), palette = palette %||% default)
  stats::setNames(pal, levels)
}

.resolve_embedding_colors <- function(x,
                                      values,
                                      color_col,
                                      group_col,
                                      use_group_colors = TRUE,
                                      palette = NULL,
                                      colors = NULL) {
  vals <- as.character(values)
  if (!length(vals)) return(NULL)

  if (isTRUE(use_group_colors) && identical(color_col, group_col)) {
    vista_cols <- .vista_group_colors(x, groups_present = vals)
    if (!is.null(vista_cols) && length(vista_cols)) {
      return(vista_cols)
    }
  }

  .resolve_palette_values(
    levels = vals,
    colors = colors,
    palette = palette,
    default = "Dark 3"
  )
}

.validate_facet_layout_value <- function(value, arg) {
  if (is.null(value)) {
    return(NULL)
  }
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || value < 1 || value != as.integer(value)) {
    cli::cli_abort("{.arg {arg}} must be NULL or a positive integer.")
  }
  as.integer(value)
}

.add_expression_facet_wrap <- function(plot,
                                       facet_formula,
                                       scales = "fixed",
                                       facet_nrow = NULL,
                                       facet_ncol = NULL) {
  facet_nrow <- .validate_facet_layout_value(facet_nrow, "facet_nrow")
  facet_ncol <- .validate_facet_layout_value(facet_ncol, "facet_ncol")

  plot + ggplot2::facet_wrap(
    facet_formula,
    scales = scales,
    nrow = facet_nrow,
    ncol = facet_ncol
  )
}

.resolve_expression_fill <- function(x,
                                     df,
                                     group_col,
                                     fill_by = NULL,
                                     default_fill = "group",
                                     allowed_special = c("group", "gene"),
                                     x_var = NULL,
                                     x_label = NULL) {
  fill_key <- fill_by %||% default_fill
  if (!is.character(fill_key) || length(fill_key) != 1L || is.na(fill_key)) {
    cli::cli_abort("{.arg fill_by} must be NULL or a single character value.")
  }

  if (identical(fill_key, "group")) {
    vals <- as.character(df[[group_col]])
    return(list(
      values = vals,
      levels = unique(vals),
      label = group_col,
      palette = .vista_group_colors(x, groups_present = vals)
    ))
  }

  if (identical(fill_key, "gene")) {
    if (!"gene" %in% colnames(df)) {
      cli::cli_abort("{.arg fill_by = 'gene'} requires a {.field gene} column in the plotting data.")
    }
    vals <- as.character(df$gene)
    return(list(
      values = vals,
      levels = unique(vals),
      label = "gene",
      palette = .resolve_palette_values(vals, default = "Dark 3")
    ))
  }

  if (identical(fill_key, "x")) {
    if (is.null(x_var) || !x_var %in% colnames(df)) {
      cli::cli_abort("{.arg fill_by = 'x'} is not available in this plotting configuration.")
    }
    vals <- as.character(df[[x_var]])
    palette <- if (identical(x_var, group_col)) {
      .vista_group_colors(x, groups_present = vals)
    } else {
      .resolve_palette_values(vals, default = "Dark 3")
    }
    return(list(
      values = vals,
      levels = unique(vals),
      label = x_label %||% x_var,
      palette = palette
    ))
  }

  if (!fill_key %in% colnames(df)) {
    allowed <- unique(c(allowed_special, colnames(df)))
    cli::cli_abort(
      "{.arg fill_by} must be one of {.val {allowed}}. Received: {.val {fill_key}}."
    )
  }

  vals <- as.character(df[[fill_key]])
  palette <- if (identical(fill_key, group_col)) {
    .vista_group_colors(x, groups_present = vals)
  } else {
    .resolve_palette_values(vals, default = "Dark 3")
  }
  list(
    values = vals,
    levels = unique(vals),
    label = fill_key,
    palette = palette
  )
}

.resolve_group_color_preference <- function(use_group_colors = TRUE,
                                            use_vista_colors = NULL) {
  if (!is.null(use_vista_colors)) {
    if (!is.logical(use_vista_colors) || length(use_vista_colors) != 1L || is.na(use_vista_colors)) {
      cli::cli_abort("{.arg use_vista_colors} must be TRUE or FALSE when supplied.")
    }
    cli::cli_warn("{.arg use_vista_colors} is deprecated; use {.arg use_group_colors} instead.")
    use_group_colors <- use_vista_colors
  }

  if (!is.logical(use_group_colors) || length(use_group_colors) != 1L || is.na(use_group_colors)) {
    cli::cli_abort("{.arg use_group_colors} must be TRUE or FALSE.")
  }

  use_group_colors
}

# Filter sample metadata to selected groups and keep a tidy tibble
.prepare_sample_metadata <- function(x, sample_group = NULL, group_column = NULL) {
  group_col <- group_column %||% .vista_group_col(x)
  meta <- as.data.frame(SummarizedExperiment::colData(x)) |>
    tibble::rownames_to_column("sample") |>
    tibble::as_tibble()
  stopifnot(group_col %in% colnames(meta))

  meta$sample <- factor(meta$sample, levels = rownames(SummarizedExperiment::colData(x)))

  if (!is.null(sample_group)) {
    meta <- dplyr::filter(meta, .data[[group_col]] %in% sample_group)
    # preserve requested group order when provided
    meta[[group_col]] <- factor(meta[[group_col]], levels = sample_group)
    meta <- dplyr::arrange(meta, .data[[group_col]], .data$sample)
  } else {
    meta[[group_col]] <- factor(meta[[group_col]], levels = unique(meta[[group_col]]))
    meta <- dplyr::arrange(meta, .data$sample)
  }

  meta$sample <- as.character(meta$sample)
  attr(meta, "group_column") <- group_col
  meta
}

# Create per-sample colors that stay consistent with the VISTA group palette.
.vista_sample_colors <- function(x, samples, group_column = NULL) {
  samples <- unique(as.character(samples))
  if (!length(samples)) return(character())

  meta <- .prepare_sample_metadata(x, sample_group = NULL, group_column = group_column)
  group_col <- attr(meta, "group_column")
  meta <- meta[match(samples, meta$sample), , drop = FALSE]

  if (anyNA(meta$sample)) {
    missing_samples <- samples[is.na(meta$sample)]
    cli::cli_abort(
      "Unknown sample(s) requested for colouring: {.val {missing_samples}}."
    )
  }

  group_vals <- as.character(meta[[group_col]])
  group_cols <- .vista_group_colors(x, groups_present = group_vals)
  if (is.null(group_cols) || !length(group_cols)) {
    pal <- colorspace::qualitative_hcl(length(samples), palette = "Dark 3")
    return(stats::setNames(pal, samples))
  }

  out <- stats::setNames(rep(NA_character_, length(samples)), samples)
  for (grp in unique(group_vals)) {
    idx <- which(group_vals == grp)
    grp_samples <- samples[idx]
    base_col <- group_cols[[grp]] %||% colorspace::qualitative_hcl(1, palette = "Dark 3")
    if (length(grp_samples) == 1) {
      out[grp_samples] <- base_col
    } else {
      ramp <- grDevices::colorRampPalette(
        c(
          colorspace::lighten(base_col, amount = 0.45),
          colorspace::darken(base_col, amount = 0.10)
        )
      )(length(grp_samples))
      names(ramp) <- grp_samples
      out[grp_samples] <- ramp
    }
  }

  out
}

# Keep only requested genes or the top-N variable genes; return matrix
.filter_genes <- function(mat, genes = NULL, top_n_genes = NULL) {
  if (!is.null(genes)) {
    if (length(genes) > 20) {
      cli::cli_abort("Maximum 20 genes can be plotted at once.")
    }
    keep <- intersect(genes, rownames(mat))
    if (!length(keep)) cli::cli_abort("None of the specified {.arg genes} were found in the data.")
    mat <- mat[keep, , drop = FALSE]
  }
  if (!is.null(top_n_genes)) {
    v <- matrixStats::rowVars(mat)
    ord <- order(v, decreasing = TRUE)
    n <- min(top_n_genes, nrow(mat))
    mat <- mat[ord[seq_len(n)], , drop = FALSE]
  }
  mat
}

# ──────────────────────────────────────────────────────────────────────────────
# PCA
# ──────────────────────────────────────────────────────────────────────────────

#' @title PCA plot
#' @description Uses normalized counts to compute principal components and plot samples,
#' optionally restricting to selected groups or genes.
#' @param x A `VISTA` object containing normalized counts.
#' @param sample_group Optional character vector of group labels (taken from the column specified by
#'   `group_column`, defaulting to the stored grouping column) used to subset samples prior to PCA. Use
#'   `NULL` to include all samples.
#' @param genes Optional character vector of gene identifiers to restrict the PCA input matrix.
#'   When `NULL`, all genes are used.
#' @param top_n_genes Optional integer selecting the top most variable genes to include. Ignored
#'   when `genes` is supplied.
#' @param label Logical; if `TRUE`, sample names are drawn next to the points.
#' @param label_size Numeric size of sample labels when `label = TRUE`.
#' @param point_size Numeric size of the plotted points.
#' @param shape_by Optional column name in `sample_info` used to map point shape. When `NULL`,
#'   shapes are not mapped.
#' @param shape_values Optional vector of shapes passed to `scale_shape_manual()` when
#'   `shape_by` is set. Use a named vector to map shapes to specific levels.
#' @param sample.seed Deprecated/unused; retained for backward compatibility.
#' @param show_clusters Logical; add normal ellipses per group when `TRUE`.
#' @param color_by Optional column name in `sample_info` used to map point colour.
#'   Defaults to the active grouping column.
#' @param use_vista_colors Deprecated alias for `use_group_colors`. When
#'   supplied, it overrides `use_group_colors`.
#' @param palette Optional qualitative palette name used when generating colours
#'   for non-group metadata levels.
#' @param colors Optional named character vector of manual colours overriding
#'   both `palette` and stored VISTA colours.
#' @param use_group_colors Logical; when `TRUE`, prefer the stored VISTA group
#'   colours when colouring by the grouping column.
#' @param group_column Optional column name in `sample_info` to use for grouping. Defaults to
#'   the stored grouping column.
#' @return A ggplot object showing the first two PCs.
#'
#' @examples
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
#'
#' # Basic PCA plot
#' get_pca_plot(vista)
#'
#' # With sample labels
#' get_pca_plot(vista, label = TRUE)
#'
#' # Using top variable genes
#' get_pca_plot(vista, top_n_genes = 100)
#'
#' # With confidence ellipses
#' get_pca_plot(vista, show_clusters = TRUE)
#'
#' @export
get_pca_plot <- function(x,
                         sample_group = NULL,
                         group_column = NULL,
                         genes = NULL,
                         top_n_genes = NULL,
                         label = FALSE,
                         label_size = 3,
                         point_size = 10,
                         shape_by = NULL,
                         shape_values = NULL,
                         sample.seed = 123,
                         show_clusters = FALSE,
                         color_by = NULL,
                         use_vista_colors = NULL,
                         palette = NULL,
                         colors = NULL,
                         use_group_colors = TRUE) {

  stopifnot(inherits(x, "VISTA"))
  use_group_colors <- .resolve_group_color_preference(
    use_group_colors = use_group_colors,
    use_vista_colors = use_vista_colors
  )

  mat <- SummarizedExperiment::assay(x)
  meta <- .prepare_sample_metadata(x, sample_group, group_column)
  group_col <- attr(meta, "group_column")
  color_col <- .resolve_plot_color_column(meta, group_col, color_by)
  if (!is.logical(label) || length(label) != 1L || is.na(label)) {
    cli::cli_abort("{.arg label} must be TRUE or FALSE.")
  }
  mat <- mat[, meta$sample, drop = FALSE]
  if (!is.null(shape_by) && !shape_by %in% colnames(meta)) {
    cli::cli_abort("Column {.val {shape_by}} not found in sample_info; cannot map shapes.")
  }

  mat <- .filter_genes(mat, genes, top_n_genes)
  mat <- mat[matrixStats::rowVars(mat) > 0, , drop = FALSE]

  pca <- stats::prcomp(t(mat), center = TRUE, scale. = TRUE)
  var_expl <- (pca$sdev^2) / sum(pca$sdev^2) * 100
  pca_df <- tibble::tibble(
    sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  ) |>
    dplyr::left_join(meta, by = "sample")

  if (!is.null(shape_by)) {
    pca_df[[shape_by]] <- factor(pca_df[[shape_by]])
  }
  pca_df[[color_col]] <- factor(pca_df[[color_col]])

  cols <- .resolve_embedding_colors(
    x = x,
    values = pca_df[[color_col]],
    color_col = color_col,
    group_col = group_col,
    use_group_colors = use_group_colors,
    palette = palette,
    colors = colors
  )

  if (is.null(shape_by)) {
    gp <- ggplot2::ggplot(
      pca_df,
      ggplot2::aes(x = PC1, y = PC2, color = .data[[color_col]], label = sample)
    )
  } else {
    gp <- ggplot2::ggplot(
      pca_df,
      ggplot2::aes(x = PC1, y = PC2, color = .data[[color_col]], shape = .data[[shape_by]], label = sample)
    )
  }

  gp <- gp +
    ggplot2::geom_point(size = point_size, alpha = 0.85) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      color = color_col,
      title = "PCA",
      x = sprintf("PC1 (%.1f%%)", var_expl[1]),
      y = sprintf("PC2 (%.1f%%)", var_expl[2])
    )

  if (!is.null(shape_by)) {
    gp <- gp + ggplot2::labs(shape = shape_by)
    if (!is.null(shape_values)) {
      shape_levels <- levels(pca_df[[shape_by]])
      if (is.null(names(shape_values))) {
        if (length(shape_values) < length(shape_levels)) {
          cli::cli_abort("Not enough {.arg shape_values} for {.arg shape_by} levels: {shape_levels}.")
        }
        shape_values <- shape_values[seq_along(shape_levels)]
        names(shape_values) <- shape_levels
      } else {
        missing_shapes <- setdiff(shape_levels, names(shape_values))
        if (length(missing_shapes) > 0) {
          cli::cli_abort("Missing {.arg shape_values} for levels: {missing_shapes}.")
        }
      }
      gp <- gp + ggplot2::scale_shape_manual(values = shape_values)
    }
  }

  if (!is.null(cols)) {
    gp <- gp + ggplot2::scale_color_manual(values = cols)
  }

  if (label) {
    gp <- gp + ggrepel::geom_text_repel(ggplot2::aes(label = sample), size = label_size)
  }

  if (show_clusters) {
    ellipse_df <- pca_df |>
      dplyr::group_by(.data[[group_col]]) |>
      dplyr::filter(dplyr::n() >= 3L) |>
      dplyr::ungroup()
    if (nrow(ellipse_df) == 0) {
      cli::cli_warn("Not enough points per group to draw ellipses (need >= 3 per group).")
    } else {
      gp <- gp + ggplot2::stat_ellipse(
        data = ellipse_df,
        ggplot2::aes(x = PC1, y = PC2, group = .data[[group_col]], color = .data[[group_col]]),
        type = "norm",
        linetype = 2,
        inherit.aes = FALSE
      )
    }
  }

  gp
}

# ──────────────────────────────────────────────────────────────────────────────
# MDS
# ──────────────────────────────────────────────────────────────────────────────

#' Generate an MDS plot for samples in a VISTA object
#'
#' Runs classical multidimensional scaling on normalized counts, optionally
#' restricting to groups or genes.
#'
#' @param x A `VISTA` object.
#' @param sample_group Optional character vector of groups to include (based on the column specified by `group_column`).
#' @param genes Optional character vector of gene identifiers to restrict the matrix.
#' @param top_n_genes Optional integer selecting the top variable genes to include.
#' @param label Logical; draw sample labels when `TRUE`.
#' @param label_size Numeric size of sample labels when `label = TRUE`.
#' @param point_size Numeric size for points.
#' @param shape_by Optional column name in `sample_info` used to map point shape. When `NULL`,
#'   shapes are not mapped.
#' @param shape_values Optional vector of shapes passed to `scale_shape_manual()` when
#'   `shape_by` is set. Use a named vector to map shapes to specific levels.
#' @param color_by Optional column name in `sample_info` used for point colour.
#'   Defaults to the active grouping column.
#' @param use_vista_colors Deprecated alias for `use_group_colors`. When
#'   supplied, it overrides `use_group_colors`.
#' @param palette Optional qualitative palette name used when generating colours
#'   for non-group metadata levels.
#' @param colors Optional named character vector of manual colours overriding
#'   both `palette` and stored VISTA colours.
#' @param use_group_colors Logical; when `TRUE`, prefer the stored VISTA group
#'   colours when colouring by the grouping column.
#' @param group_column Optional column name in `sample_info` to use for grouping/filtering.
#'
#' @export
get_mds_plot <- function(x,
                         sample_group = NULL,
                         group_column = NULL,
                         genes = NULL,
                         top_n_genes = NULL,
                         label = FALSE,
                         label_size = 3,
                         point_size = 10,
                         shape_by = NULL,
                         shape_values = NULL,
                         color_by = NULL,
                         use_vista_colors = NULL,
                         palette = NULL,
                         colors = NULL,
                         use_group_colors = TRUE) {

  stopifnot(inherits(x, "VISTA"))
  use_group_colors <- .resolve_group_color_preference(
    use_group_colors = use_group_colors,
    use_vista_colors = use_vista_colors
  )

  mat <- SummarizedExperiment::assay(x)
  meta <- .prepare_sample_metadata(x, sample_group, group_column)
  group_col <- attr(meta, "group_column")
  color_col <- .resolve_plot_color_column(meta, group_col, color_by)
  if (!is.logical(label) || length(label) != 1L || is.na(label)) {
    cli::cli_abort("{.arg label} must be TRUE or FALSE.")
  }
  mat <- mat[, meta$sample, drop = FALSE]
  if (!is.null(shape_by) && !shape_by %in% colnames(meta)) {
    cli::cli_abort("Column {.val {shape_by}} not found in sample_info; cannot map shapes.")
  }

  mat <- .filter_genes(mat, genes, top_n_genes)
  mat <- mat[matrixStats::rowVars(mat) > 0, , drop = FALSE]

  dist_matrix <- stats::dist(t(mat), method = "euclidean")
  mds <- stats::cmdscale(dist_matrix, k = 2, eig = TRUE)
  eig_vals <- mds$eig
  var_expl <- if (!is.null(eig_vals)) {
    prop <- eig_vals / sum(abs(eig_vals))
    prop[seq_len(2)] * 100
  } else {
    c(NA_real_, NA_real_)
  }
  mds_df <- tibble::tibble(
    sample = rownames(mds$points),
    Dim1 = mds$points[, 1],
    Dim2 = mds$points[, 2]
  ) |>
    dplyr::left_join(meta, by = "sample")

  if (!is.null(shape_by)) {
    mds_df[[shape_by]] <- factor(mds_df[[shape_by]])
  }
  mds_df[[color_col]] <- factor(mds_df[[color_col]])

  cols <- .resolve_embedding_colors(
    x = x,
    values = mds_df[[color_col]],
    color_col = color_col,
    group_col = group_col,
    use_group_colors = use_group_colors,
    palette = palette,
    colors = colors
  )

  if (is.null(shape_by)) {
    gp <- ggplot2::ggplot(
      mds_df,
      ggplot2::aes(x = Dim1, y = Dim2, color = .data[[color_col]], label = sample)
    )
  } else {
    gp <- ggplot2::ggplot(
      mds_df,
      ggplot2::aes(x = Dim1, y = Dim2, color = .data[[color_col]], shape = .data[[shape_by]], label = sample)
    )
  }

  gp <- gp +
    ggplot2::geom_point(size = point_size, alpha = 0.85) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      color = color_col,
      title = "MDS",
      x = if (!is.na(var_expl[1])) sprintf("Dim1 (%.1f%%)", var_expl[1]) else "Dim1",
      y = if (!is.na(var_expl[2])) sprintf("Dim2 (%.1f%%)", var_expl[2]) else "Dim2"
    )

  if (!is.null(shape_by)) {
    gp <- gp + ggplot2::labs(shape = shape_by)
    if (!is.null(shape_values)) {
      shape_levels <- levels(mds_df[[shape_by]])
      if (is.null(names(shape_values))) {
        if (length(shape_values) < length(shape_levels)) {
          cli::cli_abort("Not enough {.arg shape_values} for {.arg shape_by} levels: {shape_levels}.")
        }
        shape_values <- shape_values[seq_along(shape_levels)]
        names(shape_values) <- shape_levels
      } else {
        missing_shapes <- setdiff(shape_levels, names(shape_values))
        if (length(missing_shapes) > 0) {
          cli::cli_abort("Missing {.arg shape_values} for levels: {missing_shapes}.")
        }
      }
      gp <- gp + ggplot2::scale_shape_manual(values = shape_values)
    }
  }

  if (!is.null(cols)) {
    gp <- gp + ggplot2::scale_color_manual(values = cols)
  }
  if (label) {
    gp <- gp + ggrepel::geom_text_repel(ggplot2::aes(label = sample), size = label_size)
  }
  gp
}

# ──────────────────────────────────────────────────────────────────────────────
# UMAP
# ──────────────────────────────────────────────────────────────────────────────

#' Generate a UMAP plot for samples in a VISTA object
#'
#' Runs UMAP on normalized counts, optionally restricting to selected groups or
#' genes. UMAP is intended for exploratory sample-level structure.
#'
#' @param x A `VISTA` object.
#' @param sample_group Optional character vector of groups to include (based on
#'   `group_column`).
#' @param group_column Optional column name in `sample_info` used for
#'   filtering/grouping. Defaults to the stored grouping column.
#' @param color_by Optional column name in `sample_info` used for point color.
#'   Defaults to `group_column`.
#' @param genes Optional character vector of gene identifiers to restrict the matrix.
#' @param top_n_genes Optional integer selecting top variable genes to include.
#' @param label Logical; draw sample labels when `TRUE`.
#' @param label_size Numeric label size when `label = TRUE`.
#' @param point_size Numeric point size.
#' @param shape_by Optional column name in `sample_info` used to map point shape.
#' @param shape_values Optional vector passed to `scale_shape_manual()` when
#'   `shape_by` is set.
#' @param n_neighbors UMAP `n_neighbors` parameter.
#' @param min_dist UMAP `min_dist` parameter.
#' @param metric UMAP distance metric.
#' @param seed Integer random seed passed to UMAP.
#' @param use_vista_colors Deprecated alias for `use_group_colors`. When
#'   supplied, it overrides `use_group_colors`.
#' @param palette Optional qualitative palette name used when generating colours
#'   for non-group metadata levels.
#' @param colors Optional named character vector of manual colours overriding
#'   both `palette` and stored VISTA colours.
#' @param use_group_colors Logical; when `TRUE`, prefer the stored VISTA group
#'   colours when colouring by the grouping column.
#'
#' @return A ggplot object with UMAP1/UMAP2 coordinates.
#'
#' @examples
#' if (requireNamespace("uwot", quietly = TRUE)) {
#'   vista <- example_vista()
#'   get_umap_plot(vista, top_n_genes = 50)
#' }
#' @export
get_umap_plot <- function(x,
                          sample_group = NULL,
                          group_column = NULL,
                          color_by = NULL,
                          genes = NULL,
                          top_n_genes = NULL,
                          label = FALSE,
                          label_size = 3,
                          point_size = 10,
                          shape_by = NULL,
                          shape_values = NULL,
                          n_neighbors = 15,
                          min_dist = 0.1,
                          metric = "euclidean",
                          seed = 123,
                          use_vista_colors = NULL,
                          palette = NULL,
                          colors = NULL,
                          use_group_colors = TRUE) {

  stopifnot(inherits(x, "VISTA"))
  if (!requireNamespace("uwot", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg uwot} must be installed to compute UMAP.")
  }
  use_group_colors <- .resolve_group_color_preference(
    use_group_colors = use_group_colors,
    use_vista_colors = use_vista_colors
  )

  mat <- SummarizedExperiment::assay(x)
  meta <- .prepare_sample_metadata(x, sample_group, group_column)
  group_col <- attr(meta, "group_column")
  color_col <- .resolve_plot_color_column(meta, group_col, color_by)
  if (!is.logical(label) || length(label) != 1L || is.na(label)) {
    cli::cli_abort("{.arg label} must be TRUE or FALSE.")
  }
  mat <- mat[, meta$sample, drop = FALSE]
  if (!is.null(shape_by) && !shape_by %in% colnames(meta)) {
    cli::cli_abort("Column {.val {shape_by}} not found in sample_info; cannot map shapes.")
  }

  mat <- .filter_genes(mat, genes, top_n_genes)
  mat <- mat[matrixStats::rowVars(mat) > 0, , drop = FALSE]

  n_samples <- ncol(mat)
  if (n_samples < 3) {
    cli::cli_abort("UMAP needs at least 3 samples; found {.val {n_samples}}.")
  }
  if (n_neighbors >= n_samples) {
    adj <- max(2L, n_samples - 1L)
    cli::cli_warn(
      "{.arg n_neighbors} ({n_neighbors}) must be smaller than sample size ({n_samples}); using {adj}."
    )
    n_neighbors <- adj
  }

  um <- uwot::umap(
    t(mat),
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    metric = metric,
    n_components = 2,
    verbose = FALSE,
    ret_model = FALSE,
    seed = seed
  )

  umap_df <- tibble::tibble(
    sample = colnames(mat),
    UMAP1 = um[, 1],
    UMAP2 = um[, 2]
  ) |>
    dplyr::left_join(meta, by = "sample")

  umap_df[[color_col]] <- factor(umap_df[[color_col]])
  if (!is.null(shape_by)) {
    umap_df[[shape_by]] <- factor(umap_df[[shape_by]])
  }

  if (is.null(shape_by)) {
    gp <- ggplot2::ggplot(
      umap_df,
      ggplot2::aes(x = UMAP1, y = UMAP2, color = .data[[color_col]], label = sample)
    )
  } else {
    gp <- ggplot2::ggplot(
      umap_df,
      ggplot2::aes(x = UMAP1, y = UMAP2, color = .data[[color_col]], shape = .data[[shape_by]], label = sample)
    )
  }

  gp <- gp +
    ggplot2::geom_point(size = point_size, alpha = 0.85) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      color = color_col,
      title = "UMAP",
      x = "UMAP1",
      y = "UMAP2"
    )

  if (!is.null(shape_by)) {
    gp <- gp + ggplot2::labs(shape = shape_by)
    if (!is.null(shape_values)) {
      shape_levels <- levels(umap_df[[shape_by]])
      if (is.null(names(shape_values))) {
        if (length(shape_values) < length(shape_levels)) {
          cli::cli_abort("Not enough {.arg shape_values} for {.arg shape_by} levels: {shape_levels}.")
        }
        shape_values <- shape_values[seq_along(shape_levels)]
        names(shape_values) <- shape_levels
      } else {
        missing_shapes <- setdiff(shape_levels, names(shape_values))
        if (length(missing_shapes) > 0) {
          cli::cli_abort("Missing {.arg shape_values} for levels: {missing_shapes}.")
        }
      }
      gp <- gp + ggplot2::scale_shape_manual(values = shape_values)
    }
  }

  cols <- .resolve_embedding_colors(
    x = x,
    values = umap_df[[color_col]],
    color_col = color_col,
    group_col = group_col,
    use_group_colors = use_group_colors,
    palette = palette,
    colors = colors
  )
  if (!is.null(cols) && length(cols) > 0) {
    gp <- gp + ggplot2::scale_color_manual(values = cols)
  }

  if (label) {
    gp <- gp + ggrepel::geom_text_repel(ggplot2::aes(label = sample), size = label_size)
  }

  gp
}

# ──────────────────────────────────────────────────────────────────────────────
# Volcano
# ──────────────────────────────────────────────────────────────────────────────

#' Generate a volcano plot for a comparison in a VISTA object
#'
#' Wraps EnhancedVolcano to visualize log2FC vs p-values for a selected
#' comparison.
#'
#' @param x A `VISTA` object containing differential expression results.
#' @param sample_comparison Character scalar naming the comparison to display.
#' @param log2fc_cutoff Numeric absolute log2 fold-change threshold used to color significant points.
#' @param pval_cutoff Numeric p-value threshold used to color significant points.
#' @param label_genes Optional character vector of gene identifiers to force-label.
#' @param label_size Numeric label text size.
#' @param point_size Numeric point size.
#' @param colors Named colour vector with entries for `"Up"`, `"Down"`, and
#'   `"Other"`.
#' @param repair_genes Logical; when `TRUE`, split `gene_id` values like `ID:SYMBOL` to display the symbol.
#' @param display_id Optional ID/column name to use for plot labels. If supplied
#'   and present in `rowData(x)`, those values are used; otherwise falls back to
#'   ID mapping.
#' @param display_from Optional source ID type for mapping (used when `display_id`
#'   is not found in `rowData`).
#' @param display_orgdb Optional `OrgDb` object used for ID mapping when
#'   `display_id` is set but not found in `rowData`.
#' @param ... Additional parameters forwarded to `EnhancedVolcano::EnhancedVolcano()`.
#' @return A `ggplot2` object.
#'
#' @examples
#' vista <- example_vista()
#' comps <- names(comparisons(vista))
#' get_volcano_plot(vista, sample_comparison = comps[1])
#'
#' \donttest{
#' # Create VISTA object
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
#' # Basic volcano plot
#' comps <- names(comparisons(vista))
#' get_volcano_plot(vista, sample_comparison = comps[1])
#'
#' # With custom thresholds
#' get_volcano_plot(
#'   vista,
#'   sample_comparison = comps[1],
#'   log2fc_cutoff = 1.5,
#'   pval_cutoff = 0.01
#' )
#'
#' # Highlight specific genes
#' genes_of_interest <- rownames(vista)[seq_len(5)]
#' get_volcano_plot(
#'   vista,
#'   sample_comparison = comps[1],
#'   label_genes = genes_of_interest
#' )
#' }
#'
#' @export
get_volcano_plot <- function(x,
                             sample_comparison,
                             log2fc_cutoff = 1,
                             pval_cutoff = 0.05,
                             label_genes = NULL,
                             label_size = 3,
                             point_size = 1,
                             colors = c(Up = "#a40000", Down = "#007e2f", Other = "grey"),
                             repair_genes = TRUE,
                             display_id = NULL,
                             display_from = NULL,
                             display_orgdb = NULL,
                             ...) {
  dots <- list(...)

  stopifnot(inherits(x, "VISTA"))
  comps <- .vista_comparisons(x)
  if (!sample_comparison %in% names(comps)) {
    cli::cli_abort(c(
      "!" = "Invalid {.arg sample_comparison}.",
      "i" = "Available: {.val {names(comps)}}"
    ))
  }
  stopifnot(is.null(label_genes) || is.character(label_genes))
  color_names <- names(colors)
  colors <- as.character(colors)
  if (is.null(color_names) && length(colors) == 3L) {
    color_names <- c("Up", "Down", "Other")
  }
  names(colors) <- color_names
  if (is.null(names(colors)) || any(!c("Up", "Down", "Other") %in% names(colors))) {
    cli::cli_abort("{.arg colors} must be a named vector with entries for {.val Up}, {.val Down}, and {.val Other}.")
  }

  # Backward-compatible color aliases passed through `...`
  if ("col_up" %in% names(dots)) {
    colors[["Up"]] <- as.character(dots$col_up)[1]
    dots$col_up <- NULL
  }
  if ("col_down" %in% names(dots)) {
    colors[["Down"]] <- as.character(dots$col_down)[1]
    dots$col_down <- NULL
  }
  if ("col_other" %in% names(dots)) {
    colors[["Other"]] <- as.character(dots$col_other)[1]
    dots$col_other <- NULL
  }
  if ("col_others" %in% names(dots)) {
    colors[["Other"]] <- as.character(dots$col_others)[1]
    dots$col_others <- NULL
  }

  volcano_data <- comps[[sample_comparison]]
  gid_col <- if ("gene_id" %in% colnames(volcano_data)) "gene_id" else colnames(volcano_data)[1]

  gn <- volcano_data[[gid_col]] %||% rownames(volcano_data)
  rd <- tryCatch(SummarizedExperiment::rowData(x), error = function(e) NULL)
  if (!is.null(display_id) && !is.null(rd) && display_id %in% colnames(rd)) {
    map <- rd[[display_id]]
    names(map) <- rownames(x)
    mapped <- map[match(gn, names(map))]
    display <- ifelse(!is.na(mapped) & nzchar(mapped), mapped, gn)
  } else {
    display <- .map_gene_ids(gn, from_type = display_from, to_type = display_id, orgdb = display_orgdb)
  }
  if (repair_genes) {
    gsym <- stringr::str_replace(display, ".*:", "")
    gid <- stringr::str_replace(display, ":.*", "")
    lab <- dplyr::if_else(duplicated(gsym), gid, gsym)
  } else {
    lab <- display
  }

  vol_args <- c(
    list(
      toptable = volcano_data,
      lab = lab,
      x = "log2fc",
      y = "pvalue",
      pCutoff = pval_cutoff,
      FCcutoff = log2fc_cutoff,
      col_by_regul = TRUE,
      col_up = colors[["Up"]],
      col_down = colors[["Down"]],
      col_others = colors[["Other"]],
      selectLab = label_genes,
      labSize = label_size,
      pointSize = point_size,
      title = sample_comparison
    ),
    dots
  )
  do.call(.EnhancedVolcano2, vol_args)
}

# ──────────────────────────────────────────────────────────────────────────────
# Pairwise correlation (GGally)
# ──────────────────────────────────────────────────────────────────────────────

#' Plot pairwise correlations between samples
#'
#' Uses GGally::ggpairs on normalized expression to display correlations among
#' samples from selected groups/genes.
#'
#' @param x A `VISTA` object.
#' @param sample_group Optional character vector of groups (from the column specified by `group_column`) used to subset samples.
#' @param genes Optional character vector of gene IDs to include; defaults to all genes.
#' @param group_column Optional column name in `sample_info` defining the grouping used for filtering.
#' @param sample_order Ordering for selected samples: `"input"` keeps the
#'   current order, while `"group"` sorts by `group_column`.
#' @param value_transform Value transformation: `"log2"` (default) or `"none"`.
#' @param title Plot title.
#'
#' @export
get_pairwise_corr_plot <- function(x,
                                   sample_group = NULL,
                                   group_column = NULL,
                                   genes = NULL,
                                   sample_order = c("input", "group"),
                                   value_transform = c("log2", "none"),
                                   title = "Pairwise Sample Correlations") {
  stopifnot(inherits(x, "VISTA"))
  sample_order <- match.arg(sample_order)
  value_transform <- match.arg(value_transform)

  mat <- SummarizedExperiment::assay(x)
  meta <- .prepare_sample_metadata(x, sample_group, group_column)
  group_col <- attr(meta, "group_column")
  if (sample_order == "group" && group_col %in% colnames(meta)) {
    meta <- meta |>
      dplyr::arrange(.data[[group_col]], .data$sample)
  }
  mat <- mat[, meta$sample, drop = FALSE]

  if (!is.null(genes)) {
    missing_genes <- setdiff(genes, rownames(mat))
    if (length(missing_genes) > 0) {
      cli::cli_abort("Some genes are not found in normalized counts matrix: {.val {missing_genes}}")
    }
    mat <- mat[genes, , drop = FALSE]
  }

  if (value_transform == "log2") {
    mat <- log2(mat + 1)
  }
  GGally::ggpairs(as.data.frame(mat), title = title) +
    ggplot2::theme_minimal(base_size = 15)
}

# ──────────────────────────────────────────────────────────────────────────────
# Correlation heatmap
# ──────────────────────────────────────────────────────────────────────────────

#' Draw a sample correlation heatmap
#'
#' Plots sample-sample correlation matrix derived from normalized counts with
#' optional clustering and annotations.
#'
#' @param x A `VISTA` object.
#' @param sample_group Optional character vector of groups (referencing `group_column`) to include.
#' @param genes Optional character vector of gene IDs to limit the matrix.
#' @param corr_method Correlation method passed to `stats::cor()` (e.g., `"pearson"`).
#' @param triangle Either `"full"`, `"lower"`, or `"upper"` to control which triangle is drawn.
#' @param cluster_by Ordering strategy for samples: `"correlation"` (default),
#'   `"group"`, `"input"`, or `"none"`.
#' @param show_diagonal Logical; include the correlation diagonal when `TRUE`.
#' @param label Logical; overlay correlation coefficients as text.
#' @param show_corr_values Deprecated alias for `label`. When supplied, it
#'   overrides `label`.
#' @param label_color Color for the text labels.
#' @param col_corr_values Deprecated alias for `label_color`. When supplied, it
#'   overrides `label_color`.
#' @param label_size Numeric text size multiplier.
#' @param limits Optional numeric vector of length two giving limits for the
#'   color scale.
#' @param base_size Base theme size.
#' @param viridis_option Character viridis palette name.
#' @param viridis_direction Integer (1 or -1) controlling palette direction.
#' @param viridis_begin,viridis_end Palette endpoints between 0 and 1.
#' @param group_column Optional column name in `sample_info` defining the grouping used for filtering.
#' @aliases get_corr_heatmap
#'
#' @export
get_corr_heatmap <- function(x,
                             sample_group = NULL,
                             group_column = NULL,
                             genes = NULL,
                             corr_method = "pearson",
                             triangle = c("full", "lower", "upper"),
                             cluster_by = c("correlation", "group", "input", "none"),
                             show_diagonal = TRUE,
                             label = TRUE,
                             show_corr_values = NULL,
                             label_color = "black",
                             col_corr_values = NULL,
                             label_size = 4,
                             limits = NULL,
                             base_size = 12,
                             # NEW viridis options
                             viridis_option = "viridis",  # "magma","plasma","inferno","cividis","turbo"
                             viridis_direction = 1,       # 1 or -1
                             viridis_begin = 0,           # 0–1; raise to ~0.6 for more contrast near 1
                             viridis_end = 1) {

  stopifnot(inherits(x, "VISTA"))
  triangle <- match.arg(triangle)
  cluster_by <- match.arg(cluster_by)
  if (!is.null(show_corr_values)) {
    label <- isTRUE(show_corr_values)
  }
  if (!is.null(col_corr_values)) {
    label_color <- col_corr_values
  }

  mat <- SummarizedExperiment::assay(x)
  meta <- .prepare_sample_metadata(x, sample_group, group_column)
  group_col <- attr(meta, "group_column")
  mat <- mat[, meta$sample, drop = FALSE]

  if (!is.null(genes)) {
    missing_genes <- setdiff(genes, rownames(mat))
    if (length(missing_genes) > 0) {
      cli::cli_abort("Some genes are not found in normalized counts matrix: {.val {missing_genes}}")
    }
    mat <- mat[genes, , drop = FALSE]
  }

  mat <- log2(mat + 1)
  cor_mat <- stats::cor(mat, method = corr_method, use = "pairwise.complete.obs")

  df <- as.data.frame(as.table(cor_mat))
  colnames(df) <- c("Var1", "Var2", "value")

  if (!show_diagonal) df <- dplyr::filter(df, Var1 != Var2)
  if (triangle == "lower") {
    df <- dplyr::filter(df, as.integer(Var1) >= as.integer(Var2))
  } else if (triangle == "upper") {
    df <- dplyr::filter(df, as.integer(Var1) <= as.integer(Var2))
  }

  if (cluster_by == "correlation") {
    ord <- hclust(as.dist(1 - cor_mat))$order
    sample_levels <- rownames(cor_mat)[ord]
  } else if (cluster_by == "group" && group_col %in% colnames(meta)) {
    sample_levels <- meta |>
      dplyr::arrange(.data[[group_col]], .data$sample) |>
      dplyr::pull(.data$sample) |>
      unique()
  } else {
    sample_levels <- colnames(mat)
  }
  df$Var1 <- factor(df$Var1, levels = sample_levels)
  df$Var2 <- factor(df$Var2, levels = sample_levels)

  p <- ggplot2::ggplot(df, ggplot2::aes(Var1, Var2, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.title = ggplot2::element_blank())

  # viridis scale (sequential; works well when correlations are high)
  p <- p + ggplot2::scale_fill_viridis_c(
    option = viridis_option,
    direction = viridis_direction,
    begin = viridis_begin,
    end = viridis_end,
    limits = limits,         # e.g., c(0.8, 1)
    oob = scales::squish,
    name = NULL
  )

  if (label) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", value)),
      color = label_color,
      size = label_size / 3
    )
  }

  p
}

# ──────────────────────────────────────────────────────────────────────────────
# Expression boxplot
# ──────────────────────────────────────────────────────────────────────────────

#' Plot gene expression distributions as boxplots
#'
#' Displays per-sample or per-group distributions for selected genes using
#' normalized counts. The x-axis unit is controlled by `by`, and faceting is
#' controlled explicitly by `facet_by`.
#'
#' @param x A `VISTA` object.
#' @param genes Optional character vector of genes to display (<=20). Defaults to all genes.
#' @param sample_group Optional character vector specifying which groups (as defined by `group_column`) to include.
#' @param group_column Optional column name in `sample_info` used as the grouping variable.
#' @param log_transform Logical; apply log2(x + 1) transform before plotting.
#' @param display_id Optional column in `rowData(x)` to use for gene labels (facets).
#' @param display_from Optional source ID type for mapping (used when `display_id`
#'   is not found in `rowData`).
#' @param display_orgdb Optional `OrgDb` object used for ID mapping when
#'   `display_id` is set but not found in `rowData`.
#' @param facet_scales Facet scales argument passed to `facet_wrap()` (default `"free_y"`).
#' @param facet_nrow,facet_ncol Optional layout passed to `facet_wrap()` when faceting.
#' @param stats_group Logical; add statistical comparisons between groups when
#'   `TRUE`. Only supported when `by = "group"`.
#' @param p.label Label format for `ggpubr::stat_compare_means()`.
#' @param comparisons Optional list of specific group comparisons for `stat_compare_means()`.
#' @param pool_genes Logical; when `TRUE`, pool selected genes into one
#'   distribution per x-axis category (scenario 1).
#' @param by When `pool_genes = TRUE`, either `"group"` or `"sample"` (x-axis
#'   and optional fill). When `pool_genes = FALSE`, either `"group"` or
#'   `"gene"` (x-axis for per-gene distributions).
#' @param facet_by Faceting control: for pooled genes, `"group"` or `"none"`; for per-gene,
#'   `"none"` (default) or `"gene"`. `"auto"` selects the most readable layout.
#' @param fill_by Fill mapping. Special values are `"group"`, `"gene"`, and
#'   `"x"` (the plotted x-axis variable, when available). You may also supply a
#'   discrete column from the joined plotting data, such as a sample metadata
#'   column or `"sample"`.
#' @param sample_order Ordering used when sample names are shown on the x-axis:
#'   `"input"`, `"group"`, or `"expression"`.
#'
#' @export
get_expression_boxplot <- function(x,
                                   genes = NULL,
                                   sample_group = NULL,
                                   group_column = NULL,
                                   log_transform = TRUE,
                                   display_id = NULL,
                                   display_from = NULL,
                                   display_orgdb = NULL,
                                   facet_scales = "free_y",
                                   facet_nrow = NULL,
                                   facet_ncol = NULL,
                                   stats_group = FALSE,
                                   p.label = "p.signif",
                                   comparisons = NULL,
                                   pool_genes = FALSE,
                                   by = "group",
                                   facet_by = "auto",
                                   fill_by = NULL,
                                   sample_order = c("input", "group", "expression")) {
  stopifnot(inherits(x, "VISTA"))
  sample_order <- match.arg(sample_order)
  allowed_by <- if (pool_genes) c("group", "sample") else c("group", "gene")
  if (!is.null(by) && !by %in% allowed_by) {
    cli::cli_abort(
      "Invalid {.arg by} for {.arg pool_genes} = {pool_genes}. Use one of: {allowed_by}."
    )
  }
  allowed_facet_by <- if (pool_genes) c("auto", "group", "none") else c("auto", "gene", "none")
  if (!is.null(facet_by) && !facet_by %in% allowed_facet_by) {
    cli::cli_abort(
      "Invalid {.arg facet_by} for {.arg pool_genes} = {pool_genes}. Use one of: {allowed_facet_by}."
    )
  }
  by <- if (pool_genes) match.arg(by, c("group", "sample")) else match.arg(by, c("group", "gene"))
  facet_by <- if (pool_genes) {
    match.arg(facet_by %||% "auto", c("auto", "group", "none"))
  } else {
    match.arg(facet_by %||% "auto", c("auto", "gene", "none"))
  }
  if (identical(facet_by, "auto")) {
    facet_by <- if (pool_genes) "group" else if (!is.null(genes) && length(genes) > 1) "gene" else "none"
  }

  mat <- SummarizedExperiment::assay(x)
  meta <- .prepare_sample_metadata(x, sample_group, group_column)
  group_col <- attr(meta, "group_column")
  mat <- mat[, meta$sample, drop = FALSE]

  rd <- tryCatch(SummarizedExperiment::rowData(x), error = function(e) NULL)

  if (!is.null(genes)) {
    gene_ids <- genes
    if (!is.null(display_id) && !is.null(rd) && display_id %in% colnames(rd)) {
      map <- rd[[display_id]]
      names(map) <- rownames(x)
      mapped <- names(map)[match(genes, map)]
      mapped <- mapped[!is.na(mapped)]
      if (length(mapped)) gene_ids <- mapped
    }
    missing_genes <- setdiff(gene_ids, rownames(mat))
    if (length(missing_genes) == length(gene_ids)) {
      cli::cli_abort("None of the specified {.arg genes} were found in the data.")
    }
    keep_genes <- intersect(gene_ids, rownames(mat))
    mat <- mat[keep_genes, , drop = FALSE]
  } else if (!pool_genes && facet_by == "gene") {
    gene_vars <- matrixStats::rowVars(mat)
    if (all(is.na(gene_vars)) || all(gene_vars == 0)) {
      cli::cli_warn("All genes have zero variance; defaulting to the first 20 genes for faceting.")
      genes <- head(rownames(mat), 20)
    } else {
      top_idx <- order(gene_vars, decreasing = TRUE)[seq_len(min(20, length(gene_vars)))]
      genes <- rownames(mat)[top_idx]
      cli::cli_warn("No {.arg genes} provided; faceting by top 20 variable genes.")
    }
    mat <- mat[rownames(mat) %in% genes, , drop = FALSE]
  }

  df <- mat |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    tibble::as_tibble() |>
    tidyr::pivot_longer(-gene, names_to = "sample", values_to = "expression") |>
    dplyr::left_join(meta, by = "sample")

  # Replace gene labels with display_id when available
  if (!is.null(display_id) && !is.null(rd) && display_id %in% colnames(rd)) {
    map <- rd[[display_id]]
    names(map) <- rownames(x)
    mapped <- map[match(df$gene, names(map))]
    df$gene <- ifelse(!is.na(mapped) & nzchar(mapped), mapped, df$gene)
  }

  if (log_transform) {
    df$expression <- log2(df$expression + 1)
  }

  if (pool_genes) {
    df$gene <- "All genes"
  }

  x_var <- if (by == "group") {
    group_col
  } else if (pool_genes) {
    "sample"
  } else {
    "gene"
  }
  x_label <- if (by == "group") group_col else if (pool_genes) "sample" else "gene"

  if (pool_genes && identical(x_var, "sample")) {
    df <- .order_expression_plot_samples(df, meta, group_col, sample_order = sample_order)
  }
  fill_spec <- .resolve_expression_fill(
    x = x,
    df = df,
    group_col = group_col,
    fill_by = fill_by,
    default_fill = if (pool_genes) "x" else "group",
    allowed_special = if (pool_genes) c("x", "group", "gene") else c("group", "gene"),
    x_var = x_var,
    x_label = x_label
  )

  plt <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[x_var]],
      y = expression,
      fill = factor(fill_spec$values, levels = fill_spec$levels)
    )
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    ggplot2::labs(
      x = x_label,
      y = if (log_transform) "log2(Normalized Counts + 1)" else "Normalized Counts",
      fill = fill_spec$label
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (!is.null(fill_spec$palette)) {
    plt <- plt + ggplot2::scale_fill_manual(values = fill_spec$palette)
  }

  if (stats_group) {
    if (!requireNamespace("ggpubr", quietly = TRUE)) {
      cli::cli_abort("Package {.pkg ggpubr} must be installed for `stats_group = TRUE`.")
    }
    if (by == "group" && length(unique(df[[group_col]])) < 2) {
      cli::cli_warn("At least two groups are needed for statistical testing.")
    } else {
      plt <- plt + ggpubr::stat_compare_means(
        ggplot2::aes(group = .data[[x_var]]),
        comparisons = comparisons,
        method = "t.test",
        label = p.label,
        label.x.npc = "center"
      )
    }
  }

  if (facet_by != "none") {
    facet_var <- switch(facet_by,
                        gene = "gene",
                        group = group_col,
                        sample = "sample",
                        "none")
    if (facet_var != "none") {
      plt <- .add_expression_facet_wrap(
        plot = plt,
        facet_formula = stats::as.formula(paste0("~", facet_var)),
        scales = facet_scales,
        facet_nrow = facet_nrow,
        facet_ncol = facet_ncol
      )
    }
  }

  plt
}

# ──────────────────────────────────────────────────────────────────────────────
# Chromosome scatter (TxDb-based)
# ──────────────────────────────────────────────────────────────────────────────

#' Plot gene positions along chromosomes using a TxDb
#'
#' Retrieves gene coordinates on the fly from a user-supplied TxDb and plots
#' selected genes along chromosomes, optionally colouring by a numeric value and
#' labeling the most variable genes.
#'
#' @param x A `VISTA` object.
#' @param txdb A TxDb object (e.g., from \pkg{GenomicFeatures}).
#' @param keytype Key type in the TxDb matching \code{id_column} (default `"GENEID"`).
#' @param id_column Optional column in `rowData(x)` used to match to TxDb keys.
#'   When `NULL`, rownames(x) are used as keys.
#' @param genes Optional character vector of gene IDs to label (alternative to
#'   \code{label_top_n}). When provided, all genes are plotted but only these
#'   are labeled. Defaults to \code{NULL} (no explicit label set).
#' @param value_column Optional column in `rowData(x)` used for colouring.
#' @param comparison Optional comparison name; when supplied, uses `log2fc` from
#'   `metadata(x)$de_results[[comparison]]` for colouring. If multiple
#'   comparisons are provided, one panel per comparison is shown (log2FC
#'   clipped to +/-2).
#' @param group_value Optional group label (from `group_column`); when supplied,
#'   uses mean expression for that group for colouring (assay `norm_counts`).
#' @param label_top_n Integer; number of genes with largest |value| (or random if
#'   no value) to label. Ignored when \code{genes} is provided. Set to 0 to disable labels.
#' @param display_id Optional column in `rowData(x)` to use for point labels
#'   (fallback to gene_id/rownames).
#' @param display_from Optional source ID type for mapping when `display_id` is
#'   not present in `rowData(x)`.
#' @param display_orgdb Optional `OrgDb` object used for identifier mapping when
#'   `display_id` is not present in `rowData(x)`.
#' @param line_length Horizontal half-length (in megabases) of the tick used to
#'   mark each gene position. Default `0.02`. Increase for longer ticks.
#' @param line_width Line width of the tick marks. Default `0.6`.
#' @param filter_chrom Optional character vector of chromosomes to keep (e.g.,
#'   `c("chr1","chr2")`). When `NULL`, all chromosomes returned by the TxDb are shown.
#' @param value_label Optional legend title override for the colour scale.
#' @param scale_mode Colour scale mode: `"diverging"` (default) or
#'   `"sequential"`.
#' @param scale_low Optional low-end colour override.
#' @param scale_mid Optional midpoint colour override (used with diverging scale).
#' @param scale_high Optional high-end colour override.
#' @param scale_midpoint Numeric midpoint passed to
#'   `ggplot2::scale_color_gradient2()` when `scale_mode = "diverging"`.
#'
#' @details
#' The \code{genes} argument is only used as an explicit label set (all genes
#' are still plotted). Values in \code{genes} must match either the rownames of
#' \code{x} or the values in \code{id_column} when that is supplied. For example,
#' if \code{id_column = "ENTREZID"}, then \code{genes} should contain Entrez IDs
#' to be labeled.
#' When multiple comparisons are supplied, `value_column` and `group_value` are
#' ignored and the plot is facetted by comparison with a fixed +/-2 log2FC scale.
#'
#' @return A `ggplot2` object.
#' @keywords internal
get_chromosome_plot <- function(x,
                                txdb,
                                keytype = "GENEID",
                                id_column = NULL,
                                genes = NULL,
                                value_column = NULL,
                                comparison = NULL,
                                group_value = NULL,
                                label_top_n = 20,
                                display_id = NULL,
                                display_from = NULL,
                                display_orgdb = NULL,
                                line_length = 0.02,
                                line_width = 0.6,
                                filter_chrom = NULL,
                                value_label = NULL,
                                use_data_range = FALSE,
                                force_fc_limits = FALSE,
                                scale_mode = c("diverging", "sequential"),
                                scale_low = NULL,
                                scale_mid = NULL,
                                scale_high = NULL,
                                scale_midpoint = 0) {
  stopifnot(inherits(x, "VISTA"))
  scale_mode <- match.arg(scale_mode)
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg AnnotationDbi} is required for TxDb queries.")
  }
  if (!requireNamespace("viridis", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg viridis} is required for colouring.")
  }

  rd <- tryCatch(SummarizedExperiment::rowData(x), error = function(e) NULL)
  gene_ids <- rownames(x)
  if (!is.null(id_column) && !is.null(rd) && id_column %in% colnames(rd)) {
    map <- rd[[id_column]]
    names(map) <- rownames(x)
    mapped <- map
    gene_ids <- ifelse(!is.na(mapped) & nzchar(mapped), mapped, rownames(x))
  }

  # optional label set (do not drop other genes)
  label_set <- if (!is.null(genes)) as.character(genes) else NULL

  tx_cols <- c("TXCHROM", "TXSTART", "TXEND")
  coord_tbl <- AnnotationDbi::select(
    txdb,
    keys = unique(gene_ids),
    columns = tx_cols,
    keytype = keytype
  )
  coord_tbl <- coord_tbl[stats::complete.cases(coord_tbl[, tx_cols]), , drop = FALSE]
  if (!nrow(coord_tbl)) cli::cli_abort("No coordinates returned from TxDb for the requested genes.")

  # collapse to gene-level ranges
  coord_collapsed <- coord_tbl |>
    dplyr::group_by(.data[[keytype]]) |>
    dplyr::summarise(
      chr = dplyr::first(.data$TXCHROM),
      start = min(.data$TXSTART, na.rm = TRUE),
      end = max(.data$TXEND, na.rm = TRUE),
      .groups = "drop"
    )

  df_base <- tibble::tibble(
    gene_id = rownames(x),
    key = gene_ids
  ) |>
    dplyr::inner_join(coord_collapsed, by = c("key" = keytype))

  df_base$mid <- (df_base$start + df_base$end) / 2
  df_base$pos_mb <- df_base$mid / 1e6
  if (!is.null(filter_chrom)) {
    df_base <- dplyr::filter(df_base, .data$chr %in% filter_chrom)
    if (!nrow(df_base)) cli::cli_abort("No data remaining after filtering chromosomes.")
  }
  # order chromosomes naturally: numeric chr1..N, then X/Y/M/MT, then others
  chr_unique <- unique(df_base$chr)
  chr_lower <- tolower(chr_unique)
  num_mask <- grepl("^chr[0-9]+$", chr_lower)
  num_vals <- suppressWarnings(as.integer(sub("^chr", "", chr_lower[num_mask])))
  num_levels <- chr_unique[num_mask][order(num_vals, na.last = NA)]
  special_order <- c("chrx", "chry", "chrm", "chrmt")
  special_levels <- chr_unique[match(intersect(special_order, chr_lower), chr_lower)]
  other_levels <- setdiff(chr_unique, c(num_levels, special_levels))
  chr_levels <- c(num_levels, special_levels, other_levels)
  df_base$chr <- factor(df_base$chr, levels = chr_levels)

  # labels
  df_base$gene_label <- df_base$gene_id
  rd <- tryCatch(SummarizedExperiment::rowData(x), error = function(e) NULL)
  if (!is.null(display_id) && !is.null(rd) && display_id %in% colnames(rd)) {
    lab_map <- rd[[display_id]]
    names(lab_map) <- rownames(rd)
    lbl <- lab_map[match(df_base$gene_id, names(lab_map))]
    df_base$gene_label <- ifelse(!is.na(lbl) & nzchar(lbl), lbl, df_base$gene_label)
  } else if (!is.null(display_id) && !is.null(display_from) && !is.null(display_orgdb)) {
    lbl <- .map_gene_ids(
      ids = df_base$gene_id,
      from_type = display_from,
      to_type = display_id,
      orgdb = display_orgdb
    )
    df_base$gene_label <- ifelse(!is.na(lbl) & nzchar(lbl), lbl, df_base$gene_label)
  }

  df <- NULL
  value_label_final <- value_label %||% value_column %||% "value"
  facet_comparison <- FALSE
  facet_value_column <- FALSE
  use_data_range <- use_data_range

  # colour priority: comparison log2fc -> multiple rowData columns -> single rowData column -> group mean expression
  comps <- tryCatch(.vista_comparisons(x), error = function(e) list())
  if (!is.null(value_column) && length(value_column) > 1) {
    if (is.null(rd)) cli::cli_abort("rowData is required when supplying multiple value columns.")
    missing_cols <- setdiff(value_column, colnames(rd))
    if (length(missing_cols)) {
      cli::cli_abort("The following value columns were not found in rowData: {.val {missing_cols}}.")
    }
    if (!is.null(comparison)) {
      cli::cli_warn("`comparison` is ignored when multiple value columns are supplied.")
    }
    if (!is.null(group_value)) {
      cli::cli_warn("`group_value` is ignored when multiple value columns are supplied.")
    }
    df_list <- lapply(value_column, function(vc) {
      vals <- rd[[vc]]
      names(vals) <- rownames(rd)
      df_tmp <- df_base
      df_tmp$value <- vals[df_tmp$gene_id]
      df_tmp$value_source <- vc
      df_tmp
    })
    df <- dplyr::bind_rows(df_list)
    facet_value_column <- TRUE
    value_label_final <- value_label %||% "value"
    use_data_range <- TRUE
  } else if (!is.null(comparison) && length(comparison) > 1) {
    comp_avail <- intersect(comparison, names(comps))
    if (!length(comp_avail)) {
      cli::cli_abort("None of the requested comparisons were found in VISTA object.")
    }
    if (!is.null(value_column)) {
      cli::cli_warn("`value_column` is ignored when multiple comparisons are supplied.")
    }
    if (!is.null(group_value)) {
      cli::cli_warn("`group_value` is ignored when multiple comparisons are supplied.")
    }
    df_list <- lapply(comp_avail, function(comp_name) {
      de_tbl <- comps[[comp_name]]
      df_tmp <- df_base
      df_tmp$value <- NA_real_
      if ("gene_id" %in% names(de_tbl) && "log2fc" %in% names(de_tbl)) {
        val_map <- de_tbl$log2fc
        names(val_map) <- de_tbl$gene_id
        df_tmp$value <- val_map[df_tmp$gene_id]
      }
      df_tmp$comparison <- comp_name
      df_tmp
    })
    df <- dplyr::bind_rows(df_list)
    facet_comparison <- TRUE
    value_label_final <- "log2fc"
  } else {
    df <- df_base
    df$value <- NA_real_
    if (!is.null(comparison) && comparison %in% names(comps)) {
      de_tbl <- comps[[comparison]]
      if ("gene_id" %in% names(de_tbl) && "log2fc" %in% names(de_tbl)) {
        val_map <- de_tbl$log2fc
        names(val_map) <- de_tbl$gene_id
        df$value <- val_map[df$gene_id]
        value_label_final <- paste0(comparison, " log2fc")
      }
    } else if (!is.null(value_column) && !is.null(rd) && value_column %in% colnames(rd)) {
      vals <- rd[[value_column]]
      names(vals) <- rownames(rd)
      df$value <- vals[df$gene_id]
    } else if (!is.null(group_value)) {
      mat <- SummarizedExperiment::assay(x, "norm_counts")
      meta <- as.data.frame(SummarizedExperiment::colData(x))
      gcol <- .vista_group_col(x)
      if (!gcol %in% colnames(meta)) {
        cli::cli_warn("Group column {.val {gcol}} not found; skipping group-based colouring.")
      } else {
        samples_in_group <- rownames(meta)[meta[[gcol]] == group_value]
        if (length(samples_in_group)) {
          mat <- mat[df$gene_id, samples_in_group, drop = FALSE]
          df$value <- rowMeans(mat, na.rm = TRUE)
        } else {
          cli::cli_warn("No samples found for group {.val {group_value}}; colouring skipped.")
        }
      }
    }
  }

  # clip values only for log2fc (comparison) or when forced, if data-range is not requested
  if (!all(is.na(df$value)) && (isTRUE(force_fc_limits) || facet_comparison) && !isTRUE(use_data_range)) {
    df$value_clipped <- pmax(pmin(df$value, 2), -2)
  } else if (!all(is.na(df$value)) && !is.null(comparison) && !isTRUE(use_data_range)) {
    df$value_clipped <- pmax(pmin(df$value, 2), -2)
  } else {
    df$value_clipped <- df$value
  }

  if (all(is.na(df$value_clipped))) {
    plt <- ggplot2::ggplot(df, ggplot2::aes(y = chr)) +
      ggplot2::geom_segment(
        ggplot2::aes(x = pos_mb - line_length, xend = pos_mb + line_length, yend = chr),
        color = "grey50",
        linewidth = line_width
      ) +
      ggplot2::labs(x = "Genomic position (Mb)", y = "Chromosome") +
      ggplot2::theme_minimal()
  } else {
    limits_sym <- if ((isTRUE(force_fc_limits) || facet_comparison) && !isTRUE(use_data_range)) {
      c(-2, 2)
    } else if (isTRUE(use_data_range)) {
      rng <- range(df$value_clipped, na.rm = TRUE)
      if (is.finite(rng[1]) && is.finite(rng[2])) rng else NULL
    } else if (!is.null(comparison)) {
      c(-2, 2)
    } else {
      rng <- range(df$value_clipped, na.rm = TRUE)
      if (is.finite(rng[1]) && is.finite(rng[2])) rng else NULL
    }
    plt <- ggplot2::ggplot(df, ggplot2::aes(y = chr, color = value_clipped)) +
      ggplot2::geom_segment(
        ggplot2::aes(x = pos_mb - line_length, xend = pos_mb + line_length, yend = chr),
        linewidth = line_width
      ) +
      ggplot2::labs(x = "Genomic position (Mb)", y = "Chromosome", color = value_label_final) +
      ggplot2::theme_minimal()

    if (scale_mode == "sequential") {
      plt <- plt + ggplot2::scale_color_gradient(
        low = scale_low %||% "grey90",
        high = scale_high %||% "#2c7fb8",
        na.value = "grey70",
        limits = limits_sym
      )
    } else {
      plt <- plt + ggplot2::scale_color_gradient2(
        low = scale_low %||% "blue",
        mid = scale_mid %||% "grey80",
        high = scale_high %||% "red",
        midpoint = scale_midpoint,
        na.value = "grey70",
        limits = limits_sym
      )
    }
  }

  if (requireNamespace("scales", quietly = TRUE)) {
    plt <- plt + ggplot2::scale_x_continuous(labels = scales::label_number(suffix = " Mb"))
  }

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    if (!is.null(label_set)) {
      lab_df <- df[df$gene_id %in% label_set |
                     df$gene_label %in% label_set |
                     df$key %in% label_set, , drop = FALSE]
      if (!nrow(lab_df)) lab_df <- NULL
    } else if (label_top_n > 0) {
      if (all(is.na(df$value))) {
        lab_df <- df[sample.int(nrow(df), min(label_top_n, nrow(df))), ]
      } else {
        ord <- order(abs(df$value), decreasing = TRUE)
        lab_df <- df[ord[seq_len(min(label_top_n, length(ord)))], ]
      }
    } else {
      lab_df <- NULL
    }
    if (!is.null(lab_df)) {
      if (all(is.na(df$value_clipped))) {
        plt <- plt + ggrepel::geom_text_repel(
          data = lab_df,
          ggplot2::aes(x = pos_mb, y = chr, label = gene_label),
          size = 3,
          max.overlaps = Inf,
          color = "grey40"
        )
    } else {
      plt <- plt + ggrepel::geom_text_repel(
        data = lab_df,
        ggplot2::aes(x = pos_mb, y = chr, label = gene_label, color = value_clipped),
        size = 3,
        max.overlaps = Inf
      )
      }
    }
  }

    if (facet_comparison) {
      plt <- plt + ggplot2::facet_wrap(~comparison, ncol = 2)
    } else if (facet_value_column) {
      plt <- plt + ggplot2::facet_wrap(~value_source, ncol = 2)
    }

    plt
  }

#' Chromosome plot for expression
#'
#' Convenience wrapper around `get_chromosome_plot()` for expression-based
#' colouring (optional group mean, rowData columns, or assay columns).
#'
#' @inheritParams get_chromosome_plot
#' @param sample_group Optional character vector of groups (from `group_column`)
#'   used to subset assay columns when `value_from = "assay"`.
#' @param group_column Optional column name in `sample_info` used for
#'   `sample_group` selection and replicate summarisation.
#' @param value_from Source for `value_column` data: `"rowdata"` (default) or
#'   `"assay"`. When `"assay"`, the selected assay column is copied into
#'   `rowData` temporarily for colouring.
#' @param value_assay Assay name to pull values from when `value_from = "assay"`.
#'   Default `"norm_counts"`.
#' @param facet_value_columns Ignored (kept for compatibility); multiple
#'   `value_column`s are always arranged in a chromosome-by-column grid.
#' @param side_by_side_groups Logical; when plotting multiple `value_column`s,
#'   arrange one row per chromosome with one panel per selected sample/group.
#'   This is useful for direct side-by-side group comparison.
#' @param paginate Logical; when `TRUE`, split large multi-panel chromosome plots
#'   into pages to avoid patchwork viewport-size errors.
#' @param panels_per_page Maximum number of panels per page when `paginate = TRUE`.
#'   Set to `NULL` to disable automatic paging.
#' @param summarise_replicates Logical; when `TRUE` and `value_from = "assay"`,
#'   aggregate replicates by `group_column` before plotting.
#' @param summarise_method `"mean"` or `"median"` summarisation used when
#'   `summarise_replicates = TRUE`.
#' @param label_n Integer; number of genes with the largest absolute values to
#'   label when `genes` is not supplied. Set to 0 to disable automatic labels.
#' @param label_genes One of `"top"` (default), `"all"`, or `"none"` controlling
#'   chromosome text labels when `genes` is not supplied.
#' @param log_transform Logical; when `group_value` is used (and `value_column`
#'   is `NULL`) or when `value_from = "assay"`, apply log2(x + 1) before coloring.
#'
#' @details
#' - If multiple `value_column`s are supplied, one plot is produced per column
#'   with its own colour legend titled by that column. Chromosomes are laid out
#'   top-down in a single column: for each chromosome, plots for all
#'   `value_column`s appear sequentially. Legends for each column are shown only
#'   on the first occurrence (at the bottom) to avoid overlap. Requires
#'   \pkg{patchwork}; otherwise a list of plots is returned.
#' - Set `side_by_side_groups = TRUE` to place selected groups/samples in
#'   separate columns for each chromosome (same-row comparison).
#' - `paginate = TRUE` automatically splits large panel collections into a list
#'   of patchwork pages, which prevents "viewport too small" rendering errors in
#'   IDE plot panes.
#' - For `value_from = "assay"`, colours follow the VISTA ecosystem: group
#'   colours are used when `summarise_replicates = TRUE`, while sample colours
#'   are derived from group colours when plotting individual replicates.
#' - For `value_from = "assay"` with `value_column = NULL`, all selected
#'   assay columns are used (samples after `sample_group` filtering, or groups
#'   after replicate summarisation).
#' - Labels are kept consistent across `value_column`s: if `genes` is provided
#'   it is used for all panels; otherwise the top `label_n` (by absolute
#'   value in the first `value_column`) are used for all panels.
#' - When `value_from = "assay"`, the specified assay column is copied into
#'   `rowData` on the fly so it can be used for colouring.
#' - When `group_value` is provided (and no `value_column`), colouring is based
#'   on log2 mean expression for that group (assay `norm_counts`), using a
#'   data-driven colour scale. Ignored when `value_column` is supplied.
#'
#' @return A `ggplot2` object or a list of `ggplot2` objects when multiple
#'   `value_column`s are provided and \pkg{patchwork} is unavailable.
#' @export
get_expression_chromosome_plot <- function(x,
                                           txdb,
                                           keytype = "GENEID",
                                           id_column = NULL,
                                           genes = NULL,
                                           sample_group = NULL,
                                           group_column = NULL,
                                           value_column = NULL,
                                           value_from = c("rowdata", "assay"),
                                           value_assay = "norm_counts",
                                           facet_value_columns = FALSE,
                                           side_by_side_groups = FALSE,
                                           paginate = TRUE,
                                           panels_per_page = 24,
                                           summarise_replicates = FALSE,
                                           summarise_method = c("mean", "median"),
                                           group_value = NULL,
                                           label_n = 20,
                                           label_genes = c("top", "all", "none"),
                                           display_id = NULL,
                                           line_length = 0.02,
                                           line_width = 0.6,
                                           filter_chrom = NULL,
                                           log_transform = TRUE,
                                           value_label = "log2(mean expr)") {
  value_from <- match.arg(value_from)
  summarise_method <- match.arg(summarise_method)
  label_genes <- match.arg(label_genes)

  if (!is.null(value_column)) {
    value_column <- as.character(value_column)
  }

  assay_value_mat <- NULL
  column_color_map <- NULL

  if (value_from == "rowdata") {
    if (!is.null(sample_group)) {
      cli::cli_warn("`sample_group` is ignored when `value_from = \"rowdata\"`.")
    }
    if (!is.null(group_column)) {
      cli::cli_warn("`group_column` is ignored when `value_from = \"rowdata\"`.")
    }
    if (isTRUE(summarise_replicates)) {
      cli::cli_warn("`summarise_replicates` is ignored when `value_from = \"rowdata\"`.")
    }
    if (!is.null(value_column)) {
      rd_cols <- colnames(SummarizedExperiment::rowData(x))
      missing_cols <- setdiff(value_column, rd_cols)
      if (length(missing_cols)) {
        cli::cli_abort(
          "Value column(s) {.val {missing_cols}} not found in rowData. Set `value_from = \"assay\"` if you meant assay columns."
        )
      }
      column_color_map <- stats::setNames(
        colorspace::qualitative_hcl(length(value_column), palette = "Dark 3"),
        value_column
      )
    }
  } else {
    mat <- SummarizedExperiment::assay(x, value_assay)
    if (is.null(mat)) cli::cli_abort("Assay {.val {value_assay}} not found in object.")

    meta <- .prepare_sample_metadata(x, sample_group = sample_group, group_column = group_column)
    group_col <- attr(meta, "group_column")
    if (!nrow(meta)) cli::cli_abort("No samples available after filtering.")

    mat <- mat[, meta$sample, drop = FALSE]

    if (summarise_replicates) {
      grp_levels <- unique(as.character(meta[[group_col]]))
      summarised_list <- lapply(grp_levels, function(grp) {
        samples_grp <- meta$sample[as.character(meta[[group_col]]) == grp]
        vals <- mat[, samples_grp, drop = FALSE]
        if (summarise_method == "mean") {
          rowMeans(vals, na.rm = TRUE)
        } else {
          matrixStats::rowMedians(as.matrix(vals), na.rm = TRUE)
        }
      })
      assay_value_mat <- do.call(cbind, summarised_list)
      rownames(assay_value_mat) <- rownames(mat)
      colnames(assay_value_mat) <- grp_levels
      column_color_map <- .vista_group_colors(x, groups_present = grp_levels)
    } else {
      assay_value_mat <- as.matrix(mat)
      column_color_map <- .vista_sample_colors(x, samples = colnames(assay_value_mat), group_column = group_col)
    }

    if (is.null(value_column)) {
      value_column <- colnames(assay_value_mat)
    }
    missing_cols <- setdiff(value_column, colnames(assay_value_mat))
    if (length(missing_cols)) {
      cli::cli_abort(
        c(
          "Some {.arg value_column} entries were not found in selected assay columns.",
          "x" = "Missing: {.val {missing_cols}}",
          "i" = "Available: {.val {colnames(assay_value_mat)}}"
        )
      )
    }
    assay_value_mat <- assay_value_mat[, value_column, drop = FALSE]

    if (is.null(column_color_map) || !length(column_color_map)) {
      column_color_map <- stats::setNames(
        colorspace::qualitative_hcl(length(value_column), palette = "Dark 3"),
        value_column
      )
    }
    missing_color_names <- setdiff(value_column, names(column_color_map))
    if (length(missing_color_names)) {
      extra <- colorspace::qualitative_hcl(length(missing_color_names), palette = "Dark 3")
      names(extra) <- missing_color_names
      column_color_map <- c(column_color_map, extra)
    }
    column_color_map <- column_color_map[value_column]
    group_value <- NULL
  }

  if (!is.null(value_column)) {
    group_value <- NULL
  }

  resolve_label_setup <- function(reference_values = NULL) {
    if (!is.null(genes)) {
      return(list(label_set = as.character(genes), label_top_n = 0))
    }
    if (label_genes == "none") {
      return(list(label_set = NULL, label_top_n = 0))
    }
    if (label_genes == "all") {
      return(list(label_set = rownames(x), label_top_n = 0))
    }
    if (!is.null(reference_values) && label_n > 0) {
      ord <- order(abs(reference_values), decreasing = TRUE, na.last = NA)
      if (length(ord)) {
        ids <- names(reference_values)[ord[seq_len(min(label_n, length(ord)))]]
        return(list(label_set = ids, label_top_n = 0))
      }
    }
    list(label_set = NULL, label_top_n = label_n)
  }

  build_plot_for_value <- function(vc,
                                   filter_chr_local = filter_chrom,
                                   show_legend = TRUE,
                                   label_set = NULL,
                                   label_n_local = label_n) {
    x_local <- x
    use_value_column <- vc
    scale_mode_use <- "diverging"
    scale_low_use <- NULL
    scale_high_use <- NULL

    if (value_from == "assay") {
      vals <- assay_value_mat[, vc]
      if (log_transform) {
        vals <- log2(vals + 1)
      }
      tmp_col <- paste0(".tmp_expr_chr_", make.names(vc))
      rd_local <- as.data.frame(SummarizedExperiment::rowData(x_local))
      while (tmp_col %in% colnames(rd_local)) {
        tmp_col <- paste0(tmp_col, "_x")
      }
      rd_local[[tmp_col]] <- vals[rownames(x_local)]
      SummarizedExperiment::rowData(x_local) <- rd_local
      use_value_column <- tmp_col
      scale_mode_use <- "sequential"
      high_col <- column_color_map[[vc]] %||% "#2c7fb8"
      scale_high_use <- high_col
      scale_low_use <- colorspace::lighten(high_col, amount = 0.75)
    }

    plt <- get_chromosome_plot(
      x = x_local,
      txdb = txdb,
      keytype = keytype,
      id_column = id_column,
      genes = label_set,
      value_column = use_value_column,
      comparison = NULL,
      group_value = NULL,
      label_top_n = label_n_local,
      display_id = display_id,
      line_length = line_length,
      line_width = line_width,
      filter_chrom = filter_chr_local,
      value_label = vc,
      use_data_range = TRUE,
      force_fc_limits = FALSE,
      scale_mode = scale_mode_use,
      scale_low = scale_low_use,
      scale_high = scale_high_use
    )

    if (!isTRUE(show_legend)) {
      return(plt + ggplot2::theme(legend.position = "none"))
    }

    plt +
      ggplot2::theme(
        legend.position = "bottom",
        legend.box = "vertical",
        legend.direction = "horizontal"
      ) +
      ggplot2::guides(
        color = ggplot2::guide_colorbar(
          title.position = "top",
          barwidth = grid::unit(3.2, "cm"),
          barheight = grid::unit(0.30, "cm")
        )
      )
  }

  if (!is.null(value_column) && length(value_column) > 1) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      cli::cli_warn("Multiple value columns supplied; returning a list of plots (install {pkg patchwork} to combine automatically).")
      label_cfg <- resolve_label_setup(
        if (value_from == "assay") assay_value_mat[, value_column[1]] else NULL
      )
      return(lapply(value_column, function(vc) {
        build_plot_for_value(
          vc,
          filter_chr_local = filter_chrom,
          show_legend = TRUE,
          label_set = label_cfg$label_set,
          label_n = label_cfg$label_top_n
        )
      }))
    }

    if (!is.null(group_value)) {
      cli::cli_warn("`group_value` is ignored when multiple value columns are supplied.")
    }

    vals_first <- NULL
    if (value_from == "rowdata") {
      vals_first <- SummarizedExperiment::rowData(x)[[value_column[1]]]
      names(vals_first) <- rownames(x)
    } else {
      vals_first <- assay_value_mat[, value_column[1]]
      names(vals_first) <- rownames(x)
    }
    label_cfg <- resolve_label_setup(vals_first)

    # determine chromosome order once from TxDb
    rd <- tryCatch(SummarizedExperiment::rowData(x), error = function(e) NULL)
    gene_ids <- rownames(x)
    if (!is.null(id_column) && !is.null(rd) && id_column %in% colnames(rd)) {
      map <- rd[[id_column]]
      names(map) <- rownames(x)
      mapped <- map
      gene_ids <- ifelse(!is.na(mapped) & nzchar(mapped), mapped, rownames(x))
    }
    tx_cols <- c("TXCHROM", "TXSTART", "TXEND")
    coord_tbl <- AnnotationDbi::select(
      txdb,
      keys = unique(gene_ids),
      columns = tx_cols,
      keytype = keytype
    )
    coord_tbl <- coord_tbl[stats::complete.cases(coord_tbl[, tx_cols]), , drop = FALSE]
    coord_collapsed <- coord_tbl |>
      dplyr::group_by(.data[[keytype]]) |>
      dplyr::summarise(
        chr = dplyr::first(.data$TXCHROM),
        start = min(.data$TXSTART, na.rm = TRUE),
        end = max(.data$TXEND, na.rm = TRUE),
        .groups = "drop"
      )
    chr_unique <- unique(coord_collapsed$chr)
    chr_lower <- tolower(chr_unique)
    num_mask <- grepl("^chr[0-9]+$", chr_lower)
    num_vals <- suppressWarnings(as.integer(sub("^chr", "", chr_lower[num_mask])))
    num_levels <- chr_unique[num_mask][order(num_vals, na.last = NA)]
    special_order <- c("chrx", "chry", "chrm", "chrmt")
    special_levels <- chr_unique[match(intersect(special_order, chr_lower), chr_lower)]
    other_levels <- setdiff(chr_unique, c(num_levels, special_levels))
    chr_levels <- c(num_levels, special_levels, other_levels)
    if (!is.null(filter_chrom)) {
      chr_levels <- intersect(chr_levels, filter_chrom)
    }
    if (!length(chr_levels)) cli::cli_abort("No chromosomes to plot after filtering.")

    panels <- list()
    legend_shown <- stats::setNames(rep(FALSE, length(value_column)), value_column)
    for (chr in chr_levels) {
      for (vc in value_column) {
        show_leg <- !legend_shown[[vc]]
        legend_shown[[vc]] <- TRUE
        panels[[length(panels) + 1]] <- build_plot_for_value(
          vc,
          filter_chr_local = chr,
          show_legend = show_leg,
          label_set = label_cfg$label_set,
          label_n = label_cfg$label_top_n
        )
      }
    }
    ncol_wrap <- if (isTRUE(side_by_side_groups)) {
      length(value_column)
    } else if (isTRUE(facet_value_columns)) {
      min(2, length(value_column))
    } else {
      1
    }

    if (isTRUE(paginate) && !is.null(panels_per_page) && length(panels) > panels_per_page) {
      panels_per_page <- as.integer(panels_per_page)[1]
      if (is.na(panels_per_page) || panels_per_page < 1) {
        cli::cli_abort("{.arg panels_per_page} must be a positive integer or NULL.")
      }

      # Keep complete chromosome rows together in side-by-side mode.
      if (isTRUE(side_by_side_groups) && ncol_wrap > 1) {
        rows_per_page <- max(1L, panels_per_page %/% ncol_wrap)
        panels_per_page <- rows_per_page * ncol_wrap
      }

      idx <- split(seq_along(panels), ceiling(seq_along(panels) / panels_per_page))
      out <- lapply(idx, function(i) {
        patchwork::wrap_plots(panels[i], ncol = ncol_wrap, guides = "keep")
      })
      cli::cli_inform(
        "Returning {.val {length(out)}} plot pages; print one page at a time (for example, {.code p[[1]]})."
      )
      return(out)
    }

    return(patchwork::wrap_plots(panels, ncol = ncol_wrap, guides = "keep"))
  }

  if (!is.null(value_column) && length(value_column) == 1) {
    label_cfg <- resolve_label_setup(
      if (value_from == "assay") assay_value_mat[, value_column] else NULL
    )
    return(
      build_plot_for_value(
        value_column[[1]],
        filter_chr_local = filter_chrom,
        show_legend = TRUE,
        label_set = label_cfg$label_set,
        label_n = label_cfg$label_top_n
      )
    )
  }

  if (!is.null(value_column) && (is.null(value_label) || value_label == "log2(mean expr)")) {
    value_label <- if (length(value_column) > 1) "value" else value_column
  }

  label_cfg <- resolve_label_setup(NULL)

  # When colouring by group expression, apply log2 transform if requested
  if (!is.null(group_value) && log_transform) {
    mat <- SummarizedExperiment::assay(x, "norm_counts")
    meta <- as.data.frame(SummarizedExperiment::colData(x))
    gcol <- .vista_group_col(x)
    if (gcol %in% colnames(meta)) {
      samples_in_group <- rownames(meta)[meta[[gcol]] == group_value]
      if (length(samples_in_group)) {
        expr_vals <- rowMeans(mat[, samples_in_group, drop = FALSE], na.rm = TRUE)
        SummarizedExperiment::rowData(x)[[".tmp_group_expr"]] <- log2(expr_vals + 1)
        grp_cols <- .vista_group_colors(x, groups_present = group_value)
        high_col <- grp_cols[[group_value]] %||% "#2c7fb8"
        return(
          get_chromosome_plot(
            x = x,
            txdb = txdb,
            keytype = keytype,
            id_column = id_column,
            genes = label_cfg$label_set,
            value_column = ".tmp_group_expr",
            comparison = NULL,
            group_value = NULL,
            label_top_n = label_cfg$label_top_n,
            display_id = display_id,
            line_length = line_length,
            line_width = line_width,
            filter_chrom = filter_chrom,
            value_label = value_label,
            use_data_range = TRUE,
            force_fc_limits = FALSE,
            scale_mode = "sequential",
            scale_low = colorspace::lighten(high_col, amount = 0.75),
            scale_high = high_col
          ) +
            ggplot2::theme(
              legend.position = "bottom",
              legend.box = "vertical",
              legend.direction = "horizontal"
            )
        )
      }
    }
  }

  get_chromosome_plot(
    x = x,
    txdb = txdb,
    keytype = keytype,
    id_column = id_column,
    genes = label_cfg$label_set,
    value_column = value_column,
    comparison = NULL,
    group_value = group_value,
    label_top_n = label_cfg$label_top_n,
    display_id = display_id,
    line_length = line_length,
    line_width = line_width,
    filter_chrom = filter_chrom,
    value_label = value_label,
    use_data_range = TRUE,
    force_fc_limits = FALSE,
    scale_mode = if (!is.null(group_value)) "sequential" else "diverging"
  )
}

#' Chromosome plot for fold change
#'
#' Convenience wrapper around `get_chromosome_plot()` for fold-change colouring.
#' When multiple comparisons are supplied, panels are facetted by comparison
#' with log2FC clipped to +/-2.
#'
#' @inheritParams get_chromosome_plot
#' @param sample_comparison Optional comparison name (or vector of names) used
#'   for fold-change colouring.
#' @param label_n Integer; number of genes with the largest absolute fold-change
#'   to label when `genes` is not supplied. Set to 0 to disable labels.
#' @return A `ggplot2` object.
#' @export
get_foldchange_chromosome_plot <- function(x,
                                           txdb,
                                           keytype = "GENEID",
                                           id_column = NULL,
                                           genes = NULL,
                                           sample_comparison = NULL,
                                           value_column = NULL,
                                           label_n = 20,
                                           display_id = NULL,
                                           display_from = NULL,
                                           display_orgdb = NULL,
                                           line_length = 0.02,
                                           line_width = 0.6,
                                           filter_chrom = NULL) {
  get_chromosome_plot(
    x = x,
    txdb = txdb,
    keytype = keytype,
    id_column = id_column,
    genes = genes,
    value_column = value_column,
    comparison = sample_comparison,
    group_value = NULL,
    label_top_n = label_n,
    display_id = display_id,
    display_from = display_from,
    display_orgdb = display_orgdb,
    line_length = line_length,
    line_width = line_width,
    filter_chrom = filter_chrom,
    use_data_range = FALSE,
    force_fc_limits = TRUE
  )
}

# ──────────────────────────────────────────────────────────────────────────────
# Expression density plot
# ──────────────────────────────────────────────────────────────────────────────

#' Plot expression distributions as density curves
#'
#' Shows expression distributions pooled across the selected genes, coloured by
#' group (or sample), with optional faceting by group or sample.
#'
#' @param x A `VISTA` object.
#' @param genes Optional character vector of genes to display. Defaults to all genes.
#' @param sample_group Optional character vector specifying which groups (as defined by `group_column`) to include.
#' @param group_column Optional column name in `sample_info` used as the grouping variable.
#' @param log_transform Logical; apply log2(x + 1) transform before plotting.
#' @param facet_scales Facet scales argument passed to `facet_wrap()` (default `"free"`).
#' @param facet_nrow,facet_ncol Optional layout passed to `facet_wrap()` when faceting.
#' @param alpha Numeric transparency for density fill.
#' @param adjust Bandwidth adjustment factor passed to `geom_density()`.
#' @param color_by Either `"group"` (default) or `"sample"` to choose fill/color variable.
#' @param facet_by One of `"none"` (default), `"group"`, or `"sample"` to facet densities.
#' @param sample_order Ordering for sample-level display: `"input"`, `"group"`,
#'   or `"expression"`.
#' @param palette Optional qualitative palette name used when generating colours.
#' @param colors Optional named character vector of manual colours overriding
#'   `palette`.
#'
#' @return A `ggplot2` object.
#' @export
get_expression_density <- function(x,
                                   genes = NULL,
                                   sample_group = NULL,
                                   group_column = NULL,
                                   log_transform = TRUE,
                                   facet_scales = "free",
                                   facet_nrow = NULL,
                                   facet_ncol = NULL,
                                   alpha = 0.4,
                                   adjust = 1,
                                   color_by = c("group", "sample"),
                                   facet_by = c("none", "group", "sample"),
                                   sample_order = c("input", "group", "expression"),
                                   palette = NULL,
                                   colors = NULL) {
  stopifnot(inherits(x, "VISTA"))
  color_by <- match.arg(color_by)
  facet_by <- match.arg(facet_by)
  sample_order <- match.arg(sample_order)
  mat <- SummarizedExperiment::assay(x)
  meta <- .prepare_sample_metadata(x, sample_group, group_column)
  group_col <- attr(meta, "group_column")
  mat <- mat[, meta$sample, drop = FALSE]

  if (!is.null(genes)) {
    gene_ids <- genes
    keep <- intersect(gene_ids, rownames(mat))
    if (!length(keep)) cli::cli_abort("None of the specified {.arg genes} were found in the data.")
    mat <- mat[keep, , drop = FALSE]
  }

  df <- mat |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    tidyr::pivot_longer(-gene, names_to = "sample", values_to = "expression") |>
    dplyr::left_join(meta, by = "sample")

  if (log_transform) {
    df$expression <- log2(df$expression + 1)
    xlab <- "log2(Normalized Counts + 1)"
  } else {
    xlab <- "Normalized Counts"
  }
  df <- .order_expression_plot_samples(df, meta, group_col, sample_order = sample_order)

  fill_var <- if (color_by == "group") df[[group_col]] else df$sample
  fill_lab <- if (color_by == "group") group_col else "Sample"
  cols <- if (color_by == "group") {
    .resolve_embedding_colors(
      x = x,
      values = fill_var,
      color_col = group_col,
      group_col = group_col,
      use_group_colors = TRUE,
      palette = palette,
      colors = colors
    )
  } else {
    .resolve_palette_values(fill_var, colors = colors, palette = palette)
  }

  plt <- ggplot2::ggplot(df, ggplot2::aes(x = expression, fill = fill_var, color = fill_var)) +
    ggplot2::geom_density(alpha = alpha, adjust = adjust) +
    ggplot2::labs(x = xlab, y = "Density", fill = fill_lab, color = fill_lab) +
    ggplot2::theme_minimal()

  if (!is.null(cols)) {
    plt <- plt + ggplot2::scale_fill_manual(values = cols) + ggplot2::scale_color_manual(values = cols)
  }

  if (facet_by != "none") {
    facet_var <- if (facet_by == "group") group_col else "sample"
    plt <- .add_expression_facet_wrap(
      plot = plt,
      facet_formula = stats::as.formula(paste0("~", facet_var)),
      scales = facet_scales,
      facet_nrow = facet_nrow,
      facet_ncol = facet_ncol
    )
  }

  plt
}

# ──────────────────────────────────────────────────────────────────────────────
# Expression ridgeline (joy) plot
# ──────────────────────────────────────────────────────────────────────────────

#' Plot expression distributions as ridgelines
#'
#' Shows per-group (or per-sample) expression distributions pooled across the
#' selected genes using ridge (joy) plots. Genes are pooled; no faceting to keep
#' shapes comparable.
#'
#' @param x A `VISTA` object.
#' @param genes Optional character vector of genes to display. Defaults to all genes.
#' @param sample_group Optional character vector specifying which groups (as defined by `group_column`) to include.
#' @param group_column Optional column name in `sample_info` used as the grouping variable.
#' @param log_transform Logical; apply log2(x + 1) transform before plotting.
#' @param alpha Numeric transparency for fills.
#' @param scale Numeric scaling factor for ridges (passed to `geom_density_ridges()`).
#' @param y_by Either `"group"` (default) or `"sample"` to choose the y-axis
#'   (ridge) grouping.
#' @param color_by Either `"group"` (default) or `"sample"` to choose fill colors.
#' @param sample_order Ordering for sample-level display: `"input"`, `"group"`,
#'   or `"expression"`.
#' @param palette Optional qualitative palette name used when generating colours.
#' @param colors Optional named character vector of manual colours overriding
#'   `palette`.
#'
#' @return A `ggplot2` object.
#' @export
get_expression_joyplot <- function(x,
                                   genes = NULL,
                                   sample_group = NULL,
                                   group_column = NULL,
                                   log_transform = TRUE,
                                   alpha = 0.7,
                                   scale = 1.2,
                                   y_by = c("group", "sample"),
                                   color_by = c("group", "sample"),
                                   sample_order = c("input", "group", "expression"),
                                   palette = NULL,
                                   colors = NULL) {
  stopifnot(inherits(x, "VISTA"))
  if (!requireNamespace("ggridges", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggridges} must be installed to use `get_expression_joyplot()`.")
  }
  y_by <- match.arg(y_by)
  color_by <- match.arg(color_by)
  sample_order <- match.arg(sample_order)

  mat <- SummarizedExperiment::assay(x)
  meta <- .prepare_sample_metadata(x, sample_group, group_column)
  group_col <- attr(meta, "group_column")
  mat <- mat[, meta$sample, drop = FALSE]

  if (!is.null(genes)) {
    keep <- intersect(genes, rownames(mat))
    if (!length(keep)) cli::cli_abort("None of the specified {.arg genes} were found in the data.")
    mat <- mat[keep, , drop = FALSE]
  }

  df <- mat |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    tidyr::pivot_longer(-gene, names_to = "sample", values_to = "expression") |>
    dplyr::left_join(meta, by = "sample")

  if (log_transform) {
    df$expression <- log2(df$expression + 1)
    xlab <- "log2(Normalized Counts + 1)"
  } else {
    xlab <- "Normalized Counts"
  }
  df <- .order_expression_plot_samples(df, meta, group_col, sample_order = sample_order)

  y_var <- if (y_by == "group") group_col else "sample"
  fill_var <- if (color_by == "group") df[[group_col]] else df$sample
  fill_lab <- if (color_by == "group") group_col else "Sample"
  cols <- if (color_by == "group") {
    .resolve_embedding_colors(
      x = x,
      values = df[[group_col]],
      color_col = group_col,
      group_col = group_col,
      use_group_colors = TRUE,
      palette = palette,
      colors = colors
    )
  } else {
    .resolve_palette_values(df$sample, colors = colors, palette = palette)
  }

  plt <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = expression, y = .data[[y_var]], fill = fill_var)
  ) +
    ggridges::geom_density_ridges(alpha = alpha, scale = scale, rel_min_height = 0.01, color = "black", linewidth = 0.2) +
    ggplot2::labs(x = xlab, y = y_var, fill = fill_lab) +
    ggridges::theme_ridges() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (!is.null(cols)) {
    plt <- plt + ggplot2::scale_fill_manual(values = cols)
  }

  plt
}

# ──────────────────────────────────────────────────────────────────────────────
# Expression scatter (two samples/groups)
# ──────────────────────────────────────────────────────────────────────────────

#' Compare normalized expression between two samples or groups
#'
#' Plots gene-level expression for two selected samples or group means,
#' colours points by local density (viridis), labels the most divergent genes,
#' and reports Pearson/Spearman correlation.
#'
#' @param x A `VISTA` object.
#' @param sample_x First sample or group to plot (character scalar).
#' @param sample_y Second sample or group to plot (character scalar).
#' @param by One of `"sample"` (use individual samples) or `"group"` (average
#'   replicates within `group_column` before plotting). Default `"sample"`.
#' @param group_column Column in `sample_info` used when `by = "group"` (defaults
#'   to stored grouping column).
#' @param genes Optional character vector of genes to include; defaults to all.
#' @param log_transform Logical; apply log2(x + 1) transform. Default `TRUE`.
#' @param label_n Integer; number of most divergent genes to label
#'   (ranked by |x - y|). Set to 0 to disable labels.
#' @param label_size Numeric size for labeled genes.
#' @param method Correlation method for the subtitle; `"pearson"` (default) or
#'   `"spearman"`.
#' @param display_id Optional column in `rowData(x)` to use for point labels
#'   (fallback to gene_id/rownames when not available).
#'
#' @return A `ggplot2` object.
#' @export
get_expression_scatter <- function(x,
                                   sample_x,
                                   sample_y,
                                   by = c("sample", "group"),
                                   group_column = NULL,
                                   genes = NULL,
                                   log_transform = TRUE,
                                   label_n = 20,
                                   label_size = 3,
                                   method = c("pearson", "spearman"),
                                   display_id = NULL) {
  stopifnot(inherits(x, "VISTA"))
  by <- match.arg(by)
  method <- match.arg(method)

  mat <- SummarizedExperiment::assay(x)
  meta <- .prepare_sample_metadata(x, sample_group = NULL, group_column = group_column)
  group_col <- attr(meta, "group_column")

  if (!is.null(genes)) {
    keep <- intersect(genes, rownames(mat))
    if (!length(keep)) cli::cli_abort("None of the specified {.arg genes} were found in the data.")
    mat <- mat[keep, , drop = FALSE]
  }

  if (by == "group") {
    if (is.null(sample_x) || is.null(sample_y)) cli::cli_abort("Please provide group names for sample_x and sample_y when by = 'group'.")
    stopifnot(group_col %in% colnames(meta))
    group_means <- vapply(split(seq_len(nrow(meta)), meta[[group_col]]), function(idx) {
      rowMeans(mat[, idx, drop = FALSE], na.rm = TRUE)
    }, FUN.VALUE = numeric(nrow(mat)))
    if (is.null(dim(group_means))) group_means <- matrix(group_means, ncol = 1)
    colnames(group_means) <- names(split(seq_len(nrow(meta)), meta[[group_col]]))
    if (!all(c(sample_x, sample_y) %in% colnames(group_means))) {
      cli::cli_abort("Requested groups not found in {.field sample_info}: {.val {setdiff(c(sample_x, sample_y), colnames(group_means))}}")
    }
    xvec <- group_means[, sample_x]
    yvec <- group_means[, sample_y]
    xlab <- sprintf("%s (group mean)", sample_x)
    ylab <- sprintf("%s (group mean)", sample_y)
  } else {
    if (!all(c(sample_x, sample_y) %in% colnames(mat))) {
      cli::cli_abort("Requested samples not found in expression matrix: {.val {setdiff(c(sample_x, sample_y), colnames(mat))}}")
    }
    xvec <- mat[, sample_x]
    yvec <- mat[, sample_y]
    xlab <- sample_x
    ylab <- sample_y
  }

  df <- tibble::tibble(
    gene = rownames(mat),
    expr_x = xvec,
    expr_y = yvec
  )

  rd <- tryCatch(SummarizedExperiment::rowData(x), error = function(e) NULL)
  if (!is.null(display_id) && !is.null(rd) && display_id %in% colnames(rd)) {
    lab_map <- rd[[display_id]]
    names(lab_map) <- rownames(rd)
    labs <- lab_map[match(df$gene, names(lab_map))]
    df$gene <- ifelse(!is.na(labs) & nzchar(labs), labs, df$gene)
  }

  if (log_transform) {
    df$expr_x <- log2(df$expr_x + 1)
    df$expr_y <- log2(df$expr_y + 1)
    xlab <- paste0(xlab, " (log2+1)")
    ylab <- paste0(ylab, " (log2+1)")
  }

  df$abs_diff <- abs(df$expr_x - df$expr_y)

  cor_val <- suppressWarnings(stats::cor(df$expr_x, df$expr_y, method = method, use = "complete.obs"))
  cor_test <- suppressWarnings(stats::cor.test(df$expr_x, df$expr_y, method = method))
  r_label <- sprintf("italic(%s~r==%.3f)", tools::toTitleCase(method), cor_val)
  p_label <- format(cor_test$p.value, digits = 2, scientific = TRUE)

  if (requireNamespace("ggpointdensity", quietly = TRUE)) {
    plt <- ggplot2::ggplot(df, ggplot2::aes(x = expr_x, y = expr_y)) +
      ggpointdensity::geom_pointdensity() +
      ggplot2::scale_color_viridis_c(name = "Density") +
      ggplot2::labs(x = xlab, y = ylab) +
      ggplot2::theme_minimal()
  } else {
    plt <- ggplot2::ggplot(df, ggplot2::aes(x = expr_x, y = expr_y)) +
      ggplot2::geom_point(alpha = 0.4, size = 1.2, color = "steelblue") +
      ggplot2::labs(x = xlab, y = ylab) +
      ggplot2::theme_minimal()
  }

  xr <- range(df$expr_x, na.rm = TRUE)
  yr <- range(df$expr_y, na.rm = TRUE)
  x_pos <- xr[1] + 0.05 * diff(xr)
  y_top <- yr[2] - 0.05 * diff(yr)
  line_gap <- 0.05 * diff(yr)
  plt <- plt +
    ggplot2::annotate("text", x = x_pos, y = y_top,
                      label = r_label, hjust = 0, vjust = 1,
                      size = 3.5, color = "red", parse = TRUE) +
    ggplot2::annotate("text", x = x_pos, y = y_top - line_gap,
                      label = paste0("italic(p)==", p_label),
                      hjust = 0, vjust = 1, size = 3.2,
                      color = "red", parse = TRUE)

  if (label_n > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    top_idx <- order(df$abs_diff, decreasing = TRUE)
    top_genes <- df[top_idx[seq_len(min(label_n, nrow(df)))], ]
    plt <- plt + ggrepel::geom_text_repel(
      data = top_genes,
      ggplot2::aes(label = gene),
      size = label_size,
      max.overlaps = Inf
    )
  }

  plt
}

# ──────────────────────────────────────────────────────────────────────────────
# Expression lollipop
# ──────────────────────────────────────────────────────────────────────────────

#' Plot expression as a lollipop chart
#'
#' Displays selected genes as stem-and-dot (lollipop) plots. By default,
#' expression is summarized per group (mean across replicates). With
#' `by = "sample"`, the plot shows individual samples.
#'
#' @param x A `VISTA` object.
#' @param genes Character vector (<=15 genes) to plot.
#' @param sample_group Optional character vector of groups (from `group_column`) to include.
#' @param group_column Optional column name in `sample_info` to use for grouping samples.
#' @param by One of `"group"` (default; summarize replicates by group) or
#'   `"sample"` (show individual samples).
#' @param sample_order Ordering used when `by = "sample"`:
#'   `"input"` preserves the current sample order, `"group"` groups samples by
#'   `group_column`, and `"expression"` ranks samples by mean expression across
#'   the selected genes.
#' @param facet_by Faceting mode: `"auto"` (default; facet by gene when more
#'   than one gene is requested), `"gene"`, or `"none"`. For multiple genes,
#'   `"none"` falls back to `"gene"` to avoid unreadable combined panels.
#' @param log_transform Logical; log2-transform expression before plotting.
#' @param facet_scale Scaling option passed to `facet_wrap()` when plotting multiple genes.
#' @param facet_nrow,facet_ncol Optional layout passed to `facet_wrap()` when faceting.
#' @param point_size Numeric size of the dots.
#' @param line_size Numeric size of the stems.
#' @param label Logical; draw numeric labels above the dots.
#' @param label_digits Integer; digits to show in labels when `label = TRUE`.
#' @param display_id Optional ID/column name to use for labels/facets. If supplied
#'   and present in `rowData(x)`, those values are used; otherwise falls back to
#'   ID mapping.
#' @param display_from Optional source ID type for mapping (used when `display_id`
#'   is not found in `rowData`).
#' @param display_orgdb Optional `OrgDb` object used for ID mapping when
#'   `display_id` is set but not found in `rowData`.
#' @return A `ggplot2` object.
#'
#' @examples
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
#' genes <- rownames(vista)[seq_len(3)]
#' get_expression_lollipop(vista, genes = genes)
#' get_expression_lollipop(vista, genes = genes[seq_len(2)], by = "sample")
#'
#' @export
get_expression_lollipop <- function(x,
                                    genes,
                                    sample_group = NULL,
                                    group_column = NULL,
                                    by = c("group", "sample"),
                                    sample_order = c("input", "group", "expression"),
                                    facet_by = c("auto", "gene", "none"),
                                    log_transform = TRUE,
                                    facet_scale = "free_y",
                                    facet_nrow = NULL,
                                    facet_ncol = NULL,
                                    point_size = 6,
                                    line_size = 1.2,
                                    label = TRUE,
                                    label_digits = 1,
                                    display_id = NULL,
                                    display_from = NULL,
                                    display_orgdb = NULL) {
  stopifnot(inherits(x, "VISTA"))
  if (length(genes) > 15) cli::cli_abort("Maximum 15 genes can be plotted at once.")
  by <- match.arg(by)
  sample_order <- match.arg(sample_order)

  meta <- .prepare_sample_metadata(x, sample_group, group_column)
  group_col <- attr(meta, "group_column")
  if (nrow(meta) == 0) {
    cli::cli_abort("No samples found after applying {.arg sample_group}.")
  }

  rd <- tryCatch(SummarizedExperiment::rowData(x), error = function(e) NULL)
  mat <- SummarizedExperiment::assay(x)[, meta$sample, drop = FALSE]
  if (ncol(mat) == 0) {
    cli::cli_abort("No expression columns available after subsetting samples.")
  }

  underlying_ids <- rownames(mat)
  if (!is.null(display_id) && !is.null(rd) && display_id %in% colnames(rd)) {
    map <- rd[[display_id]]
    names(map) <- rownames(x)
    gene_ids <- names(map)[match(genes, map)]
    gene_ids <- gene_ids[!is.na(gene_ids)]
  } else {
    gene_ids <- .map_gene_ids(genes, from_type = display_from, to_type = display_from, orgdb = display_orgdb)
  }
  if (!length(gene_ids)) gene_ids <- genes

  missing_genes <- setdiff(gene_ids, underlying_ids)
  if (length(missing_genes) == length(gene_ids)) cli::cli_abort("None of the specified {.arg genes} were found.")
  keep_genes <- intersect(gene_ids, rownames(mat))
  mat <- mat[keep_genes, , drop = FALSE]

  df <- mat |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    tidyr::pivot_longer(-gene, names_to = "sample", values_to = "expression") |>
    dplyr::left_join(meta, by = "sample")

  if (!is.null(display_id) && !is.null(rd) && display_id %in% colnames(rd)) {
    map <- rd[[display_id]]
    names(map) <- rownames(x)
    mapped <- map[match(df$gene, names(map))]
    df$gene <- ifelse(!is.na(mapped) & nzchar(mapped), mapped, df$gene)
  }

  if (log_transform) {
    df$expression <- log2(df$expression + 1)
    ylab <- "log2(Normalized Counts + 1)"
  } else {
    ylab <- "Normalized Counts"
  }

  if (by == "group") {
    plot_df <- df |>
      dplyr::group_by(.data[[group_col]], .data$gene) |>
      dplyr::summarise(expression = mean(.data$expression, na.rm = TRUE), .groups = "drop")
    if (nrow(plot_df) == 0) cli::cli_abort("No data available to plot after filtering genes and samples.")
    x_var <- group_col
    x_lab <- group_col
  } else {
    plot_df <- .order_expression_plot_samples(df, meta, group_col, sample_order = sample_order)
    if (nrow(plot_df) == 0) cli::cli_abort("No data available to plot after filtering genes and samples.")
    x_var <- "sample"
    x_lab <- "Sample"
  }

  facet_mode <- .resolve_expression_plot_facet(facet_by, n_genes = length(unique(plot_df$gene)))
  cols <- .vista_group_colors(x, plot_df[[group_col]])

  plt <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data[[x_var]], y = expression, color = .data[[group_col]])
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = .data[[x_var]], y = 0, yend = expression),
      linewidth = line_size,
      lineend = "round"
    ) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::labs(x = x_lab, y = ylab, color = group_col) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (label) {
    plt <- plt + ggplot2::geom_text(
      ggplot2::aes(label = round(expression, label_digits)),
      vjust = -0.6,
      size = 4
    )
  }

  if (!is.null(cols)) {
    plt <- plt + ggplot2::scale_color_manual(values = cols)
  }

  if (facet_mode == "gene") {
    plt <- .add_expression_facet_wrap(
      plot = plt,
      facet_formula = ~gene,
      scales = facet_scale,
      facet_nrow = facet_nrow,
      facet_ncol = facet_ncol
    )
  }

  plt
}

# ──────────────────────────────────────────────────────────────────────────────
#' @title Gene expression line plot
#' @description Plots normalized expression for selected genes across samples or
#' summarized groups with optional transformations and group faceting.
#' @param x A `VISTA` object.
#' @param genes Character vector of gene identifiers to plot.
#' @param sample_group Optional character vector specifying which groups (values taken from `group_column`) to include.
#' @param group_column Optional column name in `sample_info` defining the grouping/faceting variable.
#' @param log_transform Logical; log2-transform expression values before plotting.
#' @param display_id Optional ID/column name to use for gene labels. If supplied
#'   and present in `rowData(x)`, those values are used.
#' @param display_from Optional source ID type for mapping (reserved for
#'   compatibility with other expression plotting APIs).
#' @param display_orgdb Optional `OrgDb` object used for ID mapping when
#'   `display_id` is set but not found in `rowData`.
#' @param facet_scales Scaling option passed to `facet_wrap()`.
#' @param facet_nrow,facet_ncol Optional layout passed to `facet_wrap()` when faceting.
#' @param stats_group Logical retained for API consistency. Statistical overlays
#'   are not currently added by `get_expression_lineplot()`.
#' @param p.label Label format retained for API consistency with other
#'   expression plots.
#' @param comparisons Optional list of comparisons retained for API consistency.
#' @param pool_genes Logical; when `TRUE`, average the selected genes into a
#'   single trajectory.
#' @param by Plot unit: `"sample"` (default) or `"group"` to average
#'   replicates before plotting.
#' @param facet_by Faceting mode: `"auto"` (default), `"none"`, `"group"`, or `"gene"`.
#' @param fill_by Argument retained for API consistency; ignored because line
#'   plots use colour rather than fill.
#' @param sample_order Ordering used for sample-level plots: `"input"`,
#'   `"group"`, or `"expression"`.
#' @param value_transform Deprecated compatibility alias for transformation
#'   choice; one of `"log2"`, `"zscore"`, or `"none"`.
#' @param palette Optional qualitative palette name used for gene colours.
#' @param colors Optional named character vector of manual gene colours.
#' @param line_width Line width.
#' @param point_size Point size.
#' @param base_size Base theme size.
#' @aliases get_expression_lineplot
#' @export
get_expression_lineplot <- function(x,
                                    genes = NULL,
                                    sample_group = NULL,
                                    group_column = NULL,
                                    log_transform = TRUE,
                                    display_id = NULL,
                                    display_from = NULL,
                                    display_orgdb = NULL,
                                    facet_scales = "free_y",
                                    facet_nrow = NULL,
                                    facet_ncol = NULL,
                                    stats_group = FALSE,
                                    p.label = "p.signif",
                                    comparisons = NULL,
                                    pool_genes = FALSE,
                                    by = c("sample", "group"),
                                    facet_by = c("auto", "group", "gene", "none"),
                                    fill_by = NULL,
                                    sample_order = c("input", "group", "expression"),
                                    value_transform = NULL,
                                    palette = NULL,
                                    colors = NULL,
                                    line_width = 1,
                                    point_size = 2,
                                    base_size = 12) {
  stopifnot(inherits(x, "VISTA"))
  by <- match.arg(by)
  facet_by <- match.arg(facet_by)
  sample_order <- match.arg(sample_order)

  if (!is.null(value_transform)) {
    value_transform <- match.arg(value_transform, c("log2", "zscore", "none"))
    if (value_transform == "log2") {
      log_transform <- TRUE
    } else if (value_transform == "none") {
      log_transform <- FALSE
    }
  } else {
    value_transform <- if (log_transform) "log2" else "none"
  }
  if (!is.null(fill_by)) {
    cli::cli_warn("{.arg fill_by} is ignored by {.fun get_expression_lineplot}; line plots use colour rather than fill.")
  }
  if (stats_group) {
    cli::cli_warn("{.arg stats_group} is not currently supported for {.fun get_expression_lineplot}; ignoring it.")
  }
  if (!is.null(comparisons)) {
    cli::cli_warn("{.arg comparisons} is ignored by {.fun get_expression_lineplot}.")
  }
  if (!is.null(display_from) || !is.null(display_orgdb)) {
    cli::cli_warn("{.arg display_from} and {.arg display_orgdb} are reserved for API compatibility and are not used when {.arg display_id} maps from {.field rowData(x)}.")
  }

  summarise <- identical(by, "group")
  mat <- SummarizedExperiment::assay(x)
  rd <- tryCatch(SummarizedExperiment::rowData(x), error = function(e) NULL)
  genes_use <- genes
  if (!is.null(genes_use) && !is.null(display_id) && !is.null(rd) && display_id %in% colnames(rd)) {
    map <- as.character(rd[[display_id]])
    names(map) <- rownames(x)
    mapped <- names(map)[match(genes_use, map)]
    mapped <- mapped[!is.na(mapped)]
    if (length(mapped)) genes_use <- mapped
  }
  if (is.null(genes_use)) {
    if (pool_genes) {
      genes_use <- rownames(mat)
      cli::cli_warn("No {.arg genes} provided; pooling all genes in the object.")
    } else {
      gene_vars <- matrixStats::rowVars(mat, na.rm = TRUE)
      if (all(is.na(gene_vars)) || all(gene_vars == 0)) {
        genes_use <- head(rownames(mat), 20)
      } else {
        top_idx <- order(gene_vars, decreasing = TRUE)[seq_len(min(20, length(gene_vars)))]
        genes_use <- rownames(mat)[top_idx]
      }
      cli::cli_warn("No {.arg genes} provided; using the top variable genes for plotting.")
    }
  }

  prep <- .vista_expression_long(
    x = x,
    genes = genes_use,
    sample_group = sample_group,
    group_column = group_column,
    value_transform = value_transform,
    summarise = summarise,
    sample_order = sample_order
  )
  df <- prep$df
  group_col <- prep$group_col
  ylab <- prep$ylab

  if (facet_by == "group" && summarise) {
    cli::cli_warn("{.arg facet_by = 'group'} is ignored when {.arg by = 'group'}.")
    facet_by <- "none"
  }

  if (pool_genes) {
    df <- df |>
      dplyr::group_by(.data$sample, .data[[group_col]]) |>
      dplyr::summarise(expression = mean(.data$expression, na.rm = TRUE), .groups = "drop") |>
      dplyr::mutate(gene = "All genes")
    facet_by <- "none"
  } else if (!is.null(display_id) && !is.null(rd) && display_id %in% colnames(rd)) {
    map <- as.character(rd[[display_id]])
    names(map) <- rownames(x)
    mapped <- map[match(df$gene, names(map))]
    df$gene <- ifelse(!is.na(mapped) & nzchar(mapped), mapped, df$gene)
  }

  facet_mode <- .resolve_expression_plot_facet(facet_by, n_genes = length(unique(df$gene)))
  gene_cols <- .resolve_palette_values(
    levels = unique(df$gene),
    colors = colors,
    palette = palette,
    default = "Dark 3"
  )

  plt <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = sample, y = expression, group = gene, color = gene)
  ) +
    ggplot2::geom_line(linewidth = line_width) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_color_manual(values = gene_cols, drop = FALSE) +
    ggplot2::labs(x = if (summarise) group_col else "Sample", y = ylab, color = "Gene") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (facet_mode == "group" && !summarise && group_col %in% colnames(df)) {
    plt <- .add_expression_facet_wrap(
      plot = plt,
      facet_formula = stats::as.formula(paste("~", group_col)),
      scales = facet_scales,
      facet_nrow = facet_nrow,
      facet_ncol = facet_ncol
    )
  } else if (facet_mode == "gene") {
    plt <- .add_expression_facet_wrap(
      plot = plt,
      facet_formula = ~gene,
      scales = facet_scales,
      facet_nrow = facet_nrow,
      facet_ncol = facet_ncol
    )
  }

  plt
}

# ──────────────────────────────────────────────────────────────────────────────
# Expression distribution plots: violin / density / raincloud
# ──────────────────────────────────────────────────────────────────────────────

.vista_expression_long <- function(x, genes = NULL, sample_group = NULL,
                                   group_column = NULL,
                                   value_transform = c("log2", "zscore", "none"),
                                   summarise = FALSE,
                                   sample_order = c("input", "group", "expression")) {
  value_transform <- match.arg(value_transform)
  meta <- .prepare_sample_metadata(x, sample_group, group_column)
  group_col <- attr(meta, "group_column")
  mat <- SummarizedExperiment::assay(x)[, meta$sample, drop = FALSE]
  if (!is.null(genes)) {
    missing_genes <- setdiff(genes, rownames(mat))
    if (length(missing_genes) == length(genes)) {
      cli::cli_abort("None of the specified {.arg genes} were found in the data.")
    }
    keep_genes <- intersect(genes, rownames(mat))
    mat <- mat[keep_genes, , drop = FALSE]
  }
  df <- mat |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    tidyr::pivot_longer(-gene, names_to = "sample", values_to = "expression") |>
    dplyr::left_join(meta, by = "sample")

  if (value_transform == "log2") {
    df$expression <- log2(df$expression + 1)
    ylab <- "log2(Normalized Counts + 1)"
  } else if (value_transform == "zscore") {
    df <- df |>
      dplyr::group_by(gene) |>
      dplyr::mutate(expression = scale(expression)) |>
      dplyr::ungroup()
    ylab <- "Z-score"
  } else {
    ylab <- "Normalized Counts"
  }

  if (summarise) {
    df <- df |>
      dplyr::group_by(gene, .data[[group_col]]) |>
      dplyr::summarise(expression = mean(expression, na.rm = TRUE), .groups = "drop") |>
      dplyr::mutate(sample = .data[[group_col]])
    df$sample <- factor(df$sample, levels = unique(df$sample))
    return(list(df = df, group_col = group_col, ylab = ylab))
  }

  df <- .order_expression_plot_samples(
    df = df,
    meta = meta,
    group_col = group_col,
    sample_order = sample_order
  )
  list(df = df, group_col = group_col, ylab = ylab)
}

#' @keywords internal
.resolve_expression_plot_facet <- function(facet_by = c("auto", "gene", "group", "none"),
                                           n_genes) {
  if (length(facet_by) != 1L) {
    facet_by <- "auto"
  }
  facet_by <- match.arg(facet_by)
  if (facet_by == "auto") {
    return(if (n_genes > 1L) "gene" else "none")
  }
  facet_by
}

#' @keywords internal
.order_expression_plot_samples <- function(df,
                                           meta,
                                           group_col,
                                           sample_order = c("input", "group", "expression")) {
  sample_order <- match.arg(sample_order)
  meta_ord <- meta
  meta_ord$sample <- as.character(meta_ord$sample)

  if (sample_order == "expression") {
    expr_ord <- df |>
      dplyr::group_by(.data$sample) |>
      dplyr::summarise(.expr_order = mean(.data$expression, na.rm = TRUE), .groups = "drop")
    meta_ord <- dplyr::left_join(meta_ord, expr_ord, by = "sample")
    meta_ord <- dplyr::arrange(meta_ord, dplyr::desc(.data$.expr_order), .data$sample)
  } else if (sample_order == "group") {
    meta_ord <- dplyr::arrange(
      meta_ord,
      factor(.data[[group_col]], levels = unique(meta[[group_col]])),
      .data$sample
    )
  }

  sample_levels <- unique(meta_ord$sample)
  df$sample <- factor(df$sample, levels = sample_levels)
  df
}

#' @keywords internal
.resolve_foldchange_facet_mode <- function(facet_by = c("auto", "gene", "comparison", "none"),
                                           n_genes = 1L) {
  facet_by <- match.arg(facet_by)
  if (facet_by == "auto") {
    return(if (n_genes > 1L) "gene" else "none")
  }
  facet_by
}

.resolve_foldchange_gene_ids <- function(x,
                                         genes = NULL,
                                         display_id = NULL,
                                         display_from = NULL,
                                         display_orgdb = NULL) {
  if (is.null(genes)) {
    return(NULL)
  }

  genes <- unique(stats::na.omit(as.character(genes)))
  genes <- genes[nzchar(genes)]
  if (!length(genes)) {
    return(character())
  }

  rd <- tryCatch(SummarizedExperiment::rowData(x), error = function(e) NULL)

  if (!is.null(display_id) && !is.null(rd) && display_id %in% colnames(rd)) {
    map <- as.character(rd[[display_id]])
    names(map) <- rownames(x)
    mapped <- names(map)[match(genes, map)]
    mapped <- mapped[!is.na(mapped) & nzchar(mapped)]
    if (length(mapped)) {
      return(unique(mapped))
    }
  } else if (!is.null(display_id) && !is.null(display_from) && !is.null(display_orgdb)) {
    mapped <- .map_gene_ids(
      ids = genes,
      from_type = display_id,
      to_type = display_from,
      orgdb = display_orgdb
    )
    mapped <- unique(stats::na.omit(as.character(mapped)))
    mapped <- mapped[nzchar(mapped)]
    if (length(mapped)) {
      return(mapped)
    }
  }

  genes
}

.resolve_foldchange_gene_labels <- function(x,
                                            gene_ids,
                                            display_id = NULL,
                                            display_from = NULL,
                                            display_orgdb = NULL) {
  .resolve_heatmap_row_labels(
    x = x,
    gene_ids = gene_ids,
    display_id = display_id,
    display_from = display_from,
    display_orgdb = display_orgdb,
    repair_genes = FALSE,
    prefer_symbol = FALSE
  )
}

#' Violin plot of expression values
#'
#' Mirrors the main user-facing arguments of [get_expression_boxplot()] so the
#' two geoms can be swapped with minimal code changes. Violin plots currently
#' keep group-based semantics (`by = "group"`) because a violin requires
#' replicate-level distributions within groups rather than one value per sample.
#'
#' @param x A `VISTA` object.
#' @param genes Optional character vector of gene IDs to include; defaults to
#'   all genes selected by the plotting mode.
#' @param sample_group Optional subset of groups (values of `group_column`) to keep.
#' @param group_column Grouping column in `sample_info`; defaults to the stored grouping.
#' @param log_transform Logical; log2-transform expression values before plotting.
#' @param display_id Optional ID/column name to use for labels/facets. If
#'   supplied and present in `rowData(x)`, those values are used.
#' @param display_from Optional source ID type for mapping (reserved for
#'   compatibility with [get_expression_boxplot()]).
#' @param display_orgdb Optional `OrgDb` object used for ID mapping when
#'   `display_id` is set but not found in `rowData` (reserved for compatibility
#'   with [get_expression_boxplot()]).
#' @param facet_scales Scaling option passed to `facet_wrap()`.
#' @param facet_nrow,facet_ncol Optional layout passed to `facet_wrap()` when faceting.
#' @param stats_group Logical; add statistical comparisons between groups when `TRUE`.
#' @param p.label Label format for `ggpubr::stat_compare_means()`.
#' @param comparisons Optional list of specific group comparisons for `stat_compare_means()`.
#' @param pool_genes Logical; pool all selected genes into one violin per group.
#' @param by Plot unit. Violin plots currently support only `"group"`.
#' @param facet_by Faceting mode. Uses the same argument pattern as
#'   [get_expression_boxplot()], but `pool_genes = TRUE` falls back to `"none"`
#'   because pooled violins already aggregate across genes.
#' @param fill_by Fill mapping. Uses the same values as
#'   [get_expression_boxplot()], including discrete sample metadata columns.
#' @param sample_order Ordering for sample-level display before values are
#'   grouped into violins.
#' @param value_transform Deprecated compatibility alias. `"log2"` maps to
#'   `log_transform = TRUE`, `"none"` maps to `FALSE`, and `"zscore"` applies a
#'   per-gene z-score transform.
#' @param summarise Logical retained for compatibility. Violin plots always use
#'   replicate-level values, so `summarise = TRUE` is ignored with a warning.
#'
#' @return A `ggplot2` object.
#' @export
get_expression_violinplot <- function(x,
                                      genes = NULL,
                                      sample_group = NULL,
                                      group_column = NULL,
                                      log_transform = TRUE,
                                      display_id = NULL,
                                      display_from = NULL,
                                      display_orgdb = NULL,
                                      facet_scales = "free_y",
                                      facet_nrow = NULL,
                                      facet_ncol = NULL,
                                      stats_group = FALSE,
                                      p.label = "p.signif",
                                      comparisons = NULL,
                                      pool_genes = FALSE,
                                      by = "group",
                                      facet_by = c("auto", "gene", "none"),
                                      fill_by = NULL,
                                      sample_order = c("input", "group", "expression"),
                                      value_transform = NULL,
                                      summarise = FALSE) {
  stopifnot(inherits(x, "VISTA"))
  if (!identical(by, "group")) {
    cli::cli_abort(
      "{.fun get_expression_violinplot} currently supports only {.code by = 'group'} because violins require replicate distributions within groups."
    )
  }
  sample_order <- match.arg(sample_order)
  if (!is.null(value_transform)) {
    value_transform <- match.arg(value_transform, c("log2", "zscore", "none"))
    if (value_transform == "log2") {
      log_transform <- TRUE
    } else if (value_transform == "none") {
      log_transform <- FALSE
    }
  }
  if (isTRUE(summarise)) {
    cli::cli_warn(
      "{.arg summarise = TRUE} is not meaningful for violin plots because each group becomes a single value. Using replicate-level values instead."
    )
    summarise <- FALSE
  }
  if (pool_genes && !identical(facet_by, "none") && !identical(facet_by, "auto")) {
    cli::cli_warn(
      "{.arg facet_by} is ignored when {.arg pool_genes = TRUE}; pooled violin plots use a single panel."
    )
  }

  mat <- SummarizedExperiment::assay(x)
  meta <- .prepare_sample_metadata(x, sample_group, group_column)
  group_col <- attr(meta, "group_column")
  mat <- mat[, meta$sample, drop = FALSE]

  rd <- tryCatch(SummarizedExperiment::rowData(x), error = function(e) NULL)

  if (!is.null(genes)) {
    gene_ids <- genes
    if (!is.null(display_id) && !is.null(rd) && display_id %in% colnames(rd)) {
      map <- rd[[display_id]]
      names(map) <- rownames(x)
      mapped <- names(map)[match(genes, map)]
      mapped <- mapped[!is.na(mapped)]
      if (length(mapped)) gene_ids <- mapped
    }
    missing_genes <- setdiff(gene_ids, rownames(mat))
    if (length(missing_genes) == length(gene_ids)) {
      cli::cli_abort("None of the specified {.arg genes} were found in the data.")
    }
    keep_genes <- intersect(gene_ids, rownames(mat))
    mat <- mat[keep_genes, , drop = FALSE]
  } else if (!pool_genes && identical(facet_by, "gene")) {
    gene_vars <- matrixStats::rowVars(mat)
    if (all(is.na(gene_vars)) || all(gene_vars == 0)) {
      cli::cli_warn("All genes have zero variance; defaulting to the first 20 genes for faceting.")
      genes <- head(rownames(mat), 20)
    } else {
      top_idx <- order(gene_vars, decreasing = TRUE)[seq_len(min(20, length(gene_vars)))]
      genes <- rownames(mat)[top_idx]
      cli::cli_warn("No {.arg genes} provided; faceting by top 20 variable genes.")
    }
    mat <- mat[rownames(mat) %in% genes, , drop = FALSE]
  }

  df <- mat |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    tibble::as_tibble() |>
    tidyr::pivot_longer(-gene, names_to = "sample", values_to = "expression") |>
    dplyr::left_join(meta, by = "sample")

  if (!is.null(display_id) && !is.null(rd) && display_id %in% colnames(rd)) {
    map <- rd[[display_id]]
    names(map) <- rownames(x)
    mapped <- map[match(df$gene, names(map))]
    df$gene <- ifelse(!is.na(mapped) & nzchar(mapped), mapped, df$gene)
  }

  if (!is.null(value_transform) && identical(value_transform, "zscore")) {
    df <- df |>
      dplyr::group_by(.data$gene) |>
      dplyr::mutate(expression = as.numeric(scale(.data$expression))) |>
      dplyr::ungroup()
    ylab <- "Z-score"
  } else if (log_transform) {
    df$expression <- log2(df$expression + 1)
    ylab <- "log2(Normalized Counts + 1)"
  } else {
    ylab <- "Normalized Counts"
  }

  if (pool_genes) {
    df$gene <- "All genes"
    facet_by <- "none"
  }

  facet_mode <- if (pool_genes) {
    "none"
  } else {
    .resolve_expression_plot_facet(facet_by, n_genes = length(unique(df$gene)))
  }
  fill_spec <- .resolve_expression_fill(
    x = x,
    df = df,
    group_col = group_col,
    fill_by = fill_by,
    default_fill = if (pool_genes) "x" else "group",
    allowed_special = if (pool_genes) c("x", "group", "gene") else c("group", "gene"),
    x_var = group_col,
    x_label = group_col
  )

  group_levels <- unique(as.character(df[[group_col]]))
  df[[group_col]] <- factor(as.character(df[[group_col]]), levels = group_levels)

  plt <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[group_col]],
      y = expression,
      fill = factor(fill_spec$values, levels = fill_spec$levels)
    )
  ) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.8, color = NA, scale = "width") +
    ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6, color = "gray30") +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.08, height = 0),
      size = 1.3,
      alpha = 0.65,
      shape = 21,
      stroke = 0.2,
      color = "gray20"
    ) +
    ggplot2::labs(x = group_col, y = ylab, fill = fill_spec$label) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))

  if (!is.null(fill_spec$palette)) {
    plt <- plt + ggplot2::scale_fill_manual(values = fill_spec$palette)
  }

  if (stats_group) {
    if (!requireNamespace("ggpubr", quietly = TRUE)) {
      cli::cli_abort("Package {.pkg ggpubr} must be installed for `stats_group = TRUE`.")
    }
    if (length(unique(df[[group_col]])) < 2) {
      cli::cli_warn("At least two groups are needed for statistical testing.")
    } else {
      plt <- plt + ggpubr::stat_compare_means(
        ggplot2::aes(group = .data[[group_col]]),
        comparisons = comparisons,
        method = "t.test",
        label = p.label,
        label.x.npc = "center"
      )
    }
  }

  if (facet_mode != "none") {
    plt <- .add_expression_facet_wrap(
      plot = plt,
      facet_formula = ~gene,
      scales = facet_scales,
      facet_nrow = facet_nrow,
      facet_ncol = facet_ncol
    )
  }
  plt
}

#' Raincloud plot of expression values
#'
#' Uses `ggrain::geom_rain()` to combine a half-violin, boxplot, and jittered
#' points per sample/group to show distribution, summary, and individual values.
#'
#' @inheritParams get_expression_violinplot
#' @param summarise Logical; when `TRUE`, collapse replicates within each
#'   group for each gene (one value per gene per group). This is useful for
#'   pooled multi-gene raincloud plots where each dot represents one gene-level
#'   summary in a group.
#' @param rain_side Side specification passed to `ggrain::geom_rain()`; one of
#'   `"r"`, `"l"`, `"f"`, `"f1x1"`, or `"f2x2"`.
#' @param id.long.var Optional column name passed to `ggrain::geom_rain()` as
#'   `id.long.var` to identify repeated measurements.
#' @param alpha Alpha for jittered points.
#' @param point_size Point size for jittered points.
#' @param p.label Label type passed to `ggpubr::stat_compare_means()`.
#' @param stats_group Logical; add pairwise statistical tests when `TRUE`.
#' @param stats_method Statistical method passed to `ggpubr::stat_compare_means()`.
#' @param label Logical; add text labels to points using `ggrepel`.
#' @param label_column Column name in the plotting data used for labels.
#'   Defaults to `"gene"` for expression raincloud plots.
#' @param label_size Text size for point labels.
#' @param label_max_overlaps Maximum overlaps passed to `ggrepel::geom_text_repel()`.
#' @param display_id Optional ID/column name to use for labels. If supplied and
#'   present in `rowData(x)`, those values are used; otherwise falls back to ID mapping.
#' @param display_from Optional source ID type for mapping (used when
#'   `display_id` is not found in `rowData`).
#' @param display_orgdb Optional `OrgDb` object used for ID mapping when
#'   `display_id` is set but not found in `rowData`.
#'
#' @details
#' `id.long.var` controls which repeated unit is connected by lines in
#' `ggrain::geom_rain()`.
#'
#' Recommended usage for expression raincloud plots:
#' \itemize{
#'   \item `id.long.var = NULL` (default): best for clean distribution summaries.
#'   \item `id.long.var = "gene"`: best when plotting a small number of genes and
#'   showing gene-level trajectories across x levels.
#'   \item `id.long.var = "<subject_id_column>"`: best for paired/repeated-measure
#'   designs when a subject ID exists in `sample_info`.
#'   \item `id.long.var = "sample"` or the grouping variable is usually less
#'   informative and can over-connect points.
#'   \item Point labels (`label = TRUE`) work best with `facet_by = "none"` or
#'   a small number of genes.
#' }
#'
#' For identifier display consistency with other VISTA plotting functions, set
#' `display_id` (for example, `"SYMBOL"`). When provided, `genes` can be given
#' in that ID space, and default point labels use the mapped display IDs.
#'
#' @return A `ggplot2` object.
#' @export
get_expression_raincloud <- function(x,
                                     genes = NULL,
                                     sample_group = NULL,
                                     group_column = NULL,
                                     by = "group",
                                     value_transform = c("log2", "zscore", "none"),
                                     summarise = FALSE,
                                     facet_by = c("auto", "gene", "none"),
                                     fill_by = NULL,
                                     facet_nrow = NULL,
                                     facet_ncol = NULL,
                                     sample_order = c("input", "group", "expression"),
                                     rain_side = c("r", "l", "f", "f1x1", "f2x2"),
                                     id.long.var = NULL,
                                     alpha = 0.5,
                                     point_size = 1.5,
                                     p.label = "p.signif",
                                     stats_group = FALSE,
                                     stats_method = "t.test",
                                     label = FALSE,
                                     label_column = "gene",
                                     label_size = 3,
                                     label_max_overlaps = 50,
                                     display_id = NULL,
                                     display_from = NULL,
                                     display_orgdb = NULL) {
  stopifnot(inherits(x, "VISTA"))
  if (!requireNamespace("ggrain", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggrain} must be installed to use `get_expression_raincloud()`.")
  }
  if (!identical(by, "group")) {
    cli::cli_abort(
      "{.fun get_expression_raincloud} currently supports only {.code by = 'group'} because raincloud plots require replicate distributions within groups."
    )
  }
  rain_side <- match.arg(rain_side)
  rd <- tryCatch(SummarizedExperiment::rowData(x), error = function(e) NULL)
  sample_order <- match.arg(sample_order)

  genes_use <- genes
  if (!is.null(genes) && !is.null(display_id)) {
    if (!is.null(rd) && display_id %in% colnames(rd)) {
      map <- rd[[display_id]]
      names(map) <- rownames(x)
      mapped <- names(map)[match(genes, map)]
      mapped <- unique(mapped[!is.na(mapped)])
      if (length(mapped)) genes_use <- mapped
    } else if (!is.null(display_from) && !is.null(display_orgdb)) {
      genes_use <- .map_gene_ids(
        ids = genes,
        from_type = display_id,
        to_type = display_from,
        orgdb = display_orgdb
      )
    }
  }

  prep <- .vista_expression_long(
    x, genes_use, sample_group, group_column, value_transform, summarise,
    sample_order = sample_order
  )
  df <- prep$df; group_col <- prep$group_col; ylab <- prep$ylab
  facet_mode <- .resolve_expression_plot_facet(facet_by, n_genes = length(unique(df$gene)))
  group_levels <- unique(as.character(df[[group_col]]))
  df[[group_col]] <- factor(as.character(df[[group_col]]), levels = group_levels)

  if (isTRUE(summarise) && length(unique(df$gene)) < 2L) {
    cli::cli_warn(
      "{.arg summarise = TRUE} with a single gene yields one value per group (no within-group distribution)."
    )
  }
  if (isTRUE(summarise) && identical(facet_mode, "gene")) {
    cli::cli_warn(
      "{.arg summarise = TRUE} with {.arg facet_by = 'gene'} gives one summarized value per group in each facet, so raincloud distributions are not informative. Consider {.code facet_by = 'none'} with {.code id.long.var = 'gene'}."
    )
  }

  if (!is.null(display_id)) {
    if (!is.null(rd) && display_id %in% colnames(rd)) {
      map <- rd[[display_id]]
      names(map) <- rownames(x)
      mapped <- map[match(df$gene, names(map))]
      df$gene_label <- ifelse(!is.na(mapped) & nzchar(mapped), mapped, df$gene)
    } else {
      df$gene_label <- .map_gene_ids(
        ids = df$gene,
        from_type = display_from,
        to_type = display_id,
        orgdb = display_orgdb
      )
    }
  } else {
    df$gene_label <- df$gene
  }

  label_column_effective <- label_column
  if (!is.null(display_id) && identical(label_column, "gene")) {
    label_column_effective <- "gene_label"
  }

  if (!is.null(id.long.var) && !id.long.var %in% colnames(df)) {
    cli::cli_abort("{.arg id.long.var} must be a column in the plotting data. Available columns: {.val {colnames(df)}}.")
  }
  if (!label_column_effective %in% colnames(df)) {
    cli::cli_abort("{.arg label_column} must be a column in the plotting data. Available columns: {.val {colnames(df)}}.")
  }
  fill_spec <- .resolve_expression_fill(
    x = x,
    df = df,
    group_col = group_col,
    fill_by = fill_by,
    default_fill = "group",
    allowed_special = c("group", "gene", "x"),
    x_var = group_col,
    x_label = group_col
  )

  rain_args <- list(
    rain.side = rain_side,
    point.args = list(alpha = alpha, size = point_size),
    boxplot.args = list(outlier.shape = NA, alpha = 0.7, color = "gray30"),
    violin.args = list(alpha = 0.6, color = NA)
  )
  if (!is.null(id.long.var)) {
    rain_args$id.long.var <- id.long.var
  }

  plt <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[group_col]],
      y = expression,
      fill = factor(fill_spec$values, levels = fill_spec$levels)
    )
  ) +
    do.call(ggrain::geom_rain, rain_args) +
    ggplot2::labs(x = group_col, y = ylab, fill = fill_spec$label) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))
  if (!is.null(fill_spec$palette)) {
    plt <- plt + ggplot2::scale_fill_manual(values = fill_spec$palette)
  }
  if (facet_mode != "none") {
    plt <- .add_expression_facet_wrap(
      plot = plt,
      facet_formula = ~gene,
      scales = "free_y",
      facet_nrow = facet_nrow,
      facet_ncol = facet_ncol
    )
  }

  if (label) {
    if (facet_mode != "none" && nrow(df) > 120) {
      cli::cli_warn("{.arg label = TRUE} with faceting can create crowded labels. Consider fewer genes or {.arg facet_by = 'none'}.")
    }
    label_cols <- unique(c(group_col, "expression", label_column_effective, "gene"))
    label_df <- unique(df[, label_cols, drop = FALSE])
    plt <- plt + ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(label = .data[[label_column_effective]]),
      size = label_size,
      vjust = -0.3,
      check_overlap = TRUE,
      show.legend = FALSE
    )
  }

  if (stats_group) {
    if (facet_mode != "none") {
      cli::cli_warn("Statistical comparison is only supported when {.arg facet_by = 'none'}. Skipping stats.")
    } else if (length(unique(df[[group_col]])) < 2) {
      cli::cli_warn("At least two groups are required for statistical testing.")
    } else {
      plt <- plt + ggpubr::stat_compare_means(
        ggplot2::aes(group = .data[[group_col]]),
        method = stats_method,
        label = p.label
      )
    }
  }

  plt
}

# ──────────────────────────────────────────────────────────────────────────────
# log2FC boxplot (per comparison)
# ──────────────────────────────────────────────────────────────────────────────

#' Plot fold-change distributions across comparisons
#'
#' Builds boxplots of log2 fold changes for selected genes and comparisons,
#' optionally adding statistics.
#'
#' @param x A `VISTA` object containing differential expression results.
#' @param genes Optional character vector of gene IDs to include.
#' @param sample_comparisons Optional character vector of comparison names to plot.
#' @param display_id Optional ID/column name used to interpret `genes` and,
#'   when possible, map fold-change gene identifiers to display-friendly labels.
#' @param display_from Optional source ID type for mapping when `display_id` is
#'   not present in `rowData(x)`.
#' @param display_orgdb Optional `OrgDb` object used for identifier mapping when
#'   `display_id` is not present in `rowData(x)`.
#' @param facet_scales Facet scales argument passed to `facet_wrap()` when
#'   `facet_by != "none"` (default `"free_x"`).
#' @param facet_nrow,facet_ncol Optional layout passed to `facet_wrap()` when faceting.
#' @param facet_by Faceting mode: `"auto"` (default), `"comparison"`, or `"none"`.
#' @param p.label Label type passed to `ggpubr::stat_compare_means()`.
#' @param stats_group Logical; add pairwise statistical tests when `TRUE`.
#' @param stats_method Statistical method passed to `ggpubr::stat_compare_means()`.
#'
#' @export
#' @aliases get_foldchange_boxplot
get_foldchange_boxplot <- function(x,
                                   genes = NULL,
                                   sample_comparisons = NULL,
                                   display_id = NULL,
                                   display_from = NULL,
                                   display_orgdb = NULL,
                                   facet_scales = "free_x",
                                   facet_nrow = NULL,
                                   facet_ncol = NULL,
                                   facet_by = c("auto", "comparison", "none"),
                                   p.label = "p.signif",
                                   stats_group = FALSE,
                                   stats_method = "t.test") {
  stopifnot(inherits(x, "VISTA"))

  # Build fold-change matrix from stored DE results
  comps <- .vista_comparisons(x)
  if (is.null(sample_comparisons)) {
    sample_comparisons <- names(comps)
  }
  stopifnot(all(sample_comparisons %in% names(comps)))
  facet_mode <- {
    facet_by <- match.arg(facet_by)
    if (identical(facet_by, "auto")) {
      if (length(sample_comparisons) > 1L) "comparison" else "none"
    } else {
      facet_by
    }
  }

  # wide matrix gene_id x comparisons
  all_genes <- unique(unlist(lapply(comps[sample_comparisons], \(d) d$gene_id)))
  genes_use <- .resolve_foldchange_gene_ids(
    x = x,
    genes = genes,
    display_id = display_id,
    display_from = display_from,
    display_orgdb = display_orgdb
  )
  if (!is.null(genes_use)) all_genes <- intersect(all_genes, genes_use)
  if (!length(all_genes)) cli::cli_abort("No genes found for the requested comparisons.")

  fc_mat <- vapply(sample_comparisons, function(nm) {
    dd <- comps[[nm]][, c("gene_id", "log2fc")]
    vv <- dd$log2fc[match(all_genes, dd$gene_id)]
    vv
  }, FUN.VALUE = numeric(length(all_genes)))
  colnames(fc_mat) <- sample_comparisons
  df <- tibble::tibble(gene_id = all_genes) |>
    dplyr::bind_cols(as.data.frame(fc_mat)) |>
    tidyr::pivot_longer(-gene_id, names_to = "comparison", values_to = "log2FoldChange") |>
    dplyr::filter(!is.na(log2FoldChange))

  df$gene_label <- .resolve_foldchange_gene_labels(
    x = x,
    gene_ids = df$gene_id,
    display_id = display_id,
    display_from = display_from,
    display_orgdb = display_orgdb
  )
  df$comparison <- factor(df$comparison, levels = sample_comparisons)

  pal <- .vista_comparison_colors(x, sample_comparisons)
  if (is.null(pal)) {
    pal <- colorspace::qualitative_hcl(length(sample_comparisons), palette = "Dark 3")
  }
  if (is.null(names(pal))) names(pal) <- sample_comparisons

  plt <- ggplot2::ggplot(df, ggplot2::aes(x = comparison, y = log2FoldChange, fill = comparison)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::labs(x = NULL, y = "log2 Fold Change", fill = "Comparison") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (facet_mode == "comparison") {
    plt <- .add_expression_facet_wrap(
      plot = plt,
      facet_formula = ~comparison,
      scales = facet_scales,
      facet_nrow = facet_nrow,
      facet_ncol = facet_ncol
    )
  }

  if (stats_group) {
    if (facet_mode != "none") {
      cli::cli_warn("Statistical comparison is only supported when {.arg facet_by = 'none'}. Skipping stats.")
    } else if (length(unique(df$comparison)) < 2) {
      cli::cli_warn("At least two comparisons are required for statistical testing.")
    } else {
      plt <- plt + ggpubr::stat_compare_means(
        ggplot2::aes(group = comparison),
        method = stats_method,
        label = p.label
      )
    }
  }

  plt
}

#' Raincloud plot of fold-change distributions across comparisons
#'
#' Uses `ggrain::geom_rain()` to display log2 fold-change distributions for
#' selected comparisons, with optional statistical testing across comparisons.
#'
#' @param x A `VISTA` object containing differential expression results.
#' @param genes Optional character vector of gene IDs to include.
#' @param sample_comparisons Optional character vector of comparison names to plot.
#' @param facet_scales Facet scales argument passed to `facet_wrap()` when
#'   `facet_by != "none"` (default `"free_x"`).
#' @param facet_nrow,facet_ncol Optional layout passed to `facet_wrap()` when faceting.
#' @param facet_by Faceting mode: `"auto"` (default), `"comparison"`, or `"none"`.
#' @param rain_side Side specification passed to `ggrain::geom_rain()`; one of
#'   `"r"`, `"l"`, `"f"`, `"f1x1"`, or `"f2x2"`.
#' @param id.long.var Optional column name passed to `ggrain::geom_rain()` as
#'   `id.long.var` to identify repeated measurements.
#' @param alpha Alpha for jittered points.
#' @param point_size Point size for jittered points.
#' @param p.label Label type passed to `ggpubr::stat_compare_means()`.
#' @param stats_group Logical; add pairwise statistical tests when `TRUE`.
#' @param stats_method Statistical method passed to `ggpubr::stat_compare_means()`.
#' @param label Logical; add text labels to points using `ggrepel`.
#' @param label_column Column name in the plotting data used for labels.
#'   Defaults to `"gene_id"` for fold-change raincloud plots.
#' @param label_size Text size for point labels.
#' @param label_max_overlaps Maximum overlaps passed to `ggrepel::geom_text_repel()`.
#' @param display_id Optional ID/column name to use for labels. If supplied and
#'   present in `rowData(x)`, those values are used; otherwise falls back to ID mapping.
#' @param display_from Optional source ID type for mapping (used when
#'   `display_id` is not found in `rowData`).
#' @param display_orgdb Optional `OrgDb` object used for ID mapping when
#'   `display_id` is set but not found in `rowData`.
#'
#' @details
#' `id.long.var` controls which repeated unit is connected by lines in
#' `ggrain::geom_rain()`.
#'
#' Recommended usage for fold-change raincloud plots:
#' \itemize{
#'   \item `id.long.var = NULL` (default): best for clean distribution summaries.
#'   \item `id.long.var = "gene_id"`: most useful option; connects each gene
#'   across comparisons.
#'   \item `id.long.var = "comparison"` is generally not useful because
#'   comparison is already on the x-axis.
#'   \item Continuous value columns (e.g. `log2FoldChange`) are not suitable
#'   identifiers for line connections.
#'   \item Point labels (`label = TRUE`) work best with `facet_by = "none"`
#'   unless only a small set of genes is shown.
#' }
#'
#' For identifier display consistency with other VISTA plotting functions, set
#' `display_id` (for example, `"SYMBOL"`). When provided, `genes` can be given
#' in that ID space, and default point labels use the mapped display IDs.
#'
#' @return A `ggplot2` object.
#' @export
get_foldchange_raincloud <- function(x,
                                     genes = NULL,
                                     sample_comparisons = NULL,
                                     facet_scales = "free_x",
                                     facet_nrow = NULL,
                                     facet_ncol = NULL,
                                     facet_by = c("auto", "comparison", "none"),
                                     rain_side = c("r", "l", "f", "f1x1", "f2x2"),
                                     id.long.var = NULL,
                                     alpha = 0.5,
                                     point_size = 1.5,
                                     p.label = "p.signif",
                                     stats_group = FALSE,
                                     stats_method = "t.test",
                                     label = FALSE,
                                     label_column = "gene_id",
                                     label_size = 3,
                                     label_max_overlaps = 50,
                                     display_id = NULL,
                                     display_from = NULL,
                                     display_orgdb = NULL) {
  stopifnot(inherits(x, "VISTA"))
  if (!requireNamespace("ggrain", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggrain} must be installed to use `get_foldchange_raincloud()`.")
  }
  rain_side <- match.arg(rain_side)

  comps <- .vista_comparisons(x)
  rd <- tryCatch(SummarizedExperiment::rowData(x), error = function(e) NULL)
  if (is.null(sample_comparisons)) {
    sample_comparisons <- names(comps)
  }
  stopifnot(all(sample_comparisons %in% names(comps)))
  facet_mode <- {
    facet_by <- match.arg(facet_by)
    if (identical(facet_by, "auto")) {
      if (length(sample_comparisons) > 1L) "comparison" else "none"
    } else {
      facet_by
    }
  }

  all_genes <- unique(unlist(lapply(comps[sample_comparisons], \(d) d$gene_id)))
  genes_use <- genes
  if (!is.null(genes) && !is.null(display_id)) {
    if (!is.null(rd) && display_id %in% colnames(rd)) {
      map <- rd[[display_id]]
      names(map) <- rownames(x)
      mapped <- names(map)[match(genes, map)]
      mapped <- unique(mapped[!is.na(mapped)])
      if (length(mapped)) genes_use <- mapped
    } else if (!is.null(display_from) && !is.null(display_orgdb)) {
      genes_use <- .map_gene_ids(
        ids = genes,
        from_type = display_id,
        to_type = display_from,
        orgdb = display_orgdb
      )
    }
  }
  if (!is.null(genes_use)) all_genes <- intersect(all_genes, genes_use)
  if (!length(all_genes)) cli::cli_abort("No genes found for the requested comparisons.")

  fc_mat <- vapply(sample_comparisons, function(nm) {
    dd <- comps[[nm]][, c("gene_id", "log2fc")]
    dd$log2fc[match(all_genes, dd$gene_id)]
  }, FUN.VALUE = numeric(length(all_genes)))
  colnames(fc_mat) <- sample_comparisons
  df <- tibble::tibble(gene_id = all_genes) |>
    dplyr::bind_cols(as.data.frame(fc_mat)) |>
    tidyr::pivot_longer(-gene_id, names_to = "comparison", values_to = "log2FoldChange") |>
    dplyr::filter(!is.na(log2FoldChange))

  if (!is.null(display_id)) {
    if (!is.null(rd) && display_id %in% colnames(rd)) {
      map <- rd[[display_id]]
      names(map) <- rownames(x)
      mapped <- map[match(df$gene_id, names(map))]
      df$gene_label <- ifelse(!is.na(mapped) & nzchar(mapped), mapped, df$gene_id)
    } else {
      df$gene_label <- .map_gene_ids(
        ids = df$gene_id,
        from_type = display_from,
        to_type = display_id,
        orgdb = display_orgdb
      )
    }
  } else {
    df$gene_label <- df$gene_id
  }

  label_column_effective <- label_column
  if (!is.null(display_id) && identical(label_column, "gene_id")) {
    label_column_effective <- "gene_label"
  }

  if (!is.null(id.long.var) && !id.long.var %in% colnames(df)) {
    cli::cli_abort("{.arg id.long.var} must be a column in the plotting data. Available columns: {.val {colnames(df)}}.")
  }
  if (!label_column_effective %in% colnames(df)) {
    cli::cli_abort("{.arg label_column} must be a column in the plotting data. Available columns: {.val {colnames(df)}}.")
  }

  df$comparison <- factor(df$comparison, levels = sample_comparisons)

  pal <- .vista_comparison_colors(x, sample_comparisons)
  if (is.null(pal)) {
    pal <- colorspace::qualitative_hcl(length(sample_comparisons), palette = "Dark 3")
  }
  if (is.null(names(pal))) names(pal) <- sample_comparisons

  rain_args <- list(
    rain.side = rain_side,
    point.args = list(alpha = alpha, size = point_size),
    boxplot.args = list(outlier.shape = NA, alpha = 0.8),
    violin.args = list(alpha = 0.6, color = NA)
  )
  if (!is.null(id.long.var)) {
    rain_args$id.long.var <- id.long.var
  }

  plt <- ggplot2::ggplot(df, ggplot2::aes(x = comparison, y = log2FoldChange, fill = comparison)) +
    do.call(ggrain::geom_rain, rain_args) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::labs(x = NULL, y = "log2 Fold Change", fill = "Comparison") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (facet_mode == "comparison") {
    plt <- .add_expression_facet_wrap(
      plot = plt,
      facet_formula = ~comparison,
      scales = facet_scales,
      facet_nrow = facet_nrow,
      facet_ncol = facet_ncol
    )
  }

  if (label) {
    if (!requireNamespace("ggrepel", quietly = TRUE)) {
      cli::cli_abort("Package {.pkg ggrepel} must be installed when {.arg label = TRUE}.")
    }
    if (facet_mode != "none" && nrow(df) > 120) {
      cli::cli_warn("{.arg label = TRUE} with faceting can create crowded labels. Consider fewer genes or {.arg facet_by = 'none'}.")
    }
    label_cols <- unique(c("comparison", "log2FoldChange", label_column_effective))
    label_df <- unique(df[, label_cols, drop = FALSE])
    plt <- plt + ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(label = .data[[label_column_effective]]),
      size = label_size,
      max.overlaps = label_max_overlaps,
      show.legend = FALSE
    )
  }

  if (stats_group) {
    if (facet_mode != "none") {
      cli::cli_warn("Statistical comparison is only supported when {.arg facet_by = 'none'}. Skipping stats.")
    } else if (length(unique(df$comparison)) < 2) {
      cli::cli_warn("At least two comparisons are required for statistical testing.")
    } else {
      plt <- plt + ggpubr::stat_compare_means(
        ggplot2::aes(group = comparison),
        method = stats_method,
        label = p.label
      )
    }
  }

  plt
}

# ──────────────────────────────────────────────────────────────────────────────
# Fold-change scatter between two comparisons
# ──────────────────────────────────────────────────────────────────────────────

#' Fold-change scatterplot between two comparisons
#'
#' Plots log2 fold changes from two stored comparisons against each other, with
#' points coloured by concordant/discordant regulation based on the cutoffs
#' saved in the `VISTA` object.
#'
#' @param x A `VISTA` object containing DE results.
#' @param sample_comparisons Character vector of length 2 naming the comparisons.
#' @param genes Optional character vector of gene identifiers used to subset the
#'   comparison overlap before plotting.
#' @param display_id Optional ID/column name used to interpret `genes` and to
#'   label highlighted points.
#' @param display_from Optional source ID type for mapping when `display_id` is
#'   not present in `rowData(x)`.
#' @param display_orgdb Optional `OrgDb` object used for identifier mapping when
#'   `display_id` is not present in `rowData(x)`.
#' @param label_n Integer; number of most extreme points to label (by |log2FC1| + |log2FC2|).
#' @param alpha Point transparency.
#' @param geometry Geometry used for the data layer: `"point"` or `"hex"`.
#' @param method Correlation method for the subtitle; `"pearson"` or `"spearman"`.
#' @param colors Named vector of concordance colours.
#' @param point_size Point size used when `geometry = "point"`.
#' @param label_size Text size for labeled genes.
#' @param base_size Base theme size.
#' 
#' @details Points are coloured by concordance status using fixed colours:
#' \itemize{
#'   \item Up/Up = `#1b9e77`
#'   \item Down/Down = `#7570b3`
#'   \item Up/Down = `#d95f02`
#'   \item Down/Up = `#e7298a`
#'   \item Other = `grey70`
#' }
#' Regulation is derived from the `log2fc` and `pval` cutoffs stored in
#' `cutoffs(x)` (and `p_value_type` from the same list, defaulting to `"padj"`).
#'
#' @return A `ggplot2` object.
#' @export
get_foldchange_scatter <- function(x,
                                   sample_comparisons,
                                   genes = NULL,
                                   display_id = NULL,
                                   display_from = NULL,
                                   display_orgdb = NULL,
                                   label_n = 0,
                                   alpha = 0.5,
                                   geometry = c("point", "hex"),
                                   method = c("pearson", "spearman"),
                                   colors = c(
                                     "Up/Up" = "#1b9e77",
                                     "Down/Down" = "#7570b3",
                                     "Up/Down" = "#d95f02",
                                     "Down/Up" = "#e7298a",
                                     "Other" = "grey70"
                                   ),
                                   point_size = 1.5,
                                   label_size = 3,
                                   base_size = 12) {
  stopifnot(inherits(x, "VISTA"))
  geometry <- match.arg(geometry)
  method <- match.arg(method)
  if (length(sample_comparisons) != 2) {
    cli::cli_abort("{.arg sample_comparisons} must contain exactly two comparisons.")
  }
  comps <- comparisons(x)
  missing <- setdiff(sample_comparisons, names(comps))
  if (length(missing)) cli::cli_abort("Unknown comparisons: {.val {missing}}.")

  cuts <- cutoffs(x)
  lfc_cut <- cuts$log2fc %||% cuts$lfc %||% 1
  p_cut   <- cuts$pval %||% cuts$pvalue %||% cuts$alpha %||% 0.05
  p_col_preferred <- cuts$p_value_type %||% "padj"

  clean_de <- function(df) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    if (!"gene_id" %in% names(df)) {
      rn <- rownames(df)
      if (!is.null(rn) && length(rn) == nrow(df) && all(nzchar(rn)) && !anyDuplicated(rn)) {
        df$gene_id <- rn
      } else {
        cli::cli_abort("DE table lacks a 'gene_id' column and usable rownames.")
      }
    }
    fc_col <- dplyr::coalesce(
      c(if ("log2fc" %in% names(df)) "log2fc" else NA_character_,
        if ("log2FoldChange" %in% names(df)) "log2FoldChange" else NA_character_,
        if ("logFC" %in% names(df)) "logFC" else NA_character_)
    )[1]
    if (is.na(fc_col)) cli::cli_abort("Cannot locate a log2FC column in the DE table.")
    p_col <- if (p_col_preferred %in% names(df)) {
      p_col_preferred
    } else if ("padj" %in% names(df)) {
      "padj"
    } else if ("pvalue" %in% names(df)) {
      "pvalue"
    } else if ("FDR" %in% names(df)) {
      "FDR"
    } else if ("PValue" %in% names(df)) {
      "PValue"
    } else {
      cli::cli_abort("Cannot locate a p-value column in the DE table.")
    }
    df[, c("gene_id", fc_col, p_col), drop = FALSE] |>
      stats::setNames(c("gene_id", "log2fc", "pval"))
  }

  d1 <- clean_de(comps[[sample_comparisons[1]]])
  d2 <- clean_de(comps[[sample_comparisons[2]]])
  merged <- dplyr::inner_join(d1, d2, by = "gene_id", suffix = c("1", "2"))
  genes_use <- .resolve_foldchange_gene_ids(
    x = x,
    genes = genes,
    display_id = display_id,
    display_from = display_from,
    display_orgdb = display_orgdb
  )
  if (!is.null(genes_use)) {
    merged <- dplyr::filter(merged, .data$gene_id %in% genes_use)
  }
  if (!nrow(merged)) cli::cli_abort("No overlapping genes between the two comparisons.")
  merged$display_gene <- .resolve_foldchange_gene_labels(
    x = x,
    gene_ids = merged$gene_id,
    display_id = display_id,
    display_from = display_from,
    display_orgdb = display_orgdb
  )

  merged <- merged |>
    dplyr::mutate(
      status = dplyr::case_when(
        log2fc1 >=  lfc_cut & pval1 <= p_cut & log2fc2 >=  lfc_cut & pval2 <= p_cut ~ "Up/Up",
        log2fc1 <= -lfc_cut & pval1 <= p_cut & log2fc2 <= -lfc_cut & pval2 <= p_cut ~ "Down/Down",
        log2fc1 >=  lfc_cut & pval1 <= p_cut & log2fc2 <= -lfc_cut & pval2 <= p_cut ~ "Up/Down",
        log2fc1 <= -lfc_cut & pval1 <= p_cut & log2fc2 >=  lfc_cut & pval2 <= p_cut ~ "Down/Up",
        TRUE ~ "Other"
      )
    )

  cor_val <- suppressWarnings(stats::cor(merged$log2fc1, merged$log2fc2, method = method, use = "complete.obs"))
  subtitle <- sprintf("%s correlation (n=%s): %.2f", tools::toTitleCase(method), nrow(merged), cor_val)

  status_levels <- c("Up/Up", "Down/Down", "Up/Down", "Down/Up", "Other")
  merged$status <- factor(merged$status, levels = status_levels)
  colors <- stats::setNames(as.character(colors), names(colors))
  if (is.null(names(colors)) || any(!status_levels %in% names(colors))) {
    cli::cli_abort("{.arg colors} must include named colours for {.val {status_levels}}.")
  }
  colors <- colors[status_levels]

  plt <- ggplot2::ggplot(merged, ggplot2::aes(x = log2fc1, y = log2fc2, color = status))
  if (geometry == "hex") {
    plt <- plt + ggplot2::geom_hex(show.legend = FALSE)
  } else {
    plt <- plt + ggplot2::geom_point(alpha = alpha, size = point_size)
  }
  plt <- plt +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dotted", color = "grey70") +
    ggplot2::geom_hline(yintercept = c(-lfc_cut, lfc_cut), linetype = "dotted", color = "grey70") +
    ggplot2::scale_color_manual(values = colors, drop = FALSE) +
    ggplot2::labs(
      x = paste0(sample_comparisons[1], " log2FC"),
      y = paste0(sample_comparisons[2], " log2FC"),
      color = "Regulation",
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal(base_size = base_size)

  if (label_n > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    idx <- order(abs(merged$log2fc1) + abs(merged$log2fc2), decreasing = TRUE)
    to_lab <- merged[idx[seq_len(min(label_n, length(idx)))], ]
    plt <- plt + ggrepel::geom_text_repel(
      data = to_lab,
      ggplot2::aes(label = display_gene),
      size = label_size,
      max.overlaps = Inf
    )
  }

  plt
}

# ──────────────────────────────────────────────────────────────────────────────
# Expression barplot
# ──────────────────────────────────────────────────────────────────────────────

#' Plot gene expression as barplots
#'
#' Displays selected genes as grouped summary bars or individual sample-level
#' bars. By default, expression is summarized per group with mean ± SD bars.
#' With `by = "sample"`, each sample is drawn separately.
#'
#' @param x A `VISTA` object.
#' @param genes Character vector (<=25 genes) to plot.
#' @param sample_group Optional character vector of groups (from `group_column`) to include.
#' @param log_transform Logical; log2-transform expression before plotting.
#' @param stats_group Logical; add statistical comparisons between groups when `TRUE`.
#' @param facet_scale Scaling option passed to `facet_wrap()` (deprecated; use `facet_scales`).
#' @param facet_scales Facet scales argument passed to `facet_wrap()` when faceting by gene.
#' @param facet_nrow,facet_ncol Optional layout passed to `facet_wrap()` when faceting.
#' @param p.label Label format for `ggpubr::stat_compare_means()`.
#' @param comparisons Optional list of specific group comparisons for `stat_compare_means()`.
#' @param group_column Optional column name in `sample_info` to use for grouping samples.
#' @param display_id Optional ID/column name to use for labels/facets. If supplied
#'   and present in `rowData(x)`, those values are used; otherwise falls back to
#'   ID mapping.
#' @param display_from Optional source ID type for mapping (used when `display_id`
#'   is not found in `rowData`).
#' @param display_orgdb Optional `OrgDb` object used for ID mapping when
#'   `display_id` is set but not found in `rowData`.
#' @param by One of `"group"` (default; summarize replicates by group) or
#'   `"sample"` (show one bar per sample).
#' @param sample_order Ordering used when `by = "sample"`:
#'   `"input"` preserves the current sample order, `"group"` groups samples by
#'   `group_column`, and `"expression"` ranks samples by mean expression across
#'   the selected genes.
#' @param fill_by Fill mapping for `by = "sample"` barplots. Use `"group"`
#'   (default) or any discrete column from the joined plotting data, such as a
#'   sample metadata column or `"sample"`. Group-summary barplots
#'   (`by = "group"`) only support group-based fill because each bar already
#'   represents an aggregated group mean.
#' @param facet_by Faceting mode: `"auto"` (default; facet by gene when more
#'   than one gene is requested), `"gene"`, or `"none"`. For multiple genes,
#'   `"none"` falls back to `"gene"` to preserve readability.
#' @return A `ggplot2` object.
#'
#' @examples
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
#'
#' # Plot expression for select genes
#' genes <- rownames(vista)[seq_len(3)]
#' get_expression_barplot(vista, genes = genes)
#'
#' # With statistics
#' get_expression_barplot(vista, genes = genes, stats_group = TRUE)
#'
#' # Without log transformation
#' get_expression_barplot(vista, genes = genes, log_transform = FALSE)
#'
#' @export
get_expression_barplot <- function(x,
                                   genes,
                                   sample_group = NULL,
                                   group_column = NULL,
                                   log_transform = TRUE,
                                   stats_group = FALSE,
                                   facet_scale = "free_y",
                                   facet_scales = facet_scale,
                                   facet_nrow = NULL,
                                   facet_ncol = NULL,
                                   p.label = "p.signif",
                                   comparisons = NULL,
                                   display_id = NULL,
                                   display_from = NULL,
                                   display_orgdb = NULL,
                                   by = c("group", "sample"),
                                   sample_order = c("input", "group", "expression"),
                                   fill_by = NULL,
                                   facet_by = c("auto", "gene", "none")) {
  stopifnot(inherits(x, "VISTA"))
  if (length(genes) > 25) cli::cli_abort("Maximum 25 genes can be plotted at once.")
  by <- match.arg(by)
  sample_order <- match.arg(sample_order)

  meta <- .prepare_sample_metadata(x, sample_group, group_column)
  group_col <- attr(meta, "group_column")

  if (stats_group && length(unique(meta[[group_col]])) < 2) {
    cli::cli_abort("At least two groups are needed for statistical testing.")
  }
  if (stats_group && by == "sample") {
    cli::cli_abort("{.arg stats_group = TRUE} is only supported when {.arg by = 'group'}.")
  }

  rd <- tryCatch(SummarizedExperiment::rowData(x), error = function(e) NULL)
  mat <- SummarizedExperiment::assay(x)[, meta$sample, drop = FALSE]

  # Map input genes to underlying IDs if possible
  underlying_ids <- rownames(mat)
  if (!is.null(display_id) && !is.null(rd) && display_id %in% colnames(rd)) {
    map <- rd[[display_id]]
    names(map) <- rownames(x)
    gene_ids <- names(map)[match(genes, map)]
    gene_ids <- gene_ids[!is.na(gene_ids)]
  } else if (!is.null(display_id) && !is.null(display_from) && !is.null(display_orgdb)) {
    # convert from display_id -> underlying (display_from)
    gene_ids <- .map_gene_ids(genes,
                              from_type = display_id,
                              to_type = display_from,
                              orgdb = display_orgdb)
  } else {
    gene_ids <- genes
  }
  if (!length(gene_ids)) gene_ids <- genes

  missing_genes <- setdiff(gene_ids, underlying_ids)
  if (length(missing_genes) == length(gene_ids)) cli::cli_abort("None of the specified {.arg genes} were found.")
  keep_genes <- intersect(gene_ids, rownames(mat))
  mat <- mat[keep_genes, , drop = FALSE]

  df <- mat |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    tidyr::pivot_longer(-gene, names_to = "sample", values_to = "expression") |>
    dplyr::left_join(meta, by = "sample")

  # Replace gene labels with display_id when available
  if (!is.null(display_id) && !is.null(rd) && display_id %in% colnames(rd)) {
    map <- rd[[display_id]]
    names(map) <- rownames(x)
    mapped <- map[match(df$gene, names(map))]
    df$gene <- ifelse(!is.na(mapped) & nzchar(mapped), mapped, df$gene)
  } else if (!is.null(display_id) && !is.null(display_from) && !is.null(display_orgdb)) {
    df$gene <- .map_gene_ids(df$gene,
                             from_type = display_from,
                             to_type = display_id,
                             orgdb = display_orgdb)
  }

  if (log_transform) {
    df$expression <- log2(df$expression + 1)
    ylab <- "log2(Normalized Counts + 1)"
  } else {
    ylab <- "Normalized Counts"
  }

  facet_mode <- .resolve_expression_plot_facet(facet_by, n_genes = length(unique(df$gene)))
  cols <- .vista_group_colors(x, df[[group_col]])

  if (by == "group") {
    if (!is.null(fill_by) && !identical(fill_by, "group") && !identical(fill_by, group_col)) {
      cli::cli_abort(
        "{.arg fill_by} is only supported for {.code by = 'sample'} in {.fun get_expression_barplot}. Group-summary bars must use the grouping variable."
      )
    }
    dodge <- ggplot2::position_dodge(width = 0.8)
    mean_sd_summary <- function(y) {
      m <- mean(y, na.rm = TRUE)
      s <- stats::sd(y, na.rm = TRUE)
      data.frame(y = m, ymin = m - s, ymax = m + s)
    }

    plt <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = .data[[group_col]], y = expression, fill = .data[[group_col]])
    ) +
      ggplot2::stat_summary(fun = mean, geom = "bar", position = dodge, width = 0.6) +
      ggplot2::stat_summary(fun.data = mean_sd_summary, geom = "errorbar",
                            position = dodge, width = 0.2) +
      ggplot2::labs(x = group_col, y = ylab, fill = group_col) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  } else {
    df <- .order_expression_plot_samples(df, meta, group_col, sample_order = sample_order)
    fill_spec <- .resolve_expression_fill(
      x = x,
      df = df,
      group_col = group_col,
      fill_by = fill_by,
      default_fill = "group",
      allowed_special = c("group", "sample", "x"),
      x_var = "sample",
      x_label = "Sample"
    )
    plt <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data$sample,
        y = expression,
        fill = factor(fill_spec$values, levels = fill_spec$levels)
      )
    ) +
      ggplot2::geom_col(alpha = 0.85, width = 0.7) +
      ggplot2::labs(x = "Sample", y = ylab, fill = fill_spec$label) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    cols <- fill_spec$palette
  }

  if (stats_group) {
    if (!requireNamespace("ggpubr", quietly = TRUE)) {
      cli::cli_abort("Package {.pkg ggpubr} must be installed for `stats_group = TRUE`.")
    }
    plt <- plt + ggpubr::stat_compare_means(
      ggplot2::aes(group = .data[[group_col]]),
      comparisons = comparisons,
      method = "t.test",
      label = p.label,
      label.x.npc = "center"
    )
  }

  if (!is.null(cols)) {
    plt <- plt + ggplot2::scale_fill_manual(values = cols)
  }

  if (facet_mode == "gene") {
    plt <- .add_expression_facet_wrap(
      plot = plt,
      facet_formula = ~gene,
      scales = facet_scales,
      facet_nrow = facet_nrow,
      facet_ncol = facet_ncol
    )
  }

  plt
}

#' Fold-change plotting helpers (overview)
#'
#' One-stop doc for fold-change plots:
#' \itemize{
#'   \item \code{get_foldchange_barplot()}: log2FC by comparison (bars).
#'   \item \code{get_foldchange_boxplot()}: log2FC distributions per comparison (boxes).
#'   \item \code{get_foldchange_lollipop()}: log2FC stems/dots; supports 1–2 comparisons.
#'   \item \code{get_foldchange_lineplot()}: log2FC trajectories across comparisons (optional clustering).
#' }
#'
#' Shared arguments: \code{x} (`VISTA` with DE results), \code{sample_comparisons/sample_comparison},
#' \code{genes}, \code{display_id} for label mapping, \code{sort_by} (where supported),
#' faceting controls (\code{facet_*}), and comparison colours pulled from
# ──────────────────────────────────────────────────────────────────────────────
# Fold-change lollipop
# ──────────────────────────────────────────────────────────────────────────────

#' Plot log2 fold changes as a lollipop chart (one or two comparisons)
#' 
#' Extracts log2 fold changes from stored differential expression results and
#' plots them as stems and dots, with labels and a zero reference line. You can
#' optionally provide two comparisons; in that case both comparisons are drawn
#' side-by-side, coloured by comparison.
#'
#' @param x A `VISTA` object.
#' @param sample_comparison Character vector of length 1 or 2 naming the
#'   comparison(s) to plot (must exist in `metadata(x)$de_results`).
#' @param genes Optional character vector of gene IDs to include. When `NULL`,
#'   all genes in the specified comparison(s) are shown.
#' @param sort_by How to sort genes on the y-axis: `"input"` (use supplied
#'   order), `"log2fc"` (descending log2FC; for two comparisons uses the first
#'   comparison), or `"abs_log2fc"` (descending max abs log2FC across
#'   comparisons).
#' @param palette For a single comparison, named vector of colors for `"pos"`,
#'   `"neg"`, and `"zero"` sign classes (set to `NULL` to disable color by
#'   sign). For two comparisons, a named or unnamed vector of colors with one
#'   entry per comparison (defaults to a qualitative palette).
#' @param point_size Numeric size of dots.
#' @param line_size Numeric size of stems (linewidth).
#' @param label Logical; draw numeric labels next to the dots.
#' @param label_digits Integer; digits to show in labels when `label = TRUE`.
#' @param display_id Optional column in `rowData(x)` used to interpret `genes`
#'   and to label plotted genes.
#' @param display_from Optional source ID type for mapping when `display_id` is
#'   not present in `rowData(x)`.
#' @param display_orgdb Optional `OrgDb` object used for identifier mapping when
#'   `display_id` is not present in `rowData(x)`.
#' @param dodge_width Horizontal separation between comparisons when plotting
#'   two comparisons on the same axis.
#' @param facet_scales Facet scales argument passed to `facet_wrap()` when
#'   `facet_by != "none"` (default `"free_y"`).
#' @param facet_nrow,facet_ncol Optional layout passed to `facet_wrap()` when faceting.
#' @param facet_by Faceting mode: `"auto"` (default), `"gene"`,
#'   `"comparison"`, or `"none"`.
#'
#' @return A `ggplot2` object.
#'
#' @examples
#' vista <- example_vista()
#' comp_name <- names(comparisons(vista))[1]
#' genes <- rownames(vista)[seq_len(3)]
#' get_foldchange_lollipop(vista, sample_comparison = comp_name, genes = genes)
#' @export
get_foldchange_lollipop <- function(x,
                                    sample_comparison,
                                    genes = NULL,
                                    sort_by = c("input", "log2fc", "abs_log2fc"),
                                    palette = NULL,
                                    point_size = 6,
                                    line_size = 1.2,
                                    label = TRUE,
                                    label_digits = 2,
                                    display_id = NULL,
                                    display_from = NULL,
                                    display_orgdb = NULL,
                                    dodge_width = 0.5,
                                    facet_scales = "free_y",
                                    facet_nrow = NULL,
                                    facet_ncol = NULL,
                                    facet_by = c("auto", "gene", "comparison", "none")) {
  stopifnot(inherits(x, "VISTA"))
  sort_by <- match.arg(sort_by)

  rd <- tryCatch(SummarizedExperiment::rowData(x), error = function(e) NULL)
  comps <- .vista_comparisons(x)
  if (is.null(comps) || !length(comps)) {
    cli::cli_abort("No differential expression results found in the VISTA object.")
  }
  if (missing(sample_comparison) || is.null(sample_comparison) || !length(sample_comparison)) {
    cli::cli_abort("Please provide one or two values in {.arg sample_comparison} present in the DE results.")
  }
  sample_comparison <- as.character(sample_comparison)
  if (length(sample_comparison) > 2) {
    cli::cli_abort("{.arg sample_comparison} must have length 1 or 2.")
  }

  # helper to extract and validate a comparison table
  get_comp <- function(name) {
    if (!name %in% names(comps)) {
      cli::cli_abort("Comparison {.val {name}} not found in DE results.")
    }
    tbl <- comps[[name]]
    if (!all(c("gene_id", "log2fc") %in% colnames(tbl))) {
      cli::cli_abort("DE table for {.val {name}} must contain columns {.field gene_id} and {.field log2fc}.")
    }
    tbl
  }

  # build long table
  comp_tbls <- lapply(sample_comparison, function(nm) {
    tbl <- get_comp(nm)
    genes_use <- .resolve_foldchange_gene_ids(
      x = x,
      genes = genes,
      display_id = display_id,
      display_from = display_from,
      display_orgdb = display_orgdb
    )
    if (!is.null(genes_use)) {
      tbl <- dplyr::filter(tbl, .data$gene_id %in% genes_use)
    }
    if (!nrow(tbl)) cli::cli_abort("No genes to plot for comparison {.val {nm}} after filtering.")
    dplyr::transmute(tbl, gene_id = .data$gene_id, log2fc = .data$log2fc, comparison = nm)
  })
  df <- dplyr::bind_rows(comp_tbls)
  if (!nrow(df)) cli::cli_abort("No genes to plot after filtering and merging comparisons.")

  # optional display labels
  df$gene_label <- .resolve_foldchange_gene_labels(
    x = x,
    gene_ids = df$gene_id,
    display_id = display_id,
    display_from = display_from,
    display_orgdb = display_orgdb
  )

  # sort genes
  if (sort_by == "input" && !is.null(genes)) {
    gene_levels <- unique(df$gene_id[match(genes, df$gene_id)])
  } else if (sort_by == "log2fc") {
    first_comp <- sample_comparison[1]
    ref <- dplyr::filter(df, .data$comparison == first_comp)
    gene_levels <- ref$gene_id[order(ref$log2fc, decreasing = TRUE)]
  } else if (sort_by == "abs_log2fc") {
    ref <- df |>
      dplyr::group_by(.data$gene_id) |>
      dplyr::summarise(max_abs_fc = max(abs(.data$log2fc), na.rm = TRUE), .groups = "drop")
    gene_levels <- ref$gene_id[order(ref$max_abs_fc, decreasing = TRUE)]
  } else {
    gene_levels <- unique(df$gene_id)
  }
  gene_levels <- unique(gene_levels[!is.na(gene_levels)])
  df$gene_id <- factor(df$gene_id, levels = rev(gene_levels))
  facet_mode <- .resolve_foldchange_facet_mode(
    facet_by = facet_by,
    n_genes = length(gene_levels)
  )

  single_comp <- length(sample_comparison) == 1
  x_labeller <- stats::setNames(df$gene_label, as.character(df$gene_id))

  if (facet_mode == "gene") {
    df$comparison <- factor(df$comparison, levels = sample_comparison)
    if (single_comp) {
      df$fc_dir <- dplyr::case_when(
        df$log2fc > 0 ~ "pos",
        df$log2fc < 0 ~ "neg",
        TRUE ~ "zero"
      )
      plt <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = comparison, y = log2fc, color = .data$fc_dir)
      ) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
        ggplot2::geom_segment(
          ggplot2::aes(xend = comparison, y = 0, yend = log2fc),
          linewidth = line_size,
          lineend = "round"
        ) +
        ggplot2::geom_point(size = point_size) +
        ggplot2::labs(
          x = "comparison",
          y = "log2 fold change",
          color = NULL
        ) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

      plt <- .add_expression_facet_wrap(
        plot = plt,
        facet_formula = ~gene_id,
        scales = facet_scales,
        facet_nrow = facet_nrow,
        facet_ncol = facet_ncol
      )
      plt$facet$params$labeller <- ggplot2::as_labeller(x_labeller)

      if (label) {
        plt <- plt + ggplot2::geom_text(
          ggplot2::aes(label = round(log2fc, label_digits), color = .data$fc_dir),
          vjust = -0.6,
          size = 4
        )
      }

      if (!is.null(palette)) {
        plt <- plt + ggplot2::scale_color_manual(values = palette)
      } else {
        plt <- plt + ggplot2::scale_color_manual(values = c(pos = "red4", neg = "blue4", zero = "grey50"))
      }
    } else {
      plt <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = comparison, y = log2fc, color = .data$comparison)
      ) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
        ggplot2::geom_segment(
          ggplot2::aes(xend = comparison, y = 0, yend = log2fc),
          linewidth = line_size,
          lineend = "round"
        ) +
        ggplot2::geom_point(size = point_size) +
        ggplot2::labs(
          x = "comparison",
          y = "log2 fold change",
          color = "comparison"
        ) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

      plt <- .add_expression_facet_wrap(
        plot = plt,
        facet_formula = ~gene_id,
        scales = facet_scales,
        facet_nrow = facet_nrow,
        facet_ncol = facet_ncol
      )
      plt$facet$params$labeller <- ggplot2::as_labeller(x_labeller)

      if (label) {
        plt <- plt + ggplot2::geom_text(
          ggplot2::aes(label = round(log2fc, label_digits), color = .data$comparison),
          vjust = -0.6,
          size = 4
        )
      }

      if (is.null(palette)) {
        pal_use <- .vista_comparison_colors(x, sample_comparison)
        if (is.null(pal_use)) {
          pal_use <- colorspace::qualitative_hcl(length(sample_comparison))
          names(pal_use) <- sample_comparison
        }
      } else {
        pal_use <- palette
        if (is.null(names(pal_use))) names(pal_use) <- sample_comparison
      }
      plt <- plt + ggplot2::scale_color_manual(values = pal_use)
    }

    return(plt)
  }

  base_plot <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = gene_id, y = log2fc)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
    ggplot2::scale_x_discrete(labels = x_labeller) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(panel.grid.minor.y = ggplot2::element_line(color = "grey40"))

  if (single_comp) {
    df$fc_dir <- dplyr::case_when(
      df$log2fc > 0 ~ "pos",
      df$log2fc < 0 ~ "neg",
      TRUE ~ "zero"
    )
    plt <- base_plot +
      ggplot2::geom_segment(
        ggplot2::aes(xend = gene_id, y = 0, yend = log2fc, color = fc_dir),
        linewidth = line_size,
        lineend = "round",
        data = df
      ) +
      ggplot2::geom_point(
        ggplot2::aes(color = fc_dir),
        size = point_size,
        data = df
      )
    if (label) {
      plt <- plt + ggplot2::geom_text(
        ggplot2::aes(label = round(log2fc, label_digits)),
        hjust = -0.3,
        size = 4,
        data = df
      )
    }
    plt <- plt + ggplot2::labs(
      x = "gene",
      y = "log2 fold change",
      color = NULL
    )
    if (is.null(palette)) {
      palette <- c(pos = "red4", neg = "blue4", zero = "grey50")
    }
    if (!is.null(palette)) {
      plt <- plt + ggplot2::scale_color_manual(values = palette)
    }
  } else {
    facet_comparison <- identical(facet_mode, "comparison")
    # compute fixed offsets per comparison so stems remain vertical
    comp_levels <- sample_comparison
    comp_idx <- match(df$comparison, comp_levels)
    n_comp <- length(comp_levels)
    offsets <- (comp_idx - (n_comp + 1) / 2) * dodge_width
    df$x_dodge <- as.numeric(df$gene_id) + offsets

    level_labels <- stats::setNames(
      vapply(levels(df$gene_id), function(g) df$gene_label[df$gene_id == g][1], character(1)),
      seq_along(levels(df$gene_id))
    )

    if (facet_comparison) {
      plt <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = gene_id, y = log2fc, color = .data$comparison)
      ) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
        ggplot2::geom_segment(
          ggplot2::aes(xend = gene_id, y = 0, yend = log2fc),
          linewidth = line_size,
          lineend = "round"
        ) +
        ggplot2::geom_point(size = point_size) +
        ggplot2::labs(
          x = "gene",
          y = "log2 fold change",
          color = "comparison"
        ) +
        ggplot2::coord_flip() +
        ggplot2::scale_x_discrete(labels = x_labeller) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(panel.grid.minor.y = ggplot2::element_line(color = "grey40"))

      plt <- .add_expression_facet_wrap(
        plot = plt,
        facet_formula = ~comparison,
        scales = facet_scales,
        facet_nrow = facet_nrow,
        facet_ncol = facet_ncol
      )

      if (label) {
        plt <- plt + ggplot2::geom_text(
          ggplot2::aes(label = round(log2fc, label_digits)),
          hjust = -0.3,
          size = 4
        )
      }
    } else {
      plt <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = x_dodge, y = log2fc, color = .data$comparison)
      ) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
        ggplot2::geom_segment(
          ggplot2::aes(xend = x_dodge, y = 0, yend = log2fc),
          linewidth = line_size,
          lineend = "round"
        ) +
        ggplot2::geom_point(size = point_size) +
        ggplot2::labs(
          x = "gene",
          y = "log2 fold change",
          color = "comparison"
        ) +
        ggplot2::coord_flip() +
        ggplot2::scale_x_continuous(
          breaks = seq_along(levels(df$gene_id)),
          labels = level_labels
        ) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(panel.grid.minor.y = ggplot2::element_line(color = "grey40"))

      if (label) {
        plt <- plt + ggplot2::geom_text(
          ggplot2::aes(label = round(log2fc, label_digits), color = .data$comparison),
          hjust = -0.3,
          size = 4
        )
      }
    }

    if (is.null(palette)) {
      pal_use <- .vista_comparison_colors(x, sample_comparison)
      if (is.null(pal_use)) {
        pal_use <- colorspace::qualitative_hcl(length(sample_comparison))
        names(pal_use) <- sample_comparison
      }
    } else {
      pal_use <- palette
      if (is.null(names(pal_use))) names(pal_use) <- sample_comparison
    }
    plt <- plt + ggplot2::scale_color_manual(values = pal_use)
  }

  plt
}

# ──────────────────────────────────────────────────────────────────────────────
# DEG Venn
# ──────────────────────────────────────────────────────────────────────────────

#' @title DEG Venn diagram
#' @description Visualizes overlaps between DEG sets for two to four comparisons.
#' @param x A `VISTA` object.
#' @param sample_comparisons Character vector of 2–4 comparison names to include in the Venn diagram.
#' @param regulation One of `"Up"`, `"Down"`, `"Both"`, or `"All"` selecting which genes to include.
#' @param palette Qualitative palette name passed to `colorspace::qualitative_hcl()` for fill colors.
#' @param auto_scale Logical; pass through to `ggvenn::ggvenn()` to scale circles by size.
#' @param show_percentage Logical; request percentage labels from `ggvenn::ggvenn()`.
#' @param ... Additional arguments forwarded to `ggvenn::ggvenn()`.
#' @export
get_deg_venn_diagram <- function(x,
                                 sample_comparisons,
                                 regulation = "Up",
                                 palette = "Set 2",
                                 auto_scale = FALSE,
                                 show_percentage = TRUE,
                                 ...) {
  stopifnot(inherits(x, "VISTA"))
  if (!requireNamespace("ggvenn", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggvenn} must be installed to use `get_deg_venn_diagram()`.")
  }
  if (length(sample_comparisons) < 2 || length(sample_comparisons) > 4) {
    cli::cli_abort("{.arg sample_comparisons} must contain 2 to 4 comparison names.")
  }

  gene_sets <- get_genes_by_regulation(
    x = x,
    sample_comparisons = sample_comparisons,
    regulation = regulation
  )

  fill_colors <- colorspace::qualitative_hcl(length(gene_sets), palette = palette)
  names(fill_colors) <- names(gene_sets)

  ggvenn::ggvenn(
    data = gene_sets,
    fill_color = as.character(fill_colors),
    auto_scale = auto_scale,
    show_percentage = show_percentage,
    ...
  )
}

# ──────────────────────────────────────────────────────────────────────────────
# DEG count barplot
# ──────────────────────────────────────────────────────────────────────────────

.normalize_deg_colors <- function(colors,
                                  include_other = FALSE,
                                  other_color = "grey70") {
  if (!is.character(colors) || !length(colors)) {
    cli::cli_abort("{.arg colors} must be a character vector.")
  }
  if (is.null(names(colors)) || any(!nzchar(names(colors)))) {
    min_len <- if (isTRUE(include_other)) 3 else 2
    if (length(colors) < min_len) {
      cli::cli_abort("{.arg colors} must provide at least two colors for {.val Up} and {.val Down}.")
    }
    color_names <- c("Up", "Down", "Other")
    colors <- stats::setNames(colors[seq_len(min_len)], color_names[seq_len(min_len)])
  }
  missing_keys <- setdiff(c("Up", "Down"), names(colors))
  if (length(missing_keys)) {
    cli::cli_abort("{.arg colors} must include named colors for {.val Up} and {.val Down}. Missing: {.val {missing_keys}}")
  }
  if (isTRUE(include_other)) {
    if (!"Other" %in% names(colors)) {
      colors[["Other"]] <- other_color
    }
    return(colors[c("Up", "Down", "Other")])
  }
  colors[c("Up", "Down")]
}

.collect_deg_count_data <- function(x,
                                    sample_comparisons = NULL,
                                    show_other = FALSE) {
  stopifnot(inherits(x, "VISTA"))
  if (!is.null(sample_comparisons)) sample_comparisons <- as.character(sample_comparisons)

  # pull summary from metadata (new layout)
  deg_summary <- .vista_deg_summary(x)
  if (is.null(deg_summary) || length(deg_summary) == 0) {
    cli::cli_warn("No DEG summary found in the VISTA object.")
    return(NULL)
  }
  # keep only named summaries; optionally restrict to requested comparisons
  if (is.null(names(deg_summary))) names(deg_summary) <- rep(NA_character_, length(deg_summary))
  has_name <- !is.na(names(deg_summary)) & nzchar(names(deg_summary))
  deg_summary <- deg_summary[has_name]
  if (!length(deg_summary)) {
    cli::cli_warn("No named DEG summaries found in the VISTA object.")
    return(NULL)
  }
  if (!is.null(sample_comparisons)) {
    missing <- setdiff(sample_comparisons, names(deg_summary))
    if (length(missing)) {
      cli::cli_abort("Some values in {.arg sample_comparisons} not found in DEG summary: {.val {missing}}")
    }
    deg_summary <- deg_summary[sample_comparisons]
  }

  # standardize one summary table -> columns: regulation, n
  .std_deg_sum <- function(d) {
    d <- as.data.frame(d, stringsAsFactors = FALSE)
    # regulation column
    if (!"regulation" %in% names(d)) {
      if ("Var1" %in% names(d)) names(d)[names(d) == "Var1"] <- "regulation"
      if ("reg" %in% names(d)) names(d)[names(d) == "reg"] <- "regulation"
    }
    # count column
    count_hit <- intersect(c("n", "N", "count", "Counts", "Freq", "freq", "value"), names(d))
    if (length(count_hit)) names(d)[names(d) == count_hit[1]] <- "n"

    # if still missing, try to compute
    if (!"n" %in% names(d) && "regulation" %in% names(d)) {
      d <- dplyr::count(d, .data$regulation, name = "n")
    }
    if (!"n" %in% names(d)) d$n <- NA_real_
    d
  }

  df <- purrr::imap_dfr(deg_summary, ~ .std_deg_sum(.x) |>
                          dplyr::mutate(sample_comparisons = .y))

  if (!nrow(df) || !"regulation" %in% names(df)) {
    cli::cli_warn("DEG summaries do not contain a usable {.field regulation} column.")
    return(NULL)
  }

  df <- df |>
    dplyr::filter(
      !is.na(.data$sample_comparisons),
      .data$regulation %in% c("Up", "Down")
    ) |>
    dplyr::mutate(n = suppressWarnings(as.numeric(.data$n)))

  if (!nrow(df)) {
    cli::cli_warn("No {.val Up}/{.val Down} DEG counts were found in the summary.")
    return(NULL)
  }

  df$n[!is.finite(df$n)] <- NA_real_
  df <- df |>
    dplyr::group_by(.data$sample_comparisons, .data$regulation) |>
    dplyr::summarise(n = sum(.data$n, na.rm = TRUE), .groups = "drop")

  available_comparisons <- unique(df$sample_comparisons)
  if (!is.null(sample_comparisons)) {
    missing <- setdiff(sample_comparisons, available_comparisons)
    if (length(missing)) {
      cli::cli_abort("Some values in {.arg sample_comparisons} not found in DEG summary: {.val {missing}}")
    }
    comp_levels <- sample_comparisons
  } else {
    comp_levels <- available_comparisons
  }

  df <- tidyr::complete(
    df,
    sample_comparisons = comp_levels,
    regulation = c("Up", "Down"),
    fill = list(n = 0)
  )

  if (isTRUE(show_other)) {
    comps <- .vista_comparisons(x)
    if (!is.null(sample_comparisons)) {
      comps <- comps[sample_comparisons]
    }

    other_df <- purrr::imap_dfr(comps, function(comp_tbl, comp_name) {
      comp_tbl <- as.data.frame(comp_tbl, stringsAsFactors = FALSE)
      gene_col <- intersect(c("gene_id", "gene", "geneID", "GeneID", "ID", "id"), colnames(comp_tbl))[1]
      gene_ids <- if (is.na(gene_col)) rownames(comp_tbl) else comp_tbl[[gene_col]]
      gene_ids <- unique(stats::na.omit(as.character(gene_ids)))
      gene_ids <- gene_ids[nzchar(gene_ids)]

      sig_n <- sum(df$n[df$sample_comparisons == comp_name], na.rm = TRUE)
      tibble::tibble(
        sample_comparisons = comp_name,
        regulation = "Other",
        n = max(length(gene_ids) - sig_n, 0)
      )
    })

    if (nrow(other_df)) {
      df <- dplyr::bind_rows(df, other_df)
    }
  }

  df$sample_comparisons <- factor(df$sample_comparisons, levels = comp_levels)
  df$regulation <- factor(
    df$regulation,
    levels = if (isTRUE(show_other)) c("Up", "Down", "Other") else c("Up", "Down")
  )
  df
}

.format_deg_slice_label <- function(n, pct, label = c("both", "count", "percent", "none"), label_digits = 1) {
  label <- match.arg(label)
  if (label == "none") return(rep("", length(n)))
  count_lab <- formatC(n, format = "f", digits = 0, big.mark = ",")
  pct_lab <- formatC(100 * pct, format = "f", digits = label_digits)
  pct_lab[!is.finite(pct)] <- formatC(0, format = "f", digits = label_digits)
  if (label == "both") {
    return(paste0(count_lab, " (", pct_lab, "%)"))
  }
  if (label == "count") return(count_lab)
  paste0(pct_lab, "%")
}

.plot_deg_count_circular <- function(df,
                                     label = c("both", "count", "percent", "none"),
                                     label_digits = 1,
                                     base_size = 12,
                                     colors = c(Up = "red4", Down = "blue4"),
                                     text_color = "black",
                                     donut = FALSE,
                                     facet_by = c("comparison", "none"),
                                     ncol = NULL) {
  facet_by <- match.arg(facet_by)
  label <- match.arg(label)
  colors <- .normalize_deg_colors(
    colors,
    include_other = "Other" %in% unique(as.character(df$regulation))
  )

  n_comp <- dplyr::n_distinct(df$sample_comparisons)
  if (facet_by == "none" && n_comp > 1) {
    cli::cli_abort("`facet_by = 'none'` requires exactly one comparison. Provide one {.arg sample_comparisons} value or use {.val facet_by = 'comparison'}.")
  }

  df <- df |>
    dplyr::group_by(.data$sample_comparisons) |>
    dplyr::mutate(
      total = sum(.data$n, na.rm = TRUE),
      pct = dplyr::if_else(.data$total > 0, .data$n / .data$total, 0),
      label = .format_deg_slice_label(
        .data$n,
        .data$pct,
        label = label,
        label_digits = label_digits
      )
    ) |>
    dplyr::ungroup()

  if (donut) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = 2, y = .data$n, fill = .data$regulation)) +
      ggplot2::geom_col(width = 1, color = "white") +
      ggplot2::coord_polar(theta = "y") +
      ggplot2::xlim(0.5, 2.5)
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = 1, y = .data$n, fill = .data$regulation)) +
      ggplot2::geom_col(width = 1, color = "white") +
      ggplot2::coord_polar(theta = "y")
  }

  p <- p +
    ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
    ggplot2::labs(fill = "Regulation") +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text = ggplot2::element_text(size = base_size)
    )

  if (label != "none") {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data$label),
      position = ggplot2::position_stack(vjust = 0.5),
      size = base_size / 3,
      color = text_color
    )
  }

  if (facet_by == "comparison") {
    p <- p + ggplot2::facet_wrap(~sample_comparisons, ncol = ncol)
  }

  p
}


#' Barplot of DEG counts (Up/Down) across comparisons
#' @param x A `VISTA` object containing DEG summaries.
#' @param sample_comparisons Optional character vector of comparison names to display.
#' @param label Logical; overlay numeric counts above bars when `TRUE`.
#' @param base_size Numeric base font size for the plot.
#' @param colors Named vector giving fill colors for `"Up"` and `"Down"` bars.
#' @param facet_by Either `"none"`, `"regulation"`, or `"comparison"` describing the facet variable.
#'   Use `"none"` for a single panel.
#' @param facet_scales Scale option passed to `facet_wrap()` when faceting.
#' @export
get_deg_count_barplot <- function(x,
                                  sample_comparisons = NULL,
                                  label = TRUE,
                                  base_size = 12,
                                  colors = c(Up = "red4", Down = "blue4"),
                                  facet_by = c("none", "regulation", "comparison"),
                                  facet_scales = "fixed") {
  stopifnot(inherits(x, "VISTA"))
  facet_by <- match.arg(facet_by)

  df <- .collect_deg_count_data(
    x = x,
    sample_comparisons = sample_comparisons
  )
  if (is.null(df)) return(invisible(NULL))
  colors <- .normalize_deg_colors(colors)

  # x-axis is always comparison; faceting optional
  df$sample_comparisons <- factor(df$sample_comparisons,
                                  levels = if (!is.null(sample_comparisons)) sample_comparisons else unique(df$sample_comparisons))
  df <- droplevels(df)
  aes_x <- "sample_comparisons"
  facet_formula <- switch(facet_by,
                          regulation = ~regulation,
                          comparison = ~sample_comparisons,
                          NULL)
  facet_scale_mode <- if (identical(facet_by, "comparison")) "free_x" else facet_scales

  # plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[aes_x]], y = .data$n, fill = .data$regulation)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::scale_x_discrete(drop = TRUE) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(x = NULL, y = "Gene Count", fill = "Regulation") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (label) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data$n),
      position = ggplot2::position_dodge(width = 0.9),
      vjust = -0.3,
      size = base_size / 3
    )
  }

  if (!is.null(facet_formula)) p <- p + ggplot2::facet_wrap(facet_formula, scales = facet_scale_mode)
  p
}

#' Pie chart of DEG counts (Up/Down) across comparisons
#' @param x A `VISTA` object containing DEG summaries.
#' @param sample_comparisons Optional character vector of comparison names to display.
#' @param label Label mode: `"both"` (default), `"count"`, `"percent"`, or `"none"`.
#' @param label_digits Integer number of decimals used for percentage labels.
#' @param base_size Numeric base font size for the plot.
#' @param colors Named vector giving fill colors for `"Up"` and `"Down"` slices.
#'   When `show_other = TRUE`, an `"Other"` entry may also be supplied.
#' @param show_other Logical; when `TRUE`, include non-DE genes as an `"Other"`
#'   slice using `other_color` unless overridden in `colors`.
#' @param other_color Fill colour used for the `"Other"` slice when
#'   `show_other = TRUE` and `colors` does not include `"Other"`.
#' @param text_color Colour used for pie label text.
#' @param facet_by Either `"comparison"` (default) to draw one pie per comparison,
#'   or `"none"` for a single pie.
#' @param ncol Optional number of columns when faceting.
#' @export
get_deg_count_pieplot <- function(x,
                                  sample_comparisons = NULL,
                                  label = c("both", "count", "percent", "none"),
                                  label_digits = 1,
                                  base_size = 12,
                                  colors = c(Up = "red4", Down = "blue4"),
                                  show_other = FALSE,
                                  other_color = "grey70",
                                  text_color = "black",
                                  facet_by = c("comparison", "none"),
                                  ncol = NULL) {
  stopifnot(inherits(x, "VISTA"))
  df <- .collect_deg_count_data(
    x = x,
    sample_comparisons = sample_comparisons,
    show_other = show_other
  )
  if (is.null(df)) return(invisible(NULL))
  colors <- .normalize_deg_colors(
    colors = colors,
    include_other = show_other,
    other_color = other_color
  )
  .plot_deg_count_circular(
    df = df,
    label = label,
    label_digits = label_digits,
    base_size = base_size,
    colors = colors,
    text_color = text_color,
    donut = FALSE,
    facet_by = facet_by,
    ncol = ncol
  )
}

#' Donut chart of DEG counts (Up/Down) across comparisons
#' @param x A `VISTA` object containing DEG summaries.
#' @param sample_comparisons Optional character vector of comparison names to display.
#' @param label Label mode: `"both"` (default), `"count"`, `"percent"`, or `"none"`.
#' @param label_digits Integer number of decimals used for percentage labels.
#' @param base_size Numeric base font size for the plot.
#' @param colors Named vector giving fill colors for `"Up"` and `"Down"` slices.
#'   When `show_other = TRUE`, an `"Other"` entry may also be supplied.
#' @param show_other Logical; when `TRUE`, include non-DE genes as an `"Other"`
#'   slice using `other_color` unless overridden in `colors`.
#' @param other_color Fill colour used for the `"Other"` slice when
#'   `show_other = TRUE` and `colors` does not include `"Other"`.
#' @param text_color Colour used for donut label text.
#' @param facet_by Either `"comparison"` (default) to draw one donut per comparison,
#'   or `"none"` for a single donut.
#' @param ncol Optional number of columns when faceting.
#' @export
get_deg_count_donutplot <- function(x,
                                    sample_comparisons = NULL,
                                    label = c("both", "count", "percent", "none"),
                                    label_digits = 1,
                                    base_size = 12,
                                    colors = c(Up = "red4", Down = "blue4"),
                                    show_other = FALSE,
                                    other_color = "grey70",
                                    text_color = "black",
                                    facet_by = c("comparison", "none"),
                                    ncol = NULL) {
  stopifnot(inherits(x, "VISTA"))
  df <- .collect_deg_count_data(
    x = x,
    sample_comparisons = sample_comparisons,
    show_other = show_other
  )
  if (is.null(df)) return(invisible(NULL))
  colors <- .normalize_deg_colors(
    colors = colors,
    include_other = show_other,
    other_color = other_color
  )
  .plot_deg_count_circular(
    df = df,
    label = label,
    label_digits = label_digits,
    base_size = base_size,
    colors = colors,
    text_color = text_color,
    donut = TRUE,
    facet_by = facet_by,
    ncol = ncol
  )
}


# ──────────────────────────────────────────────────────────────────────────────
# MA plot
# ──────────────────────────────────────────────────────────────────────────────

#' Generate MA plot from a VISTA object
#'
#' Create an MA plot (log2 fold change vs mean expression) for a selected
#' comparison contained in a VISTA object. Genes are coloured by their
#' regulation class and the top results can be optionally labeled with gene IDs.
#'
#' @param x A \linkS4class{VISTA} object.
#' @param sample_comparison Character scalar naming the comparison to plot. Must match
#'   one of \code{names(comparisons(x))}.
#' @param point_size Numeric point size. Default: 1.2.
#' @param alpha Numeric transparency (0-1). Default: 0.6.
#' @param colors Named character vector of colors for \code{"Up"},
#'   \code{"Down"}, and \code{"Other"} genes.
#' @param label_n Integer number of genes to label. Default: 0.
#' @param label_size Text size for labels. Default: 3.
#' @param repair_genes Logical; if \code{TRUE}, attempt to shorten gene identifiers to
#'   symbols by stripping prefixes. Default: \code{FALSE}.
#' @param display_id Optional ID/column name to use for labels. If supplied and
#'   present in `rowData(x)`, those values are used; otherwise falls back to ID
#'   mapping.
#' @param display_from Optional source ID type for mapping (used when
#'   `display_id` is not found in `rowData`).
#' @param display_orgdb Optional `OrgDb` used for ID mapping when `display_id`
#'   is set but not found in `rowData`.
#'
#' @return A \code{\link[ggplot2]{ggplot}} MA plot.
#' @export
get_ma_plot <- function(x,
                        sample_comparison,
                        point_size = 1.2,
                        alpha = 0.6,
                        colors = c(Up = "#a40000", Down = "#16317d", Other = "gray70"),
                        label_n = 0,
                        label_size = 3,
                        repair_genes = FALSE,
                        display_id = NULL,
                        display_from = NULL,
                        display_orgdb = NULL) {

  stopifnot(inherits(x, "VISTA"))

  color_names <- names(colors)
  colors <- as.character(colors)
  if (is.null(color_names) && length(colors) == 3L) {
    color_names <- c("Up", "Down", "Other")
  }
  names(colors) <- color_names
  if (is.null(names(colors)) || any(!c("Up", "Down", "Other") %in% names(colors))) {
    cli::cli_abort("{.arg colors} must be a named vector with entries for {.val Up}, {.val Down}, and {.val Other}.")
  }

  `%||%` <- function(a, b) if (!is.null(a)) a else b
  to_num <- function(v) { if (is.numeric(v)) return(v); suppressWarnings(as.numeric(gsub(",", "", as.character(v), fixed = TRUE))) }

  # --- pull DE table
  comps <- comparisons(x)
  if (!sample_comparison %in% names(comps)) cli::cli_abort("Comparison {.val {sample_comparison}} not found in comparisons(x).")
  comp_tbl <- as.data.frame(comps[[sample_comparison]], stringsAsFactors = FALSE)

  # --- standardize id / fc / p
  if (!"gene_id" %in% names(comp_tbl)) {
    id_candidates <- c("gene","geneID","GeneID","symbol","SYMBOL","ensembl","ensembl_id","ID","id")
    hit <- id_candidates[id_candidates %in% names(comp_tbl)][1]
    if (!is.na(hit)) names(comp_tbl)[names(comp_tbl) == hit] <- "gene_id"
  }
  if (!"gene_id" %in% names(comp_tbl)) {
    rn <- rownames(comp_tbl)
    if (!is.null(rn) && length(rn) == nrow(comp_tbl) && all(nzchar(rn)) && !anyDuplicated(rn)) comp_tbl$gene_id <- rn
    else cli::cli_abort("DE table lacks a {.field gene_id} column and valid rownames.")
  }
  comp_tbl$gene_id <- as.character(comp_tbl$gene_id)
  rd <- tryCatch(SummarizedExperiment::rowData(x), error = function(e) NULL)
  if (!is.null(display_id) && !is.null(rd) && display_id %in% colnames(rd)) {
    map <- rd[[display_id]]
    names(map) <- rownames(x)
    mapped <- map[match(comp_tbl$gene_id, names(map))]
    comp_tbl$display_gene <- ifelse(!is.na(mapped) & nzchar(mapped), mapped, comp_tbl$gene_id)
  } else {
    comp_tbl$display_gene <- .map_gene_ids(comp_tbl$gene_id,
                                           from_type = display_from,
                                           to_type = display_id,
                                           orgdb = display_orgdb)
  }

  if (!"log2fc" %in% names(comp_tbl)) {
    if ("log2FoldChange" %in% names(comp_tbl)) names(comp_tbl)[names(comp_tbl) == "log2FoldChange"] <- "log2fc"
    else if ("logFC" %in% names(comp_tbl))     names(comp_tbl)[names(comp_tbl) == "logFC"] <- "log2fc"
    else cli::cli_abort("DE table is missing a log2FC column.")
  }
  comp_tbl$log2fc <- to_num(comp_tbl$log2fc)

  p_col <- dplyr::first(intersect(c("pvalue","padj","qvalue","svalue","P.Value","p_val"), names(comp_tbl))) %||% NA_character_
  if (!is.na(p_col)) {
    comp_tbl[[p_col]] <- to_num(comp_tbl[[p_col]])
    comp_tbl[[p_col]][!is.finite(comp_tbl[[p_col]]) | is.na(comp_tbl[[p_col]])] <- 1
  }

  # regulation if missing
  if (!"regulation" %in% names(comp_tbl)) {
    cuts <- S4Vectors::metadata(x)$cutoffs %||% list()
    lfc_cut <- cuts$log2fc %||% 1
    p_cut   <- cuts$pvalue %||% 0.05
    if (!is.na(p_col)) {
      comp_tbl$regulation <- ifelse(comp_tbl$log2fc >=  lfc_cut & comp_tbl[[p_col]] <= p_cut, "Up",
                                    ifelse(comp_tbl$log2fc <= -lfc_cut & comp_tbl[[p_col]] <= p_cut, "Down", "Other"))
    } else {
      comp_tbl$regulation <- "Other"
    }
  }
  comp_tbl$regulation <- factor(comp_tbl$regulation, levels = c("Up","Down","Other"))

  # --- baseMean: take from rowData if present, else compute from assay
  rd_df <- as.data.frame(SummarizedExperiment::rowData(x), stringsAsFactors = FALSE)
  if (!"gene_id" %in% names(rd_df)) rd_df$gene_id <- rownames(rd_df)

  if ("baseMean" %in% names(rd_df)) {
    bm_df <- rd_df[, c("gene_id","baseMean"), drop = FALSE]
    bm_df$baseMean <- to_num(bm_df$baseMean)
  } else {
    a <- SummarizedExperiment::assay(x)
    bm_df <- data.frame(gene_id = rownames(a), baseMean = rowMeans(a, na.rm = TRUE), stringsAsFactors = FALSE)
  }

  # If DE already has a baseMean column, keep it but prefer bm_df via coalesce
  has_bm <- "baseMean" %in% names(comp_tbl)
  plot_df <- comp_tbl |>
    dplyr::left_join(bm_df, by = "gene_id", suffix = c("", ".bm"))

  if (has_bm) {
    plot_df <- plot_df |>
      dplyr::mutate(
        baseMean = dplyr::coalesce(.data$baseMean, .data$baseMean.bm)
      ) |>
      dplyr::select(-dplyr::any_of("baseMean.bm"))
  }

  plot_df <- plot_df |>
    dplyr::mutate(
      log10_baseMean = log10(.data$baseMean + 1),
      label = display_gene
    )

  if (repair_genes) plot_df$label <- sub(".*:", "", plot_df$label)

  # top labels
  label_df <- NULL
  if (label_n > 0) {
    if (!is.na(p_col)) {
      label_df <- plot_df |>
        dplyr::filter(regulation %in% c("Up","Down")) |>
        dplyr::arrange(.data[[p_col]], dplyr::desc(abs(log2fc))) |>
        dplyr::slice_head(n = label_n)
    } else {
      label_df <- plot_df |>
        dplyr::filter(regulation %in% c("Up","Down")) |>
        dplyr::arrange(dplyr::desc(abs(log2fc))) |>
        dplyr::slice_head(n = label_n)
    }
  }

  # plot
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = log10_baseMean, y = log2fc, color = regulation)) +
    ggplot2::geom_point(alpha = alpha, size = point_size) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(
      x = "log10(Mean Expression + 1)",
      y = "log2 Fold Change",
      title = paste("MA Plot:", sample_comparison),
      color = "Regulation"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right")

  if (!is.null(label_df)) {
    p <- p + ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(label = label),
      size = label_size,
      max.overlaps = Inf
    )
  }
  p
}



# ──────────────────────────────────────────────────────────────────────────────
#' @title Fold-change line plot across comparisons
#' @description Plots log2 fold-change trajectories for selected genes across multiple comparisons,
#' optionally clustering genes.
#' @param x A `VISTA` object containing differential expression results.
#' @param sample_comparisons Character vector of comparison names to include.
#' @param genes Optional character vector of gene identifiers to plot. Defaults to all genes.
#' @param display_id Optional ID/column name used to interpret `genes` and map
#'   the returned cluster table to display-friendly labels.
#' @param display_from Optional source ID type for mapping when `display_id` is
#'   not present in `rowData(x)`.
#' @param display_orgdb Optional `OrgDb` object used for identifier mapping when
#'   `display_id` is not present in `rowData(x)`.
#' @param km Optional integer specifying the number of k-means clusters to compute; `NULL` disables clustering.
#' @param facet_by Faceting mode: `"none"` (default) or `"cluster"` when
#'   k-means clustering is requested.
#' @param facet_scales Facet scales argument passed to `facet_wrap()` when
#'   `facet_by = "cluster"` (default `"fixed"`).
#' @param facet_nrow,facet_ncol Optional layout passed to `facet_wrap()` when faceting.
#' @param alpha Numeric alpha applied to individual gene lines.
#' @param palette Optional named or unnamed color vector used for cluster lines.
#' @param show_summary Logical; overlay a summary line per cluster when `TRUE`.
#' @param summary_color Color used for the summary line. When `NULL`, uses
#'   the first comparison color (if stored) for consistency across plots.
#' @param summary_linewidth Numeric line width for the summary line.
#' @param summary_fun Character string selecting `"median"` or `"mean"` for the summary statistic.
#' @param base_size Numeric base theme size.
#' @return A list with `plot` (the `ggplot2` object) and `clustered_data`
#'   (gene-to-cluster assignments).
# ──────────────────────────────────────────────────────────────────────────────

#' @export
get_foldchange_lineplot <- function(x,
                                    sample_comparisons,
                                    genes = NULL,
                                    display_id = NULL,
                                    display_from = NULL,
                                    display_orgdb = NULL,
                                    km = NULL,
                                    facet_by = c("none", "cluster"),
                                    facet_scales = "fixed",
                                    facet_nrow = NULL,
                                    facet_ncol = NULL,
                                    alpha = 0.5,
                                    palette = NULL,
                                    show_summary = TRUE,
                                    summary_color = NULL,
                                    summary_linewidth = 1,
                                    summary_fun = c("median", "mean"),
                                    base_size = 14) {

  stopifnot(inherits(x, "VISTA"))
  stopifnot(is.character(sample_comparisons))
  facet_by <- match.arg(facet_by)
  summary_fun <- match.arg(summary_fun)

  comps <- .vista_comparisons(x)
  if (!all(sample_comparisons %in% names(comps))) {
    cli::cli_abort("Some requested comparisons are not found in `metadata(x)$de_results`.")
  }

  if (is.null(summary_color)) {
    pal <- .vista_comparison_colors(x, sample_comparisons)
    summary_color <- pal[[1]] %||% "black"
  }

  comps_df <- purrr::map(comps[sample_comparisons], as.data.frame)
  res <- dplyr::bind_rows(comps_df, .id = "comparison")
  gene_id_col <- "gene_id"

  genes_use <- .resolve_foldchange_gene_ids(
    x = x,
    genes = genes,
    display_id = display_id,
    display_from = display_from,
    display_orgdb = display_orgdb
  )
  if (is.null(genes_use)) {
    genes_use <- unique(res[[gene_id_col]])
  }
  res <- dplyr::filter(res, .data[[gene_id_col]] %in% genes_use)
  if (!nrow(res)) {
    cli::cli_abort("No genes remain after filtering the requested comparisons.")
  }
  res$display_gene <- .resolve_foldchange_gene_labels(
    x = x,
    gene_ids = res[[gene_id_col]],
    display_id = display_id,
    display_from = display_from,
    display_orgdb = display_orgdb
  )
  res$comparison <- factor(res$comparison, levels = sample_comparisons)

  if (!is.null(km)) {
    # wide matrix genes x comparisons of log2fc then kmeans
    wide <- tidyr::pivot_wider(res[, c(gene_id_col, "comparison", "log2fc")],
                               names_from = comparison, values_from = log2fc)
    mat <- as.matrix(wide[,-1])
    km_fit <- stats::kmeans(scale(mat), centers = km, nstart = 10)
    cluster_df <- tibble::tibble(gene_id = wide[[gene_id_col]], cluster = factor(km_fit$cluster))
    res <- dplyr::left_join(res, cluster_df, by = gene_id_col)
  } else {
    res$cluster <- factor("All")
  }

  p <- ggplot2::ggplot(res, ggplot2::aes(x = comparison, y = log2fc, group = .data[[gene_id_col]], color = cluster)) +
    ggplot2::geom_line(alpha = alpha) +
    ggplot2::labs(x = "Comparisons", y = "Log2FC", color = "Cluster") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (!is.null(palette)) {
    if (is.null(names(palette))) {
      names(palette) <- levels(res$cluster)[seq_len(min(length(palette), nlevels(res$cluster)))]
    }
    p <- p + ggplot2::scale_color_manual(values = palette)
  }

  if (show_summary) {
    summary_fn <- switch(summary_fun,
                         median = stats::median,
                         mean = mean)
    summary_df <- res |>
      dplyr::group_by(cluster, comparison) |>
      dplyr::summarise(avg_log2FC = summary_fn(log2fc, na.rm = TRUE), .groups = "drop")

    p <- p + ggplot2::geom_line(data = summary_df,
                                ggplot2::aes(x = comparison, y = avg_log2FC, group = cluster),
                                color = summary_color,
                                linewidth = summary_linewidth)
  }

  if (facet_by == "cluster" && !is.null(km)) {
    cluster_counts <- res |>
      dplyr::distinct(display_gene, cluster) |>
      dplyr::count(cluster) |>
      dplyr::mutate(label = paste0("Cluster ", cluster, "\n(n = ", n, ")")) |>
      dplyr::select(cluster, label) |>
      tibble::deframe()

    p <- .add_expression_facet_wrap(
      plot = p,
      facet_formula = ~cluster,
      scales = facet_scales,
      facet_nrow = facet_nrow,
      facet_ncol = facet_ncol
    )
    p$facet$params$labeller <- ggplot2::labeller(cluster = cluster_counts)
  }

  list(
    plot = p,
    clustered_data = res |>
      dplyr::select(dplyr::all_of(gene_id_col), display_gene, cluster) |>
      dplyr::distinct()
  )
}


# ──────────────────────────────────────────────────────────────────────────────
# Alluvial plot of regulation transitions
# ──────────────────────────────────────────────────────────────────────────────

#' Plot alluvial diagram showing gene regulation transitions across comparisons
#'
#' @param x A VISTA object (DE results in `comparisons(x)`).
#' @param sample_comparisons Character vector of comparison names to include (>= 2).
#' @param show_other Logical; include "Other" genes. Default FALSE.
#' @return A ggplot object.
#' @export
get_deg_alluvial <- function(x,
                             sample_comparisons,
                             show_other = FALSE) {

  stopifnot(inherits(x, "VISTA"))
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE)) {
    cli::cli_abort("Packages {.pkg ComplexHeatmap} and {.pkg circlize} must be installed to draw heatmaps.")
  }
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggalluvial} must be installed to use `get_deg_alluvial()`.")
  }
  `%or%` <- function(a, b) if (!is.null(a)) a else b
  to_num <- function(v) { if (is.numeric(v)) return(v); suppressWarnings(as.numeric(gsub(",", "", as.character(v), fixed = TRUE))) }

  comps <- comparisons(x)
  if (is.null(comps) || !length(comps)) cli::cli_abort("No DE results found in comparisons(x).")
  if (!all(sample_comparisons %in% names(comps))) cli::cli_abort("Unknown comparisons: {.val {setdiff(sample_comparisons, names(comps))}}")
  if (length(sample_comparisons) < 2) cli::cli_abort("At least two comparisons required.")

  # cutoffs to derive regulation if it's missing
  cuts <- cutoffs(x) %or% list()
  lfc_cut <- cuts$log2fc %or% 1
  p_cut   <- cuts$pval %or% 0.05

  # standardize one DE table -> gene_id, log2fc, regulation
  .std_comp <- function(df) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)

    # gene_id
    if (!"gene_id" %in% names(df)) {
      id_candidates <- c("gene","geneID","GeneID","symbol","SYMBOL","ensembl","ensembl_id","ID","id")
      hit <- id_candidates[id_candidates %in% names(df)][1]
      if (!is.na(hit)) names(df)[names(df) == hit] <- "gene_id"
    }
    if (!"gene_id" %in% names(df)) {
      rn <- rownames(df)
      if (!is.null(rn) && length(rn) == nrow(df) && all(nzchar(rn)) && !anyDuplicated(rn)) {
        df$gene_id <- rn
      } else {
        cli::cli_abort("DE table lacks a {.field gene_id} column and valid rownames.")
      }
    }
    df$gene_id <- as.character(df$gene_id)

    # log2fc
    if (!"log2fc" %in% names(df)) {
      if ("log2FoldChange" %in% names(df))      names(df)[names(df) == "log2FoldChange"] <- "log2fc"
      else if ("logFC" %in% names(df))          names(df)[names(df) == "logFC"] <- "log2fc"
      else if ("LFC" %in% names(df))            names(df)[names(df) == "LFC"]   <- "log2fc"
    }
    if ("log2fc" %in% names(df)) df$log2fc <- to_num(df$log2fc)

    # regulation
    if (!"regulation" %in% names(df)) {
      p_col <- dplyr::first(intersect(c("pvalue","padj","qvalue","svalue","P.Value","p_val"), names(df)))
      if (!is.na(p_col)) {
        df[[p_col]] <- to_num(df[[p_col]])
        df[[p_col]][!is.finite(df[[p_col]]) | is.na(df[[p_col]])] <- 1
        df$regulation <- ifelse(df$log2fc >=  lfc_cut & df[[p_col]] <= p_cut, "Up",
                                ifelse(df$log2fc <= -lfc_cut & df[[p_col]] <= p_cut, "Down", "Other"))
      } else {
        df$regulation <- "Other"
      }
    }
    df$regulation <- factor(df$regulation, levels = c("Up","Down","Other"))

    df[, c("gene_id", intersect("log2fc", names(df)), "regulation"), drop = FALSE]
  }

  # 1) Extract + standardize
  full_deg_list <- purrr::map(sample_comparisons, ~ .std_comp(comps[[.x]]) |>
                                dplyr::mutate(comparison = .x))
  names(full_deg_list) <- sample_comparisons

  # 2) Intersect gene IDs across comparisons
  common_gene_ids <- purrr::map(full_deg_list, ~ .x$gene_id) |> purrr::reduce(intersect)

  # 3) Long table (+ optional removal of "Other")
  deg_long <- purrr::map2_dfr(full_deg_list, sample_comparisons, function(df, comp) {
    tmp <- dplyr::filter(df, .data$gene_id %in% common_gene_ids)
    if (!show_other) tmp <- dplyr::filter(tmp, .data$regulation %in% c("Up","Down"))
    dplyr::mutate(tmp, comparison = comp)
  })

  if (!nrow(deg_long)) cli::cli_abort("No genes remain after filtering; consider `show_other = TRUE`.")

  # 4) Wide (one column per comparison) → ensure factors & no NAs
  deg_wide <- deg_long |>
    dplyr::select(gene_id, comparison, regulation) |>
    tidyr::pivot_wider(names_from = comparison, values_from = regulation, values_fill = "Other") |>
    dplyr::distinct() |>
    dplyr::mutate(dplyr::across(tidyselect::all_of(sample_comparisons),
                                ~factor(.x, levels = c("Up","Down","Other"))))

  # 5) Back to long for plotting
  deg_long_plot <- deg_wide |>
    tidyr::pivot_longer(-gene_id, names_to = "name", values_to = "value") |>
    dplyr::mutate(
      value = forcats::fct_relevel(value, c("Up","Down","Other")),
      name  = forcats::fct_relevel(name, sample_comparisons)
    )

  # 6) Plot
  ggplot2::ggplot(
    deg_long_plot,
    ggplot2::aes(x = name, alluvium = gene_id, stratum = value, fill = value)
  ) +
    ggalluvial::geom_flow(width = 1/4) +
    ggalluvial::geom_stratum(alpha = .5, width = 1/4) +
    ggplot2::geom_text(
      ggplot2::aes(label = paste(..stratum.., "\n(n=", ..count.., ")", sep = "")),
      stat = ggalluvial::StatStratum,
      size = 4
    ) +
    ggplot2::scale_fill_manual(
      values = c(Up = "#D73027", Down = "#1A9850", Other = "#4575B4"),
      name = "Regulation",
      na.value = "grey80"
    ) +
    ggplot2::scale_x_discrete(expand = c(0.2, 0.05)) +
    ggplot2::ylab("Number of Genes") +
    ggplot2::ggtitle("Gene Regulation Transitions Across Comparisons") +
    ggplot2::theme_minimal()
}



# ──────────────────────────────────────────────────────────────────────────────
# Expression heatmap helpers
# ──────────────────────────────────────────────────────────────────────────────

.validate_heatmap_gene_input <- function(genes, arg = "genes") {
  if (is.null(genes)) {
    return(NULL)
  }
  if (!is.character(genes)) {
    cli::cli_abort("{.arg {arg}} must be a character vector of gene identifiers.")
  }
  genes <- unique(stats::na.omit(as.character(genes)))
  genes <- genes[nzchar(genes)]
  if (!length(genes)) {
    cli::cli_abort("{.arg {arg}} must contain at least one non-empty gene identifier.")
  }
  genes
}

.auto_expression_heatmap_genes <- function(x, sample_ids, top_n = 50L) {
  top_n <- as.integer(top_n)[1]
  if (!is.finite(top_n) || is.na(top_n) || top_n < 1) {
    cli::cli_abort("{.arg top_n} must be a positive integer when {.arg genes} is omitted.")
  }

  mat <- SummarizedExperiment::assay(x)[, sample_ids, drop = FALSE]
  vars <- matrixStats::rowVars(as.matrix(mat), na.rm = TRUE)
  means <- rowMeans(mat, na.rm = TRUE)
  ord <- order(vars, means, decreasing = TRUE, na.last = NA)
  rownames(mat)[utils::head(ord, n = min(top_n, nrow(mat)))]
}

.auto_foldchange_heatmap_genes <- function(x, sample_comparisons, top_n = 10L) {
  top_n <- as.integer(top_n)[1]
  if (!is.finite(top_n) || is.na(top_n) || top_n < 1) {
    cli::cli_abort("{.arg top_n} must be a positive integer when {.arg genes} is omitted.")
  }

  comps <- .vista_comparisons(x)
  picked <- character()

  for (comp_name in sample_comparisons) {
    tbl <- .tidy_de_results(comps[[comp_name]], rowname_col = "gene_id")
    if (!nrow(tbl)) {
      next
    }
    tbl <- tbl[!is.na(tbl$gene_id) & nzchar(tbl$gene_id), , drop = FALSE]
    if (!nrow(tbl)) {
      next
    }
    if ("regulation" %in% colnames(tbl)) {
      sig_tbl <- tbl[tbl$regulation %in% c("Up", "Down"), , drop = FALSE]
      if (nrow(sig_tbl)) {
        tbl <- sig_tbl
      }
    }
    ord <- order(abs(.safe_numeric(tbl$log2fc, default = NA_real_)), decreasing = TRUE, na.last = NA)
    picked <- c(picked, tbl$gene_id[utils::head(ord, n = min(top_n, length(ord)))])
  }

  picked <- unique(stats::na.omit(as.character(picked)))
  picked <- picked[nzchar(picked)]

  if (length(picked)) {
    return(picked)
  }

  fallback <- rownames(x)
  fallback <- fallback[!is.na(fallback) & nzchar(fallback)]
  utils::head(fallback, n = min(top_n, length(fallback)))
}

.resolve_heatmap_row_labels <- function(x,
                                        gene_ids,
                                        display_id = NULL,
                                        display_from = NULL,
                                        display_orgdb = NULL,
                                        repair_genes = FALSE,
                                        prefer_symbol = FALSE) {
  rd <- tryCatch(SummarizedExperiment::rowData(x), error = function(e) NULL)
  display_id_use <- display_id

  if (isTRUE(prefer_symbol) && is.null(display_id_use) && !is.null(rd) && "SYMBOL" %in% colnames(rd)) {
    display_id_use <- "SYMBOL"
  }

  labels <- gene_ids
  if (!is.null(display_id_use) && !is.null(rd) && display_id_use %in% colnames(rd)) {
    map <- rd[[display_id_use]]
    names(map) <- rownames(x)
    mapped <- map[match(gene_ids, names(map))]
    labels <- ifelse(!is.na(mapped) & nzchar(mapped), mapped, labels)
  } else if (!is.null(display_id_use)) {
    labels <- .map_gene_ids(
      gene_ids,
      from_type = display_from,
      to_type = display_id_use,
      orgdb = display_orgdb
    )
  }

  if (repair_genes || isTRUE(prefer_symbol)) {
    labels <- gsub(".*?:", "", labels)
  }

  unname(labels)
}

.resolve_heatmap_annotation_colors <- function(col_df,
                                               column_anno_palette = "Dark 3",
                                               column_anno_colors = NULL) {
  if (is.null(column_anno_colors)) {
    column_anno_colors <- list()
  }
  if (!is.list(column_anno_colors)) {
    cli::cli_abort("{.arg column_anno_colors} must be NULL or a named list of color vectors.")
  }

  anno_colors <- list()
  for (col_nm in colnames(col_df)) {
    col_levels <- unique(as.character(col_df[[col_nm]]))
    col_df[[col_nm]] <- factor(as.character(col_df[[col_nm]]), levels = col_levels)

    manual_cols <- column_anno_colors[[col_nm]]
    if (!is.null(manual_cols)) {
      manual_names <- names(manual_cols)
      manual_cols <- as.character(manual_cols)
      names(manual_cols) <- manual_names
      if (is.null(names(manual_cols)) || any(!nzchar(names(manual_cols)))) {
        cli::cli_abort(
          c(
            "{.arg column_anno_colors[['{col_nm}']]} must be a named vector.",
            "i" = "Names should match the annotation levels for {.val {col_nm}}."
          )
        )
      }
      missing_levels <- setdiff(col_levels, names(manual_cols))
      if (length(missing_levels) > 0) {
        auto_cols <- colorspace::qualitative_hcl(length(missing_levels), palette = column_anno_palette)
        names(auto_cols) <- missing_levels
        manual_cols <- c(manual_cols, auto_cols)
      }
      anno_colors[[col_nm]] <- manual_cols[col_levels]
    } else {
      level_colors <- colorspace::qualitative_hcl(length(col_levels), palette = column_anno_palette)
      names(level_colors) <- col_levels
      anno_colors[[col_nm]] <- level_colors
    }
  }

  list(col_df = col_df, anno_colors = anno_colors)
}

#' @title Expression heatmap
#' @description Summarizes expression for selected genes/groups via ComplexHeatmap
#'   with optional transformations and annotations. With only a `VISTA`
#'   object, the function will plot the top variable genes across all samples.
#' @aliases get_expression_heatmap
#' @param x A `VISTA` object.
#' @param sample_group Character vector of group labels specifying which samples
#'   to include (based on the selected grouping column).
#' @param genes Optional character vector of gene identifiers to display. When
#'   omitted, VISTA selects the top variable genes across the included samples.
#' @param top_n Integer number of genes to select automatically when
#'   `genes = NULL`. Defaults to `50`.
#' @param group_column Optional column name in `sample_info` used to interpret `sample_group`.
#' @param value_transform One of `"zscore"`, `"log2"`, or `"raw"` controlling how expression values are transformed.
#' @param summarise_replicates Logical; average replicates per group before plotting when `TRUE`.
#' @param summarise_method `"mean"` or `"median"` summarization used when `summarise_replicates = TRUE`.
#' @param convert_rowmeans Logical; subtract row means prior to plotting.
#' @param display_id Optional ID/column name to use for row labels. If supplied
#' @param display_from Optional source ID type for mapping (used when `display_id`
#' @param display_orgdb Optional `OrgDb` object used for ID mapping when
#' @param repair_genes Logical; if `TRUE`, split `gene_id` strings such as `ID:SYMBOL` to display the symbol.
#' @param show_row_names Logical; display row names (genes) beside the heatmap.
#'   When `NULL`, VISTA turns labels on automatically for auto-selected genes.
#' @param label_size Numeric font size for row names.
#' @param label_specific_rows Optional character vector of row names to highlight via `anno_mark()`.
#' @param label_specific_rows_gp `grid::gpar()` object controlling the highlighted row labels.
#' @param show_column_names Logical; draw column names when `TRUE`.
#' @param cluster_rows Logical; cluster rows when drawing the heatmap.
#' @param show_row_dend Logical; display the row dendrogram.
#' @param cluster_columns Logical; cluster columns.
#' @param kmeans_k Optional integer specifying the number of k-means clusters to compute for rows.
#' @param annotate_columns Logical or character vector. `TRUE` adds the
#'   `group_column` annotation and also includes `cluster_by` when supplied; a
#'   character vector adds multiple annotations from `sample_info`.
#' @param cluster_by Optional annotation column used to split/cluster columns.
#'   Defaults to the first annotation column when `annotate_columns` is enabled.
#' @param column_anno_palette Qualitative palette name used for the column annotation.
#' @param column_anno_colors Optional named list of annotation color vectors.
#'   Each element should be a named character vector keyed by the levels of one
#'   annotation column.
#' @param color_default Logical; use the default blue-white-red palette when `TRUE`. Set to `FALSE` to supply `col`.
#' @param col Optional `circlize::colorRamp2` function used when `color_default = FALSE`.
#' @param heatmap_name Optional legend title.
#' @param show_heatmap_legend Logical; display the heatmap legend.
#' @param return_type `"heatmap"`, `"clusters"`, or `"both"` selecting the returned object.
#' @return A `ComplexHeatmap` object, a cluster data frame, or a list containing
#'   both depending on `return_type`.
#' @examples
#' v <- example_vista()
#' if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
#'     requireNamespace("circlize", quietly = TRUE)) {
#'   hm <- get_expression_heatmap(v, return_type = "heatmap")
#'   ComplexHeatmap::draw(hm)
#' }
#' @param ... Additional arguments passed to `ComplexHeatmap::Heatmap()`.
#' @export
get_expression_heatmap <- function(x,
                                   genes = NULL,
                                   top_n = 50,
                                   sample_group = NULL,
                                   group_column = NULL,
                                   value_transform = c("zscore", "log2", "raw"),
                                   summarise_replicates = TRUE,
                                   summarise_method = c("mean", "median"),
                                   convert_rowmeans = FALSE,
                                   display_id = NULL,
                                   display_from = NULL,
                                   display_orgdb = NULL,
                                   repair_genes = FALSE,
                                   show_row_names = NULL,
                                   label_size = 10,
                                   label_specific_rows = NULL,
                                   label_specific_rows_gp = grid::gpar(fontsize = 5),
                                   show_column_names = TRUE,
                                   cluster_rows = TRUE,
                                   show_row_dend = TRUE,
                                   cluster_columns = TRUE,
                                   kmeans_k = NULL,
                                   annotate_columns = FALSE,
                                   cluster_by = NULL,
                                   column_anno_palette = "Dark 3",
                                   column_anno_colors = NULL,
                                   color_default = TRUE,
                                   col = NULL,
                                   heatmap_name = NULL,
                                   show_heatmap_legend = TRUE,
                                   return_type = c("heatmap", "clusters", "both"),
                                   ...) {

  stopifnot(inherits(x, "VISTA"))
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE)) {
    cli::cli_abort("Packages {.pkg ComplexHeatmap} and {.pkg circlize} must be installed to draw heatmaps.")
  }
  value_transform <- match.arg(value_transform)
  summarise_method <- match.arg(summarise_method)
  return_type <- match.arg(return_type)

  group_col <- group_column %||% .vista_group_col(x)

  expr <- SummarizedExperiment::assay(x) |>
    as.data.frame() |>
    tibble::rownames_to_column("gene_id")

  sample_info <- as.data.frame(SummarizedExperiment::colData(x)) |>
    tibble::rownames_to_column("sample")

  if (!(group_col %in% names(sample_info))) {
    cli::cli_abort("{.arg group_column} '{group_col}' not found in {.arg sample_info}.")
  }

  annotate_cols <- character()
  if (isTRUE(annotate_columns)) {
    annotate_cols <- unique(c(group_col, cluster_by))
    annotate_cols <- annotate_cols[!is.na(annotate_cols) & nzchar(annotate_cols)]
  } else if (is.character(annotate_columns) && length(annotate_columns) > 0) {
    annotate_cols <- unique(annotate_columns)
  } else if (!isFALSE(annotate_columns)) {
    cli::cli_abort("{.arg annotate_columns} must be TRUE, FALSE, or a non-empty character vector.")
  }

  if (length(annotate_cols) > 0) {
    missing_anno_cols <- setdiff(annotate_cols, colnames(sample_info))
    if (length(missing_anno_cols) > 0) {
      cli::cli_abort(
        c(
          "Unknown columns requested in {.arg annotate_columns}.",
          "x" = "Missing: {.val {missing_anno_cols}}"
        )
      )
    }
  }
  if (is.null(sample_group)) {
    sample_group <- unique(sample_info[[group_col]])
  }
  sample_group <- unique(as.character(sample_group))

  missing_groups <- setdiff(sample_group, unique(sample_info[[group_col]]))
  if (length(missing_groups) > 0) {
    cli::cli_abort(
      c(
        "Unknown group labels in {.arg sample_group}.",
        "x" = "Missing: {.val {missing_groups}}",
        "i" = "Available: {.val {unique(sample_info[[group_col]])}}"
      )
    )
  }

  sample_ids <- unlist(purrr::map(sample_group, function(grp) {
    sample_info$sample[sample_info[[group_col]] == grp]
  }))

  auto_genes <- is.null(genes)
  genes <- .validate_heatmap_gene_input(genes)
  if (auto_genes) {
    genes <- .auto_expression_heatmap_genes(x, sample_ids = sample_ids, top_n = top_n)
  }

  expr <- expr |>
    dplyr::select(gene_id, dplyr::all_of(sample_ids)) |>
    dplyr::filter(gene_id %in% genes)

  missing_genes <- setdiff(genes, expr$gene_id)
  if (length(missing_genes) > 0) {
    cli::cli_abort(
      c(
        "Some {.arg genes} were not found in the VISTA object.",
        "x" = "Missing: {.val {missing_genes}}"
      )
    )
  }

  mat <- expr |> tibble::column_to_rownames("gene_id") |> as.matrix()
  mat <- mat[match(genes, rownames(mat)), , drop = FALSE]

  if (is.null(show_row_names)) {
    show_row_names <- auto_genes
  }
  show_row_names <- isTRUE(show_row_names)

  rd <- tryCatch(SummarizedExperiment::rowData(x), error = function(e) NULL)
  rownames(mat) <- .resolve_heatmap_row_labels(
    x = x,
    gene_ids = rownames(mat),
    display_id = display_id,
    display_from = display_from,
    display_orgdb = display_orgdb,
    repair_genes = repair_genes,
    prefer_symbol = auto_genes
  )

  if (value_transform == "log2") {
    mat <- log2(mat + 1)
  } else if (value_transform == "zscore") {
    mat <- t(scale(t(mat)))
    mat[is.na(mat)] <- 0
  }

  if (convert_rowmeans) mat <- mat - rowMeans(mat, na.rm = TRUE)

  if (summarise_replicates) {
    groups <- sample_info |>
      dplyr::filter(sample %in% colnames(mat)) |>
      dplyr::select(sample, !!group_col)

    summarised_list <- lapply(unique(groups[[group_col]]), function(grp) {
      samples_grp <- groups$sample[groups[[group_col]] == grp]
      if (summarise_method == "mean") {
        rowMeans(mat[, samples_grp, drop = FALSE], na.rm = TRUE)
      } else {
        matrixStats::rowMedians(as.matrix(mat[, samples_grp, drop = FALSE]), na.rm = TRUE)
      }
    })

    mat <- do.call(cbind, summarised_list)
    colnames(mat) <- unique(groups[[group_col]])
  } else {
    mat <- mat[, sample_ids, drop = FALSE]
  }

  if (repair_genes && !auto_genes) {
    rownames(mat) <- gsub(".*?:", "", rownames(mat))
    if (!is.null(label_specific_rows)) label_specific_rows <- gsub(".*?:", "", label_specific_rows)
  }

  if (!is.null(label_specific_rows) && show_row_names) {
    warning("`show_row_names = TRUE` is incompatible with `label_specific_rows`; setting `show_row_names = FALSE`.")
    show_row_names <- FALSE
  }

  if (color_default) {
    rng <- range(mat, na.rm = TRUE)
    max_abs <- max(abs(rng))
    if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
    col_fun <- circlize::colorRamp2(c(-max_abs, 0, max_abs), c("blue", "white", "red"))
  } else {
    col_fun <- col
  }

  if (is.null(heatmap_name)) {
    heatmap_name <- switch(value_transform,
                           zscore = "Z-score",
                           log2 = "Log2 Expression",
                           raw = "Raw Expression")
  }

  row_anno <- NULL
  if (!is.null(label_specific_rows)) {
    label_specific_rows <- .resolve_heatmap_row_labels(
      x = x,
      gene_ids = label_specific_rows,
      display_id = display_id,
      display_from = display_from,
      display_orgdb = display_orgdb,
      repair_genes = repair_genes,
      prefer_symbol = auto_genes
    )
    match_rows <- intersect(label_specific_rows, rownames(mat))
    if (length(match_rows) > 0) {
      row_indices <- which(rownames(mat) %in% match_rows)
      row_anno <- ComplexHeatmap::rowAnnotation(
        labels = ComplexHeatmap::anno_mark(
          at = row_indices,
          labels = rownames(mat)[row_indices],
          labels_gp = label_specific_rows_gp
        )
      )
    }
  }

  col_anno <- NULL
  column_split <- NULL
  if (length(annotate_cols) > 0) {
    if (summarise_replicates) {
      groups <- sample_info |>
        dplyr::filter(sample %in% sample_ids) |>
        dplyr::select(sample, dplyr::all_of(c(group_col, annotate_cols)))

      inconsistent_labels <- character(0)

      summarise_one <- function(x) {
        ux <- unique(as.character(stats::na.omit(x)))
        if (length(ux) == 0) return(NA_character_)
        ux[1]
      }

      col_df <- vector("list", length(annotate_cols))
      names(col_df) <- annotate_cols
      for (col_nm in annotate_cols) {
        per_group_vals <- vapply(colnames(mat), function(grp) {
          vals <- groups[[col_nm]][groups[[group_col]] == grp]
          summarise_one(vals)
        }, character(1))
        has_inconsistency <- any(vapply(colnames(mat), function(grp) {
          vals <- groups[[col_nm]][groups[[group_col]] == grp]
          length(unique(as.character(stats::na.omit(vals)))) > 1
        }, logical(1)))
        if (isTRUE(has_inconsistency)) {
          inconsistent_labels <- unique(c(inconsistent_labels, col_nm))
        }
        col_df[[col_nm]] <- per_group_vals
      }
      if (length(inconsistent_labels) > 0) {
        cli::cli_warn(
          c(
            "Some annotation columns vary within groups when {.arg summarise_replicates = TRUE}.",
            "i" = "Using the first non-missing value per group for: {.val {inconsistent_labels}}"
          )
        )
      }
      col_df <- as.data.frame(col_df, stringsAsFactors = FALSE)
      rownames(col_df) <- colnames(mat)
    } else {
      col_df <- SummarizedExperiment::colData(x) |>
        as.data.frame() |>
        tibble::rownames_to_column("sample") |>
        dplyr::filter(sample %in% colnames(mat)) |>
        dplyr::select(sample, dplyr::all_of(annotate_cols)) |>
        tibble::column_to_rownames("sample")
    }

    # Ensure annotation rows align with the heatmap column order
    col_df <- col_df[colnames(mat), , drop = FALSE]

    anno_res <- .resolve_heatmap_annotation_colors(
      col_df = col_df,
      column_anno_palette = column_anno_palette,
      column_anno_colors = column_anno_colors
    )
    col_df <- anno_res$col_df
    anno_colors <- anno_res$anno_colors

    split_col <- cluster_by %||% annotate_cols[[1]]
    if (!split_col %in% colnames(col_df)) {
      cli::cli_abort(
        c(
          "{.arg cluster_by} must be one of the active annotation columns.",
          "x" = "Received: {.val {split_col}}",
          "i" = "Available: {.val {colnames(col_df)}}"
        )
      )
    }
    column_split <- col_df[[split_col]]

    col_anno <- ComplexHeatmap::HeatmapAnnotation(
      df = as.data.frame(col_df),
      col = anno_colors
    )
  } else if (!is.null(cluster_by)) {
    cli::cli_abort("{.arg cluster_by} can only be used when {.arg annotate_columns} is enabled.")
  }

  heatmap_args <- c(
    list(
      matrix = as.matrix(mat),
      name = heatmap_name,
      col = col_fun,
      cluster_rows = cluster_rows,
      show_row_dend = show_row_dend,
      show_row_names = show_row_names,
      row_names_gp = grid::gpar(fontsize = label_size),
      show_column_names = show_column_names,
      cluster_columns = cluster_columns,
      show_heatmap_legend = show_heatmap_legend,
      row_km = if (!is.null(kmeans_k)) kmeans_k else NULL,
      row_km_repeats = 1,
      top_annotation = col_anno
    ),
    list(...)
  )

  if (!is.null(column_split) && is.null(heatmap_args$column_split)) {
    heatmap_args$column_split <- column_split
  }

  ht <- do.call(ComplexHeatmap::Heatmap, heatmap_args)

  if (!is.null(row_anno)) ht <-  ht + row_anno

  if (!is.null(kmeans_k) && return_type != "heatmap") {
    drawn <- ComplexHeatmap::draw(ht)
    clusters <- ComplexHeatmap::row_order(drawn)
    row_names <- rownames(mat)
    row_clusters <- rep(NA_integer_, length(row_names))
    for (k in seq_along(clusters)) row_clusters[clusters[[k]]] <- k
    cluster_df <- tibble::tibble(gene = row_names, cluster = row_clusters)
    if (return_type == "clusters") return(cluster_df)
    return(list(heatmap = ht, clusters = cluster_df))
  }

  if (return_type == "heatmap") return(ht)
  list(heatmap = ht, clusters = NULL)
}

# ──────────────────────────────────────────────────────────────────────────────
# Fold-change barplot (per gene across comparisons)
# ──────────────────────────────────────────────────────────────────────────────

#' Plot fold-change barplots across comparisons for selected genes
#'
#' Plots log2 fold changes for selected genes across one or more comparisons.
#' `facet_by` controls per-gene or per-comparison layout explicitly.
#'
#' @param x A `VISTA` object containing differential expression results.
#' @param genes Character vector of gene IDs to plot.
#' @param sample_comparisons Optional character vector of comparison names to include; defaults to all available.
#' @param coord_flip Logical; flip axes when `TRUE`.
#' @param display_id Optional column in `rowData(x)` to use for gene labels. Input
#'   gene matching still uses `gene_id`.
#' @param display_from Optional source ID type for mapping when `display_id` is
#'   not present in `rowData(x)`.
#' @param display_orgdb Optional `OrgDb` object used for identifier mapping when
#'   `display_id` is not present in `rowData(x)`.
#' @param sort_by How to order genes when faceting: `"input"` (use supplied order),
#'   `"log2fc"` (descending log2FC of the first comparison), or `"abs_log2fc"`
#'   (descending max absolute log2FC across comparisons).
#' @param facet_scales Facet scales argument passed to `facet_wrap()` when faceting
#'   (default `"free_y"`).
#' @param facet_nrow,facet_ncol Optional layout passed to `facet_wrap()` when faceting.
#' @param facet_by Faceting mode: `"auto"` (default), `"gene"`,
#'   `"comparison"`, or `"none"`.
#' @return A `ggplot2` object.
#'
#' @examples
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
#' genes <- rownames(vista)[seq_len(3)]
#' get_foldchange_barplot(vista, genes = genes)
#' get_foldchange_barplot(vista, genes = genes, facet_by = "gene")
#'
#' @export
get_foldchange_barplot <- function(x,
                                   genes,
                                   sample_comparisons = NULL,
                                   coord_flip = FALSE,
                                   display_id = NULL,
                                   display_from = NULL,
                                   display_orgdb = NULL,
                                   sort_by = c("input", "log2fc", "abs_log2fc"),
                                   facet_scales = "free_y",
                                   facet_nrow = NULL,
                                   facet_ncol = NULL,
                                   facet_by = c("auto", "gene", "comparison", "none")) {
  stopifnot(inherits(x, "VISTA"))
  stopifnot(length(genes) > 0)
  sort_by <- match.arg(sort_by)

  rd <- tryCatch(SummarizedExperiment::rowData(x), error = function(e) NULL)
  comps <- .vista_comparisons(x)
  if (is.null(sample_comparisons)) sample_comparisons <- names(comps)
  stopifnot(all(sample_comparisons %in% names(comps)))

  # map input genes via display_id if present
  input_genes <- .resolve_foldchange_gene_ids(
    x = x,
    genes = genes,
    display_id = display_id,
    display_from = display_from,
    display_orgdb = display_orgdb
  )

  fc_list <- purrr::map(sample_comparisons, \(nm) {
    df <- comps[[nm]]
    df[df$gene_id %in% input_genes, c("gene_id", "log2fc"), drop = FALSE] |>
      dplyr::mutate(comparison = nm)
  })
  df <- dplyr::bind_rows(fc_list)
  df <- dplyr::filter(df, !is.na(log2fc))
  if (!nrow(df)) cli::cli_abort("No genes found after filtering.")

  # gene labels
  df$gene_label <- .resolve_foldchange_gene_labels(
    x = x,
    gene_ids = df$gene_id,
    display_id = display_id,
    display_from = display_from,
    display_orgdb = display_orgdb
  )

  # ordering
  if (sort_by == "input" && !is.null(genes)) {
    gene_levels <- unique(df$gene_id[match(input_genes, df$gene_id)])
  } else if (sort_by == "log2fc") {
    ref <- dplyr::filter(df, .data$comparison == sample_comparisons[1])
    gene_levels <- ref$gene_id[order(ref$log2fc, decreasing = TRUE)]
  } else { # abs_log2fc
    ref <- df |>
      dplyr::group_by(.data$gene_id) |>
      dplyr::summarise(max_abs_fc = max(abs(.data$log2fc), na.rm = TRUE), .groups = "drop")
    gene_levels <- ref$gene_id[order(ref$max_abs_fc, decreasing = TRUE)]
  }
  gene_levels <- unique(gene_levels[!is.na(gene_levels)])
  df$gene_id <- factor(df$gene_id, levels = gene_levels)
  df$comparison <- factor(df$comparison, levels = sample_comparisons)
  facet_mode <- .resolve_foldchange_facet_mode(
    facet_by = facet_by,
    n_genes = length(gene_levels)
  )

  pal <- .vista_comparison_colors(x, sample_comparisons)
  if (is.null(pal)) {
    pal <- colorspace::qualitative_hcl(length(sample_comparisons), palette = "Dark 3")
  }
  if (is.null(names(pal))) names(pal) <- sample_comparisons

  if (facet_mode == "none") {
    plt <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = gene_id, y = log2fc, fill = comparison)
    ) +
      ggplot2::geom_col(
        alpha = 0.85,
        width = 0.7,
        position = ggplot2::position_dodge(width = 0.7)
      ) +
      ggplot2::scale_fill_manual(values = pal) +
      ggplot2::labs(x = "Gene", y = "Log2 Fold Change", fill = "Comparison") +
      ggplot2::scale_x_discrete(labels = stats::setNames(df$gene_label, df$gene_id)) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal()
    if (length(sample_comparisons) == 1) {
      plt <- plt + ggplot2::guides(fill = "none")
    }
  } else if (facet_mode == "comparison") {
    plt <- ggplot2::ggplot(df, ggplot2::aes(x = gene_id, y = log2fc, fill = comparison)) +
      ggplot2::geom_col(alpha = 0.85, width = 0.7, position = ggplot2::position_dodge(width = 0.7)) +
      ggplot2::scale_fill_manual(values = pal) +
      ggplot2::labs(x = "Gene", y = "Log2 Fold Change", fill = "Comparison") +
      ggplot2::theme_minimal() +
      ggplot2::scale_x_discrete(labels = stats::setNames(df$gene_label, df$gene_id)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    plt <- .add_expression_facet_wrap(
      plot = plt,
      facet_formula = ~comparison,
      scales = facet_scales,
      facet_nrow = facet_nrow,
      facet_ncol = facet_ncol
    )
  } else {
    plt <- ggplot2::ggplot(df, ggplot2::aes(x = comparison, y = log2fc, fill = comparison)) +
      ggplot2::geom_col(alpha = 0.85, width = 0.7, position = ggplot2::position_dodge(width = 0.7)) +
      ggplot2::scale_fill_manual(values = pal) +
      ggplot2::labs(x = "Comparison", y = "Log2 Fold Change", fill = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    if (length(levels(df$gene_id)) > 1) {
      plt <- .add_expression_facet_wrap(
        plot = plt,
        facet_formula = ~gene_id,
        scales = facet_scales,
        facet_nrow = facet_nrow,
        facet_ncol = facet_ncol
      )
      plt$facet$params$labeller <- ggplot2::as_labeller(
        stats::setNames(df$gene_label, as.character(df$gene_id))
      )
    }
  }

  if (coord_flip) plt <- plt + ggplot2::coord_flip()
  plt
}

# ──────────────────────────────────────────────────────────────────────────────
# Fold-change heatmap
# ──────────────────────────────────────────────────────────────────────────────

#' @title Fold-change heatmap
#' @description Visualizes log2 fold-change matrices across comparisons with
#'   ComplexHeatmap, supporting clustering and annotations. With only a
#'   `VISTA` object, the function will plot the top DE genes across the stored
#'   comparisons.
#' @param x A `VISTA` object with stored differential expression results.
#' @param sample_comparisons Optional character vector of comparison names to
#'   include. Defaults to all available comparisons.
#' @param genes Optional character vector of gene identifiers to display. When
#'   omitted, VISTA selects the top DE genes from each comparison by absolute
#'   log2 fold-change.
#' @param top_n Integer number of genes to select per comparison when
#'   `genes = NULL`. Defaults to `10`.
#' @param display_id Optional ID/column name to use for plot labels. If supplied
#' @param display_from Optional source ID type for mapping (used when `display_id`
#' @param display_orgdb Optional `OrgDb` object used for ID mapping when
#' @param repair_genes Logical; attempt to simplify `gene_id` strings by removing prefixes.
#' @param show_row_names Logical; draw row (gene) names. When `NULL`, VISTA
#'   turns labels on automatically for auto-selected genes.
#' @param label_size Numeric font size for row names.
#' @param label_specific_rows Optional character vector of genes to highlight with `anno_mark()`.
#' @param label_specific_rows_gp `grid::gpar()` object controlling highlighted labels.
#' @param show_column_names Logical; draw column labels.
#' @param cluster_rows Logical; cluster rows.
#' @param show_row_dend Logical; display the row dendrogram.
#' @param cluster_columns Logical; cluster columns.
#' @param kmeans_k Optional integer specifying the number of k-means clusters for rows.
#' @param annotate_columns Logical; add annotation bars keyed to the sample grouping column.
#' @param column_anno_palette Qualitative palette name used for column annotations.
#' @param color_default Logical; use the default diverging palette when `TRUE`. Set to `FALSE` to supply `col`.
#' @param col Optional `circlize::colorRamp2` color function used when `color_default = FALSE`.
#' @param heatmap_name Optional legend title.
#' @param show_heatmap_legend Logical; display the heatmap legend.
#' @param return_type `"heatmap"`, `"clusters"`, or `"both"` selecting the returned value.
#' @return A `ComplexHeatmap` object, a cluster data frame, or a list containing
#'   both depending on `return_type`.
#' @examples
#' v <- example_vista()
#' if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
#'     requireNamespace("circlize", quietly = TRUE)) {
#'   hm <- get_foldchange_heatmap(v, return_type = "heatmap")
#'   ComplexHeatmap::draw(hm)
#' }
#' @param ... Additional arguments forwarded to `ComplexHeatmap::Heatmap()`.
#' @export
get_foldchange_heatmap <- function(x,
                                   sample_comparisons = NULL,
                                   genes = NULL,
                                   top_n = 10,
                                   display_id = NULL,
                                   display_from = NULL,
                                   display_orgdb = NULL,
                                   repair_genes = FALSE,
                                   show_row_names = NULL,
                                   label_size = 10,
                                   label_specific_rows = NULL,
                                   label_specific_rows_gp = grid::gpar(fontsize = 5),
                                   show_column_names = TRUE,
                                   cluster_rows = TRUE,
                                   show_row_dend = TRUE,
                                   cluster_columns = TRUE,
                                   kmeans_k = NULL,
                                   annotate_columns = FALSE,
                                   column_anno_palette = "Set2",
                                   color_default = TRUE,
                                   col = NULL,
                                   heatmap_name = NULL,
                                   show_heatmap_legend = TRUE,
                                   return_type = c("heatmap", "clusters", "both"),
                                   ...) {

  stopifnot(inherits(x, "VISTA"))
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE)) {
    cli::cli_abort("Packages {.pkg ComplexHeatmap} and {.pkg circlize} must be installed to draw heatmaps.")
  }
  return_type <- match.arg(return_type)

  comps <- .vista_comparisons(x)
  if (!length(comps)) {
    cli::cli_abort("No differential expression results are stored in this VISTA object.")
  }
  if (is.null(sample_comparisons)) {
    sample_comparisons <- names(comps)
  }
  .validate_fc_inputs(x, sample_comparisons, genes = NULL)

  auto_genes <- is.null(genes)
  genes <- .validate_heatmap_gene_input(genes)
  if (auto_genes) {
    genes <- .auto_foldchange_heatmap_genes(x, sample_comparisons = sample_comparisons, top_n = top_n)
  } else {
    .validate_fc_inputs(x, sample_comparisons, genes)
  }

  fc_mat <- get_foldchange_matrix(
    x,
    sample_comparisons = sample_comparisons,
    genes = genes
  )
  if (is.null(show_row_names)) {
    show_row_names <- auto_genes
  }
  show_row_names <- isTRUE(show_row_names)

  disp_genes <- .resolve_heatmap_row_labels(
    x = x,
    gene_ids = rownames(fc_mat),
    display_id = display_id,
    display_from = display_from,
    display_orgdb = display_orgdb,
    repair_genes = repair_genes,
    prefer_symbol = auto_genes
  )
  rownames(fc_mat) <- disp_genes

  if (repair_genes && !auto_genes) {
    rownames(fc_mat) <- gsub(".*?:", "", rownames(fc_mat))
  }
  if (!is.null(label_specific_rows) && repair_genes) {
    label_specific_rows <- gsub(".*?:", "", label_specific_rows)
  }

  if (!is.null(label_specific_rows) && show_row_names) {
    warning("`show_row_names = TRUE` is incompatible with `label_specific_rows`; setting `show_row_names = FALSE`.")
    show_row_names <- FALSE
  }

  if (color_default) {
    rng <- range(fc_mat, na.rm = TRUE)
    max_abs <- max(abs(rng))
    col_fun <- circlize::colorRamp2(c(-max_abs, 0, max_abs), c("blue", "white", "red"))
  } else {
    col_fun <- col
  }

  if (is.null(heatmap_name)) heatmap_name <- "Log2 Fold-change"

  row_anno <- NULL
  if (!is.null(label_specific_rows)) {
    label_specific_rows <- .resolve_heatmap_row_labels(
      x = x,
      gene_ids = label_specific_rows,
      display_id = display_id,
      display_from = display_from,
      display_orgdb = display_orgdb,
      repair_genes = repair_genes,
      prefer_symbol = auto_genes
    )
    match_rows <- intersect(label_specific_rows, rownames(fc_mat))
    if (length(match_rows) > 0) {
      row_indices <- which(rownames(fc_mat) %in% match_rows)
      row_anno <- ComplexHeatmap::rowAnnotation(
        labels = ComplexHeatmap::anno_mark(
          at = row_indices,
          labels = rownames(fc_mat)[row_indices],
          labels_gp = label_specific_rows_gp
        )
      )
    }
  }

  col_anno <- NULL
  if (annotate_columns) {
    group_col <- .vista_group_col(x)
    col_df <- SummarizedExperiment::colData(x) |>
      as.data.frame() |>
      tibble::rownames_to_column("sample") |>
      dplyr::filter(sample %in% colnames(fc_mat)) |>
      dplyr::select(sample, !!group_col)

    rownames(col_df) <- col_df$sample
    col_df[[group_col]] <- factor(col_df[[group_col]])
    group_levels <- levels(col_df[[group_col]])
    group_colors <- colorspace::qualitative_hcl(length(group_levels), palette = column_anno_palette)
    names(group_colors) <- group_levels

    col_anno <- ComplexHeatmap::HeatmapAnnotation(
      df = as.data.frame(col_df[, group_col, drop = FALSE]),
      col = setNames(list(group_colors), group_col)
    )
  }

  ht <- ComplexHeatmap::Heatmap(as.matrix(fc_mat),
                                name = heatmap_name,
                                col = col_fun,
                                cluster_rows = cluster_rows,
                                show_row_dend = show_row_dend,
                                show_row_names = show_row_names,
                                row_names_gp = grid::gpar(fontsize = label_size),
                                show_column_names = show_column_names,
                                cluster_columns = cluster_columns,
                                show_heatmap_legend = show_heatmap_legend,
                                row_km = if (!is.null(kmeans_k)) kmeans_k else NULL,
                                row_km_repeats = 1,
                                top_annotation = col_anno,
                                ...)

  if (!is.null(row_anno)) ht <- ht + row_anno

  if (!is.null(kmeans_k) && return_type != "heatmap") {
    drawn <- ComplexHeatmap::draw(ht)
    clusters <- ComplexHeatmap::row_order(drawn)
    row_names <- rownames(fc_mat)
    row_clusters <- rep(NA_integer_, length(row_names))
    for (k in seq_along(clusters)) row_clusters[clusters[[k]]] <- k
    cluster_df <- tibble::tibble(gene = row_names, cluster = row_clusters)
    if (return_type == "clusters") return(cluster_df)
    return(list(heatmap = ht, clusters = cluster_df))
  }

  if (return_type == "heatmap") return(ht)
  list(heatmap = ht, clusters = NULL)
}

#' @title Extract a log2 fold-change matrix
#' @description Returns a gene-by-comparison matrix of log2 fold changes stored in a VISTA object.
#'
#' @param x A VISTA object containing differential expression results.
#' @param sample_comparisons Optional character vector of comparison names. Defaults to all available comparisons.
#' @param genes Optional character vector of gene identifiers. When omitted, all genes present in `row_data(x)` are returned.
#'
#' @return A numeric matrix with genes in rows and comparisons in columns.
#' @export
get_foldchange_matrix <- function(x,
                                   sample_comparisons = NULL,
                                   genes = NULL) {

  stopifnot(inherits(x, "VISTA"))

  comps <- comparisons(x)
  if (!length(comps)) {
    cli::cli_abort("No differential expression results are stored in this VISTA object.")
  }

  if (is.null(sample_comparisons)) {
    sample_comparisons <- names(comps)
  }

  .validate_fc_inputs(x, sample_comparisons, genes)

  if (is.null(genes)) {
    genes <- .available_fc_genes(x, sample_comparisons = sample_comparisons)
  }

  norm_genes <- as.character(genes)

  extract_fc <- function(tbl) {
    df <- as.data.frame(tbl, stringsAsFactors = FALSE)

    if (!"gene_id" %in% names(df)) {
      rn <- rownames(df)
      if (!is.null(rn) && length(rn) == nrow(df) && all(nzchar(rn))) {
        df$gene_id <- rn
      } else {
        cli::cli_abort("A DE comparison table lacks a {.field gene_id} column and usable rownames.")
      }
    }

    if (!"log2fc" %in% names(df)) {
      rename_candidates <- c(log2FoldChange = "log2fc", logFC = "log2fc", LFC = "log2fc")
      for (nm in names(rename_candidates)) {
        if (nm %in% names(df)) {
          names(df)[names(df) == nm] <- rename_candidates[[nm]]
          break
        }
      }
    }

    if (!"log2fc" %in% names(df)) {
      cli::cli_abort("A DE comparison table is missing a {.field log2fc} column.")
    }

    stats::setNames(df$log2fc, df$gene_id)
  }

  fc_list <- purrr::map(sample_comparisons, ~ extract_fc(comps[[.x]]))
  names(fc_list) <- sample_comparisons

  fc_mat <- vapply(sample_comparisons, function(comp) {
    vals <- fc_list[[comp]][norm_genes]
    as.numeric(vals)
  }, FUN.VALUE = numeric(length(norm_genes)))

  rownames(fc_mat) <- norm_genes
  colnames(fc_mat) <- sample_comparisons

  fc_mat
}
