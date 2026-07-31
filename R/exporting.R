#' Save a VISTA plot object to disk
#'
#' Saves plot objects returned by VISTA plotting functions to file. Supports
#' both `ggplot`-like objects (saved via `ggplot2::ggsave()`) and
#' `ComplexHeatmap` objects (`Heatmap` / `HeatmapList`) saved via graphics
#' devices.
#'
#' @param plot A plot object. Typically `ggplot`, `patchwork`, `Heatmap`, or
#'   `HeatmapList`.
#' @param file Output file path.
#' @param width Plot width.
#' @param height Plot height.
#' @param units Units for `width` and `height`. One of `"in"`, `"cm"`, `"mm"`,
#'   or `"px"`.
#' @param dpi Resolution for raster outputs.
#' @param device Optional graphics device (e.g. `"png"`, `"pdf"`). If `NULL`,
#'   inferred from `file` extension (defaults to `"png"` when missing).
#' @param ... Additional arguments passed to `ggplot2::ggsave()` for ggplot-like
#'   objects.
#'
#' @return Invisibly, the normalized output file path.
#' @examples
#' v <- example_vista()
#' p <- get_pca_plot(v)
#' out_file <- tempfile(fileext = ".pdf")
#' save_vista_plot(p, file = out_file, width = 7, height = 5, units = "in")
#' @export
save_vista_plot <- function(plot,
                            file,
                            width = 8,
                            height = 6,
                            units = "in",
                            dpi = 300,
                            device = NULL,
                            ...) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  if (is.null(plot)) {
    cli::cli_abort("{.arg plot} cannot be NULL.")
  }
  if (!is.character(file) || length(file) != 1L || !nzchar(file)) {
    cli::cli_abort("{.arg file} must be a non-empty file path.")
  }
  if (!is.numeric(width) || length(width) != 1L || !is.finite(width) || width <= 0) {
    cli::cli_abort("{.arg width} must be a positive number.")
  }
  if (!is.numeric(height) || length(height) != 1L || !is.finite(height) || height <= 0) {
    cli::cli_abort("{.arg height} must be a positive number.")
  }
  if (!is.numeric(dpi) || length(dpi) != 1L || !is.finite(dpi) || dpi <= 0) {
    cli::cli_abort("{.arg dpi} must be a positive number.")
  }
  units <- match.arg(units, c("in", "cm", "mm", "px"))

  out_dir <- dirname(file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ext <- tolower(tools::file_ext(file))
  dev <- tolower(device %||% ext)
  if (!nzchar(dev)) dev <- "png"
  if (identical(dev, "jpg")) dev <- "jpeg"

  to_inches <- function(value, unit, dpi_val) {
    switch(
      unit,
      "in" = value,
      "cm" = value / 2.54,
      "mm" = value / 25.4,
      "px" = value / dpi_val,
      cli::cli_abort("{.arg units} must be one of {.val in}, {.val cm}, {.val mm}, {.val px}.")
    )
  }

  save_heatmap <- function(hm, dev_name) {
    if (dev_name %in% c("png", "jpeg", "tiff", "bmp")) {
      dev_fun <- switch(
        dev_name,
        png = grDevices::png,
        jpeg = grDevices::jpeg,
        tiff = grDevices::tiff,
        bmp = grDevices::bmp
      )
      dev_fun(file, width = width, height = height, units = units, res = dpi)
    } else if (identical(dev_name, "pdf")) {
      grDevices::pdf(
        file,
        width = to_inches(width, units, dpi),
        height = to_inches(height, units, dpi)
      )
    } else {
      cli::cli_abort(
        c(
          "Unsupported device {.val {dev_name}} for ComplexHeatmap objects.",
          "i" = "Use one of: {.val png}, {.val jpeg}, {.val tiff}, {.val bmp}, {.val pdf}."
        )
      )
    }
    on.exit(grDevices::dev.off(), add = TRUE)
    ComplexHeatmap::draw(hm)
  }

  if (inherits(plot, c("Heatmap", "HeatmapList"))) {
    save_heatmap(plot, dev)
  } else {
    ggplot2::ggsave(
      filename = file,
      plot = plot,
      width = width,
      height = height,
      units = units,
      dpi = dpi,
      device = device,
      ...
    )
  }

  invisible(normalizePath(file, mustWork = FALSE))
}

#' Save VISTA tabular outputs to disk
#'
#' Exports selected data components from a `VISTA` object to CSV/TSV/RDS/XLSX.
#'
#' @param x A `VISTA` object.
#' @param what Character vector specifying which object(s) to export. Supported
#'   values are `"comparison"`, `"comparisons"`, `"norm_counts"`,
#'   `"sample_info"`, `"row_data"`, `"deg_summary"`, and `"cutoffs"`.
#' @param file Output file path.
#' @param sample_comparison Optional comparison name used when `what` includes
#'   `"comparison"`. Defaults to the first comparison in `comparisons(x)`.
#' @param format Output format. One of `"csv"`, `"tsv"`, `"rds"`, `"xlsx"`.
#'   If `NULL`, inferred from `file` extension.
#' @param include_rownames Logical; include meaningful row identifiers (e.g.,
#'   gene IDs or sample names) as explicit columns where applicable.
#'
#' @return Invisibly, the normalized output file path.
#' @examples
#' v <- example_vista()
#' save_vista_data(v, what = "comparison", file = tempfile(fileext = ".csv"), format = "csv")
#' @export
save_vista_data <- function(x,
                            what = c(
                              "comparison",
                              "comparisons",
                              "norm_counts",
                              "sample_info",
                              "row_data",
                              "deg_summary",
                              "cutoffs"
                            ),
                            file,
                            sample_comparison = NULL,
                            format = NULL,
                            include_rownames = TRUE) {
  stopifnot(inherits(x, "VISTA"))

  allowed_what <- c("comparison", "comparisons", "norm_counts", "sample_info", "row_data", "deg_summary", "cutoffs")
  if (length(what) == 0) {
    cli::cli_abort("{.arg what} cannot be empty.")
  }
  what <- unique(as.character(what))
  bad_what <- setdiff(what, allowed_what)
  if (length(bad_what) > 0) {
    cli::cli_abort("Unsupported {.arg what} value(s): {.val {bad_what}}.")
  }

  if (!is.character(file) || length(file) != 1L || !nzchar(file)) {
    cli::cli_abort("{.arg file} must be a non-empty file path.")
  }
  out_dir <- dirname(file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  fmt <- format
  if (is.null(fmt) || !nzchar(fmt)) {
    fmt <- tolower(tools::file_ext(file))
  }
  if (!nzchar(fmt)) fmt <- "csv"
  fmt <- tolower(fmt)
  if (identical(fmt, "xls")) fmt <- "xlsx"
  fmt <- match.arg(fmt, c("csv", "tsv", "rds", "xlsx"))

  comps <- comparisons(x)
  if (!length(comps)) {
    cli::cli_abort("No comparisons available in this VISTA object.")
  }
  comp_name <- sample_comparison
  if (is.null(comp_name)) comp_name <- names(comps)[1]
  if (!comp_name %in% names(comps)) {
    cli::cli_abort(
      "{.arg sample_comparison} '{comp_name}' not found. Available: {.val {names(comps)}}"
    )
  }

  add_rownames <- function(df, colname) {
    if (!include_rownames) return(df)
    rn <- rownames(df)
    if (is.null(rn)) return(df)
    if (colname %in% colnames(df)) return(df)
    out <- data.frame(tmp = rn, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(out) <- colname
    cbind(out, as.data.frame(df, stringsAsFactors = FALSE))
  }

  to_df <- function(obj, key) {
    if (identical(key, "cutoffs")) {
      if (is.null(obj) || !length(obj)) {
        return(data.frame(name = character(), value = character(), stringsAsFactors = FALSE))
      }
      return(data.frame(
        name = names(obj),
        value = vapply(obj, function(v) paste(v, collapse = ", "), character(1)),
        stringsAsFactors = FALSE
      ))
    }

    if (is.matrix(obj)) obj <- as.data.frame(obj, stringsAsFactors = FALSE, check.names = FALSE)
    if (inherits(obj, "DataFrame")) obj <- as.data.frame(obj, stringsAsFactors = FALSE)
    if (!is.data.frame(obj)) {
      obj <- data.frame(value = as.character(obj), stringsAsFactors = FALSE, check.names = FALSE)
    }

    if (identical(key, "comparison")) {
      if (!"gene_id" %in% colnames(obj)) obj <- add_rownames(obj, "gene_id")
    } else if (identical(key, "norm_counts")) {
      obj <- add_rownames(obj, "gene_id")
    } else if (identical(key, "sample_info")) {
      obj <- add_rownames(obj, "sample_names")
    } else if (identical(key, "row_data")) {
      obj <- add_rownames(obj, "gene_id")
    } else if (identical(key, "deg_summary")) {
      obj <- add_rownames(obj, "sample_comparison")
    }
    rownames(obj) <- NULL
    obj
  }

  get_object <- function(key) {
    switch(
      key,
      comparison = as.data.frame(comps[[comp_name]], stringsAsFactors = FALSE),
      comparisons = lapply(comps, as.data.frame, stringsAsFactors = FALSE),
      norm_counts = norm_counts(x),
      sample_info = as.data.frame(sample_info(x), stringsAsFactors = FALSE),
      row_data = as.data.frame(row_data(x), stringsAsFactors = FALSE),
      deg_summary = as.data.frame(deg_summary(x), stringsAsFactors = FALSE),
      cutoffs = cutoffs(x)
    )
  }

  if (identical(fmt, "rds")) {
    payload <- setNames(lapply(what, get_object), what)
    if (length(payload) == 1L) payload <- payload[[1]]
    saveRDS(payload, file)
    return(invisible(normalizePath(file, mustWork = FALSE)))
  }

  if (identical(fmt, "csv") || identical(fmt, "tsv")) {
    if (length(what) != 1L) {
      cli::cli_abort(
        c(
          "{.arg what} must contain exactly one value for CSV/TSV export.",
          "i" = "Use {.arg format = \"xlsx\"} or {.arg format = \"rds\"} for multiple outputs."
        )
      )
    }
    key <- what[[1]]
    obj <- get_object(key)
    if (identical(key, "comparisons")) {
      tables <- lapply(obj, to_df, key = "comparison")
      if (length(tables)) {
        tables <- lapply(names(tables), function(nm) {
          tbl <- tables[[nm]]
          tbl$sample_comparison <- nm
          tbl
        })
        df <- dplyr::bind_rows(tables)
      } else {
        df <- data.frame(stringsAsFactors = FALSE)
      }
    } else {
      df <- to_df(obj, key)
    }
    if (identical(fmt, "csv")) {
      utils::write.csv(df, file, row.names = FALSE)
    } else {
      utils::write.table(df, file, sep = "\t", quote = TRUE, row.names = FALSE, col.names = TRUE)
    }
    return(invisible(normalizePath(file, mustWork = FALSE)))
  }

  if (!requireNamespace("writexl", quietly = TRUE)) {
    cli::cli_abort(
      c(
        "XLSX export requires the {.pkg writexl} package.",
        "i" = "Install it with {.code install.packages(\"writexl\")}, or use {.arg format = \"csv\"}/{.arg format = \"rds\"}."
      )
    )
  }

  sanitize_sheet <- function(x) {
    out <- gsub("[\\\\/?*\\[\\]:]", "_", x)
    out <- trimws(out)
    if (!nzchar(out)) out <- "sheet"
    substr(out, 1, 31)
  }

  sheet_state <- new.env(parent = emptyenv())
  sheet_state$names <- character()
  sheet_state$tables <- list()
  append_sheet <- function(name, tbl) {
    base <- sanitize_sheet(name)
    nm <- base
    idx <- 1L
    while (nm %in% sheet_state$names) {
      suffix <- paste0("_", idx)
      nm <- paste0(substr(base, 1, max(1L, 31L - nchar(suffix))), suffix)
      idx <- idx + 1L
    }
    sheet_state$names <- c(sheet_state$names, nm)
    sheet_state$tables[[nm]] <- tbl
  }

  for (key in what) {
    obj <- get_object(key)
    if (identical(key, "comparisons")) {
      for (nm in names(obj)) {
        append_sheet(paste0("comparison_", nm), to_df(obj[[nm]], "comparison"))
      }
    } else {
      append_sheet(key, to_df(obj, key))
    }
  }

  writexl::write_xlsx(sheet_state$tables, path = file)
  invisible(normalizePath(file, mustWork = FALSE))
}

#' Export a complete VISTA asset bundle
#'
#' Generates a standardized folder with selected VISTA plots, tabular outputs,
#' and a manifest describing all saved files.
#'
#' @param x A `VISTA` object.
#' @param out_dir Output directory for exported assets.
#' @param sample_comparison Optional comparison to use for comparison-specific
#'   outputs. Defaults to the first available comparison.
#' @param display_id Optional gene identifier column used in labeling for
#'   volcano/MA/heatmap plots.
#' @param include_plots Character vector of plot keys to export. Supported:
#'   `"pca"`, `"mds"`, `"corr_heatmap"`, `"deg_bar"`, `"deg_pie"`,
#'   `"deg_donut"`, `"volcano"`, `"ma"`, `"expression_heatmap"`.
#' @param include_data Character vector of data keys passed to
#'   [save_vista_data()].
#' @param plot_format Plot format (e.g. `"png"` or `"pdf"`).
#' @param width Base plot width.
#' @param height Base plot height.
#' @param heatmap_height Height used specifically for expression heatmap export.
#' @param units Plot dimension units.
#' @param dpi Raster resolution for plots.
#' @param top_n_labels Number of top genes to annotate in MA plots.
#' @param heatmap_n_genes Number of top genes used in exported expression
#'   heatmaps.
#' @param write_excel Logical; if `TRUE`, also writes a combined XLSX workbook
#'   for all requested `include_data` tables (requires \pkg{writexl}).
#' @param overwrite Logical; if `FALSE`, aborts when `out_dir` already contains
#'   files.
#'
#' @return Invisibly, a list with `out_dir`, `sample_comparison`, `manifest`,
#'   `plot_files`, and `data_files`.
#' @examples
#' v <- example_vista()
#' out_dir <- file.path(tempdir(), "vista_assets_example")
#' res <- export_vista_assets(
#'   v,
#'   out_dir = out_dir,
#'   include_plots = "pca",
#'   include_data = "comparison"
#' )
#' names(res)
#' @export
export_vista_assets <- function(x,
                                out_dir = "vista_assets",
                                sample_comparison = NULL,
                                display_id = NULL,
                                include_plots = c("pca", "mds", "corr_heatmap", "deg_bar", "volcano", "ma", "expression_heatmap"),
                                include_data = c("comparison", "norm_counts", "sample_info", "row_data", "deg_summary", "cutoffs"),
                                plot_format = "png",
                                width = 8,
                                height = 6,
                                heatmap_height = 10,
                                units = "in",
                                dpi = 300,
                                top_n_labels = 50,
                                heatmap_n_genes = 60,
                                write_excel = FALSE,
                                overwrite = TRUE) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  stopifnot(inherits(x, "VISTA"))

  if (!is.character(out_dir) || length(out_dir) != 1L || !nzchar(out_dir)) {
    cli::cli_abort("{.arg out_dir} must be a non-empty directory path.")
  }
  if (dir.exists(out_dir) && !isTRUE(overwrite)) {
    has_files <- length(list.files(out_dir, all.files = TRUE, no.. = TRUE)) > 0
    if (has_files) cli::cli_abort("{.arg out_dir} already contains files; set {.arg overwrite = TRUE} to proceed.")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  plot_format <- tolower(plot_format)
  if (identical(plot_format, "jpg")) plot_format <- "jpeg"

  allowed_plots <- c("pca", "mds", "corr_heatmap", "deg_bar", "deg_pie", "deg_donut", "volcano", "ma", "expression_heatmap")
  include_plots <- unique(as.character(include_plots))
  bad_plots <- setdiff(include_plots, allowed_plots)
  if (length(bad_plots) > 0) {
    cli::cli_abort("Unsupported {.arg include_plots} value(s): {.val {bad_plots}}.")
  }

  allowed_data <- c("comparison", "comparisons", "norm_counts", "sample_info", "row_data", "deg_summary", "cutoffs")
  include_data <- unique(as.character(include_data))
  bad_data <- setdiff(include_data, allowed_data)
  if (length(bad_data) > 0) {
    cli::cli_abort("Unsupported {.arg include_data} value(s): {.val {bad_data}}.")
  }

  comps <- comparisons(x)
  if (!length(comps)) cli::cli_abort("No comparisons available in this VISTA object.")
  comp_name <- sample_comparison %||% names(comps)[1]
  if (!comp_name %in% names(comps)) {
    cli::cli_abort(
      "{.arg sample_comparison} '{comp_name}' not found. Available: {.val {names(comps)}}"
    )
  }

  sanitize_name <- function(txt) .vista_sanitize_name(txt, fallback = "comparison")
  comp_tag <- sanitize_name(comp_name)

  plots_dir <- file.path(out_dir, "plots")
  tables_dir <- file.path(out_dir, "tables")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

  state <- new.env(parent = emptyenv())
  state$manifest <- data.frame(
    type = character(),
    key = character(),
    path = character(),
    status = character(),
    message = character(),
    stringsAsFactors = FALSE
  )
  state$plot_files <- character()
  state$data_files <- character()

  add_manifest <- function(type, key, path = NA_character_, status = "ok", message = NA_character_) {
    state$manifest <- rbind(
      state$manifest,
      data.frame(
        type = type,
        key = key,
        path = path,
        status = status,
        message = message,
        stringsAsFactors = FALSE
      )
    )
  }

  save_plot_item <- function(key, object_expr, filename, w = width, h = height) {
    result <- tryCatch(
      eval.parent(substitute(object_expr)),
      error = function(e) e
    )
    if (inherits(result, "error")) {
      add_manifest("plot", key, status = "failed", message = result$message)
      return(invisible(NULL))
    }
    out_file <- file.path(plots_dir, paste0(filename, ".", plot_format))
    path <- tryCatch(
      save_vista_plot(
        plot = result,
        file = out_file,
        width = w,
        height = h,
        units = units,
        dpi = dpi,
        device = plot_format
      ),
      error = function(e) e
    )
    if (inherits(path, "error")) {
      add_manifest("plot", key, status = "failed", message = path$message)
      return(invisible(NULL))
    }
    state$plot_files <- c(state$plot_files, path)
    add_manifest("plot", key, path = path, status = "ok")
    invisible(path)
  }

  if ("pca" %in% include_plots) {
    save_plot_item("pca", get_pca_plot(x, label = TRUE), "pca")
  }
  if ("mds" %in% include_plots) {
    save_plot_item("mds", get_mds_plot(x, label = TRUE), "mds")
  }
  if ("corr_heatmap" %in% include_plots) {
    save_plot_item("corr_heatmap", get_corr_heatmap(x), "corr_heatmap")
  }
  if ("deg_bar" %in% include_plots) {
    save_plot_item("deg_bar", get_deg_count_barplot(x, facet_by = "comparison"), "deg_count_barplot")
  }
  if ("deg_pie" %in% include_plots) {
    save_plot_item("deg_pie", get_deg_count_pieplot(x, facet_by = "comparison"), "deg_count_pieplot")
  }
  if ("deg_donut" %in% include_plots) {
    save_plot_item("deg_donut", get_deg_count_donutplot(x, facet_by = "comparison"), "deg_count_donutplot")
  }
  if ("volcano" %in% include_plots) {
    save_plot_item(
      "volcano",
      get_volcano_plot(x, sample_comparison = comp_name, display_id = display_id),
      paste0("volcano_", comp_tag)
    )
  }
  if ("ma" %in% include_plots) {
    save_plot_item(
      "ma",
      get_ma_plot(
        x,
        sample_comparison = comp_name,
        display_id = display_id,
        label_n = top_n_labels
      ),
      paste0("ma_", comp_tag)
    )
  }

  if ("expression_heatmap" %in% include_plots) {
    de_tbl <- as.data.frame(comps[[comp_name]], stringsAsFactors = FALSE, check.names = FALSE)
    if (!"gene_id" %in% colnames(de_tbl)) de_tbl$gene_id <- rownames(de_tbl)

    p_candidates <- c("padj", "p.adj", "FDR", "pvalue", "PValue")
    p_col <- p_candidates[p_candidates %in% colnames(de_tbl)]
    p_col <- if (length(p_col)) p_col[[1]] else NA_character_

    if ("regulation" %in% colnames(de_tbl)) {
      de_tbl <- de_tbl[de_tbl$regulation != "Other", , drop = FALSE]
    }
    if (!is.na(p_col)) {
      pvals <- suppressWarnings(as.numeric(de_tbl[[p_col]]))
      pvals[is.na(pvals)] <- Inf
      de_tbl <- de_tbl[order(pvals), , drop = FALSE]
    }

    hm_genes <- unique(de_tbl$gene_id)
    hm_genes <- hm_genes[!is.na(hm_genes) & nzchar(hm_genes)]
    hm_genes <- head(hm_genes, heatmap_n_genes)

    si <- as.data.frame(sample_info(x), stringsAsFactors = FALSE)
    group_col <- (S4Vectors::metadata(x)$group$column %||% colnames(si)[1])
    group_values <- unique(as.character(si[[group_col]]))
    group_values <- group_values[!is.na(group_values) & nzchar(group_values)]

    if (length(hm_genes) >= 2L && length(group_values) >= 1L) {
      save_plot_item(
        "expression_heatmap",
        get_expression_heatmap(
          x = x,
          sample_group = group_values,
          genes = hm_genes,
          display_id = display_id,
          value_transform = "zscore",
          annotate_columns = TRUE,
          summarise_replicates = FALSE,
          show_row_names = FALSE
        ),
        paste0("expression_heatmap_", comp_tag),
        w = width,
        h = heatmap_height
      )
    } else {
      add_manifest(
        "plot",
        "expression_heatmap",
        status = "skipped",
        message = "Insufficient genes or groups to export expression heatmap."
      )
    }
  }

  for (key in include_data) {
    out_file <- file.path(tables_dir, paste0(key, ".csv"))
    res <- tryCatch(
      save_vista_data(
        x = x,
        what = key,
        file = out_file,
        sample_comparison = comp_name,
        format = "csv",
        include_rownames = TRUE
      ),
      error = function(e) e
    )

    # Defensive fallback: ensure norm_counts CSV is still produced even if a
    # backend method/path issue prevents save_vista_data from writing the file.
    if (identical(key, "norm_counts") && (inherits(res, "error") || !file.exists(out_file))) {
      res <- tryCatch({
        nc <- as.data.frame(norm_counts(x), stringsAsFactors = FALSE, check.names = FALSE)
        if (!"gene_id" %in% colnames(nc)) {
          nc <- data.frame(
            gene_id = rownames(nc),
            nc,
            stringsAsFactors = FALSE,
            check.names = FALSE,
            row.names = NULL
          )
        }
        rownames(nc) <- NULL
        utils::write.csv(nc, out_file, row.names = FALSE)
        normalizePath(out_file, mustWork = FALSE)
      }, error = function(e) e)
    }

    if (inherits(res, "error") || !file.exists(out_file)) {
      msg <- if (inherits(res, "error")) res$message else "Output file was not created."
      add_manifest("data", key, status = "failed", message = msg)
    } else {
      state$data_files <- c(state$data_files, res)
      add_manifest("data", key, path = res, status = "ok")
    }
  }

  if (isTRUE(write_excel) && length(include_data) > 0) {
    xlsx_file <- file.path(tables_dir, "vista_data.xlsx")
    xlsx_res <- tryCatch(
      save_vista_data(
        x = x,
        what = include_data,
        file = xlsx_file,
        sample_comparison = comp_name,
        format = "xlsx",
        include_rownames = TRUE
      ),
      error = function(e) e
    )
    if (inherits(xlsx_res, "error")) {
      add_manifest("data", "xlsx_bundle", status = "failed", message = xlsx_res$message)
    } else {
      state$data_files <- c(state$data_files, xlsx_res)
      add_manifest("data", "xlsx_bundle", path = xlsx_res, status = "ok")
    }
  } else if (isTRUE(write_excel) && length(include_data) == 0) {
    add_manifest("data", "xlsx_bundle", status = "skipped", message = "No data keys supplied in include_data.")
  }

  manifest_file <- file.path(out_dir, "manifest.csv")
  utils::write.csv(state$manifest, manifest_file, row.names = FALSE)

  result <- list(
    out_dir = normalizePath(out_dir, mustWork = FALSE),
    sample_comparison = comp_name,
    manifest = state$manifest,
    plot_files = unique(state$plot_files),
    data_files = unique(state$data_files)
  )
  invisible(result)
}
