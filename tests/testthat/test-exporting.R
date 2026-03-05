test_that("save_vista_plot saves ggplot outputs including PDF", {
  vista <- make_small_vista()
  pca_plot <- get_pca_plot(vista)

  tmp_dir <- tempfile("vista_plot_export_")
  dir.create(tmp_dir, recursive = TRUE)

  png_file <- file.path(tmp_dir, "pca_plot.png")
  pdf_file <- file.path(tmp_dir, "pca_plot.pdf")

  expect_no_error(
    save_vista_plot(
      pca_plot,
      file = png_file,
      width = 6,
      height = 4,
      units = "in",
      dpi = 150
    )
  )
  expect_true(file.exists(png_file))

  expect_no_error(
    save_vista_plot(
      pca_plot,
      file = pdf_file,
      width = 15,
      height = 10,
      units = "cm"
    )
  )
  expect_true(file.exists(pdf_file))
})

test_that("save_vista_data exports CSV and RDS payloads", {
  vista <- make_small_vista()

  tmp_dir <- tempfile("vista_data_export_")
  dir.create(tmp_dir, recursive = TRUE)

  comparison_csv <- file.path(tmp_dir, "comparison.csv")
  comparisons_rds <- file.path(tmp_dir, "comparisons.rds")

  expect_no_error(
    save_vista_data(
      vista,
      what = "comparison",
      file = comparison_csv,
      sample_comparison = names(comparisons(vista))[1],
      format = "csv"
    )
  )
  expect_true(file.exists(comparison_csv))
  expect_true("gene_id" %in% colnames(utils::read.csv(comparison_csv, check.names = FALSE)))

  expect_no_error(
    save_vista_data(
      vista,
      what = c("comparison", "norm_counts"),
      file = comparisons_rds,
      format = "rds"
    )
  )
  expect_true(file.exists(comparisons_rds))
  payload <- readRDS(comparisons_rds)
  expect_true(is.list(payload))
  expect_true(all(c("comparison", "norm_counts") %in% names(payload)))
})

test_that("save_vista_data can export XLSX when writexl is available", {
  skip_if_not_installed("writexl")

  vista <- make_small_vista()
  tmp_dir <- tempfile("vista_xlsx_export_")
  dir.create(tmp_dir, recursive = TRUE)
  out_xlsx <- file.path(tmp_dir, "vista_tables.xlsx")

  expect_no_error(
    save_vista_data(
      vista,
      what = c("comparison", "norm_counts", "sample_info"),
      file = out_xlsx,
      format = "xlsx"
    )
  )
  expect_true(file.exists(out_xlsx))
})

test_that("export_vista_assets creates plot/data bundle and manifest", {
  vista <- make_small_vista()

  tmp_dir <- tempfile("vista_assets_bundle_")
  dir.create(tmp_dir, recursive = TRUE)

  result <- export_vista_assets(
    vista,
    out_dir = tmp_dir,
    include_plots = c("pca", "mds"),
    include_data = c("comparison", "norm_counts"),
    plot_format = "pdf",
    write_excel = FALSE
  )

  expect_true(is.list(result))
  expect_true(file.exists(file.path(tmp_dir, "manifest.csv")))
  expect_true(file.exists(file.path(tmp_dir, "plots", "pca.pdf")))
  expect_true(file.exists(file.path(tmp_dir, "plots", "mds.pdf")))
  expect_true(file.exists(file.path(tmp_dir, "tables", "comparison.csv")))
  expect_true(
    file.exists(file.path(tmp_dir, "tables", "norm_counts.csv")),
    info = paste(capture.output(print(result$manifest)), collapse = "\n")
  )
  expect_true(is.data.frame(result$manifest))
  expect_true(nrow(result$manifest) >= 4)
})

test_that("save_vista_plot validates dimensions", {
  vista <- make_small_vista()
  p <- get_pca_plot(vista)

  expect_error(
    save_vista_plot(p, file = tempfile(fileext = ".png"), width = -1),
    "width"
  )
  expect_error(
    save_vista_plot(p, file = tempfile(fileext = ".png"), height = 0),
    "height"
  )
})

test_that("save_vista_data validates format-specific rules", {
  vista <- make_small_vista()

  expect_error(
    save_vista_data(
      vista,
      what = c("comparison", "norm_counts"),
      file = tempfile(fileext = ".csv"),
      format = "csv"
    ),
    "exactly one value"
  )

  expect_error(
    save_vista_data(vista, what = "not_supported", file = tempfile(fileext = ".rds"), format = "rds"),
    "Unsupported"
  )
})

test_that("export_vista_assets validates include keys and overwrite mode", {
  vista <- make_small_vista()

  expect_error(
    export_vista_assets(vista, out_dir = tempfile("vista_assets_bad_"), include_plots = "not_a_plot"),
    "Unsupported"
  )

  tmp_dir <- tempfile("vista_assets_existing_")
  dir.create(tmp_dir, recursive = TRUE)
  writeLines("x", file.path(tmp_dir, "already_here.txt"))

  expect_error(
    export_vista_assets(vista, out_dir = tmp_dir, overwrite = FALSE),
    "already contains files"
  )
})
