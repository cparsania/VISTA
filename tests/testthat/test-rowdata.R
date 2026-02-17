test_that("set_rowdata appends annotations by key column", {
  vista <- make_small_vista()
  ann <- data.frame(
    gene_id = rownames(vista),
    SYMBOL = paste0("SYM", seq_len(nrow(vista))),
    stringsAsFactors = FALSE
  )
  vista2 <- set_rowdata(vista, annotations = ann, key_col = "gene_id", overwrite = TRUE)
  expect_true("SYMBOL" %in% colnames(rowData(vista2)))
  expect_equal(rowData(vista2)$SYMBOL[1], "SYM1")
})
