test_that("read_vista_counts standardizes matrix-like inputs", {
  data("count_data", package = "VISTA", envir = environment())

  counts_in <- count_data[seq_len(12), ]
  cnt <- read_vista_counts(
    counts_in,
    format = "matrix",
    gene_id_column = "gene_id",
    verbose = FALSE
  )

  expect_named(cnt, c("counts", "row_data", "column_geneid", "sample_names", "sample_name_map", "input_format", "report"))
  expect_identical(cnt$column_geneid, "gene_id")
  expect_identical(cnt$input_format, "matrix")
  expect_identical(cnt$counts$gene_id, counts_in$gene_id)
  expect_identical(cnt$sample_names, colnames(counts_in)[-1])
  expect_true(all(vapply(cnt$counts[-1], is.numeric, logical(1))))
})

test_that("read_vista_counts repairs file-path sample columns conservatively", {
  data("count_data", package = "VISTA", envir = environment())

  sample_cols <- colnames(count_data)[colnames(count_data) != "gene_id"][seq_len(4)]
  counts_in <- count_data[seq_len(6), c("gene_id", sample_cols), drop = FALSE]
  names(counts_in)[2:5] <- c(
    "/proj/run/03_alignment/HET_1_U/HET_1_U_star_alignAligned.sortedByCoord.out.bam",
    "/proj/run/03_alignment/HET_1_ovary/HET_1_ovary_star_alignAligned.sortedByCoord.out.bam",
    "/proj/run/03_alignment/WT_1_U/WT_1_U_star_alignAligned.sortedByCoord.out.bam",
    "/proj/run/03_alignment/WT_1_ovary/WT_1_ovary_star_alignAligned.sortedByCoord.out.bam"
  )

  cnt <- read_vista_counts(
    counts_in,
    format = "matrix",
    gene_id_column = "gene_id",
    verbose = FALSE
  )

  expect_identical(
    cnt$sample_names,
    c("HET_1_U", "HET_1_ovary", "WT_1_U", "WT_1_ovary")
  )
  expect_true(all(cnt$sample_name_map$original != cnt$sample_name_map$repaired))
})

test_that("read_vista_counts repairs generic quantification file names via parent directory", {
  data("count_data", package = "VISTA", envir = environment())

  sample_cols <- colnames(count_data)[colnames(count_data) != "gene_id"][seq_len(2)]
  counts_in <- count_data[seq_len(4), c("gene_id", sample_cols), drop = FALSE]
  names(counts_in)[2:3] <- c(
    "/proj/salmon/sample_A/quant.sf",
    "/proj/salmon/sample_B/quant.sf"
  )

  cnt <- read_vista_counts(
    counts_in,
    format = "matrix",
    gene_id_column = "gene_id",
    verbose = FALSE
  )

  expect_identical(cnt$sample_names, c("sample_A", "sample_B"))
})

test_that("read_vista_counts strips lane and read suffixes when unique", {
  data("count_data", package = "VISTA", envir = environment())

  sample_cols <- colnames(count_data)[colnames(count_data) != "gene_id"][seq_len(2)]
  counts_in <- count_data[seq_len(4), c("gene_id", sample_cols), drop = FALSE]
  names(counts_in)[2:3] <- c(
    "sampleA_S1_L001_R1_001.fastq.gz",
    "sampleB_S2_L002_R1_001.fastq.gz"
  )

  cnt <- read_vista_counts(
    counts_in,
    format = "matrix",
    gene_id_column = "gene_id",
    verbose = FALSE
  )

  expect_identical(cnt$sample_names, c("sampleA", "sampleB"))
})

test_that("read_vista_counts parses featureCounts-style tables", {
  data("count_data", package = "VISTA", envir = environment())

  fc <- data.frame(
    Geneid = count_data$gene_id[seq_len(8)],
    Chr = rep("chr1", 8),
    Start = seq(100, by = 100, length.out = 8),
    End = seq(150, by = 100, length.out = 8),
    Strand = "+",
    Length = 50,
    sampleA = unname(count_data[[2]][seq_len(8)]),
    sampleB = unname(count_data[[3]][seq_len(8)]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  tmp <- tempfile(fileext = ".txt")
  writeLines(c(
    "# Program:featureCounts",
    paste(colnames(fc), collapse = "\t"),
    apply(fc, 1, paste, collapse = "\t")
  ), con = tmp)

  cnt <- read_vista_counts(tmp, format = "auto", verbose = FALSE)

  expect_identical(cnt$input_format, "featurecounts")
  expect_identical(colnames(cnt$counts), c("gene_id", "sampleA", "sampleB"))
  expect_true(all(c("Chr", "Start", "End", "Strand", "Length") %in% colnames(cnt$row_data)))
})

test_that("read_vista_counts parses STAR per-sample files", {
  star_lines <- function(scale = 1) {
    c(
      "N_unmapped\t1\t1\t1",
      "N_multimapping\t2\t2\t2",
      "N_noFeature\t3\t3\t3",
      "N_ambiguous\t4\t4\t4",
      paste("geneA", 10 * scale, 11 * scale, 12 * scale, sep = "\t"),
      paste("geneB", 20 * scale, 21 * scale, 22 * scale, sep = "\t")
    )
  }

  f1 <- tempfile(pattern = "sampleA_", fileext = ".ReadsPerGene.out.tab")
  f2 <- tempfile(pattern = "sampleB_", fileext = ".ReadsPerGene.out.tab")
  writeLines(star_lines(1), f1)
  writeLines(star_lines(2), f2)

  cnt <- read_vista_counts(
    c(f1, f2),
    format = "star",
    count_column = "unstranded",
    verbose = FALSE
  )

  expect_identical(cnt$input_format, "star")
  expect_identical(cnt$counts$gene_id, c("geneA", "geneB"))
  expect_equal(as.numeric(cnt$counts[1, -1]), c(10, 20))
})

test_that("read_vista_counts parses HTSeq per-sample files", {
  htseq_lines <- function(scale = 1) {
    c(
      paste("geneA", 5 * scale, sep = "\t"),
      paste("geneB", 7 * scale, sep = "\t"),
      paste("__no_feature", 9 * scale, sep = "\t")
    )
  }

  f1 <- tempfile(pattern = "sampleA_", fileext = ".txt")
  f2 <- tempfile(pattern = "sampleB_", fileext = ".txt")
  writeLines(htseq_lines(1), f1)
  writeLines(htseq_lines(2), f2)

  cnt <- read_vista_counts(c(f1, f2), format = "htseq", verbose = FALSE)

  expect_identical(cnt$input_format, "htseq")
  expect_identical(cnt$counts$gene_id, c("geneA", "geneB"))
  expect_equal(as.numeric(cnt$counts[2, -1]), c(7, 14))
})

test_that("read_vista_counts parses tximport-like inputs with tx2gene", {
  tx_counts <- matrix(
    c(10, 20, 30, 40),
    nrow = 2,
    dimnames = list(c("tx1", "tx2"), c("sampleA", "sampleB"))
  )
  tx_list <- list(
    counts = tx_counts,
    abundance = tx_counts / 10,
    length = tx_counts * 2
  )
  tx2gene <- data.frame(
    transcript_id = c("tx1", "tx2"),
    gene_id = c("geneA", "geneA"),
    stringsAsFactors = FALSE
  )

  cnt <- read_vista_counts(
    tx_list,
    format = "tximport",
    tx2gene = tx2gene,
    verbose = FALSE
  )

  expect_identical(cnt$input_format, "tximport")
  expect_identical(cnt$counts$gene_id, "geneA")
  expect_equal(as.numeric(cnt$counts[1, -1]), c(30, 70))
})

test_that("read_vista_counts parses RSEM gene result files", {
  rsem_df <- data.frame(
    gene_id = c("geneA", "geneB"),
    transcript_id = c("tx1", "tx2"),
    expected_count = c(12, 24),
    TPM = c(1.2, 2.4),
    FPKM = c(0.6, 1.2),
    stringsAsFactors = FALSE
  )
  f1 <- tempfile(pattern = "sampleA_", fileext = ".genes.results")
  write.table(rsem_df, file = f1, sep = "\t", quote = FALSE, row.names = FALSE)

  cnt <- read_vista_counts(f1, format = "rsem", verbose = FALSE)

  expect_identical(cnt$input_format, "rsem")
  expect_equal(ncol(cnt$counts), 2)
  expect_identical(colnames(cnt$counts)[1], "gene_id")
  expect_true("transcript_id" %in% colnames(cnt$row_data))
})

test_that("read_vista_metadata infers sample_names and match_vista_inputs reorders samples", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  counts_in <- count_data[seq_len(10), c("gene_id", sample_metadata$sample_names[seq_len(4)]), drop = FALSE]
  sample_in <- sample_metadata[seq_len(4), c("cond_long", "sample_names"), drop = FALSE]
  sample_in <- sample_in[c(4, 2, 1, 3), , drop = FALSE]
  rownames(sample_in) <- NULL
  names(sample_in)[names(sample_in) == "sample_names"] <- "sample"

  cnt <- read_vista_counts(counts_in, format = "matrix", gene_id_column = "gene_id", verbose = FALSE)
  si <- read_vista_metadata(sample_in, verbose = FALSE)
  matched <- match_vista_inputs(cnt, si, verbose = FALSE)

  expect_identical(matched$sample_info$sample_names, colnames(matched$counts)[-1])
  expect_identical(rownames(matched$sample_info), matched$sample_info$sample_names)
})

test_that("match_vista_inputs drops unmatched samples when requested", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  counts_in <- count_data[seq_len(10), c("gene_id", sample_metadata$sample_names[seq_len(4)]), drop = FALSE]
  sample_in <- sample_metadata[seq_len(5), c("sample_names", "cond_long"), drop = FALSE]

  matched <- match_vista_inputs(
    counts = counts_in,
    sample_info = sample_in,
    column_geneid = "gene_id",
    drop_unmatched = TRUE,
    verbose = FALSE
  )

  expect_equal(length(matched$sample_names), 4)
  expect_true(all(matched$sample_names %in% sample_in$sample_names))
})

test_that("derive_vista_metadata parses split-based sample names", {
  meta <- derive_vista_metadata(
    counts = NULL,
    sample_names = c("HET_1_U", "HET_2_ovary", "WT_1_U"),
    parser = "split",
    split = "_",
    fields = c("genotype", "replicate", "tissue"),
    verbose = FALSE
  )

  expect_identical(colnames(meta), c("sample_names", "genotype", "replicate", "tissue"))
  expect_identical(meta$genotype, c("HET", "HET", "WT"))
  expect_identical(meta$tissue, c("U", "ovary", "U"))
})

test_that("derive_vista_metadata parses regex-based sample names", {
  meta <- derive_vista_metadata(
    counts = NULL,
    sample_names = c("SRR1039508", "SRR1039509"),
    parser = "regex",
    pattern = "SRR(\\d+)",
    fields = "run_id",
    verbose = FALSE
  )

  expect_identical(meta$run_id, c("1039508", "1039509"))
})

test_that("derive_vista_metadata auto-detects consistent split sample names", {
  meta <- derive_vista_metadata(
    counts = NULL,
    sample_names = c("WT-1-ovary", "WT-2-uterus"),
    parser = "auto",
    fields = c("genotype", "replicate", "tissue"),
    verbose = FALSE
  )

  expect_identical(meta$genotype, c("WT", "WT"))
  expect_identical(meta$replicate, c("1", "2"))
})

test_that("derive_vista_metadata can derive from repaired count sample names", {
  data("count_data", package = "VISTA", envir = environment())

  sample_cols <- colnames(count_data)[colnames(count_data) != "gene_id"][seq_len(2)]
  counts_in <- count_data[seq_len(4), c("gene_id", sample_cols), drop = FALSE]
  names(counts_in)[2:3] <- c(
    "/proj/run/03_alignment/HET_1_U/HET_1_U_star_alignAligned.sortedByCoord.out.bam",
    "/proj/run/03_alignment/WT_2_ovary/WT_2_ovary_star_alignAligned.sortedByCoord.out.bam"
  )

  meta <- derive_vista_metadata(
    counts = counts_in,
    column_geneid = "gene_id",
    parser = "split",
    fields = c("genotype", "replicate", "tissue"),
    verbose = FALSE
  )

  expect_identical(meta$sample_names, c("HET_1_U", "WT_2_ovary"))
  expect_identical(meta$replicate, c("1", "2"))
  expect_true(is.data.frame(attr(meta, "sample_name_map")))
})

test_that("derive_vista_metadata template adds placeholder columns", {
  meta <- derive_vista_metadata(
    counts = NULL,
    sample_names = c("sampleA", "sampleB"),
    parser = "none",
    return_type = "template",
    verbose = FALSE
  )

  expect_true(all(c("group", "batch") %in% colnames(meta)))
  expect_true(all(is.na(meta$group)))
  expect_true(all(is.na(meta$batch)))
})

# --- B6: multi-file importers must bind by gene identifier, not by position ---

test_that("STAR files with shuffled gene order are matched by identifier", {
  star_file <- function(genes, values) {
    lines <- c(
      "N_unmapped\t1\t1\t1",
      "N_multimapping\t2\t2\t2",
      "N_noFeature\t3\t3\t3",
      "N_ambiguous\t4\t4\t4",
      paste(genes, values, values, values, sep = "\t")
    )
    f <- tempfile(fileext = ".ReadsPerGene.out.tab")
    writeLines(lines, f)
    f
  }

  genes <- c("geneA", "geneB", "geneC")
  f1 <- star_file(genes, c(10, 20, 30))
  # Same genes, different row order, distinct values.
  f2 <- star_file(rev(genes), c(300, 200, 100))
  on.exit(unlink(c(f1, f2)), add = TRUE)

  cnt <- read_vista_counts(
    c(f1, f2), format = "star", sample_names = c("s1", "s2"), verbose = FALSE
  )

  expect_identical(cnt$counts$gene_id, genes)
  expect_equal(cnt$counts$s1, c(10, 20, 30))
  # Positional binding would have produced 300, 200, 100 here.
  expect_equal(cnt$counts$s2, c(100, 200, 300))
})

test_that("HTSeq files with shuffled gene order are matched by identifier", {
  htseq_file <- function(genes, values) {
    f <- tempfile(fileext = ".txt")
    writeLines(c(paste(genes, values, sep = "\t"), "__no_feature\t99"), f)
    f
  }

  genes <- c("geneA", "geneB", "geneC")
  f1 <- htseq_file(genes, c(1, 2, 3))
  f2 <- htseq_file(rev(genes), c(30, 20, 10))
  on.exit(unlink(c(f1, f2)), add = TRUE)

  cnt <- read_vista_counts(
    c(f1, f2), format = "htseq", sample_names = c("s1", "s2"), verbose = FALSE
  )

  expect_identical(cnt$counts$gene_id, genes)
  expect_equal(cnt$counts$s1, c(1, 2, 3))
  expect_equal(cnt$counts$s2, c(10, 20, 30))
})

test_that("RSEM files with shuffled gene order are matched by identifier", {
  rsem_file <- function(genes, values) {
    f <- tempfile(fileext = ".genes.results")
    utils::write.table(
      data.frame(
        gene_id = genes, transcript_id = paste0("tx_", genes),
        expected_count = values, stringsAsFactors = FALSE
      ),
      file = f, sep = "\t", quote = FALSE, row.names = FALSE
    )
    f
  }

  genes <- c("geneA", "geneB", "geneC")
  f1 <- rsem_file(genes, c(5, 6, 7))
  f2 <- rsem_file(rev(genes), c(70, 60, 50))
  on.exit(unlink(c(f1, f2)), add = TRUE)

  cnt <- read_vista_counts(
    c(f1, f2), format = "rsem", sample_names = c("s1", "s2"), verbose = FALSE
  )

  expect_identical(cnt$counts$gene_id, genes)
  expect_equal(cnt$counts$s1, c(5, 6, 7))
  expect_equal(cnt$counts$s2, c(50, 60, 70))
  # row_data is taken from the first file and must stay aligned.
  expect_identical(as.character(cnt$row_data$gene_id), genes)
})

test_that("a file that does not cover the reference gene set is an error", {
  htseq_file <- function(genes, values) {
    f <- tempfile(fileext = ".txt")
    writeLines(paste(genes, values, sep = "\t"), f)
    f
  }

  f1 <- htseq_file(c("geneA", "geneB", "geneC"), c(1, 2, 3))
  f2 <- htseq_file(c("geneA", "geneB", "geneZ"), c(1, 2, 3))
  on.exit(unlink(c(f1, f2)), add = TRUE)

  # Equal row counts, so the old positional bind accepted this silently.
  expect_error(
    read_vista_counts(
      c(f1, f2), format = "htseq", sample_names = c("s1", "s2"), verbose = FALSE
    ),
    "do not cover the reference gene set"
  )

  msg <- tryCatch(
    read_vista_counts(
      c(f1, f2), format = "htseq", sample_names = c("s1", "s2"), verbose = FALSE
    ),
    error = function(e) conditionMessage(e)
  )
  expect_match(msg, basename(f2), fixed = TRUE)
  expect_match(msg, "geneC", fixed = TRUE)
})

test_that("extra identifiers in a later file are dropped with a message", {
  htseq_file <- function(genes, values) {
    f <- tempfile(fileext = ".txt")
    writeLines(paste(genes, values, sep = "\t"), f)
    f
  }

  f1 <- htseq_file(c("geneA", "geneB"), c(1, 2))
  f2 <- htseq_file(c("geneB", "geneA", "geneEXTRA"), c(20, 10, 999))
  on.exit(unlink(c(f1, f2)), add = TRUE)

  expect_message(
    cnt <- read_vista_counts(
      c(f1, f2), format = "htseq", sample_names = c("s1", "s2"), verbose = FALSE
    ),
    "Dropped 1 identifier"
  )
  expect_identical(cnt$counts$gene_id, c("geneA", "geneB"))
  expect_equal(cnt$counts$s2, c(10, 20))
})

test_that("identically ordered files still bind positionally", {
  # Fast path: identical id vectors, including repeated identifiers, must work
  # without requiring uniqueness.
  tabs <- list(
    data.frame(id = c("g1", "g1", "g2"), n = c(1, 2, 3), stringsAsFactors = FALSE),
    data.frame(id = c("g1", "g1", "g2"), n = c(4, 5, 6), stringsAsFactors = FALSE)
  )
  out <- VISTA:::.bind_count_files(
    tabs, id_col = "id", value_col = "n", sample_names = c("a", "b")
  )
  expect_identical(out$gene_id, c("g1", "g1", "g2"))
  expect_equal(out$a, c(1, 2, 3))
  expect_equal(out$b, c(4, 5, 6))
})

test_that("duplicated identifiers with differing order are rejected", {
  tabs <- list(
    data.frame(id = c("g1", "g1", "g2"), n = c(1, 2, 3), stringsAsFactors = FALSE),
    data.frame(id = c("g2", "g1", "g1"), n = c(3, 1, 2), stringsAsFactors = FALSE)
  )
  expect_error(
    VISTA:::.bind_count_files(
      tabs, id_col = "id", value_col = "n", sample_names = c("a", "b")
    ),
    "duplicated gene identifiers"
  )
})
