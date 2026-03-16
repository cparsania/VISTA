% Generated manually: update alongside R/benchmarking.R
\name{benchmark_vista_equivalence}
\alias{benchmark_vista_equivalence}
\title{Benchmark VISTA against standalone differential-expression backends}
\usage{
benchmark_vista_equivalence(
  counts,
  sample_info,
  column_geneid,
  group_column,
  group_numerator,
  group_denominator,
  methods = c("deseq2", "edger", "limma"),
  min_counts = 10,
  min_replicates = 1,
  log2fc_cutoff = 1,
  pval_cutoff = 0.05,
  p_value_type = "padj",
  covariates = NULL,
  design_formula = NULL,
  tolerance = 1e-08,
  return_plots = FALSE
)
}
\arguments{
\item{counts}{Raw counts with a gene-id column and sample columns.}

\item{sample_info}{Sample metadata.}

\item{column_geneid}{Column name in \code{counts} containing gene identifiers.}

\item{group_column}{Grouping column in \code{sample_info}.}

\item{group_numerator}{Numerator group(s) for pairwise comparisons.}

\item{group_denominator}{Denominator group(s) for pairwise comparisons.}

\item{methods}{Subset of \code{c("deseq2", "edger", "limma")} to benchmark.}

\item{min_counts}{Minimum total counts per gene to retain.}

\item{min_replicates}{Minimum samples per group meeting filtering criteria.}

\item{log2fc_cutoff}{Absolute log2 fold-change threshold for DEG calling.}

\item{pval_cutoff}{P-value or adjusted p-value threshold.}

\item{p_value_type}{One of \code{"padj"} or \code{"pvalue"}.}

\item{covariates}{Optional covariates included in the design.}

\item{design_formula}{Optional design formula overriding automatic design construction.}

\item{tolerance}{Floating-point tolerance for numerical comparisons.}

\item{return_plots}{When \code{TRUE}, include paired VISTA/reference plots for
MA, volcano, DEG count, and PCA checks.}
}
\value{
A list with \code{valid}, \code{comparison_summary}, \code{visual_summary}, and
\code{methods}. The \code{methods} element contains the VISTA object, direct
backend outputs, structural validation, self-consistency checks, and optional
plot objects.
}
\description{
Runs VISTA and direct DESeq2, edgeR, and limma pipelines with matched
preprocessing, contrasts, and DEG thresholds, then compares DE tables,
normalized counts, DEG calls, and critical visualization payloads.
}
\examples{
\donttest{
data("count_data", package = "VISTA")
data("sample_metadata", package = "VISTA")

target_groups <- c("control", "treatment1")
sample_subset <- sample_metadata[sample_metadata$cond_long %in% target_groups, ]
count_subset <- count_data[1:150, c("gene_id", sample_subset$sample_names)]

bm <- benchmark_vista_equivalence(
  counts = count_subset,
  sample_info = sample_subset,
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control",
  methods = c("deseq2", "edger"),
  min_counts = 5,
  min_replicates = 1
)

bm$comparison_summary
}
}
