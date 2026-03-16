% Generated manually: update alongside R/benchmarking.R
\name{validate_vista_deep}
\alias{validate_vista_deep}
\title{Deep validation of VISTA differential-expression fidelity}
\usage{
validate_vista_deep(
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
  return_plots = FALSE,
  error = TRUE
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

\item{error}{When \code{TRUE}, abort if any discrepancy is detected.}
}
\value{
Invisibly returns the full benchmark report from
\code{\link{benchmark_vista_equivalence}}.
}
\description{
Combines structural validation from \code{\link{validate_vista}} with
backend-to-backend equivalence checks against standalone DESeq2, edgeR, and
limma runs.
}
\examples{
\donttest{
data("count_data", package = "VISTA")
data("sample_metadata", package = "VISTA")

target_groups <- c("control", "treatment1")
sample_subset <- sample_metadata[sample_metadata$cond_long %in% target_groups, ]
count_subset <- count_data[1:150, c("gene_id", sample_subset$sample_names)]

validate_vista_deep(
  counts = count_subset,
  sample_info = sample_subset,
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control",
  methods = "deseq2",
  min_counts = 5,
  min_replicates = 1
)
}
}
