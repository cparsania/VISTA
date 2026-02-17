#' Example RNA-seq count matrix shipped with VISTA
#'
#' A count matrix derived from the Bioconductor `airway` dataset and formatted
#' for VISTA examples. The first column stores Ensembl gene IDs; the remaining
#' columns are sample-level counts.
#'
#' @format A tibble with 63,677 rows and 9 columns.
#'   \describe{
#'     \item{gene_id}{Character column storing Ensembl gene identifiers.}
#'     \item{SRR1039508}{Numeric count values for airway sample SRR1039508.}
#'     \item{SRR1039509}{Numeric count values for airway sample SRR1039509.}
#'     \item{SRR1039512}{Numeric count values for airway sample SRR1039512.}
#'     \item{SRR1039513}{Numeric count values for airway sample SRR1039513.}
#'     \item{SRR1039516}{Numeric count values for airway sample SRR1039516.}
#'     \item{SRR1039517}{Numeric count values for airway sample SRR1039517.}
#'     \item{SRR1039520}{Numeric count values for airway sample SRR1039520.}
#'     \item{SRR1039521}{Numeric count values for airway sample SRR1039521.}
#'   }
#' @usage data(count_data)
#' @docType data
#' @name count_data
#' @keywords datasets
#' @source Derived from the `airway` Bioconductor dataset.
"count_data"

#' Sample metadata accompanying the VISTA airway example counts
#'
#' Metadata describing each column in `count_data`, including friendly sample
#' names and the experimental group assignment used throughout the vignettes and
#' tests.
#'
#' @format A tibble with 8 rows and 13 columns:
#' \describe{
#'   \item{sample_names}{Character identifiers matching the column names in
#'   `count_data`.}
#'   \item{groups}{Group labels used in examples (`"control"`, `"treatment1"`).}
#'   \item{cond_long}{Long-form condition labels (`"control"`, `"treatment1"`).}
#'   \item{cond_short}{Short-form condition labels (`"CTRL"`, `"TREAT1"`).}
#'   \item{SampleName}{Original sample label from the airway colData.}
#'   \item{dex}{Original airway treatment labels (`"untrt"`/`"trt"`).}
#'   \item{cell}{Donor/cell-line identifiers from airway.}
#'   \item{albut}{Original airway albuterol indicator.}
#'   \item{Run}{SRA run identifier (also copied into `sample_names`).}
#'   \item{avgLength}{Average transcript length metadata from airway.}
#'   \item{Experiment}{SRA experiment accession from airway.}
#'   \item{Sample}{SRA sample accession from airway.}
#'   \item{BioSample}{NCBI BioSample accession from airway.}
#' }
#' @usage data(sample_metadata)
#' @docType data
#' @name sample_metadata
#' @keywords datasets
#' @seealso [count_data]
#' @source Derived from the `airway` Bioconductor dataset.
"sample_metadata"
