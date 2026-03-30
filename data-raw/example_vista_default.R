# Regenerate the precomputed default object used by example_vista().
#
# Run from the package root:
#   /Library/Frameworks/R.framework/Resources/bin/Rscript data-raw/example_vista_default.R

pkg_dir <- normalizePath(file.path(getwd()))
load_all_ok <- requireNamespace("pkgload", quietly = TRUE)
if (!load_all_ok) {
  stop("Package 'pkgload' is required to regenerate the precomputed example object.")
}

pkgload::load_all(pkg_dir, quiet = TRUE)

vista_example_default <- VISTA:::.build_example_vista(
  n_genes = 150,
  n_per_group = 3,
  method = "deseq2"
)

save(
  vista_example_default,
  file = file.path(pkg_dir, "R", "sysdata.rda"),
  compress = "xz",
  version = 2
)
