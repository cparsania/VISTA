# The registry is the source of truth for VISTA's deprecation timelines. These
# tests make it impossible to add an alias and forget to wire up its warning --
# which is exactly how six silent aliases accumulated before 1.2.0.

test_that("the deprecation registry is well formed", {
  reg <- VISTA:::.vista_deprecations()

  expect_s3_class(reg, "data.frame")
  expect_gt(nrow(reg), 0L)
  expect_named(
    reg,
    c("fun", "old_arg", "new_arg", "deprecated_in", "defunct_in", "removed_in", "note")
  )

  expect_false(anyNA(reg))
  expect_true(all(nzchar(reg$fun)))
  expect_true(all(nzchar(reg$old_arg)))

  # A function/argument pair may only appear once.
  key <- paste(reg$fun, reg$old_arg, sep = "/")
  expect_false(anyDuplicated(key) > 0L)

  # Timelines must be ordered: deprecated < defunct < removed.
  for (col in c("deprecated_in", "defunct_in", "removed_in")) {
    expect_true(all(nzchar(reg[[col]])), info = col)
  }
  expect_true(all(
    package_version(reg$deprecated_in) < package_version(reg$defunct_in)
  ))
  expect_true(all(
    package_version(reg$defunct_in) < package_version(reg$removed_in)
  ))
})

test_that("the promised timeline is a valid Bioconductor deprecation cycle", {
  reg <- VISTA:::.vista_deprecations()
  # unclass() first: indexing a numeric_version yields another version object,
  # which has no arithmetic.
  y <- function(v) vapply(unclass(package_version(v)), function(p) p[[2]], numeric(1))

  # Bioconductor releases carry an even y; devel is odd. A cycle promised at an
  # odd y would name a release that never ships.
  for (col in c("deprecated_in", "defunct_in", "removed_in")) {
    expect_true(all(y(reg[[col]]) %% 2 == 0), info = col)
  }

  # Bioconductor requires warn -> defunct -> remove across *consecutive*
  # releases, i.e. y steps of exactly 2. A wider gap silently grants users an
  # extra cycle the warnings do not mention; a narrower one is impossible.
  expect_true(all(y(reg$defunct_in) - y(reg$deprecated_in) == 2))
  expect_true(all(y(reg$removed_in) - y(reg$defunct_in) == 2))
})

test_that("the devel version is still on track to become the promised release", {
  # Every shipped warning names `deprecated_in` as the release the rename lands
  # in. That promise only holds while this devel series is the one that becomes
  # that release: devel 1.1.z -> release 1.2.0. If devel is bumped past it
  # without shipping, the warnings already in users' hands become wrong, and the
  # fix is to update the registry rather than to let it drift.
  reg <- VISTA:::.vista_deprecations()
  target <- unique(package_version(reg$deprecated_in))
  expect_length(target, 1L)

  devel <- package_version(as.character(utils::packageVersion("VISTA")))

  # `<=` rather than "exactly one cycle before": on the release branch the
  # version *is* the promised one, which fulfils the promise rather than
  # breaking it. What must never happen is the version moving past the promised
  # release while the registry still names it -- at that point every warning
  # already in users' hands cites a release that has been and gone.
  expect_lte(
    devel, target,
    label = paste0(
      "devel ", as.character(devel), " has passed the promised ",
      as.character(target), "; ship it or re-date the registry"
    )
  )
})

test_that("every registered function is exported and still has the old formal", {
  reg <- VISTA:::.vista_deprecations()
  exported <- getNamespaceExports("VISTA")

  for (i in seq_len(nrow(reg))) {
    fun <- reg$fun[[i]]
    old <- reg$old_arg[[i]]
    new <- reg$new_arg[[i]]

    expect_true(fun %in% exported, info = sprintf("%s is not exported", fun))

    fmls <- names(formals(getExportedValue("VISTA", fun)))
    if (is.null(fmls)) {
      # S4 generic: inspect the method instead.
      fmls <- names(formals(getMethod(fun, "VISTA")))
    }

    expect_true(
      old %in% fmls || "..." %in% fmls,
      info = sprintf("%s() no longer accepts the deprecated `%s`", fun, old)
    )

    if (nzchar(new)) {
      expect_true(
        new %in% fmls || "..." %in% fmls,
        info = sprintf("%s() does not have the replacement `%s`", fun, new)
      )
    }
  }
})

test_that(".vista_deprecate_arg warns with the documented class and returns the value", {
  expect_warning(
    out <- VISTA:::.vista_deprecate_arg(
      old = "show_corr_values", value = TRUE, fun = "get_corr_heatmap"
    ),
    class = "vista_deprecated_arg"
  )
  expect_true(out)

  # It also carries base R's deprecation class, so suppressWarnings and
  # Bioconductor tooling treat it the way they treat .Deprecated().
  expect_warning(
    VISTA:::.vista_deprecate_arg(old = "col_up", value = "red", fun = "get_volcano_plot"),
    class = "deprecatedWarning"
  )

  # transform is applied to the legacy value.
  suppressWarnings(
    expect_identical(
      VISTA:::.vista_deprecate_arg(
        old = "label", value = TRUE, fun = "get_deg_count_pieplot",
        transform = function(v) if (isTRUE(v)) "both" else "none"
      ),
      "both"
    )
  )
})

test_that("the deprecation message names the replacement and the defunct release", {
  msg <- tryCatch(
    VISTA:::.vista_deprecate_arg(old = "show_corr_values", value = TRUE, fun = "get_corr_heatmap"),
    vista_deprecated_arg = function(w) conditionMessage(w)
  )
  expect_match(msg, "show_corr_values", fixed = TRUE)
  expect_match(msg, "label", fixed = TRUE)
  expect_match(msg, "1.4.0", fixed = TRUE)

  # Arguments with no replacement say so rather than pointing at "".
  msg2 <- tryCatch(
    VISTA:::.vista_deprecate_arg(old = "sample.seed", value = 1, fun = "get_pca_plot"),
    vista_deprecated_arg = function(w) conditionMessage(w)
  )
  expect_match(msg2, "no longer has any effect", fixed = TRUE)
})

test_that(".vista_defunct_arg aborts with the documented class", {
  expect_error(
    VISTA:::.vista_defunct_arg(old = "col_up", fun = "get_volcano_plot"),
    class = "vista_defunct_arg"
  )
  expect_error(
    VISTA:::.vista_defunct_arg(old = "col_up", fun = "get_volcano_plot"),
    class = "defunctError"
  )
})

test_that(".vista_check_dots accepts known names and rejects unknown ones", {
  expect_true(VISTA:::.vista_check_dots(list(), fun = "get_volcano_plot"))
  expect_true(
    VISTA:::.vista_check_dots(
      list(pointSize = 2), fun = "get_volcano_plot", allowed = "pointSize"
    )
  )

  expect_error(
    VISTA:::.vista_check_dots(list(definitely_not_real = 1), fun = "get_volcano_plot"),
    "unknown argument"
  )

  # Unnamed arguments in ... are always a mistake for these functions.
  expect_error(
    VISTA:::.vista_check_dots(stats::setNames(list(1), ""), fun = "get_volcano_plot"),
    "unnamed argument"
  )

  # Managed arguments are blocked with a distinct message.
  expect_error(
    VISTA:::.vista_check_dots(
      list(genes = "x"), fun = "get_pathway_heatmap", blocked = c("genes", "sample_group")
    ),
    "managed by"
  )
})

test_that(".vista_check_dots suggests a near-miss argument name", {
  msg <- tryCatch(
    VISTA:::.vista_check_dots(list(point_siz = 2), fun = "get_volcano_plot"),
    error = function(e) conditionMessage(e)
  )
  expect_match(msg, "did you mean", ignore.case = TRUE)
  expect_match(msg, "point_size", fixed = TRUE)
})

test_that("shared defaults are available and internally consistent", {
  d <- VISTA:::.vista_defaults()
  expect_type(d, "list")
  expect_true(all(nzchar(names(d))))
  expect_false(anyDuplicated(names(d)) > 0L)

  expect_identical(VISTA:::.vista_default("group_palette"), "Dark 2")
  expect_error(VISTA:::.vista_default("no_such_default"), "Unknown VISTA default")

  # 1.2.0 must not move any default; these pin the current values.
  expect_identical(VISTA:::.vista_default("point_size_embedding"), 10)
  expect_identical(VISTA:::.vista_default("label_size"), 3)
  expect_identical(VISTA:::.vista_default("max_genes_embedding"), 20)
})

test_that(".vista_escape_yaml produces parseable double-quoted scalars", {
  skip_if_not_installed("yaml")

  values <- c(
    "plain title",
    'has "quotes"',
    "C:\\project\\dir",
    "trailing\\",
    'mixed "q" and \\ slash'
  )

  for (v in values) {
    y <- paste0("title: \"", VISTA:::.vista_escape_yaml(v), "\"")
    parsed <- yaml::yaml.load(y)
    expect_identical(parsed$title, v, info = encodeString(v))
  }
})

test_that(".vista_sanitize_name strips path-significant characters", {
  expect_identical(
    VISTA:::.vista_sanitize_name("treatment1_VS_control"),
    "treatment1_VS_control"
  )
  # A backslash was previously whitelisted by "[^A-Za-z0-9_\\-]" and reached
  # asset filenames.
  expect_identical(VISTA:::.vista_sanitize_name("A\\B_VS_C"), "A_B_VS_C")
  expect_identical(VISTA:::.vista_sanitize_name("a/b:c*d"), "a_b_c_d")
  expect_identical(VISTA:::.vista_sanitize_name("keeps-hyphen"), "keeps-hyphen")
  expect_identical(VISTA:::.vista_sanitize_name("__trimmed__"), "trimmed")
  expect_identical(VISTA:::.vista_sanitize_name("///", fallback = "comparison"), "comparison")
  expect_false(any(grepl("[\\\\/]", VISTA:::.vista_sanitize_name(c("a\\b", "c/d")))))
})

test_that("the previously-silent aliases now warn (3a)", {
  v <- make_small_vista()
  comp <- names(comparisons(v))[[1]]

  expect_warning(
    get_corr_heatmap(v, show_corr_values = FALSE),
    class = "vista_deprecated_arg"
  )
  expect_warning(
    get_corr_heatmap(v, col_corr_values = "navy"),
    class = "vista_deprecated_arg"
  )
  expect_warning(
    get_volcano_plot(v, sample_comparison = comp, col_up = "red"),
    class = "vista_deprecated_arg"
  )
  expect_warning(
    get_pca_plot(v, sample.seed = 42),
    class = "vista_deprecated_arg"
  )
  expect_warning(
    get_expression_barplot(v, genes = rownames(v)[1:2], facet_scale = "fixed"),
    class = "vista_deprecated_arg"
  )
  expect_warning(
    get_expression_lollipop(v, genes = rownames(v)[1:2], facet_scale = "fixed"),
    class = "vista_deprecated_arg"
  )
  expect_warning(
    get_expression_violinplot(v, genes = rownames(v)[1:2], value_transform = "none"),
    class = "vista_deprecated_arg"
  )
  expect_warning(
    get_expression_lineplot(v, genes = rownames(v)[1:2], value_transform = "none"),
    class = "vista_deprecated_arg"
  )
})

test_that("deprecated aliases still produce the same result as before", {
  v <- make_small_vista()

  # legacy wins, and the outcome matches passing the new argument directly
  legacy <- suppressWarnings(get_corr_heatmap(v, label = TRUE, show_corr_values = FALSE))
  modern <- get_corr_heatmap(v, label = FALSE)
  expect_identical(length(legacy$layers), length(modern$layers))

  genes <- rownames(v)[1:2]
  a <- suppressWarnings(get_expression_barplot(v, genes = genes, facet_scale = "fixed"))
  b <- get_expression_barplot(v, genes = genes, facet_scales = "fixed")
  expect_equal(ggplot2::ggplot_build(a)$data, ggplot2::ggplot_build(b)$data)
})

test_that("gene caps are documented arguments rather than hard limits (3e)", {
  v <- make_small_vista()
  genes <- rownames(v)

  expect_error(
    get_expression_lollipop(v, genes = genes[seq_len(16)]),
    "At most 15"
  )
  expect_no_error(
    get_expression_lollipop(v, genes = genes[seq_len(16)], max_genes = 20)
  )

  expect_error(
    get_expression_barplot(v, genes = genes[seq_len(26)]),
    "At most 25"
  )
  expect_no_error(
    get_expression_barplot(v, genes = genes[seq_len(26)], max_genes = 30)
  )
})

test_that("semantic-collision renames warn and preserve behaviour (3c)", {
  v <- make_small_vista()

  # `label` was a character enum here while it is logical everywhere else.
  expect_warning(p_old <- get_deg_count_pieplot(v, label = "percent"), class = "vista_deprecated_arg")
  p_new <- get_deg_count_pieplot(v, label_type = "percent")
  expect_equal(ggplot2::ggplot_build(p_old)$data, ggplot2::ggplot_build(p_new)$data)

  # label = TRUE previously did nothing useful; it now maps to "both".
  expect_warning(p_true <- get_deg_count_donutplot(v, label = TRUE), class = "vista_deprecated_arg")
  expect_equal(
    ggplot2::ggplot_build(p_true)$data,
    ggplot2::ggplot_build(get_deg_count_donutplot(v, label_type = "both"))$data
  )

  # cluster_by named a column in the heatmaps but an ordering strategy here.
  expect_warning(c_old <- get_corr_heatmap(v, cluster_by = "input"), class = "vista_deprecated_arg")
  c_new <- get_corr_heatmap(v, order_by = "input")
  expect_equal(ggplot2::ggplot_build(c_old)$data, ggplot2::ggplot_build(c_new)$data)
})

test_that("stat_comparisons replaces the colliding comparisons argument (3c)", {
  v <- make_small_vista()
  genes <- rownames(v)[seq_len(2)]
  groups <- unique(as.character(sample_info(v)$cond_long))
  pairs <- list(groups)

  expect_warning(
    a <- get_expression_boxplot(v, genes = genes, stats_group = TRUE, comparisons = pairs),
    class = "vista_deprecated_arg"
  )
  b <- get_expression_boxplot(v, genes = genes, stats_group = TRUE, stat_comparisons = pairs)
  expect_equal(ggplot2::ggplot_build(a)$data, ggplot2::ggplot_build(b)$data)

  # The accessor of the same name is untouched.
  expect_type(comparisons(v), "list")
})

test_that("concept names are unified with warning aliases (3d)", {
  v <- make_small_vista()

  # top_n_genes -> top_n on the embeddings
  expect_warning(a <- get_pca_plot(v, top_n_genes = 50), class = "vista_deprecated_arg")
  b <- get_pca_plot(v, top_n = 50)
  expect_equal(ggplot2::ggplot_build(a)$data, ggplot2::ggplot_build(b)$data)

  expect_warning(get_mds_plot(v, top_n_genes = 50), class = "vista_deprecated_arg")

  # line_size -> linewidth on the lollipops (matches ggplot2 >= 3.4)
  genes <- rownames(v)[seq_len(2)]
  expect_warning(
    l_old <- get_expression_lollipop(v, genes = genes, line_size = 2),
    class = "vista_deprecated_arg"
  )
  l_new <- get_expression_lollipop(v, genes = genes, linewidth = 2)
  expect_equal(ggplot2::ggplot_build(l_old)$data, ggplot2::ggplot_build(l_new)$data)
})
