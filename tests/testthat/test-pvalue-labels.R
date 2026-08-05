# p-value annotation must survive ggpubr's internals changing.
#
# ggpubr >= 1.0.0 builds `label = "p.format"` into an expression calling its
# own create_p_label(), which is not visible from the environment of a mapping
# supplied by VISTA -- so the plot built fine and then died at draw time with
# `could not find function "create_p_label"`. Nothing caught it because the
# failure is in ggplot_build(), not in the call that adds the layer.
#
# These tests pin two things: the mapping never delegates label construction to
# ggpubr, and every stats-capable plot actually draws.

test_that("the p-value mapping names only public stat columns", {
  for (lab in c("p.format", "p.signif", "p.adj", "p")) {
    m <- VISTA:::.vista_pvalue_label_aes(lab)
    expect_true("label" %in% names(m), info = lab)
    txt <- paste(vapply(m, rlang::quo_text, character(1)), collapse = " ")
    # the whole point: no reference to a ggpubr internal
    expect_false(grepl("create_p_label", txt, fixed = TRUE), info = lab)
    expect_match(txt, "after_stat")
  }

  # An unrecognised label falls back to ggpubr rather than dropping silently.
  expect_null(VISTA:::.vista_pvalue_label_aes("something_else"))
  expect_null(VISTA:::.vista_pvalue_label_aes(NULL))
  expect_null(VISTA:::.vista_pvalue_label_aes(c("p", "p.adj")))
})

test_that("known labels stop passing `label` on to ggpubr; unknown ones keep it", {
  base <- ggplot2::aes(group = .data$grp)

  known <- VISTA:::.vista_compare_means_args(base, "p.format")
  expect_null(known$label)
  expect_true(all(c("group", "label") %in% names(known$mapping)))

  unknown <- VISTA:::.vista_compare_means_args(base, "..custom..")
  expect_identical(unknown$label, "..custom..")
  expect_false("label" %in% names(unknown$mapping))
  expect_true("group" %in% names(unknown$mapping))
})

test_that("stats-annotated plots build, and report the same text for each label", {
  skip_if_not_installed("ggpubr")
  v <- make_small_vista()
  genes <- rownames(v)[seq_len(6)]

  drawn_label <- function(p) {
    b <- ggplot2::ggplot_build(p)          # where the old failure surfaced
    for (ld in b$data) if ("label" %in% names(ld)) {
      return(unique(as.character(ld$label)))
    }
    character()
  }

  fmt <- get_expression_boxplot(
    v, genes = genes, stats_group = TRUE, p.label = "p.format"
  )
  sig <- get_expression_boxplot(
    v, genes = genes, stats_group = TRUE, p.label = "p.signif"
  )

  # p.format keeps ggpubr's "p = <value>" rendering
  expect_true(all(grepl("^p = ", drawn_label(fmt))))
  # p.signif is a significance code, never a "p = " string
  expect_false(any(grepl("^p = ", drawn_label(sig))))

  skip_if_not_installed("ggrain")
  rain <- get_expression_raincloud(
    v, genes = genes, value_transform = "log2", summarise = TRUE,
    facet_by = "none", stats_group = TRUE,
    stats_method = "wilcox.test", p.label = "p.format"
  )
  expect_no_error(ggplot2::ggplot_build(rain))
  expect_true(all(grepl("^p = ", drawn_label(rain))))
})
