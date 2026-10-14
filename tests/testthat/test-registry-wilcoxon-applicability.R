test_that("Wilcoxon registry exposes conditional option and parent applicability", {
  ids <- c("t_wilcoxon_signed_rank_one_sample", "t_wilcoxon_signed_rank_matched_pairs",
           "t_wilcoxon_mann_whitney")
  for (id in ids) {
    opts <- .ig_get_metadata(id)$options
    expect_identical(opts$distribution_plot$applies_when$condition$option, "method")
    expect_identical(.ig_flatten(opts$distribution_plot$applies_when$condition$`in`, "character"), "lehmann")
    expect_identical(opts$effect_size_parameter$applies_when$condition$option, "lehmann_input_method")
    expect_identical(.ig_flatten(opts$effect_size_parameter$applies_when$condition$`in`, "character"), "predefined_parent")
    parent <- opts$parent_distribution$value_applies_when
    expect_identical(names(parent), opts$parent_distribution$codes)
    expect_identical(parent$min_are$option, "method")
    expect_identical(.ig_flatten(parent$min_are$`in`, "character"), "are")
    expect_null(parent$normal)
    expect_null(parent$laplace)
    expect_null(parent$logistic)
  }
})
