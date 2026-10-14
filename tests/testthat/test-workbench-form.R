test_that("every loaded method and mode has a definition-driven editor", {
  for (method in list_methods()$id) {
    def <- .ig_get_method(method)
    expect_true(.ig_is_string(def$label$en), info = method)
    expect_true(.ig_is_string(def$label$de), info = method)
    for (mode in def$modes) {
      f <- .ig_workbench_form(method, mode)
      expect_false(any(names(f$parameters) %in% f$target), info = paste(method, mode))
      expect_true(all(names(f$options) %in% names(def$options)), info = method)
      expect_true(all(vapply(f$parameters, function(p) p$domain$length == 1 &&
                              p$domain$type %in% c("real", "integer"), logical(1))), info = method)
    }
  }
})

test_that("numeric editors reproduce canonical R specifications and reject inactive inputs", {
  f <- .ig_workbench_form("t_means_one_sample", "post_hoc", list(tails = "two"))
  s <- .ig_workbench_form_spec(f, list(d = "0,5", alpha = "0,05", n_total = "30"), "de")
  direct <- analysis_spec("t_means_one_sample", "post_hoc",
    inputs = list(d = .5, alpha = .05, n_total = 30), options = list(tails = "two"))
  expect_identical(s, direct)
  expect_error(.ig_workbench_form_spec(f, list(power = "0.8")), "inactive")
  expect_error(.ig_workbench_form_spec(f, list(d = 0.5)), "scalar text")
  expect_error(.ig_workbench_form_spec(f, list(d = "", alpha = ".05", n_total = "30")), "missing input")
  expect_error(.ig_workbench_form("t_means_one_sample", "invented"), "does not offer")
  expect_error(.ig_workbench_form("t_means_one_sample", options = list(tails = "invented")), "invalid option")
})

test_that("option changes disclose retired controls and forbidden value resets", {
  f <- .ig_workbench_form("t_wilcoxon_mann_whitney", "post_hoc",
    list(method = "are", distribution_plot = "normalized_z"))
  expect_true("distribution_plot" %in% f$dropped_options)
  expect_false("distribution_plot" %in% names(f$options))
  f <- .ig_workbench_form("t_wilcoxon_mann_whitney", "post_hoc",
    list(method = "lehmann", parent_distribution = "min_are"))
  expect_false("min_are" %in% f$option_choices$parent_distribution)
  expect_identical(f$reset_options$parent_distribution$from, "min_are")
  expect_identical(f$options$parent_distribution, "normal")
})

test_that("missing statistical choices remain explicit and block construction", {
  f <- .ig_workbench_form("boschloo_proportions_two_independent", "post_hoc")
  expect_false(f$ready)
  expect_null(f$options$variant)
  expect_setequal(f$option_choices$variant, c("minlike", "central", "less", "greater"))
  expect_error(.ig_workbench_form_spec(f, list()), "required model settings")
  expect_true(.ig_workbench_form("boschloo_proportions_two_independent", "post_hoc", list(variant = "central"))$ready)
})

test_that("default unconditional sensitivity uses its documented probability target", {
  f <- .ig_workbench_form("exact_proportions_unconditional", "sensitivity")
  expect_true(f$ready)
  expect_identical(f$target, "p1")
  expect_identical(f$options$effect_size_index, "odds_ratio")
  values <- list(p2="0.4", alpha="0.05", power="0.8", n1="30", n2="30")
  expect_identical(.ig_workbench_form_spec(f, values),
    analysis_spec("exact_proportions_unconditional", "sensitivity",
      inputs=list(p2=.4, alpha=.05, power=.8, n1=30, n2=30),options=f$options))
  expect_true(.ig_workbench_form("exact_proportions_unconditional", "sensitivity",
    list(effect_size_index = "proportion_p1"))$ready)
})
