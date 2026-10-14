test_that("Boschloo lattices preserve variants and engine conditional probabilities", {
  for (variant in c("minlike", "central", "greater", "less")) {
    s <- analysis_spec("boschloo_proportions_two_independent", "post_hoc",
      inputs = list(n1 = 4, n2 = 5, p1 = .8, p2 = .2, alpha = .05),
      options = list(variant = variant))
    r <- calculate_power(s)
    pd <- power_plot_data(r, "lattice")
    expect_identical(pd$status, "ok")
    expect_equal(nrow(pd$data), 60)
    expect_identical(pd$fixed_inputs$variant, variant)
    expect_equal(pd$fixed_inputs$n1, 4)
    expect_equal(pd$fixed_inputs$n2, 5)
    expect_equal(as.numeric(tapply(pd$data$probability, pd$data$hypothesis, sum)), c(1, 1))
    expect_equal(pd$region_probabilities$engine, pd$region_probabilities$lattice, tolerance = 1e-12)
    expect_identical(pd$result_hash, r$result_hash)
    expect_equal(power_plot_data(ig_from_json(ig_to_json(r)), "lattice")$data, pd$data)
    for (language in c("en", "de")) {
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        fig <- plot_power(r, "lattice", language)
        expect_equal(attr(fig, "inferguidR_plot_data")$data, pd$data)
        expect_s3_class(ggplot2::ggplot_build(fig), "ggplot_built")
      }
    }
  }
})

test_that("lattice uses the solved allocation and recorded numerical settings", {
  r <- calculate_power(analysis_spec("boschloo_proportions_two_independent", "a_priori",
    inputs = list(p1 = .9, p2 = .1, alpha = .05, power = .6, n_max = 20,
      allocation_ratio = 1.5, grid_points = 101, nuisance_tol = 1e-8),
    options = list(variant = "minlike")))
  pd <- power_plot_data(r, "lattice")
  expect_identical(pd$status, "ok")
  expect_equal(c(pd$fixed_inputs$n1, pd$fixed_inputs$n2), unname(r$group_sizes))
  expect_identical(pd$fixed_inputs$settings, r$method_outputs$settings)
  expect_equal(pd$region_probabilities$engine, pd$region_probabilities$lattice, tolerance = 1e-12)
})

test_that("lattice refusals are explicit and do not pretend to be distributions", {
  r <- calculate_power(analysis_spec("t_means_one_sample", "post_hoc",
    inputs = list(d = .5, alpha = .05, n_total = 30), options = list(tails = "two")))
  pd <- power_plot_data(r, "lattice")
  expect_identical(pd$status, "unavailable")
  expect_match(pd$reason, "Boschloo")
  expect_equal(nrow(pd$data), 0)
})

test_that("empty rejection regions have a valid lattice and display limits refuse before enumeration", {
  r <- calculate_power(analysis_spec("boschloo_proportions_two_independent", "post_hoc",
    inputs = list(n1 = 2, n2 = 2, p1 = .8, p2 = .2, alpha = .05), options = list(variant = "minlike")))
  pd <- power_plot_data(r, "lattice")
  expect_identical(pd$status, "ok")
  expect_false(any(pd$data$reject))
  expect_equal(pd$region_probabilities$lattice, c(0, 0, 1))
  # Exercise the display guard without running an expensive large-design calculation.
  guard <- .ig_power_boschloo_lattice(list(method_id = .bl_method_id,
    group_sizes = c(n1 = 100L, n2 = 100L)), list(status = "unavailable", data = data.frame()))
  expect_identical(guard$status, "unavailable")
  expect_match(guard$reason, "5000-table display limit")
  expect_equal(nrow(guard$data), 0)
})

test_that("the workbench renders the Boschloo lattice and preserves committed inputs", {
  skip_if_not_installed("shiny"); skip_if_not_installed("ggplot2")
  shiny::testServer(.ig_power_workspace_server, args = list(language = shiny::reactive("de")), {
    session$flushReact()
    session$setInputs(method = "boschloo_proportions_two_independent", mode = "post_hoc", research = TRUE,
      field_boschloo_proportions_two_independent_n1 = "4", field_boschloo_proportions_two_independent_n2 = "5",
      field_boschloo_proportions_two_independent_p1 = "0,8", field_boschloo_proportions_two_independent_p2 = "0,2",
      field_boschloo_proportions_two_independent_alpha = "0,05",
      field_boschloo_proportions_two_independent_tie_tolerance = "1e-7",
      field_boschloo_proportions_two_independent_grid_points = "1001",
      field_boschloo_proportions_two_independent_nuisance_tol = "1e-9",
      option_boschloo_proportions_two_independent_variant = "minlike")
    session$setInputs(calculate = 1)
    if (is.null(state())) stop(message())
    first <- state()$committed$analysis_hash
    expect_identical(current()$status$code, "ok")
    expect_type(output$distribution, "list")
    expect_match(output$distribution$src, "^data:image/png;base64,")
    session$setInputs(field_boschloo_proportions_two_independent_n1 = "6")
    expect_true(state()$dirty)
    expect_identical(state()$committed$analysis_hash, first)
    expect_equal(power_plot_data(current(), "lattice")$fixed_inputs$n1, 4)
  })
})

test_that("German Boschloo descriptions preserve each distinct ordering and direction", {
  expected <- c(minlike = "zweiseitigen Fisher", central = "getrennt maximiert", less = "p1 < p2", greater = "p1 > p2")
  for (variant in names(expected)) {
    r <- calculate_power(analysis_spec("boschloo_proportions_two_independent", "post_hoc",
      inputs = list(n1 = 2, n2 = 2, p1 = .8, p2 = .2, alpha = .05), options = list(variant = variant)))
    de <- .ig_power_display_assumptions(r, "de")
    expect_length(de, 4)
    expect_match(de[2], expected[[variant]], fixed = TRUE)
    expect_match(de[3], "Supremum", fixed = TRUE)
    expect_identical(.ig_power_display_assumptions(r, "en"), r$assumptions)
    expect_true(nzchar(.ig_boschloo_variant_labels("de")[[variant]]))
  }
})
