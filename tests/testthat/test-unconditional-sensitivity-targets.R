test_that("unconditional sensitivity retains its probability target for every effect index", {
  indices <- c("proportion_p1", "difference", "risk_ratio", "odds_ratio")
  for (method in c("exact_proportions_unconditional", "exact_proportions_unconditional_offset")) {
    offset <- grepl("offset$", method)
    for (index in indices) for (direction in c("p_h1_ge_p_h0", "p_h1_le_p_h0")) {
      inputs <- list(p2 = .4, alpha = .05, power = .6, n1 = 30L, n2 = 30L)
      if (offset) inputs <- c(inputs, switch(index,
        proportion_p1 = list(p1_h0 = .5), difference = list(difference_h0 = .1),
        risk_ratio = list(risk_ratio_h0 = 1.25), odds_ratio = list(odds_ratio_h0 = 1.5)))
      options <- list(effect_size_index = index, effect_direction = direction,
        tails = "one", computation = "exact", statistic = if (offset) "score_farrington_manning" else "z_pooled")
      spec <- analysis_spec(method, "sensitivity", inputs = inputs, options = options)
      target <- if (offset) "p1_h1" else "p1"
      expect_identical(spec$target, target)
      expect_identical(.ig_workbench_form(method, "sensitivity", options)$target, target)
      result <- calculate_power(spec)
      expect_identical(result$status$code, "ok")
      expect_equal(result$achieved$power, .6, tolerance = 1e-7)
      p1 <- result$solved[[target]]
      expect_true(if (direction == "p_h1_ge_p_h0") p1 > if (offset) .5 else .4 else p1 < if (offset) .5 else .4)
      expect_identical(ig_to_json(ig_from_json(ig_to_json(result))), ig_to_json(result))
    }
  }
})

test_that("sensitivity target availability does not allow inactive effect inputs in post hoc mode", {
  expect_error(analysis_spec("exact_proportions_unconditional", "post_hoc",
    inputs = list(p1 = .7, p2 = .4, alpha = .05, n1 = 30, n2 = 30),
    options = list(effect_size_index = "odds_ratio", tails = "one", computation = "exact", statistic = "z_pooled")),
    class = "inferguidR_invalid_input")
  expect_error(analysis_spec("exact_proportions_unconditional_offset", "post_hoc",
    inputs = list(p1_h1 = .7, risk_ratio_h0 = 1.25, p2 = .4, alpha = .05, n1 = 30, n2 = 30),
    options = list(effect_size_index = "risk_ratio", tails = "one", computation = "exact", statistic = "score_farrington_manning")),
    class = "inferguidR_invalid_input")
})
