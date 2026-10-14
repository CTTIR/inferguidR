test_that("saturated logistic enumeration does not manufacture full power", {
  for (mu in c(-5, 5)) for (tails in c("one", "two")) {
    spec <- analysis_spec("z_logistic_regression", "post_hoc",
      inputs = list(odds_ratio = 1724.76516217133, p1 = 0.5, alpha = 0.05,
                    n_total = 200, r2_other_x = 0.1, x_mu = mu, x_sigma = 0.1),
      options = list(computation = "enumeration", test_statistic = "likelihood_ratio",
                     x_distribution = "normal", tails = tails))
    result <- calculate_power(spec)
    expect_identical(result$status$code, "ok")
    expect_true(is.finite(result$quantities$ncp))
    expect_gte(result$quantities$ncp, 0)
    expect_lt(result$quantities$ncp, 1e-10)
    # In this nearly constant-response design the deviance tends to zero.
    # The enumeration convention uses size 2 alpha for one tail.
    expect_equal(result$achieved$power, if (tails == "one") 0.1 else 0.05,
                 tolerance = 1e-10)
  }
})

test_that("logistic sensitivity does not select a saturation artefact", {
  spec <- analysis_spec("z_logistic_regression", "sensitivity",
    inputs = list(p1 = 0.5, alpha = 0.05, power = 0.8, n_total = 200,
                  r2_other_x = 0.1, x_mu = 5, x_sigma = 0.1),
    options = list(computation = "enumeration", test_statistic = "likelihood_ratio",
                   x_distribution = "normal", tails = "two", effect_direction = "p2_ge_p1"))
  result <- .ig_glm_sensitivity_search(.ig_get_method(spec$method_id)$power)(
    c(spec$inputs, spec$options), spec$target)
  expect_identical(result$status, "ok")
  expect_false(result$diagnostics$target_reached)
  expect_lt(result$diagnostics$max_power, 0.2)
})
