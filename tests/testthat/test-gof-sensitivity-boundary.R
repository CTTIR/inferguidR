test_that("goodness-of-fit sensitivity admits the null effect", {
  for (n in c(1, 20, 10000)) for (df in c(1, 2.5, 20)) {
    r <- calculate_power(analysis_spec("chi2_goodness_of_fit_contingency", "sensitivity",
      list(n_total = n, df = df, alpha = .05, power = .05)))
    expect_identical(r$status$code, "ok")
    expect_identical(r$solved$w, 0)
    expect_equal(r$achieved$power, .05, tolerance = 1e-12)
  }
})

test_that("goodness-of-fit sensitivity obeys noncentrality scaling", {
  # Reference roots checked by 60-digit incomplete-gamma Poisson summation.
  for (n in c(1, 20, 10000)) {
    r <- calculate_power(analysis_spec("chi2_goodness_of_fit_contingency", "sensitivity",
      list(n_total = n, df = 1, alpha = .001, power = .99)))
    expect_identical(r$status$code, "ok")
    expect_equal(r$solved$w, sqrt(31.5492803342789 / n), tolerance = 1e-10)
    expect_equal(r$achieved$power, .99, tolerance = 1e-10)
  }
})

test_that("unresolved near-boundary targets are not represented as exact solutions", {
  r <- calculate_power(analysis_spec("chi2_goodness_of_fit_contingency", "sensitivity",
    list(n_total = 20, df = 2.5, alpha = .05, power = .05 + 1e-14)))
  expect_identical(r$status$code, "numerical_failure")
  expect_null(r$solved$w)
})
