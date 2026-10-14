test_that("public effect sizes preserve derived and supporting quantities", {
  anova <- effect_size("f_anova_oneway", "variances", var_effect = 2, var_error = 8)
  expect_identical(anova$status$code, "ok")
  expect_equal(anova$value, 0.5)
  expect_equal(anova$outputs$partial_eta_squared, 0.2)
  expect_identical(anova$details, list())

  corr <- effect_size("z_correlations_two_independent", "two_correlations", r1 = 0.2, r2 = 0.5)
  expect_identical(corr$status$code, "ok")
  expect_equal(corr$details$z1, atanh(0.2))
  expect_equal(corr$details$z2, atanh(0.5))
  expect_identical(corr$outputs, list())

  simple <- effect_size("t_means_two_independent", "common_sd", mean_1 = 1, mean_2 = 0, sd = 2)
  expect_identical(simple$outputs, list())
  expect_identical(simple$details, list())
})
