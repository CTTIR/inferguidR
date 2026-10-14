# Effect-size drawers of the G2 methods (R/effect-anova.R). Displayed values are quoted from the
# drawer captures (validation/fixtures/emulator_capture/GPW-UT-DRAW-UTL-U013..U018) and the
# documented examples DOC-057, DOC-058, DOC-059, DOC-079 (rounding-interval comparison).

g2d_ok <- function(ours, shown) {
  x <- as.numeric(shown)
  dec <- if (grepl(".", shown, fixed = TRUE)) nchar(sub("^[^.]*\\.", "", shown)) else 0
  abs(ours - x) <= 0.5 * 10^-dec + 8 * .Machine$double.eps * max(1, abs(x))
}

test_that("ANOVA drawers reproduce the captured and documented values", {
  es <- effect_size("f_anova_oneway", "variances", var_effect = 1.10059172, var_error = 36)
  expect_true(g2d_ok(es$value, "0.1748485"))                                   # U013-UO01, DOC-058
  out <- .ig_drawer_f_anova_oneway("variances", list(var_effect = 1.10059172, var_error = 36))$outputs
  expect_true(g2d_ok(out$partial_eta_squared, "0.02966507"))
  expect_true(g2d_ok(effect_size("f_anova_oneway", "eta_squared", partial_eta_squared = 0.999)$value, "31.60696"))  # U014-UE01
  expect_true(g2d_ok(effect_size("f_anova_oneway", "eta_squared", partial_eta_squared = 0.000001)$value, "0.001000001"))  # U014-UD01
  es <- effect_size("f_anova_special_main_interactions", "variances", var_effect = 100, var_error = 0.0001)
  expect_true(g2d_ok(es$value, "1000"))                                        # U015-UE01
  es <- effect_size("f_anova_special_main_interactions", "variances", var_effect = 0.85546, var_error = 1.71296)
  expect_true(g2d_ok(es$value, "0.7066856"))                                   # U015-UO01, DOC-059
  expect_true(g2d_ok(effect_size("f_ancova_main_interactions", "variances", var_effect = 0.0625, var_error = 1)$value, "0.25"))  # U017-UO01
  expect_true(g2d_ok(effect_size("f_ancova_main_interactions", "eta_squared", partial_eta_squared = 0.06)$value, "0.2526456"))  # U018-UO01
  d57 <- .ig_drawer_f_anova_oneway("means", list(means = c(2, 3, 4, 5), sizes = c(2, 3, 4, 4), sd = 6))
  expect_true(g2d_ok(d57$value, "0.1748485"))                                  # DOC-057
  expect_identical(d57$outputs$n_total, 13)
  expect_true(g2d_ok(effect_size("f_anova_oneway", "means", means = c(2, 3, 2, 5), sizes = c(5, 6, 6, 5), sd = 2)$value,
                     "0.5930904"))                                             # DOC-079
  eq <- effect_size("f_anova_oneway", "means", means = c(0, 1), n_per_group = 10, sd = 1)
  expect_equal(eq$value, 0.5)
})

test_that("drawer inputs are validated", {
  code <- function(expr) th_condition_code(expr)
  expect_identical(code(effect_size("f_anova_oneway", "variances", var_effect = 1, var_error = 0)), "invalid_input")
  expect_identical(code(effect_size("f_anova_oneway", "eta_squared", partial_eta_squared = 1)), "invalid_input")
  expect_identical(code(effect_size("f_anova_oneway", "means", means = c(1, 2), sd = 1)), "invalid_input")
  expect_identical(code(effect_size("f_anova_oneway", "means", means = c(1, 2), sizes = c(1.5, 2), sd = 1)), "invalid_input")
  expect_identical(code(effect_size("f_anova_special_main_interactions", "means", means = c(1, 2), n_per_group = 2, sd = 1)),
                   "invalid_input")
  expect_identical(code(effect_size("f_ancova_main_interactions", "variances", var_effect = 1)), "invalid_input")
})

test_that("Hotelling drawers compute the Mahalanobis distance from a covariance or SDs and correlations", {
  S <- matrix(c(4, 1.2, 1.2, 1), 2)
  es <- effect_size("f_hotelling_t2_one_group", "covariance", means = c(1, 0.5), null_means = c(0, 0), covariance = S)
  d <- c(1, 0.5)
  expect_equal(es$value, sqrt(drop(t(d) %*% solve(S) %*% d)), tolerance = 1e-14)
  es2 <- effect_size("f_hotelling_t2_two_groups", "sd_correlation", means_1 = c(1, 0.5), means_2 = c(0, 0),
                     sds = c(2, 1), correlation = matrix(c(1, 0.6, 0.6, 1), 2))
  expect_equal(es2$value, es$value, tolerance = 1e-14)
  # p = 1 reduces to |d| / sigma
  expect_equal(effect_size("f_hotelling_t2_two_groups", "covariance", means_1 = 3, means_2 = 1,
                           covariance = matrix(4))$value, 1)
  code <- function(expr) th_condition_code(expr)
  expect_identical(code(effect_size("f_hotelling_t2_one_group", "covariance", means = c(1, 0), null_means = c(0, 0),
                                    covariance = matrix(c(1, 2, 2, 1), 2))), "invalid_input")  # not positive definite
  expect_identical(code(effect_size("f_hotelling_t2_two_groups", "sd_correlation", means_1 = c(1, 0), means_2 = c(0, 0),
                                    sds = c(1, 1), correlation = matrix(c(1, 0.5, 0.4, 1), 2))), "invalid_input")  # asymmetric
})
