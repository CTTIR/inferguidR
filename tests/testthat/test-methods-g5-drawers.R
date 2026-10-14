# Group G5 effect-size drawers (R/effect-correlation-regression.R): captured
# drawer outputs (util_effect_size_drawers records and documented examples)
# and input validation.

g5d_ok <- function(value, display) {
  d <- as.numeric(display)
  m <- regmatches(display, regexec("^-?[0-9]*\\.?([0-9]*)(e([-+]?[0-9]+))?$", display))[[1]]
  half <- 0.5 * 10^(-nchar(m[2]) + if (nzchar(m[4])) as.numeric(m[4]) else 0)
  abs(value - d) <= half + 8 * .Machine$double.eps * max(1, abs(d))
}
expect_display_d <- function(value, display, record) {
  expect_true(g5d_ok(value, display), info = sprintf("%s: engine %s, G*Power %s", record, format(value, digits = 12), display))
}
drawer <- function(method, mode, ...) get(paste0(".ig_drawer_", method), asNamespace("inferguidR"))(mode, list(...))

test_that("correlation drawers take the positive root of r^2", {
  es <- effect_size("t_correlation_point_biserial", "r_squared", r_squared = .5)
  expect_identical(es$status$code, "ok")
  expect_display_d(es$value, "0.7071068", "GPW-UT-DRAW-UTL-U054-UO01")
  expect_display_d(effect_size("exact_correlation_bivariate_normal", "r_squared", r_squared = .999999)$value, "0.9999995",
                   "GPW-UT-DRAW-UTL-U001-UE01")
})

test_that("one-group slope drawer modes", {
  M <- "t_linreg_one_group_slope"
  o <- drawer(M, "rho_residual_sdx", rho = .5, sd_residual = 1, sd_x = 1)          # GPW-UT-DRAW-UTL-U055-UO01
  expect_display_d(o$outputs$sd_y, "1.154701", "U055-UO01"); expect_display_d(o$value, "0.5773503", "U055-UO01")
  o <- drawer(M, "rho_slope_sdx", rho = -.99, slope = -2, sd_x = 2)                 # U056-UE01
  expect_display_d(o$outputs$sd_y, "4.040404", "U056-UE01")
  expect_display_d(drawer(M, "rho_sdx_sdy", rho = .99, sd_x = .1, sd_y = 10)$value, "99", "U057-UE01")
  expect_display_d(drawer(M, "residual_slope_sdx", sd_residual = 1, slope = .5, sd_x = 1)$outputs$sd_y, "1.118034", "U058-UO01")
  expect_display_d(drawer(M, "residual_sdx_sdy", sd_residual = .8, sd_x = 1, sd_y = 1)$value, "0.6", "U059-UO01")
  expect_error(drawer(M, "residual_sdx_sdy", sd_residual = 2, sd_x = 1, sd_y = 1), class = "inferguidR_error")
})

test_that("two-group drawers pool the residual SD with the allocation ratio", {
  o <- drawer("t_linreg_two_groups_slopes", "sdx_sdy_slope", sd_x1 = 9.02914, sd_y1 = .66942, slope1 = -.04653,
              sd_x2 = 11.86779, sd_y2 = .68435, slope2 = -.03061, allocation_ratio = 1.571428)   # DOC-101, GPW-UT-DRAW-UTL-U061-UO01
  expect_display_d(o$outputs$rho1, "-0.6275969", "U061-UO01"); expect_display_d(o$outputs$rho2, "-0.5308293", "U061-UO01")
  expect_display_d(o$outputs$sd_residual, "0.5578413", "U061-UO01"); expect_display_d(o$value, "0.01592", "U061-UO01")
  o <- drawer("t_linreg_two_groups_slopes", "sdx_rho_slope", sd_x1 = 9, rho1 = -.6, slope1 = -.05, sd_x2 = 12, rho2 = -.5,
              slope2 = -.03, allocation_ratio = 1.5)                                               # U062-UO01
  expect_display_d(o$outputs$sd_residual, "0.6142312", "U062-UO01")
  expect_equal(o$outputs$sd_y1, 0.75)   # G*Power displays -0.05 in this field (drawer display defect, triaged)
  o <- drawer("t_linreg_two_groups_intercepts", "sdx_rho_slope", sd_x1 = 1, rho = .99, slope1 = 1, intercept1 = 0, intercept2 = 10)
  expect_display_d(o$outputs$sd_residual, "0.1424923", "GPW-UT-DRAW-UTL-U065-UE01"); expect_display_d(o$value, "10", "U065-UE01")
})

test_that("f2 drawers of the fixed-model regression tests", {
  o <- drawer("f_mrc_r2_increase", "variances", var_explained = .05, var_residual = .7)            # GPW-UT-DRAW-UTL-U041-UO01
  expect_display_d(o$outputs$partial_r2, "0.06666667", "U041-UO01"); expect_display_d(o$value, "0.07142857", "U041-UO01")
  expect_display_d(drawer("t_mrc_single_coefficient", "direct", partial_r2 = .06)$value, "0.06382979", "GPW-UT-DRAW-UTL-U067-UO01")
  expect_display_d(effect_size("f_mrc_r2_deviation_from_zero", "r_squared", rho2 = .10)$value, "0.1111111", "GPW-UT-DRAW-UTL-U039-UO02")
  B <- matrix(c(1, 0, .5, 0, 0, 1, 0, .2, .5, 0, 1, 0, 0, .2, 0, 1), 4)
  o <- drawer("f_mrc_r2_deviation_from_zero", "predictor_correlations", u = c(.3, 0, -.2, .2), B = B)
  expect_display_d(o$outputs$rho2, "0.295000", "DOC-084")
  expect_error(drawer("f_mrc_r2_deviation_from_zero", "predictor_correlations", u = c(.3, .3), B = matrix(c(1, 1.2, 1.2, 1), 2)),
               class = "inferguidR_error")
})

test_that("random-model drawer: confidence limits of rho^2 and the relative position", {
  o <- drawer("exact_mrc_random_model", "confidence_interval", n_total = 50, n_predictors = 5, r2_observed = .3, confidence = .95, position = 0)
  expect_display_d(o$outputs$ci_lower, "0.03371056", "GPW-UT-DRAW-UTL-U010-UO01")
  expect_display_d(o$outputs$ci_upper, "0.4603321", "U010-UO01")
  expect_display_d(o$outputs$bound_lower, "0.05894807", "U010-UO01")
  expect_display_d(o$outputs$bound_upper, "0.4245392", "U010-UO01")
  expect_display_d(o$value, "0.03371056", "U010-UO01")
  half <- drawer("exact_mrc_random_model", "confidence_interval", n_total = 50, n_predictors = 5, r2_observed = .3, confidence = .95, position = .5)
  expect_equal(half$value, (o$outputs$ci_lower + o$outputs$ci_upper) / 2, tolerance = 1e-14)
  B <- matrix(c(1, 0, .5, 0, 0, 1, 0, .2, .5, 0, 1, 0, 0, .2, 0, 1), 4)
  expect_display_d(drawer("exact_mrc_random_model", "predictor_correlations", u = c(.3, .1, -.2, .2), B = B)$value, "0.297083", "DOC-053")
})
