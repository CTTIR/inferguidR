# Group G5 t-test procedures: curated G*Power fixture comparisons (displayed
# values of the cited emulator-capture records, compared within their rounding
# interval), mode round trips and the two-group a priori rule.
# R/methods-correlation-regression.R; full comparisons: validation/compare/<method_id>/.

g5t_ok <- function(value, display) {
  d <- as.numeric(display)
  m <- regmatches(display, regexec("^-?[0-9]*\\.?([0-9]*)(e([-+]?[0-9]+))?$", display))[[1]]
  half <- 0.5 * 10^(-nchar(m[2]) + if (nzchar(m[4])) as.numeric(m[4]) else 0)
  abs(value - d) <= half + 8 * .Machine$double.eps * max(1, abs(d))
}
expect_display <- function(value, display, record) {
  expect_true(g5t_ok(value, display), info = sprintf("%s: engine %s, G*Power %s", record, format(value, digits = 12), display))
}
run <- function(method, mode, inputs, options = list()) calculate_power(analysis_spec(method, mode, inputs = inputs, options = options))
crit_shown <- function(r) { cr <- r$quantities$critical; if (length(cr) == 1L) cr else if (r$quantities$ncp < 0) min(cr) else max(cr) }

test_that("point biserial reproduces the captured G*Power displays in all five modes", {
  M <- "t_correlation_point_biserial"
  r <- run(M, "a_priori", list(rho = .3, alpha = .05, power = .95), list(tails = "one"))   # GPW-T-PBIS-APR-B001-O01--r01
  expect_identical(r$solved$n_total, 111L)
  expect_display(r$quantities$ncp, "3.3133098", "APR-B001-O01"); expect_display(crit_shown(r), "1.6589535", "APR-B001-O01")
  expect_display(r$achieved$power, "0.9503016", "APR-B001-O01")
  r <- run(M, "post_hoc", list(rho = .5, alpha = .05, n_total = 12), list(tails = "two"))  # GPW-T-PBIS-POS-B002-O01--r01
  expect_display(r$quantities$ncp, "2.0000000", "POS-B002-O01"); expect_display(crit_shown(r), "2.2281389", "POS-B002-O01")
  expect_display(r$achieved$power, "0.4396273", "POS-B002-O01")
  r <- run(M, "sensitivity", list(alpha = .05, power = .8, n_total = 84), list(tails = "one"))  # GPW-T-PBIS-SEN-B001-O01--r01
  expect_display(r$solved$rho, "0.2638678", "SEN-B001-O01"); expect_display(r$quantities$ncp, "2.5072484", "SEN-B001-O01")
  r <- run(M, "criterion", list(rho = .3, power = .8, n_total = 84), list(tails = "two"))   # GPW-T-PBIS-CRI-B002-O01--r01
  expect_display(r$solved$alpha, "0.0449417", "CRI-B002-O01"); expect_display(crit_shown(r), "2.0363599", "CRI-B002-O01")
  r <- run(M, "compromise", list(rho = .3, error_ratio = 1, n_total = 84), list(tails = "one"))  # GPW-T-PBIS-COM-B001-O01--r01
  expect_display(r$solved$alpha, "0.0760595", "COM-B001-O01"); expect_display(r$achieved$power, "0.9239405", "COM-B001-O01")
  expect_display(crit_shown(r), "1.4455342", "COM-B001-O01")
  r <- run(M, "post_hoc", list(rho = .999, alpha = .05, n_total = 3), list(tails = "two"))  # GPW-T-PBIS-POS-B002-E01--r01 (|ncp| > 37.62)
  expect_display(r$achieved$power, "0.9976060", "POS-B002-E01")
})

test_that("one-group slope uses the residual SD under H1, also for Slope H0 != 0", {
  M <- "t_linreg_one_group_slope"
  r <- run(M, "post_hoc", list(slope_h1 = -.0667, slope_h0 = 0, sd_x = 7.5, sd_y = 4, alpha = .05, n_total = 100), list(tails = "two"))
  expect_display(r$quantities$ncp, "-1.2605215", "GPW-T-SLOPE1-POS-B004-O01"); expect_display(crit_shown(r), "-1.9844675", "POS-B004-O01")
  expect_display(r$achieved$power, "0.2389693", "POS-B004-O01")
  r <- run(M, "post_hoc", list(slope_h1 = .5, slope_h0 = .2, sd_x = 1, sd_y = 1, alpha = .05, n_total = 50), list(tails = "two"))
  expect_display(r$quantities$ncp, "2.4494897", "SUP-G5-T-SLOPE1-POS-P01"); expect_display(r$achieved$power, "0.6701416", "POS-P01")
  r <- run(M, "post_hoc", list(slope_h1 = -.2, slope_h0 = .3, sd_x = 2, sd_y = 3, alpha = .05, n_total = 40), list(tails = "one"))
  expect_display(r$quantities$ncp, "-2.1271781", "SUP-G5-T-SLOPE1-POS-P02"); expect_display(crit_shown(r), "-1.6859545", "POS-P02")
  r <- run(M, "a_priori", list(slope_h1 = .3, sd_x = 1, sd_y = 1, alpha = .05, power = .8), list(tails = "one"))
  expect_identical(r$solved$n_total, 64L)   # GPW-T-SLOPE1-APR-B001-O01--r01
  expect_display(r$quantities$ncp, "2.5158836", "APR-B001-O01"); expect_display(r$achieved$power, "0.8005036", "APR-B001-O01")
  r <- run(M, "criterion", list(slope_h1 = .3, sd_x = 1, sd_y = 1, power = .8, n_total = 50), list(tails = "one"))
  expect_display(r$solved$alpha, "0.0868317", "GPW-T-SLOPE1-CRI-B001-O01"); expect_display(crit_shown(r), "1.3810384", "CRI-B001-O01")
  r <- run(M, "compromise", list(slope_h1 = .3, sd_x = 1, sd_y = 1, error_ratio = 1, n_total = 50), list(tails = "two"))
  expect_display(r$solved$alpha, "0.1876307", "GPW-T-SLOPE1-COM-B004-O01"); expect_display(crit_shown(r), "1.3366733", "COM-B004-O01")
})

test_that("one-group slope sensitivity solves on the side of the effect direction option", {
  M <- "t_linreg_one_group_slope"
  base <- list(sd_x = 1, sd_y = 1, alpha = .05, power = .8, n_total = 50)
  lt <- run(M, "sensitivity", c(base, slope_h0 = .2), list(tails = "one", effect_direction = "slope_lt_h0"))
  gt <- run(M, "sensitivity", c(base, slope_h0 = .2), list(tails = "one", effect_direction = "slope_gt_h0"))
  expect_display(lt$method_outputs$slope_h1, "-0.1525313", "SUP-G5-T-SLOPE1-SEN-S02")
  expect_display(gt$method_outputs$slope_h1, "0.5073807", "SUP-G5-T-SLOPE1-SEN-S01")
  expect_display(gt$quantities$ncp, "2.5222869", "SEN-S01"); expect_display(lt$quantities$ncp, "-2.5222869", "SEN-S02")
  big <- run(M, "sensitivity", list(slope_h0 = .5, sd_x = 1, sd_y = 1, alpha = .05, power = .9, n_total = 2000),
             list(tails = "two", effect_direction = "slope_gt_h0"))
  expect_display(big$method_outputs$slope_h1, "0.5600763", "SUP-G5-T-SLOPE1-SEN-S03")
  # Slope H0 = 0: the two directions are mirror images (N = 3, GPW-T-SLOPE1-SEN-B002/B003-E01).
  e <- lapply(c("slope_lt_h0", "slope_gt_h0"), function(d) run(M, "sensitivity", list(slope_h0 = 0, sd_x = 1, sd_y = 1, alpha = .05,
                                                                                         power = .95, n_total = 3), list(tails = "one", effect_direction = d)))
  expect_display(e[[2]]$method_outputs$slope_h1, "0.9905792", "SEN-B003-E01")
  expect_equal(e[[1]]$method_outputs$slope_h1, -e[[2]]$method_outputs$slope_h1, tolerance = 1e-10)
})

test_that("two-group slopes: post hoc formula and G*Power's a priori total and split", {
  M <- "t_linreg_two_groups_slopes"
  r <- run(M, "post_hoc", list(delta_slope = .5, sd_residual = 1, sd_x1 = 1, sd_x2 = 1, alpha = .05, n1 = 50, n2 = 50), list(tails = "two"))
  expect_display(r$quantities$ncp, "2.5000000", "GPW-T-SLOPE2-POS-B002-O01"); expect_display(crit_shown(r), "1.9849843", "POS-B002-O01")
  expect_display(r$achieved$power, "0.6967142", "POS-B002-O01")
  p02 <- run(M, "post_hoc", list(delta_slope = .01592, sd_residual = .5578413, sd_x1 = 9.02914, sd_x2 = 11.86779, alpha = .05, n1 = 163, n2 = 256),
             list(tails = "two"))
  expect_display(p02$quantities$ncp, "2.8121707", "SUP-G5-T-SLOPE2-POS-P02"); expect_display(p02$achieved$power, "0.8011390", "POS-P02")
  # DOC-102 / GPW-T-SLOPE2-APR-B002-O01: N 419 = 163 + 256, df 415; the engine reports the values of that
  # integer design (G*Power's post hoc P02), not the a priori display (2.8115979 / 0.8009795, n1 replaced by n2/r).
  a <- run(M, "a_priori", list(delta_slope = .01592, sd_residual = .5578413, sd_x1 = 9.02914, sd_x2 = 11.86779, alpha = .05, power = .8,
                               allocation_ratio = 1.571428), list(tails = "two"))
  expect_identical(a$solved$n_total, 419L)
  expect_equal(unname(a$group_sizes), c(163, 256))
  expect_identical(a$quantities$df, 415)
  expect_equal(a$achieved$power, p02$achieved$power, tolerance = 1e-12)
  expect_identical(a$diagnostics$n_fractional_even, 418)
  # SUP-G5-T-SLOPE2-APR-A05: N 16 (12 / 4) from the step-up check (the even-N rule alone gives 14).
  a5 <- run(M, "a_priori", list(delta_slope = 2, sd_residual = 1, sd_x1 = 1, sd_x2 = 1, alpha = .05, power = .8, allocation_ratio = .3),
            list(tails = "two"))
  expect_identical(a5$solved$n_total, 16L)
  expect_equal(unname(a5$group_sizes), c(12, 4))
  # SUP-G5-T-SLOPE2-APR-A04: G*Power's design 2 / 12 (N 14) has power 0.785 < 0.8; the engine reports N 15.
  a4 <- run(M, "a_priori", list(delta_slope = 2, sd_residual = 1, sd_x1 = 1, sd_x2 = 1, alpha = .05, power = .8, allocation_ratio = 5),
            list(tails = "one"))
  expect_identical(a4$diagnostics$n_gpower_rule, 14)
  expect_identical(a4$solved$n_total, 15L)
  expect_lt(a4$diagnostics$gpower_design_power, .8)
  expect_gte(a4$achieved$power, .8)
  r <- run(M, "sensitivity", list(sd_residual = 1, sd_x1 = 1, sd_x2 = 1, alpha = .05, power = .8, n1 = 50, n2 = 50), list(tails = "two"))
  expect_display(r$solved$delta_slope, "0.5659998", "GPW-T-SLOPE2-SEN-B002-O01")
})

test_that("two-group intercepts: post hoc formula and the even-N a priori rule", {
  M <- "t_linreg_two_groups_intercepts"
  r <- run(M, "post_hoc", list(delta_intercept = 2, sd_residual = .5, sd_x1 = .5, sd_x2 = 3, mean_x1 = -2, mean_x2 = 3, alpha = .05, n1 = 3, n2 = 3),
           list(tails = "one"))
  expect_display(r$quantities$ncp, "1.5894388", "GPW-T-INTER2-POS-B001-E01"); expect_display(crit_shown(r), "2.9199856", "POS-B001-E01")
  expect_display(r$achieved$power, "0.2900631", "POS-B001-E01")
  r <- run(M, "post_hoc", list(delta_intercept = 2, sd_residual = .5, sd_x1 = .5, sd_x2 = 3, mean_x1 = -2, mean_x2 = 3, alpha = .05, n1 = 7, n2 = 29),
           list(tails = "one"))
  expect_display(r$quantities$ncp, "2.5310693", "SUP-G5-T-INTER2-POS-P01"); expect_display(r$achieved$power, "0.7973281", "POS-P01")
  a <- run(M, "a_priori", list(delta_intercept = .5, sd_residual = 1, sd_x1 = 1, sd_x2 = 1, mean_x1 = 0, mean_x2 = 0, alpha = .05, power = .8),
           list(tails = "two"))
  expect_identical(a$solved$n_total, 128L)   # GPW-T-INTER2-APR-B002-O01--r01
  expect_display(a$quantities$ncp, "2.8284271", "APR-B002-O01"); expect_display(a$achieved$power, "0.8013621", "APR-B002-O01")
  # GPW-T-INTER2-APR-B001-E01: G*Power reports N 36 (7 / 29, power 0.797 by its own post hoc); the engine reports N 38.
  e <- run(M, "a_priori", list(delta_intercept = 2, sd_residual = .5, sd_x1 = .5, sd_x2 = 3, mean_x1 = -2, mean_x2 = 3, alpha = .05, power = .8,
                               allocation_ratio = 4), list(tails = "one"))
  expect_identical(e$diagnostics$n_gpower_rule, 36)
  expect_identical(e$solved$n_total, 38L)
  r <- run(M, "criterion", list(delta_intercept = 2, sd_residual = .5, sd_x1 = .5, sd_x2 = 3, mean_x1 = -2, mean_x2 = 3, power = .95, n1 = 3, n2 = 3),
           list(tails = "two"))
  expect_display(r$solved$alpha, "0.8295436", "GPW-T-INTER2-CRI-B002-E01"); expect_display(crit_shown(r), "0.2446420", "CRI-B002-E01")
})

test_that("single regression coefficient reproduces the captured displays", {
  M <- "t_mrc_single_coefficient"
  r <- run(M, "post_hoc", list(f2 = .15, alpha = .05, n_total = 100, n_predictors = 5), list(tails = "one"))
  expect_display(r$quantities$ncp, "3.8729833", "GPW-T-MRC1-POS-B001-O01"); expect_display(crit_shown(r), "1.6612259", "POS-B001-O01")
  expect_identical(r$quantities$df, 94); expect_display(r$achieved$power, "0.9860980", "POS-B001-O01")
  r <- run(M, "a_priori", list(f2 = 1.5, alpha = .05, power = .8, n_predictors = 1), list(tails = "two"))
  expect_identical(r$solved$n_total, 8L)     # GPW-T-MRC1-APR-B002-E01--r01
  expect_display(r$quantities$ncp, "3.4641016", "APR-B002-E01"); expect_display(r$achieved$power, "0.8212779", "APR-B002-E01")
  r <- run(M, "sensitivity", list(alpha = .05, power = .8, n_total = 100, n_predictors = 5), list(tails = "two"))
  expect_display(r$solved$f2, "0.0801237", "GPW-T-MRC1-SEN-B002-O01")
  r <- run(M, "compromise", list(f2 = .15, error_ratio = 1, n_total = 100, n_predictors = 5), list(tails = "two"))
  expect_display(r$solved$alpha, "0.0389046", "GPW-T-MRC1-COM-B002-O01"); expect_display(crit_shown(r), "2.0945120", "COM-B002-O01")
})

test_that("inverse modes recover the forward model (round trips)", {
  M <- "t_linreg_one_group_slope"
  s <- run(M, "sensitivity", list(slope_h0 = -.1, sd_x = 2, sd_y = 3, alpha = .01, power = .9, n_total = 80),
           list(tails = "two", effect_direction = "slope_lt_h0"))
  ph <- run(M, "post_hoc", list(slope_h1 = s$method_outputs$slope_h1, slope_h0 = -.1, sd_x = 2, sd_y = 3, alpha = .01, n_total = 80), list(tails = "two"))
  expect_equal(ph$achieved$power, .9, tolerance = 1e-9)
  expect_lt(s$method_outputs$slope_h1, -.1)
  c1 <- run("t_linreg_two_groups_slopes", "criterion", list(delta_slope = .3, sd_residual = 1, sd_x1 = 1, sd_x2 = 2, power = .7, n1 = 30, n2 = 45),
            list(tails = "one"))
  p1 <- run("t_linreg_two_groups_slopes", "post_hoc", list(delta_slope = .3, sd_residual = 1, sd_x1 = 1, sd_x2 = 2, alpha = c1$solved$alpha,
                                                           n1 = 30, n2 = 45), list(tails = "one"))
  expect_equal(p1$achieved$power, .7, tolerance = 1e-9)
  a <- run("t_correlation_point_biserial", "a_priori", list(rho = .25, alpha = .01, power = .9), list(tails = "two"))
  below <- run("t_correlation_point_biserial", "post_hoc", list(rho = .25, alpha = .01, n_total = a$solved$n_total - 1L), list(tails = "two"))
  expect_lt(below$achieved$power, .9)
  expect_gte(a$achieved$power, .9)
})
