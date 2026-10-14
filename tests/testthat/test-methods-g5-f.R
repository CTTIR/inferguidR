# Group G5 F-test procedures (fixed-model multiple regression): curated
# G*Power fixture comparisons (cited emulator-capture records, rounding
# interval of the displayed digits) and round trips.
# R/methods-correlation-regression.R; full comparisons: validation/compare/<method_id>/.

g5f_ok <- function(value, display) {
  d <- as.numeric(display)
  m <- regmatches(display, regexec("^-?[0-9]*\\.?([0-9]*)(e([-+]?[0-9]+))?$", display))[[1]]
  half <- 0.5 * 10^(-nchar(m[2]) + if (nzchar(m[4])) as.numeric(m[4]) else 0)
  abs(value - d) <= half + 8 * .Machine$double.eps * max(1, abs(d))
}
expect_display_f <- function(value, display, record) {
  expect_true(g5f_ok(value, display), info = sprintf("%s: engine %s, G*Power %s", record, format(value, digits = 12), display))
}
run_f <- function(method, mode, inputs) calculate_power(analysis_spec(method, mode, inputs = inputs))

test_that("R^2 deviation from zero reproduces the captured displays in all five modes", {
  M <- "f_mrc_r2_deviation_from_zero"
  r <- run_f(M, "a_priori", list(f2 = .15, alpha = .05, power = .8, n_predictors = 5))   # GPW-F-MRCR2-APR-B001-O01--r01
  expect_identical(r$solved$n_total, 92L)
  expect_display_f(r$quantities$ncp, "13.8000000", "APR-B001-O01"); expect_display_f(r$quantities$critical, "2.3205293", "APR-B001-O01")
  expect_identical(r$quantities$df, c(5, 86)); expect_display_f(r$achieved$power, "0.8041921", "APR-B001-O01")
  r <- run_f(M, "post_hoc", list(f2 = .1111111, alpha = .05, n_total = 95, n_predictors = 5))   # GPW-F-MRCR2-POS-B001-O02 (DOC-078 inputs)
  expect_display_f(r$quantities$ncp, "10.5555545", "POS-B001-O02"); expect_display_f(r$quantities$critical, "2.3168575", "POS-B001-O02")
  expect_display_f(r$achieved$power, "0.6735857", "POS-B001-O02")
  r <- run_f(M, "sensitivity", list(alpha = .05, power = .8, n_total = 100, n_predictors = 5))  # GPW-F-MRCR2-SEN-B001-O01
  expect_display_f(r$solved$f2, "0.1360036", "SEN-B001-O01"); expect_display_f(r$quantities$ncp, "13.6003596", "SEN-B001-O01")
  r <- run_f(M, "criterion", list(f2 = .15, power = .8, n_total = 100, n_predictors = 5))    # GPW-F-MRCR2-CRI-B001-O01
  expect_display_f(r$solved$alpha, "0.0341296", "CRI-B001-O01"); expect_display_f(r$quantities$critical, "2.5286630", "CRI-B001-O01")
  r <- run_f(M, "compromise", list(f2 = .15, error_ratio = 1, n_total = 100, n_predictors = 5))  # GPW-F-MRCR2-COM-B001-O01
  expect_display_f(r$solved$alpha, "0.0948576", "COM-B001-O01"); expect_display_f(r$achieved$power, "0.9051424", "COM-B001-O01")
  expect_display_f(r$quantities$critical, "1.9405714", "COM-B001-O01")
  r <- run_f(M, "post_hoc", list(f2 = 2, alpha = .05, n_total = 3, n_predictors = 1))        # GPW-F-MRCR2-POS-B001-E01 (df2 = 1)
  expect_display_f(r$quantities$critical, "161.4476", "POS-B001-E01"); expect_display_f(r$achieved$power, "0.1526964", "POS-B001-E01")
})

test_that("R^2 increase reproduces the captured displays", {
  M <- "f_mrc_r2_increase"
  r <- run_f(M, "a_priori", list(f2 = .0714286, alpha = .01, power = .8, n_tested = 4, n_predictors = 9))  # GPW-F-MRCINC-APR-B001-O01
  expect_identical(r$solved$n_total, 242L)
  expect_display_f(r$quantities$ncp, "17.2857212", "APR-B001-O01"); expect_display_f(r$quantities$critical, "3.4010085", "APR-B001-O01")
  expect_identical(r$quantities$df, c(4, 232)); expect_display_f(r$achieved$power, "0.8015724", "APR-B001-O01")
  r <- run_f(M, "post_hoc", list(f2 = .0714286, alpha = .01, n_total = 90, n_tested = 4, n_predictors = 9))  # GPW-F-MRCINC-POS-B001-O01 (DOC-086)
  expect_display_f(r$quantities$ncp, "6.4285740", "POS-B001-O01"); expect_display_f(r$quantities$critical, "3.5631096", "POS-B001-O01")
  expect_display_f(r$achieved$power, "0.2412966", "POS-B001-O01")
  r <- run_f(M, "sensitivity", list(alpha = .05, power = .8, n_total = 100, n_tested = 2, n_predictors = 5))   # GPW-F-MRCINC-SEN-B001-O01
  expect_display_f(r$solved$f2, "0.0994846", "SEN-B001-O01")
  r <- run_f(M, "criterion", list(f2 = .15, power = .8, n_total = 100, n_tested = 2, n_predictors = 5))       # GPW-F-MRCINC-CRI-B001-O01
  expect_display_f(r$solved$alpha, "0.008647099", "CRI-B001-O01"); expect_display_f(r$quantities$critical, "4.9989090", "CRI-B001-O01")
  r <- run_f(M, "a_priori", list(f2 = .018870, alpha = .05, power = .8, n_tested = 1, n_predictors = 3))       # GPW-F-MRCINC-APR-B001-O02 (DOC-093)
  expect_identical(r$solved$n_total, 418L)
})

test_that("the noncentral F shortcut for huge noncentrality is exact to double precision", {
  # The bound of .ig_g5_p_f returns 1 / 0 only when the lower tail is provably below 1e-300.
  expect_identical(.ig_g5_p_f(2.3, 5, 94, 1e9, lower_tail = FALSE), 1)
  expect_identical(.ig_g5_p_f(2.3, 5, 94, 1e9, lower_tail = TRUE), 0)
  expect_equal(.ig_g5_p_f(2.3, 5, 94, 10, lower_tail = FALSE), .ig_p_f(2.3, 5, 94, 10, lower_tail = FALSE))
  # A sensitivity scan over f2 up to 1e4 (lambda 1e9) stays fast.
  t0 <- proc.time()[["elapsed"]]
  r <- run_f("f_mrc_r2_deviation_from_zero", "sensitivity", list(alpha = .001, power = .99, n_total = 1e5, n_predictors = 100))
  expect_identical(r$status$code, "ok")
  expect_display_f(r$solved$f2, "0.001030932", "GPW-F-MRCR2-SEN-B001-D01")
  expect_lt(proc.time()[["elapsed"]] - t0, 20)
})
