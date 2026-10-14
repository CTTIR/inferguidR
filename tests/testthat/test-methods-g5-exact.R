# Group G5 exact procedures (bivariate normal correlation, random-model R^2):
# exact distributions, curated G*Power fixture comparisons (cited records),
# option branches and typed outcomes. R/methods-correlation-regression.R.

g5e_ok <- function(value, display) {
  d <- as.numeric(display)
  m <- regmatches(display, regexec("^-?[0-9]*\\.?([0-9]*)(e([-+]?[0-9]+))?$", display))[[1]]
  half <- 0.5 * 10^(-nchar(m[2]) + if (nzchar(m[4])) as.numeric(m[4]) else 0)
  abs(value - d) <= half + 8 * .Machine$double.eps * max(1, abs(d))
}
expect_display_e <- function(value, display, record) {
  expect_true(g5e_ok(value, display), info = sprintf("%s: engine %s, G*Power %s", record, format(value, digits = 12), display))
}
run_e <- function(method, mode, inputs, options = list()) calculate_power(analysis_spec(method, mode, inputs = inputs, options = options))

test_that("the exact CDF of r has the closed-form special cases and the symmetries", {
  # rho = 0: r sqrt(n - 2) / sqrt(1 - r^2) is t(n - 2)
  x <- .4; n <- 12
  expect_equal(.ig_g5_r_tail(x, 0, n, TRUE), pt(x * sqrt(n - 2) / sqrt(1 - x^2), n - 2), tolerance = 1e-13)
  # P(r <= 0) = T_{n-1}(-rho sqrt(n - 1) / sqrt(1 - rho^2))
  rho <- .45; n <- 9
  expect_equal(.ig_g5_r_tail(1e-300, rho, n, TRUE), pt(-rho / sqrt(1 - rho^2) * sqrt(n - 1), n - 1), tolerance = 1e-12)
  # P(r <= x | rho) = P(r >= -x | -rho), and the two tails add up to one
  for (cs in list(c(.3, .6, 15), c(-.2, .5, 40), c(.95, .9, 6))) {
    lo <- .ig_g5_r_tail(cs[1], cs[2], cs[3], TRUE)
    expect_equal(as.double(lo), as.double(.ig_g5_r_tail(-cs[1], -cs[2], cs[3], FALSE)), tolerance = 1e-13)
    expect_equal(as.double(lo + .ig_g5_r_tail(cs[1], cs[2], cs[3], FALSE)), 1, tolerance = 1e-13)
  }
  # m = 1: R^2 is r^2
  expect_equal(.ig_g5_r2_tail(.25, .3^2, 1, 20, TRUE), as.double(.ig_g5_r_tail(.5, .3, 20, TRUE) - .ig_g5_r_tail(-.5, .3, 20, TRUE)), tolerance = 1e-12)
  # rho^2 = 0: central beta
  expect_equal(.ig_g5_r2_tail(.1, 0, 5, 95, FALSE), pbeta(.1, 2.5, 44.5, lower.tail = FALSE), tolerance = 1e-14)
})

test_that("exact correlation reproduces the captured displays (exact and Fisher z branches)", {
  M <- "exact_correlation_bivariate_normal"
  big <- list(exact_n_threshold = 1e7)
  r <- run_e(M, "post_hoc", c(list(rho_h1 = .3, rho_h0 = .8, alpha = .05, n_total = 8), big), list(tails = "two"))   # GPW-EX-CORRBN-POS-B007-O01 (DOC-030)
  expect_display_e(r$quantities$critical[1], "0.3067474", "POS-B007-O01"); expect_display_e(r$quantities$critical[2], "0.9663957", "POS-B007-O01")
  expect_display_e(r$achieved$power, "0.4829271", "POS-B007-O01")
  r <- run_e(M, "post_hoc", list(rho_h1 = .3, rho_h0 = .8, alpha = .05, n_total = 8), list(tails = "two", computation = "fisher_z_approx"))
  expect_display_e(r$achieved$power, "0.4225990", "GPW-EX-CORRBN-POS-B010-O01")
  expect_display_e(r$method_outputs$critical_z[2], "1.9599640", "POS-B010-O01")
  r <- run_e(M, "a_priori", c(list(rho_h1 = .3, rho_h0 = 0, alpha = .05, power = .8), big), list(tails = "one"))     # GPW-EX-CORRBN-APR-B001-O01
  expect_identical(r$solved$n_total, 67L)
  expect_display_e(r$quantities$critical, "0.2026735", "APR-B001-O01"); expect_display_e(r$achieved$power, "0.8032714", "APR-B001-O01")
  r <- run_e(M, "sensitivity", c(list(rho_h0 = 0, alpha = .05, power = .95, n_total = 4), big), list(tails = "two", effect_direction = "r_ge_rho0"))
  expect_display_e(r$method_outputs$rho_h1, "0.9921287", "GPW-EX-CORRBN-SEN-B009-E01")
  expect_display_e(r$quantities$critical[2], "0.9500000", "SEN-B009-E01")
  r <- run_e(M, "sensitivity", c(list(rho_h0 = 0, alpha = .05, power = .8, n_total = 84), big), list(tails = "one", effect_direction = "r_le_rho0"))
  expect_display_e(r$method_outputs$rho_h1, "-0.2676551", "GPW-EX-CORRBN-SEN-B002-O01")
  r <- run_e(M, "criterion", c(list(rho_h1 = .3, rho_h0 = 0, power = .8, n_total = 84), big), list(tails = "two"))  # GPW-EX-CORRBN-CRI-B007-O01
  expect_display_e(r$solved$alpha, "0.0498600", "CRI-B007-O01"); expect_display_e(r$quantities$critical[2], "0.2146951", "CRI-B007-O01")
  r <- run_e(M, "compromise", c(list(rho_h1 = .95, rho_h0 = 0, error_ratio = 4, n_total = 4), big), list(tails = "one"))  # GPW-EX-CORRBN-COM-B001-E01
  expect_display_e(r$solved$alpha, "0.0495899", "COM-B001-E01"); expect_display_e(r$achieved$beta, "0.1983594", "COM-B001-E01")
  expect_display_e(r$quantities$critical, "0.9008203", "COM-B001-E01")
})

test_that("the exact-distribution threshold switches to Fisher z only above N (supplement probes)", {
  M <- "exact_correlation_bivariate_normal"
  base <- list(rho_h1 = .3, rho_h0 = 0, alpha = .05, n_total = 84)
  r50 <- run_e(M, "post_hoc", c(base, exact_n_threshold = 50), list(tails = "one"))   # SUP-G5-EX-CORRBN-POS-P01: Fisher z
  r84 <- run_e(M, "post_hoc", c(base, exact_n_threshold = 84), list(tails = "one"))   # SUP-G5-EX-CORRBN-POS-P03: exact
  expect_identical(r50$method_outputs$branch, "fisher_z")
  expect_display_e(r50$achieved$power, "0.8730282", "POS-P01"); expect_display_e(r50$method_outputs$critical_z, "1.6448536", "POS-P01")
  expect_identical(r84$method_outputs$branch, "exact")
  expect_display_e(r84$achieved$power, "0.8770141", "POS-P03"); expect_display_e(r84$method_outputs$critical_r, "0.1806951", "POS-P03")
  a <- run_e(M, "a_priori", list(rho_h1 = .3, rho_h0 = 0, alpha = .05, power = .8, exact_n_threshold = 60), list(tails = "one"))
  expect_identical(a$solved$n_total, 68L)   # SUP-G5-EX-CORRBN-APR-A01
  expect_identical(a$method_outputs$branch, "fisher_z")
})

test_that("N = 3 is computed with the exact distribution under the Fisher z option (all modes)", {
  M <- "exact_correlation_bivariate_normal"
  fz <- function(tails = "one", ...) list(tails = tails, computation = "fisher_z_approx", ...)
  r <- run_e(M, "post_hoc", list(rho_h1 = .95, rho_h0 = 0, alpha = .05, n_total = 3), fz())                       # SUP-G5-EX-CORRBN-POS-P05
  expect_identical(r$method_outputs$branch, "exact")
  expect_display_e(r$method_outputs$critical_r, "0.9876883", "SUP-G5-EX-CORRBN-POS-P05"); expect_display_e(r$achieved$power, "0.4314480", "POS-P05")
  expect_match(paste(r$assumptions, collapse = " "), "Fisher z approximation does not exist", fixed = TRUE)
  r <- run_e(M, "post_hoc", list(rho_h1 = .95, rho_h0 = 0, alpha = .05, n_total = 3), fz("two"))                  # SUP-G5-EX-CORRBN-POS-P07
  expect_display_e(r$method_outputs$critical_r[1], "-0.9969173", "SUP-G5-EX-CORRBN-POS-P07"); expect_display_e(r$achieved$power, "0.2338768", "POS-P07")
  r <- run_e(M, "post_hoc", list(rho_h1 = .95, rho_h0 = 0, alpha = .05, n_total = 4), fz())                       # SUP-G5-EX-CORRBN-POS-P06
  expect_identical(r$method_outputs$branch, "fisher_z")
  expect_display_e(r$method_outputs$critical_z, "1.6448536", "SUP-G5-EX-CORRBN-POS-P06"); expect_display_e(r$achieved$power, "0.5741411", "POS-P06")
  r <- run_e(M, "sensitivity", list(rho_h0 = 0, alpha = .05, power = .8, n_total = 3), fz(effect_direction = "r_ge_rho0"))
  expect_display_e(r$method_outputs$rho_h1, "0.9931769", "SUP-G5-EX-CORRBN-SEN-S01")
  r <- run_e(M, "criterion", list(rho_h1 = .95, rho_h0 = 0, power = .8, n_total = 3), fz())
  expect_display_e(r$solved$alpha, "0.1417935", "SUP-G5-EX-CORRBN-CRI-C01"); expect_display_e(r$method_outputs$critical_r, "0.9024137", "CRI-C01")
  r <- run_e(M, "compromise", list(rho_h1 = .95, rho_h0 = 0, error_ratio = 1, n_total = 3), fz())
  expect_display_e(r$solved$alpha, "0.1633332", "SUP-G5-EX-CORRBN-COM-M01"); expect_display_e(r$achieved$power, "0.8366668", "COM-M01")
  # A threshold below 3 selects the approximation at N = 3 as well; the exact distribution is used
  # (G*Power displays power = alpha there, SUP-G5-EX-CORRBN-POS-P08, a triaged reference defect).
  r <- run_e(M, "post_hoc", list(rho_h1 = .95, rho_h0 = 0, alpha = .05, n_total = 3, exact_n_threshold = 2), list(tails = "one"))
  expect_identical(r$method_outputs$branch, "exact")
  expect_display_e(r$achieved$power, "0.4314480", "SUP-G5-EX-CORRBN-POS-P05")
})

test_that("a priori with the Fisher z option reports the smallest N whose reported power reaches the target", {
  M <- "exact_correlation_bivariate_normal"
  apr <- function(rho, power, tails = "one") run_e(M, "a_priori", list(rho_h1 = rho, rho_h0 = 0, alpha = .05, power = power),
                                                   list(tails = tails, computation = "fisher_z_approx"))
  a <- apr(.95, .804)                                                                    # SUP-G5-EX-CORRBN-APR-A03
  expect_identical(a$solved$n_total, 5L); expect_display_e(a$achieved$power, "0.8278430", "SUP-G5-EX-CORRBN-APR-A03")
  a <- apr(.97, .8, "two")                                                               # SUP-G5-EX-CORRBN-APR-A04
  expect_identical(a$solved$n_total, 5L); expect_display_e(a$achieved$power, "0.8411004", "SUP-G5-EX-CORRBN-APR-A04")
  a <- apr(.9, .8)                                                                       # SUP-G5-EX-CORRBN-APR-A07
  expect_identical(a$solved$n_total, 6L); expect_display_e(a$achieved$power, "0.8172952", "SUP-G5-EX-CORRBN-APR-A07")
  a <- apr(.999, .8)                                                                     # SUP-G5-EX-CORRBN-APR-A06: N = 3, exact values
  expect_identical(a$solved$n_total, 3L); expect_identical(a$method_outputs$branch, "exact")
  expect_display_e(a$method_outputs$critical_r, "0.9876883", "SUP-G5-EX-CORRBN-APR-A06"); expect_display_e(a$achieved$power, "0.9614610", "APR-A06")
  # GPW-EX-CORRBN-APR-B004-E01 and SUP-G5-EX-CORRBN-APR-A02: G*Power reports N = 4, whose displayed (Fisher z)
  # power 0.5741411 is below the target; the engine reports N = 5 (same values as SUP-G5-EX-CORRBN-APR-A03).
  for (target in c(.8, .803)) {
    a <- apr(.95, target)
    expect_identical(a$solved$n_total, 5L)
    expect_gte(a$achieved$power, target)
  }
  below <- run_e(M, "post_hoc", list(rho_h1 = .95, rho_h0 = 0, alpha = .05, n_total = 4), list(tails = "one", computation = "fisher_z_approx"))
  expect_lt(below$achieved$power, .8)
  # thresholds next to the exact solution N = 67 (SUP-G5-EX-CORRBN-APR-A08, A09)
  thr <- function(x) run_e(M, "a_priori", list(rho_h1 = .3, rho_h0 = 0, alpha = .05, power = .8, exact_n_threshold = x), list(tails = "one"))
  a <- thr(66)
  expect_identical(a$solved$n_total, 68L); expect_display_e(a$achieved$power, "0.8024968", "SUP-G5-EX-CORRBN-APR-A08")
  a <- thr(67)
  expect_identical(a$solved$n_total, 67L); expect_display_e(a$achieved$power, "0.8032714", "SUP-G5-EX-CORRBN-APR-A09")
})

test_that("a priori totals are minimal by enumeration across the exact / Fisher z switch", {
  M <- "exact_correlation_bivariate_normal"
  cases <- list(
    list(ins = list(rho_h1 = .9, rho_h0 = 0, alpha = .05, power = .8), opt = list(tails = "one", computation = "fisher_z_approx")),
    list(ins = list(rho_h1 = -.97, rho_h0 = -.2, alpha = .01, power = .9), opt = list(tails = "two", computation = "fisher_z_approx")),
    list(ins = list(rho_h1 = .8, rho_h0 = .1, alpha = .05, power = .9, exact_n_threshold = 6), opt = list(tails = "two", computation = "exact")),
    list(ins = list(rho_h1 = .85, rho_h0 = 0, alpha = .1, power = .8, exact_n_threshold = 4), opt = list(tails = "one", computation = "exact")),
    list(ins = list(rho_h1 = .995, rho_h0 = .3, alpha = .05, power = .5, exact_n_threshold = 2), opt = list(tails = "one", computation = "exact")))
  for (cs in cases) {
    a <- run_e(M, "a_priori", cs$ins, cs$opt)
    at <- function(n) {
      q <- cs$ins
      q$power <- NULL
      q$n_total <- n
      run_e(M, "post_hoc", q, cs$opt)$achieved$power
    }
    N <- a$solved$n_total
    expect_gte(at(N), cs$ins$power)
    expect_equal(at(N), a$achieved$power, tolerance = 1e-12)
    if (N > 3L) expect_true(all(vapply(seq.int(3L, N - 1L), at, numeric(1)) < cs$ins$power))
  }
})

test_that("random-model R^2 reproduces the captured displays (exact and three-moment F)", {
  M <- "exact_mrc_random_model"
  r <- run_e(M, "post_hoc", list(rho2_h1 = .1, rho2_h0 = 0, alpha = .05, n_total = 95, n_predictors = 5), list(tails = "one"))
  expect_display_e(r$quantities$critical, "0.1151700", "GPW-EX-MRCRND-POS-B001-O01"); expect_display_e(r$achieved$power, "0.6626270", "POS-B001-O01")
  r <- run_e(M, "post_hoc", list(rho2_h1 = .1, rho2_h0 = 0, alpha = .05, n_total = 95, n_predictors = 5), list(tails = "one", distribution = "three_moment_f"))
  expect_display_e(r$achieved$power, "0.6626115", "GPW-EX-MRCRND-POS-B004-O01")
  lee <- run_e(M, "post_hoc", list(rho2_h1 = .4, rho2_h0 = .3, alpha = .05, n_total = 20, n_predictors = 3), list(tails = "one", distribution = "three_moment_f"))
  ex <- run_e(M, "post_hoc", list(rho2_h1 = .4, rho2_h0 = .3, alpha = .05, n_total = 20, n_predictors = 3), list(tails = "one"))
  expect_display_e(lee$quantities$critical, "0.6607579", "SUP-G5-EX-MRCRND-POS-P01"); expect_display_e(lee$achieved$power, "0.1243642", "POS-P01")
  expect_display_e(ex$quantities$critical, "0.6607293", "SUP-G5-EX-MRCRND-POS-P02"); expect_display_e(ex$achieved$power, "0.1243557", "POS-P02")
  a <- run_e(M, "a_priori", list(rho2_h1 = .05, rho2_h0 = .2, alpha = .05, power = .9, n_predictors = 5), list(tails = "one"))
  expect_identical(a$solved$n_total, 153L)  # GPW-EX-MRCRND-APR-B001-O01 (DOC-050)
  expect_display_e(a$quantities$critical, "0.1323088", "APR-B001-O01"); expect_display_e(a$achieved$power, "0.9010511", "APR-B001-O01")
  r <- run_e(M, "sensitivity", list(rho2_h0 = 0, alpha = .05, power = .8, n_total = 95, n_predictors = 5), list(tails = "one", effect_direction = "rho2_h1_ge_h0"))
  expect_display_e(r$method_outputs$rho2_h1, "0.1289187", "GPW-EX-MRCRND-SEN-B003-O01")
  r <- run_e(M, "criterion", list(rho2_h1 = .1, rho2_h0 = 0, power = .8, n_total = 95, n_predictors = 5), list(tails = "two"))
  expect_display_e(r$solved$alpha, "0.2365515", "GPW-EX-MRCRND-CRI-B007-O01")
  expect_display_e(r$quantities$critical[1], "0.0192227", "CRI-B007-O01"); expect_display_e(r$quantities$critical[2], "0.0924550", "CRI-B007-O01")
  r <- run_e(M, "compromise", list(rho2_h1 = .1, rho2_h0 = 0, error_ratio = 1, n_total = 95, n_predictors = 5), list(tails = "two"))
  expect_display_e(r$solved$alpha, "0.2148309", "GPW-EX-MRCRND-COM-B007-O01")
})

test_that("an effect direction without any admissible H1 value is infeasible, not 0", {
  # GPW-EX-MRCRND-SEN-B002-O01: G*Power displays H1 rho^2 = 0 (power alpha); no rho^2 below 0 exists.
  r <- run_e("exact_mrc_random_model", "sensitivity", list(rho2_h0 = 0, alpha = .05, power = .8, n_total = 95, n_predictors = 5),
             list(tails = "one", effect_direction = "rho2_h1_le_h0"))
  expect_identical(r$status$code, "infeasible")
  expect_null(r$solved$z_distance)
  # with H0 rho^2 > 0 the lower side has a solution
  r <- run_e("exact_mrc_random_model", "sensitivity", list(rho2_h0 = .3, alpha = .05, power = .8, n_total = 60, n_predictors = 3),
             list(tails = "two", effect_direction = "rho2_h1_le_h0"))
  expect_identical(r$status$code, "ok")
  expect_lt(r$method_outputs$rho2_h1, .3)
})
