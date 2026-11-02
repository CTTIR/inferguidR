# Tail-safe distribution helpers (R/num-dist.R) against the independent
# reference primitives of helper-reference.R (decision D18).
#
# Tolerances: probabilities use the VALIDATION_CONTRACT starting rule
# abs 1e-8 + rel 1e-6 (validation/tolerances.md section 5). Tail
# probabilities below 1e-8 are compared on the log scale with
# |log(ours) - log(ref)| <= 1e-6, i.e. a relative error of 1e-6, because the
# absolute term of the starting rule no longer tests anything there
# (tolerances.md section 5, "Small tail probabilities").

ref_reset_registry()

expect_prob_close <- function(ours, ref, info = NULL) {
  expect_true(ref_close(ours, ref, abs_tol = 1e-8, rel_tol = 1e-6),
              info = paste(info, sprintf("ours %.17g ref %.17g", ours, ref)))
}

expect_log_close <- function(ours_log, ref_log, info = NULL) {
  expect_true(is.finite(ours_log) && abs(ours_log - ref_log) <= 1e-6,
              info = paste(info, sprintf("log ours %.17g log ref %.17g", ours_log, ref_log)))
}

test_that("symmetric critical values follow the declared tail and side", {
  qupper <- function(p) qt(p, 15, lower.tail = FALSE)
  two <- .ig_critical_symmetric(0.05, "two", qupper = qupper)
  expect_identical(two, c(-qt(0.025, 15, lower.tail = FALSE), qt(0.025, 15, lower.tail = FALSE)))
  expect_identical(.ig_critical_symmetric(0.05, "one", "upper", qupper), qt(0.05, 15, lower.tail = FALSE))
  expect_identical(.ig_critical_symmetric(0.05, "one", "lower", qupper), -qt(0.05, 15, lower.tail = FALSE))
})

test_that("asymmetric critical values come from the tail they bound", {
  qtail <- function(p, lt) qchisq(p, 7, lower.tail = lt)
  expect_identical(.ig_critical_values(0.05, "two", qtail = qtail),
                   c(qchisq(0.025, 7), qchisq(0.025, 7, lower.tail = FALSE)))
  expect_identical(.ig_critical_values(0.05, "one", "upper", qtail), qchisq(0.05, 7, lower.tail = FALSE))
  expect_identical(.ig_critical_values(0.05, "one", "lower", qtail), qchisq(0.05, 7))
  # A tiny upper alpha is taken from the upper tail, not as qchisq(1 - alpha).
  up <- .ig_critical_values(1e-20, "one", "upper", qtail)
  expect_equal(pchisq(up, 7, lower.tail = FALSE), 1e-20, tolerance = 1e-10)
})

test_that("two-sided region probabilities agree with the references in both tails", {
  # Noncentral chi-square under H1 with H0 critical values of chi-square(5).
  crit <- .ig_critical_values(0.05, "two", qtail = function(p, lt) qchisq(p, 5, lower.tail = lt))
  for (ncp in c(0.5, 4, 30)) {
    reg <- .ig_two_sided_prob(function(q, lt) .ig_p_chisq(q, 5, ncp, lower_tail = lt), crit[1], crit[2])
    lo_ref <- ref_pchisq(crit[1], 5, ncp, lower_tail = TRUE)
    up_ref <- ref_pchisq(crit[2], 5, ncp, lower_tail = FALSE)
    acc_ref <- ref_pchisq(crit[2], 5, ncp, lower_tail = TRUE) - lo_ref
    info <- paste("chi-square ncp =", ncp)
    expect_prob_close(reg$lower, lo_ref, info)
    expect_prob_close(reg$upper, up_ref, info)
    expect_prob_close(reg$reject, lo_ref + up_ref, info)
    expect_prob_close(reg$accept, acc_ref, info)
    expect_equal(reg$reject + reg$accept, 1, tolerance = 1e-13)
  }
  # Noncentral F, two-sided with alpha/2 in each tail of F(3, 40).
  crit <- .ig_critical_values(0.1, "two", qtail = function(p, lt) qf(p, 3, 40, lower.tail = lt))
  reg <- .ig_two_sided_prob(function(q, lt) .ig_p_f(q, 3, 40, 6, lower_tail = lt), crit[1], crit[2])
  expect_prob_close(reg$reject, ref_pf(crit[1], 3, 40, 6, lower_tail = TRUE) +
                      ref_pf(crit[2], 3, 40, 6, lower_tail = FALSE), "F ncp = 6")
  # A tiny acceptance probability keeps its relative accuracy: beta of a
  # two-sided t test far from H0 (about 1e-17), against the reference
  # difference of two lower tails.
  tc <- qt(0.025, 30, lower.tail = FALSE)
  reg <- .ig_two_sided_prob(function(q, lt) .ig_p_t(q, 30, 11, lower_tail = lt), -tc, tc)
  ref_acc <- ref_pt(tc, 30, 11) - ref_pt(-tc, 30, 11)
  expect_lt(reg$accept, 1e-12)
  expect_log_close(log(reg$accept), log(ref_acc), "tiny beta")
  expect_error(.ig_two_sided_prob(function(q, lt) 0.5, 2, 1), "lo <= hi")
})

test_that("t-test power counts both tails when two-sided and follows the sign of ncp when one-sided", {
  df <- 40
  for (ncp in c(-2.5, -0.3, 0.3, 2.5)) {
    two <- .ig_power_t_test(df, ncp, 0.05, "two")
    tc <- qt(0.025, df, lower.tail = FALSE)
    expect_identical(two$critical, c(-tc, tc))
    upper <- ref_pt(tc, df, ncp, lower_tail = FALSE)
    lower <- ref_pt(-tc, df, ncp, lower_tail = TRUE)
    expect_prob_close(two$power, upper + lower, paste("two-sided ncp =", ncp))
    expect_equal(two$power + two$beta, 1, tolerance = 1e-14)
    one <- .ig_power_t_test(df, ncp, 0.05, "one")
    c1 <- qt(0.05, df, lower.tail = FALSE)
    expect_identical(one$critical, sign(ncp) * c1)
    expect_identical(one$side, if (ncp > 0) "upper" else "lower")
    ref_one <- if (ncp > 0) ref_pt(c1, df, ncp, lower_tail = FALSE) else ref_pt(-c1, df, ncp)
    expect_prob_close(one$power, ref_one, paste("one-sided ncp =", ncp))
  }
  # ncp = 0: the power equals alpha in either convention.
  expect_equal(.ig_power_t_test(12, 0, 0.05, "two")$power, 0.05, tolerance = 1e-14)
  expect_equal(.ig_power_t_test(12, 0, 0.05, "one")$power, 0.05, tolerance = 1e-14)
  # beta is computed directly, not as 1 - power: about 1e-20 here.
  b <- .ig_power_t_test(100, 12, 0.05, "two")$beta
  ref_b <- ref_pt(qt(0.025, 100, lower.tail = FALSE), 100, 12) - ref_pt(-qt(0.025, 100, lower.tail = FALSE), 100, 12)
  expect_lt(b, 1e-15)
  expect_log_close(log(b), log(ref_b))
})
