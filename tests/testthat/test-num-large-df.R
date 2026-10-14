# Regression cases from shared integration requests CR-G1-01/02 and G2/G4/G11.
# The tighter t tolerance protects inverse-mode rounding at large N; the
# reference integrates over the denominator with independent quadrature.
test_that("large-df noncentral t tails retain inverse-solver precision", {
  for (df in c(30000, 100000, 199998)) {
    for (lt in c(TRUE, FALSE)) {
      want <- ref_pt_gl(3.09, df, 4.47, lower_tail = lt)
      got <- .ig_p_t(3.09, df, 4.47, lower_tail = lt)
      expect_true(abs(got - want) < 2e-13, info = paste(df, lt))
    }
  }
})

test_that("saturated noncentral t tails remain probabilities", {
  for (df in c(19000, 199998)) {
    for (sgn in c(-1, 1)) {
      lo <- .ig_p_t(sgn * 3.09, df, sgn * 22.36)
      up <- .ig_p_t(sgn * 3.09, df, sgn * 22.36, lower_tail = FALSE)
      expect_true(all(is.finite(c(lo, up))) && all(c(lo, up) >= 0 & c(lo, up) <= 1))
      expect_equal(lo + up, 1, tolerance = 2e-15)
      expect_equal(.ig_p_t(sgn * 3.09, df, sgn * 22.36, log_p = TRUE),
                   ref_pt_gl(sgn * 3.09, df, sgn * 22.36, log_p = TRUE), tolerance = 1e-10)
    }
  }
})

# Central F quantiles must preserve the requested significance level even
# when either df crosses nmath qf's chi-square approximation threshold.
test_that("central F quantiles retain their tail area at large df", {
  for (ds in list(c(50, 1444522), c(1647266, 164726), c(380, 1381548.27), c(3, 400001))) {
    for (lt in c(TRUE, FALSE)) {
      for (p in c(1e-8, 0.001, 0.05, 0.9)) {
        q <- .ig_q_f(p, ds[1], ds[2], lower_tail = lt)
        lp <- stats::pf(q, ds[1], ds[2], lower.tail = lt, log.p = TRUE)
        expect_true(is.finite(q) && abs(lp - log(p)) < 1e-9,
                    info = paste(ds, p, lt, collapse = " "))
        expect_equal(.ig_q_f(log(p), ds[1], ds[2], lower_tail = lt, log_p = TRUE),
                     q, tolerance = 1e-12)
      }
    }
  }
})

test_that("central F quantiles preserve endpoints and limiting distributions", {
  expect_equal(.ig_q_f(c(0, 1), 3, 400001), c(0, Inf))
  expect_equal(.ig_q_f(c(-Inf, 0), 3, 400001, lower_tail = FALSE, log_p = TRUE), c(Inf, 0))
  expect_equal(.ig_q_f(0.001, 3, Inf, lower_tail = FALSE), qchisq(0.001, 3, lower.tail = FALSE) / 3)
  expect_equal(.ig_q_f(0.001, Inf, 5, lower_tail = FALSE), 5 / qchisq(0.001, 5))
  expect_true(is.na(.ig_q_f(NA_real_, 3, 400001)))
  expect_true(is.nan(suppressWarnings(.ig_q_f(-0.1, 3, 400001))))
})
