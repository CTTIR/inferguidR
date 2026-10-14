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

test_that("log-space helpers are accurate over the whole range", {
  expect_equal(.ig_log1mexp(-1e-20), log(1e-20), tolerance = 1e-14)
  expect_equal(.ig_log1mexp(-50), -exp(-50), tolerance = 1e-14)
  expect_identical(.ig_log1mexp(0), -Inf)
  expect_true(is.nan(.ig_log1mexp(0.1)))
  expect_true(is.na(.ig_log1mexp(NA_real_)))
  expect_equal(.ig_logspace_add(log(0.25), log(0.5)), log(0.75), tolerance = 1e-15)
  expect_equal(.ig_logspace_add(-1000, -1000), -1000 + log(2), tolerance = 1e-15)
  expect_identical(.ig_logspace_add(-Inf, -Inf), -Inf)
  expect_identical(.ig_logspace_sum(numeric(0)), -Inf)
  expect_equal(.ig_logspace_sum(c(-800, -800, -800)), -800 + log(3), tolerance = 1e-15)
})

test_that("the noncentral t agrees with the definition-based reference in every regime", {
  # Small and fractional df, the nmath series range, the regime above
  # df = 4e5 where nmath switches to a normal approximation, very large df;
  # ncp of both signs on both sides of nmath's |ncp| = 37.62 switch.
  grid <- expand.grid(q = c(-3, 0.7, 8), df = c(0.5, 2.5, 88, 4e5 + 1, 1e7),
                      ncp = c(-45, -5, 0.3, 37, 38, 120), lower = c(TRUE, FALSE))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    info <- sprintf("q=%g df=%g ncp=%g lower=%s", g$q, g$df, g$ncp, g$lower)
    ref_log <- ref_pt(g$q, g$df, g$ncp, lower_tail = g$lower, log_p = TRUE)
    ours_log <- .ig_p_t(g$q, g$df, g$ncp, lower_tail = g$lower, log_p = TRUE)
    ours <- .ig_p_t(g$q, g$df, g$ncp, lower_tail = g$lower)
    expect_equal(ours, exp(ours_log), tolerance = 1e-12, info = info)
    if (exp(ref_log) >= 1e-8) expect_prob_close(ours, exp(ref_log), info) else
      expect_log_close(ours_log, ref_log, info)
  }
})

test_that("the integral regime agrees with a reference that shares no quadrature routine", {
  # ref_pt() and the engine's integral route both call stats::integrate();
  # ref_pt_gl() uses its own Gauss-Legendre nodes. Points: |ncp| > 37.62,
  # df > 4e5 and tails below 1e-5, i.e. exactly where the engine integrates.
  pts <- list(c(45, 5, 40, 0), c(30, 5, 40, 1), c(-2, 12, -200, 0), c(3, 12, -200, 0),
              c(1.96, 4e5 + 1, 37, 0), c(2, 1e7, 0.3, 1), c(8, 30, 0.3, 0), c(-5, 30, 8, 1),
              c(1.67, 499998, 3.5, 1), c(200, 4, 1, 0))
  for (p in pts) {
    lt <- as.logical(p[4])
    info <- sprintf("q=%g df=%g ncp=%g lower=%s", p[1], p[2], p[3], lt)
    ref_log <- ref_pt_gl(p[1], p[2], p[3], lower_tail = lt, log_p = TRUE)
    ours_log <- .ig_p_t(p[1], p[2], p[3], lower_tail = lt, log_p = TRUE)
    if (exp(ref_log) >= 1e-8) expect_prob_close(exp(ours_log), exp(ref_log), info) else
      expect_log_close(ours_log, ref_log, info)
  }
  # The Gauss-Legendre reference itself against closed forms: the Cauchy
  # tail (df = 1) and the df = 2 tail 1/2 - t / (2 sqrt(2 + t^2)).
  expect_equal(ref_pt_gl(1 / tan(pi * 1e-50), 1, 0, lower_tail = FALSE), 1e-50, tolerance = 1e-10)
  expect_equal(ref_pt_gl(3, 2, 0, lower_tail = FALSE), 0.5 - 3 / (2 * sqrt(11)), tolerance = 1e-11)
})

test_that("extreme noncentrality is computed correctly where nmath approximates", {
  # nmath pnt uses the normal approximation A&S 26.7.10 for |ncp| > 37.62
  # (about 0.02 off at df = 5, ncp = 40); the helper integrates instead.
  for (q in c(30, 40, 45, 60)) {
    expect_prob_close(.ig_p_t(q, 5, 40), ref_pt(q, 5, 40), paste("q =", q))
  }
  expect_prob_close(.ig_p_t(-2, 12, -200, lower_tail = FALSE), ref_pt(-2, 12, -200, lower_tail = FALSE))
  expect_log_close(.ig_p_t(3, 12, -200, lower_tail = FALSE, log_p = TRUE),
                   ref_pt(3, 12, -200, lower_tail = FALSE, log_p = TRUE))
})

test_that("small tails are survival probabilities, not 1 - (nearly one)", {
  # Upper tail of about 1.2e-8, computed directly to ~1e-12 relative
  # accuracy; the complement 1 - P(T <= q) of the same point is less accurate.
  up_log <- .ig_p_t(200, 4, 1, lower_tail = FALSE, log_p = TRUE)
  ref_up_log <- ref_pt(200, 4, 1, lower_tail = FALSE, log_p = TRUE)
  ours_err <- abs(expm1(up_log - ref_up_log))
  naive_err <- abs((1 - .ig_p_t(200, 4, 1)) / exp(ref_up_log) - 1)
  expect_lt(ours_err, 1e-9)
  expect_gt(naive_err, ours_err)
  # A lower tail of about 2e-30, far below the resolution of 1 - p.
  lo_log <- .ig_p_t(-5, 30, 8, log_p = TRUE)
  expect_log_close(lo_log, ref_pt(-5, 30, 8, log_p = TRUE))
  expect_lt(lo_log, log(.Machine$double.eps) - 30)
  # The two tails of one point add up to one, also where both come from the
  # integral route (|ncp| > 37.62); the larger tail never exceeds one.
  for (args in list(c(1.3, 7, 0.4), c(-0.2, 30, -2), c(2, 1e6, 1.5), c(1.67, 58, 4.5e4),
                    c(1.67, 58, -4.5e4), c(-1.67, 58, 3e3))) {
    for (lt in c(TRUE, FALSE)) expect_lte(.ig_p_t(args[1], args[2], args[3], lower_tail = lt), 1)
    tot <- .ig_p_t(args[1], args[2], args[3]) + .ig_p_t(args[1], args[2], args[3], lower_tail = FALSE)
    expect_equal(tot, 1, tolerance = 1e-13)
  }
})

test_that("tails far below the double range are exact 0 or 1 by a rigorous bound", {
  # The bound is an upper bound on the log tail wherever it is finite.
  withr::local_seed(20261005)
  for (i in seq_len(150)) {
    q <- sample(c(-1, 1), 1) * 10^stats::runif(1, -2, 1.5)
    df <- 10^stats::runif(1, 0, 6)
    ncp <- sample(c(-1, 1), 1) * 10^stats::runif(1, -1, 1.8)
    lt <- stats::runif(1) < 0.5
    b <- .ig_t_tail_log_bound(q, df, ncp, lt)
    if (b < 0) {
      lp <- .ig_p_t(q, df, ncp, lower_tail = lt, log_p = TRUE)
      expect_lte(lp, b + 1e-12, label = sprintf("log tail q=%g df=%g ncp=%g lower=%s", q, df, ncp, lt))
    }
  }
  # Huge noncentrality: the probability-scale answers are exact and agree
  # with the log scale (which is still integrated).
  tc <- qt(0.025, 5e5 - 2, lower.tail = FALSE)
  expect_identical(.ig_p_t(tc, 5e5 - 2, 3.5e6), 0)
  expect_identical(.ig_p_t(tc, 5e5 - 2, 3.5e6, lower_tail = FALSE), 1)
  expect_identical(.ig_p_t(-tc, 5e5 - 2, -3.5e6, lower_tail = FALSE), 0)
  lp <- .ig_p_t(tc, 5e5 - 2, 3.5e6, log_p = TRUE)
  expect_lt(lp, -1e12)
  expect_identical(exp(lp), 0)
  # Where the bound is not decisive the result is the integral, unchanged.
  expect_gt(.ig_p_t(30, 5, 40), 0.1)
})

test_that("a quadrature roundoff message with a negligible error raises no warning", {
  # A narrow peak at w = 3e8 (the larger tail at ncp = 3.16e8, df = 4e9 - 2):
  # QUADPACK reports roundoff on pieces whose summed error is about 3e-12 of
  # the total. Such a message must not reach the user as a warning.
  tc <- qt(0.025, 4e9 - 2, lower.tail = FALSE)
  expect_silent(lp <- .ig_p_t_integral(tc, 4e9 - 2, 1e4 * sqrt(1e9), lower_tail = FALSE))
  expect_gt(lp, -1e-10)
})

test_that("special arguments of the t helper have exact answers", {
  for (ncp in c(0, 2.5, -50)) {
    expect_identical(.ig_p_t(Inf, 10, ncp), 1)
    expect_identical(.ig_p_t(-Inf, 10, ncp), 0)
    expect_identical(.ig_p_t(Inf, 10, ncp, lower_tail = FALSE), 0)
    expect_identical(.ig_p_t(Inf, 10, ncp, log_p = TRUE), 0)
    expect_identical(.ig_p_t(-Inf, 10, ncp, log_p = TRUE), -Inf)
  }
  # q = 0: P(T <= 0) = Phi(-ncp) for every df.
  expect_equal(.ig_p_t(0, 3, 1.7), pnorm(-1.7), tolerance = 1e-15)
  # df = Inf is the normal limit; ncp = 0 is the central t.
  expect_equal(.ig_p_t(1.2, Inf, 0.5), pnorm(1.2, mean = 0.5), tolerance = 1e-15)
  expect_identical(.ig_p_t(1.2, 7, 0), pt(1.2, 7))
  # Domain violations give NaN, missing values NA; arguments recycle.
  expect_true(is.nan(.ig_p_t(1, 0, 1)))
  expect_true(is.nan(.ig_p_t(1, -2, 0.5)))
  expect_true(is.nan(.ig_p_t(1, 5, Inf)))
  expect_true(is.na(.ig_p_t(NA, 5, 1)))
  expect_length(.ig_p_t(c(-1, 0, 1), 5, c(0.5, 1, 2)), 3L)
  expect_identical(.ig_p_t(numeric(0), 5, 1), numeric(0))
})

test_that("noncentral t quantiles invert the CDF, including tiny probabilities", {
  for (args in list(c(0.3, 12, 2), c(1e-20, 12, 2), c(0.9, 3, -4), c(1e-6, 4e5 + 10, 40))) {
    p <- args[1]
    x <- .ig_q_t(p, args[2], args[3])
    expect_equal(.ig_p_t(x, args[2], args[3]), p, tolerance = 1e-9)
    xu <- .ig_q_t(p, args[2], args[3], lower_tail = FALSE)
    expect_equal(.ig_p_t(xu, args[2], args[3], lower_tail = FALSE), p, tolerance = 1e-9)
  }
  x <- .ig_q_t(log(1e-200), 10, 3, log_p = TRUE)
  expect_equal(.ig_p_t(x, 10, 3, log_p = TRUE), log(1e-200), tolerance = 1e-9)
  expect_identical(.ig_q_t(0.975, 10), qt(0.975, 10))
  expect_identical(.ig_q_t(0, 10, 2), -Inf)
  expect_identical(.ig_q_t(1, 10, 2), Inf)
  # An impossible probability gives NaN with R's usual warning (as qt()).
  expect_warning(bad <- .ig_q_t(1.5, 10, 2), "NaN")
  expect_true(is.nan(bad))
})

test_that("the noncentral chi-square agrees with the convolution reference", {
  # Includes ncp >= 80, where nmath's pnchisq changes algorithm, and an upper
  # tail of 1e-23 that pchisq() returns as 0.
  pts <- list(c(3, 1, 2), c(10, 4, 7), c(150, 10, 85), c(200, 3, 120), c(0.5, 6, 30), c(400, 2, 100))
  for (p in pts) {
    for (lt in c(TRUE, FALSE)) {
      info <- sprintf("y=%g k=%g ncp=%g lower=%s", p[1], p[2], p[3], lt)
      ours <- .ig_p_chisq(p[1], p[2], p[3], lower_tail = lt)
      ref <- ref_pchisq(p[1], p[2], p[3], lower_tail = lt)
      if (ref >= 1e-8) expect_prob_close(ours, ref, info) else
        expect_log_close(log(ours), log(ref), info)
    }
  }
  expect_identical(.ig_p_chisq(-1, 3, 2), 0)
  expect_identical(.ig_p_chisq(Inf, 3, 2, lower_tail = FALSE), 0)
  expect_identical(.ig_p_chisq(0, 3, 2, lower_tail = FALSE, log_p = TRUE), 0)
  expect_true(is.nan(.ig_p_chisq(1, 3, -1)))
  q <- .ig_q_chisq(1e-12, 4, 50, lower_tail = FALSE)
  expect_equal(.ig_p_chisq(q, 4, 50, lower_tail = FALSE), 1e-12, tolerance = 1e-9)
})

test_that("the noncentral F agrees with the integral reference", {
  pts <- list(c(2, 3, 20, 5), c(0.8, 1, 10, 2), c(5, 4, 60, 30), c(1.5, 10, 200, 1))
  for (p in pts) {
    for (lt in c(TRUE, FALSE)) {
      expect_prob_close(.ig_p_f(p[1], p[2], p[3], p[4], lower_tail = lt),
                        ref_pf(p[1], p[2], p[3], p[4], lower_tail = lt),
                        sprintf("x=%g df1=%g df2=%g ncp=%g lower=%s", p[1], p[2], p[3], p[4], lt))
    }
  }
  # df2 = Inf is the scaled chi-square limit.
  expect_equal(.ig_p_f(2, 3, Inf, 4), .ig_p_chisq(6, 3, 4), tolerance = 1e-15)
  expect_identical(.ig_p_f(0, 3, 20, 5), 0)
  expect_identical(.ig_p_f(Inf, 3, 20, 5), 1)
  x <- .ig_q_f(0.95, 3, 20, 5)
  expect_equal(.ig_p_f(x, 3, 20, 5), 0.95, tolerance = 1e-10)
})

test_that("binomial tails agree with an lchoose summation in both directions", {
  pts <- list(c(3, 20, 0.3), c(15, 20, 0.3), c(0, 50, 0.01), c(45, 50, 0.99), c(900, 1000, 0.5))
  for (p in pts) {
    for (lt in c(TRUE, FALSE)) {
      ours <- .ig_p_binom(p[1], p[2], p[3], lower_tail = lt)
      ref <- ref_pbinom(p[1], p[2], p[3], lower_tail = lt)
      info <- sprintf("q=%g n=%g p=%g lower=%s", p[1], p[2], p[3], lt)
      if (ref >= 1e-8) expect_prob_close(ours, ref, info) else expect_log_close(log(ours), log(ref), info)
    }
  }
  expect_equal(sum(.ig_d_binom(0:12, 12, 0.37)), 1, tolerance = 1e-14)
})

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

test_that("noncentral F uses a checked small tail and stable complement", {
  expect_lte(.ig_p_f(1, 2, 100, 1821.4, lower_tail = FALSE), 1)
  for (q in c(1e-14, 0.1, 1, 100, 1e14)) {
    lo <- .ig_p_f(q, 2, 100, 1821.4, log_p = TRUE)
    hi <- .ig_p_f(q, 2, 100, 1821.4, lower_tail = FALSE, log_p = TRUE)
    expect_lte(lo, 0)
    expect_lte(hi, 0)
    expect_equal(exp(lo) + exp(hi), 1, tolerance = 1e-14)
  }
})

test_that("proven noncentral F underflow avoids enormous mixture allocations", {
  local_mocked_bindings(.ig_poisson_mixture_log = function(...) stop("mixture must not be called"))
  expect_identical(.ig_p_f(1, 2, 100, 1e16), 0)
  expect_identical(.ig_p_f(1, 2, 100, 1e16, lower_tail = FALSE), 1)
  expect_identical(.ig_p_f(1, 2, 100, 1e16, lower_tail = FALSE, log_p = TRUE), 0)
})

test_that("Poisson mixture term limit is enforced before constructing a window", {
  expect_true(is.nan(.ig_poisson_mixture_log(1e8, function(j) stop("oversized window evaluated"), max_terms = 10)))
  lengths <- integer()
  result <- .ig_poisson_mixture_log(10, function(j) { lengths <<- c(lengths, length(j)); j }, max_terms = 150)
  expect_true(all(lengths <= 150))
  expect_true(is.finite(result) || is.nan(result))
})

test_that("F underflow bounds handle very large denominator degrees of freedom", {
  expect_identical(.ig_p_f(1, 2, 1e30, 1821.4), 0)
  expect_lte(.ig_f_lower_log_bound(1, 0.5, 1e30, 1e10), -745.2)
  # Independent 80-digit gamma-mixture reference with a finite-denominator
  # Chernoff bracket; artifact: shared-f-tails-20261007/log-reference.json.
  expect_lt(abs(.ig_p_f(1, 2, 1e30, 1821.4, log_p = TRUE) +
                  857.69399271029176), 1e-9)
  lp <- .ig_p_f(1, 2, 100, 1821.4, log_p = TRUE)
  expect_true(is.finite(lp))
  expect_lt(lp, 0)
})
