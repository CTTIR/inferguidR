# Vertical slice: t_means_two_independent ("Means: Difference between two
# independent means (two groups)") through all five modes, one shared result
# object, the effect-size drawer, sweeps and comparisons.
#
# Evidence kinds used here, kept apart as VALIDATION_CONTRACT requires:
# * independent checks: the definition-based reference of helper-reference.R
#   (noncentral t by integrating pnorm over the chi-square density; central t
#   quantiles by root finding on that integral), which shares no numerical
#   routine with the engine (decision D18);
# * documented examples from the G*Power 3.1 manual (validation evidence kind
#   `documented_example`; printed documentation, not application captures,
#   and not G*Power validation);
# * property and inverse/forward tests.
#
# Tolerances: probabilities use the starting rule abs 1e-8 + rel 1e-6
# (validation/tolerances.md section 5); displayed manual values use the
# closed rounding interval of the last printed digit plus the default
# floating-point allowance (tolerances.md section 4, rule MR-DISPLAY);
# integers are compared exactly (MR-EXACT).

ref_reset_registry()

slice_id <- "t_means_two_independent"

slice_spec <- function(mode, tails = "two", ...) {
  analysis_spec(slice_id, mode, inputs = list(...), options = list(tails = tails))
}

slice_post_hoc <- function(d, n1, n2, alpha = 0.05, tails = "two") {
  calculate_power(slice_spec("post_hoc", tails, alpha = alpha, n1 = n1, n2 = n2, d = d))
}

slice_numbers <- function(res) {
  x <- unclass(res)
  .ig_finite_paths(x[setdiff(names(x), c("spec", "diagnostics"))], "")
}

# MR-DISPLAY: closed rounding interval of the printed value plus the default
# allowance 8 * eps * max(1, |x|) (tolerances.md section 4).
slice_display_ok <- function(ours, shown) {
  x <- as.numeric(shown)
  dec <- if (grepl(".", shown, fixed = TRUE)) nchar(sub("^[^.]*\\.", "", shown)) else 0
  half <- 0.5 * 10^-dec
  allowance <- 8 * .Machine$double.eps * max(1, abs(x))
  is.finite(ours) && ours >= x - half - allowance && ours <= x + half + allowance
}

slice_expect_display <- function(ours, shown, info) {
  expect_true(slice_display_ok(ours, shown),
              info = sprintf("%s: printed %s, ours %.10g", info, shown, ours))
}

slice_expect_prob <- function(ours, ref, info = NULL) {
  expect_true(ref_close(ours, ref, abs_tol = 1e-8, rel_tol = 1e-6),
              info = paste(info, sprintf("ours %.15g ref %.15g", ours, ref)))
}

# Registration and honest status ----------------------------------------------------

test_that("the slice method offers all five modes with an honest validation state", {
  lm <- list_methods(family = "t")
  row <- lm[lm$id == slice_id, ]
  expect_identical(row$implementation_status, "implemented")
  expect_identical(row$validation_status, "independently_checked")
  expect_identical(row$modes, "a_priori, post_hoc, sensitivity, criterion, compromise")
  vs <- validation_status(slice_id)
  expect_identical(vs$gpower_reference$status, "reference_pending")
  expect_identical(vs$gpower_reference$fixture_ids, character())
  expect_identical(vs$independent_check$status, "passed")
  expect_identical(vs$independent_check$checker, "validation/independent/t/check_t_means_two_independent.R")
  expect_identical(vs$modes_with_solver, .ig_modes)
  expect_true(all(c("DOC-001", "DOC-003", "DOC-098") %in% vs$documented_example_ids))
  expect_true(any(grepl("agreement with G\\*Power is not claimed", vs$explanation)))
  expect_false(any(grepl("gpower_validated", unlist(lm))))
  res <- slice_post_hoc(0.5, 20, 20)
  expect_identical(res$validation$method_validation_status, "independently_checked")
  expect_identical(res$validation$gpower_reference_status, "reference_pending")
})

test_that("specifications are validated against the method definition", {
  code <- function(expr) tryCatch({
    force(expr)
    NA_character_
  }, inferguidR_error = function(e) e$status$code)
  # The tail option defaults to G*Power's live initial value "One"; the
  # default is written into the specification and listed as applied.
  dflt <- analysis_spec(slice_id, "post_hoc", inputs = list(alpha = 0.05, n1 = 10, n2 = 10, d = 0.5))
  expect_identical(dflt$options$tails, "one")
  expect_identical(dflt$defaults_applied$options, "tails")
  # A target power below alpha is refused (G*Power: 'Power must be greater
  # than or equal to alpha.'); the power is at least alpha for every d.
  cond <- tryCatch(slice_spec("a_priori", alpha = 0.05, power = 0.04, d = 0.5),
                   inferguidR_invalid_input = function(e) e)
  expect_match(cond$status$reason, "below alpha")
  # Group sizes are inputs outside a priori mode; N is never an input.
  expect_identical(code(slice_spec("post_hoc", alpha = 0.05, n_total = 20, d = 0.5)), "invalid_input")
  expect_identical(code(slice_spec("a_priori", alpha = 0.05, power = 0.8, d = 0.5, n1 = 10)),
                   "invalid_input")
  expect_identical(code(slice_spec("post_hoc", alpha = 0.05, n1 = 10, n2 = 10, d = 0.5,
                                   allocation_ratio = 2)), "invalid_input")
  expect_identical(code(slice_spec("post_hoc", alpha = 0.05, n1 = 10.5, n2 = 10, d = 0.5)),
                   "invalid_input")
  # df = n1 + n2 - 2 must be at least 1.
  cond <- tryCatch(slice_spec("post_hoc", alpha = 0.05, n1 = 1, n2 = 1, d = 0.5),
                   inferguidR_invalid_input = function(e) e)
  expect_match(cond$status$reason, "n1 \\+ n2 must be at least 3")
  # The allocation ratio default is written into the specification.
  ap <- slice_spec("a_priori", alpha = 0.05, power = 0.8, d = 0.5)
  expect_identical(ap$inputs$allocation_ratio, 1)
  expect_identical(ap$defaults_applied$inputs, "allocation_ratio")
  # The G*Power-style alias of the effect size is accepted.
  al <- analysis_spec(slice_id, "post_hoc", inputs = list(alpha = 0.05, sample_size_group_1 = 10,
                                                          n2 = 10, effect_size_d = 0.5),
                      options = list(tails = "two"))
  expect_identical(al$inputs$d, 0.5)
  expect_identical(al$inputs$n1, 10L)
})

# Post hoc: forward model ---------------------------------------------------------------

test_that("post hoc power, df, ncp and critical values match the independent reference", {
  cases <- list(
    list(d = 0.5, n1 = 30, n2 = 60, tails = "two"),
    list(d = 0.5, n1 = 30, n2 = 60, tails = "one"),
    list(d = -0.8, n1 = 12, n2 = 7, tails = "one"),
    list(d = 0.1, n1 = 10, n2 = 10, tails = "two"),
    list(d = 1.2, n1 = 1, n2 = 2, tails = "two"),
    list(d = 0.3, n1 = 400, n2 = 250, tails = "one"),
    list(d = 2.5, n1 = 40, n2 = 40, tails = "two"),
    # n1 * n2 exceeds the integer range; df > 4e5 uses the integral route.
    list(d = 0.01, n1 = 250000, n2 = 250000, tails = "two")
  )
  for (cs in cases) {
    info <- sprintf("d=%g n1=%g n2=%g tails=%s", cs$d, cs$n1, cs$n2, cs$tails)
    res <- slice_post_hoc(cs$d, cs$n1, cs$n2, tails = cs$tails)
    ref <- ref_power_t2(cs$d, cs$n1, cs$n2, 0.05, cs$tails)
    expect_identical(res$status$code, "ok", info = info)
    slice_expect_prob(res$solved$power, ref$power, info)
    expect_identical(res$achieved$power, res$solved$power)
    expect_equal(res$achieved$power + res$achieved$beta, 1, tolerance = 1e-14, info = info)
    expect_identical(res$quantities$df, cs$n1 + cs$n2 - 2)
    expect_equal(res$quantities$ncp, cs$d * sqrt(cs$n1 * cs$n2 / (cs$n1 + cs$n2)), tolerance = 1e-15)
    # Quantile tolerance from the CDF tolerance and the density (tolerances.md 5).
    expect_equal(res$quantities$critical, ref$critical, tolerance = 1e-8, info = info)
    expect_identical(res$group_sizes, c(n1 = as.integer(cs$n1), n2 = as.integer(cs$n2)))
    expect_identical(res$quantities$n_total, cs$n1 + cs$n2)
    expect_identical(res$quantities$alpha_nominal, 0.05)
    expect_identical(res$achieved$alpha_attained, 0.05)
  }
})

test_that("two-sided power counts the opposite tail; one-sided tests reject on the side of d", {
  res <- slice_post_hoc(0.1, 10, 10, tails = "two")
  tc <- res$quantities$critical[2]
  ncp <- res$quantities$ncp
  upper <- ref_pt(tc, 18, ncp, lower_tail = FALSE)
  lower <- ref_pt(-tc, 18, ncp)
  expect_gt(lower, 0.01) # the opposite tail is material here
  slice_expect_prob(res$solved$power, upper + lower)
  expect_identical(res$method_outputs$statistic_distribution$rejection_side, "both")

  neg <- slice_post_hoc(-0.6, 20, 20, tails = "one")
  expect_lt(neg$quantities$critical, 0)
  expect_identical(neg$method_outputs$statistic_distribution$rejection_side, "lower")
  slice_expect_prob(neg$solved$power, ref_pt(neg$quantities$critical, 38, neg$quantities$ncp))
  expect_gt(neg$solved$power, 0.5)

  for (tails in c("one", "two")) {
    zero <- slice_post_hoc(0, 15, 25, tails = tails)
    expect_equal(zero$solved$power, 0.05, tolerance = 1e-13, info = tails)
  }
  # A one-sided alpha above one half puts the critical value on the other
  # side of zero (rejection still on the side of d).
  for (d in c(0.4, -0.4)) {
    big_alpha <- slice_post_hoc(d, 15, 9, alpha = 0.7, tails = "one")
    ref <- ref_power_t2(d, 15, 9, 0.7, "one")
    expect_equal(big_alpha$quantities$critical, ref$critical, tolerance = 1e-8)
    expect_identical(sign(big_alpha$quantities$critical), -sign(d))
    slice_expect_prob(big_alpha$solved$power, ref$power, paste("alpha 0.7, d =", d))
  }
  # Two-sided power at alpha = one-sided power at alpha/2 plus the opposite tail.
  two <- slice_post_hoc(0.4, 25, 25, alpha = 0.05, tails = "two")
  one <- slice_post_hoc(0.4, 25, 25, alpha = 0.025, tails = "one")
  expect_identical(one$quantities$critical, two$quantities$critical[2])
  slice_expect_prob(two$solved$power - one$solved$power, ref_pt(-two$quantities$critical[2], 48, two$quantities$ncp))
})

test_that("sign reversal with group exchange leaves the power unchanged", {
  for (tails in c("one", "two")) {
    a <- slice_post_hoc(0.45, 18, 31, tails = tails)
    b <- slice_post_hoc(-0.45, 31, 18, tails = tails)
    expect_equal(b$solved$power, a$solved$power, tolerance = 1e-14, info = tails)
    expect_equal(b$achieved$beta, a$achieved$beta, tolerance = 1e-13, info = tails)
    expect_identical(b$quantities$ncp, -a$quantities$ncp)
    expect_identical(b$group_sizes, c(n1 = 31L, n2 = 18L))
    if (tails == "one") {
      expect_identical(b$quantities$critical, -a$quantities$critical)
    } else {
      expect_identical(b$quantities$critical, a$quantities$critical)
    }
    # The noncentrality is symmetric in n1 and n2, so relabelling the groups
    # without changing d leaves the result unchanged as well.
    c2 <- slice_post_hoc(0.45, 31, 18, tails = tails)
    expect_identical(c2$solved$power, a$solved$power)
  }
})

# A priori ---------------------------------------------------------------------------------

test_that("a priori: smallest even N, predecessor fails by the independent reference", {
  for (cs in list(list(d = 0.5, power = 0.95, tails = "one"), list(d = 2, power = 0.95, tails = "one"),
                  list(d = 0.8, power = 0.9, tails = "two"))) {
    info <- sprintf("d=%g power=%g tails=%s", cs$d, cs$power, cs$tails)
    res <- solve_design(slice_spec("a_priori", cs$tails, alpha = 0.05, power = cs$power, d = cs$d))
    N <- res$solved$n_total
    expect_identical(res$status$code, "ok", info = info)
    expect_type(N, "integer")
    expect_identical(N %% 2L, 0L)
    expect_identical(res$group_sizes, c(n1 = N %/% 2L, n2 = N %/% 2L))
    expect_gte(res$achieved$power, cs$power)
    # Independent evidence of minimality: the reference power at N and N - 2.
    expect_gte(ref_power_t2(cs$d, N / 2, N / 2, 0.05, cs$tails)$power, cs$power)
    expect_lt(ref_power_t2(cs$d, N / 2 - 1, N / 2 - 1, 0.05, cs$tails)$power, cs$power)
    dg <- res$diagnostics
    expect_identical(dg$predecessor_n_total, N - 2)
    expect_lt(dg$predecessor_power, cs$power)
    expect_match(dg$minimality, "^proven")
    expect_identical(dg$enumerated_range, c(4, N))
    expect_false(dg$nonmonotone_detected)
    expect_match(dg$allocation_rule, "n1 = n2 = N/2")
    expect_identical(dg$target_power, cs$power)
    # The continuous relaxation lies below the integer solution.
    expect_lte(dg$continuous_relaxation$n_total, N)
    expect_gt(dg$continuous_relaxation$n_total, N - 2)
  }
})

test_that("a priori with a large solution: bisection with predecessor and window checks", {
  res <- solve_design(slice_spec("a_priori", "two", alpha = 0.05, power = 0.9, d = 0.05))
  N <- res$solved$n_total
  expect_identical(res$diagnostics$search, "bisection")
  expect_match(res$diagnostics$minimality, "predecessor fails")
  expect_identical(N, as.integer(2 * ceiling(res$diagnostics$continuous_relaxation$n_total / 2)))
  pred <- slice_post_hoc(0.05, N / 2 - 1, N / 2 - 1)
  expect_lt(pred$solved$power, 0.9)
  expect_gte(slice_post_hoc(0.05, N / 2, N / 2)$solved$power, 0.9)
  expect_lt(res$diagnostics$evaluations, 150)
})

test_that("allocation rule: even N, n1 = floor(N/(1 + r) + 1/2), no empty group", {
  p <- list(allocation_ratio = 1)
  expect_null(.ig_design_t_means_two_independent(p, 7))
  expect_identical(.ig_design_t_means_two_independent(p, 8), c(n1 = 4, n2 = 4))
  expect_identical(.ig_design_t_means_two_independent(list(allocation_ratio = 2), 9, continuous = TRUE),
                   c(n1 = 3, n2 = 6))
  # Round half up at an exact tie: 134 / 4 = 33.5 -> 34, 170 / 4 = 42.5 -> 43.
  expect_identical(.ig_design_t_means_two_independent(list(allocation_ratio = 3), 134), c(n1 = 34, n2 = 100))
  expect_identical(.ig_design_t_means_two_independent(list(allocation_ratio = 3), 170), c(n1 = 43, n2 = 127))
  expect_identical(.ig_design_t_means_two_independent(list(allocation_ratio = 0.1), 887520),
                   c(n1 = 806836, n2 = 80684))
  expect_null(.ig_design_t_means_two_independent(list(allocation_ratio = 3), 135))
  # An empty group is not an admissible design (floor(4/11 + 1/2) = 0).
  expect_null(.ig_design_t_means_two_independent(list(allocation_ratio = 10), 4))
  expect_identical(.ig_design_t_means_two_independent(list(allocation_ratio = 10), 6), c(n1 = 1, n2 = 5))
  # Both group sizes are non-decreasing in N, so the rounded-split power is
  # non-decreasing over even N (used by the a priori search).
  for (r in c(0.1, 0.37, 2.5, 10)) {
    gs <- vapply(seq(4, 400, by = 2), function(N) {
      g <- .ig_design_t_means_two_independent(list(allocation_ratio = r), N)
      if (is.null(g)) c(0, 0) else unname(g)
    }, numeric(2))
    expect_true(all(diff(gs[1, ]) >= 0) && all(diff(gs[2, ]) >= 0), info = r)
  }
  # Unequal groups are fully supported when the group sizes are given.
  uneq <- slice_post_hoc(0.5, 30, 60)
  expect_identical(uneq$status$code, "ok")
  slice_expect_prob(uneq$solved$power, ref_power_t2(0.5, 30, 60, 0.05, "two")$power)
})

test_that("a priori with unequal allocation follows G*Power's fractional criterion", {
  # N* is the smallest even N whose power with fractional group sizes
  # N/(1 + r), rN/(1 + r) reaches the target; checked here with the
  # independent reference (ref_power_t2 accepts fractional group sizes).
  for (cs in list(list(d = 0.5, r = 3, tails = "one", power = 0.8), list(d = 0.8, r = 0.2, tails = "one", power = 0.8),
                  list(d = 1.5, r = 2.5, tails = "two", power = 0.95))) {
    info <- sprintf("d=%g r=%g tails=%s", cs$d, cs$r, cs$tails)
    res <- solve_design(slice_spec("a_priori", cs$tails, alpha = 0.05, power = cs$power, d = cs$d,
                                   allocation_ratio = cs$r))
    expect_identical(res$status$code, "ok", info = info)
    N <- res$solved$n_total
    expect_identical(N %% 2L, 0L)
    frac <- function(M) ref_power_t2(cs$d, M / (1 + cs$r), cs$r * M / (1 + cs$r), 0.05, cs$tails)$power
    expect_gte(frac(N), cs$power)
    expect_lt(frac(N - 2), cs$power)
    n1 <- floor(N / (1 + cs$r) + 0.5)
    expect_identical(res$group_sizes, c(n1 = as.integer(n1), n2 = as.integer(N - n1)), info = info)
    expect_gte(res$achieved$power, cs$power)
    expect_false(isTRUE(res$diagnostics$adjusted_from_gpower_rule))
    expect_true(res$diagnostics$gpower_allocation_rule$reaches_target)
    expect_match(res$diagnostics$allocation_rule, "floor\\(N/\\(1 \\+ r\\) \\+ 1/2\\)")
  }
  # The fractional criterion can choose a larger N than the integer one; the
  # integer-criterion N is reported for transparency (r = 0.2: 70 vs 72; the
  # G*Power capture of this request is SUP-RSC-T-MEAN2-APR-A03--r01).
  res <- solve_design(slice_spec("a_priori", "one", alpha = 0.05, power = 0.8, d = 0.8, allocation_ratio = 0.2))
  expect_identical(res$solved$n_total, 72L)
  expect_identical(res$diagnostics$integer_criterion_n_total, 70)
})

test_that("unequal allocation: designs that miss the target or leave a group empty are not reported", {
  # r = 10, d = 1.5, one tail (capture GPW-T-MEAN2-APR-B001-E01--r01): the
  # rule as computed by the engine gives N = 36 (3/33) whose
  # integer-design power is below 0.8; inferguidR reports the next even N
  # whose rounded split reaches the target and records the G*Power design.
  res <- solve_design(slice_spec("a_priori", "one", alpha = 0.05, power = 0.8, d = 1.5, allocation_ratio = 10))
  expect_identical(res$status$code, "ok")
  N <- res$solved$n_total
  gp <- res$diagnostics$gpower_allocation_rule
  expect_identical(gp$n_total, 36)
  expect_identical(c(gp$n1, gp$n2), c(3, 33))
  expect_false(gp$reaches_target)
  expect_lt(gp$integer_design_power, 0.8)
  expect_lt(ref_power_t2(1.5, 3, 33, 0.05, "one")$power, 0.8)
  expect_true(res$diagnostics$adjusted_from_gpower_rule)
  expect_match(res$diagnostics$adjustment_note, "gives N = 36 \\(n1 = 3, n2 = 33\\)")
  expect_gt(N, 36L)
  expect_gte(res$achieved$power, 0.8)
  expect_gte(ref_power_t2(1.5, res$group_sizes[["n1"]], res$group_sizes[["n2"]], 0.05, "one")$power, 0.8)
  pred <- .ig_design_t_means_two_independent(list(allocation_ratio = 10), N - 2)
  expect_lt(ref_power_t2(1.5, pred[["n1"]], pred[["n2"]], 0.05, "one")$power, 0.8)
  # r = 10, d = 10 (capture SUP-RSC-T-MEAN2-APR-A09--r01): N = 4 leaves group
  # 1 empty.
  res <- solve_design(slice_spec("a_priori", "one", alpha = 0.05, power = 0.8, d = 10, allocation_ratio = 10))
  expect_identical(res$status$code, "ok")
  expect_true(res$diagnostics$gpower_allocation_rule$empty_group)
  expect_identical(res$solved$n_total, 6L)
  expect_identical(res$group_sizes, c(n1 = 1L, n2 = 5L))
  expect_match(res$diagnostics$adjustment_note, "leaves group 1 empty")
  # An extreme ratio whose groups stay empty up to the integer range is
  # infeasible without numbers.
  inf <- solve_design(slice_spec("a_priori", "one", alpha = 0.05, power = 0.8, d = 1, allocation_ratio = 1e12))
  expect_identical(inf$status$code, "infeasible")
  expect_identical(slice_numbers(inf), character())
  # No effect: the power does not change with N.
  flat <- solve_design(slice_spec("a_priori", "two", alpha = 0.05, power = 0.8, d = 0, allocation_ratio = 2))
  expect_identical(flat$status$code, "infeasible")
  expect_match(flat$status$reason, "power does not change with N")
})

test_that("unreachable targets are infeasible and report no number", {
  flat <- solve_design(slice_spec("a_priori", "one", alpha = 0.05, power = 0.8, d = 0))
  expect_identical(flat$status$code, "infeasible")
  expect_match(flat$status$reason, "power does not change with N")
  expect_identical(slice_numbers(flat), character())
  tiny <- solve_design(slice_spec("a_priori", "two", alpha = 0.05, power = 0.8, d = 1e-8))
  expect_identical(tiny$status$code, "infeasible")
  expect_match(tiny$status$reason, "not reached by any admissible design with N <= 2147483646")
  expect_identical(slice_numbers(tiny), character())
  # A target power below alpha is refused at specification time.
  low <- tryCatch(slice_spec("sensitivity", "two", alpha = 0.05, power = 0.03, n1 = 10, n2 = 10),
                  inferguidR_invalid_input = function(e) e)
  expect_s3_class(low, "inferguidR_invalid_input")
  expect_match(low$status$reason, "target power 0.03 is below alpha 0.05")
})

# Inverse modes -----------------------------------------------------------------------------

test_that("every mode inverts the forward model, for both tail options", {
  for (tails in c("one", "two")) {
    fwd <- slice_post_hoc(0.6, 25, 35, alpha = 0.05, tails = tails)
    P <- fwd$solved$power

    se <- calculate_power(slice_spec("sensitivity", tails, alpha = 0.05, power = P, n1 = 25, n2 = 35))
    expect_identical(se$status$code, "ok")
    expect_equal(se$solved$d, 0.6, tolerance = 1e-8, info = tails)
    expect_lt(abs(se$diagnostics$residual), 1e-9)
    expect_equal(se$diagnostics$other_branch_roots, -se$solved$d, tolerance = 1e-8)
    expect_equal(se$achieved$power, P, tolerance = 1e-9)

    cr <- calculate_power(slice_spec("criterion", tails, power = P, n1 = 25, n2 = 35, d = 0.6))
    expect_identical(cr$status$code, "ok")
    expect_equal(cr$solved$alpha, 0.05, tolerance = 1e-8, info = tails)
    expect_identical(cr$achieved$alpha_attained, cr$solved$alpha)
    tc <- qt(if (tails == "two") cr$solved$alpha / 2 else cr$solved$alpha, 58, lower.tail = FALSE)
    expect_equal(abs(cr$quantities$critical), rep(tc, if (tails == "two") 2 else 1), tolerance = 1e-12)

    q <- fwd$achieved$beta / 0.05
    co <- calculate_power(slice_spec("compromise", tails, error_ratio = q, n1 = 25, n2 = 35, d = 0.6))
    expect_identical(co$status$code, "ok")
    expect_equal(co$solved$alpha, 0.05, tolerance = 1e-8, info = tails)
    expect_equal(co$solved$power, P, tolerance = 1e-8)
    expect_equal(co$achieved$error_ratio, q, tolerance = 1e-9)

    eq <- slice_post_hoc(0.6, 30, 30, alpha = 0.05, tails = tails)
    ap <- solve_design(slice_spec("a_priori", tails, alpha = 0.05, power = eq$solved$power, d = 0.6))
    expect_identical(ap$solved$n_total, 60L, info = tails)
    expect_identical(ap$achieved$power, eq$solved$power)
  }
})

test_that("compromise: beta = q * alpha with nominal and attained alpha equal (continuous test)", {
  for (q in c(1, 4, 0.25)) {
    co <- calculate_power(slice_spec("compromise", "two", error_ratio = q, n1 = 30, n2 = 60, d = 0.5))
    a <- co$solved$alpha
    expect_equal(co$achieved$beta, q * a, tolerance = 1e-10, info = q)
    expect_identical(co$achieved$alpha_attained, a)
    expect_identical(co$quantities$alpha_nominal, a)
    expect_equal(co$quantities$critical[2], qt(a / 2, 88, lower.tail = FALSE), tolerance = 1e-12)
    check <- slice_post_hoc(0.5, 30, 60, alpha = a)
    expect_equal(check$solved$power, co$solved$power, tolerance = 1e-14)
    slice_expect_prob(co$solved$power, ref_power_t2(0.5, 30, 60, a, "two")$power, paste("q =", q))
  }
})

test_that("extreme but admissible requests give a checked answer or an honest typed failure", {
  # Compromise with beta/alpha = 1e-6: alpha close to one, where the logit
  # scale keeps the objective well conditioned.
  co <- calculate_power(slice_spec("compromise", "two", error_ratio = 1e-6, n1 = 20, n2 = 20, d = 0.5))
  expect_identical(co$status$code, "ok")
  expect_gt(co$solved$alpha, 0.99999)
  expect_equal(co$achieved$beta, 1e-6 * co$solved$alpha, tolerance = 1e-8)
  slice_expect_prob(co$solved$power, ref_power_t2(0.5, 20, 20, co$solved$alpha, "two")$power, "q = 1e-6")
  # Compromise with q = 1e-9 needs a beta of about 1e-9 formed from tail
  # probabilities near one half; whether double precision resolves it may
  # depend on the platform's math library. Either outcome must be honest.
  co9 <- calculate_power(slice_spec("compromise", "two", error_ratio = 1e-9, n1 = 20, n2 = 20, d = 0.5))
  expect_true(co9$status$code %in% c("ok", "numerical_failure"))
  if (identical(co9$status$code, "ok")) {
    expect_equal(co9$achieved$beta, 1e-9 * co9$solved$alpha, tolerance = 1e-6)
  } else {
    expect_match(co9$status$reason, "without a continuous trend")
    expect_identical(slice_numbers(co9), character())
  }
  # Criterion with a tiny solution alpha, and with alpha close to one.
  cr <- calculate_power(slice_spec("criterion", "one", power = 0.8, n1 = 500, n2 = 500, d = 2))
  expect_identical(cr$status$code, "ok")
  expect_lt(cr$solved$alpha, 1e-100)
  slice_expect_prob(cr$achieved$power, ref_power_t2(2, 500, 500, cr$solved$alpha, "one")$power, "tiny alpha")
  cr <- calculate_power(slice_spec("criterion", "two", power = 0.9999999, n1 = 20, n2 = 20, d = 0.01))
  expect_identical(cr$status$code, "ok")
  expect_gt(cr$solved$alpha, 0.999999)
  # A target power equal to alpha is the minimum of the power curve (d = 0):
  # touched, not crossed, for both tail options.
  for (tails in c("one", "two")) {
    se <- calculate_power(slice_spec("sensitivity", tails, alpha = 0.05, power = 0.05, n1 = 20, n2 = 20))
    expect_identical(se$status$code, "infeasible", info = tails)
    expect_match(se$status$reason, "touched there but not crossed", info = tails)
    expect_identical(slice_numbers(se), character())
  }
  # Very unequal and very large designs.
  big <- calculate_power(slice_spec("sensitivity", "two", alpha = 0.05, power = 0.8, n1 = 250000, n2 = 250000))
  expect_identical(big$status$code, "ok")
  slice_expect_prob(ref_power_t2(big$solved$d, 250000, 250000, 0.05, "two")$power, 0.8, "large design")
  # The scan reaches ncp of several million, where the tails are far below
  # the double range; no quadrature message may leak into the result.
  expect_identical(big$warnings, character())
  lop <- slice_post_hoc(0.5, 1e6, 1)
  slice_expect_prob(lop$solved$power, ref_power_t2(0.5, 1e6, 1, 0.05, "two")$power, "1e6 vs 1")
})

# Documented examples (documented_example evidence; not G*Power validation) ----------

# Values as transcribed in validation/documented-examples.csv (and the
# fixture records validation/fixtures/documented_example/<id>--r01.json);
# the last test checks that this table still equals the ledger.
slice_doc_apriori <- list(
  "DOC-001" = list(d = "0.5", alpha = "0.05", power = "0.95",
                   out = c(ncp = "3.3166248", critical_t = "1.6536580", df = "174", n1 = "88",
                           n2 = "88", n_total = "176", actual_power = "0.9514254")),
  "DOC-003" = list(d = "2.0000000", alpha = "0.05", power = "0.95",
                   out = c(ncp = "3.7416574", critical_t = "1.7822876", df = "12", n1 = "7",
                           n2 = "7", n_total = "14", actual_power = "0.9695815"))
)

slice_doc_xy <- data.frame(
  id = sprintf("DOC-%03d", 5:17),
  power = c("0.600000", "0.610000", "0.620000", "0.630000", "0.640000", "0.650000", "0.660000",
            "0.670000", "0.680000", "0.690000", "0.700000", "0.710000", "0.720000"),
  d010 = c("1442.62", "1482.33", "1522.88", "1564.30", "1606.66", "1650.00", "1694.38", "1739.86",
           "1786.51", "1834.40", "1883.62", "1934.25", "1986.40"),
  d015 = c("641.920", "659.568", "677.588", "696.000", "714.825", "734.086", "753.810", "774.023",
           "794.757", "816.043", "837.918", "860.423", "883.600"),
  d020 = c("361.676", "371.604", "381.740", "392.096", "402.685", "413.520", "424.614", "435.984",
           "447.646", "459.620", "471.925", "484.583", "497.621"),
  stringsAsFactors = FALSE
)

test_that("documented examples DOC-001 and DOC-003 (a priori) are reproduced", {
  for (id in names(slice_doc_apriori)) {
    ex <- slice_doc_apriori[[id]]
    info <- paste("documented_example", id)
    res <- solve_design(slice_spec("a_priori", "one", alpha = as.numeric(ex$alpha),
                                   power = as.numeric(ex$power), d = as.numeric(ex$d)))
    expect_identical(res$status$code, "ok", info = info)
    # MR-EXACT: integers.
    expect_identical(res$solved$n_total, as.integer(ex$out[["n_total"]]), info = info)
    expect_identical(unname(res$group_sizes), as.integer(ex$out[c("n1", "n2")]), info = info)
    expect_identical(res$quantities$df, as.numeric(ex$out[["df"]]), info = info)
    # MR-DISPLAY: printed reals.
    slice_expect_display(res$quantities$ncp, ex$out[["ncp"]], info)
    slice_expect_display(res$quantities$critical, ex$out[["critical_t"]], info)
    slice_expect_display(res$achieved$power, ex$out[["actual_power"]], info)
  }
})

test_that("documented examples DOC-005 to DOC-017 (X-Y plot table) are reproduced", {
  # The manual's X-Y plot table lists the total sample size at which the
  # power curve reaches the target with fractional group sizes N/2. The a
  # priori solver reports exactly this quantity as its continuous relaxation
  # (a diagnostic, never as the solved design).
  for (i in seq_len(nrow(slice_doc_xy))) {
    row <- slice_doc_xy[i, ]
    for (dd in c(d010 = 0.1, d015 = 0.15, d020 = 0.2)) {
      col <- names(which(c(d010 = 0.1, d015 = 0.15, d020 = 0.2) == dd))
      res <- solve_design(slice_spec("a_priori", "one", alpha = 0.05, power = as.numeric(row$power), d = dd))
      slice_expect_display(res$diagnostics$continuous_relaxation$n_total, row[[col]],
                           sprintf("documented_example %s d = %g", row$id, dd))
      expect_gte(res$solved$n_total, res$diagnostics$continuous_relaxation$n_total)
    }
  }
})

test_that("documented example DOC-098 (post hoc) and its recorded manual misprint", {
  res <- slice_post_hoc(0.51639778, 82, 82, alpha = 0.05, tails = "one")
  ref <- ref_power_t2(0.51639778, 82, 82, 0.05, "one")$power
  slice_expect_prob(res$solved$power, ref, "DOC-098 independent reference")
  # The manual prints 0.930308 but states that the power is exactly that of
  # the preceding example (0.950308). Every engine gives 0.9503082; the
  # printed value is logged as a probable misprint (specs/engine-decisions.md
  # 5.2; validation/discrepancy-log.csv) and is not used as a fixture value.
  expect_false(slice_display_ok(res$solved$power, "0.930308"))
  slice_expect_display(res$solved$power, "0.950308", "documented_example DOC-098 (DOC-097 value)")
})

test_that("documented examples DOC-002 and DOC-099 (effect size) are reproduced", {
  es <- effect_size(slice_id, "two_sds", mean_1 = 0, mean_2 = 1, sd_1 = 0.5, sd_2 = 0.5)
  expect_identical(es$status$code, "ok")
  # The drawer image shows d = 2; d = (mu1 - mu2) / sigma is -2 (manual
  # 21.1). Which sign G*Power displays is a capture item.
  expect_identical(es$magnitude, 2)
  expect_identical(es$value, -2)
  slice_expect_display(es$magnitude, "2", "documented_example DOC-002")
  hand <- effect_size(slice_id, "point_biserial_r", r = 0.25, n1 = 64, n2 = 100)
  slice_expect_display(hand$value, "0.52930772", "documented_example DOC-099")
  # With equal groups the conversion reduces to 2 r / sqrt(1 - r^2), the d of DOC-098.
  eq <- effect_size(slice_id, "point_biserial_r", r = 0.25, n1 = 82, n2 = 82)
  slice_expect_display(eq$value, "0.51639778", "DOC-098 input d")
})

test_that("the embedded documented values equal the transcription ledger (development checkout)", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
  csv <- file.path(root, "validation", "documented-examples.csv")
  skip_if_not(file.exists(csv), "validation ledger not available (installed package)")
  ex <- utils::read.csv(csv, stringsAsFactors = FALSE, encoding = "UTF-8")
  ex <- ex[ex$method_id == slice_id, ]
  led <- function(id, field) jsonlite::fromJSON(ex[ex$example_id == id, field])
  labels <- c(ncp = "Noncentrality parameter \u03b4", critical_t = "Critical t", df = "Df",
              n1 = "Sample size group 1", n2 = "Sample size group 2", n_total = "Total sample size",
              actual_power = "Actual power")
  for (id in names(slice_doc_apriori)) {
    out <- led(id, "outputs_json")
    inp <- led(id, "inputs_json")
    expect_identical(unlist(out[labels]), stats::setNames(slice_doc_apriori[[id]]$out, labels), info = id)
    expect_identical(inp[["Effect size d"]], slice_doc_apriori[[id]]$d)
    expect_identical(inp[["Tail(s)"]], "One")
    expect_identical(inp[["Allocation ratio N2/N1"]], "1")
  }
  for (i in seq_len(nrow(slice_doc_xy))) {
    id <- slice_doc_xy$id[i]
    out <- led(id, "outputs_json")
    expect_identical(unname(unlist(out)), unname(unlist(slice_doc_xy[i, c("d010", "d015", "d020")])), info = id)
    expect_identical(led(id, "inputs_json")[["Power (1-\u03b2 err prob)"]], slice_doc_xy$power[i])
  }
  expect_identical(led("DOC-098", "outputs_json")[[1]], "0.930308")
  expect_identical(led("DOC-098", "inputs_json")[["Effect size d"]], "0.51639778")
  expect_identical(led("DOC-002", "outputs_json")[["Effect size d"]], "2")
  expect_identical(led("DOC-099", "outputs_json")[["d"]], "0.52930772")
  # Fixture records carry the same printed strings.
  fx_dir <- file.path(root, "validation", "fixtures", "documented_example")
  for (id in c(names(slice_doc_apriori), slice_doc_xy$id)) {
    f <- file.path(fx_dir, paste0(id, "--r01.json"))
    expect_true(file.exists(f), info = id)
    rec <- jsonlite::fromJSON(f, simplifyVector = FALSE)
    expect_identical(rec$evidence_kind, "documented_example")
    shown <- vapply(rec$outputs, `[[`, "", "displayed")
    names(shown) <- vapply(rec$outputs, `[[`, "", "label_as_displayed")
    expect_identical(shown, unlist(led(id, "outputs_json"))[names(shown)], info = id)
  }
})

# Shared result object ------------------------------------------------------------------------

test_that("results of every mode share one object type and round-trip through JSON", {
  results <- list(
    solve_design(slice_spec("a_priori", "two", alpha = 0.05, power = 0.8, d = 0.5)),
    slice_post_hoc(0.5, 30, 60),
    calculate_power(slice_spec("sensitivity", "one", alpha = 0.05, power = 0.8, n1 = 30, n2 = 60)),
    calculate_power(slice_spec("criterion", "two", power = 0.8, n1 = 30, n2 = 60, d = 0.5)),
    calculate_power(slice_spec("compromise", "one", error_ratio = 2, n1 = 30, n2 = 60, d = 0.5)),
    solve_design(slice_spec("a_priori", alpha = 0.05, power = 0.8, d = 0.5, allocation_ratio = 3))
  )
  for (res in results) {
    expect_s3_class(res, "power_result")
    expect_identical(ig_from_json(ig_to_json(res)), res, info = res$mode)
    expect_match(res$r_call, "^inferguidR::(calculate_power|solve_design)\\(inferguidR::analysis_spec\\(")
    if (requireNamespace("jsonvalidate", quietly = TRUE)) {
      chk <- .ig_validate_schema(ig_to_json(res), "power_result")
      expect_true(chk$valid, info = paste(res$mode, paste(chk$errors, collapse = "; ")))
    }
    expect_identical(res$warnings, character())
  }
  # The reproducible call recreates the identical specification. The string
  # is generated by the package itself in this test (not imported data), so
  # parsing and evaluating it here is the point of the check.
  call_text <- sub("^inferguidR::solve_design\\((.*)\\)$", "\\1", results[[1]]$r_call)
  rebuilt <- eval(parse(text = sub("^inferguidR::", "", call_text)))
  expect_identical(rebuilt$spec_hash, results[[1]]$spec_hash)
})

# Effect-size drawer -------------------------------------------------------------------------

test_that("effect_size() computes d with its formula, source and caveats", {
  cs <- effect_size(slice_id, "common_sd", mean_1 = 10, mean_2 = 7, sd = 4)
  expect_identical(cs$value, 0.75)
  expect_identical(cs$effect, "d")
  expect_match(cs$formula, "\\(mean_1 - mean_2\\) / sd")
  expect_identical(cs$warnings, character())
  two <- effect_size(slice_id, "two_sds", inputs = list(mean_1 = 1, mean_2 = 0, sd_1 = 1, sd_2 = 2))
  expect_equal(two$value, 1 / sqrt(2.5), tolerance = 1e-15)
  expect_match(two$warnings, "equal group sizes")
  # format.ig_effect_size() is called directly so that the check does not
  # depend on S3 registration in a development load.
  expect_true(any(grepl("d = 0.632", format.ig_effect_size(two), fixed = TRUE)))
  cond <- tryCatch(effect_size(slice_id, "two_sds", mean_1 = 1, mean_2 = 0, sd_1 = -1, sd_2 = 1, extra = 2),
                   inferguidR_invalid_input = function(e) e)
  problems <- unlist(cond$status$details$problems)
  expect_true(any(grepl("'sd_1' must be positive", problems)))
  expect_true(any(grepl("input 'extra' is not used", problems)))
  expect_s3_class(tryCatch(effect_size(slice_id, "nope", r = 1), error = function(e) e),
                  "inferguidR_invalid_input")
  expect_s3_class(tryCatch(effect_size(slice_id, "point_biserial_r", r = 1, n1 = 5, n2 = 5),
                           error = function(e) e), "inferguidR_invalid_input")
  expect_s3_class(tryCatch(effect_size(slice_id, "common_sd", mean_1 = 1, inputs = list(mean_1 = 2),
                                       mean_2 = 0, sd = 1), error = function(e) e),
                  "inferguidR_invalid_input")
  # The generic tests have no drawer in G*Power (registry: drawer null).
  none <- effect_size("t_generic", "any", delta = 2)
  expect_identical(none$status$code, "unsupported")
  expect_true(is.na(none$value))
})

# Sweeps and comparisons ----------------------------------------------------------------------

test_that("sweep_design() keeps one full specification and result per value, gaps included", {
  base <- slice_spec("a_priori", "two", alpha = 0.05, power = 0.8, d = 0.5)
  sw <- sweep_design(base, "d", c(0.6, 0, 0.5, 0.4))
  expect_identical(sw$value, c(0.6, 0, 0.5, 0.4))
  expect_identical(sw$status, c("ok", "infeasible", "ok", "ok"))
  expect_true(all(sw$admissible))
  expect_true(is.na(sw$solved_n_total[2]))
  expect_match(sw$reason[2], "power does not change with N")
  expect_identical(sw$solved_n_total[c(1, 3, 4)], c(90, 128, 200))
  expect_identical(sw$spec_hash[3], base$spec_hash)
  expect_identical(sw$result[[3]]$result_hash, calculate_power(base)$result_hash)
  expect_identical(attr(sw, "free"), "d")
  expect_identical(attr(sw, "base_spec_hash"), base$spec_hash)
  expect_false("d" %in% names(attr(sw, "fixed_inputs")))

  ph <- slice_spec("post_hoc", "two", alpha = 0.05, n1 = 20, n2 = 20, d = 0.5)
  sw <- sweep_design(ph, "n1", c(10, 10.5, 30))
  expect_identical(sw$admissible, c(TRUE, FALSE, TRUE))
  expect_identical(sw$status[2], "invalid_input")
  expect_null(sw$result[[2]])
  expect_true(is.na(sw$achieved_power[2]))
  expect_lt(sw$achieved_power[1], sw$achieved_power[3])
  expect_identical(sw$n_total, c(30, NA, 50))

  sw <- sweep_design(base, "allocation_ratio", c(1, 2))
  expect_identical(sw$status, c("ok", "ok"))
  expect_identical(sw$solved_n_total, c(128, 144))
  expect_s3_class(tryCatch(sweep_design(base, "n_total", 1:3), error = function(e) e),
                  "inferguidR_invalid_input")
  expect_s3_class(tryCatch(sweep_design(base, "n1", 1:3), error = function(e) e),
                  "inferguidR_invalid_input")
  expect_s3_class(tryCatch(sweep_design(base, "d", character()), error = function(e) e),
                  "inferguidR_invalid_input")
})

test_that("compare_designs() tabulates specifications and results side by side", {
  base <- slice_spec("post_hoc", "two", alpha = 0.05, n1 = 30, n2 = 30, d = 0.5)
  ap <- solve_design(slice_spec("a_priori", "two", alpha = 0.05, power = 0.8, d = 0.5))
  pend <- slice_spec("a_priori", "two", alpha = 0.05, power = 0.8, d = 0, allocation_ratio = 2)
  cmp <- compare_designs(equal = base, unequal = update(base, inputs = list(n1 = 20, n2 = 40)), ap,
                         designs = list(pending = pend))
  expect_identical(cmp$design, c("equal", "unequal", "design_3", "pending"))
  expect_identical(cmp$status, c("ok", "ok", "ok", "infeasible"))
  expect_identical(cmp$n_total, c(60, 60, 128, NA))
  expect_gt(cmp$achieved_power[1], cmp$achieved_power[2])
  expect_identical(cmp$result_hash[3], ap$result_hash)
  expect_identical(cmp$critical[[1]], calculate_power(base)$quantities$critical)
  expect_null(cmp$solved[[4]])
  expect_true(is.na(cmp$achieved_power[4]))
  expect_s3_class(tryCatch(compare_designs(list(1)), error = function(e) e), "inferguidR_invalid_input")
  expect_s3_class(tryCatch(compare_designs(), error = function(e) e), "inferguidR_invalid_input")
})
