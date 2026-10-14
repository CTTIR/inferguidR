# G11 method definitions: registration, specifications, conventions, typed
# outcomes, allocation rule, properties and drawers for
# f_variances_two_groups, chi2_variance_one_sample and
# chi2_goodness_of_fit_contingency (R/methods-variances.R, R/effect-variances.R).
#
# Conventions tested here come from the certified live inventory and the
# supplement captures named in specs/method-notes/<id>.md; numerical
# correctness rests on test-methods-g11-reference.R and
# validation/independent/g11/.

ref_reset_registry()

g11_ids <- c("f_variances_two_groups", "chi2_variance_one_sample", "chi2_goodness_of_fit_contingency")

g11_spec <- function(method, mode, inputs, options = list()) analysis_spec(method, mode, inputs = inputs, options = options)
g11_calc <- function(...) calculate_power(g11_spec(...))
g11_numbers <- function(res) {
  x <- unclass(res)
  .ig_finite_paths(x[c("solved", "achieved", "quantities", "group_sizes", "method_outputs")], "")
}
g11_code <- function(expr) tryCatch({
  force(expr)
  "ok"
}, inferguidR_error = function(e) e$status$code)

test_that("the three methods are implemented in all five modes with an honest validation state", {
  lm <- list_methods()
  for (id in g11_ids) {
    row <- lm[lm$id == id, ]
    expect_identical(row$implementation_status, "implemented", info = id)
    expect_identical(row$validation_status, "independently_checked", info = id)
    expect_identical(row$modes, "a_priori, post_hoc, sensitivity, criterion, compromise", info = id)
    vs <- validation_status(id)
    expect_identical(vs$gpower_reference$status, "reference_pending", info = id)
    expect_identical(vs$independent_check$checker, "validation/independent/g11/check_g11.R", info = id)
    expect_identical(vs$modes_with_solver, .ig_modes, info = id)
  }
  expect_identical(.ig_get_method("f_variances_two_groups")$family, "F")
  expect_identical(.ig_get_method("chi2_goodness_of_fit_contingency")$family, "chi2")
})

test_that("parameters and options follow the certified live inventory", {
  # Live defaults: Tail(s) One, Allocation ratio N2/N1 1, Effect direction var1>=var0 (sensitivity only).
  sp <- g11_spec("f_variances_two_groups", "a_priori", list(ratio = 1.5, alpha = 0.05, power = 0.8))
  expect_identical(sp$options$tails, "one")
  expect_identical(sp$inputs$allocation_ratio, 1)
  expect_true(all(c("allocation_ratio") %in% sp$defaults_applied$inputs))
  expect_identical(sp$defaults_applied$options, "tails")
  se <- g11_spec("f_variances_two_groups", "sensitivity", list(alpha = 0.05, power = 0.8, n1 = 50, n2 = 50))
  expect_identical(se$options$effect_direction, "var1_ge_var0")
  expect_identical(g11_code(g11_spec("f_variances_two_groups", "post_hoc", list(ratio = 1.5, alpha = 0.05, n1 = 50, n2 = 50),
                                     list(effect_direction = "var1_le_var0"))), "invalid_input")
  expect_identical(g11_code(g11_spec("f_variances_two_groups", "post_hoc",
                                     list(ratio = 1.5, alpha = 0.05, n1 = 50, n2 = 50, allocation_ratio = 2))), "invalid_input")
  # Live labels map onto the parameters (registry inputs_by_mode_live).
  live <- .ig_live_labels("f_variances_two_groups", "sensitivity")
  expect_true(all(c("Tail(s)", "Effect direction", "Sample size group 1", "Sample size group 2") %in% live$inputs))
  expect_identical(.ig_get_method("f_variances_two_groups")$options$effect_direction$labels[["var1_le_var0"]], "var1\u2264var0")
  expect_identical(g11_spec("chi2_variance_one_sample", "post_hoc", list(ratio = 2, alpha = 0.05, n_total = 10))$options$tails, "one")
  expect_identical(names(.ig_get_method("chi2_goodness_of_fit_contingency")$options), character())
})

test_that("domains and consistency checks refuse inputs as G*Power does, with typed invalid_input", {
  f <- "f_variances_two_groups"
  expect_identical(g11_code(g11_spec(f, "post_hoc", list(ratio = 1.5, alpha = 1.5, n1 = 50, n2 = 50))), "invalid_input")
  expect_identical(g11_code(g11_spec(f, "post_hoc", list(ratio = -1, alpha = 0.05, n1 = 50, n2 = 50))), "invalid_input")
  expect_identical(g11_code(g11_spec(f, "post_hoc", list(ratio = 1.5, alpha = 0.05, n1 = 1, n2 = 50))), "invalid_input")
  expect_identical(g11_code(g11_spec(f, "compromise", list(ratio = 1.5, error_ratio = 0, n1 = 50, n2 = 50))), "invalid_input")
  # Power below alpha (G*Power: "Power must be greater than or equal to alpha.").
  for (mode in c("a_priori", "sensitivity")) {
    ins <- list(alpha = 0.05, power = 0.04)
    ins <- c(ins, if (mode == "a_priori") list(ratio = 1.5) else list(n1 = 50, n2 = 50))
    expect_identical(g11_code(g11_spec(f, mode, ins)), "invalid_input", info = mode)
  }
  expect_identical(g11_code(g11_spec("chi2_variance_one_sample", "post_hoc", list(ratio = 1.5, alpha = 0.05, n_total = 1))),
                   "invalid_input")
  g <- "chi2_goodness_of_fit_contingency"
  expect_identical(g11_code(g11_spec(g, "post_hoc", list(w = 0.3, alpha = 0.05, n_total = 150, df = 0))), "invalid_input")
  expect_identical(g11_code(g11_spec(g, "post_hoc", list(w = -0.3, alpha = 0.05, n_total = 150, df = 5))), "invalid_input")
  expect_identical(g11_code(g11_spec(g, "post_hoc", list(w = 0.3, alpha = 0.05, n_total = 0, df = 5))), "invalid_input")
  # Real degrees of freedom >= 1 are admissible (G*Power accepts Df = 2.5).
  expect_identical(g11_calc(g, "post_hoc", list(w = 0.3, alpha = 0.05, n_total = 150, df = 2.5))$status$code, "ok")
})

test_that("zero effects in a priori mode are infeasible and report no number", {
  for (cs in list(list("f_variances_two_groups", list(ratio = 1, alpha = 0.05, power = 0.8), list(tails = "one")),
                  list("chi2_variance_one_sample", list(ratio = 1, alpha = 0.05, power = 0.8), list(tails = "two")),
                  list("chi2_goodness_of_fit_contingency", list(w = 0, alpha = 0.05, power = 0.8, df = 5), list()))) {
    res <- solve_design(g11_spec(cs[[1]], "a_priori", cs[[2]], cs[[3]]))
    expect_identical(res$status$code, "infeasible", info = cs[[1]])
    expect_match(res$status$reason, "not reached by any admissible design", info = cs[[1]])
    expect_identical(g11_numbers(res), character(), info = cs[[1]])
  }
})

test_that("goodness-of-fit sensitivity admits the documented zero-effect boundary", {
  se <- g11_calc("chi2_goodness_of_fit_contingency", "sensitivity", list(alpha = 0.05, power = 0.05, n_total = 150, df = 5))
  expect_identical(se$status$code, "ok")
  expect_identical(se$solved$w, 0)
  expect_equal(se$quantities$ncp, 0)
  expect_equal(se$achieved$power, .05, tolerance = 1e-12)
})

test_that("tail conventions: one tail on the side of the ratio, ratio 1 upper for F and lower for chi2", {
  f <- function(ratio, tails) g11_calc("f_variances_two_groups", "post_hoc", list(ratio = ratio, alpha = 0.05, n1 = 20, n2 = 30),
                                       list(tails = tails))
  r1 <- f(1, "one")
  expect_equal(r1$quantities$critical, stats::qf(0.05, 19, 29, lower.tail = FALSE), tolerance = 1e-12)
  expect_equal(r1$achieved$power, 0.05, tolerance = 1e-12)
  expect_identical(r1$method_outputs$statistic_distribution$rejection_side, "upper")
  rl <- f(0.5, "one")
  expect_equal(rl$quantities$critical, stats::qf(0.05, 19, 29), tolerance = 1e-12)
  expect_identical(rl$method_outputs$statistic_distribution$rejection_side, "lower")
  r2 <- f(1.5, "two")
  expect_equal(r2$quantities$critical, c(stats::qf(0.025, 19, 29), stats::qf(0.025, 19, 29, lower.tail = FALSE)), tolerance = 1e-12)
  expect_equal(r2$achieved$power + r2$achieved$beta, 1, tolerance = 1e-14)
  c1 <- g11_calc("chi2_variance_one_sample", "post_hoc", list(ratio = 1, alpha = 0.05, n_total = 20), list(tails = "one"))
  expect_equal(c1$quantities$critical, stats::qchisq(0.05, 19), tolerance = 1e-12)
  expect_identical(c1$method_outputs$statistic_distribution$rejection_side, "lower")
  c2 <- g11_calc("chi2_variance_one_sample", "post_hoc", list(ratio = 1.5, alpha = 0.05, n_total = 20), list(tails = "one"))
  expect_equal(c2$quantities$critical, stats::qchisq(0.05, 19, lower.tail = FALSE), tolerance = 1e-12)
  # Ratio 0 puts the statistic at 0: power 1 in the lower rejection region.
  z <- g11_calc("chi2_variance_one_sample", "post_hoc", list(ratio = 0, alpha = 0.05, n_total = 20), list(tails = "one"))
  expect_identical(c(z$achieved$power, z$achieved$beta), c(1, 0))
})

test_that("exchanging the groups and inverting the ratio leaves the F-test power unchanged", {
  for (tails in c("one", "two")) for (ratio in c(0.3, 1.7)) {
    a <- g11_calc("f_variances_two_groups", "post_hoc", list(ratio = ratio, alpha = 0.05, n1 = 12, n2 = 40), list(tails = tails))
    b <- g11_calc("f_variances_two_groups", "post_hoc", list(ratio = 1 / ratio, alpha = 0.05, n1 = 40, n2 = 12), list(tails = tails))
    expect_equal(a$achieved$power, b$achieved$power, tolerance = 1e-12, info = paste(tails, ratio))
    expect_equal(sort(a$quantities$critical), sort(1 / b$quantities$critical), tolerance = 1e-12, info = paste(tails, ratio))
  }
})

test_that("critical values are exact where stats::qf switches to its large-df approximation", {
  res <- g11_calc("f_variances_two_groups", "post_hoc", list(ratio = 1.02, alpha = 0.001, n1 = 1647267, n2 = 164727),
                  list(tails = "one"))
  cr <- res$quantities$critical
  expect_equal(stats::pf(cr, 1647266, 164726, lower.tail = FALSE), 0.001, tolerance = 1e-10)
  # qf() returns df2 / qchisq(., df2) for df1 > 4e5; its upper tail is 0.0016 here.
  expect_gt(abs(stats::pf(stats::qf(0.001, 1647266, 164726, lower.tail = FALSE), 1647266, 164726, lower.tail = FALSE) - 0.001), 5e-4)
})

test_that("a priori allocation: ratio 1 gives equal groups; other ratios adjust rounded designs to reach target", {
  f <- "f_variances_two_groups"
  eq <- solve_design(g11_spec(f, "a_priori", list(ratio = 1.5, alpha = 0.05, power = 0.8), list(tails = "one")))
  expect_identical(eq$status$code, "ok")
  expect_identical(eq$group_sizes[["n1"]], eq$group_sizes[["n2"]])
  n <- eq$group_sizes[["n1"]]
  prev <- g11_calc(f, "post_hoc", list(ratio = 1.5, alpha = 0.05, n1 = n - 1, n2 = n - 1), list(tails = "one"))
  expect_lt(prev$achieved$power, 0.8)
  # Rule for r != 1: N* = smallest N whose fractional design reaches the target;
  # n1 = floor(N*/(1 + r) + 1/2), n2 = N* - n1.
  fw <- .ig_power_f_variances_two_groups
  rule <- function(ratio, alpha, power, r, tails) {
    fp <- function(N) fw(list(ratio = ratio, alpha = alpha, tails = tails, n1 = N / (1 + r), n2 = r * N / (1 + r)))$power
    N <- 4
    while (fp(N) < power) N <- N + 1
    n1 <- floor(N / (1 + r) + 0.5)
    c(n1 = n1, n2 = N - n1)
  }
  for (cs in list(list(1.5, 0.05, 0.8, 1.000001, "one"), list(0.5, 0.05, 0.8, 2, "one"), list(0.7, 0.05, 0.8, 1.5, "two"))) {
    res <- solve_design(g11_spec(f, "a_priori", list(ratio = cs[[1]], alpha = cs[[2]], power = cs[[3]], allocation_ratio = cs[[4]]),
                                 list(tails = cs[[5]])))
    info <- paste(unlist(cs), collapse = " ")
    expect_identical(res$status$code, "ok", info = info)
    expected <- rule(cs[[1]], cs[[2]], cs[[3]], cs[[4]], cs[[5]])
    expect_identical(names(res$group_sizes), names(expected), info = info)
    expect_identical(as.double(res$group_sizes), as.double(expected), info = info)
    expect_gte(res$achieved$power, cs[[3]])
    expect_match(res$diagnostics$allocation_rule, "fractional group sizes", info = info)
  }
  # Rounding G*Power's N=427 design misses the requested power.
  dec <- solve_design(g11_spec(f, "a_priori", list(ratio = 1.5, alpha = 0.05, power = 0.8, allocation_ratio = 2),
                               list(tails = "two")))
  expect_identical(dec$status$code, "ok")
  expect_gte(dec$achieved$power, 0.8)
  expect_true(dec$diagnostics$adjusted_from_gpower_rule)
  expect_equal(dec$diagnostics$n_total_fractional_rule, 427)
  expect_match(paste(dec$warnings, collapse = " "), "below.*target")
  # Independent central-F calculation verifies the first attaining rounded split.
  ref <- function(N) {
    n1 <- floor(N / 3 + 0.5); n2 <- N - n1
    pf(qf(0.025, n1-1, n2-1) / 1.5, n1-1, n2-1) +
      pf(qf(0.025, n1-1, n2-1, lower.tail=FALSE) / 1.5, n1-1, n2-1, lower.tail=FALSE)
  }
  N <- 427
  while (ref(N) < 0.8) N <- N + 1
  expect_equal(dec$solved$n_total, N)
  expect_lt(ref(N - 1), 0.8)
  expect_equal(dec$achieved$power, ref(N), tolerance = 1e-10)
  bounded <- .ig_apriori_f_variances_two_groups(
    list(ratio = 1.5, alpha = 0.05, power = 0.8, allocation_ratio = 2, tails = "two"), 4, 427)
  expect_identical(bounded$status, "infeasible")
  expect_match(bounded$reason, "sample-size bound")
  # A split with a group below 2 subjects moves to the next admissible total.
  expect_null(.ig_design_f_variances_two_groups(list(allocation_ratio = 10), 12))
  adj <- solve_design(g11_spec(f, "a_priori", list(ratio = 1e4, alpha = 0.05, power = 0.8, allocation_ratio = 10),
                               list(tails = "one")))
  expect_identical(adj$status$code, "ok")
  expect_true(all(adj$group_sizes >= 2))
  expect_match(adj$diagnostics$split_adjustment, "fewer than 2")
})

test_that("a priori search for unequal groups starts at the first total with a computable criterion", {
  f <- "f_variances_two_groups"
  fw <- .ig_power_f_variances_two_groups
  # Allocation 0.2: the first total with both fractional group sizes above 1 is N = 6 in double arithmetic
  # ((1 + r) / r = 5.999...), with n2 = 1 + 2e-16 and a critical value far outside the double range.
  expect_true(is.na(fw(list(ratio = 2, alpha = 0.05, tails = "one", n1 = 6 / 1.2, n2 = 0.2 * 6 / 1.2))$power))
  # Allocation 0.3334 and 2.9999: df 1.5e-4 and 2.5e-5 at the first total (N = 4).
  for (cs in list(c(0.2, 2), c(0.3334, 1.8), c(2.9999, 1.8))) {
    r <- cs[1]
    p <- list(ratio = cs[2], alpha = 0.05, tails = "one")
    res <- solve_design(g11_spec(f, "a_priori", list(ratio = cs[2], alpha = 0.05, power = 0.8, allocation_ratio = r), list(tails = "one")))
    info <- paste("allocation", r)
    expect_identical(res$status$code, "ok", info = info)
    expect_match(res$diagnostics$fractional_search_start, "not computable", info = info)
    expect_null(res$diagnostics$split_adjustment, info = info)
    expect_identical(res$warnings, character(), info = info)
    # The rule, evaluated here from the first total with a computable criterion.
    fp <- function(N) fw(c(p, list(n1 = N / (1 + r), n2 = r * N / (1 + r))))$power
    N <- 4
    repeat {
      v <- if (min(N / (1 + r), r * N / (1 + r)) > 1) suppressWarnings(fp(N)) else NA_real_
      if (!is.na(v) && v >= 0.8) break
      N <- N + 1
    }
    n1 <- floor(N / (1 + r) + 0.5)
    expect_identical(as.double(res$group_sizes[c("n1", "n2")]), c(n1, N - n1), info = info)
    expect_identical(as.double(res$solved$n_total), N, info = info)
    expect_gt(N, 50)
    expect_gte(res$achieved$power, 0.8)
  }
  # No allocation ratio between 0.01 and 10 takes a design with a df next to 0 for the solution.
  for (r in c(0.1667, 0.2001, 0.25, 0.3333, 0.5, 0.7, 1.9999, 3.9999, 7)) {
    res <- solve_design(g11_spec(f, "a_priori", list(ratio = 2, alpha = 0.05, power = 0.8, allocation_ratio = r), list(tails = "one")))
    expect_identical(res$status$code, "ok", info = paste("allocation", r))
    expect_null(res$diagnostics$split_adjustment, info = paste("allocation", r))
    if (identical(res$status$code, "ok")) expect_gt(res$solved$n_total, 50)
  }
})

test_that("a critical value outside the double range is a typed numerical failure, not a number", {
  # The upper 1e-200 quantile of F(1, 1) is about 4e399 and the lower one about 2.5e-400; the lower 1e-200
  # quantile of chi2(1) is about 1.6e-400. No double is such a quantile.
  expect_true(is.nan(.ig_g11_qf(1e-200, 1, 1, FALSE)))
  expect_true(is.nan(.ig_g11_qf(1e-200, 1, 1, TRUE)))
  expect_true(is.nan(.ig_g11_qchisq(1e-200, 1, TRUE)))
  expect_true(is.nan(.ig_g11_qf(0.05, 2, 1e-3, FALSE)))
  for (cs in list(list("f_variances_two_groups", list(ratio = 2, alpha = 1e-200, n1 = 2, n2 = 2), list(tails = "one")),
                  list("f_variances_two_groups", list(ratio = 0.5, alpha = 1e-200, n1 = 2, n2 = 2), list(tails = "one")),
                  list("f_variances_two_groups", list(ratio = 2, alpha = 1e-200, n1 = 2, n2 = 2), list(tails = "two")),
                  list("chi2_variance_one_sample", list(ratio = 0.5, alpha = 1e-200, n_total = 2), list(tails = "one")),
                  list("chi2_variance_one_sample", list(ratio = 2, alpha = 1e-200, n_total = 2), list(tails = "two")))) {
    res <- g11_calc(cs[[1]], "post_hoc", cs[[2]], cs[[3]])
    info <- paste(cs[[1]], cs[[2]]$ratio, cs[[3]]$tails)
    expect_identical(res$status$code, "numerical_failure", info = info)
    expect_match(res$status$reason, "outside the range of double-precision numbers", info = info)
    expect_identical(g11_numbers(res), character(), info = info)
  }
  # Representable critical values keep their relative accuracy in the far tails. Closed forms with 1 df:
  # P(F(1, 1) > x) = (2 / pi) atan(1 / sqrt(x)) and P(chi2(1) <= x) = erf(sqrt(x / 2)), so the power of the
  # one-tailed tests at ratio 2 and 1/2 is alpha * sqrt(2) up to a relative error of order alpha^2.
  for (a in c(1e-20, 1e-100, 1e-150)) {
    up <- g11_calc("f_variances_two_groups", "post_hoc", list(ratio = 2, alpha = a, n1 = 2, n2 = 2), list(tails = "one"))
    lo <- g11_calc("f_variances_two_groups", "post_hoc", list(ratio = 0.5, alpha = a, n1 = 2, n2 = 2), list(tails = "one"))
    ch <- g11_calc("chi2_variance_one_sample", "post_hoc", list(ratio = 0.5, alpha = a, n_total = 2), list(tails = "one"))
    for (res in list(up, lo, ch)) {
      expect_identical(res$status$code, "ok")
      expect_equal(res$achieved$power / a, sqrt(2), tolerance = 1e-10, info = paste(res$method_id, a))
    }
  }
  # The upper critical value of the chi-square tests is representable for every admissible alpha.
  expect_identical(g11_calc("chi2_variance_one_sample", "post_hoc", list(ratio = 2, alpha = 1e-200, n_total = 2),
                            list(tails = "one"))$status$code, "ok")
  expect_identical(g11_calc("chi2_goodness_of_fit_contingency", "post_hoc", list(w = 0.3, alpha = 1e-200, n_total = 50, df = 1))$status$code, "ok")
})

test_that("sensitivity: the effect direction selects the root below or above 1", {
  f <- "f_variances_two_groups"
  for (tails in c("one", "two")) {
    lo <- g11_calc(f, "sensitivity", list(alpha = 0.05, power = 0.8, n1 = 50, n2 = 50), list(tails = tails, effect_direction = "var1_le_var0"))
    hi <- g11_calc(f, "sensitivity", list(alpha = 0.05, power = 0.8, n1 = 50, n2 = 50), list(tails = tails, effect_direction = "var1_ge_var0"))
    expect_identical(c(lo$status$code, hi$status$code), c("ok", "ok"), info = tails)
    expect_lt(lo$solved$ratio, 1)
    expect_gt(hi$solved$ratio, 1)
    expect_equal(lo$achieved$power, 0.8, tolerance = 1e-9)
    expect_equal(hi$achieved$power, 0.8, tolerance = 1e-9)
  }
  # F(1,1) has CDF 2/pi * atan(sqrt(x)). At alpha=.05 and power=.95,
  # the one-tail variance ratio is cot(pi*.025)^4, beyond the old 1e4 scan.
  far <- g11_calc(f, "sensitivity", list(alpha = 0.05, power = 0.95, n1 = 2, n2 = 2), list(tails = "one"))
  expect_identical(far$status$code, "ok")
  expect_equal(far$solved$ratio, 1 / tan(pi * 0.025)^4, tolerance = 1e-9)
  expect_equal(far$achieved$power, 0.95, tolerance = 1e-9)
  mirror <- g11_calc(f, "sensitivity", list(alpha = 0.05, power = 0.95, n1 = 2, n2 = 2),
                     list(tails = "one", effect_direction = "var1_le_var0"))
  expect_identical(mirror$status$code, "ok")
  expect_equal(mirror$solved$ratio * far$solved$ratio, 1, tolerance = 1e-8)

})

test_that("every mode inverts the forward model", {
  cases <- list(
    list("f_variances_two_groups", list(ratio = 1.8, n1 = 25, n2 = 40), "ratio", list(tails = "two")),
    list("chi2_variance_one_sample", list(ratio = 0.6, n_total = 40), "ratio", list(tails = "one")),
    list("chi2_goodness_of_fit_contingency", list(w = 0.25, n_total = 120, df = 3), "w", list()))
  for (cs in cases) {
    m <- cs[[1]]
    base <- cs[[2]]
    ph <- g11_calc(m, "post_hoc", c(base, list(alpha = 0.05)), cs[[4]])
    pw <- ph$achieved$power
    cr <- g11_calc(m, "criterion", c(base, list(power = pw)), cs[[4]])
    expect_equal(cr$solved$alpha, 0.05, tolerance = 1e-8, info = m)
    eff <- base[[cs[[3]]]]
    sens_in <- c(base[setdiff(names(base), cs[[3]])], list(alpha = 0.05, power = pw))
    opts <- cs[[4]]
    if (m == "f_variances_two_groups") opts$effect_direction <- if (eff >= 1) "var1_ge_var0" else "var1_le_var0"
    se <- g11_calc(m, "sensitivity", sens_in, opts)
    if (m != "chi2_variance_one_sample") expect_equal(se$solved[[cs[[3]]]], eff, tolerance = 1e-7, info = m)
    co <- g11_calc(m, "compromise", c(base, list(error_ratio = ph$achieved$beta / 0.05)), cs[[4]])
    expect_equal(co$solved$alpha, 0.05, tolerance = 1e-7, info = m)
    expect_equal(co$achieved$beta, co$spec$inputs$error_ratio * co$solved$alpha, tolerance = 1e-9, info = m)
    expect_identical(co$achieved$alpha_attained, co$solved$alpha)
  }
  # chi2 one-sample sensitivity reports the root above 1 (as G*Power); 0.6 is the other branch.
  se <- g11_calc("chi2_variance_one_sample", "sensitivity", list(n_total = 40, alpha = 0.05,
                 power = g11_calc("chi2_variance_one_sample", "post_hoc", list(ratio = 0.6, n_total = 40, alpha = 0.05),
                                  list(tails = "one"))$achieved$power), list(tails = "one"))
  expect_gt(se$solved$ratio, 1)
  expect_equal(se$diagnostics$other_branch_roots, 0.6, tolerance = 1e-7)
})

test_that("goodness of fit: huge noncentrality gives power exactly 1 without a long mixture", {
  t0 <- proc.time()[["elapsed"]]
  res <- g11_calc("chi2_goodness_of_fit_contingency", "post_hoc", list(w = 100, alpha = 0.05, n_total = 1e5, df = 5))
  expect_identical(c(res$achieved$power, res$achieved$beta), c(1, 0))
  expect_lt(proc.time()[["elapsed"]] - t0, 5)
  ap <- solve_design(g11_spec("chi2_goodness_of_fit_contingency", "a_priori", list(w = 5, alpha = 0.05, power = 0.8, df = 1)))
  expect_identical(ap$solved$n_total, 1L)
})

test_that("results of every mode round-trip through JSON and carry their assumptions", {
  results <- list(
    solve_design(g11_spec("f_variances_two_groups", "a_priori", list(ratio = 0.5, alpha = 0.05, power = 0.8, allocation_ratio = 2),
                          list(tails = "one"))),
    g11_calc("f_variances_two_groups", "sensitivity", list(alpha = 0.05, power = 0.8, n1 = 30, n2 = 30), list(effect_direction = "var1_le_var0")),
    g11_calc("chi2_variance_one_sample", "criterion", list(ratio = 1.5, power = 0.8, n_total = 50), list(tails = "two")),
    g11_calc("chi2_variance_one_sample", "compromise", list(ratio = 1.5, error_ratio = 1, n_total = 50)),
    g11_calc("chi2_goodness_of_fit_contingency", "post_hoc", list(w = 0.3, alpha = 0.05, n_total = 150, df = 5))
  )
  for (res in results) {
    expect_identical(res$status$code, "ok", info = res$method_id)
    expect_identical(ig_from_json(ig_to_json(res)), res, info = paste(res$method_id, res$mode))
    expect_true(length(res$assumptions) >= 3L)
    if (requireNamespace("jsonvalidate", quietly = TRUE)) {
      chk <- .ig_validate_schema(ig_to_json(res), "power_result")
      expect_true(chk$valid, info = paste(res$method_id, res$mode, paste(chk$errors, collapse = "; ")))
    }
  }
})

test_that("effect-size drawers validate their inputs and report formula and source", {
  es <- effect_size("f_variances_two_groups", "from_variances", var0 = 1.5, var1 = 1)
  expect_identical(es$effect, "ratio")
  expect_match(es$formula, "var1 / var0")
  expect_identical(g11_code(effect_size("f_variances_two_groups", "from_variances", var0 = 0, var1 = 1)), "invalid_input")
  expect_identical(g11_code(effect_size("chi2_variance_one_sample", "from_variances", var0 = 1)), "invalid_input")
  expect_identical(g11_code(effect_size("chi2_variance_one_sample", "nope", var0 = 1, var1 = 2)), "invalid_input")
  g <- "chi2_goodness_of_fit_contingency"
  expect_identical(g11_code(effect_size(g, "cell_probabilities", p0 = c(0.5, 0.6), p1 = c(0.5, 0.5))), "invalid_input")
  expect_identical(g11_code(effect_size(g, "cell_probabilities", p0 = c(0, 1), p1 = c(0.5, 0.5))), "invalid_input")
  expect_identical(g11_code(effect_size(g, "cell_probabilities", p0 = c(0.5, 0.5), p1 = c(0.2, 0.3, 0.5))), "invalid_input")
  tab <- effect_size(g, "contingency_table", p1 = matrix(c(0.3, 0.2, 0.1, 0.4), 2))
  expect_identical(tab$status$code, "ok")
  expect_match(tab$warnings, "Df = 1")
  indep <- effect_size(g, "contingency_table", p1 = outer(c(0.4, 0.6), c(0.3, 0.7)))
  expect_lt(indep$value, 1e-12)
})

test_that("F variance sensitivity selects either direction without invalid scan points", {
  for (ratio in c(0.5, 2)) {
    post <- g11_calc("f_variances_two_groups", "post_hoc",
      list(ratio = ratio, n1 = 30, n2 = 30, alpha = 0.05), list(tails = "two"))
    res <- g11_calc("f_variances_two_groups", "sensitivity",
      list(n1 = 30, n2 = 30, alpha = 0.05, power = post$achieved$power),
      list(tails = "two", effect_direction = if (ratio < 1) "var1_le_var0" else "var1_ge_var0"))
    expect_identical(res$status$code, "ok")
    expect_equal(res$solved$ratio, ratio, tolerance = 1e-8)
    expect_equal(res$diagnostics$failed_evaluations, 0)
  }
})
