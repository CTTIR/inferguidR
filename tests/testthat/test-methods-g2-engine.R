# G2 engines (R/methods-anova.R): fixed-effects ANOVA (one-way, special effects), ANCOVA and
# Hotelling's T^2 (one group, two groups) through all five modes.
#
# Evidence kinds kept apart (VALIDATION_CONTRACT): property, boundary and inverse/forward tests
# here; the independent reference (noncentral F by integrating the convolution chi-square over the
# denominator density, helper-reference.R) in test-methods-g2-reference.R; G*Power fixture
# comparisons in test-methods-g2-fixtures.R.

ref_reset_registry()

g2_ids <- c("f_anova_oneway", "f_anova_special_main_interactions", "f_ancova_main_interactions",
            "f_hotelling_t2_one_group", "f_hotelling_t2_two_groups")

g2_calc <- function(id, mode, ...) calculate_power(analysis_spec(id, mode, inputs = list(...)))

test_that("the five G2 methods are registered with all modes and an honest validation state", {
  lm <- list_methods(family = "F")
  for (id in g2_ids) {
    row <- lm[lm$id == id, ]
    expect_identical(row$implementation_status, "implemented", info = id)
    expect_identical(row$validation_status, "independently_checked", info = id)
    expect_identical(row$modes, "a_priori, post_hoc, sensitivity, criterion, compromise", info = id)
    vs <- validation_status(id)
    expect_identical(vs$gpower_reference$status, "reference_pending", info = id)
    expect_identical(vs$independent_check$checker, "validation/independent/g2/check_g2.R", info = id)
    def <- .ig_get_method(id)
    expect_true(is.function(def$a_priori_search), info = id)
    expect_identical(def$effect_direction, "increasing", info = id)
  }
  expect_false(any(grepl("gpower_validated", unlist(lm))))
})

test_that("forward models map the effect, df and noncentrality as defined", {
  r <- g2_calc("f_anova_oneway", "post_hoc", f = 0.25, n_total = 120L, n_groups = 4L, alpha = 0.05)
  expect_equal(r$quantities$ncp, 0.25^2 * 120)
  expect_equal(r$quantities$df, c(3, 116))
  expect_equal(r$method_outputs$numerator_df, 3)
  expect_equal(unname(r$group_sizes), rep(30L, 4))
  r <- g2_calc("f_anova_oneway", "post_hoc", f = 0.25, n_total = 22L, n_groups = 4L, alpha = 0.05)
  expect_null(r$group_sizes)
  expect_true(any(grepl("average group size", r$assumptions)))
  r <- g2_calc("f_anova_special_main_interactions", "post_hoc", f = 0.3, n_total = 108L, numerator_df = 4L,
               n_groups = 36L, alpha = 0.05)
  expect_equal(r$quantities$df, c(4, 72))
  expect_equal(r$quantities$ncp, 0.09 * 108)
  r <- g2_calc("f_ancova_main_interactions", "post_hoc", f = 0.25, n_total = 120L, numerator_df = 2L,
               n_groups = 6L, n_covariates = 1L, alpha = 0.05)
  expect_equal(r$quantities$df, c(2, 113))
  r <- g2_calc("f_hotelling_t2_one_group", "post_hoc", delta = 0.5, n_total = 40L, n_response = 3L, alpha = 0.05)
  expect_equal(r$quantities$df, c(3, 37))
  expect_equal(r$quantities$ncp, 10)
  r <- g2_calc("f_hotelling_t2_two_groups", "post_hoc", delta = 0.5, n1 = 30L, n2 = 20L, n_response = 3L, alpha = 0.05)
  expect_equal(r$quantities$df, c(3, 46))
  expect_equal(r$quantities$ncp, 30 * 20 / 50 * 0.25)
  expect_equal(r$group_sizes, c(n1 = 30L, n2 = 20L))
  expect_identical(r$method_outputs$statistic_distribution$h1$family, "F")
  expect_identical(r$method_outputs$statistic_distribution$rejection_side, "upper")
})

test_that("the critical F avoids stats::qf's chi-square approximation for df2 > 4e5", {
  # stats::qf uses qchisq(p, df1) / df1 when df2 > 4e5 (nmath qf.c); the engine solves
  # P(F > c) = alpha through the beta quantile. Found with GPW-F-HOT2-APR-B001-D01.
  crit <- .ig_g2_crit_f(0.001, 50, 1444522)
  expect_equal(stats::pf(crit, 50, 1444522, lower.tail = FALSE), 0.001, tolerance = 1e-10)
  approx <- stats::qf(0.001, 50, 1444522, lower.tail = FALSE)
  expect_gt(abs(approx / crit - 1), 1e-5)
  for (a in c(0.05, 1e-6, 0.9)) for (d in list(c(1, 1), c(3, 116), c(99, 99900), c(2, 3e6), c(1, 2e9))) {
    c0 <- .ig_g2_crit_f(a, d[1], d[2])
    expect_equal(stats::pf(c0, d[1], d[2], lower.tail = FALSE, log.p = TRUE), log(a), tolerance = 1e-10,
                 info = paste(a, d[1], d[2]))
    if (d[2] <= 4e5) expect_equal(c0, stats::qf(a, d[1], d[2], lower.tail = FALSE), tolerance = 1e-10)
  }
  # extreme alpha (criterion and compromise search range)
  for (a in c(1e-300, 1 - 1e-12)) {
    c0 <- .ig_g2_crit_f(a, 3, 50)
    expect_true(is.finite(c0) && c0 > 0)
    expect_equal(stats::pf(c0, 3, 50, lower.tail = FALSE, log.p = TRUE), log(a), tolerance = 1e-9)
  }
})

test_that("power and beta are both computed directly; huge noncentrality is exact by the bound", {
  t1 <- .ig_g2_f_test(7.5, 2, 113, 0.05)
  expect_equal(t1$power + t1$beta, 1, tolerance = 1e-14)
  t0 <- .ig_g2_f_test(0, 3, 50, 0.01)
  expect_equal(t0$power, 0.01, tolerance = 1e-9)
  tb <- .ig_g2_f_test(400, 3, 100, 0.05)
  expect_false(tb$bounded)
  expect_true(tb$beta > 0 && tb$beta < 1e-30)
  big <- .ig_g2_f_test(1e13, 49, 99950, 1e-4)
  expect_true(big$bounded)
  expect_identical(c(big$power, big$beta), c(1, 0))
  # the bound never claims beta = 0 where it is representable
  near <- .ig_g2_f_test(2500, 2, 1, 0.05)
  expect_false(near$bounded)
  expect_gt(near$beta, 0)
})

test_that("invalid and inconsistent inputs are refused with typed outcomes", {
  code <- function(expr) th_condition_code(expr)
  expect_identical(code(analysis_spec("f_anova_oneway", "a_priori",
                                      inputs = list(f = 0, alpha = 0.05, power = 0.8, n_groups = 4L))), "invalid_input")
  expect_identical(code(analysis_spec("f_anova_oneway", "a_priori",
                                      inputs = list(f = 0.25, alpha = 0.05, power = 0.04, n_groups = 4L))), "invalid_input")
  expect_identical(code(analysis_spec("f_anova_oneway", "sensitivity",
                                      inputs = list(alpha = 0.05, power = 0.04, n_total = 120L, n_groups = 4L))), "invalid_input")
  expect_identical(code(analysis_spec("f_anova_oneway", "post_hoc",
                                      inputs = list(f = -0.25, alpha = 0.05, n_total = 120L, n_groups = 4L))), "invalid_input")
  expect_identical(code(analysis_spec("f_anova_oneway", "post_hoc",
                                      inputs = list(f = 0.25, alpha = 0.05, n_total = 4L, n_groups = 4L))), "invalid_input")
  expect_identical(code(analysis_spec("f_anova_special_main_interactions", "post_hoc",
                                      inputs = list(f = 0.5, alpha = 0.05, n_total = 6L, numerator_df = 2L, n_groups = 6L))), "invalid_input")
  expect_identical(code(analysis_spec("f_ancova_main_interactions", "post_hoc",
                                      inputs = list(f = 0.5, alpha = 0.05, n_total = 7L, numerator_df = 2L, n_groups = 6L,
                                                    n_covariates = 1L))), "invalid_input")
  expect_identical(code(analysis_spec("f_hotelling_t2_one_group", "post_hoc",
                                      inputs = list(delta = 1, alpha = 0.05, n_total = 2L, n_response = 2L))), "invalid_input")
  expect_identical(code(analysis_spec("f_hotelling_t2_two_groups", "post_hoc",
                                      inputs = list(delta = 1, alpha = 0.05, n1 = 2L, n2 = 1L, n_response = 2L))), "invalid_input")
  expect_identical(code(analysis_spec("f_hotelling_t2_two_groups", "a_priori",
                                      inputs = list(delta = 0, alpha = 0.05, power = 0.8, n_response = 3L))), "invalid_input")
  # f = 0 is admissible outside a priori mode: power = alpha
  r <- g2_calc("f_anova_oneway", "post_hoc", f = 0, n_total = 120L, n_groups = 4L, alpha = 0.05)
  expect_equal(r$achieved$power, 0.05, tolerance = 1e-9)
  # boundary designs (df2 = 1) are computed
  expect_identical(g2_calc("f_ancova_main_interactions", "post_hoc", f = 0.5, n_total = 8L, numerator_df = 2L,
                           n_groups = 6L, n_covariates = 1L, alpha = 0.05)$status$code, "ok")
  expect_identical(g2_calc("f_hotelling_t2_two_groups", "post_hoc", delta = 2.5, n1 = 2L, n2 = 2L, n_response = 2L,
                           alpha = 0.05)$status$code, "ok")
  # sensitivity with target power equal to alpha: touched but not crossed
  r <- g2_calc("f_anova_oneway", "sensitivity", alpha = 0.05, power = 0.05, n_total = 120L, n_groups = 4L)
  expect_identical(r$status$code, "infeasible")
})

# Smallest admissible design by plain enumeration with the forward model (no solver code).
g2_scan <- function(fw, p, from, step, target) {
  n <- from
  while (fw(c(p, list(n_total = n)))$power < target) n <- n + step
  n
}

test_that("one-way a priori designs are the smallest multiple of k with at least 2 per group", {
  for (cs in list(list(f = 0.25, k = 10L, power = 0.95), list(f = 1.5, k = 2L, power = 0.8),
                  list(f = 0.6, k = 7L, power = 0.8), list(f = 50, k = 3L, power = 0.8))) {
    r <- g2_calc("f_anova_oneway", "a_priori", f = cs$f, n_groups = cs$k, alpha = 0.05, power = cs$power)
    expect_identical(r$status$code, "ok")
    n <- r$solved$n_total
    expect_equal(n %% cs$k, 0)
    p <- list(f = cs$f, n_groups = cs$k, alpha = 0.05)
    expect_equal(n, g2_scan(.ig_power_f_anova_oneway, p, 2 * cs$k, cs$k, cs$power))
    expect_true(all(r$group_sizes == n / cs$k))
    expect_identical(r$diagnostics$step, as.double(cs$k))
  }
  # N = 386 reaches the target for f 0.25, k 10, but is not a multiple of k (DOC-055)
  p386 <- .ig_power_f_anova_oneway(list(f = 0.25, n_groups = 10, n_total = 386, alpha = 0.05))$power
  expect_gt(p386, 0.95)
})

test_that("single-N a priori designs are minimal (special effects, ANCOVA, Hotelling one group)", {
  r <- g2_calc("f_anova_special_main_interactions", "a_priori", f = 0.1, numerator_df = 8L, n_groups = 30L,
               alpha = 0.05, power = 0.95)
  expect_equal(r$solved$n_total, g2_scan(.ig_power_f_anova_special, list(f = 0.1, numerator_df = 8, n_groups = 30, alpha = 0.05),
                                         31, 1, 0.95))
  r <- g2_calc("f_anova_special_main_interactions", "a_priori", f = 50, numerator_df = 2L, n_groups = 3L,
               alpha = 0.05, power = 0.8)
  expect_equal(r$solved$n_total, 4)
  r <- g2_calc("f_ancova_main_interactions", "a_priori", f = 1.5, numerator_df = 1L, n_groups = 2L,
               n_covariates = 2L, alpha = 0.05, power = 0.8)
  expect_equal(r$solved$n_total, 8)
  expect_equal(r$quantities$df, c(1, 4))
  r <- g2_calc("f_hotelling_t2_one_group", "a_priori", delta = 0.5, n_response = 3L, alpha = 0.05, power = 0.8)
  expect_equal(r$solved$n_total, g2_scan(.ig_power_f_hotelling_one, list(delta = 0.5, n_response = 3, alpha = 0.05), 4, 1, 0.8))
  r <- g2_calc("f_hotelling_t2_one_group", "a_priori", delta = 50, n_response = 2L, alpha = 0.05, power = 0.8)
  expect_equal(r$solved$n_total, 3)
})

test_that("Hotelling two-group a priori follows G*Power's split and never reports a design below the target", {
  fw <- function(n1, n2, d, p, a = 0.05) .ig_power_f_hotelling_two(list(n1 = n1, n2 = n2, delta = d, n_response = p, alpha = a))$power
  # r = 1: even N (179 would reach the target with 90/89)
  r <- g2_calc("f_hotelling_t2_two_groups", "a_priori", delta = 0.5, n_response = 3L, allocation_ratio = 1,
               alpha = 0.05, power = 0.8)
  expect_equal(r$group_sizes, c(n1 = 90L, n2 = 90L))
  expect_gt(fw(90, 89, 0.5, 3), 0.8)
  # r != 1: any N, n1 = round(N / (1 + r)) half up, n2 = N - n1 (tie 36.5 / 109.5 -> 37 / 109)
  r <- g2_calc("f_hotelling_t2_two_groups", "a_priori", delta = 0.6, n_response = 2L, allocation_ratio = 3,
               alpha = 0.05, power = 0.8)
  expect_equal(r$group_sizes, c(n1 = 37L, n2 = 109L))
  expect_identical(r$diagnostics$gpower_rule$n_total, 146)
  # G*Power's design 52 / 157 misses the target; the engine reports the next design of the same split
  r <- g2_calc("f_hotelling_t2_two_groups", "a_priori", delta = 0.5, n_response = 2L, allocation_ratio = 3,
               alpha = 0.05, power = 0.8)
  expect_identical(c(r$diagnostics$gpower_rule$n1, r$diagnostics$gpower_rule$n2), c(52, 157))
  expect_lt(r$diagnostics$gpower_rule$power_integer_design, 0.8)
  expect_equal(r$group_sizes, c(n1 = 53L, n2 = 157L))
  expect_gte(r$achieved$power, 0.8)
  expect_match(r$diagnostics$design_choice, "misses the target")
  # an empty group in G*Power's split is never reported
  r <- g2_calc("f_hotelling_t2_two_groups", "a_priori", delta = 10, n_response = 2L, allocation_ratio = 100,
               alpha = 0.05, power = 0.8)
  expect_identical(r$diagnostics$gpower_rule$n1, 0)
  expect_true(all(r$group_sizes >= 1))
  expect_match(r$diagnostics$design_choice, "empty group")
  # the design rule reproduces every reported split
  for (rr in c(0.2, 1, 1.01, 3, 7)) for (n in c(5, 46, 147, 1001)) {
    p <- list(allocation_ratio = rr)
    g <- .ig_design_f_hotelling_two(p, n)
    if (!is.null(g)) expect_equal(.ig_g2_split_gpower(n, rr), g)
  }
  expect_null(.ig_design_f_hotelling_two(list(allocation_ratio = 1), 7))
})

test_that("inverse modes recover the forward model", {
  cases <- list(
    list(id = "f_anova_oneway", eff = "f", x = list(n_total = 120L, n_groups = 4L)),
    list(id = "f_anova_special_main_interactions", eff = "f", x = list(n_total = 120L, numerator_df = 2L, n_groups = 6L)),
    list(id = "f_ancova_main_interactions", eff = "f", x = list(n_total = 120L, numerator_df = 2L, n_groups = 6L, n_covariates = 1L)),
    list(id = "f_hotelling_t2_one_group", eff = "delta", x = list(n_total = 40L, n_response = 3L)),
    list(id = "f_hotelling_t2_two_groups", eff = "delta", x = list(n1 = 30L, n2 = 30L, n_response = 3L))
  )
  for (cs in cases) {
    s <- calculate_power(analysis_spec(cs$id, "sensitivity", inputs = c(cs$x, list(alpha = 0.05, power = 0.8))))
    expect_identical(s$status$code, "ok", info = cs$id)
    e <- s$solved[[cs$eff]]
    ph <- calculate_power(analysis_spec(cs$id, "post_hoc", inputs = c(cs$x, list(alpha = 0.05), stats::setNames(list(e), cs$eff))))
    expect_equal(ph$achieved$power, 0.8, tolerance = 1e-9, info = cs$id)
    cr <- calculate_power(analysis_spec(cs$id, "criterion", inputs = c(cs$x, list(power = 0.8), stats::setNames(list(0.3), cs$eff))))
    expect_identical(cr$status$code, "ok", info = cs$id)
    ph <- calculate_power(analysis_spec(cs$id, "post_hoc", inputs = c(cs$x, list(alpha = cr$solved$alpha), stats::setNames(list(0.3), cs$eff))))
    expect_equal(ph$achieved$power, 0.8, tolerance = 1e-9, info = cs$id)
    co <- calculate_power(analysis_spec(cs$id, "compromise", inputs = c(cs$x, list(error_ratio = 2), stats::setNames(list(0.3), cs$eff))))
    expect_identical(co$status$code, "ok", info = cs$id)
    expect_equal(co$achieved$beta / co$solved$alpha, 2, tolerance = 1e-8, info = cs$id)
  }
})

test_that("sensitivity scans over very large required effects stay fast and exact", {
  t0 <- proc.time()[["elapsed"]]
  r <- g2_calc("f_hotelling_t2_one_group", "sensitivity", alpha = 0.05, power = 0.95, n_total = 3L, n_response = 2L)
  expect_identical(r$status$code, "ok")
  expect_equal(r$achieved$power, 0.95, tolerance = 1e-9)
  r <- g2_calc("f_anova_oneway", "sensitivity", alpha = 0.001, power = 0.99, n_total = 100000L, n_groups = 50L)
  expect_identical(r$status$code, "ok")
  expect_lt(proc.time()[["elapsed"]] - t0, 20)
})

test_that("identities with the two-sided two-sample t test hold", {
  for (cs in list(c(25, 35, 0.6), c(10, 10, 1.2))) {
    h <- g2_calc("f_hotelling_t2_two_groups", "post_hoc", delta = cs[3], n1 = as.integer(cs[1]), n2 = as.integer(cs[2]),
                 n_response = 1L, alpha = 0.05)
    t2 <- calculate_power(analysis_spec("t_means_two_independent", "post_hoc",
                                        inputs = list(d = cs[3], n1 = as.integer(cs[1]), n2 = as.integer(cs[2]), alpha = 0.05),
                                        options = list(tails = "two")))
    expect_equal(h$achieved$power, t2$achieved$power, tolerance = 1e-9)
  }
  o <- g2_calc("f_anova_oneway", "post_hoc", f = 0.4, n_total = 64L, n_groups = 2L, alpha = 0.05)
  t2 <- calculate_power(analysis_spec("t_means_two_independent", "post_hoc",
                                      inputs = list(d = 0.8, n1 = 32L, n2 = 32L, alpha = 0.05), options = list(tails = "two")))
  expect_equal(o$achieved$power, t2$achieved$power, tolerance = 1e-9)
})
