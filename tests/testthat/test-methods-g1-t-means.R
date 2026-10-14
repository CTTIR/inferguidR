# t_means_one_sample and t_means_matched_pairs (group G1) ------------------------------------
#
# Evidence kinds, kept apart (VALIDATION_CONTRACT):
# * independent checks against the definition-based references of
#   helper-reference.R (noncentral t by integrating pnorm over the
#   chi-square density, central t quantiles by root finding on it; decision
#   D18), on selected points (the grid runs in validation/independent/g1/);
# * documented examples of the G*Power manual (documented_example evidence,
#   not application captures);
# * property, inverse/forward and typed-outcome tests.
# The G*Power emulator-capture comparisons are in test-methods-g1-fixtures.R.
#
# Tolerances: probabilities abs 1e-8 + rel 1e-6 (validation/tolerances.md
# section 5); quantiles 1e-8; integers exact; printed manual values by the
# closed rounding interval of the last printed digit (MR-DISPLAY).

ref_reset_registry()

g1_t_spec <- function(method, mode, tails = "two", ...) {
  analysis_spec(method, mode, inputs = list(...), options = list(tails = tails))
}

g1_expect_prob <- function(ours, ref, info = NULL) {
  expect_true(ref_close(ours, ref, abs_tol = 1e-8, rel_tol = 1e-6),
              info = paste(info, sprintf("ours %.15g ref %.15g", ours, ref)))
}

g1_display <- function(ours, shown) {
  x <- as.numeric(shown)
  dec <- if (grepl(".", shown, fixed = TRUE)) nchar(sub("^[^.]*\\.", "", shown)) else 0
  half <- 0.5 * 10^-dec
  is.finite(ours) && abs(ours - x) <= half + 8 * .Machine$double.eps * max(1, abs(x))
}

# Reference power of the one-sample t test (df = N - 1, ncp = effect sqrt(N)).
g1_ref_power_1 <- function(effect, n, alpha, tails) {
  df <- n - 1
  ncp <- effect * sqrt(n)
  if (tails == "two") {
    tc <- ref_qt_upper(alpha / 2, df)
    list(power = ref_pt(tc, df, ncp, lower_tail = FALSE) + ref_pt(-tc, df, ncp), critical = c(-tc, tc))
  } else {
    tc <- ref_qt_upper(alpha, df)
    if (effect >= 0) list(power = ref_pt(tc, df, ncp, lower_tail = FALSE), critical = tc) else
      list(power = ref_pt(-tc, df, ncp), critical = -tc)
  }
}

g1_numbers <- function(res) {
  x <- unclass(res)
  .ig_finite_paths(x[setdiff(names(x), c("spec", "diagnostics"))], "")
}

test_that("the one-sample and matched-pairs methods are registered with honest states", {
  for (id in c("t_means_one_sample", "t_means_matched_pairs")) {
    lm <- list_methods(family = "t")
    row <- lm[lm$id == id, ]
    expect_identical(row$implementation_status, "implemented", info = id)
    expect_identical(row$validation_status, "independently_checked", info = id)
    vs <- validation_status(id)
    expect_identical(vs$gpower_reference$status, "reference_pending", info = id)
    expect_identical(vs$independent_check$checker, "validation/independent/g1/check_g1.R")
    expect_identical(vs$modes_with_solver, .ig_modes)
    def <- .ig_get_method(id)
    meta <- .ig_get_metadata(id)
    expect_identical(def$modes, meta$modes)
    # Options exactly as the certified inventory lists them (codes and default).
    expect_identical(def$options$tails$values, meta$options$tails$codes)
    expect_identical(def$options$tails$default, meta$options$tails$default)
  }
})

test_that("specifications: defaults, aliases, the N domain and the power/alpha constraint", {
  sp <- analysis_spec("t_means_one_sample", "post_hoc", inputs = list(d = 0.5, alpha = 0.05, total_sample_size = 34))
  expect_identical(sp$options$tails, "one")
  expect_identical(sp$defaults_applied$options, "tails")
  expect_identical(sp$inputs$n_total, 34L)
  sp2 <- analysis_spec("t_means_matched_pairs", "post_hoc", inputs = list(effect_size_dz = 0.5, alpha = 0.05,
                                                                         n_pairs = 34))
  expect_identical(sp2$inputs$dz, 0.5)
  code <- function(expr) tryCatch({
    force(expr)
    NA_character_
  }, inferguidR_error = function(e) e$status$code)
  # N = 1 has df = 0; N = 2 (df = 1) is admissible (G*Power's input fields
  # refuse it, its a priori search returns it; method notes).
  expect_identical(code(g1_t_spec("t_means_one_sample", "post_hoc", d = 0.5, alpha = 0.05, n_total = 1)),
                   "invalid_input")
  expect_true(is.na(code(g1_t_spec("t_means_one_sample", "post_hoc", d = 0.5, alpha = 0.05, n_total = 2))))
  expect_identical(code(g1_t_spec("t_means_matched_pairs", "a_priori", dz = 0.5, alpha = 0.05, power = 0.8,
                                  n_total = 10)), "invalid_input")
  cond <- tryCatch(g1_t_spec("t_means_one_sample", "a_priori", d = 0.5, alpha = 0.05, power = 0.04),
                   inferguidR_invalid_input = function(e) e)
  expect_match(cond$status$reason, "below alpha")
  # power = alpha is allowed (G*Power: 'greater than or equal to alpha').
  ok <- solve_design(g1_t_spec("t_means_one_sample", "a_priori", d = 0.5, alpha = 0.05, power = 0.05))
  expect_identical(ok$status$code, "ok")
  expect_identical(ok$solved$n_total, 2L)
})

test_that("post hoc power, df, ncp and critical values match the independent reference", {
  cases <- list(
    list(e = 0.5, n = 34, tails = "one"), list(e = 0.5, n = 34, tails = "two"),
    list(e = -0.25, n = 25, tails = "two"), list(e = -0.7, n = 12, tails = "one"),
    list(e = 1.2, n = 2, tails = "two"), list(e = 0.05, n = 5000, tails = "one"),
    list(e = 3, n = 300, tails = "two")
  )
  for (id in c("t_means_one_sample", "t_means_matched_pairs")) {
    eff <- if (id == "t_means_one_sample") "d" else "dz"
    for (cs in cases) {
      info <- sprintf("%s e=%g n=%g %s", id, cs$e, cs$n, cs$tails)
      args <- c(list(method = id, mode = "post_hoc", inputs = stats::setNames(list(cs$e, 0.05, cs$n),
                                                                             c(eff, "alpha", "n_total")),
                     options = list(tails = cs$tails)))
      res <- calculate_power(do.call(analysis_spec, args))
      ref <- g1_ref_power_1(cs$e, cs$n, 0.05, cs$tails)
      expect_identical(res$status$code, "ok", info = info)
      g1_expect_prob(res$achieved$power, ref$power, info)
      expect_equal(res$achieved$power + res$achieved$beta, 1, tolerance = 1e-14, info = info)
      expect_identical(res$quantities$df, cs$n - 1)
      expect_equal(res$quantities$ncp, cs$e * sqrt(cs$n), tolerance = 1e-15)
      expect_equal(res$quantities$critical, ref$critical, tolerance = 1e-8, info = info)
      expect_identical(res$quantities$n_total, cs$n)
    }
  }
})

test_that("documented examples DOC-120, DOC-121 (one sample) and DOC-104, DOC-106 (pairs)", {
  # documented_example evidence (printed manual, older build), not G*Power
  # validation; values as in validation/fixtures/documented_example/.
  r120 <- solve_design(g1_t_spec("t_means_one_sample", "a_priori", "one", d = 0.625, alpha = 0.05, power = 0.95))
  expect_identical(r120$solved$n_total, 30L)
  expect_true(g1_display(r120$quantities$ncp, "3.423266"))
  expect_true(g1_display(r120$quantities$critical, "1.699127"))
  expect_true(g1_display(r120$achieved$power, "0.955144"))
  r121 <- solve_design(g1_t_spec("t_means_one_sample", "a_priori", "two", d = 0.1, alpha = 0.01, power = 0.9))
  expect_identical(r121$solved$n_total, 1492L)
  expect_true(g1_display(r121$quantities$ncp, "3.862642"))
  expect_true(g1_display(r121$achieved$power, "0.900169"))
  r104 <- calculate_power(g1_t_spec("t_means_matched_pairs", "post_hoc", "two", dz = 0.421637, alpha = 0.05,
                                    n_total = 50))
  expect_true(g1_display(r104$achieved$power, "0.832114"))
  expect_true(g1_display(r104$quantities$critical[2], "2.009575"))
  r106 <- calculate_power(g1_t_spec("t_means_matched_pairs", "post_hoc", "two", dz = 0.2828427, alpha = 0.05,
                                    n_total = 50))
  expect_true(g1_display(r106$achieved$power, "0.500352"))
  # Drawer DOC-103 / DOC-119.
  expect_true(g1_display(effect_size("t_means_matched_pairs", "from_group_parameters", mean_1 = 0, mean_2 = 0.4,
                                     sd_1 = 1, sd_2 = 1, rho = 0.55)$magnitude, "0.421637"))
  expect_true(g1_display(effect_size("t_means_one_sample", "from_means", mean_h0 = 10, mean_h1 = 15, sd = 8)$value,
                         "0.625"))
})

test_that("a priori: smallest N from N = 2, predecessor fails by the independent reference", {
  for (cs in list(list(e = 0.625, power = 0.95, tails = "one"), list(e = 0.3, power = 0.8, tails = "two"),
                  list(e = 10, power = 0.8, tails = "one"), list(e = -0.4, power = 0.9, tails = "one"))) {
    for (id in c("t_means_one_sample", "t_means_matched_pairs")) {
      eff <- if (id == "t_means_one_sample") "d" else "dz"
      info <- sprintf("%s e=%g power=%g %s", id, cs$e, cs$power, cs$tails)
      res <- solve_design(do.call(analysis_spec, list(id, "a_priori",
                                                      inputs = stats::setNames(list(cs$e, 0.05, cs$power),
                                                                               c(eff, "alpha", "power")),
                                                      options = list(tails = cs$tails))))
      N <- res$solved$n_total
      expect_identical(res$status$code, "ok", info = info)
      expect_gte(g1_ref_power_1(cs$e, N, 0.05, cs$tails)$power, cs$power)
      if (N > 2) expect_lt(g1_ref_power_1(cs$e, N - 1, 0.05, cs$tails)$power, cs$power)
      expect_identical(res$diagnostics$n_lower, 2)
      expect_match(res$diagnostics$allocation_rule, "N is the design")
    }
  }
  # Large effects reach the target at N = 2 (df = 1), as G*Power's search does.
  big <- solve_design(g1_t_spec("t_means_one_sample", "a_priori", "one", d = 10, alpha = 0.05, power = 0.8))
  expect_identical(big$solved$n_total, 2L)
  expect_identical(big$quantities$df, 1)
})

test_that("inverse modes recover the forward model for both tails and both signs", {
  for (id in c("t_means_one_sample", "t_means_matched_pairs")) {
    eff <- if (id == "t_means_one_sample") "d" else "dz"
    for (tails in c("one", "two")) {
      for (e in c(0.45, -0.45)) {
        info <- sprintf("%s %s e=%g", id, tails, e)
        mk <- function(mode, ...) do.call(analysis_spec, list(id, mode, inputs = list(...), options = list(tails = tails)))
        fwd <- calculate_power(do.call(analysis_spec, list(id, "post_hoc",
                                                          inputs = stats::setNames(list(e, 0.05, 40), c(eff, "alpha", "n_total")),
                                                          options = list(tails = tails))))
        P <- fwd$achieved$power
        se <- calculate_power(mk("sensitivity", alpha = 0.05, power = P, n_total = 40))
        expect_identical(se$status$code, "ok", info = info)
        expect_equal(se$solved[[eff]], abs(e), tolerance = 1e-8, info = info)
        expect_equal(se$diagnostics$other_branch_roots, -abs(e), tolerance = 1e-8, info = info)
        crit <- calculate_power(do.call(analysis_spec, list(id, "criterion",
                                                           inputs = stats::setNames(list(e, P, 40), c(eff, "power", "n_total")),
                                                           options = list(tails = tails))))
        expect_equal(crit$solved$alpha, 0.05, tolerance = 1e-8, info = info)
        # The critical value has the sign of the effect (one tail) on the
        # original scale (G*Power displays the mirrored value in criterion mode).
        if (tails == "one") expect_identical(sign(crit$quantities$critical), sign(e))
        q <- fwd$achieved$beta / 0.05
        co <- calculate_power(do.call(analysis_spec, list(id, "compromise",
                                                         inputs = stats::setNames(list(e, q, 40), c(eff, "error_ratio", "n_total")),
                                                         options = list(tails = tails))))
        expect_equal(co$solved$alpha, 0.05, tolerance = 1e-8, info = info)
        expect_equal(co$achieved$beta, q * co$solved$alpha, tolerance = 1e-9, info = info)
      }
    }
  }
})

test_that("sign reversal leaves power unchanged and mirrors the one-sided critical value", {
  for (tails in c("one", "two")) {
    a <- calculate_power(g1_t_spec("t_means_matched_pairs", "post_hoc", tails, dz = 0.37, alpha = 0.05, n_total = 21))
    b <- calculate_power(g1_t_spec("t_means_matched_pairs", "post_hoc", tails, dz = -0.37, alpha = 0.05, n_total = 21))
    expect_identical(b$achieved$power, a$achieved$power)
    expect_identical(b$achieved$beta, a$achieved$beta)
    expect_identical(b$quantities$critical, if (tails == "one") -a$quantities$critical else a$quantities$critical)
    expect_identical(b$method_outputs$statistic_distribution$rejection_side, if (tails == "one") "lower" else "both")
  }
  # A one-sample test equals the matched-pairs test with the same effect.
  one <- calculate_power(g1_t_spec("t_means_one_sample", "post_hoc", "two", d = 0.37, alpha = 0.05, n_total = 21))
  pair <- calculate_power(g1_t_spec("t_means_matched_pairs", "post_hoc", "two", dz = 0.37, alpha = 0.05, n_total = 21))
  expect_identical(one$achieved$power, pair$achieved$power)
})

test_that("power near one never exceeds one; large designs stay accurate", {
  # nmath's lower tail for t = 3.09, df = 199998, ncp = 22.36 is -3.8e-11;
  # the complement rule keeps the power inside [0, 1].
  p <- .ig_t_power_g1(199998, -22.36, 0.001, "one")
  expect_lte(p$power, 1)
  expect_identical(p$power + p$beta, 1)
  res <- calculate_power(g1_t_spec("t_means_one_sample", "sensitivity", "one", alpha = 0.001, power = 0.99,
                                   n_total = 100000))
  expect_identical(res$status$code, "ok")
  g1_expect_prob(ref_pt_gl(res$quantities$critical, 99999, res$quantities$ncp, lower_tail = FALSE), 0.99,
                 "sensitivity at N = 1e5")
})

test_that("unreachable targets and degenerate requests are typed outcomes without numbers", {
  flat <- solve_design(g1_t_spec("t_means_one_sample", "a_priori", "two", d = 0, alpha = 0.05, power = 0.8))
  expect_identical(flat$status$code, "infeasible")
  expect_match(flat$status$reason, "power does not change with N")
  expect_identical(g1_numbers(flat), character())
  touch <- calculate_power(g1_t_spec("t_means_matched_pairs", "sensitivity", "one", alpha = 0.05, power = 0.05,
                                     n_total = 20))
  expect_identical(touch$status$code, "infeasible")
  expect_identical(g1_numbers(touch), character())
})

test_that("effect-size drawers: formulas, signs and input checks", {
  es <- effect_size("t_means_one_sample", "from_means", mean_h0 = 15, mean_h1 = 10, sd = 8)
  expect_identical(es$value, -0.625)
  expect_identical(es$magnitude, 0.625)
  dz <- effect_size("t_means_matched_pairs", "from_differences", mean_diff = -5, sd_diff = 0.1)
  expect_identical(dz$value, -50)
  gp <- effect_size("t_means_matched_pairs", "from_group_parameters", mean_1 = 0, mean_2 = 5, sd_1 = 1, sd_2 = 3,
                    rho = 0.99)
  expect_equal(gp$value, -5 / sqrt(1 + 9 - 2 * 0.99 * 3), tolerance = 1e-15)
  expect_match(gp$source, "displays \\|dz\\|")
  bad <- function(...) tryCatch(effect_size(...), inferguidR_invalid_input = function(e) e)
  expect_s3_class(bad("t_means_matched_pairs", "from_group_parameters", mean_1 = 0, mean_2 = 1, sd_1 = 1, sd_2 = 1,
                      rho = 1), "inferguidR_invalid_input")
  expect_s3_class(bad("t_means_one_sample", "from_means", mean_h0 = 0, mean_h1 = 1, sd = 0), "inferguidR_invalid_input")
  expect_s3_class(bad("t_means_one_sample", "nope", mean_h0 = 0), "inferguidR_invalid_input")
  expect_s3_class(bad("t_means_matched_pairs", "from_differences", mean_diff = 1, sd_diff = 1, extra = 2),
                  "inferguidR_invalid_input")
})
