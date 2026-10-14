# Wilcoxon procedures (group G10): definitions, forward model, modes, typed outcomes ------
#
# Methods: t_wilcoxon_signed_rank_one_sample, t_wilcoxon_signed_rank_matched_pairs,
# t_wilcoxon_mann_whitney (R/methods-wilcoxon.R, R/effect-wilcoxon.R; method
# notes specs/method-notes/<method_id>.md).
#
# Evidence kinds are kept apart (VALIDATION_CONTRACT):
# * independent checks: closed forms written here (normal, Laplace and
#   logistic difference CDFs; the Lehmann power recomputed from the published
#   formulas) on selected points; the full check with its own quadrature,
#   noncentral t reference and Monte Carlo is validation/independent/g10/;
# * documented examples of the G*Power manual (printed, not application
#   captures); where the replay capture of G*Power 3.1.9.7 differs from the
#   print, the capture is cited by its record id;
# * emulator captures (record ids cited) for the labelled differences and the
#   reference defects documented in validation/compare/<method_id>/summary.md;
# * property and inverse/forward tests.
#
# Tolerances: MR-DISPLAY (closed rounding interval of the printed digits plus
# 8 eps max(1, |x|)) for displayed values; 1e-12 for moments against closed
# forms; 1e-10 for the Lehmann power recomputed from the same moments; solver
# residuals 1e-8 (validation/tolerances.md sections 4 and 5).

ref_reset_registry()

g10_ids <- c("t_wilcoxon_signed_rank_one_sample", "t_wilcoxon_signed_rank_matched_pairs", "t_wilcoxon_mann_whitney")

g10_shown_ok <- function(ours, shown) {
  x <- as.numeric(shown)
  dec <- if (grepl(".", shown, fixed = TRUE)) nchar(sub("^[^.]*\\.", "", shown)) else 0
  is.finite(ours) && abs(ours - x) <= 0.5 * 10^-dec + 8 * .Machine$double.eps * max(1, abs(x))
}
g10_expect_shown <- function(ours, shown, info = "") {
  expect_true(g10_shown_ok(ours, shown), info = sprintf("%s: shown %s, ours %.10g", info, shown, ours))
}
g10_calc <- function(id, mode, inputs, options = list(), target = NULL) {
  calculate_power(analysis_spec(id, mode, inputs = inputs, options = options, target = target))
}
g10_post <- function(id, inputs, options = list()) g10_calc(id, "post_hoc", inputs, options)
lehmann <- function(...) c(list(method = "lehmann"), list(...))

# Closed forms used as independent references on selected points.
g10_laplace_g1 <- function(t) {
  b <- 1 / sqrt(2)
  u <- abs(t)
  up <- 1 - 0.5 * exp(-u / b) * (1 + u / (2 * b))
  if (t >= 0) up else 1 - up
}
g10_logistic_g1 <- function(t) {
  x <- abs(t) / (sqrt(3) / pi)
  cx <- if (x < 0.1) 0.5 + x / 6 - x^3 / 180 + x^5 / 5040 - x^7 / 151200 else (-expm1(-x) - x * exp(-x)) / expm1(-x)^2
  if (t >= 0) cx else 1 - cx
}

# Registration ------------------------------------------------------------------------

test_that("the three Wilcoxon methods are registered with all five modes and an honest status", {
  for (id in g10_ids) {
    def <- .ig_get_method(id)
    expect_s3_class(def, "ig_method")
    expect_identical(def$family, "t")
    expect_identical(def$modes, c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"))
    expect_identical(def$status, "independently_checked")
    expect_identical(def$checker, "validation/independent/g10/check_wilcoxon.R")
    expect_identical(.ig_implementation_status(id), "implemented")
    vs <- validation_status(id)
    expect_identical(vs$gpower_reference$status, "reference_pending")
    # Options and codes as in the certified registry; 'output_moments' is a
    # presentation setting and not part of a specification.
    meta <- .ig_get_metadata(id)
    reg <- setdiff(names(meta$options), "output_moments")
    expect_setequal(names(def$options), reg)
    for (o in reg) expect_setequal(def$options[[o]]$values, meta$options[[o]]$codes)
    for (o in reg) expect_identical(def$options[[o]]$default, meta$options[[o]]$default, info = paste(id, o))
    expect_identical(def$targets$sensitivity, c(.ig_wilcoxon_info(id)$effect, "p1"))
  }
})

test_that("the A.R.E. constants equal 12 sigma^2 (int f^2)^2 of the standardised parents", {
  are <- .ig_wilcoxon_are
  expect_equal(are[["normal"]], 12 * (1 / (2 * sqrt(pi)))^2, tolerance = 1e-15)
  b <- 1 / sqrt(2)
  expect_equal(are[["laplace"]], 12 * (1 / (4 * b))^2, tolerance = 1e-15)
  s <- sqrt(3) / pi
  expect_equal(are[["logistic"]], 12 * (1 / (6 * s))^2, tolerance = 1e-15)
  expect_identical(are[["min_are"]], 108 / 125)
})

# Moments and the Lehmann forward model ------------------------------------------------------

test_that("moments agree with closed forms", {
  for (D in c(0, 0.1, 0.5, 1.7)) {
    m <- .ig_wilcoxon_moments("signed_rank", "normal", D)
    expect_equal(m[["p1"]], pnorm(D), tolerance = 1e-12)
    expect_equal(m[["p2"]], pnorm(sqrt(2) * D), tolerance = 1e-12)
    w <- .ig_wilcoxon_moments("rank_sum", "normal", D)
    expect_equal(w[["p1"]], pnorm(D / sqrt(2)), tolerance = 1e-12)
    expect_identical(w[["p2"]], w[["p3"]])
    for (parent in c("laplace", "logistic")) {
      sr <- .ig_wilcoxon_moments("signed_rank", parent, D)
      rs <- .ig_wilcoxon_moments("rank_sum", parent, D)
      g1 <- if (parent == "laplace") g10_laplace_g1 else g10_logistic_g1
      expect_lt(abs(rs[["p1"]] - g1(D)), 1e-12)
      expect_lt(abs(sr[["p2"]] - g1(2 * D)), 1e-12)
      # identity of the two kinds: p2 of the signed-rank test at D is p1 of
      # the rank-sum test at 2 D, p3 likewise
      rs2 <- .ig_wilcoxon_moments("rank_sum", parent, 2 * D)
      expect_lt(abs(sr[["p3"]] - rs2[["p3"]]), 1e-12)
    }
  }
  # no effect: p1 = 1/2, p2 = 1/2, p3 = 1/3 (signed rank); p1 = 1/2, p2 = p3 = 1/3 (rank sum)
  expect_equal(unname(.ig_wilcoxon_moments("signed_rank", "logistic", 0)), c(0.5, 0.5, 1 / 3), tolerance = 1e-12)
  expect_equal(unname(.ig_wilcoxon_moments("rank_sum", "laplace", 0)), c(0.5, 1 / 3, 1 / 3), tolerance = 1e-12)
})

test_that("the Lehmann power equals the published formula recomputed from the moments", {
  N <- 30
  m <- .ig_wilcoxon_moments("signed_rank", "laplace", 0.4)
  e1 <- N * m[[1]] + choose(N, 2) * m[[2]]
  v1 <- N * m[[1]] * (1 - m[[1]]) + choose(N, 2) * (2 * (m[[1]] - m[[2]])^2 + 3 * m[[2]] * (1 - m[[2]])) +
    6 * choose(N, 3) * (m[[3]] - m[[2]]^2)
  c0 <- N * (N + 1) / 4 + qnorm(0.975) * sqrt(N * (N + 1) * (2 * N + 1) / 24)
  l0 <- N * (N + 1) / 2 - c0
  ref <- pnorm((c0 + 0.5 - e1) / sqrt(v1), lower.tail = FALSE) + pnorm((l0 - 0.5 - e1) / sqrt(v1))
  r <- g10_post(g10_ids[1], list(d = 0.4, alpha = 0.05, n_total = N),
                lehmann(parent_distribution = "laplace", tails = "two", continuity_correction = "on"))
  expect_lt(abs(r$achieved$power - ref), 1e-10)
  expect_lt(abs(r$achieved$power + r$achieved$beta - 1), 1e-14)
  # rank sum with unequal groups: p2 belongs to the group-2 size (supplement captures M01/M02)
  mm <- list(method = "lehmann", lehmann_input_method = "manual_moments")
  a <- g10_post(g10_ids[3], list(p1 = 0.62, p2 = 0.47, p3 = 0.42, alpha = 0.05, n1 = 5, n2 = 40), mm)
  b <- g10_post(g10_ids[3], list(p1 = 0.62, p2 = 0.47, p3 = 0.42, alpha = 0.05, n1 = 40, n2 = 5), mm)
  g10_expect_shown(a$achieved$power, "0.2146992", "SUP-G10-T-WMW-POS-M01--r01")
  g10_expect_shown(b$achieved$power, "0.1386597", "SUP-G10-T-WMW-POS-M02--r01")
})

test_that("documented examples of the manual are reproduced", {
  r <- g10_post(g10_ids[1], list(d = 0.1, alpha = 0.05, n_total = 649), lehmann())   # DOC-122
  g10_expect_shown(r$achieved$power, "0.8000589", "DOC-122")
  for (x in list(c("moment_p1", "0.5398278"), c("moment_p2", "0.5562315"), c("moment_p3", "0.3913924"))) {
    g10_expect_shown(r$method_outputs[[x[1]]], x[2], paste("DOC-122", x[1]))
  }
  r <- g10_post(g10_ids[1], list(d = 0.8, alpha = 0.05, n_total = 11), lehmann(parent_distribution = "laplace"))  # DOC-123
  g10_expect_shown(r$achieved$power, "0.8308341", "DOC-123")
  r <- g10_post(g10_ids[2], list(dz = 1.13842, alpha = 0.05, n_total = 10),       # DOC-126
                lehmann(parent_distribution = "laplace", tails = "two", continuity_correction = "on"))
  g10_expect_shown(r$achieved$power, "0.8534958", "DOC-126")
  g10_expect_shown(r$method_outputs$moment_p3, "0.9122280", "DOC-126 p3")
  r <- g10_post(g10_ids[3], list(d = 0.375, alpha = 0.05, n1 = 67, n2 = 134),     # DOC-129
                lehmann(parent_distribution = "laplace", tails = "two", continuity_correction = "on"))
  g10_expect_shown(r$achieved$power, "0.8472103", "DOC-129")
  g10_expect_shown(r$method_outputs$moment_p1, "0.6277817", "DOC-129 p1")
  r <- g10_post(g10_ids[3], list(d = 0.375, alpha = 0.05, n1 = 67, n2 = 134),     # DOC-130
                lehmann(tails = "two", continuity_correction = "on"))
  g10_expect_shown(r$achieved$power, "0.68", "DOC-130")
  r <- g10_post(g10_ids[3], list(p1 = 0.6, relative_sd = 1.15, alpha = 0.05, n1 = 67, n2 = 134),  # DOC-131
                lehmann(parent_distribution = "logistic", effect_size_parameter = "p1", continuity_correction = "on"))
  g10_expect_shown(r$achieved$power, "0.6461276", "DOC-131")
  g10_expect_shown(r$method_outputs$effect_size_d, "0.2911838", "DOC-131 d")
  g10_expect_shown(r$method_outputs$moment_p2, "0.4251206", "DOC-131 p2")
  r <- g10_post(g10_ids[3], list(d = 0.2911838, alpha = 0.05, n1 = 67, n2 = 134),  # DOC-133
                lehmann(parent_distribution = "laplace", continuity_correction = "on"))
  g10_expect_shown(r$achieved$power, "0.7551564", "DOC-133")
  g10_expect_shown(r$method_outputs$moment_p1, "0.6005708", "DOC-133 p1")
  # DOC-124/DOC-125 (A.R.E.): the print differs from G*Power 3.1.9.7; the
  # replay captures GPW-T-WSR1-POS-B003-O01 and -B004-O01 are reproduced.
  r <- g10_post(g10_ids[1], list(d = 0.8, alpha = 0.05, n_total = 11), list(parent_distribution = "laplace"))
  g10_expect_shown(r$achieved$power, "0.9275986", "GPW-T-WSR1-POS-B003-O01 (DOC-124 prints 0.923)")
  g10_expect_shown(r$quantities$df, "15.5000000", "GPW-T-WSR1-POS-B003-O01 df")
  r <- g10_post(g10_ids[1], list(d = 0.8, alpha = 0.05, n_total = 11), list(parent_distribution = "min_are"))
  g10_expect_shown(r$achieved$power, "0.7319798", "GPW-T-WSR1-POS-B004-O01 (DOC-125 prints 0.7319)")
})

test_that("conventions found in the captures hold", {
  # Lehmann with a predefined parent uses |d| (GPW-T-WMW-POS-B009-E02: d = -0.5)
  pos <- g10_post(g10_ids[3], list(d = 0.5, alpha = 0.05, n1 = 5, n2 = 40), lehmann())
  neg <- g10_post(g10_ids[3], list(d = -0.5, alpha = 0.05, n1 = 5, n2 = 40), lehmann())
  expect_identical(pos$achieved$power, neg$achieved$power)
  expect_identical(neg$method_outputs$moment_p1, pos$method_outputs$moment_p1)
  # A.R.E. one tail: rejection on the side of d (GPW-T-WMW-POS-B001-E01)
  r <- g10_post(g10_ids[3], list(d = -0.5, alpha = 0.05, n1 = 5, n2 = 40), list())
  g10_expect_shown(r$quantities$critical, "-1.6829047", "GPW-T-WMW-POS-B001-E01 critical t")
  g10_expect_shown(r$achieved$power, "0.2637667", "GPW-T-WMW-POS-B001-E01 power")
  g10_expect_shown(r$quantities$df, "40.9718346", "GPW-T-WMW-POS-B001-E01 df")
  # manual moments are used as given: p1 below 1/2 gives power below alpha
  # with one tail and counts the lower region with two (SUP-G10-T-WMW-POS-M07/M08)
  mm <- list(method = "lehmann", lehmann_input_method = "manual_moments")
  one <- g10_post(g10_ids[3], list(p1 = 0.45, p2 = 0.25, p3 = 0.22, alpha = 0.05, n1 = 50, n2 = 50), mm)
  two <- g10_post(g10_ids[3], list(p1 = 0.45, p2 = 0.25, p3 = 0.22, alpha = 0.05, n1 = 50, n2 = 50), c(mm, tails = "two"))
  g10_expect_shown(one$achieved$power, "4.335365e-05", "M07")
  expect_true(g10_shown_ok(one$achieved$power * 1e5, "4.335365"))
  g10_expect_shown(two$achieved$power, "0.0427456", "M08")
  # continuity correction lowers the power; the rank-sum critical values move by 0.5
  a <- g10_post(g10_ids[3], list(d = 0.5, alpha = 0.05, n1 = 50, n2 = 50), lehmann(distribution_plot = "normal_rank_sum"))
  b <- g10_post(g10_ids[3], list(d = 0.5, alpha = 0.05, n1 = 50, n2 = 50),
                lehmann(distribution_plot = "normal_rank_sum", continuity_correction = "on"))
  expect_lt(b$achieved$power, a$achieved$power)
  expect_equal(b$quantities$critical - a$quantities$critical, 0.5, tolerance = 1e-12)
  g10_expect_shown(a$quantities$critical, "1488.598", "inventory variant: lower = upper critical T, one tail")
  # the distribution plot changes only the reported critical values
  z <- g10_post(g10_ids[3], list(d = 0.5, alpha = 0.05, n1 = 50, n2 = 50), lehmann())
  expect_identical(z$achieved$power, a$achieved$power)
  g10_expect_shown(z$quantities$critical, "1.6448536", "critical z")
  # p1 effect-size parameter: d = shift(p1) / relative SD; the sign is kept
  r <- g10_post(g10_ids[1], list(p1 = 0.3, relative_sd = 2, alpha = 0.05, n_total = 40),
                lehmann(effect_size_parameter = "p1"))
  expect_equal(r$method_outputs$effect_size_d, qnorm(0.3) / 2, tolerance = 1e-13)
  expect_equal(r$method_outputs$effect_size_used, abs(qnorm(0.3) / 2), tolerance = 1e-13)
})

# Modes ------------------------------------------------------------------------------------

test_that("every mode solves consistently for each power method", {
  for (id in g10_ids) {
    eff <- .ig_wilcoxon_info(id)$effect
    two <- id == "t_wilcoxon_mann_whitney"
    nfix <- if (two) list(n1 = 30, n2 = 45) else list(n_total = 40)
    for (opt in list(list(parent_distribution = "logistic", tails = "two"),
                     lehmann(parent_distribution = "laplace", continuity_correction = "on"))) {
      e <- stats::setNames(list(0.45), eff)
      ap <- g10_calc(id, "a_priori", c(e, list(alpha = 0.05, power = 0.8), if (two) list(allocation_ratio = 1.5)), opt)
      expect_identical(ap$status$code, "ok")
      expect_gte(ap$achieved$power, 0.8)
      expect_lt(ap$diagnostics$predecessor_power, 0.8)
      if (two) expect_equal(unname(ap$group_sizes[["n1"]]), floor(ap$solved$n_total / 2.5 + 0.5))
      se <- g10_calc(id, "sensitivity", c(list(alpha = 0.05, power = 0.9), nfix), opt)
      expect_identical(se$status$code, "ok")
      back <- g10_post(id, c(stats::setNames(list(se$solved[[eff]]), eff), list(alpha = 0.05), nfix), opt)
      expect_lt(abs(back$achieved$power - 0.9), 1e-8)
      cr <- g10_calc(id, "criterion", c(e, list(power = 0.85), nfix), opt)
      expect_identical(cr$status$code, "ok")
      back <- g10_post(id, c(e, list(alpha = cr$solved$alpha), nfix), opt)
      expect_lt(abs(back$achieved$power - 0.85), 1e-8)
      co <- g10_calc(id, "compromise", c(e, list(error_ratio = 2), nfix), opt)
      expect_identical(co$status$code, "ok")
      expect_lt(abs(co$achieved$beta / (2 * co$solved$alpha) - 1), 1e-6)
    }
  }
})

test_that("the p1 branch solves p1 in the sensitivity mode, self-consistently", {
  id <- g10_ids[3]
  opt <- lehmann(effect_size_parameter = "p1", parent_distribution = "normal")
  se <- g10_calc(id, "sensitivity", list(relative_sd = 1, alpha = 0.05, power = 0.8, n1 = 50, n2 = 50), opt, target = "p1")
  expect_identical(se$status$code, "ok")
  # G*Power solves the same shift (GPW-T-WMW-SEN-B011-O01 displays p2 0.4867778 and the d branch
  # GPW-T-WMW-SEN-B010-O01 gives d 0.5135913) ...
  g10_expect_shown(se$method_outputs$moment_p2, "0.4867778", "GPW-T-WMW-SEN-B011-O01 p2")
  expect_lt(abs(se$method_outputs$effect_size_d - 0.51359136), 1e-7)
  # ... but displays p1 0.6962311, whose post hoc power is 0.9732339
  # (SUP-G10-T-WMW-POS-P01); the engine's solved p1 reproduces the target.
  back <- g10_post(id, list(p1 = se$solved$p1, relative_sd = 1, alpha = 0.05, n1 = 50, n2 = 50), opt)
  expect_lt(abs(back$achieved$power - 0.8), 1e-8)
  shown <- g10_post(id, list(p1 = 0.6962311, relative_sd = 1, alpha = 0.05, n1 = 50, n2 = 50), opt)
  g10_expect_shown(shown$achieved$power, "0.9732339", "SUP-G10-T-WMW-POS-P01--r01")
})

test_that("the a priori design is the smallest admissible design reaching the target (labelled difference)", {
  # GPW-T-WMW-APR-B002-E01 (A.R.E. logistic, d 1.2, one tail, ratio 3): G*Power
  # chooses N = 24 by the power of fractional group sizes; N = 22 (6/16)
  # already reaches the target with integer groups.
  r <- g10_calc(g10_ids[3], "a_priori", list(d = 1.2, alpha = 0.05, power = 0.8, allocation_ratio = 3),
                list(parent_distribution = "logistic"))
  expect_equal(r$solved$n_total, 22)
  expect_equal(unname(r$group_sizes), c(6, 16))
  expect_gte(r$achieved$power, 0.8)
  expect_lt(r$diagnostics$predecessor_power, 0.8)
  # with equal allocation the two rules coincide (inventory variant mainTail(s)=Two: 110 / 110)
  r <- g10_calc(g10_ids[3], "a_priori", list(d = 0.5, alpha = 0.05, power = 0.95), list(tails = "two"))
  expect_equal(unname(r$group_sizes), c(110, 110))
  g10_expect_shown(r$quantities$df, "208.0845", "inventory variant df")
})

# Typed outcomes ---------------------------------------------------------------------------

test_that("inadmissible or inconsistent inputs are refused with a typed status", {
  inv <- function(expr) {
    e <- tryCatch(expr, inferguidR_error = function(e) e)
    if (inherits(e, "inferguidR_error")) e$status$code else e$status$code
  }
  # G*Power: "The minimal A.R.E. of the Wilcoxon test versus the t test is 0.864."
  expect_identical(inv(g10_post(g10_ids[1], list(are = 0.2, d = 0.5, alpha = 0.05, n_total = 50),
                                list(are_specification = "manual_are"))), "invalid_input")
  # Total sample size below 3 (SUP-G10-T-WSR1-POS-D02); n1 + n2 = 3 is admissible (D03)
  expect_identical(inv(g10_post(g10_ids[1], list(d = 0.5, alpha = 0.05, n_total = 2))), "invalid_input")
  expect_identical(inv(g10_post(g10_ids[3], list(d = 0.5, alpha = 0.05, n1 = 1, n2 = 1))), "invalid_input")
  r <- g10_post(g10_ids[3], list(d = 0.5, alpha = 0.05, n1 = 1, n2 = 2))
  g10_expect_shown(r$achieved$power, "0.0771681", "SUP-G10-T-WMW-POS-D03--r01")
  # power below alpha (a priori and sensitivity)
  expect_identical(inv(g10_calc(g10_ids[2], "a_priori", list(dz = 0.5, alpha = 0.05, power = 0.04))), "invalid_input")
  # min ARE belongs to the A.R.E. method
  expect_identical(inv(g10_post(g10_ids[3], list(d = 0.5, alpha = 0.05, n1 = 20, n2 = 20),
                                lehmann(parent_distribution = "min_are"))), "invalid_input")
  # no sensitivity with manual moments (G*Power refuses too)
  expect_identical(inv(g10_calc(g10_ids[3], "sensitivity", list(p1 = 0.6, p2 = 0.45, p3 = 0.4, alpha = 0.05,
                                                                 power = 0.8, n1 = 20, n2 = 20),
                                list(method = "lehmann", lehmann_input_method = "manual_moments"))), "invalid_input")
  # the sensitivity target must match the effect-size parameter
  expect_identical(inv(g10_calc(g10_ids[3], "sensitivity", list(relative_sd = 1, p1 = 0.6, alpha = 0.05, power = 0.8,
                                                                 n1 = 20, n2 = 20), lehmann(effect_size_parameter = "p1"))),
                   "invalid_input")
  expect_identical(inv(g10_calc(g10_ids[3], "sensitivity", list(d = 0.5, alpha = 0.05, power = 0.8, n1 = 20, n2 = 20),
                                lehmann(), target = "p1")), "invalid_input")
  # moments that no distribution has (G*Power's initial p1 0.5, p2 = p3 = 0.65)
  expect_identical(inv(g10_post(g10_ids[3], list(p1 = 0.5, p2 = 0.65, p3 = 0.65, alpha = 0.05, n1 = 50, n2 = 50),
                                list(method = "lehmann", lehmann_input_method = "manual_moments"))), "invalid_input")
  expect_identical(inv(g10_post(g10_ids[1], list(p1 = 0.7, p2 = 0.6, p3 = 0.65, alpha = 0.05, n_total = 30),
                                list(method = "lehmann", lehmann_input_method = "manual_moments"))), "invalid_input")
  # options that do not apply are refused
  expect_identical(inv(g10_post(g10_ids[3], list(d = 0.5, alpha = 0.05, n1 = 20, n2 = 20),
                                list(continuity_correction = "on"))), "invalid_input")
  expect_identical(inv(g10_post(g10_ids[3], list(p1 = 0.6, p2 = 0.45, p3 = 0.4, alpha = 0.05, n1 = 20, n2 = 20),
                                list(method = "lehmann", lehmann_input_method = "manual_moments",
                                     effect_size_parameter = "p1"))), "invalid_input")
})

test_that("unreachable targets are infeasible without numbers", {
  r <- g10_calc(g10_ids[3], "a_priori", list(d = 0, alpha = 0.05, power = 0.8), list())
  expect_identical(r$status$code, "infeasible")
  expect_null(r$solved$n_total)
  # GPW-T-WSR1-SEN-B032-E01: N = 5, two tails, continuity correction; the
  # critical rank sum exceeds its maximum 15 and the approximate power stays
  # below 0.23, so the target 0.95 is unattainable (G*Power displays d = 5.03).
  r <- g10_calc(g10_ids[1], "sensitivity", list(alpha = 0.05, power = 0.95, n_total = 5),
                lehmann(tails = "two", continuity_correction = "on"))
  expect_identical(r$status$code, "infeasible")
  expect_identical(.ig_finite_paths(unclass(r)[c("solved", "achieved")], ""), character())
  # a huge shift with N = 5: the statistic is degenerate at its maximum under H1
  r <- g10_post(g10_ids[1], list(d = 50, alpha = 0.05, n_total = 5), lehmann(tails = "two", continuity_correction = "on"))
  expect_identical(r$status$code, "ok")
  expect_identical(r$achieved$power, 0)
})

# Effect-size drawers -----------------------------------------------------------------------

test_that("the drawers follow the captured sign conventions", {
  expect_equal(effect_size(g10_ids[1], "means_sd", mean_h0 = 10, mean_h1 = 12, sd = 4)$value, 0.5)  # U049-UO01
  expect_equal(effect_size(g10_ids[2], "differences", mean_diff = -5, sd_diff = 0.1)$value, -50)    # U050-UE01
  es <- effect_size(g10_ids[2], "group_parameters", mean_1 = 6500, mean_2 = 5600, sd_1 = 1500, sd_2 = 1300, r = 0.85)
  g10_expect_shown(es$magnitude, "1.13842", "DOC-127")
  es <- effect_size(g10_ids[2], "group_parameters", mean_1 = 0, mean_2 = 5, sd_1 = 1, sd_2 = 3, r = 0.99)
  g10_expect_shown(es$magnitude, "2.481458", "U051-UE01")
  expect_lt(es$value, 0)
  g10_expect_shown(effect_size(g10_ids[3], "common_sd", mean_1 = 0.3, mean_2 = 0.27, sd = 0.08)$magnitude, "0.375", "DOC-128")
  g10_expect_shown(effect_size(g10_ids[3], "two_sds", mean_1 = 0, mean_2 = 10, sd_1 = 0.1, sd_2 = 10)$magnitude,
                   "1.414143", "U053-UE01")
  expect_error(effect_size(g10_ids[3], "common_sd", mean_1 = 0, mean_2 = 1, sd = 0), class = "inferguidR_invalid_input")
  expect_error(effect_size(g10_ids[2], "group_parameters", mean_1 = 0, mean_2 = 1, sd_1 = 1, sd_2 = 1, r = 1),
               class = "inferguidR_invalid_input")
  expect_error(effect_size(g10_ids[1], "common_sd", mean_1 = 0, mean_2 = 1, sd = 1), class = "inferguidR_invalid_input")
})

 test_that("p1 sensitivity selects its active target by default", {
  opt <- lehmann(effect_size_parameter = "p1", parent_distribution = "normal")
  inp <- list(relative_sd = 1, alpha = 0.05, power = 0.8, n1 = 50, n2 = 50)
  implicit <- analysis_spec(g10_ids[3], "sensitivity", inputs = inp, options = opt)
  explicit <- analysis_spec(g10_ids[3], "sensitivity", inputs = inp, options = opt, target = "p1")
  expect_identical(implicit, explicit)
})
