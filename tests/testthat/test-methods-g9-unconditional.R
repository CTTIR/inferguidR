# Unconditional exact tests of two proportions (group G9) -----------------------------------
#
# exact_proportions_unconditional and exact_proportions_unconditional_offset
# (R/methods-unconditional.R, R/effect-unconditional.R; method notes in
# specs/method-notes/). The slow, broad independent check is
# validation/independent/g9/check_unconditional.R; the full G*Power fixture
# comparison is validation/compare/<method_id>/. These tests are a fast subset.

g9_m1 <- "exact_proportions_unconditional"
g9_m2 <- "exact_proportions_unconditional_offset"
g9_s1 <- c("z_pooled", "z_pooled_cc", "z_unpooled", "z_unpooled_cc", "mantel_haenszel", "likelihood_ratio",
           "t_test")

# Test-local references (written from the 2 x 2 counts; no engine function).
g9_pmf <- function(n, p) exp(lchoose(n, 0:n) + (0:n) * log(p) + (n - 0:n) * log1p(-p))

g9_ref_stat <- function(stat, x1, x2, n1, n2) {
  N <- n1 + n2
  m <- x1 + x2
  cr <- x1 * n2 - x2 * n1
  fin <- function(num, den) {
    out <- num / den
    out[num == 0] <- 0
    i <- num != 0 & den == 0
    out[i] <- sign(num[i]) * Inf
    out
  }
  pooled <- sqrt(n1 * n2 * m * (N - m) / N)
  unp <- sqrt(n2^2 * x1 * (n1 - x1) / n1 + n1^2 * x2 * (n2 - x2) / n2)
  cc <- sign(cr) * pmax(abs(cr) - N / 2, 0)
  switch(stat,
         z_pooled = fin(cr, pooled),
         z_pooled_cc = fin(cc, pooled),
         z_unpooled = fin(cr, unp),
         z_unpooled_cc = fin(cc, unp),
         mantel_haenszel = sign(cr) * sqrt(ifelse(m == 0 | m == N, 0, (N - 1) * cr^2 / (n1 * n2 * m * (N - m)))),
         likelihood_ratio = {
           e <- function(o, ex) ifelse(o == 0, 0, o * log(o / ex))
           sign(cr) * sqrt(pmax(2 * (e(x1, n1 * m / N) + e(n1 - x1, n1 * (N - m) / N) + e(x2, n2 * m / N) +
                                       e(n2 - x2, n2 * (N - m) / N)), 0))
         },
         t_test = fin(cr / (n1 * n2), sqrt((x1 * (n1 - x1) / n1 + x2 * (n2 - x2) / n2) / (N - 2) * (1 / n1 + 1 / n2))))
}

g9_ref_power <- function(stat, n1, n2, p1, p2, alpha, tails) {
  g <- expand.grid(x1 = 0:n1, x2 = 0:n2)
  z <- g9_ref_stat(stat, g$x1, g$x2, n1, n2)
  a <- if (tails == "two") alpha / 2 else alpha
  c <- if (stat == "t_test") stats::qt(a, n1 + n2 - 2, lower.tail = FALSE) else stats::qnorm(a, lower.tail = FALSE)
  dir <- if (p1 >= p2) 1 else -1
  rej <- if (tails == "two") z >= c | z <= -c else if (dir > 0) z >= c else z <= -c
  P1 <- as.vector(outer(g9_pmf(n1, p1), g9_pmf(n2, p2)))
  P0 <- as.vector(outer(g9_pmf(n1, p2), g9_pmf(n2, p2)))
  c(power = sum(P1[rej]), beta = sum(P1[!rej]), alpha = sum(P0[rej]))
}

g9_display_ok <- function(ours, shown) {
  x <- as.numeric(shown)
  mant <- sub("[eE].*$", "", shown)
  ex <- if (grepl("[eE]", shown)) as.numeric(sub("^.*[eE]", "", shown)) else 0
  dec <- if (grepl(".", mant, fixed = TRUE)) nchar(sub("^[^.]*\\.", "", mant)) else 0
  half <- 0.5 * 10^(ex - dec)
  allowance <- 8 * .Machine$double.eps * max(1, abs(x))
  is.finite(ours) && ours >= x - half - allowance && ours <= x + half + allowance
}

g9_spec <- function(method, mode, inputs, options = list(), target = NULL) {
  analysis_spec(method, mode, inputs = inputs, options = options, target = target)
}

# Registration ------------------------------------------------------------------------------

test_that("both methods are registered with the certified options and an honest status", {
  lm <- list_methods(family = "exact")
  for (id in c(g9_m1, g9_m2)) {
    row <- lm[lm$id == id, ]
    expect_identical(row$implementation_status, "implemented", info = id)
    expect_identical(row$validation_status, "independently_checked", info = id)
    expect_identical(row$modes, "a_priori, post_hoc, sensitivity, criterion, compromise", info = id)
    vs <- validation_status(id)
    expect_identical(vs$gpower_reference$status, "reference_pending", info = id)
    expect_identical(vs$independent_check$checker, "validation/independent/g9/check_unconditional.R", info = id)
    def <- .ig_get_method(id)
    meta <- .ig_get_metadata(id)
    # Every registry option exists in the definition with the registry codes and default.
    for (o in names(meta$options)) {
      expect_identical(def$options[[o]]$values, meta$options[[o]]$codes, info = paste(id, o))
      expect_identical(def$options[[o]]$default, meta$options[[o]]$default, info = paste(id, o))
    }
    expect_setequal(names(def$options), names(meta$options))
  }
  expect_identical(.ig_get_method(g9_m1)$options$effect_size_index$default, "odds_ratio")
  expect_identical(.ig_get_method(g9_m2)$options$effect_size_index$default, "proportion_p1")
  expect_identical(.ig_get_method(g9_m1)$targets$sensitivity, "p1")
  expect_identical(.ig_get_method(g9_m2)$targets$sensitivity, "p1_h1")
})

# Statistics and regions ------------------------------------------------------------------------

test_that("the inequality statistics equal a count-based reference at every table", {
  for (d in list(c(4L, 6L), c(7L, 3L), c(1L, 5L), c(9L, 9L))) {
    g <- expand.grid(x1 = 0:d[1], x2 = 0:d[2])
    for (st in g9_s1) {
      if (st == "t_test" && sum(d) < 3) next
      eng <- .ig_unc_stat(st, g$x1, g$x2, d[1], d[2], list(scale = "difference", value = 0))
      ref <- g9_ref_stat(st, g$x1, g$x2, d[1], d[2])
      expect_identical(is.infinite(eng), is.infinite(ref), info = paste(st, d[1], d[2]))
      fin <- is.finite(ref)
      expect_lt(max(abs(eng[fin] - ref[fin])), 1e-12)
    }
  }
})

test_that("zero-variance tables follow the documented conventions", {
  h0 <- list(scale = "difference", value = 0)
  # p1hat = 1, p2hat = 0: infinitely extreme for the unpooled and t statistics.
  expect_identical(.ig_unc_stat("z_unpooled", 3, 0, 3, 4, h0), Inf)
  expect_identical(.ig_unc_stat("z_unpooled", 0, 4, 3, 4, h0), -Inf)
  expect_identical(.ig_unc_stat("t_test", 3, 0, 3, 4, h0), Inf)
  expect_identical(.ig_unc_stat("z_unpooled_cc", 3, 0, 3, 4, h0), Inf)
  # 0 / 0 (equal proportions with zero variance) is 0.
  for (st in g9_s1) {
    expect_identical(.ig_unc_stat(st, c(0, 3), c(0, 4), 3, 4, h0), c(0, 0), info = st)
  }
  # Offset: a non-zero numerator over a zero standard error.
  expect_identical(.ig_unc_stat("z_pooled", 0, 0, 3, 3, list(scale = "difference", value = 0.1)), -Inf)
})

test_that("the threshold search gives the region of the full enumeration", {
  set.seed(20261005)
  h0 <- list(scale = "difference", value = 0)
  for (i in 1:40) {
    st <- sample(g9_s1, 1)
    n1 <- sample(2:25, 1)
    n2 <- sample(2:25, 1)
    crit <- stats::runif(1, -0.5, 3)
    side <- sample(c("upper", "lower", "both"), 1)
    if (side == "both") crit <- abs(crit) + 0.1
    thr <- .ig_unc_region_threshold(st, n1, n2, h0, crit, side)
    full <- .ig_unc_region_full(st, n1, n2, h0, crit, side)
    expect_identical(thr$tables_rejected, full$tables_rejected, info = paste(st, n1, n2, crit, side))
    pr <- list(c(stats::runif(1, 0.05, 0.95), stats::runif(1, 0.05, 0.95)))
    expect_equal(.ig_unc_region_probs(thr, n1, n2, pr), .ig_unc_region_probs(full, n1, n2, pr), tolerance = 1e-13)
  }
})

test_that("exact power, beta and actual alpha equal the enumeration reference", {
  for (st in g9_s1) for (tl in c("one", "two")) {
    for (cs in list(c(5, 8, 0.65, 0.4, 0.05), c(12, 7, 0.25, 0.4, 0.01), c(3, 3, 0.9, 0.5, 0.05))) {
      ref <- g9_ref_power(st, cs[1], cs[2], cs[3], cs[4], cs[5], tl)
      res <- calculate_power(g9_spec(g9_m1, "post_hoc",
                                     list(p1 = cs[3], p2 = cs[4], alpha = cs[5], n1 = cs[1], n2 = cs[2]),
                                     list(statistic = st, tails = tl, effect_size_index = "proportion_p1")))
      info <- paste(st, tl, paste(cs, collapse = " "))
      expect_identical(res$status$code, "ok", info = info)
      expect_lt(abs(res$achieved$power - ref[["power"]]), 1e-13)
      expect_lt(abs(res$achieved$beta - ref[["beta"]]), 1e-13)
      expect_lt(abs(res$achieved$alpha_attained - ref[["alpha"]]), 1e-13)
    }
  }
})

test_that("restricted MLEs maximise the constrained likelihood", {
  ll <- function(p1, p2, x1, x2, n1, n2) x1 * log(p1) + (n1 - x1) * log1p(-p1) + x2 * log(p2) + (n2 - x2) * log1p(-p2)
  for (tab in list(c(3, 7, 10, 12), c(9, 2, 10, 8), c(5, 5, 20, 20))) {
    x1 <- tab[1]; x2 <- tab[2]; n1 <- tab[3]; n2 <- tab[4]
    r <- .ig_unc_rmle_diff(x1, x2, n1, n2, 0.15)
    opt <- stats::optimize(function(p) ll(p + 0.15, p, x1, x2, n1, n2), c(0, 0.85), maximum = TRUE, tol = 1e-12)
    expect_lt(abs(r[1, 2] - opt$maximum), 1e-6)
    expect_equal(r[1, 1] - r[1, 2], 0.15, tolerance = 1e-12)
    r <- .ig_unc_rmle_ratio(x1, x2, n1, n2, 1.3)
    opt <- stats::optimize(function(p) ll(1.3 * p, p, x1, x2, n1, n2), c(0, 1 / 1.3), maximum = TRUE, tol = 1e-12)
    expect_lt(abs(r[1, 2] - opt$maximum), 1e-6)
    r <- .ig_unc_rmle_or(x1, x2, n1, n2, 2.5)
    p1f <- function(p) 2.5 * p / (1 - p + 2.5 * p)
    opt <- stats::optimize(function(p) ll(p1f(p), p, x1, x2, n1, n2), c(0, 1), maximum = TRUE, tol = 1e-12)
    expect_lt(abs(r[1, 2] - opt$maximum), 1e-6)
  }
})

test_that("the Gart-Nam statistic solves z = S + g (S^2 - 1) / 6", {
  g <- expand.grid(x1 = 0:6, x2 = 0:6)
  for (h0 in list(list(scale = "difference", value = 0.1), list(scale = "ratio", value = 1.2))) {
    S <- .ig_unc_stat("score_gart_nam", g$x1, g$x2, 6, 6, h0)
    z <- .ig_unc_stat("score_farrington_manning", g$x1, g$x2, 6, 6, h0)
    r <- if (h0$scale == "difference") .ig_unc_rmle_diff(g$x1, g$x2, 6, 6, 0.1) else .ig_unc_rmle_ratio(g$x1, g$x2, 6, 6, 1.2)
    phi <- if (h0$scale == "ratio") 1.2 else 1
    v <- r[, 1] * (1 - r[, 1]) / 6 + phi^2 * r[, 2] * (1 - r[, 2]) / 6
    gam <- (r[, 1] * (1 - r[, 1]) * (1 - 2 * r[, 1]) / 36 - phi^3 * r[, 2] * (1 - r[, 2]) * (1 - 2 * r[, 2]) / 36) / v^1.5
    ok <- is.finite(S) & is.finite(z) & v > 0
    expect_lt(max(abs(S[ok] + gam[ok] * (S[ok]^2 - 1) / 6 - z[ok])), 1e-10)
  }
})

# Modes -------------------------------------------------------------------------------------

test_that("a priori returns the smallest design of the allocation sequence reaching the target", {
  sp <- g9_spec(g9_m1, "a_priori", list(p1 = 0.65, p2 = 0.4, alpha = 0.05, power = 0.8),
                list(effect_size_index = "proportion_p1"))
  res <- solve_design(sp)
  expect_identical(res$status$code, "ok")
  expect_identical(unname(res$group_sizes), c(46L, 46L))
  expect_match(res$diagnostics$minimality, "^proven")
  # Every earlier design fails (exact power is not monotone in N: 37/37 has
  # less power than 36/36).
  pw <- vapply(1:45, function(n) g9_ref_power("z_pooled", n, n, 0.65, 0.4, 0.05, "one")[["power"]], numeric(1))
  expect_true(all(pw < 0.8))
  expect_gt(pw[36], pw[37])
  # Unequal allocation: n2 = floor(r n1 + 0.5).
  res <- solve_design(g9_spec(g9_m1, "a_priori", list(p1 = 0.95, p2 = 0.05, alpha = 0.05, power = 0.8,
                                                     allocation_ratio = 4),
                              list(effect_size_index = "proportion_p1", statistic = "mantel_haenszel")))
  expect_identical(unname(res$group_sizes), c(2L, 8L))
  # No effect: infeasible without numbers.
  res <- solve_design(g9_spec(g9_m1, "a_priori", list(p1 = 0.4, p2 = 0.4, alpha = 0.05, power = 0.8),
                              list(effect_size_index = "proportion_p1")))
  expect_identical(res$status$code, "infeasible")
  expect_null(res$solved$n_total)
})

test_that("sensitivity solves p1 on the side of the effect direction", {
  for (dirn in c("p_h1_ge_p_h0", "p_h1_le_p_h0")) for (tl in c("one", "two")) {
    res <- calculate_power(g9_spec(g9_m1, "sensitivity", list(p2 = 0.4, alpha = 0.05, power = 0.8, n1 = 40, n2 = 40),
                                   list(effect_size_index = "proportion_p1", effect_direction = dirn, tails = tl)))
    expect_identical(res$status$code, "ok")
    p1 <- res$solved$p1
    expect_true(if (dirn == "p_h1_ge_p_h0") p1 > 0.4 else p1 < 0.4)
    expect_lt(abs(g9_ref_power("z_pooled", 40, 40, p1, 0.4, 0.05, tl)[["power"]] - 0.8), 1e-9)
    # the displayed transformations of the solved p1
    expect_equal(res$method_outputs$odds_ratio, (p1 / (1 - p1)) / (0.4 / 0.6))
  }
  # Unreachable target (n = 3 per group): infeasible, not a boundary value.
  res <- calculate_power(g9_spec(g9_m1, "sensitivity", list(p2 = 0.05, alpha = 0.05, power = 0.95, n1 = 3, n2 = 3),
                                 list(effect_size_index = "proportion_p1", effect_direction = "p_h1_le_p_h0")))
  expect_identical(res$status$code, "infeasible")
})

test_that("criterion returns the smallest nominal alpha whose rule reaches the target", {
  res <- calculate_power(g9_spec(g9_m1, "criterion", list(p1 = 0.65, p2 = 0.4, power = 0.8, n1 = 40, n2 = 40),
                                 list(effect_size_index = "proportion_p1")))
  expect_identical(res$status$code, "ok")
  a <- res$solved$alpha
  expect_gte(res$achieved$power, 0.8)
  below <- g9_ref_power("z_pooled", 40, 40, 0.65, 0.4, a * (1 - 1e-9), "one")[["power"]]
  expect_lt(below, 0.8)
  expect_equal(res$achieved$power, g9_ref_power("z_pooled", 40, 40, 0.65, 0.4, a, "one")[["power"]], tolerance = 1e-12)
})

test_that("compromise: beta = q alpha for the large-sample computation; exact is reference_pending", {
  res <- calculate_power(g9_spec(g9_m1, "compromise", list(p1 = 0.65, p2 = 0.4, error_ratio = 1, n1 = 40, n2 = 40),
                                 list(effect_size_index = "proportion_p1", computation = "large_sample")))
  expect_identical(res$status$code, "ok")
  expect_equal(res$achieved$beta, res$solved$alpha, tolerance = 1e-8)
  res <- calculate_power(g9_spec(g9_m1, "compromise", list(p1 = 0.65, p2 = 0.4, error_ratio = 1, n1 = 40, n2 = 40),
                                 list(effect_size_index = "proportion_p1")))
  expect_identical(res$status$code, "reference_pending")
  expect_null(res$solved$alpha)
})

test_that("unsupported branches and inconsistent inputs are typed outcomes without numbers", {
  # Statistic not offered for the index (G*Power's list depends on it).
  expect_identical(th_condition_code(g9_spec(g9_m2, "post_hoc",
    list(risk_ratio_h1 = 1.6, risk_ratio_h0 = 1.2, p2 = 0.4, alpha = 0.05, n1 = 40, n2 = 40),
    list(effect_size_index = "risk_ratio"))), "invalid_input")
  # p1 outside (0, 1) from a difference.
  expect_identical(th_condition_code(g9_spec(g9_m2, "post_hoc",
    list(difference_h1 = 0.55, difference_h0 = 0, p2 = 0.5, alpha = 0.05, n1 = 40, n2 = 40),
    list(effect_size_index = "difference"))), "invalid_input")
  # Target power not above alpha in a priori mode.
  expect_identical(th_condition_code(g9_spec(g9_m1, "a_priori", list(p1 = 0.65, p2 = 0.4, alpha = 0.05, power = 0.04),
                                             list(effect_size_index = "proportion_p1"))), "invalid_input")
  # G*Power's admissible range of the offset test: N >= 3.
  expect_identical(th_condition_code(g9_spec(g9_m2, "post_hoc",
    list(p1_h1 = 0.65, p1_h0 = 0.5, p2 = 0.4, alpha = 0.05, n1 = 1, n2 = 1))), "invalid_input")
  # Gart-Nam with the large-sample computation at a zero effect: G*Power's form
  # divides the effect by the statistic at the expected table and has no value.
  res <- calculate_power(g9_spec(g9_m2, "post_hoc", list(p1_h1 = 0.5, p1_h0 = 0.5, p2 = 0.4, alpha = 0.05, n1 = 40, n2 = 40),
                                 list(statistic = "score_gart_nam", computation = "large_sample")))
  expect_identical(res$status$code, "unsupported")
  expect_null(res$achieved$power)
})

test_that("the large-sample Gart-Nam power uses the statistic at the expected H1 table", {
  gn <- function(inputs, options = list()) {
    calculate_power(g9_spec(g9_m2, "post_hoc", inputs,
                            c(list(statistic = "score_gart_nam", computation = "large_sample"), options)))
  }
  # Test-local reference: S_E solves z = S + g (S^2 - 1) / 6 at the expected
  # table (uniroot on the increasing branch); power = Phi((d - c d / S_E) / sd1).
  ref <- function(scale, h, p1, p2, n1, n2, alpha, two) {
    phi <- if (scale == "ratio") h else 1
    d <- if (scale == "ratio") p1 - h * p2 else p1 - p2 - h
    r <- if (scale == "ratio") .ig_unc_rmle_ratio(n1 * p1, n2 * p2, n1, n2, h) else .ig_unc_rmle_diff(n1 * p1, n2 * p2, n1, n2, h)
    v <- r[1, 1] * (1 - r[1, 1]) / n1 + phi^2 * r[1, 2] * (1 - r[1, 2]) / n2
    g <- (r[1, 1] * (1 - r[1, 1]) * (1 - 2 * r[1, 1]) / n1^2 - phi^3 * r[1, 2] * (1 - r[1, 2]) * (1 - 2 * r[1, 2]) / n2^2) / v^1.5
    z <- d / sqrt(v)
    S <- stats::uniroot(function(s) s + g * (s^2 - 1) / 6 - z, c(z - 1, z + 1), tol = 1e-14)$root
    sd1 <- sqrt(p1 * (1 - p1) / n1 + phi^2 * p2 * (1 - p2) / n2)
    c <- stats::qnorm(if (two) alpha / 2 else alpha, lower.tail = FALSE)
    dir <- sign(d)
    stats::pnorm((dir * d - c * d / S) / sd1) + if (two) stats::pnorm((-dir * d - c * d / S) / sd1) else 0
  }
  res <- gn(list(p1_h1 = 0.9, p1_h0 = 0.5, p2 = 0.4, alpha = 0.01, n1 = 30, n2 = 60))
  expect_identical(res$status$code, "ok")
  expect_equal(res$achieved$power, ref("difference", 0.1, 0.9, 0.4, 30, 60, 0.01, FALSE), tolerance = 1e-10)
  expect_identical(res$achieved$alpha_attained, 0.01)
  res <- gn(list(p1_h1 = 0.35, p1_h0 = 0.5, p2 = 0.4, alpha = 0.05, n1 = 25, n2 = 40), list(tails = "two"))
  expect_equal(res$achieved$power, ref("difference", 0.1, 0.35, 0.4, 25, 40, 0.05, TRUE), tolerance = 1e-10)
  res <- gn(list(risk_ratio_h1 = 2, risk_ratio_h0 = 1.2, p2 = 0.4, alpha = 0.05, n1 = 100, n2 = 50),
            list(effect_size_index = "risk_ratio"))
  expect_equal(res$achieved$power, ref("ratio", 1.2, 0.8, 0.4, 100, 50, 0.05, FALSE), tolerance = 1e-10)
  # Without skewness (null difference 0, equal groups) the Gart-Nam statistic
  # is the Farrington-Manning statistic, and so is its large-sample power.
  a <- gn(list(p1_h1 = 0.55, p1_h0 = 0.4, p2 = 0.4, alpha = 0.05, n1 = 50, n2 = 50))
  b <- calculate_power(g9_spec(g9_m2, "post_hoc", list(p1_h1 = 0.55, p1_h0 = 0.4, p2 = 0.4, alpha = 0.05, n1 = 50, n2 = 50),
                               list(statistic = "score_farrington_manning", computation = "large_sample")))
  expect_equal(a$achieved$power, b$achieved$power, tolerance = 1e-12)
  # The other modes run through the generic solvers.
  res <- solve_design(g9_spec(g9_m2, "a_priori", list(p1_h1 = 0.65, p1_h0 = 0.5, p2 = 0.4, alpha = 0.05, power = 0.8),
                              list(statistic = "score_gart_nam", computation = "large_sample")))
  expect_identical(res$status$code, "ok")
  expect_gte(res$achieved$power, 0.8)
  res <- calculate_power(g9_spec(g9_m2, "sensitivity", list(p1_h0 = 0.5, p2 = 0.4, alpha = 0.05, power = 0.8, n1 = 40, n2 = 40),
                                 list(statistic = "score_gart_nam", computation = "large_sample")))
  expect_identical(res$status$code, "ok")
  expect_lt(abs(res$achieved$power - 0.8), 1e-9)
  expect_gt(res$solved$p1_h1, 0.5)
})

test_that("an exhaustive a priori search beyond the budget is unsupported, not unverified", {
  testthat::local_mocked_bindings(.ig_unc_budget = function() 1e5)
  res <- solve_design(g9_spec(g9_m2, "a_priori",
    list(difference_h1 = 0.05, difference_h0 = 0, p2 = 0.4, alpha = 0.001, power = 0.99, allocation_ratio = 0.25),
    list(effect_size_index = "difference", statistic = "z_unpooled", tails = "two")))
  expect_identical(res$status$code, "unsupported")
  expect_match(res$status$reason, "budget")
})

test_that("effect-size conversions", {
  es <- effect_size(g9_m1, "odds_ratio", p1 = 0.65, p2 = 0.5)
  expect_identical(es$status$code, "ok")
  expect_equal(es$value, (0.65 / 0.35))
  es <- effect_size(g9_m1, "p1_from_odds_ratio", odds_ratio = 1.857143, p2 = 0.5)
  expect_equal(es$value, 0.65, tolerance = 1e-6)
  es <- effect_size(g9_m2, "p1_from_risk_ratio", risk_ratio = 1.2, p2 = 0.4)
  expect_equal(es$value, 0.48)
  expect_identical(th_condition_code(effect_size(g9_m1, "p1_from_difference", difference = 0.7, p2 = 0.5)), "invalid_input")
})

test_that("results of every mode round-trip through JSON", {
  results <- list(
    calculate_power(g9_spec(g9_m1, "post_hoc", list(p1 = 0.65, p2 = 0.4, alpha = 0.05, n1 = 10, n2 = 12),
                            list(effect_size_index = "proportion_p1"))),
    solve_design(g9_spec(g9_m1, "a_priori", list(odds_ratio = 3, p2 = 0.4, alpha = 0.05, power = 0.8))),
    calculate_power(g9_spec(g9_m2, "sensitivity", list(p1_h0 = 0.5, p2 = 0.4, alpha = 0.05, power = 0.8, n1 = 20, n2 = 20),
                            list(statistic = "score_farrington_manning"))),
    calculate_power(g9_spec(g9_m2, "criterion", list(p1_h1 = 0.8, p1_h0 = 0.5, p2 = 0.4, power = 0.8, n1 = 20, n2 = 20))))
  for (res in results) {
    expect_identical(res$status$code, "ok", info = res$mode)
    expect_identical(ig_from_json(ig_to_json(res)), res, info = res$mode)
  }
})
