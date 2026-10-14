# Generic tests: t_generic, z_generic, f_generic, chi2_generic (group G1) -------------------
#
# Independent checks on selected points (decision D18; the grid runs in
# validation/independent/g1/):
# * noncentral t: ref_pt() / ref_qt_upper() of helper-reference.R;
# * normal: Phi from the incomplete gamma function, Q(x) = pgamma(x^2/2, 1/2,
#   upper) / 2 (nmath pgamma, not pnorm's Cody algorithm), quantiles by
#   bisection on it;
# * noncentral chi-square: ref_pchisq() (convolution (Z + sqrt(lambda))^2 +
#   chi2(k - 1), central pchisq only; the engine sums Poisson mixtures);
# * noncentral F: ref_pf() (the convolution integrated over the denominator
#   chi-square density; the engine sums Poisson mixtures of pbeta).
# The checks evaluate the reference at the engine's critical values and
# solutions (tail area alpha at the critical value, target power at the
# solved effect or alpha), so no reference root finding is needed for F and
# chi-square.
#
# Tolerances: probabilities abs 1e-8 + rel 1e-6 (tolerances.md section 5);
# small tail probabilities (< 1e-8) relative 1e-6.

ref_reset_registry()

g1_gen_prob <- function(ours, ref, info = NULL) {
  ok <- if (ref >= 1e-8) ref_close(ours, ref, 1e-8, 1e-6) else ours > 0 && abs(log(ours) - log(ref)) <= 1e-6
  expect_true(ok, info = paste(info, sprintf("ours %.15g ref %.15g", ours, ref)))
}

# Standard normal tails without pnorm(): Q(x) = P(Z > x).
g1_ref_q <- function(x) {
  if (x >= 0) 0.5 * stats::pgamma(x^2 / 2, 0.5, lower.tail = FALSE) else 1 - 0.5 * stats::pgamma(x^2 / 2, 0.5, lower.tail = FALSE)
}
g1_ref_p <- function(x) g1_ref_q(-x)
g1_ref_qnorm_upper <- function(p) {
  lo <- -40
  hi <- 40
  for (i in 1:200) {
    m <- (lo + hi) / 2
    if (g1_ref_q(m) > p) lo <- m else hi <- m
  }
  (lo + hi) / 2
}

g1_gen <- function(method, mode, tails = NULL, ...) {
  calculate_power(analysis_spec(method, mode, inputs = list(...),
                                options = if (is.null(tails)) list() else list(tails = tails)))
}

g1_gen_numbers <- function(res) {
  x <- unclass(res)
  .ig_finite_paths(x[setdiff(names(x), c("spec", "diagnostics"))], "")
}

test_that("the generic tests are registered with the inventory's modes and options", {
  for (id in c("t_generic", "z_generic", "f_generic", "chi2_generic")) {
    def <- .ig_get_method(id)
    meta <- .ig_get_metadata(id)
    expect_identical(def$modes, meta$modes, info = id)
    expect_identical(def$modes, c("post_hoc", "sensitivity", "criterion", "compromise"))
    expect_identical(names(def$options), names(meta$options), info = id)
    if (length(meta$options)) expect_identical(def$options$tails$default, meta$options$tails$default)
    expect_identical(def$status, "independently_checked")
    expect_identical(.ig_implementation_status(id), "implemented")
    es <- effect_size(id, "any", x = 1)
    expect_identical(es$status$code, "unsupported", info = id)
  }
  expect_s3_class(tryCatch(analysis_spec("t_generic", "a_priori", inputs = list(delta = 2, alpha = 0.05, power = 0.8)),
                           error = function(e) e), "inferguidR_invalid_input")
  # The H1 SD of the generic z test defaults to G*Power's initial value 1.
  sp <- analysis_spec("z_generic", "post_hoc", inputs = list(mu = 2.5, alpha = 0.05))
  expect_identical(sp$inputs$sigma, 1)
  expect_true("sigma" %in% sp$defaults_applied$inputs)
})

test_that("generic t: post hoc against the reference, non-integer df, both signs", {
  for (cs in list(list(delta = 2.5, df = 30, tails = "one"), list(delta = -1.25, df = 24, tails = "two"),
                  list(delta = -1, df = 1, tails = "one"), list(delta = 2.5, df = 0.5, tails = "one"),
                  list(delta = 40, df = 1e5, tails = "two"))) {
    info <- sprintf("delta=%g df=%g %s", cs$delta, cs$df, cs$tails)
    res <- g1_gen("t_generic", "post_hoc", cs$tails, delta = cs$delta, alpha = 0.05, df = cs$df)
    expect_identical(res$status$code, "ok", info = info)
    a <- if (cs$tails == "two") 0.025 else 0.05
    tc <- ref_qt_upper(a, cs$df)
    expect_equal(abs(res$quantities$critical), rep(tc, if (cs$tails == "two") 2 else 1), tolerance = 1e-8)
    ref_pw <- if (cs$tails == "two") {
      ref_pt(tc, cs$df, cs$delta, lower_tail = FALSE) + ref_pt(-tc, cs$df, cs$delta)
    } else if (cs$delta >= 0) ref_pt(tc, cs$df, cs$delta, lower_tail = FALSE) else ref_pt(-tc, cs$df, cs$delta)
    # beta directly, as a difference of same-side tails (no cancellation).
    ref_beta <- if (cs$tails == "two") {
      ref_pt(tc, cs$df, cs$delta) - ref_pt(-tc, cs$df, cs$delta)
    } else if (cs$delta >= 0) ref_pt(tc, cs$df, cs$delta) else ref_pt(-tc, cs$df, cs$delta, lower_tail = FALSE)
    g1_gen_prob(res$achieved$power, ref_pw, info)
    g1_gen_prob(res$achieved$beta, ref_beta, info)
  }
  # Ledger entry DOC-134 (manual 25.3; alpha not printed, so no fixture
  # record): delta = -1.25, df = 24, two tails. The manual prints 0.224525;
  # the reference and G*Power 3.1.9.7 (captures GPW-T-GENT-POS-B002-O01 and
  # GPW-T-MEAN1-POS-B002-O01 show 0.2245240) give 0.22452403, so the printed
  # last digit is off by one (manual defect, method notes).
  res <- g1_gen("t_generic", "post_hoc", "two", delta = -1.25, alpha = 0.05, df = 24)
  ref <- ref_pt(ref_qt_upper(0.025, 24), 24, -1.25, lower_tail = FALSE) + ref_pt(-ref_qt_upper(0.025, 24), 24, -1.25)
  g1_gen_prob(res$achieved$power, ref, "DOC-134")
  expect_lte(abs(res$achieved$power - 0.224524), 5e-7)
  expect_gt(abs(res$achieved$power - 0.224525), 5e-7)
})

test_that("generic z: power, beta and critical values from the definition", {
  for (cs in list(list(mu = 2.5, sigma = 1, tails = "one"), list(mu = -0.5, sigma = 0.2, tails = "one"),
                  list(mu = -0.5, sigma = 0.2, tails = "two"), list(mu = 0.001, sigma = 5, tails = "two"),
                  list(mu = 6, sigma = 0.2, tails = "one"))) {
    info <- sprintf("mu=%g sigma=%g %s", cs$mu, cs$sigma, cs$tails)
    res <- g1_gen("z_generic", "post_hoc", cs$tails, mu = cs$mu, sigma = cs$sigma, alpha = 0.05)
    expect_identical(res$status$code, "ok", info = info)
    if (cs$tails == "two") {
      zc <- g1_ref_qnorm_upper(0.025)
      expect_equal(res$quantities$critical, c(-zc, zc), tolerance = 1e-10)
      up <- g1_ref_q((zc - cs$mu) / cs$sigma)
      lo <- g1_ref_p((-zc - cs$mu) / cs$sigma)
      g1_gen_prob(res$achieved$power, up + lo, info)
    } else {
      zc <- g1_ref_qnorm_upper(0.05)
      side <- if (cs$mu >= 0) 1 else -1
      expect_equal(res$quantities$critical, side * zc, tolerance = 1e-10)
      g1_gen_prob(res$achieved$power, g1_ref_q((zc - abs(cs$mu)) / cs$sigma), info)
      g1_gen_prob(res$achieved$beta, g1_ref_p((zc - abs(cs$mu)) / cs$sigma), info)
    }
  }
  # With sigma < 1 the power at mu = 0 is below alpha (one tail:
  # Phi(-z / sigma)); a target between that minimum and alpha is solvable.
  zero <- g1_gen("z_generic", "post_hoc", "one", mu = 0, sigma = 0.5, alpha = 0.05)
  g1_gen_prob(zero$achieved$power, g1_ref_q(g1_ref_qnorm_upper(0.05) / 0.5), "mu = 0")
  low <- g1_gen("z_generic", "sensitivity", "one", sigma = 0.5, alpha = 0.05, power = 0.04)
  expect_identical(low$status$code, "ok")
  g1_gen_prob(g1_ref_q((g1_ref_qnorm_upper(0.05) - low$solved$mu) / 0.5), 0.04, "sensitivity below alpha")
})

test_that("generic F and chi-square: tail area alpha at the critical value, power by convolution", {
  f <- g1_gen("f_generic", "post_hoc", lambda = 12, alpha = 0.05, df1 = 3, df2 = 60)
  expect_identical(f$status$code, "ok")
  g1_gen_prob(ref_pf(f$quantities$critical, 3, 60, 0, lower_tail = FALSE), 0.05, "F critical")
  g1_gen_prob(f$achieved$power, ref_pf(f$quantities$critical, 3, 60, 12, lower_tail = FALSE), "F power")
  expect_identical(f$quantities$df, c(3, 60))
  ff <- g1_gen("f_generic", "post_hoc", lambda = 12, alpha = 0.05, df1 = 2.5, df2 = 7.5)
  g1_gen_prob(ff$achieved$power, ref_pf(ff$quantities$critical, 2.5, 7.5, 12, lower_tail = FALSE), "F df 2.5/7.5")
  for (cs in list(list(lambda = 10, df = 5), list(lambda = 0.5, df = 1), list(lambda = 400, df = 1000))) {
    info <- sprintf("chi2 lambda=%g df=%g", cs$lambda, cs$df)
    x <- g1_gen("chi2_generic", "post_hoc", lambda = cs$lambda, alpha = 0.05, df = cs$df)
    g1_gen_prob(ref_pchisq(x$quantities$critical, cs$df, 0, lower_tail = FALSE), 0.05, info)
    g1_gen_prob(x$achieved$power, ref_pchisq(x$quantities$critical, cs$df, cs$lambda, lower_tail = FALSE), info)
    g1_gen_prob(x$achieved$beta, ref_pchisq(x$quantities$critical, cs$df, cs$lambda, lower_tail = TRUE), info)
  }
  # lambda = 0: the power is alpha.
  z <- g1_gen("chi2_generic", "post_hoc", lambda = 0, alpha = 0.05, df = 5)
  expect_equal(z$achieved$power, 0.05, tolerance = 1e-14)
})

test_that("generic F: beta stays the complement of the power for a huge critical F", {
  # Found by the review cross-check (validation/independent/g1/review-mpmath):
  # with alpha 1e-5 and df2 = 0.7 the critical F is about 5e13, and the lower
  # tail of .ig_p_f() lost its accuracy (beta 0.99990775452 instead of
  # 0.99990800573, power + beta = 1 - 2.5e-7; df2 = 0.5: beta 1 instead of
  # 1 - 1.27e-4). The reference values are from the arbitrary-precision
  # mpmath reference of that check (no G*Power number).
  f <- g1_gen("f_generic", "post_hoc", lambda = 300, alpha = 1e-5, df1 = 1, df2 = 0.7)
  expect_identical(f$status$code, "ok")
  expect_equal(f$achieved$power + f$achieved$beta, 1, tolerance = 1e-15)
  expect_equal(f$achieved$power, 9.1994273610289155e-5, tolerance = 1e-10)
  expect_equal(f$achieved$beta, 0.99990800572638971, tolerance = 1e-13)
  for (cs in list(c(0.5, 0.5, 3000), c(20, 1, 3000), c(2, 0.5, 30))) {
    g <- g1_gen("f_generic", "post_hoc", lambda = cs[3], alpha = 1e-5, df1 = cs[1], df2 = cs[2])
    expect_equal(g$achieved$power + g$achieved$beta, 1, tolerance = 1e-15, info = paste(cs, collapse = "/"))
  }
})

test_that("inverse modes: sensitivity, criterion and compromise solutions reach their targets", {
  # Sensitivity solutions above lambda = 100, which G*Power 3.1.9.7 displays
  # as 0 (validation/compare/chi2_generic, f_generic), are found.
  s1 <- g1_gen("chi2_generic", "sensitivity", alpha = 0.001, power = 0.99, df = 1000)
  expect_identical(s1$status$code, "ok")
  expect_gt(s1$solved$lambda, 100)
  g1_gen_prob(ref_pchisq(s1$quantities$critical, 1000, s1$solved$lambda, lower_tail = FALSE), 0.99, "chi2 sens")
  s2 <- g1_gen("f_generic", "sensitivity", alpha = 0.05, power = 0.575, df1 = 1, df2 = 1)
  expect_identical(s2$status$code, "ok")
  expect_gt(s2$solved$lambda, 100)
  expect_equal(s2$achieved$power, 0.575, tolerance = 1e-9)
  # Criterion and compromise for every generic test.
  for (cs in list(list(m = "t_generic", tails = "two", a = list(delta = 2.5, df = 30)),
                  list(m = "t_generic", tails = "one", a = list(delta = -1, df = 1)),
                  list(m = "z_generic", tails = "two", a = list(mu = 2.5, sigma = 1)),
                  list(m = "f_generic", tails = NULL, a = list(lambda = 12, df1 = 3, df2 = 60)),
                  list(m = "chi2_generic", tails = NULL, a = list(lambda = 10, df = 5)))) {
    info <- cs$m
    cr <- do.call(g1_gen, c(list(cs$m, "criterion", cs$tails), cs$a, list(power = 0.8)))
    expect_identical(cr$status$code, "ok", info = info)
    fw <- do.call(g1_gen, c(list(cs$m, "post_hoc", cs$tails), cs$a, list(alpha = cr$solved$alpha)))
    expect_equal(fw$achieved$power, 0.8, tolerance = 1e-9, info = info)
    co <- do.call(g1_gen, c(list(cs$m, "compromise", cs$tails), cs$a, list(error_ratio = 2)))
    expect_identical(co$status$code, "ok", info = info)
    expect_equal(co$achieved$beta, 2 * co$solved$alpha, tolerance = 1e-9, info = info)
    se <- do.call(g1_gen, c(list(cs$m, "sensitivity", cs$tails), cs$a[setdiff(names(cs$a), c("delta", "mu", "lambda"))],
                            list(alpha = 0.05, power = 0.8)))
    expect_identical(se$status$code, "ok", info = info)
    expect_equal(se$achieved$power, 0.8, tolerance = 1e-9, info = info)
  }
})

test_that("domains and typed outcomes of the generic tests", {
  inv <- function(...) tryCatch(analysis_spec(...), inferguidR_invalid_input = function(e) e)
  expect_s3_class(inv("t_generic", "post_hoc", inputs = list(delta = 2, alpha = 0.05, df = 0)), "inferguidR_invalid_input")
  expect_s3_class(inv("f_generic", "post_hoc", inputs = list(lambda = -1, alpha = 0.05, df1 = 3, df2 = 60)),
                  "inferguidR_invalid_input")
  expect_s3_class(inv("chi2_generic", "post_hoc", inputs = list(lambda = 1, alpha = 0.05, df = 0)),
                  "inferguidR_invalid_input")
  expect_s3_class(inv("z_generic", "post_hoc", inputs = list(mu = 1, sigma = 0, alpha = 0.05)), "inferguidR_invalid_input")
  # A target below the power at no effect is infeasible (power = alpha at
  # lambda = 0 for F and chi-square).
  low <- g1_gen("chi2_generic", "sensitivity", alpha = 0.05, power = 0.04, df = 5)
  expect_identical(low$status$code, "infeasible")
  expect_identical(g1_gen_numbers(low), character())
  # Criterion with a solution alpha within 1e-12 of one is infeasible.
  far <- g1_gen("z_generic", "criterion", "one", mu = 0.001, sigma = 5, power = 0.9999)
  expect_identical(far$status$code, "infeasible")
  expect_identical(g1_gen_numbers(far), character())
})
