# Repeated-measures ANOVA (group G3): f_rm_anova_between, f_rm_anova_within,
# f_rm_anova_within_between. Formulas and evidence: specs/method-notes/<id>.md.
# Curated G*Power comparisons: displayed strings of emulator_capture records
# (G*Power 3.1.9.7 under Wine), compared within the rounding interval of the
# last displayed digit (validation/tolerances.md); record ids are cited.

g3_display_ok <- function(value, shown) {
  x <- as.numeric(shown)
  dec <- if (grepl("e", shown, ignore.case = TRUE)) {
    mant <- sub("[eE].*$", "", shown)
    nchar(sub("^[^.]*\\.?", "", mant)) - as.integer(sub("^.*[eE]", "", shown))
  } else nchar(sub("^[^.]*\\.?", "", shown))
  abs(value - x) <= 0.5 * 10^(-dec) + 8 * .Machine$double.eps * max(1, abs(x))
}

g3_calc <- function(id, mode, inputs, options = list()) {
  calculate_power(analysis_spec(id, mode, inputs = inputs, options = options))
}

test_that("the three RM-ANOVA definitions are registered with the certified options", {
  for (id in c("f_rm_anova_between", "f_rm_anova_within", "f_rm_anova_within_between")) {
    def <- .ig_get_method(id)
    expect_s3_class(def, "ig_method")
    expect_identical(def$family, "F")
    expect_identical(def$modes, c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"))
    expect_identical(def$options$effect_size_specification$values,
                     .ig_get_metadata(id)$options$effect_size_specification$codes)
    expect_identical(def$options$effect_size_specification$default, "gpower30")
    expect_identical(.ig_implementation_status(id), "implemented")
    expect_identical(def$status, "independently_checked")
  }
  expect_null(.ig_get_method("f_rm_anova_between")$params$epsilon)
  expect_identical(.ig_get_method("f_rm_anova_within")$params$epsilon$default, 1)
})

test_that("forward models follow the stated df and noncentrality for every convention", {
  f <- 0.3; N <- 90; k <- 3; m <- 5; rho <- 0.4; eps <- 0.6
  lam <- function(id, sp) {
    inp <- list(f = f, alpha = 0.05, n_total = N, k = k, m = m)
    if (sp == "gpower30") inp$rho <- rho
    if (id != "f_rm_anova_between") inp$epsilon <- eps
    r <- g3_calc(id, "post_hoc", inp, list(effect_size_specification = sp))
    c(r$quantities$ncp, r$quantities$df)
  }
  expect_equal(lam("f_rm_anova_between", "gpower30"), c(f^2 * N * m / (1 + (m - 1) * rho), k - 1, N - k), ignore_attr = TRUE)
  expect_equal(lam("f_rm_anova_between", "cohen1988"), c(f^2 * N, k - 1, N - k), ignore_attr = TRUE)
  expect_equal(lam("f_rm_anova_between", "spss"), c(f^2 * (N - k), k - 1, N - k), ignore_attr = TRUE)
  expect_equal(lam("f_rm_anova_within", "gpower30"),
               c(f^2 * N * m * eps / (1 - rho), (m - 1) * eps, (N - k) * (m - 1) * eps), ignore_attr = TRUE)
  expect_equal(lam("f_rm_anova_within", "gpower30_implicit_rho"), c(f^2 * N * eps, (m - 1) * eps, (N - k) * (m - 1) * eps),
               ignore_attr = TRUE)
  expect_equal(lam("f_rm_anova_within_between", "spss"),
               c(f^2 * (N - k) * (m - 1) * eps, (k - 1) * (m - 1) * eps, (N - k) * (m - 1) * eps), ignore_attr = TRUE)
  expect_identical(lam("f_rm_anova_within", "cohen1988"), lam("f_rm_anova_within", "gpower30_implicit_rho"))
})

test_that("post hoc agrees with G*Power captures (every convention)", {
  cases <- list(
    # GPW-F-RMBET-POS-B001-O01--r01
    list("f_rm_anova_between", list(f = 0.25, alpha = 0.05, n_total = 60, k = 3, m = 3, rho = 0.5), "gpower30",
         c(ncp = "5.6250000", crit = "3.1588427", power = "0.5316231")),
    # GPW-F-RMBET-POS-B003-E01--r01
    list("f_rm_anova_between", list(f = 1, alpha = 0.05, n_total = 4, k = 2, m = 2), "spss",
         c(ncp = "2.0000000", crit = "18.5128205", power = "0.1382528")),
    # GPW-F-RMWIT-POS-B001-D01--r01
    list("f_rm_anova_within", list(f = 0.02, alpha = 0.0001, n_total = 100000, k = 20, m = 20, rho = -0.04,
                                   epsilon = 0.0526316), "gpower30",
         c(ncp = "40.4858462", crit = "15.1379223", power = "0.9932823")),
    # GPW-F-RMWIT-POS-B002-D01--r01
    list("f_rm_anova_within", list(f = 0.02, alpha = 0.0001, n_total = 100000, k = 20, m = 20, epsilon = 0.0526316),
         "gpower30_implicit_rho", c(ncp = "2.1052640", crit = "15.1379223", power = "0.007349897")),
    # SUP-G3-RMWIT-POS-X04--r01 (rho = -0.9, epsilon = 0.5)
    list("f_rm_anova_within", list(f = 0.25, alpha = 0.05, n_total = 60, k = 2, m = 4, rho = -0.9, epsilon = 0.5),
         "gpower30", c(ncp = "3.9473684", crit = "3.4270477", power = "0.4410632")),
    # SUP-G3-RMWIT-POS-X03--r01 (a single group)
    list("f_rm_anova_within", list(f = 0.25, alpha = 0.05, n_total = 60, k = 1, m = 4, rho = 0.5, epsilon = 1),
         "gpower30", c(ncp = "30.0000000", crit = "2.6556472", power = "0.9979918")),
    # GPW-F-RMWXB-POS-B003-D01--r01
    list("f_rm_anova_within_between", list(f = 0.02, alpha = 0.0001, n_total = 100000, k = 20, m = 20,
                                           epsilon = 0.0526316), "spss",
         c(ncp = "39.9920160", crit = "2.6738983", power = "0.7044601"))
  )
  for (cs in cases) {
    r <- g3_calc(cs[[1]], "post_hoc", cs[[2]], list(effect_size_specification = cs[[3]]))
    expect_identical(r$status$code, "ok")
    expect_true(g3_display_ok(r$quantities$ncp, cs[[4]][["ncp"]]), info = cs[[1]])
    expect_true(g3_display_ok(r$quantities$critical, cs[[4]][["crit"]]), info = cs[[1]])
    expect_true(g3_display_ok(r$achieved$power, cs[[4]][["power"]]), info = cs[[1]])
  }
})

test_that("a priori uses equal groups (N a multiple of k) as G*Power does", {
  # SUP-G3-RMBET-APR-K01--r01: the smallest N over all integers is 193 (power 0.9009), G*Power 196
  r <- g3_calc("f_rm_anova_between", "a_priori", list(f = 0.2, alpha = 0.05, power = 0.9, k = 4, m = 3, rho = 0.3))
  expect_equal(r$solved$n_total, 196)
  expect_true(g3_display_ok(r$achieved$power, "0.9050457"))
  expect_equal(unname(r$group_sizes), rep(49, 4))
  # SUP-G3-RMWIT-APR-K01--r01 (96, not 94) and GPW-F-RMWIT-APR-B001-O01--r01 (24)
  r <- g3_calc("f_rm_anova_within", "a_priori",
               list(f = 0.15, alpha = 0.05, power = 0.9, k = 3, m = 4, rho = 0.5, epsilon = 0.8))
  expect_equal(r$solved$n_total, 96)
  expect_true(g3_display_ok(r$quantities$df[2], "223.2"))
  expect_true(g3_display_ok(r$achieved$power, "0.9069420"))
  # SUP-G3-RMWXB-APR-K02--r01 (implicit rho, 284 not 281)
  r <- g3_calc("f_rm_anova_within_between", "a_priori", list(f = 0.25, alpha = 0.05, power = 0.8, k = 4, m = 3, epsilon = 0.7),
               list(effect_size_specification = "gpower30_implicit_rho"))
  expect_equal(r$solved$n_total, 284)
  expect_true(g3_display_ok(r$achieved$power, "0.8054354"))
  # GPW-F-RMWIT-APR-B002-D01--r01: denominator df beyond 4e5 (critical F of the finite df)
  r <- g3_calc("f_rm_anova_within", "a_priori", list(f = 0.02, alpha = 0.001, power = 0.99, k = 20, m = 20,
                                                     epsilon = 0.0526316),
               list(effect_size_specification = "gpower30_implicit_rho"))
  expect_equal(r$solved$n_total, 1498600)
  expect_true(g3_display_ok(r$quantities$critical, "10.8276060"))
  expect_true(g3_display_ok(r$achieved$power, "0.9900002"))
})

test_that("sensitivity, criterion and compromise agree with G*Power captures", {
  # GPW-F-RMBET-SEN-B003-E01--r01
  r <- g3_calc("f_rm_anova_between", "sensitivity", list(alpha = 0.05, power = 0.95, n_total = 4, k = 2, m = 2),
               list(effect_size_specification = "spss"))
  expect_true(g3_display_ok(r$solved$f, "5.4953957"))
  expect_true(g3_display_ok(r$quantities$ncp, "60.3987483"))
  # GPW-F-RMWIT-SEN-B001-D01--r01
  r <- g3_calc("f_rm_anova_within", "sensitivity", list(alpha = 0.001, power = 0.99, n_total = 100000, k = 20, m = 20,
                                                        rho = -0.04, epsilon = 0.0526316))
  expect_true(g3_display_ok(r$solved$f, "0.0176557"))
  # GPW-F-RMWIT-CRI-B001-O01--r01
  r <- g3_calc("f_rm_anova_within", "criterion", list(f = 0.25, power = 0.8, n_total = 60, k = 2, m = 4, rho = 0.5))
  expect_true(g3_display_ok(r$solved$alpha, "7.228082e-05"))
  expect_true(g3_display_ok(r$quantities$critical, "7.7128407"))
  # GPW-F-RMWIT-CRI-B002-D01--r01: alpha near one (the power is computed as 1 - beta there)
  r <- g3_calc("f_rm_anova_within", "criterion", list(f = 0.02, power = 0.9999, n_total = 100000, k = 20, m = 20,
                                                      epsilon = 0.0526316),
               list(effect_size_specification = "gpower30_implicit_rho"))
  expect_identical(r$status$code, "ok")
  expect_true(g3_display_ok(r$solved$alpha, "0.9997135"))
  expect_true(g3_display_ok(r$quantities$critical, "1.289524e-07"))
  # GPW-F-RMWXB-COM-B001-D01--r01
  r <- g3_calc("f_rm_anova_within_between", "compromise",
               list(f = 0.02, error_ratio = 0.01, n_total = 100000, k = 20, m = 20, rho = -0.04, epsilon = 0.0526316))
  expect_true(g3_display_ok(r$solved$alpha, "0.1578034"))
  expect_true(g3_display_ok(r$achieved$beta, "0.001578034"))
  expect_true(g3_display_ok(r$quantities$critical, "1.3202515"))
})

test_that("inadmissible inputs are typed refusals, as G*Power's input errors", {
  refuse <- function(id, mode, inputs, options = list()) {
    out <- tryCatch(g3_calc(id, mode, inputs, options), inferguidR_error = function(e) e$status$code)
    if (is.character(out)) out else out$status$code
  }
  w <- list(f = 0.25, alpha = 0.05, n_total = 60, k = 2, m = 4, rho = 0.5, epsilon = 1)
  # SUP-G3-RMWIT-POS-X01/X05: epsilon outside [1/(m - 1), 1]
  expect_identical(refuse("f_rm_anova_within", "post_hoc", modifyList(w, list(epsilon = 0.2))), "invalid_input")
  expect_identical(refuse("f_rm_anova_within", "post_hoc", modifyList(w, list(epsilon = 1.2))), "invalid_input")
  # SUP-G3-RMWIT-POS-X02, SUP-G3-RMBET-POS-X04: rho = 1 ('outside the admissible range (-0.999999, 0.999999)')
  expect_identical(refuse("f_rm_anova_within", "post_hoc", modifyList(w, list(rho = 1))), "invalid_input")
  b <- list(f = 0.25, alpha = 0.05, n_total = 60, k = 3, m = 3, rho = 0.5)
  expect_identical(refuse("f_rm_anova_between", "post_hoc", modifyList(b, list(rho = 1))), "invalid_input")
  # SUP-G3-RMBET-POS-X02/X03: 'There must be at least 2 groups / measurements.'
  expect_identical(refuse("f_rm_anova_between", "post_hoc", modifyList(b, list(k = 1))), "invalid_input")
  expect_identical(refuse("f_rm_anova_between", "post_hoc", modifyList(b, list(m = 1))), "invalid_input")
  # SUP-G3-RMBET-POS-X01: 1 + (m - 1) rho < 0 (G*Power: internal math error)
  expect_identical(refuse("f_rm_anova_between", "post_hoc", modifyList(b, list(m = 4, rho = -0.5))), "invalid_input")
  # GPW-F-RMBET-POS-B001-I03--r01: N not above k; GPW-F-RMBET-APR-B001-I02: f = 0
  expect_identical(refuse("f_rm_anova_between", "post_hoc", modifyList(b, list(n_total = 3))), "invalid_input")
  expect_identical(refuse("f_rm_anova_between", "a_priori",
                          list(f = 0, alpha = 0.05, power = 0.8, k = 3, m = 3, rho = 0.5)), "invalid_input")
  # GPW-F-RMWIT-APR-B001-I03--r01: power below alpha
  expect_identical(refuse("f_rm_anova_within", "a_priori",
                          list(f = 0.25, alpha = 0.05, power = 0.04, k = 2, m = 4, rho = 0.5)), "invalid_input")
  # the correlation is not an input of the conventions without it
  expect_identical(refuse("f_rm_anova_between", "post_hoc", b, list(effect_size_specification = "spss")), "invalid_input")
})

test_that("the critical value stays exact where stats::qf() switches to the chi-square limit", {
  # stats::qf() uses qchisq(p, df1) / df1 for df2 > 4e5; GPW-F-RMMWB-SEN-B010-D01--r01 shows 1.5566587
  expect_true(g3_display_ok(.ig_rm_qf_upper(0.001, 81, 436735.2863), "1.5566587"))
  expect_false(g3_display_ok(stats::qf(0.001, 81, 436735.2863, lower.tail = FALSE), "1.5566587"))
  for (a in c(1e-300, 1e-12, 0.05, 0.5, 0.9, 1 - 1e-9)) {
    x <- .ig_rm_qf_upper(a, 3, 40)
    tail <- if (a <= 0.5) stats::pf(x, 3, 40, lower.tail = FALSE) else stats::pf(x, 3, 40)
    expect_equal(tail, if (a <= 0.5) a else 1 - a, tolerance = 1e-12)
  }
})

test_that("huge noncentralities return power 1 without summing a huge Poisson mixture", {
  t0 <- proc.time()[["elapsed"]]
  out <- .ig_rm_f_test(3, 500, 1e10, 0.05)
  expect_identical(out$power, 1)
  expect_identical(out$beta, 0)
  expect_lt(proc.time()[["elapsed"]] - t0, 2)
  # the bound is not used where beta is representable
  mid <- .ig_rm_f_test(3, 500, 30, 0.05)
  expect_gt(mid$beta, 0)
  expect_equal(mid$power + mid$beta, 1, tolerance = 1e-15)
})
