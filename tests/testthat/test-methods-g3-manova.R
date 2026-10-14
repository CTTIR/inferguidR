# Repeated-measures MANOVA (group G3): f_rm_manova_between, f_rm_manova_within,
# f_rm_manova_within_between. Formulas and evidence: specs/method-notes/<id>.md.
# Curated G*Power comparisons (emulator_capture records, rounding interval of
# the last displayed digit); record ids are cited.

g3m_display_ok <- function(value, shown) {
  x <- as.numeric(shown)
  dec <- if (grepl("e", shown, ignore.case = TRUE)) {
    mant <- sub("[eE].*$", "", shown)
    nchar(sub("^[^.]*\\.?", "", mant)) - as.integer(sub("^.*[eE]", "", shown))
  } else nchar(sub("^[^.]*\\.?", "", shown))
  abs(value - x) <= 0.5 * 10^(-dec) + 8 * .Machine$double.eps * max(1, abs(x))
}

g3m_calc <- function(id, mode, inputs, options = list()) {
  calculate_power(analysis_spec(id, mode, inputs = inputs, options = options))
}

test_that("the three RM-MANOVA definitions are registered with the reconciled options", {
  for (id in c("f_rm_manova_between", "f_rm_manova_within", "f_rm_manova_within_between")) {
    def <- .ig_get_method(id)
    meta <- .ig_get_metadata(id)
    expect_identical(def$options$approximation$values, meta$options$approximation$codes)
    expect_identical(def$options$approximation$default, "pillai_v")
    expect_identical(def$options$algorithm$values, meta$options$algorithm$codes)
    expect_identical(def$options$algorithm$default, "obrien_shieh_1999")
    expect_identical(.ig_implementation_status(id), "implemented")
  }
  # D25: 'Use mean correlation' is user-reachable in both states for between/within
  # (validation/raw/emulator/supplement/g3/d25/), fixed unchecked for within-between.
  expect_identical(.ig_get_method("f_rm_manova_between")$options$use_mean_correlation$values, c("on", "off"))
  expect_identical(.ig_get_method("f_rm_manova_within")$options$use_mean_correlation$default, "on")
  expect_identical(.ig_get_method("f_rm_manova_within_between")$options$use_mean_correlation$values, "off")
})

test_that("the effect conversion with the correlation matches the captures with rho != 0", {
  # SUP-G3-RMMB-POS-R01/R02 and SUP-G3-RMMW-POS-R01/R02 (Pillai V, O'Brien-Shieh)
  cases <- list(
    list("f_rm_manova_between", 0.5, c(ncp = "5.6250000", power = "0.5316231", crit = "0.0857143")),
    list("f_rm_manova_between", -0.3, c(ncp = "28.1250000", power = "0.9978431", crit = "0.3191489")),
    list("f_rm_manova_within", 0.5, c(ncp = "22.5000000", power = "0.9897392", crit = "0.2727273")),
    list("f_rm_manova_within", -0.3, c(ncp = "8.6538462", power = "0.7288763", crit = "0.1260504"))
  )
  for (cs in cases) {
    r <- g3m_calc(cs[[1]], "post_hoc", list(f = 0.25, alpha = 0.05, n_total = 60, k = 3, m = 3, rho = cs[[2]]))
    expect_true(g3m_display_ok(r$quantities$ncp, cs[[3]][["ncp"]]), info = cs[[1]])
    expect_true(g3m_display_ok(r$achieved$power, cs[[3]][["power"]]), info = cs[[1]])
    expect_true(g3m_display_ok(r$method_outputs$pillai_v, cs[[3]][["crit"]]), info = cs[[1]])
  }
})

test_that("'Use mean correlation' unchecked takes the MANOVA effect size directly (D25)", {
  # SUP-G3-RMMB-POS-U01--r01 and SUP-G3-RMMW-POS-U01--r01
  r <- g3m_calc("f_rm_manova_between", "post_hoc", list(f = 0.25, alpha = 0.05, n_total = 60, k = 3, m = 3),
                list(use_mean_correlation = "off"))
  expect_true(g3m_display_ok(r$quantities$ncp, "3.7500000"))
  expect_true(g3m_display_ok(r$achieved$power, "0.3744311"))
  r <- g3m_calc("f_rm_manova_within", "post_hoc", list(f = 0.25, alpha = 0.05, n_total = 60, k = 3, m = 3),
                list(use_mean_correlation = "off"))
  expect_true(g3m_display_ok(r$achieved$power, "0.3741177"))
  # SUP-G3-RMMB-POS-U03--r01: Wilks U with Muller-Peterson shows the criterion (0.9411765)
  r <- g3m_calc("f_rm_manova_between", "post_hoc", list(f = 0.25, alpha = 0.05, n_total = 60, k = 3, m = 3),
                list(use_mean_correlation = "off", approximation = "wilks_u", algorithm = "muller_peterson_1984"))
  expect_true(g3m_display_ok(r$quantities$ncp, "3.5625000"))
  expect_true(g3m_display_ok(r$method_outputs$wilks_u, "0.9411765"))
  # SUP-G3-RMMB-APR-U05--r01: a priori 159 (multiple of k = 3; smallest integer N 158)
  r <- g3m_calc("f_rm_manova_between", "a_priori", list(f = 0.25, alpha = 0.05, power = 0.8, k = 3, m = 3),
                list(use_mean_correlation = "off"))
  expect_equal(r$solved$n_total, 159)
  expect_true(g3m_display_ok(r$achieved$power, "0.8048873"))
  # campaign records of the restored branch: GPW-F-RMMB-COM-B002-O01 (Wilks, OS) and -B004-O01 (Wilks, MP)
  r <- g3m_calc("f_rm_manova_between", "compromise", list(f = 0.25, error_ratio = 1, n_total = 60, k = 3, m = 3),
                list(use_mean_correlation = "off", approximation = "wilks_u"))
  expect_true(g3m_display_ok(r$solved$alpha, "0.2723162"))
  expect_true(g3m_display_ok(r$quantities$critical, "1.3309333"))
  r <- g3m_calc("f_rm_manova_between", "compromise", list(f = 0.25, error_ratio = 1, n_total = 60, k = 3, m = 3),
                list(use_mean_correlation = "off", approximation = "wilks_u", algorithm = "muller_peterson_1984"))
  expect_true(g3m_display_ok(r$solved$alpha, "0.2801999"))
  expect_true(g3m_display_ok(r$method_outputs$wilks_u, "0.9411765"))
  # rho is not an input when the box is unchecked
  expect_error(analysis_spec("f_rm_manova_between", "post_hoc",
                             inputs = list(f = 0.25, alpha = 0.05, n_total = 60, k = 3, m = 3, rho = 0.5),
                             options = list(use_mean_correlation = "off")), class = "inferguidR_invalid_input")
})

test_that("between/within: every approximation and algorithm agrees with the captures", {
  # GPW-F-RMMB-APR-B001-O01 (Wilks, OS), -B003-O01 (Wilks, MP; criterion field not compared, see the notes)
  r <- g3m_calc("f_rm_manova_between", "a_priori", list(f = 0.25, alpha = 0.05, power = 0.8, k = 3, m = 3, rho = 0),
                list(approximation = "wilks_u"))
  expect_equal(r$solved$n_total, 57)
  expect_true(g3m_display_ok(r$achieved$power, "0.8196637"))
  expect_true(g3m_display_ok(r$method_outputs$wilks_u, "0.8421053"))
  r <- g3m_calc("f_rm_manova_between", "a_priori", list(f = 0.25, alpha = 0.05, power = 0.8, k = 3, m = 3, rho = 0),
                list(approximation = "wilks_u", algorithm = "muller_peterson_1984"))
  expect_equal(r$solved$n_total, 60)
  expect_true(g3m_display_ok(r$quantities$ncp, "10.6875000"))
  # GPW-F-RMMW-SEN-B009-D01--r01 (T2, OS) and GPW-F-RMMW-COM-B005-O01--r01 (T1, OS)
  r <- g3m_calc("f_rm_manova_within", "sensitivity", list(alpha = 0.001, power = 0.99, n_total = 100000, k = 10, m = 10,
                                                          rho = 0), list(approximation = "hotelling_lawley_t2"))
  expect_true(g3m_display_ok(r$solved$f, "0.006944651"))
  expect_true(g3m_display_ok(r$method_outputs$hotelling_lawley_t2, "0.0004822817"))
  r <- g3m_calc("f_rm_manova_within", "compromise", list(f = 0.25, error_ratio = 1, n_total = 60, k = 3, m = 3, rho = 0),
                list(approximation = "hotelling_lawley_t1"))
  expect_true(g3m_display_ok(r$solved$alpha, "0.0955069"))
  expect_true(g3m_display_ok(r$quantities$critical, "2.4498642"))
  # SUP-G3-RMMW-POS-X01--r01: a single group, df2 = N - k - m + 2
  r <- g3m_calc("f_rm_manova_within", "post_hoc", list(f = 0.25, alpha = 0.05, n_total = 20, k = 1, m = 3, rho = 0.5))
  expect_equal(unname(r$quantities$df), c(2, 18))
  expect_true(g3m_display_ok(r$achieved$power, "0.6088495"))
  # SUP-G3-RMMB-CRI-D01--r01: criterion mode with m = 4, rho = 0 (where G*Power's defaults coincide)
  r <- g3m_calc("f_rm_manova_between", "criterion", list(f = 0.25, power = 0.8, n_total = 60, k = 3, m = 4, rho = 0))
  expect_true(g3m_display_ok(r$solved$alpha, "0.0102222"))
  expect_true(g3m_display_ok(r$quantities$critical, "4.9722909"))
})

test_that("within-between: the four F approximations agree with the captures", {
  base <- list(alpha = 0.05, power = 0.8, k = 3, m = 3)
  # GPW-F-RMMWB-APR-B002-O01 (Wilks, OS), -B006-O01 (T1, OS), -B010-O01 (T2, OS), -B014-O01 (Pillai, OS), -B012-O01 (T2, MP)
  exp_n <- c(wilks_u = 98, hotelling_lawley_t1 = 98, hotelling_lawley_t2 = 164, pillai_v = 98)
  exp_pw <- c(wilks_u = "0.8003705", hotelling_lawley_t1 = "0.8002489", hotelling_lawley_t2 = "0.8008683",
              pillai_v = "0.8004895")
  for (ap in names(exp_n)) {
    r <- g3m_calc("f_rm_manova_within_between", "a_priori", c(base, f = 0.25), list(approximation = ap))
    expect_equal(r$solved$n_total, exp_n[[ap]], info = ap)
    expect_true(g3m_display_ok(r$achieved$power, exp_pw[[ap]]), info = ap)
  }
  r <- g3m_calc("f_rm_manova_within_between", "a_priori", c(base, f = 0.25),
                list(approximation = "hotelling_lawley_t2", algorithm = "muller_peterson_1984"))
  expect_equal(r$solved$n_total, 168)
  expect_true(g3m_display_ok(r$quantities$df[2], "195.7647"))
  expect_true(g3m_display_ok(r$method_outputs$hotelling_lawley_t2, "0.0747549"))
  # GPW-F-RMMWB-POS-B002-D01--r01: Wilks with p = q = 9 (Rao's g = 6.46)
  r <- g3m_calc("f_rm_manova_within_between", "post_hoc", list(f = 0.02, alpha = 0.0001, n_total = 100000, k = 10, m = 10),
                list(approximation = "wilks_u"))
  expect_true(g3m_display_ok(r$quantities$ncp, "258.5013"))
  expect_true(g3m_display_ok(r$quantities$df[2], "646145.9"))
  expect_true(g3m_display_ok(r$quantities$critical, "1.6923560"))
  # GPW-F-RMMWB-SEN-B010-D01--r01: T2 sensitivity with df2 beyond 4e5
  r <- g3m_calc("f_rm_manova_within_between", "sensitivity", list(alpha = 0.001, power = 0.99, n_total = 100000, k = 10, m = 10),
                list(approximation = "hotelling_lawley_t2"))
  expect_true(g3m_display_ok(r$solved$f, "0.0147683"))
  expect_true(g3m_display_ok(r$quantities$critical, "1.5566587"))
  # GPW-F-RMMWB-CRI-B016-E01--r01 (Pillai, MP)
  r <- g3m_calc("f_rm_manova_within_between", "criterion", list(f = 1, power = 0.95, n_total = 6, k = 2, m = 2),
                list(algorithm = "muller_peterson_1984"))
  expect_true(g3m_display_ok(r$solved$alpha, "0.6822235"))
})

test_that("the engine keeps the criterion and the entered design where G*Power's display does not", {
  # GPW-F-RMMB-SEN-B003-E01--r01: G*Power shows 'Wilks U = 1.7423003' (the effect size f; a Wilks
  # U above 1 is impossible). The engine reports the criterion of the solved effect.
  r <- g3m_calc("f_rm_manova_between", "sensitivity", list(alpha = 0.05, power = 0.95, n_total = 6, k = 2, m = 2, rho = 0),
                list(approximation = "wilks_u", algorithm = "muller_peterson_1984"))
  expect_true(g3m_display_ok(r$solved$f, "1.7423003"))
  expect_lt(r$method_outputs$wilks_u, 1)
  expect_equal(r$method_outputs$wilks_u, 1 / (1 + 2 * r$solved$f^2))
  # SUP-G3-RMMB-CRI-D02--r01: G*Power's criterion-mode lambda uses m = 4 (lambda 15); the engine uses m = 6
  r <- g3m_calc("f_rm_manova_between", "criterion", list(f = 0.25, power = 0.8, n_total = 60, k = 3, m = 6, rho = 0))
  expect_equal(r$quantities$ncp, 0.0625 * 60 * 6)
})

test_that("Hotelling-Lawley T2 with s > 1 is computed down to error df p + 2 as G*Power does", {
  t2 <- function(mode, inputs) {
    g3m_calc("f_rm_manova_within_between", mode, c(inputs, list(k = 3, m = 4)), list(approximation = "hotelling_lawley_t2"))
  }
  # REV-G3-RMMWB-POS-E02--r01: nu_e = p + 3 (N = 9), McKeon's df2 at its limit 4
  r <- t2("post_hoc", list(f = 0.42, alpha = 0.05, n_total = 9))
  expect_equal(unname(r$quantities$df[2]), 4)
  expect_true(g3m_display_ok(r$quantities$ncp, "1.5876000"))
  expect_true(g3m_display_ok(r$quantities$critical, "6.1631323"))
  expect_true(g3m_display_ok(r$achieved$power, "0.0732819"))
  expect_true(g3m_display_ok(r$method_outputs$hotelling_lawley_t2, "0.1764000"))
  expect_true(any(grepl("no finite variance", r$assumptions, fixed = TRUE)))
  # REV-G3-RMMWB-POS-E05--r01: nu_e = p + 2 (N = 8), continuation of the formula
  r <- t2("post_hoc", list(f = 0.42, alpha = 0.05, n_total = 8))
  expect_true(g3m_display_ok(unname(r$quantities$df[2]), "2.8571429"))
  expect_true(g3m_display_ok(r$quantities$ncp, "1.2096000"))
  expect_true(g3m_display_ok(r$quantities$critical, "9.6415979"))
  expect_true(g3m_display_ok(r$achieved$power, "0.0633836"))
  expect_true(g3m_display_ok(r$method_outputs$hotelling_lawley_t2, "0.1512000"))
  # REV-G3-RMMWB-POS-E30--r01 and REV-G3-RMMWB-APR-E07--r01: N = 8 misses the target (0.6664361),
  # N = 9 is G*Power's a priori answer (actual power 0.9314582)
  r <- t2("post_hoc", list(f = 3, alpha = 0.05, n_total = 8))
  expect_true(g3m_display_ok(r$achieved$power, "0.6664361"))
  r <- t2("a_priori", list(f = 3, alpha = 0.05, power = 0.8))
  expect_identical(as.numeric(r$solved$n_total), 9)
  expect_true(g3m_display_ok(r$achieved$power, "0.9314582"))
  expect_true(g3m_display_ok(r$method_outputs$hotelling_lawley_t2, "9.0000000"))
  # no note above the moment condition
  r <- t2("post_hoc", list(f = 0.42, alpha = 0.05, n_total = 10))
  expect_false(any(grepl("no finite variance", r$assumptions, fixed = TRUE)))
})

test_that("inadmissible designs and inputs are typed refusals", {
  refuse <- function(id, mode, inputs, options = list()) {
    out <- tryCatch(g3m_calc(id, mode, inputs, options), inferguidR_error = function(e) e$status$code)
    if (is.character(out)) out else out$status$code
  }
  # SUP-G3-RMMB-POS-X01: 'There must be at least 2 measurements.'
  expect_identical(refuse("f_rm_manova_between", "post_hoc",
                          list(f = 0.25, alpha = 0.05, n_total = 60, k = 3, m = 1, rho = 0.5)), "invalid_input")
  # GPW-F-RMMW-COM-B001-I03: N = 0
  expect_identical(refuse("f_rm_manova_within", "compromise",
                          list(f = 0.25, error_ratio = 1, n_total = 0, k = 3, m = 3, rho = 0)), "invalid_input")
  # McKeon's T2 with s > 1 needs nu_e >= p + 2 (at nu_e = p + 1 its multiplier divides by zero;
  # G*Power does not return there: review probe REV-G3-RMMWB-POS-E06, no record)
  expect_identical(refuse("f_rm_manova_within_between", "post_hoc", list(f = 0.42, alpha = 0.05, n_total = 7, k = 3, m = 4),
                          list(approximation = "hotelling_lawley_t2")), "invalid_input")
  # REV-G3-RMMWB-POS-E11: T1 with denominator df 0; REV-G3-RMMW-POS-E04: within test with denominator df 0
  expect_identical(refuse("f_rm_manova_within_between", "post_hoc", list(f = 0.42, alpha = 0.05, n_total = 6, k = 3, m = 4),
                          list(approximation = "hotelling_lawley_t1")), "invalid_input")
  expect_identical(refuse("f_rm_manova_within", "post_hoc",
                          list(f = 0.42, alpha = 0.05, n_total = 7, k = 4, m = 5, rho = 0)), "invalid_input")
  # GPW-F-RMMWB-COM-B001-I01: beta/alpha ratio 0
  expect_identical(refuse("f_rm_manova_within_between", "compromise",
                          list(f = 0.25, error_ratio = 0, n_total = 60, k = 3, m = 3)), "invalid_input")
})
