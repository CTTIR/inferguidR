# Tests for the MANOVA procedures f_manova_global and f_manova_special
# (group G4; R/methods-manova.R, R/effect-manova.R). Method notes:
# specs/method-notes/f_manova_global.md and f_manova_special.md. The slow
# independent check is validation/independent/g4/check_f_manova.R; the full
# G*Power comparison is validation/compare/f_manova_global/ and
# validation/compare/f_manova_special/. Displayed G*Power strings below are
# copied from the emulator-capture fixture records named by their record ids
# (validation/fixtures/emulator_capture/, emulator_capture_supplement/g4/ and
# emulator_capture_supplement/review-G4-manova/); when the fixtures are
# present the copies are checked against them.

g4_ids <- c("f_manova_global", "f_manova_special")
g4_approx <- c("wilks_u", "hotelling_lawley_t1", "hotelling_lawley_t2", "pillai_v")
g4_alg <- c("obrien_shieh_1999", "muller_peterson_1984")

g4_run <- function(method, mode, inputs, approximation = "pillai_v", algorithm = "obrien_shieh_1999") {
  calculate_power(analysis_spec(method, mode, inputs = inputs,
                                options = list(approximation = approximation, algorithm = algorithm)))
}

# Closed rounding interval of a displayed G*Power string (validation/tolerances.md
# section 4): half a unit of the last displayed digit plus 8 eps max(1, |x|).
g4_within <- function(display, value) {
  x <- as.numeric(display)
  m <- regmatches(display, regexec("^-?[0-9]*\\.?([0-9]*)(e([-+]?[0-9]+))?$", display))[[1]]
  dec <- nchar(m[2])
  hw <- if (nzchar(m[3])) 0.5 * 10^(as.numeric(m[4]) - dec) else 0.5 * 10^(-dec)
  isTRUE(abs(value - x) <= hw + 8 * .Machine$double.eps * max(1, abs(x)))
}

# Engine values keyed by G*Power output labels (as validation/compare/f_manova_global/mapping_common.R).
g4_outputs <- function(res) {
  q <- res$quantities
  mo <- res$method_outputs
  lab <- c(pillai_v = "Pillai V", wilks_u = "Wilks U", hotelling_lawley_t1 = "Hotelling-Lawley T1",
           hotelling_lawley_t2 = "Hotelling-Lawley T2")[[mo$criterion]]
  eff <- c(pillai_v = "Effect size f\u00b2(V)", wilks_u = "Effect size f\u00b2(U)",
           hotelling_lawley_t1 = "Effect size f\u00b2(T1)", hotelling_lawley_t2 = "Effect size f\u00b2(T2)")[[mo$criterion]]
  out <- list("Noncentrality parameter \u03bb" = q$ncp, "Critical F" = q$critical, "Numerator df" = q$df[1],
              "Denominator df" = q$df[2], "Total sample size" = res$solved$n_total, "Actual power" = res$achieved$power,
              "Power (1-\u03b2 err prob)" = res$achieved$power, "\u03b1 err prob" = res$solved$alpha,
              "\u03b2 err prob" = res$achieved$beta)
  out[[lab]] <- mo$criterion_value
  out[[eff]] <- res$solved$f2
  out
}

# Curated fixture records, one or more per analysis type (record id = name).
g4_fixtures <- list(
  "GPW-F-MANG-APR-B007-O01--r01" = list(method = "f_manova_global", mode = "a_priori", options = c(approximation = "pillai_v", algorithm = "obrien_shieh_1999"),
    inputs = list(f2 = 0.0625, alpha = 0.05, power = 0.80, n_groups = 3, n_responses = 3),
    outputs = c("Noncentrality parameter \u03bb" = "14.2500000", "Critical F" = "2.1399548", "Numerator df" = "6.0000000", "Denominator df" = "220", "Total sample size" = "114", "Actual power" = "0.8076177", "Pillai V" = "0.1176471")),
  "GPW-F-MANG-APR-B006-O01--r01" = list(method = "f_manova_global", mode = "a_priori", options = c(approximation = "hotelling_lawley_t2", algorithm = "muller_peterson_1984"),
    inputs = list(f2 = 0.0625, alpha = 0.05, power = 0.80, n_groups = 3, n_responses = 3),
    outputs = c("Noncentrality parameter \u03bb" = "14.2227626", "Critical F" = "2.1385678", "Numerator df" = "6.0000000", "Denominator df" = "227.5642", "Total sample size" = "177", "Actual power" = "0.8071704", "Hotelling-Lawley T2" = "0.0829280")),
  "GPW-F-MANS-APR-B001-O01--r01" = list(method = "f_manova_special", mode = "a_priori", options = c(approximation = "wilks_u", algorithm = "obrien_shieh_1999"),
    inputs = list(f2 = 0.0625, alpha = 0.05, power = 0.80, n_groups = 3, n_predictors = 2, n_responses = 3),
    outputs = c("Noncentrality parameter \u03bb" = "14.1250000", "Critical F" = "2.1407279", "Numerator df" = "6.0000000", "Denominator df" = "216", "Total sample size" = "113", "Actual power" = "0.8033002", "Wilks U" = "0.8858131")),
  "SUP-G4-MANG-APR-A22--r01" = list(method = "f_manova_global", mode = "a_priori", options = c(approximation = "wilks_u", algorithm = "muller_peterson_1984"),
    inputs = list(f2 = 0.1, alpha = 0.05, power = 0.90, n_groups = 4, n_responses = 3),
    outputs = c("Noncentrality parameter \u03bb" = "20.9452008", "Critical F" = "1.9247842", "Numerator df" = "9.0000000", "Denominator df" = "209.452", "Total sample size" = "92", "Actual power" = "0.9050730", "Wilks U" = "0.7929779")),
  "GPW-F-MANG-APR-B003-D01--r01" = list(method = "f_manova_global", mode = "a_priori", options = c(approximation = "hotelling_lawley_t1", algorithm = "obrien_shieh_1999"),
    inputs = list(f2 = 0.0004, alpha = 0.001, power = 0.99, n_groups = 20, n_responses = 20),
    outputs = c("Noncentrality parameter \u03bb" = "176.472", "Critical F" = "1.2393918", "Numerator df" = "380", "Denominator df" = "440403", "Total sample size" = "23220", "Actual power" = "0.9900578", "Hotelling-Lawley T1" = "0.0076")),
  "GPW-F-MANG-POS-B003-O01--r01" = list(method = "f_manova_global", mode = "post_hoc", options = c(approximation = "hotelling_lawley_t1", algorithm = "obrien_shieh_1999"),
    inputs = list(f2 = 0.0625, alpha = 0.05, n_total = 90, n_groups = 3, n_responses = 3),
    outputs = c("Noncentrality parameter \u03bb" = "11.2500000", "Critical F" = "2.1529110", "Numerator df" = "6.0000000", "Denominator df" = "168", "Power (1-\u03b2 err prob)" = "0.6857889", "Hotelling-Lawley T1" = "0.1250000")),
  "GPW-F-MANS-POS-B008-O01--r01" = list(method = "f_manova_special", mode = "post_hoc", options = c(approximation = "pillai_v", algorithm = "muller_peterson_1984"),
    inputs = list(f2 = 0.0625, alpha = 0.05, n_total = 90, n_groups = 3, n_predictors = 2, n_responses = 3),
    outputs = c("Noncentrality parameter \u03bb" = "10.7500000", "Critical F" = "2.1516329", "Numerator df" = "6.0000000", "Denominator df" = "172", "Power (1-\u03b2 err prob)" = "0.6625896", "Pillai V" = "0.1176471")),
  "SUP-G4-MANS-POS-B07--r01" = list(method = "f_manova_special", mode = "post_hoc", options = c(approximation = "wilks_u", algorithm = "obrien_shieh_1999"),
    inputs = list(f2 = 0.05, alpha = 0.05, n_total = 80, n_groups = 8, n_predictors = 3, n_responses = 4),
    outputs = c("Noncentrality parameter \u03bb" = "10.5830052", "Critical F" = "1.8054384", "Numerator df" = "12.0000000", "Denominator df" = "182.8483", "Power (1-\u03b2 err prob)" = "0.5051334", "Wilks U" = "0.8788978")),
  "GPW-F-MANG-POS-B005-D01--r01" = list(method = "f_manova_global", mode = "post_hoc", options = c(approximation = "hotelling_lawley_t2", algorithm = "obrien_shieh_1999"),
    inputs = list(f2 = 0.0004, alpha = 0.0001, n_total = 100000, n_groups = 20, n_responses = 20),
    outputs = c("Noncentrality parameter \u03bb" = "381.9605", "Critical F" = "1.2926411", "Numerator df" = "380", "Denominator df" = "954511.7", "Power (1-\u03b2 err prob)" = "1.0000000", "Hotelling-Lawley T2" = "0.003819605")),
  "SUP-G4-MANG-POS-A16--r01" = list(method = "f_manova_global", mode = "post_hoc", options = c(approximation = "wilks_u", algorithm = "obrien_shieh_1999"),
    inputs = list(f2 = 0.0625, alpha = 0.05, n_total = 90, n_groups = 3, n_responses = 1),
    outputs = c("Noncentrality parameter \u03bb" = "5.6250000", "Critical F" = "3.1012958", "Numerator df" = "2.0000000", "Denominator df" = "87.0000000", "Power (1-\u03b2 err prob)" = "0.5396133", "Wilks U" = "0.9411765")),
  "GPW-F-MANG-SEN-B001-O01--r01" = list(method = "f_manova_global", mode = "sensitivity", options = c(approximation = "wilks_u", algorithm = "obrien_shieh_1999"),
    inputs = list(alpha = 0.05, power = 0.80, n_total = 90, n_groups = 3, n_responses = 3),
    outputs = c("Noncentrality parameter \u03bb" = "14.1351189", "Critical F" = "2.1522643", "Numerator df" = "6.0000000", "Denominator df" = "170", "Effect size f\u00b2(U)" = "0.0785284", "Wilks U" = "0.8596800")),
  "GPW-F-MANS-SEN-B006-O01--r01" = list(method = "f_manova_special", mode = "sensitivity", options = c(approximation = "hotelling_lawley_t2", algorithm = "muller_peterson_1984"),
    inputs = list(alpha = 0.05, power = 0.80, n_total = 90, n_groups = 3, n_predictors = 2, n_responses = 3),
    outputs = c("Noncentrality parameter \u03bb" = "14.4077047", "Critical F" = "2.1808830", "Numerator df" = "6.0000000", "Denominator df" = "111.5731", "Effect size f\u00b2(T2)" = "0.1291324", "Hotelling-Lawley T2" = "0.1704752")),
  "SUP-G4-MANG-SEN-A19--r01" = list(method = "f_manova_global", mode = "sensitivity", options = c(approximation = "wilks_u", algorithm = "obrien_shieh_1999"),
    inputs = list(alpha = 0.05, power = 0.80, n_total = 60, n_groups = 4, n_responses = 3),
    outputs = c("Noncentrality parameter \u03bb" = "16.6707074", "Critical F" = "1.9517291", "Numerator df" = "9.0000000", "Denominator df" = "131.5724", "Effect size f\u00b2(U)" = "0.1141640", "Wilks U" = "0.7686669")),
  "GPW-F-MANG-CRI-B007-O01--r01" = list(method = "f_manova_global", mode = "criterion", options = c(approximation = "pillai_v", algorithm = "obrien_shieh_1999"),
    inputs = list(f2 = 0.0625, power = 0.80, n_total = 90, n_groups = 3, n_responses = 3),
    outputs = c("Noncentrality parameter \u03bb" = "11.2500000", "Critical F" = "1.7798447", "Numerator df" = "6.0000000", "Denominator df" = "172", "\u03b1 err prob" = "0.1058373", "Pillai V" = "0.1176471")),
  "GPW-F-MANS-CRI-B004-O01--r01" = list(method = "f_manova_special", mode = "criterion", options = c(approximation = "hotelling_lawley_t1", algorithm = "muller_peterson_1984"),
    inputs = list(f2 = 0.0625, power = 0.80, n_total = 90, n_groups = 3, n_predictors = 2, n_responses = 3),
    outputs = c("Noncentrality parameter \u03bb" = "10.5000000", "Critical F" = "1.6844941", "Numerator df" = "6.0000000", "Denominator df" = "168", "\u03b1 err prob" = "0.1276126", "Hotelling-Lawley T1" = "0.1250000")),
  "GPW-F-MANG-CRI-B001-E01--r01" = list(method = "f_manova_global", mode = "criterion", options = c(approximation = "wilks_u", algorithm = "obrien_shieh_1999"),
    inputs = list(f2 = 1, power = 0.95, n_total = 4, n_groups = 2, n_responses = 1),
    outputs = c("Noncentrality parameter \u03bb" = "4.0000000", "Critical F" = "0.2070060", "Numerator df" = "1.0000000", "Denominator df" = "2.0000000", "\u03b1 err prob" = "0.6937404", "Wilks U" = "0.5000000")),
  "GPW-F-MANG-COM-B001-O01--r01" = list(method = "f_manova_global", mode = "compromise", options = c(approximation = "wilks_u", algorithm = "obrien_shieh_1999"),
    inputs = list(f2 = 0.0625, error_ratio = 1, n_total = 90, n_groups = 3, n_responses = 3),
    outputs = c("Noncentrality parameter \u03bb" = "11.2500000", "Critical F" = "1.5971541", "Numerator df" = "6.0000000", "Denominator df" = "170", "\u03b1 err prob" = "0.1507343", "\u03b2 err prob" = "0.1507343", "Power (1-\u03b2 err prob)" = "0.8492657", "Wilks U" = "0.8858131")),
  "GPW-F-MANS-COM-B005-O01--r01" = list(method = "f_manova_special", mode = "compromise", options = c(approximation = "hotelling_lawley_t2", algorithm = "obrien_shieh_1999"),
    inputs = list(f2 = 0.0625, error_ratio = 1, n_total = 90, n_groups = 3, n_predictors = 2, n_responses = 3),
    outputs = c("Noncentrality parameter \u03bb" = "7.4258893", "Critical F" = "1.3812265", "Numerator df" = "6.0000000", "Denominator df" = "111.5731", "\u03b1 err prob" = "0.2283352", "\u03b2 err prob" = "0.2283352", "Power (1-\u03b2 err prob)" = "0.7716648", "Hotelling-Lawley T2" = "0.0825099")),
  "GPW-F-MANG-COM-B008-E01--r01" = list(method = "f_manova_global", mode = "compromise", options = c(approximation = "pillai_v", algorithm = "muller_peterson_1984"),
    inputs = list(f2 = 1, error_ratio = 4, n_total = 4, n_groups = 2, n_responses = 1),
    outputs = c("Noncentrality parameter \u03bb" = "2.0000000", "Critical F" = "4.8899509", "Numerator df" = "1.0000000", "Denominator df" = "2.0000000", "\u03b1 err prob" = "0.1575499", "\u03b2 err prob" = "0.6301996", "Power (1-\u03b2 err prob)" = "0.3698004", "Pillai V" = "0.5000000")),
  "GPW-F-MANG-POS-B006-E01--r01" = list(method = "f_manova_global", mode = "post_hoc", options = c(approximation = "hotelling_lawley_t2", algorithm = "muller_peterson_1984"),
    inputs = list(f2 = 1, alpha = 0.05, n_total = 4, n_groups = 2, n_responses = 1),
    outputs = c("Noncentrality parameter \u03bb" = "2.0000000", "Critical F" = "18.5128205", "Numerator df" = "1.0000000", "Denominator df" = "2.0000000", "Power (1-\u03b2 err prob)" = "0.1382528", "Hotelling-Lawley T2" = "-nan(ind)")),
  # A priori for Hotelling-Lawley T2 with Muller-Peterson where a design with error df < p has a
  # higher formula power than its successors (review captures; see the a priori domain test below).
  "REV-G4-MANS-APR-S17--r01" = list(method = "f_manova_special", mode = "a_priori", options = c(approximation = "hotelling_lawley_t2", algorithm = "muller_peterson_1984"),
    inputs = list(f2 = 2, alpha = 0.05, power = 0.3, n_groups = 3, n_predictors = 5, n_responses = 8),
    outputs = c("Noncentrality parameter \u03bb" = "32.0000000", "Critical F" = "2.1507110", "Numerator df" = "40.0000000", "Denominator df" = "16.0000000", "Total sample size" = "19", "Actual power" = "0.3605771", "Hotelling-Lawley T2" = "4.0000000")),
  "REV-G4-MANG-APR-T12--r01" = list(method = "f_manova_global", mode = "a_priori", options = c(approximation = "hotelling_lawley_t2", algorithm = "muller_peterson_1984"),
    inputs = list(f2 = 1, alpha = 0.05, power = 0.8, n_groups = 4, n_responses = 9),
    outputs = c("Noncentrality parameter \u03bb" = "37.5209125", "Critical F" = "1.7851425", "Numerator df" = "27.0000000", "Denominator df" = "37.5209125", "Total sample size" = "32", "Actual power" = "0.8183523", "Hotelling-Lawley T2" = "1.9733840"))
)

# Definitions -------------------------------------------------------------------------------

test_that("both MANOVA procedures are registered with the certified options and modes", {
  expect_identical(.ig_registry_problems()[grepl("manova", .ig_registry_problems())], character())
  for (id in g4_ids) {
    def <- .ig_get_method(id)
    meta <- .ig_get_metadata(id)
    expect_s3_class(def, "ig_method")
    expect_identical(def$family, "F")
    expect_setequal(def$modes, meta$modes)
    expect_identical(def$status, "independently_checked")
    root <- th_pkg_root()
    expect_true(is.null(root) || file.exists(file.path(root, def$checker)))
    for (o in c("approximation", "algorithm")) {
      expect_setequal(def$options[[o]]$values, meta$options[[o]]$codes)
      expect_identical(def$options[[o]]$default, meta$options[[o]]$default)
    }
    expect_identical(def$options$approximation$default, "pillai_v")
    expect_identical(def$options$algorithm$default, "obrien_shieh_1999")
    pn <- c("f2", "alpha", "power", "error_ratio", "n_total", "n_groups", "n_responses",
            if (id == "f_manova_special") "n_predictors")
    expect_setequal(names(def$params), pn)
    expect_identical(def$targets$sensitivity, "f2")
  }
  # Live input labels of every analysis type map onto the parameters.
  map <- c("\u03b1 err prob" = "alpha", "Power (1-\u03b2 err prob)" = "power", "\u03b2/\u03b1 ratio" = "error_ratio",
           "Total sample size" = "n_total", "Number of groups" = "n_groups", "Number of predictors" = "n_predictors",
           "Response variables" = "n_responses", "Effect size f\u00b2(V)" = "f2")
  for (id in g4_ids) for (m in .ig_get_method(id)$modes) {
    def <- .ig_get_method(id)
    live <- .ig_live_labels(id, m)$inputs
    expect_true(all(live %in% names(map)), info = paste(id, m))
    used <- names(def$params)[vapply(def$params, function(p) is.null(p$modes) || m %in% p$modes, logical(1))]
    expect_setequal(unname(map[live]), setdiff(used, def$targets[[m]]))
  }
})

# Closed forms ------------------------------------------------------------------------------

test_that("degrees of freedom, noncentrality and criterion values follow the four approximations", {
  # Written in the SAS / Muller-Peterson (m, n) forms, independently of the engine's algebra.
  ref <- function(approx, ne, p, q, f2) {
    s <- min(p, q); m <- (abs(p - q) - 1) / 2; n <- (ne - p - 1) / 2
    switch(approx,
           pillai_v = list(df1 = s * (2 * m + s + 1), df2 = s * (2 * n + s + 1), g = s, crit = s * f2 / (1 + f2)),
           hotelling_lawley_t1 = list(df1 = p * q, df2 = 2 * (s * n + 1), g = s, crit = s * f2),
           wilks_u = {
             t <- if (p^2 + q^2 - 5 > 0) sqrt((p^2 * q^2 - 4) / (p^2 + q^2 - 5)) else 1
             list(df1 = p * q, df2 = (ne - (p - q + 1) / 2) * t - (p * q - 2) / 2, g = t, crit = (1 + f2)^(-t))
           },
           hotelling_lawley_t2 = {
             b <- (p + 2 * n) * (q + 2 * n) / (2 * (2 * n + 1) * (n - 1))
             cc <- (2 + (p * q + 2) / (b - 1)) / (2 * n)
             list(df1 = p * q, df2 = 4 + (p * q + 2) / (b - 1), g = cc, crit = cc * f2)
           })
  }
  for (id in g4_ids) for (ap in g4_approx) for (al in g4_alg) {
    for (d in list(c(N = 60, k = 4, p = 3, q = 3), c(N = 97, k = 6, p = 2, q = 5), c(N = 40, k = 3, p = 5, q = 2))) {
      inp <- list(f2 = 0.08, alpha = 0.05, n_total = d[["N"]], n_groups = d[["k"]], n_responses = d[["p"]])
      q <- d[["k"]] - 1
      if (id == "f_manova_special") { inp$n_predictors <- d[["q"]]; q <- d[["q"]] }
      r <- g4_run(id, "post_hoc", inp, ap, al)
      rf <- ref(ap, d[["N"]] - d[["k"]], d[["p"]], q, 0.08)
      info <- paste(id, ap, al, paste(d, collapse = "/"))
      expect_equal(r$quantities$df, c(rf$df1, rf$df2), tolerance = 1e-13, info = info)
      lam <- if (al == "obrien_shieh_1999") d[["N"]] * rf$g * 0.08 else rf$df2 * 0.08
      expect_equal(r$quantities$ncp, lam, tolerance = 1e-13, info = info)
      expect_equal(r$method_outputs$criterion_value, rf$crit, tolerance = 1e-13, info = info)
      expect_identical(r$method_outputs$criterion, ap)
    }
  }
})

test_that("for s = 1 Pillai V with O'Brien-Shieh noncentrality is Hotelling's exact two-sample T^2", {
  mu1 <- c(0, 0, 0); mu2 <- c(0.5, 0.2, -0.1); Sigma <- diag(3); n <- 25
  B <- 0.5 * tcrossprod(mu1 - (mu1 + mu2) / 2) + 0.5 * tcrossprod(mu2 - (mu1 + mu2) / 2)
  V <- sum(diag(B %*% solve(B + Sigma)))
  lam_exact <- n * n / (2 * n) * drop(crossprod(mu1 - mu2, solve(Sigma, mu1 - mu2)))
  r <- g4_run("f_manova_global", "post_hoc", list(f2 = V / (1 - V), alpha = 0.05, n_total = 2 * n, n_groups = 2, n_responses = 3))
  expect_equal(r$quantities$ncp, lam_exact, tolerance = 1e-12)
  expect_identical(r$quantities$df, c(3, 2 * n - 2 - 3 + 1))
  # p = 1: every approximation is the univariate F(q, N - k), lambda = N f2.
  for (ap in g4_approx) {
    r1 <- g4_run("f_manova_global", "post_hoc", list(f2 = 0.1, alpha = 0.05, n_total = 63, n_groups = 7, n_responses = 1), ap)
    expect_equal(r1$quantities$df, c(6, 56), tolerance = 1e-14, info = ap)
    expect_equal(r1$quantities$ncp, 6.3, tolerance = 1e-14, info = ap)
  }
})

# Probabilities against an independent reference -----------------------------------------------

test_that("power and beta agree with an integral reference that shares no routine with the engine", {
  # Route A of validation/independent/g4/ref_f_manova.R in compact form: outer
  # integral over the chi-square(df2) density, inner Poisson mixture of central
  # pchisq() (the engine sums Poisson mixtures of pbeta() and never integrates).
  ref_tail <- function(c, df1, df2, lam, lower) {
    inner <- function(y) {
      j <- 0:400
      w <- stats::dpois(j, lam / 2, log = TRUE) + stats::pchisq(y, df1 + 2 * j, lower.tail = lower, log.p = TRUE)
      exp(max(w)) * sum(exp(w - max(w)))
    }
    f <- function(v) vapply(v, function(vi) inner(c * df1 * vi / df2), numeric(1)) * stats::dchisq(v, df2)
    lo <- stats::qchisq(1e-15, df2); hi <- stats::qchisq(1e-15, df2, lower.tail = FALSE); md <- df2
    stats::integrate(f, lo, md, rel.tol = 1e-12, abs.tol = 0)$value + stats::integrate(f, md, hi, rel.tol = 1e-12, abs.tol = 0)$value
  }
  pts <- list(list("f_manova_global", "wilks_u", "obrien_shieh_1999", list(f2 = 0.0625, alpha = 0.05, n_total = 60, n_groups = 4, n_responses = 3)),
              list("f_manova_special", "hotelling_lawley_t2", "muller_peterson_1984",
                   list(f2 = 0.05, alpha = 0.01, n_total = 80, n_groups = 8, n_predictors = 3, n_responses = 4)),
              list("f_manova_global", "pillai_v", "obrien_shieh_1999", list(f2 = 1, alpha = 0.05, n_total = 12, n_groups = 2, n_responses = 1)))
  for (pt in pts) {
    r <- g4_run(pt[[1]], "post_hoc", pt[[4]], pt[[2]], pt[[3]])
    crit <- r$quantities$critical
    expect_equal(stats::pf(crit, r$quantities$df[1], r$quantities$df[2], lower.tail = FALSE), pt[[4]]$alpha, tolerance = 1e-10)
    pw <- ref_tail(crit, r$quantities$df[1], r$quantities$df[2], r$quantities$ncp, FALSE)
    bt <- ref_tail(crit, r$quantities$df[1], r$quantities$df[2], r$quantities$ncp, TRUE)
    expect_lt(abs(r$achieved$power - pw), 1e-9)
    expect_lt(abs(r$achieved$beta / bt - 1), 1e-7)
    expect_equal(r$achieved$power + r$achieved$beta, 1, tolerance = 1e-14)
  }
})

test_that("critical F stays exact where stats::qf() switches to the chi-square limit (df2 > 4e5)", {
  for (d in list(c(380, 954511.7), c(380, 1381548.27), c(3, 400001))) {
    crit <- .ig_manova_q_f_upper(1e-4, d[1], d[2])
    expect_equal(stats::pf(crit, d[1], d[2], lower.tail = FALSE), 1e-4, tolerance = 1e-10)
    # nmath qf() returns qchisq(p, df1) / df1 here, whose tail is off by more than 1e-4 relative.
    expect_gt(abs(stats::pf(stats::qf(1e-4, d[1], d[2], lower.tail = FALSE), d[1], d[2], lower.tail = FALSE) / 1e-4 - 1), 1e-4)
  }
  expect_identical(.ig_manova_q_f_upper(0.05, 6, 172), stats::qf(0.05, 6, 172, lower.tail = FALSE))
})

test_that("the Chernoff bound of the acceptance probability is rigorous and skips only underflowing tails", {
  for (x in list(c(2, 6, 172, 11.25), c(1.5, 380, 1381548, 552.811), c(1.24, 380, 1381548, 552.811), c(4, 9, 60, 200))) {
    lb <- .ig_manova_log_beta_bound(x[1], x[2], x[3], x[4])
    exact <- .ig_p_f(x[1], x[2], x[3], x[4], lower_tail = TRUE, log_p = TRUE)
    expect_gte(lb, exact - 1e-9)
  }
  tl <- .ig_manova_tails(1.2392914, 380, 1899601, 1e7)
  expect_identical(c(tl$power, tl$beta), c(1, 0))
  expect_true(tl$beta_bound)
  expect_false(.ig_manova_tails(2, 6, 172, 11.25)$beta_bound)
})

# Curated G*Power fixture subset ----------------------------------------------------------------

test_that("curated emulator-capture fixtures are reproduced in every analysis type", {
  for (id in names(g4_fixtures)) {
    fx <- g4_fixtures[[id]]
    inp <- lapply(fx$inputs, function(v) v)
    for (nm in intersect(names(inp), c("n_total", "n_groups", "n_predictors", "n_responses"))) inp[[nm]] <- as.integer(inp[[nm]])
    r <- g4_run(fx$method, fx$mode, inp, fx$options[["approximation"]], fx$options[["algorithm"]])
    expect_identical(r$status$code, "ok", info = id)
    eo <- g4_outputs(r)
    for (lab in names(fx$outputs)) {
      shown <- fx$outputs[[lab]]
      if (shown == "-nan(ind)") {
        expect_true(is.na(eo[[lab]]), info = paste(id, lab))
        expect_match(r$method_outputs$criterion_value_note, "McKeon", info = id)
        next
      }
      expect_true(g4_within(shown, eo[[lab]]), info = sprintf("%s %s: G*Power %s, engine %.12g", id, lab, shown, eo[[lab]]))
    }
  }
})

test_that("the curated strings equal the fixture records (development checkout)", {
  root <- th_pkg_root()
  skip_if(is.null(root), "fixture records are not part of the installed package")
  fx_dirs <- file.path(root, "validation", "fixtures", c("emulator_capture", file.path("emulator_capture_supplement", "g4"),
                                                         file.path("emulator_capture_supplement", "review-G4-manova")))
  skip_if_not(all(dir.exists(fx_dirs)), "fixture records are not part of the installed package")
  for (id in names(g4_fixtures)) {
    f <- file.path(fx_dirs, paste0(id, ".json"))
    f <- f[file.exists(f)]
    expect_length(f, 1L)
    rec <- jsonlite::fromJSON(f, simplifyVector = FALSE)
    shown <- vapply(rec$outputs, `[[`, "", "displayed")
    names(shown) <- vapply(rec$outputs, `[[`, "", "label_as_displayed")
    shown <- shown[!vapply(rec$outputs, function(o) startsWith(o$name, "message"), logical(1))]
    expect_identical(shown, g4_fixtures[[id]]$outputs, info = id)
    expect_identical(rec$method$method_id, g4_fixtures[[id]]$method)
    expect_identical(rec$method$analysis_type, g4_fixtures[[id]]$mode)
  }
})

# Typed outcomes ----------------------------------------------------------------------------------

test_that("G*Power's refusals are typed refusals without numbers", {
  code <- function(expr) th_condition_code(expr)
  # GPW-F-MANG-APR-B001-I03 ('Power must be greater than or equal to alpha.')
  expect_identical(code(analysis_spec("f_manova_global", "a_priori", inputs = list(f2 = 0.0625, alpha = 0.05, power = 0.04,
                                                                                  n_groups = 3, n_responses = 3))), "invalid_input")
  # GPW-F-MANG-APR-B001-I02 (effect size 0 outside the admissible range)
  expect_identical(code(analysis_spec("f_manova_global", "a_priori", inputs = list(f2 = 0, alpha = 0.05, power = 0.8,
                                                                                  n_groups = 3, n_responses = 3))), "invalid_input")
  # GPW-F-MANG-POS-B001-I03 (N = 0) and SUP-G4-MANG-POS-A05 (one group)
  expect_identical(code(analysis_spec("f_manova_global", "post_hoc", inputs = list(f2 = 0.0625, alpha = 0.05, n_total = 0,
                                                                                  n_groups = 3, n_responses = 3))), "invalid_input")
  expect_identical(code(analysis_spec("f_manova_global", "post_hoc", inputs = list(f2 = 0.0625, alpha = 0.05, n_total = 20,
                                                                                  n_groups = 1, n_responses = 2))), "invalid_input")
  expect_identical(code(analysis_spec("f_manova_special", "post_hoc", inputs = list(f2 = 0.0625, alpha = 0.05, n_total = 6,
                                                                                   n_groups = 6, n_predictors = 3, n_responses = 2))), "invalid_input")
  # SUP-G4-MANG-POS-A07: Pillai df2 = 0 (G*Power: 'Internal math error').
  cond <- tryCatch(analysis_spec("f_manova_global", "post_hoc", inputs = list(f2 = 0.5, alpha = 0.05, n_total = 6,
                                                                              n_groups = 2, n_responses = 5)),
                   inferguidR_error = function(e) e)
  expect_s3_class(cond, "inferguidR_invalid_input")
  expect_match(conditionMessage(cond), "denominator df is 0")
  # McKeon's scale c = 0/0 at error df p + 1: O'Brien-Shieh is refused (G*Power does not finish:
  # case GPW-F-MANG-POS-B005-E01 timed out), Muller-Peterson computes the power and leaves T2
  # undefined with a reason (GPW-F-MANG-POS-B006-E01 shows '-nan(ind)').
  expect_identical(code(analysis_spec("f_manova_global", "post_hoc", inputs = list(f2 = 1, alpha = 0.05, n_total = 4,
                                                                                  n_groups = 2, n_responses = 1),
                                      options = list(approximation = "hotelling_lawley_t2"))), "invalid_input")
  r <- g4_run("f_manova_global", "post_hoc", list(f2 = 1, alpha = 0.05, n_total = 4L, n_groups = 2L, n_responses = 1L),
              "hotelling_lawley_t2", "muller_peterson_1984")
  expect_identical(r$status$code, "ok")
  expect_true(is.na(r$method_outputs$criterion_value))
  expect_true(g4_within("0.1382528", r$achieved$power))
})

test_that("the a priori search starts at the smallest admissible design", {
  # SUP-G4-MANG-APR-A04 / SUP-G4-MANS-APR-B04: G*Power reports 'Internal math error' (Pillai's df2 is
  # <= 0 at the smallest N); the engine finds the smallest admissible design. G*Power's failure is not a
  # start rule: it returns designs although the approximation is undefined at N = k + 1
  # (REV-G4-MANS-APR-S07: N = 12; REV-G4-MANG-APR-R08: N = 18), which the engine reproduces.
  r <- solve_design(analysis_spec("f_manova_global", "a_priori", inputs = list(f2 = 50, alpha = 0.05, power = 0.8,
                                                                              n_groups = 3, n_responses = 5)))
  expect_identical(r$solved$n_total, 9L)
  expect_identical(r$diagnostics$first_admissible_n_total, 9)
  expect_match(r$diagnostics$start_rule, "smallest N > k at which the approximation is defined")
  r <- solve_design(analysis_spec("f_manova_special", "a_priori", inputs = list(f2 = 50, alpha = 0.05, power = 0.8,
                                                                               n_groups = 3, n_predictors = 2, n_responses = 4)))
  expect_identical(r$solved$n_total, 6L)
  # SUP-G4-MANG-APR-A03 / SUP-G4-MANS-APR-B03: where N = k + 1 is admissible both agree (6 and 4).
  r <- solve_design(analysis_spec("f_manova_global", "a_priori", inputs = list(f2 = 50, alpha = 0.05, power = 0.8,
                                                                              n_groups = 3, n_responses = 2)))
  expect_identical(r$solved$n_total, 6L)
  r <- solve_design(analysis_spec("f_manova_special", "a_priori", inputs = list(f2 = 50, alpha = 0.05, power = 0.8,
                                                                               n_groups = 3, n_predictors = 2, n_responses = 2)))
  expect_identical(r$solved$n_total, 4L)
  r <- solve_design(analysis_spec("f_manova_special", "a_priori", inputs = list(f2 = 3, alpha = 0.05, power = 0.95, n_groups = 4,
                                                                               n_predictors = 3, n_responses = 5),
                                  options = list(algorithm = "muller_peterson_1984")))
  expect_identical(r$solved$n_total, 12L)
  expect_true(g4_within("0.9662805", r$achieved$power))
  r <- solve_design(analysis_spec("f_manova_global", "a_priori", inputs = list(f2 = 0.9, alpha = 0.05, power = 0.8, n_groups = 3,
                                                                              n_responses = 6),
                                  options = list(approximation = "hotelling_lawley_t1")))
  expect_identical(r$solved$n_total, 18L)
  expect_true(g4_within("0.8573913", r$achieved$power))
  # Global effects: multiples of k only; the predecessor (N - k) misses the target.
  r <- solve_design(analysis_spec("f_manova_global", "a_priori", inputs = list(f2 = 0.0625, alpha = 0.05, power = 0.8,
                                                                              n_groups = 3, n_responses = 3)))
  expect_identical(r$solved$n_total %% 3L, 0L)
  expect_identical(r$diagnostics$predecessor_n_total, r$solved$n_total - 3)
  expect_lt(r$diagnostics$predecessor_power, 0.8)
  # An unreachable target is a typed outcome.
  r <- solve_design(analysis_spec("f_manova_global", "a_priori", inputs = list(f2 = 1e-15, alpha = 0.05, power = 0.99,
                                                                              n_groups = 3, n_responses = 3)))
  expect_identical(r$status$code, "infeasible")
  expect_match(r$status$reason, "multiples of k = 3")
  expect_false(any(grepl("N <= 7158278", r$status$reason, fixed = TRUE)))
})

test_that("a priori designs of Hotelling-Lawley T2 lie in the domain of McKeon's approximation", {
  # Below error df p + 2 McKeon's scale c is 0/0 or negative. With Muller-Peterson the formula still
  # returns a power there, and at some error df < p it exceeds the power of every design up to a
  # much larger N: special effects, k = 3, 5 predictors, p = 8, f2 = 2 gives 0.3606 at N = 10 (error
  # df 7), less than 0.3 from N = 11 to 18 and 0.3606 again at N = 19. G*Power returns N = 19
  # (REV-G4-MANS-APR-S17) and, for global effects with k = 4, p = 9, f2 = 1, N = 32 where N = 12 has
  # formula power 0.896 (REV-G4-MANG-APR-T12); both are in the curated fixture test above.
  mp <- list(approximation = "hotelling_lawley_t2", algorithm = "muller_peterson_1984")
  inp <- list(f2 = 2, alpha = 0.05, n_groups = 3L, n_predictors = 5L, n_responses = 8L)
  pw <- function(N) calculate_power(analysis_spec("f_manova_special", "post_hoc", inputs = c(inp, list(n_total = N)), options = mp))
  expect_gt(pw(10L)$achieved$power, 0.36)
  expect_true(is.na(pw(10L)$method_outputs$criterion_value))
  expect_true(all(vapply(11:18, function(N) pw(N)$achieved$power, numeric(1)) < 0.3))
  r <- solve_design(analysis_spec("f_manova_special", "a_priori", inputs = c(inp, list(power = 0.3)), options = mp))
  expect_identical(r$solved$n_total, 19L)
  expect_identical(r$diagnostics$first_admissible_n_total, 13)        # error df 10 = p + 2
  expect_match(r$diagnostics$start_rule, "error df >= p \\+ 2")
  expect_gt(r$method_outputs$approximation_g, 0)
  # The design rule (a priori candidates) excludes error df < p + 2 for T2 only.
  def <- .ig_get_method("f_manova_special")
  p_t2 <- c(inp, list(alpha = 0.05, power = 0.3), mp)
  expect_null(def$design(p_t2, 10))
  expect_null(def$design(p_t2, 12))                                    # error df p + 1: c = 0/0
  expect_identical(def$design(p_t2, 13), c(n_total = 13))
  p_v <- p_t2
  p_v$approximation <- "pillai_v"
  expect_identical(def$design(p_v, 10), c(n_total = 10))
  # Fixed N is evaluated as G*Power does (REV-G4-MANS-POS-S15: power 0.9988374 at error df 5 < p = 6,
  # where G*Power shows T2 = -42.0000000; the engine reports the criterion as undefined with the reason).
  r <- calculate_power(analysis_spec("f_manova_special", "post_hoc",
                                     inputs = list(f2 = 3, alpha = 0.05, n_total = 10L, n_groups = 5L, n_predictors = 4L,
                                                   n_responses = 6L), options = mp))
  expect_true(g4_within("0.9988374", r$achieved$power))
  expect_true(is.na(r$method_outputs$criterion_value))
  expect_match(r$method_outputs$criterion_value_note, "c f2 = -42", fixed = TRUE)
  expect_true(any(grepl("singular", r$assumptions)))
  # The a priori answer for that problem is the design G*Power's own fixed-N numbers confirm
  # (REV-G4-MANS-POS-T03: 0.9453248 at N = 22; REV-G4-MANS-POS-T04: 0.8891200 at N = 21); G*Power's
  # a priori search ends in 'Internal math error' there (REV-G4-MANS-APR-S16).
  r <- solve_design(analysis_spec("f_manova_special", "a_priori",
                                  inputs = list(f2 = 3, alpha = 0.05, power = 0.9, n_groups = 5L, n_predictors = 4L,
                                                n_responses = 6L), options = mp))
  expect_identical(r$solved$n_total, 22L)
  expect_true(g4_within("0.9453248", r$achieved$power))
  expect_true(g4_within("0.8891200", r$diagnostics$predecessor_power))
  # O'Brien-Shieh admits the same designs (0 < c < Inf), so both algorithms share the candidates.
  for (N in 9:16) for (pq in list(c(4, 2), c(6, 4), c(8, 5), c(1, 1))) {
    tm <- .ig_manova_terms(N, 3, pq[1], pq[2], "hotelling_lawley_t2")
    os_ok <- is.null(.ig_manova_problem(tm, "hotelling_lawley_t2", "obrien_shieh_1999"))
    expect_identical(os_ok, (N - 3) >= pq[1] + 2, info = paste(N, pq[1], pq[2]))
  }
})

test_that("a total that is not a multiple of k keeps N in the error df (global effects)", {
  # SUP-G4-MANG-POS-A01: N = 100, k = 3: lambda 12.5, df2 194, power 0.8097158 (G*Power shows its
  # Pillai V field as -4.787243e+307, a display defect; the engine reports s f2 / (1 + f2)).
  r <- g4_run("f_manova_global", "post_hoc", list(f2 = 0.0625, alpha = 0.05, n_total = 100L, n_groups = 3L, n_responses = 2L))
  expect_true(g4_within("12.5000000", r$quantities$ncp))
  expect_identical(r$quantities$df[2], 194)
  expect_true(g4_within("0.8097158", r$achieved$power))
  expect_equal(r$method_outputs$criterion_value, 2 * 0.0625 / 1.0625)
  expect_match(r$method_outputs$group_size_note, "average group size")
})

test_that("inverse modes meet their targets and results round-trip through JSON", {
  res <- list(
    solve_design(analysis_spec("f_manova_special", "a_priori", inputs = list(f2 = 0.05, alpha = 0.05, power = 0.8, n_groups = 8,
                                                                            n_predictors = 3, n_responses = 4),
                               options = list(approximation = "hotelling_lawley_t2", algorithm = "muller_peterson_1984"))),
    g4_run("f_manova_global", "post_hoc", list(f2 = 0.1, alpha = 0.05, n_total = 60, n_groups = 4, n_responses = 3), "wilks_u"),
    g4_run("f_manova_global", "sensitivity", list(alpha = 0.05, power = 0.8, n_total = 60, n_groups = 4, n_responses = 3), "wilks_u"),
    g4_run("f_manova_special", "criterion", list(f2 = 0.05, power = 0.8, n_total = 80, n_groups = 8, n_predictors = 3, n_responses = 4)),
    g4_run("f_manova_special", "compromise", list(f2 = 0.05, error_ratio = 1, n_total = 80, n_groups = 8, n_predictors = 3,
                                                  n_responses = 4), "hotelling_lawley_t1")
  )
  for (r in res) {
    expect_identical(r$status$code, "ok", info = r$mode)
    expect_identical(ig_from_json(ig_to_json(r)), r, info = r$mode)
  }
  expect_lt(abs(res[[3]]$achieved$power - 0.8), 1e-9)
  expect_lt(abs(res[[4]]$achieved$power - 0.8), 1e-9)
  expect_lt(abs(res[[5]]$achieved$beta / res[[5]]$solved$alpha - 1), 1e-8)
  # SUP-G4-MANG-SEN-A19, SUP-G4-MANS-POS-B07, SUP-G4-MANS-APR-B18, SUP-G4-MANS-CRI-B16, SUP-G4-MANS-COM-B17.
  expect_true(g4_within("0.1141640", res[[3]]$solved$f2))
  expect_identical(res[[1]]$solved$n_total, 218L)
  expect_true(g4_within("0.1758687", res[[4]]$solved$alpha))
  expect_true(g4_within("0.1886154", res[[5]]$solved$alpha))
})

# Effect-size drawers ----------------------------------------------------------------------------

test_that("the drawers invert the criterion of the main analysis", {
  # Global effects: SUP-G4-MANG-DRW-D01 (Pillai V 0.4 -> 0.25), D02 (Wilks U 0.4, k = 4, p = 3 ->
  # 0.4571688: Rao's t, not s), D03 (T1 0.4 -> 0.1333333), D05 / D06 (T2 0.4 at N = 60 / 200).
  g <- "f_manova_global"
  expect_true(g4_within("0.25", effect_size(g, "pillai_v", pillai_v = 0.4, n_groups = 3, n_responses = 2)$value))
  expect_true(g4_within("0.4571688", effect_size(g, "wilks_u", wilks_u = 0.4, n_groups = 4, n_responses = 3)$value))
  expect_true(g4_within("0.1333333", effect_size(g, "hotelling_lawley_t1", hotelling_lawley_t1 = 0.4, n_groups = 4, n_responses = 3)$value))
  expect_true(g4_within("0.2608696", effect_size(g, "hotelling_lawley_t2", hotelling_lawley_t2 = 0.4, n_total = 60, n_groups = 4,
                                                 n_responses = 3)$value))
  expect_true(g4_within("0.2562648", effect_size(g, "hotelling_lawley_t2", hotelling_lawley_t2 = 0.4, n_total = 200, n_groups = 4,
                                                 n_responses = 3)$value))
  # Round trip through the forward model for every approximation of both procedures.
  for (id in g4_ids) for (ap in g4_approx) {
    x <- if (ap == "wilks_u") 0.7 else 0.3
    args <- list(method = id, input_mode = ap, n_responses = 4)
    args[[ap]] <- x
    if (id == "f_manova_global") args$n_groups <- 5 else args$n_predictors <- 3
    if (ap == "hotelling_lawley_t2") { args$n_total <- 90; args$n_groups <- 5 }
    es <- do.call(effect_size, args)
    inp <- list(f2 = es$value, alpha = 0.05, n_total = 90, n_groups = 5, n_responses = 4)
    if (id == "f_manova_special") inp$n_predictors <- 3
    expect_equal(g4_run(id, "post_hoc", inp, ap)$method_outputs$criterion_value, x, tolerance = 1e-13, info = paste(id, ap))
  }
  # Special effects: the engine uses q = number of predictors, as the main analysis does. G*Power's
  # drawer uses q - 1 (SUP-G4-MANS-DRW-E03: T1 0.4, p = 4, 3 predictors -> 0.2 = 0.4 / 2; the main
  # analysis at f2 = 0.2 reports T1 = 0.6, SUP-G4-MANS-POS-B19), a G*Power defect.
  es <- effect_size("f_manova_special", "hotelling_lawley_t1", hotelling_lawley_t1 = 0.4, n_responses = 4, n_predictors = 3)
  expect_equal(es$value, 0.4 / 3)
  r <- g4_run("f_manova_special", "post_hoc", list(f2 = 0.2, alpha = 0.05, n_total = 80L, n_groups = 8L, n_predictors = 3L,
                                                   n_responses = 4L), "hotelling_lawley_t1")
  expect_true(g4_within("0.6000000", r$method_outputs$criterion_value))
  # Inadmissible criterion values are refused (G*Power shows '-nan(ind)', SUP-G4-MANG-DRW-D07).
  expect_s3_class(tryCatch(effect_size(g, "pillai_v", pillai_v = 1.2, n_groups = 2, n_responses = 3), error = function(e) e),
                  "inferguidR_invalid_input")
  expect_s3_class(tryCatch(effect_size(g, "hotelling_lawley_t2", hotelling_lawley_t2 = 0.4, n_groups = 4, n_responses = 3),
                           error = function(e) e), "inferguidR_invalid_input")
})
