# Curated G*Power fixture comparisons for the G11 methods (layer 1 evidence:
# G*Power 3.1.9.7 under Wine, emulator_capture, D22/D23; documented examples
# are secondary evidence). The full comparison of every fixture runs from
# validation/compare/<method_id>/ (validation/tools/compare_fixtures.R); this
# file keeps a fast subset with at least one record per method and analysis
# type. Every displayed string below is copied from the fixture record named
# next to it (validation/fixtures/...); the last test checks that they still
# equal the records in a development checkout.
#
# Match rules (validation/tolerances.md): MR-EXACT for integers, MR-DISPLAY
# (closed rounding interval of the last displayed digit plus 8 eps max(1, |x|))
# for reals. Critical values are compared only where G*Power displays the
# critical values of the test (its criterion-mode display of the upper value is
# a confirmed reference defect, validation/compare/*/summary.md).

ref_reset_registry()

g11_display_ok <- function(ours, shown) {
  if (!is.finite(ours)) return(FALSE)
  s <- trimws(shown)
  parts <- strsplit(tolower(s), "e", fixed = TRUE)[[1]]
  mant <- parts[1]
  ex <- if (length(parts) == 2L) as.numeric(parts[2]) else 0
  dec <- if (grepl(".", mant, fixed = TRUE)) nchar(sub("^[^.]*\\.", "", mant)) else 0
  x <- as.numeric(s)
  half <- 0.5 * 10^(ex - dec)
  ours >= x - half - 8 * .Machine$double.eps * max(1, abs(x)) && ours <= x + half + 8 * .Machine$double.eps * max(1, abs(x))
}

g11_engine_outputs <- function(res) {
  q <- res$quantities
  crit <- q$critical
  lo <- if (length(crit)) crit[1] else NA_real_
  hi <- if (length(crit)) crit[length(crit)] else NA_real_
  df <- q$df
  gs <- res$group_sizes
  list("Lower critical F" = lo, "Upper critical F" = hi, "Lower critical \u03c7\u00b2" = lo, "Upper critical \u03c7\u00b2" = hi,
       "Critical \u03c7\u00b2" = lo, "Numerator df" = df[1], "Denominator df" = df[length(df)], "Df" = df[1],
       "Noncentrality parameter \u03bb" = q$ncp,
       "Sample size group 1" = if (!is.null(gs)) gs[["n1"]] else NA_real_,
       "Sample size group 2" = if (!is.null(gs)) gs[["n2"]] else NA_real_,
       "Total sample size" = res$solved$n_total %||% NA_real_, "Actual power" = res$achieved$power,
       "Power (1-\u03b2 err prob)" = res$achieved$power, "\u03b1 err prob" = res$solved$alpha %||% NA_real_,
       "\u03b2 err prob" = res$achieved$beta, "Ratio var1/var0" = res$solved$ratio %||% NA_real_,
       "Effect size w" = res$solved$w %||% NA_real_)
}

g11_integer_labels <- c("Numerator df", "Denominator df", "Df", "Sample size group 1", "Sample size group 2", "Total sample size")

F2 <- "f_variances_two_groups"
V1 <- "chi2_variance_one_sample"
GF <- "chi2_goodness_of_fit_contingency"

# id, method, mode, engine inputs, options, displayed outputs (label = string).
g11_records <- list(
  list("GPW-F-VAR2-APR-B001-O01--r01", F2, "a_priori", list(ratio = 1.5, alpha = 0.05, power = 0.80, allocation_ratio = 1), list(tails = "one"),
       c("Lower critical F" = "1.3069118", "Upper critical F" = "1.3069118", "Numerator df" = "152", "Denominator df" = "152",
         "Sample size group 1" = "153", "Sample size group 2" = "153", "Actual power" = "0.8017014")),
  list("GPW-F-VAR2-APR-B001-D01--r01", F2, "a_priori", list(ratio = 1.02, alpha = 0.001, power = 0.99, allocation_ratio = 0.1), list(tails = "one"),
       c("Lower critical F" = "1.0113786", "Numerator df" = "1647266", "Denominator df" = "164726",
         "Sample size group 1" = "1647267", "Sample size group 2" = "164727", "Actual power" = "0.9900001")),
  list("SUP-G11-F-VAR2-APR-A03--r01", F2, "a_priori", list(ratio = 0.5, alpha = 0.05, power = 0.80, allocation_ratio = 2), list(tails = "one"),
       c("Lower critical F" = "0.6228799", "Numerator df" = "40", "Denominator df" = "82", "Sample size group 1" = "41",
         "Sample size group 2" = "83", "Actual power" = "0.8004010")),
  # Review captures: allocation ratios whose first fractional design of G*Power's rule has a df next to 0
  # (0.3334: 1.5e-4; 2.9999: 2.5e-5), and a rounding tie that is not exact in binary (156 / 1.6 = 97.5 -> 98).
  list("REV-G11-F-VAR2-APR-R08--r01", F2, "a_priori", list(ratio = 1.8, alpha = 0.05, power = 0.80, allocation_ratio = 0.3334), list(tails = "one"),
       c("Lower critical F" = "1.4973107", "Numerator df" = "151", "Denominator df" = "50", "Sample size group 1" = "152",
         "Sample size group 2" = "51", "Actual power" = "0.8013214")),
  list("REV-G11-F-VAR2-APR-R09--r01", F2, "a_priori", list(ratio = 1.8, alpha = 0.05, power = 0.80, allocation_ratio = 2.9999), list(tails = "one"),
       c("Lower critical F" = "1.4507060", "Numerator df" = "47", "Denominator df" = "142", "Sample size group 1" = "48",
         "Sample size group 2" = "143", "Actual power" = "0.8015503")),
  list("REV-G11-F-VAR2-APR-R06--r01", F2, "a_priori", list(ratio = 0.5, alpha = 0.05, power = 0.90, allocation_ratio = 0.6), list(tails = "one"),
       c("Lower critical F" = "0.6842280", "Numerator df" = "97", "Denominator df" = "57", "Sample size group 1" = "98",
         "Sample size group 2" = "58", "Actual power" = "0.9004318")),
  list("GPW-F-VAR2-POS-B004-O01--r01", F2, "post_hoc", list(ratio = 1.5, alpha = 0.05, n1 = 50L, n2 = 50L), list(tails = "two"),
       c("Lower critical F" = "0.5674762", "Upper critical F" = "1.7621885", "Numerator df" = "49", "Denominator df" = "49",
         "Power (1-\u03b2 err prob)" = "0.2879654")),
  list("GPW-F-VAR2-POS-B001-E01--r01", F2, "post_hoc", list(ratio = 0.5, alpha = 0.05, n1 = 2L, n2 = 2L), list(tails = "one"),
       c("Lower critical F" = "0.006193959", "Power (1-\u03b2 err prob)" = "0.0705660")),
  list("SUP-G11-F-VAR2-POS-P01--r01", F2, "post_hoc", list(ratio = 1, alpha = 0.05, n1 = 20L, n2 = 30L), list(tails = "one"),
       c("Lower critical F" = "1.9581455", "Upper critical F" = "1.9581455", "Power (1-\u03b2 err prob)" = "0.0500000")),
  list("GPW-F-VAR2-SEN-B002-O01--r01", F2, "sensitivity", list(alpha = 0.05, power = 0.80, n1 = 50L, n2 = 50L),
       list(tails = "one", effect_direction = "var1_le_var0"),
       c("Lower critical F" = "0.6221655", "Ratio var1/var0" = "0.4884418")),
  list("GPW-F-VAR2-SEN-B006-O01--r01", F2, "sensitivity", list(alpha = 0.05, power = 0.80, n1 = 50L, n2 = 50L),
       list(tails = "two", effect_direction = "var1_ge_var0"),
       c("Lower critical F" = "0.5674762", "Upper critical F" = "1.7621885", "Ratio var1/var0" = "2.2446293")),
  list("GPW-F-VAR2-CRI-B001-O01--r01", F2, "criterion", list(ratio = 1.5, power = 0.80, n1 = 50L, n2 = 50L), list(tails = "one"),
       c("Lower critical F" = "1.1776009", "\u03b1 err prob" = "0.2846947")),
  list("GPW-F-VAR2-CRI-B004-O01--r01", F2, "criterion", list(ratio = 1.5, power = 0.80, n1 = 50L, n2 = 50L), list(tails = "two"),
       c("Lower critical F" = "0.8317710", "\u03b1 err prob" = "0.5215471")),
  list("GPW-F-VAR2-COM-B004-O01--r01", F2, "compromise", list(ratio = 1.5, error_ratio = 1, n1 = 50L, n2 = 50L), list(tails = "two"),
       c("Lower critical F" = "0.7537550", "Upper critical F" = "1.3266910", "\u03b1 err prob" = "0.3257306",
         "\u03b2 err prob" = "0.3257306", "Power (1-\u03b2 err prob)" = "0.6742694")),
  list("GPW-CHI2-VAR1-APR-B002-O01--r01", V1, "a_priori", list(ratio = 0.6666667, alpha = 0.05, power = 0.80), list(tails = "two"),
       c("Lower critical \u03c7\u00b2" = "75.9456997", "Upper critical \u03c7\u00b2" = "131.8375", "Df" = "102",
         "Total sample size" = "103", "Actual power" = "0.8024134")),
  list("GPW-CHI2-VAR1-APR-B001-D01--r01", V1, "a_priori", list(ratio = 1.02, alpha = 0.001, power = 0.99), list(tails = "one"),
       c("Lower critical \u03c7\u00b2" = "151193.5", "Df" = "149498", "Total sample size" = "149499", "Actual power" = "0.9900002")),
  list("GPW-CHI2-VAR1-POS-B002-O01--r01", V1, "post_hoc", list(ratio = 1.5, alpha = 0.05, n_total = 50L), list(tails = "two"),
       c("Lower critical \u03c7\u00b2" = "31.5549165", "Upper critical \u03c7\u00b2" = "70.2224136", "Df" = "49",
         "Power (1-\u03b2 err prob)" = "0.5623138")),
  list("SUP-G11-CHI2-VAR1-POS-P01--r01", V1, "post_hoc", list(ratio = 1, alpha = 0.05, n_total = 20L), list(tails = "one"),
       c("Lower critical \u03c7\u00b2" = "10.1170131", "Power (1-\u03b2 err prob)" = "0.0500000")),
  list("GPW-CHI2-VAR1-SEN-B002-O01--r01", V1, "sensitivity", list(alpha = 0.05, power = 0.80, n_total = 50L), list(tails = "two"),
       c("Lower critical \u03c7\u00b2" = "31.5549165", "Upper critical \u03c7\u00b2" = "70.2224136", "Ratio var1/var0" = "1.7323910")),
  list("GPW-CHI2-VAR1-SEN-B001-E01--r01", V1, "sensitivity", list(alpha = 0.05, power = 0.95, n_total = 2L), list(tails = "one"),
       c("Lower critical \u03c7\u00b2" = "3.8414588", "Ratio var1/var0" = "976.9385")),
  list("GPW-CHI2-VAR1-CRI-B001-O01--r01", V1, "criterion", list(ratio = 1.5, power = 0.80, n_total = 50L), list(tails = "one"),
       c("Lower critical \u03c7\u00b2" = "60.8016261", "\u03b1 err prob" = "0.1201716")),
  list("GPW-CHI2-VAR1-COM-B002-E01--r01", V1, "compromise", list(ratio = 4, error_ratio = 4, n_total = 2L), list(tails = "two"),
       c("Lower critical \u03c7\u00b2" = "0.008609574", "Upper critical \u03c7\u00b2" = "3.1935822", "\u03b1 err prob" = "0.1478557",
         "\u03b2 err prob" = "0.5914228", "Power (1-\u03b2 err prob)" = "0.4085772")),
  list("GPW-CHI2-GOF-APR-B001-O01--r01", GF, "a_priori", list(w = 0.3, alpha = 0.05, power = 0.80, df = 5), list(),
       c("Noncentrality parameter \u03bb" = "12.8700000", "Critical \u03c7\u00b2" = "11.0704977", "Total sample size" = "143",
         "Actual power" = "0.8015133")),
  list("GPW-CHI2-GOF-APR-B001-D01--r01", GF, "a_priori", list(w = 0.01, alpha = 0.001, power = 0.99, df = 200), list(),
       c("Noncentrality parameter \u03bb" = "135.0496", "Critical \u03c7\u00b2" = "267.5405", "Total sample size" = "1350496",
         "Actual power" = "0.9900001")),
  list("GPW-CHI2-GOF-POS-B001-D01--r01", GF, "post_hoc", list(w = 0.01, alpha = 0.0001, n_total = 100000L, df = 200), list(),
       c("Noncentrality parameter \u03bb" = "10.0000000", "Critical \u03c7\u00b2" = "283.0603", "Power (1-\u03b2 err prob)" = "0.0007405085")),
  list("SUP-G11-CHI2-GOF-POS-P01--r01", GF, "post_hoc", list(w = 0.3, alpha = 0.05, n_total = 150L, df = 2.5), list(),
       c("Noncentrality parameter \u03bb" = "13.5000000", "Critical \u03c7\u00b2" = "6.9280761", "Power (1-\u03b2 err prob)" = "0.9009149")),
  list("GPW-CHI2-GOF-SEN-B001-O01--r01", GF, "sensitivity", list(alpha = 0.05, power = 0.80, n_total = 150L, df = 5), list(),
       c("Noncentrality parameter \u03bb" = "12.8276066", "Critical \u03c7\u00b2" = "11.0704977", "Effect size w" = "0.2924335")),
  list("GPW-CHI2-GOF-CRI-B001-E01--r01", GF, "criterion", list(w = 0.9, power = 0.95, n_total = 2L, df = 1), list(),
       c("Noncentrality parameter \u03bb" = "1.6200000", "Critical \u03c7\u00b2" = "0.0197632", "\u03b1 err prob" = "0.8882006")),
  list("GPW-CHI2-GOF-COM-B001-D01--r01", GF, "compromise", list(w = 0.01, error_ratio = 0.01, n_total = 100000L, df = 200), list(),
       c("Noncentrality parameter \u03bb" = "10.0000000", "Critical \u03c7\u00b2" = "164.0937", "\u03b1 err prob" = "0.9700658",
         "\u03b2 err prob" = "0.009700658", "Power (1-\u03b2 err prob)" = "0.9902993"))
)

# Drawer records: id, method, var0, var1, displayed ratio.
g11_drawer_records <- list(
  list("DOC-096--r01", F2, 1.5, 1, "0.6666667"),
  list("GPW-UT-DRAW-UTL-U043-UE01--r02", F2, 0.000001, 1000000, "1e+12"),
  list("DOC-136--r01", V1, 1.5, 1, "0.6666667"),
  list("GPW-UT-DRAW-UTL-U074-UD01--r02", V1, 1, 1.000001, "1.000001")
)

# G*Power refusals (TD-0010): id, method, mode, inputs, options, accepted engine statuses.
g11_refusals <- list(
  list("GPW-F-VAR2-APR-B001-I03--r01", F2, "a_priori", list(ratio = 1.5, alpha = 0.05, power = 0.04, allocation_ratio = 1),
       list(tails = "one"), c("invalid_input", "infeasible")),
  list("GPW-F-VAR2-POS-B001-I03--r01", F2, "post_hoc", list(ratio = 1.5, alpha = 0.05, n1 = 0L, n2 = 50L), list(tails = "one"),
       "invalid_input"),
  list("GPW-CHI2-VAR1-APR-B001-I02--r01", V1, "a_priori", list(ratio = 1, alpha = 0.05, power = 0.80), list(tails = "one"),
       c("invalid_input", "infeasible")),
  list("GPW-CHI2-GOF-APR-B001-I02--r01", GF, "a_priori", list(w = 0, alpha = 0.05, power = 0.80, df = 5), list(),
       c("invalid_input", "infeasible")),
  list("SUP-G11-CHI2-GOF-POS-P02--r01", GF, "post_hoc", list(w = 0.3, alpha = 0.05, n_total = 150L, df = 0), list(), "invalid_input")
)

test_that("curated G*Power records are reproduced in every analysis type of the three methods", {
  seen <- character()
  for (r in g11_records) {
    info <- r[[1]]
    res <- calculate_power(analysis_spec(r[[2]], r[[3]], inputs = r[[4]], options = r[[5]]))
    expect_identical(res$status$code, "ok", info = info)
    eng <- g11_engine_outputs(res)
    for (lab in names(r[[6]])) {
      v <- eng[[lab]]
      if (lab %in% g11_integer_labels) {
        expect_identical(as.double(v), as.double(r[[6]][[lab]]), info = paste(info, lab))
      } else {
        expect_true(g11_display_ok(v, r[[6]][[lab]]), info = sprintf("%s %s: shown %s, engine %.12g", info, lab, r[[6]][[lab]], v))
      }
    }
    seen <- c(seen, paste(r[[2]], r[[3]]))
  }
  for (m in c(F2, V1, GF)) for (mode in .ig_modes) expect_true(paste(m, mode) %in% seen, info = paste(m, mode))
})

test_that("curated drawer records are reproduced", {
  for (r in g11_drawer_records) {
    es <- effect_size(r[[2]], "from_variances", var0 = r[[3]], var1 = r[[4]])
    expect_true(g11_display_ok(es$value, r[[5]]), info = r[[1]])
  }
})

test_that("G*Power refusals are refused by the engine without numbers (TD-0010)", {
  for (r in g11_refusals) {
    out <- tryCatch(calculate_power(analysis_spec(r[[2]], r[[3]], inputs = r[[4]], options = r[[5]])),
                    inferguidR_error = function(e) e)
    code <- if (inherits(out, "inferguidR_error")) out$status$code else out$status$code
    expect_true(code %in% r[[6]], info = paste(r[[1]], code))
    if (inherits(out, "power_result")) {
      x <- unclass(out)
      expect_identical(.ig_finite_paths(x[c("solved", "achieved", "quantities", "group_sizes")], ""), character(), info = r[[1]])
    }
  }
})

test_that("the embedded displayed strings equal the fixture records (development checkout)", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
  fx <- file.path(root, "validation", "fixtures")
  skip_if_not(dir.exists(fx), "validation fixtures not available (installed package)")
  find_record <- function(id) {
    cand <- c(file.path(fx, c("emulator_capture", "documented_example"), paste0(id, ".json")),
              file.path(fx, "emulator_capture_supplement", c("g11", "review-G11-variances-gof"), paste0(id, ".json")))
    cand[file.exists(cand)][1]
  }
  shown <- function(rec) {
    outs <- c(rec$outputs, rec$effect_size_drawer$outputs)
    stats::setNames(vapply(outs, `[[`, "", "displayed"), vapply(outs, `[[`, "", "label_as_displayed"))
  }
  for (r in g11_records) {
    f <- find_record(r[[1]])
    expect_false(is.na(f), info = r[[1]])
    rec <- jsonlite::fromJSON(f, simplifyVector = FALSE)
    expect_identical(rec$method$method_id, r[[2]], info = r[[1]])
    expect_identical(rec$method$analysis_type, r[[3]], info = r[[1]])
    s <- shown(rec)
    expect_identical(s[names(r[[6]])], r[[6]], info = r[[1]])
    # Every engine input appears among the record's displayed inputs.
    nums <- suppressWarnings(as.numeric(vapply(rec$inputs, function(i) i$displayed_in_reference, "")))
    for (nm in names(r[[4]])) {
      v <- as.double(r[[4]][[nm]])
      expect_true(any(abs(nums - v) <= 1e-12 * max(1, abs(v)), na.rm = TRUE), info = paste(r[[1]], nm))
    }
  }
  for (r in g11_drawer_records) {
    rec <- jsonlite::fromJSON(find_record(r[[1]]), simplifyVector = FALSE)
    expect_true(r[[5]] %in% shown(rec), info = r[[1]])
  }
  for (r in g11_refusals) {
    rec <- jsonlite::fromJSON(find_record(r[[1]]), simplifyVector = FALSE)
    expect_true(any(names(shown(rec)) == "Input error"), info = r[[1]])
  }
})
