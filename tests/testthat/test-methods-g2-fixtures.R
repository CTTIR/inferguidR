# Curated G*Power fixture comparisons for the G2 methods, one or more per mode and method.
#
# Every displayed value below is quoted from an emulator_capture fixture (G*Power 3.1.9.7 under
# Wine, D22/D23) or a documented example, named by its record id in validation/fixtures/
# (emulator_capture/, emulator_capture_supplement/g2/, documented_example/). Match rules are the
# records' own: integers exact (MR-EXACT), reals within the closed rounding interval of the
# displayed digits plus 8 eps max(1, |x|) (MR-DISPLAY). The complete comparison of all records is
# validation/compare/<method_id>/summary.md (validation/tools/compare_fixtures.R).

ref_reset_registry()

g2f_ok <- function(ours, shown) {
  x <- as.numeric(shown)
  dec <- if (grepl(".", shown, fixed = TRUE)) nchar(sub("^[^.]*\\.", "", shown)) else 0
  half <- 0.5 * 10^-dec
  is.finite(ours) && ours >= x - half - 8 * .Machine$double.eps * max(1, abs(x)) &&
    ours <= x + half + 8 * .Machine$double.eps * max(1, abs(x))
}

# Compare a result with displayed strings keyed by short names.
g2f_check <- function(res, shown, record) {
  expect_identical(res$status$code, "ok", info = record)
  q <- res$quantities
  gsv <- function(nm) if (!is.null(res$group_sizes) && nm %in% names(res$group_sizes)) res$group_sizes[[nm]] else NA
  vals <- list(lambda = q$ncp, crit = q$critical, df1 = q$df[1], df2 = q$df[2],
               n = res$solved$n_total %||% q$n_total, n1 = gsv("n1"), n2 = gsv("n2"),
               power = res$achieved$power, alpha = res$solved$alpha, beta = res$achieved$beta,
               f = res$solved$f, delta = res$solved$delta)
  for (nm in names(shown)) {
    if (nm %in% c("df1", "df2", "n", "n1", "n2")) {
      expect_equal(as.double(vals[[nm]]), as.double(shown[[nm]]), info = paste(record, nm))
    } else {
      expect_true(g2f_ok(vals[[nm]], shown[[nm]]), info = sprintf("%s %s: G*Power %s, engine %.10g", record, nm,
                                                                  shown[[nm]], vals[[nm]]))
    }
  }
}

g2f <- function(id, mode, ...) calculate_power(analysis_spec(id, mode, inputs = list(...)))

test_that("one-way ANOVA reproduces the captures in all five modes", {
  g2f_check(g2f("f_anova_oneway", "a_priori", f = 0.25, alpha = 0.05, power = 0.95, n_groups = 10L),
            list(lambda = "24.3750000", crit = "1.9045377", df1 = "9", df2 = "380", n = "390", power = "0.9523634"),
            "GPW-F-ANOVA1-APR-B001-O01--r01")
  g2f_check(g2f("f_anova_oneway", "a_priori", f = 0.02, alpha = 0.001, power = 0.99, n_groups = 50L),
            list(lambda = "79.7200000", crit = "1.7420159", df2 = "199250", n = "199300", power = "0.9900071"),
            "GPW-F-ANOVA1-APR-B001-D01--r01")
  g2f_check(g2f("f_anova_oneway", "a_priori", f = 0.6, alpha = 0.05, power = 0.8, n_groups = 7L),
            list(lambda = "17.6400000", crit = "2.3239938", df2 = "42", n = "49", power = "0.8493675"),
            "SUP-G2-F-ANOVA1-APR-B02--r01")
  g2f_check(g2f("f_anova_oneway", "post_hoc", f = 0.5930904, alpha = 0.05, n_total = 22L, n_groups = 4L),
            list(lambda = "7.7386369", crit = "3.1599076", df1 = "3", df2 = "18", power = "0.5360105"),
            "GPW-F-ANOVA1-POS-B001-O01--r01")
  g2f_check(g2f("f_anova_oneway", "post_hoc", f = 0.02, alpha = 0.0001, n_total = 100000L, n_groups = 50L),
            list(lambda = "40.0000000", crit = "1.9310050", df2 = "99950", power = "0.3456502"),
            "GPW-F-ANOVA1-POS-B001-D01--r01")
  g2f_check(g2f("f_anova_oneway", "sensitivity", alpha = 0.05, power = 0.8, n_total = 120L, n_groups = 4L),
            list(lambda = "11.2765814", crit = "2.6828094", f = "0.3065477"), "GPW-F-ANOVA1-SEN-B001-O01--r01")
  g2f_check(g2f("f_anova_oneway", "sensitivity", alpha = 0.05, power = 0.95, n_total = 3L, n_groups = 2L),
            list(lambda = "624.0359", crit = "161.4476", f = "14.4226201"), "GPW-F-ANOVA1-SEN-B001-E01--r01")
  g2f_check(g2f("f_anova_oneway", "criterion", f = 0.25, power = 0.8, n_total = 120L, n_groups = 4L),
            list(lambda = "7.5000000", crit = "1.7651720", alpha = "0.1577003"), "GPW-F-ANOVA1-CRI-B001-O01--r01")
  g2f_check(g2f("f_anova_oneway", "criterion", f = 0.02, power = 0.9999, n_total = 100000L, n_groups = 50L),
            list(crit = "0.8236806", alpha = "0.8055042"), "GPW-F-ANOVA1-CRI-B001-D01--r01")
  g2f_check(g2f("f_anova_oneway", "compromise", f = 0.25, error_ratio = 1, n_total = 200L, n_groups = 10L),
            list(lambda = "12.5000000", crit = "1.4762105", alpha = "0.1591943", beta = "0.1591943", power = "0.8408057"),
            "GPW-F-ANOVA1-COM-B001-O01--r01")
  g2f_check(g2f("f_anova_oneway", "compromise", f = 0.02, error_ratio = 0.01, n_total = 100000L, n_groups = 50L),
            list(crit = "1.0508134", alpha = "0.3766118", beta = "0.003766118", power = "0.9962339"),
            "GPW-F-ANOVA1-COM-B001-D01--r01")
  g2f_check(g2f("f_anova_oneway", "compromise", f = 1.5, error_ratio = 4, n_total = 6L, n_groups = 2L),
            list(lambda = "13.5000000", crit = "7.4891492", alpha = "0.0520859", beta = "0.2083435", power = "0.7916565"),
            "SUP-G2-F-ANOVA1-COM-B02--r01")
})

test_that("special-effects ANOVA reproduces the captures in all five modes", {
  id <- "f_anova_special_main_interactions"
  g2f_check(g2f(id, "a_priori", f = 0.1, alpha = 0.05, power = 0.95, numerator_df = 8L, n_groups = 30L),
            list(lambda = "22.8300000", crit = "1.9425075", df2 = "2253", n = "2283", power = "0.9500780"),
            "GPW-F-ANOVASP-APR-B001-O01--r01")
  g2f_check(g2f(id, "a_priori", f = 0.02, alpha = 0.001, power = 0.99, numerator_df = 99L, n_groups = 100L),
            list(lambda = "102.6576", crit = "1.4974259", df2 = "256544", n = "256644", power = "0.9900004"),
            "GPW-F-ANOVASP-APR-B001-D01--r01")
  g2f_check(g2f(id, "post_hoc", f = 0.2450722, alpha = 0.05, n_total = 108L, numerator_df = 4L, n_groups = 36L),
            list(lambda = "6.4865214", crit = "2.4989186", df2 = "72", power = "0.4756346"), "GPW-F-ANOVASP-POS-B001-O02--r01")
  g2f_check(g2f(id, "post_hoc", f = 0.25, alpha = 0.05, n_total = 120L, numerator_df = 6L, n_groups = 24L),
            list(lambda = "7.5000000", crit = "2.1945162", df2 = "96", power = "0.4701291"), "SUP-G2-F-ANOVASP-POS-B02--r01")
  g2f_check(g2f(id, "sensitivity", alpha = 0.001, power = 0.99, n_total = 100000L, numerator_df = 99L, n_groups = 100L),
            list(lambda = "102.7037", crit = "1.4976602", f = "0.0320474"), "GPW-F-ANOVASP-SEN-B001-D01--r01")
  g2f_check(g2f(id, "criterion", f = 0.25, power = 0.8, n_total = 120L, numerator_df = 2L, n_groups = 6L),
            list(crit = "2.2210137", df2 = "114", alpha = "0.1131733"), "GPW-F-ANOVASP-CRI-B001-O01--r01")
  g2f_check(g2f(id, "compromise", f = 1.5, error_ratio = 4, n_total = 3L, numerator_df = 1L, n_groups = 2L),
            list(lambda = "6.7500000", crit = "19.5654966", alpha = "0.1415450", beta = "0.5661801", power = "0.4338199"),
            "GPW-F-ANOVASP-COM-B001-E01--r01")
})

test_that("ANCOVA reproduces the captures in all five modes", {
  id <- "f_ancova_main_interactions"
  g2f_check(g2f(id, "a_priori", f = 0.25, alpha = 0.05, power = 0.8, numerator_df = 2L, n_groups = 6L, n_covariates = 1L),
            list(lambda = "9.8750000", crit = "3.0559594", df2 = "151", n = "158", power = "0.8019839"),
            "GPW-F-ANCOVA-APR-B001-O01--r01")
  g2f_check(g2f(id, "a_priori", f = 0.8, alpha = 0.05, power = 0.8, numerator_df = 1L, n_groups = 2L, n_covariates = 2L),
            list(lambda = "9.6000000", crit = "4.8443357", df2 = "11", n = "15", power = "0.8050166"),
            "SUP-G2-F-ANCOVA-APR-B02--r01")
  g2f_check(g2f(id, "post_hoc", f = 0.25, alpha = 0.05, n_total = 120L, numerator_df = 2L, n_groups = 6L, n_covariates = 1L),
            list(lambda = "7.5000000", crit = "3.0765743", df2 = "113", power = "0.6752487"), "GPW-F-ANCOVA-POS-B001-O01--r01")
  g2f_check(g2f(id, "post_hoc", f = 1.5, alpha = 0.05, n_total = 8L, numerator_df = 1L, n_groups = 2L, n_covariates = 2L),
            list(lambda = "18.0000000", crit = "7.7086474", df2 = "4", power = "0.8802120"), "SUP-G2-F-ANCOVA-POS-B03--r01")
  g2f_check(g2f(id, "sensitivity", alpha = 0.05, power = 0.8, n_total = 120L, numerator_df = 2L, n_groups = 6L, n_covariates = 1L),
            list(lambda = "9.8947514", crit = "3.0765743", f = "0.2871520"), "GPW-F-ANCOVA-SEN-B001-O01--r01")
  g2f_check(g2f(id, "criterion", f = 1.5, power = 0.95, n_total = 4L, numerator_df = 1L, n_groups = 2L, n_covariates = 0L),
            list(lambda = "9.0000000", crit = "1.4992924", df2 = "2", alpha = "0.3454346"), "GPW-F-ANCOVA-CRI-B001-E01--r01")
  g2f_check(g2f(id, "compromise", f = 0.02, error_ratio = 0.01, n_total = 100000L, numerator_df = 99L, n_groups = 100L,
                n_covariates = 10L),
            list(crit = "0.9658324", df2 = "99890", alpha = "0.5775379", beta = "0.005775379", power = "0.9942246"),
            "GPW-F-ANCOVA-COM-B001-D01--r01")
  # G*Power refuses this a priori problem (reference defect, see the comparison summary); its own
  # post hoc values (SUP-G2-F-ANCOVA-POS-B03/B04) show that N = 8 is the answer.
  r <- g2f(id, "a_priori", f = 1.5, alpha = 0.05, power = 0.8, numerator_df = 1L, n_groups = 2L, n_covariates = 2L)
  expect_equal(r$solved$n_total, 8)
  r7 <- g2f(id, "post_hoc", f = 1.5, alpha = 0.05, n_total = 7L, numerator_df = 1L, n_groups = 2L, n_covariates = 2L)
  expect_true(g2f_ok(r7$achieved$power, "0.7490785"))
})

test_that("Hotelling one group reproduces the captures in all five modes", {
  id <- "f_hotelling_t2_one_group"
  g2f_check(g2f(id, "a_priori", delta = 0.5, alpha = 0.05, power = 0.8, n_response = 3L),
            list(lambda = "12.0000000", crit = "2.8115435", df1 = "3", df2 = "45", n = "48", power = "0.8038167"),
            "GPW-F-HOT1-APR-B001-O01--r01")
  g2f_check(g2f(id, "a_priori", delta = 0.02, alpha = 0.001, power = 0.99, n_response = 50L),
            list(lambda = "80.2692000", crit = "1.7333833", df2 = "200623", n = "200673", power = "0.9900003"),
            "GPW-F-HOT1-APR-B001-D01--r01")
  g2f_check(g2f(id, "post_hoc", delta = 2.5, alpha = 0.05, n_total = 3L, n_response = 2L),
            list(lambda = "18.7500000", crit = "199.5", df2 = "1", power = "0.1759758"), "GPW-F-HOT1-POS-B001-E01--r01")
  g2f_check(g2f(id, "sensitivity", alpha = 0.05, power = 0.95, n_total = 3L, n_response = 2L),
            list(lambda = "1535.582", crit = "199.5", delta = "22.6243376"), "GPW-F-HOT1-SEN-B001-E01--r01")
  g2f_check(g2f(id, "criterion", delta = 0.5, power = 0.8, n_total = 40L, n_response = 3L),
            list(crit = "2.3406399", df2 = "37", alpha = "0.0891064"), "GPW-F-HOT1-CRI-B001-O01--r01")
  g2f_check(g2f(id, "compromise", delta = 0.5, error_ratio = 1, n_total = 40L, n_response = 3L),
            list(lambda = "10.0000000", crit = "1.9474633", alpha = "0.1388679", beta = "0.1388679", power = "0.8611321"),
            "GPW-F-HOT1-COM-B001-O01--r01")
})

test_that("Hotelling two groups reproduces the captures in all five modes", {
  id <- "f_hotelling_t2_two_groups"
  g2f_check(g2f(id, "a_priori", delta = 0.5, alpha = 0.05, power = 0.8, allocation_ratio = 1, n_response = 3L),
            list(lambda = "11.2500000", crit = "2.6559389", df2 = "176", n1 = "90", n2 = "90", power = "0.8039869"),
            "GPW-F-HOT2-APR-B001-O01--r01")
  g2f_check(g2f(id, "a_priori", delta = 2.5, alpha = 0.05, power = 0.8, allocation_ratio = 4, n_response = 2L),
            list(lambda = "14.7321429", crit = "3.9822980", df2 = "11", n1 = "3", n2 = "11", power = "0.8543089"),
            "GPW-F-HOT2-APR-B001-E01--r01")
  g2f_check(g2f(id, "a_priori", delta = 0.6, alpha = 0.05, power = 0.8, allocation_ratio = 3, n_response = 2L),
            list(lambda = "9.9443836", crit = "3.0593761", df2 = "143", n1 = "37", n2 = "109", power = "0.8044430"),
            "SUP-G2-F-HOT2-APR-A01--r01")
  g2f_check(g2f(id, "a_priori", delta = 0.5, alpha = 0.05, power = 0.8, allocation_ratio = 1.01, n_response = 3L),
            list(lambda = "11.1871508", crit = "2.6562339", df2 = "175", n1 = "89", n2 = "90", power = "0.8014937"),
            "SUP-G2-F-HOT2-APR-A08--r01")
  # GPW-F-HOT2-APR-B001-D01: G*Power's design (n1 1203811, n2 240762, actual power 0.9900000, below
  # the target 0.99 at full precision) is reproduced in the diagnostics; the engine reports the next
  # design of the same split that reaches the target (documented difference).
  r <- g2f(id, "a_priori", delta = 0.02, alpha = 0.001, power = 0.99, allocation_ratio = 0.2, n_response = 50L)
  expect_identical(c(r$diagnostics$gpower_rule$n1, r$diagnostics$gpower_rule$n2), c(1203811, 240762))
  expect_true(g2f_ok(r$diagnostics$gpower_rule$power_integer_design, "0.9900000"))
  expect_lt(r$diagnostics$gpower_rule$power_integer_design, 0.99)
  expect_gte(r$achieved$power, 0.99)
  g2f_check(g2f(id, "post_hoc", delta = 0.02, alpha = 0.001, n1 = 1203811L, n2 = 240762L, n_response = 50L),
            list(lambda = "80.2540111", crit = "1.7332395", df2 = "1444522", power = "0.9900000"),
            "GPW-F-HOT2-APR-B001-D01--r01 (G*Power's design as a post hoc input)")
  g2f_check(g2f(id, "post_hoc", delta = 0.02, alpha = 0.0001, n1 = 50000L, n2 = 50000L, n_response = 50L),
            list(lambda = "10.0000000", crit = "1.9198356", df2 = "99949", power = "0.003977601"),
            "GPW-F-HOT2-POS-B001-D01--r01")
  g2f_check(g2f(id, "sensitivity", alpha = 0.05, power = 0.8, n1 = 30L, n2 = 30L, n_response = 3L),
            list(lambda = "11.6927766", crit = "2.7694309", delta = "0.8829034"), "GPW-F-HOT2-SEN-B001-O01--r01")
  g2f_check(g2f(id, "criterion", delta = 0.5, power = 0.8, n1 = 30L, n2 = 30L, n_response = 3L),
            list(crit = "0.9386598", df2 = "56", alpha = "0.4281773"), "GPW-F-HOT2-CRI-B001-O01--r01")
  g2f_check(g2f(id, "compromise", delta = 2.5, error_ratio = 4, n1 = 2L, n2 = 2L, n_response = 2L),
            list(lambda = "6.2500000", crit = "18.0828911", alpha = "0.1640319", beta = "0.6561277", power = "0.3438723"),
            "GPW-F-HOT2-COM-B001-E01--r01")
})

test_that("G*Power refusals are matched by typed refusals without numbers", {
  code <- function(expr) th_condition_code(expr)
  # GPW-F-ANOVA1-APR-B001-I03 'Power must be greater than or equal to alpha.'
  expect_identical(code(analysis_spec("f_anova_oneway", "a_priori", inputs = list(f = 0.25, alpha = 0.05, power = 0.04,
                                                                                    n_groups = 4L))), "invalid_input")
  # GPW-F-HOT2-APR-B001-I02 'Effect size f² is outside the admissible range (1e-06, 50).'
  expect_identical(code(analysis_spec("f_hotelling_t2_two_groups", "a_priori",
                                      inputs = list(delta = 0, alpha = 0.05, power = 0.8, allocation_ratio = 1,
                                                    n_response = 3L))), "invalid_input")
  # SUP-G2-F-ANCOVA-POS-B01 'The number of covariate must lie between 0 and 0.'
  expect_identical(code(analysis_spec("f_ancova_main_interactions", "post_hoc",
                                      inputs = list(f = 0.5, alpha = 0.05, n_total = 7L, numerator_df = 2L, n_groups = 6L,
                                                    n_covariates = 1L))), "invalid_input")
  # GPW-F-HOT1-POS-B001-I03 'Sample size must be greater than the number of response variables (3) + 1.'
  expect_identical(code(analysis_spec("f_hotelling_t2_one_group", "post_hoc",
                                      inputs = list(delta = 0.5, alpha = 0.05, n_total = 0L, n_response = 3L))), "invalid_input")
})
