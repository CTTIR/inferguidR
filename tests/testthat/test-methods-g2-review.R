# G2 review captures (adversarial review 2026-10-05 and 2026-10-06, specs/reviews/G2-anova-ancova-hotelling.md).
#
# Every displayed value is quoted from an emulator_capture fixture in
# validation/fixtures/emulator_capture_supplement/review-G2-anova-ancova-hotelling/ (G*Power 3.1.9.7
# under Wine, D22/D23), named by its record id (R: first round, S: second round). Match rule as in
# test-methods-g2-fixtures.R: integers exact, reals within the closed rounding interval of the
# displayed digits plus 8 eps max(1, |x|); the last displayed digit of a value in exponent notation
# ("1.754034e-06") is scaled by its exponent.

ref_reset_registry()

g2r_ok <- function(ours, shown) {
  x <- as.numeric(shown)
  mant <- sub("[eE].*$", "", shown)
  ex <- if (grepl("[eE]", shown)) as.numeric(sub("^.*[eE]", "", shown)) else 0
  dec <- if (grepl(".", mant, fixed = TRUE)) nchar(sub("^[^.]*\\.", "", mant)) else 0
  is.finite(ours) && abs(ours - x) <= 0.5 * 10^(ex - dec) + 8 * .Machine$double.eps * max(1, abs(x))
}
g2r <- function(id, mode, ...) calculate_power(analysis_spec(id, mode, inputs = list(...)))

test_that("review captures of the ANOVA, ANCOVA and one-group Hotelling engines agree", {
  r <- g2r("f_anova_oneway", "a_priori", f = 0.3, alpha = 0.05, power = 0.9, n_groups = 5L)  # REV-G2-F-ANOVA1-APR-R01
  expect_equal(r$solved$n_total, 180)
  expect_true(g2r_ok(r$achieved$power, "0.9073603"))
  expect_true(g2r_ok(r$quantities$critical, "2.4232862"))
  r <- g2r("f_anova_oneway", "compromise", f = 0.8, error_ratio = 1, n_total = 7L, n_groups = 6L)  # REV-G2-F-ANOVA1-COM-R04
  expect_true(g2r_ok(r$solved$alpha, "0.4316708"))
  expect_true(g2r_ok(r$quantities$critical, "2.6847311"))
  r <- g2r("f_anova_special_main_interactions", "a_priori", f = 0.5, alpha = 0.05, power = 0.8,
           numerator_df = 20L, n_groups = 5L)                                               # REV-G2-F-ANOVASP-APR-R03
  expect_equal(r$solved$n_total, 99)
  expect_true(g2r_ok(r$achieved$power, "0.8054383"))
  r <- g2r("f_anova_special_main_interactions", "compromise", f = 1, error_ratio = 1, n_total = 5L,
           numerator_df = 3L, n_groups = 4L)                                                # REV-G2-F-ANOVASP-COM-R01
  expect_true(g2r_ok(r$solved$alpha, "0.3942563"))
  r <- g2r("f_ancova_main_interactions", "post_hoc", f = 1, alpha = 0.05, n_total = 9L, numerator_df = 2L,
           n_groups = 3L, n_covariates = 5L)                                                # REV-G2-F-ANCOVA-POS-R01
  expect_true(g2r_ok(r$achieved$power, "0.1260480"))
  expect_true(g2r_ok(r$quantities$critical, "199.5"))
  r <- g2r("f_hotelling_t2_one_group", "a_priori", delta = 0.5, alpha = 0.05, power = 0.8, n_response = 100L)  # REV-G2-F-HOT1-APR-R04
  expect_equal(r$solved$n_total, 233)
  expect_true(g2r_ok(r$achieved$power, "0.8015421"))
  r <- g2r("f_hotelling_t2_one_group", "sensitivity", alpha = 0.05, power = 0.8, n_total = 9L, n_response = 8L)  # REV-G2-F-HOT1-SEN-R02
  expect_true(g2r_ok(r$solved$delta, "18.6586829"))
})

test_that("G*Power refuses feasible ANCOVA a priori problems; the engine answers them", {
  # REV-G2-F-ANCOVA-APR-R03 / -R04: 'Cannot compute result! ...'. Smallest N recomputed in the review.
  r <- g2r("f_ancova_main_interactions", "a_priori", f = 1.2, alpha = 0.05, power = 0.8, numerator_df = 1L,
           n_groups = 2L, n_covariates = 2L)
  expect_equal(r$solved$n_total, 9)
  r <- g2r("f_ancova_main_interactions", "a_priori", f = 2, alpha = 0.05, power = 0.8, numerator_df = 1L,
           n_groups = 2L, n_covariates = 2L)
  expect_equal(r$solved$n_total, 7)
})

test_that("Hotelling two groups: the allocation rule predicts G*Power's a priori designs", {
  id <- "f_hotelling_t2_two_groups"
  # REV-G2-F-HOT2-APR-R05: G*Power 54 / 82, actual power 0.8996072 (below the target 0.90)
  r <- g2r(id, "a_priori", delta = 0.7, alpha = 0.05, power = 0.9, allocation_ratio = 1.5, n_response = 4L)
  expect_identical(c(r$diagnostics$gpower_rule$n1, r$diagnostics$gpower_rule$n2), c(54, 82))
  expect_true(g2r_ok(r$diagnostics$gpower_rule$power_integer_design, "0.8996072"))
  expect_equal(r$group_sizes, c(n1 = 55L, n2 = 82L))
  expect_gte(r$achieved$power, 0.9)
  # REV-G2-F-HOT2-APR-R08: r = 1.0001 follows the r != 1 rule (89 / 90)
  r <- g2r(id, "a_priori", delta = 0.5, alpha = 0.05, power = 0.8, allocation_ratio = 1.0001, n_response = 3L)
  expect_equal(r$group_sizes, c(n1 = 89L, n2 = 90L))
  expect_true(g2r_ok(r$achieved$power, "0.8014937"))
  # REV-G2-F-HOT2-APR-R04: G*Power's N* = 155 (62 / 93) is reported although 62 / 92 also reaches 0.80
  r <- g2r(id, "a_priori", delta = 0.62, alpha = 0.05, power = 0.8, allocation_ratio = 1.5, n_response = 6L)
  expect_equal(r$group_sizes, c(n1 = 62L, n2 = 93L))
  expect_true(g2r_ok(r$achieved$power, "0.8028609"))
  p92 <- .ig_power_f_hotelling_two(list(n1 = 62, n2 = 92, delta = 0.62, n_response = 6, alpha = 0.05))$power
  expect_gt(p92, 0.8)
})

test_that("Hotelling two groups with unequal groups: engine uses n1 n2 / N in criterion and compromise", {
  id <- "f_hotelling_t2_two_groups"
  # REV-G2-F-HOT2-POS-R03: G*Power's own post hoc at the alpha its criterion mode solved (0.2521328)
  r <- g2r(id, "post_hoc", delta = 0.9, alpha = 0.2521328, n1 = 25L, n2 = 15L, n_response = 4L)
  expect_true(g2r_ok(r$quantities$ncp, "7.5937500"))
  expect_true(g2r_ok(r$achieved$power, "0.8314709"))
  # REV-G2-F-HOT2-CRI-R01: G*Power shows lambda 8.1 = N Delta^2 / 4 (reference defect); the engine keeps
  # lambda = n1 n2 / N Delta^2 and its alpha reproduces the target power in post hoc mode.
  cr <- g2r(id, "criterion", delta = 0.9, power = 0.85, n1 = 25L, n2 = 15L, n_response = 4L)
  expect_equal(cr$quantities$ncp, 25 * 15 / 40 * 0.81)
  ph <- g2r(id, "post_hoc", delta = 0.9, alpha = cr$solved$alpha, n1 = 25L, n2 = 15L, n_response = 4L)
  expect_equal(ph$achieved$power, 0.85, tolerance = 1e-9)
  # REV-G2-F-HOT2-COM-R01: G*Power displays lambda 6.912 but solves alpha with 7.2 (beta / alpha 2.108);
  # the engine's beta / alpha equals the requested 2.
  co <- g2r(id, "compromise", delta = 1.2, error_ratio = 2, n1 = 8L, n2 = 12L, n_response = 3L)
  expect_true(g2r_ok(co$quantities$ncp, "6.9120000"))
  expect_equal(co$achieved$beta / co$solved$alpha, 2, tolerance = 1e-8)
  # REV-G2-F-HOT2-POS-R04: G*Power's post hoc at its compromise alpha gives its displayed compromise power
  r <- g2r(id, "post_hoc", delta = 1.2, alpha = 0.1427415, n1 = 8L, n2 = 12L, n_response = 3L)
  expect_true(g2r_ok(r$achieved$power, "0.6990409"))
})

# Second round (2026-10-06) -------------------------------------------------------------------------

test_that("second round: every method agrees with G*Power for df2 above 4e5", {
  # stats::qf switches to a chi-square approximation there; the campaign had one such record (HOT2 D01).
  r <- g2r("f_anova_oneway", "post_hoc", f = 0.004, alpha = 0.001, n_total = 1000000L, n_groups = 5L)  # REV-G2-F-ANOVA1-POS-S01
  expect_equal(r$quantities$df, c(4, 999995))
  expect_true(g2r_ok(r$quantities$critical, "4.6167447"))
  expect_true(g2r_ok(r$achieved$power, "0.5270298"))
  r <- g2r("f_anova_oneway", "a_priori", f = 0.003, alpha = 0.01, power = 0.9, n_groups = 4L)          # REV-G2-F-ANOVA1-APR-S01
  expect_equal(r$solved$n_total, 2138612)
  expect_true(g2r_ok(r$quantities$critical, "3.7816314"))
  expect_true(g2r_ok(r$achieved$power, "0.9000006"))
  r <- g2r("f_anova_special_main_interactions", "a_priori", f = 0.004, alpha = 0.05, power = 0.95,
           numerator_df = 3L, n_groups = 8L)                                                         # REV-G2-F-ANOVASP-APR-S01
  expect_equal(r$solved$n_total, 1073123)
  expect_true(g2r_ok(r$quantities$critical, "2.6049176"))
  expect_true(g2r_ok(r$achieved$power, "0.9500001"))
  r <- g2r("f_ancova_main_interactions", "criterion", f = 0.01, power = 0.95, n_total = 450000L,
           numerator_df = 2L, n_groups = 3L, n_covariates = 2L)                                      # REV-G2-F-ANCOVA-CRI-S01
  expect_true(g2r_ok(r$solved$alpha, "1.754034e-06"))
  expect_true(g2r_ok(r$quantities$critical, "13.2539829"))
  r <- g2r("f_hotelling_t2_one_group", "a_priori", delta = 0.005, alpha = 0.05, power = 0.9, n_response = 3L)  # REV-G2-F-HOT1-APR-S01
  expect_equal(r$solved$n_total, 566864)
  expect_true(g2r_ok(r$quantities$critical, "2.6049250"))
  expect_true(g2r_ok(r$achieved$power, "0.9000003"))
  r <- g2r("f_hotelling_t2_two_groups", "post_hoc", delta = 0.01, alpha = 0.0001, n1 = 400000L, n2 = 600000L,
           n_response = 30L)                                                                         # REV-G2-F-HOT2-POS-S01
  expect_true(g2r_ok(r$quantities$critical, "2.2544657"))
  expect_true(g2r_ok(r$achieved$power, "0.1388395"))
  r <- g2r("f_hotelling_t2_two_groups", "a_priori", delta = 0.004, alpha = 0.05, power = 0.8, allocation_ratio = 1,
           n_response = 2L)                                                                          # REV-G2-F-HOT2-APR-S08
  expect_equal(r$group_sizes, c(n1 = 1204338L, n2 = 1204338L))
  expect_true(g2r_ok(r$achieved$power, "0.8000001"))
})

test_that("second round: a priori with target power equal to alpha gives the smallest design", {
  r <- g2r("f_anova_oneway", "a_priori", f = 0.25, alpha = 0.05, power = 0.05, n_groups = 4L)          # REV-G2-F-ANOVA1-APR-S02
  expect_equal(r$solved$n_total, 8)
  expect_true(g2r_ok(r$achieved$power, "0.0636852"))
  r <- g2r("f_anova_special_main_interactions", "a_priori", f = 0.3, alpha = 0.05, power = 0.05,
           numerator_df = 2L, n_groups = 6L)                                                         # REV-G2-F-ANOVASP-APR-S02
  expect_equal(r$solved$n_total, 7)
  expect_true(g2r_ok(r$achieved$power, "0.0575592"))
  r <- g2r("f_hotelling_t2_one_group", "a_priori", delta = 0.5, alpha = 0.05, power = 0.05, n_response = 3L)  # REV-G2-F-HOT1-APR-S02
  expect_equal(r$solved$n_total, 4)
  expect_true(g2r_ok(r$quantities$critical, "215.7073"))
  expect_true(g2r_ok(r$achieved$power, "0.0579257"))
  # REV-G2-F-ANCOVA-APR-S02: G*Power refuses the same problem for the ANCOVA ('Cannot compute result!', one
  # covariate); the engine answers the smallest design N = k + c + 1.
  r <- g2r("f_ancova_main_interactions", "a_priori", f = 0.3, alpha = 0.05, power = 0.05, numerator_df = 2L,
           n_groups = 3L, n_covariates = 1L)
  expect_equal(r$solved$n_total, 5)
  expect_gte(r$achieved$power, 0.05)
  # REV-G2-F-ANOVA1-SEN-S01: in sensitivity mode G*Power answers this target with f = 0, lambda 0 and the
  # critical F below; f = 0 has power alpha exactly.
  r <- g2r("f_anova_oneway", "post_hoc", f = 0, alpha = 0.05, n_total = 60L, n_groups = 3L)
  expect_true(g2r_ok(r$quantities$critical, "3.1588427"))
  expect_equal(r$achieved$power, 0.05, tolerance = 1e-12)
})

test_that("second round: G*Power's input limits", {
  code <- function(expr) th_condition_code(expr)
  # 'There must be at least 2 groups.' REV-G2-F-ANOVA1-POS-S05, -ANOVASP-POS-S02, -ANCOVA-POS-S02
  expect_identical(code(analysis_spec("f_anova_oneway", "post_hoc",
                                      inputs = list(f = 0.3, alpha = 0.05, n_total = 20L, n_groups = 1L))), "invalid_input")
  expect_identical(code(analysis_spec("f_anova_special_main_interactions", "post_hoc",
                                      inputs = list(f = 0.3, alpha = 0.05, n_total = 30L, numerator_df = 1L, n_groups = 1L))),
                   "invalid_input")
  expect_identical(code(analysis_spec("f_ancova_main_interactions", "post_hoc",
                                      inputs = list(f = 0.3, alpha = 0.05, n_total = 30L, numerator_df = 1L, n_groups = 1L,
                                                    n_covariates = 1L))), "invalid_input")
  # 'Degrees of Freedom must be greater or equal one.' REV-G2-F-ANOVASP-POS-S03
  expect_identical(code(analysis_spec("f_anova_special_main_interactions", "post_hoc",
                                      inputs = list(f = 0.3, alpha = 0.05, n_total = 30L, numerator_df = 0L, n_groups = 3L))),
                   "invalid_input")
  # 'There must be at least one response variable.' REV-G2-F-HOT1-POS-S04
  expect_identical(code(analysis_spec("f_hotelling_t2_one_group", "post_hoc",
                                      inputs = list(delta = 0.5, alpha = 0.05, n_total = 10L, n_response = 0L))), "invalid_input")
  # limits G*Power accepts: f = 50, alpha = 1e-07, beta / alpha = 100 and 0.01
  r <- g2r("f_anova_oneway", "post_hoc", f = 50, alpha = 0.05, n_total = 10L, n_groups = 2L)           # REV-G2-F-ANOVA1-POS-S02
  expect_true(g2r_ok(r$achieved$power, "1.0000000"))
  r <- g2r("f_anova_oneway", "post_hoc", f = 0.3, alpha = 1e-7, n_total = 90L, n_groups = 3L)          # REV-G2-F-ANOVA1-POS-S04
  expect_true(g2r_ok(r$quantities$critical, "19.5099010"))
  expect_true(g2r_ok(r$achieved$power, "0.001569288"))
  r <- g2r("f_anova_oneway", "compromise", f = 0.3, error_ratio = 100, n_total = 60L, n_groups = 3L)   # REV-G2-F-ANOVA1-COM-S01
  expect_true(g2r_ok(r$solved$alpha, "0.007594647"))
  expect_true(g2r_ok(r$achieved$beta, "0.7594647"))
  r <- g2r("f_anova_oneway", "compromise", f = 0.3, error_ratio = 0.01, n_total = 60L, n_groups = 3L)  # REV-G2-F-ANOVA1-COM-S02
  expect_true(g2r_ok(r$solved$alpha, "0.8874206"))
  expect_true(g2r_ok(r$achieved$beta, "0.008874206"))
  # numerator df above N and 50 covariates are computed
  r <- g2r("f_anova_special_main_interactions", "post_hoc", f = 0.5, alpha = 0.05, n_total = 12L, numerator_df = 40L,
           n_groups = 5L)                                                                            # REV-G2-F-ANOVASP-POS-S05
  expect_true(g2r_ok(r$achieved$power, "0.0603220"))
  r <- g2r("f_ancova_main_interactions", "post_hoc", f = 0.3, alpha = 0.05, n_total = 60L, numerator_df = 2L,
           n_groups = 3L, n_covariates = 50L)                                                        # REV-G2-F-ANCOVA-POS-S03
  expect_equal(r$quantities$df, c(2, 7))
  expect_true(g2r_ok(r$achieved$power, "0.3701362"))
})

test_that("second round: the two-group split is taken in double arithmetic at binary ties", {
  id <- "f_hotelling_t2_two_groups"
  ap <- function(delta, ratio) g2r(id, "a_priori", delta = delta, alpha = 0.05, power = 0.8, allocation_ratio = ratio,
                                   n_response = 2L)
  # N* / (1 + r) is exactly m + 1/2 but one ulp below it as a double: G*Power shows the rounded-down n1.
  r <- ap(3.17, 0.12)                                                                               # REV-G2-F-HOT2-APR-S01
  expect_equal(r$group_sizes, c(n1 = 12L, n2 = 2L))
  expect_true(g2r_ok(r$achieved$power, "0.9054819"))
  r <- ap(1.52, 0.68)                                                                               # REV-G2-F-HOT2-APR-S02
  expect_equal(r$group_sizes, c(n1 = 12L, n2 = 9L))
  expect_true(g2r_ok(r$achieved$power, "0.8157990"))
  r <- ap(1.27, 1.24)                                                                               # REV-G2-F-HOT2-APR-S03
  expect_equal(r$group_sizes, c(n1 = 12L, n2 = 16L))
  expect_true(g2r_ok(r$achieved$power, "0.8065011"))
  r <- ap(1.62, 0.12)                                                                               # REV-G2-F-HOT2-APR-S05
  expect_equal(r$group_sizes, c(n1 = 37L, n2 = 5L))
  expect_true(g2r_ok(r$achieved$power, "0.8419284"))
  # REV-G2-F-HOT2-APR-S04 and -S06: the rounded-down design misses the target (7 / 26, 7 / 59); it is kept in
  # the diagnostics and the engine reports the next design.
  r <- ap(1.37, 3.4)
  expect_identical(c(r$diagnostics$gpower_rule$n1, r$diagnostics$gpower_rule$n2), c(7, 26))
  expect_true(g2r_ok(r$diagnostics$gpower_rule$power_integer_design, "0.7873393"))
  expect_equal(r$group_sizes, c(n1 = 8L, n2 = 26L))
  r <- ap(1.24, 7.8)
  expect_identical(c(r$diagnostics$gpower_rule$n1, r$diagnostics$gpower_rule$n2), c(7, 59))
  expect_true(g2r_ok(r$diagnostics$gpower_rule$power_integer_design, "0.7789251"))
  expect_equal(r$group_sizes, c(n1 = 8L, n2 = 59L))
  # REV-G2-F-HOT2-APR-S07: r = 0.01 leaves group 2 empty in G*Power (14 / 0, power = alpha)
  r <- ap(10, 0.01)
  expect_identical(c(r$diagnostics$gpower_rule$n1, r$diagnostics$gpower_rule$n2), c(14, 0))
  expect_true(all(r$group_sizes >= 1))
  expect_gte(r$achieved$power, 0.8)
})

test_that("second round: equal groups agree in criterion and compromise mode, unequal groups show the defect", {
  id <- "f_hotelling_t2_two_groups"
  cr <- g2r(id, "criterion", delta = 0.7, power = 0.9, n1 = 22L, n2 = 22L, n_response = 5L)            # REV-G2-F-HOT2-CRI-S01
  expect_true(g2r_ok(cr$solved$alpha, "0.5706073"))
  expect_true(g2r_ok(cr$quantities$critical, "0.7796437"))
  # REV-G2-F-HOT2-CRI-S02: for 21 / 23 G*Power displays the 22 / 22 values (lambda 5.39 = N Delta^2 / 4)
  cu <- g2r(id, "criterion", delta = 0.7, power = 0.9, n1 = 21L, n2 = 23L, n_response = 5L)
  expect_equal(cu$quantities$ncp, 21 * 23 / 44 * 0.49)
  expect_false(g2r_ok(cu$solved$alpha, "0.5706073"))
  co <- g2r(id, "compromise", delta = 0.7, error_ratio = 3, n1 = 22L, n2 = 22L, n_response = 5L)       # REV-G2-F-HOT2-COM-S01
  expect_true(g2r_ok(co$solved$alpha, "0.1463953"))
  expect_true(g2r_ok(co$achieved$beta, "0.4391859"))
  # REV-G2-F-HOT2-COM-S02 (15 / 25, q 100): G*Power's alpha 0.008613912 and beta 0.8727750 have the ratio
  # 101.3; its beta and power are those of lambda = n1 n2 / N Delta^2 at that alpha.
  cu <- g2r(id, "compromise", delta = 0.6, error_ratio = 100, n1 = 15L, n2 = 25L, n_response = 2L)
  expect_equal(cu$achieved$beta / cu$solved$alpha, 100, tolerance = 1e-8)
  ph <- g2r(id, "post_hoc", delta = 0.6, alpha = 0.008613912, n1 = 15L, n2 = 25L, n_response = 2L)
  expect_true(g2r_ok(ph$quantities$ncp, "3.3750000"))
  expect_true(g2r_ok(ph$quantities$critical, "5.4211754"))
  expect_true(g2r_ok(ph$achieved$power, "0.1272250"))
  expect_true(g2r_ok(ph$achieved$beta, "0.8727750"))
})

test_that("the two-group a priori search stays inside its search limit", {
  # N* = 209 (52 / 157) misses the target; the next design is 210. With the limit at 209 there is none.
  p <- list(delta = 0.5, alpha = 0.05, power = 0.8, allocation_ratio = 3, n_response = 2)
  at_limit <- .ig_g2_search_hotelling_two(p, 4, 209)
  expect_identical(at_limit$status, "infeasible")
  expect_null(at_limit$n)
  expect_equal(.ig_g2_search_hotelling_two(p, 4, 210)$n, 210)
})

test_that("the two-group a priori diagnostics claim a smallest design only for equal allocation", {
  id <- "f_hotelling_t2_two_groups"
  # REV-G2-F-HOT2-APR-R04: 62 / 93 is G*Power's design although 62 / 92 reaches the target as well
  r <- g2r(id, "a_priori", delta = 0.62, alpha = 0.05, power = 0.8, allocation_ratio = 1.5, n_response = 6L)
  expect_equal(r$group_sizes, c(n1 = 62L, n2 = 93L))
  expect_false(grepl("^proven", r$diagnostics$minimality))
  expect_match(r$diagnostics$minimality, "not claimed to be the smallest design")
  expect_match(r$diagnostics$search_scale, "fractional group sizes")
  # SUP-G2-F-HOT2-APR-A02: the engine moves from G*Power's 52 / 157 to the next design
  r <- g2r(id, "a_priori", delta = 0.5, alpha = 0.05, power = 0.8, allocation_ratio = 3, n_response = 2L)
  expect_match(r$diagnostics$minimality, "first total above G\\*Power's N = 209")
  # r = 1: fractional and integer design coincide, the enumeration proves minimality over even totals
  r <- g2r(id, "a_priori", delta = 0.5, alpha = 0.05, power = 0.8, allocation_ratio = 1, n_response = 3L)
  expect_match(r$diagnostics$minimality, "^proven")
  expect_null(r$diagnostics$search_scale)
})
