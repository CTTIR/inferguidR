# Benchmark cases for the power reference engines (development only).
#
# Inputs come from the original function's help-page example or from a G*Power
# manual example (the documented_example records of the reference corpus, cited
# by record id). A case derived from either (for example the other tail) says
# so. No input is made up for a benchmark.
#
# A case is list(case_id, source, mode, inputs, options, fields).

power_case <- function(case_id, source, mode, inputs, options = list(), fields = NULL) {
  if (is.null(fields)) fields <- if (identical(mode, "a_priori")) c("status", "n_total", "power") else c("status", "power")
  list(case_id = case_id, source = source, mode = mode, inputs = inputs, options = options, fields = fields)
}
both_tails <- function(case_id, source, mode, inputs, given = "two") {
  other <- if (identical(given, "two")) "one" else "two"
  list(power_case(paste0(case_id, "-", given), source, mode, inputs, list(tails = given)),
    power_case(paste0(case_id, "-", other), paste(source, "with the other tail (derived)"), mode, inputs,
      list(tails = other)))
}
manual <- function(id) paste("G*Power 3.1 manual example", id)

eta_to_f <- function(eta2) sqrt(eta2 / (1 - eta2))
pwrss_mixed <- "pwrss ?power.f.mixed.anova example"

power_benchmark_cases <- list(
  f_anova_oneway = function() list(
    power_case("doc055", manual("DOC-055"), "a_priori", list(f = 0.25, alpha = 0.05, power = 0.95, n_groups = 10)),
    power_case("doc080", manual("DOC-080"), "post_hoc", list(f = 0.5930904, alpha = 0.05, n_total = 22, n_groups = 4)),
    power_case("webpower-help", "WebPower ?wp.anova example wp.anova(f = 0.25, k = 4, n = 100)", "post_hoc",
      list(f = 0.25, alpha = 0.05, n_total = 100, n_groups = 4)),
    power_case("webpower-help-n101", "WebPower ?wp.anova example with n = 101, not a multiple of k (derived)",
      "post_hoc", list(f = 0.25, alpha = 0.05, n_total = 101, n_groups = 4))),
  f_anova_special_main_interactions = function() list(
    power_case("doc061", manual("DOC-061"), "post_hoc",
      list(f = 0.7066856, alpha = 0.05, n_total = 108, numerator_df = 2, n_groups = 36)),
    power_case("doc066", manual("DOC-066"), "a_priori",
      list(f = 0.1, alpha = 0.05, power = 0.95, numerator_df = 8, n_groups = 30)),
    power_case("doc069", manual("DOC-069"), "post_hoc",
      list(f = 0.475164, alpha = 0.05, n_total = 20, numerator_df = 1, n_groups = 4)),
    power_case("webpower-help", "WebPower ?wp.kanova example wp.kanova(n = 360, ndf = 4, f = 0.3, ng = 18)",
      "post_hoc", list(f = 0.3, alpha = 0.05, n_total = 360, numerator_df = 4, n_groups = 18))),
  f_ancova_main_interactions = function() list(
    power_case("pwrss-help", paste("pwrss ?power.f.ancova example (eta^2 = 0.059, two levels, power 0.8)",
      "with one covariate and f = sqrt(eta^2 / (1 - eta^2)) (derived)"), "a_priori",
      list(f = eta_to_f(0.059), alpha = 0.05, power = 0.8, numerator_df = 1, n_groups = 2, n_covariates = 1)),
    power_case("pwrss-help-posthoc", "the same at N = 128 (derived)", "post_hoc",
      list(f = eta_to_f(0.059), alpha = 0.05, n_total = 128, numerator_df = 1, n_groups = 2, n_covariates = 1))),
  f_rm_anova_between = function() list(
    power_case("pwrss-help", paste(pwrss_mixed, "(eta^2 = 0.038, 2 between x 2 within, rho 0.5), f from eta^2 (derived)"),
      "a_priori", list(f = eta_to_f(0.038), alpha = 0.05, power = 0.8, k = 2, m = 2, rho = 0.5),
      list(effect_size_specification = "gpower30"))),
  f_rm_anova_within = function() list(
    power_case("pwrss-help", paste(pwrss_mixed, "(eta^2 = 0.022, 1 group x 2 within, rho 0.5), f from eta^2 (derived)"),
      "a_priori", list(f = eta_to_f(0.022), alpha = 0.05, power = 0.8, k = 1, m = 2, rho = 0.5, epsilon = 1),
      list(effect_size_specification = "gpower30")),
    power_case("pwrss-help-2x2", paste(pwrss_mixed, "(eta^2 = 0.01, 2 x 2, rho 0.5, within), f from eta^2 (derived)"),
      "a_priori", list(f = eta_to_f(0.01), alpha = 0.05, power = 0.8, k = 2, m = 2, rho = 0.5, epsilon = 1),
      list(effect_size_specification = "gpower30"))),
  f_rm_anova_within_between = function() list(
    power_case("pwrss-help", paste(pwrss_mixed, "(eta^2 = 0.01, 2 x 2, rho 0.5, interaction), f from eta^2 (derived)"),
      "a_priori", list(f = eta_to_f(0.01), alpha = 0.05, power = 0.8, k = 2, m = 2, rho = 0.5, epsilon = 1),
      list(effect_size_specification = "gpower30"))),
  f_mrc_r2_deviation_from_zero = function() list(
    power_case("doc078", manual("DOC-078"), "post_hoc", list(f2 = 0.1111111, alpha = 0.05, n_total = 95, n_predictors = 5)),
    power_case("pwr-help", "pwr ?pwr.f2.test example pwr.f2.test(u = 5, v = 89, f2 = 0.1 / 0.9)", "post_hoc",
      list(f2 = 0.1 / 0.9, alpha = 0.05, n_total = 95, n_predictors = 5)),
    power_case("doc078-apriori", paste(manual("DOC-078"), "a priori for power 0.8 (derived)"), "a_priori",
      list(f2 = 0.1111111, alpha = 0.05, power = 0.8, n_predictors = 5))),
  f_mrc_r2_increase = function() list(
    power_case("doc086", manual("DOC-086"), "post_hoc",
      list(f2 = 0.0714286, alpha = 0.01, n_total = 90, n_tested = 4, n_predictors = 9)),
    power_case("doc087", manual("DOC-087"), "a_priori",
      list(f2 = 0.0714286, alpha = 0.01, power = 0.8, n_tested = 4, n_predictors = 9)),
    power_case("doc089", manual("DOC-089"), "post_hoc",
      list(f2 = 0.075, alpha = 0.01, n_total = 200, n_tested = 3, n_predictors = 12))),
  f_generic = function() list(
    power_case("pwrss-help", "pwrss ?power.f.test example power.f.test(ncp = 1, df1 = 4, df2 = 100)", "post_hoc",
      list(lambda = 1, alpha = 0.05, df1 = 4, df2 = 100))),
  t_means_one_sample = function() c(
    both_tails("doc120", manual("DOC-120"), "a_priori", list(d = 0.625, alpha = 0.05, power = 0.95), "one"),
    list(power_case("doc121-two", manual("DOC-121"), "a_priori", list(d = 0.1, alpha = 0.01, power = 0.9),
      list(tails = "two"))),
    both_tails("doc120-posthoc", paste(manual("DOC-120"), "at its solved N = 30 (derived)"), "post_hoc",
      list(d = 0.625, alpha = 0.05, n_total = 30), "one")),
  t_means_matched_pairs = function() c(
    both_tails("doc104", manual("DOC-104"), "post_hoc", list(dz = 0.421637, alpha = 0.05, n_total = 50)),
    list(power_case("doc106-two", manual("DOC-106"), "post_hoc", list(dz = 0.2828427, alpha = 0.05, n_total = 50),
      list(tails = "two")))),
  t_means_two_independent = function() c(
    both_tails("stats-help-n20", "stats ?power.t.test example power.t.test(n = 20, delta = 1)", "post_hoc",
      list(d = 1, alpha = 0.05, n1 = 20, n2 = 20)),
    both_tails("stats-help-power90", "stats ?power.t.test example power.t.test(power = .90, delta = 1)", "a_priori",
      list(d = 1, alpha = 0.05, power = 0.9, allocation_ratio = 1)),
    both_tails("pwr-help-90-60", "pwr ?pwr.t2n.test example pwr.t2n.test(d = 0.6, n1 = 90, n2 = 60)", "post_hoc",
      list(d = 0.6, alpha = 0.05, n1 = 90, n2 = 60), "one"),
    list(power_case("doc001-one", manual("DOC-001"), "a_priori",
        list(d = 0.5, alpha = 0.05, power = 0.95, allocation_ratio = 1), list(tails = "one")),
      power_case("doc003-one", manual("DOC-003"), "a_priori",
        list(d = 2, alpha = 0.05, power = 0.95, allocation_ratio = 1), list(tails = "one")))),
  t_correlation_point_biserial = function() c(
    list(power_case("doc004-one", manual("DOC-004"), "a_priori", list(rho = 0.3, alpha = 0.05, power = 0.95),
        list(tails = "one")),
      power_case("doc097-one", manual("DOC-097"), "a_priori", list(rho = 0.25, alpha = 0.05, power = 0.95),
        list(tails = "one"))),
    both_tails("doc083", manual("DOC-083"), "post_hoc", list(rho = 0.5, alpha = 0.05, n_total = 12))),
  t_linreg_one_group_slope = function()
    both_tails("doc073", manual("DOC-073"), "post_hoc",
      list(slope_h1 = -0.0667, slope_h0 = 0, sd_x = 7.5, sd_y = 4, alpha = 0.05, n_total = 100)),
  t_mrc_single_coefficient = function() list(
    power_case("webpower-help-two", paste("WebPower ?wp.regression example wp.regression(n = 100, p1 = 3, f2 = 0.1)",
      "as the test of one of the three coefficients (derived)"), "post_hoc",
      list(f2 = 0.1, alpha = 0.05, n_total = 100, n_predictors = 3), list(tails = "two")),
    power_case("webpower-help-apriori-two", paste("WebPower ?wp.regression example (one of three coefficients),",
      "a priori for power 0.8 (derived)"), "a_priori", list(f2 = 0.1, alpha = 0.05, power = 0.8, n_predictors = 3),
      list(tails = "two"))),
  t_generic = function()
    both_tails("pwrss-help-ncp196", "pwrss ?power.t.test example power.t.test(ncp = 1.96, df = 100)", "post_hoc",
      list(delta = 1.96, alpha = 0.05, df = 100))
)
