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

glm_options <- function(procedure, vc, x, tails = "two", mode_input = NULL) c(list(tails = tails,
  computation = "large_sample", test_procedure = procedure), if (identical(procedure, "demidenko"))
  list(variance_correction = vc), list(x_distribution = x), if (!is.null(mode_input))
  list(effect_size_input_mode = mode_input))

mess_binom <- "MESS ?power_binom_test example power_binom_test(n = 50, p0 = .50, pa = .75)"
pwrss_oneprop <- "pwrss ?power.exact.oneprop example (prob = 0.45, null.prob = 0.50, n = 500)"
exact_help <- "Exact ?power.exact.test example (p1 = 0.15, p2 = 0.60, n1 = 15, n2 = 30)"

pwrss_np <- "pwrss ?power.np.wilcoxon example power.np.wilcoxon(d = 0.25, power = 0.80)"
wilcoxon_cases <- function(effect, design_inputs, posthoc_inputs) c(
  both_tails("pwrss-help", pwrss_np, "a_priori", c(stats::setNames(list(0.25), effect), list(alpha = 0.05, power = 0.8),
    design_inputs)),
  lapply(c("laplace", "logistic"), function(parent) power_case(paste0("pwrss-help-", parent),
    paste(pwrss_np, "with a", parent, "parent, post hoc (derived)"), "post_hoc",
    c(stats::setNames(list(0.25), effect), list(alpha = 0.05), posthoc_inputs),
    list(tails = "two", parent_distribution = parent))))

power_benchmark_cases <- list(
  t_wilcoxon_signed_rank_one_sample = function() wilcoxon_cases("d", list(), list(n_total = 120)),
  t_wilcoxon_signed_rank_matched_pairs = function() wilcoxon_cases("dz", list(), list(n_total = 120)),
  t_wilcoxon_mann_whitney = function() wilcoxon_cases("d", list(allocation_ratio = 1), list(n1 = 200, n2 = 300)),
  exact_proportion_one_sample_binomial = function() c(
    both_tails("mess-help", mess_binom, "post_hoc", list(g = 0.25, p0 = 0.5, alpha = 0.05, n_total = 50)),
    list(power_case("mess-help-minor-tail", paste(mess_binom, "with the minor-tail balancing rule (derived)"), "post_hoc",
        list(g = 0.25, p0 = 0.5, alpha = 0.05, n_total = 50),
        list(tails = "two", alpha_balancing = "minor_tail_half_rest_major"))),
    both_tails("pwrss-help", pwrss_oneprop, "post_hoc", list(g = -0.05, p0 = 0.5, alpha = 0.05, n_total = 500))),
  exact_proportion_sign_test = function()
    both_tails("mess-help", paste(mess_binom, "as a sign test (derived)"), "post_hoc",
      list(g = 0.25, alpha = 0.05, n_total = 50)),
  exact_generic_binomial = function() c(
    both_tails("mess-help", mess_binom, "post_hoc", list(p2 = 0.75, p1 = 0.5, alpha = 0.05, n_total = 50)),
    both_tails("pwrss-help", pwrss_oneprop, "post_hoc", list(p2 = 0.45, p1 = 0.5, alpha = 0.05, n_total = 500))),
  exact_proportions_mcnemar = function() list(
    power_case("mess-help-two", paste("MESS ?power_mcnemar_test example (paid = .125, psi = 2), i.e. pd = 0.375 and",
      "OR = 2, post hoc at n = 120 (derived)"), "post_hoc", list(or = 2, pd = 0.375, alpha = 0.05, n_total = 120),
      list(tails = "two"))),
  exact_proportions_fisher = function()
    both_tails("exact-help", paste(exact_help, "with method = \"Fisher\""), "post_hoc",
      list(p1 = 0.15, p2 = 0.6, alpha = 0.05, n1 = 15, n2 = 30)),
  exact_correlation_bivariate_normal = function() list(
    power_case("pwrss-help-exact", "pwrss ?power.exact.onecor example (rho = 0.20) post hoc at n = 200 (derived)",
      "post_hoc", list(rho_h1 = 0.2, rho_h0 = 0, alpha = 0.05, n_total = 200), list(tails = "two", computation = "exact")),
    power_case("pwrss-help-exact-one", "the same, one-tailed (derived)", "post_hoc",
      list(rho_h1 = 0.2, rho_h0 = 0, alpha = 0.05, n_total = 200), list(tails = "one", computation = "exact")),
    power_case("pwrss-help-fisher-z", "the same with the Fisher z approximation (derived)", "post_hoc",
      list(rho_h1 = 0.2, rho_h0 = 0, alpha = 0.05, n_total = 200), list(tails = "two", computation = "fisher_z_approx"))),
  chi2_goodness_of_fit_contingency = function() list(
    power_case("pwr-help-1", "pwr ?pwr.chisq.test example pwr.chisq.test(w = 0.289, df = 3, N = 100)", "post_hoc",
      list(w = 0.289, df = 3, alpha = 0.05, n_total = 100)),
    power_case("pwr-help-2", "pwr ?pwr.chisq.test example pwr.chisq.test(w = 0.346, df = 2, N = 140, sig.level = 0.01)",
      "post_hoc", list(w = 0.346, df = 2, alpha = 0.01, n_total = 140)),
    power_case("pwr-help-1-apriori", "the first example a priori for power 0.8 (derived)", "a_priori",
      list(w = 0.289, df = 3, alpha = 0.05, power = 0.8))),
  chi2_generic = function() list(
    power_case("pwrss-help", "pwrss ?power.chisq.test example power.chisq.test(ncp = 20, df = 100)", "post_hoc",
      list(lambda = 20, df = 100, alpha = 0.05))),
  boschloo_proportions_two_independent = function() lapply(c("minlike", "less", "greater"), function(v)
    power_case(paste0("exact-help-", v), paste(exact_help, "with method = \"Boschloo\",", v), "post_hoc",
      list(p1 = 0.15, p2 = 0.6, alpha = 0.05, n1 = 15, n2 = 30), list(variant = v))),
  z_proportions_two_independent = function() c(
    both_tails("stats-help-n50", "stats ?power.prop.test example power.prop.test(n = 50, p1 = .50, p2 = .75)",
      "post_hoc", list(p1 = 0.5, p2 = 0.75, alpha = 0.05, n1 = 50, n2 = 50)),
    list(power_case("stats-help-power90", "stats ?power.prop.test example power.prop.test(p1 = .50, p2 = .75, power = .90)",
        "a_priori", list(p1 = 0.5, p2 = 0.75, alpha = 0.05, power = 0.9, allocation_ratio = 1), list(tails = "two")),
      power_case("mess-help-ratio2", paste("MESS ?power_prop_test example (p1 = .65, p2 = .85, ratio = 2) at n1 = 50",
        "(derived)"), "post_hoc", list(p1 = 0.65, p2 = 0.85, alpha = 0.05, n1 = 50, n2 = 100), list(tails = "two")),
      power_case("pwr-help-h030", paste("pwr ?pwr.2p2n.test example pwr.2p2n.test(h = 0.30, n1 = 80, n2 = 245,",
        "alternative = \"greater\") with p1 = 0.5 (derived)"), "post_hoc",
        list(p1 = 0.5, h = 0.3, alpha = 0.05, n1 = 80, n2 = 245),
        list(tails = "one", arcsin_transform = "on", input = "cohens_h")))),
  z_correlations_two_independent = function() list(
    power_case("doc140", manual("DOC-140"), "post_hoc", list(q = -0.4028126, alpha = 0.05, n1 = 260, n2 = 51),
      list(tails = "two")),
    power_case("doc140-one", paste(manual("DOC-140"), "one-tailed (derived)"), "post_hoc",
      list(q = -0.4028126, alpha = 0.05, n1 = 260, n2 = 51), list(tails = "one")),
    power_case("doc141", manual("DOC-141"), "a_priori", list(q = -0.4028126, alpha = 0.05, power = 0.726352,
      allocation_ratio = 1), list(tails = "two"))),
  z_correlations_two_dependent_common_index = function() list(
    power_case("doc146", manual("DOC-146"), "a_priori", list(rho_ac = 0.2, rho_ab = 0.4, rho_bc = 0.5, alpha = 0.05,
      power = 0.8), list(tails = "one"))),
  z_correlations_two_dependent_no_common_index = function() list(
    power_case("doc143", manual("DOC-143"), "a_priori", list(rho_cd = 0.2, rho_ab = 0.1, rho_ac = 0.5, rho_ad = 0.4,
      rho_bc = -0.4, rho_bd = 0.8, alpha = 0.05, power = 0.8), list(tails = "one")),
    power_case("doc144", manual("DOC-144"), "a_priori", list(rho_cd = 0.2, rho_ab = 0.1, rho_ac = 0, rho_ad = 0,
      rho_bc = 0, rho_bd = 0, alpha = 0.05, power = 0.8), list(tails = "one"))),
  z_logistic_regression = function() list(
    power_case("doc150-hsieh", manual("DOC-150"), "a_priori", list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05,
      power = 0.95, r2_other_x = 0, x_mu = 0, x_sigma = 1), glm_options("hsieh", "on", "normal")),

    power_case("doc152-demidenko", manual("DOC-152"), "a_priori", list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05,
      power = 0.95, r2_other_x = 0, x_mu = 0, x_sigma = 1), glm_options("demidenko", "off", "normal")),

    power_case("doc156-binary", manual("DOC-156"), "a_priori", list(p2 = 0.1, p1 = 0.05, alpha = 0.05, power = 0.95,
      r2_other_x = 0, x_pi = 0.5), glm_options("demidenko", "off", "binomial", mode_input = "two_probabilities")),
    power_case("doc152-posthoc-one", paste(manual("DOC-152"), "post hoc at N = 300, one-tailed (derived)"), "post_hoc",
      list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, n_total = 300, r2_other_x = 0, x_mu = 0, x_sigma = 1),
      glm_options("demidenko", "off", "normal", tails = "one"))),
  z_poisson_regression = function() list(
    power_case("doc163-demidenko", manual("DOC-163"), "a_priori", list(exp_beta1 = 1.3, base_rate = 0.85,
      mean_exposure = 1, alpha = 0.05, power = 0.95, r2_other_x = 0, x_pi = 0.5),
      glm_options("demidenko", "off", "binomial", tails = "one")),
    power_case("webpower-help", paste("WebPower ?wp.poisson example wp.poisson(n = 4406, exp0 = 2.798, exp1 = 0.8938,",
      "family = \"Bernoulli\", parameter = 0.53)"), "post_hoc", list(exp_beta1 = 0.8938, base_rate = 2.798,
      mean_exposure = 1, alpha = 0.05, n_total = 4406, r2_other_x = 0, x_pi = 0.53),
      glm_options("demidenko", "off", "binomial"))),
  z_generic = function()
    both_tails("pwrss-help", "pwrss ?power.z.test example power.z.test(mean = 1.96)", "post_hoc",
      list(mu = 1.96, sigma = 1, alpha = 0.05)),
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
