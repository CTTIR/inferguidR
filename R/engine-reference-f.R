# Reference forward models for the F-test procedures (noncentrality f^2 N as
# in G*Power). Originals are called with every setting explicit; pwrss takes
# eta^2 = f^2 / (1 + f^2).

.ig_ref_quiet <- function(expr) suppressMessages(suppressWarnings(expr))
.ig_ref_eta2 <- function(f) f^2 / (1 + f^2)

.ig_ref_f_anova_oneway <- function(p) {
  if (p$n_total <= p$n_groups) .ig_ref_outside("WebPower::wp.anova needs more observations than groups")
  fit <- WebPower::wp.anova(k = p$n_groups, n = p$n_total, f = p$f, alpha = p$alpha)
  .ig_ref_result(fit$power, "WebPower::wp.anova", n_total = p$n_total)
}

.ig_ref_f_anova_special <- function(p) {
  if (p$n_total <= p$n_groups) .ig_ref_outside("WebPower::wp.kanova needs more observations than groups")
  fit <- WebPower::wp.kanova(n = p$n_total, ndf = p$numerator_df, f = p$f, ng = p$n_groups, alpha = p$alpha)
  .ig_ref_result(fit$power, "WebPower::wp.kanova", n_total = p$n_total)
}

.ig_ref_f_ancova <- function(p) {
  if (p$n_total <= p$n_groups + p$n_covariates)
    .ig_ref_outside("pwrss::power.f.ancova needs more observations than groups plus covariates")
  fit <- .ig_ref_quiet(pwrss::power.f.ancova(eta.squared = .ig_ref_eta2(p$f), factor.levels = p$n_groups,
    k.covariates = p$n_covariates, n.total = p$n_total, alpha = p$alpha, ceil.n = FALSE, verbose = FALSE))
  .ig_ref_result(fit$power, "pwrss::power.f.ancova(eta.squared = f^2 / (1 + f^2))", n_total = p$n_total)
}

.ig_ref_f_rm_anova <- function(effect) function(p) {
  if (p$n_total <= p$k) .ig_ref_outside("pwrss::power.f.mixed.anova needs more subjects than groups")
  fit <- .ig_ref_quiet(pwrss::power.f.mixed.anova(eta.squared = .ig_ref_eta2(p$f), factor.levels = c(p$k, p$m),
    factor.type = c("between", "within"), rho.within = p$rho, epsilon = if (is.null(p$epsilon)) 1 else p$epsilon,
    n.total = p$n_total, alpha = p$alpha, effect = effect, ceil.n = FALSE, verbose = FALSE))
  .ig_ref_result(fit$power, sprintf("pwrss::power.f.mixed.anova(effect = \"%s\")", effect), n_total = p$n_total)
}

.ig_ref_f_mrc_r2_deviation_from_zero <- function(p) {
  v <- p$n_total - p$n_predictors - 1
  if (v < 1) .ig_ref_outside("pwr::pwr.f2.test needs N > number of predictors + 1")
  fit <- pwr::pwr.f2.test(u = p$n_predictors, v = v, f2 = p$f2, sig.level = p$alpha)
  .ig_ref_result(fit$power, "pwr::pwr.f2.test(u = p, v = N - p - 1)", n_total = p$n_total)
}

.ig_ref_f_mrc_r2_increase <- function(p) {
  if (p$n_total <= p$n_predictors + 1) .ig_ref_outside("WebPower::wp.regression needs N > number of predictors + 1")
  fit <- WebPower::wp.regression(n = p$n_total, p1 = p$n_predictors, p2 = p$n_predictors - p$n_tested, f2 = p$f2,
    alpha = p$alpha)
  .ig_ref_result(fit$power, "WebPower::wp.regression(p1 = p, p2 = p - q)", n_total = p$n_total)
}

.ig_ref_f_generic <- function(p) {
  fit <- .ig_ref_quiet(pwrss::power.f.test(ncp = p$lambda, df1 = p$df1, df2 = p$df2, alpha = p$alpha,
    plot = FALSE, verbose = FALSE))
  .ig_ref_result(fit$power, "pwrss::power.f.test")
}

.ig_ref_gpower30 <- function(options, inputs) identical(options$effect_size_specification, "gpower30")

.ig_engine_info_f <- function() list(
  f_anova_oneway = list(source = "WebPower takes the total sample size (any N)",
    branches = list(.ig_ref_branch("WebPower", "wp.anova", "0.9.4", .ig_ref_f_anova_oneway, "WebPower::wp.anova",
      n_min = function(inputs) inputs$n_groups + 1,
      search = function(forward) .ig_g2_search_single(forward, n_min_fun = function(p) 2 * as.double(p$n_groups),
        step_fun = function(p) as.double(p$n_groups),
        grid_text = "N = k n with n >= 2 subjects in each of the k groups (equal group sizes)")))),
  f_anova_special_main_interactions = list(source = "WebPower; not pwr.f2.test, whose noncentrality differs",
    branches = list(.ig_ref_branch("WebPower", "wp.kanova", "0.9.4", .ig_ref_f_anova_special, "WebPower::wp.kanova",
      n_min = function(inputs) inputs$n_groups + 1,
      search = function(forward) .ig_g2_search_single(forward, n_min_fun = function(p) as.double(p$n_groups) + 1,
        step_fun = function(p) 1, grid_text = "every integer N > k (df2 = N - k >= 1)")))),
  f_ancova_main_interactions = list(source = "pwrss one-factor ANCOVA; other numerator df are not verified",
    branches = list(.ig_ref_branch("pwrss", "power.f.ancova", "1.3.1", .ig_ref_f_ancova, "pwrss::power.f.ancova",
      when = function(options, inputs) isTRUE(inputs$numerator_df == inputs$n_groups - 1),
      n_min = function(inputs) inputs$n_groups + inputs$n_covariates + 1,
      search = function(forward) .ig_g2_search_single(forward,
        n_min_fun = function(p) as.double(p$n_groups) + as.double(p$n_covariates) + 1, step_fun = function(p) 1,
        grid_text = "every integer N > k + c (df2 = N - k - c >= 1)")))),
  f_rm_anova_between = list(source = "pwrss mixed ANOVA, G*Power 3.0 effect-size specification only",
    branches = list(.ig_ref_branch("pwrss", "power.f.mixed.anova", "1.3.1", .ig_ref_f_rm_anova("between"),
      "pwrss::power.f.mixed.anova", when = .ig_ref_gpower30, n_min = function(inputs) inputs$k + 1))),
  f_rm_anova_within = list(source = "pwrss mixed ANOVA, G*Power 3.0 effect-size specification only",
    branches = list(.ig_ref_branch("pwrss", "power.f.mixed.anova", "1.3.1", .ig_ref_f_rm_anova("within"),
      "pwrss::power.f.mixed.anova", when = .ig_ref_gpower30, n_min = function(inputs) inputs$k + 1))),
  f_rm_anova_within_between = list(source = "pwrss mixed ANOVA, G*Power 3.0 effect-size specification only",
    branches = list(.ig_ref_branch("pwrss", "power.f.mixed.anova", "1.3.1", .ig_ref_f_rm_anova("interaction"),
      "pwrss::power.f.mixed.anova", when = .ig_ref_gpower30, n_min = function(inputs) inputs$k + 1))),
  f_mrc_r2_deviation_from_zero = list(source = "pwr noncentrality f^2 (u + v + 1) = f^2 N",
    branches = list(.ig_ref_branch("pwr", "pwr.f2.test", "1.3.0", .ig_ref_f_mrc_r2_deviation_from_zero,
      "pwr::pwr.f2.test", n_min = function(inputs) inputs$n_predictors + 2))),
  f_mrc_r2_increase = list(source = "WebPower; not pwr.f2.test, whose noncentrality differs",
    branches = list(.ig_ref_branch("WebPower", "wp.regression", "0.9.4", .ig_ref_f_mrc_r2_increase,
      "WebPower::wp.regression", n_min = function(inputs) inputs$n_predictors + 2))),
  f_generic = list(source = "pwrss generic noncentral F",
    branches = list(.ig_ref_branch("pwrss", "power.f.test", "1.3.1", .ig_ref_f_generic, "pwrss::power.f.test"))))
