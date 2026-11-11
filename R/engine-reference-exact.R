# Reference forward models for the exact and chi-squared procedures. The exact
# tests are discrete: their power is not monotone in N and the base engines
# search a priori designs with dedicated rules, so their reference engines
# answer post hoc analyses only (generic_search = FALSE).

.ig_ref_binomial <- function(p0_of, prob_of) function(p) {
  fit <- .ig_ref_quiet(pwrss::power.exact.oneprop(prob = prob_of(p), null.prob = p0_of(p), n = p$n_total,
    alpha = p$alpha, alternative = .ig_ref_t_one_two(p$tails), verbose = FALSE))
  .ig_ref_result(fit$power, "pwrss::power.exact.oneprop", n_total = p$n_total, discrete = TRUE)
}

.ig_ref_mcnemar <- function(p) {
  # MESS takes the ratio psi >= 1 and the probability of the less frequent
  # discordant cell, pd / (1 + psi); the two-tailed test is symmetric in OR and 1 / OR.
  psi <- max(p$or, 1 / p$or)
  fit <- MESS::power_mcnemar_test(n = p$n_total, paid = p$pd / (1 + psi), psi = psi, sig.level = p$alpha,
    alternative = .ig_ref_t_one_two(p$tails), method = "cond.exact")
  .ig_ref_result(fit$power, "MESS::power_mcnemar_test(method = \"cond.exact\")", n_total = p$n_total, discrete = TRUE)
}

.ig_ref_fisher <- function(p) {
  fit <- Exact::power.exact.test(p1 = p$p1, p2 = p$p2, n1 = p$n1, n2 = p$n2, alpha = p$alpha, method = "fisher",
    alternative = if (identical(p$tails, "two")) "two.sided" else if (p$p1 < p$p2) "less" else "greater")
  .ig_ref_result(fit$power, "Exact::power.exact.test(method = \"fisher\")", n_total = p$n1 + p$n2,
    group_sizes = c(n1 = p$n1, n2 = p$n2), discrete = TRUE)
}

.ig_ref_correlation <- function(p) {
  exact <- identical(p$computation, "exact")
  fit <- .ig_ref_quiet((if (exact) pwrss::power.exact.onecor else pwrss::power.z.onecor)(rho = p$rho_h1,
    null.rho = p$rho_h0, n = p$n_total, alpha = p$alpha, alternative = .ig_ref_t_one_two(p$tails), verbose = FALSE))
  .ig_ref_result(fit$power, if (exact) "pwrss::power.exact.onecor" else "pwrss::power.z.onecor", n_total = p$n_total)
}

.ig_ref_chi2_gof <- function(p) {
  if (p$n_total <= p$df) .ig_ref_outside("pwr::pwr.chisq.test needs N > df")
  fit <- pwr::pwr.chisq.test(w = p$w, N = p$n_total, df = p$df, sig.level = p$alpha)
  .ig_ref_result(fit$power, "pwr::pwr.chisq.test", n_total = p$n_total)
}

.ig_ref_chi2_generic <- function(p) {
  fit <- .ig_ref_quiet(pwrss::power.chisq.test(ncp = p$lambda, df = p$df, alpha = p$alpha, plot = FALSE,
    verbose = FALSE))
  .ig_ref_result(fit$power, "pwrss::power.chisq.test")
}

.ig_ref_boschloo <- function(p) {
  args <- list(p1 = p$p1, p2 = p$p2, n1 = p$n1, n2 = p$n2, alpha = p$alpha, method = "boschloo")
  args <- c(args, if (identical(p$variant, "minlike")) list(alternative = "two.sided", tsmethod = "square") else
    list(alternative = p$variant))
  fit <- do.call(Exact::power.exact.test, args)
  .ig_ref_result(fit$power, "Exact::power.exact.test(method = \"boschloo\")", n_total = p$n1 + p$n2,
    group_sizes = c(n1 = p$n1, n2 = p$n2), discrete = TRUE)
}

.ig_ref_binomial_branches <- function(p0_of, prob_of) list(
  .ig_ref_branch("pwrss", "power.exact.oneprop", "1.3.1", .ig_ref_binomial(p0_of, prob_of),
    "pwrss::power.exact.oneprop", generic_search = FALSE,
    when = function(options, inputs) !identical(options$tails, "two") || identical(options$alpha_balancing, "half_each_side")))

.ig_engine_info_exact <- function() list(
  exact_proportion_one_sample_binomial = list(source = "one-tailed and half-each-side two-sided: pwrss; the other balancing rules have no matching original",
    branches = .ig_ref_binomial_branches(function(p) p$p0, function(p) p$p0 + p$g)),
  exact_proportion_sign_test = list(source = "one-tailed and half-each-side two-sided: pwrss; the other balancing rules have no matching original",
    branches = .ig_ref_binomial_branches(function(p) 0.5, function(p) 0.5 + p$g)),
  exact_generic_binomial = list(source = "one-tailed and half-each-side two-sided: pwrss; the other balancing rules have no matching original",
    branches = .ig_ref_binomial_branches(function(p) p$p1, function(p) p$p2)),
  exact_proportions_mcnemar = list(source = paste("fast approximation, two-tailed: MESS conditional exact; one-tailed",
      "(MESS is 2.8e-5 off) and the exact unconditional computation (exact2x2, about 1e-7 off) have no matching original"),
    branches = list(.ig_ref_branch("MESS", "power_mcnemar_test", "0.6.0", .ig_ref_mcnemar, "MESS::power_mcnemar_test",
      generic_search = FALSE, when = function(options, inputs) identical(options$computation, "fast_approx") &&
        identical(options$tails, "two") && identical(options$alpha_balancing, "half_each_side")))),
  exact_proportions_fisher = list(source = "Exact; the symmetric Fisher two-sided rule and one-tailed tests only",
    branches = list(.ig_ref_branch("Exact", "power.exact.test", "3.3", .ig_ref_fisher, "Exact::power.exact.test",
      generic_search = FALSE, when = function(options, inputs) identical(options$computation, "exact") &&
        !is.null(inputs$n1) && (!identical(options$tails, "two") ||
          identical(options$two_sided_balancing, "symmetric_fisher_statistic"))))),
  exact_correlation_bivariate_normal = list(source = "pwrss exact distribution of r, or Fisher z",
    branches = list(.ig_ref_branch("pwrss", "power.exact.onecor", "1.3.1", .ig_ref_correlation,
      "pwrss::power.exact.onecor", generic_search = FALSE,
      when = function(options, inputs) !is.null(inputs$rho_h1)))),
  chi2_goodness_of_fit_contingency = list(source = "pwr noncentral chi-squared",
    branches = list(.ig_ref_branch("pwr", "pwr.chisq.test", "1.3.0", .ig_ref_chi2_gof, "pwr::pwr.chisq.test",
      n_min = function(inputs) inputs$df + 1))),
  chi2_generic = list(source = "pwrss generic noncentral chi-squared",
    branches = list(.ig_ref_branch("pwrss", "power.chisq.test", "1.3.1", .ig_ref_chi2_generic,
      "pwrss::power.chisq.test"))),
  boschloo_proportions_two_independent = list(source = "Exact; Exact's central variant differs at equal proportions",
    branches = list(.ig_ref_branch("Exact", "power.exact.test", "3.3", .ig_ref_boschloo, "Exact::power.exact.test",
      generic_search = FALSE, when = function(options, inputs) options$variant %in% c("minlike", "less", "greater") &&
        !is.null(inputs$n1)))))
