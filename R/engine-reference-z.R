# Reference forward models for the z-test procedures. Branches cover the
# option combinations an original computes with the base definition; other
# combinations are unsupported. One-tailed tests are on the side of the effect.

.ig_ref_side <- function(tails, effect) if (identical(tails, "two")) "two.sided" else if (effect >= 0) "greater" else "less"

.ig_ref_z_proportions <- function(original) function(p) {
  n1 <- as.double(p$n1)
  n2 <- as.double(p$n2)
  p1 <- p$p1
  p2 <- if (!is.null(p$p2)) p$p2 else sin(asin(sqrt(p1)) + p$h / 2)^2
  if (n1 < 2 || n2 < 2) .ig_ref_outside("the originals need at least two observations per group")
  power <- switch(original,
    stats = stats::power.prop.test(n = n1, p1 = p1, p2 = p2, sig.level = p$alpha,
      alternative = .ig_ref_t_one_two(p$tails), strict = TRUE)$power,
    MESS = MESS::power_prop_test(n = n1, p1 = p1, p2 = p2, ratio = n2 / n1, sig.level = p$alpha,
      alternative = "two.sided")$power,
    pwr = pwr::pwr.2p2n.test(h = abs(2 * asin(sqrt(p2)) - 2 * asin(sqrt(p1))), n1 = n1, n2 = n2, sig.level = p$alpha,
      alternative = if (identical(p$tails, "two")) "two.sided" else "greater")$power)
  .ig_ref_result(power, c(stats = "stats::power.prop.test(strict = TRUE)", MESS = "MESS::power_prop_test",
    pwr = "pwr::pwr.2p2n.test")[[original]], n_total = n1 + n2, group_sizes = c(n1 = n1, n2 = n2))
}

.ig_ref_z_correlations_two_independent <- function(p) {
  n1 <- as.double(p$n1)
  n2 <- as.double(p$n2)
  if (n1 < 4 || n2 < 4) .ig_ref_outside("pwrss::power.z.twocors needs at least four observations per group")
  # The test depends on q = atanh(r1) - atanh(r2) only: r1 = tanh(|q|), r2 = 0.
  fit <- .ig_ref_quiet(pwrss::power.z.twocors(rho1 = tanh(abs(p$q)), rho2 = 0, n2 = n2, n.ratio = n1 / n2,
    alpha = p$alpha, alternative = .ig_ref_t_one_two(p$tails), ceil.n = FALSE, verbose = FALSE))
  .ig_ref_result(fit$power, "pwrss::power.z.twocors", n_total = n1 + n2, group_sizes = c(n1 = n1, n2 = n2))
}

.ig_ref_z_steiger <- function(common) function(p) {
  if (p$n_total < 4) .ig_ref_outside("pwrss::power.z.twocors.steiger needs at least four observations")
  args <- if (common) list(rho12 = p$rho_ab, rho13 = p$rho_ac, rho23 = p$rho_bc) else
    list(rho12 = p$rho_ab, rho34 = p$rho_cd, rho13 = p$rho_ac, rho14 = p$rho_ad, rho23 = p$rho_bc, rho24 = p$rho_bd)
  fit <- .ig_ref_quiet(do.call(pwrss::power.z.twocors.steiger, c(args, list(n = p$n_total, alpha = p$alpha,
    alternative = .ig_ref_t_one_two(p$tails), pooled = FALSE, common.index = common, ceil.n = FALSE,
    verbose = FALSE))))
  .ig_ref_result(fit$power, "pwrss::power.z.twocors.steiger(pooled = FALSE)", n_total = p$n_total)
}

.ig_ref_logistic_p2 <- function(p) if (!is.null(p$p2)) p$p2 else p$p1 * p$odds_ratio / (1 - p$p1 + p$p1 * p$odds_ratio)

.ig_ref_z_logistic <- function(original) function(p) {
  if (p$n_total < 10) .ig_ref_outside("the logistic originals need at least ten observations")
  p2 <- .ig_ref_logistic_p2(p)
  odds <- p2 / (1 - p2) / (p$p1 / (1 - p$p1))
  binary <- identical(p$x_distribution, "binomial")
  power <- switch(original,
    WebPower = .ig_ref_quiet(WebPower::wp.logistic(n = p$n_total, p0 = p$p1, p1 = p2, alpha = p$alpha,
      alternative = .ig_ref_side(p$tails, log(odds)), family = if (binary) "Bernoulli" else "normal",
      parameter = if (binary) p$x_pi else c(p$x_mu, p$x_sigma)))$power,
    powerMediation = powerMediation::powerLogisticCon(n = p$n_total, p1 = p$p1, OR = odds, alpha = p$alpha))
  .ig_ref_result(power, c(WebPower = "WebPower::wp.logistic", powerMediation = "powerMediation::powerLogisticCon")[[original]], n_total = p$n_total)
}

.ig_ref_z_poisson <- function(original) function(p) {
  if (p$n_total < 10) .ig_ref_outside("the Poisson originals need at least ten observations")
  binary <- identical(p$x_distribution, "binomial")
  power <- switch(original,
    WebPower = .ig_ref_quiet(WebPower::wp.poisson(n = p$n_total, exp0 = p$base_rate * p$mean_exposure,
      exp1 = p$exp_beta1, alpha = p$alpha, alternative = .ig_ref_side(p$tails, log(p$exp_beta1)),
      family = if (binary) "Bernoulli" else "normal", parameter = if (binary) p$x_pi else c(p$x_mu, p$x_sigma)))$power)
  .ig_ref_result(power, "WebPower::wp.poisson",
    n_total = p$n_total)
}

.ig_ref_z_generic <- function(p) {
  fit <- .ig_ref_quiet(pwrss::power.z.test(mean = abs(p$mu), sd = p$sigma, alpha = p$alpha,
    alternative = .ig_ref_t_one_two(p$tails), plot = FALSE, verbose = FALSE))
  .ig_ref_result(fit$power, "pwrss::power.z.test")
}

# Option tests for the GLM branches.
.ig_ref_glm_large_wald <- function(options) identical(options$computation, "large_sample") &&
  (is.null(options$test_statistic) || identical(options$test_statistic, "wald"))
.ig_ref_std_normal <- function(options, inputs) identical(options$x_distribution, "normal") &&
  isTRUE(inputs$x_mu == 0) && isTRUE(inputs$x_sigma == 1)
.ig_ref_normal_or_binary <- function(options) options$x_distribution %in% c("normal", "binomial")
.ig_ref_no_other_x <- function(inputs) isTRUE(inputs$r2_other_x == 0)

.ig_engine_info_z <- function() list(
  z_proportions_two_independent = list(source = paste("equal groups: stats; unequal groups two-sided: MESS; arcsine:",
      "pwr. Unequal one-sided and continuity correction have no matching original"),
    branches = list(
      .ig_ref_branch("pwr", "pwr.2p2n.test", "1.3.0", .ig_ref_z_proportions("pwr"), "pwr::pwr.2p2n.test",
        when = function(options, inputs) identical(options$arcsin_transform, "on") &&
          identical(options$continuity_correction, "off"), n_min = 4),
      .ig_ref_branch("stats", "power.prop.test", "4.1.0", .ig_ref_z_proportions("stats"), "stats::power.prop.test",
        when = function(options, inputs) identical(options$continuity_correction, "off") &&
          (isTRUE(inputs$n1 == inputs$n2) || isTRUE(inputs$allocation_ratio == 1)), n_min = 4),
      .ig_ref_branch("MESS", "power_prop_test", "0.6.0", .ig_ref_z_proportions("MESS"), "MESS::power_prop_test",
        when = function(options, inputs) identical(options$continuity_correction, "off") &&
          identical(options$tails, "two") && !is.null(inputs$n1), generic_search = FALSE))),
  z_correlations_two_independent = list(source = "pwrss Fisher-z test of two independent correlations",
    branches = list(.ig_ref_branch("pwrss", "power.z.twocors", "1.3.1", .ig_ref_z_correlations_two_independent,
      "pwrss::power.z.twocors", when = function(options, inputs) is.null(inputs$allocation_ratio) ||
        isTRUE(inputs$allocation_ratio == 1), n_min = 8))),
  z_correlations_two_dependent_common_index = list(source = "pwrss Steiger test with pooled = FALSE",
    branches = list(.ig_ref_branch("pwrss", "power.z.twocors.steiger", "1.3.1", .ig_ref_z_steiger(TRUE),
      "pwrss::power.z.twocors.steiger", n_min = 4))),
  z_correlations_two_dependent_no_common_index = list(source = "pwrss Steiger test with pooled = FALSE",
    branches = list(.ig_ref_branch("pwrss", "power.z.twocors.steiger", "1.3.1", .ig_ref_z_steiger(FALSE),
      "pwrss::power.z.twocors.steiger", n_min = 4))),
  z_logistic_regression = list(source = paste("large-sample Wald test: Demidenko without variance correction",
      "(WebPower), Hsieh (powerMediation). Demidenko with variance correction (pwrss, 2.5e-7 off), enumeration and",
      "likelihood ratio have no matching original"),
    branches = list(
      .ig_ref_branch("WebPower", "wp.logistic", "0.9.4", .ig_ref_z_logistic("WebPower"), "WebPower::wp.logistic",
        when = function(options, inputs) .ig_ref_glm_large_wald(options) && identical(options$test_procedure, "demidenko") &&
          identical(options$variance_correction, "off") && .ig_ref_normal_or_binary(options) && .ig_ref_no_other_x(inputs),
        n_min = 10),
      .ig_ref_branch("powerMediation", "powerLogisticCon", "0.3.4", .ig_ref_z_logistic("powerMediation"),
        "powerMediation::powerLogisticCon", when = function(options, inputs) .ig_ref_glm_large_wald(options) &&
          identical(options$test_procedure, "hsieh") && identical(options$tails, "two") &&
          .ig_ref_std_normal(options, inputs) && .ig_ref_no_other_x(inputs), n_min = 10))),
  z_poisson_regression = list(source = paste("large-sample Wald test: Demidenko (WebPower). Signorini",
      "(powerMediation, 6.2e-6 off) and other combinations have no matching original"),
    branches = list(
      .ig_ref_branch("WebPower", "wp.poisson", "0.9.4", .ig_ref_z_poisson("WebPower"), "WebPower::wp.poisson",
        when = function(options, inputs) .ig_ref_glm_large_wald(options) && identical(options$test_procedure, "demidenko") &&
          (identical(options$variance_correction, "off") || identical(options$x_distribution, "normal")) &&
          .ig_ref_normal_or_binary(options) && .ig_ref_no_other_x(inputs), n_min = 10))),
  z_generic = list(source = "pwrss generic z",
    branches = list(.ig_ref_branch("pwrss", "power.z.test", "1.3.1", .ig_ref_z_generic, "pwrss::power.z.test"))))
