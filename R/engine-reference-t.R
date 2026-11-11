# Reference forward models for the t-test procedures. Each returns the power
# computed by the original with every setting stated explicitly, and echoes the
# design quantities the solvers read (n_total, group_sizes). A one-tailed test
# is on the side of the effect's sign, as in G*Power and the base engines, so
# the originals receive |effect| and their upper-tailed alternative.

.ig_ref_t_one_two <- function(tails) if (identical(tails, "two")) "two.sided" else "one.sided"

.ig_ref_result <- function(power, call, n_total = NULL, group_sizes = NULL) {
  c(list(power = as.double(power)),
    if (!is.null(n_total)) list(n_total = as.double(n_total)),
    if (!is.null(group_sizes)) list(group_sizes = group_sizes),
    list(assumptions = sprintf("power computed by %s (reference engine)", call)))
}

.ig_ref_t_means_one_sample <- function(p) {
  fit <- stats::power.t.test(n = p$n_total, delta = abs(p$d), sd = 1, sig.level = p$alpha, type = "one.sample",
    alternative = .ig_ref_t_one_two(p$tails), strict = TRUE)
  .ig_ref_result(fit$power, "stats::power.t.test(type = \"one.sample\", strict = TRUE)", n_total = p$n_total)
}

.ig_ref_t_means_matched_pairs <- function(p) {
  fit <- stats::power.t.test(n = p$n_total, delta = abs(p$dz), sd = 1, sig.level = p$alpha, type = "paired",
    alternative = .ig_ref_t_one_two(p$tails), strict = TRUE)
  .ig_ref_result(fit$power, "stats::power.t.test(type = \"paired\", strict = TRUE)", n_total = p$n_total)
}

.ig_ref_t_means_two_independent <- function(p) {
  n1 <- as.double(p$n1)
  n2 <- as.double(p$n2)
  if (n1 < 2 || n2 < 2) .ig_ref_outside("pwr::pwr.t2n.test needs at least two observations per group")
  fit <- pwr::pwr.t2n.test(n1 = n1, n2 = n2, d = abs(p$d), sig.level = p$alpha,
    alternative = if (identical(p$tails, "two")) "two.sided" else "greater")
  .ig_ref_result(fit$power, "pwr::pwr.t2n.test", n_total = n1 + n2, group_sizes = c(n1 = n1, n2 = n2))
}

.ig_ref_t_correlation_point_biserial <- function(p) {
  if (p$n_total < 5) .ig_ref_outside("pwrss::power.t.regression needs at least 3 degrees of freedom (N >= 5)")
  fit <- pwrss::power.t.regression(beta = abs(p$rho), sd.predictor = 1, sd.outcome = 1, r.squared = p$rho^2,
    k.total = 1, n = p$n_total, alpha = p$alpha, alternative = .ig_ref_t_one_two(p$tails), verbose = FALSE)
  .ig_ref_result(fit$power, "pwrss::power.t.regression(beta = rho, r.squared = rho^2)", n_total = p$n_total)
}

.ig_ref_t_linreg_one_group_slope <- function(p) {
  if (p$n_total < 5) .ig_ref_outside("pwrss::power.t.regression needs at least 3 degrees of freedom (N >= 5)")
  k <- p$sd_x / p$sd_y
  # The H1 slope as the base engine derives it when only a distance is given.
  b1 <- if (!is.null(p$slope_h1)) p$slope_h1 else
    tanh(atanh(p$slope_h0 * k) + .ig_g5_effect_direction_sign(p, "slope_lt_h0") * p$z_distance) / k
  fit <- pwrss::power.t.regression(beta = b1, null.beta = p$slope_h0, sd.predictor = p$sd_x,
    sd.outcome = p$sd_y, r.squared = (b1 * k)^2, k.total = 1, n = p$n_total, alpha = p$alpha,
    alternative = .ig_ref_t_one_two(p$tails), verbose = FALSE)
  .ig_ref_result(fit$power, "pwrss::power.t.regression", n_total = p$n_total)
}

.ig_ref_t_mrc_single_coefficient <- function(p) {
  fit <- WebPower::wp.regression(n = p$n_total, p1 = p$n_predictors, p2 = p$n_predictors - 1, f2 = p$f2,
    alpha = p$alpha)
  .ig_ref_result(fit$power, "WebPower::wp.regression(p1 = p, p2 = p - 1)", n_total = p$n_total)
}

.ig_ref_t_generic <- function(p) {
  fit <- pwrss::power.t.test(ncp = abs(p$delta), df = p$df, alpha = p$alpha,
    alternative = .ig_ref_t_one_two(p$tails), plot = FALSE, verbose = FALSE)
  .ig_ref_result(fit$power, "pwrss::power.t.test")
}

.ig_ref_branch <- function(package, functions, min_version, power, call, generic_search = TRUE,
                           when = function(options, inputs) TRUE, n_min = NULL, search = NULL)
  list(when = when, package = package, functions = functions, min_version = min_version, power = power,
       generic_search = generic_search, call = call, n_min = n_min, search = search)

.ig_engine_info_t <- function() list(
  t_means_one_sample = list(source = "stats computes the same noncentral t power",
    branches = list(.ig_ref_branch("stats", "power.t.test", "4.1.0", .ig_ref_t_means_one_sample,
      "stats::power.t.test"))),
  t_means_matched_pairs = list(source = "stats computes the same noncentral t power",
    branches = list(.ig_ref_branch("stats", "power.t.test", "4.1.0", .ig_ref_t_means_matched_pairs,
      "stats::power.t.test"))),
  t_means_two_independent = list(source = paste("pwr covers unequal group sizes; an unbalanced a priori search",
      "brackets designs with fewer than two observations per group, which the original cannot evaluate"),
    branches = list(
      .ig_ref_branch("pwr", "pwr.t2n.test", "1.3.0", .ig_ref_t_means_two_independent, "pwr::pwr.t2n.test",
        when = function(options, inputs) is.null(inputs$allocation_ratio) || isTRUE(inputs$allocation_ratio == 1),
        n_min = 4),
      .ig_ref_branch("pwr", "pwr.t2n.test", "1.3.0", .ig_ref_t_means_two_independent, "pwr::pwr.t2n.test",
        generic_search = FALSE))),
  t_correlation_point_biserial = list(source = "pwrss regression t test (not pwr.r.test, which uses Fisher z)",
    branches = list(.ig_ref_branch("pwrss", "power.t.regression", "1.3.1", .ig_ref_t_correlation_point_biserial,
      "pwrss::power.t.regression", n_min = 5))),
  t_linreg_one_group_slope = list(source = "pwrss regression t test with the residual SD",
    branches = list(.ig_ref_branch("pwrss", "power.t.regression", "1.3.1", .ig_ref_t_linreg_one_group_slope,
      "pwrss::power.t.regression", n_min = 5))),
  t_mrc_single_coefficient = list(source = "WebPower F test of one coefficient; two-tailed only",
    branches = list(.ig_ref_branch("WebPower", "wp.regression", "0.9.4", .ig_ref_t_mrc_single_coefficient,
      "WebPower::wp.regression", when = function(options, inputs) identical(options$tails, "two")))),
  t_generic = list(source = "pwrss generic noncentral t",
    branches = list(.ig_ref_branch("pwrss", "power.t.test", "1.3.1", .ig_ref_t_generic, "pwrss::power.t.test"))))
