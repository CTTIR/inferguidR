# Reference forward models for the Wilcoxon procedures: pwrss's A.R.E. method
# (method = "guenther") with a predefined normal, Laplace or logistic parent.
# The Lehmann method, manual or minimum ARE, effect size p1 and the continuity
# correction have no matching original.

.ig_ref_wilcoxon <- function(design) function(p) {
  two_groups <- identical(design, "independent")
  n1 <- as.double(if (two_groups) p$n1 else p$n_total)
  n2 <- as.double(if (two_groups) p$n2 else p$n_total)
  if (min(n1, n2) < (if (two_groups) 3 else 5))
    .ig_ref_outside("pwrss::power.np.wilcoxon needs n >= 5 (one sample, pairs) or >= 3 per group")
  effect <- abs(if (identical(design, "paired")) p$dz else p$d)
  # pwrss rejects designs whose adjusted degrees of freedom fall below 3.
  fit <- tryCatch(.ig_ref_quiet(pwrss::power.np.wilcoxon(d = effect, n2 = n2, n.ratio = n1 / n2, alpha = p$alpha,
    alternative = .ig_ref_t_one_two(p$tails), design = design, distribution = p$parent_distribution,
    method = "guenther", ceil.n = FALSE, verbose = FALSE)),
    error = function(e) .ig_ref_outside(paste("pwrss::power.np.wilcoxon:", conditionMessage(e))))
  .ig_ref_result(fit$power, sprintf("pwrss::power.np.wilcoxon(design = \"%s\", method = \"guenther\")", design),
    n_total = if (two_groups) n1 + n2 else n1, group_sizes = if (two_groups) c(n1 = n1, n2 = n2))
}

.ig_ref_wilcoxon_when <- function(options, inputs) identical(options$method, "are") &&
  identical(options$are_specification, "predefined_parent") &&
  options$parent_distribution %in% c("normal", "laplace", "logistic") &&
  (is.null(options$effect_size_parameter) || identical(options$effect_size_parameter, "d")) &&
  (is.null(options$continuity_correction) || identical(options$continuity_correction, "off"))

.ig_engine_info_wilcoxon <- function() list(
  t_wilcoxon_signed_rank_one_sample = list(source = "pwrss A.R.E. (Guenther) with a predefined parent",
    branches = list(.ig_ref_branch("pwrss", "power.np.wilcoxon", "1.3.1", .ig_ref_wilcoxon("one.sample"),
      "pwrss::power.np.wilcoxon", when = .ig_ref_wilcoxon_when, n_min = 6))),
  t_wilcoxon_signed_rank_matched_pairs = list(source = "pwrss A.R.E. (Guenther) with a predefined parent",
    branches = list(.ig_ref_branch("pwrss", "power.np.wilcoxon", "1.3.1", .ig_ref_wilcoxon("paired"),
      "pwrss::power.np.wilcoxon", when = .ig_ref_wilcoxon_when, n_min = 6))),
  t_wilcoxon_mann_whitney = list(source = "pwrss A.R.E. (Guenther) with a predefined parent; balanced a priori only",
    branches = list(.ig_ref_branch("pwrss", "power.np.wilcoxon", "1.3.1", .ig_ref_wilcoxon("independent"),
      "pwrss::power.np.wilcoxon", when = function(options, inputs) .ig_ref_wilcoxon_when(options, inputs) &&
        (is.null(inputs$allocation_ratio) || isTRUE(inputs$allocation_ratio == 1)), n_min = 8))))
