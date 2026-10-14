# Method definitions: t tests for means (group G1) -------------------------------------
#
# Ownership: R/methods-t.R belongs to group G1 (phase 06). It holds the three
# t tests for means of the G*Power inventory:
#
#   t_means_one_sample        Means: Difference from constant (one sample case)
#   t_means_matched_pairs     Means: Difference between two dependent means (matched pairs)
#   t_means_two_independent   Means: Difference between two independent means (two groups)
#                             (the phase 05 vertical slice, extended here by the a priori
#                             rule for unequal allocation)
#
# and the shared noncentral-t power helper .ig_t_power_g1(), which
# R/methods-generic.R reuses for t_generic. Method notes:
# specs/method-notes/<method_id>.md.

.ig_define_methods_t_g1 <- function() {
  list(.ig_method_t_means_one_sample(), .ig_method_t_means_matched_pairs(),
       .ig_method_t_means_two_independent())
}

# Shared pieces --------------------------------------------------------------------------

# Power of a t test with H0 distribution t(df) and H1 distribution t(df, ncp)
# in G*Power's convention (manual 25.5): two tails put alpha/2 in each tail
# and count both rejection tails; one tail rejects on the side of the sign of
# ncp (the upper tail for ncp = 0).
#
# The probabilities are always evaluated at |ncp|. For a negative ncp the
# test is the mirror image of the test at -ncp (T -> -T), so power, beta and
# the absolute critical values are identical and only the side changes. This
# is an exact identity, used because it keeps every evaluation on the upper
# side of a positive noncentrality: nmath's pt() returns lower tails of a
# negative ncp as 1 - (upper tail) and can exceed 1 by a few 1e-11 (for
# example ncp = -22.36, df = 199998, t = -3.09: 1 + 3.8e-11), which the
# solvers reject as a power outside [0, 1].
.ig_t_power_g1 <- function(df, ncp, alpha, tails) {
  core <- .ig_power_t_test(df, abs(ncp), alpha, tails)
  if (ncp < 0 && !identical(tails, "two")) {
    core$critical <- -core$critical
    core$side <- "lower"
  }
  # Power and beta are complementary. The smaller one is the accurately
  # computed probability (small tails are integrated, R/num-dist.R); the larger
  # one is taken as its complement, so that a power near one never exceeds one
  # through the absolute error of nmath's series (for example ncp = 22.36,
  # df = 199998, t = 3.09: nmath's lower tail is -3.8e-11, its upper tail
  # 1 + 3.8e-11).
  if (is.finite(core$power) && is.finite(core$beta)) {
    if (core$beta < core$power) core$power <- 1 - core$beta else core$beta <- 1 - core$power
  }
  core
}

.ig_t_tail_note_g1 <- function(tails, effect = "d") {
  if (identical(tails, "two")) {
    "two-sided test: alpha/2 in each tail; the power counts both rejection tails (G*Power convention, manual 25.5)"
  } else {
    sprintf(paste("one-sided test: the rejection region lies on the side of the sign of %s (G*Power manual 25.5);",
                  "%s = 0 uses the upper tail, where the power equals alpha"), effect, effect)
  }
}

# A target power below alpha cannot be a planning target of these tests: the
# power is at least alpha for every effect size (it equals alpha at no
# effect), so every design "reaches" it. G*Power refuses it with the message
# 'Power must be greater than or equal to alpha.' (a priori and sensitivity).
.ig_t_power_alpha_constraint_g1 <- function(p) {
  if (!is.null(p$power) && !is.null(p$alpha) && is.numeric(p$power) && is.numeric(p$alpha) &&
      p$power < p$alpha) {
    return(sprintf(paste("the target power %s is below alpha %s; the power of this test is at least alpha for",
                         "every effect size (G*Power: 'Power must be greater than or equal to alpha.')"),
                   .ig_format_double(p$power), .ig_format_double(p$alpha)))
  }
  character()
}

# Forward model of the one-sample t test, shared by t_means_one_sample and
# t_means_matched_pairs (manual 19.5 and 20.5): df = N - 1,
# delta = effect * sqrt(N).
.ig_power_t_single_g1 <- function(n, effect, alpha, tails, effect_name, assumptions) {
  n <- as.double(n)
  df <- n - 1
  ncp <- effect * sqrt(n)
  core <- .ig_t_power_g1(df, ncp, alpha, tails)
  list(
    power = core$power,
    beta = core$beta,
    alpha_attained = alpha,
    critical = core$critical,
    df = df,
    ncp = ncp,
    n_total = n,
    statistic_distribution = list(
      h0 = list(family = "t", df = df, ncp = 0),
      h1 = list(family = "t", df = df, ncp = ncp),
      rejection_side = core$side
    ),
    assumptions = c(assumptions, .ig_t_tail_note_g1(tails, effect_name))
  )
}

.ig_t_tails_option_g1 <- function() {
  ig_option("tails", values = c("one", "two"), default = "one", labels = c(one = "One", two = "Two"),
            label = list(en = "Tail(s)", de = "Seiten"))
}

.ig_t_engine_g1 <- function() {
  list(package = "stats",
       fun = "pt/qt (nmath AS 243); inferguidR integral for |ncp| > 37.62, df > 2e4 and tails < 1e-5",
       version = NA_character_)
}

# Means: Difference from constant (one sample case) ------------------------------------------
#
# G*Power: Tests > Means > One group: Difference from constant; label
# "t tests - Means: Difference from constant (one sample case)" (manual
# ch. 20). H0: mu = mu0; d = (mu - mu0) / sigma;
#   H0 distribution t(N - 1), H1 distribution t(N - 1, d sqrt(N)).
# Single-sample design: N is the design. G*Power's a priori search starts at
# N = 2 (df = 1; supplement capture SUP-G1-T-MEAN1-APR-A01), while its input
# fields refuse N < 3 ('Total sample size is outside the admissible range (3,
# 2147483646)'); inferguidR admits N >= 2 everywhere (method notes).

.ig_method_t_means_one_sample <- function() {
  ig_method(
    id = "t_means_one_sample",
    family = "t",
    label = list(en = "Means: Difference from constant (one sample case)",
                 de = "Mittelwerte: Differenz zu einer Konstanten (eine Stichprobe)"),
    params = ig_params(
      ig_param("d", role = "effect", domain = ig_domain_real(),
               label = list(en = "Effect size d", de = "Effektst\u00e4rke d"),
               description = "Standardised difference d = (mu - mu0) / sigma (manual 20.1).",
               aliases = "effect_size_d"),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      ig_param("n_total", role = "design", domain = ig_domain_integer(lower = 2), unit = "subjects",
               label = list(en = "Total sample size", de = "Gesamtstichprobenumfang"),
               description = paste("Sample size N (df = N - 1). Solved in a priori mode; an input in the other",
                                   "modes. G*Power's input fields require N >= 3, its a priori search starts at",
                                   "N = 2; inferguidR admits N >= 2 (df >= 1) in every mode."),
               aliases = "total_sample_size")
    ),
    options = ig_options(.ig_t_tails_option_g1()),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = .ig_power_t_means_one_sample,
    design = NULL,
    discrete = FALSE,
    effect_direction = "increasing",
    engine = .ig_t_engine_g1(),
    checker = "validation/independent/g1/check_g1.R",
    reference = list(
      gpower_path = "Tests > Means > One group: Difference from constant",
      gpower_label = "t tests - Means: Difference from constant (one sample case)",
      manual_pages = "ch. 20 (PDF pp. 50-51); tail convention 25.5 (PDF p. 63)",
      allocation_rule = "single sample: N is the design (a priori search from N = 2, as G*Power)",
      method_notes = "specs/method-notes/t_means_one_sample.md"
    ),
    status = "independently_checked",
    constraints = .ig_t_power_alpha_constraint_g1
  )
}

.ig_power_t_means_one_sample <- function(p) {
  .ig_power_t_single_g1(p$n_total, p$d, p$alpha, p$tails, "d", c(
    "one sample from a normal population with unknown variance (one-sample t test)",
    "d = (mu - mu0) / sigma; df = N - 1; ncp = d * sqrt(N) (manual 20.5)"))
}

# Means: Difference between two dependent means (matched pairs) ------------------------------
#
# G*Power: Tests > Means > Two dependent groups (matched pairs); label "t
# tests - Means: Difference between two dependent means (matched pairs)"
# (manual ch. 19). The test is the one-sample t test of the pair differences
# z = x - y: dz = mu_z / sigma_z with sigma_z = sqrt(sx^2 + sy^2 - 2 rho sx sy);
# H1 distribution t(N - 1, dz sqrt(N)). "Total sample size" is the number of
# pairs (manual 19.3: "110 pairs, that is, 220 subjects").

.ig_method_t_means_matched_pairs <- function() {
  ig_method(
    id = "t_means_matched_pairs",
    family = "t",
    label = list(en = "Means: Difference between two dependent means (matched pairs)",
                 de = "Mittelwerte: Differenz zweier abh\u00e4ngiger Mittelwerte (verbundene Paare)"),
    params = ig_params(
      ig_param("dz", role = "effect", domain = ig_domain_real(),
               label = list(en = "Effect size dz", de = "Effektst\u00e4rke dz"),
               description = paste("Standardised mean difference of the pairs dz = mu_z / sigma_z with",
                                   "sigma_z = sqrt(sx^2 + sy^2 - 2 rho sx sy) (manual 19.1)."),
               aliases = "effect_size_dz"),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      ig_param("n_total", role = "design", domain = ig_domain_integer(lower = 2), unit = "pairs",
               label = list(en = "Total sample size (number of pairs)", de = "Anzahl der Paare"),
               description = paste("Number of pairs N (df = N - 1); 2N subjects are measured. Solved in a priori",
                                   "mode; G*Power's input fields require N >= 3, its a priori search starts at",
                                   "N = 2; inferguidR admits N >= 2 in every mode."),
               aliases = c("total_sample_size", "n_pairs"))
    ),
    options = ig_options(.ig_t_tails_option_g1()),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = .ig_power_t_means_matched_pairs,
    design = NULL,
    discrete = FALSE,
    effect_direction = "increasing",
    engine = .ig_t_engine_g1(),
    checker = "validation/independent/g1/check_g1.R",
    reference = list(
      gpower_path = "Tests > Means > Two dependent groups (matched pairs)",
      gpower_label = "t tests - Means: Difference between two dependent means (matched pairs)",
      manual_pages = "ch. 19 (PDF pp. 48-49); tail convention 25.5 (PDF p. 63)",
      allocation_rule = "N pairs: N is the design (a priori search from N = 2 pairs, as G*Power)",
      method_notes = "specs/method-notes/t_means_matched_pairs.md"
    ),
    status = "independently_checked",
    constraints = .ig_t_power_alpha_constraint_g1
  )
}

.ig_power_t_means_matched_pairs <- function(p) {
  .ig_power_t_single_g1(p$n_total, p$dz, p$alpha, p$tails, "dz", c(
    "N pairs whose differences z = x - y are normal with unknown variance (one-sample t test of the differences)",
    "dz = mu_z / sigma_z; df = N - 1; ncp = dz * sqrt(N), N = number of pairs (manual 19.5)"))
}

# Means: Difference between two independent means (two groups) ------------------------------
#
# G*Power: Tests > Means > Two independent groups; distribution-based label
# "t tests - Means: Difference between two independent means (two groups)"
# (manual ch. 21; examples in 1.2.3). Pooled-variance Student t test:
#
#   H0 distribution t(N - 2), H1 distribution t(N - 2, delta),
#   N = n1 + n2,  delta = d * sqrt(n1 n2 / (n1 + n2))     (manual 21.5)
#
# Tail convention (manual 25.5, specs/engine-decisions.md 5.1): two-sided
# tests put alpha/2 in each tail and count both rejection tails; a one-sided
# critical value has the sign of delta.
#
# A priori allocation (allocation ratio r = N2/N1; validation/compare/
# t_means_two_independent/allocation-rule.md, 17 a priori captures):
#
#   1. N* = the smallest even N >= 4 whose power with fractional group sizes
#      N/(1 + r) and rN/(1 + r) reaches the target (G*Power 3.1.9.7's rule);
#   2. n1 = floor(N*/(1 + r) + 1/2), n2 = N* - n1, and every reported quantity
#      is evaluated at this integer design;
#   3. where G*Power's design misses the target power at the integer design
#      (captures APR-B001-E01, SUP-RSC A04) or leaves a group empty
#      (SUP-RSC A09), inferguidR reports the smallest even N > N* whose
#      rounded split has both groups non-empty and reaches the target, and
#      records G*Power's design in diagnostics$gpower_allocation_rule. This
#      is a deliberate, documented difference (an a priori answer of
#      inferguidR always reaches the target power).
#
# With r = 1 the fractional and integer criteria coincide (n1 = n2 = N/2) and
# the generic integer search is used.

.ig_method_t_means_two_independent <- function() {
  ig_method(
    id = "t_means_two_independent",
    family = "t",
    label = list(en = "Means: Difference between two independent means (two groups)",
                 de = "Mittelwerte: Differenz zweier unabh\u00e4ngiger Mittelwerte (zwei Gruppen)"),
    params = ig_params(
      ig_param("d", role = "effect", domain = ig_domain_real(),
               label = list(en = "Effect size d", de = "Effektst\u00e4rke d"),
               description = paste("Standardised mean difference d = (mu1 - mu2) / sigma with the common",
                                   "within-group standard deviation sigma (manual 21.1)."),
               aliases = "effect_size_d"),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      ig_param("n_total", role = "design", domain = ig_domain_integer(lower = 4), unit = "subjects",
               label = list(en = "Total sample size", de = "Gesamtstichprobenumfang"),
               description = "Solved in a priori mode (n1 + n2).",
               modes = "a_priori"),
      ig_param("n1", role = "design", domain = ig_domain_integer(lower = 1), unit = "subjects",
               label = list(en = "Sample size group 1", de = "Stichprobenumfang Gruppe 1"),
               aliases = "sample_size_group_1",
               modes = c("post_hoc", "sensitivity", "criterion", "compromise")),
      ig_param("n2", role = "design", domain = ig_domain_integer(lower = 1), unit = "subjects",
               label = list(en = "Sample size group 2", de = "Stichprobenumfang Gruppe 2"),
               aliases = "sample_size_group_2",
               modes = c("post_hoc", "sensitivity", "criterion", "compromise")),
      ig_param("allocation_ratio", role = "design",
               domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE),
               default = 1,
               label = list(en = "Allocation ratio N2/N1", de = "Allokationsverh\u00e4ltnis N2/N1"),
               description = paste("Group-size ratio n2/n1 of the a priori design (G*Power's initial field",
                                   "value 1, live capture). Unequal ratios follow G*Power's rule: the smallest",
                                   "even N whose fractional-allocation power reaches the target, split as",
                                   "n1 = floor(N/(1 + r) + 1/2), n2 = N - n1 (method notes)."),
               aliases = "allocation_ratio_n2_n1",
               modes = "a_priori")
    ),
    options = ig_options(.ig_t_tails_option_g1()),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = .ig_power_t_means_two_independent,
    design = .ig_design_t_means_two_independent,
    discrete = FALSE,
    effect_direction = "increasing",
    engine = .ig_t_engine_g1(),
    checker = "validation/independent/t/check_t_means_two_independent.R",
    reference = list(
      gpower_path = "Tests > Means > Two independent groups",
      gpower_label = "t tests - Means: Difference between two independent means (two groups)",
      manual_pages = "ch. 21 (PDF p. 52); examples 1.2.3-1.2.4 (PDF pp. 4-6); tail convention 25.5 (PDF p. 63)",
      allocation_rule = paste(
        "allocation ratio r = N2/N1: N = the smallest even N >= 4 whose power with fractional group sizes",
        "N/(1 + r), rN/(1 + r) reaches the target (G*Power 3.1.9.7), split as n1 = floor(N/(1 + r) + 1/2),",
        "n2 = N - n1; where that design misses the target or leaves a group empty, the next even N whose split",
        "reaches the target (validation/compare/t_means_two_independent/allocation-rule.md); r = 1: n1 = n2 = N/2"),
      method_notes = "specs/method-notes/t_means_two_independent.md"
    ),
    status = "independently_checked",
    constraints = function(p) {
      out <- .ig_t_power_alpha_constraint_g1(p)
      if (!is.null(p$n1) && !is.null(p$n2) && as.double(p$n1) + as.double(p$n2) < 3) {
        out <- c(out, "n1 + n2 must be at least 3 (df = n1 + n2 - 2 >= 1)")
      }
      out
    },
    a_priori_search = .ig_a_priori_t_means_two_independent
  )
}

# Integer split of a total N under allocation ratio r = N2/N1 (G*Power 3.1.9.7,
# allocation-rule.md): even N only; n1 = floor(N/(1 + r) + 1/2) (round half
# up), n2 = N - n1. NULL for an odd N or when a group would be empty.
# `continuous = TRUE` returns the fractional split used by the fractional
# criterion, the continuous relaxation and G*Power's X-Y plot tables
# (documented examples DOC-005 to DOC-017).
.ig_design_t_means_two_independent <- function(p, n_total, continuous = FALSE) {
  r <- as.double(p$allocation_ratio)
  if (continuous) return(c(n1 = n_total / (1 + r), n2 = r * n_total / (1 + r)))
  if (n_total %% 2 != 0) return(NULL)
  n1 <- if (r == 1) n_total / 2 else floor(n_total / (1 + r) + 0.5)
  n2 <- n_total - n1
  if (n1 < 1 || n2 < 1) return(NULL)
  c(n1 = n1, n2 = n2)
}

# Forward model.
.ig_power_t_means_two_independent <- function(p) {
  # Group sizes arrive as integers (their domain type); products such as
  # n1 * n2 overflow the integer range, so all arithmetic is in doubles.
  n1 <- as.double(p$n1)
  n2 <- as.double(p$n2)
  df <- n1 + n2 - 2
  ncp <- p$d * sqrt(n1 * n2 / (n1 + n2))
  core <- .ig_t_power_g1(df, ncp, p$alpha, p$tails)
  list(
    power = core$power,
    beta = core$beta,
    alpha_attained = p$alpha,
    critical = core$critical,
    df = df,
    ncp = ncp,
    n_total = n1 + n2,
    group_sizes = c(n1 = n1, n2 = n2),
    statistic_distribution = list(
      h0 = list(family = "t", df = df, ncp = 0),
      h1 = list(family = "t", df = df, ncp = ncp),
      rejection_side = core$side
    ),
    assumptions = c(
      "two independent samples from normal populations with equal variances (pooled-variance Student t test; G*Power offers no Welch variant)",
      "d = (mu1 - mu2) / sigma; df = n1 + n2 - 2; ncp = d * sqrt(n1 * n2 / (n1 + n2)) (manual 21.5)",
      .ig_t_tail_note_g1(p$tails, "d")
    )
  )
}

# Smallest even N in [lo, hi] with ok(N) TRUE, for a predicate that is
# monotone in N (FALSE ... FALSE TRUE ... TRUE over even N); NA when ok(hi)
# is FALSE. `start` (an even N, optional) seeds the doubling.
.ig_smallest_even_g1 <- function(ok, lo, hi, start = NA_real_) {
  even_up <- function(x) 2 * ceiling(x / 2)
  lo <- even_up(lo)
  hi <- 2 * floor(hi / 2)
  if (lo > hi) return(NA_real_)
  if (ok(lo)) return(lo)
  L <- lo
  H <- if (is.finite(start)) min(max(even_up(start), lo + 2), hi) else min(2 * lo, hi)
  while (!ok(H)) {
    if (H >= hi) return(NA_real_)
    L <- H
    H <- min(even_up(2 * H), hi)
  }
  while (H - L > 2) {
    M <- L + 2 * floor((H - L) / 4)
    if (ok(M)) H <- M else L <- M
  }
  H
}

# A priori search of t_means_two_independent (ig_method(a_priori_search =)).
.ig_a_priori_t_means_two_independent <- function(p, n_lo, n_max) {
  r <- as.double(p$allocation_ratio)
  target <- p$power
  at <- function(n1, n2) {
    q <- p
    q$n1 <- n1
    q$n2 <- n2
    .ig_power_t_means_two_independent(q)$power
  }
  frac_power <- function(N) at(N / (1 + r), r * N / (1 + r))

  # Continuous relaxation (fractional group sizes): a diagnostic, and the
  # start of the searches. Same construction as the generic solver.
  relax <- tryCatch({
    f <- function(u) frac_power(exp(u)) - target
    br <- .ig_bracket_monotone(f, log(max(n_lo, 2)), log(n_lo), log(n_max), decreasing = FALSE)
    if (isTRUE(br$found)) {
      rr <- if (isTRUE(br$exact)) list(x = br$a, kind = "crossing") else .ig_refine_root(f, br$a, br$b, br$fa, br$fb)
      if (rr$kind %in% c("crossing", "steep_crossing")) exp(rr$x) else NA_real_
    } else NA_real_
  }, error = function(e) NA_real_)
  relax_diag <- if (is.finite(relax)) list(continuous_relaxation = list(
    n_total = relax,
    note = "root of the power curve with fractional group sizes; not an admissible design, used only to start the integer search"
  )) else list()

  if (r == 1) {
    # Equal allocation: fractional and integer designs coincide; the generic
    # integer search (bisection, predecessor and minimality evidence).
    design_at <- function(N) .ig_design_t_means_two_independent(p, N)
    search <- .ig_integer_search(
      power_at = function(N) {
        gs <- design_at(N)
        at(gs[["n1"]], gs[["n2"]])
      },
      admissible = function(N) !is.null(design_at(N)),
      n_lo = n_lo, n_max = n_max, target = target, discrete = FALSE, n_start = relax)
    search$diagnostics <- c(search$diagnostics, relax_diag)
    return(search)
  }

  # Unequal allocation: G*Power's fractional criterion over even N.
  lo <- max(4, n_lo)
  n_star <- .ig_smallest_even_g1(function(N) frac_power(N) >= target, lo, n_max, start = relax)
  if (is.na(n_star)) {
    p_lim <- frac_power(2 * floor(n_max / 2))
    flat <- abs(p_lim - frac_power(lo)) <= 1e-14
    return(list(status = "infeasible", reason = sprintf(paste(
      "the target power %s is not reached with fractional group sizes for any even N <= %s",
      "(power there: %s)%s"), format(target), format(2 * floor(n_max / 2), scientific = FALSE),
      format(signif(p_lim, 7)), if (flat) "; power does not change with N" else ""),
      diagnostics = c(list(allocation_ratio = r), relax_diag)))
  }
  split <- function(N) {
    n1 <- floor(N / (1 + r) + 0.5)
    c(n1 = n1, n2 = N - n1)
  }
  int_ok <- function(N) {
    s <- split(N)
    s[["n1"]] >= 1 && s[["n2"]] >= 1 && at(s[["n1"]], s[["n2"]]) >= target
  }
  gs_star <- split(n_star)
  empty <- gs_star[["n1"]] < 1 || gs_star[["n2"]] < 1
  power_star <- if (empty) NA_real_ else at(gs_star[["n1"]], gs_star[["n2"]])
  gp_rule <- list(
    n_total = n_star, n1 = unname(gs_star[["n1"]]), n2 = unname(gs_star[["n2"]]),
    fractional_power = frac_power(n_star),
    fractional_power_predecessor = if (n_star - 2 >= lo) frac_power(n_star - 2) else NA_real_,
    integer_design_power = power_star,
    empty_group = empty,
    reaches_target = !empty && power_star >= target
  )
  diag <- c(list(
    search = paste("G*Power 3.1.9.7 rule for unequal allocation: smallest even N >= 4 whose fractional-allocation",
                   "power reaches the target (even-N bisection; that power increases with N), split as",
                   "n1 = floor(N / (1 + r) + 1/2), n2 = N - n1"),
    allocation_ratio = r,
    n_lower = lo, n_upper_limit = n_max
  ), relax_diag)
  warnings <- character()
  if (gp_rule$reaches_target) {
    n <- n_star
    diag$minimality <- if (n_star - 2 >= lo) {
      sprintf(paste("proven for the fractional criterion: N = %s is the smallest even N whose fractional-allocation",
                    "power reaches the target (fractional power at N - 2 = %s)"),
              format(n_star, scientific = FALSE), format(signif(gp_rule$fractional_power_predecessor, 10)))
    } else {
      "proven for the fractional criterion: N = 4 is the smallest admissible even N"
    }
  } else {
    # The rounded split is non-decreasing in N in both groups, so its power is
    # non-decreasing over even N and the bisection is exact.
    n <- .ig_smallest_even_g1(int_ok, n_star + 2, n_max)
    if (is.na(n)) {
      return(list(status = "infeasible", reason = sprintf(paste(
        "G*Power's allocation rule gives N = %s (n1 = %s, n2 = %s), which %s, and no even N up to %s",
        "has a rounded split with both groups non-empty that reaches the target power %s"),
        format(n_star, scientific = FALSE), format(gs_star[["n1"]]), format(gs_star[["n2"]]),
        if (empty) "leaves a group empty" else "misses the target power",
        format(2 * floor(n_max / 2), scientific = FALSE), format(target)),
        diagnostics = c(diag, list(gpower_allocation_rule = gp_rule))))
    }
    reason <- if (empty) {
      sprintf("leaves group %s empty (G*Power 3.1.9.7 then displays an actual power equal to alpha without a message)",
              if (gs_star[["n1"]] < 1) "1" else "2")
    } else {
      sprintf("has power %s at its integer design, below the target %s (G*Power 3.1.9.7 reports it without a warning)",
              format(signif(power_star, 7)), format(target))
    }
    warnings <- sprintf(paste(
      "G*Power's allocation rule for N2/N1 = %s gives N = %s (n1 = %s, n2 = %s), which %s; inferguidR reports",
      "the smallest larger even N whose rounded split reaches the target (see diagnostics$gpower_allocation_rule)"),
      .ig_format_double(r), format(n_star, scientific = FALSE), format(gs_star[["n1"]]), format(gs_star[["n2"]]),
      reason)
    diag$minimality <- sprintf(paste(
      "N = %s is the smallest even N above G*Power's N = %s whose rounded split has both groups non-empty and",
      "reaches the target (rounded-split power is non-decreasing over even N)"),
      format(n, scientific = FALSE), format(n_star, scientific = FALSE))
    diag$adjusted_from_gpower_rule <- TRUE
    diag$adjustment_note <- warnings
  }
  # Transparency: the smallest even N whose rounded split reaches the target
  # (integer-design criterion), which can lie below G*Power's N.
  n_int <- .ig_smallest_even_g1(int_ok, lo, n, start = relax)
  if (!is.na(n_int) && n_int != n) {
    diag$integer_criterion_n_total <- n_int
    diag$integer_criterion_note <- paste(
      "smallest even N whose rounded split reaches the target; G*Power's rule chooses N by the fractional-allocation",
      "power instead, which can give a larger N")
  }
  pred <- n - 2
  if (pred >= lo) {
    diag$predecessor_n_total <- pred
    gp <- split(pred)
    diag$predecessor_power <- if (gp[["n1"]] >= 1 && gp[["n2"]] >= 1) at(gp[["n1"]], gp[["n2"]]) else NA_real_
  }
  diag$gpower_allocation_rule <- gp_rule
  list(status = "ok", n = n, power = at(split(n)[["n1"]], split(n)[["n2"]]), diagnostics = diag,
       warnings = warnings)
}
