# Method definitions: variance tests and goodness of fit (group G11) ------------------
#
# Ownership: R/methods-variances.R belongs to the G11 task. It defines three
# procedures of the reference inventory:
#
#   f_variances_two_groups            "Variance: Test of equality (two sample case)" (F tests)
#   chi2_variance_one_sample          "Variance: Difference from constant (one sample case)" (chi2 tests)
#   chi2_goodness_of_fit_contingency  "Goodness-of-fit tests: Contingency tables" (chi2 tests)
#
# Providers (collected by R/core-registry.R): .ig_define_methods_F_g11() and
# .ig_define_methods_chi2_g11(). Effect-size drawers: R/effect-variances.R.
# Mathematical contracts, conventions and evidence: specs/method-notes/<id>.md.
# Fixture comparisons: validation/compare/<id>/; independent checks:
# validation/independent/g11/.

.ig_define_methods_F_g11 <- function() {
  list(.ig_method_f_variances_two_groups())
}

.ig_define_methods_chi2_g11 <- function() {
  list(.ig_method_chi2_variance_one_sample(), .ig_method_chi2_goodness_of_fit_contingency())
}

# Shared helpers --------------------------------------------------------------------------

# Quantile of a central distribution from its R quantile function, verified by
# the round trip through the CDF and refined by root finding when the round
# trip misses. nmath's qf() returns df2 / qchisq(., df2) as soon as df1 > 4e5
# (and qchisq(.) / df1 for df2 > 4e5), a large-df approximation: at df1 =
# 1647266, df2 = 164726 the upper 0.001 quantile from qf() has upper-tail
# probability 0.0016 under pf() (TOMS 708, exact to double precision). G*Power
# displays the exact quantile there (campaign record GPW-F-VAR2-APR-B001-D01).
# The refinement solves log P(X <= x) = log p (or the upper tail) on the log(x)
# scale, so tiny quantiles keep their relative accuracy.
#
# A refined value is returned only after it has been verified through the tail
# probability. A quantile outside the range of normal doubles cannot be
# represented: the search then ends at the range limit (the largest double, 0
# or a subnormal number), which is not the quantile. Examples: the upper 1e-200
# quantile of F(1, 1) is about 4e399, the lower one about 2.5e-400, and the
# upper 0.05 quantile of F(2, 0.001) exceeds 1e1900. NaN is returned there, and
# the forward models report NA power with a reason instead of a power computed
# from a wrong critical value.
#   q0       start value (the R quantile)
#   logtail  function(x, lower_tail) returning the log tail probability
.ig_g11_quantile <- function(p, q0, logtail, lower_tail) {
  if (!is.finite(p) || p <= 0 || p >= 1) return(NaN)
  lp <- log(p)
  h <- function(u) logtail(exp(u), lower_tail) - lp
  verified <- function(u) {
    x <- exp(u)
    r <- logtail(x, lower_tail) - lp
    if (is.finite(x) && x >= .Machine$double.xmin && is.finite(r) && abs(r) <= 1e-6) x else NaN
  }
  if (is.finite(q0) && q0 > 0) {
    r0 <- h(log(q0))
    if (is.finite(r0) && abs(r0) <= 1e-12) return(q0)
    u0 <- log(q0)
  } else {
    u0 <- 0
    r0 <- h(u0)
  }
  # h is increasing in u for the lower tail and decreasing for the upper tail.
  sgn <- if (lower_tail) 1 else -1
  g <- function(u) sgn * h(u)
  a <- u0
  b <- u0
  ga <- sgn * r0
  gb <- ga
  step <- 1e-3
  for (i in seq_len(200L)) {
    if (!is.na(ga) && !is.na(gb) && ga <= 0 && gb >= 0) break
    if (is.na(ga) || ga > 0) {
      a <- a - step
      ga <- g(a)
    }
    if (is.na(gb) || gb < 0) {
      b <- b + step
      gb <- g(b)
    }
    step <- step * 2
  }
  if (!(!is.na(ga) && !is.na(gb) && ga <= 0 && gb >= 0)) return(NaN)
  if (ga == 0) return(verified(a))
  if (gb == 0) return(verified(b))
  if (is.finite(ga) && is.finite(gb)) {
    return(verified(stats::uniroot(g, c(a, b), f.lower = ga, f.upper = gb, tol = 1e-15, maxiter = 2000L)$root))
  }
  # An infinite end value (exp(u) beyond the double range at a bracket end):
  # plain bisection. It converges to the quantile when the quantile is
  # representable and to the range limit otherwise; verified() tells them apart.
  for (i in seq_len(400L)) {
    m <- (a + b) / 2
    gm <- g(m)
    if (is.na(gm)) return(NaN)
    if (gm <= 0) a <- m else b <- m
    if (b - a <= 1e-15 * max(1, abs(a))) break
  }
  verified((a + b) / 2)
}

.ig_g11_qf <- function(p, df1, df2, lower_tail) {
  .ig_g11_quantile(p, stats::qf(p, df1, df2, lower.tail = lower_tail),
                   function(x, lt) stats::pf(x, df1, df2, lower.tail = lt, log.p = TRUE), lower_tail)
}

.ig_g11_qchisq <- function(p, df, lower_tail) {
  .ig_g11_quantile(p, stats::qchisq(p, df, lower.tail = lower_tail),
                   function(x, lt) stats::pchisq(x, df, lower.tail = lt, log.p = TRUE), lower_tail)
}

# Power of a test whose statistic is S = ratio * X under H1 and S = X under H0,
# X a central distribution with CDF `ptail(q, lower_tail)` and quantile
# function `qtail(p, lower_tail)` (both on the X scale). Two tails: alpha/2 in
# each tail, both rejection tails counted. One tail: the rejection side is
# `side` ("upper" or "lower"). beta is the non-rejection probability computed
# directly (a difference of same-side tails for two tails), not 1 - power.
.ig_g11_scaled_power <- function(ratio, alpha, tails, side, ptail, qtail) {
  # P(ratio * X <= q) and P(ratio * X > q); ratio = 0 puts all mass at 0.
  ph1 <- function(q, lt) {
    if (ratio == 0) return(as.double(if (lt) q >= 0 else q < 0))
    ptail(q / ratio, lt)
  }
  if (identical(tails, "two")) {
    crit <- c(qtail(alpha / 2, TRUE), qtail(alpha / 2, FALSE))
    if (anyNA(crit) || crit[1] > crit[2]) {
      return(list(power = NA_real_, beta = NA_real_, critical = crit, side = "both"))
    }
    region <- .ig_two_sided_prob(ph1, crit[1], crit[2])
    return(list(power = region$reject, beta = region$accept, critical = crit, side = "both"))
  }
  upper <- identical(side, "upper")
  crit <- qtail(alpha, !upper)
  if (is.na(crit)) return(list(power = NA_real_, beta = NA_real_, critical = NA_real_, side = side))
  list(power = ph1(crit, !upper), beta = ph1(crit, upper), critical = crit, side = side)
}

# Reason given for an NA power: a critical value that could not be computed.
.ig_g11_no_critical_note <- function(h0) {
  paste0("the critical values of ", h0, " could not be computed: a critical value lies outside the range of ",
         "double-precision numbers (alpha too small for these degrees of freedom)")
}

# G*Power refuses a target power below alpha in a priori and sensitivity mode
# ("Power must be greater than or equal to alpha."); the constraint applies
# whenever both are inputs.
.ig_g11_constraints_power_alpha <- function(p) {
  if (!is.null(p$alpha) && !is.null(p$power) && p$power < p$alpha) {
    sprintf(paste("the target power %s is below alpha %s; every design of this test has power of at",
                  "least about alpha, so no sample size or effect size is defined by it"),
            .ig_format_double(p$power), .ig_format_double(p$alpha))
  } else {
    character()
  }
}

# Variance: Test of equality (two sample case) ------------------------------------------
#
# G*Power: Tests > Variances > Two groups; "F tests - Variance: Test of
# equality (two sample case)" (manual ch. 15). Groups 1 and 2 with variances
# var1 and var0, sample sizes n1 and n2. Test statistic F = s1^2 / s2^2
# (numerator = group 1, as the live captures show: Numerator df = n1 - 1).
#
#   H0 distribution F(n1 - 1, n2 - 1); H1 distribution ratio * F(n1 - 1, n2 - 1),
#   ratio = var1 / var0                                      (manual 15.5)
#
# Two tails: alpha/2 in each tail, both rejection tails counted. One tail: the
# rejection region is in the upper tail for ratio >= 1 and in the lower tail
# for ratio < 1 (live capture SUP-G11-F-VAR2-POS-P01: ratio 1 uses the upper
# tail). In sensitivity mode the option "Effect direction" (var1 <= var0 or
# var1 >= var0) selects on which side of 1 the solved ratio lies.
#
# A priori allocation (live captures, specs/method-notes/f_variances_two_groups.md):
#   r = N2/N1 = 1: equal groups n1 = n2 = n, the smallest n reaching the target;
#   r != 1: N* = the smallest integer N whose design with fractional group sizes
#     N/(1 + r), rN/(1 + r) reaches the target power; n1 = floor(N*/(1 + r) + 1/2),
#     n2 = N* - n1. Critical values, df and "Actual power" are those of the
#     integer design, whose power can be below the target (G*Power reports it
#     without a message). inferguidR reports the design when it reaches the
#     target; otherwise it enumerates subsequent admissible rounded designs
#     and reports the first reaching the target, with an explicit warning.

.ig_method_f_variances_two_groups <- function() {
  ig_method(
    id = "f_variances_two_groups",
    family = "F",
    label = list(en = "Variance: Test of equality (two sample case)",
                 de = "Varianz: Test auf Gleichheit (zwei Stichproben)"),
    params = ig_params(
      ig_param("ratio", role = "effect",
               domain = ig_domain_real(0, Inf, lower_inclusive = TRUE, upper_inclusive = FALSE,
                                       scan = c(1e-12, 1e12)),
               label = list(en = "Ratio var1/var0", de = "Verh\u00e4ltnis var1/var0"),
               description = paste("Variance ratio var1 / var0 of group 1 to group 2 under H1 (manual 15.1);",
                                   "the test statistic s1^2 / s2^2 follows ratio * F(n1 - 1, n2 - 1)."),
               aliases = "ratio_var1_var0"),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      ig_param("n_total", role = "design", domain = ig_domain_integer(lower = 4), unit = "subjects",
               label = list(en = "Total sample size", de = "Gesamtstichprobenumfang"),
               description = "Solved in a priori mode (n1 + n2).",
               modes = "a_priori"),
      ig_param("n1", role = "design", domain = ig_domain_integer(lower = 2), unit = "subjects",
               label = list(en = "Sample size group 1", de = "Stichprobenumfang Gruppe 1"),
               aliases = "sample_size_group_1",
               modes = c("post_hoc", "sensitivity", "criterion", "compromise")),
      ig_param("n2", role = "design", domain = ig_domain_integer(lower = 2), unit = "subjects",
               label = list(en = "Sample size group 2", de = "Stichprobenumfang Gruppe 2"),
               aliases = "sample_size_group_2",
               modes = c("post_hoc", "sensitivity", "criterion", "compromise")),
      ig_param("allocation_ratio", role = "design",
               domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE),
               default = 1,
               label = list(en = "Allocation ratio N2/N1", de = "Allokationsverh\u00e4ltnis N2/N1"),
               description = paste("Group-size ratio n2/n1 of the a priori design. Default 1 = G*Power's",
                                   "initial field value (certified live inventory)."),
               aliases = "allocation_ratio_n2_n1",
               modes = "a_priori")
    ),
    options = ig_options(
      ig_option("tails", values = c("one", "two"), default = "one", labels = c(one = "One", two = "Two"),
                label = list(en = "Tail(s)", de = "Seiten")),
      ig_option("effect_direction", values = c("var1_le_var0", "var1_ge_var0"), default = "var1_ge_var0",
                labels = c(var1_le_var0 = "var1\u2264var0", var1_ge_var0 = "var1\u2265var0"),
                applies_when = list(mode = "sensitivity"),
                label = list(en = "Effect direction", de = "Effektrichtung"))
    ),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = .ig_power_f_variances_two_groups,
    design = .ig_design_f_variances_two_groups,
    discrete = FALSE,
    # Select the requested crossing while evaluating the full power curve.
    effect_direction = function(p) {
      if (identical(p$effect_direction, "var1_le_var0")) "decreasing" else "increasing"
    },
    engine = list(package = "stats",
                  fun = paste("pf (nmath, TOMS 708 pbeta); critical values from qf, verified through pf",
                              "and refined by root finding where qf's large-df approximation misses"),
                  version = NA_character_),
    checker = "validation/independent/g11/check_g11.R",
    reference = list(
      gpower_path = "Tests > Variances > Two groups",
      gpower_label = "F tests - Variance: Test of equality (two sample case)",
      manual_pages = "ch. 15 (PDF p. 39)",
      allocation_rule = paste(
        "allocation ratio 1: equal groups, the smallest n per group reaching the target power;",
        "other ratios (G*Power 3.1.9.7, live captures): N* = smallest integer N whose design with",
        "fractional group sizes N/(1 + r), rN/(1 + r) reaches the target, n1 = floor(N*/(1 + r) + 1/2),",
        "n2 = N* - n1 (both at least 2); if rounding misses the target, enumerate subsequent",
        "admissible totals and report the first reaching the target, with a warning"),
      method_notes = "specs/method-notes/f_variances_two_groups.md"
    ),
    status = "independently_checked",
    constraints = .ig_g11_constraints_power_alpha,
    a_priori_search = .ig_apriori_f_variances_two_groups
  )
}

# Allocation rule (see the header). `continuous = TRUE` gives the fractional
# group sizes of G*Power's search criterion.
.ig_design_f_variances_two_groups <- function(p, n_total, continuous = FALSE) {
  r <- as.double(p$allocation_ratio)
  N <- as.double(n_total)
  if (continuous) return(c(n1 = N / (1 + r), n2 = r * N / (1 + r)))
  if (r == 1) {
    if (N %% 2 != 0 || N < 4) return(NULL)
    return(c(n1 = N / 2, n2 = N / 2))
  }
  n1 <- floor(N / (1 + r) + 0.5)
  n2 <- N - n1
  if (n1 < 2 || n2 < 2) return(NULL)
  c(n1 = n1, n2 = n2)
}

.ig_power_f_variances_two_groups <- function(p) {
  n1 <- as.double(p$n1)
  n2 <- as.double(p$n2)
  ratio <- as.double(p$ratio)
  d1 <- n1 - 1
  d2 <- n2 - 1
  base <- list(alpha_attained = p$alpha, df = c(d1, d2), n_total = n1 + n2, group_sizes = c(n1 = n1, n2 = n2),
               outputs = list(numerator_df = d1, denominator_df = d2))
  side <- if (ratio >= 1) "upper" else "lower"
  core <- .ig_g11_scaled_power(
    ratio, p$alpha, p$tails, side,
    ptail = function(q, lt) stats::pf(q, d1, d2, lower.tail = lt),
    qtail = function(pr, lt) .ig_g11_qf(pr, d1, d2, lt))
  notes <- list()
  if (is.na(core$power)) {
    notes$power <- notes$beta <- notes$critical <- .ig_g11_no_critical_note("F(n1 - 1, n2 - 1)")
  }
  tail_note <- if (identical(p$tails, "two")) {
    "two-sided test: alpha/2 in each tail of F(n1 - 1, n2 - 1); the power counts both rejection tails"
  } else {
    paste("one-sided test: upper rejection tail for ratio >= 1, lower tail for ratio < 1",
          "(G*Power 3.1.9.7; ratio = 1 uses the upper tail, where the power equals alpha)")
  }
  c(list(
    power = core$power,
    beta = core$beta,
    critical = core$critical,
    notes = notes,
    statistic_distribution = list(
      h0 = list(family = "F", df1 = d1, df2 = d2),
      h1 = list(family = "scaled_F", df1 = d1, df2 = d2, scale = ratio,
                note = "the statistic s1^2 / s2^2 is ratio * F(df1, df2) under H1"),
      rejection_side = core$side
    ),
    assumptions = c(
      "two independent samples from normal populations; F = s1^2 / s2^2 with group 1 in the numerator",
      "H0: var1 = var0; H1: var1 / var0 = ratio; the statistic follows ratio * F(n1 - 1, n2 - 1) (manual 15.5)",
      tail_note)
  ), base)
}

# A priori search (G*Power 3.1.9.7 allocation rule; see the header).
.ig_apriori_f_variances_two_groups <- function(p, n_lo, n_max) {
  r <- as.double(p$allocation_ratio)
  target <- p$power
  fw <- function(n1, n2) .ig_power_f_variances_two_groups(c(p, list(n1 = n1, n2 = n2)))$power
  if (r == 1) {
    s <- .ig_integer_search(power_at = function(N) fw(N / 2, N / 2),
                            admissible = function(N) N %% 2 == 0 && N >= 4,
                            n_lo = n_lo, n_max = n_max, target = target, discrete = FALSE)
    s$diagnostics$allocation_rule <- "allocation ratio 1: n1 = n2 = N/2, the smallest even N whose power reaches the target"
    return(s)
  }
  # Smallest N at which both fractional group sizes exceed 1 (df > 0) ...
  n_frac_lo <- max(n_lo, floor(max(1 + r, (1 + r) / r)) + 1)
  fp <- function(N) fw(N / (1 + r), r * N / (1 + r))
  # ... and at which the criterion can be evaluated. A fractional group size
  # barely above 1 gives a df close to 0 (r = 0.2, N = 6: n2 = 1 + 2e-16 in
  # double arithmetic; r = 0.3334, N = 4: df 1.5e-4), where the critical value
  # lies outside the double range and the forward model returns NA. The power
  # of such a design tends to alpha as its df tends to 0, so these totals are
  # not the rule's N* for a target above alpha; the search starts at the first
  # total with a computable criterion (df grows with N, so that total is found
  # by doubling and bisection).
  # nmath's accuracy warnings at such a df belong to totals that are skipped;
  # the search evaluates the totals it uses again.
  uncomputable <- function(N) is.na(suppressWarnings(fp(N)))
  n_frac_first <- n_frac_lo
  if (uncomputable(n_frac_lo)) {
    lo <- n_frac_lo
    step <- 1
    repeat {
      hi <- min(lo + step, n_max)
      if (!uncomputable(hi)) break
      if (hi >= n_max) {
        return(list(status = "numerical_failure", reason = paste(
          "the power of the design with fractional group sizes N/(1 + r), rN/(1 + r) (G*Power's search criterion",
          "for unequal groups) could not be computed for any total sample size:",
          .ig_g11_no_critical_note("F(n1 - 1, n2 - 1)")), diagnostics = list()))
      }
      lo <- hi
      step <- 2 * step
    }
    while (hi - lo > 1) {
      mid <- floor((lo + hi) / 2)
      if (uncomputable(mid)) lo <- mid else hi <- mid
    }
    n_frac_lo <- hi
  }
  s <- .ig_integer_search(power_at = fp, admissible = function(N) N >= n_frac_lo,
                          n_lo = n_frac_lo, n_max = n_max, target = target, discrete = FALSE)
  if (!identical(s$status, "ok")) {
    s$diagnostics$criterion <- "power of the design with fractional group sizes N/(1 + r), rN/(1 + r)"
    return(s)
  }
  n_star <- s$n
  N <- n_star
  split <- NULL
  steps <- 0
  repeat {
    split <- .ig_design_f_variances_two_groups(p, N)
    if (!is.null(split)) break
    steps <- steps + 1
    N <- N + 1
    if (N > n_max || steps > 1000) {
      return(list(status = "infeasible", reason = sprintf(
        "no design with both groups of at least 2 subjects within 1000 totals above N = %s under the split n1 = floor(N/(1 + r) + 1/2)",
        format(n_star, scientific = FALSE)), diagnostics = s$diagnostics))
    }
  }
  pw <- fw(split[["n1"]], split[["n2"]])
  diag <- list(
    allocation_rule = paste("G*Power 3.1.9.7 rule for allocation ratio r != 1: N* = smallest integer N whose",
                            "design with fractional group sizes N/(1 + r), rN/(1 + r) reaches the target;",
                            "n1 = floor(N*/(1 + r) + 1/2), n2 = N* - n1"),
    criterion = "power of the design with fractional group sizes (G*Power's search criterion)",
    fractional_search = s$diagnostics,
    n_total_fractional_rule = n_star,
    fractional_power_at_n = fp(n_star),
    fractional_power_at_predecessor = if (n_star - 1 >= n_frac_lo) fp(n_star - 1) else NA_real_,
    integer_power = pw,
    minimality = paste("minimal over the rule's design sequence: the predecessor total has fractional power",
                       "below the target; smaller integer designs are not considered by this rule")
  )
  if (n_frac_lo > n_frac_first) {
    diag$fractional_search_start <- sprintf(paste(
      "the criterion is not computable for N = %s to %s (a fractional group size barely above 1: df close to 0,",
      "critical value outside the double range; the power there is close to alpha); the search starts at N = %s"),
      format(n_frac_first, scientific = FALSE), format(n_frac_lo - 1, scientific = FALSE),
      format(n_frac_lo, scientific = FALSE))
  }
  if (N > n_star) {
    diag$split_adjustment <- sprintf(
      "the split of N* = %s leaves a group with fewer than 2 subjects; the next total with an admissible split, N = %s, is used (G*Power 3.1.9.7 reports no result here)",
      format(n_star, scientific = FALSE), format(N, scientific = FALSE))
  }
  warnings <- character()
  if (is.na(pw)) {
    return(list(status = "numerical_failure", reason = "the rounded allocation has no computable power", diagnostics = diag))
  }
  if (pw < target) {
    initial_N <- N
    initial_power <- pw
    # Enumerate the subsequent rounded splits: this verifies minimality over
    # the declared sequence without assuming monotonicity after rounding.
    repeat {
      predecessor_N <- N
      predecessor_power <- pw
      if (N >= n_max) return(list(status = "infeasible", reason =
        "no subsequent admissible rounded design reaches the target within the sample-size bound", diagnostics = diag))
      N <- N + 1
      split <- .ig_design_f_variances_two_groups(p, N)
      if (is.null(split)) next
      pw <- fw(split[["n1"]], split[["n2"]])
      if (is.na(pw)) return(list(status = "numerical_failure", reason =
        "power could not be evaluated while checking subsequent rounded designs", diagnostics = diag))
      if (pw >= target) break
    }
    warnings <- sprintf(paste(
      "The reference allocation at N = %s has power %s below the target %s after rounding;",
      "the first subsequent admissible rounded design reaching the target, N = %s, is reported."),
      format(initial_N), format(initial_power, digits = 10), format(target), format(N))
    diag$adjusted_from_gpower_rule <- TRUE
    diag$unadjusted_integer_power <- initial_power
    diag$integer_power <- pw
    diag$adjustment_predecessor <- list(n_total = predecessor_N, power = predecessor_power)
    diag$adjustment_enumerated_range <- c(initial_N, N)
    diag$minimality <- "first target-reaching admissible rounded design at or above the reference fractional-rule total; subsequent totals checked by enumeration"
  }
  list(status = "ok", n = N, power = pw, diagnostics = diag, warnings = warnings)
}

# Variance: Difference from constant (one sample case) ----------------------------------
#
# G*Power: Tests > Variances > One group; "chi2 tests - Variance: Difference
# from constant (one sample case)" (manual ch. 26). Sample of size N from a
# normal population with variance var1, tested against the constant var0.
#
#   statistic (N - 1) s^2 / var0; H0 distribution chi2(N - 1),
#   H1 distribution ratio * chi2(N - 1), ratio = var1 / var0     (manual 26.5)
#
# Two tails: alpha/2 in each tail. One tail: upper rejection tail for
# ratio > 1, lower tail for ratio <= 1 (live capture SUP-G11-CHI2-VAR1-POS-P01:
# ratio 1 uses the lower tail, unlike the two-sample F test). Sensitivity mode
# has no effect-direction option; G*Power reports the root above 1.

.ig_method_chi2_variance_one_sample <- function() {
  ig_method(
    id = "chi2_variance_one_sample",
    family = "chi2",
    label = list(en = "Variance: Difference from constant (one sample case)",
                 de = "Varianz: Abweichung von einer Konstanten (eine Stichprobe)"),
    params = ig_params(
      ig_param("ratio", role = "effect",
               domain = ig_domain_real(0, Inf, lower_inclusive = TRUE, upper_inclusive = FALSE),
               label = list(en = "Ratio var1/var0", de = "Verh\u00e4ltnis var1/var0"),
               description = paste("Ratio var1 / var0 of the population variance under H1 to the constant",
                                   "var0 (manual 26.1)."),
               aliases = "ratio_var1_var0"),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      ig_param("n_total", role = "design", domain = ig_domain_integer(lower = 2), unit = "subjects",
               label = list(en = "Total sample size", de = "Gesamtstichprobenumfang"),
               aliases = "total_sample_size")
    ),
    options = ig_options(
      ig_option("tails", values = c("one", "two"), default = "one", labels = c(one = "One", two = "Two"),
                label = list(en = "Tail(s)", de = "Seiten"))
    ),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = .ig_power_chi2_variance_one_sample,
    sensitivity_search = .ig_sensitivity_chi2_variance_one_sample,
    design = NULL,
    discrete = FALSE,
    # Power falls on (0, 1] and rises on (1, Inf); "increasing" selects the
    # root above 1, as G*Power's sensitivity mode does.
    effect_direction = "increasing",
    engine = list(package = "stats",
                  fun = "pchisq (nmath pgamma); critical values from qchisq, verified through pchisq",
                  version = NA_character_),
    checker = "validation/independent/g11/check_g11.R",
    reference = list(
      gpower_path = "Tests > Variances > One group",
      gpower_label = "chi2 tests - Variance: Difference from constant (one sample case)",
      manual_pages = "ch. 26 (PDF p. 64)",
      method_notes = "specs/method-notes/chi2_variance_one_sample.md"
    ),
    status = "independently_checked",
    constraints = .ig_g11_constraints_power_alpha
  )
}

.ig_power_chi2_variance_one_sample <- function(p) {
  N <- as.double(p$n_total)
  ratio <- as.double(p$ratio)
  df <- N - 1
  side <- if (ratio > 1) "upper" else "lower"
  core <- .ig_g11_scaled_power(
    ratio, p$alpha, p$tails, side,
    ptail = function(q, lt) stats::pchisq(q, df, lower.tail = lt),
    qtail = function(pr, lt) .ig_g11_qchisq(pr, df, lt))
  notes <- list()
  if (is.na(core$power)) notes$power <- notes$beta <- notes$critical <- .ig_g11_no_critical_note("chi2(N - 1)")
  tail_note <- if (identical(p$tails, "two")) {
    "two-sided test: alpha/2 in each tail of chi2(N - 1); the power counts both rejection tails"
  } else {
    paste("one-sided test: upper rejection tail for ratio > 1, lower tail for ratio <= 1",
          "(G*Power 3.1.9.7; ratio = 1 uses the lower tail, where the power equals alpha)")
  }
  list(
    power = core$power,
    beta = core$beta,
    alpha_attained = p$alpha,
    critical = core$critical,
    df = df,
    n_total = N,
    notes = notes,
    statistic_distribution = list(
      h0 = list(family = "chisq", df = df, ncp = 0),
      h1 = list(family = "scaled_chisq", df = df, scale = ratio,
                note = "the statistic (N - 1) s^2 / var0 is ratio * chi2(df) under H1"),
      rejection_side = core$side
    ),
    assumptions = c(
      "one sample from a normal population; statistic (N - 1) s^2 / var0",
      "H0: var1 = var0; H1: var1 / var0 = ratio; the statistic follows ratio * chi2(N - 1) (manual 26.5)",
      tail_note)
  )
}

# Goodness-of-fit tests: Contingency tables -----------------------------------------------
#
# G*Power: Tests > Proportions > Multigroup: Goodness-of-Fit; "chi2 tests -
# Goodness-of-fit tests: Contingency tables" (no manual chapter). Pearson's
# chi-square test with Df degrees of freedom; effect size w (Cohen 1988, ch. 7):
# w = sqrt(sum((p1_i - p0_i)^2 / p0_i)) over the cells (for a contingency table
# with p0 the product of the margins of p1, the same formula gives the
# independence test with Df = (rows - 1)(columns - 1)).
#
#   H0 distribution chi2(Df); H1 distribution noncentral chi2(Df, lambda),
#   lambda = N w^2; upper rejection tail only (no tail option).
#
# Df is a real number >= 1 in G*Power (live capture SUP-G11-CHI2-GOF-POS-P01
# accepts Df = 2.5; POS-P02 refuses Df = 0).

.ig_method_chi2_goodness_of_fit_contingency <- function() {
  ig_method(
    id = "chi2_goodness_of_fit_contingency",
    family = "chi2",
    label = list(en = "Goodness-of-fit tests: Contingency tables",
                 de = "Anpassungstests: Kontingenztafeln"),
    params = ig_params(
      ig_param("w", role = "effect",
               domain = ig_domain_real(0, Inf, lower_inclusive = TRUE, upper_inclusive = FALSE),
               label = list(en = "Effect size w", de = "Effektst\u00e4rke w"),
               description = "Cohen's w = sqrt(sum((p1 - p0)^2 / p0)) over the cells (Cohen 1988, ch. 7).",
               aliases = "effect_size_w"),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      ig_param("df", role = "design",
               domain = ig_domain_real(1, Inf, lower_inclusive = TRUE, upper_inclusive = FALSE),
               label = list(en = "Df", de = "Freiheitsgrade"),
               description = paste("Degrees of freedom of the chi-square test (cells - 1 for a goodness-of-fit",
                                   "test, (rows - 1)(columns - 1) for an independence test); a real number >= 1",
                                   "as in G*Power."),
               aliases = "degrees_of_freedom"),
      ig_param("n_total", role = "design", domain = ig_domain_integer(lower = 1), unit = "subjects",
               label = list(en = "Total sample size", de = "Gesamtstichprobenumfang"),
               aliases = "total_sample_size")
    ),
    options = ig_options(),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = .ig_power_chi2_goodness_of_fit_contingency,
    sensitivity_search = .ig_sensitivity_chi2_goodness_of_fit,
    design = NULL,
    discrete = FALSE,
    effect_direction = "increasing",
    engine = list(package = "stats",
                  fun = paste("pchisq/qchisq; noncentral chi-square as a log-scale Poisson mixture of central",
                              "pchisq terms (R/num-dist.R .ig_p_chisq) with a Chernoff bound for tails below",
                              "the double range"),
                  version = NA_character_),
    checker = "validation/independent/g11/check_g11.R",
    reference = list(
      gpower_path = "Tests > Proportions > Multigroup: Goodness-of-Fit",
      gpower_label = "chi2 tests - Goodness-of-fit tests: Contingency tables",
      manual_pages = "no manual chapter (APP-ONLY); Cohen (1988) ch. 7 for w",
      method_notes = "specs/method-notes/chi2_goodness_of_fit_contingency.md"
    ),
    status = "independently_checked",
    constraints = .ig_g11_constraints_power_alpha
  )
}

# Rigorous upper bound on log P(X <= c) for X ~ chi2(k, lambda) (Chernoff):
# for s > 0, P(X <= c) <= exp(s c) E[exp(-s X)] with
# E[exp(-s X)] = (1 + 2 s)^(-k/2) exp(-lambda s / (1 + 2 s)). The exponent is
# minimised over log(s); 0 is returned when no s gives a bound below 1.
.ig_g11_nchisq_lower_log_bound <- function(c, k, lambda) {
  if (c >= k + lambda) return(0)
  f <- function(ls) {
    s <- exp(ls)
    s * c - (k / 2) * log1p(2 * s) - lambda * s / (1 + 2 * s)
  }
  o <- stats::optimize(f, c(-40, 40))
  min(0, o$objective)
}

.ig_power_chi2_goodness_of_fit_contingency <- function(p) {
  N <- as.double(p$n_total)
  w <- as.double(p$w)
  df <- as.double(p$df)
  lambda <- N * w^2
  crit <- .ig_g11_qchisq(p$alpha, df, FALSE)
  notes <- list()
  power <- beta <- NA_real_
  if (!is.finite(crit) || !is.finite(lambda)) {
    notes$power <- notes$beta <- "the critical value or the noncentrality parameter is not finite"
  } else if (lambda == 0) {
    power <- stats::pchisq(crit, df, lower.tail = FALSE)
    beta <- stats::pchisq(crit, df, lower.tail = TRUE)
  } else if (.ig_g11_nchisq_lower_log_bound(crit, df, lambda) < -800) {
    # beta is below exp(-800), far under the smallest double: exactly 0 and 1.
    power <- 1
    beta <- 0
  } else if (lambda > 1e10) {
    notes$power <- notes$beta <- paste("noncentral chi-square with lambda > 1e10 near the centre of its",
                                       "distribution is outside the computational range of the engine")
  } else {
    power <- .ig_p_chisq(crit, df, lambda, lower_tail = FALSE)
    beta <- .ig_p_chisq(crit, df, lambda, lower_tail = TRUE)
  }
  list(
    power = power,
    beta = beta,
    alpha_attained = p$alpha,
    critical = crit,
    df = df,
    ncp = lambda,
    n_total = N,
    notes = notes,
    statistic_distribution = list(
      h0 = list(family = "chisq", df = df, ncp = 0),
      h1 = list(family = "chisq", df = df, ncp = lambda),
      rejection_side = "upper"
    ),
    assumptions = c(
      "Pearson chi-square test of N independent observations classified into cells",
      "H1 distribution noncentral chi2(Df, lambda) with lambda = N w^2 (large-sample approximation of the test statistic)",
      "upper rejection tail only")
  )
}
