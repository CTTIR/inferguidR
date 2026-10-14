# Method definitions: proportions (exact binomial family, McNemar, Fisher, z) -------
#
# Ownership: R/methods-proportions.R belongs to group G8 (proportions, exact
# and z). It defines six registry procedures through two provider functions:
#
#   .ig_define_methods_exact_g8()  exact_proportion_one_sample_binomial,
#                                  exact_proportion_sign_test,
#                                  exact_generic_binomial,
#                                  exact_proportions_mcnemar,
#                                  exact_proportions_fisher
#   .ig_define_methods_z_g8()      z_proportions_two_independent
#
# Method notes (mathematical contracts, G*Power conventions, documented
# deviations): specs/method-notes/<method_id>.md. Effect-size drawers:
# R/effect-proportions.R.
#
# Conventions shared by the exact tests (method notes, section "Rejection
# regions"):
#
# * A rejection region of a test on a binomial count X ~ Bin(N, p) is written
#   (lower, upper): reject when X <= lower or X >= upper. An empty lower tail
#   is lower = -1, an empty upper tail upper = N + 1. Every region has
#   attained size <= nominal alpha (a level-alpha test); no region is ever
#   widened beyond its alpha share (G*Power 3.1.9.7's post hoc routine does
#   that for tiny N; documented as a reference defect).
# * One tail: the region lies on the side of the effect (upper when the H1
#   proportion is >= the H0 proportion).
# * Two tails: G*Power's three alpha-balancing options (manual 4.2) as
#   specified in the method notes of exact_proportion_one_sample_binomial.
# * Power and beta are computed directly from binomial tails (beta never as
#   1 - power), attained alpha from the same region under H0.

.ig_define_methods_exact_g8 <- function() {
  list(
    .ig_method_exact_proportion_one_sample_binomial(),
    .ig_method_exact_proportion_sign_test(),
    .ig_method_exact_generic_binomial(),
    .ig_method_exact_proportions_mcnemar(),
    .ig_method_exact_proportions_fisher()
  )
}

.ig_define_methods_z_g8 <- function() {
  list(.ig_method_z_proportions_two_independent())
}

# Shared constants ---------------------------------------------------------------------

# G*Power's admissible proportion range (eps = 1e-6, manual 4.1).
.ig_prop_eps <- 1e-6

.ig_prop_tails_option <- function() {
  ig_option("tails", values = c("one", "two"), default = "one",
            labels = c(one = "One", two = "Two"),
            label = list(en = "Tail(s)", de = "Seiten"))
}

.ig_prop_balancing_option <- function() {
  ig_option("alpha_balancing",
            values = c("half_each_side", "minor_tail_half_rest_major", "half_each_side_then_increase"),
            default = "half_each_side",
            labels = c(half_each_side = "Assign alpha/2 to both sides",
                       minor_tail_half_rest_major = "Assign to minor tail alpha/2, then rest to major tail",
                       half_each_side_then_increase = "Assign alpha/2 to both sides, then increase to minimize the difference of alpha1 + alpha2 to alpha"),
            applies_when = list(option = "tails", "in" = "two"),
            label = list(en = "Alpha balancing in two-sided tests", de = "Alpha-Aufteilung bei zweiseitigen Tests"))
}

# Inputs that are each admissible but jointly empty: a target power below
# alpha (G*Power: "Power must be greater than or equal to alpha.").
.ig_prop_power_alpha_constraint <- function(p) {
  if (!is.null(p$alpha) && !is.null(p$power) && is.numeric(p$power) && is.numeric(p$alpha) &&
      p$power < p$alpha) {
    return(sprintf(paste("the target power %s is below alpha %s; no level-alpha test is planned for a power",
                         "below its significance level (G*Power refuses the same input)"),
                   .ig_format_double(p$power), .ig_format_double(p$alpha)))
  }
  character()
}

# Binomial tails --------------------------------------------------------------------------

# P(X <= k) and P(X >= k) for X ~ Bin(N, p), vectorised over k, exact 0/1
# outside the support.
.ig_binom_cdf <- function(k, N, p) {
  out <- stats::pbinom(k, N, p)
  out[k < 0] <- 0
  out[k >= N] <- 1
  out
}

.ig_binom_sf <- function(k, N, p) {
  out <- stats::pbinom(k - 1, N, p, lower.tail = FALSE)
  out[k <= 0] <- 1
  out[k > N] <- 0
  out
}

# Largest k in -1..N with P(X <= k) <= a (-1: empty lower tail).
.ig_binom_lower_crit <- function(N, p, a) {
  if (!(a > 0)) return(-1)
  if (a >= 1) return(N)
  k <- stats::qbinom(a, N, p)
  while (k >= 0 && .ig_binom_cdf(k, N, p) > a) k <- k - 1
  while (k + 1 <= N && .ig_binom_cdf(k + 1, N, p) <= a) k <- k + 1
  k
}

# Smallest k in 0..N+1 with P(X >= k) <= a (N + 1: empty upper tail).
.ig_binom_upper_crit <- function(N, p, a) {
  if (!(a > 0)) return(N + 1)
  if (a >= 1) return(0)
  k <- stats::qbinom(a, N, p, lower.tail = FALSE) + 1
  while (k <= N && .ig_binom_sf(k, N, p) > a) k <- k + 1
  while (k - 1 >= 0 && .ig_binom_sf(k - 1, N, p) <= a) k <- k - 1
  k
}

# Probability of the region (lower, upper) under Bin(N, p) and of its
# complement (computed directly, without 1 - x cancellation).
.ig_binom_region_prob <- function(region, N, p) {
  .ig_binom_cdf(region[[1]], N, p) + .ig_binom_sf(region[[2]], N, p)
}

.ig_binom_region_accept <- function(region, N, p) {
  lo <- region[[1]]
  up <- region[[2]]
  if (up - lo <= 1) return(0)
  a <- max(lo + 1, 0)
  b <- min(up - 1, N)
  if (b < a) return(0)
  # P(a <= X <= b): from whichever tail representation keeps both terms small.
  low_form <- .ig_binom_cdf(b, N, p) - .ig_binom_cdf(a - 1, N, p)
  up_form <- .ig_binom_sf(a, N, p) - .ig_binom_sf(b + 1, N, p)
  v <- if (.ig_binom_cdf(b, N, p) <= 0.5) low_form else up_form
  min(max(v, 0), 1)
}

# Rejection region of the binomial test of H0: p = p0 at level alpha for the
# H1 proportion p1 (which fixes the one-sided direction and the option-2
# choice). Returns list(lower, upper, side, alpha_lower, alpha_upper, rule).
.ig_binom_region <- function(N, p0, alpha, tails, balancing, p1) {
  up_dir <- p1 >= p0
  if (identical(tails, "one")) {
    if (up_dir) {
      k <- .ig_binom_upper_crit(N, p0, alpha)
      return(list(lower = -1, upper = k, side = "upper",
                  rule = "one tail on the side of the effect (upper: H1 proportion >= H0 proportion)"))
    }
    k <- .ig_binom_lower_crit(N, p0, alpha)
    return(list(lower = k, upper = N + 1, side = "lower",
                rule = "one tail on the side of the effect (lower: H1 proportion < H0 proportion)"))
  }
  balancing <- balancing %||% "half_each_side"
  kl <- .ig_binom_lower_crit(N, p0, alpha / 2)
  ku <- .ig_binom_upper_crit(N, p0, alpha / 2)
  if (identical(balancing, "half_each_side")) {
    return(list(lower = kl, upper = ku, side = "both", rule = "alpha/2 assigned to each tail (option 1)"))
  }
  if (identical(balancing, "minor_tail_half_rest_major")) {
    # Both assignments of the minor tail are level-alpha regions; G*Power
    # 3.1.9.7 uses the one with the larger power under H1 (probe grid of
    # 800 post hoc captures; method notes). Ties keep the lower tail minor.
    A <- c(kl, .ig_binom_upper_crit(N, p0, alpha - .ig_binom_cdf(kl, N, p0)))
    B <- c(.ig_binom_lower_crit(N, p0, alpha - .ig_binom_sf(ku, N, p0)), ku)
    # larger power = smaller non-rejection probability (computed directly, so
    # that powers near 1 are compared without cancellation)
    useA <- .ig_binom_region_accept(A, N, p1) <= .ig_binom_region_accept(B, N, p1)
    r <- if (useA) A else B
    return(list(lower = r[1], upper = r[2], side = "both",
                minor_tail = if (useA) "lower" else "upper",
                rule = paste("alpha/2 to the minor tail, alpha minus its attained size to the major tail (option 2);",
                             "the minor tail is the one whose assignment gives the larger power")))
  }
  # Option 3 (manual 4.2): start from option 1, then repeatedly widen the tail
  # with the smaller incremental alpha while the attained total stays <= alpha
  # (manual 4.6: the attained alpha never exceeds the nominal one). Equal
  # increments widen the tail on the side of the effect.
  lo <- kl
  up <- ku
  repeat {
    tot <- .ig_binom_cdf(lo, N, p0) + .ig_binom_sf(up, N, p0)
    il <- if (lo + 1 < up) stats::dbinom(lo + 1, N, p0) else Inf
    iu <- if (up - 1 > lo) stats::dbinom(up - 1, N, p0) else Inf
    step <- min(il, iu)
    if (!is.finite(step) || tot + step > alpha) break
    # increments equal up to rounding (for example both tails at pi0 = 1/2)
    # count as a tie, resolved towards the side of the effect
    tie <- is.finite(il) && is.finite(iu) && abs(il - iu) <= 1e-10 * max(il, iu)
    if ((!tie && il < iu) || (tie && !up_dir)) lo <- lo + 1 else up <- up - 1
  }
  list(lower = lo, upper = up, side = "both",
       rule = paste("alpha/2 to each tail, then the tail with the smaller increment is widened while the",
                    "attained total stays <= alpha (option 3, manual 4.2 and 4.6)"))
}

# Forward evaluation of a binomial region test.
.ig_binom_test_power <- function(N, p0, p1, alpha, tails, balancing) {
  r <- .ig_binom_region(N, p0, alpha, tails, balancing, p1)
  region <- c(r$lower, r$upper)
  att <- .ig_binom_region_prob(region, N, p0)
  pw <- .ig_binom_region_prob(region, N, p1)
  beta <- .ig_binom_region_accept(region, N, p1)
  list(region = region, rule = r$rule, side = r$side, alpha_attained = min(att, 1), power = min(pw, 1),
       beta = beta, alpha_lower = .ig_binom_cdf(region[1], N, p0), alpha_upper = .ig_binom_sf(region[2], N, p0),
       minor_tail = r$minor_tail)
}

# Power of the randomised one-sided most powerful level-alpha test of
# p = p0 against the simple alternative p1 (Neyman-Pearson). No level-alpha
# test with N observations has more power at p1, and the MP power is
# non-decreasing in N (an MP test with N + 1 observations can ignore one), so
# the smallest N with MP power >= target is a lower bound for the a priori
# solution of every option.
.ig_binom_mp_power <- function(N, p0, p1, alpha) {
  if (p1 == p0) return(alpha)
  if (p1 > p0) {
    k <- .ig_binom_upper_crit(N, p0, alpha)          # S(k) <= alpha < S(k - 1)
    if (k <= 0) return(1)
    s_k <- .ig_binom_sf(k, N, p0)
    m <- stats::dbinom(k - 1, N, p0)
    gam <- if (m > 0) min(max((alpha - s_k) / m, 0), 1) else 0
    return(min(1, .ig_binom_sf(k, N, p1) + gam * stats::dbinom(k - 1, N, p1)))
  }
  k <- .ig_binom_lower_crit(N, p0, alpha)            # F(k) <= alpha < F(k + 1)
  if (k >= N) return(1)
  f_k <- .ig_binom_cdf(k, N, p0)
  m <- stats::dbinom(k + 1, N, p0)
  gam <- if (m > 0) min(max((alpha - f_k) / m, 0), 1) else 0
  min(1, .ig_binom_cdf(k, N, p1) + gam * stats::dbinom(k + 1, N, p1))
}

# Smallest N >= n_lo where a non-decreasing function of N reaches the target
# (doubling, then bisection); NA if not reached up to n_max.
.ig_monotone_first_n <- function(f, target, n_lo, n_max) {
  if (f(n_lo) >= target) return(n_lo)
  lo <- n_lo
  hi <- n_lo
  repeat {
    hi <- min(max(2 * hi, hi + 1), n_max)
    if (f(hi) >= target) break
    lo <- hi
    if (hi >= n_max) return(NA_real_)
  }
  while (hi - lo > 1) {
    mid <- floor((lo + hi) / 2)
    if (f(mid) >= target) hi <- mid else lo <- mid
  }
  hi
}

# A priori search for a test whose power is bounded by a monotone function
# `bound(N)` (proven upper bound of the power at N): every N below the
# smallest N with bound(N) >= target fails; from there every N is enumerated
# upwards (power is not monotone in N for discrete tests, manual 4.3).
.ig_bounded_enumeration <- function(power_at, bound, target, n_lo, n_max, bound_text,
                                    max_enumeration = 2e6) {
  nb <- .ig_monotone_first_n(bound, target, n_lo, n_max)
  if (is.na(nb)) {
    return(list(status = "infeasible", reason = sprintf(paste(
      "the target power %s is not reached for any N <= %s: even %s stays below it"),
      format(target), format(n_max, scientific = FALSE), bound_text),
      diagnostics = list(bound = bound_text, n_upper_limit = n_max)))
  }
  N <- nb
  count <- 0
  seq_n <- numeric()
  seq_p <- numeric()
  while (N <= n_max) {
    pw <- power_at(N)
    count <- count + 1
    seq_n <- c(utils::tail(seq_n, 9L), N)
    seq_p <- c(utils::tail(seq_p, 9L), pw)
    if (pw >= target) {
      return(list(status = "ok", n = N, power = pw, diagnostics = list(
        search = "bounded enumeration",
        minimality = sprintf(paste("proven: every N below %s fails because %s is below the target there (a proven",
                                   "upper bound of the power, non-decreasing in N); every N from %s upwards was",
                                   "evaluated"), format(nb, scientific = FALSE), bound_text,
                             format(nb, scientific = FALSE)),
        monotonicity = "not assumed (discrete test; power is not monotone in N)",
        bound = bound_text, n_bound = nb, evaluations = count,
        power_sequence = list(n_total = seq_n, power = seq_p))))
    }
    if (count >= max_enumeration) break
    N <- N + 1
  }
  list(status = "infeasible", reason = sprintf(
    "the target power %s is not reached by any N from %s to %s (enumeration limit)", format(target),
    format(nb, scientific = FALSE), format(N, scientific = FALSE)),
    diagnostics = list(bound = bound_text, n_bound = nb, evaluations = count))
}

# Binomial forward model shared by the one-sample test, the sign test and the
# generic binomial test.
.ig_binom_forward <- function(N, p0, p1, alpha, tails, balancing, effect_note) {
  N <- as.double(N)
  bad <- !(is.finite(p1) && p1 > 0 && p1 < 1)
  if (bad) {
    return(list(power = NA_real_, beta = NA_real_, alpha_attained = NA_real_,
                notes = list(power = sprintf("the H1 proportion %s lies outside (0, 1)", format(p1)),
                             beta = "undefined: the H1 proportion lies outside (0, 1)",
                             alpha_attained = "undefined: the H1 proportion lies outside (0, 1)")))
  }
  t <- .ig_binom_test_power(N, p0, p1, alpha, tails, balancing)
  crit <- if (identical(tails, "one")) {
    if (t$side == "upper") t$region[2] else t$region[1]
  } else t$region
  list(
    power = t$power,
    beta = t$beta,
    alpha_attained = t$alpha_attained,
    critical = crit,
    n_total = N,
    outputs = list(
      lower_critical_n = if (identical(tails, "one")) crit else t$region[1],
      upper_critical_n = if (identical(tails, "one")) crit else t$region[2],
      alpha_attained_lower = t$alpha_lower,
      alpha_attained_upper = t$alpha_upper,
      proportion_h0 = p0,
      proportion_h1 = p1
    ),
    statistic_distribution = list(
      h0 = list(family = "binomial", size = N, prob = p0),
      h1 = list(family = "binomial", size = N, prob = p1),
      rejection_region = list(lower = t$region[1], upper = t$region[2]),
      rejection_side = t$side
    ),
    assumptions = c(
      sprintf("X ~ Bin(N, p): H0 proportion %s, H1 proportion %s (%s)", format(p0), format(p1), effect_note),
      paste("reject H0 when X <= lower critical N or X >= upper critical N (", t$rule, "); the attained size",
            "never exceeds alpha; empty tails are reported as lower = -1 and upper = N + 1"),
      if (identical(tails, "one")) "one-sided: both critical-N outputs show the single critical count, as G*Power displays it"
    )
  )
}

# A priori search shared by the binomial tests.
.ig_binom_a_priori_search <- function(p0_of, p1_of) {
  function(p, n_lo, n_max) {
    p0 <- p0_of(p)
    p1 <- p1_of(p)
    target <- p$power
    if (abs(p1 - p0) <= .ig_prop_eps) {
      return(list(status = "infeasible", reason = sprintf(paste(
        "the H1 proportion %s equals the H0 proportion %s up to eps = 1e-6: the power equals the attained size",
        "(<= alpha) for every N, so the target power %s is never reached (G*Power requires |effect| > 1e-6 in",
        "a priori analyses)"), format(p1), format(p0), format(target)),
        diagnostics = list(rule = "no effect")))
    }
    fwd <- function(N) .ig_binom_test_power(N, p0, p1, p$alpha, p$tails, p$alpha_balancing)$power
    .ig_bounded_enumeration(
      power_at = fwd,
      bound = function(N) .ig_binom_power_bound(N, p0, p1, p$alpha, p$tails, p$alpha_balancing),
      target = target, n_lo = max(n_lo, 1), n_max = n_max,
      bound_text = paste("the most powerful level-alpha power (randomised one-sided Neyman-Pearson test; two tails",
                         "with alpha/2 per side: MP at alpha/2 plus the exact power of the opposite tail)"))
  }
}

# One-sample binomial test ----------------------------------------------------------------
#
# G*Power: Tests > Proportions > One group: Difference from constant
# ("Exact - Proportion: Difference from constant (binomial test, one sample
# case)", manual ch. 4). H0: pi = pi0, H1: pi = pi0 + g.

.ig_method_exact_proportion_one_sample_binomial <- function() {
  id <- "exact_proportion_one_sample_binomial"
  ig_method(
    id = id,
    family = "exact",
    label = list(en = "Proportion: Difference from constant (binomial test, one sample case)",
                 de = "Anteil: Abweichung von einer Konstanten (Binomialtest, eine Stichprobe)"),
    params = ig_params(
      ig_param("g", role = "effect",
               domain = ig_domain_real(-1, 1, lower_inclusive = FALSE, upper_inclusive = FALSE),
               label = list(en = "Effect size g", de = "Effektst\u00e4rke g"),
               description = paste("g = pi - pi0, the H1 proportion minus the constant proportion; pi0 + g must lie in",
                                   "(1e-6, 1 - 1e-6) (manual 4.1)."),
               aliases = "effect_size_g"),
      ig_param("p0", role = "nuisance",
               domain = ig_domain_real(.ig_prop_eps, 1 - .ig_prop_eps, lower_inclusive = FALSE, upper_inclusive = FALSE),
               default = 0.5,
               label = list(en = "Constant proportion", de = "Konstanter Anteil"),
               description = "H0 proportion pi0. Default 0.5 = G*Power's initial field value (certified inventory).",
               aliases = "constant_proportion"),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      ig_param_n_total(lower = 1)
    ),
    options = ig_options(.ig_prop_tails_option(), .ig_prop_balancing_option()),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = function(p) {
      # The static domain of g is (-1, 1); the admissible range of pi0 + g
      # depends on pi0 and is enforced by `constraints` for every
      # specification. Only the sensitivity scan evaluates g beyond it: there
      # the H1 proportion is held at the nearest admissible value, a monotone
      # continuation of the power curve that only brackets the root (the
      # reported effect always lies inside the admissible range).
      p1 <- min(max(p$p0 + p$g, .ig_prop_eps), 1 - .ig_prop_eps)
      .ig_binom_forward(p$n_total, p$p0, p1, p$alpha, p$tails, p$alpha_balancing,
                        "H1 proportion = constant proportion + g")
    },
    design = NULL,
    discrete = TRUE,
    effect_direction = "increasing",
    engine = list(package = "stats", fun = "pbinom/qbinom (binomial tails; inferguidR rejection regions)",
                  version = NA_character_),
    checker = "validation/independent/g8/check_binomial.R",
    reference = list(
      gpower_path = "Tests > Proportions > One group: Difference From Constant",
      gpower_label = "Exact - Proportion: Difference from constant (binomial test, one sample case)",
      manual_pages = "ch. 4 (PDF pp. 11-13); options 4.2; a priori rule 4.3",
      allocation_rule = "single sample: N is the design",
      method_notes = "specs/method-notes/exact_proportion_one_sample_binomial.md"
    ),
    status = "independently_checked",
    constraints = function(p) {
      out <- .ig_prop_power_alpha_constraint(p)
      if (!is.null(p$g) && !is.null(p$p0)) {
        p1 <- p$p0 + p$g
        if (!(p1 > .ig_prop_eps && p1 < 1 - .ig_prop_eps)) {
          out <- c(out, sprintf(paste("the effect size g = %s with constant proportion %s gives an H1 proportion %s",
                                      "outside the admissible range (1e-6, 1 - 1e-6)"),
                                .ig_format_double(p$g), .ig_format_double(p$p0), .ig_format_double(p1)))
        }
      }
      out
    },
    a_priori_search = .ig_binom_a_priori_search(function(p) p$p0, function(p) p$p0 + p$g)
  )
}

# Sign test ---------------------------------------------------------------------------------
#
# G*Power: Tests > Proportions > Sign test (manual ch. 8): the binomial test
# with pi0 = 0.5 and effect g = pi - 0.5.

.ig_method_exact_proportion_sign_test <- function() {
  id <- "exact_proportion_sign_test"
  ig_method(
    id = id,
    family = "exact",
    label = list(en = "Proportion: Sign test (binomial test)", de = "Anteil: Vorzeichentest (Binomialtest)"),
    params = ig_params(
      ig_param("g", role = "effect",
               domain = ig_domain_real(-0.5 + .ig_prop_eps, 0.5 - .ig_prop_eps, lower_inclusive = FALSE,
                                       upper_inclusive = FALSE),
               label = list(en = "Effect size g", de = "Effektst\u00e4rke g"),
               description = "g = pi - 0.5 (manual 8.1); admissible range (-0.499999, 0.499999) as in G*Power.",
               aliases = "effect_size_g"),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      ig_param_n_total(lower = 1)
    ),
    options = ig_options(.ig_prop_tails_option(), .ig_prop_balancing_option()),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = function(p) {
      .ig_binom_forward(p$n_total, 0.5, 0.5 + p$g, p$alpha, p$tails, p$alpha_balancing,
                        "sign test: H0 proportion 0.5, H1 proportion 0.5 + g")
    },
    design = NULL,
    discrete = TRUE,
    effect_direction = "increasing",
    engine = list(package = "stats", fun = "pbinom/qbinom (binomial tails; inferguidR rejection regions)",
                  version = NA_character_),
    checker = "validation/independent/g8/check_binomial.R",
    reference = list(
      gpower_path = "Tests > Proportions > Sign test",
      gpower_label = "Exact - Proportion: Sign test (binomial test)",
      manual_pages = "ch. 8 (PDF p. 22); options and conventions as ch. 4",
      allocation_rule = "single sample: N is the design",
      method_notes = "specs/method-notes/exact_proportion_sign_test.md"
    ),
    status = "independently_checked",
    constraints = .ig_prop_power_alpha_constraint,
    a_priori_search = .ig_binom_a_priori_search(function(p) 0.5, function(p) 0.5 + p$g)
  )
}

# Generic binomial test -----------------------------------------------------------------------
#
# G*Power: Tests > Generic binomial test (manual ch. 9, a stub): the
# binomial test of H0: pi = p1 against H1: pi = p2, both proportions entered
# directly ("Proportion p1" = H0, "Proportion p2" = H1; certified inventory).

.ig_method_exact_generic_binomial <- function() {
  id <- "exact_generic_binomial"
  rng <- ig_domain_real(.ig_prop_eps, 1 - .ig_prop_eps, lower_inclusive = FALSE, upper_inclusive = FALSE)
  ig_method(
    id = id,
    family = "exact",
    label = list(en = "Generic binomial test", de = "Generischer Binomialtest"),
    params = ig_params(
      ig_param("p2", role = "effect", domain = rng,
               label = list(en = "Proportion p2", de = "Anteil p2"),
               description = "H1 proportion (G*Power's 'Proportion p2'); admissible range (1e-6, 1 - 1e-6).",
               aliases = "proportion_p2"),
      ig_param("p1", role = "nuisance", domain = rng, default = 0.5,
               label = list(en = "Proportion p1", de = "Anteil p1"),
               description = paste("H0 proportion (G*Power's 'Proportion p1'). Default 0.5 = G*Power's initial field",
                                   "value (certified inventory)."),
               aliases = "proportion_p1"),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      ig_param_n_total(lower = 1)
    ),
    options = ig_options(.ig_prop_tails_option(), .ig_prop_balancing_option()),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = function(p) {
      .ig_binom_forward(p$n_total, p$p1, p$p2, p$alpha, p$tails, p$alpha_balancing,
                        "generic binomial test: H0 proportion p1, H1 proportion p2")
    },
    design = NULL,
    discrete = TRUE,
    effect_direction = "increasing",
    engine = list(package = "stats", fun = "pbinom/qbinom (binomial tails; inferguidR rejection regions)",
                  version = NA_character_),
    checker = "validation/independent/g8/check_binomial.R",
    reference = list(
      gpower_path = "Tests > Generic binomial test",
      gpower_label = "Exact - Generic binomial test",
      manual_pages = "ch. 9 (PDF p. 23; the chapter documents the options only)",
      allocation_rule = "single sample: N is the design",
      method_notes = "specs/method-notes/exact_generic_binomial.md"
    ),
    status = "independently_checked",
    constraints = .ig_prop_power_alpha_constraint,
    a_priori_search = .ig_binom_a_priori_search(function(p) p$p1, function(p) p$p2)
  )
}

# Two independent proportions: z test -----------------------------------------------------
#
# G*Power: Tests > Proportions > Two independent groups: inequality, z-test
# ("z tests - Proportions: Difference between two independent proportions";
# no manual chapter, behaviour from the certified inventory and the
# captures; specs/method-notes/z_proportions_two_independent.md).
#
# * Plain test: D = p2 - p1; H0 SD s0 = sqrt(pbar (1 - pbar) (1/n1 + 1/n2))
#   with the pooled pbar = (n1 p1 + n2 p2) / (n1 + n2); H1 SD
#   s1 = sqrt(p1 (1 - p1) / n1 + p2 (1 - p2) / n2); power
#   Phi((|D| - z s0) / s1) (+ the opposite tail for two tails).
# * Arcsine transform: h = 2 asin(sqrt(p2)) - 2 asin(sqrt(p1)) (Cohen),
#   se = sqrt(1/n1 + 1/n2), power Phi(|h| / se - z) (+ opposite tail).
# * Continuity correction (Fleiss, Tytun & Ury 1980): G*Power evaluates the
#   power of the uncorrected statistic at the effective group sizes of the
#   FTU relation n' = (n/4) (1 + sqrt(1 + 2 (r + 1) / (n r |D|)))^2, i.e.
#   n1_eff = (n1 - (r + 1) / (2 r |D|))^2 / n1, n2_eff = r n1_eff
#   (r = n2 / n1), each rounded up to an integer. For n1 <= (r + 1) / (2 r |D|)
#   the FTU relation has no positive solution; G*Power then uses the spurious
#   root of the squared equation (power above the uncorrected one: a
#   reference defect), inferguidR uses the smallest admissible effective size
#   (1 per group).
# * "Plot normal distribution (instead of z)" reports the critical value on
#   the difference scale (z * s0; arcsine: z * se) instead of the z scale.
# * Critical value sign: the side of p2 - p1, as G*Power displays it.

.ig_z2p_effective_n <- function(n1, n2, D) {
  if (!(D > 0)) return(c(n1, n2))
  r <- n2 / n1
  a4 <- (r + 1) / (2 * r * D)
  if (n1 <= a4) return(c(1, max(1, ceiling(r - 1e-12))))
  e1 <- (n1 - a4)^2 / n1
  c(max(1, ceiling(e1 - 1e-9 * e1)), max(1, ceiling(r * e1 - 1e-9 * r * e1)))
}

.ig_z2p_core <- function(p1, p2, n1, n2, alpha, tails, arcsine) {
  zc <- if (identical(tails, "two")) stats::qnorm(alpha / 2, lower.tail = FALSE) else stats::qnorm(alpha, lower.tail = FALSE)
  if (arcsine) {
    eff <- 2 * asin(sqrt(p2)) - 2 * asin(sqrt(p1))
    s0 <- sqrt(1 / n1 + 1 / n2)
    s1 <- s0
  } else {
    eff <- p2 - p1
    pbar <- (n1 * p1 + n2 * p2) / (n1 + n2)
    s0 <- sqrt(pbar * (1 - pbar) * (1 / n1 + 1 / n2))
    s1 <- sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
  }
  a <- (abs(eff) - zc * s0) / s1
  if (identical(tails, "two")) {
    b <- (-abs(eff) - zc * s0) / s1
    power <- stats::pnorm(a) + stats::pnorm(b)
    beta <- max(stats::pnorm(-a) - stats::pnorm(b), 0)
  } else {
    power <- stats::pnorm(a)
    beta <- stats::pnorm(a, lower.tail = FALSE)
  }
  list(power = min(power, 1), beta = beta, zc = zc, eff = eff, s0 = s0, s1 = s1)
}

.ig_z2p_p2 <- function(p) {
  if (identical(p$input, "cohens_h")) {
    th <- asin(sqrt(p$p1)) + p$h / 2
    # Outside (0, pi/2) sin^2 would fold back; the sensitivity scan only needs
    # a monotone continuation there (constraints reject such h otherwise).
    th <- min(max(th, asin(sqrt(.ig_prop_eps))), asin(sqrt(1 - .ig_prop_eps)))
    return(sin(th)^2)
  }
  p$p2
}

.ig_power_z_proportions_two_independent <- function(p) {
  n1 <- as.double(p$n1)
  n2 <- as.double(p$n2)
  p1 <- p$p1
  p2 <- .ig_z2p_p2(p)
  cc <- identical(p$continuity_correction, "on")
  as <- identical(p$arcsin_transform, "on")
  D <- abs(p2 - p1)
  ne <- if (cc) .ig_z2p_effective_n(n1, n2, D) else c(n1, n2)
  below_ftu <- cc && D > 0 && n1 <= (n2 / n1 + 1) / (2 * (n2 / n1) * D)
  core <- .ig_z2p_core(p1, p2, ne[1], ne[2], p$alpha, p$tails, as)
  sgn <- if (p2 >= p1) 1 else -1
  crit_z <- sgn * core$zc
  crit_d <- sgn * core$zc * core$s0
  shown <- if (identical(p$distribution_plot, "on")) crit_d else crit_z
  critical <- if (identical(p$tails, "two")) c(-abs(shown), abs(shown)) else shown
  list(
    power = core$power,
    beta = core$beta,
    alpha_attained = p$alpha,
    critical = critical,
    n_total = n1 + n2,
    group_sizes = c(n1 = n1, n2 = n2),
    outputs = list(
      critical_effect_side = shown,
      critical_z = crit_z,
      critical_difference_scale = crit_d,
      proportion_p1 = p1,
      proportion_p2 = p2,
      cohens_h = 2 * asin(sqrt(p2)) - 2 * asin(sqrt(p1)),
      effective_n1 = ne[1],
      effective_n2 = ne[2]
    ),
    statistic_distribution = list(
      scale = if (as) "arcsine difference h" else "proportion difference p2 - p1",
      h0 = list(family = "normal", mean = 0, sd = core$s0),
      h1 = list(family = "normal", mean = core$eff, sd = core$s1),
      rejection_side = if (identical(p$tails, "two")) "both" else if (sgn > 0) "upper" else "lower"
    ),
    assumptions = c(
      if (as) "arcsine-transformed difference h = 2 asin(sqrt(p2)) - 2 asin(sqrt(p1)) with SE sqrt(1/n1 + 1/n2) (Cohen 1988)"
      else "normal approximation: pooled-proportion SD under H0, unpooled SD under H1",
      if (cc) paste("continuity correction (Fleiss, Tytun & Ury 1980) as G*Power implements it: the uncorrected power at the",
                    "effective group sizes of the FTU relation, rounded up; below the FTU minimum the effective size is 1",
                    "per group (G*Power uses a spurious root there; method notes)"),
      if (below_ftu) paste("group 1 is not larger than the Fleiss-Tytun-Ury minimum (r + 1) / (2 r |p2 - p1|): the corrected",
                           "test is evaluated at the smallest effective group sizes and the approximation is not reliable here"),
      if (identical(p$tails, "two")) "two-sided: alpha/2 in each tail; the power counts both rejection tails"
      else "one-sided: the rejection region lies on the side of p2 - p1",
      if (identical(p$distribution_plot, "on")) "critical value reported on the difference scale ('Plot normal distribution')"
    )
  )
}

# A priori allocation (G*Power's split rule, also found for the t test):
# n1 = floor(N / (1 + r) + 1/2), n2 = N - n1; with r = 1 only even N (equal
# groups). Designs with an empty group are not admissible. The total N is
# chosen by .ig_two_group_a_priori_search() below.
.ig_design_z_proportions_two_independent <- function(p, n_total, continuous = FALSE) {
  r <- p$allocation_ratio
  if (continuous) return(c(n1 = n_total / (1 + r), n2 = r * n_total / (1 + r)))
  if (r == 1) {
    if (n_total %% 2 != 0) return(NULL)
    return(c(n1 = n_total / 2, n2 = n_total / 2))
  }
  n1 <- floor(n_total / (1 + r) + 0.5)
  n2 <- n_total - n1
  if (n1 < 1 || n2 < 1) return(NULL)
  c(n1 = n1, n2 = n2)
}

# A priori search of the z test of two proportions. `at(n1, n2)` is the power
# at (possibly fractional) group sizes.
#
# * Allocation ratio 1: equal groups, the smallest n per group reaching the
#   target (generic integer search over even N).
# * Other ratios: G*Power 3.1.9.7 takes the smallest N >= 2 whose power with
#   the fractional group sizes N/(1+r), rN/(1+r) reaches the target and splits
#   it as n1 = floor(N/(1+r) + 1/2), n2 = N - n1 (captures
#   SUP-G8-ZPROP2-APR-A01..A08 and the campaign a priori records). That
#   design is reported when both groups are non-empty and its own power
#   reaches the target; otherwise the smallest larger N whose rounded split
#   does (G*Power reports the design below the target without a message).
#   G*Power's design is always kept in diagnostics$gpower_allocation_rule.
.ig_two_group_a_priori_search <- function(p, n_lo, n_max, at, design, max_adjust = 1e5) {
  r <- as.double(p$allocation_ratio)
  target <- p$power
  int_power <- function(N) {
    gs <- design(p, N)
    if (is.null(gs)) return(NA_real_)
    at(gs[["n1"]], gs[["n2"]])
  }
  if (r == 1) {
    return(.ig_integer_search(power_at = int_power, admissible = function(N) !is.null(design(p, N)),
                              n_lo = n_lo, n_max = n_max, target = target, discrete = FALSE))
  }
  frac_power <- function(N) at(N / (1 + r), r * N / (1 + r))
  lo <- max(2, n_lo)
  if (abs(frac_power(n_max) - frac_power(lo)) <= 1e-14 && frac_power(lo) < target) {
    return(list(status = "infeasible", reason = sprintf(
      "the power does not change with N (%s for N = %s and N = %s), so the target power %s is never reached",
      format(signif(frac_power(lo), 7)), format(lo), format(n_max, scientific = FALSE), format(target)),
      diagnostics = list(allocation_ratio = r)))
  }
  n_star <- .ig_monotone_first_n(frac_power, target, lo, n_max)
  if (is.na(n_star)) {
    return(list(status = "infeasible", reason = sprintf(paste(
      "the target power %s is not reached with fractional group sizes for any N <= %s (power there: %s)"),
      format(target), format(n_max, scientific = FALSE), format(signif(frac_power(n_max), 7))),
      diagnostics = list(allocation_ratio = r)))
  }
  split <- function(N) {
    n1 <- floor(N / (1 + r) + 0.5)
    c(n1 = n1, n2 = N - n1)
  }
  int_ok <- function(N) {
    v <- int_power(N)
    !is.na(v) && v >= target
  }
  gs_star <- split(n_star)
  empty <- gs_star[["n1"]] < 1 || gs_star[["n2"]] < 1
  power_star <- if (empty) NA_real_ else at(gs_star[["n1"]], gs_star[["n2"]])
  gp_rule <- list(
    n_total = n_star, n1 = unname(gs_star[["n1"]]), n2 = unname(gs_star[["n2"]]),
    fractional_power = frac_power(n_star),
    fractional_power_predecessor = if (n_star - 1 >= lo) frac_power(n_star - 1) else NA_real_,
    integer_design_power = power_star, empty_group = empty,
    reaches_target = !empty && power_star >= target)
  diag <- list(
    search = paste("G*Power 3.1.9.7 rule for unequal allocation: smallest N >= 2 whose fractional-allocation power",
                   "reaches the target (doubling and bisection; that power is non-decreasing in N), split as",
                   "n1 = floor(N / (1 + r) + 1/2), n2 = N - n1"),
    allocation_ratio = r, n_lower = lo, n_upper_limit = n_max)
  warnings <- character()
  if (gp_rule$reaches_target) {
    n <- n_star
    diag$minimality <- sprintf(paste(
      "proven for the fractional criterion: N = %s is the smallest N whose fractional-allocation power reaches",
      "the target%s"), format(n_star, scientific = FALSE),
      if (is.na(gp_rule$fractional_power_predecessor)) "" else
        sprintf(" (fractional power at N - 1: %s)", format(signif(gp_rule$fractional_power_predecessor, 10))))
  } else {
    n <- NA_real_
    N <- n_star
    for (k in seq_len(max_adjust)) {
      N <- N + 1
      if (N > n_max) break
      if (int_ok(N)) {
        n <- N
        break
      }
    }
    why <- if (empty) "leaves a group empty" else sprintf(
      "has power %s at its integer design, below the target %s (G*Power 3.1.9.7 reports it without a message)",
      format(signif(power_star, 7)), format(target))
    if (is.na(n)) {
      return(list(status = "infeasible", reason = sprintf(paste(
        "G*Power's allocation rule gives N = %s (n1 = %s, n2 = %s), which %s, and none of the next %s totals has a",
        "rounded split that reaches the target power"), format(n_star, scientific = FALSE), format(gs_star[["n1"]]),
        format(gs_star[["n2"]]), why, format(max_adjust, scientific = FALSE)),
        diagnostics = c(diag, list(gpower_allocation_rule = gp_rule))))
    }
    warnings <- sprintf(paste(
      "G*Power's allocation rule for N2/N1 = %s gives N = %s (n1 = %s, n2 = %s), which %s; inferguidR reports the",
      "smallest larger N whose rounded split reaches the target (see diagnostics$gpower_allocation_rule)"),
      .ig_format_double(r), format(n_star, scientific = FALSE), format(gs_star[["n1"]]), format(gs_star[["n2"]]), why)
    diag$minimality <- sprintf(paste(
      "N = %s is the smallest N above G*Power's N = %s whose rounded split has both groups non-empty and reaches",
      "the target (every N in between was evaluated)"), format(n, scientific = FALSE),
      format(n_star, scientific = FALSE))
    diag$adjusted_from_gpower_rule <- TRUE
    diag$adjustment_note <- warnings
  }
  diag$gpower_allocation_rule <- gp_rule
  list(status = "ok", n = n, power = int_power(n), diagnostics = diag, warnings = warnings)
}

.ig_z2p_a_priori_search <- function(p, n_lo, n_max) {
  at <- function(n1, n2) {
    q <- p
    q$n1 <- n1
    q$n2 <- n2
    .ig_power_z_proportions_two_independent(q)$power
  }
  .ig_two_group_a_priori_search(p, n_lo, n_max, at, .ig_design_z_proportions_two_independent)
}

.ig_method_z_proportions_two_independent <- function() {
  rng <- ig_domain_real(.ig_prop_eps, 1 - .ig_prop_eps, lower_inclusive = FALSE, upper_inclusive = FALSE)
  onoff <- function(name, label, lab_en) {
    ig_option(name, values = c("on", "off"), default = "off", labels = c(on = "checked", off = "unchecked"),
              label = list(en = lab_en))
  }
  ig_method(
    id = "z_proportions_two_independent",
    family = "z",
    label = list(en = "Proportions: Difference between two independent proportions",
                 de = "Anteile: Differenz zweier unabh\u00e4ngiger Anteile"),
    params = ig_params(
      ig_param("p2", role = "effect", domain = rng,
               label = list(en = "Proportion p2", de = "Anteil p2"),
               description = "Proportion of group 2 (H1); the effect is p2 - p1.",
               aliases = "proportion_p2",
               applies_when = list(option = "input", "in" = "proportions")),
      ig_param("h", role = "effect",
               domain = ig_domain_real(-pi, pi, lower_inclusive = FALSE, upper_inclusive = FALSE),
               label = list(en = "Effect size h", de = "Effektst\u00e4rke h"),
               description = "Cohen's h = 2 asin(sqrt(p2)) - 2 asin(sqrt(p1)) (option input = cohens_h).",
               aliases = "effect_size_h",
               applies_when = list(option = "input", "in" = "cohens_h")),
      ig_param("p1", role = "nuisance", domain = rng, default = 0.5,
               label = list(en = "Proportion p1", de = "Anteil p1"),
               description = "Proportion of group 1. Default 0.5 = G*Power's initial field value (certified inventory).",
               aliases = "proportion_p1"),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      ig_param("n_total", role = "design", domain = ig_domain_integer(lower = 2), unit = "subjects",
               label = list(en = "Total sample size", de = "Gesamtstichprobenumfang"),
               description = "Solved in a priori mode (n1 + n2).", modes = "a_priori"),
      ig_param("n1", role = "design", domain = ig_domain_integer(lower = 1), unit = "subjects",
               label = list(en = "Sample size group 1", de = "Stichprobenumfang Gruppe 1"),
               aliases = "sample_size_group_1", modes = c("post_hoc", "sensitivity", "criterion", "compromise")),
      ig_param("n2", role = "design", domain = ig_domain_integer(lower = 1), unit = "subjects",
               label = list(en = "Sample size group 2", de = "Stichprobenumfang Gruppe 2"),
               aliases = "sample_size_group_2", modes = c("post_hoc", "sensitivity", "criterion", "compromise")),
      ig_param("allocation_ratio", role = "design",
               domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE), default = 1,
               label = list(en = "Allocation ratio N2/N1", de = "Allokationsverh\u00e4ltnis N2/N1"),
               description = "Group-size ratio n2/n1 of the a priori design; default 1 = G*Power's initial field value.",
               aliases = "allocation_ratio_n2_n1", modes = "a_priori")
    ),
    options = ig_options(
      .ig_prop_tails_option(),
      onoff("continuity_correction", "Use continuity correction", "Use continuity correction"),
      onoff("arcsin_transform", "Use arcsin transform", "Use arcsin transform"),
      ig_option("input", values = c("proportions", "cohens_h"), default = "proportions",
                labels = c(proportions = "Use two proportions p1 / p2", cohens_h = "Use Cohen's effect size index h"),
                label = list(en = "Input")),
      onoff("distribution_plot", "Plot normal distribution (instead of z)", "Plot normal distribution (instead of z)")
    ),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = .ig_power_z_proportions_two_independent,
    design = .ig_design_z_proportions_two_independent,
    discrete = FALSE,
    effect_direction = "increasing",
    engine = list(package = "stats", fun = "pnorm/qnorm (normal approximation; FTU continuity correction)",
                  version = NA_character_),
    checker = "validation/independent/g8/check_z_proportions.R",
    reference = list(
      gpower_path = "Tests > Proportions > Two independent groups: inequality, z-test",
      gpower_label = "z tests - Proportions: Difference between two independent proportions",
      manual_pages = "no manual chapter (application-only procedure)",
      allocation_rule = paste("allocation ratio 1: n1 = n2 = N/2 (N even), the smallest design reaching the target;",
                              "otherwise G*Power's rule: the smallest N whose power with fractional group sizes",
                              "N/(1+r), rN/(1+r) reaches the target, split as n1 = floor(N/(1+r) + 1/2), n2 = N - n1;",
                              "where that integer design misses the target, the smallest larger N whose split reaches it"),
      method_notes = "specs/method-notes/z_proportions_two_independent.md"
    ),
    status = "independently_checked",
    a_priori_search = .ig_z2p_a_priori_search,
    targets = list(sensitivity = c("p2", "h")),
    constraints = function(p) {
      out <- .ig_prop_power_alpha_constraint(p)
      if (!is.null(p$h) && !is.null(p$p1)) {
        th <- asin(sqrt(p$p1)) + p$h / 2
        if (!(th > asin(sqrt(.ig_prop_eps)) && th < asin(sqrt(1 - .ig_prop_eps)))) {
          out <- c(out, sprintf("Cohen's h = %s with p1 = %s implies a proportion p2 outside (1e-6, 1 - 1e-6)",
                                .ig_format_double(p$h), .ig_format_double(p$p1)))
        }
      }
      out
    }
  )
}

# Vectorised binomial rejection regions ------------------------------------------------------
#
# The same rules as .ig_binom_region() for a vector of sample sizes n (scalar
# p0, alpha and H1 proportion p1); used by McNemar's test, where one region is
# needed per number of discordant pairs. validation/independent/g8 checks
# that both implementations agree.

.ig_binom_lower_crit_vec <- function(n, p, a) {
  if (!(a > 0)) return(rep(-1, length(n)))
  if (a >= 1) return(as.double(n))
  k <- stats::qbinom(a, n, p)
  repeat {
    bad <- k >= 0 & .ig_binom_cdf(k, n, p) > a
    if (!any(bad)) break
    k[bad] <- k[bad] - 1
  }
  repeat {
    inc <- k + 1 <= n & .ig_binom_cdf(k + 1, n, p) <= a
    if (!any(inc)) break
    k[inc] <- k[inc] + 1
  }
  k
}

.ig_binom_upper_crit_vec <- function(n, p, a) {
  if (!(a > 0)) return(as.double(n) + 1)
  if (a >= 1) return(rep(0, length(n)))
  k <- stats::qbinom(a, n, p, lower.tail = FALSE) + 1
  repeat {
    bad <- k <= n & .ig_binom_sf(k, n, p) > a
    if (!any(bad)) break
    k[bad] <- k[bad] + 1
  }
  repeat {
    dec <- k - 1 >= 0 & .ig_binom_sf(k - 1, n, p) <= a
    if (!any(dec)) break
    k[dec] <- k[dec] - 1
  }
  k
}

.ig_binom_region_vec <- function(n, p0, alpha, tails, balancing, p1) {
  n <- as.double(n)
  up_dir <- p1 >= p0
  if (identical(tails, "one")) {
    if (up_dir) return(list(lower = rep(-1, length(n)), upper = .ig_binom_upper_crit_vec(n, p0, alpha)))
    return(list(lower = .ig_binom_lower_crit_vec(n, p0, alpha), upper = n + 1))
  }
  balancing <- balancing %||% "half_each_side"
  kl <- .ig_binom_lower_crit_vec(n, p0, alpha / 2)
  ku <- .ig_binom_upper_crit_vec(n, p0, alpha / 2)
  if (identical(balancing, "half_each_side")) return(list(lower = kl, upper = ku))
  if (identical(balancing, "minor_tail_half_rest_major")) {
    a_up <- vapply(alpha - .ig_binom_cdf(kl, n, p0), function(x) x, 0)
    a_lo <- alpha - .ig_binom_sf(ku, n, p0)
    Au <- vapply(seq_along(n), function(j) .ig_binom_upper_crit(n[j], p0, a_up[j]), 0)
    Bl <- vapply(seq_along(n), function(j) .ig_binom_lower_crit(n[j], p0, a_lo[j]), 0)
    acc_a <- vapply(seq_along(n), function(j) .ig_binom_region_accept(c(kl[j], Au[j]), n[j], p1), 0)
    acc_b <- vapply(seq_along(n), function(j) .ig_binom_region_accept(c(Bl[j], ku[j]), n[j], p1), 0)
    useA <- acc_a <= acc_b
    return(list(lower = ifelse(useA, kl, Bl), upper = ifelse(useA, Au, ku)))
  }
  lo <- kl
  up <- ku
  repeat {
    tot <- .ig_binom_cdf(lo, n, p0) + .ig_binom_sf(up, n, p0)
    il <- ifelse(lo + 1 < up, stats::dbinom(lo + 1, n, p0), Inf)
    iu <- ifelse(up - 1 > lo, stats::dbinom(up - 1, n, p0), Inf)
    step <- pmin(il, iu)
    go <- is.finite(step) & tot + step <= alpha
    if (!any(go)) break
    tie <- is.finite(il) & is.finite(iu) & abs(il - iu) <= 1e-10 * pmax(il, iu)
    widen_lo <- go & ((!tie & il < iu) | (tie & !up_dir))
    widen_up <- go & !widen_lo
    lo[widen_lo] <- lo[widen_lo] + 1
    up[widen_up] <- up[widen_up] - 1
  }
  list(lower = lo, upper = up)
}

# McNemar test ------------------------------------------------------------------------------
#
# G*Power: Tests > Proportions > Two dependent groups: McNemar ("Exact -
# Proportions: Inequality, two dependent groups (McNemar)", manual ch. 5).
# Paired binary responses; discordant cells pi12 (yes/no) and pi21 (no/yes),
# pi_D = pi12 + pi21, odds ratio OR = pi12 / pi21. Conditional on the number of
# discordant pairs N_D = i the test is the binomial test of theta = 1/2 for
# the count of (yes/no) pairs, with theta = OR / (1 + OR) under H1.
#
# * Exact (unconditional) power: sum over i of P(N_D = i) Pow(i), with
#   N_D ~ Bin(N, pi_D) (the exponent order printed in manual 5.5 is swapped;
#   registry note) and the conditional binomial test of the alpha-balancing
#   option; attained alpha = sum of P(N_D = i) times the conditional attained
#   size. Used when N < the exact-computation limit (default 1e6).
# * Fast approximation (manual 5.2.2: "assume number of discordant pairs as
#   constant"): one binomial test with N_D = ceiling(N pi_D) (probe grid of 60
#   post hoc captures; method notes).

.ig_mcnemar_terms <- function(N, pd) {
  N <- as.double(N)
  if (N <= 2000) return(list(i = 0:N, w = stats::dbinom(0:N, N, pd), omitted = 0))
  lo <- stats::qbinom(1e-18, N, pd)
  hi <- stats::qbinom(1e-18, N, pd, lower.tail = FALSE)
  i <- seq(lo, hi)
  omitted <- (if (lo > 0) stats::pbinom(lo - 1, N, pd) else 0) + stats::pbinom(hi, N, pd, lower.tail = FALSE)
  list(i = i, w = stats::dbinom(i, N, pd), omitted = omitted)
}

.ig_mcnemar_conditional <- function(i, alpha, tails, balancing, theta) {
  pos <- i > 0
  out <- list(power = numeric(length(i)), alpha = numeric(length(i)), accept = rep(1, length(i)))
  if (!any(pos)) return(out)
  ii <- as.double(i[pos])
  r <- .ig_binom_region_vec(ii, 0.5, alpha, tails, balancing, theta)
  out$power[pos] <- .ig_binom_cdf(r$lower, ii, theta) + .ig_binom_sf(r$upper, ii, theta)
  out$alpha[pos] <- .ig_binom_cdf(r$lower, ii, 0.5) + .ig_binom_sf(r$upper, ii, 0.5)
  # P(lower < X < upper) under H1 from the tail representation with small terms.
  a <- pmax(r$lower + 1, 0)
  b <- pmin(r$upper - 1, ii)
  low_form <- .ig_binom_cdf(b, ii, theta) - .ig_binom_cdf(a - 1, ii, theta)
  up_form <- .ig_binom_sf(a, ii, theta) - .ig_binom_sf(b + 1, ii, theta)
  acc <- ifelse(.ig_binom_cdf(b, ii, theta) <= 0.5, low_form, up_form)
  acc[b < a] <- 0
  out$accept[pos] <- pmin(pmax(acc, 0), 1)
  out
}

.ig_mcnemar_is_exact <- function(p) {
  identical(p$computation, "exact_unconditional") && as.double(p$n_total) < p$exact_n_limit
}

# N_D = ceiling(N pi_D); a product that is an integer up to representation
# error (for example 30 * 0.1) is that integer.
.ig_mcnemar_discordant <- function(N, pd) {
  x <- as.double(N) * pd
  max(1, ceiling(x - 8 * .Machine$double.eps * x))
}

.ig_power_exact_proportions_mcnemar <- function(p) {
  N <- as.double(p$n_total)
  or <- p$or
  pd <- p$pd
  theta <- or / (1 + or)
  p12 <- pd * or / (1 + or)
  p21 <- pd / (1 + or)
  common_out <- list(proportion_p12 = p12, proportion_p21 = p21, conditional_h1_proportion = theta)
  dist <- list(h0 = list(family = "binomial (conditional on the discordant pairs)", prob = 0.5),
               h1 = list(family = "binomial (conditional on the discordant pairs)", prob = theta))
  if (.ig_mcnemar_is_exact(p)) {
    tm <- .ig_mcnemar_terms(N, pd)
    cnd <- .ig_mcnemar_conditional(tm$i, p$alpha, p$tails, p$alpha_balancing, theta)
    power <- sum(tm$w * cnd$power)
    beta <- sum(tm$w * cnd$accept)
    att <- sum(tm$w * cnd$alpha)
    dist$discordant_pairs <- list(family = "binomial", size = N, prob = pd)
    return(list(
      power = min(power, 1), beta = min(beta, 1), alpha_attained = min(att, 1),
      n_total = N,
      outputs = c(common_out, list(computation = 1, truncated_mass = tm$omitted)),
      statistic_distribution = dist,
      assumptions = c(
        "exact unconditional power: sum over the number of discordant pairs N_D ~ Bin(N, pi_D) of the conditional binomial test power (theta = 1/2 under H0, OR/(1 + OR) under H1)",
        sprintf("conditional rejection regions per N_D follow the alpha-balancing rule; attained alpha is the unconditional size; omitted tail mass of N_D %s", format(signif(tm$omitted, 3))),
        "critical values are not reported: the exact computation uses one rejection region per number of discordant pairs"
      )
    ))
  }
  nd <- .ig_mcnemar_discordant(N, pd)
  t <- .ig_binom_test_power(nd, 0.5, theta, p$alpha, p$tails, p$alpha_balancing)
  crit <- if (identical(p$tails, "one")) {
    if (t$side == "upper") t$region[2] else t$region[1]
  } else t$region
  dist$discordant_pairs <- list(fixed = nd)
  list(
    power = t$power, beta = t$beta, alpha_attained = t$alpha_attained,
    critical = crit,
    n_total = N,
    outputs = c(common_out, list(
      lower_critical_n = if (identical(p$tails, "one")) crit else t$region[1],
      upper_critical_n = if (identical(p$tails, "one")) crit else t$region[2],
      discordant_pairs = nd, computation = 0)),
    statistic_distribution = dist,
    assumptions = c(
      sprintf("fast approximation: the number of discordant pairs is taken as fixed, N_D = ceiling(N pi_D) = %s", format(nd)),
      paste("binomial test of theta = 1/2 on the (yes/no) count among the N_D discordant pairs (", t$rule, ")"),
      if (identical(p$computation, "exact_unconditional")) "N is at or above the exact-computation limit, so the approximation is used (manual 5.2.2)"
    )
  )
}

# Upper bound of the conditional power at N_D = i used by the a priori
# searches: the randomised one-sided most powerful level-alpha test (Neyman-
# Pearson); for two tails with alpha/2 per side (option 1) the tighter
# MP power at level alpha/2 on the side of the effect plus the exact power of
# the opposite tail (each tail of the region has size <= alpha/2).
.ig_binom_power_bound <- function(n, p0, p1, alpha, tails, balancing) {
  if (identical(tails, "two") && identical(balancing %||% "half_each_side", "half_each_side")) {
    opp <- if (p1 >= p0) .ig_binom_cdf(.ig_binom_lower_crit(n, p0, alpha / 2), n, p1) else
      .ig_binom_sf(.ig_binom_upper_crit(n, p0, alpha / 2), n, p1)
    return(min(1, .ig_binom_mp_power(n, p0, p1, alpha / 2) + opp))
  }
  .ig_binom_mp_power(n, p0, p1, alpha)
}

# A priori search of McNemar's test. Conditional terms depend only on N_D = i
# and are cached across N. The bound is the conditional bound averaged over
# N_D (non-decreasing in N: the conditional bound increases with i and N_D
# increases stochastically with N) for the exact computation, and the bound
# at N_D = ceiling(N pi_D) for the fast approximation (the maximum of both
# where the computation switches at the exact-computation limit).
.ig_mcnemar_a_priori_search <- function(p, n_lo, n_max) {
  theta <- p$or / (1 + p$or)
  if (abs(theta - 0.5) <= .ig_prop_eps) {
    return(list(status = "infeasible", reason = paste(
      "with odds ratio 1 the power equals the attained size (<= alpha) for every N, so the target power is never",
      "reached (G*Power excludes OR = 1 in a priori analyses)"), diagnostics = list(rule = "no effect")))
  }
  cache <- new.env(parent = emptyenv())
  cache$n <- 0
  cache$power <- numeric()
  cache$bound <- numeric()
  grow <- function(imax) {
    if (imax <= cache$n) return(invisible())
    i <- seq(cache$n + 1, imax)
    cnd <- .ig_mcnemar_conditional(i, p$alpha, p$tails, p$alpha_balancing, theta)
    cache$power <- c(cache$power, cnd$power)
    cache$bound <- c(cache$bound, vapply(i, function(k) .ig_binom_power_bound(k, 0.5, theta, p$alpha, p$tails,
                                                                              p$alpha_balancing), 0))
    cache$n <- imax
  }
  at <- function(i, what) {
    grow(max(i))
    v <- numeric(length(i))
    pos <- i > 0
    v[pos] <- cache[[what]][i[pos]]
    v
  }
  exact <- identical(p$computation, "exact_unconditional")
  avg <- function(N, what) {
    tm <- .ig_mcnemar_terms(N, p$pd)
    sum(tm$w * at(tm$i, what)) + if (identical(what, "bound")) tm$omitted else 0
  }
  fwd <- function(N) {
    if (exact && N < p$exact_n_limit) return(avg(N, "power"))
    q <- p
    q$n_total <- N
    .ig_power_exact_proportions_mcnemar(q)$power
  }
  bound <- function(N) {
    b_fast <- at(.ig_mcnemar_discordant(N, p$pd), "bound")
    if (!exact) return(b_fast)
    max(b_fast, avg(N, "bound"))
  }
  out <- .ig_bounded_enumeration(fwd, bound, p$power, max(n_lo, 1), n_max,
                                 bound_text = paste("the most powerful conditional level-alpha power (Neyman-Pearson; two",
                                                    "tails with alpha/2 per side: MP at alpha/2 plus the opposite tail),",
                                                    "averaged over the number of discordant pairs for the exact computation"))
  if (identical(out$status, "ok") && !(exact && out$n < p$exact_n_limit)) {
    # Fast approximation: the power depends on N only through N_D =
    # ceiling(N pi_D). G*Power 3.1.9.7 reports N = ceiling(n_D / pi_D) for the
    # smallest sufficient n_D (with the power of n_D), although its post hoc
    # analysis uses ceiling(N pi_D) at that N; kept for transparency.
    nd <- .ig_mcnemar_discordant(out$n, p$pd)
    x <- nd / p$pd
    out$diagnostics$gpower_fast_rule <- list(
      discordant_pairs = nd, n_total = ceiling(x - 8 * .Machine$double.eps * x),
      note = paste("G*Power 3.1.9.7 reports N = ceiling(n_D / pi_D) for the smallest sufficient number of discordant pairs",
                   "n_D; inferguidR reports the smallest N with ceiling(N pi_D) = n_D, the convention of the post hoc analysis"))
  }
  out
}

.ig_method_exact_proportions_mcnemar <- function() {
  ig_method(
    id = "exact_proportions_mcnemar",
    family = "exact",
    label = list(en = "Proportions: Inequality, two dependent groups (McNemar)",
                 de = "Anteile: Ungleichheit, zwei abh\u00e4ngige Gruppen (McNemar)"),
    params = ig_params(
      ig_param("or", role = "effect",
               domain = ig_domain_real(1e-6, 1e6, lower_inclusive = FALSE, upper_inclusive = FALSE),
               label = list(en = "Odds ratio", de = "Odds Ratio"),
               description = "OR = pi12 / pi21 (manual 5.1); admissible range (1e-6, 1e6) as in G*Power.",
               aliases = "odds_ratio"),
      ig_param("pd", role = "nuisance",
               domain = ig_domain_real(.ig_prop_eps, 1 - .ig_prop_eps, lower_inclusive = FALSE, upper_inclusive = FALSE),
               default = 0.3,
               label = list(en = "Prop discordant pairs", de = "Anteil diskordanter Paare"),
               description = "pi_D = pi12 + pi21. Default 0.3 = G*Power's initial field value (certified inventory).",
               aliases = "prop_discordant_pairs"),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      ig_param_n_total(lower = 1),
      ig_param("exact_n_limit", role = "control", domain = ig_domain_real(1, Inf, upper_inclusive = FALSE),
               default = 1e6,
               label = list(en = "Exact computation if N <", de = "Exakte Berechnung falls N <"),
               description = paste("The exact (unconditional) power is used when N is below this limit, the fast",
                                   "approximation otherwise (manual 5.2.2). Default 1e6 = G*Power's initial field",
                                   "value (Options dialog, certified capture)."),
               applies_when = list(option = "computation", "in" = "exact_unconditional"))
    ),
    options = ig_options(
      .ig_prop_tails_option(), .ig_prop_balancing_option(),
      ig_option("computation", values = c("exact_unconditional", "fast_approx"), default = "fast_approx",
                labels = c(exact_unconditional = "exact (unconditional) power if N < limit",
                           fast_approx = "faster approximation (assume number of discordant pairs as constant)"),
                label = list(en = "Computation", de = "Berechnung"))
    ),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = .ig_power_exact_proportions_mcnemar,
    design = NULL,
    discrete = TRUE,
    effect_direction = "increasing",
    engine = list(package = "stats", fun = "dbinom/pbinom/qbinom (sum over discordant pairs; inferguidR regions)",
                  version = NA_character_),
    checker = "validation/independent/g8/check_mcnemar.R",
    reference = list(
      gpower_path = "Tests > Proportions > Two dependent groups: McNemar",
      gpower_label = "Exact - Proportions: Inequality, two dependent groups (McNemar)",
      manual_pages = "ch. 5 (PDF pp. 14-16); options 5.2; implementation notes 5.5",
      allocation_rule = "single sample of pairs: N is the design",
      method_notes = "specs/method-notes/exact_proportions_mcnemar.md"
    ),
    status = "independently_checked",
    constraints = .ig_prop_power_alpha_constraint,
    a_priori_search = .ig_mcnemar_a_priori_search
  )
}

# Fisher's exact test (unconditional power of the conditional test) ------------------------------
#
# G*Power: Tests > Proportions > Two independent groups: inequality, Fisher's
# exact test (manual ch. 6). For each total number of successes m the test
# conditions on m: x1 | m ~ hypergeometric under H0; the rejection set R_m is
# built from the H0 conditional probabilities; the power is the probability
# of the union of the R_m under independent Bin(n1, p1) and Bin(n2, p2)
# (manual 6.6.1). Rules (method notes; probe grid of 120 post hoc captures):
#
# * one tail: the tail of x1 on the side of p1 - p2 (upper when p1 > p2),
#   the largest tail with conditional size <= alpha;
# * two tails, "symmetrical with respect to the statistic": reject the tables
#   whose two-sided conditional p-value (sum of the H0 probabilities of the
#   tables at least as extreme by the statistic) is <= alpha; statistics:
#   Fisher (table probability; ties within a relative 1e-7), Pearson chi-square
#   (|x1 - E x1|) and the likelihood ratio G^2 of manual 6.6.1;
# * two tails, "maximize alpha" (undocumented option; rule from the captures):
#   the two tails of the Fisher-statistic region, then the adjacent table
#   with the smaller probability is added while the conditional size stays
#   <= alpha (equal increments: the tail away from the effect first; G*Power
#   resolves such ties by floating-point rounding);
# * "Actual alpha": the unconditional size at the common proportion p2
#   (G*Power's convention, confirmed on the probe grid).
# * Large-sample approximation (manual 6.6.2): continuity-corrected z test
#   with pooled H0 SD s0 and unpooled H1 SD s1; with c = (1/n1 + 1/n2)/2 the
#   test rejects when |p2^ - p1^| - c > z s0, so the power is
#   Phi((|D| - c - z s0)/s1) on the side of the effect plus, for two tails,
#   Phi((-|D| - c - z s0)/s1) (the correction moves both tails outwards).

.ig_fisher_up_dir <- function(p1, p2) p1 > p2

# Rejection indicator of one margin (x = support of x1, h = H0 conditional
# probabilities), the effect direction up (x1 large extreme) and the rule.
.ig_fisher_margin_region <- function(x, h, m, n1, n2, alpha, tails, rule, up) {
  K <- length(h)
  lim <- alpha * (1 + 1e-12)
  if (identical(tails, "one")) {
    pv <- if (up) rev(cumsum(rev(h))) else cumsum(h)
    return(pv <= lim)
  }
  fisher_region <- function() {
    # p(x) = sum of h(y) over h(y) <= h(x) (1 + 1e-7), as in stats::fisher.test.
    o <- order(h)
    hs <- h[o]
    cs <- cumsum(hs)
    idx <- findInterval(h * (1 + 1e-7), hs)
    pv <- cs[pmax(idx, 1)]
    pv[idx == 0] <- 0
    pv <= lim
  }
  if (identical(rule, "symmetric_fisher_statistic")) return(fisher_region())
  if (identical(rule, "maximize_alpha")) {
    # Start from the two tails of the Fisher-statistic region, then add the
    # adjacent table with the smaller probability while the conditional size
    # stays <= alpha; equal increments: the tail away from the effect first.
    rej <- fisher_region()
    if (all(rej)) return(rej)
    a <- if (rej[1]) which(!rej)[1] - 1 else 0
    b <- if (rej[K]) K + 2 - which(!rev(rej))[1] else K + 1
    cl <- cumsum(h)
    cu <- rev(cumsum(rev(h)))
    repeat {
      tot <- (if (a > 0) cl[a] else 0) + (if (b <= K) cu[b] else 0)
      il <- if (a + 1 < b) h[a + 1] else Inf
      iu <- if (b - 1 > a) h[b - 1] else Inf
      st <- min(il, iu)
      if (!is.finite(st) || tot + st > lim) break
      tie <- is.finite(il) && is.finite(iu) && abs(il - iu) <= 1e-10 * max(il, iu)
      if ((!tie && il < iu) || (tie && up)) a <- a + 1 else b <- b - 1
    }
    return(seq_len(K) <= a | seq_len(K) >= b)
  }
  # Ordering statistic: larger = more extreme (manual 6.6.1).
  stat <- switch(rule,
    symmetric_pearson_chi2_statistic = abs(x - n1 * m / (n1 + n2)),
    symmetric_lr_statistic = {
      N <- n1 + n2
      e <- c(n1 * m / N, n1 * (N - m) / N, n2 * m / N, n2 * (N - m) / N)
      xl <- function(o, ee) ifelse(o > 0, o * log(o / ee), 0)
      2 * (xl(x, e[1]) + xl(n1 - x, e[2]) + xl(m - x, e[3]) + xl(n2 - m + x, e[4]))
    })
  o <- order(-stat)
  cs <- cumsum(h[o])
  # number of tables with stat >= stat(x) (relative tie tolerance 1e-7):
  # K minus the number of values below the threshold
  thr <- stat - 1e-7 * abs(stat)
  cnt <- K - findInterval(thr, sort(stat), left.open = TRUE)
  pv <- cs[pmax(cnt, 1)]
  pv <= lim
}

# Margins with non-negligible probability under H1 or under the common
# proportion p2 (attained alpha); the omitted mass bounds the truncation error.
.ig_fisher_margins <- function(n1, n2, p1, p2) {
  N <- n1 + n2
  if (N <= 400) return(list(m = 0:N, omitted = 0))
  mu <- n1 * p1 + n2 * p2
  sd1 <- sqrt(n1 * p1 * (1 - p1) + n2 * p2 * (1 - p2))
  lo <- min(stats::qbinom(1e-18, N, p2), floor(mu - 12 * sd1 - 5))
  hi <- max(stats::qbinom(1e-18, N, p2, lower.tail = FALSE), ceiling(mu + 12 * sd1 + 5))
  lo <- max(lo, 0)
  hi <- min(hi, N)
  omitted_h0 <- (if (lo > 0) stats::pbinom(lo - 1, N, p2) else 0) + stats::pbinom(hi, N, p2, lower.tail = FALSE)
  list(m = lo:hi, omitted = omitted_h0)
}

# Exact unconditional power, beta and attained alpha (at the common
# proportion p2). `bound = TRUE` also returns the unconditional power of the
# randomised one-sided conditional most powerful test (the UMPU test) in the
# direction of the effect, an upper bound of the power of every conditional
# level-alpha test (used by the a priori search). With `stop_below` (a target
# power) the margins are processed from the centre of the H1 distribution of
# the total number of successes outwards and the evaluation stops as soon as
# the power accumulated so far plus the H1 mass not yet processed is below
# the target: the design then provably misses it (`reaches = FALSE`; no
# power is returned).
.ig_fisher_exact <- function(n1, n2, p1, p2, alpha, tails, rule, bound = FALSE, stop_below = NULL) {
  n1 <- as.double(n1)
  n2 <- as.double(n2)
  N <- n1 + n2
  up <- .ig_fisher_up_dir(p1, p2)
  mg <- .ig_fisher_margins(n1, n2, p1, p2)
  ms <- mg$m
  if (!is.null(stop_below)) ms <- ms[order(abs(ms - (n1 * p1 + n2 * p2)))]
  power <- 0
  beta <- 0
  att <- 0
  bnd <- 0
  h1_mass <- 0
  for (m in ms) {
    s_lo <- max(0, m - n2)
    s_hi <- min(n1, m)
    if (s_hi - s_lo > 200) {
      # Window covering the central H0 conditional range and the H1 range of
      # x1 (each beyond 1e-20 in both tails); tables outside have H0
      # conditional probability < 1e-20 (they are rejected by every rule, their
      # mass is negligible in the p-values) and H1 probability < 1e-20.
      xl <- min(stats::qhyper(1e-20, n1, n2, m), stats::qbinom(1e-20, n1, p1))
      xh <- max(stats::qhyper(1e-20, n1, n2, m, lower.tail = FALSE), stats::qbinom(1e-20, n1, p1, lower.tail = FALSE))
      x <- seq(max(s_lo, xl), min(s_hi, xh))
      h <- stats::dhyper(x, n1, n2, m)
    } else {
      x <- seq(s_lo, s_hi)
      h <- stats::dhyper(x, n1, n2, m)
      h <- h / sum(h)
    }
    rej <- .ig_fisher_margin_region(x, h, m, n1, n2, alpha, tails, rule, up)
    lt <- stats::dbinom(x, n1, p1, log = TRUE) + stats::dbinom(m - x, n2, p2, log = TRUE)
    pt <- exp(lt)
    power <- power + sum(pt[rej])
    beta <- beta + sum(pt[!rej])
    h1_mass <- h1_mass + sum(pt)
    if (!is.null(stop_below) && power + (1 - h1_mass) + 1e-12 < stop_below) {
      return(list(power = NA_real_, reaches = FALSE, power_upper = power + (1 - h1_mass)))
    }
    att <- att + stats::dbinom(m, N, p2) * sum(h[rej])
    if (bound) {
      # conditional MP (randomised) test of psi = 1 against the H1 odds ratio
      ord <- if (up) rev(seq_along(x)) else seq_along(x)
      cs <- cumsum(h[ord])
      full <- cs <= alpha * (1 + 1e-12)
      k <- sum(full)
      g1 <- if (k < length(x)) (alpha - (if (k > 0) cs[k] else 0)) / h[ord][k + 1] else 0
      pm <- sum(pt)
      cond <- if (pm > 0) pt[ord] / pm else rep(0, length(x))
      mp <- (if (k > 0) sum(cond[seq_len(k)]) else 0) + (if (k < length(x)) max(min(g1, 1), 0) * cond[k + 1] else 0)
      bnd <- bnd + pm * mp
    }
  }
  omitted_h1 <- max(0, 1 - h1_mass)
  list(power = min(power, 1), beta = min(beta + 0, 1), alpha_attained = min(att, 1),
       bound = if (bound) min(1, bnd + omitted_h1) else NA_real_,
       omitted = max(mg$omitted, omitted_h1), margins = length(mg$m), reaches = TRUE)
}

.ig_fisher_large_sample <- function(n1, n2, p1, p2, alpha, tails) {
  n1 <- as.double(n1)
  n2 <- as.double(n2)
  D <- p2 - p1
  pbar <- (n1 * p1 + n2 * p2) / (n1 + n2)
  s0 <- sqrt(pbar * (1 - pbar) * (1 / n1 + 1 / n2))
  s1 <- sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
  cc <- (1 / n1 + 1 / n2) / 2
  sgn <- if (D >= 0) 1 else -1
  mu <- sgn * (abs(D) - cc) / s0      # H1 mean of Z on the side of the effect (correction toward 0)
  sdz <- s1 / s0
  zc <- if (identical(tails, "two")) stats::qnorm(alpha / 2, lower.tail = FALSE) else stats::qnorm(alpha, lower.tail = FALSE)
  # The corrected test rejects when |p2^ - p1^| - cc > zc s0. With the estimated
  # difference ~ N(D, s1^2) under H1 the tail on the side of the effect has
  # probability Phi((|D| - cc - zc s0) / s1); the opposite tail is corrected
  # away from the effect as well: Phi((-|D| - cc - zc s0) / s1).
  a_eff <- (abs(D) - cc - zc * s0) / s1
  if (identical(tails, "two")) {
    a_opp <- (-abs(D) - cc - zc * s0) / s1
    power <- stats::pnorm(a_eff) + stats::pnorm(a_opp)
    beta <- max(stats::pnorm(a_eff, lower.tail = FALSE) - stats::pnorm(a_opp), 0)
  } else {
    power <- stats::pnorm(a_eff)
    beta <- stats::pnorm(a_eff, lower.tail = FALSE)
  }
  list(power = min(power, 1), beta = beta, alpha_attained = alpha, zc = zc, mu = mu, sdz = sdz, s0 = s0, sgn = sgn)
}

.ig_fisher_is_exact <- function(p, N) identical(p$computation, "exact") && N < p$exact_n_limit

# Sensitivity branch: with the option effect_direction the scan for p1 is
# restricted to p1 <= p2 or p1 >= p2 by holding p1 at p2 on the other side (a
# flat continuation with power = attained size, below every admissible
# target), so that the root is unique.
.ig_fisher_p1 <- function(p) {
  p1 <- p$p1
  if (identical(p$effect_direction, "p1_le_p2")) p1 <- min(p1, p$p2)
  if (identical(p$effect_direction, "p1_ge_p2")) p1 <- max(p1, p$p2)
  p1
}

.ig_power_exact_proportions_fisher <- function(p) {
  n1 <- as.double(p$n1)
  n2 <- as.double(p$n2)
  N <- n1 + n2
  p1 <- .ig_fisher_p1(p)
  p2 <- p$p2
  rule <- p$two_sided_balancing %||% "symmetric_fisher_statistic"
  if (.ig_fisher_is_exact(p, N)) {
    ex <- .ig_fisher_exact(n1, n2, p1, p2, p$alpha, p$tails, rule)
    return(list(
      power = ex$power, beta = ex$beta, alpha_attained = ex$alpha_attained,
      n_total = N, group_sizes = c(n1 = n1, n2 = n2),
      outputs = list(proportion_p1 = p1, proportion_p2 = p2, computation = 1, truncated_mass = ex$omitted),
      statistic_distribution = list(
        h0 = list(family = "hypergeometric (conditional on the total number of successes)"),
        h1 = list(family = "product binomial", size = c(n1, n2), prob = c(p1, p2))),
      assumptions = c(
        "exact unconditional power of the conditional test (manual 6.6.1): rejection sets per total number of successes from the H0 hypergeometric distribution",
        if (identical(p$tails, "two")) sprintf("two-sided rule: %s", rule)
        else "one-sided: the tail of x1 on the side of p1 - p2",
        "attained alpha: unconditional size at the common proportion p2 (G*Power's convention)",
        sprintf("omitted probability mass of the truncated margins: %s", format(signif(ex$omitted, 3)))
      )
    ))
  }
  ls <- .ig_fisher_large_sample(n1, n2, p1, p2, p$alpha, p$tails)
  crit <- if (identical(p$tails, "two")) c(-ls$zc, ls$zc) else ls$sgn * ls$zc
  list(
    power = ls$power, beta = ls$beta, alpha_attained = p$alpha,
    critical = crit,
    n_total = N, group_sizes = c(n1 = n1, n2 = n2),
    outputs = list(proportion_p1 = p1, proportion_p2 = p2, computation = 0),
    statistic_distribution = list(h0 = list(family = "normal", mean = 0, sd = 1),
                                  h1 = list(family = "normal", mean = ls$mu, sd = ls$sdz)),
    assumptions = c(
      "large-sample approximation (manual 6.6.2): continuity-corrected z test with pooled H0 variance and unpooled H1 variance; the correction (1/n1 + 1/n2)/2 moves each rejection boundary away from 0",
      if (identical(p$computation, "exact")) "N is at or above the exact-computation limit, so the approximation is used"
    )
  )
}

# Allocation rule of G*Power 3.1.9.7 for this procedure (captures
# SUP-G8-FISHER-APR-A01..A10 and the campaign a priori records): the design
# is indexed by the size of group 1, n2 = floor(r n1 + 1/2) (round half up),
# N = n1 + n2. A total N is admissible when it has that form (for r = 1: the
# even N); it then determines n1 uniquely.
.ig_fisher_n2 <- function(n1, r) floor(r * n1 + 0.5)

.ig_design_exact_proportions_fisher <- function(p, n_total, continuous = FALSE) {
  r <- as.double(p$allocation_ratio)
  if (continuous) return(c(n1 = n_total / (1 + r), n2 = r * n_total / (1 + r)))
  c0 <- n_total / (1 + r)
  for (n1 in unique(c(floor(c0), ceiling(c0), floor(c0) - 1, ceiling(c0) + 1))) {
    if (n1 < 1) next
    n2 <- .ig_fisher_n2(n1, r)
    if (n2 >= 1 && n1 + n2 == n_total) return(c(n1 = n1, n2 = n2))
  }
  NULL
}

# A priori search over the size of group 1 (the designs n1, n2 = round(r n1)).
# Exact computation: bounded enumeration, with the unconditional power of the
# randomised conditional one-sided UMPU test as the bound (an upper bound of
# every conditional level-alpha test; non-decreasing in n1 because both group
# sizes are non-decreasing in n1 and a smaller design's test is a level-alpha
# conditional test of the larger one that ignores the additional
# observations). Every design from the first one whose bound reaches the
# target is evaluated; designs that provably miss the target are abandoned
# early (.ig_fisher_exact(stop_below =)). Large-sample approximation: the
# power is non-decreasing in n1 (doubling and bisection; G*Power's criterion
# with a fractional second group, see below).
.ig_fisher_a_priori_search <- function(p, n_lo, n_max) {
  if (abs(p$p1 - p$p2) <= 0) {
    return(list(status = "infeasible", reason = paste(
      "with p1 = p2 the power equals the attained size (<= alpha) for every N, so the target power is never reached"),
      diagnostics = list(rule = "no effect")))
  }
  r <- as.double(p$allocation_ratio)
  target <- p$power
  rule <- p$two_sided_balancing %||% "symmetric_fisher_statistic"
  n_of <- function(n1) n1 + .ig_fisher_n2(n1, r)
  n1_min <- max(1, ceiling(0.5 / r - 1e-9))
  while (.ig_fisher_n2(n1_min, r) < 1 || n_of(n1_min) < n_lo) n1_min <- n1_min + 1
  n1_max <- floor(n_max / (1 + r))
  while (n1_max >= n1_min && n_of(n1_max) > n_max) n1_max <- n1_max - 1
  if (n1_max < n1_min) return(list(status = "infeasible", reason = "no admissible design", diagnostics = list()))
  ls_power <- function(n1) .ig_fisher_large_sample(n1, .ig_fisher_n2(n1, r), p$p1, p$p2, p$alpha, p$tails)$power
  # Large-sample branch: G*Power 3.1.9.7 takes the smallest total N* whose
  # power with the fractional group sizes N/(1+r), rN/(1+r) reaches the target
  # and displays the smallest design n1, n2 = floor(r n1 + 1/2) with
  # n1 + n2 >= N*; that integer design can miss the target (capture
  # SUP-G8-FISHER-APR-A09). It is reported when its own power reaches the
  # target, otherwise the smallest larger n1 whose integer design does;
  # G*Power's design is kept in diagnostics$gpower_allocation_rule.
  ls_frac <- function(N) .ig_fisher_large_sample(N / (1 + r), r * N / (1 + r), p$p1, p$p2, p$alpha, p$tails)$power
  gp_rule <- NULL
  ls_search <- function(from) {
    if (from > n1_max) return(NA_real_)
    n_star <- .ig_monotone_first_n(ls_frac, target, max(2, n_of(from)), n_max)
    if (is.na(n_star)) return(NA_real_)
    n1s <- max(from, floor(n_star / (1 + r)) - 1)
    while (n1s <= n1_max && n_of(n1s) < n_star) n1s <- n1s + 1
    if (n1s > n1_max) return(NA_real_)
    pw <- ls_power(n1s)
    gp_rule <<- list(n_total = n_of(n1s), n1 = n1s, n2 = .ig_fisher_n2(n1s, r), fractional_n_total = n_star,
                     fractional_power = ls_frac(n_star), integer_design_power = pw, reaches_target = pw >= target)
    if (pw >= target) return(n1s)
    if (n1s + 1 > n1_max) return(NA_real_)
    .ig_monotone_first_n(ls_power, target, n1s + 1, n1_max)
  }
  ls_result <- function(n1, search) {
    adjusted <- !is.null(gp_rule) && !isTRUE(gp_rule$reaches_target)
    note <- if (adjusted) sprintf(paste(
      "G*Power's rule gives n1 = %s, n2 = %s, whose power %s is below the target %s (G*Power 3.1.9.7 reports it",
      "without a message); inferguidR reports the smallest larger design that reaches the target"),
      format(gp_rule$n1, scientific = FALSE), format(gp_rule$n2, scientific = FALSE),
      format(signif(gp_rule$integer_design_power, 7)), format(target)) else character()
    list(status = "ok", n = n_of(n1), power = ls_power(n1), warnings = note, diagnostics = c(list(
      search = search,
      minimality = paste0("G*Power's criterion: the smallest design n1, n2 = floor(r n1 + 1/2) whose total is at least the ",
                          "smallest N whose fractional-allocation power reaches the target (that power is non-decreasing ",
                          "in N)", if (adjusted) "; adjusted upwards to the smallest n1 whose integer design reaches the target" else ""),
      gpower_allocation_rule = gp_rule),
      if (adjusted) list(adjusted_from_gpower_rule = TRUE, adjustment_note = note)))
  }
  alloc <- "designs n1, n2 = floor(r n1 + 1/2) (G*Power's allocation rule for this procedure)"
  if (identical(p$computation, "exact")) {
    # exact regime: the n1 whose total is below the exact-computation limit
    n1_lim <- n1_min
    if (.ig_fisher_is_exact(p, n_of(n1_min))) {
      n1_lim <- min(n1_max, ceiling(p$exact_n_limit / (1 + r)) + 1)
      while (n1_lim >= n1_min && !.ig_fisher_is_exact(p, n_of(n1_lim))) n1_lim <- n1_lim - 1
      bound_at <- function(n1) .ig_fisher_exact(n1, .ig_fisher_n2(n1, r), p$p1, p$p2, p$alpha, p$tails, rule, bound = TRUE)$bound
      nb <- .ig_monotone_first_n(bound_at, target, n1_min, n1_lim)
      if (!is.na(nb)) {
        count <- 0
        for (n1 in seq(nb, n1_lim)) {
          e <- .ig_fisher_exact(n1, .ig_fisher_n2(n1, r), p$p1, p$p2, p$alpha, p$tails, rule, stop_below = target)
          count <- count + 1
          if (isTRUE(e$reaches) && e$power >= target) {
            return(list(status = "ok", n = n_of(n1), power = e$power, diagnostics = list(
              search = paste("bounded enumeration over the", alloc),
              minimality = sprintf(paste(
                "proven: designs with n1 < %s fail because the unconditional power of the randomised conditional",
                "one-sided UMPU test (an upper bound, non-decreasing in n1) is below the target there; every design",
                "from n1 = %s upwards was evaluated"), format(nb, scientific = FALSE), format(nb, scientific = FALSE)),
              monotonicity = "not assumed (discrete test; power is not monotone in the sample size)",
              n1_bound = nb, evaluations = count)))
          }
        }
      }
      n1_lim <- n1_lim + 1
    }
    # designs at or above the exact-computation limit use the approximation
    n1 <- ls_search(n1_lim)
    if (is.na(n1)) {
      return(list(status = "infeasible", reason = sprintf(
        "the target power %s is not reached by any admissible design up to N = %s", format(target),
        format(n_max, scientific = FALSE)), diagnostics = list(allocation = alloc)))
    }
    return(ls_result(n1, paste("exact computation below the exact-computation limit (no design reaches the target",
                               "there); large-sample approximation above it: doubling and bisection over the", alloc)))
  }
  n1 <- ls_search(n1_min)
  if (is.na(n1)) {
    return(list(status = "infeasible", reason = sprintf(
      "the target power %s is not reached by any admissible design up to N = %s (power there: %s)", format(target),
      format(n_max, scientific = FALSE), format(signif(ls_power(n1_max), 7))), diagnostics = list(allocation = alloc)))
  }
  ls_result(n1, paste("doubling and bisection over the", alloc))
}

.ig_method_exact_proportions_fisher <- function() {
  prop <- ig_domain_real(0, 1, lower_inclusive = FALSE, upper_inclusive = FALSE)
  ig_method(
    id = "exact_proportions_fisher",
    family = "exact",
    label = list(en = "Proportions: Inequality, two independent groups (Fisher's exact test)",
                 de = "Anteile: Ungleichheit, zwei unabh\u00e4ngige Gruppen (exakter Test nach Fisher)"),
    params = ig_params(
      ig_param("p1", role = "effect", domain = prop,
               label = list(en = "Proportion p1", de = "Anteil p1"),
               description = "Success proportion of group 1 (H1); the effect is p1 versus p2 (manual 6.2).",
               aliases = "proportion_p1"),
      ig_param("p2", role = "nuisance", domain = prop, default = 0.6,
               label = list(en = "Proportion p2", de = "Anteil p2"),
               description = "Success proportion of group 2. Default 0.6 = G*Power's initial field value (certified inventory).",
               aliases = "proportion_p2"),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      ig_param("n_total", role = "design", domain = ig_domain_integer(lower = 2), unit = "subjects",
               label = list(en = "Total sample size", de = "Gesamtstichprobenumfang"),
               description = "Solved in a priori mode (n1 + n2).", modes = "a_priori"),
      ig_param("n1", role = "design", domain = ig_domain_integer(lower = 1), unit = "subjects",
               label = list(en = "Sample size group 1", de = "Stichprobenumfang Gruppe 1"),
               aliases = "sample_size_group_1", modes = c("post_hoc", "sensitivity", "criterion", "compromise")),
      ig_param("n2", role = "design", domain = ig_domain_integer(lower = 1), unit = "subjects",
               label = list(en = "Sample size group 2", de = "Stichprobenumfang Gruppe 2"),
               aliases = "sample_size_group_2", modes = c("post_hoc", "sensitivity", "criterion", "compromise")),
      ig_param("allocation_ratio", role = "design",
               domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE), default = 1,
               label = list(en = "Allocation ratio N2/N1", de = "Allokationsverh\u00e4ltnis N2/N1"),
               description = "Group-size ratio n2/n1 of the a priori design; default 1 = G*Power's initial field value.",
               aliases = "allocation_ratio_n2_n1", modes = "a_priori"),
      ig_param("exact_n_limit", role = "control", domain = ig_domain_real(2, Inf, upper_inclusive = FALSE),
               default = 10000,
               label = list(en = "Use exact distribution if N <", de = "Exakte Verteilung falls N <"),
               description = paste("The exact computation is used when n1 + n2 is below this limit, the large-sample",
                                   "approximation otherwise. Default 10000 = G*Power's initial field value (Options",
                                   "dialog, certified capture)."),
               applies_when = list(option = "computation", "in" = "exact"))
    ),
    options = ig_options(
      .ig_prop_tails_option(),
      ig_option("computation", values = c("exact", "large_sample_z"), default = "exact",
                labels = c(exact = "Use exact distribution if N < limit", large_sample_z = "Use large sample approximation (z-Test)"),
                label = list(en = "Computation", de = "Berechnung")),
      ig_option("two_sided_balancing",
                values = c("symmetric_fisher_statistic", "symmetric_lr_statistic", "symmetric_pearson_chi2_statistic",
                           "maximize_alpha"),
                default = "symmetric_fisher_statistic",
                labels = c(symmetric_fisher_statistic = "Symmetrical with respect to the Fisher statistic",
                           symmetric_lr_statistic = "Symmetrical with respect to the likelihood ratio statistic",
                           symmetric_pearson_chi2_statistic = "Symmetrical with respect to the Pearson chi square statistic",
                           maximize_alpha = "Maximize alpha"),
                applies_when = list(all = list(list(option = "tails", "in" = "two"),
                                               list(option = "computation", "in" = "exact"))),
                label = list(en = "Alpha-balancing in two-sided exact tests")),
      ig_option("effect_direction", values = c("p1_le_p2", "p1_ge_p2"), default = "p1_le_p2",
                labels = c(p1_le_p2 = "p1 <= p2", p1_ge_p2 = "p1 >= p2"),
                applies_when = list(mode = "sensitivity"),
                label = list(en = "Effect direction", de = "Effektrichtung"))
    ),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = .ig_power_exact_proportions_fisher,
    design = .ig_design_exact_proportions_fisher,
    discrete = TRUE,
    effect_direction = "unknown",
    engine = list(package = "stats", fun = "dhyper/dbinom (exact enumeration over margins); pnorm (large-sample)",
                  version = NA_character_),
    checker = "validation/independent/g8/check_fisher.R",
    reference = list(
      gpower_path = "Tests > Proportions > Two independent groups: inequality, Fisher's exact test",
      gpower_label = "Exact - Proportions: Inequality, two independent groups (Fisher's exact test)",
      manual_pages = "ch. 6 (PDF p. 17); implementation notes 6.6.1 (exact), 6.6.2 (large sample)",
      allocation_rule = paste("designs are indexed by the size of group 1: n2 = floor(r n1 + 1/2), N = n1 + n2",
                              "(G*Power's rule for this procedure; r = 1: equal groups); the smallest such design",
                              "reaching the target"),
      method_notes = "specs/method-notes/exact_proportions_fisher.md"
    ),
    status = "independently_checked",
    constraints = .ig_prop_power_alpha_constraint,
    a_priori_search = .ig_fisher_a_priori_search
  )
}
