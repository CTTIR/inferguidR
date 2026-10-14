# Method definitions: MANOVA global effects and special effects (F family) ---------------
#
# Ownership: R/methods-manova.R belongs to group G4 (blueprint phases 06-10 for
# f_manova_global and f_manova_special). Effect-size drawers: R/effect-manova.R.
# Method notes (mathematical contract, evidence, G*Power defects):
# specs/method-notes/f_manova_global.md, specs/method-notes/f_manova_special.md.
#
# G*Power 3.1.9.7: Tests > Means > Multivariate: "MANOVA: Global effects" and
# "MANOVA: Special effects and interactions" (F tests). The manual has no
# chapter on either procedure. The model below was determined from the
# certified live inventory (specs/gpower-capabilities.json), the emulator
# captures (validation/fixtures/emulator_capture/GPW-F-MANG-*, GPW-F-MANS-*;
# supplement captures under validation/fixtures/emulator_capture_supplement/g4/)
# and the primary sources of the four F approximations (Rao 1951; Pillai and
# Samson 1959; McKeon 1974; Pillai and Mijares 1959) and of the two
# noncentrality algorithms (Muller and Peterson 1984; O'Brien and Shieh 1999).
#
# Notation: N total sample size, k number of groups, p response variables,
# q hypothesis df per response (global effects: q = k - 1; special effects:
# q = number of predictors), error df nu_e = N - k, s = min(p, q),
# df1 = p q. Every approximation writes its criterion as F = f2 df2 / df1, so
# the effect size f2 (G*Power's f²(V), f²(U), f²(T1), f²(T2)) is defined
# relative to the chosen approximation:
#
#   Pillai V (Pillai & Mijares 1959)      df2 = s (nu_e - p + s)            g = s   V = s f2 / (1 + f2)
#   Hotelling-Lawley T1 (Pillai & Samson) df2 = s (nu_e - p - 1) + 2        g = s   T = s f2
#   Wilks U (Rao 1951)                    df2 = t (nu_e - (p - q + 1) / 2) - (p q - 2) / 2
#                                         t = sqrt((p^2 q^2 - 4) / (p^2 + q^2 - 5)), t = 1 when
#                                         p^2 + q^2 - 5 <= 0                g = t   U = (1 + f2)^(-t)
#   Hotelling-Lawley T2 (McKeon 1974)     b = (nu_e + q - p - 1)(nu_e - 1) / ((nu_e - p - 3)(nu_e - p))
#                                         df2 = 4 + (p q + 2) / (b - 1)
#                                         g = (df2 - 2) / (nu_e - p - 1)   T = g f2
#
# Noncentrality: Muller & Peterson (1984) lambda = df2 f2; O'Brien & Shieh
# (1999) lambda = N g f2 (the error df in the leading term g nu_e of df2 is
# replaced by N). The test rejects when F > F_crit, the upper-alpha quantile
# of F(df1, df2); power = P(F(df1, df2, lambda) > F_crit). There is no tail
# option.
#
# Design: global effects use equal groups, so admissible a priori designs are
# multiples of k (G*Power reports N as a multiple of k; given a total that is
# not a multiple it notes that "the analysis was based on the average group
# size", which leaves nu_e = N - k unchanged). Special effects use N alone.
# An a priori design is admissible when N > k and the approximation is defined
# (df2 > 0; for O'Brien-Shieh also 0 < g < Inf). For Hotelling-Lawley T2 an
# a priori design must also lie in the domain of McKeon's approximation,
# n = (nu_e - p - 1) / 2 > 0, i.e. nu_e >= p + 2, where the scale c is positive
# and finite (the same designs O'Brien-Shieh admits). Below it the formula
# still returns numbers under Muller-Peterson (c is 0/0 at nu_e = p + 1 and
# negative for nu_e <= p), and for nu_e < p their power can exceed that of
# every design up to much larger N; such a design is not an answer to "the
# smallest N" (captures REV-G4-MANS-APR-S17, REV-G4-MANG-APR-T12: G*Power
# returns the design with nu_e >= p + 2). inferguidR returns the smallest
# admissible design that reaches the target. G*Power's own a priori search
# ends in "Internal math error" for a number of inputs whose solution lies
# near the smallest designs (method notes, "A priori search").
#
# Numerical note: stats::qf() returns qchisq(p, df1) / df1 whenever df2 > 4e5
# (nmath qf.c), which moves the critical F in the fifth significant digit for
# the large designs of the fixtures (df2 about 1e6). The critical value is
# therefore found by inverting stats::pf() there (.ig_manova_q_f_upper());
# a change request for R/num-dist.R (.ig_q_f) is pending.

.ig_define_methods_F_g4 <- function() {
  list(.ig_method_f_manova("global"), .ig_method_f_manova("special"))
}

.ig_manova_approximations <- c(
  wilks_u = "Wilks U (Rao, 1951)",
  hotelling_lawley_t1 = "Hotelling-Lawley T1 (Pillai and Samson, 1959)",
  hotelling_lawley_t2 = "Hotelling-Lawley T2 (McKeon, 1974)",
  pillai_v = "Pillai V (Pillai and Mijares, 1959)"
)

.ig_manova_algorithms <- c(
  obrien_shieh_1999 = "O'Brien and Shieh (1999), recommended",
  muller_peterson_1984 = "Muller and Peterson (1984)"
)

# G*Power's label of the criterion output (the effect-size field is labelled
# f²(V), f²(U), f²(T1), f²(T2) accordingly).
.ig_manova_criterion_labels <- c(pillai_v = "Pillai V", wilks_u = "Wilks U",
                                 hotelling_lawley_t1 = "Hotelling-Lawley T1",
                                 hotelling_lawley_t2 = "Hotelling-Lawley T2")

.ig_method_f_manova <- function(effect) {
  global <- identical(effect, "global")
  id <- if (global) "f_manova_global" else "f_manova_special"
  label <- if (global) {
    list(en = "MANOVA: Global effects", de = "MANOVA: Globale Effekte")
  } else {
    list(en = "MANOVA: Special effects and interactions",
         de = "MANOVA: Spezielle Effekte und Interaktionen")
  }
  q_text <- if (global) "q = k - 1 (the global group effect)" else "q = number of predictors"
  params <- list(
    ig_param("f2", role = "effect",
             domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE),
             label = list(en = "Effect size f\u00b2", de = "Effektst\u00e4rke f\u00b2"),
             description = paste(
               "Effect size f2 of the chosen approximation (G*Power: f\u00b2(V), f\u00b2(U), f\u00b2(T1) or",
               "f\u00b2(T2)): the F statistic of the approximation is f2 df2 / df1. The same number means",
               "a different population effect under each approximation (Pillai V = s f2 / (1 + f2),",
               "Wilks U = (1 + f2)^(-t), Hotelling-Lawley T1 = s f2, T2 = c f2)."),
             aliases = "effect_size_f2"),
    ig_param_alpha(),
    ig_param_power(),
    ig_param_error_ratio(),
    ig_param("n_total", role = "design", domain = ig_domain_integer(lower = 3), unit = "subjects",
             label = list(en = "Total sample size", de = "Gesamtstichprobenumfang"),
             description = paste(
               "Total sample size N (> number of groups); error df nu_e = N - k.",
               if (global) "Solved in a priori mode as a multiple of k (equal groups)." else
                 "Solved in a priori mode (any integer N)."),
             aliases = "total_sample_size"),
    ig_param("n_groups", role = "design", domain = ig_domain_integer(lower = 2), unit = "groups",
             label = list(en = "Number of groups", de = "Anzahl der Gruppen"),
             description = paste("Number of groups (cells) k of the design; error df nu_e = N - k.",
                                 if (global) "Hypothesis df per response q = k - 1." else ""),
             aliases = "number_of_groups"),
    if (!global) {
      ig_param("n_predictors", role = "design", domain = ig_domain_integer(lower = 1),
               label = list(en = "Number of predictors", de = "Anzahl der Pr\u00e4diktoren"),
               description = paste("Hypothesis (numerator) df q of the tested special effect or",
                                   "interaction per response variable; df1 = p q."),
               aliases = "number_of_predictors")
    },
    ig_param("n_responses", role = "design", domain = ig_domain_integer(lower = 1),
             label = list(en = "Response variables", de = "Antwortvariablen"),
             description = "Number of response (dependent) variables p.",
             aliases = "response_variables")
  )
  ig_method(
    id = id,
    family = "F",
    label = label,
    params = do.call(ig_params, Filter(Negate(is.null), params)),
    options = ig_options(
      ig_option("approximation", values = names(.ig_manova_approximations), default = "pillai_v",
                labels = .ig_manova_approximations,
                label = list(en = "Approximation (F-Transformation)", de = "Approximation (F-Transformation)")),
      ig_option("algorithm", values = names(.ig_manova_algorithms), default = "obrien_shieh_1999",
                labels = .ig_manova_algorithms,
                label = list(en = "Algorithm (noncentrality)", de = "Algorithmus (Nichtzentralit\u00e4t)"))
    ),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = function(p) .ig_power_f_manova(p, effect),
    design = function(p, n_total, continuous = FALSE) .ig_design_f_manova(p, n_total, effect, continuous),
    discrete = FALSE,
    effect_direction = "increasing",
    engine = list(package = "stats",
                  fun = paste("pf (noncentral: Poisson mixture of pbeta, R/num-dist.R .ig_p_f);",
                              "critical F by qf, or by inverting pf when df2 > 4e5"),
                  version = NA_character_),
    checker = "validation/independent/g4/check_f_manova.R",
    reference = list(
      gpower_path = paste("Tests > Means > Multivariate:", label$en),
      gpower_label = paste("F tests -", label$en),
      manual_pages = "no manual chapter (procedure documented only in the application)",
      allocation_rule = if (global) {
        paste("equal groups: admissible a priori designs are multiples of k with N > k at which the",
              "approximation is defined (df2 > 0; O'Brien-Shieh also needs 0 < g < Inf; Hotelling-Lawley T2",
              "also needs error df >= p + 2)")
      } else {
        paste("any N > k at which the approximation is defined (df2 > 0; O'Brien-Shieh also needs",
              "0 < g < Inf; Hotelling-Lawley T2 also needs error df >= p + 2); G*Power uses N alone",
              "(error df N - k)")
      },
      hypothesis_df = q_text,
      sources = paste("Rao (1951); Pillai and Samson (1959); McKeon (1974); Pillai and Mijares (1959);",
                      "Muller and Peterson (1984); O'Brien and Shieh (1999)"),
      method_notes = sprintf("specs/method-notes/%s.md", id)
    ),
    status = "independently_checked",
    constraints = function(p) .ig_constraints_f_manova(p, effect),
    a_priori_search = function(p, n_lo, n_max) .ig_a_priori_f_manova(p, n_lo, n_max, effect)
  )
}

# Hypothesis df per response.
.ig_manova_q <- function(p, effect) {
  if (identical(effect, "global")) as.double(p$n_groups) - 1 else as.double(p$n_predictors)
}

# Degrees of freedom and the O'Brien-Shieh multiplier g of one approximation
# (see the header). Plain IEEE arithmetic as in the formulas: a zero
# denominator gives an infinite b (McKeon's df2 then tends to 4).
.ig_manova_terms <- function(N, k, p, q, approximation) {
  ne <- N - k
  s <- min(p, q)
  out <- list(error_df = ne, s = s, q = q, df1 = p * q)
  if (identical(approximation, "pillai_v")) {
    out$df2 <- s * (ne - p + s)
    out$g <- s
  } else if (identical(approximation, "hotelling_lawley_t1")) {
    out$df2 <- s * (ne - p - 1) + 2
    out$g <- s
  } else if (identical(approximation, "wilks_u")) {
    t <- .ig_manova_rao_t(p, q)
    out$rao_t <- t
    out$df2 <- t * (ne - (p - q + 1) / 2) - (p * q - 2) / 2
    out$g <- t
  } else if (identical(approximation, "hotelling_lawley_t2")) {
    b <- ((ne + q - p - 1) * (ne - 1)) / ((ne - p - 3) * (ne - p))
    out$mckeon_b <- b
    out$df2 <- 4 + (p * q + 2) / (b - 1)
    out$g <- (out$df2 - 2) / (ne - p - 1)
  } else {
    stop("unknown MANOVA approximation: ", approximation, call. = FALSE)
  }
  out
}

# Rao's (1951) exponent t of Wilks U (often written s in Rao's formula).
.ig_manova_rao_t <- function(p, q) {
  d <- p^2 + q^2 - 5
  if (d > 0) sqrt((p^2 * q^2 - 4) / d) else 1
}

# Reason why the approximation is undefined for a design, or NULL.
.ig_manova_problem <- function(tm, approximation, algorithm) {
  name <- .ig_manova_criterion_labels[[approximation]]
  if (!(tm$error_df >= 1)) {
    return(sprintf("the total sample size must exceed the number of groups (error df N - k = %s)",
                   .ig_format_double(tm$error_df)))
  }
  if (!is.finite(tm$df2) || tm$df2 <= 0) {
    return(sprintf(paste("the %s F approximation is undefined for this design: its denominator df is %s",
                         "(error df %s, %s response variable(s), hypothesis df %s); increase N"),
                   name, .ig_format_double(tm$df2), .ig_format_double(tm$error_df),
                   .ig_format_double(tm$df1 / tm$q), .ig_format_double(tm$q)))
  }
  if (identical(algorithm, "obrien_shieh_1999") && (!is.finite(tm$g) || tm$g <= 0)) {
    return(sprintf(paste("the O'Brien-Shieh noncentrality N * g * f2 is undefined for this design: the %s",
                         "multiplier g = (df2 - 2) / (nu_e - p - 1) is %s (error df %s); increase N or use",
                         "the Muller-Peterson algorithm"),
                   name, .ig_format_double(tm$g), .ig_format_double(tm$error_df)))
  }
  NULL
}

.ig_manova_design_problem <- function(p, N, effect) {
  k <- as.double(p$n_groups)
  if (!(N > k)) {
    return(sprintf("the total sample size (%s) must be greater than the number of groups (%s)",
                   .ig_format_double(N), .ig_format_double(k)))
  }
  tm <- .ig_manova_terms(N, k, as.double(p$n_responses), .ig_manova_q(p, effect), p$approximation)
  .ig_manova_problem(tm, p$approximation, p$algorithm)
}

# Reason why N is not a candidate of the a priori search, or NULL. On top of
# .ig_manova_design_problem(): McKeon's approximation of Hotelling-Lawley T2
# is derived for n = (nu_e - p - 1) / 2 > 0, i.e. nu_e >= p + 2, where its
# scale c = (df2 - 2) / (nu_e - p - 1) is positive and finite. For nu_e <= p + 1
# the formula still gives df2 > 0 for many designs (c is 0/0 at nu_e = p + 1
# and negative for nu_e <= p), and with Muller-Peterson noncentrality its
# power at some nu_e < p exceeds that of all designs up to a much larger N
# (special effects, k = 3, q = 5, p = 8, f2 = 2: 0.36 at N = 10, below 0.3
# from N = 11 to 18, 0.36 again at N = 19). Such a design is no solution of
# the a priori problem; G*Power returns N = 19 there (REV-G4-MANS-APR-S17;
# global effects: REV-G4-MANG-APR-T12). O'Brien-Shieh admits exactly the
# same designs through 0 < g < Inf. Fixed N (post hoc and the other analysis
# types) is not affected: there the formula is evaluated as G*Power does.
.ig_manova_a_priori_problem <- function(p, N, effect) {
  pr <- .ig_manova_design_problem(p, N, effect)
  if (!is.null(pr)) return(pr)
  if (identical(p$approximation, "hotelling_lawley_t2")) {
    ne <- N - as.double(p$n_groups)
    r <- as.double(p$n_responses)
    if (!(ne >= r + 2)) {
      return(sprintf(paste("error df %s is below p + 2 = %s: outside the domain of McKeon's approximation of",
                           "Hotelling-Lawley T2 (its scale c is not positive there)"),
                     .ig_format_double(ne), .ig_format_double(r + 2)))
    }
  }
  NULL
}

.ig_constraints_f_manova <- function(p, effect) {
  out <- character()
  if (!is.null(p$power) && !is.null(p$alpha) && p$power < p$alpha) {
    out <- c(out, sprintf("the target power (%s) must not be below alpha (%s)",
                          .ig_format_double(p$power), .ig_format_double(p$alpha)))
  }
  if (!is.null(p$n_total)) {
    pr <- .ig_manova_design_problem(p, as.double(p$n_total), effect)
    if (!is.null(pr)) out <- c(out, pr)
  }
  out
}

# Allocation rule. Global effects: equal groups of N / k (N a multiple of k).
# Special effects: N alone (G*Power uses no group sizes; the error df is
# N - k). NULL when N is not an admissible a priori design (the rule is used
# by the a priori solver only). `continuous = TRUE` gives the fractional
# relaxation used only to start the a priori search.
.ig_design_f_manova <- function(p, n_total, effect, continuous = FALSE) {
  k <- as.double(p$n_groups)
  global <- identical(effect, "global")
  if (continuous) {
    return(if (global) stats::setNames(rep(n_total / k, k), paste0("n", seq_len(k))) else c(n_total = n_total))
  }
  if (global && (n_total %% k) != 0) return(NULL)
  if (!is.null(.ig_manova_a_priori_problem(p, n_total, effect))) return(NULL)
  if (global) stats::setNames(rep(n_total / k, k), paste0("n", seq_len(k))) else c(n_total = n_total)
}

# Upper-alpha quantile of the central F distribution. stats::qf() is used
# where it is accurate; for df2 > 4e5 (or df1 > 4e5) nmath replaces the F
# quantile by a chi-square quantile, so the root of the upper tail of
# stats::pf() (incomplete beta, no such switch) is taken instead.
.ig_manova_q_f_upper <- function(alpha, df1, df2) {
  x0 <- stats::qf(alpha, df1, df2, lower.tail = FALSE)
  if (df1 <= 4e5 && df2 <= 4e5) return(x0)
  .ig_invert_cdf(function(x, lt) stats::pf(x, df1, df2, lower.tail = lt, log.p = TRUE),
                 alpha, FALSE, FALSE, support = c(0, Inf),
                 start = if (is.finite(x0) && x0 > 0) x0 else 1)
}

# Rejection and acceptance probabilities under H1, P(F' > crit) and
# P(F' <= crit) for F' ~ F(df1, df2, ncp). The smaller tail is computed
# (Poisson mixture of incomplete beta tails, .ig_p_f) and the larger one is its
# complement, so both stay in [0, 1], add up to one and the small one keeps its
# relative accuracy (the mixture of the larger tail can exceed one by a few
# units in the twelfth digit). Which tail is smaller is decided by the bulk of
# F' (mean of the noncentral chi-square numerator over df1); when the
# acceptance probability is below the double range by a rigorous Chernoff
# bound (.ig_manova_log_beta_bound) it is 0 and the power is 1 without summing
# a mixture of up to millions of terms (sensitivity scans reach lambda ~ 1e10).
.ig_manova_tails <- function(crit, df1, df2, ncp) {
  if (.ig_manova_log_beta_bound(crit, df1, df2, ncp) < -745.2) {
    return(list(power = 1, beta = 0, beta_bound = TRUE))
  }
  if (crit * df1 >= df1 + ncp) {
    power <- .ig_p_f(crit, df1, df2, ncp, lower_tail = FALSE)
    beta <- 1 - power
    if (power > 0.5) {
      beta <- .ig_p_f(crit, df1, df2, ncp, lower_tail = TRUE)
      power <- 1 - beta
    }
  } else {
    beta <- .ig_p_f(crit, df1, df2, ncp, lower_tail = TRUE)
    power <- 1 - beta
    if (beta > 0.5) {
      power <- .ig_p_f(crit, df1, df2, ncp, lower_tail = FALSE)
      beta <- 1 - power
    }
  }
  list(power = min(max(power, 0), 1), beta = min(max(beta, 0), 1), beta_bound = FALSE)
}

# Rigorous upper bound of log P(F' <= crit), F' = (X1 / df1) / (X2 / df2),
# X1 ~ chi2(df1, ncp), X2 ~ chi2(df2). For any d > 0,
#   P(F' <= c) <= P(X1 <= c df1 (1 + d)) + P(X2 >= df2 (1 + d)),
# with the Chernoff bounds
#   log P(X2 >= n (1 + d)) <= -(n / 2) (d - log(1 + d)),
#   log P(X1 <= a) <= theta a - ncp theta / (1 + 2 theta) - (df1 / 2) log(1 + 2 theta),
# minimised at 1 + 2 theta = (df1 + sqrt(df1^2 + 4 a ncp)) / (2 a) (trivial,
# 0, when that is not above 1, i.e. a >= df1 + ncp). d is chosen so that the
# X2 term is at most exp(-750). Used only to recognise acceptance probabilities below
# the double range; returns 0 (no information) otherwise.
.ig_manova_log_beta_bound <- function(crit, df1, df2, ncp) {
  .ig_f_lower_log_bound(crit, df1, df2, ncp)
}

# Criterion value implied by f2 under the chosen approximation.
.ig_manova_criterion <- function(f2, tm, approximation) {
  switch(approximation,
         pillai_v = tm$s * f2 / (1 + f2),
         hotelling_lawley_t1 = tm$s * f2,
         wilks_u = (1 + f2)^(-tm$rao_t),
         hotelling_lawley_t2 = if (is.finite(tm$g) && tm$g > 0) tm$g * f2 else NA_real_)
}

# Forward model.
.ig_power_f_manova <- function(p, effect) {
  N <- as.double(p$n_total)
  k <- as.double(p$n_groups)
  r <- as.double(p$n_responses)
  q <- .ig_manova_q(p, effect)
  approx <- p$approximation
  alg <- p$algorithm
  f2 <- as.double(p$f2)
  tm <- .ig_manova_terms(N, k, r, q, approx)
  problem <- .ig_manova_problem(tm, approx, alg)
  if (!is.null(problem)) .ig_abort("invalid_input", problem)
  df1 <- tm$df1
  df2 <- tm$df2
  ncp <- if (identical(alg, "obrien_shieh_1999")) N * tm$g * f2 else df2 * f2
  crit <- .ig_manova_q_f_upper(p$alpha, df1, df2)
  tails <- .ig_manova_tails(crit, df1, df2, ncp)
  power <- tails$power
  beta <- tails$beta
  cv <- .ig_manova_criterion(f2, tm, approx)
  outputs <- list(
    criterion = approx,
    criterion_label = .ig_manova_criterion_labels[[approx]],
    criterion_value = cv,
    f2 = f2,
    s = tm$s,
    hypothesis_df = q,
    error_df = tm$error_df,
    noncentrality_multiplier = if (identical(alg, "obrien_shieh_1999")) N * tm$g else df2,
    approximation_g = if (is.finite(tm$g)) tm$g else NA_real_,
    group_size = N / k
  )
  if (!is.null(tm$rao_t)) outputs$rao_t <- tm$rao_t
  if (!is.null(tm$mckeon_b)) outputs$mckeon_b <- if (is.finite(tm$mckeon_b)) tm$mckeon_b else NA_real_
  notes <- list()
  if (is.na(cv)) {
    # G*Power's display there: '-nan(ind)' for c = 0/0 (error df p + 1,
    # GPW-F-MANG-POS-B006-E01) and the product c f2 for a negative c (error
    # df <= p: -4.2 in REV-G4-MANG-POS-R16, -42 in REV-G4-MANS-POS-S15).
    outputs$criterion_value_note <- sprintf(paste(
      "Hotelling-Lawley T2 = c f2 is undefined: McKeon's scale c = (df2 - 2) / (nu_e - p - 1) is %s at",
      "error df %s, outside the domain nu_e >= p + 2 of the approximation (G*Power displays %s); the",
      "Muller-Peterson power does not use c"),
      .ig_format_double(tm$g), .ig_format_double(tm$error_df),
      if (is.finite(tm$g)) sprintf("the product c f2 = %s, which is not the value of a trace",
                                   .ig_format_double(signif(tm$g * f2, 7))) else "'-nan(ind)'")
  }
  if (identical(effect, "global") && (N %% k) != 0) {
    outputs$group_size_note <- paste(
      "N is not a multiple of k: the average group size N / k is reported; the error df N - k and the",
      "noncentrality use N (G*Power notes 'The analysis was based on the average group size')")
  }
  assumptions <- c(
    "p response variables, multivariate normal with a covariance matrix common to all groups",
    sprintf("%s F approximation: df1 = p q = %s, df2 = %s; q = %s", outputs$criterion_label,
            .ig_format_double(df1), .ig_format_double(signif(df2, 15)),
            if (identical(effect, "global")) "k - 1" else "number of predictors"),
    if (identical(alg, "obrien_shieh_1999")) {
      "noncentrality (O'Brien and Shieh 1999): lambda = N g f2 with g the leading coefficient of df2 in the error df"
    } else {
      "noncentrality (Muller and Peterson 1984): lambda = df2 f2"
    },
    "one-sided: H0 is rejected for large F (upper tail), as in G*Power",
    if (tm$error_df < r) {
      sprintf(paste("error df %s is below the number of response variables %s: the error SSCP matrix of a",
                    "sample is singular; the approximation is used outside the range it was derived for"),
              .ig_format_double(tm$error_df), .ig_format_double(r))
    }
  )
  list(
    power = power,
    beta = beta,
    alpha_attained = p$alpha,
    critical = crit,
    df = c(df1, df2),
    ncp = ncp,
    n_total = N,
    statistic_distribution = list(
      h0 = list(family = "F", df1 = df1, df2 = df2, ncp = 0),
      h1 = list(family = "F", df1 = df1, df2 = df2, ncp = ncp),
      rejection_side = "upper"
    ),
    outputs = outputs,
    notes = notes,
    assumptions = assumptions
  )
}

# A priori search. Global effects search the per-group size m = N / k (so
# that every multiple of k is a candidate, whatever k is); special effects
# search N. The search starts at the smallest admissible design (N > k with
# the approximation defined; Hotelling-Lawley T2: error df >= p + 2, see
# .ig_manova_a_priori_problem()) and uses the generic integer search of
# R/num-root.R (bracketing, bisection, minimality evidence).
.ig_a_priori_f_manova <- function(p, n_lo, n_max, effect) {
  global <- identical(effect, "global")
  k <- as.double(p$n_groups)
  unit <- if (global) k else 1
  at <- function(N) {
    pp <- p
    pp$n_total <- N
    pp
  }
  valid <- function(N) is.null(.ig_manova_a_priori_problem(p, N, effect))
  # First admissible design: every approximation is defined once
  # nu_e > p + 3 (b > 1 for McKeon's T2); scan up to there.
  m <- ceiling(max(n_lo, k + 1) / unit)
  m_scan_end <- ceiling((k + as.double(p$n_responses) + 8) / unit) + 1
  first <- NA_real_
  while (m <= max(m_scan_end, m)) {
    if (valid(m * unit)) {
      first <- m
      break
    }
    if (m >= m_scan_end) break
    m <- m + 1
  }
  m_max <- floor(n_max / unit)
  if (is.na(first) || first > m_max) {
    return(list(status = "infeasible", reason = sprintf(
      "no admissible design with N <= %s: the %s approximation stays undefined", format(n_max, scientific = FALSE),
      .ig_manova_criterion_labels[[p$approximation]]),
      diagnostics = list(n_upper_limit = n_max)))
  }
  power_at <- function(mm) .ig_power_f_manova(at(mm * unit), effect)$power
  # Continuous relaxation (fractional N; diagnostic and start value only).
  relax <- tryCatch({
    f <- function(u) .ig_power_f_manova(at(exp(u)), effect)$power - p$power
    lo <- log(first * unit)
    if (f(lo) >= 0) NA_real_ else {
      br <- .ig_bracket_monotone(f, lo, lo, log(n_max), decreasing = FALSE)
      if (isTRUE(br$found)) {
        r <- if (isTRUE(br$exact)) list(x = br$a, kind = "crossing") else .ig_refine_root(f, br$a, br$b, br$fa, br$fb)
        if (r$kind %in% c("crossing", "steep_crossing")) exp(r$x) else NA_real_
      } else NA_real_
    }
  }, inferguidR_error = function(e) NA_real_, error = function(e) NA_real_)
  res <- .ig_integer_search(power_at = power_at, admissible = function(mm) valid(mm * unit),
                            n_lo = first, n_max = m_max, target = p$power, discrete = FALSE,
                            n_start = if (is.finite(relax)) relax / unit else NA_real_)
  d <- res$diagnostics
  scale_n <- function(x) if (is.null(x)) NULL else x * unit
  for (nm in c("n_lower", "n_upper_limit", "bracket", "enumerated_range", "predecessor_n_total")) {
    if (!is.null(d[[nm]])) d[[nm]] <- scale_n(d[[nm]])
  }
  if (!is.null(d$power_sequence)) d$power_sequence$n_total <- scale_n(d$power_sequence$n_total)
  if (is.finite(relax)) {
    d$continuous_relaxation <- list(
      n_total = relax,
      note = "root of the power curve with fractional N; not an admissible design, used only to start the integer search")
  }
  d$first_admissible_n_total <- first * unit
  d$search_unit <- if (global) {
    sprintf("multiples of k = %s (equal group sizes N / k)", .ig_format_double(k))
  } else "N"
  d$start_rule <- paste("the search starts at the smallest N > k at which the approximation is defined",
                        "(Hotelling-Lawley T2: error df >= p + 2) and returns the smallest admissible design",
                        "reaching the target; G*Power 3.1.9.7 reports 'Internal math error' for a number of",
                        "inputs whose solution lies near the smallest designs")
  res$diagnostics <- d
  if (identical(res$status, "ok")) res$n <- res$n * unit
  if (!identical(res$status, "ok") && .ig_is_string(res$reason) && global) {
    res$reason <- paste0(.ig_manova_scale_reason(res$reason, unit),
                         sprintf(" (searched over multiples of k = %s)", .ig_format_double(k)))
  }
  res
}

# The generic integer search reports sizes in its own unit; for global effects
# that is m = N / k. Rewrite "N <= m", "N = m" and "m <= N" as totals.
.ig_manova_scale_reason <- function(txt, unit) {
  rescale <- function(s) {
    num <- regmatches(s, regexpr("[0-9]+", s))
    sub(num, format(as.numeric(num) * unit, scientific = FALSE), s, fixed = TRUE)
  }
  for (pattern in c("N (<=|=) [0-9]+", "[0-9]+ <= N\\b")) {
    mt <- gregexpr(pattern, txt, perl = TRUE)
    regmatches(txt, mt) <- lapply(regmatches(txt, mt), function(v) vapply(v, rescale, "", USE.NAMES = FALSE))
  }
  txt
}
