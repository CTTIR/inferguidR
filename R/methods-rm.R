# Method definitions: repeated-measures ANOVA and MANOVA (group G3) ----------------------
#
# Ownership: group G3 (repeated measures). Six G*Power procedures of the F
# family, registered through the provider .ig_define_methods_F_g3_rm():
#
#   f_rm_anova_between          ANOVA: Repeated measures, between factors
#   f_rm_anova_within           ANOVA: Repeated measures, within factors
#   f_rm_anova_within_between   ANOVA: Repeated measures, within-between interaction
#   f_rm_manova_between         MANOVA: Repeated measures, between factors
#   f_rm_manova_within          MANOVA: Repeated measures, within factors
#   f_rm_manova_within_between  MANOVA: Repeated measures, within-between interaction
#
# The G*Power manual has no chapter for these procedures. The mathematical
# contracts were written from the certified live inventory and the
# emulator captures (G*Power 3.1.9.7 under Wine) and are stated, with the
# evidence for every formula, in specs/method-notes/<method_id>.md.
#
# Design: k groups of n subjects each (N = k n), m repeated measurements per
# subject, compound-symmetric covariance (common SD, correlation rho between
# any two measurements); for the within-subject tests the nonsphericity
# correction epsilon multiplies numerator df, denominator df and the
# noncentrality as in the Geisser-Greenhouse correction.
#
# RM-ANOVA (univariate approach). With f the effect size of the selected
# convention ("Effect size specification" option):
#
#   between:         df1 = k - 1,              df2 = N - k
#   within:          df1 = (m - 1) eps,        df2 = (N - k)(m - 1) eps
#   within-between:  df1 = (k - 1)(m - 1) eps, df2 = (N - k)(m - 1) eps
#
#   lambda, "as in GPower 3.0" (f with the correlation as an input):
#     between  f^2 N m / (1 + (m - 1) rho)
#     within / within-between  f^2 N m eps / (1 - rho)
#   lambda, "with implicit rho" and "Cohen (1988)" (f(V)): f^2 N [eps]
#   lambda, "as in SPSS" (f(U), partial eta^2 / (1 - partial eta^2)): f^2 df2
#
# RM-MANOVA (multivariate approach): the general MANOVA F approximations with
# p response contrasts, q hypothesis df and error df nu_e = N - k:
#
#   between:         p = 1,     q = k - 1   (test of the subject means)
#   within:          p = m - 1, q = 1       (test of the m - 1 contrasts)
#   within-between:  p = m - 1, q = k - 1
#
# s = min(p, q). With "Use mean correlation in effect size calculation"
# checked (between and within only) the entered f is converted with the
# correlation: f^2 m / (1 + (m - 1) rho) (between) and f^2 m / (1 - rho)
# (within); unchecked, f is the MANOVA effect size of the approximation
# (f(V), f(U), f(T1), f(T2)) itself. Approximations (df1 = p q):
#
#   Pillai V (Pillai & Mijares 1959):   df2 = s (nu_e - p + s);  V = s f^2 / (1 + f^2)
#   Wilks U (Rao 1951):                 df2 = g (nu_e - (p - q + 1)/2) - (p q - 2)/2,
#                                       g = sqrt((p^2 q^2 - 4) / (p^2 + q^2 - 5)) (1 when p^2 + q^2 <= 5);
#                                       U = (1 + f^2)^(-g)
#   Hotelling-Lawley T1 (Pillai & Samson 1959): df2 = s (nu_e - p - 1) + 2;  T1 = s f^2
#   Hotelling-Lawley T2 (McKeon 1974):  df2 = 4 + (p q + 2) / (B - 1),
#                                       B = (nu_e + q - p - 1)(nu_e - 1) / ((nu_e - p - 3)(nu_e - p));
#                                       T2 = f^2 (df2 - 2) / (nu_e - p - 1);
#                                       admitted for nu_e >= p + 2 (df2 = 4 at nu_e = p + 3,
#                                       the formula's continuation at nu_e = p + 2), as G*Power
#
#   lambda, O'Brien & Shieh (1999):  N c f^2 with c = s (Pillai, T1), g (Wilks),
#                                    (df2 - 2)/(nu_e - p - 1) (T2)
#   lambda, Muller & Peterson (1984): df2 f^2
#
# A priori designs: equal group sizes, N a multiple of k (n >= 2 per group),
# for the three RM-ANOVA tests and the RM-MANOVA between/within tests; any N
# with valid df for the RM-MANOVA within-between test (G*Power's own MANOVA
# rule; evidence in the method notes).
#
# The power is that of the upper-tail F test: critical value
# F_crit = F^-1(1 - alpha; df1, df2), power = P(F'(df1, df2, lambda) > F_crit),
# beta computed directly as the lower tail (R/num-dist.R, .ig_p_f()).

.ig_define_methods_F_g3_rm <- function() {
  list(.ig_method_rm_anova("between"), .ig_method_rm_anova("within"),
       .ig_method_rm_anova("within_between"), .ig_method_rm_manova("between"),
       .ig_method_rm_manova("within"), .ig_method_rm_manova("within_between"))
}

# Shared pieces -----------------------------------------------------------------------

.ig_rm_ids <- list(
  anova = c(between = "f_rm_anova_between", within = "f_rm_anova_within",
            within_between = "f_rm_anova_within_between"),
  manova = c(between = "f_rm_manova_between", within = "f_rm_manova_within",
             within_between = "f_rm_manova_within_between")
)

.ig_rm_labels <- list(
  anova = c(between = "ANOVA: Repeated measures, between factors",
            within = "ANOVA: Repeated measures, within factors",
            within_between = "ANOVA: Repeated measures, within-between interaction"),
  manova = c(between = "MANOVA: Repeated measures, between factors",
             within = "MANOVA: Repeated measures, within factors",
             within_between = "MANOVA: Repeated measures, within-between interaction")
)

.ig_rm_labels_de <- list(
  anova = c(between = "ANOVA: Messwiederholung, Zwischensubjektfaktoren",
            within = "ANOVA: Messwiederholung, Innersubjektfaktoren",
            within_between = "ANOVA: Messwiederholung, Wechselwirkung inner x zwischen"),
  manova = c(between = "MANOVA: Messwiederholung, Zwischensubjektfaktoren",
             within = "MANOVA: Messwiederholung, Innersubjektfaktoren",
             within_between = "MANOVA: Messwiederholung, Wechselwirkung inner x zwischen")
)

.ig_rm_paths <- list(
  anova = c(between = "Tests > Means > Repeated measures: Between factors, ANOVA-approach",
            within = "Tests > Means > Repeated measures: Within factors, ANOVA-approach",
            within_between = "Tests > Means > Repeated measures: Within-between interaction, ANOVA-approach"),
  manova = c(between = "Tests > Means > Repeated measures: Between factors, MANOVA-approach",
             within = "Tests > Means > Repeated measures: Within factors, MANOVA-approach",
             within_between = "Tests > Means > Repeated measures: Within-between interaction, MANOVA-approach")
)

.ig_rm_spec_labels <- c(gpower30 = "as in GPower 3.0",
                        gpower30_implicit_rho = "as in GPower 3.0 with implicit rho",
                        spss = "as in SPSS", cohen1988 = "as in Cohen (1988) - recommended")

.ig_rm_approx_labels <- c(wilks_u = "Wilks U (Rao, 1951)",
                          hotelling_lawley_t1 = "Hotelling-Lawley T1 (Pillai and Samson, 1959)",
                          hotelling_lawley_t2 = "Hotelling-Lawley T2 (McKeon, 1974)",
                          pillai_v = "Pillai V (Pillai and Mijares, 1959)")

.ig_rm_criterion_labels <- c(wilks_u = "Wilks U", hotelling_lawley_t1 = "Hotelling-Lawley T1",
                             hotelling_lawley_t2 = "Hotelling-Lawley T2", pillai_v = "Pillai V")

.ig_rm_effect_symbol <- c(wilks_u = "f(U)", hotelling_lawley_t1 = "f(T1)",
                          hotelling_lawley_t2 = "f(T2)", pillai_v = "f(V)")

# Upper critical value of the central F distribution.
#
# stats::qf() replaces F(df1, df2) by chi-square(df1) / df1 when df2 > 4e5
# (nmath qf.c), which moves the critical value in the 5th to 7th digit at the
# denominator df of large designs (for example 1.556574794 instead of
# 1.556658706 at df = (81, 436735.3), alpha = 0.001; G*Power displays
# 1.5566587). The value is therefore taken from the beta quantile of the
# smaller tail and polished by Newton steps on stats::pf(), which is exact
# for finite df (pbeta, TOMS 708).
.ig_rm_qf_upper <- function(alpha, df1, df2) {
  # The smaller tail is solved: the upper tail for alpha <= 0.5, the lower
  # tail (probability 1 - alpha, exact in double precision for alpha >= 0.5)
  # otherwise, so the target keeps its relative accuracy at both ends.
  upper <- alpha <= 0.5
  target <- if (upper) alpha else 1 - alpha
  x <- stats::qbeta(target, df1 / 2, df2 / 2, lower.tail = !upper)
  crit <- if (is.finite(x) && x <= 0.5) {
    df2 / df1 * x / (1 - x)
  } else {
    y <- stats::qbeta(target, df2 / 2, df1 / 2, lower.tail = upper)
    df2 / df1 * (1 - y) / y
  }
  if (!is.finite(crit) || crit <= 0) return(crit)
  lt <- log(target)
  for (i in seq_len(8L)) {
    lp <- stats::pf(crit, df1, df2, lower.tail = !upper, log.p = TRUE)
    ld <- stats::df(crit, df1, df2, log = TRUE)
    # d log P / d crit = -+ density / P (upper / lower tail).
    step <- (lp - lt) * exp(lp - ld) * (if (upper) 1 else -1)
    if (!is.finite(step)) break
    new <- crit + step
    if (new <= 0) new <- crit / 2
    done <- abs(new - crit) <= 4 * .Machine$double.eps * crit
    crit <- new
    if (done) break
  }
  crit
}

# Shared Chernoff bound for beta, including fractional positive df1.
.ig_rm_beta_log_bound <- function(crit, df1, df2, ncp) {
  .ig_f_lower_log_bound(crit, df1, df2, ncp)
}

# Upper-tail F test: critical value, power and beta (computed directly).
# When beta is provably below the double range (huge noncentrality, for
# example during a sensitivity scan), power = 1 and beta = 0 are returned
# without summing a Poisson mixture of billions of terms.
.ig_rm_f_test <- function(df1, df2, ncp, alpha) {
  crit <- .ig_rm_qf_upper(alpha, df1, df2)
  if (!is.finite(crit)) {
    stop(sprintf("no finite critical F for alpha = %s, df = (%s, %s)", format(alpha), format(df1),
                 format(df2)), call. = FALSE)
  }
  if (.ig_rm_beta_log_bound(crit, df1, df2, ncp) < -745.2) {
    return(list(critical = crit, power = 1, beta = 0))
  }
  # The smaller tail is computed directly and the larger one as its
  # complement: the beta argument of the upper tail, df2 / (df1 F + df2),
  # carries no relative accuracy in its complement when F_crit is tiny
  # (alpha near 1), which made the power noisy at the 1e-9 level there.
  lower <- .ig_p_f(crit, df1, df2, ncp, lower_tail = TRUE)
  if (is.finite(lower) && lower <= 0.5) {
    return(list(critical = crit, power = 1 - lower, beta = lower))
  }
  upper <- .ig_p_f(crit, df1, df2, ncp, lower_tail = FALSE)
  list(critical = crit, power = upper, beta = 1 - upper)
}

.ig_rm_is_whole <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x == round(x)

# Equal-group design: N must be a multiple of k with at least `n_min` per group.
.ig_rm_design_equal <- function(k, n_total, n_min = 2, continuous = FALSE) {
  k <- as.double(k)
  if (continuous) return(stats::setNames(rep(n_total / k, k), paste0("n", seq_len(k))))
  if (n_total %% k != 0 || n_total / k < n_min) return(NULL)
  stats::setNames(rep(n_total / k, k), paste0("n", seq_len(k)))
}

.ig_rm_common_params <- function(effect, k_lower, m_lower, f_label, f_description, with_epsilon,
                                 rho_applies = NULL) {
  ps <- list(
    ig_param("f", role = "effect", domain = ig_domain_real(0, Inf, lower_inclusive = FALSE),
             label = list(en = f_label, de = sub("Effect size", "Effektst\u00e4rke", f_label)),
             description = f_description, aliases = "effect_size_f"),
    ig_param_alpha(),
    ig_param_power(),
    ig_param_error_ratio(),
    ig_param("n_total", role = "design", domain = ig_domain_integer(lower = 2), unit = "subjects",
             label = list(en = "Total sample size", de = "Gesamtstichprobenumfang"),
             description = "Total number of subjects N (k groups); solved in a priori mode."),
    ig_param("k", role = "design", domain = ig_domain_integer(lower = k_lower), unit = "groups",
             label = list(en = "Number of groups", de = "Anzahl der Gruppen"),
             description = "Number of independent groups (levels of the between-subjects factor).",
             aliases = "number_of_groups"),
    ig_param("m", role = "design", domain = ig_domain_integer(lower = m_lower), unit = "measurements",
             label = list(en = "Number of measurements", de = "Anzahl der Messungen"),
             description = "Number of repeated measurements per subject (levels of the within-subjects factor).",
             aliases = "number_of_measurements"),
    ig_param("rho", role = "nuisance",
             domain = ig_domain_real(-1, 1, lower_inclusive = FALSE, upper_inclusive = FALSE),
             label = list(en = "Corr among rep measures", de = "Korrelation der Messwiederholungen"),
             description = paste("Common correlation between any two of the m measurements (compound",
                                 "symmetry). An input only where the effect-size convention uses it."),
             aliases = "corr_among_rep_measures", applies_when = rho_applies)
  )
  if (with_epsilon) {
    ps[[length(ps) + 1L]] <- ig_param(
      "epsilon", role = "nuisance",
      domain = ig_domain_real(0, 1, lower_inclusive = FALSE, upper_inclusive = TRUE), default = 1,
      label = list(en = "Nonsphericity correction \u03b5", de = "Nonsph\u00e4rizit\u00e4tskorrektur \u03b5"),
      description = paste("Nonsphericity correction epsilon in [1/(m - 1), 1] (Greenhouse-Geisser type):",
                          "multiplies both df and the noncentrality. The default 1 (sphericity) is",
                          "G*Power's initial field value (live capture 2026-10-05)."),
      aliases = "nonsphericity_epsilon")
  }
  ig_params(ps)
}

# Cross-parameter checks shared by both approaches.
.ig_rm_constraints_common <- function(p, effect) {
  out <- character()
  k <- p$k
  m <- p$m
  if (!is.null(p$n_total) && !is.null(k) && as.double(p$n_total) <= as.double(k)) {
    out <- c(out, sprintf("the total sample size N = %s must be greater than the number of groups k = %s",
                          format(p$n_total), format(k)))
  }
  if (!is.null(p$epsilon) && !is.null(m) && m >= 2 && p$epsilon < 1 / (as.double(m) - 1) * (1 - 1e-12)) {
    out <- c(out, sprintf(paste("the nonsphericity correction epsilon = %s is below its lower bound",
                                "1/(m - 1) = %s"), format(p$epsilon), format(1 / (m - 1), digits = 7)))
  }
  if (!is.null(p$rho) && !is.null(m)) {
    if (effect == "between") {
      if (1 + (as.double(m) - 1) * p$rho <= 0) {
        out <- c(out, sprintf(paste("the correlation rho = %s gives 1 + (m - 1) rho = %s <= 0 for m = %s",
                                    "(the variance of the subject mean must be positive); rho must exceed",
                                    "-1/(m - 1)"), format(p$rho), format(1 + (m - 1) * p$rho), format(m)))
      }
    } else if (p$rho >= 1) {
      out <- c(out, sprintf(paste("the correlation rho = %s leaves no within-subject error variance",
                                  "(1 - rho must be positive)"), format(p$rho)))
    }
  }
  out
}

# RM-ANOVA ------------------------------------------------------------------------------

.ig_rm_anova_df <- function(effect, n_total, k, m, epsilon) {
  N <- as.double(n_total)
  k <- as.double(k)
  m <- as.double(m)
  eps <- if (is.null(epsilon)) 1 else as.double(epsilon)
  switch(effect,
         between = c(df1 = k - 1, df2 = N - k),
         within = c(df1 = (m - 1) * eps, df2 = (N - k) * (m - 1) * eps),
         within_between = c(df1 = (k - 1) * (m - 1) * eps, df2 = (N - k) * (m - 1) * eps))
}

.ig_rm_anova_ncp <- function(effect, spec, f, n_total, m, rho, epsilon, df2) {
  f2 <- f^2
  N <- as.double(n_total)
  m <- as.double(m)
  eps <- if (is.null(epsilon)) 1 else as.double(epsilon)
  switch(spec,
         gpower30 = if (effect == "between") f2 * N * m / (1 + (m - 1) * rho) else f2 * N * m * eps / (1 - rho),
         gpower30_implicit_rho = ,
         cohen1988 = if (effect == "between") f2 * N else f2 * N * eps,
         spss = f2 * df2)
}

.ig_power_rm_anova <- function(effect) {
  force(effect)
  function(p) {
    spec <- p$effect_size_specification
    N <- as.double(p$n_total)
    k <- as.double(p$k)
    df <- .ig_rm_anova_df(effect, N, k, p$m, p$epsilon)
    ncp <- .ig_rm_anova_ncp(effect, spec, p$f, N, p$m, p$rho, p$epsilon, df[["df2"]])
    core <- .ig_rm_f_test(df[["df1"]], df[["df2"]], ncp, p$alpha)
    gs <- if (N %% k == 0) stats::setNames(rep(N / k, k), paste0("n", seq_len(k))) else NULL
    lam_txt <- switch(spec,
      gpower30 = if (effect == "between") "lambda = f^2 N m / (1 + (m - 1) rho)" else
        "lambda = f^2 N m epsilon / (1 - rho)",
      gpower30_implicit_rho = , cohen1988 = if (effect == "between") "lambda = f(V)^2 N" else
        "lambda = f(V)^2 N epsilon",
      spss = "lambda = f(U)^2 df2 (f(U)^2 = partial eta^2 / (1 - partial eta^2))")
    df_txt <- switch(effect,
      between = "df1 = k - 1, df2 = N - k",
      within = "df1 = (m - 1) epsilon, df2 = (N - k)(m - 1) epsilon",
      within_between = "df1 = (k - 1)(m - 1) epsilon, df2 = (N - k)(m - 1) epsilon")
    list(
      power = core$power, beta = core$beta, alpha_attained = p$alpha,
      critical = core$critical, df = df, ncp = ncp, n_total = N, group_sizes = gs,
      outputs = list(n_per_group = N / k),
      statistic_distribution = list(
        h0 = list(family = "F", df1 = df[["df1"]], df2 = df[["df2"]], ncp = 0),
        h1 = list(family = "F", df1 = df[["df1"]], df2 = df[["df2"]], ncp = ncp),
        rejection_side = "upper"),
      assumptions = c(
        sprintf("%s: k independent groups, m repeated measurements per subject, multivariate normal errors with compound symmetry (common variance, correlation rho)",
                .ig_rm_labels$anova[[effect]]),
        sprintf("effect size convention '%s'; %s; %s", .ig_rm_spec_labels[[spec]], lam_txt, df_txt),
        if (effect != "between") "the nonsphericity correction epsilon multiplies df1, df2 and lambda (Greenhouse-Geisser type correction, as G*Power)",
        "upper-tail F test: F_crit = F^-1(1 - alpha; df1, df2), power = P(F'(df1, df2, lambda) > F_crit)"
      )
    )
  }
}

.ig_method_rm_anova <- function(effect) {
  id <- .ig_rm_ids$anova[[effect]]
  with_eps <- effect != "between"
  rho_cond <- list(option = "effect_size_specification", "in" = "gpower30")
  f_desc <- paste(
    "Effect size in the convention selected by 'effect_size_specification': f with the correlation",
    "as a separate input ('as in GPower 3.0'; G*Power's main-window label 'Effect size f'), f(V)",
    "('with implicit rho' and 'Cohen (1988)'; label 'Effect size f(V)'), or f(U) = sqrt(partial",
    "eta^2 / (1 - partial eta^2)) ('as in SPSS'; label 'Effect size f(U)'). See the method notes.")
  ig_method(
    id = id, family = "F",
    label = list(en = .ig_rm_labels$anova[[effect]], de = .ig_rm_labels_de$anova[[effect]]),
    params = .ig_rm_common_params(effect, k_lower = if (effect == "within") 1 else 2, m_lower = 2,
                                  f_label = "Effect size f", f_description = f_desc,
                                  with_epsilon = with_eps, rho_applies = rho_cond),
    options = ig_options(
      ig_option("effect_size_specification", values = names(.ig_rm_spec_labels), default = "gpower30",
                labels = .ig_rm_spec_labels,
                label = list(en = "Effect size specification", de = "Spezifikation der Effektst\u00e4rke"))
    ),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = .ig_power_rm_anova(effect),
    design = function(p, n_total, continuous = FALSE) .ig_rm_design_equal(p$k, n_total, 2, continuous),
    discrete = FALSE,
    effect_direction = "increasing",
    engine = list(package = "stats", fun = "qf; noncentral F as a Poisson mixture of pbeta (R/num-dist.R .ig_p_f)",
                  version = NA_character_),
    checker = "validation/independent/g3/check_rm.R",
    reference = list(
      gpower_path = .ig_rm_paths$anova[[effect]],
      gpower_label = paste("F tests -", .ig_rm_labels$anova[[effect]]),
      manual_pages = "none (no manual chapter; formulas from the certified live inventory and emulator captures)",
      allocation_rule = paste("a priori: equal groups, N a multiple of k with at least 2 subjects per group;",
                              "the smallest such N reaching the target power (G*Power's rule, confirmed by",
                              "emulator captures; method notes)"),
      method_notes = sprintf("specs/method-notes/%s.md", id)
    ),
    status = "independently_checked",
    constraints = function(p) {
      out <- .ig_rm_constraints_common(p, effect)
      if (!is.null(p$power) && !is.null(p$alpha) && p$power < p$alpha) {
        out <- c(out, sprintf(paste("the target power %s is below alpha %s; the power of this test is at",
                                    "least alpha for every design"), format(p$power), format(p$alpha)))
      }
      out
    }
  )
}

# RM-MANOVA -----------------------------------------------------------------------------

.ig_rm_manova_dims <- function(effect, k, m) {
  k <- as.double(k)
  m <- as.double(m)
  switch(effect,
         between = c(p = 1, q = k - 1),
         within = c(p = m - 1, q = 1),
         within_between = c(p = m - 1, q = k - 1))
}

# Rao's g for Wilks' U (1 when p^2 + q^2 <= 5, where the exact F transform
# has g = 1).
.ig_rm_rao_g <- function(p, q) {
  den <- p^2 + q^2 - 5
  if (den <= 0) return(1)
  sqrt((p^2 * q^2 - 4) / den)
}

# The F approximation of one MANOVA criterion. Returns df1, df2, the
# multiplier c of the O'Brien-Shieh noncentrality (lambda = N c f^2), the
# population criterion value and Rao's g.
.ig_rm_manova_approx <- function(p, q, nu_e, f2, approximation) {
  s <- min(p, q)
  df1 <- p * q
  g <- NA_real_
  if (approximation == "pillai_v") {
    df2 <- s * (nu_e - p + s)
    mult <- s
    crit <- s * f2 / (1 + f2)
  } else if (approximation == "wilks_u") {
    g <- .ig_rm_rao_g(p, q)
    df2 <- g * (nu_e - (p - q + 1) / 2) - (p * q - 2) / 2
    mult <- g
    crit <- exp(-g * log1p(f2))
  } else if (approximation == "hotelling_lawley_t1") {
    df2 <- s * (nu_e - p - 1) + 2
    mult <- s
    crit <- s * f2
  } else {
    if (s == 1) {
      # McKeon's df2 and multiplier reduce exactly to these for s = 1.
      df2 <- if (q == 1) nu_e - p + 1 else nu_e
      mult <- 1
    } else {
      # df2 = 4 + (p q + 2) / (B - 1) with B - 1 written without cancellation:
      # B - 1 = ((p + q + 1) nu_e - (p^2 + 2p + q - 1)) / ((nu_e - p - 3)(nu_e - p)).
      # The numerator is positive for nu_e >= p + 2. At nu_e = p + 3 B is infinite and
      # df2 takes its limit 4; at nu_e = p + 2 (B < 0) the formula continues to
      # 2 < df2 < 4 with a positive multiplier. G*Power computes both designs
      # (REV-G3-RMMWB-POS-E02: df2 4.0000000; -E05: df2 2.8571429) and its a priori
      # search admits them (REV-G3-RMMWB-APR-E07, -E39). For nu_e <= p + 1 the
      # multiplier (df2 - 2) / (nu_e - p - 1) is undefined or not positive: no valid df.
      num <- (p + q + 1) * nu_e - (p^2 + 2 * p + q - 1)
      df2 <- if (nu_e - p - 2 >= 0 && num > 0) 4 + (df1 + 2) * (nu_e - p - 3) * (nu_e - p) / num else NA_real_
      mult <- (df2 - 2) / (nu_e - p - 1)
    }
    crit <- f2 * mult
  }
  list(df1 = df1, df2 = df2, mult = mult, criterion = crit, s = s, g = g)
}

.ig_rm_manova_f2 <- function(effect, p) {
  f2 <- p$f^2
  if (identical(p$use_mean_correlation, "on")) {
    m <- as.double(p$m)
    f2 <- if (effect == "between") f2 * m / (1 + (m - 1) * p$rho) else f2 * m / (1 - p$rho)
  }
  f2
}

.ig_rm_manova_eval <- function(effect, p) {
  N <- as.double(p$n_total)
  k <- as.double(p$k)
  dims <- .ig_rm_manova_dims(effect, k, p$m)
  f2 <- .ig_rm_manova_f2(effect, p)
  ap <- .ig_rm_manova_approx(dims[["p"]], dims[["q"]], N - k, f2, p$approximation)
  ncp <- if (identical(p$algorithm, "obrien_shieh_1999")) N * ap$mult * f2 else ap$df2 * f2
  c(ap, list(ncp = ncp, f2 = f2, p_dim = dims[["p"]], q_dim = dims[["q"]], nu_e = N - k))
}

# Smallest problems first: the df of the approximation must be positive and
# finite at this design.
.ig_rm_manova_df_ok <- function(ev) {
  is.finite(ev$df2) && ev$df2 > 0 && is.finite(ev$mult) && ev$mult > 0 && ev$df1 > 0
}

.ig_power_rm_manova <- function(effect) {
  force(effect)
  function(p) {
    ev <- .ig_rm_manova_eval(effect, p)
    if (!.ig_rm_manova_df_ok(ev)) {
      stop(sprintf("the %s approximation has no valid df at N = %s (df2 = %s)", p$approximation,
                   format(p$n_total), format(ev$df2)), call. = FALSE)
    }
    core <- .ig_rm_f_test(ev$df1, ev$df2, ev$ncp, p$alpha)
    N <- as.double(p$n_total)
    k <- as.double(p$k)
    gs <- if (N %% k == 0) stats::setNames(rep(N / k, k), paste0("n", seq_len(k))) else NULL
    outputs <- list(n_per_group = N / k, s = ev$s, p_response_contrasts = ev$p_dim,
                    q_hypothesis_df = ev$q_dim, f2_effect = ev$f2)
    outputs[[p$approximation]] <- ev$criterion
    if (!is.na(ev$g)) outputs$rao_g <- ev$g
    conv <- if (identical(p$use_mean_correlation, "on")) {
      if (effect == "between") "the entered f is converted with the correlation: f(V)^2 = f^2 m / (1 + (m - 1) rho)" else
        "the entered f is converted with the correlation: f(V)^2 = f^2 m / (1 - rho)"
    } else {
      sprintf("the entered effect size is the MANOVA effect size %s of the selected approximation",
              .ig_rm_effect_symbol[[p$approximation]])
    }
    list(
      power = core$power, beta = core$beta, alpha_attained = p$alpha,
      critical = core$critical, df = c(df1 = ev$df1, df2 = ev$df2), ncp = ev$ncp, n_total = N,
      group_sizes = gs, outputs = outputs,
      statistic_distribution = list(
        h0 = list(family = "F", df1 = ev$df1, df2 = ev$df2, ncp = 0),
        h1 = list(family = "F", df1 = ev$df1, df2 = ev$df2, ncp = ev$ncp),
        rejection_side = "upper"),
      assumptions = c(
        sprintf("%s: k independent groups, m repeated measurements per subject, multivariate normal errors; p = %s response contrasts, q = %s hypothesis df, error df N - k",
                .ig_rm_labels$manova[[effect]], format(ev$p_dim), format(ev$q_dim)),
        conv,
        sprintf("F approximation: %s; noncentrality: %s", .ig_rm_approx_labels[[p$approximation]],
                if (identical(p$algorithm, "obrien_shieh_1999")) "O'Brien and Shieh (1999), lambda = N c f^2" else
                  "Muller and Peterson (1984), lambda = df2 f^2"),
        if (identical(p$approximation, "hotelling_lawley_t2") && ev$s > 1 && ev$nu_e - ev$p_dim - 3 <= 0) {
          sprintf(paste("error df %s is not above p + 3 = %s: the Hotelling-Lawley trace has no finite variance there,",
                        "so McKeon's two-moment match is not defined; df2 = %s is the continuation of his formula",
                        "(its limit 4 at error df p + 3), as G*Power computes it"),
                  format(ev$nu_e), format(ev$p_dim + 3), format(ev$df2, digits = 7))
        },
        "upper-tail F test: F_crit = F^-1(1 - alpha; df1, df2), power = P(F'(df1, df2, lambda) > F_crit)"
      )
    )
  }
}

.ig_method_rm_manova <- function(effect) {
  id <- .ig_rm_ids$manova[[effect]]
  wb <- effect == "within_between"
  umc_values <- if (wb) "off" else c("on", "off")
  umc_default <- if (wb) "off" else "on"
  rho_cond <- list(option = "use_mean_correlation", "in" = "on")
  f_desc <- if (wb) {
    paste("MANOVA effect size of the selected approximation: f(V) (Pillai V, f^2 = V/(s - V)), f(U)",
          "(Wilks U, f^2 = U^(-1/g) - 1), f(T1) (Hotelling-Lawley T1, f^2 = T1/s) or f(T2)",
          "(Hotelling-Lawley T2); G*Power relabels the field accordingly.")
  } else {
    paste("With 'Use mean correlation' checked (default): the effect size f of the univariate",
          "convention, converted with the correlation rho (G*Power label 'Effect size f'); unchecked:",
          "the MANOVA effect size of the approximation, entered directly (label 'Effect size f(V)').")
  }
  ig_method(
    id = id, family = "F",
    label = list(en = .ig_rm_labels$manova[[effect]], de = .ig_rm_labels_de$manova[[effect]]),
    params = .ig_rm_common_params(effect, k_lower = if (effect == "within") 1 else 2, m_lower = 2,
                                  f_label = if (wb) "Effect size f(V)" else "Effect size f",
                                  f_description = f_desc, with_epsilon = FALSE, rho_applies = rho_cond),
    options = ig_options(
      ig_option("approximation", values = names(.ig_rm_approx_labels), default = "pillai_v",
                labels = .ig_rm_approx_labels,
                label = list(en = "Approximation (F-Transformation)", de = "Approximation (F-Transformation)")),
      ig_option("algorithm", values = c("obrien_shieh_1999", "muller_peterson_1984"),
                default = "obrien_shieh_1999",
                labels = c(obrien_shieh_1999 = "O'Brien and Shieh (1999), recommended",
                           muller_peterson_1984 = "Muller and Peterson (1984)"),
                label = list(en = "Algorithm", de = "Algorithmus")),
      ig_option("use_mean_correlation", values = umc_values, default = umc_default,
                labels = c(on = "checked", off = "unchecked")[umc_values],
                label = list(en = "Use mean correlation in effect size calculation",
                             de = "Mittlere Korrelation in der Effektst\u00e4rkeberechnung verwenden"))
    ),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = .ig_power_rm_manova(effect),
    design = function(p, n_total, continuous = FALSE) {
      if (wb) {
        if (!continuous) {
          ev <- .ig_rm_manova_eval(effect, c(p, list(n_total = n_total)))
          if (!.ig_rm_manova_df_ok(ev)) return(NULL)
        }
        return(c(n_total = n_total))
      }
      gs <- .ig_rm_design_equal(p$k, n_total, 2, continuous)
      if (is.null(gs) || continuous) return(gs)
      ev <- .ig_rm_manova_eval(effect, c(p, list(n_total = n_total)))
      if (!.ig_rm_manova_df_ok(ev)) return(NULL)
      gs
    },
    discrete = FALSE,
    effect_direction = "increasing",
    engine = list(package = "stats", fun = "qf; noncentral F as a Poisson mixture of pbeta (R/num-dist.R .ig_p_f)",
                  version = NA_character_),
    checker = "validation/independent/g3/check_rm.R",
    reference = list(
      gpower_path = .ig_rm_paths$manova[[effect]],
      gpower_label = paste("F tests -", .ig_rm_labels$manova[[effect]]),
      manual_pages = "none (no manual chapter; formulas from the certified live inventory and emulator captures)",
      allocation_rule = if (wb) {
        "a priori: the smallest N with valid df of the approximation (no group-size constraint; G*Power's MANOVA rule, emulator captures)"
      } else {
        "a priori: equal groups, N a multiple of k with at least 2 subjects per group and valid df; the smallest such N reaching the target power (emulator captures)"
      },
      method_notes = sprintf("specs/method-notes/%s.md", id)
    ),
    status = "independently_checked",
    constraints = function(p) {
      out <- .ig_rm_constraints_common(p, effect)
      if (!is.null(p$power) && !is.null(p$alpha) && p$power < p$alpha) {
        out <- c(out, sprintf(paste("the target power %s is below alpha %s; the power of this test is at",
                                    "least alpha for every design"), format(p$power), format(p$alpha)))
      }
      if (!length(out) && !is.null(p$n_total) && !is.null(p$k) && !is.null(p$m) && !is.null(p$approximation)) {
        q <- p
        if (is.null(q$f)) q$f <- 1
        if (identical(q$use_mean_correlation, "on") && is.null(q$rho)) q$rho <- 0
        ev <- .ig_rm_manova_eval(effect, q)
        if (!.ig_rm_manova_df_ok(ev)) {
          out <- c(out, sprintf(paste("N = %s is too small for the %s approximation with k = %s and m = %s",
                                      "(denominator df %s); a larger N is needed"),
                                format(p$n_total), .ig_rm_approx_labels[[p$approximation]], format(p$k),
                                format(p$m), format(ev$df2, digits = 7)))
        }
      }
      out
    }
  )
}
