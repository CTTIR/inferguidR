# Method definitions: fixed-effects ANOVA, ANCOVA and Hotelling T^2 (F family, group G2) -------
#
# Ownership: R/methods-anova.R belongs to group G2 (blueprint phases 06-10 for the fixed-effects
# ANOVA, ANCOVA and Hotelling T^2 procedures). Effect-size drawers: R/effect-anova.R. Method
# notes: specs/method-notes/<method_id>.md. Independent checks: validation/independent/g2/.
#
# All five procedures are upper-tailed F tests. Under H0 the statistic follows F(df1, df2); under
# H1 it follows the noncentral F(df1, df2, lambda):
#
#   method                             lambda                 df1      df2
#   f_anova_oneway                     f^2 N                  k - 1    N - k
#   f_anova_special_main_interactions  f^2 N                  u        N - k
#   f_ancova_main_interactions         f^2 N                  u        N - k - c
#   f_hotelling_t2_one_group           N Delta^2              p        N - p
#   f_hotelling_t2_two_groups          n1 n2 / N * Delta^2    p        N - p - 1
#
# (k groups or cells, u numerator df of the tested effect, c covariates, p response variables;
# manual sections 10.5 and 11.5 for the ANOVA rows; the ANCOVA and Hotelling rows have no manual
# chapter and were established from the emulator captures, see the method notes.)
#
# Numerics. The critical value solves P(F(df1, df2) > c) = alpha through the beta quantile
# (.ig_g2_crit_f()): stats::qf() replaces the F distribution by chi-square(df1) / df1 when
# df2 > 4e5, which moves the critical value by 2.3e-5 at df1 = 50, df2 = 1444522, alpha = 0.001
# (its upper tail is then 0.00100027 instead of 0.001; found in the Hotelling two-group capture
# GPW-F-HOT2-APR-B001-D01). The noncentral tails come from .ig_p_f() (R/num-dist.R: Poisson
# mixture of central pbeta terms of the requested tail, summed on the log scale). Both tails are
# summed; the smaller one is reported as summed, so that it keeps its relative accuracy, and the
# larger one as its complement (.ig_g2_f_test()). Where a rigorous bound shows beta below the double range the
# tails are exact 1 and 0 (.ig_g2_beta_negligible()); this keeps sensitivity scans over very large
# effects from summing millions of Poisson terms.
#
# A priori searches (method-specific `a_priori_search`, checked by the generic solver):
#
# * one-way ANOVA: smallest N that is a multiple of k with at least 2 subjects per group (searched
#   over the per-group size); DOC-055 and the captures show N = 390 for k = 10 although N = 386
#   already reaches the target.
# * special effects, ANCOVA, Hotelling one group: smallest integer N with df2 >= 1.
# * Hotelling two groups: G*Power's allocation rule found in the captures (method notes). For
#   allocation ratio r = 1 the designs are n1 = n2 = N / 2 (smallest even N reaching the target).
#   For r != 1, N* is the smallest N (any parity) whose power with fractional group sizes
#   N / (1 + r), r N / (1 + r) reaches the target, n1 = round(N* / (1 + r)) half up and
#   n2 = N* - n1. If that design has an empty group or misses the target with its integer group
#   sizes, inferguidR moves to the next N under the same split whose design has both groups and
#   reaches the target, and reports G*Power's design in the diagnostics (a documented difference:
#   G*Power reports such designs, captures GPW-F-HOT2-APR-B001-D01 and SUP-G2-F-HOT2-APR-A02).

.ig_define_methods_F_g2 <- function() {
  list(
    .ig_method_f_anova_oneway(),
    .ig_method_f_anova_special_main_interactions(),
    .ig_method_f_ancova_main_interactions(),
    .ig_method_f_hotelling_t2_one_group(),
    .ig_method_f_hotelling_t2_two_groups()
  )
}

# Shared numerics --------------------------------------------------------------------------------

.ig_g2_engine <- list(
  package = "stats",
  fun = paste("central F quantile from qbeta (no chi-square approximation for df2 > 4e5);",
              "noncentral F tails from the inferguidR Poisson mixture of pbeta terms (.ig_p_f)"),
  version = NA_character_
)

# Upper critical value of the central F distribution: P(F(df1, df2) > c) = alpha. With
# B = df1 F / (df1 F + df2) ~ Beta(df1 / 2, df2 / 2), c = (df2 / df1) b / (1 - b) for the upper
# alpha quantile b of B; when b > 1/2 the complement w = 1 - b comes from the lower alpha quantile
# of Beta(df2 / 2, df1 / 2), so neither form loses precision by cancellation. The result is checked
# against the central F tail and refined by root finding on the log scale if it is off by more
# than 1e-10 (relative) in alpha.
.ig_g2_crit_f <- function(alpha, df1, df2) {
  b <- suppressWarnings(stats::qbeta(alpha, df1 / 2, df2 / 2, lower.tail = FALSE))
  crit <- if (is.finite(b) && b <= 0.5) {
    df2 / df1 * b / (1 - b)
  } else {
    w <- suppressWarnings(stats::qbeta(alpha, df2 / 2, df1 / 2))
    df2 / df1 * (1 - w) / w
  }
  la <- log(alpha)
  ok <- is.finite(crit) && crit > 0 && {
    lp <- stats::pf(crit, df1, df2, lower.tail = FALSE, log.p = TRUE)
    is.finite(lp) && abs(lp - la) <= 1e-10 * max(1, abs(la))
  }
  if (!ok) {
    crit <- .ig_invert_cdf(function(x, lt) stats::pf(x, df1, df2, lower.tail = lt, log.p = TRUE),
                           alpha, lower_tail = FALSE, log_p = FALSE, support = c(0, Inf),
                           start = if (is.finite(crit) && crit > 0) crit else 1)
  }
  crit
}

# TRUE when beta = P(F' <= crit), F' ~ F(df1, df2, lambda), is provably below the double range.
# With F' = (X1 / df1) / (X2 / df2), X1 ~ chi-square(df1, lambda) and X2 ~ chi-square(df2):
#   P(F' <= c) <= P(X1 <= c df1 s / df2) + P(X2 > s)                       for any s > 0,
#   P(X1 <= x) <= P((Z + sqrt(lambda))^2 <= x) <= Phi(sqrt(x) - sqrt(lambda)) for x < lambda.
# s is the upper exp(-800) quantile of chi-square(df2); both terms below exp(-800) give beta = 0
# in double precision. Only central pchisq/qchisq and pnorm are used.
.ig_g2_beta_negligible <- function(crit, df1, df2, lambda) {
  if (!(lambda > 1000)) return(FALSE)
  s <- stats::qchisq(-800, df2, lower.tail = FALSE, log.p = TRUE)
  z <- sqrt(crit * df1 * s / df2) - sqrt(lambda)
  is.finite(z) && z < 0 && stats::pnorm(z, log.p = TRUE) < -800
}

# Power and beta of the upper-tailed F test with noncentrality lambda.
.ig_g2_f_test <- function(lambda, df1, df2, alpha) {
  crit <- .ig_g2_crit_f(alpha, df1, df2)
  if (.ig_g2_beta_negligible(crit, df1, df2, lambda)) {
    return(list(power = 1, beta = 0, critical = crit, bounded = TRUE))
  }
  power <- .ig_p_f(crit, df1, df2, lambda, lower_tail = FALSE)
  beta <- .ig_p_f(crit, df1, df2, lambda, lower_tail = TRUE)
  # The smaller tail keeps its directly summed value (relative accuracy); the larger one is its
  # complement. .ig_p_f forms a tail close to 1 from a beta argument that is itself within about
  # 1e-14 of 1 when the critical value is extreme, and that argument has lost its last digits: at
  # df (100, 1), lambda 1e6, alpha 1e-7 the directly summed beta was 0.99999159 where
  # 1 - power = 0.99998997 (power + beta - 1 = 1.6e-6; review 2026-10-06). The complement of the
  # small tail has no such loss.
  if (is.finite(power) && is.finite(beta)) {
    if (power <= beta) beta <- 1 - power else power <- 1 - beta
  }
  list(power = power, beta = beta, critical = crit, bounded = FALSE)
}

# Forward-model result shared by the five methods.
.ig_g2_forward <- function(p, lambda, df1, df2, n_total, group_sizes = NULL, assumptions = character()) {
  core <- .ig_g2_f_test(lambda, df1, df2, p$alpha)
  notes <- list()
  for (nm in c("power", "beta")) {
    if (!is.finite(core[[nm]])) {
      core[[nm]] <- NA_real_
      notes[[nm]] <- sprintf(paste("the noncentral F tail could not be summed (df1 = %s, df2 = %s,",
                                   "lambda = %s)"), format(df1), format(df2), format(signif(lambda, 8)))
    }
  }
  list(
    power = core$power,
    beta = core$beta,
    alpha_attained = p$alpha,
    critical = core$critical,
    df = c(df1, df2),
    ncp = lambda,
    n_total = n_total,
    group_sizes = group_sizes,
    statistic_distribution = list(
      h0 = list(family = "F", df1 = df1, df2 = df2, ncp = 0),
      h1 = list(family = "F", df1 = df1, df2 = df2, ncp = lambda),
      rejection_side = "upper"
    ),
    outputs = list(numerator_df = df1, denominator_df = df2),
    notes = notes,
    assumptions = c(assumptions,
                    "upper-tailed F test: reject H0 when F > F_crit with P(F(df1, df2) > F_crit) = alpha",
                    if (isTRUE(core$bounded)) paste(
                      "beta is below the double-precision range (rigorous bound), so power = 1 and",
                      "beta = 0 are exact to double precision"))
  )
}

# The analysis mode, read from the inputs a specification carries (constraints receive inputs
# only): compromise has error_ratio, sensitivity lacks the effect, criterion lacks alpha, post hoc
# lacks power, a priori has effect, alpha and power.
.ig_g2_mode <- function(p, effect) {
  if (!is.null(p$error_ratio)) return("compromise")
  if (is.null(p[[effect]])) return("sensitivity")
  if (is.null(p$alpha)) return("criterion")
  if (is.null(p$power)) return("post_hoc")
  "a_priori"
}

# Constraints shared by the five methods; `design_check(p)` adds the design-specific ones.
.ig_g2_constraints <- function(effect, effect_label, design_check) {
  force(effect)
  force(effect_label)
  force(design_check)
  function(p) {
    mode <- .ig_g2_mode(p, effect)
    out <- character()
    if (mode == "a_priori" && isTRUE(p[[effect]] <= 0)) {
      out <- c(out, sprintf(paste(
        "in a priori mode the effect size %s must be positive: with %s = 0 the power equals alpha for",
        "every sample size, so no N reaches a target above alpha (G*Power's admissible range is",
        "1e-06 to 50)"), effect_label, effect))
    }
    if (mode %in% c("a_priori", "sensitivity") && !is.null(p$alpha) && !is.null(p$power) &&
        p$power < p$alpha) {
      out <- c(out, sprintf(paste(
        "the target power %s is below alpha %s; the power of this F test is never below alpha",
        "(G*Power: 'Power must be greater than or equal to alpha')"),
        .ig_format_double(p$power), .ig_format_double(p$alpha)))
    }
    c(out, design_check(p, mode))
  }
}

# A priori search over single-N designs N = n_min + step (m - 1), m = 1, 2, ... (step = k for the
# one-way ANOVA, whose designs are multiples of k; 1 otherwise), with the generic integer search
# (R/num-root.R) applied to the index m. `forward(p)` is the method's forward model.
.ig_g2_search_single <- function(forward, n_min_fun, step_fun, grid_text) {
  force(forward)
  force(n_min_fun)
  force(step_fun)
  force(grid_text)
  function(p, n_lo, n_max) {
    n_min <- as.double(n_min_fun(p))
    step <- as.double(step_fun(p))
    while (n_min < n_lo) n_min <- n_min + step
    if (n_min > n_max) {
      return(list(status = "infeasible", reason = sprintf(
        "the smallest admissible design (N = %s) exceeds the largest searched sample size %s",
        format(n_min, scientific = FALSE), format(n_max, scientific = FALSE)), diagnostics = list()))
    }
    m_max <- floor((n_max - n_min) / step) + 1
    n_of <- function(m) n_min + step * (m - 1)
    power_at <- function(m) {
      q <- p
      q$n_total <- n_of(m)
      forward(q)$power
    }
    s <- .ig_integer_search(power_at = power_at, admissible = function(m) TRUE, n_lo = 1, n_max = m_max,
                            target = p$power, discrete = FALSE)
    d <- s$diagnostics
    for (nm in c("n_lower", "n_upper_limit", "predecessor_n_total")) {
      if (is.numeric(d[[nm]])) d[[nm]] <- n_of(d[[nm]])
    }
    for (nm in c("bracket", "enumerated_range")) {
      if (is.numeric(d[[nm]])) d[[nm]] <- n_of(d[[nm]])
    }
    if (is.list(d$power_sequence) && is.numeric(d$power_sequence$n_total)) {
      d$power_sequence$n_total <- n_of(d$power_sequence$n_total)
    }
    d$admissible_designs <- grid_text
    d$search_index <- "the integer search runs over m with N = n_lower + step * (m - 1)"
    d$step <- step
    if (!identical(s$status, "ok")) {
      reason <- if (identical(s$status, "infeasible")) {
        sprintf("the target power %s is not reached by any admissible design with N <= %s%s",
                format(p$power), format(n_of(m_max), scientific = FALSE),
                if (is.numeric(d$power_at_limit)) sprintf(" (power there: %s)", format(signif(d$power_at_limit, 7))) else "")
      } else s$reason
      return(list(status = s$status, reason = reason, diagnostics = d))
    }
    list(status = "ok", n = n_of(s$n), power = s$power, diagnostics = d)
  }
}

.ig_g2_effect_f <- function() {
  ig_param("f", role = "effect", domain = ig_domain_real(0, Inf),
           label = list(en = "Effect size f", de = "Effektst\u00e4rke f"),
           description = paste("Cohen's f = sigma_m / sigma of the tested effect (f^2 = eta^2 / (1 - eta^2)",
                               "with the partial eta^2 of that effect). Must be positive in a priori mode."),
           aliases = "effect_size_f")
}

.ig_g2_n_total <- function(modes = NULL) {
  ig_param("n_total", role = "design", domain = ig_domain_integer(lower = 2), unit = "subjects",
           label = list(en = "Total sample size", de = "Gesamtstichprobenumfang"),
           description = "Total number of subjects N (solved in a priori mode).",
           aliases = "total_sample_size", modes = modes)
}

.ig_g2_reference <- function(path, label, manual, rule, id) {
  list(gpower_path = path, gpower_label = label, manual_pages = manual, allocation_rule = rule,
       method_notes = sprintf("specs/method-notes/%s.md", id))
}

.ig_g2_count_text <- function(x) format(x, scientific = FALSE, trim = TRUE)

# ANOVA: Fixed effects, omnibus, one-way -----------------------------------------------------------
#
# Manual ch. 10: H0 mu_1 = ... = mu_k; lambda = f^2 N, df1 = k - 1, df2 = N - k (10.5). For N not
# a multiple of k G*Power warns that "the analysis was based on the average group size" and keeps
# lambda = f^2 N and df2 = N - k (captures GPW-F-ANOVA1-POS-B001-O01 and -E01).

.ig_method_f_anova_oneway <- function() {
  id <- "f_anova_oneway"
  ig_method(
    id = id,
    family = "F",
    label = list(en = "ANOVA: Fixed effects, omnibus, one-way",
                 de = "ANOVA: Feste Effekte, global, einfaktoriell"),
    params = ig_params(
      .ig_g2_effect_f(),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      .ig_g2_n_total(),
      ig_param("n_groups", role = "design", domain = ig_domain_integer(lower = 2), unit = "groups",
               label = list(en = "Number of groups", de = "Anzahl der Gruppen"),
               description = "Number of groups k (numerator df k - 1).",
               aliases = c("number_of_groups", "k"))
    ),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = .ig_power_f_anova_oneway,
    design = NULL,
    discrete = FALSE,
    effect_direction = "increasing",
    engine = .ig_g2_engine,
    checker = "validation/independent/g2/check_g2.R",
    reference = .ig_g2_reference(
      "Tests > Means > Many groups: ANOVA: One-way (one independent variable)",
      "F tests - ANOVA: Fixed effects, omnibus, one-way",
      "ch. 10 (PDF pp. 24-25); lambda = f^2 N, df (k - 1, N - k) in 10.5",
      paste("equal group sizes: a priori N is the smallest multiple of k with at least 2 subjects per",
            "group whose power reaches the target (DOC-055, captures GPW-F-ANOVA1-APR-*); other modes",
            "accept any N > k (lambda = f^2 N, as G*Power's 'average group size')"),
      id),
    status = "independently_checked",
    constraints = .ig_g2_constraints("f", "f", function(p, mode) {
      if (!is.null(p$n_total) && !is.null(p$n_groups) && p$n_total <= p$n_groups) {
        sprintf("the total sample size N = %s must exceed the number of groups k = %s (df2 = N - k >= 1)",
                .ig_g2_count_text(p$n_total), .ig_g2_count_text(p$n_groups))
      } else character()
    }),
    a_priori_search = .ig_g2_search_single(
      function(p) .ig_power_f_anova_oneway(p),
      n_min_fun = function(p) 2 * as.double(p$n_groups),
      step_fun = function(p) as.double(p$n_groups),
      grid_text = "N = k n with n >= 2 subjects in each of the k groups (equal group sizes)")
  )
}

.ig_power_f_anova_oneway <- function(p) {
  n <- as.double(p$n_total)
  k <- as.double(p$n_groups)
  equal <- n %% k == 0
  gs <- if (equal && k <= 1e5) {
    stats::setNames(rep(n / k, k), sprintf("n_group_%0*d", nchar(format(k, scientific = FALSE)), seq_len(k)))
  } else NULL
  .ig_g2_forward(
    p, lambda = p$f^2 * n, df1 = k - 1, df2 = n - k, n_total = n, group_sizes = gs,
    assumptions = c(
      "k independent groups from normal populations with a common variance; fixed effects",
      "f = sigma_m / sigma with sigma_m the (size-weighted) SD of the group means; lambda = f^2 N, df1 = k - 1, df2 = N - k (manual 10.5)",
      if (!equal) sprintf(paste(
        "N = %s is not a multiple of k = %s: as in G*Power, the analysis is based on the average group",
        "size N / k (lambda = f^2 N, df2 = N - k); no integer group sizes are reported"),
        .ig_g2_count_text(n), .ig_g2_count_text(k))
    ))
}

# ANOVA: Fixed effects, special, main effects and interactions -------------------------------------
#
# Manual ch. 11: lambda = f^2 N, df1 = numerator df of the tested effect, df2 = N - k with k the
# total number of groups (cells) of the design (11.5). No restriction of N to multiples of k:
# DOC-066 prints N = 2283 for k = 30. k >= 2: G*Power refuses a single group ("There must be at
# least 2 groups.", REV-G2-F-ANOVASP-POS-S02; the ANCOVA procedure likewise, REV-G2-F-ANCOVA-POS-S02).

.ig_method_f_anova_special_main_interactions <- function() {
  id <- "f_anova_special_main_interactions"
  ig_method(
    id = id,
    family = "F",
    label = list(en = "ANOVA: Fixed effects, special, main effects and interactions",
                 de = "ANOVA: Feste Effekte, speziell, Haupteffekte und Interaktionen"),
    params = ig_params(
      .ig_g2_effect_f(),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      .ig_g2_n_total(),
      ig_param("numerator_df", role = "design", domain = ig_domain_integer(lower = 1), unit = "df",
               label = list(en = "Numerator df", de = "Z\u00e4hler-df"),
               description = "Numerator degrees of freedom u of the tested main effect, interaction or contrast.",
               aliases = "df1"),
      ig_param("n_groups", role = "design", domain = ig_domain_integer(lower = 2), unit = "groups",
               label = list(en = "Number of groups", de = "Anzahl der Gruppen"),
               description = "Total number of groups (cells) k of the design, at least 2; df2 = N - k.",
               aliases = c("number_of_groups", "k"))
    ),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = .ig_power_f_anova_special,
    design = NULL,
    discrete = FALSE,
    effect_direction = "increasing",
    engine = .ig_g2_engine,
    checker = "validation/independent/g2/check_g2.R",
    reference = .ig_g2_reference(
      "Tests > Means > Many groups: ANOVA: Main effects and interactions (two or more independent variables)",
      "F tests - ANOVA: Fixed effects, special, main effects and interactions",
      "ch. 11 (PDF pp. 26-30); lambda = f^2 N, df (u, N - k) in 11.5",
      "a priori N is the smallest integer N > k whose power reaches the target (DOC-066: N = 2283 for k = 30)",
      id),
    status = "independently_checked",
    constraints = .ig_g2_constraints("f", "f", function(p, mode) {
      if (!is.null(p$n_total) && !is.null(p$n_groups) && p$n_total <= p$n_groups) {
        sprintf("the total sample size N = %s must exceed the number of groups k = %s (df2 = N - k >= 1)",
                .ig_g2_count_text(p$n_total), .ig_g2_count_text(p$n_groups))
      } else character()
    }),
    a_priori_search = .ig_g2_search_single(
      function(p) .ig_power_f_anova_special(p),
      n_min_fun = function(p) as.double(p$n_groups) + 1,
      step_fun = function(p) 1,
      grid_text = "every integer N > k (df2 = N - k >= 1)")
  )
}

.ig_power_f_anova_special <- function(p) {
  n <- as.double(p$n_total)
  k <- as.double(p$n_groups)
  u <- as.double(p$numerator_df)
  .ig_g2_forward(
    p, lambda = p$f^2 * n, df1 = u, df2 = n - k, n_total = n,
    assumptions = c(
      "fixed-effects factorial design with k cells, normal errors with a common variance",
      "f = sqrt(partial eta^2 / (1 - partial eta^2)) of the tested effect; lambda = f^2 N, df1 = u, df2 = N - k (manual 11.5)"
    ))
}

# ANCOVA: Fixed effects, main effects and interactions ---------------------------------------------
#
# No manual chapter. The captures show lambda = f^2 N, df1 = u and df2 = N - k - c with c the number
# of covariates (GPW-F-ANCOVA-POS-B001-O01: N 120, k 6, c 1 gives df2 113 and lambda 7.5); the
# critical F and power follow from F(u, N - k - c, f^2 N).

.ig_method_f_ancova_main_interactions <- function() {
  id <- "f_ancova_main_interactions"
  ig_method(
    id = id,
    family = "F",
    label = list(en = "ANCOVA: Fixed effects, main effects and interactions",
                 de = "ANCOVA: Feste Effekte, Haupteffekte und Interaktionen"),
    params = ig_params(
      .ig_g2_effect_f(),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      .ig_g2_n_total(),
      ig_param("numerator_df", role = "design", domain = ig_domain_integer(lower = 1), unit = "df",
               label = list(en = "Numerator df", de = "Z\u00e4hler-df"),
               description = "Numerator degrees of freedom u of the tested (covariate-adjusted) effect.",
               aliases = "df1"),
      ig_param("n_groups", role = "design", domain = ig_domain_integer(lower = 2), unit = "groups",
               label = list(en = "Number of groups", de = "Anzahl der Gruppen"),
               description = "Total number of groups (cells) k of the design, at least 2.",
               aliases = c("number_of_groups", "k")),
      ig_param("n_covariates", role = "design", domain = ig_domain_integer(lower = 0), unit = "covariates",
               label = list(en = "Number of covariates", de = "Anzahl der Kovariaten"),
               description = "Number of covariates c; each costs one error degree of freedom (df2 = N - k - c).",
               aliases = "number_of_covariates")
    ),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = .ig_power_f_ancova,
    design = NULL,
    discrete = FALSE,
    effect_direction = "increasing",
    engine = .ig_g2_engine,
    checker = "validation/independent/g2/check_g2.R",
    reference = .ig_g2_reference(
      "Tests > Means > Many groups: ANCOVA: Fixed effects, main effects and interactions",
      "F tests - ANCOVA: Fixed effects, main effects and interactions",
      "no manual chapter; model established from the emulator captures (method notes)",
      "a priori N is the smallest integer N > k + c whose power reaches the target",
      id),
    status = "independently_checked",
    constraints = .ig_g2_constraints("f", "f", function(p, mode) {
      if (!is.null(p$n_total) && !is.null(p$n_groups) && !is.null(p$n_covariates) &&
          p$n_total <= as.double(p$n_groups) + as.double(p$n_covariates)) {
        sprintf(paste("the total sample size N = %s must exceed the number of groups plus covariates",
                      "k + c = %s (df2 = N - k - c >= 1)"),
                .ig_g2_count_text(p$n_total), .ig_g2_count_text(as.double(p$n_groups) + as.double(p$n_covariates)))
      } else character()
    }),
    a_priori_search = .ig_g2_search_single(
      function(p) .ig_power_f_ancova(p),
      n_min_fun = function(p) as.double(p$n_groups) + as.double(p$n_covariates) + 1,
      step_fun = function(p) 1,
      grid_text = "every integer N > k + c (df2 = N - k - c >= 1)")
  )
}

.ig_power_f_ancova <- function(p) {
  n <- as.double(p$n_total)
  k <- as.double(p$n_groups)
  cv <- as.double(p$n_covariates)
  u <- as.double(p$numerator_df)
  .ig_g2_forward(
    p, lambda = p$f^2 * n, df1 = u, df2 = n - k - cv, n_total = n,
    assumptions = c(
      "fixed-effects ANCOVA with k cells and c covariates (fixed covariate values, common within-cell regression slopes), normal errors with a common variance",
      "f = sqrt(partial eta^2 / (1 - partial eta^2)) of the covariate-adjusted effect; lambda = f^2 N, df1 = u, df2 = N - k - c",
      "no G*Power manual chapter: the model was established from the emulator captures (method notes)"
    ))
}

# Hotelling's T^2: one group mean vector ------------------------------------------------------------
#
# No manual chapter. T^2 = N (xbar - mu0)' S^-1 (xbar - mu0); (N - p) / (p (N - 1)) T^2 ~
# F(p, N - p, lambda = N Delta^2) with Delta^2 = (mu - mu0)' Sigma^-1 (mu - mu0) (Anderson 2003,
# 5.2). The captures confirm df (p, N - p) and lambda = N Delta^2 (GPW-F-HOT1-POS-B001-O01).

.ig_method_f_hotelling_t2_one_group <- function() {
  id <- "f_hotelling_t2_one_group"
  ig_method(
    id = id,
    family = "F",
    label = list(en = "Hotellings T\u00b2: One group mean vector",
                 de = "Hotellings T\u00b2: Mittelwertvektor einer Gruppe"),
    params = ig_params(
      .ig_g2_effect_delta(),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      .ig_g2_n_total(),
      .ig_g2_n_response()
    ),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = .ig_power_f_hotelling_one,
    design = NULL,
    discrete = FALSE,
    effect_direction = "increasing",
    engine = .ig_g2_engine,
    checker = "validation/independent/g2/check_g2.R",
    reference = .ig_g2_reference(
      "Tests > Means > Multivariate: Hotelling's T2: One group mean vector",
      "F tests - Hotellings T\u00b2: One group mean vector",
      "no manual chapter; model established from the emulator captures (method notes)",
      "a priori N is the smallest integer N > p whose power reaches the target",
      id),
    status = "independently_checked",
    constraints = .ig_g2_constraints("delta", "Delta", function(p, mode) {
      if (!is.null(p$n_total) && !is.null(p$n_response) && p$n_total <= p$n_response) {
        sprintf(paste("the sample size N = %s must exceed the number of response variables p = %s",
                      "(df2 = N - p >= 1)"), .ig_g2_count_text(p$n_total), .ig_g2_count_text(p$n_response))
      } else character()
    }),
    a_priori_search = .ig_g2_search_single(
      function(p) .ig_power_f_hotelling_one(p),
      n_min_fun = function(p) as.double(p$n_response) + 1,
      step_fun = function(p) 1,
      grid_text = "every integer N > p (df2 = N - p >= 1)")
  )
}

.ig_g2_effect_delta <- function() {
  ig_param("delta", role = "effect", domain = ig_domain_real(0, Inf),
           label = list(en = "Effect size \u0394", de = "Effektst\u00e4rke \u0394"),
           description = paste("Mahalanobis distance Delta = sqrt(d' Sigma^-1 d) of the mean difference d",
                               "(mu - mu0 for one group, mu1 - mu2 for two groups) with the common covariance",
                               "matrix Sigma. Must be positive in a priori mode."),
           aliases = "effect_size_delta")
}

.ig_g2_n_response <- function() {
  ig_param("n_response", role = "design", domain = ig_domain_integer(lower = 1), unit = "variables",
           label = list(en = "Response variables", de = "Zielvariablen"),
           description = "Number of response variables p (numerator df).",
           aliases = c("response_variables", "number_of_response_variables"))
}

.ig_power_f_hotelling_one <- function(p) {
  n <- as.double(p$n_total)
  pv <- as.double(p$n_response)
  .ig_g2_forward(
    p, lambda = n * p$delta^2, df1 = pv, df2 = n - pv, n_total = n,
    assumptions = c(
      "one sample from a p-variate normal distribution with unknown covariance matrix",
      "F = (N - p) / (p (N - 1)) T^2 ~ F(p, N - p, lambda = N Delta^2) (Anderson 2003, 5.2)",
      "no G*Power manual chapter: the model was confirmed with the emulator captures (method notes)"
    ))
}

# Hotelling's T^2: two group mean vectors -----------------------------------------------------------
#
# No manual chapter. T^2 = n1 n2 / N (xbar1 - xbar2)' S^-1 (xbar1 - xbar2) with the pooled S;
# (N - p - 1) / (p (N - 2)) T^2 ~ F(p, N - p - 1, lambda = n1 n2 / N Delta^2) (Anderson 2003, 5.3).
# The captures confirm df2 = N - p - 1 and lambda = n1 n2 / N Delta^2 (GPW-F-HOT2-POS-B001-O01).

.ig_method_f_hotelling_t2_two_groups <- function() {
  id <- "f_hotelling_t2_two_groups"
  ig_method(
    id = id,
    family = "F",
    label = list(en = "Hotellings T\u00b2: Two group mean vectors",
                 de = "Hotellings T\u00b2: Mittelwertvektoren zweier Gruppen"),
    params = ig_params(
      .ig_g2_effect_delta(),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      ig_param("n_total", role = "design", domain = ig_domain_integer(lower = 3), unit = "subjects",
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
               description = paste("Group-size ratio n2/n1 of the a priori design (G*Power's initial value",
                                   "is 1, the default here)."),
               aliases = "allocation_ratio_n2_n1",
               modes = "a_priori"),
      .ig_g2_n_response()
    ),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = .ig_power_f_hotelling_two,
    design = .ig_design_f_hotelling_two,
    discrete = FALSE,
    effect_direction = "increasing",
    engine = .ig_g2_engine,
    checker = "validation/independent/g2/check_g2.R",
    reference = .ig_g2_reference(
      "Tests > Means > Multivariate: Hotelling's T2: Two group mean vectors",
      "F tests - Hotellings T\u00b2: Two group mean vectors",
      "no manual chapter; model and allocation rule established from the emulator captures (method notes)",
      paste("G*Power's rule (captures): N* = smallest N (even N when r = 1) whose power with fractional",
            "group sizes N/(1+r), rN/(1+r) reaches the target; n1 = round(N*/(1+r)) half up, n2 = N* - n1.",
            "If that design has an empty group or misses the target with its integer group sizes, the next",
            "N under the same split whose design reaches the target is reported (G*Power's design is given",
            "in diagnostics$gpower_rule)"),
      id),
    status = "independently_checked",
    constraints = .ig_g2_constraints("delta", "Delta", function(p, mode) {
      if (!is.null(p$n1) && !is.null(p$n2) && !is.null(p$n_response) &&
          as.double(p$n1) + as.double(p$n2) < as.double(p$n_response) + 2) {
        sprintf(paste("n1 + n2 = %s must be at least the number of response variables plus 2 = %s",
                      "(df2 = n1 + n2 - p - 1 >= 1)"),
                .ig_g2_count_text(as.double(p$n1) + as.double(p$n2)),
                .ig_g2_count_text(as.double(p$n_response) + 2))
      } else character()
    }),
    a_priori_search = .ig_g2_search_hotelling_two
  )
}

# Integer split of a total N (G*Power's rule, captures GPW-F-HOT2-APR-* and SUP-G2-F-HOT2-APR-*):
# n1 = round(N / (1 + r)) half up, n2 = N - n1; with r = 1 only even totals are designs
# (n1 = n2 = N / 2). Used by the a priori search and by the solver's admissibility check of the
# design the search reports. `continuous = TRUE` returns the fractional split.
.ig_design_f_hotelling_two <- function(p, n_total, continuous = FALSE) {
  r <- as.double(p$allocation_ratio)
  if (continuous) return(c(n1 = n_total / (1 + r), n2 = r * n_total / (1 + r)))
  if (r == 1 && n_total %% 2 != 0) return(NULL)
  n1 <- floor(n_total / (1 + r) + 0.5)
  n2 <- n_total - n1
  if (n1 < 1 || n2 < 1) return(NULL)
  c(n1 = n1, n2 = n2)
}

.ig_g2_split_gpower <- function(n, r) {
  n1 <- floor(n / (1 + r) + 0.5)
  c(n1 = n1, n2 = n - n1)
}

# A priori search of the Hotelling two-group test: N* by the power of the fractional design (for
# r = 1 over even N, where the fractional and the integer design coincide), G*Power's integer
# split of N*, and if that design has an empty group or misses the target, the next N under the
# same split whose design has two groups and reaches the target.
.ig_g2_search_hotelling_two <- function(p, n_lo, n_max) {
  r <- as.double(p$allocation_ratio)
  pv <- as.double(p$n_response)
  target <- p$power
  equal <- r == 1
  at <- function(n1, n2) {
    q <- p
    q$n1 <- n1
    q$n2 <- n2
    .ig_power_f_hotelling_two(q)$power
  }
  frac_power <- function(n) at(n / (1 + r), r * n / (1 + r))
  # Smallest total with df2 = N - p - 1 >= 1.
  n_min <- max(n_lo, pv + 2)
  s <- .ig_integer_search(power_at = frac_power,
                          admissible = if (equal) function(n) n %% 2 == 0 else function(n) TRUE,
                          n_lo = n_min, n_max = n_max, target = target, discrete = FALSE)
  if (!identical(s$status, "ok")) {
    return(list(status = s$status, reason = paste(s$reason, "(power with fractional group sizes)"),
                diagnostics = s$diagnostics))
  }
  n_star <- s$n
  gp <- .ig_g2_split_gpower(n_star, r)
  gp_ok <- all(gp >= 1)
  gp_power <- if (gp_ok) at(gp[["n1"]], gp[["n2"]]) else NA_real_
  rule <- list(
    n_total = n_star,
    power_fractional = s$power,
    n1 = gp[["n1"]], n2 = gp[["n2"]],
    power_integer_design = gp_power,
    note = paste("G*Power's a priori design: N* = smallest N (even N when r = 1) whose power with",
                 "fractional group sizes N/(1+r), rN/(1+r) reaches the target; n1 = round(N*/(1+r))",
                 "half up, n2 = N* - n1")
  )
  diag <- c(s$diagnostics, list(criterion = "power with fractional group sizes (G*Power rule)",
                                gpower_rule = rule))
  # The integer search above ran on the power with fractional group sizes. For r = 1 that is the
  # power of the design itself (n1 = n2 = N / 2), so its minimality statement stands. For r != 1 it
  # is G*Power's criterion: a smaller total can reach the target with its integer split (capture
  # REV-G2-F-HOT2-APR-R04: 62 / 93 is reported, 62 / 92 also reaches 0.80), so the statement of the
  # search is replaced by what was shown.
  if (!equal) {
    diag$search_scale <- paste("bracket, predecessor and enumerated range refer to the power with",
                               "fractional group sizes, not to integer designs")
  }
  chosen <- NULL
  warn <- character()
  if (gp_ok && gp_power >= target) {
    chosen <- gp
    diag$design_choice <- "G*Power's design reaches the target with its integer group sizes and is reported"
    if (!equal) {
      diag$minimality <- sprintf(paste(
        "shown for the fractional criterion only: N = %s is the smallest total whose power with fractional",
        "group sizes reaches the target (G*Power's rule). A smaller total can reach the target with its",
        "integer split; such designs are not searched, so N is not claimed to be the smallest design"),
        .ig_g2_count_text(n_star))
    }
  } else {
    # Next total under the same split whose design has two non-empty groups and reaches the target.
    # seq_len() keeps the scan inside the search limit when N* is the limit itself (seq() would
    # count down and leave it).
    n_last <- min(n_max, n_star + 1e5)
    for (n in n_star + seq_len(max(0, n_last - n_star))) {
      if (equal && n %% 2 != 0) next
      g <- .ig_g2_split_gpower(n, r)
      if (any(g < 1)) next
      if (at(g[["n1"]], g[["n2"]]) >= target) {
        chosen <- g
        break
      }
    }
    if (is.null(chosen) && n_last >= n_max) {
      return(list(status = "infeasible", reason = sprintf(paste(
        "G*Power's design for this problem (n1 = %s, n2 = %s) %s, and no later design under the same",
        "allocation rule reaches the target up to the largest searched sample size %s"),
        .ig_g2_count_text(gp[["n1"]]), .ig_g2_count_text(gp[["n2"]]),
        if (!gp_ok) "has an empty group" else "misses the target power",
        .ig_g2_count_text(n_max)), diagnostics = diag))
    }
    if (is.null(chosen)) {
      return(list(status = "numerical_failure", reason = paste(
        "no design under G*Power's allocation rule with two non-empty groups reaches the target within",
        "1e5 totals above the fractional solution"), diagnostics = diag))
    }
    diag$design_choice <- sprintf(paste(
      "G*Power's design (n1 = %s, n2 = %s) %s; inferguidR reports the next design under the same",
      "allocation rule with two non-empty groups whose integer group sizes reach the target"),
      .ig_g2_count_text(gp[["n1"]]), .ig_g2_count_text(gp[["n2"]]),
      if (!gp_ok) "has an empty group" else sprintf("misses the target (power %s < %s)",
                                                   format(gp_power, digits = 10), format(target)))
    diag$minimality <- sprintf(paste(
      "N = %s is the first total above G*Power's N = %s whose integer split has two non-empty groups and",
      "reaches the target; totals below G*Power's N were not searched with the integer split, so N is not",
      "claimed to be the smallest design"),
      .ig_g2_count_text(sum(chosen)), .ig_g2_count_text(n_star))
    warn <- paste0(diag$design_choice, " (diagnostics$gpower_rule)")
  }
  total <- sum(chosen)
  list(status = "ok", n = total, power = at(chosen[["n1"]], chosen[["n2"]]), diagnostics = diag,
       warnings = warn)
}

.ig_power_f_hotelling_two <- function(p) {
  n1 <- as.double(p$n1)
  n2 <- as.double(p$n2)
  n <- n1 + n2
  pv <- as.double(p$n_response)
  whole <- n1 == round(n1) && n2 == round(n2)
  .ig_g2_forward(
    p, lambda = n1 * n2 / n * p$delta^2, df1 = pv, df2 = n - pv - 1, n_total = n,
    group_sizes = if (whole) c(n1 = n1, n2 = n2) else NULL,
    assumptions = c(
      "two independent samples from p-variate normal distributions with a common covariance matrix",
      "F = (N - p - 1) / (p (N - 2)) T^2 ~ F(p, N - p - 1, lambda = n1 n2 / N Delta^2) (Anderson 2003, 5.3)",
      "no G*Power manual chapter: the model was confirmed with the emulator captures (method notes)"
    ))
}
