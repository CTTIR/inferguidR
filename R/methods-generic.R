# Method definitions: generic tests (group G1) -------------------------------------------
#
# Ownership: R/methods-generic.R belongs to group G1 (phase 06, decision D09).
# G*Power's four generic tests take the noncentrality parameter and the
# degrees of freedom directly; they offer no a priori analysis (manual ch. 25
# for the generic t test; the z, F and chi-square generic tests have no
# manual chapter and are specified from the certified live inventory and the
# emulator captures):
#
#   t_generic     H0 t(df),            H1 t(df, delta)
#   z_generic     H0 N(0, 1),          H1 N(mu, sigma)
#   f_generic     H0 F(df1, df2),      H1 F(df1, df2, lambda)
#   chi2_generic  H0 chi2(df),         H1 chi2(df, lambda)
#
# None of them has an effect-size drawer (registry: drawer null). Method
# notes: specs/method-notes/<method_id>.md.

.ig_define_methods_generic_g1 <- function() {
  list(.ig_method_t_generic(), .ig_method_z_generic(), .ig_method_f_generic(),
       .ig_method_chi2_generic())
}

.ig_generic_modes_g1 <- c("post_hoc", "sensitivity", "criterion", "compromise")

# Generic t test ------------------------------------------------------------------------
#
# G*Power: Tests > Generic > Generic t; label "t tests - Generic t test"
# (manual ch. 25). Inputs: Tail(s), noncentrality parameter delta, df (real,
# > 0; G*Power refuses df <= 0.01). Tail convention 25.5 as for every t test
# (.ig_t_power_g1() in R/methods-t.R).

.ig_method_t_generic <- function() {
  ig_method(
    id = "t_generic",
    family = "t",
    label = list(en = "Generic t test", de = "Generischer t-Test"),
    params = ig_params(
      ig_param("delta", role = "effect", domain = ig_domain_real(),
               label = list(en = "Noncentrality parameter \u03b4", de = "Nichtzentralit\u00e4tsparameter \u03b4"),
               description = "Noncentrality parameter of the H1 distribution t(df, delta) (manual 25.1).",
               aliases = c("noncentrality_parameter", "ncp")),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      ig_param("df", role = "design",
               domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE),
               unit = "degrees of freedom",
               label = list(en = "Df", de = "Freiheitsgrade"),
               description = paste("Degrees of freedom of the t distribution (any positive real number;",
                                   "G*Power refuses df <= 0.01)."),
               aliases = "degrees_of_freedom")
    ),
    options = ig_options(.ig_t_tails_option_g1()),
    modes = .ig_generic_modes_g1,
    power = .ig_power_t_generic,
    design = NULL,
    discrete = FALSE,
    effect_direction = "increasing",
    engine = .ig_t_engine_g1(),
    checker = "validation/independent/g1/check_g1.R",
    reference = list(
      gpower_path = "Tests > Generic > Generic t",
      gpower_label = "t tests - Generic t test",
      manual_pages = "ch. 25 (PDF p. 63), tail convention 25.5",
      method_notes = "specs/method-notes/t_generic.md"
    ),
    status = "independently_checked"
  )
}

.ig_power_t_generic <- function(p) {
  df <- as.double(p$df)
  ncp <- as.double(p$delta)
  core <- .ig_t_power_g1(df, ncp, p$alpha, p$tails)
  list(
    power = core$power,
    beta = core$beta,
    alpha_attained = p$alpha,
    critical = core$critical,
    df = df,
    ncp = ncp,
    statistic_distribution = list(
      h0 = list(family = "t", df = df, ncp = 0),
      h1 = list(family = "t", df = df, ncp = ncp),
      rejection_side = core$side
    ),
    assumptions = c(
      "H0: central t(df); H1: noncentral t(df, delta), delta and df given directly (manual 25.1)",
      .ig_t_tail_note_g1(p$tails, "delta")
    )
  )
}

# Generic z test ------------------------------------------------------------------------
#
# G*Power: Tests > Generic > Generic z; label "z tests - Generic z test" (no
# manual chapter). H0: Z ~ N(0, 1); H1: Z ~ N(mu, sigma), mu the
# "Noncentrality parameter mu", sigma the "Noncentral dist. SD sigma"
# (G*Power refuses sigma < 1e-8). Critical values are standard normal
# quantiles: one tail, z_c = sign(mu) * z(1 - alpha) (upper tail for
# mu = 0); two tails, -z(1 - alpha/2) and z(1 - alpha/2). The power is the
# H1 probability of the rejection region, so with sigma != 1 the power at
# mu = 0 differs from alpha (Phi(-z(1 - alpha)/sigma) for one tail).

.ig_method_z_generic <- function() {
  ig_method(
    id = "z_generic",
    family = "z",
    label = list(en = "Generic z test", de = "Generischer z-Test"),
    params = ig_params(
      ig_param("mu", role = "effect", domain = ig_domain_real(),
               label = list(en = "Noncentrality parameter \u03bc", de = "Nichtzentralit\u00e4tsparameter \u03bc"),
               description = "Mean of the H1 distribution N(mu, sigma) of the test statistic.",
               aliases = c("noncentrality_parameter", "ncp")),
      ig_param("sigma", role = "nuisance",
               domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE),
               default = 1,
               label = list(en = "Noncentral dist. SD \u03c3", de = "SD der nichtzentralen Verteilung \u03c3"),
               description = paste("Standard deviation of the H1 distribution N(mu, sigma) of the test statistic",
                                   "(G*Power's initial field value 1, live capture; G*Power refuses sigma < 1e-8)."),
               aliases = c("noncentral_sd", "h1_sd")),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio()
    ),
    options = ig_options(.ig_t_tails_option_g1()),
    modes = .ig_generic_modes_g1,
    power = .ig_power_z_generic,
    design = NULL,
    discrete = FALSE,
    effect_direction = "increasing",
    engine = list(package = "stats", fun = "pnorm/qnorm", version = NA_character_),
    checker = "validation/independent/g1/check_g1.R",
    reference = list(
      gpower_path = "Tests > Generic > Generic z",
      gpower_label = "z tests - Generic z test",
      manual_pages = "no manual chapter (live inventory and emulator captures)",
      method_notes = "specs/method-notes/z_generic.md"
    ),
    status = "independently_checked"
  )
}

.ig_power_z_generic <- function(p) {
  mu <- as.double(p$mu)
  sigma <- as.double(p$sigma)
  alpha <- p$alpha
  # Tail of X ~ N(mu, sigma): P(X <= q) or P(X > q).
  ptail <- function(q, lower_tail) .ig_p_norm((q - mu) / sigma, lower_tail = lower_tail)
  if (identical(p$tails, "two")) {
    c2 <- .ig_q_norm(alpha / 2, lower_tail = FALSE)
    crit <- c(-c2, c2)
    region <- .ig_two_sided_prob(ptail, crit[1], crit[2])
    power <- region$reject
    beta <- region$accept
    side <- "both"
  } else {
    c1 <- .ig_q_norm(alpha, lower_tail = FALSE)
    # Evaluated at |mu| (mirror image of the test at -mu), as for the t tests.
    m <- abs(mu)
    power <- .ig_p_norm((c1 - m) / sigma, lower_tail = FALSE)
    beta <- .ig_p_norm((c1 - m) / sigma, lower_tail = TRUE)
    side <- if (mu < 0) "lower" else "upper"
    crit <- if (mu < 0) -c1 else c1
  }
  list(
    power = power,
    beta = beta,
    alpha_attained = alpha,
    critical = crit,
    ncp = mu,
    statistic_distribution = list(
      h0 = list(family = "normal", mean = 0, sd = 1),
      h1 = list(family = "normal", mean = mu, sd = sigma),
      rejection_side = side
    ),
    assumptions = c(
      "H0: Z ~ N(0, 1); H1: Z ~ N(mu, sigma), mu and sigma given directly (G*Power generic z test; no manual chapter)",
      if (identical(p$tails, "two")) {
        "two-sided test: alpha/2 in each tail of N(0, 1); the power counts both rejection tails under H1"
      } else {
        "one-sided test: the rejection region lies on the side of the sign of mu (upper tail for mu = 0)"
      }
    )
  )
}

# Generic F test ------------------------------------------------------------------------
#
# G*Power: Tests > Generic > Generic F; label "F tests - Generic F test" (no
# manual chapter). Upper-tail test: F_c = F^-1(1 - alpha; df1, df2), power
# P(F(df1, df2, lambda) > F_c). df1, df2 real (G*Power refuses df < 1);
# lambda >= 0. Noncentral F tails: Poisson mixtures of central pbeta terms of
# the requested tail (R/num-dist.R .ig_p_f()).

.ig_method_f_generic <- function() {
  ig_method(
    id = "f_generic",
    family = "F",
    label = list(en = "Generic F test", de = "Generischer F-Test"),
    params = ig_params(
      ig_param("lambda", role = "effect", domain = ig_domain_real(0, Inf, upper_inclusive = FALSE),
               label = list(en = "Noncentrality parameter \u03bb", de = "Nichtzentralit\u00e4tsparameter \u03bb"),
               description = "Noncentrality parameter of the H1 distribution F(df1, df2, lambda).",
               aliases = c("noncentrality_parameter", "ncp")),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      ig_param("df1", role = "design",
               domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE),
               unit = "degrees of freedom",
               label = list(en = "Numerator df", de = "Z\u00e4hler-Freiheitsgrade"),
               description = "Numerator degrees of freedom (positive real; G*Power refuses df < 1).",
               aliases = "numerator_df"),
      ig_param("df2", role = "design",
               domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE),
               unit = "degrees of freedom",
               label = list(en = "Denominator df", de = "Nenner-Freiheitsgrade"),
               description = "Denominator degrees of freedom (positive real; G*Power refuses df < 1).",
               aliases = "denominator_df")
    ),
    modes = .ig_generic_modes_g1,
    power = .ig_power_f_generic,
    design = NULL,
    discrete = FALSE,
    effect_direction = "increasing",
    engine = list(package = "stats",
                  fun = "qf; noncentral tails as Poisson mixtures of pbeta (inferguidR .ig_p_f)",
                  version = NA_character_),
    checker = "validation/independent/g1/check_g1.R",
    reference = list(
      gpower_path = "Tests > Generic > Generic F",
      gpower_label = "F tests - Generic F test",
      manual_pages = "no manual chapter (live inventory and emulator captures)",
      method_notes = "specs/method-notes/f_generic.md"
    ),
    status = "independently_checked"
  )
}

.ig_power_f_generic <- function(p) {
  df1 <- as.double(p$df1)
  df2 <- as.double(p$df2)
  lambda <- as.double(p$lambda)
  crit <- .ig_q_f(p$alpha, df1, df2, lower_tail = FALSE)
  power <- .ig_p_f(crit, df1, df2, lambda, lower_tail = FALSE)
  beta <- .ig_p_f(crit, df1, df2, lambda, lower_tail = TRUE)
  # Power and beta are complementary; the smaller one is the accurately
  # computed tail and the larger one is taken as its complement (as for the t
  # tests, .ig_t_power_g1()). .ig_p_f() evaluates the lower tail from the beta
  # argument y = df1 F / (df1 F + df2), whose complement loses its relative
  # accuracy when y is close to 1 (huge critical F, e.g. alpha 1e-5 with
  # df2 = 0.7): the lower tail was off by up to 1.3e-4 (power + beta =
  # 1.000127), while the upper tail, formed from the complement directly,
  # was accurate (review-G1-means-generic, mpmath cross-check; the root-cause
  # fix in R/num-dist.R is a change request).
  if (is.finite(power) && is.finite(beta)) {
    if (beta < power) power <- 1 - beta else beta <- 1 - power
  }
  list(
    power = power,
    beta = beta,
    alpha_attained = p$alpha,
    critical = crit,
    df = c(df1, df2),
    ncp = lambda,
    statistic_distribution = list(
      h0 = list(family = "F", df1 = df1, df2 = df2, ncp = 0),
      h1 = list(family = "F", df1 = df1, df2 = df2, ncp = lambda),
      rejection_side = "upper"
    ),
    assumptions = c(
      "H0: central F(df1, df2); H1: noncentral F(df1, df2, lambda), lambda and both df given directly",
      "upper-tail test: reject when F > F^-1(1 - alpha; df1, df2)"
    )
  )
}

# Generic chi-square test ---------------------------------------------------------------
#
# G*Power: Tests > Generic > Generic chi-square; label "chi-square tests -
# Generic chi-square test" (no manual chapter). Upper-tail test:
# c = chi2^-1(1 - alpha; df), power P(chi2(df, lambda) > c). df real
# (G*Power refuses df <= 1e-6); lambda >= 0.

.ig_method_chi2_generic <- function() {
  ig_method(
    id = "chi2_generic",
    family = "chi2",
    label = list(en = "Generic \u03c7\u00b2 test", de = "Generischer \u03c7\u00b2-Test"),
    params = ig_params(
      ig_param("lambda", role = "effect", domain = ig_domain_real(0, Inf, upper_inclusive = FALSE),
               label = list(en = "Noncentrality parameter \u03bb", de = "Nichtzentralit\u00e4tsparameter \u03bb"),
               description = "Noncentrality parameter of the H1 distribution chi2(df, lambda).",
               aliases = c("noncentrality_parameter", "ncp")),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      ig_param("df", role = "design",
               domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE),
               unit = "degrees of freedom",
               label = list(en = "Df", de = "Freiheitsgrade"),
               description = "Degrees of freedom (positive real; G*Power refuses df <= 1e-6).",
               aliases = "degrees_of_freedom")
    ),
    modes = .ig_generic_modes_g1,
    power = .ig_power_chi2_generic,
    design = NULL,
    discrete = FALSE,
    effect_direction = "increasing",
    engine = list(package = "stats",
                  fun = "qchisq; noncentral tails as Poisson mixtures of pchisq (inferguidR .ig_p_chisq)",
                  version = NA_character_),
    checker = "validation/independent/g1/check_g1.R",
    reference = list(
      gpower_path = "Tests > Generic > Generic chi-square",
      gpower_label = "\u03c7\u00b2 tests - Generic \u03c7\u00b2 test",
      manual_pages = "no manual chapter (live inventory and emulator captures)",
      method_notes = "specs/method-notes/chi2_generic.md"
    ),
    status = "independently_checked"
  )
}

.ig_power_chi2_generic <- function(p) {
  df <- as.double(p$df)
  lambda <- as.double(p$lambda)
  crit <- .ig_q_chisq(p$alpha, df, lower_tail = FALSE)
  power <- .ig_p_chisq(crit, df, lambda, lower_tail = FALSE)
  beta <- .ig_p_chisq(crit, df, lambda, lower_tail = TRUE)
  list(
    power = power,
    beta = beta,
    alpha_attained = p$alpha,
    critical = crit,
    df = df,
    ncp = lambda,
    statistic_distribution = list(
      h0 = list(family = "chisq", df = df, ncp = 0),
      h1 = list(family = "chisq", df = df, ncp = lambda),
      rejection_side = "upper"
    ),
    assumptions = c(
      "H0: central chi2(df); H1: noncentral chi2(df, lambda), lambda and df given directly",
      "upper-tail test: reject when X > chi2^-1(1 - alpha; df)"
    )
  )
}
