# Method definitions: correlation and regression (group G5) -------------------------
#
# Ownership: this file belongs to group G5 (blueprint phase 08, correlation and
# regression). It defines nine G*Power 3.1.9.7 procedures in three registry
# families through three provider functions:
#
#   .ig_define_methods_t_g5()      t_correlation_point_biserial,
#                                  t_linreg_one_group_slope,
#                                  t_linreg_two_groups_slopes,
#                                  t_linreg_two_groups_intercepts,
#                                  t_mrc_single_coefficient
#   .ig_define_methods_F_g5()      f_mrc_r2_deviation_from_zero, f_mrc_r2_increase
#   .ig_define_methods_exact_g5()  exact_correlation_bivariate_normal,
#                                  exact_mrc_random_model
#
# Method notes (mathematical contracts, conventions, evidence):
# specs/method-notes/<method_id>.md. Effect-size drawers:
# R/effect-correlation-regression.R. Independent checks:
# validation/independent/g5/. G*Power fixture comparisons:
# validation/compare/<method_id>/.
#
# Conventions shared by the nine procedures (all taken from the certified live
# inventory and the emulator captures; details in the method notes):
#
# * t tests: two tails put alpha/2 in each tail and count both rejection
#   tails; one tail rejects on the side of the noncentrality parameter
#   (.ig_power_t_test(), manual 25.5). F tests are one-sided (upper tail).
# * The exact tests (correlation, random-model R^2) reject in the H0
#   distribution of the statistic with alpha/2 per tail (two tails) or alpha
#   on the side of the effect (one tail).
# * G*Power refuses a target power below alpha ("Power must be greater than
#   or equal to alpha.") and an effect equal to the null value in a priori
#   mode; both are invalid inputs here (constraints).
# * Sensitivity mode with an "Effect direction" option (one-group slope,
#   exact correlation, random-model R^2) solves a Fisher-z distance from the
#   null value (`z_distance`) on the side the option selects; the solved H1
#   value is reported in method_outputs. The generic sensitivity scan is
#   centred at zero, so the natural parameter (whose null value can be
#   anywhere in its domain) could miss roots next to the null value.

.ig_define_methods_t_g5 <- function() {
  list(.ig_method_t_correlation_point_biserial(),
       .ig_method_t_linreg_one_group_slope(),
       .ig_method_t_linreg_two_groups_slopes(),
       .ig_method_t_linreg_two_groups_intercepts(),
       .ig_method_t_mrc_single_coefficient())
}

.ig_define_methods_F_g5 <- function() {
  list(.ig_method_f_mrc_r2_deviation_from_zero(),
       .ig_method_f_mrc_r2_increase())
}

.ig_define_methods_exact_g5 <- function() {
  list(.ig_method_exact_correlation_bivariate_normal(),
       .ig_method_exact_mrc_random_model())
}

# Shared building blocks ----------------------------------------------------------------

.ig_g5_modes <- c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise")

.ig_g5_option_tails <- function() {
  ig_option("tails", values = c("one", "two"), default = "one", labels = c(one = "One", two = "Two"),
            label = list(en = "Tail(s)", de = "Seiten"))
}

.ig_g5_n_total <- function(lower, description) {
  ig_param("n_total", role = "design", domain = ig_domain_integer(lower = lower), unit = "subjects",
           label = list(en = "Total sample size", de = "Gesamtstichprobenumfang"),
           description = description)
}

.ig_g5_count <- function(name, label, description, lower = 1) {
  ig_param(name, role = "design", domain = ig_domain_integer(lower = lower), unit = "predictors",
           label = list(en = label), description = description)
}

.ig_g5_positive <- function(name, label, description, role = "nuisance") {
  ig_param(name, role = role, domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE),
           label = list(en = label), description = description)
}

# z_distance: the sensitivity solve target of the procedures with an effect
# direction option (see the file header).
.ig_g5_z_distance <- function(description) {
  ig_param("z_distance", role = "effect", domain = ig_domain_real(0, Inf, upper_inclusive = FALSE),
           label = list(en = "Fisher-z distance from H0 (sensitivity)"), description = description,
           modes = "sensitivity")
}

.ig_g5_effect_direction_sign <- function(p, negative_code) {
  if (identical(p$effect_direction, negative_code)) -1 else 1
}

# The specification holds `power` and `alpha` together only in a priori and
# sensitivity mode; an effect value is present in every mode but sensitivity.
.ig_g5_power_alpha_problem <- function(p) {
  if (!is.null(p$alpha) && !is.null(p$power) && p$power < p$alpha) {
    sprintf("the target power %s is below alpha %s; it must be at least alpha (G*Power: 'Power must be greater than or equal to alpha.')",
            format(p$power), format(p$alpha))
  } else {
    character()
  }
}

.ig_g5_is_a_priori <- function(p, effect) {
  !is.null(p$alpha) && !is.null(p$power) && !is.null(p[[effect]]) && is.null(p$n_total) && is.null(p$n1)
}

.ig_g5_reference <- function(path, label, manual, notes, extra = list()) {
  c(list(gpower_path = path, gpower_label = label, manual_pages = manual,
         method_notes = notes, reference_status = "reference_pending (emulator captures compared in validation/compare/)"),
    extra)
}

.ig_g5_t_assumption <- function(tails) {
  if (identical(tails, "two")) {
    "two-sided t test: alpha/2 in each tail; the power counts both rejection tails (G*Power convention, manual 25.5)"
  } else {
    "one-sided t test: the rejection region lies on the side of the noncentrality parameter (manual 25.5); ncp = 0 uses the upper tail"
  }
}

# Common output list of a t-based forward model.
.ig_g5_t_result <- function(df, ncp, alpha, tails, n_total, group_sizes = NULL, outputs = list(),
                            assumptions = character()) {
  core <- .ig_power_t_test(df, ncp, alpha, tails)
  # nmath's noncentral t forms some upper tails as 1 - p and can return
  # 1 + 2e-11 at df ~ 1e5 (for example df = 99996, ncp = 15.8, t = 3.09; the
  # same defect is triaged for t_means_two_independent). Noise of at most
  # 1e-9 outside [0, 1] is clamped here; anything larger is left to the
  # solver's range check. Root cause: R/num-dist.R (change request).
  clamp <- function(v) if (!is.na(v) && v > 1 && v <= 1 + 1e-9) 1 else if (!is.na(v) && v < 0 && v >= -1e-9) 0 else v
  core$power <- clamp(core$power)
  core$beta <- clamp(core$beta)
  out <- list(
    power = core$power, beta = core$beta, alpha_attained = alpha, critical = core$critical,
    df = df, ncp = ncp, n_total = n_total,
    statistic_distribution = list(h0 = list(family = "t", df = df, ncp = 0),
                                  h1 = list(family = "t", df = df, ncp = ncp),
                                  rejection_side = core$side),
    outputs = outputs,
    assumptions = c(assumptions, .ig_g5_t_assumption(tails))
  )
  if (!is.null(group_sizes)) out$group_sizes <- group_sizes
  out
}

# Noncentral F tails with a rigorous shortcut for very large noncentrality.
# .ig_p_f sums a Poisson mixture of about 10 sqrt(lambda) terms, which takes
# seconds for lambda ~ 1e9 (values reached by the sensitivity scan). With
# F = (X1/df1) / (X2/df2), X1 = (Z + sqrt(lambda))^2 + chi^2_(df1-1) and any s > 0:
#   P(F <= q) <= P(X1 <= q df1 s / df2) + P(X2 > s) <= Phi(sqrt(q df1 s / df2) - sqrt(lambda)) + P(X2 > s).
# When this bound is below 1e-300 the lower tail is 0 and the upper tail 1 to
# double precision; otherwise .ig_p_f is used. Root cause and general fix:
# R/num-dist.R (change request).
.ig_g5_p_f <- function(q, df1, df2, lambda, lower_tail) {
  if (lambda > 0 && q > 0 && is.finite(q)) {
    s <- df2 + 60 * sqrt(2 * df2) + 1500
    x <- q * df1 * s / df2
    if (sqrt(lambda) > sqrt(x) + 38) {
      lb <- .ig_logspace_add(stats::pnorm(sqrt(x) - sqrt(lambda), log.p = TRUE),
                             stats::pchisq(s, df2, lower.tail = FALSE, log.p = TRUE))
      if (lb < log(1e-300)) return(if (lower_tail) 0 else 1)
    }
  }
  .ig_p_f(q, df1, df2, lambda, lower_tail = lower_tail)
}

# Common output list of an F-based forward model (upper-tail test).
.ig_g5_f_result <- function(df1, df2, lambda, alpha, n_total, outputs = list(), assumptions = character()) {
  crit <- .ig_q_f(alpha, df1, df2, lower_tail = FALSE)
  list(
    power = .ig_g5_p_f(crit, df1, df2, lambda, lower_tail = FALSE),
    beta = .ig_g5_p_f(crit, df1, df2, lambda, lower_tail = TRUE),
    alpha_attained = alpha, critical = crit, df = c(df1, df2), ncp = lambda, n_total = n_total,
    statistic_distribution = list(h0 = list(family = "F", df1 = df1, df2 = df2, ncp = 0),
                                  h1 = list(family = "F", df1 = df1, df2 = df2, ncp = lambda),
                                  rejection_side = "upper"),
    outputs = c(list(numerator_df = df1, denominator_df = df2), outputs),
    assumptions = c(assumptions, "one-sided F test (upper tail), as in G*Power")
  )
}

# Single-sample a priori design: every N with at least `min_n` subjects
# (NULL below, so the generic integer search skips those totals).
.ig_g5_single_design <- function(min_n_fun) {
  function(p, n_total) {
    if (n_total < min_n_fun(p)) return(NULL)
    c(n_total = n_total)
  }
}

# =====================================================================================
# t family
# =====================================================================================

# Correlation: Point biserial model ---------------------------------------------------
#
# Manual ch. 16 (PDF pp. 40-41): t test of H0 rho = 0 for the point-biserial
# correlation, delta = sqrt(rho^2 N / (1 - rho^2)), df = N - 2 (16.5).
# G*Power's effect field is the magnitude |rho| (admissible range
# (0, 0.999999); 1e-6 in a priori mode).

.ig_method_t_correlation_point_biserial <- function() {
  ig_method(
    id = "t_correlation_point_biserial", family = "t",
    label = list(en = "Correlation: Point biserial model", de = "Korrelation: Punktbiseriales Modell"),
    params = ig_params(
      ig_param("rho", role = "effect", domain = ig_domain_real(0, 1, upper_inclusive = FALSE),
               label = list(en = "Effect size |\u03c1|"),
               description = "Magnitude of the point-biserial correlation under H1 (H0: rho = 0).",
               aliases = "effect_size_abs_rho"),
      ig_param_alpha(), ig_param_power(), ig_param_error_ratio(),
      .ig_g5_n_total(3, "Total sample size N (df = N - 2).")
    ),
    options = ig_options(.ig_g5_option_tails()),
    modes = .ig_g5_modes,
    power = .ig_power_t_correlation_point_biserial,
    design = NULL, discrete = FALSE, effect_direction = "increasing",
    engine = list(package = "stats", fun = "pt/qt via .ig_power_t_test()", version = NA_character_),
    checker = "validation/independent/g5/check_g5_t.R",
    reference = .ig_g5_reference("Tests > Correlation and regression > Correlation: Point biserial model",
                                 "t tests - Correlation: Point biserial model",
                                 "ch. 16 (PDF pp. 40-41)", "specs/method-notes/t_correlation_point_biserial.md"),
    status = "independently_checked",
    constraints = function(p) {
      c(.ig_g5_power_alpha_problem(p),
        if (.ig_g5_is_a_priori(p, "rho") && p$rho == 0)
          "in a priori mode |rho| must be greater than 0 (with rho = 0 the power equals alpha for every N)")
    }
  )
}

.ig_power_t_correlation_point_biserial <- function(p) {
  N <- as.double(p$n_total)
  rho <- p$rho
  ncp <- sqrt(N) * rho / sqrt(1 - rho^2)
  .ig_g5_t_result(N - 2, ncp, p$alpha, p$tails, N,
                  assumptions = c("binary predictor and normally distributed outcome within groups (point-biserial model)",
                                  "df = N - 2; ncp = sqrt(N rho^2 / (1 - rho^2)) (manual 16.5)"))
}

# Linear bivariate regression: One group, size of slope ---------------------------------
#
# Manual ch. 12 (PDF pp. 31-32). H0 b = b0. The noncentrality uses the
# residual SD under H1, sigma_y sqrt(1 - rho1^2) with rho1 = b1 sigma_x /
# sigma_y: delta = sqrt(N) sigma_x (b1 - b0) / (sigma_y sqrt(1 - rho1^2)),
# df = N - 2. The manual prints sigma_y in place of the residual SD (12.5); the
# printed example DOC-073 and the captures (including Slope H0 != 0,
# SUP-G5-T-SLOPE1-POS-P01/P02) follow the residual SD.

.ig_method_t_linreg_one_group_slope <- function() {
  ig_method(
    id = "t_linreg_one_group_slope", family = "t",
    label = list(en = "Linear bivariate regression: One group, size of slope", de = "Lineare bivariate Regression: Steigung in einer Gruppe"),
    params = ig_params(
      ig_param("slope_h1", role = "effect", domain = ig_domain_real(), label = list(en = "Slope H1"),
               description = "Slope b under H1; |slope_h1 sd_x / sd_y| < 1.",
               modes = c("a_priori", "post_hoc", "criterion", "compromise")),
      .ig_g5_z_distance(paste(
        "Sensitivity solve target: |atanh(rho1) - atanh(rho0)| with rho = slope sd_x / sd_y; the effect",
        "direction option gives the side. The solved Slope H1 is method_outputs$slope_h1.")),
      ig_param("slope_h0", role = "nuisance", domain = ig_domain_real(), default = 0,
               label = list(en = "Slope H0"),
               description = "Slope b0 under H0. Default 0 (G*Power's live field default)."),
      .ig_g5_positive("sd_x", "Std dev \u03c3_x", "SD of the X values (normalised by N, manual 12.1)."),
      .ig_g5_positive("sd_y", "Std dev \u03c3_y", "SD of the Y values."),
      ig_param_alpha(), ig_param_power(), ig_param_error_ratio(),
      .ig_g5_n_total(3, "Total sample size N (df = N - 2).")
    ),
    options = ig_options(
      .ig_g5_option_tails(),
      ig_option("effect_direction", values = c("slope_lt_h0", "slope_gt_h0"), default = "slope_gt_h0",
                labels = c(slope_lt_h0 = "Slope < H0", slope_gt_h0 = "Slope > H0"),
                applies_when = list(mode = "sensitivity"), label = list(en = "Effect direction"))
    ),
    modes = .ig_g5_modes,
    targets = list(sensitivity = "z_distance"),
    power = .ig_power_t_linreg_one_group_slope,
    design = NULL, discrete = FALSE, effect_direction = "increasing",
    engine = list(package = "stats", fun = "pt/qt via .ig_power_t_test()", version = NA_character_),
    checker = "validation/independent/g5/check_g5_t.R",
    reference = .ig_g5_reference("Tests > Correlation and regression > Linear bivariate regression: One group, size of slope",
                                 "t tests - Linear bivariate regression: One group, size of slope",
                                 "ch. 12 (PDF pp. 31-32)", "specs/method-notes/t_linreg_one_group_slope.md"),
    status = "independently_checked",
    constraints = function(p) {
      out <- .ig_g5_power_alpha_problem(p)
      k <- p$sd_x / p$sd_y
      if (!is.null(p$slope_h1) && abs(p$slope_h1 * k) >= 1) {
        out <- c(out, "|slope_h1 * sd_x / sd_y| must be below 1 (it is the correlation implied under H1)")
      }
      if (is.null(p$slope_h1) && abs(p$slope_h0 * k) >= 1) {
        out <- c(out, "in sensitivity mode |slope_h0 * sd_x / sd_y| must be below 1")
      }
      if (.ig_g5_is_a_priori(p, "slope_h1") && p$slope_h1 == p$slope_h0) {
        out <- c(out, "in a priori mode slope_h1 must differ from slope_h0 (otherwise the power equals alpha for every N)")
      }
      out
    }
  )
}

.ig_power_t_linreg_one_group_slope <- function(p) {
  N <- as.double(p$n_total)
  k <- p$sd_x / p$sd_y
  b1 <- if (!is.null(p$slope_h1)) p$slope_h1 else {
    rho0 <- p$slope_h0 * k
    tanh(atanh(rho0) + .ig_g5_effect_direction_sign(p, "slope_lt_h0") * p$z_distance) / k
  }
  rho1 <- b1 * k
  if (abs(rho1) >= 1) {
    # Only reachable from a sensitivity distance so large that tanh() rounds to +-1.
    side_power <- 1
    return(list(power = side_power, beta = 0, alpha_attained = p$alpha,
                critical = .ig_power_t_test(N - 2, 0, p$alpha, p$tails)$critical,
                df = N - 2, ncp = sign(rho1) * Inf, n_total = N,
                outputs = list(slope_h1 = b1),
                notes = list(ncp = "infinite: the implied correlation is +-1 to double precision"),
                assumptions = "limit |rho1| -> 1: the residual SD is 0 and the power is 1"))
  }
  ncp <- sqrt(N) * p$sd_x * (b1 - p$slope_h0) / (p$sd_y * sqrt(1 - rho1^2))
  .ig_g5_t_result(N - 2, ncp, p$alpha, p$tails, N,
                  outputs = list(slope_h1 = b1, residual_sd_h1 = p$sd_y * sqrt(1 - rho1^2), rho_h1 = rho1),
                  assumptions = c("X fixed with SD sd_x (normalised by N); normal residuals with constant SD",
                                  "ncp = sqrt(N) sd_x (b1 - b0) / (sd_y sqrt(1 - rho1^2)), rho1 = b1 sd_x / sd_y (residual SD under H1); df = N - 2"))
}

# Linear bivariate regression: Two groups (slopes and intercepts) ----------------------
#
# Slopes: manual ch. 17 (= 18, a duplicate), PDF pp. 42-47; intercepts: no
# manual chapter (app only). With integer group sizes n1, n2:
#
#   slopes:      delta = |D| / (sigma sqrt(1/(n1 sx1^2) + 1/(n2 sx2^2)))
#   intercepts:  delta = |D| / (sigma sqrt((1 + mx1^2/sx1^2)/n1 + (1 + mx2^2/sx2^2)/n2))
#   df = n1 + n2 - 4
#
# (SDs of X normalised by n, manual 17.3). Post hoc captures with unequal
# groups follow these formulas (SUP-G5-T-SLOPE2-POS-P01/P02,
# SUP-G5-T-INTER2-POS-P01). In a priori mode G*Power displays a
# noncentrality computed with n1 replaced by n2 / r (the nominal allocation
# ratio r = N2/N1), which is not the noncentrality of the displayed design;
# inferguidR reports the values of the integer design (method notes,
# discrepancy triage in validation/compare/).
#
# A priori search (validation/compare/t_linreg_two_groups_slopes/allocation-rule.md):
#   1. N0 = smallest even N >= 6 whose power with fractional sizes
#      N/(1+r), rN/(1+r) and df = N - 4 reaches the target;
#   2. slopes only: N = smallest N >= N0 whose design n1 = floor(N/(1+r) + 1/2),
#      n2 = N - n1 reaches the target under G*Power's displayed
#      noncentrality (n2 with n1 replaced by n2/r, df = N - 4); intercepts: N = N0;
#   3. the design is n1 = floor(N/(1+r) + 1/2), n2 = N - n1 (each group >= 2).
# If the power of that integer design is below the target, the next totals
# are evaluated until one reaches it (G*Power reports such designs; inferguidR
# does not); the G*Power total is kept in the diagnostics.

.ig_g5_two_group_params <- function(effect_name, effect_label, effect_description, intercepts) {
  base <- list(
    ig_param(effect_name, role = "effect", domain = ig_domain_real(0, Inf, upper_inclusive = FALSE),
             label = list(en = effect_label), description = effect_description),
    ig_param_alpha(), ig_param_power(), ig_param_error_ratio(),
    ig_param("n_total", role = "design", domain = ig_domain_integer(lower = 5), unit = "subjects",
             label = list(en = "Total sample size", de = "Gesamtstichprobenumfang"),
             description = "Solved in a priori mode (n1 + n2; df = N - 4).", modes = "a_priori"),
    ig_param("n1", role = "design", domain = ig_domain_integer(lower = 2), unit = "subjects",
             label = list(en = "Sample size group 1"), aliases = "sample_size_group_1",
             description = "Group 1 size (at least 2 observations per regression line).",
             modes = c("post_hoc", "sensitivity", "criterion", "compromise")),
    ig_param("n2", role = "design", domain = ig_domain_integer(lower = 2), unit = "subjects",
             label = list(en = "Sample size group 2"), aliases = "sample_size_group_2",
             description = "Group 2 size (at least 2 observations per regression line).",
             modes = c("post_hoc", "sensitivity", "criterion", "compromise")),
    ig_param("allocation_ratio", role = "design",
             domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE), default = 1,
             label = list(en = "Allocation ratio N2/N1"), aliases = "allocation_ratio_n2_n1",
             description = "Group-size ratio n2/n1 of the a priori design (G*Power's live field default 1).",
             modes = "a_priori"),
    .ig_g5_positive("sd_residual", "Std dev residual \u03c3", "Common SD of the residuals (pooled over the groups)."),
    .ig_g5_positive("sd_x1", "Std dev \u03c3_x1", "SD of X in group 1 (normalised by n1)."),
    .ig_g5_positive("sd_x2", "Std dev \u03c3_x2", "SD of X in group 2 (normalised by n2).")
  )
  if (intercepts) {
    base <- c(base, list(
      ig_param("mean_x1", role = "nuisance", domain = ig_domain_real(), label = list(en = "Mean \u03bc_x1"),
               description = "Mean of X in group 1."),
      ig_param("mean_x2", role = "nuisance", domain = ig_domain_real(), label = list(en = "Mean \u03bc_x2"),
               description = "Mean of X in group 2.")))
  }
  do.call(ig_params, base)
}

.ig_g5_two_group_constraints <- function(effect) {
  function(p) {
    out <- .ig_g5_power_alpha_problem(p)
    if (!is.null(p$n1) && !is.null(p$n2) && as.double(p$n1) + as.double(p$n2) < 5) {
      out <- c(out, "n1 + n2 must be at least 5 (df = n1 + n2 - 4 >= 1)")
    }
    if (!is.null(p$alpha) && !is.null(p$power) && !is.null(p[[effect]]) && is.null(p$n1) && p[[effect]] == 0) {
      out <- c(out, sprintf("in a priori mode %s must be greater than 0 (otherwise the power equals alpha for every N)", effect))
    }
    out
  }
}

# Variance factor v(n1, n2) with delta = |D| / (sigma sqrt(v)).
.ig_g5_slopes_var <- function(p, n1, n2) 1 / (n1 * p$sd_x1^2) + 1 / (n2 * p$sd_x2^2)
.ig_g5_intercepts_var <- function(p, n1, n2) {
  (1 + p$mean_x1^2 / p$sd_x1^2) / n1 + (1 + p$mean_x2^2 / p$sd_x2^2) / n2
}

.ig_g5_two_group_power <- function(var_fun, effect, kind) {
  function(p) {
    n1 <- as.double(p$n1)
    n2 <- as.double(p$n2)
    df <- n1 + n2 - 4
    ncp <- p[[effect]] / (p$sd_residual * sqrt(var_fun(p, n1, n2)))
    formula <- if (kind == "slopes") {
      "ncp = |D slope| / (sigma sqrt(1/(n1 sx1^2) + 1/(n2 sx2^2))), df = n1 + n2 - 4 (Dupont & Plummer 1998 with a noncentral t, manual 17.5)"
    } else {
      "ncp = |D intercept| / (sigma sqrt((1 + mx1^2/sx1^2)/n1 + (1 + mx2^2/sx2^2)/n2)), df = n1 + n2 - 4 (separate regressions; no manual chapter)"
    }
    .ig_g5_t_result(df, ncp, p$alpha, p$tails, n1 + n2, group_sizes = c(n1 = n1, n2 = n2),
                    assumptions = c("X fixed in each group (SDs normalised by n); normal residuals with a common SD", formula))
  }
}

# G*Power's integer split and the a priori design rule.
.ig_g5_two_group_split <- function(N, r) {
  n1 <- floor(N / (1 + r) + 0.5)
  c(n1 = n1, n2 = N - n1)
}

.ig_g5_two_group_design <- function(p, n_total, continuous = FALSE) {
  r <- p$allocation_ratio
  if (continuous) return(c(n1 = n_total / (1 + r), n2 = r * n_total / (1 + r)))
  gs <- .ig_g5_two_group_split(n_total, r)
  if (n_total < 5 || any(gs < 2)) return(NULL)
  gs
}

.ig_g5_two_group_a_priori_search <- function(var_fun, effect, kind) {
  function(p, n_lo, n_max) {
    r <- p$allocation_ratio
    target <- p$power
    pw <- function(n1, n2, df) {
      ncp <- p[[effect]] / (p$sd_residual * sqrt(var_fun(p, n1, n2)))
      .ig_power_t_test(df, ncp, p$alpha, p$tails)$power
    }
    frac_ok <- function(N) isTRUE(pw(N / (1 + r), r * N / (1 + r), N - 4) >= target)
    evals <- 0L
    # Step 1: smallest even N >= 6 with fractional power >= target (power is
    # increasing in N here, so doubling and bisection over even N are exact).
    lo <- 4
    hi <- 6
    while (!frac_ok(hi)) {
      evals <- evals + 1L
      lo <- hi
      hi <- hi * 2
      if (hi > n_max) {
        return(list(status = "infeasible", reason = sprintf(
          "the target power %s is not reached by G*Power's fractional criterion for N <= %s",
          format(target), format(n_max, scientific = FALSE)), diagnostics = list(evaluations = evals)))
      }
    }
    while (hi - lo > 2) {
      mid <- lo + 2 * floor((hi - lo) / 4)
      evals <- evals + 1L
      if (frac_ok(mid)) hi <- mid else lo <- mid
    }
    n_frac <- hi
    n_gp <- n_frac
    # Step 2 (slopes): G*Power's check of its displayed design.
    if (kind == "slopes") {
      disp_ok <- function(N) {
        gs <- .ig_g5_two_group_split(N, r)
        gs[["n2"]] >= 1 && isTRUE(pw(gs[["n2"]] / r, gs[["n2"]], N - 4) >= target)
      }
      while (!disp_ok(n_gp)) {
        n_gp <- n_gp + 1
        evals <- evals + 1L
        if (n_gp > n_max) {
          return(list(status = "infeasible", reason = "G*Power's design check does not terminate below the search limit",
                      diagnostics = list(n_fractional_even = n_frac, evaluations = evals)))
        }
      }
    }
    diag <- list(
      rule = if (kind == "slopes") {
        "G*Power 3.1.9.7: smallest even N with fractional power >= target, then the smallest N whose rounded design passes G*Power's displayed-noncentrality check"
      } else {
        "G*Power 3.1.9.7: smallest even N with fractional power >= target (no check of the rounded design)"
      },
      n_fractional_even = n_frac, n_gpower_rule = n_gp,
      split = "n1 = floor(N / (1 + r) + 1/2), n2 = N - n1", evaluations = evals)
    gs <- .ig_g5_two_group_split(n_gp, r)
    if (any(gs < 2)) {
      return(list(status = "infeasible", reason = sprintf(paste(
        "G*Power's allocation rule gives n1 = %s, n2 = %s at N = %s; each regression line needs at least 2",
        "observations (G*Power reports such a design); use an allocation ratio closer to 1"),
        format(gs[["n1"]]), format(gs[["n2"]]), format(n_gp, scientific = FALSE)), diagnostics = diag))
    }
    N <- n_gp
    p_int <- pw(gs[["n1"]], gs[["n2"]], N - 4)
    diag$gpower_design_power <- p_int
    if (p_int < target) {
      # The G*Power design misses the target; report the next design that
      # reaches it (the difference is listed in the diagnostics).
      repeat {
        N <- N + 1
        evals <- evals + 1L
        if (N > n_max) {
          return(list(status = "infeasible", reason = "no design reaching the target above the G*Power total",
                      diagnostics = diag))
        }
        g2 <- .ig_g5_two_group_split(N, r)
        if (any(g2 < 2)) next
        p_int <- pw(g2[["n1"]], g2[["n2"]], N - 4)
        if (p_int >= target) break
      }
      diag$note <- sprintf(paste(
        "G*Power's total N = %s (n1 = %s, n2 = %s) has power %s below the target %s for that integer design;",
        "the next total reaching the target is reported"), format(n_gp, scientific = FALSE),
        format(gs[["n1"]]), format(gs[["n2"]]), format(signif(diag$gpower_design_power, 7)), format(target))
      diag$evaluations <- evals
    }
    list(status = "ok", n = N, power = p_int, diagnostics = diag)
  }
}

.ig_g5_two_group_method <- function(kind) {
  slopes <- kind == "slopes"
  id <- if (slopes) "t_linreg_two_groups_slopes" else "t_linreg_two_groups_intercepts"
  effect <- if (slopes) "delta_slope" else "delta_intercept"
  var_fun <- if (slopes) .ig_g5_slopes_var else .ig_g5_intercepts_var
  label <- if (slopes) "Linear bivariate regression: Two groups, difference between slopes" else
    "Linear bivariate regression: Two groups, difference between intercepts"
  ig_method(
    id = id, family = "t", label = list(en = label, de = if (slopes) "Lineare bivariate Regression: Unterschied zweier Steigungen" else "Lineare bivariate Regression: Unterschied zweier Achsenabschnitte"),
    params = .ig_g5_two_group_params(
      effect, if (slopes) "|\u0394 slope|" else "|\u0394 intercept|",
      if (slopes) "Absolute difference |b1 - b2| between the slopes (H0: b1 = b2)." else
        "Absolute difference |a1 - a2| between the intercepts (H0: a1 = a2).",
      intercepts = !slopes),
    options = ig_options(.ig_g5_option_tails()),
    modes = .ig_g5_modes,
    power = .ig_g5_two_group_power(var_fun, effect, kind),
    design = .ig_g5_two_group_design,
    discrete = FALSE, effect_direction = "increasing",
    engine = list(package = "stats", fun = "pt/qt via .ig_power_t_test()", version = NA_character_),
    checker = "validation/independent/g5/check_g5_t.R",
    reference = .ig_g5_reference(
      paste("Tests > Correlation and regression >", label), paste("t tests -", label),
      if (slopes) "ch. 17 (= 18), PDF pp. 42-47" else "no manual chapter (app only)",
      sprintf("specs/method-notes/%s.md", id),
      list(allocation_rule = paste(
        if (slopes) "G*Power 3.1.9.7 slope rule:" else "G*Power 3.1.9.7 intercept rule:",
        "smallest even N with fractional power >= target",
        if (slopes) "then the smallest N whose rounded design passes G*Power's displayed-noncentrality check;" else ";",
        "n1 = floor(N/(1+r) + 1/2), n2 = N - n1; reported values are those of the integer design",
        "(validation/compare/t_linreg_two_groups_slopes/allocation-rule.md)"))),
    status = "independently_checked",
    constraints = .ig_g5_two_group_constraints(effect),
    a_priori_search = .ig_g5_two_group_a_priori_search(var_fun, effect, kind)
  )
}

.ig_method_t_linreg_two_groups_slopes <- function() .ig_g5_two_group_method("slopes")
.ig_method_t_linreg_two_groups_intercepts <- function() .ig_g5_two_group_method("intercepts")

# Linear multiple regression: Fixed model, single regression coefficient --------------
#
# No manual chapter (app only). t test of one coefficient in a model with p
# predictors: delta = sqrt(f^2 N), df = N - p - 1 (f^2 the partial effect of the
# predictor, drawer "Variance explained by predictor / Residual variance").

.ig_method_t_mrc_single_coefficient <- function() {
  ig_method(
    id = "t_mrc_single_coefficient", family = "t",
    label = list(en = "Linear multiple regression: Fixed model, single regression coefficient", de = "Multiple lineare Regression: Festes Modell, einzelner Regressionskoeffizient"),
    params = ig_params(
      ig_param("f2", role = "effect", domain = ig_domain_real(0, Inf, upper_inclusive = FALSE),
               label = list(en = "Effect size f\u00b2"),
               description = "Partial f^2 of the tested predictor (variance it explains / residual variance)."),
      ig_param_alpha(), ig_param_power(), ig_param_error_ratio(),
      .ig_g5_n_total(3, "Total sample size N (df = N - p - 1, so N >= p + 2)."),
      .ig_g5_count("n_predictors", "Number of predictors", "Total number p of predictors in the model.")
    ),
    options = ig_options(.ig_g5_option_tails()),
    modes = .ig_g5_modes,
    power = .ig_power_t_mrc_single_coefficient,
    design = .ig_g5_single_design(function(p) p$n_predictors + 2),
    discrete = FALSE, effect_direction = "increasing",
    engine = list(package = "stats", fun = "pt/qt via .ig_power_t_test()", version = NA_character_),
    checker = "validation/independent/g5/check_g5_t.R",
    reference = .ig_g5_reference(
      "Tests > Correlation and regression > Linear multiple regression: Fixed model, single regression coefficient",
      "t tests - Linear multiple regression: Fixed model, single regression coefficient",
      "no manual chapter (app only)", "specs/method-notes/t_mrc_single_coefficient.md"),
    status = "independently_checked",
    constraints = function(p) {
      out <- .ig_g5_power_alpha_problem(p)
      if (!is.null(p$n_total) && p$n_total < p$n_predictors + 2) {
        out <- c(out, sprintf("N must be at least the number of predictors + 2 (%s)", format(p$n_predictors + 2)))
      }
      if (.ig_g5_is_a_priori(p, "f2") && p$f2 == 0) {
        out <- c(out, "in a priori mode f2 must be greater than 0 (otherwise the power equals alpha for every N)")
      }
      out
    }
  )
}

.ig_power_t_mrc_single_coefficient <- function(p) {
  N <- as.double(p$n_total)
  .ig_g5_t_result(N - p$n_predictors - 1, sqrt(p$f2 * N), p$alpha, p$tails, N,
                  assumptions = c("fixed predictors; normal residuals with constant variance",
                                  "ncp = sqrt(f2 N), df = N - p - 1 (t test of one regression coefficient)"))
}

# =====================================================================================
# F family
# =====================================================================================
#
# Manual ch. 13 (R^2 deviation from zero) and ch. 14 (R^2 increase): fixed
# predictors, lambda = f^2 N (13.5, 14.5); df1 = m (resp. q tested
# predictors), df2 = N - m - 1 with m the total number of predictors.

.ig_method_f_mrc_r2_deviation_from_zero <- function() {
  ig_method(
    id = "f_mrc_r2_deviation_from_zero", family = "F",
    label = list(en = "Linear multiple regression: Fixed model, R\u00b2 deviation from zero", de = "Multiple lineare Regression: Festes Modell, Abweichung von R\u00b2 von null"),
    params = ig_params(
      ig_param("f2", role = "effect", domain = ig_domain_real(0, Inf, upper_inclusive = FALSE),
               label = list(en = "Effect size f\u00b2"), description = "f^2 = R^2 / (1 - R^2) (manual 13.1)."),
      ig_param_alpha(), ig_param_power(), ig_param_error_ratio(),
      .ig_g5_n_total(3, "Total sample size N (df2 = N - m - 1, so N >= m + 2)."),
      .ig_g5_count("n_predictors", "Number of predictors", "Number m of predictors (df1 = m).")
    ),
    modes = .ig_g5_modes,
    power = function(p) {
      N <- as.double(p$n_total)
      .ig_g5_f_result(p$n_predictors, N - p$n_predictors - 1, p$f2 * N, p$alpha, N,
                      assumptions = c("fixed predictors; normal residuals", "lambda = f2 N, df1 = m, df2 = N - m - 1 (manual 13.5)"))
    },
    design = .ig_g5_single_design(function(p) p$n_predictors + 2),
    discrete = FALSE, effect_direction = "increasing",
    engine = list(package = "stats + inferguidR", fun = "qf; noncentral F tail as Poisson mixture of pbeta (.ig_p_f)",
                  version = NA_character_),
    checker = "validation/independent/g5/check_g5_f.R",
    reference = .ig_g5_reference(
      "Tests > Correlation and regression > Linear multiple regression: Fixed model, R\u00b2 deviation from zero",
      "F tests - Linear multiple regression: Fixed model, R\u00b2 deviation from zero",
      "ch. 13 (PDF pp. 33-35)", "specs/method-notes/f_mrc_r2_deviation_from_zero.md"),
    status = "independently_checked",
    constraints = function(p) {
      out <- .ig_g5_power_alpha_problem(p)
      if (!is.null(p$n_total) && p$n_total < p$n_predictors + 2) {
        out <- c(out, sprintf("N must be at least the number of predictors + 2 (%s)", format(p$n_predictors + 2)))
      }
      if (.ig_g5_is_a_priori(p, "f2") && p$f2 == 0) {
        out <- c(out, "in a priori mode f2 must be greater than 0 (otherwise the power equals alpha for every N)")
      }
      out
    }
  )
}

.ig_method_f_mrc_r2_increase <- function() {
  ig_method(
    id = "f_mrc_r2_increase", family = "F",
    label = list(en = "Linear multiple regression: Fixed model, R\u00b2 increase", de = "Multiple lineare Regression: Festes Modell, R\u00b2-Zuwachs"),
    params = ig_params(
      ig_param("f2", role = "effect", domain = ig_domain_real(0, Inf, upper_inclusive = FALSE),
               label = list(en = "Effect size f\u00b2"),
               description = "f^2 of the tested predictors: (R^2_AB - R^2_A) / (1 - R^2_AB) (manual 14.1)."),
      ig_param_alpha(), ig_param_power(), ig_param_error_ratio(),
      .ig_g5_n_total(3, "Total sample size N (df2 = N - m - 1, so N >= m + 2)."),
      .ig_g5_count("n_tested", "Number of tested predictors", "Number q of tested predictors (df1 = q)."),
      .ig_g5_count("n_predictors", "Total number of predictors", "Total number m of predictors (df2 = N - m - 1).")
    ),
    modes = .ig_g5_modes,
    power = function(p) {
      N <- as.double(p$n_total)
      .ig_g5_f_result(p$n_tested, N - p$n_predictors - 1, p$f2 * N, p$alpha, N,
                      assumptions = c("fixed predictors; normal residuals", "lambda = f2 N, df1 = q, df2 = N - m - 1 (manual 14.5)"))
    },
    design = .ig_g5_single_design(function(p) p$n_predictors + 2),
    discrete = FALSE, effect_direction = "increasing",
    engine = list(package = "stats + inferguidR", fun = "qf; noncentral F tail as Poisson mixture of pbeta (.ig_p_f)",
                  version = NA_character_),
    checker = "validation/independent/g5/check_g5_f.R",
    reference = .ig_g5_reference(
      "Tests > Correlation and regression > Linear multiple regression: Fixed model, R\u00b2 increase",
      "F tests - Linear multiple regression: Fixed model, R\u00b2 increase",
      "ch. 14 (PDF pp. 36-38)", "specs/method-notes/f_mrc_r2_increase.md"),
    status = "independently_checked",
    constraints = function(p) {
      out <- .ig_g5_power_alpha_problem(p)
      if (p$n_tested > p$n_predictors) {
        out <- c(out, "the number of tested predictors must not exceed the total number of predictors")
      }
      if (!is.null(p$n_total) && p$n_total < p$n_predictors + 2) {
        out <- c(out, sprintf("N must be at least the total number of predictors + 2 (%s)", format(p$n_predictors + 2)))
      }
      if (.ig_g5_is_a_priori(p, "f2") && p$f2 == 0) {
        out <- c(out, "in a priori mode f2 must be greater than 0 (otherwise the power equals alpha for every N)")
      }
      out
    }
  )
}

# =====================================================================================
# Exact distributions: sample correlation r and squared multiple correlation R^2
# =====================================================================================
#
# Beta mixtures. Both distributions are mixtures of regularized incomplete
# beta functions with negative-binomial weights; the series are summed on the
# log scale over a window of terms that is widened until the terms at both
# edges are below exp(-46) of the largest term. The incomplete betas of one
# series differ by integer steps in their first shape, so one pbeta() call
# anchors the series and the others follow from the exact increments
#   I_y(s, b) - I_y(s + 1, b) = y^s (1 - y)^b / (s B(s, b)),
# accumulated in the direction in which they add up (no cancellation):
# I_y(s, b) (decreasing in s) from the largest shape down, its complement
# I_{1-y}(b, s) (increasing in s) from the smallest shape up.

.ig_g5_max_terms <- 5e5

# log sum_j w_j I_j for j >= 0, log w_j = logw(j); I_j = I_y(j + s0, b)
# ("lower") or I_{1-y}(b, j + s0) ("upper"). Returns NA when more than
# .ig_g5_max_terms terms would be needed.
.ig_g5_beta_mixture_log <- function(logw, s0, b, y, yc, tail, center, spread) {
  if (y <= 0 || yc <= 0) {
    zero_lower <- y <= 0
    if ((tail == "lower") == zero_lower) return(-Inf)
    # I is 1 for every term: the sum of the weights.
    half <- ceiling(14 * spread) + 40
    if (2 * half > .ig_g5_max_terms) return(NA_real_)
    j <- seq.int(max(0, floor(center) - half), floor(center) + half)
    return(.ig_logspace_sum(logw(j)))
  }
  ly <- log(y)
  lyc <- log(yc)
  half <- ceiling(14 * spread) + 40
  if (2 * half > .ig_g5_max_terms) return(NA_real_)
  lo <- max(0, floor(center) - half)
  hi <- floor(center) + half
  repeat {
    j <- seq.int(lo, hi)
    s <- j + s0
    k <- length(s)
    sd <- s[-k]
    ld <- if (k > 1L) sd * ly + b * lyc - log(sd) - lbeta(sd, b) else numeric(0)
    if (tail == "lower") {
      lv <- c(ld, stats::pbeta(y, s[k], b, log.p = TRUE))
      M <- max(lv)
      li <- log(rev(cumsum(rev(exp(lv - M))))) + M
    } else {
      lv <- c(stats::pbeta(yc, b, s[1L], log.p = TRUE), ld)
      M <- max(lv)
      li <- log(cumsum(exp(lv - M))) + M
    }
    lt <- logw(j) + li
    fin <- is.finite(lt)
    if (!any(fin)) return(-Inf)
    mx <- max(lt[fin])
    grow_left <- lo > 0 && is.finite(lt[1L]) && lt[1L] - mx > -46
    grow_right <- is.finite(lt[k]) && lt[k] - mx > -46
    if (!grow_left && !grow_right) return(.ig_logspace_sum(lt))
    if (hi - lo > .ig_g5_max_terms) return(NA_real_)
    step <- max(half, hi - lo)
    if (grow_left) lo <- max(0, lo - step)
    if (grow_right) hi <- hi + step
  }
}

# Sample correlation r of n bivariate normal pairs with population
# correlation rho. For x >= 0 (y = x^2, b = (n - 2)/2, a = (n - 1)/2,
# c = rho / sqrt(1 - rho^2)):
#
#   P(r <= x) = T_{n-1}(-c sqrt(n - 1)) + 1/2 sum_j [P_j I_y(j + 1/2, b) + Q_j I_y(j + 1, b)]
#   P(r >  x) = 1/2 sum_j [P_j I_{1-y}(b, j + 1/2) + Q_j I_{1-y}(b, j + 1)]
#
#   P_j = Gamma(a + j) / (Gamma(a) j!) rho^(2j) (1 - rho^2)^a           (negative binomial)
#   Q_j = Gamma(a + j + 1/2) / (Gamma(a) Gamma(j + 3/2)) rho^(2j+1) (1 - rho^2)^a
#
# This is the series of the noncentral t distribution (Lenth 1989, AS 243) of
# t = r sqrt(n - 2) / sqrt(1 - r^2) given the X sum of squares W ~ chi^2(n - 1)
# (noncentrality c sqrt(W)), averaged over W term by term; P(r <= 0) is the
# central t probability T_{n-1}(-c sqrt(n - 1)). Negative x uses
# P(r <= x | rho) = P(r >= -x | -rho). For rho >= 0 every term is positive;
# for rho < 0 the Q_j are negative and the upper tail at x >= 0 is a
# difference, whose relative accuracy is returned in attr "rel_err".
.ig_g5_r_tail <- function(x, rho, n, lower_tail) {
  if (is.na(x) || is.na(rho) || is.na(n)) return(NA_real_)
  if (x <= -1) return(if (lower_tail) 0 else 1)
  if (x >= 1) return(if (lower_tail) 1 else 0)
  if (abs(rho) >= 1) {
    below <- as.double(sign(rho) <= x)
    return(if (lower_tail) below else 1 - below)
  }
  if (x < 0) return(.ig_g5_r_tail(-x, -rho, n, !lower_tail))
  a <- (n - 1) / 2
  b <- (n - 2) / 2
  y <- x * x
  yc <- (1 - x) * (1 + x)
  if (rho == 0) {
    return(if (lower_tail) 0.5 + 0.5 * stats::pbeta(y, 0.5, b) else 0.5 * stats::pbeta(yc, b, 0.5))
  }
  r2 <- rho * rho
  lr <- log(abs(rho))
  l1r <- log1p(-r2)
  logP <- function(j) stats::dnbinom(j, size = a, prob = 1 - r2, log = TRUE)
  logQ <- function(j) (2 * j + 1) * lr + a * l1r - log(j + 0.5) - lbeta(a, j + 0.5)
  center <- a * r2 / (1 - r2)
  spread <- sqrt(a * r2) / (1 - r2) + 1
  series <- function(tail) {
    e <- .ig_g5_beta_mixture_log(logP, 0.5, b, y, yc, tail, center, spread)
    o <- .ig_g5_beta_mixture_log(logQ, 1, b, y, yc, tail, center, spread)
    c(e, o)
  }
  T0 <- stats::pt(-rho / sqrt(1 - r2) * sqrt(n - 1), n - 1)
  if (rho > 0) {
    # Both series are positive; the smaller tail is summed directly and the
    # larger one is its complement (a sum near 1 would carry rounding above 1).
    up <- series("upper")
    if (anyNA(up)) return(NA_real_)
    U <- 0.5 * (exp(up[1]) + exp(up[2]))
    if (U <= 0.5) return(if (lower_tail) 1 - U else U)
    lo <- series("lower")
    if (anyNA(lo)) return(NA_real_)
    F <- min(T0 + 0.5 * (exp(lo[1]) + exp(lo[2])), 0.5)
    return(if (lower_tail) F else 1 - F)
  }
  # rho < 0, x >= 0: P(r <= x) >= P(r <= 0) >= 1/2.
  if (lower_tail) {
    lo <- series("lower")
    if (anyNA(lo)) return(NA_real_)
    return(min(1, T0 + 0.5 * (exp(lo[1]) - exp(lo[2]))))
  }
  up <- series("upper")
  if (anyNA(up)) return(NA_real_)
  scale <- 0.5 * (exp(up[1]) + exp(up[2]))
  U <- 0.5 * (exp(up[1]) - exp(up[2]))
  structure(max(U, 0), rel_err = if (U > 0) 1e-13 * scale / U else Inf)
}

# Squared multiple correlation R^2 of m random predictors, N observations,
# population rho2 (Lee 1971; Benton & Krishnamoorthy 2003):
#   P(R^2 <= x) = sum_k NB(k; (N - 1)/2, rho2) I_x(m/2 + k, (N - m - 1)/2),
# every term positive; the upper tail uses the complements I_{1-x}((N-m-1)/2, m/2 + k).
.ig_g5_r2_tail <- function(x, rho2, m, N, lower_tail) {
  if (is.na(x) || is.na(rho2)) return(NA_real_)
  if (x <= 0) return(if (lower_tail) 0 else 1)
  if (x >= 1) return(if (lower_tail) 1 else 0)
  b <- (N - m - 1) / 2
  if (rho2 == 0) return(stats::pbeta(x, m / 2, b, lower.tail = lower_tail))
  if (rho2 >= 1) return(if (lower_tail) 0 else 1)
  a <- (N - 1) / 2
  logw <- function(k) stats::dnbinom(k, size = a, prob = 1 - rho2, log = TRUE)
  center <- a * rho2 / (1 - rho2)
  spread <- sqrt(a * rho2) / (1 - rho2) + 1
  l <- .ig_g5_beta_mixture_log(logw, m / 2, b, x, 1 - x, if (lower_tail) "lower" else "upper", center, spread)
  if (is.na(l)) NA_real_ else exp(l)
}

# Three-moment noncentral F approximation of the R^2 distribution (Lee 1971,
# p. 123; G*Power option "3-Moment F-Approximation"). Given the predictors,
# R^2 / (1 - R^2) = chi^2_m(phi W) / chi^2_(N-m-1) with phi = rho2 / (1 - rho2)
# and W ~ chi^2_(N-1). The numerator's first three cumulants
#   kappa_r = 2^(r-1) (r-1)! (m + n (gamma^r - 1)),  n = N - 1, gamma = 1 / (1 - rho2),
# are matched by g chi^2_nu(lambda):
#   g = (A2 - sqrt(A2^2 - A1 A3)) / A1,  lambda = (A2 - g A1) / g^2,  nu = A1 / g - lambda,
# with A_r = m + n (gamma^r - 1). Then
#   P(R^2 <= x) ~ P(F(nu, N - m - 1, lambda) <= x / (1 - x) (N - m - 1) / (g nu)).
# For rho2 = 0 it is the exact central beta distribution (g = 1, lambda = 0, nu = m).
.ig_g5_r2_lee_parameters <- function(rho2, m, N) {
  n <- N - 1
  gam <- 1 / (1 - rho2)
  A <- m + n * (gam^(1:3) - 1)
  disc <- max(A[2]^2 - A[1] * A[3], 0)
  g <- (A[2] - sqrt(disc)) / A[1]
  lambda <- max((A[2] - g * A[1]) / g^2, 0)
  nu <- A[1] / g - lambda
  list(g = g, lambda = lambda, nu = nu, df2 = N - m - 1)
}

.ig_g5_r2_tail_lee <- function(x, rho2, m, N, lower_tail) {
  if (x <= 0) return(if (lower_tail) 0 else 1)
  if (x >= 1) return(if (lower_tail) 1 else 0)
  if (rho2 == 0) return(stats::pbeta(x, m / 2, (N - m - 1) / 2, lower.tail = lower_tail))
  if (rho2 >= 1) return(if (lower_tail) 0 else 1)
  lp <- .ig_g5_r2_lee_parameters(rho2, m, N)
  q <- x / (1 - x) * lp$df2 / (lp$g * lp$nu)
  .ig_g5_p_f(q, lp$nu, lp$df2, lp$lambda, lower_tail = lower_tail)
}

# Quantiles by root finding on the log tail probability in a transformed
# scale (Fisher z for r, logit for R^2), with a small memo (critical values
# are reused across the evaluations of one solve).
.ig_g5_quantile_cache <- new.env(parent = emptyenv())

.ig_g5_cached <- function(key, compute) {
  hit <- .ig_g5_quantile_cache[[key]]
  if (!is.null(hit)) return(hit)
  if (length(ls(.ig_g5_quantile_cache, all.names = TRUE)) > 20000L) {
    rm(list = ls(.ig_g5_quantile_cache, all.names = TRUE), envir = .ig_g5_quantile_cache)
  }
  val <- compute()
  assign(key, val, envir = .ig_g5_quantile_cache)
  val
}

# Solve log(tail(g(u))) = log(p) for u, where tail is monotone in u
# (increasing for lower_tail, decreasing otherwise); u0 a start, step a scale.
.ig_g5_invert <- function(logtail, lp, lower_tail, u0, step, lower = -Inf, upper = Inf) {
  h <- function(u) {
    v <- logtail(u) - lp
    if (lower_tail) v else -v
  }
  a <- u0
  fa <- h(a)
  if (is.na(fa)) return(NA_real_)
  b <- a
  fb <- fa
  dir <- if (fa < 0) 1 else -1
  for (i in seq_len(200L)) {
    b <- min(max(a + dir * step, lower), upper)
    fb <- h(b)
    if (is.na(fb)) return(NA_real_)
    if (sign(fb) != sign(fa) || fb == 0) break
    if (b == lower || b == upper) return(NA_real_)
    a <- b
    fa <- fb
    step <- step * 2
  }
  if (fa == 0) return(a)
  if (fb == 0) return(b)
  if (sign(fb) == sign(fa)) return(NA_real_)
  stats::uniroot(h, sort(c(a, b)), f.lower = if (a < b) fa else fb, f.upper = if (a < b) fb else fa,
                 tol = 1e-14 * max(1, abs(a), abs(b)), maxiter = 1000L)$root
}

# Quantile of r: x with P(r <= x) = p (lower_tail) or P(r > x) = p.
.ig_g5_r_quantile <- function(p, rho, n, lower_tail) {
  if (rho == 0) {
    t <- stats::qt(p, n - 2, lower.tail = lower_tail)
    return(t / sqrt(n - 2 + t^2))
  }
  key <- sprintf("r|%.17g|%.17g|%.17g|%d", p, rho, n, lower_tail)
  .ig_g5_cached(key, function() {
    lp <- log(p)
    logtail <- function(z) {
      v <- .ig_g5_r_tail(tanh(z), rho, n, lower_tail)
      if (is.na(v)) NA_real_ else log(v)
    }
    sc <- 1 / sqrt(max(n - 3, 0.5))
    z0 <- atanh(rho) + stats::qnorm(p, lower.tail = lower_tail) * sc
    z <- .ig_g5_invert(logtail, lp, lower_tail, z0, 0.25 * sc, lower = -30, upper = 30)
    if (is.na(z)) return(NA_real_)
    x <- tanh(z)
    v <- .ig_g5_r_tail(x, rho, n, lower_tail)
    re <- attr(v, "rel_err")
    if (!is.null(re) && re > 1e-9) return(NA_real_)
    x
  })
}

# Quantile of R^2 (exact or three-moment approximation).
.ig_g5_r2_quantile <- function(p, rho2, m, N, lower_tail, approx = FALSE) {
  if (rho2 == 0) return(stats::qbeta(p, m / 2, (N - m - 1) / 2, lower.tail = lower_tail))
  key <- sprintf("r2|%.17g|%.17g|%.17g|%.17g|%d|%d", p, rho2, m, N, lower_tail, approx)
  .ig_g5_cached(key, function() {
    tail_fun <- if (approx) .ig_g5_r2_tail_lee else .ig_g5_r2_tail
    lp <- log(p)
    logtail <- function(u) {
      v <- tail_fun(stats::plogis(u), rho2, m, N, lower_tail)
      if (is.na(v)) NA_real_ else log(v)
    }
    # Start from the mean of R^2 under rho2 (about rho2 + (1 - rho2) m / (N - 1)).
    mu <- min(max(rho2 + (1 - rho2) * m / (N - 1), 1e-6), 1 - 1e-6)
    u <- .ig_g5_invert(logtail, lp, lower_tail, stats::qlogis(mu), 0.25, lower = -700, upper = 37)
    if (is.na(u)) NA_real_ else stats::plogis(u)
  })
}

# A priori search of the two exact procedures. Their forward models invert an
# exact distribution at every N, so the generic search's full enumeration of
# all totals below 4000 is slow. This search uses the generic policy for
# large solutions at every N: doubling and bisection, then the predecessor and
# a window of 50 totals below the solution must fail; any design in the window
# that reaches the target triggers a full enumeration from the smallest total
# (so non-monotone power, for example at the exact/Fisher-z switch, is not
# missed inside the window). Power is assumed non-decreasing in N outside it
# (a continuous test whose statistic concentrates at the H1 value).
.ig_g5_exact_a_priori_search <- function(forward, min_n) {
  function(p, n_lo, n_max) {
    first <- max(n_lo, min_n(p))
    target <- p$power
    cache <- new.env(parent = emptyenv())
    evals <- 0L
    P <- function(N) {
      key <- format(N, scientific = FALSE)
      v <- cache[[key]]
      if (!is.null(v)) return(v)
      q <- p
      q$n_total <- N
      v <- forward(q)$power
      evals <<- evals + 1L
      if (!is.numeric(v) || length(v) != 1L || is.na(v)) {
        .ig_abort("numerical_failure", sprintf("the power could not be computed for N = %s", key))
      }
      assign(key, v, envir = cache)
      v
    }
    diag <- function(...) c(list(...), list(evaluations = evals, n_lower = first, window = 50L))
    if (first > n_max) {
      return(list(status = "infeasible", reason = "no admissible total sample size in the search range",
                  diagnostics = diag()))
    }
    if (P(first) >= target) {
      return(list(status = "ok", n = first, power = P(first),
                  diagnostics = diag(search = "smallest design", minimality = "proven: the smallest admissible design reaches the target")))
    }
    L <- first
    H <- max(first + 1, 2 * first)
    repeat {
      H <- min(H, n_max)
      if (P(H) >= target) break
      if (H >= n_max) {
        return(list(status = "infeasible", reason = sprintf(
          "the target power %s is not reached by any N <= %s (power there: %s)", format(target),
          format(n_max, scientific = FALSE), format(signif(P(H), 7))), diagnostics = diag()))
      }
      L <- H
      H <- 2 * H
    }
    while (H - L > 1) {
      M <- floor((L + H) / 2)
      if (P(M) >= target) H <- M else L <- M
    }
    nonmono <- FALSE
    W <- L - 1
    while (W >= max(first, L - 50)) {
      if (P(W) >= target) {
        nonmono <- TRUE
        break
      }
      W <- W - 1
    }
    if (nonmono) {
      N <- first
      while (P(N) < target) N <- N + 1
      return(list(status = "ok", n = N, power = P(N), diagnostics = diag(
        search = "enumeration after non-monotone power was found in the window",
        minimality = "proven: every total from the smallest one up to the solution was evaluated",
        nonmonotone_detected = TRUE, bracket = c(L, H))))
    }
    list(status = "ok", n = H, power = P(H), diagnostics = diag(
      search = "doubling and bisection over N",
      minimality = "predecessor fails; the 50 totals below it also fail (power assumed non-decreasing in N below the window)",
      predecessor_n_total = L, predecessor_power = P(L), bracket = c(L, H), nonmonotone_detected = FALSE))
  }
}

# Two-sided region probabilities that tolerate an NA tail (a series beyond
# the term cap): the result is then NA instead of an error.
.ig_g5_two_sided <- function(ptail, lo, hi) {
  lower <- as.double(ptail(lo, TRUE))
  upper <- as.double(ptail(hi, FALSE))
  below_hi <- as.double(ptail(hi, TRUE))
  if (anyNA(c(lower, upper, below_hi))) return(list(reject = NA_real_, accept = NA_real_))
  accept <- if (below_hi <= 0.5) below_hi - lower else as.double(ptail(lo, FALSE)) - upper
  list(reject = lower + upper, accept = max(accept, 0))
}

# =====================================================================================
# exact_correlation_bivariate_normal
# =====================================================================================
#
# Manual ch. 3 (PDF pp. 9-10). H0 rho = rho0 for a bivariate normal
# population; the exact distribution of r under rho0 (critical values) and
# rho1 (power), or the Fisher z approximation
# z = (atanh(r) - atanh(rho0)) sqrt(N - 3) ~ N((atanh(rho1) - atanh(rho0)) sqrt(N - 3), 1).
# Options dialog: "use exact distribution if N <" x (default 10000) or "use
# large sample approximation (Fisher Z)". The captures show that the exact
# branch is used for N <= x and the approximation for N > x
# (SUP-G5-EX-CORRBN-POS-P01..P04, APR-A01), and that with the Fisher z option
# N = 3, where the approximation does not exist (variance 1 / (N - 3)), is
# computed with the exact distribution (SUP-G5-EX-CORRBN-POS-P05/P07, SEN-S01,
# CRI-C01, COM-M01). inferguidR uses the exact distribution at N = 3 under
# either option.
#
# A priori: inferguidR reports the smallest N whose reported power (the
# branch of that N) reaches the target; with a threshold G*Power does the
# same (SUP-G5-EX-CORRBN-APR-A01, A08-A10). With the Fisher z option
# G*Power 3.1.9.7 judges N = 4 by the exact distribution during its search
# but displays the Fisher z power of N = 4, so it can report N = 4 with a
# displayed power below the target (GPW-EX-CORRBN-APR-B004-E01,
# SUP-G5-EX-CORRBN-APR-A02/A05; controls A03, A04, A06, A07; method note).
# That design is not copied.

.ig_method_exact_correlation_bivariate_normal <- function() {
  ig_method(
    id = "exact_correlation_bivariate_normal", family = "exact",
    label = list(en = "Correlation: Bivariate normal model", de = "Korrelation: Bivariates Normalverteilungsmodell"),
    params = ig_params(
      ig_param("rho_h1", role = "effect", domain = ig_domain_real(-1, 1, FALSE, FALSE),
               label = list(en = "Correlation \u03c1 H1"), description = "Population correlation under H1.",
               modes = c("a_priori", "post_hoc", "criterion", "compromise")),
      .ig_g5_z_distance(paste(
        "Sensitivity solve target: |atanh(rho1) - atanh(rho0)|; the effect direction option gives the side.",
        "The solved Correlation rho H1 is method_outputs$rho_h1.")),
      ig_param("rho_h0", role = "nuisance", domain = ig_domain_real(-1, 1, FALSE, FALSE), default = 0,
               label = list(en = "Correlation \u03c1 H0"),
               description = "Population correlation under H0. Default 0 (G*Power's live field default)."),
      ig_param_alpha(), ig_param_power(), ig_param_error_ratio(),
      .ig_g5_n_total(3, "Total sample size N (number of pairs)."),
      ig_param("exact_n_threshold", role = "control",
               domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE), default = 10000,
               label = list(en = "use exact distribution if N <"),
               description = paste("Options-dialog field: the exact distribution is used for N <= this value and the",
                                   "Fisher z approximation above it (live default 10000; the captures show N = x uses",
                                   "the exact branch). N = 3 always uses the exact distribution."),
               applies_when = list(option = "computation", "in" = "exact"))
    ),
    options = ig_options(
      .ig_g5_option_tails(),
      ig_option("computation", values = c("exact", "fisher_z_approx"), default = "exact",
                labels = c(exact = "use exact distribution if N < x", fisher_z_approx = "use large sample approximation (Fisher Z)"),
                label = list(en = "Computation")),
      ig_option("effect_direction", values = c("r_le_rho0", "r_ge_rho0"), default = "r_ge_rho0",
                labels = c(r_le_rho0 = "r \u2264 \u03c1", r_ge_rho0 = "r \u2265 \u03c1"),
                applies_when = list(mode = "sensitivity"), label = list(en = "Effect direction"))
    ),
    modes = .ig_g5_modes,
    targets = list(sensitivity = "z_distance"),
    power = .ig_power_exact_correlation_bivariate_normal,
    design = NULL, discrete = FALSE, effect_direction = "increasing",
    a_priori_search = .ig_g5_exact_a_priori_search(.ig_power_exact_correlation_bivariate_normal, function(p) 3),
    engine = list(package = "inferguidR",
                  fun = "exact CDF of r as a negative-binomial mixture of pbeta (.ig_g5_r_tail); pnorm for Fisher z",
                  version = NA_character_),
    checker = "validation/independent/g5/check_g5_exact.R",
    reference = .ig_g5_reference("Tests > Correlation and regression > Correlation: Bivariate normal model",
                                 "Exact - Correlation: Bivariate normal model",
                                 "ch. 3 (PDF pp. 9-10)", "specs/method-notes/exact_correlation_bivariate_normal.md"),
    status = "independently_checked",
    constraints = function(p) {
      c(.ig_g5_power_alpha_problem(p),
        if (.ig_g5_is_a_priori(p, "rho_h1") && p$rho_h1 == p$rho_h0)
          "in a priori mode rho_h1 must differ from rho_h0 (G*Power: 'In a priori calculations the effect size must not be zero')")
    }
  )
}

.ig_power_exact_correlation_bivariate_normal <- function(p) {
  N <- as.double(p$n_total)
  rho0 <- p$rho_h0
  sens <- is.null(p$rho_h1)
  rho1 <- if (!sens) p$rho_h1 else {
    tanh(atanh(rho0) + .ig_g5_effect_direction_sign(p, "r_le_rho0") * p$z_distance)
  }
  tails <- p$tails
  side <- if (identical(tails, "two")) "both" else if (sens) {
    if (identical(p$effect_direction, "r_le_rho0")) "lower" else "upper"
  } else if (rho1 >= rho0) "upper" else "lower"
  fisher_selected <- identical(p$computation, "fisher_z_approx") || N > p$exact_n_threshold
  # The Fisher z statistic has variance 1 / (N - 3), so the approximation does
  # not exist for N = 3. With the Fisher z option G*Power 3.1.9.7 then
  # computes with the exact distribution in every analysis type
  # (SUP-G5-EX-CORRBN-POS-P05/P07, SEN-S01, CRI-C01, COM-M01); so does
  # inferguidR, visibly (outputs$branch and the assumptions). The same holds
  # here when the approximation is selected by a threshold below 3, where
  # G*Power displays the degenerate value power = alpha (SUP-G5-EX-CORRBN-POS-P08).
  fisher <- fisher_selected && N > 3
  alpha <- p$alpha
  outputs <- list(rho_h1 = rho1, branch = if (fisher) "fisher_z" else "exact")
  notes <- list()
  if (fisher) {
    s <- sqrt(N - 3)
    mu1 <- (atanh(rho1) - atanh(rho0)) * s
    if (!is.finite(mu1)) mu1 <- sign(rho1 - rho0) * Inf
    zc <- switch(side,
                 both = c(-1, 1) * stats::qnorm(alpha / 2, lower.tail = FALSE),
                 upper = stats::qnorm(alpha, lower.tail = FALSE),
                 lower = stats::qnorm(alpha))
    if (side == "both") {
      reg <- .ig_two_sided_prob(function(q, lt) stats::pnorm(q, mu1, lower.tail = lt), zc[1], zc[2])
      power <- reg$reject
      beta <- reg$accept
    } else if (side == "upper") {
      power <- stats::pnorm(zc, mu1, lower.tail = FALSE)
      beta <- stats::pnorm(zc, mu1)
    } else {
      power <- stats::pnorm(zc, mu1)
      beta <- stats::pnorm(zc, mu1, lower.tail = FALSE)
    }
    outputs$critical_z <- zc
    outputs$critical_r <- tanh(atanh(rho0) + zc / s)
    critical <- zc
    sd_hint <- list(h0 = list(family = "normal", mean = 0, sd = 1), h1 = list(family = "normal", mean = mu1, sd = 1))
    assumption <- sprintf(paste("Fisher z approximation: z = (atanh(r) - atanh(rho0)) sqrt(N - 3) ~ N(mu1, 1),",
                                "mu1 = (atanh(rho1) - atanh(rho0)) sqrt(N - 3)%s"),
                          if (identical(p$computation, "fisher_z_approx")) "" else
                            sprintf(" (N = %s exceeds the exact-distribution threshold %s)",
                                    format(N, scientific = FALSE), format(p$exact_n_threshold)))
  } else {
    q0 <- function(prob, lt) .ig_g5_r_quantile(prob, rho0, N, lt)
    rc <- switch(side, both = c(q0(alpha / 2, TRUE), q0(alpha / 2, FALSE)),
                 upper = q0(alpha, FALSE), lower = q0(alpha, TRUE))
    if (anyNA(rc)) {
      return(list(power = NA_real_, alpha_attained = alpha,
                  notes = list(power = paste("the critical value of the exact r distribution could not be computed to the",
                                             "required relative accuracy (extreme tail or |rho| too close to 1)")),
                  outputs = outputs))
    }
    ptail <- function(q, lt) .ig_g5_r_tail(q, rho1, N, lt)
    if (side == "both") {
      reg <- .ig_g5_two_sided(ptail, rc[1], rc[2])
      power <- reg$reject
      beta <- reg$accept
    } else if (side == "upper") {
      power <- ptail(rc, FALSE)
      beta <- ptail(rc, TRUE)
    } else {
      power <- ptail(rc, TRUE)
      beta <- ptail(rc, FALSE)
    }
    power <- as.double(power)
    beta <- as.double(beta)
    if (is.na(power) || is.na(beta)) {
      return(list(power = NA_real_, alpha_attained = alpha,
                  notes = list(power = paste("the exact distribution of r under H1 needs more than",
                                             format(.ig_g5_max_terms), "series terms (|rho1| too close to 1)")),
                  outputs = outputs))
    }
    outputs$critical_r <- rc
    critical <- rc
    sd_hint <- list(h0 = list(family = "sample_correlation", rho = rho0, n = N),
                    h1 = list(family = "sample_correlation", rho = rho1, n = N))
    assumption <- "exact distribution of the sample correlation of N bivariate normal pairs under rho0 and rho1"
    if (fisher_selected) {
      assumption <- c(assumption, paste(
        "N = 3: the Fisher z approximation does not exist (its variance is 1 / (N - 3)); the exact distribution",
        "is used instead",
        if (identical(p$computation, "fisher_z_approx")) "(as in G*Power 3.1.9.7)" else
          "(G*Power 3.1.9.7 displays the degenerate value power = alpha for a threshold below 3)"))
    }
  }
  sd_hint$rejection_side <- side
  list(power = power, beta = beta, alpha_attained = alpha, critical = critical, n_total = N,
       statistic_distribution = sd_hint, outputs = outputs, notes = notes,
       assumptions = c("bivariate normal population; H0 rho = rho0", assumption,
                       if (side == "both") "two-sided: alpha/2 in each tail of the H0 distribution" else
                         sprintf("one-sided: alpha in the %s tail of the H0 distribution (side of rho1 - rho0)", side)))
}

# =====================================================================================
# exact_mrc_random_model
# =====================================================================================
#
# Manual ch. 7 (PDF pp. 18-21). Random (multivariate normal) predictors; H0
# rho^2 = rho0^2 tested with the exact distribution of R^2 (Lee 1971; Benton &
# Krishnamoorthy 2003), or Lee's three-moment noncentral F approximation
# (option), which the captures show is used for both the H0 critical values
# and the power (SUP-G5-EX-MRCRND-POS-P01 vs P02).

.ig_method_exact_mrc_random_model <- function() {
  ig_method(
    id = "exact_mrc_random_model", family = "exact",
    label = list(en = "Linear multiple regression: Random model", de = "Multiple lineare Regression: Zufallsmodell"),
    params = ig_params(
      ig_param("rho2_h1", role = "effect", domain = ig_domain_real(0, 1, upper_inclusive = FALSE),
               label = list(en = "H1 \u03c1\u00b2"), description = "Population squared multiple correlation under H1.",
               modes = c("a_priori", "post_hoc", "criterion", "compromise")),
      .ig_g5_z_distance(paste(
        "Sensitivity solve target: |atanh(sqrt(rho2_h1)) - atanh(sqrt(rho2_h0))|; the effect direction option",
        "gives the side (on the lower side the distance stops at H1 rho^2 = 0). The solved H1 rho^2 is",
        "method_outputs$rho2_h1.")),
      ig_param("rho2_h0", role = "nuisance", domain = ig_domain_real(0, 1, upper_inclusive = FALSE), default = 0,
               label = list(en = "H0 \u03c1\u00b2"),
               description = "Population squared multiple correlation under H0. Default 0 (G*Power's live field default)."),
      ig_param_alpha(), ig_param_power(), ig_param_error_ratio(),
      .ig_g5_n_total(3, "Total sample size N (N >= m + 2)."),
      .ig_g5_count("n_predictors", "Number of predictors", "Number m of random predictors.")
    ),
    options = ig_options(
      .ig_g5_option_tails(),
      ig_option("distribution", values = c("exact", "three_moment_f"), default = "exact",
                labels = c(exact = "Exact distribution", three_moment_f = "3-Moment F-Approximation"),
                label = list(en = "Distribution")),
      ig_option("effect_direction", values = c("rho2_h1_le_h0", "rho2_h1_ge_h0"), default = "rho2_h1_ge_h0",
                labels = c(rho2_h1_le_h0 = "\u03c1\u00b2 H1\u2264H0", rho2_h1_ge_h0 = "\u03c1\u00b2 H1\u2265H0"),
                applies_when = list(mode = "sensitivity"), label = list(en = "Effect direction"))
    ),
    modes = .ig_g5_modes,
    targets = list(sensitivity = "z_distance"),
    power = .ig_power_exact_mrc_random_model,
    design = .ig_g5_single_design(function(p) p$n_predictors + 2),
    discrete = FALSE, effect_direction = "increasing",
    a_priori_search = .ig_g5_exact_a_priori_search(.ig_power_exact_mrc_random_model, function(p) p$n_predictors + 2),
    engine = list(package = "inferguidR",
                  fun = "exact CDF of R^2 as a negative-binomial mixture of pbeta (.ig_g5_r2_tail); Lee three-moment F (.ig_p_f)",
                  version = NA_character_),
    checker = "validation/independent/g5/check_g5_exact.R",
    reference = .ig_g5_reference("Tests > Correlation and regression > Linear multiple regression: Random model",
                                 "Exact - Linear multiple regression: Random model",
                                 "ch. 7 (PDF pp. 18-21)", "specs/method-notes/exact_mrc_random_model.md"),
    status = "independently_checked",
    constraints = function(p) {
      out <- .ig_g5_power_alpha_problem(p)
      if (!is.null(p$n_total) && p$n_total < p$n_predictors + 2) {
        out <- c(out, sprintf("N must be at least the number of predictors + 2 (%s)", format(p$n_predictors + 2)))
      }
      if (.ig_g5_is_a_priori(p, "rho2_h1") && p$rho2_h1 == p$rho2_h0) {
        out <- c(out, "in a priori mode rho2_h1 must differ from rho2_h0 (otherwise the power equals alpha for every N)")
      }
      out
    }
  )
}

.ig_power_exact_mrc_random_model <- function(p) {
  N <- as.double(p$n_total)
  m <- p$n_predictors
  r20 <- p$rho2_h0
  sens <- is.null(p$rho2_h1)
  notes <- list()
  if (!sens) {
    r21 <- p$rho2_h1
  } else {
    z0 <- atanh(sqrt(r20))
    if (identical(p$effect_direction, "rho2_h1_le_h0")) {
      z1 <- z0 - p$z_distance
      if (z1 < 0) {
        z1 <- 0
        notes$rho2_h1 <- "distance beyond H1 rho^2 = 0: the lower side ends at rho^2 = 0"
      }
    } else {
      z1 <- z0 + p$z_distance
    }
    r21 <- tanh(z1)^2
  }
  tails <- p$tails
  side <- if (identical(tails, "two")) "both" else if (sens) {
    if (identical(p$effect_direction, "rho2_h1_le_h0")) "lower" else "upper"
  } else if (r21 >= r20) "upper" else "lower"
  approx <- identical(p$distribution, "three_moment_f")
  alpha <- p$alpha
  q0 <- function(prob, lt) .ig_g5_r2_quantile(prob, r20, m, N, lt, approx)
  rc <- switch(side, both = c(q0(alpha / 2, TRUE), q0(alpha / 2, FALSE)),
               upper = q0(alpha, FALSE), lower = q0(alpha, TRUE))
  outputs <- list(rho2_h1 = r21, f2_h1 = r21 / (1 - r21), distribution = if (approx) "three_moment_f" else "exact")
  if (anyNA(rc)) {
    return(list(power = NA_real_, alpha_attained = alpha, outputs = outputs,
                notes = list(power = "the critical R^2 under H0 could not be computed")))
  }
  tail1 <- if (approx) .ig_g5_r2_tail_lee else .ig_g5_r2_tail
  ptail <- function(q, lt) tail1(q, r21, m, N, lt)
  if (side == "both") {
    reg <- .ig_g5_two_sided(ptail, rc[1], rc[2])
    power <- reg$reject
    beta <- reg$accept
  } else if (side == "upper") {
    power <- ptail(rc, FALSE)
    beta <- ptail(rc, TRUE)
  } else {
    power <- ptail(rc, TRUE)
    beta <- ptail(rc, FALSE)
  }
  if (is.na(power) || is.na(beta)) {
    return(list(power = NA_real_, alpha_attained = alpha, outputs = outputs,
                notes = list(power = paste("the exact R^2 distribution under H1 needs more than",
                                           format(.ig_g5_max_terms), "series terms (rho^2 too close to 1)"))))
  }
  outputs$critical_r2 <- rc
  list(power = power, beta = beta, alpha_attained = alpha, critical = rc, n_total = N,
       df = c(m, N - m - 1),
       statistic_distribution = list(h0 = list(family = "squared_multiple_correlation", rho2 = r20, m = m, n = N),
                                     h1 = list(family = "squared_multiple_correlation", rho2 = r21, m = m, n = N),
                                     rejection_side = side, approximation = if (approx) "lee_three_moment_f" else "none"),
       outputs = outputs, notes = notes,
       assumptions = c("random predictors: (Y, X1..Xm) multivariate normal (unconditional model, Gatsonis & Sampson 1989)",
                       if (approx) "Lee (1971) three-moment noncentral F approximation of the R^2 distribution under H0 and H1" else
                         "exact distribution of R^2 under rho0^2 and rho1^2",
                       if (side == "both") "two-sided: alpha/2 in each tail of the H0 distribution" else
                         sprintf("one-sided: alpha in the %s tail of the H0 distribution (side of rho1^2 - rho0^2)", side)))
}
