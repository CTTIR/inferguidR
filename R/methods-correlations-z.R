# Method definitions: correlation z tests (group G7) ------------------------------------
#
# Ownership: this file belongs to group G7 (blueprint phase 08, correlation
# procedures of the z family). It defines four procedures of the certified
# G*Power 3.1.9.7 inventory:
#
#   z_correlations_two_independent               Fisher-z difference of two independent r's
#   z_correlations_two_dependent_common_index    Steiger (1980) Z1* (rho_ab vs rho_ac)
#   z_correlations_two_dependent_no_common_index Steiger (1980) Z2* (rho_ab vs rho_cd)
#   z_tetrachoric_correlation                    Wald test of a tetrachoric correlation
#                                                (Brown & Benedetti 1977 / Bonett & Price 2005)
#
# Method notes (mathematical contract, G*Power conventions, evidence):
# specs/method-notes/<method_id>.md. Effect-size drawers:
# R/effect-correlations-z.R. Independent checks: validation/independent/g7/.
#
# Conventions shared by the four procedures (from the manual and the
# emulator captures; method notes give the record ids):
#
# * The test statistic Z is standard normal under H0 and N(mu, sd1) under H1.
#   It is oriented so that a positive effect (H1 value above the H0 value)
#   gives a positive mu; G*Power displays the critical z with the sign of the
#   effect, which is the critical value on the rejection side of a one-tailed
#   test (manual 25.5 for t; the same convention is shown by every capture
#   of these four tests).
# * One tail: alpha in the tail on the side of the effect; a zero effect uses
#   the lower tail, as G*Power displays it (critical z -1.6448536 for q = 0,
#   rho_ac = rho_ab, rho_cd = rho_ab and rho H1 = rho H0; captures
#   GPW-Z-CORR2-APR-B001-I02, SUP-RG7-Z-CORR2-POS-R13, SUP-RG7-Z-DEPCI-POS-R07,
#   SUP-RG7-Z-DEPNCI-POS-R07, SUP-RG7-Z-TETRA-POS-R17). Two tails: alpha/2 in
#   each tail and both rejection tails count in the power (G*Power convention,
#   manual 25.5).
# * Sensitivity mode of the three procedures with an "Effect direction"
#   option: the option fixes the rejection side of a one-tailed test and the
#   half-line in which the H1 correlation is searched (manual 28.3.3: "two
#   possible solutions, namely one for rho_ac <= rho_ab and one for rho_ac >=
#   rho_ab"). See .ig_g7_sensitivity_domain() for how this is passed to the
#   generic solver.

.ig_define_methods_z_g7 <- function() {
  list(
    .ig_method_z_correlations_two_independent(),
    .ig_method_z_correlations_two_dependent_common_index(),
    .ig_method_z_correlations_two_dependent_no_common_index(),
    .ig_method_z_tetrachoric_correlation()
  )
}

# Shared helpers ----------------------------------------------------------------------------

# G*Power requires the H0 and H1 correlations of the dependent-correlation
# and tetrachoric procedures in [-0.999999, 0.999999] (manual 28.1; live
# messages 'Effect size r_cd is outside the admissible range (-0.999999,
# 0.999999).', 'Corr rho_ab is outside the admissible range (-0.999999,
# 0.999999).', 'The H0 correlation must lie in the interval [-0.999999,
# 0.999999].': SUP-RG7-Z-DEPCI-POS-R15/R16, SUP-RG7-Z-DEPNCI-POS-R13,
# SUP-RG7-Z-TETRA-POS-R19/R23). The bound is the domain of those parameters
# and the end of the sensitivity search.
.ig_g7_rho_eps <- 1e-6

# Power of a z test whose statistic is N(0, 1) under H0 and N(mu, sd1) under
# H1. `side` (+1 or -1) is the rejection side of a one-tailed test. beta is
# computed directly as the non-rejection probability (no 1 - power), with
# the tail orientation that avoids cancellation.
.ig_g7_z_power <- function(mu, sd1, alpha, tails, side) {
  if (!is.finite(mu) || !is.finite(sd1) || sd1 < 0) {
    return(list(power = NA_real_, beta = NA_real_, critical = NA_real_, side = side))
  }
  if (identical(tails, "two")) {
    zc <- stats::qnorm(alpha / 2, lower.tail = FALSE)
    if (sd1 == 0) {
      reject <- abs(mu) > zc
      return(list(power = as.double(reject), beta = as.double(!reject), critical = c(-zc, zc), side = side))
    }
    a <- (-zc - mu) / sd1
    b <- (zc - mu) / sd1
    power <- stats::pnorm(b, lower.tail = FALSE) + stats::pnorm(a)
    beta <- if (a + b > 0) {
      stats::pnorm(a, lower.tail = FALSE) - stats::pnorm(b, lower.tail = FALSE)
    } else {
      stats::pnorm(b) - stats::pnorm(a)
    }
    return(list(power = min(power, 1), beta = max(beta, 0), critical = c(-zc, zc), side = side))
  }
  zc <- stats::qnorm(alpha, lower.tail = FALSE)
  if (side >= 0) {
    if (sd1 == 0) return(list(power = as.double(mu > zc), beta = as.double(mu <= zc), critical = zc, side = 1))
    u <- (zc - mu) / sd1
    list(power = stats::pnorm(u, lower.tail = FALSE), beta = stats::pnorm(u), critical = zc, side = 1)
  } else {
    if (sd1 == 0) return(list(power = as.double(mu < -zc), beta = as.double(mu >= -zc), critical = -zc, side = -1))
    u <- (-zc - mu) / sd1
    list(power = stats::pnorm(u), beta = stats::pnorm(u, lower.tail = FALSE), critical = -zc, side = -1)
  }
}

# Rejection side of a one-tailed test: the effect direction option in
# sensitivity mode (codes *_ge_* / *_le_*), otherwise the sign of the effect.
# No effect uses the lower tail: G*Power 3.1.9.7 displays the critical z
# -1.6448536 for a zero effect in all four procedures (the power is alpha on
# either side).
.ig_g7_side <- function(effect, p) {
  d <- p$effect_direction
  if (!is.null(d)) return(if (grepl("_le_", d, fixed = TRUE)) -1 else 1)
  if (!is.na(effect) && effect > 0) 1 else -1
}

# Sensitivity scans use the requested half-interval and, for dependent
# correlations, its intersection with positive semi-definite matrices.
.ig_g7_sensitivity_domain <- function(p, kind = NULL) {
  lim <- 1 - .ig_g7_rho_eps
  if (is.null(kind)) {
    null <- p$rho_h0
    iv <- c(lower = -lim, upper = lim)
  } else {
    sp <- .ig_g7_dep_spec(kind)
    null <- p$rho_ab
    iv <- .ig_g7_psd_interval(function(v) sp$build(p, v), null)
  }
  seg <- if (grepl("_le_", p$effect_direction, fixed = TRUE))
    c(iv[["lower"]], null) else c(null, iv[["upper"]])
  ig_domain_real(seg[1], seg[2])
}

.ig_g7_sensitivity_direction <- function(p) {
  if (grepl("_le_", p$effect_direction, fixed = TRUE)) "decreasing" else "increasing"
}

.ig_g7_tail_note <- function(tails, what) {
  if (identical(tails, "two")) {
    "two-sided test: alpha/2 in each tail; the power counts both rejection tails (G*Power convention, manual 25.5)"
  } else {
    sprintf(paste("one-sided test: the rejection region lies on the side of the effect (%s), or on the side",
                  "given by the effect direction option in sensitivity mode; no effect uses the lower tail"), what)
  }
}

# G*Power refuses a target power below alpha in a priori and sensitivity
# mode ("Power must be greater than or equal to alpha"); both are given only
# in those modes.
.ig_g7_power_alpha_constraint <- function(p) {
  if (!is.null(p$power) && !is.null(p$alpha) && is.numeric(p$power) && is.numeric(p$alpha) &&
      p$power < p$alpha) {
    sprintf("the target power %s is below alpha %s; it must be at least alpha",
            .ig_format_double(p$power), .ig_format_double(p$alpha))
  } else {
    character()
  }
}

# Smallest eigenvalue of a symmetric matrix (PSD check of correlation
# matrices). Values above -1e-12 count as positive semi-definite (rounding of
# exactly singular matrices such as rho_bc = 1 with rho_ac = rho_ab).
.ig_g7_psd_tol <- 1e-12

.ig_g7_min_eigen <- function(m) {
  min(eigen(m, symmetric = TRUE, only.values = TRUE)$values)
}

.ig_g7_criterion_alpha_note <- function(alpha, tails) {
  if (!identical(tails, "two") && is.finite(alpha) && alpha > 0.5) {
    paste("alpha above 0.5 in a one-tailed test: the critical z lies on the opposite side of zero from",
          "the effect (critical z = sign(effect) * qnorm(1 - alpha))")
  } else {
    NULL
  }
}

# ===========================================================================================
# Correlations: Two independent Pearson r's (manual ch. 27)
# ===========================================================================================
#
#   q = z1 - z2 (Fisher z, z = atanh(r)), q in [-10, 10]
#   H0: N(0, 1); H1: N(q / s, 1), s = sqrt(1/(n1 - 3) + 1/(n2 - 3))      (manual 27.5)
#
# A priori allocation (emulator captures; method notes, "Allocation"): the
# total N is the smallest N whose power with fractional group sizes
# N/(1 + r), rN/(1 + r) reaches the target (r = allocation ratio N2/N1); the
# integer split is n1 = floor(N/(1 + r) + 1/2), n2 = N - n1. With r = 1 the
# groups are equal (N even).

.ig_method_z_correlations_two_independent <- function() {
  ig_method(
    id = "z_correlations_two_independent",
    family = "z",
    label = list(en = "Correlations: Two independent Pearson r's",
                 de = "Korrelationen: Zwei unabh\u00e4ngige Pearson-r"),
    params = ig_params(
      ig_param("q", role = "effect", domain = ig_domain_real(-10, 10),
               label = list(en = "Effect size q", de = "Effektst\u00e4rke q"),
               description = paste("Cohen's q = atanh(r1) - atanh(r2), the difference of the Fisher-z transformed",
                                   "correlations of the two populations (manual 27.1); G*Power requires q in [-10, 10]."),
               aliases = "effect_size_q"),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      ig_param("n_total", role = "design", domain = ig_domain_integer(lower = 8), unit = "subjects",
               label = list(en = "Total sample size", de = "Gesamtstichprobenumfang"),
               description = "Solved in a priori mode (n1 + n2).",
               modes = "a_priori"),
      ig_param("n1", role = "design", domain = ig_domain_integer(lower = 3), unit = "subjects",
               label = list(en = "Sample size group 1", de = "Stichprobenumfang Gruppe 1"),
               description = paste("At least 3 (as in G*Power). With n = 3 the Fisher-z variance 1/(n - 3) is",
                                   "infinite and the power equals alpha."),
               aliases = "sample_size_group_1",
               modes = c("post_hoc", "sensitivity", "criterion", "compromise")),
      ig_param("n2", role = "design", domain = ig_domain_integer(lower = 3), unit = "subjects",
               label = list(en = "Sample size group 2", de = "Stichprobenumfang Gruppe 2"),
               description = paste("At least 3 (as in G*Power). With n = 3 the Fisher-z variance 1/(n - 3) is",
                                   "infinite and the power equals alpha."),
               aliases = "sample_size_group_2",
               modes = c("post_hoc", "sensitivity", "criterion", "compromise")),
      ig_param("allocation_ratio", role = "design",
               domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE),
               default = 1,
               label = list(en = "Allocation ratio N2/N1", de = "Allokationsverh\u00e4ltnis N2/N1"),
               description = paste("Group-size ratio n2/n1 of the a priori design; default 1 (the initial value of",
                                   "the G*Power field, live capture 2026-10-05)."),
               aliases = "allocation_ratio_n2_n1",
               modes = "a_priori")
    ),
    options = ig_options(
      ig_option("tails", values = c("one", "two"), default = "one",
                labels = c(one = "One", two = "Two"),
                label = list(en = "Tail(s)", de = "Seiten"))
    ),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = .ig_power_z_correlations_two_independent,
    design = .ig_design_z_correlations_two_independent,
    discrete = FALSE,
    effect_direction = "increasing",
    engine = list(package = "stats", fun = "pnorm/qnorm (closed form, manual 27.5)", version = NA_character_),
    checker = "validation/independent/g7/check_z_correlations_two_independent.R",
    reference = list(
      gpower_path = "Tests > Correlation and regression > Correlations: Two independent Pearson r's",
      gpower_label = "z tests - Correlations: Two independent Pearson r's",
      manual_pages = "ch. 27 (PDF p. 65); tail convention 25.5",
      allocation_rule = paste(
        "N = smallest total whose power with fractional group sizes N/(1+r), rN/(1+r) reaches the target;",
        "r = 1: n1 = n2 = N/2 (N even); r != 1: n1 = floor(N/(1+r) + 1/2), n2 = N - n1 (G*Power 3.1.9.7,",
        "emulator captures); where that split misses the target, the smallest larger N whose split reaches it",
        "(diagnostics$adjusted_from_gpower_rule)"),
      method_notes = "specs/method-notes/z_correlations_two_independent.md"
    ),
    status = "independently_checked",
    constraints = .ig_g7_power_alpha_constraint,
    a_priori_search = .ig_a_priori_z_correlations_two_independent
  )
}

# Fisher-z standard error of z1 - z2 and the H1 mean for (possibly
# fractional) group sizes.
.ig_z2_core <- function(q, n1, n2) {
  s <- sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
  list(s = s, mu = q / s)
}

.ig_power_z_correlations_two_independent <- function(p) {
  n1 <- as.double(p$n1)
  n2 <- as.double(p$n2)
  core <- .ig_z2_core(p$q, n1, n2)
  side <- .ig_g7_side(p$q, p)
  zp <- .ig_g7_z_power(core$mu, 1, p$alpha, p$tails, side)
  notes <- list()
  if (is.na(zp$power)) {
    notes$power <- notes$beta <- notes$critical <- notes$ncp <-
      "group sizes below 3 leave the Fisher-z variance 1/(n - 3) undefined"
  }
  assumptions <- character()
  if (min(n1, n2) == 3) {
    assumptions <- paste("a group of size 3 has infinite Fisher-z variance 1/(n - 3): the test carries no",
                         "information and its power equals alpha (limit n -> 3, as in G*Power)")
  }
  list(
    power = zp$power,
    beta = zp$beta,
    alpha_attained = p$alpha,
    critical = zp$critical,
    ncp = core$mu,
    n_total = n1 + n2,
    group_sizes = c(n1 = n1, n2 = n2),
    statistic_distribution = list(
      h0 = list(family = "normal", mean = 0, sd = 1),
      h1 = list(family = "normal", mean = core$mu, sd = 1),
      rejection_side = if (identical(p$tails, "two")) "both" else if (zp$side > 0) "upper" else "lower"
    ),
    outputs = list(se_z_difference = core$s),
    notes = notes,
    assumptions = c(
      "two independent samples from bivariate normal populations; Fisher z test of rho1 = rho2",
      "q = atanh(r1) - atanh(r2); H1: N(q / s, 1) with s = sqrt(1/(n1 - 3) + 1/(n2 - 3)) (manual 27.5)",
      .ig_g7_tail_note(p$tails, "the sign of q"),
      assumptions
    )
  )
}

# Allocation rule. Equal allocation (r = 1): N even, n1 = n2 = N/2. Other
# ratios: n1 = floor(N/(1 + r) + 1/2), n2 = N - n1 (G*Power 3.1.9.7). A design
# with a group below 4 is not admissible. `continuous = TRUE` gives the
# fractional split used by G*Power's search criterion.
.ig_design_z_correlations_two_independent <- function(p, n_total, continuous = FALSE) {
  r <- p$allocation_ratio
  if (continuous) return(c(n1 = n_total / (1 + r), n2 = r * n_total / (1 + r)))
  if (r == 1) {
    if (n_total %% 2 != 0 || n_total / 2 < 4) return(NULL)
    return(c(n1 = n_total / 2, n2 = n_total / 2))
  }
  n1 <- floor(n_total / (1 + r) + 0.5)
  n2 <- n_total - n1
  if (n1 < 4 || n2 < 4) return(NULL)
  c(n1 = n1, n2 = n2)
}

# A priori search (G*Power 3.1.9.7 rule, method notes "Allocation"): the
# smallest admissible N whose power with fractional group sizes reaches the
# target. The fractional power is strictly increasing in N for q != 0 (s
# decreases), so bracketing plus bisection finds N*, and the predecessor
# proves minimality. The integer design at N* is returned when its power
# reaches the target. Otherwise G*Power reports a design below the target
# (as for the t test, validation/compare/t_means_two_independent/
# allocation-rule.md); the engine then reports the smallest larger N whose
# split reaches the target, records G*Power's design in
# diagnostics$gpower_allocation_rule and flags diagnostics$adjusted_from_gpower_rule
# (the convention of group G1 for the t test).
.ig_a_priori_z_correlations_two_independent <- function(p, n_lo, n_max) {
  r <- p$allocation_ratio
  target <- p$power
  step <- if (r == 1) 2 else 1
  fw <- .ig_power_z_correlations_two_independent
  frac_power <- function(N) {
    gs <- .ig_design_z_correlations_two_independent(p, N, continuous = TRUE)
    if (any(gs <= 3)) return(NA_real_)
    fw(c(p, list(n1 = gs[["n1"]], n2 = gs[["n2"]])))$power
  }
  int_power <- function(N) {
    gs <- .ig_design_z_correlations_two_independent(p, N)
    if (is.null(gs)) return(NA_real_)
    fw(c(p, as.list(gs)))$power
  }
  admissible <- function(N) !is.null(.ig_design_z_correlations_two_independent(p, N))
  # Smallest admissible N (both groups at least 4).
  first <- max(n_lo, 8)
  while (first <= n_max && !admissible(first)) first <- first + 1
  diag <- list(criterion = "power with fractional group sizes N/(1+r), rN/(1+r) (G*Power 3.1.9.7 rule)",
               allocation_ratio = r, n_lower = first, n_upper_limit = n_max, step = step)
  if (first > n_max) {
    return(list(status = "infeasible", reason = "no admissible design with both groups of at least 4",
                diagnostics = diag))
  }
  if (p$q == 0) {
    return(list(status = "infeasible", reason = paste(
      "with q = 0 the power equals alpha for every N, so the target power is not reached by any design;",
      "power does not change with N"), diagnostics = c(diag, list(flat = TRUE))))
  }
  p_first <- frac_power(first)
  if (is.na(p_first)) {
    stop("the fractional power is undefined at the smallest admissible design", call. = FALSE)
  }
  if (p_first < target) {
    # Bracket by doubling, then bisect over the admissible lattice.
    lo <- first
    hi <- NA_real_
    cand <- first
    repeat {
      cand <- min(max(2 * cand, cand + step), n_max)
      if (step == 2 && cand %% 2 != first %% 2) cand <- cand - 1
      pc <- frac_power(cand)
      if (pc >= target) {
        hi <- cand
        break
      }
      if (cand >= n_max - step + 1) break
      lo <- cand
    }
    if (is.na(hi)) {
      p_lim <- frac_power(lo)
      flat <- abs(p_lim - p_first) <= 1e-14
      return(list(status = "infeasible", reason = sprintf(
        "the target power %s is not reached by any design with N <= %s (power at N = %s: %s)%s",
        format(target), format(n_max, scientific = FALSE), format(lo, scientific = FALSE),
        format(signif(p_lim, 7)), if (flat) "; power does not change with N" else ""),
        diagnostics = c(diag, list(power_at_limit = p_lim, flat = flat))))
    }
    while (hi - lo > step) {
      mid <- lo + step * floor((hi - lo) / (2 * step))
      if (mid <= lo) mid <- lo + step
      if (frac_power(mid) >= target) hi <- mid else lo <- mid
    }
    n_star <- hi
  } else {
    n_star <- first
  }
  gs <- .ig_design_z_correlations_two_independent(p, n_star)
  p_int <- int_power(n_star)
  diag$fractional_power <- frac_power(n_star)
  if (n_star > first) diag$predecessor <- list(n_total = n_star - step, fractional_power = frac_power(n_star - step))
  diag$minimality <- if (n_star > first) {
    "proven: the fractional power is strictly increasing in N and the predecessor fails"
  } else {
    "proven: the smallest admissible design reaches the target"
  }
  # Continuous root of the fractional power (diagnostic only).
  diag$continuous_relaxation <- tryCatch({
    f <- function(N) frac_power(N) - target
    lo_c <- max(n_star - step, first)
    if (!is.na(f(lo_c)) && f(lo_c) < 0) {
      list(n_total = stats::uniroot(f, c(lo_c, n_star), tol = 1e-10)$root,
           note = "root of the power curve with fractional group sizes; not an admissible design")
    } else NULL
  }, error = function(e) NULL)
  diag$gpower_allocation_rule <- list(n_total = n_star, n1 = gs[["n1"]], n2 = gs[["n2"]], power = p_int,
                                      reaches_target = isTRUE(p_int >= target))
  # Smallest admissible N whose integer design reaches the target (integer-
  # design criterion; transparency only, G*Power chooses N by the fractional
  # power). Both group sizes of the split are non-decreasing in N, so the
  # integer-design power is non-decreasing and a bisection finds it.
  smallest_int <- function(from) {
    if (from > n_max) return(NA_real_)
    v <- int_power(from)
    if (!is.na(v) && v >= target) return(from)
    lo <- from
    hi <- from
    repeat {
      hi <- min(max(2 * hi, hi + step), n_max)
      if (step == 2 && hi %% 2 != from %% 2) hi <- hi - 1
      v <- int_power(hi)
      if (!is.na(v) && v >= target) break
      if (hi >= n_max - step + 1) return(NA_real_)
      lo <- hi
    }
    while (hi - lo > step) {
      mid <- lo + step * floor((hi - lo) / (2 * step))
      if (mid <= lo) mid <- lo + step
      v <- int_power(mid)
      if (!is.na(v) && v >= target) hi <- mid else lo <- mid
    }
    hi
  }
  if (isTRUE(p_int >= target)) {
    n_int <- smallest_int(first)
    if (!is.na(n_int) && n_int != n_star) {
      diag$integer_criterion_n_total <- n_int
      diag$integer_criterion_note <- paste(
        "smallest N whose integer design reaches the target; G*Power's rule chooses N by the power with",
        "fractional group sizes instead")
    }
    return(list(status = "ok", n = n_star, power = p_int, diagnostics = diag))
  }
  # G*Power's design misses the target (it reports it without a warning, as
  # for the t test). As group G1 does for the t test (R/methods-t.R), the
  # engine reports the smallest larger total whose split reaches the target
  # and records G*Power's design in the diagnostics.
  n_adj <- smallest_int(n_star + step)
  if (is.na(n_adj)) {
    return(list(status = "infeasible", reason = sprintf(paste(
      "G*Power's allocation rule gives N = %s (n1 = %s, n2 = %s), whose power %s misses the target %s, and no",
      "larger N up to %s has a split that reaches it"), format(n_star, scientific = FALSE), format(gs[["n1"]]),
      format(gs[["n2"]]), format(signif(p_int, 7)), format(target), format(n_max, scientific = FALSE)),
      diagnostics = diag))
  }
  g2 <- .ig_design_z_correlations_two_independent(p, n_adj)
  note <- sprintf(paste(
    "G*Power's allocation rule for N2/N1 = %s gives N = %s (n1 = %s, n2 = %s), whose power %s is below the target",
    "%s (G*Power 3.1.9.7 reports it without a warning); inferguidR reports the smallest larger N whose split",
    "reaches the target, N = %s (n1 = %s, n2 = %s)"),
    .ig_format_double(r), format(n_star, scientific = FALSE), format(gs[["n1"]]), format(gs[["n2"]]),
    format(signif(p_int, 7)), format(target), format(n_adj, scientific = FALSE), format(g2[["n1"]]), format(g2[["n2"]]))
  diag$adjusted_from_gpower_rule <- TRUE
  diag$adjustment_note <- note
  diag$minimality <- sprintf(paste(
    "N = %s is the smallest N above G*Power's N = %s whose split reaches the target (the integer-design power is",
    "non-decreasing in N)"), format(n_adj, scientific = FALSE), format(n_star, scientific = FALSE))
  n_int <- smallest_int(first)
  if (!is.na(n_int) && n_int != n_adj) diag$integer_criterion_n_total <- n_int
  list(status = "ok", n = n_adj, power = int_power(n_adj), diagnostics = diag, warnings = note)
}

# ===========================================================================================
# Correlations: Two dependent Pearson r's (manual ch. 28)
# ===========================================================================================
#
# Steiger (1980) Z1* (common index) and Z2* (no common index) with the Fisher
# z transform. With c0 the covariance of the transformed correlations under
# H0 (H1 correlation set equal to rho_ab) and c1 under H1 (manual 28.5.3):
#
#   s0 = sqrt((2 - 2 c0) / (N - 3)),  Z = (z_H1 - z_ab) / s0 ~ N(mu, sd1) under H1
#   mu = (z_H1 - z_ab) / s0,  sd1 = sqrt((2 - 2 c1) / (N - 3)) / s0 = sqrt((1 - c1) / (1 - c0))
#
# Manual 28.5.3 prints c1 = c_ab;ac for the general case; the covariance of
# the two tested correlations, c_ab;cd at the H1 matrix, is used (the live
# captures, for example GPW-Z-DEPNCI-POS-B001-O01, reproduce only with it).

# Steiger (1980) eq. 3 / manual eq. 13 divided by (1 - rho_ab^2)(1 - rho_ac^2).
.ig_g7_c_common <- function(rab, rac, rbc) {
  (rbc * (1 - rab^2 - rac^2) - 0.5 * rab * rac * (1 - rab^2 - rac^2 - rbc^2)) / ((1 - rab^2) * (1 - rac^2))
}

# Steiger (1980) eq. 2 / manual eq. 12 divided by (1 - rho_ab^2)(1 - rho_cd^2).
.ig_g7_c_no_common <- function(rab, rac, rad, rbc, rbd, rcd) {
  0.5 * ((rac - rab * rbc) * (rbd - rbc * rcd) + (rad - rac * rcd) * (rbc - rab * rac) +
           (rac - rad * rcd) * (rbd - rab * rad) + (rad - rab * rbd) * (rbc - rbd * rcd)) /
    ((1 - rab^2) * (1 - rcd^2))
}

# Correlation matrices (manual 28: C2 for a common index, C1 otherwise) with
# the H1 correlation x.
.ig_g7_cor_common <- function(p, x) {
  matrix(c(1, p$rho_ab, x,
           p$rho_ab, 1, p$rho_bc,
           x, p$rho_bc, 1), 3, 3)
}

.ig_g7_cor_no_common <- function(p, x) {
  matrix(c(1, p$rho_ab, p$rho_ac, p$rho_ad,
           p$rho_ab, 1, p$rho_bc, p$rho_bd,
           p$rho_ac, p$rho_bc, 1, x,
           p$rho_ad, p$rho_bd, x, 1), 4, 4)
}

# Interval of H1 correlations x with a positive semi-definite matrix (the
# smallest eigenvalue is concave in x, so the set is an interval containing
# the H0 value when the H0 matrix is valid), intersected with
# [-1 + eps, 1 - eps].
.ig_g7_psd_interval <- function(build, x0) {
  lam <- function(x) .ig_g7_min_eigen(build(x)) + .ig_g7_psd_tol
  lim <- 1 - .ig_g7_rho_eps
  end <- function(to) {
    if (lam(to) >= 0) return(to)
    stats::uniroot(lam, sort(c(x0, to)), tol = 1e-14)$root
  }
  c(lower = end(-lim), upper = end(lim))
}

.ig_g7_dep_spec <- function(kind) {
  if (kind == "common") {
    list(h1 = "rho_ac", build = .ig_g7_cor_common,
         cov = function(p, x) .ig_g7_c_common(p$rho_ab, x, p$rho_bc),
         h1_label = "rho_ac", ge_code = "rho_ac_ge_rho_ab")
  } else {
    list(h1 = "rho_cd", build = .ig_g7_cor_no_common,
         cov = function(p, x) .ig_g7_c_no_common(p$rho_ab, p$rho_ac, p$rho_ad, p$rho_bc, p$rho_bd, x),
         h1_label = "rho_cd", ge_code = "rho_cd_ge_rho_ab")
  }
}

.ig_power_z_dependent <- function(p, kind) {
  sp <- .ig_g7_dep_spec(kind)
  n <- as.double(p$n_total)
  x_in <- p[[sp$h1]]
  x <- x_in
  notes <- list()
  seg <- NULL
  if (!is.null(p$effect_direction)) {
    build <- function(v) sp$build(p, v)
    iv <- .ig_g7_psd_interval(build, p$rho_ab)
    seg <- if (grepl("_le_", p$effect_direction, fixed = TRUE)) c(iv[["lower"]], p$rho_ab) else c(p$rho_ab, iv[["upper"]])
  }
  c0 <- sp$cov(p, p$rho_ab)
  c1 <- sp$cov(p, x)
  z_ab <- atanh(p$rho_ab)
  z_h1 <- atanh(x)
  s0 <- sqrt((2 - 2 * c0) / (n - 3))
  mu <- (z_h1 - z_ab) / s0
  sd1 <- sqrt(max(1 - c1, 0) / (1 - c0))
  side <- .ig_g7_side(z_h1 - z_ab, p)
  zp <- if (c0 < 1) .ig_g7_z_power(mu, sd1, p$alpha, p$tails, side) else
    list(power = NA_real_, beta = NA_real_, critical = NA_real_, side = side)
  if (is.na(zp$power)) {
    notes$power <- notes$beta <- notes$critical <- notes$ncp <- paste(
      "the H0 covariance of the two transformed correlations is 1 (for example rho_bc = 1), so the",
      "standard error of their difference is 0 and the test statistic is undefined")
  }
  outputs <- list(c0 = c0, c1 = c1, se_h0 = s0, sd_h1 = sd1, z_h0 = z_ab, z_h1 = z_h1)
  if (!is.null(seg)) {
    outputs$admissible_segment <- seg
    outputs$evaluated_at <- x
  }
  list(
    power = zp$power,
    beta = zp$beta,
    alpha_attained = p$alpha,
    critical = zp$critical,
    ncp = mu,
    n_total = n,
    statistic_distribution = list(
      h0 = list(family = "normal", mean = 0, sd = 1),
      h1 = list(family = "normal", mean = mu, sd = sd1),
      rejection_side = if (identical(p$tails, "two")) "both" else if (zp$side > 0) "upper" else "lower"
    ),
    outputs = outputs,
    notes = notes,
    assumptions = c(
      sprintf("multivariate normal data; Steiger (1980) %s test of rho_ab = %s with Fisher-z transformed sample correlations",
              if (kind == "common") "Z1*" else "Z2*", sp$h1_label),
      paste("H0 covariance c0 at the H1 correlation set equal to rho_ab, H1 covariance c1 at the H1 matrix;",
            "s0 = sqrt((2 - 2 c0)/(N - 3)), H1: N((z_H1 - z_ab)/s0, sqrt((1 - c1)/(1 - c0))) (manual 28.5.3)"),
      .ig_g7_tail_note(p$tails, sprintf("the sign of %s - rho_ab", sp$h1_label))
    )
  )
}

.ig_power_z_correlations_two_dependent_common_index <- function(p) .ig_power_z_dependent(p, "common")

.ig_power_z_correlations_two_dependent_no_common_index <- function(p) .ig_power_z_dependent(p, "no_common")

# Cross-parameter checks: valid correlation matrices under H0 (H1
# correlation = rho_ab) and, when the H1 correlation is an input, under H1
# (manual 28: "The correlation matrix is not valid, that is, not positive
# semi-definite"); target power not below alpha.
.ig_g7_dep_constraints <- function(kind) {
  force(kind)
  function(p) {
    sp <- .ig_g7_dep_spec(kind)
    need <- if (kind == "common") c("rho_ab", "rho_bc") else c("rho_ab", "rho_ac", "rho_ad", "rho_bc", "rho_bd")
    out <- .ig_g7_power_alpha_constraint(p)
    if (!all(vapply(need, function(nm) is.numeric(p[[nm]]), logical(1)))) return(out)
    m0 <- .ig_g7_min_eigen(sp$build(p, p$rho_ab))
    if (m0 < -.ig_g7_psd_tol) {
      out <- c(out, sprintf(paste("the correlation matrix under H0 (%s = rho_ab) is not positive semi-definite",
                                  "(smallest eigenvalue %s)"), sp$h1_label, format(signif(m0, 4))))
    }
    x <- p[[sp$h1]]
    if (is.numeric(x)) {
      m1 <- .ig_g7_min_eigen(sp$build(p, x))
      if (m1 < -.ig_g7_psd_tol) {
        out <- c(out, sprintf("the correlation matrix under H1 is not positive semi-definite (smallest eigenvalue %s)",
                              format(signif(m1, 4))))
      }
    }
    out
  }
}

.ig_g7_rho_param <- function(name, role, label, description, aliases = character(), open = TRUE) {
  # `open`: the H0 and H1 correlations, which G*Power admits in
  # [-0.999999, 0.999999] only; the other correlations may be +-1.
  ig_param(name, role = role,
           domain = if (open) ig_domain_real(-(1 - .ig_g7_rho_eps), 1 - .ig_g7_rho_eps) else
             ig_domain_real(-1, 1),
           label = label, description = description, aliases = aliases)
}

.ig_method_z_correlations_two_dependent_common_index <- function() {
  ig_method(
    id = "z_correlations_two_dependent_common_index",
    family = "z",
    label = list(en = "Correlations: Two dependent Pearson r's (common index)",
                 de = "Korrelationen: Zwei abh\u00e4ngige Pearson-r (gemeinsamer Index)"),
    params = ig_params(
      .ig_g7_rho_param("rho_ac", "effect", list(en = "H1 corr \u03c1_ac", de = "H1-Korrelation \u03c1_ac"),
                       "Correlation rho_ac assumed under H1 (the effect size of this procedure, manual 28.1).",
                       aliases = "h1_corr_rho_ac"),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      ig_param("n_total", role = "design", domain = ig_domain_integer(lower = 4), unit = "subjects",
               label = list(en = "Sample size", de = "Stichprobenumfang"),
               description = "Number of cases (length of the data vectors a, b, c); at least 4 (N - 3 > 0).",
               aliases = "sample_size"),
      .ig_g7_rho_param("rho_ab", "nuisance", list(en = "H0 corr \u03c1_ab", de = "H0-Korrelation \u03c1_ab"),
                       "Correlation rho_ab, the value of rho_ac under H0.", aliases = "h0_corr_rho_ab"),
      .ig_g7_rho_param("rho_bc", "nuisance", list(en = "Corr \u03c1_bc", de = "Korrelation \u03c1_bc"),
                       "Correlation between b and c (the dependence of the two tested correlations).",
                       aliases = "corr_rho_bc", open = FALSE)
    ),
    options = ig_options(
      ig_option("tails", values = c("one", "two"), default = "one",
                labels = c(one = "One", two = "Two"),
                label = list(en = "Tail(s)", de = "Seiten")),
      ig_option("effect_direction", values = c("rho_ac_le_rho_ab", "rho_ac_ge_rho_ab"),
                default = "rho_ac_le_rho_ab",
                labels = c(rho_ac_le_rho_ab = "\u03c1_ac \u2264 \u03c1_ab", rho_ac_ge_rho_ab = "\u03c1_ac \u2265 \u03c1_ab"),
                applies_when = list(mode = "sensitivity"),
                label = list(en = "Effect direction", de = "Effektrichtung"))
    ),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = .ig_power_z_correlations_two_dependent_common_index,
    design = NULL,
    discrete = FALSE,
    effect_direction = .ig_g7_sensitivity_direction,
    sensitivity_domain = function(p) .ig_g7_sensitivity_domain(p, "common"),
    engine = list(package = "stats", fun = "pnorm/qnorm; Steiger (1980) covariance in closed form",
                  version = NA_character_),
    checker = "validation/independent/g7/check_z_correlations_two_dependent.R",
    reference = list(
      gpower_path = "Tests > Correlation and regression > Correlations: Two dependent Pearson r's (common index)",
      gpower_label = "z tests - Correlations: Two dependent Pearson r's (common index)",
      manual_pages = "ch. 28 (PDF pp. 66-69); Steiger (1980) eqs. 3, 9-12",
      method_notes = "specs/method-notes/z_correlations_two_dependent_common_index.md"
    ),
    status = "independently_checked",
    constraints = .ig_g7_dep_constraints("common")
  )
}

.ig_method_z_correlations_two_dependent_no_common_index <- function() {
  ig_method(
    id = "z_correlations_two_dependent_no_common_index",
    family = "z",
    label = list(en = "Correlations: Two dependent Pearson r's (no common index)",
                 de = "Korrelationen: Zwei abh\u00e4ngige Pearson-r (kein gemeinsamer Index)"),
    params = ig_params(
      .ig_g7_rho_param("rho_cd", "effect", list(en = "H1 corr \u03c1_cd", de = "H1-Korrelation \u03c1_cd"),
                       "Correlation rho_cd assumed under H1 (the effect size of this procedure, manual 28.1).",
                       aliases = "h1_corr_rho_cd"),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      ig_param("n_total", role = "design", domain = ig_domain_integer(lower = 4), unit = "subjects",
               label = list(en = "Sample size", de = "Stichprobenumfang"),
               description = "Number of cases (length of the data vectors a, b, c, d); at least 4 (N - 3 > 0).",
               aliases = "sample_size"),
      .ig_g7_rho_param("rho_ab", "nuisance", list(en = "H0 corr \u03c1_ab", de = "H0-Korrelation \u03c1_ab"),
                       "Correlation rho_ab, the value of rho_cd under H0.", aliases = "h0_corr_rho_ab"),
      .ig_g7_rho_param("rho_ac", "nuisance", list(en = "Corr \u03c1_ac", de = "Korrelation \u03c1_ac"),
                       "Correlation between a and c.", aliases = "corr_rho_ac", open = FALSE),
      .ig_g7_rho_param("rho_ad", "nuisance", list(en = "Corr \u03c1_ad", de = "Korrelation \u03c1_ad"),
                       "Correlation between a and d.", aliases = "corr_rho_ad", open = FALSE),
      .ig_g7_rho_param("rho_bc", "nuisance", list(en = "Corr \u03c1_bc", de = "Korrelation \u03c1_bc"),
                       "Correlation between b and c.", aliases = "corr_rho_bc", open = FALSE),
      .ig_g7_rho_param("rho_bd", "nuisance", list(en = "Corr \u03c1_bd", de = "Korrelation \u03c1_bd"),
                       "Correlation between b and d.", aliases = "corr_rho_bd", open = FALSE)
    ),
    options = ig_options(
      ig_option("tails", values = c("one", "two"), default = "one",
                labels = c(one = "One", two = "Two"),
                label = list(en = "Tail(s)", de = "Seiten")),
      ig_option("effect_direction", values = c("rho_cd_le_rho_ab", "rho_cd_ge_rho_ab"),
                default = "rho_cd_ge_rho_ab",
                labels = c(rho_cd_le_rho_ab = "\u03c1_cd \u2264 \u03c1_ab", rho_cd_ge_rho_ab = "\u03c1_cd \u2265 \u03c1_ab"),
                applies_when = list(mode = "sensitivity"),
                label = list(en = "Effect direction", de = "Effektrichtung"))
    ),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = .ig_power_z_correlations_two_dependent_no_common_index,
    design = NULL,
    discrete = FALSE,
    effect_direction = .ig_g7_sensitivity_direction,
    sensitivity_domain = function(p) .ig_g7_sensitivity_domain(p, "no_common"),
    engine = list(package = "stats", fun = "pnorm/qnorm; Steiger (1980) covariance in closed form",
                  version = NA_character_),
    checker = "validation/independent/g7/check_z_correlations_two_dependent.R",
    reference = list(
      gpower_path = "Tests > Correlation and regression > Correlations: Two dependent Pearson r's (no common index)",
      gpower_label = "z tests - Correlations: Two dependent Pearson r's (no common index)",
      manual_pages = "ch. 28 (PDF pp. 66-69); Steiger (1980) eqs. 2, 9-12",
      method_notes = "specs/method-notes/z_correlations_two_dependent_no_common_index.md"
    ),
    status = "independently_checked",
    constraints = .ig_g7_dep_constraints("no_common")
  )
}

# ===========================================================================================
# Correlation: Tetrachoric model (manual ch. 31)
# ===========================================================================================
#
# Latent standard bivariate normal (X, Y) with correlation rho, dichotomised
# at h and k with P(X > h) = marginal prop x and P(Y > k) = marginal prop y
# (manual 31.1). Table orientation of the manual: rows Y = 0, 1; columns
# X = 0, 1; a = p11 = P(X <= h, Y <= k), b = p12, c = p21, d = p22.
#
# Wald test W = (r - rho0) / se0(r) (manual 31.0.2); H0 N(0, 1), H1 N(m1, s1)
# with m1 = (rho - rho0) / s_r0, s1 = s_r1 / s_r0 (manual 31.5.3), where s_r
# is the asymptotic standard error at the cell probabilities of the model:
#
# * exact mode (Brown & Benedetti 1977): Pearson's (1901) asymptotic standard
#   error of the tetrachoric r (manual 31.5.1). The engine evaluates it as the
#   variance of the influence function of r (cells weighted by p_ij), which is
#   algebraically the manual's formula with its first term corrected to
#   (a + d)(b + c)/4 and Phi_1, Phi_2 read as normal CDFs, and has no
#   cancellation when a cell probability is tiny. The H0 and H1 correlations
#   are the inputs ("In the exact mode it is not necessary to estimate r in
#   order to calculate power, because it is already given in the input").
# * approximation mode (Bonett & Price 2005, manual 31.5.2): rho* =
#   cos(pi / (1 + omega^c)), s(r*) = k sqrt(sum 1/(N p_ij)) at the model cell
#   probabilities; m1 = (rho*_1 - rho*_0) / s*_0, s1 = s*_1 / s*_0.

.ig_method_z_tetrachoric_correlation <- function() {
  prob <- function(name, label, description, aliases) {
    ig_param(name, role = "nuisance",
             domain = ig_domain_real(0, 1, lower_inclusive = FALSE, upper_inclusive = FALSE),
             label = label, description = description, aliases = aliases)
  }
  ig_method(
    id = "z_tetrachoric_correlation",
    family = "z",
    label = list(en = "Correlation: Tetrachoric model", de = "Korrelation: Tetrachorisches Modell"),
    params = ig_params(
      .ig_g7_rho_param("rho_h1", "effect", list(en = "H1 corr \u03c1", de = "H1-Korrelation \u03c1"),
                       "Tetrachoric correlation assumed under H1 (the effect size of this procedure, manual 31.1).",
                       aliases = "h1_corr_rho"),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      ig_param("n_total", role = "design", domain = ig_domain_integer(lower = 2), unit = "subjects",
               label = list(en = "Total sample size", de = "Gesamtstichprobenumfang"),
               description = "Number of observations in the 2 x 2 table; at least 2 (as in G*Power).",
               aliases = "total_sample_size"),
      .ig_g7_rho_param("rho_h0", "nuisance", list(en = "H0 corr \u03c1", de = "H0-Korrelation \u03c1"),
                       "Tetrachoric correlation rho0 under H0.", aliases = "h0_corr_rho"),
      prob("marginal_x", list(en = "Marginal prop x", de = "Randwahrscheinlichkeit x"),
           "Marginal probability P(X = 1) = P(X > h) = p12 + p22 (manual 31.1).", "marginal_prop_x"),
      prob("marginal_y", list(en = "Marginal prop y", de = "Randwahrscheinlichkeit y"),
           "Marginal probability P(Y = 1) = P(Y > k) = p21 + p22 (manual 31.1).", "marginal_prop_y")
    ),
    options = ig_options(
      ig_option("tails", values = c("one", "two"), default = "one",
                labels = c(one = "One", two = "Two"),
                label = list(en = "Tail(s)", de = "Seiten")),
      ig_option("computation_of_correlation", values = c("exact_brown_benedetti", "approx_bonett_price"),
                default = "exact_brown_benedetti",
                labels = c(exact_brown_benedetti = "Use exact computation (Brown and Benedetti, 1977)",
                           approx_bonett_price = "Use approximation (Bonett and Price, 2005)"),
                label = list(en = "Computation of correlation", de = "Berechnung der Korrelation")),
      ig_option("effect_direction", values = c("rho_h1_le_h0", "rho_h1_ge_h0"),
                default = "rho_h1_ge_h0",
                labels = c(rho_h1_le_h0 = "\u03c1 H1 \u2264 \u03c1 H0", rho_h1_ge_h0 = "\u03c1 H1 \u2265 \u03c1 H0"),
                applies_when = list(mode = "sensitivity"),
                label = list(en = "Effect direction", de = "Effektrichtung"))
    ),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = .ig_power_z_tetrachoric_correlation,
    design = NULL,
    discrete = FALSE,
    effect_direction = .ig_g7_sensitivity_direction,
    sensitivity_domain = function(p) .ig_g7_sensitivity_domain(p, NULL),
    engine = list(package = "stats",
                  fun = paste("pnorm/qnorm/dnorm; bivariate normal cell probabilities by stats::integrate of",
                              "phi(x) Phi((k - rho x)/sqrt(1 - rho^2)); asymptotic SE in influence-function form"),
                  version = NA_character_),
    checker = "validation/independent/g7/check_z_tetrachoric_correlation.R",
    reference = list(
      gpower_path = "Tests > Correlation and regression > Correlation: Tetrachoric model",
      gpower_label = "z tests - Correlation: Tetrachoric model",
      manual_pages = "ch. 31 (PDF pp. 80-83); manual defects in specs/method-notes/z_tetrachoric_correlation.md",
      method_notes = "specs/method-notes/z_tetrachoric_correlation.md"
    ),
    status = "independently_checked",
    constraints = .ig_g7_power_alpha_constraint
  )
}

# Bivariate normal CDF P(X <= h, Y <= k) for correlation rho in (-1, 1), by
# one-dimensional integration of phi(x) Phi((k - rho x) / sqrt(1 - rho^2))
# over x <= h (integrating over the variable with the smaller marginal). The
# integrand has a step of width sqrt(1 - rho^2)/|rho| at x = k / rho; the range
# is split there. Each cell of a 2 x 2 table is computed as its own orthant
# probability (never as a difference), so small cells keep their relative
# accuracy.
.ig_g7_bvn <- function(h, k, rho) {
  if (h == -Inf || k == -Inf) return(0)
  if (h == Inf) return(stats::pnorm(k))
  if (k == Inf) return(stats::pnorm(h))
  if (rho == 0) return(stats::pnorm(h) * stats::pnorm(k))
  if (stats::pnorm(k) < stats::pnorm(h)) {
    tmp <- h
    h <- k
    k <- tmp
  }
  s <- sqrt((1 - rho) * (1 + rho))
  f <- function(x) stats::dnorm(x) * stats::pnorm((k - rho * x) / s)
  brk <- k / rho
  ends <- if (is.finite(brk) && brk < h) c(-Inf, brk, h) else c(-Inf, h)
  total <- 0
  for (i in seq_len(length(ends) - 1L)) {
    v <- tryCatch(stats::integrate(f, ends[i], ends[i + 1L], rel.tol = 1e-12, abs.tol = 0,
                                   subdivisions = 1000L)$value,
                  error = function(e) stats::integrate(f, ends[i], ends[i + 1L], rel.tol = 1e-9, abs.tol = 0,
                                                       subdivisions = 1000L)$value)
    total <- total + v
  }
  total
}

# Cell probabilities of the dichotomised bivariate normal: thresholds
# h = qnorm(1 - px) (column X), k = qnorm(1 - py) (row Y).
.ig_g7_tetra_cells <- function(rho, px, py) {
  h <- stats::qnorm(px, lower.tail = FALSE)
  k <- stats::qnorm(py, lower.tail = FALSE)
  if (rho == 0) {
    pr <- c(a = (1 - px) * (1 - py), b = px * (1 - py), c = (1 - px) * py, d = px * py)
  } else {
    pr <- c(a = .ig_g7_bvn(h, k, rho), b = .ig_g7_bvn(-h, k, -rho),
            c = .ig_g7_bvn(h, -k, -rho), d = .ig_g7_bvn(-h, -k, rho))
  }
  list(p = pr, h = h, k = k)
}

.ig_g7_phi2 <- function(h, k, rho) {
  s2 <- (1 - rho) * (1 + rho)
  exp(-(h^2 - 2 * rho * h * k + k^2) / (2 * s2)) / (2 * pi * sqrt(s2))
}

# Asymptotic standard error of the tetrachoric r (Pearson 1901; Brown &
# Benedetti 1977) at cell probabilities `pr` and correlation `rho`. With
# A = Phi((k - rho h)/s) and B = Phi((h - rho k)/s) (s = sqrt(1 - rho^2)), the
# influence function of r is g / phi2(h, k; rho) with g = I_a - A I_{X=0} -
# B I_{Y=0}, so Var(r) = sum p (g - gbar)^2 / (N phi2^2). This is the
# variance written out in manual 31.5.1 (first term (a + d)(b + c)/4).
.ig_g7_tetra_se_exact <- function(pr, h, k, rho, n) {
  s <- sqrt((1 - rho) * (1 + rho))
  ua <- (k - rho * h) / s
  ub <- (h - rho * k) / s
  A <- stats::pnorm(ua)
  B <- stats::pnorm(ub)
  # 1 - A - B without cancellation.
  ga <- if (A < B) stats::pnorm(ub, lower.tail = FALSE) - A else stats::pnorm(ua, lower.tail = FALSE) - B
  g <- c(a = ga, b = -B, c = -A, d = 0)
  gbar <- sum(pr * g)
  v <- sum(pr * (g - gbar)^2)
  dens <- .ig_g7_phi2(h, k, rho)
  sqrt(v / n) / dens
}

# Bonett & Price (2005) approximation at cell probabilities `pr` (manual
# 31.5.2): rho* = cos(pi/(1 + omega^c)), omega = ad/(bc),
# c = (1 - |p1. - p.1|/5 - (1/2 - p_min)^2)/2, s(ln omega) = sqrt(sum 1/(N p)),
# s(r*) = k s(ln omega) with k = pi c w sin(pi/(1 + w))/(1 + w)^2, w = omega^c.
.ig_g7_tetra_bp <- function(pr, px, py, n) {
  pm <- min(px, py, 1 - px, 1 - py)
  cexp <- (1 - abs(px - py) / 5 - (0.5 - pm)^2) / 2
  omega <- (pr[["a"]] * pr[["d"]]) / (pr[["b"]] * pr[["c"]])
  w <- omega^cexp
  sl <- sqrt(sum(1 / pr) / n)
  r <- cos(pi / (1 + w))
  kk <- pi * cexp * w * sin(pi / (1 + w)) / (1 + w)^2
  list(r = r, se = kk * sl, se_log_omega = sl, omega = omega, c = cexp)
}

.ig_power_z_tetrachoric_correlation <- function(p) {
  n <- as.double(p$n_total)
  px <- p$marginal_x
  py <- p$marginal_y
  r0 <- p$rho_h0
  r1_in <- p$rho_h1
  r1 <- r1_in
  seg <- NULL
  if (!is.null(p$effect_direction)) {
    lim <- 1 - .ig_g7_rho_eps
    seg <- if (grepl("_le_", p$effect_direction, fixed = TRUE)) c(-lim, r0) else c(r0, lim)
  }
  exact <- !identical(p$computation_of_correlation, "approx_bonett_price")
  c0 <- .ig_g7_tetra_cells(r0, px, py)
  c1 <- if (r1 == r0) c0 else .ig_g7_tetra_cells(r1, px, py)
  outputs <- list(computation = if (exact) "exact_brown_benedetti" else "approx_bonett_price",
                  thresholds = c(h = c0$h, k = c0$k), cells_h0 = c0$p, cells_h1 = c1$p)
  if (exact) {
    est0 <- r0
    est1 <- r1
    se0 <- .ig_g7_tetra_se_exact(c0$p, c0$h, c0$k, r0, n)
    se1 <- .ig_g7_tetra_se_exact(c1$p, c1$h, c1$k, r1, n)
  } else {
    b0 <- .ig_g7_tetra_bp(c0$p, px, py, n)
    b1 <- if (r1 == r0) b0 else .ig_g7_tetra_bp(c1$p, px, py, n)
    est0 <- b0$r
    est1 <- b1$r
    se0 <- b0$se
    se1 <- b1$se
    outputs$omega <- c(h0 = b0$omega, h1 = b1$omega)
    outputs$bp_exponent_c <- b0$c
  }
  mu <- (est1 - est0) / se0
  sd1 <- se1 / se0
  side <- .ig_g7_side(est1 - est0, p)
  zp <- .ig_g7_z_power(mu, sd1, p$alpha, p$tails, side)
  notes <- list()
  if (is.na(zp$power)) {
    notes$power <- notes$beta <- notes$critical <- notes$ncp <- paste(
      "the standard error of r is not finite at these cell probabilities (a cell probability of zero or a",
      "correlation too close to +-1 for the margins)")
  }
  crit_r <- if (is.na(zp$power)) NA_real_ else est0 + zp$critical * se0
  outputs$rho_h1_used <- est1
  outputs$rho_h0_used <- est0
  outputs$se_h0 <- se0
  outputs$se_h1 <- se1
  outputs$critical_r <- crit_r
  if (!exact && identical(p$tails, "two") && !is.na(zp$power)) {
    # The Bonett-Price interval of r* under H0 at the two-sided critical z
    # (manual 31.5.2, confidence interval formula with omega and s(ln omega) of
    # the H0 table). G*Power 3.1.9.7 displays these bounds as "Critical r* lwr /
    # upr" in two-tailed approximation mode; they are not the rejection region
    # of the Wald test whose power is computed (that region is critical_r).
    zc <- zp$critical[2]
    lw <- exp(log(b0$omega) - zc * b0$se_log_omega)
    up <- exp(log(b0$omega) + zc * b0$se_log_omega)
    outputs$bonett_price_interval_h0 <- c(cos(pi / (1 + lw^b0$c)), cos(pi / (1 + up^b0$c)))
  }
  if (!is.null(seg)) {
    outputs$admissible_segment <- seg
    outputs$evaluated_at <- r1
  }
  list(
    power = zp$power,
    beta = zp$beta,
    alpha_attained = p$alpha,
    critical = zp$critical,
    ncp = mu,
    n_total = n,
    statistic_distribution = list(
      h0 = list(family = "normal", mean = 0, sd = 1),
      h1 = list(family = "normal", mean = mu, sd = sd1),
      rejection_side = if (identical(p$tails, "two")) "both" else if (zp$side > 0) "upper" else "lower"
    ),
    outputs = outputs,
    notes = notes,
    assumptions = c(
      "2 x 2 table from dichotomising a standard bivariate normal (X, Y) at fixed thresholds (tetrachoric model)",
      if (exact) {
        paste("exact mode: Wald test with Pearson's / Brown and Benedetti's asymptotic standard error at the",
              "model cell probabilities of the H0 and H1 correlations (manual 31.5.1, first term corrected)")
      } else {
        paste("approximation mode: Bonett and Price (2005) r* = cos(pi/(1 + omega^c)) and s(r*) at the model",
              "cell probabilities (manual 31.5.2); m1 = (r*1 - r*0)/s*0, s1 = s*1/s*0")
      },
      .ig_g7_tail_note(p$tails, "the sign of rho H1 - rho H0")
    )
  )
}
