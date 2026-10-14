# Method definitions: Wilcoxon tests (t family, group G10) ----------------------------
#
# Ownership: group G10 (blueprint phases 06-10 for the three Wilcoxon
# procedures of the G*Power "t tests" family). Method notes:
#   specs/method-notes/t_wilcoxon_signed_rank_one_sample.md
#   specs/method-notes/t_wilcoxon_signed_rank_matched_pairs.md
#   specs/method-notes/t_wilcoxon_mann_whitney.md
#
# G*Power 3.1.9.7 offers two large-sample power methods for each test
# (Options dialog), both reproduced here:
#
# * A.R.E. method ("Asymptotic Relative Efficiency (based on t test)"): the
#   power of the t test at the effective sample size N k, with k the
#   asymptotic relative efficiency of the Wilcoxon test against the t test
#   for the parent distribution (normal 3/pi, Laplace 3/2, logistic pi^2/9,
#   the lower bound 0.864, or a value specified manually, at least 0.864):
#     signed-rank:  t(N k - 1, d sqrt(N k))
#     rank-sum:     t((n1 + n2) k - 2, d sqrt(n1 n2 k / (n1 + n2)))
#   The degrees of freedom are not whole numbers. The tail convention is
#   that of G*Power's t tests (manual 25.5): two tails count both rejection
#   regions; one tail rejects on the side of the sign of d.
#
# * Lehmann method ("Asymptotic Power (Lehmann method based on z test)"):
#   the normal approximation of the rank statistic with its exact mean and
#   variance under the shift alternative (Lehmann 1975; Shieh, Jan & Randles
#   2007), written with the moments p1, p2, p3 of the parent distribution
#   (or moments given manually). Signed-rank statistic V (sum of the ranks of
#   the positive values, equivalently the number of positive Walsh averages):
#     E0 = N(N + 1)/4, V0 = N(N + 1)(2N + 1)/24,
#     E1 = N(N - 1) p2 / 2 + N p1,
#     V1 = N(N - 1)(N - 2)(p3 - p2^2) + N(N - 1)[2(p1 - p2)^2 + 3 p2 (1 - p2)]/2
#          + N p1 (1 - p1),
#   with p1 = P(X > 0), p2 = P(X1 + X2 > 0), p3 = P(X1 + X2 > 0, X1 + X3 > 0).
#   Rank-sum statistic U (Mann-Whitney count, group sizes m = n1, n = n2):
#     E0 = mn/2, V0 = mn(m + n + 1)/12,
#     E1 = mn p1,
#     V1 = mn[p1(1 - p1) + (n - 1)(p2 - p1^2) + (m - 1)(p3 - p1^2)],
#   with p1 = P(X < Y), p2 = P(X < Y, X < Y'), p3 = P(X < Y, X' < Y)
#   (X from group 1, Y from group 2). The test rejects when the statistic
#   exceeds c = E0 + z sigma0 (plus 0.5 with continuity correction); two
#   tails add the lower region E0 - z sigma0 (minus 0.5). The power is the
#   probability of the rejection region under N(E1, V1).
#
#   Three corrections of the printed manual formulas (ch. 22-24) were needed
#   to reproduce G*Power (specs/engine-decisions.md G3; confirmed on 3,100
#   emulator captures): the moments are P(X > 0)-type quantities, the
#   variance uses (p3 - p2^2), and the continuity correction is applied
#   conservatively (c + 0.5).
#
# Conventions found in the captures (validation/compare/<method_id>/):
# * Lehmann with a predefined parent uses |d|: the moments and the power of
#   d = -0.5 equal those of d = 0.5 (both rejection regions of the
#   symmetric parent agree). Manually specified moments are used as given.
# * Effect-size parameter p1 (Lehmann, predefined parent): the shift is the
#   one whose p1 under the parent equals the input, and d = shift /
#   "Relative SD"; all moments and the power are then those of d.
# * A priori (rank-sum): even total N, n1 = floor(N/(1 + r) + 1/2),
#   n2 = N - n1 (as for G*Power's two-group t test). The smallest such
#   design whose power reaches the target is reported; G*Power chooses N by
#   the power of the fractional group sizes instead, which can give a
#   different N (documented in the method note; not copied).

.ig_define_methods_t_g10 <- function() {
  list(.ig_method_wilcoxon("t_wilcoxon_signed_rank_one_sample"),
       .ig_method_wilcoxon("t_wilcoxon_signed_rank_matched_pairs"),
       .ig_method_wilcoxon("t_wilcoxon_mann_whitney"))
}

# Asymptotic relative efficiencies of the Wilcoxon tests against the t test
# for the predefined parents (Lehmann 1975; manual ch. 22-24). "min_are" is
# the lower bound 108/125 = 0.864 over all continuous symmetric parents
# (Hodges & Lehmann 1956).
.ig_wilcoxon_are <- c(normal = 3 / pi, laplace = 1.5, logistic = pi^2 / 9, min_are = 0.864)

# Smallest A.R.E. accepted for a manually specified value (G*Power: "The
# minimal A.R.E. of the Wilcoxon test versus the t test is 0.864.").
.ig_wilcoxon_min_are <- 0.864

# Per-method identity.
.ig_wilcoxon_info <- function(id) {
  switch(id,
    t_wilcoxon_signed_rank_one_sample = list(
      id = id, kind = "signed_rank", effect = "d",
      label = list(en = "Means: Wilcoxon signed-rank test (one sample case)",
                   de = "Mittelwerte: Wilcoxon-Vorzeichen-Rang-Test (Einstichprobenfall)"),
      effect_label = list(en = "Effect size d", de = "Effektst\u00e4rke d"),
      effect_description = paste("Standardised shift d = (mu - mu0) / sigma of a distribution symmetric about",
                                 "mu (manual ch. 22; analogous to the one-sample t test)."),
      gpower_path = "Tests > Means > One group: Wilcoxon (non-parametric)",
      gpower_label = "t tests - Means: Wilcoxon signed-rank test (one sample case)",
      manual_pages = "ch. 22 (PDF p. 53); formulas ch. 23.5 (PDF pp. 57-58); A.R.E. table ch. 23",
      notes = "specs/method-notes/t_wilcoxon_signed_rank_one_sample.md",
      checker = "validation/independent/g10/check_wilcoxon.R"),
    t_wilcoxon_signed_rank_matched_pairs = list(
      id = id, kind = "signed_rank", effect = "dz",
      label = list(en = "Means: Wilcoxon signed-rank test (matched pairs)",
                   de = "Mittelwerte: Wilcoxon-Vorzeichen-Rang-Test (verbundene Paare)"),
      effect_label = list(en = "Effect size dz", de = "Effektst\u00e4rke dz"),
      effect_description = paste("Standardised mean difference dz = mu_z / sigma_z of the paired differences",
                                 "z = x - y (shift model; manual ch. 23)."),
      gpower_path = "Tests > Means > Two dependent groups (matched pairs): Wilcoxon (non-parametric)",
      gpower_label = "t tests - Means: Wilcoxon signed-rank test (matched pairs)",
      manual_pages = "ch. 23 (PDF pp. 56-58)",
      notes = "specs/method-notes/t_wilcoxon_signed_rank_matched_pairs.md",
      checker = "validation/independent/g10/check_wilcoxon.R"),
    t_wilcoxon_mann_whitney = list(
      id = id, kind = "rank_sum", effect = "d",
      label = list(en = "Means: Wilcoxon-Mann-Whitney test (two groups)",
                   de = "Mittelwerte: Wilcoxon-Mann-Whitney-Test (zwei Gruppen)"),
      effect_label = list(en = "Effect size d", de = "Effektst\u00e4rke d"),
      effect_description = paste("Standardised shift d = (mu1 - mu2) / sigma between two distributions of the",
                                 "same shape (shift model; manual ch. 24)."),
      gpower_path = "Tests > Means > Two independent groups: Wilcoxon (non-parametric)",
      gpower_label = "t tests - Means: Wilcoxon-Mann-Whitney test (two groups)",
      manual_pages = "ch. 24 (PDF pp. 59-62)",
      notes = "specs/method-notes/t_wilcoxon_mann_whitney.md",
      checker = "validation/independent/g10/check_wilcoxon.R"),
    stop("unknown Wilcoxon method ", id, call. = FALSE))
}

# Applicability conditions (registry grammar; an option that does not apply
# is absent, so a condition on it is false).
.ig_wilcoxon_conditions <- function(effect) {
  list(
    effect = list(any = list(list(option = "method", "in" = "are"),
                             list(option = "effect_size_parameter", "in" = "d"))),
    p1 = list(any = list(list(option = "effect_size_parameter", "in" = "p1"),
                         list(option = "lehmann_input_method", "in" = "manual_moments"))),
    relative_sd = list(option = "effect_size_parameter", "in" = "p1"),
    moments = list(option = "lehmann_input_method", "in" = "manual_moments"),
    are = list(option = "are_specification", "in" = "manual_are")
  )
}

.ig_method_wilcoxon <- function(id) {
  info <- .ig_wilcoxon_info(id)
  cond <- .ig_wilcoxon_conditions(info$effect)
  two_groups <- identical(info$kind, "rank_sum")
  eff <- info$effect
  sample_params <- if (two_groups) {
    list(
      ig_param("n_total", role = "design", domain = ig_domain_integer(lower = 4), unit = "subjects",
               label = list(en = "Total sample size", de = "Gesamtstichprobenumfang"),
               description = "Solved in a priori mode (n1 + n2; even, split by the allocation ratio).",
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
                                   "value 1). Admissible designs: even N, n1 = floor(N/(1 + r) + 1/2),",
                                   "n2 = N - n1, both at least 1."),
               aliases = "allocation_ratio_n2_n1",
               modes = "a_priori")
    )
  } else {
    list(ig_param("n_total", role = "design", domain = ig_domain_integer(lower = 3), unit = "subjects",
                  label = list(en = "Total sample size", de = "Gesamtstichprobenumfang"),
                  description = paste("Number of observations (or pairs). G*Power admits N >= 3",
                                      "(live probe: N = 2 is refused, N = 3 is computed)."),
                  aliases = "total_sample_size"))
  }
  ig_method(
    id = id,
    family = "t",
    label = info$label,
    params = ig_params(
      ig_param(eff, role = "effect", domain = ig_domain_real(), label = info$effect_label,
               description = paste(info$effect_description,
                                   "Used by the A.R.E. method and by the Lehmann method with",
                                   "effect_size_parameter = 'd' (the Lehmann method uses |d|)."),
               aliases = paste0("effect_size_", eff), applies_when = cond$effect),
      ig_param("p1", role = "effect",
               domain = ig_domain_real(0, 1, lower_inclusive = FALSE, upper_inclusive = FALSE),
               label = list(en = "Moment p1", de = "Moment p1"),
               description = paste(if (two_groups) "p1 = P(X < Y) (X group 1, Y group 2)" else "p1 = P(X > 0)",
                                   "under the shift alternative. Lehmann method only: the effect-size",
                                   "parameter when effect_size_parameter = 'p1' (the shift with this p1",
                                   "under the parent, divided by relative_sd, is d), or the first of the",
                                   "manually specified moments."),
               aliases = "moment_p1", applies_when = cond$p1),
      ig_param("relative_sd", role = "nuisance",
               domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE),
               default = 1, label = list(en = "Relative SD", de = "Relative Standardabweichung"),
               description = paste("G*Power's 'Relative SD' of the p1 effect-size parameter: d is the shift",
                                   "that gives p1 under the standardised parent, divided by this value",
                                   "(reproduces the manual example 24.3 Ex. 2 and every capture). The",
                                   "default 1 is G*Power's initial field value."),
               applies_when = cond$relative_sd),
      ig_param("p2", role = "nuisance",
               domain = ig_domain_real(0, 1, lower_inclusive = FALSE, upper_inclusive = FALSE),
               label = list(en = "Moment p2", de = "Moment p2"),
               description = if (two_groups) "p2 = P(X < Y, X < Y') (manually specified moments)." else
                 "p2 = P(X1 + X2 > 0) (manually specified moments).",
               aliases = "moment_p2", applies_when = cond$moments),
      ig_param("p3", role = "nuisance",
               domain = ig_domain_real(0, 1, lower_inclusive = FALSE, upper_inclusive = FALSE),
               label = list(en = "Moment p3", de = "Moment p3"),
               description = if (two_groups) "p3 = P(X < Y, X' < Y) (manually specified moments)." else
                 "p3 = P(X1 + X2 > 0, X1 + X3 > 0) (manually specified moments).",
               aliases = "moment_p3", applies_when = cond$moments),
      ig_param("are", role = "nuisance",
               domain = ig_domain_real(.ig_wilcoxon_min_are, Inf, upper_inclusive = FALSE),
               label = list(en = "A.R.E.", de = "A.R.E."),
               description = paste("Asymptotic relative efficiency of the Wilcoxon test against the t test,",
                                   "specified manually (A.R.E. method). Values below 0.864, the minimum",
                                   "over continuous symmetric parents, are refused (as in G*Power)."),
               aliases = "are_value", applies_when = cond$are),
      ig_param_alpha(),
      ig_param_power(),
      ig_param_error_ratio(),
      sample_params
    ),
    options = ig_options(
      ig_option("tails", values = c("one", "two"), default = "one",
                labels = c(one = "One", two = "Two"), label = list(en = "Tail(s)", de = "Seiten")),
      ig_option("method", values = c("are", "lehmann"), default = "are",
                labels = c(are = "Asymptotic Relative Efficiency (based on t test)",
                           lehmann = "Asymptotic Power (Lehmann method based on z test)"),
                label = list(en = "Power method", de = "Power-Methode")),
      ig_option("are_specification", values = c("predefined_parent", "manual_are"),
                default = "predefined_parent",
                labels = c(predefined_parent = "Use predefined parent distributions",
                           manual_are = "Specify A.R.E. manually"),
                applies_when = list(option = "method", "in" = "are"),
                label = list(en = "A.R.E. specification", de = "A.R.E.-Angabe")),
      ig_option("lehmann_input_method", values = c("predefined_parent", "manual_moments"),
                default = "predefined_parent",
                labels = c(predefined_parent = "Use predefined parent distributions",
                           manual_moments = "Specify moments p1,p2,p3 manually"),
                applies_when = list(option = "method", "in" = "lehmann"),
                label = list(en = "Input method", de = "Eingabemethode")),
      ig_option("effect_size_parameter", values = c("d", "p1"), default = "d",
                labels = c(d = "d = \u0394/\u03c3", p1 = "p1 = Prob(X < Y)"),
                applies_when = list(option = "lehmann_input_method", "in" = "predefined_parent"),
                label = list(en = "Effect size parameter", de = "Effektst\u00e4rkeparameter")),
      ig_option("continuity_correction", values = c("on", "off"), default = "off",
                labels = c(on = "checked", off = "unchecked"),
                applies_when = list(option = "method", "in" = "lehmann"),
                label = list(en = "Use continuity correction", de = "Stetigkeitskorrektur")),
      ig_option("distribution_plot", values = c("normalized_z", "normal_rank_sum"), default = "normalized_z",
                labels = c(normalized_z = "Normalized z-Distribution",
                           normal_rank_sum = "Normal approximation of rank sum"),
                applies_when = list(option = "method", "in" = "lehmann"),
                label = list(en = "Distribution plots", de = "Verteilungsgrafik")),
      ig_option("parent_distribution", values = c("normal", "laplace", "logistic", "min_are"),
                default = "normal",
                labels = c(normal = "Normal", laplace = "Laplace", logistic = "Logistic",
                           min_are = "min ARE (A.R.E. method only)"),
                applies_when = list(any = list(list(option = "are_specification", "in" = "predefined_parent"),
                                               list(option = "lehmann_input_method", "in" = "predefined_parent"))),
                label = list(en = "Parent distribution", de = "Ausgangsverteilung"))
    ),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = function(p) .ig_power_wilcoxon(p, info),
    design = if (two_groups) .ig_design_wilcoxon_mann_whitney else NULL,
    discrete = FALSE,
    effect_direction = "increasing",
    engine = list(package = "stats",
                  fun = paste("A.R.E.: noncentral t via .ig_t_power_g1 / .ig_power_t_test (pt/qt, inferguidR integral route);",
                              "Lehmann: pnorm/qnorm with moments by stats::integrate"),
                  version = NA_character_),
    checker = info$checker,
    reference = list(
      gpower_path = info$gpower_path,
      gpower_label = info$gpower_label,
      manual_pages = info$manual_pages,
      allocation_rule = if (two_groups) paste(
        "even N; n1 = floor(N/(1 + r) + 1/2), n2 = N - n1 (r = allocation ratio N2/N1; G*Power's split);",
        "the smallest such design whose power reaches the target (G*Power selects N by the power of the",
        "fractional group sizes, which can give a different N; method note)") else NULL,
      method_notes = info$notes
    ),
    status = "independently_checked",
    targets = list(sensitivity = c(eff, "p1")),
    constraints = function(p, mode, target) .ig_wilcoxon_constraints(p, info, mode, target)
  )
}

# Cross-input checks at specification time ----------------------------------------------

.ig_wilcoxon_constraints <- function(p, info, mode, target) {
  out <- character()
  eff <- info$effect
  if (identical(info$kind, "rank_sum") && !is.null(p$n1) && !is.null(p$n2) &&
      as.double(p$n1) + as.double(p$n2) < 3) {
    out <- c(out, "n1 + n2 must be at least 3 (G*Power's admissible total sample size)")
  }
  if (!is.null(p$alpha) && !is.null(p$power) && p$power < p$alpha) {
    out <- c(out, "power must be greater than or equal to alpha (as in G*Power)")
  }
  lehmann <- identical(p$method, "lehmann")
  if (lehmann && identical(p$parent_distribution, "min_are")) {
    out <- c(out, "the parent distribution 'min_are' is available for the A.R.E. method only")
  }
  if (identical(mode, "sensitivity")) {
    if (lehmann && identical(p$lehmann_input_method, "manual_moments")) {
      out <- c(out, "a sensitivity analysis is not possible if the moments are specified manually (as in G*Power)")
    } else if (lehmann && identical(p$effect_size_parameter, "p1")) {
      if (!identical(target, "p1")) out <- c(out, "with effect_size_parameter = 'p1' the sensitivity target is p1; use target = \"p1\"")
    } else if (!identical(target, eff)) {
      out <- c(out, sprintf("the sensitivity target p1 applies only with method = 'lehmann' and effect_size_parameter = 'p1'; use target = \"%s\"", eff))
    }
  }
  if (lehmann && identical(p$lehmann_input_method, "manual_moments") &&
      !is.null(p$p1) && !is.null(p$p2) && !is.null(p$p3)) {
    out <- c(out, .ig_wilcoxon_moment_problems(info$kind, p$p1, p$p2, p$p3))
  }
  out
}

# Necessary conditions for manually given moments to belong to a
# distribution. Rank-sum: {X < Y, X < Y'} and {X < Y, X' < Y} are sub-events
# of {X < Y}, so p2, p3 <= p1; and p2 = E[P(Y > X | X)^2] >= p1^2, likewise
# p3 >= p1^2 (Jensen). Signed-rank: p3 <= p2 and p3 = E[P(X1 + X2 > 0 |
# X1)^2] >= p2^2. A set violating them is refused (G*Power computes a
# number from such inputs; its own initial values p1 = 0.5, p2 = p3 = 0.65
# of the rank-sum test violate p2, p3 <= p1).
.ig_wilcoxon_moment_problems <- function(kind, p1, p2, p3) {
  tol <- 1e-12
  bad <- if (identical(kind, "rank_sum")) {
    c("p2 <= p1" = p2 > p1 + tol, "p3 <= p1" = p3 > p1 + tol,
      "p2 >= p1^2" = p2 < p1^2 - tol, "p3 >= p1^2" = p3 < p1^2 - tol)
  } else {
    c("p3 <= p2" = p3 > p2 + tol, "p3 >= p2^2" = p3 < p2^2 - tol)
  }
  if (!any(bad)) return(character())
  sprintf(paste("the manually specified moments violate %s, which every distribution satisfies (%s), so they",
                "cannot describe a shift alternative"), paste(names(bad)[bad], collapse = " and "),
          if (identical(kind, "rank_sum")) "p1 = P(X < Y), p2 = P(X < Y, X < Y'), p3 = P(X < Y, X' < Y)" else
            "p2 = P(X1 + X2 > 0), p3 = P(X1 + X2 > 0, X1 + X3 > 0)")
}

# Allocation rule (rank-sum test) ---------------------------------------------------------

.ig_design_wilcoxon_mann_whitney <- function(p, n_total, continuous = FALSE) {
  r <- p$allocation_ratio
  if (continuous) return(c(n1 = n_total / (1 + r), n2 = r * n_total / (1 + r)))
  if (n_total %% 2 != 0) return(NULL)
  n1 <- floor(n_total / (1 + r) + 0.5)
  n2 <- n_total - n1
  if (n1 < 1 || n2 < 1) return(NULL)
  c(n1 = n1, n2 = n2)
}

# Parent distributions (standardised to variance 1) ---------------------------------------

.ig_wilcoxon_parent <- function(parent) {
  switch(parent,
    normal = list(
      density = function(x) stats::dnorm(x),
      cdf = function(x, lower_tail = TRUE) stats::pnorm(x, lower.tail = lower_tail),
      quantile = function(p) stats::qnorm(p),
      kinks = FALSE),
    laplace = {
      b <- 1 / sqrt(2)
      list(
        density = function(x) exp(-abs(x) / b) / (2 * b),
        cdf = function(x, lower_tail = TRUE) {
          y <- if (lower_tail) x else -x
          ifelse(y < 0, 0.5 * exp(y / b), 1 - 0.5 * exp(-y / b))
        },
        quantile = function(p) ifelse(p < 0.5, b * log(2 * p), -b * log(2 * (1 - p))),
        kinks = TRUE)
    },
    logistic = {
      s <- sqrt(3) / pi
      list(
        density = function(x) stats::dlogis(x, 0, s),
        cdf = function(x, lower_tail = TRUE) stats::plogis(x, 0, s, lower.tail = lower_tail),
        quantile = function(p) stats::qlogis(p, 0, s),
        kinks = FALSE)
    },
    .ig_abort("invalid_input", sprintf("parent distribution '%s' has no moments (Lehmann method)", parent)))
}

# G_k(t) = integral of f(u) F(u + t)^k du, k = 1, 2, for the standardised
# parent (f density, F CDF). G_1 is the CDF of the difference of two
# independent parent variables. For t >= 0 the complement
# 1 - G_k(t) = integral of f(u) S(u + t) (1 + [k = 2] F(u + t)) du is
# integrated (S = 1 - F), so values near one keep their accuracy; the
# Laplace parent is split at its kinks (u = 0 and u = -t).
.ig_wilcoxon_gk <- function(t, k, parent) {
  P <- .ig_wilcoxon_parent(parent)
  if (t >= 0) {
    g <- if (k == 1) {
      function(u) P$density(u) * P$cdf(u + t, lower_tail = FALSE)
    } else {
      function(u) {
        s <- P$cdf(u + t, lower_tail = FALSE)
        P$density(u) * s * (2 - s)
      }
    }
  } else {
    g <- if (k == 1) {
      function(u) P$density(u) * P$cdf(u + t)
    } else {
      function(u) P$density(u) * P$cdf(u + t)^2
    }
  }
  cuts <- sort(unique(c(-Inf, 0, -t, Inf)))
  total <- 0
  for (j in seq_len(length(cuts) - 1L)) {
    total <- total + .ig_wilcoxon_integrate(g, cuts[j], cuts[j + 1L])
  }
  if (t >= 0) 1 - total else total
}

.ig_wilcoxon_integrate <- function(g, lower, upper) {
  if (lower == upper) return(0)
  r <- tryCatch(stats::integrate(g, lower, upper, rel.tol = 1e-12, abs.tol = 1e-15, subdivisions = 2000L),
                error = function(e) NULL)
  if (is.null(r)) {
    r <- tryCatch(stats::integrate(g, lower, upper, rel.tol = 1e-10, abs.tol = 1e-14, subdivisions = 2000L),
                  error = function(e) {
                    .ig_abort("numerical_failure", paste("the moment integral of the parent distribution",
                                                         "could not be evaluated:", conditionMessage(e)))
                  })
  }
  r$value
}

# Moments p1, p2, p3 of the shift alternative with shift `delta` (in units
# of the parent SD), cached by (kind, parent, delta).
#   signed-rank: p1 = F(delta), p2 = G_1(2 delta), p3 = G_2(2 delta)
#   rank-sum:    p1 = G_1(delta), p2 = p3 = G_2(delta) (symmetric parent)
.ig_wilcoxon_cache <- new.env(parent = emptyenv())

.ig_wilcoxon_moments <- function(kind, parent, delta) {
  key <- sprintf("%s|%s|%a", kind, parent, delta)
  hit <- .ig_wilcoxon_cache[[key]]
  if (!is.null(hit)) return(hit)
  if (length(ls(.ig_wilcoxon_cache)) > 20000L) rm(list = ls(.ig_wilcoxon_cache), envir = .ig_wilcoxon_cache)
  P <- .ig_wilcoxon_parent(parent)
  m <- if (identical(kind, "signed_rank")) {
    c(p1 = P$cdf(delta), p2 = .ig_wilcoxon_gk(2 * delta, 1, parent), p3 = .ig_wilcoxon_gk(2 * delta, 2, parent))
  } else {
    p23 <- .ig_wilcoxon_gk(delta, 2, parent)
    c(p1 = .ig_wilcoxon_gk(delta, 1, parent), p2 = p23, p3 = p23)
  }
  assign(key, m, envir = .ig_wilcoxon_cache)
  m
}

# Shift whose p1 under the parent equals `p1`: the parent quantile for the
# signed-rank tests (p1 = F(delta)); for the rank-sum test sqrt(2) qnorm(p1)
# for the normal parent and the root of G_1(delta) = p1 otherwise.
.ig_wilcoxon_shift_from_p1 <- function(kind, parent, p1) {
  P <- .ig_wilcoxon_parent(parent)
  if (identical(kind, "signed_rank")) return(P$quantile(p1))
  if (identical(parent, "normal")) return(sqrt(2) * stats::qnorm(p1))
  f <- function(x) .ig_wilcoxon_gk(x, 1, parent) - p1
  lo <- -1
  hi <- 1
  while (f(lo) > 0 && lo > -1e3) lo <- 2 * lo
  while (f(hi) < 0 && hi < 1e3) hi <- 2 * hi
  stats::uniroot(f, c(lo, hi), tol = 1e-15, maxiter = 1000L)$root
}

# Forward model ----------------------------------------------------------------------------

.ig_power_wilcoxon <- function(p, info) {
  kind <- info$kind
  eff <- info$effect
  if (identical(kind, "rank_sum")) {
    n1 <- as.double(p$n1)
    n2 <- as.double(p$n2)
    n <- n1 + n2
  } else {
    n <- as.double(p$n_total)
  }
  tail_note <- if (identical(p$tails, "two")) {
    "two-sided test: alpha/2 in each tail; the power counts both rejection regions (G*Power convention, manual 25.5)"
  } else {
    "one-sided test"
  }
  if (identical(p$method, "are")) {
    manual <- identical(p$are_specification, "manual_are")
    k <- if (manual) as.double(p$are) else .ig_wilcoxon_are[[p$parent_distribution]]
    d <- p[[eff]]
    if (is.null(d)) .ig_abort("invalid_input", sprintf("the A.R.E. method needs the effect size %s", eff))
    if (identical(kind, "rank_sum")) {
      df <- n * k - 2
      ncp <- d * sqrt(n1 * n2 * k / n)
    } else {
      df <- n * k - 1
      ncp <- d * sqrt(n * k)
    }
    if (!(df > 0)) {
      .ig_abort("invalid_input", sprintf("the effective degrees of freedom N k - %d = %s are not positive",
                                         if (identical(kind, "rank_sum")) 2L else 1L, .ig_format_double(df)))
    }
    # .ig_t_power_g1() (R/methods-t.R) is the t family's wrapper of
    # .ig_power_t_test(): probabilities at |ncp| (exact mirror identity) and
    # the larger of power and beta as the complement of the smaller one, so a
    # power near one never exceeds one through nmath's absolute error (at
    # df near 1e5 the plain helper returned 1 + 4e-11 and the sensitivity scan
    # stopped; review 2026-10-06).
    core <- .ig_t_power_g1(df, ncp, p$alpha, p$tails)
    return(list(
      power = core$power,
      beta = core$beta,
      alpha_attained = p$alpha,
      critical = core$critical,
      df = df,
      ncp = ncp,
      n_total = n,
      group_sizes = if (identical(kind, "rank_sum")) c(n1 = n1, n2 = n2) else NULL,
      outputs = list(are = k, effective_n = n * k, power_method = "are"),
      statistic_distribution = list(
        h0 = list(family = "t", df = df, ncp = 0),
        h1 = list(family = "t", df = df, ncp = ncp),
        rejection_side = core$side
      ),
      assumptions = c(
        sprintf("A.R.E. method: power of the t test at the effective sample size N k with k = %s (%s)",
                format(signif(k, 10)), if (manual) "specified manually" else
                  sprintf("%s parent", p$parent_distribution)),
        if (identical(kind, "rank_sum")) "df = k (n1 + n2) - 2; ncp = d sqrt(n1 n2 k / (n1 + n2))" else
          "df = N k - 1; ncp = d sqrt(N k)",
        "large-sample approximation: for small N it can overstate the power of the exact Wilcoxon test",
        if (identical(p$tails, "two")) tail_note else
          paste("one-sided test: the rejection region lies on the side of the sign of d (G*Power manual 25.5);",
                "d = 0 uses the upper tail")
      )
    ))
  }

  # Lehmann method.
  cc <- identical(p$continuity_correction, "on")
  a <- if (cc) 0.5 else 0
  manual <- identical(p$lehmann_input_method, "manual_moments")
  d_used <- NA_real_
  shift <- NA_real_
  if (manual) {
    m <- c(p1 = as.double(p$p1), p2 = as.double(p$p2), p3 = as.double(p$p3))
  } else {
    parent <- p$parent_distribution
    if (identical(p$effect_size_parameter, "p1")) {
      shift <- .ig_wilcoxon_shift_from_p1(kind, parent, p$p1)
      d_used <- shift / p$relative_sd
    } else {
      d_used <- p[[eff]]
      if (is.null(d_used)) .ig_abort("invalid_input", sprintf("the Lehmann method needs the effect size %s", eff))
    }
    m <- .ig_wilcoxon_moments(kind, parent, abs(d_used))
  }
  p1 <- m[["p1"]]
  p2 <- m[["p2"]]
  p3 <- m[["p3"]]
  if (identical(kind, "rank_sum")) {
    e0 <- n1 * n2 / 2
    v0 <- n1 * n2 * (n + 1) / 12
    e1 <- n1 * n2 * p1
    v1 <- n1 * n2 * (p1 * (1 - p1) + (n2 - 1) * (p2 - p1^2) + (n1 - 1) * (p3 - p1^2))
  } else {
    e0 <- n * (n + 1) / 4
    v0 <- n * (n + 1) * (2 * n + 1) / 24
    e1 <- n * (n - 1) * p2 / 2 + n * p1
    v1 <- n * (n - 1) * (n - 2) * (p3 - p2^2) + n * (n - 1) * (2 * (p1 - p2)^2 + 3 * p2 * (1 - p2)) / 2 +
      n * p1 * (1 - p1)
  }
  degenerate <- FALSE
  if (!(v1 > 0)) {
    if (manual) {
      .ig_abort("invalid_input", sprintf(paste(
        "the manually specified moments give a non-positive variance %s of the rank statistic under H1",
        "at this sample size; they cannot belong to a distribution"), .ig_format_double(v1)))
    }
    # Predefined parent with a very large shift: the moments equal one to
    # double precision and the H1 distribution of the statistic is
    # degenerate at its maximum.
    degenerate <- TRUE
    v1 <- 0
  }
  s0 <- sqrt(v0)
  s1 <- sqrt(v1)
  z <- .ig_q_norm(if (identical(p$tails, "two")) p$alpha / 2 else p$alpha, lower_tail = FALSE)
  c_up <- e0 + z * s0 + a
  c_lo <- e0 - z * s0 - a
  ptail <- function(q, lt) {
    if (degenerate) return(if (lt) as.double(e1 <= q) else as.double(e1 > q))
    .ig_p_norm(q, e1, s1, lower_tail = lt)
  }
  if (identical(p$tails, "two")) {
    region <- .ig_two_sided_prob(ptail, c_lo, c_up)
    power <- region$reject
    beta <- region$accept
    crit_z <- c(-z, z)
    crit_t <- c(c_lo, c_up)
  } else {
    power <- ptail(c_up, FALSE)
    beta <- ptail(c_up, TRUE)
    crit_z <- z
    crit_t <- c_up
  }
  plot_t <- identical(p$distribution_plot, "normal_rank_sum")
  stat_name <- if (identical(kind, "rank_sum")) "Mann-Whitney U" else "signed-rank V"
  outputs <- list(
    power_method = "lehmann",
    moment_p1 = p1, moment_p2 = p2, moment_p3 = p3,
    critical_z = crit_z,
    critical_rank_sum = crit_t,
    critical_rank_sum_lower = if (identical(p$tails, "two")) c_lo else c_up,
    critical_rank_sum_upper = c_up,
    mean_h0 = e0, sd_h0 = s0, mean_h1 = e1, sd_h1 = s1
  )
  if (!manual) {
    outputs$effect_size_used <- abs(d_used)
    outputs$effect_size_d <- d_used
    if (!is.na(shift)) outputs$shift_from_p1 <- shift
  }
  list(
    power = power,
    beta = beta,
    alpha_attained = p$alpha,
    critical = if (plot_t) crit_t else crit_z,
    n_total = n,
    group_sizes = if (identical(kind, "rank_sum")) c(n1 = n1, n2 = n2) else NULL,
    outputs = outputs,
    statistic_distribution = if (plot_t) list(
      scale = "rank statistic",
      h0 = list(family = "normal", mean = e0, sd = s0),
      h1 = list(family = "normal", mean = e1, sd = s1),
      rejection_side = if (identical(p$tails, "two")) "both" else "upper",
      continuity_correction = a
    ) else list(
      scale = "normalized z",
      h0 = list(family = "normal", mean = 0, sd = 1),
      h1 = list(family = "normal", mean = (e1 - e0) / s0, sd = s1 / s0),
      rejection_side = if (identical(p$tails, "two")) "both" else "upper",
      continuity_shift = a / s0
    ),
    assumptions = c(
      sprintf(paste("Lehmann method: normal approximation of the %s statistic with mean %s and variance %s",
                    "under H1 from the moments p1, p2, p3%s"), stat_name,
              if (identical(kind, "rank_sum")) "n1 n2 p1" else "N(N - 1) p2 / 2 + N p1",
              "(Lehmann 1975; corrected manual formulas, specs/engine-decisions.md G3)",
              if (manual) " (specified manually)" else sprintf(" of the %s parent", p$parent_distribution)),
      if (!manual) paste("the moments are those of the shift |d| (G*Power uses the magnitude of d with a",
                         "predefined parent; for a symmetric parent both rejection regions give the same power)"),
      if (!manual && identical(p$effect_size_parameter, "p1"))
        "d = shift with the given p1 under the parent / relative SD (G*Power's 'Relative SD')",
      if (cc) "continuity correction: the critical value of the rank statistic is moved by 0.5 away from its null mean (conservative)",
      "the significance level is nominal; the size of the exact rank test differs",
      if (degenerate) "the moments equal one to double precision; the statistic is degenerate at its maximum under H1",
      tail_note
    )
  )
}
