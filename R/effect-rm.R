# Effect-size drawers: repeated-measures ANOVA and MANOVA (group G3) ---------------------
#
# The "Determine =>" dialogs of the six repeated-measures procedures
# (certified live inventory, effect_size_drawer_live). Reached through
# effect_size() (R/core-api.R), which looks up `.ig_drawer_<method_id>()`.
#
# Input modes (live drawer controls in brackets):
#
#   direct          partial eta^2 -> f = sqrt(eta^2 / (1 - eta^2))   [group 'Direct']
#   from_variances  variance explained by the effect, variance within ->
#                   f = sqrt(var_effect / var_within), partial eta^2 =
#                   var_effect / (var_effect + var_within)            [group 'From variances']
#   from_means      group means, common SD within groups, optional group
#                   sizes -> f = sigma_m / sigma with sigma_m^2 =
#                   sum_i w_i (mu_i - mu_w)^2, w_i = n_i / N           [combo 'Effect size from means';
#                                                                      between factors only]
#   convert         the same effect between the effect-size conventions
#                   (not a G*Power drawer mode: the conversion implied by
#                   equal noncentralities, see the method notes)
#   from_criterion  MANOVA criterion (Pillai V, Wilks U, Hotelling-Lawley
#                   T1/T2) -> f of the approximation                  [combo 'Effectsize from criterion';
#                                                                      MANOVA within-between]
#   from_means_cov  group x measurement means and the covariance matrix ->
#                   population criterion and f                        [combo 'Effectsize from mean and
#                                                                      cov'; MANOVA within-between]
#
# 'direct' and 'from_variances' are compared with emulator captures
# (validation/compare/<method_id>/). 'from_means' follows the ANOVA
# definition of f (manual ch. 10.1); 'from_criterion' and 'from_means_cov'
# use the population definitions of the MANOVA effect sizes; the drawer
# grids of G*Power could not be driven by the capture harness, so these
# modes carry no G*Power comparison (method notes).

.ig_rm_drawer_modes <- list(
  f_rm_anova_between = c("from_means", "from_variances", "direct", "convert"),
  f_rm_anova_within = c("from_variances", "direct", "convert"),
  f_rm_anova_within_between = c("from_variances", "direct", "convert"),
  f_rm_manova_between = c("from_means", "from_variances", "direct", "convert"),
  f_rm_manova_within = c("from_variances", "direct", "convert"),
  f_rm_manova_within_between = c("from_criterion", "from_means_cov")
)

.ig_drawer_f_rm_anova_between <- function(input_mode, inputs) .ig_rm_drawer("f_rm_anova_between", input_mode, inputs)
.ig_drawer_f_rm_anova_within <- function(input_mode, inputs) .ig_rm_drawer("f_rm_anova_within", input_mode, inputs)
.ig_drawer_f_rm_anova_within_between <- function(input_mode, inputs) {
  .ig_rm_drawer("f_rm_anova_within_between", input_mode, inputs)
}
.ig_drawer_f_rm_manova_between <- function(input_mode, inputs) .ig_rm_drawer("f_rm_manova_between", input_mode, inputs)
.ig_drawer_f_rm_manova_within <- function(input_mode, inputs) .ig_rm_drawer("f_rm_manova_within", input_mode, inputs)
.ig_drawer_f_rm_manova_within_between <- function(input_mode, inputs) {
  .ig_rm_drawer("f_rm_manova_within_between", input_mode, inputs)
}

.ig_rm_drawer <- function(method_id, input_mode, inputs) {
  modes <- .ig_rm_drawer_modes[[method_id]]
  if (!.ig_is_string(input_mode) || !(input_mode %in% modes)) {
    .ig_abort_invalid(sprintf("`input_mode` must be one of: %s", paste(modes, collapse = ", ")),
                      "invalid effect-size input")
  }
  approach <- if (startsWith(method_id, "f_rm_manova")) "manova" else "anova"
  effect <- sub("^f_rm_(an|man)ova_", "", method_id)
  switch(input_mode,
         direct = .ig_rm_drawer_direct(inputs),
         from_variances = .ig_rm_drawer_from_variances(inputs),
         from_means = .ig_rm_drawer_from_means(inputs),
         convert = .ig_rm_drawer_convert(inputs, approach, effect),
         from_criterion = .ig_rm_drawer_from_criterion(inputs),
         from_means_cov = .ig_rm_drawer_from_means_cov(inputs))
}

# Input checks -----------------------------------------------------------------------

# Validate a named list of drawer inputs: `required` and `optional` names,
# numeric scalars unless listed in `vectors`, `matrices` or `codes`.
.ig_rm_drawer_inputs <- function(inputs, required, optional = character(), vectors = character(),
                                 matrices = character(), codes = list()) {
  problems <- character()
  unknown <- setdiff(names(inputs), c(required, optional))
  if (length(unknown)) {
    problems <- c(problems, sprintf("input '%s' is not used in this input mode (inputs: %s)", unknown,
                                    paste(c(required, optional), collapse = ", ")))
  }
  vals <- list()
  for (nm in c(required, intersect(optional, names(inputs)))) {
    v <- inputs[[nm]]
    if (is.null(v)) {
      problems <- c(problems, sprintf("missing input '%s'", nm))
      next
    }
    if (nm %in% names(codes)) {
      if (!.ig_is_string(v) || !(v %in% codes[[nm]])) {
        problems <- c(problems, sprintf("'%s' must be one of: %s", nm, paste(codes[[nm]], collapse = ", ")))
      } else {
        vals[[nm]] <- v
      }
      next
    }
    ok <- is.numeric(v) && length(v) >= 1L && all(is.finite(v))
    if (nm %in% matrices) {
      ok <- ok && is.matrix(v)
    } else if (!(nm %in% vectors)) {
      ok <- ok && length(v) == 1L
    }
    if (!ok) {
      problems <- c(problems, sprintf("'%s' must be %s", nm,
                                      if (nm %in% matrices) "a finite numeric matrix" else if (nm %in% vectors)
                                        "a finite numeric vector" else "a single finite number"))
      next
    }
    vals[[nm]] <- if (is.matrix(v)) v else as.double(v)
  }
  list(values = vals, problems = problems)
}

.ig_rm_drawer_result <- function(f, formula, source, inputs, warnings = character(), extra = list()) {
  c(list(effect = "f", value = f, magnitude = abs(f), formula = formula, source = source,
         inputs = inputs, warnings = warnings), extra)
}

# Modes ------------------------------------------------------------------------------

.ig_rm_drawer_direct <- function(inputs) {
  chk <- .ig_rm_drawer_inputs(inputs, "partial_eta2")
  e <- chk$values$partial_eta2
  if (!is.null(e) && (e <= 0 || e >= 1)) chk$problems <- c(chk$problems, "'partial_eta2' must lie in (0, 1)")
  if (length(chk$problems)) .ig_abort_invalid(chk$problems, "invalid effect-size input")
  .ig_rm_drawer_result(sqrt(e / (1 - e)), "f = sqrt(partial_eta2 / (1 - partial_eta2))",
                       "G*Power drawer group 'Direct' (partial eta^2 to f)", chk$values,
                       extra = list(partial_eta2 = e))
}

.ig_rm_drawer_from_variances <- function(inputs) {
  chk <- .ig_rm_drawer_inputs(inputs, c("variance_effect", "variance_within"))
  v <- chk$values
  if (!is.null(v$variance_effect) && v$variance_effect < 0) {
    chk$problems <- c(chk$problems, "'variance_effect' must not be negative")
  }
  if (!is.null(v$variance_within) && v$variance_within <= 0) {
    chk$problems <- c(chk$problems, "'variance_within' must be positive")
  }
  if (length(chk$problems)) .ig_abort_invalid(chk$problems, "invalid effect-size input")
  .ig_rm_drawer_result(sqrt(v$variance_effect / v$variance_within),
                       "f = sqrt(variance_effect / variance_within); partial eta^2 = variance_effect / (variance_effect + variance_within)",
                       "G*Power drawer group 'From variances'", v,
                       extra = list(partial_eta2 = v$variance_effect / (v$variance_effect + v$variance_within)))
}

.ig_rm_drawer_from_means <- function(inputs) {
  chk <- .ig_rm_drawer_inputs(inputs, c("means", "sd"), optional = "n", vectors = c("means", "n"))
  v <- chk$values
  if (!is.null(v$means) && length(v$means) < 2L) chk$problems <- c(chk$problems, "'means' needs at least 2 groups")
  if (!is.null(v$sd) && v$sd <= 0) chk$problems <- c(chk$problems, "'sd' must be positive")
  if (!is.null(v[["n"]]) && !is.null(v$means)) {
    if (length(v[["n"]]) != length(v$means)) chk$problems <- c(chk$problems, "'n' must have one size per group")
    else if (any(v[["n"]] <= 0)) chk$problems <- c(chk$problems, "'n' must be positive")
  }
  if (length(chk$problems)) .ig_abort_invalid(chk$problems, "invalid effect-size input")
  w <- if (is.null(v[["n"]])) rep(1 / length(v$means), length(v$means)) else v[["n"]] / sum(v[["n"]])
  mu <- sum(w * v$means)
  sigma_m <- sqrt(sum(w * (v$means - mu)^2))
  warn <- character()
  if (!is.null(v[["n"]]) && any(v[["n"]] != v[["n"]][1])) {
    warn <- paste("unequal group sizes weight the means here, while the a priori design of this test uses",
                  "equal groups (N a multiple of k)")
  }
  .ig_rm_drawer_result(sigma_m / v$sd, "f = sigma_m / sd, sigma_m^2 = sum_i w_i (mean_i - weighted mean)^2, w_i = n_i / N",
                       "G*Power drawer 'Effect size from means' (ANOVA definition of f, manual ch. 10.1); the f of the 'as in GPower 3.0' convention",
                       v, warn, extra = list(sigma_m = sigma_m))
}

# f in one convention -> f in another, through equal noncentralities
# (lambda / N per subject). RM-ANOVA conventions are the codes of
# 'effect_size_specification'; RM-MANOVA between/within use 'on' (f with the
# correlation) and 'off' (f(V)).
.ig_rm_drawer_convert <- function(inputs, approach, effect) {
  codes <- if (approach == "anova") names(.ig_rm_spec_labels) else c("on", "off")
  chk <- .ig_rm_drawer_inputs(inputs, c("f", "from", "to"), optional = c("measurements", "rho", "n_total", "groups"),
                              codes = list(from = codes, to = codes))
  v <- chk$values
  v[["m"]] <- v$measurements
  v[["k"]] <- v$groups
  need <- character()
  uses <- c(v$from, v$to)
  if (any(uses %in% c("gpower30", "on"))) need <- c(need, "measurements", "rho")
  if (any(uses == "spss")) need <- c(need, "measurements", "n_total", "groups")
  miss <- setdiff(unique(need), names(v))
  if (length(miss)) chk$problems <- c(chk$problems, sprintf("missing input '%s' for this conversion", miss))
  if (!is.null(v$f) && v$f < 0) chk$problems <- c(chk$problems, "'f' must not be negative")
  if (!is.null(v$rho) && (v$rho <= -1 || v$rho >= 1)) chk$problems <- c(chk$problems, "'rho' must lie in (-1, 1)")
  if (!is.null(v[["m"]]) && (v[["m"]] < 2 || v[["m"]] != round(v[["m"]]))) chk$problems <- c(chk$problems, "'measurements' must be a whole number >= 2")
  if (!is.null(v[["k"]]) && (v[["k"]] < 1 || v[["k"]] != round(v[["k"]]))) chk$problems <- c(chk$problems, "'groups' must be a whole number >= 1")
  if (!is.null(v$n_total) && !is.null(v[["k"]]) && v$n_total <= v[["k"]]) {
    chk$problems <- c(chk$problems, "'n_total' must exceed 'groups'")
  }
  if (!length(chk$problems) && !is.null(v$rho) && !is.null(v[["m"]])) {
    if (effect == "between" && 1 + (v[["m"]] - 1) * v$rho <= 0) chk$problems <- c(chk$problems, "1 + (m - 1) rho must be positive")
  }
  if (length(chk$problems)) .ig_abort_invalid(chk$problems, "invalid effect-size input")
  # factor c(code) with f(V)^2 = c * f_code^2
  fac <- function(code) {
    switch(code,
           gpower30 = , on = if (effect == "between") v[["m"]] / (1 + (v[["m"]] - 1) * v$rho) else v[["m"]] / (1 - v$rho),
           gpower30_implicit_rho = , cohen1988 = , off = 1,
           spss = if (effect == "between") (v$n_total - v[["k"]]) / v$n_total else
             (v$n_total - v[["k"]]) * (v[["m"]] - 1) / v$n_total)
  }
  out <- v$f * sqrt(fac(v$from) / fac(v$to))
  warn <- if (any(uses == "spss")) "the SPSS convention depends on N and k (f(U)^2 = lambda / df2), so the converted value holds for this design only" else character()
  v[["m"]] <- NULL
  v[["k"]] <- NULL
  .ig_rm_drawer_result(out, sprintf("f_to = f_from * sqrt(c(from) / c(to)), f(V)^2 = c * f^2 (c = %s for '%s', %s for '%s')",
                                    format(fac(v$from), digits = 10), v$from, format(fac(v$to), digits = 10), v$to),
                       "equal noncentrality per subject under both conventions (method notes; not a G*Power drawer mode)",
                       v, warn)
}

.ig_rm_criterion_to_f2 <- function(value, approximation, p, q, n_total = NULL, k = NULL) {
  s <- min(p, q)
  switch(approximation,
         pillai_v = value / (s - value),
         wilks_u = value^(-1 / .ig_rm_rao_g(p, q)) - 1,
         hotelling_lawley_t1 = value / s,
         hotelling_lawley_t2 = {
           ap <- .ig_rm_manova_approx(p, q, n_total - k, 1, "hotelling_lawley_t2")
           value / ap$mult
         })
}

.ig_rm_drawer_from_criterion <- function(inputs) {
  chk <- .ig_rm_drawer_inputs(inputs, c("criterion", "approximation", "groups", "measurements"), optional = "n_total",
                              codes = list(approximation = names(.ig_rm_approx_labels)))
  v <- chk$values
  if (!length(chk$problems)) {
    v[["k"]] <- v$groups
    v[["m"]] <- v$measurements
    if (v[["k"]] < 2 || v[["k"]] != round(v[["k"]])) chk$problems <- c(chk$problems, "'groups' must be a whole number >= 2")
    if (v[["m"]] < 2 || v[["m"]] != round(v[["m"]])) chk$problems <- c(chk$problems, "'measurements' must be a whole number >= 2")
  }
  if (!length(chk$problems)) {
    p <- v[["m"]] - 1
    q <- v[["k"]] - 1
    s <- min(p, q)
    x <- v$criterion
    bad <- switch(v$approximation,
                  pillai_v = if (x <= 0 || x >= s) sprintf("Pillai V must lie in (0, s) = (0, %s)", s),
                  wilks_u = if (x <= 0 || x >= 1) "Wilks U must lie in (0, 1)",
                  if (x <= 0) "the Hotelling-Lawley trace must be positive")
    if (!is.null(bad)) chk$problems <- c(chk$problems, bad)
    if (v$approximation == "hotelling_lawley_t2" && is.null(v$n_total)) {
      chk$problems <- c(chk$problems, "'n_total' is needed for Hotelling-Lawley T2 (its effect size depends on N)")
    }
  }
  if (length(chk$problems)) .ig_abort_invalid(chk$problems, "invalid effect-size input")
  p <- v[["m"]] - 1
  q <- v[["k"]] - 1
  f2 <- .ig_rm_criterion_to_f2(v$criterion, v$approximation, p, q, v$n_total, v[["k"]])
  if (!is.finite(f2) || f2 < 0) {
    .ig_abort_invalid("the criterion gives no admissible effect size for this design (check N for T2)",
                      "invalid effect-size input")
  }
  v[["k"]] <- NULL
  v[["m"]] <- NULL
  .ig_rm_drawer_result(sqrt(f2), .ig_rm_f2_formula(v$approximation),
                       "MANOVA effect size of the approximation (G*Power drawer 'Effectsize from criterion')", v,
                       extra = list(s = min(p, q)))
}

.ig_rm_f2_formula <- function(approximation) {
  switch(approximation,
         pillai_v = "f(V)^2 = V / (s - V), s = min(m - 1, k - 1)",
         wilks_u = "f(U)^2 = U^(-1/g) - 1, g = sqrt((p^2 q^2 - 4) / (p^2 + q^2 - 5)) (1 when p^2 + q^2 <= 5)",
         hotelling_lawley_t1 = "f(T1)^2 = T / s",
         hotelling_lawley_t2 = "f(T2)^2 = T (nu_e - p - 1) / (df2 - 2) with McKeon's df2 (depends on N)")
}

# Population criterion of the group x measurement interaction from the
# k x m matrix of means, the m x m covariance matrix (or SDs with a
# correlation matrix or a common correlation) and group weights.
.ig_rm_drawer_from_means_cov <- function(inputs) {
  chk <- .ig_rm_drawer_inputs(inputs, c("means", "approximation"),
                              optional = c("sigma", "sd", "corr", "rho", "n", "n_total"),
                              vectors = c("sd", "n"), matrices = c("means", "sigma", "corr"),
                              codes = list(approximation = names(.ig_rm_approx_labels)))
  v <- chk$values
  if (length(chk$problems)) .ig_abort_invalid(chk$problems, "invalid effect-size input")
  M <- v$means
  k <- nrow(M)
  m <- ncol(M)
  probs <- character()
  if (k < 2 || m < 2) probs <- c(probs, "'means' must be a k x m matrix with k >= 2 groups and m >= 2 measurements")
  given <- c(sigma = !is.null(v$sigma), sd = !is.null(v$sd))
  if (sum(given) != 1L) probs <- c(probs, "give either 'sigma' (covariance matrix) or 'sd' with 'corr' or 'rho'")
  if (given[["sd"]] && (is.null(v$corr) == is.null(v$rho))) probs <- c(probs, "with 'sd' give exactly one of 'corr' and 'rho'")
  if (!is.null(v[["n"]]) && (length(v[["n"]]) != k || any(v[["n"]] <= 0))) probs <- c(probs, "'n' must give one positive size per group")
  if (length(probs)) .ig_abort_invalid(probs, "invalid effect-size input")
  S <- if (given[["sigma"]]) v$sigma else {
    sd <- if (length(v$sd) == 1L) rep(v$sd, m) else v$sd
    R <- if (!is.null(v$corr)) v$corr else {
      r <- matrix(v$rho, m, m)
      diag(r) <- 1
      r
    }
    if (length(sd) != m || any(sd <= 0)) .ig_abort_invalid("'sd' must give one positive SD per measurement (or one common SD)",
                                                          "invalid effect-size input")
    if (!all(dim(R) == c(m, m)) || any(abs(diag(R) - 1) > 1e-12) || any(abs(R) > 1)) {
      .ig_abort_invalid("'corr' must be an m x m correlation matrix", "invalid effect-size input")
    }
    outer(sd, sd) * R
  }
  if (!all(dim(S) == c(m, m)) || max(abs(S - t(S))) > 1e-12 * max(abs(S))) {
    .ig_abort_invalid("the covariance matrix must be a symmetric m x m matrix", "invalid effect-size input")
  }
  ev <- eigen((S + t(S)) / 2, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-12 * max(abs(ev))) {
    .ig_abort_invalid(sprintf("the covariance matrix must be positive definite (smallest eigenvalue %s)",
                              format(min(ev), digits = 4)), "invalid effect-size input")
  }
  w <- if (is.null(v[["n"]])) rep(1 / k, k) else v[["n"]] / sum(v[["n"]])
  # Orthonormal contrasts orthogonal to the unit vector (the within-subject
  # contrasts); the criterion does not depend on the basis chosen.
  C <- qr.Q(qr(cbind(1, diag(m))))[, -1, drop = FALSE][, seq_len(m - 1), drop = FALSE]
  mu <- colSums(w * M)
  D <- sweep(M, 2, mu)
  B <- crossprod(D * sqrt(w))
  Bc <- crossprod(C, B %*% C)
  Sc <- crossprod(C, S %*% C)
  theta <- Re(eigen(solve(Sc, Bc), only.values = TRUE)$values)
  theta <- pmax(theta, 0)
  p <- m - 1
  q <- k - 1
  s <- min(p, q)
  crit <- switch(v$approximation,
                 pillai_v = sum(theta / (1 + theta)),
                 wilks_u = prod(1 / (1 + theta)),
                 hotelling_lawley_t1 = sum(theta),
                 hotelling_lawley_t2 = sum(theta))
  if (v$approximation == "hotelling_lawley_t2" && is.null(v$n_total)) {
    .ig_abort_invalid("'n_total' is needed for Hotelling-Lawley T2 (its effect size depends on N)",
                      "invalid effect-size input")
  }
  f2 <- if (v$approximation == "hotelling_lawley_t2") {
    .ig_rm_criterion_to_f2(crit, "hotelling_lawley_t2", p, q, v$n_total, k)
  } else {
    .ig_rm_criterion_to_f2(crit, v$approximation, p, q)
  }
  warn <- character()
  if (v$approximation == "hotelling_lawley_t2") {
    warn <- "for Hotelling-Lawley T2 the population trace is converted with N-dependent McKeon df (as G*Power's T2 output)"
  }
  .ig_rm_drawer_result(sqrt(f2), paste(.ig_rm_f2_formula(v$approximation),
                                       "; criterion from the eigenvalues of (C' Sigma C)^-1 C' B C, B = sum_i w_i (mu_i - mu)(mu_i - mu)'"),
                       "population MANOVA criterion of the interaction (O'Brien and Shieh 1999); not compared with G*Power (grid input)",
                       c(v, list(sigma_used = S)), warn,
                       extra = list(criterion = crit, eigenvalues = theta, s = s))
}
