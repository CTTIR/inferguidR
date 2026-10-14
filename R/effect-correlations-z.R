# Effect-size drawers: correlation z tests (group G7) -------------------------------------
#
# G*Power's "Determine" dialogs of the four correlation z tests (manual
# 27.1, 28.1, 31.1), reached through effect_size() (R/core-api.R), which looks
# up `.ig_drawer_<method_id>()`. Method notes: specs/method-notes/<method_id>.md.
#
# * z_correlations_two_independent: q = atanh(r1) - atanh(r2) from two
#   correlations (drawer_Z_2Corr / RT_DIALOG 168).
# * z_correlations_two_dependent_common_index and ..._no_common_index: the H1
#   correlation from the H0 correlation and Cohen's q,
#   rho_H1 = tanh(atanh(rho_ab) - q) (manual 28.1: r2 = (a - 1)/(a + 1),
#   a = exp(-2q + ln((1 + r1)/(1 - r1)))).
# * z_tetrachoric_correlation: H1 corr rho from four cell probabilities
#   (exact rho by inverting the bivariate normal CDF, plus the Bonett-Price
#   rho*), or from a confidence interval computed from observed frequencies
#   (manual 31.1: p_ij = (f_ij + 0.5)/(N + 2); interval r -+ z s(r)).
#
# Every drawer returns the effect (`value`, signed) and, where G*Power shows
# more fields, a named list `details` with them.

.ig_g7_drawer_num <- function(inputs, required, optional = character(), input_mode) {
  problems <- character()
  unknown <- setdiff(names(inputs), c(required, optional))
  if (length(unknown)) {
    problems <- c(problems, sprintf("input '%s' is not used in input mode '%s' (inputs: %s)", unknown, input_mode,
                                    paste(c(required, optional), collapse = ", ")))
  }
  vals <- list()
  for (nm in c(required, optional)) {
    v <- inputs[[nm]]
    if (is.null(v)) {
      if (nm %in% required) problems <- c(problems, sprintf("missing input '%s'", nm))
      next
    }
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v)) {
      problems <- c(problems, sprintf("'%s' must be a single finite number", nm))
      next
    }
    vals[[nm]] <- as.double(v)
  }
  list(values = vals, problems = problems)
}

.ig_g7_check_mode <- function(input_mode, modes) {
  if (!.ig_is_string(input_mode) || !(input_mode %in% modes)) {
    .ig_abort_invalid(sprintf("`input_mode` must be one of: %s", paste(modes, collapse = ", ")),
                      "invalid effect-size input")
  }
}

.ig_g7_check_rho <- function(vals, names) {
  out <- character()
  for (nm in names) {
    v <- vals[[nm]]
    if (!is.null(v) && abs(v) >= 1) out <- c(out, sprintf("'%s' must lie in (-1, 1)", nm))
  }
  out
}

# Two independent Pearson r's: q from two correlations ------------------------------------

.ig_drawer_z_correlations_two_independent <- function(input_mode, inputs) {
  .ig_g7_check_mode(input_mode, "two_correlations")
  chk <- .ig_g7_drawer_num(inputs, c("r1", "r2"), input_mode = input_mode)
  problems <- c(chk$problems, .ig_g7_check_rho(chk$values, c("r1", "r2")))
  if (length(problems)) .ig_abort_invalid(problems, "invalid effect-size input")
  v <- chk$values
  q <- atanh(v$r1) - atanh(v$r2)
  warnings <- character()
  if (abs(q) > 10) {
    warnings <- sprintf("q = %s lies outside [-10, 10], the range G*Power accepts for this test (manual 27.1)",
                        format(signif(q, 7)))
  }
  list(effect = "q", value = q, formula = "q = atanh(r1) - atanh(r2) (difference of Fisher z values)",
       source = "manual 27.1; G*Power drawer 'Correlation coefficient r1, r2'",
       inputs = v, warnings = warnings, details = list(z1 = atanh(v$r1), z2 = atanh(v$r2)))
}

# Two dependent Pearson r's: H1 correlation from rho_ab and Cohen's q ------------------------

.ig_g7_drawer_dependent <- function(input_mode, inputs, h1) {
  .ig_g7_check_mode(input_mode, "cohens_q")
  chk <- .ig_g7_drawer_num(inputs, c("rho_ab", "q"), input_mode = input_mode)
  problems <- c(chk$problems, .ig_g7_check_rho(chk$values, "rho_ab"))
  if (length(problems)) .ig_abort_invalid(problems, "invalid effect-size input")
  v <- chk$values
  rho <- tanh(atanh(v$rho_ab) - v$q)
  list(effect = h1, value = rho,
       formula = sprintf("%s = tanh(atanh(rho_ab) - q), i.e. q = atanh(rho_ab) - atanh(%s)", h1, h1),
       source = "manual 28.1 (r2 = (a - 1)/(a + 1), a = exp(-2q + ln((1 + r1)/(1 - r1)))); G*Power drawer 'from Cohen's q'",
       inputs = v, warnings = character(), details = list())
}

.ig_drawer_z_correlations_two_dependent_common_index <- function(input_mode, inputs) {
  .ig_g7_drawer_dependent(input_mode, inputs, "rho_ac")
}

.ig_drawer_z_correlations_two_dependent_no_common_index <- function(input_mode, inputs) {
  .ig_g7_drawer_dependent(input_mode, inputs, "rho_cd")
}

# Tetrachoric correlation --------------------------------------------------------------------

# Exact tetrachoric correlation of a 2 x 2 probability table (manual
# orientation: a = p11 = P(X = 0, Y = 0), b = p12, c = p21, d = p22): the rho
# with Phi2(h, k; rho) = a, h = qnorm(a + c), k = qnorm(a + b). Phi2 is strictly
# increasing in rho (derivative phi2(h, k; rho)), so the root is unique.
.ig_g7_tetra_invert <- function(pr) {
  h <- stats::qnorm(pr[["a"]] + pr[["c"]])
  k <- stats::qnorm(pr[["a"]] + pr[["b"]])
  f <- function(r) .ig_g7_bvn(h, k, r) - pr[["a"]]
  lim <- 1 - 1e-12
  fl <- f(-lim)
  fu <- f(lim)
  if (fl > 0 || fu < 0) {
    return(list(rho = NA_real_, h = h, k = k,
                note = "the table is outside the range of the tetrachoric model (no rho in (-1, 1) gives p11)"))
  }
  r <- stats::uniroot(f, c(-lim, lim), f.lower = fl, f.upper = fu, tol = 1e-15, maxiter = 1000L)$root
  list(rho = r, h = h, k = k, note = NULL)
}

.ig_drawer_z_tetrachoric_correlation <- function(input_mode, inputs) {
  .ig_g7_check_mode(input_mode, c("cell_probabilities", "ci_from_frequencies"))
  comp <- inputs$computation %||% "exact_brown_benedetti"
  inputs$computation <- NULL
  if (!.ig_is_string(comp) || !(comp %in% c("exact_brown_benedetti", "approx_bonett_price"))) {
    .ig_abort_invalid("`computation` must be 'exact_brown_benedetti' or 'approx_bonett_price'",
                      "invalid effect-size input")
  }
  if (input_mode == "cell_probabilities") {
    cells <- c("p11", "p12", "p21", "p22")
    chk <- .ig_g7_drawer_num(inputs, character(), optional = cells, input_mode = input_mode)
    problems <- chk$problems
    v <- chk$values
    given <- intersect(cells, names(v))
    if (length(given) < 3L) problems <- c(problems, "give at least three of p11, p12, p21, p22")
    if (length(problems)) .ig_abort_invalid(problems, "invalid effect-size input")
    filled <- character()
    if (length(given) == 3L) {
      miss <- setdiff(cells, given)
      v[[miss]] <- 1 - sum(unlist(v[given]))
      filled <- miss
    }
    pv <- unlist(v[cells])
    if (any(pv <= 0 | pv >= 1)) problems <- c(problems, "every cell probability must lie in (0, 1)")
    if (abs(sum(pv) - 1) > 1e-9) {
      problems <- c(problems, sprintf("the four cell probabilities must sum to 1 (sum %s)", format(sum(pv), digits = 15)))
    }
    if (length(problems)) .ig_abort_invalid(problems, "invalid effect-size input")
    pr <- c(a = pv[["p11"]], b = pv[["p12"]], c = pv[["p21"]], d = pv[["p22"]])
    inv <- .ig_g7_tetra_invert(pr)
    px <- pr[["b"]] + pr[["d"]]
    py <- pr[["c"]] + pr[["d"]]
    bp <- .ig_g7_tetra_bp(pr, px, py, 1)
    details <- list(correlation_rho = inv$rho, approx_correlation_rho_star = bp$r,
                    marginal_prob_x = px, marginal_prob_y = py, filled_cell = filled)
    return(list(effect = "rho_h1", value = inv$rho,
                formula = paste("rho with Phi2(qnorm(p11 + p21), qnorm(p11 + p12); rho) = p11 (exact tetrachoric",
                                "correlation); rho* = cos(pi/(1 + omega^c)) (Bonett and Price 2005)"),
                source = "manual 31.1 (effect size dialog, cell probabilities); the exact rho is used as H1 corr rho",
                inputs = v, warnings = if (!is.null(inv$note)) inv$note else character(), details = details))
  }
  # Confidence interval from observed frequencies.
  chk <- .ig_g7_drawer_num(inputs, c("f11", "f12", "f21", "f22", "position", "confidence"),
                           input_mode = input_mode)
  problems <- chk$problems
  v <- chk$values
  for (nm in c("f11", "f12", "f21", "f22")) {
    if (!is.null(v[[nm]]) && v[[nm]] < 0) problems <- c(problems, sprintf("'%s' must be non-negative", nm))
  }
  if (!is.null(v$position) && (v$position < 0 || v$position > 1)) problems <- c(problems, "'position' must lie in [0, 1]")
  if (!is.null(v$confidence) && (v$confidence <= 0 || v$confidence >= 1)) {
    problems <- c(problems, "'confidence' must lie in (0, 1)")
  }
  if (length(problems)) .ig_abort_invalid(problems, "invalid effect-size input")
  if (identical(comp, "approx_bonett_price")) {
    .ig_abort("reference_pending", paste(
      "the confidence interval from observed frequencies in approximation mode is reference_pending: the manual",
      "gives two different cell estimators ((f + 0.5)/(N + 2) in 31.1, (f + 0.5)/N in 31.5.2) and no G*Power",
      "capture of the approximation-mode dialog exists; use computation = 'exact_brown_benedetti'"))
  }
  f <- c(a = v$f11, b = v$f12, c = v$f21, d = v$f22)
  n <- sum(f)
  if (n <= 0) .ig_abort_invalid("the frequencies must not all be zero", "invalid effect-size input")
  pr <- (f + 0.5) / (n + 2)
  inv <- .ig_g7_tetra_invert(pr)
  if (is.na(inv$rho)) .ig_abort("numerical_failure", inv$note)
  se <- .ig_g7_tetra_se_exact(pr, inv$h, inv$k, inv$rho, n)
  z <- stats::qnorm((1 + v$confidence) / 2)
  lwr <- inv$rho - z * se
  upr <- inv$rho + z * se
  value <- lwr + (upr - lwr) * v$position
  warnings <- character()
  if (lwr < -1 || upr > 1) {
    warnings <- c(warnings, paste("the interval r -+ z s(r) extends beyond [-1, 1] (it is not truncated, as in",
                                  "G*Power); a value outside (-1, 1) cannot be used as H1 corr rho"))
  }
  if (min(f) < 5) {
    warnings <- c(warnings, "a cell frequency below 5: the asymptotic standard error can be poor (manual 31.0.1, 31.5.1)")
  }
  details <- list(total_sample_size = n, correlation_r = inv$rho, est_marginal_prob_x = pr[["b"]] + pr[["d"]],
                  est_marginal_prob_y = pr[["c"]] + pr[["d"]], ci_lower = lwr, ci_upper = upr,
                  std_error_r = se, cell_estimates = pr)
  list(effect = "rho_h1", value = value,
       formula = paste("p_ij = (f_ij + 0.5)/(N + 2); r: exact tetrachoric r of the estimated table; s(r): asymptotic",
                       "standard error (Pearson 1901, Brown and Benedetti 1977) with N = sum f; H1 corr rho =",
                       "L + (U - L) * position with [L, U] = r -+ qnorm((1 + confidence)/2) s(r)"),
       source = "manual 31.1 (effect size dialog, C.I. from observed frequencies, exact mode)",
       inputs = v, warnings = warnings, details = details)
}
