# Effect-size drawers: correlation and regression (group G5) ---------------------------
#
# Drawers of the nine G5 procedures (R/methods-correlation-regression.R),
# reached through effect_size() as `.ig_drawer_<method_id>(input_mode,
# inputs)`. Each returns the effect parameter of its method as `value` plus
# every quantity G*Power's drawer displays in `outputs` (named list), with
# the formula and its source. Input modes and labels follow the certified
# live inventory (registry effect_size$drawer) and the drawer captures filed
# under util_effect_size_drawers.
#
# | method | input modes |
# |---|---|
# | t_correlation_point_biserial | r_squared |
# | exact_correlation_bivariate_normal | r_squared |
# | t_linreg_one_group_slope | rho_residual_sdx, rho_slope_sdx, rho_sdx_sdy, residual_slope_sdx, residual_sdx_sdy |
# | t_linreg_two_groups_slopes | sdx_sdy_rho, sdx_sdy_slope, sdx_rho_slope |
# | t_linreg_two_groups_intercepts | sdx_sdy_rho, sdx_sdy_slope, sdx_rho_slope |
# | t_mrc_single_coefficient, f_mrc_r2_increase | variances, direct |
# | f_mrc_r2_deviation_from_zero | r_squared, predictor_correlations |
# | exact_mrc_random_model | confidence_interval, predictor_correlations |
#
# Relations (manual 12.1, 17.1): sigma_y = b sigma_x / rho = sigma_res /
# sqrt(1 - rho^2); the pooled residual SD of two groups with allocation ratio
# r = n2/n1 is sqrt((sigma_r1^2 + r sigma_r2^2) / (1 + r)) (manual 17.3, eq. 3).

# Validation of drawer inputs: names, numbers, ranges. `spec` is a named list
# of checks: "real", "positive", "probability" (0 <= x < 1), "corr" (-1 < x < 1),
# "unit_open" (0 < x < 1), "count" (whole number >= 1), "vector", "matrix".
.ig_g5_drawer_inputs <- function(input_mode, modes, required, inputs) {
  if (!.ig_is_string(input_mode) || !(input_mode %in% modes)) {
    .ig_abort_invalid(sprintf("`input_mode` must be one of: %s", paste(modes, collapse = ", ")),
                      "invalid effect-size input")
  }
  problems <- character()
  unknown <- setdiff(names(inputs), names(required))
  if (length(unknown)) {
    problems <- c(problems, sprintf("input '%s' is not used in input mode '%s' (inputs: %s)", unknown, input_mode,
                                    paste(names(required), collapse = ", ")))
  }
  vals <- list()
  for (nm in names(required)) {
    v <- inputs[[nm]]
    kind <- required[[nm]]
    if (is.null(v)) {
      problems <- c(problems, sprintf("missing input '%s'", nm))
      next
    }
    if (kind == "vector") {
      if (!is.numeric(v) || !length(v) || any(!is.finite(v))) problems <- c(problems, sprintf("'%s' must be a finite numeric vector", nm))
      vals[[nm]] <- as.double(v)
      next
    }
    if (kind == "matrix") {
      if (!is.numeric(v) || any(!is.finite(v))) {
        problems <- c(problems, sprintf("'%s' must be a finite numeric matrix", nm))
      }
      vals[[nm]] <- v
      next
    }
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v)) {
      problems <- c(problems, sprintf("'%s' must be a single finite number", nm))
      next
    }
    v <- as.double(v)
    bad <- switch(kind,
                  real = FALSE,
                  positive = v <= 0,
                  probability = v < 0 || v >= 1,
                  corr = v <= -1 || v >= 1,
                  unit_open = v <= 0 || v >= 1,
                  unit_closed = v < 0 || v > 1,
                  count = v < 1 || v != round(v),
                  FALSE)
    if (bad) {
      problems <- c(problems, sprintf("'%s' is outside its range (%s)", nm, switch(kind,
        positive = "> 0", probability = "[0, 1)", corr = "(-1, 1)", unit_open = "(0, 1)",
        unit_closed = "[0, 1]", count = "whole number >= 1", "")))
    }
    vals[[nm]] <- v
  }
  if (length(problems)) .ig_abort_invalid(problems, "invalid effect-size input")
  vals
}

.ig_g5_drawer_result <- function(effect, value, outputs, formula, source, inputs, warnings = character()) {
  list(effect = effect, value = value, outputs = outputs, formula = formula, source = source,
       inputs = inputs, warnings = warnings)
}

# Squared correlation -> correlation (point biserial and bivariate normal).
.ig_g5_drawer_r_squared <- function(input_mode, inputs, effect, label) {
  v <- .ig_g5_drawer_inputs(input_mode, "r_squared", list(r_squared = "probability"), inputs)
  r <- sqrt(v$r_squared)
  .ig_g5_drawer_result(effect, r, stats::setNames(list(r), effect), "rho = sqrt(r_squared)",
                       sprintf("G*Power drawer 'from coefficient of determination' (%s); the positive root", label), v)
}

.ig_drawer_t_correlation_point_biserial <- function(input_mode, inputs) {
  .ig_g5_drawer_r_squared(input_mode, inputs, "rho", "manual 16, drawer_T_Corr_PBM / RT_DIALOG 156")
}

.ig_drawer_exact_correlation_bivariate_normal <- function(input_mode, inputs) {
  .ig_g5_drawer_r_squared(input_mode, inputs, "rho_h1", "manual 3, Fig. 5")
}

# One group, size of slope (manual 12.1, eqs. 1-2).
.ig_drawer_t_linreg_one_group_slope <- function(input_mode, inputs) {
  modes <- c("rho_residual_sdx", "rho_slope_sdx", "rho_sdx_sdy", "residual_slope_sdx", "residual_sdx_sdy")
  req <- switch(if (.ig_is_string(input_mode) && input_mode %in% modes) input_mode else "rho_residual_sdx",
                rho_residual_sdx = list(rho = "corr", sd_residual = "positive", sd_x = "positive"),
                rho_slope_sdx = list(rho = "corr", slope = "real", sd_x = "positive"),
                rho_sdx_sdy = list(rho = "corr", sd_x = "positive", sd_y = "positive"),
                residual_slope_sdx = list(sd_residual = "positive", slope = "real", sd_x = "positive"),
                residual_sdx_sdy = list(sd_residual = "positive", sd_x = "positive", sd_y = "positive"))
  v <- .ig_g5_drawer_inputs(input_mode, modes, req, inputs)
  src <- "manual 12.1: sigma_y = b sigma_x / rho = sigma / sqrt(1 - rho^2); drawer_T_LinRegress / RT_DIALOG 229"
  if (input_mode == "rho_residual_sdx") {
    sd_y <- v$sd_residual / sqrt(1 - v$rho^2)
    slope <- v$rho * sd_y / v$sd_x
    return(.ig_g5_drawer_result("slope_h1", slope, list(sd_y = sd_y, slope_h1 = slope),
                                "sd_y = sd_residual / sqrt(1 - rho^2); slope = rho sd_y / sd_x", src, v))
  }
  if (input_mode == "rho_slope_sdx") {
    if (v$rho == 0) .ig_abort_invalid("'rho' must not be 0 in input mode rho_slope_sdx (sd_y = slope sd_x / rho)",
                                      "invalid effect-size input")
    sd_y <- v$slope * v$sd_x / v$rho
    if (sd_y <= 0) .ig_abort_invalid("slope and rho must have the same sign (sd_y = slope sd_x / rho > 0)",
                                     "invalid effect-size input")
    return(.ig_g5_drawer_result("slope_h1", v$slope, list(sd_y = sd_y, slope_h1 = v$slope),
                                "sd_y = slope sd_x / rho", src, v))
  }
  if (input_mode == "rho_sdx_sdy") {
    slope <- v$rho * v$sd_y / v$sd_x
    return(.ig_g5_drawer_result("slope_h1", slope, list(slope_h1 = slope), "slope = rho sd_y / sd_x", src, v))
  }
  if (input_mode == "residual_slope_sdx") {
    sd_y <- sqrt(v$sd_residual^2 + v$slope^2 * v$sd_x^2)
    return(.ig_g5_drawer_result("slope_h1", v$slope, list(sd_y = sd_y, slope_h1 = v$slope),
                                "sd_y = sqrt(sd_residual^2 + slope^2 sd_x^2)", src, v))
  }
  if (v$sd_residual >= v$sd_y) {
    .ig_abort_invalid("sd_residual must be below sd_y (sd_residual = sd_y sqrt(1 - rho^2))", "invalid effect-size input")
  }
  slope <- sqrt(v$sd_y^2 - v$sd_residual^2) / v$sd_x
  .ig_g5_drawer_result("slope_h1", slope, list(slope_h1 = slope),
                       "slope = sqrt(sd_y^2 - sd_residual^2) / sd_x (positive root)", src, v)
}

# Two groups: slopes (manual 17.3) and intercepts (app only; same relations).
.ig_g5_pooled_residual <- function(sr1, sr2, r) sqrt((sr1^2 + r * sr2^2) / (1 + r))

.ig_drawer_t_linreg_two_groups_slopes <- function(input_mode, inputs) {
  modes <- c("sdx_sdy_rho", "sdx_sdy_slope", "sdx_rho_slope")
  m <- if (.ig_is_string(input_mode) && input_mode %in% modes) input_mode else "sdx_sdy_rho"
  common <- list(sd_x1 = "positive", sd_x2 = "positive", allocation_ratio = "positive")
  req <- c(common, switch(m,
    sdx_sdy_rho = list(sd_y1 = "positive", rho1 = "corr", sd_y2 = "positive", rho2 = "corr"),
    sdx_sdy_slope = list(sd_y1 = "positive", slope1 = "real", sd_y2 = "positive", slope2 = "real"),
    sdx_rho_slope = list(rho1 = "corr", slope1 = "real", rho2 = "corr", slope2 = "real")))
  v <- .ig_g5_drawer_inputs(input_mode, modes, req, inputs)
  if (m == "sdx_sdy_rho") {
    b1 <- v$rho1 * v$sd_y1 / v$sd_x1
    b2 <- v$rho2 * v$sd_y2 / v$sd_x2
    sy1 <- v$sd_y1
    sy2 <- v$sd_y2
    r1 <- v$rho1
    r2 <- v$rho2
  } else if (m == "sdx_sdy_slope") {
    b1 <- v$slope1
    b2 <- v$slope2
    sy1 <- v$sd_y1
    sy2 <- v$sd_y2
    r1 <- b1 * v$sd_x1 / sy1
    r2 <- b2 * v$sd_x2 / sy2
    if (abs(r1) >= 1 || abs(r2) >= 1) {
      .ig_abort_invalid("|slope sd_x / sd_y| must be below 1 in each group", "invalid effect-size input")
    }
  } else {
    if (v$rho1 == 0 || v$rho2 == 0) .ig_abort_invalid("rho1 and rho2 must not be 0 in input mode sdx_rho_slope",
                                                       "invalid effect-size input")
    b1 <- v$slope1
    b2 <- v$slope2
    r1 <- v$rho1
    r2 <- v$rho2
    sy1 <- b1 * v$sd_x1 / r1
    sy2 <- b2 * v$sd_x2 / r2
    if (sy1 <= 0 || sy2 <= 0) .ig_abort_invalid("slope and rho must have the same sign in each group",
                                                "invalid effect-size input")
  }
  sr <- .ig_g5_pooled_residual(sy1 * sqrt(1 - r1^2), sy2 * sqrt(1 - r2^2), v$allocation_ratio)
  d <- abs(b1 - b2)
  outputs <- switch(m,
    sdx_sdy_rho = list(slope1 = b1, slope2 = b2, sd_residual = sr, delta_slope = d),
    sdx_sdy_slope = list(rho1 = r1, rho2 = r2, sd_residual = sr, delta_slope = d),
    sdx_rho_slope = list(sd_y1 = sy1, sd_y2 = sy2, sd_residual = sr, delta_slope = d))
  .ig_g5_drawer_result("delta_slope", d, outputs,
                       "b_i = rho_i sd_yi / sd_xi; sd_ri = sd_yi sqrt(1 - rho_i^2); sd_residual = sqrt((sd_r1^2 + r sd_r2^2) / (1 + r)); delta = |b1 - b2|",
                       "manual 17.3 (eqs. 3-5); drawer_T_LinRegress2G / RT_DIALOG 230", v)
}

.ig_drawer_t_linreg_two_groups_intercepts <- function(input_mode, inputs) {
  modes <- c("sdx_sdy_rho", "sdx_sdy_slope", "sdx_rho_slope")
  m <- if (.ig_is_string(input_mode) && input_mode %in% modes) input_mode else "sdx_sdy_rho"
  common <- list(sd_x1 = "positive", intercept1 = "real", intercept2 = "real")
  req <- c(common, switch(m,
    sdx_sdy_rho = list(sd_y1 = "positive", rho = "corr"),
    sdx_sdy_slope = list(sd_y1 = "positive", slope1 = "real"),
    sdx_rho_slope = list(rho = "corr", slope1 = "real")))
  v <- .ig_g5_drawer_inputs(input_mode, modes, req, inputs)
  if (m == "sdx_sdy_rho") {
    sy <- v$sd_y1
    r <- v$rho
    b <- r * sy / v$sd_x1
  } else if (m == "sdx_sdy_slope") {
    sy <- v$sd_y1
    b <- v$slope1
    r <- b * v$sd_x1 / sy
    if (abs(r) >= 1) .ig_abort_invalid("|slope1 sd_x1 / sd_y1| must be below 1", "invalid effect-size input")
  } else {
    if (v$rho == 0) .ig_abort_invalid("rho must not be 0 in input mode sdx_rho_slope", "invalid effect-size input")
    r <- v$rho
    b <- v$slope1
    sy <- b * v$sd_x1 / r
    if (sy <= 0) .ig_abort_invalid("slope1 and rho must have the same sign", "invalid effect-size input")
  }
  sr <- sy * sqrt(1 - r^2)
  d <- abs(v$intercept2 - v$intercept1)
  outputs <- list(sd_residual = sr, delta_intercept = d, slope1 = b, rho = r, sd_y1 = sy)
  .ig_g5_drawer_result("delta_intercept", d, outputs,
                       "sd_residual = sd_y1 sqrt(1 - rho^2) with sd_y1 = slope1 sd_x1 / rho; delta = |a2 - a1|",
                       "G*Power drawer drawer_T_LinRegress2GI / RT_DIALOG 231 (no manual chapter); relations of manual 12.1", v)
}

# f^2 from variances or partial R^2 (single coefficient and R^2 increase).
.ig_g5_drawer_f2_partial <- function(input_mode, inputs, explained_label, source) {
  modes <- c("variances", "direct")
  m <- if (.ig_is_string(input_mode) && input_mode %in% modes) input_mode else "variances"
  req <- switch(m, variances = list(var_explained = "positive", var_residual = "positive"),
                direct = list(partial_r2 = "unit_open"))
  v <- .ig_g5_drawer_inputs(input_mode, modes, req, inputs)
  if (m == "variances") {
    f2 <- v$var_explained / v$var_residual
    pr2 <- v$var_explained / (v$var_explained + v$var_residual)
    return(.ig_g5_drawer_result("f2", f2, list(partial_r2 = pr2, f2 = f2),
                                sprintf("f2 = %s / residual variance; partial R^2 = f2 / (1 + f2)", explained_label), source, v))
  }
  f2 <- v$partial_r2 / (1 - v$partial_r2)
  .ig_g5_drawer_result("f2", f2, list(f2 = f2), "f2 = partial R^2 / (1 - partial R^2)", source, v)
}

.ig_drawer_t_mrc_single_coefficient <- function(input_mode, inputs) {
  .ig_g5_drawer_f2_partial(input_mode, inputs, "variance explained by the predictor",
                           "G*Power drawer drawer_T_MRC_1Pred / RT_DIALOG 234 (no manual chapter)")
}

.ig_drawer_f_mrc_r2_increase <- function(input_mode, inputs) {
  .ig_g5_drawer_f2_partial(input_mode, inputs, "variance explained by the special effect",
                           "manual 14.1-14.2 (f2 = V_S / V_E); drawer_F_MRC_Special / RT_DIALOG 162/240")
}

# rho^2 = u B^-1 u' from predictor correlations (manual 7.1, 13.1).
.ig_g5_rho2_from_correlations <- function(u, B) {
  B <- as.matrix(B)
  m <- length(u)
  if (!identical(dim(B), c(m, m))) {
    .ig_abort_invalid(sprintf("'B' must be a %d x %d matrix (one row and column per predictor)", m, m),
                      "invalid effect-size input")
  }
  problems <- character()
  if (any(abs(u) > 1)) problems <- c(problems, "every correlation in 'u' must lie in [-1, 1]")
  if (any(abs(B) > 1) || any(abs(diag(B) - 1) > 1e-12) || max(abs(B - t(B))) > 1e-12) {
    problems <- c(problems, "'B' must be a symmetric correlation matrix (unit diagonal, entries in [-1, 1])")
  }
  if (length(problems)) .ig_abort_invalid(problems, "invalid effect-size input")
  ev <- eigen(B, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-12 * max(1, max(ev))) {
    .ig_abort_invalid(sprintf("'B' is not positive definite (smallest eigenvalue %s)", format(signif(min(ev), 4))),
                      "invalid effect-size input")
  }
  rho2 <- as.double(crossprod(u, solve(B, u)))
  if (rho2 < 0 || rho2 > 1) {
    .ig_abort_invalid(sprintf("u B^-1 u' = %s lies outside [0, 1]; the correlations are inconsistent", format(rho2)),
                      "invalid effect-size input")
  }
  rho2
}

.ig_drawer_f_mrc_r2_deviation_from_zero <- function(input_mode, inputs) {
  modes <- c("r_squared", "predictor_correlations")
  m <- if (.ig_is_string(input_mode) && input_mode %in% modes) input_mode else "r_squared"
  req <- switch(m, r_squared = list(rho2 = "probability"), predictor_correlations = list(u = "vector", B = "matrix"))
  v <- .ig_g5_drawer_inputs(input_mode, modes, req, inputs)
  rho2 <- if (m == "r_squared") v$rho2 else .ig_g5_rho2_from_correlations(v$u, v$B)
  if (rho2 >= 1) .ig_abort_invalid("rho^2 must be below 1", "invalid effect-size input")
  f2 <- rho2 / (1 - rho2)
  .ig_g5_drawer_result("f2", f2, list(rho2 = rho2, f2 = f2),
                       if (m == "r_squared") "f2 = rho2 / (1 - rho2)" else "rho2 = u B^-1 u'; f2 = rho2 / (1 - rho2)",
                       "manual 13.1-13.2; drawer_F_MRC_Omni / RT_DIALOG 161 + 241", v)
}

# Random model: H1 rho^2 from a confidence interval of an observed R^2, or from
# predictor correlations (manual 7.1). The confidence limits invert the exact
# R^2 distribution in rho^2 (P(R^2 <= R2obs | rho^2) is decreasing in rho^2):
#   two-sided:  P(R^2 >= R2obs | L) = alpha/2,  P(R^2 <= R2obs | U) = alpha/2
#   one-sided:  (L1, 1) with P(R^2 >= R2obs | L1) = alpha, (0, U1) with P(R^2 <= R2obs | U1) = alpha
# (a limit is 0 when the equation has no root in [0, 1)). H1 rho^2 =
# L + position (U - L) with position in [0, 1] ("Rel C.I. pos to use").
.ig_drawer_exact_mrc_random_model <- function(input_mode, inputs) {
  modes <- c("confidence_interval", "predictor_correlations")
  m <- if (.ig_is_string(input_mode) && input_mode %in% modes) input_mode else "confidence_interval"
  if (m == "predictor_correlations") {
    v <- .ig_g5_drawer_inputs(input_mode, modes, list(u = "vector", B = "matrix"), inputs)
    rho2 <- .ig_g5_rho2_from_correlations(v$u, v$B)
    return(.ig_g5_drawer_result("rho2_h1", rho2, list(rho2_h1 = rho2), "rho2 = u B^-1 u'",
                                "manual 7.1 (effect size from predictor correlations)", v))
  }
  v <- .ig_g5_drawer_inputs(input_mode, modes,
                            list(n_total = "count", n_predictors = "count", r2_observed = "unit_open",
                                 confidence = "unit_open", position = "unit_closed"), inputs)
  N <- v$n_total
  k <- v$n_predictors
  if (N < k + 2) .ig_abort_invalid("n_total must be at least n_predictors + 2", "invalid effect-size input")
  alpha <- 1 - v$confidence
  obs <- v$r2_observed
  limit <- function(prob, upper_tail) {
    # upper_tail: P(R^2 >= obs | rho2) = prob (increasing in rho2) -> lower
    # limit; otherwise P(R^2 <= obs | rho2) = prob (decreasing) -> upper limit.
    # g(rho2) = +-(log tail - log prob) is increasing; a limit is 0 when g(0) >= 0.
    g <- function(rho2) {
      t <- .ig_g5_r2_tail(obs, rho2, k, N, lower_tail = !upper_tail)
      if (is.na(t)) return(NA_real_)
      v <- log(max(t, .Machine$double.xmin)) - log(prob)
      if (upper_tail) v else -v
    }
    g0 <- g(0)
    if (is.na(g0)) .ig_abort("numerical_failure", "the R^2 distribution could not be evaluated at rho^2 = 0")
    if (g0 >= 0) return(0)
    lo <- 0
    glo <- g0
    for (hi in 1 - 10^-(seq(0.3, 12, by = 0.3))) {
      ghi <- g(hi)
      if (is.na(ghi)) {
        .ig_abort("numerical_failure", sprintf(
          "the confidence limit lies above rho^2 = %s, where the exact R^2 distribution is not computed", format(lo)))
      }
      if (ghi >= 0) return(stats::uniroot(g, c(lo, hi), f.lower = glo, f.upper = ghi, tol = 1e-13)$root)
      lo <- hi
      glo <- ghi
    }
    .ig_abort("numerical_failure", "no confidence limit below rho^2 = 1 - 1e-12")
  }
  ci_lo <- limit(alpha / 2, TRUE)
  ci_hi <- limit(alpha / 2, FALSE)
  b_lo <- limit(alpha, TRUE)
  b_hi <- limit(alpha, FALSE)
  rho2 <- ci_lo + v$position * (ci_hi - ci_lo)
  .ig_g5_drawer_result("rho2_h1", rho2,
                       list(ci_lower = ci_lo, ci_upper = ci_hi, bound_lower = b_lo, bound_upper = b_hi, rho2_h1 = rho2),
                       "confidence limits for rho^2 by inverting the exact R^2 distribution; rho2 = L + position (U - L)",
                       "manual 7.1 and 7.3.2 (effect size from C.I.; DOC-051, DOC-052)", v)
}
