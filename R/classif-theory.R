# Binormal theory model ----------------------------------------------------------------
#
# Non-events S0 ~ N(0, 1), events S1 ~ N(d, sigma^2) with
# d = qnorm(AUC) * sqrt(1 + sigma^2), so that P(S1 > S0) = AUC exactly. A
# cutoff c is a real number on the score scale standardised to the non-event
# distribution; it is preserved when AUC or sigma change. Tails are analytic:
#
#   TPR(c) = P(S1 >= c) = pnorm((c - d) / sigma, lower.tail = FALSE)
#   FPR(c) = P(S0 >= c) = pnorm(c, lower.tail = FALSE)
#
# (>= and > coincide for continuous distributions). The model N scales
# expected counts only; rates do not depend on it. Theoretical quantities have
# no empirical confidence interval. Grids are used for display and to bracket
# optima; every reported optimum is refined numerically or analytically.

#' Binormal score model
#'
#' Parameters of the binormal model used for theoretical classification
#' scenarios: non-events \eqn{S_0 \sim N(0, 1)}, events
#' \eqn{S_1 \sim N(d, \sigma^2)} with \eqn{d = \Phi^{-1}(AUC)\sqrt{1 + \sigma^2}}.
#'
#' @param auc Area under the ROC curve in (0, 1).
#' @param sigma Ratio of the event to the non-event standard deviation (> 0).
#' @return A list with `auc`, `sigma`, `d`, the class distributions and
#'   `auc_check` (the AUC recomputed as \eqn{\Phi(d / \sqrt{1 + \sigma^2})}).
#' @examples
#' binormal_model(0.8, 1.5)
#' @concept classification
#' @seealso
#' * [classification_spec()]
#' * [evaluate_cutoff()]
#' * [plot_classification()]
#' @export
binormal_model <- function(auc, sigma = 1) {
  .ig_cls_check_binormal(auc, sigma)
  d <- stats::qnorm(auc) * sqrt(1 + sigma^2)
  list(auc = auc, sigma = sigma, d = d,
       non_event = c(mean = 0, sd = 1), event = c(mean = d, sd = sigma),
       auc_check = stats::pnorm(d / sqrt(1 + sigma^2)),
       cutoff_scale = "score standardised to the non-event distribution (z)")
}

.ig_cls_check_binormal <- function(auc, sigma) {
  num1 <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!num1(auc) || auc <= 0 || auc >= 1) .ig_abort("invalid_input", "`auc` must be a single number in (0, 1)")
  if (!num1(sigma) || sigma <= 0) .ig_abort("invalid_input", "`sigma` must be a single number > 0")
  invisible(TRUE)
}

#' Analytic rates of the binormal model at given cutoffs
#'
#' @param cutoff Numeric vector of real-valued cutoffs (alert if score >=
#'   cutoff) on the score scale standardised to the non-event distribution.
#' @inheritParams binormal_model
#' @return A data frame with `cutoff`, `tpr` and `fpr`, computed from the
#'   normal upper tails (no discretisation).
#' @examples
#' binormal_rates(c(0, 1, 2), auc = 0.8, sigma = 1)
#' @concept classification
#' @seealso
#' * [classification_spec()]
#' * [evaluate_cutoff()]
#' * [plot_classification()]
#' @export
binormal_rates <- function(cutoff, auc, sigma = 1) {
  .ig_cls_check_binormal(auc, sigma)
  if (!is.numeric(cutoff) || !length(cutoff) || anyNA(cutoff)) {
    .ig_abort("invalid_input", "`cutoff` must be numeric without NA")
  }
  d <- stats::qnorm(auc) * sqrt(1 + sigma^2)
  data.frame(cutoff = cutoff,
             tpr = stats::pnorm((cutoff - d) / sigma, lower.tail = FALSE),
             fpr = stats::pnorm(cutoff, lower.tail = FALSE))
}

.ig_cls_theory_d <- function(th) stats::qnorm(th$auc) * sqrt(1 + th$sigma^2)

# Expected counts at real-valued cutoffs (vectorised).
.ig_cls_theory_counts <- function(cutoff, th, pi) {
  d <- .ig_cls_theory_d(th)
  tpr <- stats::pnorm((cutoff - d) / th$sigma, lower.tail = FALSE)
  fpr <- stats::pnorm(cutoff, lower.tail = FALSE)
  n <- th$n_total
  P <- pi * n
  N <- (1 - pi) * n
  data.frame(tp = P * tpr, fn = P * (1 - tpr), fp = N * fpr, tn = N * (1 - fpr),
             tpr = tpr, fpr = fpr)
}

.ig_cls_theory_range <- function(th) {
  d <- .ig_cls_theory_d(th)
  c(min(-4, d - 4 * th$sigma), max(4, d + 4 * th$sigma))
}

# Display grid for the theoretical sweep (expected counts).
.ig_cls_theory_grid <- function(spec, grid = 401) {
  th <- spec$theory
  if (is.numeric(grid) && length(grid) == 1L) {
    if (!is.finite(grid) || grid < 2 || grid != round(grid)) .ig_abort("invalid_input", "`grid` must be a whole number >= 2 or a vector of cutoffs")
    r <- .ig_cls_theory_range(th)
    cuts <- seq(r[1], r[2], length.out = grid)
  } else if (is.numeric(grid) && all(is.finite(grid))) {
    cuts <- sort(unique(as.double(grid)))
  } else {
    .ig_abort("invalid_input", "`grid` must be a whole number >= 2 or a vector of finite cutoffs")
  }
  ec <- .ig_cls_theory_counts(cuts, th, spec$prevalence$value)
  data.frame(set_index = NA_real_, cutoff = cuts, cutoff_kind = "display_grid",
             boundary_score = NA_real_, alerts = ec$tp + ec$fp,
             tp = ec$tp, fn = ec$fn, fp = ec$fp, tn = ec$tn,
             comparator = spec$models[[1]]$comparator, count_kind = "expected",
             stringsAsFactors = FALSE)
}

# Theoretical curve summaries: AUC (model parameter), analytic partial AUC
# and average precision at the scenario prevalence by numerical integration.
.ig_cls_theory_curve_summaries <- function(spec, partial) {
  th <- spec$theory
  d <- .ig_cls_theory_d(th)
  s <- th$sigma
  pi <- spec$prevalence$value
  # Integrals over pieces (split points inside the range); NA when
  # stats::integrate does not reach the requested accuracy, which the rows
  # below report as "not_computed" with reason "integration_failed".
  integ <- function(f, lo, hi, at = numeric()) {
    pts <- sort(unique(c(lo, at[at > lo & at < hi], hi)))
    tryCatch(sum(vapply(seq_len(length(pts) - 1L), function(k) {
      stats::integrate(f, pts[k], pts[k + 1L], rel.tol = 1e-11, abs.tol = 0, subdivisions = 1000L)$value
    }, double(1))), error = function(e) NA_real_)
  }
  failed_row <- function(id, variant = NA_character_, setting = NA_character_) {
    .ig_cls_summary_row(id, NA, "not_computed", "integration_failed", variant = variant,
                        estimate_kind = "theoretical", setting = setting)
  }
  rows <- list(.ig_cls_summary_row("roc_auc", th$auc, estimate_kind = "theoretical",
                                   setting = "model parameter; equals P(S1 > S0)"))
  if (!is.null(partial)) {
    lo <- partial$range[1]; hi <- partial$range[2]
    setting <- sprintf("%s in [%s, %s]; analytic integral", toupper(partial$focus), format(lo), format(hi))
    if (partial$focus == "fpr") {
      # integral of TPR d(FPR) with FPR = 1 - Phi(c): substitute c; split at
      # c = d, where TPR changes fastest (a near step for small sigma)
      c_lo <- stats::qnorm(hi, lower.tail = FALSE)
      c_hi <- stats::qnorm(lo, lower.tail = FALSE)
      raw <- integ(function(z) stats::pnorm((z - d) / s, lower.tail = FALSE) * stats::dnorm(z), c_lo, c_hi,
                   at = d)
      min_area <- (hi^2 - lo^2) / 2
    } else {
      # integral of (1 - FPR) d(TPR) with TPR = 1 - Phi((c - d) / s), written
      # over the standardised event score u = (c - d) / s
      u_lo <- stats::qnorm(hi, lower.tail = FALSE)
      u_hi <- stats::qnorm(lo, lower.tail = FALSE)
      raw <- integ(function(u) stats::pnorm(d + s * u) * stats::dnorm(u), u_lo, u_hi, at = -d / s)
      min_area <- (hi - lo) - (hi^2 - lo^2) / 2
    }
    max_area <- hi - lo
    if (is.na(raw)) {
      rows <- c(rows, lapply(c("raw", "normalised", "mcclish"), function(v) failed_row("partial_auc", v, setting)))
    } else {
      mcclish <- (1 + (raw - min_area) / (max_area - min_area)) / 2
      rows <- c(rows, list(
        .ig_cls_summary_row("partial_auc", raw, variant = "raw", estimate_kind = "theoretical", setting = setting),
        .ig_cls_summary_row("partial_auc", raw / max_area, variant = "normalised", estimate_kind = "theoretical",
                            setting = setting),
        .ig_cls_summary_row("partial_auc", mcclish, if (raw >= min_area) "ok" else "undefined",
                            if (raw >= min_area) NA_character_ else "mcclish_below_diagonal",
                            variant = "mcclish", estimate_kind = "theoretical", setting = setting)))
    }
  }
  if (pi > 0 && pi < 1) {
    # AP = integral of PPV d(TPR), written over the standardised event score
    # u = (c - d) / s (weight dnorm(u), whatever sigma is); PPV in log space so
    # that both tails stay accurate.
    lodds <- log((1 - pi) / pi)
    ppv <- function(z) {
      lt <- stats::pnorm((z - d) / s, lower.tail = FALSE, log.p = TRUE)
      lf <- stats::pnorm(z, lower.tail = FALSE, log.p = TRUE)
      stats::plogis(lt - lf - lodds)
    }
    ap <- integ(function(u) ppv(d + s * u) * stats::dnorm(u), -Inf, Inf, at = c(-d / s, 0))
    setting <- sprintf("integral of PPV over recall at assumed prevalence %s", format(pi))
    rows[[length(rows) + 1L]] <- if (is.na(ap)) failed_row("average_precision", setting = setting) else
      .ig_cls_summary_row("average_precision", ap, estimate_kind = "theoretical", setting = setting)
  } else {
    rows[[length(rows) + 1L]] <- .ig_cls_summary_row(
      "average_precision", NA, "undefined", if (pi == 0) "no_events" else "no_non_events",
      estimate_kind = "theoretical")
  }
  rows[[length(rows) + 1L]] <- .ig_cls_summary_row(
    "pr_area_trapezoid", NA, "not_applicable", "theoretical_continuous_curve",
    estimate_kind = "theoretical")
  .ig_cls_summary_frame(rows)
}

# Optimisation helpers -------------------------------------------------------------------

# Closed-form maximiser of Youden's J for the binormal model. Stationary
# points solve phi(c) = phi((c - d) / sigma) / sigma.
.ig_cls_theory_youden <- function(th) {
  d <- .ig_cls_theory_d(th)
  s <- th$sigma
  J <- function(z) stats::pnorm((z - d) / s, lower.tail = FALSE) - stats::pnorm(z, lower.tail = FALSE)
  if (abs(s - 1) < 1e-12) {
    roots <- if (d != 0) d / 2 else numeric()
  } else {
    A <- 1 / s^2 - 1
    B <- -2 * d / s^2
    C <- d^2 / s^2 + 2 * log(s)
    disc <- B^2 - 4 * A * C
    roots <- if (disc < 0) numeric() else (-B + c(-1, 1) * sqrt(disc)) / (2 * A)
  }
  if (!length(roots)) return(list(cutoff = NA_real_, value = NA_real_, roots = roots))
  vals <- J(roots)
  best <- which.max(vals)
  list(cutoff = roots[best], value = vals[best], roots = roots)
}
