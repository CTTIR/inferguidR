# ROC and precision-recall summaries, own DeLong ---------------------------------------
#
# Conventions (method notes: specs/method-notes/classification.md):
#
# * ROC AUC = P(S1 > S0) + 0.5 P(S1 = S0) over all event/non-event pairs, in
#   the declared direction (never auto-flipped). It equals the trapezoidal area
#   of the empirical ROC through all decision sets, where a tied score group
#   is a straight (diagonal) segment. Frequency weights act as replication.
# * Partial AUC integrates the same piecewise-linear ROC over a stated FPR
#   range (TPR as the integrand) or TPR range ((1 - FPR) as the integrand);
#   raw, normalised (divided by the range width) and McClish-standardised
#   values are reported separately.
# * Average precision AP = sum_k (R_k - R_(k-1)) P_k over the decision sets
#   (step sum); the trapezoidal PR area sum_k (R_k - R_(k-1)) (P_k + P_(k-1)) / 2
#   with P_0 := 1 is a different quantity and is reported separately.
# * DeLong variance from placement values (ties count 1/2):
#   Var = S10 / n_events + S01 / n_non_events, S = sample variances of the
#   placement values (denominator n - 1). Interval AUC +- z sqrt(Var), clipped
#   to [0, 1]; for a fixed score set only.

# AUC, placement values and DeLong variance of one score column.
.ig_cls_auc_stats <- function(score, y, w, direction) {
  g <- .ig_cls_score_groups(score, y, w, direction)
  a <- rev(g$pos)          # ascending oriented score
  b <- rev(g$neg)
  A <- g$n_pos
  B <- g$n_neg
  out <- list(n_pos = A, n_neg = B, auc = NA_real_, var = NA_real_, status = "ok", reason = NA_character_)
  if (A == 0) return(utils::modifyList(out, list(status = "undefined", reason = "no_events")))
  if (B == 0) return(utils::modifyList(out, list(status = "undefined", reason = "no_non_events")))
  b_below <- cumsum(b) - b
  a_above <- A - cumsum(a)
  v10 <- (b_below + 0.5 * b) / B   # placement of an event with this score
  v01 <- (a_above + 0.5 * a) / A   # placement of a non-event with this score
  auc <- sum(a * v10) / A
  var <- NA_real_
  if (A >= 2 && B >= 2) {
    s10 <- sum(a * (v10 - auc)^2) / (A - 1)
    s01 <- sum(b * (v01 - auc)^2) / (B - 1)
    var <- s10 / A + s01 / B
  }
  out$auc <- auc
  out$var <- var
  out$groups <- list(u = rev(g$u), v10 = v10, v01 = v01)
  out
}

# ROC points (FPR, TPR) of the decision sets, from alert none to alert all.
.ig_cls_roc_points <- function(ds) {
  P <- ds$tp[1] + ds$fn[1]
  N <- ds$fp[1] + ds$tn[1]
  list(fpr = ds$fp / N, tpr = ds$tp / P)
}

# Area under a piecewise-linear curve through (x, y) over [lo, hi] in x;
# segments with equal x contribute nothing. Used for the ROC (x = FPR) and,
# with swapped roles, for the TPR focus.
.ig_cls_segment_area <- function(x, y, lo, hi) {
  area <- 0
  for (i in seq_len(length(x) - 1L)) {
    x0 <- x[i]; x1 <- x[i + 1L]
    if (x1 <= x0) next
    a <- max(x0, lo); b <- min(x1, hi)
    if (b <= a) next
    ya <- y[i] + (y[i + 1L] - y[i]) * (a - x0) / (x1 - x0)
    yb <- y[i] + (y[i + 1L] - y[i]) * (b - x0) / (x1 - x0)
    area <- area + (b - a) * (ya + yb) / 2
  }
  area
}

.ig_cls_check_partial <- function(partial_auc) {
  if (is.null(partial_auc)) return(NULL)
  ok <- is.list(partial_auc) && setequal(names(partial_auc), c("focus", "range")) &&
    .ig_is_string(partial_auc$focus) && partial_auc$focus %in% c("fpr", "tpr") &&
    is.numeric(partial_auc$range) && length(partial_auc$range) == 2L &&
    all(is.finite(partial_auc$range)) && partial_auc$range[1] >= 0 &&
    partial_auc$range[2] <= 1 && partial_auc$range[1] < partial_auc$range[2]
  if (!isTRUE(ok)) {
    .ig_abort("invalid_input",
              "partial_auc must be list(focus = 'fpr' or 'tpr', range = c(lo, hi)) with 0 <= lo < hi <= 1")
  }
  list(focus = partial_auc$focus, range = as.double(partial_auc$range))
}

# Partial AUC values (raw, normalised, McClish) from ROC points.
.ig_cls_partial_values <- function(fpr, tpr, partial) {
  lo <- partial$range[1]; hi <- partial$range[2]
  if (partial$focus == "fpr") {
    raw <- .ig_cls_segment_area(fpr, tpr, lo, hi)
    min_area <- (hi^2 - lo^2) / 2               # chance diagonal over the range
  } else {
    raw <- .ig_cls_segment_area(tpr, 1 - fpr, lo, hi)
    min_area <- (hi - lo) - (hi^2 - lo^2) / 2   # integral of (1 - t) over [lo, hi]
  }
  max_area <- hi - lo
  mcclish <- (1 + (raw - min_area) / (max_area - min_area)) / 2
  list(raw = raw, normalised = raw / max_area, mcclish = mcclish,
       mcclish_ok = raw >= min_area, min_area = min_area, max_area = max_area)
}

# AP (step sum) and trapezoidal PR area from decision sets.
.ig_cls_pr_values <- function(ds) {
  P <- ds$tp[1] + ds$fn[1]
  if (P == 0) return(list(ap = NA_real_, trap = NA_real_, status = "undefined", reason = "no_events"))
  r <- ds$tp / P
  prec <- ifelse(ds$alerts > 0, ds$tp / ds$alerts, NA_real_)
  k <- seq_len(nrow(ds))[-1L]
  dr <- r[k] - r[k - 1L]
  ap <- sum(ifelse(dr > 0, dr * prec[k], 0))
  p_prev <- c(1, prec[k][-length(k)])   # P_0 := 1 at recall 0
  trap <- sum(ifelse(dr > 0, dr * (prec[k] + p_prev) / 2, 0))
  list(ap = ap, trap = trap, status = "ok", reason = NA_character_)
}

# DeLong interval record for a single AUC.
.ig_cls_delong_interval <- function(st, level, clustered) {
  z <- stats::qnorm((1 - level) / 2, lower.tail = FALSE)
  rec <- list(metric_id = "roc_auc", variant = NA_character_, method = "delong", level = level,
              lower = NA_real_, upper = NA_real_, status = "ok", reason_code = NA_character_,
              note = NA_character_)
  bad <- function(status, code) utils::modifyList(rec, list(status = status, reason_code = code))
  if (st$status != "ok") return(bad("undefined", st$reason))
  if (clustered) return(bad("unsupported", "clustered_units"))
  if (is.na(st$var)) return(bad("undefined", "single_case_or_control"))
  if (st$var <= 0) return(bad("undefined", "delong_degenerate"))
  lo <- st$auc - z * sqrt(st$var)
  hi <- st$auc + z * sqrt(st$var)
  clipped <- lo < 0 || hi > 1
  rec$lower <- max(0, lo)
  rec$upper <- min(1, hi)
  rec$note <- if (clipped) "normal-approximation limits clipped to [0, 1]" else NA_character_
  rec
}

.ig_cls_interval_frame <- function(recs) {
  if (!length(recs)) {
    return(data.frame(metric_id = character(), variant = character(), method = character(),
                      level = double(), lower = double(), upper = double(), status = character(),
                      reason_code = character(), reason = character(), note = character(),
                      stringsAsFactors = FALSE))
  }
  df <- data.frame(
    metric_id = vapply(recs, function(r) r$metric_id, ""),
    variant = vapply(recs, function(r) as.character(r$variant %||% NA_character_), ""),
    method = vapply(recs, function(r) r$method, ""),
    level = vapply(recs, function(r) as.double(r$level), 0),
    lower = vapply(recs, function(r) as.double(r$lower), 0),
    upper = vapply(recs, function(r) as.double(r$upper), 0),
    status = vapply(recs, function(r) r$status, ""),
    reason_code = vapply(recs, function(r) as.character(r$reason_code %||% NA_character_), ""),
    note = vapply(recs, function(r) as.character(r$note %||% NA_character_), ""),
    stringsAsFactors = FALSE)
  df$reason <- .ig_cls_reason_text(df$reason_code)
  df[, c("metric_id", "variant", "method", "level", "lower", "upper", "status", "reason_code",
         "reason", "note")]
}

.ig_cls_summary_row <- function(metric_id, value, status = "ok", reason_code = NA_character_,
                                variant = NA_character_, estimate_kind = "observed",
                                setting = NA_character_) {
  if (status != "ok") value <- NA_real_
  list(metric_id = metric_id, variant = variant, value = as.double(value), status = status,
       reason_code = reason_code, estimate_kind = estimate_kind, setting = setting)
}

.ig_cls_summary_frame <- function(rows) {
  df <- data.frame(
    metric_id = vapply(rows, `[[`, "", "metric_id"),
    variant = vapply(rows, function(r) as.character(r$variant), ""),
    value = vapply(rows, function(r) r$value, 0),
    status = vapply(rows, `[[`, "", "status"),
    reason_code = vapply(rows, function(r) as.character(r$reason_code), ""),
    estimate_kind = vapply(rows, `[[`, "", "estimate_kind"),
    setting = vapply(rows, function(r) as.character(r$setting), ""),
    stringsAsFactors = FALSE)
  df$reason <- .ig_cls_reason_text(df$reason_code)
  df
}

# Curve summaries of one empirical model: summaries table and intervals.
.ig_cls_curve_summaries <- function(score, y, w, direction, ds, partial, level, auc_method,
                                    clustered, estimate_kind) {
  st <- .ig_cls_auc_stats(score, y, w, direction)
  rows <- list(.ig_cls_summary_row("roc_auc", st$auc, st$status, st$reason,
                                   estimate_kind = estimate_kind,
                                   setting = "Mann-Whitney, ties count 1/2"))
  intervals <- list()
  if (identical(auc_method, "delong")) intervals <- list(.ig_cls_delong_interval(st, level, clustered))
  if (!is.null(partial)) {
    setting <- sprintf("%s in [%s, %s]", toupper(partial$focus), format(partial$range[1]),
                       format(partial$range[2]))
    if (st$status != "ok") {
      for (v in c("raw", "normalised", "mcclish")) {
        rows[[length(rows) + 1L]] <- .ig_cls_summary_row("partial_auc", NA, st$status, st$reason,
                                                         variant = v, estimate_kind = estimate_kind,
                                                         setting = setting)
      }
    } else {
      pts <- .ig_cls_roc_points(ds)
      pv <- .ig_cls_partial_values(pts$fpr, pts$tpr, partial)
      rows[[length(rows) + 1L]] <- .ig_cls_summary_row("partial_auc", pv$raw, variant = "raw",
                                                       estimate_kind = estimate_kind, setting = setting)
      rows[[length(rows) + 1L]] <- .ig_cls_summary_row("partial_auc", pv$normalised, variant = "normalised",
                                                       estimate_kind = estimate_kind, setting = setting)
      rows[[length(rows) + 1L]] <- .ig_cls_summary_row(
        "partial_auc", pv$mcclish, if (pv$mcclish_ok) "ok" else "undefined",
        if (pv$mcclish_ok) NA_character_ else "mcclish_below_diagonal",
        variant = "mcclish", estimate_kind = estimate_kind, setting = setting)
    }
  }
  pr <- .ig_cls_pr_values(ds)
  rows[[length(rows) + 1L]] <- .ig_cls_summary_row("average_precision", pr$ap, pr$status, pr$reason,
                                                   estimate_kind = estimate_kind,
                                                   setting = "step sum over decision sets; observed prevalence")
  rows[[length(rows) + 1L]] <- .ig_cls_summary_row("pr_area_trapezoid", pr$trap, pr$status, pr$reason,
                                                   estimate_kind = estimate_kind,
                                                   setting = "trapezoid, precision at recall 0 set to 1; observed prevalence")
  list(summaries = .ig_cls_summary_frame(rows), intervals = intervals, auc_stats = st)
}

# Paired DeLong comparison of two score columns on the same units.
.ig_cls_paired_delong <- function(score_a, score_b, y, w, dir_a, dir_b, level) {
  if (is.null(w)) w <- rep(1, length(y))
  keep <- w > 0
  sa <- score_a[keep]; sb <- score_b[keep]; y <- y[keep]; w <- w[keep]
  st_a <- .ig_cls_auc_stats(sa, y, w, dir_a)
  st_b <- .ig_cls_auc_stats(sb, y, w, dir_b)
  out <- list(auc_a = st_a$auc, auc_b = st_b$auc, difference = NA_real_, se = NA_real_,
              lower = NA_real_, upper = NA_real_, z = NA_real_, p_value = NA_real_,
              status = "ok", reason_code = NA_character_)
  if (st_a$status != "ok") return(utils::modifyList(out, list(status = "undefined", reason_code = st_a$reason)))
  out$difference <- st_a$auc - st_b$auc
  A <- st_a$n_pos; B <- st_a$n_neg
  if (A < 2 || B < 2) return(utils::modifyList(out, list(status = "undefined", reason_code = "single_case_or_control")))
  place <- function(st, s, dir) {
    o <- .ig_cls_orientation(dir) * s
    idx <- match(o, st$groups$u)
    list(v10 = st$groups$v10[idx], v01 = st$groups$v01[idx])
  }
  pa <- place(st_a, sa, dir_a)
  pb <- place(st_b, sb, dir_b)
  ev <- y == 1L
  ne <- y == 0L
  cov_w <- function(x1, x2, m1, m2, ww, n) sum(ww * (x1 - m1) * (x2 - m2)) / (n - 1)
  s10 <- matrix(c(cov_w(pa$v10[ev], pa$v10[ev], st_a$auc, st_a$auc, w[ev], A),
                  cov_w(pa$v10[ev], pb$v10[ev], st_a$auc, st_b$auc, w[ev], A),
                  cov_w(pb$v10[ev], pa$v10[ev], st_b$auc, st_a$auc, w[ev], A),
                  cov_w(pb$v10[ev], pb$v10[ev], st_b$auc, st_b$auc, w[ev], A)), 2)
  s01 <- matrix(c(cov_w(pa$v01[ne], pa$v01[ne], st_a$auc, st_a$auc, w[ne], B),
                  cov_w(pa$v01[ne], pb$v01[ne], st_a$auc, st_b$auc, w[ne], B),
                  cov_w(pb$v01[ne], pa$v01[ne], st_b$auc, st_a$auc, w[ne], B),
                  cov_w(pb$v01[ne], pb$v01[ne], st_b$auc, st_b$auc, w[ne], B)), 2)
  S <- s10 / A + s01 / B
  contrast <- c(1, -1)
  v <- as.double(t(contrast) %*% S %*% contrast)
  diff <- st_a$auc - st_b$auc
  if (!is.finite(v) || v <= 0) {
    return(utils::modifyList(out, list(status = "undefined", reason_code = "delong_degenerate")))
  }
  z <- stats::qnorm((1 - level) / 2, lower.tail = FALSE)
  out$se <- sqrt(v)
  out$lower <- diff - z * sqrt(v)
  out$upper <- diff + z * sqrt(v)
  out$z <- diff / sqrt(v)
  out$p_value <- 2 * stats::pnorm(-abs(out$z))
  out
}

#' ROC AUC and partial AUC of scored models
#'
#' Computes the empirical ROC AUC (ties count 1/2, declared direction, never
#' auto-flipped), an optional partial AUC over a stated FPR or TPR range, and
#' the DeLong confidence interval (own implementation; for a fixed set of
#' scores, not for a selected cutoff). For a theoretical model the AUC is the
#' model parameter and has no sampling interval.
#'
#' @param spec A `classification_spec` (scored data with attached data, or
#'   theory).
#' @param partial_auc Optional `list(focus = "fpr" | "tpr", range = c(lo, hi))`.
#'   pROC's `partial.auc = c(1, 0.8), partial.auc.focus = "specificity"`
#'   corresponds to `list(focus = "fpr", range = c(0, 0.2))`.
#' @param level Confidence level (default: the spec's uncertainty level).
#' @return A data frame with one row per model and summary (`metric_id`,
#'   `variant`, `value`, `status`, `reason_code`, `reason`, `estimate_kind`,
#'   `setting`) and the DeLong interval columns `ci_lower`, `ci_upper`,
#'   `ci_status`, `ci_reason`.
#' @examples
#' d <- data.frame(y = rep(c(1, 0), c(6, 8)),
#'                 s = c(9, 7, 7, 6, 5, 3, 7, 5, 4, 4, 3, 2, 2, 1))
#' spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "s", cutoff = 5)
#' classification_auc(spec, partial_auc = list(focus = "fpr", range = c(0, 0.2)))
#' @family classification
#' @export
classification_auc <- function(spec, partial_auc = NULL, level = NULL) {
  .ig_cls_check_spec(spec)
  partial <- .ig_cls_check_partial(partial_auc)
  level <- level %||% spec$uncertainty$level
  if (identical(spec$input_kind, "counts_2x2")) {
    .ig_abort("unsupported", .ig_cls_reason_text("aggregate_counts_only"))
  }
  if (identical(spec$input_kind, "theoretical_binormal")) {
    cs <- .ig_cls_theory_curve_summaries(spec, partial)
    out <- cbind(model_id = "binormal", cs, ci_lower = NA_real_, ci_upper = NA_real_,
                 ci_status = "not_applicable", ci_reason = .ig_cls_reason_text("theoretical_no_interval"),
                 stringsAsFactors = FALSE)
    return(out[out$metric_id %in% c("roc_auc", "partial_auc"), , drop = FALSE])
  }
  prep <- .ig_cls_spec_data(spec)
  blocks <- lapply(seq_along(spec$models), function(i) {
    m <- spec$models[[i]]
    ds <- .ig_cls_decision_sets(prep$scores[[i]], prep$y, prep$w, m$direction, m$comparator)
    cs <- .ig_cls_curve_summaries(prep$scores[[i]], prep$y, prep$w, m$direction, ds, partial, level,
                                  "delong", !is.null(prep$cluster),
                                  if (is.null(prep$w)) "observed" else "frequency_weighted")
    s <- cs$summaries[cs$summaries$metric_id %in% c("roc_auc", "partial_auc"), , drop = FALSE]
    iv <- cs$intervals[[1]]
    s$ci_lower <- ifelse(s$metric_id == "roc_auc", iv$lower, NA_real_)
    s$ci_upper <- ifelse(s$metric_id == "roc_auc", iv$upper, NA_real_)
    s$ci_status <- ifelse(s$metric_id == "roc_auc", iv$status, "not_computed")
    s$ci_reason <- ifelse(s$metric_id == "roc_auc", .ig_cls_reason_text(iv$reason_code),
                          .ig_cls_reason_text("no_interval_method"))
    cbind(model_id = m$model_id, s, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, blocks)
  rownames(out) <- NULL
  out
}

#' Precision-recall summaries of scored models
#'
#' Returns average precision (step sum over the decision sets) and the
#' trapezoidal PR area (precision at recall 0 set to 1) as two distinct
#' quantities. Both refer to the observed prevalence; with `prevalence`, the
#' average precision is also transported to that prevalence under fixed
#' class-conditional score distributions.
#'
#' @param spec A scored-data `classification_spec` with attached data.
#' @param prevalence Optional scenario prevalence in (0, 1).
#' @return A data frame with `model_id`, `metric_id`, `value`, `status`,
#'   `reason_code`, `reason`, `estimate_kind` and `setting`.
#' @examples
#' d <- data.frame(y = c(1, 0, 1, 1, 0, 0, 1, 0), s = c(.9, .8, .7, .6, .55, .5, .4, .3))
#' spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "s", cutoff = .5)
#' pr_summary(spec, prevalence = 0.1)
#' @family classification
#' @export
pr_summary <- function(spec, prevalence = NULL) {
  .ig_cls_check_spec(spec)
  if (!identical(spec$input_kind, "scored_data")) {
    .ig_abort("unsupported", if (identical(spec$input_kind, "counts_2x2"))
      .ig_cls_reason_text("aggregate_counts_only") else
        "pr_summary() needs scored data; theoretical summaries are part of evaluate_cutoff()")
  }
  if (!is.null(prevalence) && !(is.numeric(prevalence) && length(prevalence) == 1L &&
                                is.finite(prevalence) && prevalence > 0 && prevalence < 1)) {
    .ig_abort("invalid_input", "prevalence must be a single number in (0, 1)")
  }
  prep <- .ig_cls_spec_data(spec)
  blocks <- lapply(seq_along(spec$models), function(i) {
    m <- spec$models[[i]]
    ds <- .ig_cls_decision_sets(prep$scores[[i]], prep$y, prep$w, m$direction, m$comparator)
    pr <- .ig_cls_pr_values(ds)
    kind <- if (is.null(prep$w)) "observed" else "frequency_weighted"
    rows <- list(
      .ig_cls_summary_row("average_precision", pr$ap, pr$status, pr$reason, estimate_kind = kind,
                          setting = "step sum over decision sets; observed prevalence"),
      .ig_cls_summary_row("pr_area_trapezoid", pr$trap, pr$status, pr$reason, estimate_kind = kind,
                          setting = "trapezoid, precision at recall 0 set to 1; observed prevalence"))
    if (!is.null(prevalence)) {
      sc <- .ig_cls_scenario_ap(ds, prevalence)
      rows[[3]] <- .ig_cls_summary_row("average_precision", sc$value, sc$status, sc$reason,
                                       variant = "scenario", estimate_kind = "expected_scenario",
                                       setting = sprintf("step sum at assumed prevalence %s", format(prevalence)))
    }
    cbind(model_id = m$model_id, .ig_cls_summary_frame(rows), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, blocks)
  rownames(out) <- NULL
  out
}

# Average precision transported to prevalence pi (expected counts with the
# observed TPR/FPR of every decision set).
.ig_cls_scenario_ap <- function(ds, pi) {
  P <- ds$tp[1] + ds$fn[1]
  N <- ds$fp[1] + ds$tn[1]
  if (P == 0) return(list(value = NA_real_, status = "undefined", reason = "no_events"))
  if (N == 0) return(list(value = NA_real_, status = "undefined", reason = "no_non_events"))
  exp_ds <- data.frame(tp = pi * ds$tp / P, fp = (1 - pi) * ds$fp / N)
  exp_ds$fn <- pi - exp_ds$tp
  exp_ds$tn <- (1 - pi) - exp_ds$fp
  exp_ds$alerts <- exp_ds$tp + exp_ds$fp
  pr <- .ig_cls_pr_values(exp_ds)
  list(value = pr$ap, status = pr$status, reason = pr$reason)
}

#' Compare scored models on the same units
#'
#' Pairwise comparison of the ROC AUC of models evaluated on the same units
#' (paired design), using the DeLong covariance of the placement values (own
#' implementation). Exclusions of missing values are matched across models by
#' [classification_spec()].
#'
#' @param spec A scored-data `classification_spec` with at least two models
#'   and attached data.
#' @param level Confidence level (default: the spec's uncertainty level).
#' @return A data frame with one row per model pair: `model_a`, `model_b`,
#'   `auc_a`, `auc_b`, `difference`, `se`, `lower`, `upper`, `z`, `p_value`,
#'   `status`, `reason_code`, `reason`, `method`.
#' @examples
#' d <- data.frame(y = rep(c(1, 0), each = 6),
#'                 a = c(8, 7, 7, 5, 4, 2, 6, 3, 3, 2, 1, 1),
#'                 b = c(5, 6, 4, 4, 3, 3, 5, 4, 2, 2, 2, 1))
#' spec <- classification_spec(d, truth = "y", positive_label = 1, scores = c("a", "b"),
#'                             cutoff = 4)
#' compare_classifiers(spec)
#' @family classification
#' @export
compare_classifiers <- function(spec, level = NULL) {
  .ig_cls_check_spec(spec)
  if (!identical(spec$input_kind, "scored_data") || length(spec$models) < 2L) {
    .ig_abort("invalid_input", "compare_classifiers() needs a scored-data specification with at least two models")
  }
  level <- level %||% spec$uncertainty$level
  prep <- .ig_cls_spec_data(spec)
  .ig_cls_comparisons(spec, prep, level)
}

.ig_cls_comparisons <- function(spec, prep, level) {
  k <- length(spec$models)
  rows <- list()
  for (i in seq_len(k - 1L)) for (j in (i + 1L):k) {
    ma <- spec$models[[i]]; mb <- spec$models[[j]]
    r <- .ig_cls_paired_delong(prep$scores[[i]], prep$scores[[j]], prep$y, prep$w,
                               ma$direction, mb$direction, level)
    if (!is.null(prep$cluster) && r$status == "ok") {
      # The point difference stands; the DeLong inference assumes independent units.
      r <- utils::modifyList(r, list(se = NA_real_, lower = NA_real_, upper = NA_real_,
                                     z = NA_real_, p_value = NA_real_,
                                     status = "unsupported", reason_code = "clustered_units"))
    }
    rows[[length(rows) + 1L]] <- data.frame(
      model_a = ma$model_id, model_b = mb$model_id,
      auc_a = r$auc_a, auc_b = r$auc_b, difference = r$difference, se = r$se,
      lower = r$lower, upper = r$upper, z = r$z, p_value = r$p_value, level = level,
      inference_status = r$status, reason_code = r$reason_code,
      reason = .ig_cls_reason_text(r$reason_code), method = "delong_paired",
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
