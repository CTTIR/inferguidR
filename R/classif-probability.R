# Probability views: Brier score, calibration and decision-curve net benefit --------------
#
# These quantities need scores declared as predicted probabilities of the
# positive class (score_scale = "probability"); for raw scores they are
# reported as not applicable. Discrimination says nothing about calibration,
# so calibration is assessed separately:
#
# * Brier score: weighted mean of (p - y)^2.
# * Calibration in the large: mean(p) - mean(y).
# * Calibration intercept: a in logit P(y = 1) = a + logit(p) (logit(p) as
#   offset); calibration slope: b in logit P(y = 1) = a + b logit(p). Both by
#   maximum likelihood (stats::glm.fit, binomial family, frequency weights).
# * Grouped calibration table: units grouped by quantiles of p. The grouping
#   is a display diagnostic; it never enters a decision or a metric.
# * Net benefit at an explicit decision threshold probability pt with the
#   harm weight w = pt / (1 - pt): NB = TP/N - FP/N * w for the rule "alert if
#   p >= pt" (or "> pt" under the model's comparator); references alert all
#   (pi - (1 - pi) w) and alert none (0). A generic score cutoff is not a
#   threshold probability, and net benefit presumes calibrated probabilities.

.ig_cls_check_thresholds <- function(threshold_probability) {
  if (is.null(threshold_probability)) return(NULL)
  ok <- is.numeric(threshold_probability) && length(threshold_probability) >= 1L &&
    all(is.finite(threshold_probability)) && all(threshold_probability > 0 & threshold_probability < 1)
  if (!ok) .ig_abort("invalid_input", "`threshold_probability` must be numeric values in (0, 1)")
  sort(unique(as.double(threshold_probability)))
}

.ig_cls_check_bins <- function(bins) {
  if (!is.numeric(bins) || length(bins) != 1L || !is.finite(bins) || bins < 2 || bins != round(bins)) {
    .ig_abort("invalid_input", "`calibration_bins` must be a whole number >= 2")
  }
  as.double(bins)
}

# Reason why the probability views of a model are not applicable, or NULL
# when they apply: they need probabilities of the positive class, so higher
# scores must mean positive.
.ig_cls_probability_unavailable <- function(score_scale, direction) {
  if (!identical(score_scale, "probability")) return("not_probability_scale")
  if (!identical(direction, "higher_is_positive")) return("probability_direction_reversed")
  NULL
}

# Summary rows (Brier, calibration in the large, intercept, slope).
.ig_cls_probability_rows <- function(p, y, w, score_scale, estimate_kind, direction = "higher_is_positive") {
  ids <- c("brier_score", "calibration_in_the_large", "calibration_intercept", "calibration_slope")
  why_not <- .ig_cls_probability_unavailable(score_scale, direction)
  if (!is.null(why_not)) {
    return(lapply(ids, function(id) .ig_cls_summary_row(id, NA, "not_applicable", why_not,
                                                       estimate_kind = estimate_kind)))
  }
  if (is.null(w)) w <- rep(1, length(y))
  keep <- w > 0
  p <- p[keep]; y <- y[keep]; w <- w[keep]
  sw <- sum(w)
  rows <- list(
    .ig_cls_summary_row("brier_score", sum(w * (p - y)^2) / sw, estimate_kind = estimate_kind,
                        setting = "weighted mean of (p - y)^2"),
    .ig_cls_summary_row("calibration_in_the_large", sum(w * p) / sw - sum(w * y) / sw,
                        estimate_kind = estimate_kind, setting = "mean predicted minus observed event rate")
  )
  fit <- .ig_cls_calibration_fit(p, y, w)
  for (id in c("calibration_intercept", "calibration_slope")) {
    f <- fit[[id]]
    rows[[length(rows) + 1L]] <- .ig_cls_summary_row(id, f$value, f$status, f$reason_code,
                                                     estimate_kind = estimate_kind, setting = f$setting)
  }
  rows
}

# Logistic recalibration fits.
.ig_cls_calibration_fit <- function(p, y, w) {
  bad <- function(code) list(value = NA_real_, status = "undefined", reason_code = code,
                             setting = NA_character_)
  out <- list(calibration_intercept = NULL, calibration_slope = NULL)
  if (sum(w[y == 1L]) == 0) return(list(calibration_intercept = bad("no_events"), calibration_slope = bad("no_events")))
  if (sum(w[y == 0L]) == 0) return(list(calibration_intercept = bad("no_non_events"), calibration_slope = bad("no_non_events")))
  lp <- stats::qlogis(p)
  if (any(!is.finite(lp))) {
    return(list(calibration_intercept = bad("probability_endpoints"),
                calibration_slope = bad("probability_endpoints")))
  }
  fit_one <- function(x, offset) {
    # tight convergence so that the reported coefficients are accurate to
    # far more digits than are ever displayed
    f <- tryCatch(suppressWarnings(stats::glm.fit(x, y, weights = w, offset = offset,
                                                  family = stats::binomial(),
                                                  control = list(epsilon = 1e-12, maxit = 100))),
                  error = function(e) NULL)
    if (is.null(f) || !isTRUE(f$converged) || isTRUE(f$boundary) || any(!is.finite(f$coefficients))) NULL else f
  }
  n <- length(y)
  f_int <- fit_one(matrix(1, n, 1L), lp)
  f_slope <- fit_one(cbind(1, lp), rep(0, n))
  out$calibration_intercept <- if (is.null(f_int)) bad("calibration_fit_failed") else
    list(value = unname(f_int$coefficients[1]), status = "ok", reason_code = NA_character_,
         setting = "logistic model logit P(y = 1) = a + logit(p), maximum likelihood")
  out$calibration_slope <- if (is.null(f_slope)) bad("calibration_fit_failed") else
    list(value = unname(f_slope$coefficients[2]), status = "ok", reason_code = NA_character_,
         setting = "logistic model logit P(y = 1) = a + b logit(p), maximum likelihood")
  out
}

# Grouped calibration table (quantile groups of p; display diagnostic).
.ig_cls_calibration_table <- function(p, y, w, bins) {
  if (is.null(w)) w <- rep(1, length(y))
  keep <- w > 0
  p <- p[keep]; y <- y[keep]; w <- w[keep]
  o <- order(p)
  cw <- cumsum(w[o]) / sum(w)
  grp <- integer(length(p))
  grp[o] <- pmin(bins, ceiling(cw * bins - 1e-12))
  grp[grp < 1L] <- 1L
  # units with equal p stay in one group (the group of their first member)
  key <- match(p, unique(p))
  first <- tapply(grp, key, min)
  grp <- as.integer(first[key])
  ug <- sort(unique(grp))
  rows <- lapply(seq_along(ug), function(j) {
    g <- grp == ug[j]
    list(group = as.double(j), n_units = sum(w[g]), events = sum(w[g & y == 1L]),
         mean_predicted = sum(w[g] * p[g]) / sum(w[g]),
         observed_rate = sum(w[g & y == 1L]) / sum(w[g]),
         min_predicted = min(p[g]), max_predicted = max(p[g]))
  })
  .ig_cls_rows_to_cols(rows)
}

# Net benefit table for one probability score.
.ig_cls_net_benefit_table <- function(p, y, w, direction, comparator, thresholds, estimate_kind) {
  if (is.null(w)) w <- rep(1, length(y))
  n <- sum(w)
  P <- sum(w[y == 1L])
  N <- n - P
  rows <- lapply(thresholds, function(pt) {
    hw <- pt / (1 - pt)
    alert <- if (identical(comparator, ">")) p > pt else p >= pt
    tp <- sum(w[alert & y == 1L])
    fp <- sum(w[alert & y == 0L])
    list(threshold_probability = pt, harm_weight = hw, tp = tp, fp = fp, n_total = n,
         net_benefit = tp / n - fp / n * hw, net_benefit_alert_all = P / n - N / n * hw,
         net_benefit_alert_none = 0, rule = sprintf("alert if p %s %s", if (identical(comparator, ">")) ">" else ">=",
                                                    format(pt)),
         estimate_kind = estimate_kind)
  })
  .ig_cls_rows_to_cols(rows)
}

.ig_cls_probability_model_check <- function(spec, m) {
  if (!identical(spec$input_kind, "scored_data")) {
    .ig_abort("unsupported", "probability views need scored data with predicted probabilities")
  }
  why_not <- .ig_cls_probability_unavailable(m$score_scale, m$direction)
  if (!is.null(why_not)) {
    .ig_abort("unsupported", sprintf("model '%s': %s", m$model_id, .ig_cls_reason_text(why_not)))
  }
}

#' Calibration of predicted probabilities
#'
#' Brier score, calibration in the large, calibration intercept and slope
#' (logistic recalibration) and a grouped calibration table for scored models
#' whose scores are declared as predicted probabilities of the positive class
#' (`score_scale = "probability"`). Discrimination (AUC) does not establish
#' calibration; this function assesses it separately. The groups of the table
#' (quantiles of the predicted probability) are a display diagnostic only.
#'
#' @param spec A scored-data `classification_spec` with attached data.
#' @param bins Number of quantile groups for the table (default 10); units with
#'   equal predicted probability stay in one group.
#' @return A list with `summaries` (one row per model and quantity: `model_id`,
#'   `metric_id`, `value`, `status`, `reason_code`, `reason`, `setting`) and
#'   `table` (grouped observed versus predicted, one block per model).
#' @examples
#' ex <- classification_examples("scored_synthetic")
#' calibration_summary(ex$spec)$summaries
#' @family classification
#' @export
calibration_summary <- function(spec, bins = 10) {
  .ig_cls_check_spec(spec)
  bins <- .ig_cls_check_bins(bins)
  for (m in spec$models) .ig_cls_probability_model_check(spec, m)
  prep <- .ig_cls_spec_data(spec)
  sums <- list()
  tabs <- list()
  for (i in seq_along(spec$models)) {
    m <- spec$models[[i]]
    .ig_cls_probability_model_check(spec, m)
    kind <- if (is.null(prep$w)) "observed" else "frequency_weighted"
    rows <- .ig_cls_probability_rows(prep$scores[[i]], prep$y, prep$w, m$score_scale, kind)
    sums[[i]] <- cbind(model_id = m$model_id, .ig_cls_summary_frame(rows), stringsAsFactors = FALSE)
    tb <- as.data.frame(.ig_cls_calibration_table(prep$scores[[i]], prep$y, prep$w, bins))
    tabs[[i]] <- cbind(model_id = m$model_id, tb, stringsAsFactors = FALSE)
  }
  s <- do.call(rbind, sums)
  t <- do.call(rbind, tabs)
  rownames(s) <- NULL
  rownames(t) <- NULL
  list(summaries = s, table = t)
}

#' Decision-curve net benefit
#'
#' Net benefit of the rule "alert if predicted probability >= pt" at explicit
#' decision threshold probabilities `pt`, with the harm weight
#' \eqn{w = p_t / (1 - p_t)} (the odds at the threshold: one false alert is
#' weighed as \eqn{w} missed true alerts):
#' \deqn{NB = TP/N - (FP/N) w,}
#' together with the references "alert all" (\eqn{\pi - (1-\pi) w}) and
#' "alert none" (0). Scores must be predicted probabilities of the positive
#' class; a generic score cutoff is not a threshold probability. Net benefit
#' presumes calibrated probabilities, see [calibration_summary()].
#'
#' @param spec A scored-data `classification_spec` with attached data and
#'   probability scores.
#' @param threshold_probability Numeric threshold probabilities in (0, 1). No
#'   default: thresholds express a declared harm-benefit trade-off.
#' @return A data frame with `model_id`, `threshold_probability`,
#'   `harm_weight`, `tp`, `fp`, `n_total`, `net_benefit`,
#'   `net_benefit_alert_all`, `net_benefit_alert_none`, `rule` and
#'   `estimate_kind`.
#' @examples
#' ex <- classification_examples("scored_synthetic")
#' net_benefit(ex$spec, threshold_probability = c(0.1, 0.2, 0.3))
#' @family classification
#' @export
net_benefit <- function(spec, threshold_probability) {
  .ig_cls_check_spec(spec)
  if (missing(threshold_probability)) {
    .ig_abort("invalid_input", .ig_cls_reason_text("threshold_probability_missing"))
  }
  thr <- .ig_cls_check_thresholds(threshold_probability)
  for (m in spec$models) .ig_cls_probability_model_check(spec, m)
  prep <- .ig_cls_spec_data(spec)
  out <- list()
  for (i in seq_along(spec$models)) {
    m <- spec$models[[i]]
    .ig_cls_probability_model_check(spec, m)
    kind <- if (is.null(prep$w)) "observed" else "frequency_weighted"
    tb <- as.data.frame(.ig_cls_net_benefit_table(prep$scores[[i]], prep$y, prep$w, m$direction,
                                                  m$comparator, thr, kind), stringsAsFactors = FALSE)
    out[[i]] <- cbind(model_id = m$model_id, tb, stringsAsFactors = FALSE)
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}
