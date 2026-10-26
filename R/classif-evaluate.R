# Evaluation of an operating rule ----------------------------------------------------------
#
# evaluate_cutoff() applies the operating rule of a classification_spec to
# each model and returns a classification_result:
#
# * scored data: exact decision sets from the original scores; the operating
#   point is the decision set selected by the manual cutoff (counts are also
#   recomputed directly from the scores and must agree) or by the selection
#   objective (exhaustive search); curve summaries, intervals, probability
#   views and an optional prevalence scenario;
# * aggregate 2x2 counts: the one fixed operating point; curves, sweeps and
#   selection are unavailable with a reason;
# * theoretical binormal model: expected counts at a real-valued cutoff
#   (analytic tails), the model's AUC, partial AUC and average precision, and
#   a display grid for curves; no sampling intervals.

.ig_cls_model_head <- function(m) {
  list(model_id = m$model_id, name = m$name, version = m$version, score_meaning = m$score_meaning,
       score_scale = m$score_scale, direction = m$direction, comparator = m$comparator,
       tie_rule = m$tie_rule, evaluation_split = m$evaluation_split)
}

.ig_cls_rule_text <- function(cutoff, cutoff_kind, comparator) {
  if (identical(cutoff_kind, "alert_none")) return("alert no unit")
  if (identical(cutoff_kind, "alert_all")) return("alert every unit")
  sprintf("alert if score %s %s", comparator, .ig_format_double(cutoff))
}

.ig_cls_op_point <- function(cutoff, cutoff_kind, counts, count_kind, comparator, set_index = NULL,
                             boundary_score = NULL) {
  tp <- counts[["tp"]]; fn <- counts[["fn"]]; fp <- counts[["fp"]]; tn <- counts[["tn"]]
  list(cutoff = if (is.finite(cutoff)) cutoff else NA_real_, cutoff_kind = cutoff_kind,
       set_index = set_index, boundary_score = boundary_score,
       tp = tp, fn = fn, fp = fp, tn = tn,
       n_alerts = tp + fp, n_non_alerts = fn + tn, n_events = tp + fn, n_non_events = fp + tn,
       n_total = tp + fn + fp + tn, count_kind = count_kind, comparator = comparator,
       rule_text = .ig_cls_rule_text(cutoff, cutoff_kind, comparator))
}

.ig_cls_view <- function(available, reason_code = NA_character_) {
  list(available = available, reason_code = if (available) NA_character_ else reason_code)
}

.ig_cls_interval_cols <- function(rows) {
  if (!length(rows)) return(NULL)
  cols <- .ig_cls_rows_to_cols(rows)
  cols$reason <- .ig_cls_reason_text(cols$reason_code)
  cols
}

.ig_cls_summary_cols <- function(rows) {
  if (!length(rows)) return(NULL)
  .ig_cls_cols(.ig_cls_summary_frame(rows))
}

# A model record without an operating point (selection infeasible or an
# internal consistency failure).
.ig_cls_failed_model <- function(rec, code, reason, selection = NULL) {
  rec$status <- .ig_cls_status(code, reason)
  rec$selection <- selection
  rec$operating_point <- NULL
  rec$metrics <- NULL
  rec$scenario <- NULL
  if (!is.null(rec$intervals)) {
    keep <- rec$intervals$conditioning != "fixed_cutoff"
    rec$intervals <- if (any(keep)) lapply(rec$intervals, function(v) v[keep]) else NULL
  }
  if (!is.null(rec$probability)) rec$probability$net_benefit <- NULL
  rec
}

# Scored data -------------------------------------------------------------------------------

.ig_cls_eval_scored <- function(spec, settings) {
  prep <- .ig_cls_spec_data(spec)
  models <- lapply(seq_along(spec$models), function(i) .ig_cls_eval_scored_model(spec, prep, i, settings))
  comparisons <- NULL
  if (length(spec$models) >= 2L) {
    comparisons <- .ig_cls_cols(.ig_cls_comparisons(spec, prep, spec$uncertainty$level))
  }
  diagnostics <- list(n_input_rows = prep$n_input, n_excluded_rows = prep$n_excluded,
                      excluded_by = .ig_named_or_empty(prep$excluded_by))
  list(models = models, comparisons = comparisons, diagnostics = diagnostics)
}

.ig_cls_eval_scored_model <- function(spec, prep, i, settings) {
  m <- spec$models[[i]]
  op <- spec$operating_rule
  unc <- spec$uncertainty
  score <- prep$scores[[i]]
  count_kind <- if (is.null(prep$w)) "observed" else "frequency_weighted"
  clustered <- !is.null(prep$cluster)
  case_control <- identical(spec$prevalence$sampling, "case_control")
  scen <- spec$prevalence
  n_all <- prep$n_events + prep$n_non_events
  ds <- .ig_cls_decision_sets(score, prep$y, prep$w, m$direction, m$comparator)
  rec <- .ig_cls_model_head(m)
  rec$count_kind <- count_kind

  # Cutoff-independent summaries and views.
  cs <- .ig_cls_curve_summaries(score, prep$y, prep$w, m$direction, ds, settings$partial_auc,
                                unc$level, unc$auc_method, clustered, count_kind)
  sum_rows <- lapply(seq_len(nrow(cs$summaries)), function(r) as.list(cs$summaries[r, c(
    "metric_id", "variant", "value", "status", "reason_code", "estimate_kind", "setting")]))
  if (identical(scen$kind, "scenario")) {
    sap <- .ig_cls_scenario_ap(ds, scen$value)
    sum_rows[[length(sum_rows) + 1L]] <- .ig_cls_summary_row(
      "average_precision", sap$value, sap$status, sap$reason, variant = "scenario",
      estimate_kind = "expected_scenario",
      setting = sprintf("step sum at %s prevalence %s; fixed class-conditional distributions", scen$source,
                        format(scen$value)))
  }
  sum_rows <- c(sum_rows, .ig_cls_probability_rows(score, prep$y, prep$w, m$score_scale, count_kind,
                                                    m$direction))
  rec$summaries <- .ig_cls_summary_cols(sum_rows)
  auc_rows <- lapply(cs$intervals, function(iv) {
    .ig_cls_interval_row(iv$metric_id, iv$method, iv$level, iv$lower, iv$upper, iv$status,
                         iv$reason_code, conditioning = "fixed_scores", note = iv$note)
  })
  # Probability views need probabilities of the positive class (higher means
  # positive), as in calibration_summary() and net_benefit().
  why_not_prob <- .ig_cls_probability_unavailable(m$score_scale, m$direction)
  is_prob <- is.null(why_not_prob)
  rec$views <- list(
    roc = .ig_cls_view(prep$n_events > 0 && prep$n_non_events > 0,
                       if (prep$n_events == 0) "no_events" else "no_non_events"),
    pr = .ig_cls_view(prep$n_events > 0, "no_events"),
    sweep = .ig_cls_view(TRUE),
    score_distribution = .ig_cls_view(TRUE),
    calibration = .ig_cls_view(is_prob, why_not_prob),
    net_benefit = .ig_cls_view(is_prob && !is.null(settings$threshold_probability),
                               if (!is_prob) why_not_prob else "threshold_probability_missing"))
  rec$curve <- if (isTRUE(settings$include_curves)) .ig_cls_cols(ds[, names(.ig_cls_table_columns$curve)]) else NULL
  if (is_prob) {
    rec$probability <- list(
      calibration = .ig_cls_calibration_table(score, prep$y, prep$w, settings$calibration_bins),
      net_benefit = if (is.null(settings$threshold_probability)) NULL else
        .ig_cls_net_benefit_table(score, prep$y, prep$w, m$direction, m$comparator,
                                  settings$threshold_probability, count_kind))
  }
  rec$intervals <- .ig_cls_interval_cols(auc_rows)

  # Selection basis: observed cells or expected cells at the scenario prevalence.
  # Without events or without non-events the scenario transports nothing: the
  # basis cells are NA, every metric on them is undefined ("counts_undefined")
  # and a selection is infeasible for that stated reason.
  basis_reason <- NULL
  basis <- if (identical(op$prevalence_basis, "scenario")) {
    tr <- .ig_cls_transport(ds$tp, ds$fn, ds$fp, ds$tn, scen$value, n_all)
    if (tr$status[1] != "ok") basis_reason <- tr$reason[1]
    tr$counts
  } else ds[, c("tp", "fn", "fp", "tn")]

  direct <- function(cutoff) .ig_cls_counts_at(score, prep$y, prep$w, m$direction, m$comparator, cutoff)
  agrees <- function(cnt, row) {
    !is.na(row) && cnt[["set_index"]] == ds$set_index[row] &&
      all(unlist(ds[row, c("tp", "fn", "fp", "tn")]) == cnt[c("tp", "fn", "fp", "tn")])
  }
  if (identical(op$selection, "manual")) {
    cnt <- direct(op$score_cutoff)
    row <- match(cnt[["set_index"]], ds$set_index)
    if (!agrees(cnt, row)) {
      return(.ig_cls_failed_model(rec, "numerical_failure",
                                  "direct counts at the cutoff disagree with the enumerated decision sets"))
    }
    sel <- .ig_cls_manual_selection(op, basis[row, , drop = FALSE])
    cutoff <- op$score_cutoff
    kind <- "manual"
  } else {
    sel <- .ig_cls_select_sets(ds, basis, op)
    if (!identical(sel$status, "ok")) {
      if (!is.null(basis_reason)) {
        sel$reason_code <- basis_reason
        sel$reason <- .ig_cls_selection_reason(basis_reason)
      }
      sel$row <- NULL
      return(.ig_cls_failed_model(rec, "infeasible", sel$reason, sel))
    }
    row <- sel$row
    cutoff <- ds$cutoff[row]
    kind <- ds$cutoff_kind[row]
    if (!is.na(cutoff) && !agrees(direct(cutoff), row)) {
      return(.ig_cls_failed_model(rec, "numerical_failure",
                                  "direct counts at the selected cutoff disagree with the enumerated decision sets"))
    }
  }
  sel$row <- NULL
  rec$selection <- sel
  counts <- c(tp = ds$tp[row], fn = ds$fn[row], fp = ds$fp[row], tn = ds$tn[row])
  rec$operating_point <- .ig_cls_op_point(if (is.na(cutoff)) NA_real_ else cutoff, kind, counts, count_kind,
                                          m$comparator, set_index = ds$set_index[row],
                                          boundary_score = ds$boundary_score[row])
  rec$metrics <- .ig_cls_result_metric_cols(counts, count_kind, list(case_control = case_control))
  note <- if (!identical(op$selection, "manual")) {
    "conditional on the cutoff selected on these data; the selection uncertainty is not included"
  } else NA_character_
  iv <- .ig_cls_proportion_intervals(counts, unc$proportion_methods, unc$level,
                                     list(count_kind = count_kind, clustered = clustered,
                                          case_control = case_control, note = note))
  if (!is.null(unc$bootstrap_metrics)) {
    iv <- c(iv, .ig_cls_bootstrap(score, prep$y, prep$w, prep$cluster, m$direction, m$comparator,
                                  cutoff, kind, unc, settings$partial_auc, m$score_scale, note))
  }
  rec$intervals <- .ig_cls_interval_cols(c(iv, auc_rows))
  if (identical(scen$kind, "scenario")) {
    rec$scenario <- .ig_cls_scenario_record(counts, scen$value, n_all, scen$source)
  }
  rec$status <- .ig_cls_status("ok")
  rec
}

# Aggregate counts --------------------------------------------------------------------------

.ig_cls_eval_counts <- function(spec, settings) {
  m <- spec$models[[1]]
  op <- spec$operating_rule
  unc <- spec$uncertainty
  counts <- .ig_cls_check_counts(spec$counts, "observed")
  case_control <- identical(spec$prevalence$sampling, "case_control")
  rec <- .ig_cls_model_head(m)
  rec$count_kind <- "observed"
  kind <- if (identical(op$cutoff_source, "binary_test_convention")) "binary_test_convention" else
    "fixed_operating_point"
  rec$operating_point <- .ig_cls_op_point(op$score_cutoff, kind, counts, "observed", m$comparator)
  rec$selection <- .ig_cls_manual_selection(op, as.data.frame(as.list(counts)))
  rec$metrics <- .ig_cls_result_metric_cols(counts, "observed", list(case_control = case_control))
  iv <- .ig_cls_proportion_intervals(counts, unc$proportion_methods, unc$level,
                                     list(count_kind = "observed", clustered = FALSE,
                                          case_control = case_control))
  if (identical(unc$auc_method, "delong")) {
    iv[[length(iv) + 1L]] <- .ig_cls_interval_row("roc_auc", "delong", unc$level, status = "not_applicable",
                                                  reason_code = "aggregate_counts_only",
                                                  conditioning = "fixed_scores")
  }
  rec$intervals <- .ig_cls_interval_cols(iv)
  na_rows <- list(.ig_cls_summary_row("roc_auc", NA, "not_applicable", "aggregate_counts_only"))
  if (!is.null(settings$partial_auc)) {
    na_rows[[2]] <- .ig_cls_summary_row("partial_auc", NA, "not_applicable", "aggregate_counts_only",
                                        variant = "raw")
  }
  na_rows <- c(na_rows, list(
    .ig_cls_summary_row("average_precision", NA, "not_applicable", "aggregate_counts_only"),
    .ig_cls_summary_row("pr_area_trapezoid", NA, "not_applicable", "aggregate_counts_only")))
  rec$summaries <- .ig_cls_summary_cols(na_rows)
  if (identical(spec$prevalence$kind, "scenario")) {
    rec$scenario <- .ig_cls_scenario_record(counts, spec$prevalence$value, sum(counts),
                                            spec$prevalence$source)
  }
  no <- .ig_cls_view(FALSE, "aggregate_counts_only")
  rec$views <- list(roc = no, pr = no, sweep = no, score_distribution = no,
                    calibration = no, net_benefit = no)
  rec$status <- .ig_cls_status("ok")
  list(models = list(rec), comparisons = NULL, diagnostics = list())
}

# Theoretical binormal model ------------------------------------------------------------------

.ig_cls_eval_theory <- function(spec, settings) {
  m <- spec$models[[1]]
  op <- spec$operating_rule
  unc <- spec$uncertainty
  th <- spec$theory
  pi <- spec$prevalence$value
  rec <- .ig_cls_model_head(m)
  rec$count_kind <- "expected"
  rec$summaries <- .ig_cls_cols(.ig_cls_theory_curve_summaries(spec, settings$partial_auc))
  ok_view <- .ig_cls_view(TRUE)
  rec$views <- list(roc = ok_view, pr = .ig_cls_view(pi > 0, "no_events"), sweep = ok_view,
                    score_distribution = ok_view,
                    calibration = .ig_cls_view(FALSE, "not_probability_scale"),
                    net_benefit = .ig_cls_view(FALSE, "not_probability_scale"))
  auc_iv <- if (identical(unc$auc_method, "delong")) list(.ig_cls_interval_row(
    "roc_auc", "delong", unc$level, status = "not_applicable", reason_code = "theoretical_no_interval",
    conditioning = "fixed_scores")) else list()
  rec$intervals <- .ig_cls_interval_cols(auc_iv)
  if (identical(op$selection, "manual")) {
    cutoff <- op$score_cutoff
    ec <- .ig_cls_theory_counts(cutoff, th, pi)
    sel <- .ig_cls_manual_selection(op, ec)
  } else {
    sel <- .ig_cls_theory_select(th, pi, op)
    if (!identical(sel$status, "ok")) {
      rec$curve <- if (isTRUE(settings$include_curves)) .ig_cls_theory_curve(spec, NULL) else NULL
      return(list(models = list(.ig_cls_failed_model(rec, "infeasible", sel$reason, sel)),
                  comparisons = NULL, diagnostics = list()))
    }
    cutoff <- sel$cutoff
    ec <- .ig_cls_theory_counts(cutoff, th, pi)
  }
  kind <- if (is.finite(cutoff)) "real" else if (cutoff > 0) "alert_none" else "alert_all"
  counts <- c(tp = ec$tp, fn = ec$fn, fp = ec$fp, tn = ec$tn)
  rec$selection <- sel
  rec$operating_point <- .ig_cls_op_point(cutoff, kind, counts, "expected", m$comparator)
  rec$operating_point$cutoff_scale <- "score standardised to the non-event distribution (z)"
  rec$metrics <- .ig_cls_result_metric_cols(counts, "expected_theoretical")
  iv <- .ig_cls_proportion_intervals(counts, unc$proportion_methods, unc$level,
                                     list(count_kind = "expected", expected_reason = "theoretical_no_interval"))
  rec$intervals <- .ig_cls_interval_cols(c(iv, auc_iv))
  rec$curve <- if (isTRUE(settings$include_curves)) .ig_cls_theory_curve(spec, cutoff) else NULL
  rec$status <- .ig_cls_status("ok")
  list(models = list(rec), comparisons = NULL, diagnostics = list())
}

# Display grid of the theoretical model with the limiting rules and the
# operating point (exact, not snapped to the grid).
.ig_cls_theory_curve <- function(spec, cutoff) {
  g <- .ig_cls_theory_grid(spec, 401)
  g <- g[, names(.ig_cls_table_columns$curve)]
  th <- spec$theory
  pi <- spec$prevalence$value
  extra <- list()
  for (z in c(Inf, -Inf, if (!is.null(cutoff) && is.finite(cutoff)) cutoff)) {
    ec <- .ig_cls_theory_counts(z, th, pi)
    extra[[length(extra) + 1L]] <- data.frame(
      set_index = NA_real_, cutoff = if (is.finite(z)) z else NA_real_,
      cutoff_kind = if (!is.finite(z)) (if (z > 0) "alert_none" else "alert_all") else "operating_point",
      boundary_score = NA_real_, alerts = ec$tp + ec$fp, tp = ec$tp, fn = ec$fn, fp = ec$fp, tn = ec$tn,
      stringsAsFactors = FALSE)
  }
  out <- rbind(g, do.call(rbind, extra))
  # order from alert none (highest cutoff) to alert all, as for decision sets
  key <- ifelse(out$cutoff_kind == "alert_none", Inf, ifelse(out$cutoff_kind == "alert_all", -Inf, out$cutoff))
  out <- out[order(-key), ]
  .ig_cls_cols(out)
}

# Transfer of a cutoff chosen on other data -----------------------------------------------------

# Check only supplied identities; disjoint identifiers do not prove design independence.
.ig_cls_transfer_overlap <- function(spec, src) {
  source_spec <- attr(src, .ig_cls_spec_attr, exact = TRUE)
  target_data <- .ig_cls_spec_data(spec)
  if (!identical(.ig_cls_fingerprint(target_data), spec$data_fingerprint))
    .ig_abort("invalid_input", "the attached evaluation data no longer match their fingerprint")
  source_data <- NULL
  if (!is.null(source_spec)) {
    .ig_cls_check_spec(source_spec)
    if (!identical(source_spec$spec_hash, src$spec_hash) ||
        !identical(source_spec$data_fingerprint, src$data_fingerprint))
      .ig_abort("invalid_input", "the attached tuning specification does not match the source result")
    source_data <- attr(source_spec, .ig_cls_data_attr, exact = TRUE)
    if (!is.null(source_data) && !identical(.ig_cls_fingerprint(source_data), source_spec$data_fingerprint))
      .ig_abort("invalid_input", "the attached tuning data no longer match their fingerprint")
  }
  check <- function(field, label) {
    a <- source_data[[field]]; b <- target_data[[field]]
    if (is.null(a) || is.null(b)) return("not_available")
    if (any(a %in% b)) .ig_abort("invalid_input", paste("tuning and evaluation data share", label,
      "- this cannot be treated as held-out evaluation"))
    "no_overlap_detected"
  }
  list(unit_identity = check("id", "unit identifiers"),
       cluster_identity = check("cluster", "cluster identifiers"),
       scope = "Only provided identifiers in a common namespace are compared. Absence of overlap does not prove independence or exclude model-training leakage.")
}

.ig_cls_transfer <- function(spec, src, model_id) {
  if (!identical(spec$input_kind, "scored_data")) {
    .ig_abort("invalid_input", "a cutoff from another result can only be applied to scored data")
  }
  if (!identical(src$input_kind, "scored_data")) {
    .ig_abort("invalid_input", "the source result must come from scored data (a theoretical cutoff is on a standardised scale)")
  }
  if (!identical(.ig_cls_status_code(src$status), "ok")) {
    .ig_abort("invalid_input", "the source result has no operating point (its status is not ok)")
  }
  ids <- vapply(src$models, function(m) m$model_id, "")
  if (is.null(model_id)) {
    if (length(ids) != 1L) .ig_abort("invalid_input", "the source result has several models; give `model_id`")
    model_id <- ids
  }
  j <- match(model_id, ids)
  if (!.ig_is_string(model_id) || is.na(j)) .ig_abort("invalid_input", "`model_id` is not a model of the source result")
  sm <- src$models[[j]]
  op <- sm$operating_point
  if (is.null(op$cutoff) || is.na(op$cutoff)) {
    .ig_abort("invalid_input", sprintf(
      "the source rule (%s) has no numeric cutoff that would carry over to other data", op$rule_text))
  }
  for (m in spec$models) {
    if (!identical(m$comparator, sm$comparator) || !identical(m$direction, sm$direction) ||
        !identical(m$score_scale, sm$score_scale)) {
      .ig_abort("invalid_input", sprintf(
        "model '%s' uses a different rule (comparator, direction or score scale) than the source model '%s'",
        m$model_id, model_id))
    }
  }
  if (identical(src$data_fingerprint, spec$data_fingerprint)) {
    .ig_abort("invalid_input", "the source result was computed on the same data; this would be an in-sample evaluation, not a held-out one")
  }
  overlap_check <- .ig_cls_transfer_overlap(spec, src)
  sel_data <- if (identical(src$evaluation$selection_data, "predeclared_rule") &&
                  identical(src$evaluation$selection, "manual")) "predeclared_rule" else "held_out_tuning"
  rule <- cutoff_rule(cutoff = op$cutoff, selection_data = sel_data)
  record <- list(source_result_hash = .ig_hash_value(.ig_unclass_frozen(src)),
                 source_spec_hash = src$spec_hash, source_model_id = model_id,
                 source_procedure = sm$selection$procedure %||% NA_character_,
                 source_evaluation_label = src$evaluation$label, cutoff = op$cutoff, overlap_check = overlap_check)
  list(spec = update(spec, rule = rule), record = record)
}

#' Evaluate an operating rule (cutoff) of a classification specification
#'
#' Applies the operating rule of `spec` to every model and returns an
#' immutable `classification_result`. For scored data the counts come from the
#' original scores under the declared comparator (never from histogram bins);
#' every distinct decision set is enumerated, so ROC/PR summaries and sweeps are
#' exact. A selection rule in `spec` (for example `"max_youden"`) is carried
#' out on the evaluated data and labelled in-sample; see [select_cutoff()].
#'
#' The result holds, per model: the operating point (cutoff, comparator,
#' TP/FN/FP/TN and totals with their count kind), every operating-point metric
#' of the [metric_dictionary()] with status, reason, denominator and unit,
#' Wilson and Clopper-Pearson intervals (conditional on the fixed cutoff), the
#' ROC AUC with its own DeLong interval, optional partial AUC, average
#' precision and the trapezoidal PR area as distinct quantities, optional
#' bootstrap intervals, calibration and net benefit for probability scores, a
#' prevalence scenario when one is declared, the decision sets (or the
#' theoretical display grid) and the availability of each view. It also states
#' how the cutoff was chosen relative to the data and lists the limitations
#' that apply.
#'
#' @param spec A `classification_spec`.
#' @param cutoff Optional override: a number (a manual cutoff for all models),
#'   or a `classification_result` computed on *other* (tuning) data whose
#'   selected cutoff is then evaluated here as a held-out rule.
#'   Shared supplied unit or cluster identifiers are refused when both live
#'   specifications retain their data. Imported source results cannot supply
#'   this identity check. Identity absence is recorded; disjoint identifiers
#'   never prove independence or rule out model-training leakage. Identifiers
#'   must share the same meaning across splits; none are exported by this check.
#' @param model_id With a result as `cutoff`: the model whose cutoff is used
#'   (needed when the result has several models).
#' @param partial_auc Optional partial AUC range,
#'   `list(focus = "fpr" | "tpr", range = c(lo, hi))`.
#' @param threshold_probability Optional decision threshold probabilities for
#'   net benefit (probability scores only); see [net_benefit()].
#' @param calibration_bins Number of quantile groups of the calibration table.
#' @param include_curves Store the decision sets (or the theoretical display
#'   grid) in the result (needed for sweeps and curve plots).
#' @return A `classification_result`; see [classification_table()] for its
#'   tables. Its status is `"ok"`, or `"infeasible"` when a selection rule has
#'   no admissible operating point (no number is reported then). The
#'   specification it was evaluated with (which differs from `spec` after a
#'   cutoff override or a transferred rule) is returned by
#'   [classification_evaluated_spec()].
#' @examples
#' # Fixed 2x2 table
#' res <- evaluate_cutoff(classification_spec(counts = c(tp = 40, fn = 10, fp = 20, tn = 130)))
#' res
#'
#' # Scored data, manual cutoff, partial AUC over FPR 0-0.2
#' ex <- classification_examples("scored_synthetic")
#' r2 <- evaluate_cutoff(ex$spec, cutoff = 0.3, partial_auc = list(focus = "fpr", range = c(0, 0.2)))
#' classification_table(r2, "summaries")
#'
#' # Rule chosen on tuning data, evaluated on test data
#' tune <- classification_examples("scored_synthetic", part = "tuning", models = "model_a")$spec
#' test <- classification_examples("scored_synthetic", part = "test", models = "model_a")$spec
#' chosen <- select_cutoff(tune, "max_youden")
#' held_out <- evaluate_cutoff(test, cutoff = chosen)
#' held_out$evaluation$label
#' @seealso [select_cutoff()], [classification_table()], [cutoff_sweep()],
#'   [prevalence_scenario()]
#' @concept classification
#' @seealso
#' * [classification_spec()]
#' * [plot_classification()]
#' @export
evaluate_cutoff <- function(spec, cutoff = NULL, model_id = NULL, partial_auc = NULL,
                            threshold_probability = NULL, calibration_bins = 10,
                            include_curves = TRUE) {
  .ig_cls_check_spec(spec)
  partial <- .ig_cls_check_partial(partial_auc)
  thr <- .ig_cls_check_thresholds(threshold_probability)
  bins <- .ig_cls_check_bins(calibration_bins)
  if (!.ig_is_flag(include_curves)) .ig_abort("invalid_input", "`include_curves` must be TRUE or FALSE")
  if (!is.null(thr) && !identical(spec$input_kind, "scored_data")) {
    .ig_abort("invalid_input", "net benefit needs scored data with predicted probabilities")
  }
  transfer <- NULL
  if (!is.null(cutoff)) {
    if (inherits(cutoff, "classification_result")) {
      tr <- .ig_cls_transfer(spec, cutoff, model_id)
      spec <- tr$spec
      transfer <- tr$record
    } else if (is.numeric(cutoff) && length(cutoff) == 1L && is.finite(cutoff)) {
      if (!is.null(model_id)) .ig_abort("invalid_input", "`model_id` is only used with a result as `cutoff`")
      if (identical(spec$input_kind, "counts_2x2")) {
        .ig_abort("invalid_input", "the cutoff of aggregate 2x2 counts is fixed by the counts")
      }
      spec <- update(spec, cutoff = cutoff)
    } else {
      .ig_abort("invalid_input", "`cutoff` must be a single finite number or a classification_result")
    }
  } else if (!is.null(model_id)) {
    .ig_abort("invalid_input", "`model_id` is only used with a result as `cutoff`")
  }
  settings <- list(partial_auc = partial, threshold_probability = thr, calibration_bins = bins,
                   include_curves = include_curves)
  out <- switch(spec$input_kind,
                scored_data = .ig_cls_eval_scored(spec, settings),
                counts_2x2 = .ig_cls_eval_counts(spec, settings),
                theoretical_binormal = .ig_cls_eval_theory(spec, settings))
  args <- character()
  if (!is.null(transfer)) args <- c(args, "cutoff = <tuning result>")
  if (!is.null(partial)) args <- c(args, paste0("partial_auc = ", deparse1(partial)))
  if (!is.null(thr)) args <- c(args, paste0("threshold_probability = ", deparse1(thr)))
  if (bins != 10) args <- c(args, paste0("calibration_bins = ", bins))
  if (!include_curves) args <- c(args, "include_curves = FALSE")
  r_call <- sprintf("inferguidR::evaluate_cutoff(spec%s)  # spec: classification_spec %s",
                    if (length(args)) paste0(", ", paste(args, collapse = ", ")) else "",
                    spec$spec_hash)
  res <- .ig_new_classification_result(spec, out$models, out$comparisons, settings, transfer,
                                       out$diagnostics, r_call)
  # The evaluated specification (which differs from the input when a cutoff
  # override, a transferred rule or select_cutoff() changed the rule) travels
  # with the result as an attribute; it is never serialised.
  attr(res, .ig_cls_spec_attr) <- spec
  res
}

.ig_cls_spec_attr <- "ig_cls_evaluated_spec"

#' Specification behind a classification result
#'
#' [evaluate_cutoff()] with a cutoff override or a transferred rule, and
#' [select_cutoff()], evaluate an updated specification. This accessor returns
#' that evaluated specification (with its data attached), which is the one an
#' [analysis_envelope()] or a report must pair with the result.
#'
#' The specification travels with the result in the R session only. A result
#' read back from JSON has none; its envelope already holds the matching
#' specification.
#'
#' @param x A `classification_result` computed in this session.
#' @return The evaluated `classification_spec`.
#' @examples
#' ex <- classification_examples("scored_synthetic")
#' res <- select_cutoff(ex$spec, "max_youden")
#' spec <- classification_evaluated_spec(res)
#' spec$operating_rule$selection
#' analysis_envelope(spec, res)
#' @concept classification
#' @seealso
#' * [classification_spec()]
#' * [evaluate_cutoff()]
#' * [plot_classification()]
#' @export
classification_evaluated_spec <- function(x) {
  if (!inherits(x, "classification_result")) {
    .ig_abort("invalid_input", "`x` must be a classification_result (see evaluate_cutoff())")
  }
  spec <- attr(x, .ig_cls_spec_attr, exact = TRUE)
  if (is.null(spec)) {
    .ig_abort("invalid_input", paste(
      "this result carries no evaluated specification (results read from JSON do not);",
      "use the specification stored in its analysis envelope"))
  }
  if (!identical(spec$spec_hash, x$spec_hash)) {
    .ig_abort("invalid_input", "the attached specification does not match the result (spec_hash differs)")
  }
  spec
}
