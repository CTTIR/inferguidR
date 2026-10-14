# Classification results ------------------------------------------------------------------
#
# A classification_result is an immutable plain-data record (lists, atomic
# vectors and column lists; no data frames, no classed members), so that it
# serialises to canonical JSON and comes back from it with the same hash. It
# refers to its specification by spec_hash and is paired with it in an
# analysis_envelope. Statuses are stored as plain lists {code, reason,
# details}. Tables are stored as column lists and returned as data frames by
# classification_table(); reason texts are stored in English and translated
# by the accessor.
#
# Guard rails enforced by the constructor:
#
# * a model whose status is not "ok" carries no operating point, no
#   operating-point metrics, no scenario and no fixed-cutoff interval;
# * in every metric, summary and interval table a value is finite exactly
#   when the row status is "ok", and every other row has a reason code;
# * G*Power validation is never claimed (reference status "not_applicable").

.ig_cls_result_schema_version <- "inferguidR.classification_result/0.1.0-draft"

# Plain status record.
.ig_cls_status <- function(code, reason = "") {
  .ig_status_to_list(ig_status(code, reason))
}

.ig_cls_status_code <- function(st) {
  if (is.null(st)) return(NA_character_)
  if (is_ig_status(st)) return(st$code)
  if (is.list(st) && .ig_is_string(st$code)) return(st$code)
  NA_character_
}

# A list of row lists (same names) as a column list.
.ig_cls_rows_to_cols <- function(rows) {
  if (!length(rows)) return(NULL)
  nms <- names(rows[[1]])
  out <- lapply(nms, function(nm) {
    vals <- lapply(rows, function(r) {
      v <- r[[nm]]
      if (is.null(v) || length(v) == 0L) NA else v
    })
    unlist(vals, use.names = FALSE)
  })
  names(out) <- nms
  out
}

# A data frame as a column list (factors to character, integers to double).
.ig_cls_cols <- function(df) {
  if (is.null(df) || !nrow(df)) return(NULL)
  out <- lapply(as.list(df), function(col) {
    if (is.factor(col)) col <- as.character(col)
    if (is.integer(col)) col <- as.double(col)
    unname(col)
  })
  names(out) <- names(df)
  out
}

# Column order and storage type of the stored tables.
.ig_cls_table_columns <- list(
  metrics = c(metric_id = "character", value = "double", status = "character",
              reason_code = "character", reason = "character", denominator = "character",
              denominator_value = "double", unit = "character", estimate_kind = "character",
              prevalence_dependent = "logical", note = "character"),
  intervals = c(metric_id = "character", variant = "character", method = "character",
                level = "double", lower = "double", upper = "double", status = "character",
                reason_code = "character", reason = "character", conditioning = "character",
                note = "character"),
  summaries = c(metric_id = "character", variant = "character", value = "double",
                status = "character", reason_code = "character", reason = "character",
                estimate_kind = "character", setting = "character"),
  curve = c(set_index = "double", cutoff = "double", cutoff_kind = "character",
            boundary_score = "double", alerts = "double", tp = "double", fn = "double",
            fp = "double", tn = "double"),
  bounds = c(constraint = "character", source = "character", metric_id = "character",
             sense = "character", bound = "double", n_satisfied = "double", n_points = "double",
             best_value = "double", best_status = "character", satisfied_at_choice = "logical"),
  tied = c(set_index = "double", cutoff = "double", cutoff_kind = "character", alerts = "double",
           objective_value = "double"),
  calibration = c(group = "double", n_units = "double", events = "double",
                  mean_predicted = "double", observed_rate = "double", min_predicted = "double",
                  max_predicted = "double"),
  net_benefit = c(threshold_probability = "double", harm_weight = "double", tp = "double",
                  fp = "double", n_total = "double", net_benefit = "double",
                  net_benefit_alert_all = "double", net_benefit_alert_none = "double",
                  rule = "character", estimate_kind = "character"),
  operating_point = c(cutoff = "double", cutoff_kind = "character", set_index = "double",
                      boundary_score = "double", tp = "double", fn = "double", fp = "double",
                      tn = "double", n_alerts = "double", n_non_alerts = "double",
                      n_events = "double", n_non_events = "double", n_total = "double",
                      count_kind = "character", comparator = "character", rule_text = "character",
                      cutoff_scale = "character"),
  selection = c(procedure = "character", status = "character", reason_code = "character",
                reason = "character", objective_metric = "character",
                objective_direction = "character", objective_value = "double",
                prevalence_basis = "character", selection_data = "character",
                tie_break = "character", n_tied = "double", n_candidates = "double",
                n_defined = "double", n_feasible = "double", n_admissible = "double",
                constraints_met = "logical", method = "character"),
  comparisons = c(model_a = "character", model_b = "character", auc_a = "double",
                  auc_b = "double", difference = "double", se = "double", lower = "double",
                  upper = "double", z = "double", p_value = "double", level = "double",
                  inference_status = "character", reason_code = "character",
                  reason = "character", method = "character"),
  limitations = c(code = "character", text = "character")
)

# Column list (possibly rebuilt from JSON) as a typed data frame.
.ig_cls_frame <- function(cols, table, lang = "en") {
  spec <- .ig_cls_table_columns[[table]]
  if (is.null(cols)) {
    if (is.null(spec)) return(data.frame())
    return(as.data.frame(lapply(spec, function(tp) vector(tp, 0L)), stringsAsFactors = FALSE))
  }
  # Canonical JSON writes a length-1 column as a scalar and NA as null, so a
  # one-row table read back from JSON can hold NULL for a missing cell.
  n <- max(1L, vapply(cols, length, 1L))
  conv <- function(v, tp) {
    if (is.null(v)) v <- NA
    if (is.list(v)) v <- unlist(lapply(v, function(e) if (is.null(e)) NA else e))
    if (length(v) == 1L && n > 1L) v <- rep(v, n)
    switch(tp, character = as.character(v), double = as.double(v), logical = as.logical(v), v)
  }
  nms <- if (is.null(spec)) names(cols) else c(intersect(names(spec), names(cols)),
                                               setdiff(names(cols), names(spec)))
  out <- lapply(nms, function(nm) {
    tp <- if (!is.null(spec) && nm %in% names(spec)) spec[[nm]] else NA
    v <- cols[[nm]]
    if (is.na(tp)) {
      if (is.null(v)) v <- NA
      if (is.list(v)) v <- unlist(lapply(v, function(e) if (is.null(e)) NA else e))
      if (length(v) == 1L && n > 1L) v <- rep(v, n)
      v
    } else conv(v, tp)
  })
  names(out) <- nms
  df <- as.data.frame(out, stringsAsFactors = FALSE)
  if (identical(lang, "de")) df <- .ig_cls_translate_frame(df, table)
  df
}

# Translate reason texts and denominator labels of a table to German.
.ig_cls_translate_frame <- function(df, table) {
  if ("reason_code" %in% names(df) && "reason" %in% names(df)) {
    df$reason <- .ig_cls_reason_text(df$reason_code, "de")
  }
  if (table == "metrics" && "denominator" %in% names(df)) {
    dict <- .ig_cls_dictionary()$metrics
    df$denominator <- vapply(df$metric_id, function(id) {
      d <- dict[[id]]$denominator
      if (is.null(d)) NA_character_ else sprintf("%s (%s)", d$expression, d$de %||% d$en)
    }, character(1), USE.NAMES = FALSE)
  }
  df
}

# Metric table of one operating point as a column list, with the
# prevalence-dependence flag and notes.
.ig_cls_result_metric_cols <- function(counts, estimate_kind, context = list()) {
  tab <- .ig_cls_metric_table(counts, estimate_kind, "en")
  dict <- .ig_cls_dictionary()$metrics
  tab$prevalence_dependent <- vapply(tab$metric_id, function(id) {
    !identical(dict[[id]]$prevalence_dependence, "invariant_given_fixed_class_conditional_distributions")
  }, logical(1), USE.NAMES = FALSE)
  tab$note <- NA_character_
  if (isTRUE(context$case_control)) {
    dep <- tab$prevalence_dependent & tab$unit != "count"
    tab$note[dep] <- "sample-specific: outcome-selected sampling fixes the sample prevalence; see a prevalence scenario for population values"
  }
  if (estimate_kind %in% c("expected_theoretical", "expected_scenario")) {
    cnt <- tab$unit == "count"
    tab$note[cnt] <- "expected count, not observed units"
  }
  .ig_cls_cols(tab)
}

# Check a stored table: finite values exactly for ok rows, reasons otherwise.
.ig_cls_check_table <- function(cols, value_cols, where) {
  if (is.null(cols)) return(invisible(TRUE))
  st <- cols$status
  for (vc in value_cols) {
    v <- cols[[vc]]
    if (is.null(v)) next
    bad_ok <- st == "ok" & !is.finite(v)
    bad_num <- st != "ok" & !is.na(v)
    if (any(bad_ok) || any(bad_num)) {
      stop(sprintf("internal error: %s: values in '%s' do not match their statuses", where, vc), call. = FALSE)
    }
  }
  if (any(st != "ok" & is.na(cols$reason_code))) {
    stop(sprintf("internal error: %s: a row that is not ok has no reason code", where), call. = FALSE)
  }
  invisible(TRUE)
}

.ig_cls_guard_model <- function(m) {
  code <- .ig_cls_status_code(m$status)
  .ig_cls_check_table(m$metrics, "value", paste0(m$model_id, "$metrics"))
  .ig_cls_check_table(m$summaries, "value", paste0(m$model_id, "$summaries"))
  .ig_cls_check_table(m$intervals, c("lower", "upper"), paste0(m$model_id, "$intervals"))
  if (!is.null(m$scenario$metrics)) .ig_cls_check_table(m$scenario$metrics, "value", paste0(m$model_id, "$scenario"))
  if (!identical(code, "ok")) {
    leak <- c(.ig_finite_paths(m$operating_point, "operating_point"),
              .ig_finite_paths(m$metrics, "metrics"),
              .ig_finite_paths(m$scenario, "scenario"),
              .ig_finite_paths(m$probability$net_benefit, "net_benefit"))
    if (!is.null(m$intervals) && any(m$intervals$conditioning == "fixed_cutoff" & m$intervals$status == "ok")) {
      leak <- c(leak, "intervals")
    }
    if (length(leak)) {
      stop(sprintf("internal error: model '%s' with status '%s' reports operating-point numbers (%s)",
                   m$model_id, code, paste(leak, collapse = ", ")), call. = FALSE)
    }
  }
  invisible(TRUE)
}

# Validation state of the classifier components. "independently_checked"
# means compared with validation/independent/classifier/reference.R, which
# shares no numerical routine with the engine (decision D18); the run and its
# match rules are recorded in validation/independent/classifier/summary.md,
# and the brute-force core of the comparison is repeated by
# tests/testthat/test-classif-independent.R. A second, cross-language check
# (validation/independent/classifier/crosscheck/: Python with exact rational
# arithmetic and mpmath, no shared code or library) covers the same
# components. pROC and yardstick are comparison packages only
# (tests/testthat/test-classif-compare.R). The bootstrap has no independent
# check, only property tests, and stays "development".
.ig_cls_validation_record <- function() {
  list(
    gpower_reference_status = "not_applicable",
    note = paste(
      "G*Power is not a reference for classifier metrics, ROC/PR summaries, calibration or",
      "prevalence scenarios; no G*Power validation is claimed for any classification result."),
    independent_checker = "validation/independent/classifier/check_classifier.R",
    checker_summary = "validation/independent/classifier/summary.md",
    cross_check_summary = "validation/independent/classifier/crosscheck/summary.md",
    comparison_packages = c("pROC", "yardstick"),
    components = list(
      confusion_counts = "independently_checked",
      operating_point_metrics = "independently_checked",
      roc_auc = "independently_checked",
      partial_auc = "independently_checked",
      average_precision = "independently_checked",
      pr_area_trapezoid = "independently_checked",
      delong = "independently_checked",
      delong_paired = "independently_checked",
      wilson = "independently_checked",
      clopper_pearson = "independently_checked",
      binormal_theory = "independently_checked",
      prevalence_scenario = "independently_checked",
      ppv_gap_max_prevalence = "independently_checked",
      cutoff_selection = "independently_checked",
      calibration = "independently_checked",
      net_benefit = "independently_checked",
      bootstrap = "development"
    )
  )
}

.ig_cls_limitation_texts <- c(
  theoretical_model = "Values are expected values under an assumed binormal model: no data, no sampling uncertainty, not an empirical performance estimate.",
  aggregate_counts_only = "Only aggregate 2x2 counts are available: no ROC or precision-recall curve, no cutoff sweep and no cutoff selection.",
  in_sample_selection = "The cutoff was selected on the evaluated data. Its performance is descriptive and optimistic for these data; it is not validated operational performance.",
  same_data_manual_cutoff = "The manual cutoff is treated as chosen with access to these data (default); declare a predeclared rule or held-out tuning if that is true.",
  split_independence_unverified = "Held-out status is a design declaration. Supplied identities are checked for overlap when available; absence of overlap does not prove independence or rule out model-training leakage. Imported source results cannot provide their private identities.",
  apparent_performance = "The model is declared as developed on the evaluated data (evaluation split 'same_data'): all metrics, including the AUC, are apparent and optimistic.",
  resampling_split = "The evaluation split is declared as resampling; a selected cutoff is only honestly evaluated when the selection is repeated inside every resampling loop, which this result does not do.",
  fixed_cutoff_intervals = "Intervals of operating-point metrics are conditional on the fixed cutoff; for a data-driven cutoff they do not include the selection uncertainty.",
  post_selection_not_computed = "No post-selection (optimism-corrected or nested-resampling) interval is computed.",
  case_control_sampling = "Outcome-selected (case-control) sampling: the sample prevalence is fixed by design, so PPV, NPV, accuracy and other prevalence-dependent values are sample-specific; use a prevalence scenario for population values.",
  prevalence_scenario = "Scenario values assume fixed class-conditional score distributions; they are not observed measurements and not a claim about real population shift.",
  clustered_units = "Clustered units: binomial and DeLong intervals assume independent units and are not reported; use the cluster bootstrap.",
  frequency_weights = "Frequency weights are treated as counts of replicated independent units.",
  single_class = "Only one actual class is present: rates and summaries that need both classes are undefined.",
  net_benefit_calibration = "Net benefit presumes calibrated predicted probabilities; check the calibration summary.",
  missing_excluded = "Rows with missing values were excluded for all models alike."
)

.ig_cls_evaluation_record <- function(spec, transfer) {
  op <- spec$operating_rule
  kind <- spec$input_kind
  splits <- unique(vapply(spec$models, function(m) m$evaluation_split, ""))
  label <- if (kind == "theoretical_binormal") "theoretical_model" else
    if (kind == "counts_2x2") "fixed_operating_point" else
      if (!identical(op$selection, "manual")) "in_sample_selection" else
        switch(op$selection_data, same_data = "same_data_manual_cutoff",
               held_out_tuning = "held_out_rule", predeclared_rule = "predeclared_rule")
  statement <- switch(label,
    theoretical_model = if (identical(op$selection, "manual"))
      "Expected values of the assumed binormal model at a stated cutoff; no outcome data are involved." else
        "Optimum of the assumed binormal model; no outcome data are involved, so this is a property of the model, not an estimate.",
    fixed_operating_point = switch(op$selection_data,
      predeclared_rule = "Counts of a rule fixed before these data were seen.",
      held_out_tuning = "Counts of a rule chosen on separate tuning data.",
      "Counts of one operating point; the cutoff behind them may have been chosen on the same data, so treat the values as descriptive."),
    in_sample_selection = "The cutoff was chosen on the evaluated data; the values describe these data and are not validated operational performance.",
    same_data_manual_cutoff = "A manual cutoff evaluated on data that may have informed it; treat the values as descriptive.",
    held_out_rule = if (any(splits %in% c("held_out_test", "external")))
      "A rule declared to be chosen on separate tuning data, evaluated on declared independent test or external data. These declarations do not establish independence." else
        "A rule declared to be chosen on separate tuning data, evaluated on other data. These declarations do not establish independence.",
    predeclared_rule = "A rule fixed before these data were seen, evaluated on them.")
  list(label = label, statement = statement, selection = op$selection,
       selection_data = op$selection_data, evaluation_split = splits, transfer = transfer)
}

.ig_cls_limitations <- function(spec, models, evaluation) {
  op <- spec$operating_rule
  codes <- character()
  kind <- spec$input_kind
  if (kind == "theoretical_binormal") codes <- c(codes, "theoretical_model")
  if (kind == "counts_2x2") codes <- c(codes, "aggregate_counts_only")
  if (evaluation$label == "held_out_rule") codes <- c(codes, "split_independence_unverified")
  if (evaluation$label == "in_sample_selection") codes <- c(codes, "in_sample_selection")
  if (evaluation$label == "same_data_manual_cutoff") codes <- c(codes, "same_data_manual_cutoff")
  if (kind != "theoretical_binormal" && "same_data" %in% evaluation$evaluation_split) codes <- c(codes, "apparent_performance")
  if ("resampling" %in% evaluation$evaluation_split) codes <- c(codes, "resampling_split")
  data_driven <- kind == "scored_data" && (!identical(op$selection, "manual") ||
                                             identical(op$selection_data, "same_data"))
  if (kind != "theoretical_binormal") codes <- c(codes, "fixed_cutoff_intervals")
  if (data_driven) codes <- c(codes, "post_selection_not_computed")
  if (identical(spec$prevalence$sampling, "case_control")) codes <- c(codes, "case_control_sampling")
  if (identical(spec$prevalence$kind, "scenario") && kind != "theoretical_binormal") codes <- c(codes, "prevalence_scenario")
  m1 <- spec$models[[1]]
  if (!is.null(m1$cluster_column)) codes <- c(codes, "clustered_units")
  if (!is.null(m1$weight_column)) codes <- c(codes, "frequency_weights")
  single <- any(vapply(models, function(m) {
    op <- m$operating_point
    !is.null(op) && (op$n_events == 0 || op$n_non_events == 0)
  }, logical(1)))
  if (single) codes <- c(codes, "single_class")
  if (any(vapply(models, function(m) !is.null(m$probability$net_benefit), logical(1)))) {
    codes <- c(codes, "net_benefit_calibration")
  }
  if (!is.null(spec$models[[1]]$missing_value_policy) && identical(m1$missing_value_policy, "exclude")) {
    codes <- c(codes, "missing_excluded")
  }
  codes <- unique(codes)
  list(code = codes, text = unname(.ig_cls_limitation_texts[codes]))
}

# Internal constructor.
.ig_new_classification_result <- function(spec, models, comparisons = NULL, settings = list(),
                                          transfer = NULL, diagnostics = list(), r_call = NULL) {
  if (!inherits(spec, "classification_spec")) stop("`spec` must be a classification_spec", call. = FALSE)
  for (m in models) .ig_cls_guard_model(m)
  codes <- vapply(models, function(m) .ig_cls_status_code(m$status), "")
  status <- if (all(codes == "ok")) .ig_cls_status("ok") else {
    bad <- which(codes != "ok")
    first <- codes[bad[1]]
    .ig_cls_status(first, paste(vapply(models[bad], function(m) {
      sprintf("model '%s': %s", m$model_id, m$status$reason)
    }, ""), collapse = "; "))
  }
  evaluation <- .ig_cls_evaluation_record(spec, transfer)
  x <- list(
    schema_version = .ig_cls_result_schema_version,
    object_type = "classification_result",
    status = status,
    spec_hash = spec$spec_hash,
    data_fingerprint = spec$data_fingerprint,
    input_kind = spec$input_kind,
    evaluation = evaluation,
    prevalence = spec$prevalence,
    theory = if (is.null(spec$theory)) NULL else
      c(spec$theory[c("auc", "sigma", "n_total", "model")], list(d = .ig_cls_theory_d(spec$theory))),
    uncertainty = spec$uncertainty,
    settings = .ig_named_or_empty(settings),
    models = unname(models),
    comparisons = comparisons,
    limitations = .ig_cls_limitations(spec, models, evaluation),
    diagnostics = .ig_named_or_empty(diagnostics),
    validation = .ig_cls_validation_record(),
    r_call = r_call %||% "inferguidR::evaluate_cutoff(spec)"
  )
  x <- .ig_doubles(.ig_sort_names(x))
  structure(x[order(names(x), method = "radix")], class = c("classification_result", "ig_frozen"))
}

#' Tables of a classification result
#'
#' Returns one table of a `classification_result` as a data frame, with one
#' block per model where the table is per model. Reason texts and denominator
#' labels are translated with `lang`.
#'
#' @param x A `classification_result` (also one rebuilt from JSON).
#' @param table One of `"metrics"` (all operating-point metrics with status,
#'   reason, denominator, unit and estimate kind), `"intervals"`,
#'   `"summaries"` (ROC AUC, partial AUC, average precision, trapezoidal PR
#'   area, Brier score and calibration), `"operating_point"`, `"curve"`
#'   (decision sets or theoretical display grid), `"selection"`, `"bounds"`
#'   (constraint feasibility), `"tied"` (tied operating points),
#'   `"scenario"` (prevalence-scenario metrics), `"calibration"`,
#'   `"net_benefit"`, `"comparisons"` or `"limitations"`.
#' @param model_id Optional model id(s) to keep.
#' @param lang `"en"` or `"de"`.
#' @return A data frame.
#' @examples
#' res <- evaluate_cutoff(classification_spec(counts = c(tp = 40, fn = 10, fp = 20, tn = 130)))
#' classification_table(res, "intervals")[1:4, c("metric_id", "method", "lower", "upper")]
#' @family classification
#' @export
classification_table <- function(x, table = c("metrics", "intervals", "summaries", "operating_point",
                                              "curve", "selection", "bounds", "tied", "scenario",
                                              "calibration", "net_benefit", "comparisons",
                                              "limitations"),
                                 model_id = NULL, lang = c("en", "de")) {
  if (!inherits(x, "classification_result")) {
    .ig_abort("invalid_input", "`x` must be a classification_result (see evaluate_cutoff())")
  }
  table <- match.arg(table)
  lang <- match.arg(lang)
  if (table == "comparisons") return(.ig_cls_frame(x$comparisons, "comparisons", lang))
  if (table == "limitations") return(.ig_cls_frame(x$limitations, "limitations", lang))
  blocks <- list()
  for (m in x$models) {
    if (!is.null(model_id) && !m$model_id %in% model_id) next
    cols <- switch(table,
      metrics = m$metrics, intervals = m$intervals, summaries = m$summaries, curve = m$curve,
      operating_point = m$operating_point, bounds = m$selection$bounds, tied = m$selection$tied,
      scenario = m$scenario$metrics, calibration = m$probability$calibration,
      net_benefit = m$probability$net_benefit,
      selection = if (is.null(m$selection)) NULL else .ig_cls_selection_row(m$selection))
    if (is.null(cols)) next
    tab_kind <- switch(table, scenario = "metrics", operating_point = "operating_point",
                       selection = "selection", table)
    df <- .ig_cls_frame(cols, tab_kind, lang)
    if (table == "scenario") {
      df <- cbind(scenario_prevalence = m$scenario$prevalence, scenario_n_total = m$scenario$n_total,
                  scenario_source = m$scenario$source, df, stringsAsFactors = FALSE)
    }
    blocks[[length(blocks) + 1L]] <- cbind(model_id = m$model_id, df, stringsAsFactors = FALSE)
  }
  if (!length(blocks)) {
    empty <- .ig_cls_frame(NULL, table, lang)
    return(cbind(model_id = character(0), empty, stringsAsFactors = FALSE)[0, , drop = FALSE])
  }
  out <- do.call(rbind, blocks)
  rownames(out) <- NULL
  out
}

# Scalar fields of a selection record as one row.
.ig_cls_selection_row <- function(sel) {
  keep <- c("procedure", "status", "reason_code", "objective_metric", "objective_direction",
            "objective_value", "prevalence_basis", "selection_data", "tie_break", "n_tied",
            "n_candidates", "n_defined", "n_feasible", "n_admissible", "constraints_met", "method")
  out <- list()
  for (k in keep) {
    v <- sel[[k]]
    out[[k]] <- if (is.null(v) || length(v) != 1L) NA else v
  }
  out$reason <- if (is.na(out$reason_code)) NA_character_ else
    sel$reason %||% .ig_cls_reason_text(out$reason_code)
  out
}

#' @export
format.classification_result <- function(x, ...) {
  st <- x$status
  lines <- c(sprintf("<classification_result> %s", x$input_kind),
             sprintf("  status:     %s", if (identical(st$code, "ok")) "ok" else paste0(st$code, ": ", st$reason)),
             sprintf("  evaluation: %s", x$evaluation$statement))
  fmt <- function(v) if (is.null(v) || is.na(v)) "NA" else format(signif(v, 4))
  for (m in x$models) {
    lines <- c(lines, sprintf("  model %s (%s)", m$model_id, .ig_cls_status_code(m$status)))
    op <- m$operating_point
    if (!is.null(op)) {
      lines <- c(lines,
                 sprintf("    rule:   %s", op$rule_text),
                 sprintf("    counts: TP %s, FN %s, FP %s, TN %s (%s)", fmt(op$tp), fmt(op$fn), fmt(op$fp),
                         fmt(op$tn), op$count_kind))
      mt <- .ig_cls_frame(m$metrics, "metrics")
      show <- mt[mt$metric_id %in% c("tpr", "fpr", "ppv", "npv"), ]
      txt <- vapply(seq_len(nrow(show)), function(i) {
        if (show$status[i] == "ok") sprintf("%s %s", toupper(show$metric_id[i]), format(signif(show$value[i], 4))) else
          sprintf("%s %s (%s)", toupper(show$metric_id[i]), show$status[i], show$reason_code[i])
      }, "")
      lines <- c(lines, paste0("    ", paste(txt, collapse = "; ")))
    } else if (!is.null(m$selection) && !identical(m$selection$status, "ok")) {
      lines <- c(lines, sprintf("    selection: %s", m$selection$reason %||% m$selection$reason_code))
    }
    sm <- .ig_cls_frame(m$summaries, "summaries")
    auc <- sm[sm$metric_id == "roc_auc", , drop = FALSE]
    if (nrow(auc)) {
      lines <- c(lines, sprintf("    ROC AUC: %s", if (auc$status[1] == "ok") format(signif(auc$value[1], 4)) else
        paste(auc$status[1], auc$reason_code[1])))
    }
  }
  if (length(x$limitations$code)) {
    lines <- c(lines, paste0("  limitation: ", x$limitations$code))
  }
  c(lines, sprintf("  spec hash:  %s", x$spec_hash))
}

#' @export
print.classification_result <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}
