# Portable reports contain explicit aggregates, never attached classifier rows.
.ig_report_plain <- function(x) {
  if (is.data.frame(x)) return(lapply(seq_len(nrow(x)), function(i) .ig_report_plain(as.list(x[i, , drop = FALSE]))))
  if (is.list(x)) return(lapply(x, .ig_report_plain))
  if (is.numeric(x)) {
    one <- function(v) if (is.na(v)) NULL else if (is.infinite(v))
      list(value = NULL, status = if (v > 0) "positive_infinity" else "negative_infinity") else v
    if (length(x) == 1L) return(one(x))
    return(lapply(x, one))
  }
  if (is.atomic(x)) { x[is.na(x)] <- NA; return(unname(x)) }
  x
}

.ig_report_model <- function(m) {
  keep <- c("model_id", "name", "version", "score_meaning", "score_scale", "direction",
    "comparator", "tie_rule", "evaluation_split", "count_kind", "status", "operating_point",
    "metrics", "intervals", "summaries", "selection", "scenario")
  out <- m[intersect(keep, names(m))]
  out$operating_point$boundary_score <- NULL
  out$operating_point$set_index <- NULL
  out
}

#' Freeze a portable analysis report
#'
#' Creates an immutable report document from matching specifications and
#' results. All export formats use this same frozen content. The portable
#' document omits attached classifier rows and score-derived curve tables.
#' Aggregate-cell mosaics and prevalence scenarios retain their plotted values;
#' undefined conditional rates remain explicitly unavailable.
#' Selected ROC and precision-recall points use only retained aggregate cells;
#' undefined denominators remain unavailable. Prevalence-scenario precision is
#' labelled separately and is not a new observed measurement. These selected
#' points do not reconstruct full curves or imply an AUC from one cutoff.
#' Original private scored data must be reattached for recomputation.
#' User-supplied titles and metadata are included as literal text.
#' Observed Boschloo reports retain independent-group counts and numerical
#' optimisation bounds; those bounds are not sampling confidence intervals.
#'
#' Comparison scope accepts separate analysis envelopes or multiple classifier
#' models in one completed envelope. Paired ROC-AUC interval figures use only
#' retained comparison aggregates, with marginal intervals and no multiplicity
#' adjustment. Missing intervals are not replaced or clipped.
#'
#' @param envelopes An analysis_envelope or a nonempty list of envelopes.
#' @param recommendation Optional advisor recommendation.
#' @param report_kind One of single_analysis, model_cutoff, comparison,
#'   session or incomplete_settings.
#' @param title Report title; defaults to an analysis report title.
#' @param locale Label language, en or de.
#' @param display_timezone Display timezone recorded with the UTC timestamp.
#' @param sweeps Optional list of entries with analysis (one-based envelope index) and sweep (a sweep_design() result). At most20sweeps and1000rows per sweep are retained; requested values must be finite. Each sweep must match its completed frozen base. Failed rows require report_kind = "incomplete_settings". Specifications, results, failures and individually evaluated points are retained without interpolation. Documents with sweeps use schema0.2; existing schema0.1 documents remain readable.
#' @param resources Optional list of entries with analysis (one-based index) and resource_id. Currently inferguidR_boschloo_minlike is available only for the exact minlike observed test or power variant. Frozen titles, scope, file hashes and source references use a versioned bundled catalog. No resource code is executed and no website is requested. Attached resources use report schema0.3.
#' @return A report_document inheriting from report_snapshot, suitable for
#'   export_report() and safe ig_to_json()/ig_from_json() round trips.
#' @export
create_report_snapshot <- function(envelopes, recommendation = NULL,
                                   report_kind = "single_analysis", title = NULL,
                                   locale = "en", display_timezone = "UTC", sweeps = list(), resources = list()) {
  if (inherits(envelopes, "analysis_envelope")) envelopes <- list(envelopes)
  if (!is.list(envelopes) || !length(envelopes) ||
      !all(vapply(envelopes, inherits, logical(1), "analysis_envelope")))
    .ig_abort("invalid_input", "envelopes must contain analysis_envelope objects")
  envelopes <- lapply(envelopes, function(e) {
    rebuilt <- analysis_envelope(e$spec, e$result)
    if (!identical(rebuilt$analysis_hash, e$analysis_hash))
      .ig_abort("invalid_input", "the envelope was modified after creation")
    rebuilt
  })
  if (is.null(title)) title <- if (locale == "de") "Analysebericht" else "Analysis report"
  if (!.ig_is_string(display_timezone) || !display_timezone %in% c("UTC", OlsonNames()))
    .ig_abort("invalid_input", "display_timezone must be UTC or an Olson timezone")
  frozen <- .ig_new_report_snapshot(envelopes, recommendation, report_kind, title,
    locale, display_timezone)
  specs <- results <- metrics <- uncertainty <- visualisations <- validation <- reproduction <- list()
  for (i in seq_along(envelopes)) {
    e <- envelopes[[i]]; r <- e$result
    specs[[i]] <- list(analysis_kind = e$analysis_kind, analysis_hash = e$analysis_hash,
      specification = if (e$analysis_kind == "power") .ig_spec_to_list(e$spec) else .ig_unclass_frozen(e$spec))
    if (e$analysis_kind == "power") {
      results[[i]] <- if (is.null(r)) list(status = list(code = "not_computed", reason = "Settings only")) else .ig_result_to_list(r)
      metrics[[i]] <- list(analysis = i, achieved = r$achieved, solved = r$solved, quantities = r$quantities)
      uncertainty[[i]] <- list(analysis = i, status = "not_applicable", reason = "Prospective power calculation; no sampling interval for the design inputs.")
      validation[[i]] <- list(analysis = i, status = "research_preview", reason = "Complete numerical release certification remains pending.", evidence = r$validation)
      reproduction[[i]] <- list(analysis = i, r_code = if (is.null(r)) .ig_spec_call(e$spec) else r$r_call,
        private_data_required = FALSE)
      if (!is.null(r) && r$status$code == "ok") {
        pd <- power_plot_data(r)
        visualisations[[length(visualisations) + 1L]] <- list(analysis = i, plot_type = "distribution",
          status = pd$status, reason = pd$reason, axes = list(x = "test statistic", y = pd$ordinate),
          fixed_parameters = e$spec$inputs, critical = pd$critical, data = pd$data,
          caption = paste(c(pd$notes, "Displayed densities or probability masses describe the fixed design."), collapse = " "))
      }
    } else if (e$analysis_kind == "observed_test") {
      if (!e$spec$variant_id %in% c("boschloo_minlike", "boschloo_central", "boschloo_lower", "boschloo_upper") ||
          (!is.null(r) && !identical(r$method, "boschloo")))
        .ig_abort("unsupported", "Only Boschloo observed-test reports are implemented")
      results[[i]] <- if (is.null(r)) list(status = list(code = "not_computed", reason = "Settings only")) else
        .ig_unclass_frozen(r)[setdiff(names(r), "spec")]
      metrics[[i]] <- list(analysis = i, p_value = r$p_value, p_upper = r$p_upper,
        alpha = r$alpha, decision = r$decision, variant = e$spec$options$variant)
      uncertainty[[i]] <- list(analysis = i, status = "not_computed",
        reason = "No sampling confidence interval is computed. p_value and p_upper bound numerical optimisation error, not sampling uncertainty.",
        numerical_bounds = list(lower = r$p_value, upper = r$p_upper),
        evidence = r$sides)
      validation[[i]] <- list(analysis = i, status = "research_preview",
        reason = "Complete numerical release certification remains pending.",
        computation_status = r$status, settings = e$spec$options)
      reproduction[[i]] <- list(analysis = i, r_code = r$r_call,
        private_data_required = FALSE, note = "Independent binomial groups; counts and event orientation are retained.")
      d <- e$spec$inputs
      visualisations[[length(visualisations) + 1L]] <- list(analysis = i,
        plot_type = "two_group_counts", status = "ok",
        axes = list(x = "independent group", y = "observed count"),
        data = list(group_labels = d$group_labels, event_label = d$event_label,
          event = c(d$x1, d$x2), non_event = c(d$n1 - d$x1, d$n2 - d$x2)),
        caption = "Observed event and non-event counts in independent groups. These are not classifier confusion cells.")
    } else {
      safe <- if (is.null(r)) list(status = list(code = "not_computed", reason = "Settings only")) else
        r[intersect(c("status", "input_kind", "prevalence", "evaluation", "theory", "settings",
                      "diagnostics", "limitations", "comparisons", "validation", "data_fingerprint"), names(r))]
      if (!is.null(r)) safe$models <- lapply(r$models, .ig_report_model)
      results[[i]] <- safe
      metrics[[i]] <- list(analysis = i, models = lapply(safe$models, function(m)
        list(model_id = m$model_id, metrics = m$metrics, summaries = m$summaries, scenario = m$scenario)))
      uncertainty[[i]] <- list(analysis = i, settings = e$spec$uncertainty,
        models = lapply(safe$models, function(m) list(model_id = m$model_id, intervals = m$intervals)))
      validation[[i]] <- list(analysis = i, evidence = r$validation)
      reproduction[[i]] <- list(analysis = i, r_code = r$r_call,
        private_data_required = identical(e$spec$input_kind, "scored_data"),
        note = "Reattach the original private input using its fingerprint before recomputing scored-data results.")
      for (m in safe$models) if (!is.null(m$operating_point))
        visualisations[[length(visualisations) + 1L]] <- list(analysis = i, model_id = m$model_id,
          plot_type = "confusion_counts", status = "ok", axes = list(x = "cell", y = m$count_kind),
          data = m$operating_point[c("tp", "fn", "fp", "tn")],
          caption = "Selected-cutoff cells. AUC summarises the full curve; prevalence-scenario values are expected, not observed measurements.")
      for (m in safe$models) visualisations <- c(visualisations,
        .ig_report_classifier_aggregate_figures(m, i, locale))
      for (m in safe$models) visualisations <- c(visualisations,
        .ig_report_classifier_operating_figures(m, i, locale))
      visualisations <- c(visualisations, .ig_report_classifier_comparison_figures(safe, i, locale))
    }
  }
  portable_sweeps <- .ig_report_sweeps(sweeps,envelopes)
  sweeps_complete <- all(vapply(portable_sweeps,function(sw) all(vapply(sw$rows,function(row) identical(row$status,"ok"),TRUE)),TRUE))
  if(!sweeps_complete && report_kind!="incomplete_settings") .ig_abort("invalid_input","Sweeps contain failed rows; use report_kind = incomplete_settings to retain them explicitly.")
  visualisations <- c(visualisations,.ig_report_sweep_figures(portable_sweeps,locale))
  portable_resources <- .ig_report_resources(resources, specs, locale)
  doc <- list(schema_version = if(length(resources)) "inferguidR.report_document/0.3.0" else if(length(portable_sweeps)) "inferguidR.report_document/0.2.0" else "inferguidR.report_document/0.1.0", object_type = "report_document",
    snapshot_id = frozen$snapshot_id, analysis_hash = frozen$analysis_hash, created_at = frozen$created_at,
    display_timezone = display_timezone, locale = locale, report_kind = report_kind,
    completeness = if(sweeps_complete) frozen$completeness else "incomplete", title = title,
    recommendation = if (is.null(recommendation)) NULL else .ig_recommendation_to_list(recommendation),
    specifications = specs, results = results, metrics = metrics, uncertainty = uncertainty,
    comparisons = lapply(results, function(r) r$comparisons), sweeps = portable_sweeps,
    visualisations = visualisations, validation = validation, provenance = frozen$runtime,
    reproducibility = reproduction,
    resources = portable_resources,
    privacy = list(raw_rows_included = FALSE, score_curve_tables_included = FALSE,
      note = "Explicit titles, model metadata and selected cutoffs are included. Review these before sharing."))
  doc <- .ig_json_parse(.ig_canonical_json(.ig_report_plain(doc)))
  doc$snapshot_hash <- .ig_hash_value(doc)
  structure(doc, class = c("report_document", "report_snapshot", "ig_frozen"))
}

.ig_report_document_to_list <- function(x) unclass(x)

.ig_report_document_from_list <- function(x) {
  if (!identical(x$object_type, "report_document") ||
      (!.ig_is_string(x$schema_version) || !x$schema_version %in% c("inferguidR.report_document/0.1.0","inferguidR.report_document/0.2.0","inferguidR.report_document/0.3.0")))
    .ig_abort("invalid_input", "unsupported report document version")
  required <- c("schema_version", "object_type", "snapshot_id", "snapshot_hash", "analysis_hash",
    "created_at", "display_timezone", "locale", "report_kind", "completeness", "title", "recommendation",
    "specifications", "results", "metrics", "uncertainty", "comparisons", "sweeps", "visualisations",
    "validation", "provenance", "reproducibility", "resources", "privacy")
  if (!setequal(names(x), required) || anyDuplicated(names(x)))
    .ig_abort("invalid_input", "missing or unexpected report fields")
  scan <- function(value) {
    if (!is.list(value)) return(invisible(NULL))
    if (any(names(value) %in% c("ig_cls_data", "boundary_score", "scores", "raw_rows", "curve")))
      .ig_abort("invalid_input", "private score-derived content is not allowed in portable reports")
    if (!is.null(value$metric_id) && !is.null(value$value)) {
      for (i in seq_along(value$value)) {
        v <- value$value[[i]]
        if (!is.null(v) && (!is.numeric(v) || length(v) != 1L || !is.finite(v)))
          .ig_abort("invalid_input", "report metric values must be finite numbers or null")
      }
    }
    invisible(lapply(value, scan))
  }
  scan(x)
  if (!length(x$specifications) || length(x$specifications) != length(x$results))
    .ig_abort("invalid_input", "report specifications and results must be paired")
  for (entry in x$specifications) {
    if (identical(entry$analysis_kind, "power")) .ig_spec_from_list(entry$specification) else
      if (entry$analysis_kind %in% c("classification", "observed_test")) .ig_aux_spec_from_list(entry$specification, entry$analysis_kind) else
        .ig_abort("invalid_input", "unsupported report analysis kind")
  }
  expected <- x$snapshot_hash; x$snapshot_hash <- NULL
  if (!.ig_is_string(expected) || !identical(expected, .ig_hash_value(x)))
    .ig_abort("invalid_input", "report content hash mismatch")
  if (!.ig_is_string(x$locale) || !x$locale %in% c("en", "de") ||
      !.ig_is_string(x$report_kind) || !x$report_kind %in% .ig_report_kinds ||
      !.ig_is_utc_timestamp(x$created_at) || !isFALSE(x$privacy$raw_rows_included))
    .ig_abort("invalid_input", "invalid report metadata")
  .ig_report_validate_content(x)
  x$snapshot_hash <- expected
  structure(x, class = c("report_document", "report_snapshot", "ig_frozen"))
}

#' @export
format.report_document <- function(x, ...) {
  c(sprintf("<report_document> %s", x$snapshot_id), x$title,
    paste(x$report_kind, x$completeness), paste("Analysis hash:", x$analysis_hash))
}
