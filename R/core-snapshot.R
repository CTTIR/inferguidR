# Report snapshots (draft; phase 13A finalises exporters and the public API) -------------
#
# A report_snapshot freezes one or more analysis envelopes, an optional
# recommendation, plot data and texts at one moment. Every exporter (PDF,
# JSON, TXT) reads only the snapshot, so input changes during rendering
# cannot change an export. Identity:
#
# * `analysis_hash`: deterministic hash of the analysed content (envelope
#   analysis hashes, recommendation hash, report kind); independent of time;
# * `snapshot_id`: identifier of this freeze (content plus creation time);
# * `snapshot_hash`: integrity hash over the whole frozen snapshot.
#
# No rendering happens here.

.ig_report_kinds <- c("single_analysis", "model_cutoff", "comparison", "session",
                      "incomplete_settings")

.ig_utc_now <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
}

.ig_is_utc_timestamp <- function(x) {
  .ig_is_string(x) && grepl("^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}(\\.\\d{1,6})?Z$", x)
}

# Internal constructor (phase 13A exposes create_report_snapshot()).
.ig_new_report_snapshot <- function(envelopes, recommendation = NULL,
                                    report_kind = "single_analysis", title = NULL,
                                    locale = "en", display_timezone = NULL,
                                    plot_data = list(), texts = list(),
                                    created_at = NULL, runtime = NULL,
                                    snapshot_id = NULL) {
  if (inherits(envelopes, "analysis_envelope")) envelopes <- list(envelopes)
  problems <- character()
  if (!is.list(envelopes) || !length(envelopes) ||
      !all(vapply(envelopes, inherits, logical(1), what = "analysis_envelope"))) {
    .ig_abort_invalid("`envelopes` must be a non-empty list of analysis_envelope objects",
                      "invalid report snapshot")
  }
  if (!.ig_is_string(report_kind) || !(report_kind %in% .ig_report_kinds)) {
    problems <- c(problems, sprintf("report_kind must be one of: %s", paste(.ig_report_kinds, collapse = ", ")))
  }
  if (!.ig_is_string(locale) || !(locale %in% c("en", "de"))) problems <- c(problems, "locale must be 'en' or 'de'")
  if (!is.null(title) && !.ig_is_string(title)) problems <- c(problems, "title must be a string or NULL")
  if (!is.null(recommendation) && !inherits(recommendation, "recommendation")) {
    problems <- c(problems, "recommendation must be a recommendation object or NULL")
  }
  if (!is.list(plot_data) || (length(plot_data) && !is.null(names(plot_data)))) {
    problems <- c(problems, "plot_data must be an unnamed list of plot-data records")
  }
  if (!is.list(texts) || (length(texts) && !.ig_is_named(texts)) ||
      !all(vapply(texts, .ig_is_string, logical(1)))) {
    problems <- c(problems, "texts must be a named list of strings")
  }
  if (is.null(created_at)) created_at <- .ig_utc_now()
  if (!.ig_is_utc_timestamp(created_at)) {
    problems <- c(problems, "created_at must be an ISO-8601 UTC timestamp (YYYY-MM-DDTHH:MM:SS[.sss]Z)")
  }
  if (is.null(display_timezone)) {
    display_timezone <- Sys.timezone()
    if (is.na(display_timezone) || !nzchar(display_timezone)) display_timezone <- "UTC"
  }
  if (!.ig_is_string(display_timezone)) problems <- c(problems, "display_timezone must be a string")

  # Completeness: a final-result report needs computed, ok results.
  status_of <- function(env) {
    if (is.null(env$result)) return("missing")
    st <- env$result$status
    if (is_ig_status(st)) st$code else if (is.list(st)) as.character(st$code) else "unknown"
  }
  codes <- vapply(envelopes, status_of, "")
  complete <- all(codes == "ok")
  if (.ig_is_string(report_kind)) {
    if (report_kind != "incomplete_settings" && !complete) {
      problems <- c(problems, sprintf(
        "a '%s' report needs completed results with status 'ok' (found: %s); use report_kind 'incomplete_settings' for a labelled settings-only report",
        report_kind, paste(unique(codes[codes != "ok"]), collapse = ", ")))
    }
    if (report_kind == "single_analysis" && length(envelopes) != 1L) {
      problems <- c(problems, "a 'single_analysis' report holds exactly one envelope")
    }
    multiple_models <- length(envelopes) == 1L &&
      identical(envelopes[[1]]$analysis_kind, "classification") && length(envelopes[[1]]$spec$models) >= 2L
    if (report_kind == "comparison" && length(envelopes) < 2L && !multiple_models) {
      problems <- c(problems, "a 'comparison' report needs two envelopes or multiple classifier models in one envelope")
    }
    if (report_kind == "model_cutoff" &&
        !any(vapply(envelopes, function(e) identical(e$analysis_kind, "classification"), logical(1)))) {
      problems <- c(problems, "a 'model_cutoff' report needs a classification envelope")
    }
  }
  if (length(problems)) .ig_abort_invalid(problems, "invalid report snapshot")

  analysis_hash <- .ig_hash_value(list(
    report_kind = report_kind,
    envelopes = I(vapply(envelopes, `[[`, "", "analysis_hash")),
    recommendation = if (is.null(recommendation)) NULL else recommendation$recommendation_hash
  ))
  if (is.null(runtime)) runtime <- c(.ig_runtime(), list(platform = R.version$platform))
  runtime <- .ig_sort_names(runtime)
  plot_data <- .ig_sort_names(plot_data)
  texts <- .ig_sort_names(texts)
  computed_id <- paste0("snap-", substr(.ig_hash_value(list(
    analysis_hash = analysis_hash, created_at = created_at, locale = locale, title = title
  )), 1L, 24L))
  if (!is.null(snapshot_id) && !identical(snapshot_id, computed_id)) {
    .ig_abort("invalid_input", "snapshot_id does not match the snapshot content and creation time")
  }
  x <- list(
    schema_version = .ig_schema_versions[["report_snapshot"]],
    object_type = "report_snapshot",
    snapshot_id = computed_id,
    analysis_hash = analysis_hash,
    created_at = created_at,
    display_timezone = display_timezone,
    locale = locale,
    report_kind = report_kind,
    completeness = if (complete) "complete" else "incomplete",
    title = title,
    envelopes = envelopes,
    recommendation = recommendation,
    plot_data = plot_data,
    texts = .ig_named_or_empty(texts),
    runtime = runtime,
    snapshot_hash = NA_character_
  )
  x$snapshot_hash <- .ig_snapshot_hash(x)
  structure(x, class = c("report_snapshot", "ig_frozen"))
}

.ig_snapshot_hash <- function(x) {
  lst <- .ig_snapshot_to_list(x)
  lst$snapshot_hash <- NULL
  .ig_hash_value(lst)
}

.ig_snapshot_to_list <- function(x) {
  x <- unclass(x)
  list(
    schema_version = x$schema_version, object_type = x$object_type,
    snapshot_id = x$snapshot_id, analysis_hash = x$analysis_hash,
    created_at = x$created_at, display_timezone = x$display_timezone,
    locale = x$locale, report_kind = x$report_kind, completeness = x$completeness,
    title = x$title,
    envelopes = lapply(x$envelopes, .ig_envelope_to_list),
    recommendation = if (is.null(x$recommendation)) NULL else .ig_recommendation_to_list(x$recommendation),
    plot_data = .ig_to_plain(x$plot_data),
    texts = .ig_named_or_empty(x$texts),
    runtime = x$runtime,
    snapshot_hash = x$snapshot_hash
  )
}

.ig_snapshot_from_list <- function(x) {
  need <- c("schema_version", "object_type", "snapshot_id", "analysis_hash", "created_at",
            "envelopes", "snapshot_hash")
  miss <- setdiff(need, names(x))
  if (length(miss)) .ig_abort_invalid(sprintf("missing field '%s'", miss), "malformed report_snapshot")
  envelopes <- lapply(x$envelopes, .ig_envelope_from_list)
  rec <- if (is.null(x$recommendation)) NULL else .ig_recommendation_from_list(x$recommendation)
  texts <- x$texts %||% list()
  snap <- .ig_new_report_snapshot(envelopes, recommendation = rec,
                                  report_kind = x$report_kind, title = x$title,
                                  locale = x$locale, display_timezone = x$display_timezone,
                                  plot_data = .ig_simplify_parsed(x$plot_data %||% list()),
                                  texts = texts, created_at = x$created_at,
                                  runtime = x$runtime, snapshot_id = x$snapshot_id)
  if (!identical(snap$analysis_hash, x$analysis_hash) ||
      !identical(snap$snapshot_hash, x$snapshot_hash)) {
    .ig_abort("invalid_input", "snapshot hashes do not match the content (the document was modified)")
  }
  snap
}

#' @export
format.report_snapshot <- function(x, ...) {
  c(sprintf("<report_snapshot> %s", x$snapshot_id),
    sprintf("  kind:          %s (%s)", x$report_kind, x$completeness),
    sprintf("  created (UTC): %s", x$created_at),
    sprintf("  envelopes:     %d", length(x$envelopes)),
    sprintf("  analysis hash: %s", x$analysis_hash))
}

#' @export
print.report_snapshot <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}
