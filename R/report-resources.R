# Versioned teaching metadata is retained independently of mutable web pages.
# Keep historical catalogs when revising a lesson so old reports remain readable.
.ig_report_resource_catalog <- function(version = "2026-10-13.1") {
  if (!identical(version, "2026-10-13.1"))
    .ig_abort("invalid_input", "Unknown report resource catalog version.")
  path <- system.file("metadata", paste0("report-resources-", version, ".json"), package = "inferguidR")
  if (!nzchar(path)) .ig_abort("unsupported", "Report resource catalog is unavailable.")
  jsonlite::read_json(path, simplifyVector = FALSE)
}

.ig_report_resources <- function(requests, specifications, locale, version = "2026-10-13.1") {
  empty <- list(status = "not_attached", reason = "No verified topic-specific resource bundle was supplied.")
  if (!is.list(requests) || length(requests) > 20L || (length(requests) && !is.null(names(requests))))
    .ig_abort("invalid_input", "resources must be a list of at most 20 entries.")
  if (!length(requests)) return(empty)
  catalog <- .ig_report_resource_catalog(version)
  entries <- lapply(requests, function(request) {
    if (!is.list(request) || anyDuplicated(names(request)) ||
        !setequal(names(request), c("analysis", "resource_id")))
      .ig_abort("invalid_input", "Each resource needs only analysis and resource_id.")
    i <- request$analysis
    if (!is.numeric(i) || length(i) != 1L || !is.finite(i) || i != floor(i) || i < 1L || i > length(specifications))
      .ig_abort("invalid_input", "Invalid resource analysis index.")
    if (!.ig_is_string(request$resource_id)) .ig_abort("invalid_input", "Invalid resource identifier.")
    hits <- Filter(function(entry) identical(entry$resource_id, request$resource_id), catalog$entries)
    if (length(hits) != 1L) .ig_abort("invalid_input", "Unknown verified report resource.")
    spec <- specifications[[i]]
    valid <- .ig_report_resource_matches(spec)
    if (!valid) .ig_abort("invalid_input", "Resource does not teach the selected analysis variant.")
    entry <- hits[[1L]]
    list(analysis = as.integer(i), resource_id = entry$resource_id,
      title = entry[[paste0("title_", locale)]], locale = locale, match = entry$match,
      analysis_variant_id = entry$analysis_variant_id,
      scope = entry[[paste0("scope_", locale)]], source = entry[[paste0("source_", locale)]],
      access = entry[[paste0("access_", locale)]], files = entry$files,
      references = lapply(entry$references, function(reference)
        list(title = reference[[paste0("title_", locale)]], url = reference$url, role = reference$role)))
  })
  keys <- vapply(entries, function(e) paste(e$analysis, e$resource_id), "")
  if (anyDuplicated(keys)) .ig_abort("invalid_input", "Duplicate report resource.")
  list(status = "attached", catalog_version = version, entries = entries,
    evidence_note = if (locale == "de")
      "Lehrmaterial; kein numerischer Validierungsnachweis. Externe Webseiten koennen sich aendern." else
      "Teaching material, not numerical validation evidence. External web pages may change.")
}

.ig_report_validate_resources <- function(x) {
  r <- x$resources
  if (!is.list(r) || anyDuplicated(names(r)) || !.ig_is_string(x$locale) || !x$locale %in% c("en", "de"))
    .ig_abort("invalid_input", "Invalid resource declaration or locale.")
  if (identical(r$status, "attached")) {
    if (!x$schema_version %in% c("inferguidR.report_document/0.3.0","inferguidR.report_document/0.4.0") ||
        !is.list(r$entries) || !length(r$entries) || !is.null(names(r$entries)))
      .ig_abort("invalid_input", "Invalid version or entries for attached report resources.")
    requests <- lapply(r$entries, function(e) {
      if (!is.list(e)) .ig_abort("invalid_input", "Invalid resource record.")
      list(analysis = e$analysis, resource_id = e$resource_id)
    })
    expected <- .ig_report_resources(requests, x$specifications, x$locale, r$catalog_version)
  } else expected <- .ig_report_resources(list(), x$specifications, x$locale)
  if (!identical(.ig_canonical_json(r), .ig_canonical_json(expected)))
    .ig_abort("invalid_input", "Unsupported or altered report resource declaration.")
  invisible(TRUE)
}

.ig_report_resource_matches <- function(spec) {
  (identical(spec$analysis_kind, "observed_test") && identical(spec$specification$variant_id, "boschloo_minlike")) ||
    (identical(spec$analysis_kind, "power") &&
     identical(spec$specification$method_id, "boschloo_proportions_two_independent") &&
     identical(spec$specification$options$variant, "minlike"))
}

.ig_report_envelope_resources <- function(envelope) {
  if (!inherits(envelope, "analysis_envelope")) return(list())
  spec <- list(analysis_kind = envelope$analysis_kind, specification = unclass(envelope$spec))
  if (!.ig_report_resource_matches(spec)) return(list())
  list(list(analysis = 1L, resource_id = "inferguidR_boschloo_minlike"))
}
