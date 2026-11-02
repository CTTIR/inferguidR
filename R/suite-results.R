# Public shared views leave historical power objects and their identities intact.
#' @export
as_inference.power_result <- function(x, ...) {
  status <- x$status$code
  rows <- list()
  for (section in c("solved", "achieved", "quantities")) {
    for (name in names(x[[section]])) {
      value <- x[[section]][[name]]
      if (!is.numeric(value)) next
      for (i in seq_along(value)) {
        id <- paste(section, name, i, sep = ":")
        rows[[length(rows) + 1L]] <- data.frame(row_id = id, component = section,
          estimand = name, estimate = as.double(value[i]), status = status,
          reason = if (status == "ok") NA_character_ else x$status$reason,
          stringsAsFactors = FALSE)
      }
    }
  }
  if (!length(rows)) rows <- list(data.frame(row_id = "planning:status", component = "planning",
    estimand = NA_character_, estimate = NA_real_, status = status, reason = x$status$reason,
    stringsAsFactors = FALSE))
  infercalcR::new_inference_result(do.call(rbind, rows), list(method_id = x$method_id,
    variant_id = x$spec$options$variant %||% x$method_id, result_id = x$result_hash,
    status = x$status, mode = x$mode, assumptions = x$assumptions, warnings = x$warnings,
    validation = x$validation, provenance = x$provenance, runtime = x$runtime,
    population = "prospective study assumptions", legacy_schema = x$schema_version,
    legacy_identity = x$result_hash), kind = "planning", payload = unclass(x))
}
#' @export
result_data.power_result <- function(x, ...) infercalcR::result_data(infercalcR::as_inference(x))
#' @export
result_metadata.power_result <- function(x, ...) infercalcR::result_metadata(infercalcR::as_inference(x))

#' Observed-method relations to planning engines
#' @export
planning_methods <- function(method_id = NULL) {
  file <- system.file("suite", "planning-relations.csv", package = "inferguidR")
  out <- utils::read.csv(file, stringsAsFactors = FALSE, na.strings = "")
  if (!is.null(method_id)) {
    if (!is.character(method_id) || anyNA(method_id) || any(!method_id %in% out$method_id))
      .ig_abort("invalid_input", "method_id must identify an inventoried observed method")
    out <- out[match(method_id, out$method_id), , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}
