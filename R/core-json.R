
# A parsed JSON array of scalars (list) as an atomic vector; NULL stays NULL.
.ig_flatten <- function(x, mode = NULL) {
  if (is.null(x)) return(NULL)
  if (is.list(x)) {
    if (length(x) == 0L) {
      return(if (is.null(mode)) character() else vector(mode, 0L))
    }
    if (any(vapply(x, function(e) is.list(e) || length(e) != 1L, logical(1)))) {
      stop("expected an array of scalars", call. = FALSE)
    }
    x <- unlist(x, use.names = FALSE)
  }
  if (!is.null(mode)) x <- .ig_coerce(x, mode)
  x
}

.ig_coerce <- function(x, mode) {
  switch(mode,
         character = as.character(x),
         double = as.double(x),
         numeric = as.double(x),
         integer = as.integer(x),
         logical = as.logical(x),
         x)
}

# JSON Schema validation -----------------------------------------------------

# Bundled schemas (inst/schemas/<name>.schema.json). `common` holds shared
# definitions referenced by the others.
.ig_schema_names <- function() {
  c("power_spec", "power_result", "analysis_envelope", "report_snapshot", "report_document",
    "recommendation", "registry", "common")
}

.ig_schema_path <- function(name) {
  if (!(name %in% .ig_schema_names())) stop("unknown schema: ", name, call. = FALSE)
  path <- system.file("schemas", paste0(name, ".schema.json"), package = "inferguidR")
  if (!nzchar(path)) stop("bundled schema not found: ", name, call. = FALSE)
  path
}

# Validate JSON text against a bundled schema. Returns a list with
# `checked` (FALSE when 'jsonvalidate' is unavailable), `valid` and `errors`.
.ig_validate_schema <- function(json, name) {
  if (!requireNamespace("jsonvalidate", quietly = TRUE)) {
    return(list(checked = FALSE, valid = NA, errors = character(),
                note = "package 'jsonvalidate' is not installed; schema validation skipped"))
  }
  cache <- .ig_state$schema_validators
  validator <- cache[[name]]
  if (is.null(validator)) {
    validator <- jsonvalidate::json_validator(.ig_schema_path(name), engine = "ajv")
    cache[[name]] <- validator
  }
  if (inherits(json, "json")) json <- unclass(json)
  ok <- validator(json, verbose = TRUE, greedy = TRUE)
  errors <- character()
  if (!isTRUE(ok)) {
    err <- attr(ok, "errors")
    if (is.data.frame(err) && nrow(err)) {
      where <- if (!is.null(err$instancePath)) err$instancePath else err$dataPath
      where[is.na(where) | !nzchar(where)] <- "<root>"
      errors <- unique(paste0(where, ": ", err$message))
    } else {
      errors <- "document does not match the schema"
    }
  }
  list(checked = TRUE, valid = isTRUE(ok), errors = errors, note = "")
}

# Public JSON API -------------------------------------------------------------

#' Canonical JSON for inferguidR objects
#'
#' Serialises a specification, result, envelope, report snapshot or
#' recommendation (or any plain list of data) to canonical JSON: sorted keys,
#' shortest round-trip numbers and explicit markers for non-finite values.
#' The same object always produces the same text, so the text can be hashed.
#'
#' @param x An inferguidR object or a plain list/vector of data.
#' @param pretty If `TRUE`, indent the output (two spaces). Pretty output
#'   contains the same tokens; hashes are always computed on the compact form.
#' @return A character string of class `json`.
#' @examples
#' ig_to_json(list(b = 1, a = c(0.1 + 0.2, Inf)))
#' @seealso [ig_from_json()], [ig_hash()]
#' @family serialisation
#' @export
ig_to_json <- function(x, pretty = FALSE) {
  txt <- .ig_canonical_json(x)
  if (isTRUE(pretty)) txt <- as.character(jsonlite::prettify(txt, indent = 2))
  structure(txt, class = "json")
}

#' Read inferguidR objects from JSON safely
#'
#' Parses JSON produced by [ig_to_json()] and rebuilds the object through the
#' same validators that construct it in R. Imported text is treated as data:
#' nothing is evaluated, URLs are never fetched, duplicate keys are rejected,
#' and unknown object types or schema versions are refused with a typed
#' `invalid_input` error. When the 'jsonvalidate' package is installed and
#' `validate` is `TRUE`, the document is also checked against the bundled JSON
#' Schema.
#'
#' @param json A JSON string (or `NULL` when `file` is given).
#' @param file Path of a local JSON file (alternative to `json`).
#' @param type Optional expected object type (`"power_spec"`,
#'   `"power_result"`, `"analysis_envelope"`, `"report_snapshot"`,
#'   `"recommendation"`).
#' @param validate Validate against the JSON Schema when possible. Defaults
#'   to the JSON-schema validation option described in [ig_config()].
#' @return The rebuilt object. Errors are conditions of class
#'   `inferguidR_invalid_input` carrying an [ig_status()] in `$status`.
#' @examples
#' spec <- analysis_spec("t_means_one_sample", "post_hoc",
#'   inputs = list(alpha = 0.05, n_total = 20, d = 0.5),
#'   options = list(tails = "two"))
#' identical(ig_from_json(ig_to_json(spec)), spec)
#' @seealso [ig_to_json()]
#' @family serialisation
#' @export
ig_from_json <- function(json = NULL, file = NULL, type = NULL,
                         validate = getOption("inferguidR.validate_json_schema", TRUE)) {
  if (is.null(json) == is.null(file)) {
    .ig_abort("invalid_input", "supply exactly one of `json` or `file`")
  }
  if (!is.null(file)) {
    if (!.ig_is_string(file) || grepl("^[A-Za-z][A-Za-z0-9+.-]*://", file)) {
      .ig_abort("invalid_input", "`file` must be the path of a local file")
    }
    if (!file.exists(file)) .ig_abort("invalid_input", paste("file not found:", file))
    size <- file.size(file)
    if (is.na(size) || size > .ig_json_max_bytes) {
      .ig_abort("invalid_input", sprintf("file is larger than %d bytes", .ig_json_max_bytes))
    }
    json <- readChar(file, size, useBytes = TRUE)
    Encoding(json) <- "UTF-8"
  }
  if (inherits(json, "json")) json <- unclass(json)
  parsed <- .ig_json_parse(json)
  if (!is.list(parsed) || is.null(names(parsed))) {
    .ig_abort("invalid_input", "JSON document must be an object")
  }
  object_type <- parsed$object_type
  if (!.ig_is_string(object_type) || !(object_type %in% names(.ig_object_types()))) {
    .ig_abort("invalid_input",
              sprintf("unknown or missing object_type; expected one of: %s",
                      paste(names(.ig_object_types()), collapse = ", ")))
  }
  if (!is.null(type) && !identical(type, object_type)) {
    .ig_abort("invalid_input",
              sprintf("expected a %s but the document is a %s", type, object_type))
  }
  info <- .ig_object_types()[[object_type]]
  if (!identical(parsed$schema_version, info$schema_version) && !(identical(object_type,"report_document") && parsed$schema_version %in% c("inferguidR.report_document/0.2.0","inferguidR.report_document/0.3.0","inferguidR.report_document/0.4.0"))) {
    .ig_abort("invalid_input",
              sprintf("unsupported schema_version '%s' for %s; this version of inferguidR reads '%s'",
                      paste(format(parsed$schema_version), collapse = ""),
                      object_type, info$schema_version))
  }
  if (isTRUE(validate)) {
    sv <- .ig_validate_schema(paste(json, collapse = "\n"), info$schema)
    if (isTRUE(sv$checked) && !isTRUE(sv$valid)) {
      .ig_abort_invalid(sv$errors, sprintf("%s does not match its JSON Schema", object_type))
    }
  }
  info$from_list(parsed)
}

#' Content hash of an inferguidR object
#'
#' Returns the SHA-256 identity hash of an object: the specification hash of
#' a `power_spec`, the result hash of a `power_result`, the analysis hash of
#' an envelope or report snapshot, and the recommendation hash of a
#' recommendation. For plain data it is the SHA-256 of the canonical JSON.
#'
#' @param x An inferguidR object or plain data.
#' @return A 64-character lower-case hexadecimal string.
#' @examples
#' ig_hash(list(a = 1, b = "x")) == ig_hash(list(b = "x", a = 1))
#' @family serialisation
#' @export
ig_hash <- function(x) {
  if (inherits(x, "power_spec")) return(x$spec_hash)
  if (inherits(x, "power_result")) return(x$result_hash)
  if (inherits(x, "analysis_envelope")) return(x$analysis_hash)
  if (inherits(x, "report_snapshot")) return(x$analysis_hash)
  if (inherits(x, "recommendation")) return(x$recommendation_hash)
  .ig_hash_value(x)
}

# Object type table: schema version, schema name, list converters.
.ig_object_types <- function() {
  list(
    power_spec = list(schema_version = .ig_schema_versions[["power_spec"]],
                      schema = "power_spec", class = "power_spec",
                      to_list = .ig_spec_to_list, from_list = .ig_spec_from_list),
    power_result = list(schema_version = .ig_schema_versions[["power_result"]],
                        schema = "power_result", class = "power_result",
                        to_list = .ig_result_to_list, from_list = .ig_result_from_list),
    analysis_envelope = list(schema_version = .ig_schema_versions[["analysis_envelope"]],
                             schema = "analysis_envelope", class = "analysis_envelope",
                             to_list = .ig_envelope_to_list, from_list = .ig_envelope_from_list),
    report_document = list(schema_version = "inferguidR.report_document/0.1.0",
                           schema = "report_document", class = "report_document",
                           to_list = .ig_report_document_to_list, from_list = .ig_report_document_from_list),
    report_snapshot = list(schema_version = .ig_schema_versions[["report_snapshot"]],
                           schema = "report_snapshot", class = "report_snapshot",
                           to_list = .ig_snapshot_to_list, from_list = .ig_snapshot_from_list),
    recommendation = list(schema_version = .ig_schema_versions[["recommendation"]],
                          schema = "recommendation", class = "recommendation",
                          to_list = .ig_recommendation_to_list,
                          from_list = .ig_recommendation_from_list)
  )
}

# Schema versions of the serialised objects.
.ig_schema_versions <- c(
  power_spec = "inferguidR.power_spec/0.1.0",
  power_result = "inferguidR.power_result/0.1.0",
  analysis_envelope = "inferguidR.analysis_envelope/0.1.0",
  report_snapshot = "inferguidR.report_snapshot/0.1.0-draft",
  recommendation = "inferguidR.recommendation/0.1.0",
  classification_spec = "inferguidR.classification_spec/0.1.0-draft",
  observed_test_spec = "inferguidR.observed_test_spec/0.1.0-draft"
)

# Convert ig objects (possibly nested inside plain lists) to plain lists.
.ig_to_plain <- function(x) {
  if (is.object(x) && !inherits(x, "AsIs")) {
    if (inherits(x, "ig_status")) return(.ig_status_to_list(x))
    for (tp in .ig_object_types()) {
      if (inherits(x, tp$class)) return(tp$to_list(x))
    }
    if (inherits(x, .ig_frozen_aux_classes)) return(.ig_unclass_frozen(x))
  }
  if (is.list(x) && !is.data.frame(x) && length(x)) {
    nms <- names(x)
    asis <- inherits(x, "AsIs")
    out <- lapply(unclass(x), .ig_to_plain)
    names(out) <- nms
    if (asis) out <- I(out)
    return(out)
  }
  x
}

# Spec/result classes owned by later phases that the envelope carries as
# frozen plain data.
.ig_frozen_aux_classes <- c("classification_spec", "observed_test_spec",
                            "classification_result", "observed_test_result")
