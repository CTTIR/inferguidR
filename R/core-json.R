# Canonical JSON, hashing and safe import ------------------------------------
#
# Canonical form (used for every hash and every export):
#
# * object keys are sorted in C-locale (byte) order at every level;
# * no insignificant whitespace;
# * numbers use the shortest decimal representation that reads back to the
#   identical IEEE double (15, 16 or 17 significant digits); negative zero is
#   written as 0. jsonlite's `digits = NA` writes at most 15 significant
#   digits, which does not round-trip (0.1 + 0.2 is written as 0.3), so the
#   number text is produced here and jsonlite is used for string escaping,
#   structural parsing and pretty-printing; numerical tokens use R conversion;
# * missing values (NA) are `null`; the object that holds them must say why;
# * non-finite numbers are written as the marker object
#   {"nonfinite": "positive_infinity" | "negative_infinity" | "not_a_number",
#    "value": null} because JSON has no Inf or NaN tokens;
# * length-1 atomic vectors are scalars unless wrapped in I(); named atomic
#   vectors and named lists are objects; unnamed lists are arrays.
#
# Hashes are SHA-256 (digest) of the UTF-8 bytes of the canonical text.
# Imported JSON is parsed as data only: no R code is ever evaluated.

.ig_nonfinite_codes <- c("positive_infinity", "negative_infinity", "not_a_number")

# Maximum accepted size of an imported JSON document (bytes).
.ig_json_max_bytes <- 64 * 1024^2

# Shortest round-trip decimal text for doubles (non-finite values are
# returned as R prints them; the JSON writer encodes them separately).
.ig_format_double <- function(x) {
  out <- as.character(x)
  fin <- which(is.finite(x))
  if (!length(fin)) return(out)
  xf <- x[fin]
  txt <- sprintf("%.15g", xf)
  bad <- as.numeric(txt) != xf
  if (any(bad)) {
    txt[bad] <- sprintf("%.16g", xf[bad])
    bad <- as.numeric(txt) != xf
    if (any(bad)) txt[bad] <- sprintf("%.17g", xf[bad])
  }
  txt[xf == 0] <- "0"
  out[fin] <- txt
  out
}

.ig_nonfinite_marker <- function(code) {
  sprintf('{"nonfinite":"%s","value":null}', code)
}

# JSON text for the elements of an atomic vector (one string per element).
.ig_json_atomic_elements <- function(x, path) {
  if (is.factor(x)) x <- as.character(x)
  if (is.logical(x)) {
    out <- ifelse(x, "true", "false")
    out[is.na(x)] <- "null"
    return(out)
  }
  if (is.integer(x)) {
    out <- as.character(x)
    out[is.na(x)] <- "null"
    return(out)
  }
  if (is.double(x)) {
    fin <- is.finite(x)
    if (all(fin)) return(.ig_format_double(x))
    out <- character(length(x))
    out[fin] <- .ig_format_double(x[fin])
    out[is.na(x) & !is.nan(x)] <- "null"
    out[is.nan(x)] <- .ig_nonfinite_marker("not_a_number")
    out[!is.na(x) & x == Inf] <- .ig_nonfinite_marker("positive_infinity")
    out[!is.na(x) & x == -Inf] <- .ig_nonfinite_marker("negative_infinity")
    return(out)
  }
  if (is.character(x)) {
    x <- enc2utf8(x)
    # Printable ASCII needs no escaping except quotes, backslashes and HTML
    # closing-tag prefixes. Keep all other strings on the established encoder.
    safe <- !is.na(x) & !grepl('[^ -~]|["\\\\<]', x, useBytes = TRUE)
    out <- rep("null",length(x))
    if(any(safe)) out[safe] <- paste0('"',x[safe],'"')
    escaped <- which(!is.na(x) & !safe)
    if(length(escaped)) out[escaped] <- vapply(x[escaped],function(s)
      as.character(jsonlite::toJSON(s,auto_unbox=TRUE)),character(1),USE.NAMES=FALSE)
    return(out)
  }
  stop(sprintf("cannot serialise a value of type '%s' at %s", typeof(x),
               if (nzchar(path)) path else "<root>"), call. = FALSE)
}

.ig_json_object <- function(keys, values, path) {
  if (anyNA(keys) || !all(nzchar(keys))) {
    stop(sprintf("partially named list at %s cannot be serialised",
                 if (nzchar(path)) path else "<root>"), call. = FALSE)
  }
  if (anyDuplicated(keys)) {
    stop(sprintf("duplicate key '%s' at %s", keys[anyDuplicated(keys)],
                 if (nzchar(path)) path else "<root>"), call. = FALSE)
  }
  ord <- order(keys, method = "radix")
  key_txt <- .ig_json_atomic_elements(keys[ord], path)
  paste0("{", paste0(key_txt, ":", values[ord], collapse = ","), "}")
}

# Recursive canonical serialiser.
.ig_json_value <- function(x, path = "") {
  if (is.null(x)) return("null")
  if (inherits(x, "ig_status")) x <- .ig_status_to_list(x)
  if (is.function(x) || is.environment(x) || inherits(x, "formula")) {
    stop(sprintf("cannot serialise a %s at %s", class(x)[1],
                 if (nzchar(path)) path else "<root>"), call. = FALSE)
  }
  if (is.data.frame(x)) x <- as.list(x)
  if (is.list(x)) {
    nms <- names(x)
    if (length(x) == 0L) return(if (is.null(nms)) "[]" else "{}")
    if (is.null(nms)) {
      parts <- vapply(seq_along(x), function(i) {
        .ig_json_value(x[[i]], paste0(path, "[", i, "]"))
      }, character(1))
      return(paste0("[", paste(parts, collapse = ","), "]"))
    }
    parts <- vapply(seq_along(x), function(i) {
      .ig_json_value(x[[i]], paste0(path, "/", nms[i]))
    }, character(1))
    return(.ig_json_object(nms, parts, path))
  }
  if (!is.atomic(x)) {
    stop(sprintf("cannot serialise a value of type '%s' at %s", typeof(x),
                 if (nzchar(path)) path else "<root>"), call. = FALSE)
  }
  elements <- .ig_json_atomic_elements(unname(x), path)
  if (!is.null(names(x))) {
    return(.ig_json_object(names(x), elements, path))
  }
  if (length(x) == 1L && !inherits(x, "AsIs")) return(elements)
  paste0("[", paste(elements, collapse = ","), "]")
}

# Canonical JSON text of an R value (plain lists/vectors or ig objects).
.ig_canonical_json <- function(x) {
  enc2utf8(.ig_json_value(.ig_to_plain(x)))
}

.ig_sha256 <- function(txt) {
  digest::digest(enc2utf8(txt), algo = "sha256", serialize = FALSE)
}

# SHA-256 of the canonical JSON of a plain value.
.ig_hash_value <- function(x) .ig_sha256(.ig_canonical_json(x))

# Parsing ------------------------------------------------------------------

# Parse JSON text into plain lists (no simplification), reject duplicate
# keys and decode non-finite markers. Never evaluates anything.
.ig_json_parse <- function(txt) {
  if (inherits(txt, "json")) txt <- unclass(txt)
  if (!is.character(txt) || length(txt) < 1L || anyNA(txt)) {
    .ig_abort("invalid_input", "JSON input must be a character string")
  }
  txt <- paste(txt, collapse = "\n")
  if (nchar(txt, type = "bytes") > .ig_json_max_bytes) {
    .ig_abort("invalid_input", sprintf("JSON input exceeds %d bytes", .ig_json_max_bytes))
  }
  parsed <- tryCatch(
    jsonlite::parse_json(txt, simplifyVector = FALSE),
    error = function(e) {
      .ig_abort("invalid_input", paste("input is not valid JSON:", conditionMessage(e)))
    }
  )
  dup <- .ig_find_duplicate_keys(parsed)
  if (length(dup)) {
    .ig_abort_invalid(sprintf("duplicate key '%s'", dup),
                      "JSON objects must not repeat keys")
  }
  .ig_decode_nonfinite(.ig_json_restore_numbers(parsed, txt))
}


# jsonlite structural parsing can round a valid decimal one ULP away from
# R's conversion used by the canonical writer. Read numerical tokens using
# that same conversion, preserving string contents and the parsed structure.
.ig_json_restore_numbers <- function(parsed, txt) {
  pattern <- '"(?:[^"\\\\]++|\\\\.)*+"(*SKIP)(*F)|-?(?:0|[1-9][0-9]*)(?:\\.[0-9]+)?(?:[eE][+-]?[0-9]+)?'
  positions <- gregexpr(pattern, txt, perl = TRUE)
  tokens <- regmatches(txt, positions)[[1L]]
  values <- as.numeric(tokens)
  cursor <- 0L
  restore <- function(x) {
    if (is.numeric(x)) {
      cursor <<- cursor + 1L
      if (length(x) != 1L || cursor > length(values))
        .ig_abort("invalid_input", "JSON numeric token structure mismatch")
      value <- values[[cursor]]
      if (is.integer(x) && is.finite(value) && value == x) x else value
    } else if (is.list(x)) {
      for (i in seq_along(x)) x[i] <- list(restore(x[[i]]))
      x
    } else x
  }
  out <- restore(parsed)
  if (cursor != length(values))
    .ig_abort("invalid_input", "JSON numeric token count mismatch")
  out
}

.ig_find_duplicate_keys <- function(x, path = "") {
  if (!is.list(x)) return(character())
  out <- character()
  nms <- names(x)
  if (!is.null(nms)) {
    d <- unique(nms[duplicated(nms)])
    if (length(d)) out <- c(out, paste0(path, "/", d))
  }
  for (i in seq_along(x)) {
    sub <- if (is.null(nms)) paste0(path, "[", i, "]") else paste0(path, "/", nms[i])
    out <- c(out, .ig_find_duplicate_keys(x[[i]], sub))
  }
  out
}

.ig_is_nonfinite_marker <- function(x) {
  is.list(x) && length(x) == 2L && setequal(names(x), c("nonfinite", "value")) &&
    is.null(x$value) && .ig_is_string(x$nonfinite) &&
    x$nonfinite %in% .ig_nonfinite_codes
}

.ig_decode_nonfinite <- function(x) {
  if (.ig_is_nonfinite_marker(x)) {
    return(switch(x$nonfinite,
                  positive_infinity = Inf,
                  negative_infinity = -Inf,
                  not_a_number = NaN))
  }
  if (is.list(x)) {
    for (i in seq_along(x)) {
      val <- .ig_decode_nonfinite(x[[i]])
      if (!is.null(val)) x[[i]] <- val
    }
  }
  x
}

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
