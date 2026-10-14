# Analysis envelopes ------------------------------------------------------------------
#
# An analysis_envelope is the discriminated container shared by all analysis
# kinds: `power` (power_spec / power_result), `classification`
# (classification_spec / classification_result, phase 10A) and
# `observed_test` (observed_test_spec / observed_test_result, phase 09A).
# The core defines the object shapes and validates them; the user-facing
# constructors and the computations of the later kinds belong to their
# phases. `analysis_hash` is the SHA-256 of (kind, spec hash, result hash);
# provenance (package and R versions) is kept outside the hash.

.ig_analysis_kinds <- c("power", "classification", "observed_test")

#' Wrap a specification and its result in an analysis envelope
#'
#' The envelope is the unit that reports, comparisons and the application
#' pass around. It pairs one immutable specification with its result (or no
#' result yet) and records shared provenance.
#'
#' @param spec A `power_spec` (or, from later phases, a `classification_spec`
#'   or `observed_test_spec`).
#' @param result The matching result, or `NULL` when not computed yet.
#' @param kind Optional analysis kind (`"power"`, `"classification"`,
#'   `"observed_test"`); inferred from `spec` when omitted.
#' @return An immutable object of class `analysis_envelope`.
#' @examples
#' spec <- analysis_spec("t_means_one_sample", "post_hoc",
#'   inputs = list(alpha = 0.05, n_total = 20, d = 0.5),
#'   options = list(tails = "two"))
#' analysis_envelope(spec, calculate_power(spec))
#' @export
analysis_envelope <- function(spec, result = NULL, kind = NULL) {
  .ig_new_envelope(spec, result, kind, provenance = NULL)
}

.ig_new_envelope <- function(spec, result, kind, provenance) {
  inferred <- if (inherits(spec, "power_spec")) "power" else
    if (inherits(spec, "classification_spec")) "classification" else
      if (inherits(spec, "observed_test_spec")) "observed_test" else NA_character_
  if (is.na(inferred)) {
    .ig_abort("invalid_input",
              "`spec` must be a power_spec, classification_spec or observed_test_spec")
  }
  if (!is.null(kind) && (!.ig_is_string(kind) || !identical(kind, inferred))) {
    .ig_abort("invalid_input", sprintf("analysis kind '%s' does not match a %s",
                                       paste(format(kind), collapse = ""), class(spec)[1]))
  }
  kind <- inferred
  result_hash <- NULL
  if (kind == "power") {
    .ig_check_spec_integrity(spec)
    if (!is.null(result)) {
      if (!inherits(result, "power_result")) {
        .ig_abort("invalid_input", "the result of a power analysis must be a power_result")
      }
      if (!identical(result$spec_hash, spec$spec_hash)) {
        .ig_abort("invalid_input", "the result belongs to a different specification (spec_hash differs)")
      }
      if (!identical(.ig_result_hash(result), result$result_hash)) {
        .ig_abort("invalid_input", "the result was modified after creation (hash mismatch)")
      }
      result_hash <- result$result_hash
    }
  } else {
    .ig_check_aux_spec(spec, kind)
    if (!is.null(result)) {
      cls <- paste0(kind, "_result")
      if (!inherits(result, cls) || !.ig_is_string(result$spec_hash) ||
          !is.list(result$status)) {
        .ig_abort("invalid_input", sprintf("the result must be a %s with spec_hash and status", cls))
      }
      if (!identical(result$spec_hash, spec$spec_hash)) {
        .ig_abort("invalid_input", "the result belongs to a different specification (spec_hash differs)")
      }
      result_hash <- .ig_hash_value(.ig_unclass_frozen(result))
    }
  }
  if (is.null(provenance)) provenance <- .ig_runtime()
  provenance <- .ig_sort_names(provenance)
  analysis_hash <- .ig_hash_value(list(analysis_kind = kind, spec_hash = spec$spec_hash,
                                       result_hash = result_hash))
  structure(list(
    schema_version = .ig_schema_versions[["analysis_envelope"]],
    object_type = "analysis_envelope",
    analysis_kind = kind,
    spec = spec,
    result = result,
    spec_hash = spec$spec_hash,
    result_hash = result_hash,
    analysis_hash = analysis_hash,
    provenance = provenance
  ), class = c("analysis_envelope", "ig_frozen"))
}

.ig_envelope_to_list <- function(x) {
  spec <- if (x$analysis_kind == "power") .ig_spec_to_list(x$spec) else .ig_unclass_frozen(x$spec)
  result <- if (is.null(x$result)) NULL else if (x$analysis_kind == "power") {
    .ig_result_to_list(x$result)
  } else {
    .ig_unclass_frozen(x$result)
  }
  list(schema_version = x$schema_version, object_type = x$object_type,
       analysis_kind = x$analysis_kind, spec = spec, result = result,
       spec_hash = x$spec_hash, result_hash = x$result_hash,
       analysis_hash = x$analysis_hash, provenance = x$provenance)
}

.ig_envelope_from_list <- function(x) {
  need <- c("schema_version", "object_type", "analysis_kind", "spec", "spec_hash", "analysis_hash")
  miss <- setdiff(need, names(x))
  if (length(miss)) .ig_abort_invalid(sprintf("missing field '%s'", miss), "malformed analysis_envelope")
  kind <- x$analysis_kind
  if (!.ig_is_string(kind) || !(kind %in% .ig_analysis_kinds)) {
    .ig_abort("invalid_input", "unknown analysis_kind")
  }
  if (kind == "power") {
    spec <- .ig_spec_from_list(x$spec)
    result <- if (is.null(x$result)) NULL else .ig_result_from_list(x$result)
  } else {
    # Later phases provide full rebuilders; until then the shape is checked
    # and the content hash must match.
    spec <- .ig_aux_spec_from_list(x$spec, kind)
    result <- if (is.null(x$result)) NULL else
      structure(.ig_simplify_parsed(x$result), class = c(paste0(kind, "_result"), "ig_frozen"))
  }
  env <- .ig_new_envelope(spec, result, kind, provenance = x$provenance)
  if (!identical(env$analysis_hash, x$analysis_hash)) {
    .ig_abort("invalid_input", "analysis_hash does not match the envelope content (the document was modified)")
  }
  env
}

#' @export
format.analysis_envelope <- function(x, ...) {
  status <- if (is.null(x$result)) "no result yet" else if (is.list(x$result$status)) {
    format(.ig_status_from_list(unclass(x$result$status)))
  } else {
    "unknown"
  }
  c(sprintf("<analysis_envelope> %s", x$analysis_kind),
    sprintf("  spec:          %s", if (x$analysis_kind == "power") paste(x$spec$method_id, x$spec$mode, sep = " | ") else x$spec$object_type),
    sprintf("  result:        %s", status),
    sprintf("  analysis hash: %s", x$analysis_hash))
}

#' @export
print.analysis_envelope <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}

# Classification specification shape (draft; phase 10A owns the constructor
# and the computations) --------------------------------------------------------------------

.ig_classification_vocab <- list(
  input_kind = c("theoretical_binormal", "scored_data", "counts_2x2"),
  direction = c("higher_is_positive", "lower_is_positive"),
  comparator = c(">=", ">", "<=", "<"),
  score_scale = c("probability", "raw_score"),
  evaluation_split = c("same_data", "held_out_tuning", "held_out_test", "external",
                       "resampling", "not_applicable"),
  selection = c("manual", "max_youden", "max_f1", "max_mcc", "constraint", "capacity",
                "cost_utility"),
  selection_data = c("same_data", "held_out_tuning", "predeclared_rule"),
  prevalence_kind = c("observed", "scenario")
)

# Low-level constructor of a classification_spec. Validates the shape only:
# model metadata, operating rule (score cutoff, comparator, selection),
# prevalence, counts or theory parameters, uncertainty settings. The term
# `score_cutoff` is a classifier threshold, never a significance critical
# value.
.ig_new_classification_spec <- function(input_kind, models, operating_rule,
                                        prevalence = list(kind = "observed", value = NULL, source = NULL),
                                        counts = NULL, theory = NULL, data_fingerprint = NULL,
                                        uncertainty = list(method = NULL, level = NULL,
                                                           seed = NULL, replicates = NULL)) {
  x <- list(schema_version = .ig_schema_versions[["classification_spec"]],
            object_type = "classification_spec",
            input_kind = input_kind, models = models, operating_rule = operating_rule,
            prevalence = prevalence, counts = counts, theory = theory,
            data_fingerprint = data_fingerprint, uncertainty = uncertainty)
  problems <- .ig_classification_spec_problems(x)
  if (length(problems)) .ig_abort_invalid(problems, "invalid classification specification")
  x <- .ig_doubles(.ig_sort_names(x))
  x$spec_hash <- .ig_hash_value(x)
  structure(x[order(names(x), method = "radix")], class = c("classification_spec", "ig_frozen"))
}

.ig_classification_spec_problems <- function(x) {
  v <- .ig_classification_vocab
  p <- character()
  enum_ok <- function(val, allowed) .ig_is_string(val) && val %in% allowed
  if (!enum_ok(x$input_kind, v$input_kind)) {
    p <- c(p, sprintf("input_kind must be one of: %s", paste(v$input_kind, collapse = ", ")))
  }
  if (!is.list(x$models) || !length(x$models) || !is.null(names(x$models))) {
    p <- c(p, "models must be a non-empty list of model records")
  } else {
    ids <- character()
    for (i in seq_along(x$models)) {
      m <- x$models[[i]]
      where <- sprintf("models[%d]", i)
      for (f in c("model_id", "name", "score_meaning", "positive_label", "tie_rule",
                  "missing_value_policy")) {
        if (!.ig_is_string(m[[f]])) p <- c(p, sprintf("%s$%s must be a string", where, f))
      }
      if (!enum_ok(m$direction, v$direction)) p <- c(p, sprintf("%s$direction must be one of: %s", where, paste(v$direction, collapse = ", ")))
      if (!enum_ok(m$comparator, v$comparator)) p <- c(p, sprintf("%s$comparator must be one of: %s", where, paste(v$comparator, collapse = ", ")))
      if (!enum_ok(m$score_scale, v$score_scale)) p <- c(p, sprintf("%s$score_scale must be one of: %s", where, paste(v$score_scale, collapse = ", ")))
      if (!enum_ok(m$evaluation_split, v$evaluation_split)) p <- c(p, sprintf("%s$evaluation_split must be one of: %s", where, paste(v$evaluation_split, collapse = ", ")))
      if (!is.null(m$source_hash) && !(.ig_is_string(m$source_hash) && grepl("^[0-9a-f]{64}$", m$source_hash))) {
        p <- c(p, sprintf("%s$source_hash must be a SHA-256 hex string or NULL", where))
      }
      ids <- c(ids, if (.ig_is_string(m$model_id)) m$model_id else NA_character_)
    }
    if (anyDuplicated(stats::na.omit(ids))) p <- c(p, "model_id values must be unique")
  }
  r <- x$operating_rule
  if (!is.list(r)) {
    p <- c(p, "operating_rule must be a list")
  } else {
    if (!is.null(r$score_cutoff) && !(is.numeric(r$score_cutoff) && length(r$score_cutoff) == 1L &&
                                       is.finite(r$score_cutoff))) {
      p <- c(p, "operating_rule$score_cutoff must be a finite number or NULL")
    }
    if (!enum_ok(r$selection, v$selection)) p <- c(p, sprintf("operating_rule$selection must be one of: %s", paste(v$selection, collapse = ", ")))
    if (!enum_ok(r$selection_data, v$selection_data)) p <- c(p, sprintf("operating_rule$selection_data must be one of: %s", paste(v$selection_data, collapse = ", ")))
    if (identical(r$selection, "manual") && is.null(r$score_cutoff)) {
      p <- c(p, "a manual selection needs operating_rule$score_cutoff")
    }
  }
  pr <- x$prevalence
  if (!is.list(pr) || !enum_ok(pr$kind, v$prevalence_kind)) {
    p <- c(p, "prevalence$kind must be 'observed' or 'scenario'")
  } else if (!is.null(pr$value) && !(is.numeric(pr$value) && length(pr$value) == 1L &&
                                       is.finite(pr$value) && pr$value >= 0 && pr$value <= 1)) {
    p <- c(p, "prevalence$value must be a probability or NULL")
  } else if (identical(pr$kind, "scenario") && is.null(pr$value)) {
    p <- c(p, "a prevalence scenario needs prevalence$value")
  }
  if (identical(x$input_kind, "counts_2x2")) {
    cnt <- x$counts
    ok <- is.list(cnt) && all(c("tp", "fn", "fp", "tn") %in% names(cnt)) &&
      all(vapply(cnt[c("tp", "fn", "fp", "tn")], function(z) {
        is.numeric(z) && length(z) == 1L && is.finite(z) && z >= 0 && z == round(z)
      }, logical(1)))
    if (!ok) p <- c(p, "counts_2x2 input needs counts tp, fn, fp, tn as non-negative whole numbers")
  }
  if (identical(x$input_kind, "theoretical_binormal")) {
    th <- x$theory
    ok <- is.list(th) && is.numeric(th$auc) && length(th$auc) == 1L && is.finite(th$auc) &&
      th$auc > 0 && th$auc < 1 && is.numeric(th$sigma) && length(th$sigma) == 1L &&
      is.finite(th$sigma) && th$sigma > 0
    if (!ok) p <- c(p, "theoretical_binormal input needs theory$auc in (0, 1) and theory$sigma > 0")
  }
  if (identical(x$input_kind, "scored_data") && !.ig_is_string(x$data_fingerprint)) {
    p <- c(p, "scored_data input needs a data_fingerprint (hash of the local data; no raw rows)")
  }
  p
}

# Observed-test specification shape (draft; phase 09A owns the constructor).
.ig_new_observed_test_spec <- function(variant_id, inputs = list(), options = list(),
                                       data_fingerprint = NULL) {
  p <- character()
  if (!.ig_is_string(variant_id) || !grepl("^[a-z0-9]+(_[a-z0-9]+)*$", variant_id)) {
    p <- c(p, "variant_id must be a lower-case identifier")
  }
  if (!is.list(inputs) || (length(inputs) && !.ig_is_named(inputs))) p <- c(p, "inputs must be a named list")
  if (!is.list(options) || (length(options) && !.ig_is_named(options))) p <- c(p, "options must be a named list")
  if (!is.null(data_fingerprint) && !.ig_is_string(data_fingerprint)) p <- c(p, "data_fingerprint must be a string or NULL")
  if (length(p)) .ig_abort_invalid(p, "invalid observed-test specification")
  x <- list(schema_version = .ig_schema_versions[["observed_test_spec"]],
            object_type = "observed_test_spec", variant_id = variant_id,
            inputs = .ig_named_or_empty(.ig_sort_names(inputs)),
            options = .ig_named_or_empty(.ig_sort_names(options)),
            data_fingerprint = data_fingerprint)
  x <- .ig_doubles(.ig_sort_names(x))
  x$spec_hash <- .ig_hash_value(x)
  structure(x[order(names(x), method = "radix")], class = c("observed_test_spec", "ig_frozen"))
}

.ig_check_aux_spec <- function(spec, kind) {
  body <- .ig_unclass_frozen(spec)
  stored <- body$spec_hash
  body$spec_hash <- NULL
  if (!.ig_is_string(stored) || !identical(.ig_hash_value(body), stored)) {
    .ig_abort("invalid_input", sprintf("the %s was modified after creation (hash mismatch)", class(spec)[1]))
  }
  if (kind == "classification") {
    problems <- .ig_classification_spec_problems(body)
    if (length(problems)) .ig_abort_invalid(problems, "invalid classification specification")
  }
  invisible(TRUE)
}

.ig_aux_spec_from_list <- function(x, kind) {
  if (!is.list(x)) .ig_abort("invalid_input", "spec must be an object")
  cls <- if (kind == "classification") "classification_spec" else "observed_test_spec"
  if (!identical(x$object_type, cls) ||
      !identical(x$schema_version, .ig_schema_versions[[cls]])) {
    .ig_abort("invalid_input", sprintf("not a supported %s document", cls))
  }
  # Parsed numbers come back as integers where whole; the hash is computed on
  # canonical JSON, which does not distinguish them. Empty objects become
  # empty lists, as in the constructors.
  body <- .ig_doubles(.ig_sort_names(.ig_simplify_parsed(x)))
  spec <- structure(body, class = c(cls, "ig_frozen"))
  .ig_check_aux_spec(spec, kind)
  spec
}
