# Saving and restoring advisor sessions -----------------------------------------
#
# A saved state records the session schema version, the rules id, version
# and content hash, the language, the raw answer trace, the override, the
# ignored deep-link answers, the last invalidation, and the outcome and
# status it produced, sealed by a SHA-256 state hash over the canonical JSON.
#
# Restoring (rules.json engine.restore_policy_en):
# * the document is parsed as data only (no evaluation, duplicate keys
#   rejected) and must have the exact session schema version;
# * a different rules version raises `migration_required` (a condition of
#   class inferguidR_migration_required): the state is never re-interpreted
#   silently; the rule change log is attached for an explicit migration;
# * with the same rules version the trace is replayed and must reproduce the
#   saved outcome and status exactly.

.ig_advisor_state_type <- "advisor_session_state"

.ig_advisor_state_body <- function(session, rs) {
  rec_status <- .ig_advisor_session_status(session, rs)
  list(
    object_type = .ig_advisor_state_type,
    schema_version = .ig_advisor_session_schema,
    rules = session$rules,
    language = session$language,
    answers = lapply(session$answers, function(a) list(question_id = a$question_id, option_id = a$option_id)),
    override = session$override,
    ignored_answers = I(as.character(session$ignored_answers)),
    invalidated = list(question_ids = I(as.character(session$invalidated$question_ids)),
                       override = session$invalidated$override),
    expected = list(node = session$state$node, status = rec_status)
  )
}

# Final status of a session (same status rules as the recommendation).
.ig_advisor_session_status <- function(session, rs) {
  walk <- session$state
  walk$pending <- walk$pending_question_id
  .ig_advisor_status(rs, .ig_advisor_outcome_of(rs, walk), walk, session$override)$status
}

#' Save and restore an advisor session
#'
#' `advisor_save()` writes the session state as canonical JSON: session
#' schema version, rules id, version and content hash, language, raw answer
#' trace, override, ignored deep-link answers and the outcome and status the
#' state produced, sealed by a SHA-256 state hash.
#'
#' `advisor_restore()` reads such a document as data (nothing is evaluated).
#' With the same rules version the answers are replayed and must reproduce
#' the saved outcome and status. A state saved under a different rules
#' version is not re-interpreted: a condition of class
#' `inferguidR_migration_required` (status `invalid_input`) is raised, with the
#' rule change log in `$status$details$change_log`, so that an explicit
#' migration can map the answers and report changes.
#'
#' @inheritParams advisor_answer
#' @param file Path of a local file to write (`advisor_save()`) or read
#'   (`advisor_restore()`).
#' @param json A JSON string produced by `advisor_save()`.
#' @return `advisor_save()` returns the JSON text (class `json`), invisibly
#'   when `file` is given. `advisor_restore()` returns an `advisor_session`.
#' @examples
#' s <- advisor_answer(advisor_session(), "compare")
#' js <- advisor_save(s)
#' identical(advisor_restore(js), s)
#' @family advisor
#' @export
advisor_save <- function(session, file = NULL, rules = NULL) {
  .ig_advisor_check_session(session)
  rs <- .ig_advisor_session_rules(session, rules)
  body <- .ig_advisor_state_body(session, rs)
  body$state_hash <- .ig_hash_value(body)
  json <- ig_to_json(body)
  if (!is.null(file)) {
    if (!.ig_is_string(file) || grepl("^[A-Za-z][A-Za-z0-9+.-]*://", file)) {
      .ig_abort("invalid_input", "`file` must be a single local path")
    }
    if (dir.exists(file)) .ig_abort("invalid_input", paste("`file` is a directory, not a file:", file))
    fail <- function(e) {
      .ig_abort("invalid_input", sprintf("the session could not be written to %s (%s)", file, conditionMessage(e)))
    }
    con <- tryCatch(file(file, open = "wb"), error = fail, warning = fail)
    on.exit(close(con), add = TRUE)
    writeBin(charToRaw(enc2utf8(unclass(json))), con)
    return(invisible(json))
  }
  json
}

#' @rdname advisor_save
#' @export
advisor_restore <- function(json = NULL, file = NULL, rules = NULL) {
  if (is.null(json) == is.null(file)) .ig_abort("invalid_input", "supply exactly one of `json` or `file`")
  if (!is.null(file)) json <- .ig_advisor_read_text(file)
  if (inherits(json, "json")) json <- unclass(json)
  doc <- .ig_json_parse(json)
  if (!is.list(doc) || is.null(names(doc)) || !identical(doc$object_type, .ig_advisor_state_type)) {
    .ig_abort("invalid_input", "not a saved advisor session (object_type must be 'advisor_session_state')")
  }
  if (!identical(doc$schema_version, .ig_advisor_session_schema)) {
    .ig_abort("invalid_input", sprintf(
      "unsupported advisor session schema_version '%s'; this version of inferguidR reads '%s'",
      paste(format(doc$schema_version), collapse = ""), .ig_advisor_session_schema))
  }
  expected_keys <- c("answers", "expected", "ignored_answers", "invalidated", "language", "object_type",
                     "override", "rules", "schema_version", "state_hash")
  unknown <- setdiff(names(doc), expected_keys)
  if (length(unknown)) {
    .ig_abort("invalid_input", sprintf("unknown field '%s' in the saved session", unknown[1]))
  }
  absent <- setdiff(expected_keys, names(doc))
  if (length(absent)) {
    .ig_abort("invalid_input", sprintf("the saved session lacks the field(s): %s", paste(absent, collapse = ", ")))
  }
  body <- doc
  body$state_hash <- NULL
  if (!.ig_is_string(doc$state_hash) || !identical(.ig_hash_value(body), doc$state_hash)) {
    .ig_abort("invalid_input", "state_hash does not match the content (the saved session was modified)")
  }
  rs <- if (is.null(rules)) advisor_rules() else rules
  if (!inherits(rs, "advisor_rules")) .ig_abort("invalid_input", "`rules` must be an advisor_rules object or NULL")
  ref <- doc$rules
  if (!is.list(ref) || !.ig_is_string(ref$rules_id) || !.ig_is_string(ref$rules_version) ||
      !.ig_is_string(ref$rules_sha256)) {
    .ig_abort("invalid_input", "the saved session lacks rules_id, rules_version or rules_sha256")
  }
  .ig_advisor_check_rules_identity(ref, rs)
  .ig_advisor_check_state_shape(doc)
  chr <- function(z) as.character(.ig_flatten(z, "character"))
  inv <- doc$invalidated
  session <- tryCatch(
    .ig_advisor_new_session(rs, doc$answers, doc$language, override = doc$override,
                            ignored_answers = chr(doc$ignored_answers),
                            invalidated = list(question_ids = chr(inv$question_ids), override = inv$override)),
    inferguidR_override_refused = function(e) {
      .ig_abort("invalid_input", paste("the saved override is not allowed by these rules:", conditionMessage(e)))
    })
  exp <- doc$expected
  got <- list(node = session$state$node, status = .ig_advisor_session_status(session, rs))
  if (!identical(exp$node, got$node) || !identical(exp$status, got$status)) {
    .ig_abort("invalid_input", sprintf(
      "replaying the saved answers gives %s (%s), not the saved %s (%s)",
      got$node, got$status, exp$node, exp$status))
  }
  session
}

# Shape of the fields of a saved state (after the hash and rules checks), so
# that a resealed but malformed document is refused with a typed error.
.ig_advisor_check_state_shape <- function(doc) {
  string_array <- function(z) is.list(z) && is.null(names(z)) && all(vapply(z, .ig_is_string, logical(1)))
  string_or_null <- function(z) is.null(z) || .ig_is_string(z)
  p <- character()
  ans_ok <- is.list(doc$answers) && is.null(names(doc$answers)) &&
    all(vapply(doc$answers, function(a) is.list(a) && .ig_is_string(a$question_id) && .ig_is_string(a$option_id),
               logical(1)))
  if (!ans_ok) p <- c(p, "answers must be a list of {question_id, option_id} records")
  if (!.ig_is_string(doc$language)) p <- c(p, "language must be a string")
  if (!string_or_null(doc$override)) p <- c(p, "override must be an analysis variant id or null")
  if (!string_array(doc$ignored_answers)) p <- c(p, "ignored_answers must be an array of question ids")
  inv <- doc$invalidated
  if (!is.list(inv) || !string_array(inv$question_ids) || !string_or_null(inv$override)) {
    p <- c(p, "invalidated must hold question_ids (array of question ids) and override (id or null)")
  }
  exp <- doc$expected
  if (!is.list(exp) || !.ig_is_string(exp$node) || !.ig_is_string(exp$status)) {
    p <- c(p, "expected must hold the node and status the state produced")
  }
  if (length(p)) .ig_abort_invalid(p, "malformed saved advisor session")
  invisible(TRUE)
}

.ig_advisor_read_text <- function(file) {
  if (!.ig_is_string(file) || grepl("^[A-Za-z][A-Za-z0-9+.-]*://", file)) {
    .ig_abort("invalid_input", "`file` must be the path of a local file")
  }
  .ig_advisor_read_local(file, "saved session file")
}
