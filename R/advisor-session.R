# Advisor sessions ------------------------------------------------------------------
#
# An advisor_session is an immutable value: the ordered raw answers, an
# optional method override and the session language, tied to one rules
# version. Every operation returns a new session; the walk (current question
# or outcome, effective answers, transition trace) is recomputed from the
# answers by the rule interpreter, so a session can never hold a state the
# rules do not produce.

.ig_advisor_session_schema <- "inferguidR.advisor_session/0.1.0"
.ig_advisor_languages <- c("en", "de")

# Typed advisor condition: status `code`, optional extra condition class
# `inferguidR_<subclass>` (e.g. override_refused, migration_required).
.ig_advisor_abort <- function(code, reason, details = list(), subclass = NULL) {
  status <- ig_status(code, reason, details)
  cond <- structure(
    class = c(if (!is.null(subclass)) paste0("inferguidR_", subclass),
              paste0("inferguidR_", code), "inferguidR_error", "error", "condition"),
    list(message = reason, call = NULL, status = status)
  )
  stop(cond)
}

.ig_advisor_check_language <- function(language) {
  if (!.ig_is_string(language) || !(language %in% .ig_advisor_languages)) {
    .ig_abort("invalid_input", sprintf("language must be one of: %s",
                                       paste(.ig_advisor_languages, collapse = ", ")))
  }
  language
}

.ig_advisor_check_session <- function(session) {
  if (!inherits(session, "advisor_session")) {
    .ig_abort("invalid_input", "`session` must be an advisor_session (see advisor_session())")
  }
  invisible(session)
}

# Resolve the rule set for a session and check that it is the one the session
# was built with (same id, version and content hash).
.ig_advisor_session_rules <- function(session, rules = NULL) {
  rs <- if (is.null(rules)) advisor_rules() else rules
  if (!inherits(rs, "advisor_rules")) {
    .ig_abort("invalid_input", "`rules` must be an advisor_rules object (see advisor_rules()) or NULL")
  }
  if (is.null(session)) return(rs)
  .ig_advisor_check_rules_identity(session$rules, rs)
  rs
}

# Compare a recorded rules reference with a loaded rule set.
.ig_advisor_check_rules_identity <- function(ref, rs) {
  v <- rs$versions
  if (!identical(ref$rules_id, v$rules_id)) {
    .ig_abort("invalid_input", sprintf("the session belongs to rule set '%s', not '%s'",
                                       format(ref$rules_id), v$rules_id))
  }
  if (!identical(ref$rules_version, v$rules_version)) {
    .ig_advisor_migration_abort(ref$rules_version, rs)
  }
  if (!identical(ref$rules_sha256, v$rules_sha256)) {
    .ig_abort("invalid_input", sprintf(
      "rules version %s has different content from the one the session was created with (content hash mismatch)",
      v$rules_version))
  }
  invisible(TRUE)
}

.ig_advisor_migration_abort <- function(saved_version, rs) {
  log <- rs$docs$rules$change_log %||% list()
  entries <- lapply(log, function(e) list(rules_version = e$rules_version %||% "",
                                          summary_en = e$summary_en %||% "",
                                          migration_en = e$migration_en %||% ""))
  .ig_advisor_abort(
    "invalid_input",
    sprintf(paste0("migration_required: the saved state uses rules version %s, this installation has %s; ",
                   "it is not re-interpreted silently and needs an explicit migration"),
            format(saved_version), rs$versions$rules_version),
    list(reason_code = "migration_required",
         saved_rules_version = format(saved_version),
         current_rules_version = rs$versions$rules_version,
         change_log = entries),
    subclass = "migration_required")
}

# Check an override variant against the outcome a walk reached. Returns NULL
# when allowed, otherwise a refusal record.
.ig_advisor_override_refusal <- function(rs, walk, variant_id) {
  if (!.ig_is_string(variant_id)) {
    return(list(refusal = "unknown_variant", reason_en = "analysis_variant_id must be a single string",
                reason_de = "analysis_variant_id muss eine einzelne Zeichenkette sein"))
  }
  if (is.null(rs$variants[[variant_id]])) {
    return(list(refusal = "unknown_variant",
                reason_en = sprintf("'%s' is not an analysis variant of rules %s", variant_id,
                                    rs$versions$rules_version),
                reason_de = sprintf("'%s' ist keine Analysevariante der Regeln %s", variant_id,
                                    rs$versions$rules_version)))
  }
  outcome <- .ig_advisor_outcome_of(rs, walk)
  if (!is.null(walk$pending) || !(outcome$status %in% c("recommended", "conditional"))) {
    return(list(refusal = "no_method",
                reason_en = "the walkthrough has not reached a method outcome, so there is no variant to replace",
                reason_de = "Der Ablauf hat kein Methodenergebnis erreicht; es gibt keine Variante, die ersetzt werden kann"))
  }
  if (identical(variant_id, "cochran_armitage_trend") &&
      !identical(unname(walk$effective["q.outcome_type"]), "binary")) {
    return(list(refusal = "incompatible_outcome_type",
                reason_en = "The ordered-proportion trend test requires a recorded binary outcome.",
                reason_de = "Der Trendtest fuer Anteile erfordert eine dokumentierte binaere Zielgroesse."))
  }
  if (identical(variant_id, "variance_f_test_two") &&
      !identical(unname(walk$effective["q.n_groups"]), "two")) {
    return(list(refusal = "incompatible_group_count",
                reason_en = "The F variance-ratio test requires exactly two independent groups.",
                reason_de = "Der F-Test des Varianzverhaeltnisses erfordert genau zwei unabhaengige Gruppen."))
  }
  # The outcome's own variant needs no override (the session keeps none).
  if (identical(variant_id, outcome$analysis_variant_id)) return(NULL)
  for (r in outcome$rejected_routes %||% list()) {
    if (identical(r$analysis_variant_id, variant_id)) {
      return(list(refusal = "rejected_route", reason_en = r$reason_en, reason_de = r$reason_de))
    }
  }
  for (a in outcome$alternatives %||% list()) {
    if (identical(a$analysis_variant_id, variant_id)) {
      if (isTRUE(a$override_allowed)) return(NULL)
      return(list(refusal = "override_not_allowed", reason_en = a$why_not_preferred_en,
                  reason_de = a$why_not_preferred_de))
    }
  }
  list(refusal = "not_listed",
       reason_en = sprintf("'%s' is not listed as an alternative of %s and is not compatible with the recorded design",
                           variant_id, outcome$outcome_id),
       reason_de = sprintf("'%s' ist keine Alternative zu %s und passt nicht zum erfassten Design",
                           variant_id, outcome$outcome_id))
}

# Internal constructor: replays `answers` and builds the session value.
.ig_advisor_new_session <- function(rs, answers, language, override = NULL,
                                    ignored_answers = character(),
                                    invalidated = list(question_ids = character(), override = NULL)) {
  .ig_advisor_check_language(language)
  answers <- lapply(answers, function(a) list(question_id = a$question_id, option_id = a$option_id))
  walk <- .ig_advisor_walk(rs, answers)
  if (!is.null(override)) {
    ref <- .ig_advisor_override_refusal(rs, walk, override)
    if (!is.null(ref)) {
      .ig_advisor_abort("invalid_input", sprintf("override '%s' refused (%s): %s", override, ref$refusal,
                                                 ref$reason_en),
                        list(refusal = ref$refusal, analysis_variant_id = override), subclass = "override_refused")
    }
    if (identical(override, .ig_advisor_outcome_of(rs, walk)$analysis_variant_id)) override <- NULL
  }
  v <- rs$versions
  structure(list(
    schema_version = .ig_advisor_session_schema,
    object_type = "advisor_session",
    rules = list(rules_id = v$rules_id, rules_version = v$rules_version, rules_sha256 = v$rules_sha256),
    language = language,
    answers = answers,
    override = override,
    ignored_answers = as.character(ignored_answers),
    invalidated = list(question_ids = as.character(invalidated$question_ids),
                       override = invalidated$override),
    state = list(node = walk$node, outcome_id = walk$outcome_id,
                 pending_question_id = walk$pending, steps = walk$steps,
                 effective = walk$effective, raw = walk$raw)
  ), class = c("advisor_session", "ig_frozen"))
}

#' Start a guided advisor walkthrough
#'
#' Creates an empty advisor session at the first question of the rule graph.
#' Sessions are immutable: [advisor_answer()], [advisor_back()],
#' [advisor_override()] and [advisor_set_language()] return new sessions.
#' Only questions that apply to the answers given so far are ever asked;
#' 'I do not know' is an answer on every question and never invents a
#' default (its policy either stops with `needs_information`, routes without
#' an extra assumption, or routes with a recorded assumption and a downgrade
#' to `conditional`).
#'
#' @param language Language of question texts: `"en"` or `"de"`. Decisions do
#'   not depend on the language.
#' @param rules An `advisor_rules` object, or `NULL` for the bundled rules.
#' @return An object of class `advisor_session`.
#' @examples
#' s <- advisor_session()
#' s <- advisor_answer(s, "compare")
#' advisor_question(s)$question_id
#' @seealso [advisor_answer()], [advisor_recommendation()], [advisor_save()]
#' @family advisor
#' @export
advisor_session <- function(language = "en", rules = NULL) {
  rs <- .ig_advisor_session_rules(NULL, rules)
  .ig_advisor_new_session(rs, list(), language)
}

#' Answer, change or retract advisor questions
#'
#' `advisor_answer()` answers the current question, or changes the answer to
#' a question answered earlier. Changing an earlier answer is backtracking:
#' every answer recorded after it, and any override, are invalidated (listed
#' in `session$invalidated`) and the walkthrough continues from the changed
#' question. Giving the same answer again changes nothing. A question that is
#' neither current nor answered is refused, so only applicable follow-ups can
#' be answered.
#'
#' `advisor_back()` retracts the last answer, or every answer from question
#' `to` onwards, and makes that question current again.
#'
#' @param session An `advisor_session`.
#' @param answer An option id of the question (see [advisor_question()]),
#'   for example `"dont_know"`.
#' @param question_id The question to answer. `NULL` (default) answers the
#'   current question.
#' @param to Question id to go back to (it must have been answered). `NULL`
#'   retracts only the last answer.
#' @param rules An `advisor_rules` object, or `NULL` for the bundled rules. It
#'   must be the rule set the session was created with.
#' @return A new `advisor_session`. Invalid answers raise a condition of class
#'   `inferguidR_invalid_input`.
#' @examples
#' s <- advisor_session()
#' s <- advisor_answer(s, "compare")
#' s <- advisor_answer(s, "analysis")
#' s <- advisor_answer(s, "describe", question_id = "q.aim") # backtracking
#' s$invalidated$question_ids
#' @family advisor
#' @export
advisor_answer <- function(session, answer, question_id = NULL, rules = NULL) {
  .ig_advisor_check_session(session)
  rs <- .ig_advisor_session_rules(session, rules)
  st <- session$state
  qid <- question_id %||% st$pending_question_id
  if (is.null(qid)) {
    .ig_abort("invalid_input", paste0("the walkthrough has reached an outcome; give `question_id` to change ",
                                      "an earlier answer"))
  }
  if (!.ig_is_string(qid) || is.null(rs$questions[[qid]])) {
    .ig_abort("invalid_input", sprintf("unknown question id '%s'", paste(format(qid), collapse = "")))
  }
  if (!.ig_is_string(answer) || is.null(rs$options[[qid]][[answer]])) {
    .ig_abort("invalid_input", sprintf("'%s' is not an option of question '%s'; options: %s",
                                       paste(format(answer), collapse = ""), qid,
                                       paste(names(rs$options[[qid]]), collapse = ", ")))
  }
  answered <- vapply(session$answers, `[[`, "", "question_id")
  pos <- match(qid, answered)
  invalidated <- list(question_ids = character(), override = NULL)
  if (!is.na(pos)) {
    if (identical(session$answers[[pos]]$option_id, answer)) return(session)
    invalidated <- list(question_ids = answered[seq_along(answered) > pos], override = session$override)
    new_answers <- c(session$answers[seq_len(pos - 1L)],
                     list(list(question_id = qid, option_id = answer)))
  } else if (identical(qid, st$pending_question_id)) {
    new_answers <- c(session$answers, list(list(question_id = qid, option_id = answer)))
  } else {
    .ig_abort("invalid_input", sprintf(
      "question '%s' does not apply at this point of the walkthrough (current question: %s)",
      qid, st$pending_question_id %||% "none, an outcome was reached"),
      list(current_question_id = st$pending_question_id %||% "",
           answered_question_ids = as.list(answered)))
  }
  .ig_advisor_new_session(rs, new_answers, session$language, override = NULL,
                          ignored_answers = session$ignored_answers, invalidated = invalidated)
}

#' @rdname advisor_answer
#' @export
advisor_back <- function(session, to = NULL, rules = NULL) {
  .ig_advisor_check_session(session)
  rs <- .ig_advisor_session_rules(session, rules)
  answered <- vapply(session$answers, `[[`, "", "question_id")
  if (!length(answered)) .ig_abort("invalid_input", "no answer to go back to")
  if (!is.null(to) && !.ig_is_string(to)) {
    .ig_abort("invalid_input", "`to` must be NULL or a single question id")
  }
  pos <- if (is.null(to)) length(answered) else match(to, answered)
  if (is.na(pos)) {
    .ig_abort("invalid_input", sprintf("question '%s' has not been answered in this session",
                                       paste(format(to), collapse = "")))
  }
  .ig_advisor_new_session(rs, session$answers[seq_len(pos - 1L)], session$language,
                          ignored_answers = session$ignored_answers,
                          invalidated = list(question_ids = answered[pos:length(answered)],
                                             override = session$override))
}

#' Replace the recommended variant by a listed alternative
#'
#' A user may replace the variant of the outcome only by a variant listed in
#' its alternatives with `override_allowed = TRUE`. The override is recorded
#' with the alternative's why-not text and the status becomes at most
#' `conditional`. Variants listed as rejected routes are refused with the
#' stated reason; any other variant is refused as not compatible with the
#' recorded design (`refusal = "not_listed"`). Listed alternatives must also
#' satisfy their design guards: the F variance-ratio test requires two groups,
#' and the ordered-proportion trend test requires a recorded binary outcome.
#'
#' @param session An `advisor_session` that has reached a method outcome.
#' @param analysis_variant_id The variant to use instead.
#' @inheritParams advisor_answer
#' @return A new `advisor_session`. A refusal is a condition of classes
#'   `inferguidR_override_refused` and `inferguidR_invalid_input`; its
#'   `$status$details$refusal` is one of `"rejected_route"`, `"not_listed"`,
#'   `"override_not_allowed"`, `"no_method"`, `"unknown_variant"`,
#'   `"incompatible_group_count"` or `"incompatible_outcome_type"`.
#' @examples
#' s <- advisor_replay(c(q.aim = "compare", q.stage = "analysis",
#'   q.outcome_type = "continuous", q.unit = "individuals",
#'   q.structure = "independent", q.independence_check = "no",
#'   q.n_factors = "one", q.n_groups = "two", q.covariates = "no",
#'   q.estimand_continuous = "mean", q.variance = "not_justified",
#'   q.direction = "two_sided", q.missingness = "none"))
#' s2 <- advisor_override(s, "t_means_two_independent_student")
#' advisor_recommendation(s2)$status
#' @family advisor
#' @export
advisor_override <- function(session, analysis_variant_id, rules = NULL) {
  .ig_advisor_check_session(session)
  rs <- .ig_advisor_session_rules(session, rules)
  walk <- .ig_advisor_walk(rs, session$answers)
  ref <- .ig_advisor_override_refusal(rs, walk, analysis_variant_id)
  if (!is.null(ref)) {
    .ig_advisor_abort(
      "invalid_input",
      sprintf("override '%s' refused (%s): %s", paste(format(analysis_variant_id), collapse = ""),
              ref$refusal, ref$reason_en),
      list(refusal = ref$refusal, analysis_variant_id = paste(format(analysis_variant_id), collapse = ""),
           outcome_id = walk$outcome_id %||% "", reason_en = ref$reason_en %||% "",
           reason_de = ref$reason_de %||% ""),
      subclass = "override_refused")
  }
  .ig_advisor_new_session(rs, session$answers, session$language, override = analysis_variant_id,
                          ignored_answers = session$ignored_answers)
}

#' @rdname advisor_override
#' @export
advisor_clear_override <- function(session, rules = NULL) {
  .ig_advisor_check_session(session)
  rs <- .ig_advisor_session_rules(session, rules)
  .ig_advisor_new_session(rs, session$answers, session$language,
                          ignored_answers = session$ignored_answers,
                          invalidated = list(question_ids = character(), override = session$override))
}

#' Change the language of an advisor session
#'
#' Switches question texts between English and German. Answers, outcome and
#' status are unchanged.
#'
#' @inheritParams advisor_answer
#' @param language `"en"` or `"de"`.
#' @return A new `advisor_session`.
#' @examples
#' advisor_question(advisor_set_language(advisor_session(), "de"))$text
#' @family advisor
#' @export
advisor_set_language <- function(session, language, rules = NULL) {
  .ig_advisor_check_session(session)
  rs <- .ig_advisor_session_rules(session, rules)
  .ig_advisor_new_session(rs, session$answers, language, override = session$override,
                          ignored_answers = session$ignored_answers, invalidated = session$invalidated)
}

#' Replay an answer map (deep link)
#'
#' Rebuilds a session from a map of question ids to option ids, for example
#' from a deep link. The map is replayed through the graph in transition
#' order: an unknown question or option id makes the whole map invalid
#' (`invalid_input`), answers to questions that the path does not reach are
#' reported in `session$ignored_answers` and do not influence the result, and
#' a map that stops at an applicable but unanswered question leaves that
#' question current (the recommendation is then `needs_information` with
#' that pending question).
#'
#' @param answers A named character vector or named list: question id ->
#'   option id.
#' @inheritParams advisor_session
#' @return An `advisor_session`.
#' @examples
#' s <- advisor_replay(c(q.aim = "compare", q.stage = "analysis",
#'   q.outcome_type = "continuous", q.variance = "justified"))
#' s$ignored_answers
#' advisor_question(s)$question_id
#' @family advisor
#' @export
advisor_replay <- function(answers, language = "en", rules = NULL) {
  rs <- .ig_advisor_session_rules(NULL, rules)
  if (is.list(answers)) {
    ok <- vapply(answers, .ig_is_string, logical(1))
    if (!all(ok)) .ig_abort("invalid_input", "every answer in the map must be a single option id")
    answers <- vapply(answers, identity, "")
  }
  if (length(answers) && (!is.character(answers) || !.ig_is_named(answers) || anyNA(answers))) {
    .ig_abort("invalid_input", "`answers` must be a named character vector (question id -> option id)")
  }
  problems <- character()
  dup <- unique(names(answers)[duplicated(names(answers))])
  if (length(dup)) problems <- c(problems, sprintf("question '%s' is given more than once", dup))
  for (qid in unique(names(answers))) {
    if (is.null(rs$questions[[qid]])) {
      problems <- c(problems, sprintf("unknown question id '%s'", qid))
    } else if (is.null(rs$options[[qid]][[answers[[qid]]]])) {
      problems <- c(problems, sprintf("'%s' is not an option of question '%s'", answers[[qid]], qid))
    }
  }
  if (length(problems)) .ig_abort_invalid(problems, "invalid answer map")
  seq_answers <- list()
  repeat {
    walk <- .ig_advisor_walk(rs, seq_answers)
    nxt <- walk$pending
    if (is.null(nxt) || !(nxt %in% names(answers))) break
    seq_answers[[length(seq_answers) + 1L]] <- list(question_id = nxt, option_id = answers[[nxt]])
  }
  used <- vapply(seq_answers, `[[`, "", "question_id")
  .ig_advisor_new_session(rs, seq_answers, language,
                          ignored_answers = setdiff(names(answers), used))
}

#' Show an advisor question
#'
#' Returns the text of a question in the session language: the question,
#' why it matters, the options (with the worked explanation of 'I do not
#' know' and the assumption it records, if any), the confirmation rules that
#' apply when data are supplied, and the recorded answer.
#'
#' @inheritParams advisor_answer
#' @param question_id A question id; `NULL` (default) gives the current
#'   question, or `NULL` when an outcome has been reached.
#' @param language `"en"` or `"de"`; `NULL` uses the session language.
#' @return An object of class `advisor_question`, or `NULL`.
#' @examples
#' q <- advisor_question(advisor_session(language = "de"))
#' q$text
#' q$options[[length(q$options)]]$explanation
#' @family advisor
#' @export
advisor_question <- function(session, question_id = NULL, language = NULL, rules = NULL) {
  .ig_advisor_check_session(session)
  rs <- .ig_advisor_session_rules(session, rules)
  lang <- .ig_advisor_check_language(language %||% session$language)
  qid <- question_id %||% session$state$pending_question_id
  if (is.null(qid)) return(NULL)
  # Only a single id string is looked up: `[[` with a number would silently
  # select a question by position, and with a vector would index recursively.
  q <- if (.ig_is_string(qid)) rs$questions[[qid]]
  if (is.null(q)) .ig_abort("invalid_input", sprintf("unknown question id '%s'", paste(format(qid), collapse = " ")))
  txt <- function(x, field) x[[paste0(field, "_", lang)]]
  options <- lapply(q$options, function(o) {
    list(option_id = o$option_id, label = txt(o, "label"),
         dont_know = identical(o$option_id, "dont_know"),
         dont_know_policy = o$dont_know_policy, routes_as = o$routes_as,
         explanation = txt(o, "explanation"), assumption = txt(o, "assumption"))
  })
  answered <- vapply(session$answers, `[[`, "", "question_id")
  conf <- lapply(q$confirmation_rules %||% list(), function(cr) {
    list(rule_id = cr$rule_id, applies_when = cr$applies_when, action = txt(cr, "action"))
  })
  structure(list(
    question_id = qid, language = lang,
    text = txt(q, "text"), why_it_matters = txt(q, "why_it_matters"),
    options = options,
    confirmation_rules = conf,
    is_current = identical(qid, session$state$pending_question_id),
    on_path = qid %in% answered,
    answer = .ig_advisor_lookup(session$state$raw, qid),
    effective_answer = .ig_advisor_lookup(session$state$effective, qid)
  ), class = "advisor_question")
}

# Value of a named character vector, or NULL when the name is absent.
.ig_advisor_lookup <- function(x, key) {
  if (!(key %in% names(x))) return(NULL)
  unname(x[[key]])
}

#' @export
format.advisor_question <- function(x, ...) {
  opts <- vapply(x$options, function(o) sprintf("  [%s] %s", o$option_id, o$label), "")
  c(sprintf("<advisor_question> %s (%s)", x$question_id, x$language), paste0("  ", x$text), opts)
}

#' @export
print.advisor_question <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}

#' @export
format.advisor_session <- function(x, ...) {
  st <- x$state
  ans <- vapply(st$steps, function(s) {
    if (identical(s$raw_answer, s$effective_answer)) sprintf("%s = %s", s$question_id, s$raw_answer)
    else sprintf("%s = %s (routes as %s)", s$question_id, s$raw_answer, s$effective_answer)
  }, "")
  c(sprintf("<advisor_session> rules %s %s, language %s", x$rules$rules_id, x$rules$rules_version, x$language),
    if (length(ans)) paste0("  ", ans) else "  (no answers yet)",
    if (!is.null(st$pending_question_id)) sprintf("  current question: %s", st$pending_question_id)
    else sprintf("  outcome: %s", st$outcome_id),
    if (!is.null(x$override)) sprintf("  override: %s", x$override),
    if (length(x$ignored_answers)) sprintf("  ignored answers: %s", paste(x$ignored_answers, collapse = ", ")))
}

#' @export
print.advisor_session <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}
