# Advisor recommendation object -------------------------------------------------------
#
# Shape per the DIGITAL_ADVISOR contract: answered questions, rule trace,
# estimand/design, preferred method and variant, rationale, assumptions
# (met/unknown/not met), alternatives with reasons, missing information,
# R implementation, planning compatibility, learning links, visual preset,
# references and rule versions. A recommendation references specifications
# by hash; it never contains numerical results. The rule engine that fills
# it is phase 12A; the core validates the shape.

.ig_recommendation_vocab <- list(
  status = c("recommended", "conditional", "needs_information", "outside_supported_scope"),
  planning_compatibility = c("exact_gpower_match", "gpower_comparator_different_method",
                             "independently_validated_extension", "not_yet_validated",
                             "not_applicable"),
  planning_availability = c("planned_with_gpower_reference_candidate",
                            "planned_independent_validation", "simulation_only_if_implemented",
                            "unavailable", "not_applicable", "available"),
  assumption_state = c("met", "unknown", "not_met", "not_applicable")
)

# Internal constructor. `planning$specs` may hold power_spec objects; they are
# stored as their hashes in `planning$spec_refs`.
.ig_new_recommendation <- function(status, outcome_id, rules, answers, rule_trace,
                                   rationale, estimand = NULL, design = list(),
                                   preferred = NULL, assumptions = list(),
                                   alternatives = list(), missing_information = character(),
                                   r_implementation = list(template_id = NULL,
                                                           template_status = "not_applicable"),
                                   planning = list(planning_compatibility = "not_applicable",
                                                   planning_availability = "not_applicable"),
                                   learning_links = character(), visual_preset = NULL,
                                   references = character(), next_actions = character()) {
  if (is.list(planning) && !is.null(planning$specs)) {
    specs <- planning$specs
    if (inherits(specs, "power_spec")) specs <- list(specs)
    ok <- vapply(specs, inherits, logical(1), what = "power_spec")
    if (!all(ok)) .ig_abort_invalid("planning$specs must contain power_spec objects only")
    planning$spec_refs <- unique(c(as.character(unlist(planning$spec_refs)),
                                   vapply(specs, `[[`, "", "spec_hash")))
    planning$specs <- NULL
  }
  x <- list(
    schema_version = .ig_schema_versions[["recommendation"]],
    object_type = "recommendation",
    status = status,
    outcome_id = outcome_id,
    rules = rules,
    answers = answers,
    rule_trace = rule_trace,
    estimand = estimand,
    design = .ig_named_or_empty(design),
    preferred = preferred,
    rationale = rationale,
    assumptions = assumptions,
    alternatives = alternatives,
    missing_information = as.character(missing_information),
    r_implementation = r_implementation,
    planning = list(
      planning_compatibility = planning$planning_compatibility,
      planning_availability = planning$planning_availability,
      gpower_method_id = planning$gpower_method_id,
      spec_refs = as.character(unlist(planning$spec_refs)),
      evidence_ids = as.character(unlist(planning$evidence_ids)),
      note_en = planning$note_en,
      note_de = planning$note_de
    ),
    learning_links = as.character(learning_links),
    visual_preset = visual_preset,
    references = as.character(references),
    next_actions = as.character(next_actions)
  )
  problems <- .ig_recommendation_problems(x)
  if (length(problems)) .ig_abort_invalid(problems, "invalid recommendation")
  # Canonical member order inside the free-form records (as in canonical
  # JSON), so that a recommendation read back from JSON is identical.
  for (f in c("rules", "answers", "rule_trace", "estimand", "preferred", "rationale",
              "assumptions", "alternatives", "r_implementation")) {
    x[f] <- list(.ig_sort_names(x[[f]]))
  }
  x$design <- .ig_named_or_empty(.ig_sort_names(x$design))
  x$recommendation_hash <- .ig_hash_value(.ig_recommendation_body(x))
  structure(x, class = c("recommendation", "ig_frozen"))
}

.ig_recommendation_problems <- function(x) {
  v <- .ig_recommendation_vocab
  p <- character()
  str_or_null <- function(z) is.null(z) || .ig_is_string(z)
  if (!.ig_is_string(x$status) || !(x$status %in% v$status)) {
    p <- c(p, sprintf("status must be one of: %s", paste(v$status, collapse = ", ")))
  }
  if (!.ig_is_string(x$outcome_id)) p <- c(p, "outcome_id must be a string")
  r <- x$rules
  if (!is.list(r) || !.ig_is_string(r$rules_id) || !.ig_is_string(r$rules_version)) {
    p <- c(p, "rules must give rules_id and rules_version")
  } else {
    for (f in c("uzh_tree_version", "descrtab2_tree_version", "analysis_variants_version")) {
      if (!str_or_null(r[[f]])) p <- c(p, sprintf("rules$%s must be a string or NULL", f))
    }
  }
  if (!is.list(x$answers) || !is.null(names(x$answers))) {
    p <- c(p, "answers must be an unnamed list of answer records")
  } else {
    for (i in seq_along(x$answers)) {
      a <- x$answers[[i]]
      if (!is.list(a) || !.ig_is_string(a$question_id) || !is.character(a$option_ids) ||
          anyNA(a$option_ids)) {
        p <- c(p, sprintf("answers[%d] needs question_id and option_ids (character)", i))
      }
    }
  }
  if (!is.list(x$rule_trace) || !length(x$rule_trace) || !is.null(names(x$rule_trace))) {
    p <- c(p, "rule_trace must be a non-empty list of steps (every recommendation is explainable)")
  } else {
    for (i in seq_along(x$rule_trace)) {
      s <- x$rule_trace[[i]]
      if (!is.list(s) || !.ig_is_string(s$question_id) || !.ig_is_string(s$transition_id) ||
          !.ig_is_string(s$to)) {
        p <- c(p, sprintf("rule_trace[%d] needs question_id, transition_id and to", i))
      }
    }
  }
  if (!is.list(x$rationale) || !.ig_is_string(x$rationale$en)) p <- c(p, "rationale$en must be a string")
  if (!is.null(x$estimand) && (!is.list(x$estimand) || !.ig_is_string(x$estimand$en))) {
    p <- c(p, "estimand must be NULL or a list with en")
  }
  pref <- x$preferred
  if (!is.null(pref) && (!is.list(pref) || !.ig_is_string(pref$analysis_variant_id))) {
    p <- c(p, "preferred must be NULL or a list with analysis_variant_id")
  }
  for (i in seq_along(x$assumptions)) {
    a <- x$assumptions[[i]]
    if (!is.list(a) || !.ig_is_string(a$assumption_id) || !.ig_is_string(a$state) ||
        !(a$state %in% v$assumption_state) || !.ig_is_string(a$text_en)) {
      p <- c(p, sprintf("assumptions[%d] needs assumption_id, state (%s) and text_en", i,
                        paste(v$assumption_state, collapse = "/")))
    }
  }
  for (i in seq_along(x$alternatives)) {
    a <- x$alternatives[[i]]
    if (!is.list(a) || !.ig_is_string(a$analysis_variant_id) || !.ig_is_string(a$why_not_preferred_en)) {
      p <- c(p, sprintf("alternatives[%d] needs analysis_variant_id and why_not_preferred_en", i))
    }
  }
  ri <- x$r_implementation
  if (!is.list(ri) || !str_or_null(ri$template_id) || !.ig_is_string(ri$template_status)) {
    p <- c(p, "r_implementation needs template_id (string or NULL) and template_status")
  }
  pl <- x$planning
  if (!.ig_is_string(pl$planning_compatibility) || !(pl$planning_compatibility %in% v$planning_compatibility)) {
    p <- c(p, sprintf("planning_compatibility must be one of: %s", paste(v$planning_compatibility, collapse = ", ")))
  }
  if (!.ig_is_string(pl$planning_availability) || !(pl$planning_availability %in% v$planning_availability)) {
    p <- c(p, sprintf("planning_availability must be one of: %s", paste(v$planning_availability, collapse = ", ")))
  }
  if (identical(pl$planning_compatibility, "exact_gpower_match") && !length(pl$evidence_ids)) {
    p <- c(p, "planning_compatibility 'exact_gpower_match' requires evidence_ids of authentic G*Power captures")
  }
  if (!is.null(pl$gpower_method_id)) {
    if (!.ig_is_string(pl$gpower_method_id) || !(pl$gpower_method_id %in% .ig_list_method_ids())) {
      p <- c(p, "planning$gpower_method_id must be a registry method id or NULL")
    } else if (identical(.ig_get_metadata(pl$gpower_method_id)$reference_status, "not_applicable")) {
      # Extensions such as Boschloo's test are registry methods, but they are
      # not G*Power procedures and cannot be named as one.
      p <- c(p, sprintf(
        "planning$gpower_method_id '%s' is an inferguidR extension, not a G*Power procedure",
        pl$gpower_method_id))
    }
  }
  if (!str_or_null(pl$note_en)) p <- c(p, "planning$note_en must be a string or NULL")
  if (!str_or_null(pl$note_de)) p <- c(p, "planning$note_de must be a string or NULL")
  if (length(pl$spec_refs) && !all(grepl("^[0-9a-f]{64}$", pl$spec_refs))) {
    p <- c(p, "planning$spec_refs must be SHA-256 specification hashes")
  }
  if (!str_or_null(x$visual_preset)) p <- c(p, "visual_preset must be a string or NULL")
  # Status consistency
  if (.ig_is_string(x$status)) {
    if (x$status %in% c("recommended", "conditional") && is.null(pref)) {
      p <- c(p, sprintf("a '%s' recommendation needs a preferred method", x$status))
    }
    if (x$status %in% c("needs_information", "outside_supported_scope") && !is.null(pref)) {
      p <- c(p, sprintf("a '%s' outcome must not name a preferred method", x$status))
    }
    if (x$status == "needs_information" && !length(x$missing_information)) {
      p <- c(p, "a 'needs_information' outcome must list the missing information")
    }
    if (x$status == "outside_supported_scope" &&
        !identical(pl$planning_compatibility, "not_applicable")) {
      p <- c(p, "an 'outside_supported_scope' outcome has planning_compatibility 'not_applicable'")
    }
  }
  if (.ig_contains_result(x)) {
    p <- c(p, "a recommendation references specifications, never numerical results")
  }
  p
}

.ig_contains_result <- function(x) {
  if (inherits(x, c("power_result", "classification_result", "observed_test_result"))) return(TRUE)
  if (is.list(x)) {
    if (.ig_is_string(x$object_type) &&
        x$object_type %in% c("power_result", "classification_result", "observed_test_result")) {
      return(TRUE)
    }
    return(any(vapply(x, .ig_contains_result, logical(1))))
  }
  FALSE
}

.ig_recommendation_body <- function(x) {
  lst <- .ig_recommendation_to_list(x)
  lst$recommendation_hash <- NULL
  lst
}

.ig_recommendation_to_list <- function(x) {
  x <- unclass(x)
  list(
    schema_version = x$schema_version, object_type = x$object_type,
    status = x$status, outcome_id = x$outcome_id, rules = x$rules,
    answers = lapply(x$answers, function(a) {
      a$option_ids <- I(as.character(a$option_ids))
      a
    }),
    rule_trace = x$rule_trace, estimand = x$estimand,
    design = .ig_named_or_empty(x$design), preferred = x$preferred,
    rationale = x$rationale, assumptions = x$assumptions,
    alternatives = x$alternatives,
    missing_information = I(x$missing_information),
    r_implementation = x$r_implementation,
    planning = list(planning_compatibility = x$planning$planning_compatibility,
                    planning_availability = x$planning$planning_availability,
                    gpower_method_id = x$planning$gpower_method_id,
                    spec_refs = I(x$planning$spec_refs),
                    evidence_ids = I(x$planning$evidence_ids),
                    note_en = x$planning$note_en,
                    note_de = x$planning$note_de),
    learning_links = I(x$learning_links), visual_preset = x$visual_preset,
    references = I(x$references), next_actions = I(x$next_actions),
    recommendation_hash = x$recommendation_hash
  )
}

.ig_recommendation_from_list <- function(x) {
  if (!identical(x$object_type, "recommendation") ||
      !identical(x$schema_version, .ig_schema_versions[["recommendation"]])) {
    .ig_abort("invalid_input", "not a supported recommendation document")
  }
  chr <- function(z) .ig_flatten(z, "character") %||% character()
  answers <- lapply(x$answers %||% list(), function(a) {
    a$option_ids <- chr(a$option_ids)
    a
  })
  pl <- x$planning %||% list()
  rec <- .ig_new_recommendation(
    status = x$status, outcome_id = x$outcome_id, rules = x$rules, answers = answers,
    rule_trace = x$rule_trace %||% list(), rationale = x$rationale, estimand = x$estimand,
    design = x$design %||% list(), preferred = x$preferred,
    assumptions = x$assumptions %||% list(), alternatives = x$alternatives %||% list(),
    missing_information = chr(x$missing_information),
    r_implementation = x$r_implementation,
    planning = list(planning_compatibility = pl$planning_compatibility,
                    planning_availability = pl$planning_availability,
                    gpower_method_id = pl$gpower_method_id,
                    spec_refs = chr(pl$spec_refs), evidence_ids = chr(pl$evidence_ids),
                    note_en = pl$note_en, note_de = pl$note_de),
    learning_links = chr(x$learning_links), visual_preset = x$visual_preset,
    references = chr(x$references), next_actions = chr(x$next_actions)
  )
  if (!identical(rec$recommendation_hash, x$recommendation_hash)) {
    .ig_abort("invalid_input", "recommendation_hash does not match the content (the document was modified)")
  }
  rec
}

#' @export
format.recommendation <- function(x, ...) {
  c(sprintf("<recommendation> %s (%s)", x$outcome_id, x$status),
    if (!is.null(x$preferred)) sprintf("  preferred: %s", x$preferred$analysis_variant_id),
    sprintf("  rationale: %s", x$rationale$en),
    sprintf("  planning:  %s / %s", x$planning$planning_compatibility,
            x$planning$planning_availability),
    sprintf("  rules:     %s %s", x$rules$rules_id, x$rules$rules_version))
}

#' @export
print.recommendation <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}
