# Advisor recommendation ------------------------------------------------------------
#
# Builds the core `recommendation` object (R/core-recommendation.R) from a
# session. Mapping of the advisor record onto the core shape:
#
# * answers       raw answer (option_ids), effective answer, dont_know policy;
# * rule_trace    one step per answered question: transition id, order, kind
#                 and target; a pending walk ends with an
#                 engine.incomplete_policy step to out.needs_information;
# * rules         rule provenance: rules id/version/schema/status/hash, the
#                 source-tree, variant and resource versions, content hashes of
#                 every input, the base status and the status rules and
#                 modifiers that were applied, and the source-tree nodes and
#                 crosswalk rows behind the outcome;
# * design        effective answers by question, stage, destination, pending
#                 question and ignored deep-link answers;
# * preferred     variant (name, family, settings, analysis inputs), override
#                 record, rejected routes with reasons (why-not), planning
#                 inputs, what the G*Power manual documents for the variant
#                 (provenance, not parity) and the learning-resource lookup
#                 with explicit gaps;
# * assumptions   method assumptions with the variant they were stated for
#                 (state "unknown": the advisor does not verify them against
#                 data), assumptions recorded by 'I do not know' answers and
#                 status modifiers, and with an override the override record
#                 plus the overriding variant's own assumptions (or an
#                 explicit record that the rules state none);
# * alternatives  why-not texts and override_allowed;
# * r_implementation template id with an explicit availability status (no
#                 code and no function names until a reviewed template exists);
# * planning      compatibility and availability exactly as the rules state
#                 (no advisor route currently claims exact_gpower_match).
#
# The recommendation references specifications only; it never contains a
# numerical result.

.ig_advisor_match_rank <- c(exact = 1L, related = 2L, broader = 3L)

# Installed template file for a template id ("" when not shipped).
.ig_advisor_template_file <- function(template_id) {
  if (!.ig_is_string(template_id) || !grepl("^tpl\\.[a-z0-9_]+\\.v[0-9]+$", template_id)) return("")
  system.file("advisor", "templates", paste0(template_id, ".R"), package = "inferguidR")
}

# R implementation record: the template id and whether it is shipped. The
# r_functions call sketches of analysis-variants.json are drafts for the
# templates and are deliberately not exported, not even as function names:
# some name engines the project has rejected (for example exact2x2::boschloo,
# decision D12), so only a reviewed template may say which functions to call.
.ig_advisor_template_record <- function(template_id, rules_status) {
  if (is.null(template_id)) {
    return(list(template_id = NULL, template_status = "not_applicable", rules_template_status = rules_status,
                note_en = "No method is recommended, so there is no R template.",
                note_de = "Es wird keine Methode empfohlen, daher gibt es keine R-Vorlage."))
  }
  available <- nzchar(.ig_advisor_template_file(template_id))
  list(template_id = template_id,
       template_status = if (available) "available" else "template_not_yet_available",
       rules_template_status = rules_status,
       note_en = if (available) "The R template is bundled with the package." else
         "The R template for this variant is not yet available; no substitute code is generated.",
       note_de = if (available) "Die R-Vorlage ist im Paket enthalten." else
         "Die R-Vorlage f\u00fcr diese Variante ist noch nicht verf\u00fcgbar; es wird kein Ersatzcode erzeugt.")
}

# What the G*Power manual documents for a variant (analysis-variants.json
# gpower_candidate). This is documentary provenance only: it never implies
# parity, which needs matching authentic application evidence and qualification.
.ig_advisor_manual_relation <- function(variant) {
  gc <- variant$gpower_candidate
  if (!is.list(gc) || !.ig_is_string(gc$relation)) return(NULL)
  list(relation = gc$relation, manual_section = gc$manual_section,
       manual_page = if (is.numeric(gc$manual_page) && length(gc$manual_page) == 1L) as.integer(gc$manual_page),
       basis_en = gc$basis,
       parity_note_en = "Documented in the G*Power manual only; no parity is claimed without authentic application captures.")
}

# Learning resources for an analysis variant, looked up through the variant's
# learning_method_ids in the resource registry's crosswalk. Each role (course
# lab, method tutorial, overview) gets a status: exact, related_only,
# broader_only or gap. Links are teaching evidence only.
.ig_advisor_learning <- function(rs, variant_id) {
  v <- rs$variants[[variant_id]]
  mids <- as.character(unlist(v$learning_method_ids))
  role_of <- function(type) {
    if (identical(type, "course_lab")) "course" else if (identical(type, "tutorial")) "tutorial" else "overview"
  }
  links <- list()
  for (r in rs$resources) {
    hits <- Filter(function(cw) .ig_is_string(cw$method_id) && cw$method_id %in% mids, r$crosswalk %||% list())
    if (!length(hits)) next
    m <- vapply(hits, function(cw) cw$match %||% "broader", "")
    rank <- .ig_advisor_match_rank[m]
    rank[is.na(rank)] <- 99L
    best <- m[which.min(rank)]
    issues <- as.character(unlist(r$known_issues))
    # `match` is the best match over the variant's method ids; `method_matches`
    # says which method id each match is for, so that a resource matched
    # exactly on an auxiliary id (for example an effect-size lesson) is not
    # read as an exact lesson on the variant itself.
    hit_ids <- vapply(hits, `[[`, "", "method_id")
    hit_order <- order(match(hit_ids, mids), rank, method = "radix")
    links[[length(links) + 1L]] <- list(
      resource_id = r$resource_id, role = role_of(r$resource_type), resource_type = r$resource_type,
      match = best, method_ids = as.list(unique(hit_ids)),
      method_matches = lapply(hit_order, function(k) list(method_id = hit_ids[[k]], match = m[[k]])),
      title = r$title, published_url = r$published_url, repo_url = r$repo_url,
      known_issues = as.list(issues), caution = length(issues) > 0L
    )
  }
  role_block <- function(role) {
    ls <- Filter(function(l) identical(l$role, role), links)
    if (length(ls)) {
      key <- vapply(ls, function(l) sprintf("%02d %s", .ig_advisor_match_rank[[l$match]] %||% 99L,
                                            l$resource_id), "")
      ls <- ls[order(key, method = "radix")]
    }
    matches <- vapply(ls, `[[`, "", "match")
    status <- if (!length(ls)) "gap" else if ("exact" %in% matches) "exact" else
      if ("related" %in% matches) "related_only" else "broader_only"
    label <- switch(role, course = c("CTTIR course lab", "CTTIR-Kurslabor"),
                    tutorial = c("CTTIR method tutorial", "CTTIR-Methodentutorial"),
                    overview = c("CTTIR overview resource", "CTTIR-\u00dcbersichtsressource"))
    note <- switch(status,
                   exact = c("", ""),
                   related_only = c(sprintf("No %s teaches exactly this variant; the linked material is related only.", label[1]),
                                    sprintf("Kein %s behandelt genau diese Variante; das verlinkte Material ist nur verwandt.", label[2])),
                   broader_only = c(sprintf("No %s teaches this variant; the linked material is broader.", label[1]),
                                    sprintf("Kein %s behandelt diese Variante; das verlinkte Material ist allgemeiner.", label[2])),
                   gap = c(sprintf("Gap: no %s for this variant in the resource registry; none is substituted.", label[1]),
                           sprintf("L\u00fccke: kein %s f\u00fcr diese Variante im Ressourcenregister; es wird nichts ersetzt.", label[2])))
    list(status = status, note_en = note[1], note_de = note[2], links = ls)
  }
  list(learning_method_ids = as.list(mids),
       resources_version = rs$versions$learning_resources_version,
       course = role_block("course"), tutorial = role_block("tutorial"), overview = role_block("overview"),
       evidence_note_en = "Learning links are teaching material, never numerical validation evidence.",
       evidence_note_de = "Lernlinks sind Lehrmaterial, nie Belege f\u00fcr eine numerische Validierung.")
}

# Assumptions of an overriding variant. Outcomes, not variants, carry
# assumption texts, so they are taken from the outcome(s) whose primary
# variant is the overriding one (ids already listed are not repeated). Where
# no outcome has that variant as its primary, the rules state no assumptions
# for it, and an explicit record says so instead of leaving the list silently
# incomplete.
.ig_advisor_override_assumptions <- function(rs, variant_id, have) {
  homes <- Filter(function(o) identical(o$analysis_variant_id, variant_id), rs$outcomes)
  out <- list()
  for (o in homes) {
    for (a in o$assumptions %||% list()) {
      if (a$assumption_id %in% have) next
      have <- c(have, a$assumption_id)
      out[[length(out) + 1L]] <- list(
        assumption_id = a$assumption_id, state = "unknown", text_en = a$text_en, text_de = a$text_de,
        source = "override_variant", question_id = NULL, stated_for = variant_id,
        from_outcome_id = o$outcome_id)
    }
  }
  if (length(homes)) return(out)
  v <- rs$variants[[variant_id]]
  list(list(
    assumption_id = paste0("override.", variant_id, ".assumptions_not_stated"), state = "unknown",
    text_en = sprintf(paste0("Rules %s state no method-specific assumptions for %s; ",
                             "they must be checked separately before it is used."),
                      rs$versions$rules_version, v$name_en),
    text_de = sprintf(paste0("Die Regeln %s nennen keine methodenspezifischen Annahmen f\u00fcr %s; ",
                             "sie m\u00fcssen vor der Verwendung gesondert gepr\u00fcft werden."),
                      rs$versions$rules_version, v$name_de),
    source = "override_variant", question_id = NULL, stated_for = variant_id))
}

# Rule trace records for the core object.
.ig_advisor_trace <- function(walk) {
  steps <- lapply(walk$steps, function(s) {
    list(step = s$step, question_id = s$question_id, raw_answer = s$raw_answer,
         effective_answer = s$effective_answer, dont_know_policy = s$dont_know_policy,
         transition_id = s$transition_id, transition_order = s$transition_order,
         kind = s$kind, to = s$to)
  })
  if (!is.null(walk$pending_question_id)) {
    steps[[length(steps) + 1L]] <- list(
      step = length(steps) + 1L, question_id = walk$pending_question_id, raw_answer = NULL,
      effective_answer = NULL, dont_know_policy = NULL, transition_id = "engine.incomplete_policy",
      transition_order = NULL, kind = "pending", to = "out.needs_information")
  }
  steps
}

#' Recommendation for an advisor session
#'
#' Builds the immutable, versioned recommendation for the current state of a
#' session: answered questions (raw and effective), the rule trace, outcome,
#' status (`recommended`, `conditional`, `needs_information` or
#' `outside_supported_scope`), estimand, rationale, assumptions (including
#' those recorded by 'I do not know' answers and status modifiers, and with
#' an override the overriding variant's own assumptions or an explicit record
#' that the rules state none), alternatives with why-not texts and rejected
#' routes, missing information,
#' the R template id with an explicit availability status, planning
#' compatibility and availability, learning resources with explicit gaps,
#' source-tree nodes of both trees and the rules, variants and resources
#' versions.
#'
#' Statuses follow the rule set: a `recommended` outcome is downgraded to
#' `conditional` by an 'I do not know' answer that routes with an assumption,
#' by an applicable status modifier or by an override; statuses are never
#' upgraded. An unfinished walkthrough gives `needs_information` with the
#' pending question. Planning compatibility is never `exact_gpower_match`
#' without matching authentic 'G*Power' evidence and a qualified mapping;
#' no advisor route currently claims this status.
#'
#' @inheritParams advisor_answer
#' @return A `recommendation` object (see [ig_to_json()] for export).
#' @examples
#' s <- advisor_replay(c(q.aim = "compare", q.stage = "planning",
#'   q.outcome_type = "binary", q.unit = "individuals",
#'   q.structure = "independent", q.independence_check = "no",
#'   q.n_factors = "one", q.n_groups = "two", q.covariates = "no",
#'   q.binary_margins = "group_sizes_fixed", q.direction = "two_sided",
#'   q.missingness = "none"))
#' rec <- advisor_recommendation(s)
#' rec$preferred$analysis_variant_id
#' rec$planning$planning_compatibility
#' @family advisor
#' @export
advisor_recommendation <- function(session, rules = NULL) {
  .ig_advisor_check_session(session)
  rs <- .ig_advisor_session_rules(session, rules)
  walk <- session$state
  walk$pending <- walk$pending_question_id
  outcome <- .ig_advisor_outcome_of(rs, walk)
  override <- session$override
  st <- .ig_advisor_status(rs, outcome, walk, override)
  status <- st$status
  has_method <- status %in% c("recommended", "conditional")
  primary <- outcome$analysis_variant_id
  vid <- if (has_method) override %||% primary else NULL
  overridden <- !is.null(override) && has_method
  v <- if (!is.null(vid)) rs$variants[[vid]] else NULL
  alt_of <- function(id) {
    for (a in outcome$alternatives %||% list()) if (identical(a$analysis_variant_id, id)) return(a)
    NULL
  }
  ov_alt <- if (overridden) alt_of(override) else NULL
  ver <- rs$versions

  # Answers and design
  answers <- lapply(walk$steps, function(s) {
    list(question_id = s$question_id, option_ids = s$raw_answer, effective_option_id = s$effective_answer,
         dont_know_policy = s$dont_know_policy, step = s$step)
  })
  design <- list(
    answers = .ig_named_or_empty(as.list(walk$effective)),
    stage = .ig_advisor_lookup(walk$effective, "q.stage"),
    destination = outcome$destination,
    pending_question_id = walk$pending_question_id,
    ignored_answers = as.list(session$ignored_answers)
  )

  # Rule provenance
  rules <- list(
    rules_id = ver$rules_id, rules_version = ver$rules_version,
    rules_schema_version = ver$rules_schema_version, rules_status = ver$rules_status,
    rules_sha256 = ver$rules_sha256,
    uzh_tree_version = ver$uzh_tree_version, descrtab2_tree_version = ver$descrtab2_tree_version,
    analysis_variants_version = ver$analysis_variants_version,
    learning_resources_version = ver$learning_resources_version,
    data_sha256 = ver$data_sha256,
    base_status = st$base_status,
    status_rules_applied = as.list(st$rules_applied),
    status_modifiers_applied = as.list(st$modifiers_applied),
    source_nodes = list(uzh = as.list(as.character(unlist(outcome$source_nodes$uzh))),
                        descrtab2 = as.list(as.character(unlist(outcome$source_nodes$descrtab2)))),
    source_nodes_note = outcome$source_nodes_note,
    crosswalk_ids = as.list(as.character(unlist(outcome$crosswalk_ids)))
  )

  # Preferred method and variant
  preferred <- NULL
  if (has_method) {
    pl <- outcome$planning
    preferred <- list(
      analysis_variant_id = vid, outcome_variant_id = primary,
      name_en = v$name_en, name_de = v$name_de, family = v$family,
      outcome_title_en = outcome$title_en, outcome_title_de = outcome$title_de,
      overridden = overridden,
      override = if (overridden) list(
        analysis_variant_id = override, replaces = primary,
        why_not_preferred_en = ov_alt$why_not_preferred_en,
        why_not_preferred_de = ov_alt$why_not_preferred_de,
        status_cap = "conditional") else NULL,
      variant_settings = if (overridden) .ig_empty_object() else .ig_named_or_empty(outcome$variant_settings),
      analysis_inputs = as.list(as.character(unlist(outcome$analysis_inputs))),
      rejected_routes = lapply(outcome$rejected_routes %||% list(), function(r) {
        list(analysis_variant_id = r$analysis_variant_id, reason_en = r$reason_en, reason_de = r$reason_de)
      }),
      planning_inputs = if (overridden) list() else as.list(as.character(unlist(pl$planning_inputs))),
      planning_note_de = NULL, # set below from the planning record
      planning_compatibility_target = if (overridden) NULL else pl$planning_compatibility_target,
      gpower_manual_relation = .ig_advisor_manual_relation(v),
      learning = .ig_advisor_learning(rs, vid)
    )
  }

  # Assumptions: those the outcome states for its variant, then those recorded
  # on the path, then (with an override) the override itself and the
  # assumptions of the overriding variant.
  assumptions <- lapply(outcome$assumptions %||% list(), function(a) {
    list(assumption_id = a$assumption_id, state = "unknown", text_en = a$text_en, text_de = a$text_de,
         source = "outcome", question_id = NULL, stated_for = primary)
  })
  assumptions <- c(assumptions, st$assumptions_made)
  if (overridden) {
    assumptions[[length(assumptions) + 1L]] <- list(
      assumption_id = paste0("override.", override), state = "unknown",
      text_en = ov_alt$why_not_preferred_en, text_de = ov_alt$why_not_preferred_de,
      source = "override", question_id = NULL)
    have <- vapply(assumptions, `[[`, "", "assumption_id")
    assumptions <- c(assumptions, .ig_advisor_override_assumptions(rs, override, have))
  }

  # Alternatives (with an override the replaced variant becomes an alternative)
  alt_record <- function(a) {
    av <- rs$variants[[a$analysis_variant_id]]
    list(analysis_variant_id = a$analysis_variant_id,
         why_not_preferred_en = a$why_not_preferred_en, why_not_preferred_de = a$why_not_preferred_de,
         override_allowed = isTRUE(a$override_allowed) &&
           is.null(.ig_advisor_override_refusal(rs, walk, a$analysis_variant_id)), name_en = av$name_en, name_de = av$name_de)
  }
  alts <- outcome$alternatives %||% list()
  if (overridden) {
    alts <- Filter(function(a) !identical(a$analysis_variant_id, override), alts)
    alts <- c(list(list(
      analysis_variant_id = primary,
      why_not_preferred_en = "Recommended by the rules for the recorded design; replaced by a user override.",
      why_not_preferred_de = paste0("Von den Regeln f\u00fcr das erfasste Design empfohlen; ",
                                    "durch eine Auswahl der Nutzerin oder des Nutzers ersetzt."),
      override_allowed = TRUE)), alts)
  }
  alternatives <- lapply(alts, alt_record)

  # Missing information (question ids)
  missing <- if (!is.null(walk$pending_question_id)) {
    walk$pending_question_id
  } else if (identical(status, "needs_information")) {
    walk$steps[[length(walk$steps)]]$question_id
  } else {
    vapply(Filter(function(s) identical(s$dont_know_policy, "route_with_assumption") &&
                    identical(s$raw_answer, "dont_know"), walk$steps), `[[`, "", "question_id")
  }

  # Estimand and rationale (an override states what changed)
  estimand <- if (overridden) {
    list(en = v$estimand_en, de = NULL)
  } else if (.ig_is_string(outcome$estimand_en)) {
    list(en = outcome$estimand_en, de = outcome$estimand_de)
  }
  rationale <- list(en = outcome$rationale_en, de = outcome$rationale_de)
  if (overridden) {
    pv <- rs$variants[[primary]]
    rationale <- list(
      en = sprintf("%s User override: %s replaces %s. The rules did not prefer it because: %s",
                   outcome$rationale_en, v$name_en, pv$name_en, ov_alt$why_not_preferred_en),
      de = sprintf("%s Auswahl durch die Nutzerin oder den Nutzer: %s ersetzt %s. Die Regeln bevorzugten sie nicht, weil: %s",
                   outcome$rationale_de, v$name_de, pv$name_de, ov_alt$why_not_preferred_de))
  }

  # Planning: preserve the rules unless an explicit selected-variant mapping applies.
  # The core planning record keeps both notes (note_en, note_de); method
  # outcomes also carry the German note in preferred$planning_note_de.
  pl <- outcome$planning
  planning <- list(planning_compatibility = pl$planning_compatibility,
                   planning_availability = pl$planning_availability,
                   gpower_method_id = NULL, spec_refs = character(), evidence_ids = character(),
                   note_en = pl$note_en, note_de = pl$note_de)
  if (overridden && !identical(pl$planning_compatibility, "not_applicable")) {
    planning$planning_compatibility <- "not_yet_validated"
    planning$planning_availability <- "unavailable"
    planning$note_en <- sprintf(paste0(
      "Planning for %s (a user override) is not defined by rules %s; no sample size is offered. ",
      "The planning stated for %s (%s / %s) does not apply to the overriding variant."),
      override, ver$rules_version, primary, pl$planning_compatibility, pl$planning_availability)
    planning$note_de <- sprintf(paste0(
      "Die Planung f\u00fcr %s (Auswahl durch die Nutzerin oder den Nutzer) ist in den Regeln %s nicht ",
      "festgelegt; es wird keine Fallzahl angeboten. Die f\u00fcr %s angegebene Planung (%s / %s) gilt ",
      "nicht f\u00fcr die gew\u00e4hlte Variante."),
      override, ver$rules_version, primary, pl$planning_compatibility, pl$planning_availability)
  } else if (overridden) {
    # Planning stays not_applicable, but the rules' note was written for the
    # replaced variant; say so instead of presenting it as the override's.
    planning$note_en <- sprintf("Rules %s define no planning for %s (a user override). Note stated for %s: %s",
                                ver$rules_version, override, primary, pl$note_en)
    planning$note_de <- sprintf(paste0("Die Regeln %s legen f\u00fcr %s (Auswahl durch die Nutzerin oder den ",
                                       "Nutzer) keine Planung fest. Hinweis zu %s: %s"),
                                ver$rules_version, override, primary, pl$note_de)
  }
  # Explicit mappings apply to the selected variant, including minlike overrides.
  # Unknown answers never become an available design merely through a rule default.
  if (has_method && identical(pl$planning_availability, "available") &&
      .ig_advisor_boschloo_outcome(outcome$outcome_id, primary) &&
      !is.null(.ig_advisor_boschloo_variant(vid))) {
    planning$planning_compatibility <- "independently_validated_extension"
    planning$planning_availability <- if (length(missing)) "unavailable" else "available"
    planning$note_en <- paste0(vid, ": ", .ig_advisor_boschloo_planning_note("en"))
    planning$note_de <- paste0(vid, ": ", .ig_advisor_boschloo_planning_note("de"))
    preferred$planning_inputs <- as.list(c("p1", "p2", "alpha", "power or n1/n2", "allocation_ratio and allocation policy for a priori planning"))
  }
  if (has_method) preferred$planning_note_de <- planning$note_de

  template <- .ig_advisor_template_record(
    if (has_method) paste0("tpl.", vid, ".v1") else NULL,
    if (has_method && !overridden) outcome$r_template_status else if (has_method) NULL else outcome$r_template_status)

  learning_links <- if (has_method) {
    unique(unlist(lapply(c("course", "tutorial", "overview"), function(role) {
      vapply(preferred$learning[[role]]$links, `[[`, "", "resource_id")
    })))
  } else character()

  .ig_new_recommendation(
    status = status, outcome_id = outcome$outcome_id, rules = rules, answers = answers,
    rule_trace = .ig_advisor_trace(walk), rationale = rationale, estimand = estimand,
    design = design, preferred = preferred, assumptions = assumptions,
    alternatives = alternatives, missing_information = as.character(missing),
    r_implementation = template, planning = planning,
    learning_links = as.character(learning_links), visual_preset = NULL,
    references = as.character(unlist(outcome$references)),
    next_actions = as.character(unlist(outcome$next_actions))
  )
}
