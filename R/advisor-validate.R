# Rule-set validation -------------------------------------------------------------
#
# advisor_validate_rules() makes the checks of the phase 01A draft walker
# permanent: structure, bilingual texts, 'I do not know' handling, forbidden
# routing inputs, honest planning claims, cross-file references, an
# exhaustive enumeration of every answer path, and design-consistency
# invariants evaluated on every path. Problems are returned as a table of
# (check, message) rows so that callers and tests can see which rule fired.

# Checks whose failure makes compiling or walking the graph meaningless.
.ig_advisor_blocking_checks <- c("schema", "unique_ids", "options", "condition", "transitions")

# Terms that must not appear in option ids/labels (acceptance scenario AS21):
# routing must never depend on a normality, variance or sphericity pre-test
# or on an observed p-value.
.ig_advisor_forbidden_terms <- c("shapiro", "kolmogorov", "lilliefors", "levene", "bartlett",
                                 "mauchly", "normality test", "pre-test", "p_value",
                                 "significant", "normalit\u00e4tstest", "vortest", "signifikant")

.ig_advisor_dont_know_policies <- c("needs_information", "route_without_assumption",
                                    "route_with_assumption")

.ig_advisor_problem_table <- function(check = character(), message = character()) {
  data.frame(check = as.character(check), message = as.character(message),
             stringsAsFactors = FALSE)
}

# Structural checks (blocking). Returns a problem table.
.ig_advisor_structure_problems <- function(rs) {
  P <- list(check = character(), message = character())
  add <- function(check, msg) {
    if (length(msg)) {
      P$check <<- c(P$check, rep(check, length(msg)))
      P$message <<- c(P$message, msg)
    }
  }
  r <- rs$docs$rules
  if (!identical(r$schema_version, .ig_advisor_rules_schema)) {
    add("schema", sprintf("rules schema_version must be '%s'", .ig_advisor_rules_schema))
  }
  if (!.ig_is_string(r$rules_id) || !.ig_is_string(r$rules_version)) {
    add("schema", "rules_id and rules_version must be strings")
  }
  if (!identical(rs$docs$variants$schema_version, .ig_advisor_variants_schema)) {
    add("schema", sprintf("analysis-variants schema_version must be '%s'", .ig_advisor_variants_schema))
  }
  for (f in c("questions", "transitions", "outcomes")) {
    if (!is.list(r[[f]]) || !length(r[[f]]) || !is.null(names(r[[f]]))) {
      add("schema", sprintf("rules.%s must be a non-empty array", f))
    } else if (!all(vapply(r[[f]], .ig_is_named, logical(1)))) {
      add("schema", sprintf("every entry of rules.%s must be a JSON object", f))
    }
  }
  if (!is.list(r$engine) || !.ig_is_named(r$engine)) add("schema", "rules.engine must be a JSON object")
  if (length(P$check)) return(.ig_advisor_problem_table(P$check, P$message))
  start <- r$engine$start_question
  if (!.ig_is_string(start) || is.null(rs$questions[[start]])) {
    add("schema", "engine.start_question must name an existing question")
  }
  if (is.null(rs$outcomes[["out.needs_information"]])) {
    add("schema", "the outcome out.needs_information is required (guards and pending walks use it)")
  }
  # Unique identifiers
  ids <- function(items, key) vapply(items, function(z) if (.ig_is_string(z[[key]])) z[[key]] else NA_character_, "")
  qids <- ids(r$questions, "question_id")
  tids <- ids(r$transitions, "transition_id")
  oids <- ids(r$outcomes, "outcome_id")
  for (spec in list(list("question", qids), list("transition", tids), list("outcome", oids))) {
    v <- spec[[2]]
    if (anyNA(v)) add("unique_ids", sprintf("every %s needs a string id", spec[[1]]))
    d <- unique(v[duplicated(v) & !is.na(v)])
    if (length(d)) add("unique_ids", sprintf("duplicate %s id '%s'", spec[[1]], d))
  }
  both <- intersect(qids, oids)
  if (length(both)) add("unique_ids", sprintf("id '%s' is both a question and an outcome", both))
  if (length(rs$variants) != length(rs$docs$variants$variants %||% list()) ||
      anyDuplicated(names(rs$variants))) {
    add("unique_ids", "analysis variant ids must be present and unique")
  }
  # Options
  for (q in rs$questions) {
    opts <- q$options
    if (!is.list(opts) || !length(opts) || !is.null(names(opts))) {
      add("options", sprintf("question '%s' has no options array", q$question_id))
      next
    }
    if (!all(vapply(opts, .ig_is_named, logical(1)))) {
      add("options", sprintf("question '%s': every option must be a JSON object", q$question_id))
      next
    }
    oid <- ids(opts, "option_id")
    if (anyNA(oid)) add("options", sprintf("question '%s' has an option without option_id", q$question_id))
    d <- unique(oid[duplicated(oid) & !is.na(oid)])
    if (length(d)) add("options", sprintf("question '%s' repeats option id '%s'", q$question_id, d))
  }
  # Conditions
  for (q in rs$questions) {
    add("condition", .ig_advisor_condition_problems(q$applies_when, rs$questions,
                                                    paste0(q$question_id, "/applies_when")))
  }
  for (t in r$transitions) {
    add("condition", .ig_advisor_condition_problems(t$when, rs$questions,
                                                    paste0(t$transition_id %||% "?", "/when")))
  }
  for (o in rs$outcomes) {
    for (m in o$status_modifiers %||% list()) {
      add("condition", .ig_advisor_condition_problems(
        m$when, rs$questions, sprintf("%s/%s/when", o$outcome_id, m$modifier_id %||% "?")))
    }
  }
  # Transitions
  for (t in r$transitions) {
    tid <- t$transition_id %||% "?"
    if (!.ig_is_string(t$from) || is.null(rs$questions[[t$from]])) {
      add("transitions", sprintf("%s: 'from' must name an existing question", tid))
    }
    if (!.ig_is_string(t$to) || (is.null(rs$questions[[t$to]]) && is.null(rs$outcomes[[t$to]]))) {
      add("transitions", sprintf("%s: unknown target '%s'", tid, format(t$to)))
    }
    if (!.ig_is_string(t$kind) || !(t$kind %in% c("route", "guard"))) {
      add("transitions", sprintf("%s: kind must be 'route' or 'guard'", tid))
    }
  }
  for (qid in names(rs$questions)) {
    ts <- rs$transitions[[qid]]
    if (!length(ts)) {
      add("transitions", sprintf("question '%s' has no transitions (dead end)", qid))
      next
    }
    ord <- vapply(ts, function(t) if (is.numeric(t$order) && length(t$order) == 1L) as.numeric(t$order) else NA_real_,
                  numeric(1))
    if (anyNA(ord) || !identical(as.numeric(sort(ord)), as.numeric(seq_along(ts)))) {
      add("transitions", sprintf("question '%s': transition orders must be 1..%d without gaps or repeats",
                                 qid, length(ts)))
    }
    last <- ts[[length(ts)]]
    if (!identical(last$kind, "guard") || !is.list(last$when) || !isTRUE(last$when$always) ||
        length(last$when) != 1L || !identical(last$to, "out.needs_information")) {
      add("transitions", sprintf(
        "question '%s': the last transition must be an unconditional guard to out.needs_information", qid))
    }
    early <- vapply(ts[-length(ts)], function(t) identical(t$kind, "guard"), logical(1))
    if (any(early)) add("transitions", sprintf("question '%s': a guard transition must come last", qid))
  }
  .ig_advisor_problem_table(P$check, P$message)
}

# Content checks (non-blocking). Returns a problem table.
.ig_advisor_content_problems <- function(rs) {
  P <- list(check = character(), message = character())
  add <- function(check, msg) {
    if (length(msg)) {
      P$check <<- c(P$check, rep(check, length(msg)))
      P$message <<- c(P$message, msg)
    }
  }
  r <- rs$docs$rules
  vocab <- .ig_recommendation_vocab
  # Bilingual texts: every *_en string in questions and outcomes needs a *_de string.
  for (q in rs$questions) add("bilingual", .ig_advisor_bilingual_problems(q, q$question_id))
  for (o in rs$outcomes) add("bilingual", .ig_advisor_bilingual_problems(o, o$outcome_id))
  for (q in rs$questions) {
    for (f in c("text_en", "text_de", "why_it_matters_en", "why_it_matters_de")) {
      if (!.ig_is_string(q[[f]])) add("bilingual", sprintf("%s: %s is missing", q$question_id, f))
    }
    for (o in q$options) {
      for (f in c("label_en", "label_de")) {
        if (!.ig_is_string(o[[f]])) add("bilingual", sprintf("%s:%s: %s is missing", q$question_id, o$option_id, f))
      }
    }
  }
  for (o in rs$outcomes) {
    for (f in c("title_en", "title_de", "rationale_en", "rationale_de")) {
      if (!.ig_is_string(o[[f]])) add("bilingual", sprintf("%s: %s is missing", o$outcome_id, f))
    }
  }
  # 'I do not know' option on every question
  for (q in rs$questions) {
    dk <- rs$options[[q$question_id]][["dont_know"]]
    if (is.null(dk)) {
      add("dont_know_option", sprintf("question '%s' has no 'dont_know' option", q$question_id))
      next
    }
    pol <- dk$dont_know_policy
    if (!.ig_is_string(pol) || !(pol %in% .ig_advisor_dont_know_policies)) {
      add("dont_know_option", sprintf("question '%s': dont_know_policy must be one of %s", q$question_id,
                                      paste(.ig_advisor_dont_know_policies, collapse = ", ")))
    }
    if (!.ig_is_string(dk$explanation_en) || !.ig_is_string(dk$explanation_de)) {
      add("dont_know_option", sprintf("question '%s': the dont_know option needs explanation_en and explanation_de",
                                      q$question_id))
    }
    if (.ig_is_string(pol) && pol %in% c("route_without_assumption", "route_with_assumption")) {
      ra <- dk$routes_as
      if (!.ig_is_string(ra) || identical(ra, "dont_know") || is.null(rs$options[[q$question_id]][[ra]])) {
        add("dont_know_option", sprintf("question '%s': routes_as must name another option", q$question_id))
      }
      if (pol == "route_with_assumption" &&
          (!.ig_is_string(dk$assumption_en) || !.ig_is_string(dk$assumption_de))) {
        add("dont_know_option", sprintf("question '%s': route_with_assumption needs assumption_en and assumption_de",
                                        q$question_id))
      }
    }
  }
  # Routing inputs: design/aim only; no pre-test or significance wording
  for (q in rs$questions) {
    if (!identical(q$input_kind, "design_or_aim")) {
      add("input_kind", sprintf("question '%s': input_kind must be 'design_or_aim'", q$question_id))
    }
  }
  if (!length(r$engine$forbidden_routing_inputs)) {
    add("input_kind", "engine.forbidden_routing_inputs must be declared")
  }
  pattern <- paste(.ig_advisor_forbidden_terms, collapse = "|")
  for (q in rs$questions) {
    for (o in q$options) {
      txt <- paste(o$option_id, o$label_en %||% "", o$label_de %||% "")
      hit <- regmatches(tolower(txt), regexpr(pattern, tolower(txt)))
      if (length(hit)) {
        add("forbidden_terms", sprintf("%s:%s: option text contains '%s' (pre-test or significance routing)",
                                       q$question_id, o$option_id, hit))
      }
    }
  }
  # Outcomes
  status_vocab <- vocab$status
  for (o in rs$outcomes) {
    id <- o$outcome_id
    st <- o$status
    if (!.ig_is_string(st) || !(st %in% status_vocab)) {
      add("outcome_status", sprintf("%s: status must be one of %s", id, paste(status_vocab, collapse = ", ")))
      next
    }
    vid <- o$analysis_variant_id
    has_method <- st %in% c("recommended", "conditional")
    if (has_method && !.ig_is_string(vid)) {
      add("outcome_status", sprintf("%s: a '%s' outcome needs an analysis_variant_id", id, st))
    }
    if (!has_method && !is.null(vid)) {
      add("outcome_status", sprintf("%s: a '%s' outcome must not name an analysis variant", id, st))
    }
    if (!has_method && !length(o$next_actions)) {
      add("next_actions", sprintf("%s: a '%s' outcome needs at least one next action (no dead end)", id, st))
    }
    pl <- o$planning
    if (!is.list(pl) || !.ig_is_string(pl$planning_compatibility) ||
        !(pl$planning_compatibility %in% vocab$planning_compatibility)) {
      add("planning", sprintf("%s: planning_compatibility outside the vocabulary", id))
    } else if (identical(pl$planning_compatibility, "exact_gpower_match") ||
               (identical(pl$planning_compatibility, "independently_validated_extension") &&
                !.ig_advisor_boschloo_outcome(id, vid))) {
      add("planning", sprintf("%s: claims '%s' without a validated mapping",
                              id, pl$planning_compatibility))
    }
    if (is.list(pl) && identical(pl$planning_availability, "available") &&
        (!.ig_advisor_boschloo_outcome(id, vid) ||
         !identical(pl$planning_compatibility, "independently_validated_extension")))
      add("planning", sprintf("%s: no implemented same-test planning mapping", id))
    if (is.list(pl) && (!.ig_is_string(pl$planning_availability) ||
                        !(pl$planning_availability %in% vocab$planning_availability))) {
      add("planning", sprintf("%s: planning_availability outside the vocabulary", id))
    }
    if (identical(st, "outside_supported_scope") && is.list(pl) &&
        !identical(pl$planning_compatibility, "not_applicable")) {
      add("planning", sprintf("%s: an outside_supported_scope outcome has planning_compatibility not_applicable", id))
    }
    if (identical(o$destination, "classification_workspace") && is.list(pl) &&
        (!identical(pl$planning_compatibility, "not_applicable") ||
         !identical(pl$planning_availability, "not_applicable"))) {
      add("planning", sprintf("%s: classifier outcomes never carry test-power planning", id))
    }
    # Template id: exactly the template of the primary variant (code export
    # cannot select another test).
    if (.ig_is_string(vid)) {
      if (!identical(o$r_template_id, paste0("tpl.", vid, ".v1"))) {
        add("template", sprintf("%s: r_template_id must be 'tpl.%s.v1'", id, vid))
      }
      if (!.ig_is_string(o$r_template_status) || identical(o$r_template_status, "not_applicable")) {
        add("template", sprintf("%s: a method outcome needs a template status other than not_applicable", id))
      }
    } else if (!is.null(o$r_template_id) || !identical(o$r_template_status, "not_applicable")) {
      add("template", sprintf("%s: an outcome without a method has no template (status not_applicable)", id))
    }
    # Variant references
    refs <- c(if (.ig_is_string(vid)) vid,
              vapply(o$alternatives %||% list(), function(a) a$analysis_variant_id %||% NA_character_, ""),
              vapply(o$rejected_routes %||% list(), function(a) a$analysis_variant_id %||% NA_character_, ""))
    bad <- refs[is.na(refs) | !(refs %in% names(rs$variants))]
    if (length(bad)) add("variant_refs", sprintf("%s: unknown analysis variant '%s'", id, bad))
    for (a in o$alternatives %||% list()) {
      if (!.ig_is_flag(a$override_allowed)) {
        add("variant_refs", sprintf("%s: alternative '%s' needs a logical override_allowed", id,
                                    format(a$analysis_variant_id)))
      }
    }
    # Literature and source references
    bad <- setdiff(unlist(o$references), names(rs$references))
    if (length(bad)) add("references", sprintf("%s: unknown reference id '%s'", id, bad))
    sn <- o$source_nodes
    bad_u <- setdiff(unlist(sn$uzh), rs$uzh_nodes)
    bad_d <- setdiff(unlist(sn$descrtab2), rs$descrtab2_nodes)
    if (length(bad_u)) add("source_nodes", sprintf("%s: unknown UZH node '%s'", id, bad_u))
    if (length(bad_d)) add("source_nodes", sprintf("%s: unknown DescrTab2 node '%s'", id, bad_d))
    if (!is.list(sn) || is.null(names(sn)) || !all(c("uzh", "descrtab2") %in% names(sn))) {
      add("source_nodes", sprintf("%s: source_nodes must list 'uzh' and 'descrtab2' node ids", id))
    }
    bad <- setdiff(unlist(o$crosswalk_ids), rs$docs$crosswalk$crosswalk_id)
    if (length(bad)) add("crosswalk_refs", sprintf("%s: unknown crosswalk id '%s'", id, bad))
    for (m in o$status_modifiers %||% list()) {
      if (!identical(m$effect, "conditional")) {
        add("modifiers", sprintf("%s/%s: the only modifier effect is 'conditional'", id, m$modifier_id %||% "?"))
      }
      # A modifier that sends the user to a mixed model "listed as an
      # alternative" is only valid where that alternative is listed.
      if (.ig_is_string(m$assumption_en) &&
          grepl("mixed model listed as an alternative", m$assumption_en, fixed = TRUE)) {
        alts <- vapply(o$alternatives %||% list(), function(a) a$analysis_variant_id %||% "", "")
        if (!("lmm_repeated_measures" %in% alts)) {
          add("modifiers", sprintf("%s/%s: points to a mixed-model alternative that the outcome does not list",
                                   id, m$modifier_id))
        }
      }
    }
  }
  # Required assumptions and refused routes (contract invariants)
  for (inv in .ig_advisor_static_invariants()) {
    o <- rs$outcomes[[inv$outcome_id]]
    if (is.null(o)) next # a removed outcome is reported by the reachability/summary checks
    if (!is.null(inv$assumption_ids)) {
      have <- vapply(o$assumptions %||% list(), function(a) a$assumption_id %||% "", "")
      miss <- setdiff(inv$assumption_ids, have)
      if (length(miss)) add("required_assumptions", sprintf("%s: assumption '%s' is required (%s)",
                                                            inv$outcome_id, miss, inv$reason))
    }
    if (!is.null(inv$rejected_ids)) {
      have <- vapply(o$rejected_routes %||% list(), function(a) a$analysis_variant_id %||% "", "")
      miss <- setdiff(inv$rejected_ids, have)
      if (length(miss)) add("rejected_routes", sprintf("%s: variant '%s' must be a rejected route (%s)",
                                                       inv$outcome_id, miss, inv$reason))
    }
  }
  # Analysis variants: usage and learning ids
  usage <- list()
  for (o in rs$outcomes) {
    if (.ig_is_string(o$analysis_variant_id)) usage[[length(usage) + 1L]] <- c(o$analysis_variant_id, o$outcome_id, "primary")
    for (a in o$alternatives %||% list()) usage[[length(usage) + 1L]] <- c(a$analysis_variant_id %||% "", o$outcome_id, "alternative")
    for (a in o$rejected_routes %||% list()) usage[[length(usage) + 1L]] <- c(a$analysis_variant_id %||% "", o$outcome_id, "rejected_route")
  }
  usage_key <- vapply(usage, paste, "", collapse = "|")
  learning_ids <- unique(unlist(lapply(rs$resources, function(x) {
    vapply(x$crosswalk %||% list(), function(c) c$method_id %||% "", "")
  })))
  for (v in rs$variants) {
    vid <- v$analysis_variant_id
    ub <- vapply(v$used_by %||% list(), function(u) paste(vid, u$outcome_id %||% "", u$role %||% "", sep = "|"), "")
    mine <- usage_key[startsWith(usage_key, paste0(vid, "|"))]
    if (!setequal(ub, mine)) {
      add("variant_usage", sprintf("variant '%s': used_by does not match its use in rules.json", vid))
    }
    bad <- setdiff(unlist(v$learning_method_ids), learning_ids)
    if (length(bad)) add("learning_refs", sprintf("variant '%s': learning_method_id '%s' is not in the resource registry",
                                                  vid, bad))
  }
  # Crosswalk: reasons and coverage of every source node used
  cw <- rs$docs$crosswalk
  need <- c("crosswalk_id", "uzh_node_id", "descrtab2_node_id", "difference_reason")
  if (!all(need %in% names(cw))) {
    add("crosswalk_refs", sprintf("source crosswalk needs columns %s", paste(need, collapse = ", ")))
  } else {
    empty <- cw$crosswalk_id[!nzchar(trimws(cw$difference_reason))]
    if (length(empty)) add("crosswalk_reason", sprintf("crosswalk row '%s' has no difference_reason", empty))
    used_u <- unique(unlist(lapply(rs$outcomes, function(o) o$source_nodes$uzh)))
    used_d <- unique(unlist(lapply(rs$outcomes, function(o) o$source_nodes$descrtab2)))
    miss <- c(setdiff(used_u, cw$uzh_node_id), setdiff(used_d, cw$descrtab2_node_id))
    if (length(miss)) add("crosswalk_refs", sprintf("source node '%s' has no crosswalk row", miss))
  }
  # Documented counts
  s <- r$summary
  status_counts <- table(vapply(rs$outcomes, function(o) o$status %||% NA_character_, ""))
  expect <- list(
    n_questions = length(rs$questions),
    n_transitions = length(r$transitions),
    n_guard_transitions = sum(vapply(r$transitions, function(t) identical(t$kind, "guard"), logical(1))),
    n_outcomes = length(rs$outcomes),
    n_analysis_variants = length(rs$variants),
    n_crosswalk_rows = nrow(cw),
    n_source_rule_dispositions = length(r$source_rule_dispositions)
  )
  for (nm in names(expect)) {
    if (!identical(as.integer(s[[nm]]), as.integer(expect[[nm]]))) {
      add("summary_counts", sprintf("summary.%s is %s but the rule set has %d", nm,
                                    format(s[[nm]] %||% "missing"), expect[[nm]]))
    }
  }
  documented <- unlist(s$outcomes_by_status)
  if (is.null(documented) ||
      !setequal(names(documented), names(status_counts)) ||
      !all(as.integer(documented[names(status_counts)]) == as.integer(status_counts))) {
    add("summary_counts", "summary.outcomes_by_status does not match the outcomes")
  }
  .ig_advisor_problem_table(P$check, P$message)
}

# Every *_en string of a record (recursively) needs a non-empty *_de sibling.
.ig_advisor_bilingual_problems <- function(x, where) {
  out <- character()
  if (!is.list(x)) return(out)
  nms <- names(x)
  if (!is.null(nms)) {
    for (nm in nms[grepl("_en$", nms)]) {
      if (.ig_is_string(x[[nm]])) {
        de <- sub("_en$", "_de", nm)
        if (!.ig_is_string(x[[de]])) out <- c(out, sprintf("%s: %s has no %s text", where, nm, de))
      }
    }
  }
  for (i in seq_along(x)) {
    sub_where <- if (!is.null(nms) && nzchar(nms[i])) paste0(where, "/", nms[i]) else paste0(where, "[", i, "]")
    out <- c(out, .ig_advisor_bilingual_problems(x[[i]], sub_where))
  }
  out
}

# Static contract invariants (BOSCHLOO, DIGITAL_ADVISOR, rules-review R03).
.ig_advisor_static_invariants <- function() {
  ties <- "Boschloo ordering ties must be handled explicitly (BOSCHLOO.md, a.ordering_ties)"
  rank <- "the rank-test null hypothesis is identical distributions (rules-review R03)"
  paired <- "paired binary data never use an independent-groups test (DIGITAL_ADVISOR acceptance scenario)"
  list(
    list(outcome_id = "out.two_group.boschloo_central", assumption_ids = "a.ordering_ties", reason = ties),
    list(outcome_id = "out.two_group.boschloo_one_sided", assumption_ids = "a.ordering_ties", reason = ties),
    list(outcome_id = "out.two_group.wmw", assumption_ids = "a.rank_null_identical", reason = rank),
    list(outcome_id = "out.kgroup.kruskal", assumption_ids = "a.rank_null_identical", reason = rank),
    list(outcome_id = "out.paired.mcnemar",
         rejected_ids = c("boschloo_two_sided_central", "boschloo_two_sided_minlike", "boschloo_one_sided",
                          "fisher_exact_2x2", "prop_two_independent_pearson_chisq"),
         reason = paired)
  )
}

# Design-consistency invariants evaluated on every enumerated answer path.
# `outcomes`/`exclude` select outcomes by regular expression and `require`
# must hold for their effective answers; `when` + `allowed` restrict which
# outcomes a design may reach. Conditions use the rule grammar and the same
# allowlisted compiler.
.ig_advisor_path_invariants <- function() {
  q <- function(id, ...) list(q = id, "in" = list(...))
  all_of <- function(...) list(all = list(...))
  any_of <- function(...) list(any = list(...))
  ot <- function(...) q("q.outcome_type", ...)
  indep <- all_of(q("q.structure", "independent"), q("q.independence_check", "no"))
  one_record <- all_of(q("q.unit", "individuals"), q("q.records_per_individual", "one"))
  list(
    list(id = "independent_groups", outcomes = "^out\\.(two_group|kgroup|factorial|spread)\\.|^out\\.ancova$",
         require = indep,
         text = "independent-group methods need structure = independent and a confirmed independence check"),
    list(id = "contingency_tables", outcomes = "^out\\.rxc\\.",
         require = any_of(indep, all_of(q("q.aim", "associate"), one_record)),
         text = "contingency-table tests need independent groups or one record per individual"),
    list(id = "paired", outcomes = "^out\\.paired\\.",
         require = all_of(q("q.structure", "paired"), q("q.pairing_id", "yes")),
         text = "paired methods need structure = paired and a confirmed pairing id"),
    list(id = "repeated", outcomes = "^out\\.repeated\\.",
         require = all_of(q("q.structure", "repeated"), q("q.pairing_id", "yes")),
         text = "repeated-measures methods need structure = repeated and a confirmed id"),
    list(id = "one_sample", outcomes = "^out\\.one_sample\\.",
         require = all_of(q("q.structure", "one_sample"), q("q.reference_value", "yes")),
         text = "one-sample tests need a justified reference value"),
    list(id = "clustered", outcomes = "^out\\.clustered\\.",
         require = any_of(q("q.unit", "clusters"), q("q.records_per_individual", "several")),
         text = "cluster models only for clusters or several records per individual"),
    list(id = "one_record_per_individual", outcomes = "^out\\.(assoc|predict)\\.",
         require = one_record,
         text = "correlation and prediction outcomes need one record per individual"),
    list(id = "regression_independence", outcomes = "^out\\.model\\.",
         require = any_of(indep, one_record),
         text = "regression outcomes need independent units"),
    list(id = "groups_over_time", outcomes = "^out\\.", exclude = "^out\\.group_time_design$",
         require = list(not = q("q.structure", "groups_over_time")),
         text = "groups measured over time never reach an analysis outcome"),
    list(id = "pooled_variance_justified", outcomes = "^out\\.(two_group\\.student|kgroup\\.anova_f)$",
         require = q("q.variance", "justified"),
         text = "pooled-variance tests only when equal variances are justified in advance"),
    list(id = "type_continuous",
         outcomes = paste0("^out\\.(two_group\\.(welch|student)|kgroup\\.(welch_anova|anova_f)|factorial\\.anova|",
                           "ancova|spread\\..+|paired\\.(t|wilcoxon)|repeated\\.(rm_anova|lmm)|one_sample\\.t|",
                           "clustered\\.lmm|model\\.linear_.+|predict\\.linear)$"),
         require = ot("continuous"), text = "method needs a continuous outcome"),
    list(id = "type_continuous_or_ordinal",
         outcomes = "^out\\.(two_group\\.wmw|kgroup\\.kruskal|repeated\\.friedman)$",
         require = ot("continuous", "ordinal"), text = "rank method needs a continuous or ordinal outcome"),
    list(id = "type_binary",
         outcomes = paste0("^out\\.(two_group\\.(boschloo_central|boschloo_one_sided|fisher_2x2)|paired\\.mcnemar|",
                           "repeated\\.cochran_q|one_sample\\.binomial|clustered\\.binary|model\\.logistic|",
                           "predict\\.logistic)$"),
         require = ot("binary"), text = "method needs a binary outcome"),
    list(id = "type_ordinal", outcomes = "^out\\.(paired\\.sign|one_sample\\.sign|model\\.ordinal)$",
         require = ot("ordinal"), text = "method needs an ordinal outcome"),
    list(id = "type_nominal", outcomes = "^out\\.(paired\\.nominal|one_sample\\.gof|model\\.multinomial)$",
         require = ot("nominal"), text = "method needs a nominal outcome"),
    list(id = "type_count", outcomes = "^out\\.model\\.count$", require = ot("count"),
         text = "method needs a count outcome"),
    list(id = "type_categorical_table", outcomes = "^out\\.rxc\\.",
         require = any_of(ot("binary", "nominal"), q("q.assoc_types", "both_categorical")),
         text = "contingency-table tests need categorical variables"),
    list(id = "boschloo_not_replaced",
         when = all_of(indep, ot("binary"), q("q.n_groups", "two"), q("q.covariates", "no"),
                       q("q.binary_margins", "group_sizes_fixed")),
         allowed = "^out\\.two_group\\.boschloo_(central|one_sided)$",
         text = "two independent binomial samples are analysed with Boschloo's test"),
    list(id = "welch_not_replaced",
         when = all_of(indep, ot("continuous"), q("q.n_groups", "two"), q("q.covariates", "no"),
                       q("q.estimand_continuous", "mean"), q("q.variance", "not_justified")),
         allowed = "^out\\.two_group\\.welch$",
         text = "without justified equal variances two means are compared with Welch's test")
  )
}

.ig_advisor_compile_invariants <- function(invs) {
  lapply(invs, function(inv) {
    inv$require_fn <- if (!is.null(inv$require)) .ig_advisor_compile_condition(inv$require)
    inv$when_fn <- if (!is.null(inv$when)) .ig_advisor_compile_condition(inv$when)
    inv
  })
}

# Exhaustive depth-first enumeration of every answer path (all options of
# every question, including 'dont_know'). Shares no code with the session
# walk beyond the compiled conditions and the effective-answer rule.
.ig_advisor_enumerate <- function(rs, invariants = .ig_advisor_path_invariants(), max_examples = 5L) {
  invs <- .ig_advisor_compile_invariants(invariants)
  acc <- new.env(parent = emptyenv())
  acc$n_paths <- 0L
  acc$max_len <- 0L
  acc$outcome_hits <- stats::setNames(integer(length(rs$outcomes)), names(rs$outcomes))
  acc$question_hits <- stats::setNames(integer(length(rs$questions)), names(rs$questions))
  all_tids <- unlist(lapply(rs$transitions, function(ts) vapply(ts, `[[`, "", "transition_id")))
  acc$transition_hits <- stats::setNames(integer(length(all_tids)), all_tids)
  acc$combos <- character()
  acc$issues <- list(guard_reached = character(), applies_when = character(), cycle = character(),
                     design_consistency = character())
  acc$issue_counts <- c(guard_reached = 0L, applies_when = 0L, cycle = 0L, design_consistency = 0L)
  record <- function(type, msg) {
    acc$issue_counts[[type]] <- acc$issue_counts[[type]] + 1L
    if (length(acc$issues[[type]]) < max_examples) acc$issues[[type]] <- c(acc$issues[[type]], msg)
  }
  acc$examples <- list()
  # Which outcome-selecting invariants apply to an outcome, and whether an
  # outcome is allowed by a `when` invariant (cached per outcome id).
  applies_to <- lapply(rs$outcomes, function(o) {
    oid <- o$outcome_id
    Filter(function(inv) !is.null(inv$outcomes) && grepl(inv$outcomes, oid) &&
             (is.null(inv$exclude) || !grepl(inv$exclude, oid)), invs)
  })
  when_invs <- Filter(function(inv) is.null(inv$outcomes) && !is.null(inv$when_fn), invs)
  allowed_by <- lapply(rs$outcomes, function(o) {
    vapply(when_invs, function(inv) grepl(inv$allowed, o$outcome_id), logical(1))
  })
  path_txt <- function(raw) paste(sprintf("%s=%s", names(raw), raw), collapse = " ")
  finish <- function(oid, eff, raw, dk_assume, depth) {
    acc$n_paths <- acc$n_paths + 1L
    acc$max_len <- max(acc$max_len, depth)
    acc$outcome_hits[[oid]] <- acc$outcome_hits[[oid]] + 1L
    o <- rs$outcomes[[oid]]
    status <- o$status
    downgradable <- identical(status, "recommended")
    if (downgradable && dk_assume) status <- "conditional"
    mods <- rs$compiled$modifiers[[oid]]
    for (j in seq_along(mods)) {
      if (mods[[j]](eff) && identical(status, "recommended")) status <- "conditional"
    }
    combo <- paste(oid, status, sep = "|")
    if (is.null(acc$examples[[combo]])) {
      acc$combos <- c(acc$combos, combo)
      acc$examples[[combo]] <- raw
    }
    for (inv in applies_to[[oid]]) {
      if (!inv$require_fn(eff)) {
        record("design_consistency", sprintf("[%s] %s reached by: %s (%s)", inv$id, oid, path_txt(raw), inv$text))
      }
    }
    ok <- allowed_by[[oid]]
    for (j in seq_along(when_invs)) {
      if (!ok[j] && when_invs[[j]]$when_fn(eff)) {
        record("design_consistency", sprintf("[%s] %s reached by: %s (%s)", when_invs[[j]]$id, oid,
                                             path_txt(raw), when_invs[[j]]$text))
      }
    }
  }
  dfs <- function(node, eff, raw, dk_assume, visited, depth) {
    acc$question_hits[[node]] <- acc$question_hits[[node]] + 1L
    if (!rs$compiled$applies[[node]](eff)) {
      record("applies_when", sprintf("%s reached while its applies_when is false, by: %s", node, path_txt(raw)))
    }
    ts <- rs$transitions[[node]]
    fs <- rs$compiled$transitions[[node]]
    for (opt in rs$options[[node]]) {
      e <- .ig_advisor_effective(opt)
      eff2 <- eff
      eff2[[node]] <- e
      raw2 <- raw
      raw2[[node]] <- opt$option_id
      dk2 <- dk_assume || (identical(opt$option_id, "dont_know") &&
                             identical(opt$dont_know_policy, "route_with_assumption"))
      k <- NA_integer_
      for (i in seq_along(fs)) if (fs[[i]](eff2)) { k <- i; break }
      tr <- ts[[k]]
      acc$transition_hits[[tr$transition_id]] <- acc$transition_hits[[tr$transition_id]] + 1L
      if (identical(tr$kind, "guard")) {
        record("guard_reached", sprintf("guard %s reached by: %s", tr$transition_id, path_txt(raw2)))
      }
      if (!is.null(rs$outcomes[[tr$to]])) {
        finish(tr$to, eff2, raw2, dk2, depth)
      } else if (tr$to %in% c(visited, node)) {
        record("cycle", sprintf("cycle at %s by: %s", tr$to, path_txt(raw2)))
      } else {
        dfs(tr$to, eff2, raw2, dk2, c(visited, node), depth + 1L)
      }
    }
  }
  empty <- stats::setNames(character(), character())
  dfs(rs$docs$rules$engine$start_question, empty, empty, FALSE, character(), 1L)
  guards <- unlist(lapply(rs$transitions, function(ts) {
    vapply(Filter(function(t) identical(t$kind, "guard"), ts), `[[`, "", "transition_id")
  }), use.names = FALSE)
  list(
    answer_paths = acc$n_paths,
    max_questions_on_a_path = acc$max_len,
    outcome_status_combinations = length(acc$combos),
    combos = sort(acc$combos),
    # One raw answer path (question id -> option id) per outcome-status combination.
    examples = acc$examples[sort(acc$combos)],
    outcome_hits = acc$outcome_hits,
    question_hits = acc$question_hits,
    transition_hits = acc$transition_hits,
    guard_transitions = guards,
    unreachable_outcomes = names(acc$outcome_hits)[acc$outcome_hits == 0L],
    unreachable_questions = names(acc$question_hits)[acc$question_hits == 0L],
    unused_transitions = setdiff(names(acc$transition_hits)[acc$transition_hits == 0L], guards),
    issue_counts = acc$issue_counts,
    issues = acc$issues
  )
}

#' Validate an advisor rule set
#'
#' Checks a rule graph and its companion files before they are used:
#'
#' * structure: unique ids, allowlisted conditions that reference existing
#'   questions and options, contiguous transition orders ending in an
#'   unconditional guard, existing targets;
#' * content: English and German texts for every question, option and
#'   outcome; an 'I do not know' option with explanation and policy on every
#'   question; routing only on design and aim (no normality, variance or
#'   sphericity pre-test and no p-value wording); no planning claim of
#'   `exact_gpower_match` or `independently_validated_extension`; template ids
#'   that match the primary variant; resolvable variant, reference,
#'   source-node, crosswalk and learning ids; documented summary counts;
#' * walk: an exhaustive enumeration of every answer path, which must never
#'   reach a guard, must respect `applies_when`, must not cycle, must reach
#'   every question, outcome and non-guard transition, and must match the
#'   documented path count; design-consistency invariants are evaluated on
#'   every path (for example, paired data never reach an independent-groups
#'   test, and two independent binomial samples are always analysed with
#'   Boschloo's test).
#'
#' @inheritParams advisor_rules
#' @param walk If `TRUE` (default), enumerate all answer paths. The walk is
#'   skipped when structural problems make it meaningless.
#' @return An object of class `advisor_rules_validation`: a list with `valid`
#'   (logical), `problems` (data frame with columns `check` and `message`),
#'   `walk` (enumeration summary, or `NULL`) and `versions`. Files that cannot
#'   be read or parsed (not JSON, duplicate keys, remote paths, unexpected
#'   shapes) give `valid = FALSE` with a `schema` problem; the function does
#'   not raise errors for defective rule sets.
#' @examples
#' \donttest{
#' v <- advisor_validate_rules()
#' v$valid
#' v$walk$answer_paths
#' }
#' @seealso [advisor_rules()]
#' @family advisor
#' @export
advisor_validate_rules <- function(rules = NULL, variants = NULL, uzh_tree = NULL,
                                   descrtab2_tree = NULL, crosswalk = NULL, resources = NULL,
                                   walk = TRUE) {
  paths <- list(rules = rules, variants = variants, uzh_tree = uzh_tree,
                descrtab2_tree = descrtab2_tree, crosswalk = crosswalk, resources = resources)
  rs <- tryCatch(.ig_advisor_load_checked(paths), inferguidR_error = function(e) e)
  if (inherits(rs, "inferguidR_error")) {
    return(structure(list(valid = FALSE,
                          problems = .ig_advisor_problem_table("schema", conditionMessage(rs)),
                          walk = NULL, versions = NULL),
                     class = "advisor_rules_validation"))
  }
  problems <- .ig_advisor_guarded(.ig_advisor_structure_problems(rs), "schema",
                                  "the structure could not be checked")
  blocking <- nrow(problems) > 0L
  problems <- rbind(problems, .ig_advisor_guarded(.ig_advisor_content_problems(rs), "schema",
                                                  "the content could not be checked"))
  w <- NULL
  if (isTRUE(walk) && !blocking) {
    rs <- .ig_advisor_compile(rs)
    w <- tryCatch(.ig_advisor_enumerate(rs), error = function(e) e)
    if (inherits(w, "error")) {
      problems <- rbind(problems, .ig_advisor_problem_table(
        "walk_failed", paste("the answer paths could not be enumerated:", conditionMessage(w))))
      rownames(problems) <- NULL
      return(structure(list(valid = FALSE, problems = problems, walk = NULL, versions = rs$versions),
                       class = "advisor_rules_validation"))
    }
    add <- function(check, msg) {
      if (length(msg)) problems <<- rbind(problems, .ig_advisor_problem_table(rep(check, length(msg)), msg))
    }
    for (type in names(w$issue_counts)) {
      if (w$issue_counts[[type]] > 0L) {
        add(type, sprintf("%d path(s); e.g. %s", w$issue_counts[[type]], w$issues[[type]]))
      }
    }
    add("unreachable_outcome", sprintf("outcome '%s' is not reached by any answer path", w$unreachable_outcomes))
    add("unreachable_question", sprintf("question '%s' is not reached by any answer path", w$unreachable_questions))
    add("unused_transition", sprintf("transition '%s' never fires", w$unused_transitions))
    dc <- rs$docs$rules$summary$draft_check
    documented <- list(answer_paths = dc$answer_paths_enumerated,
                       max_questions_on_a_path = dc$max_questions_on_a_path,
                       outcome_status_combinations = dc$outcome_status_combinations_reachable)
    for (nm in names(documented)) {
      if (!identical(as.integer(documented[[nm]]), as.integer(w[[nm]]))) {
        add("summary_counts", sprintf("walk found %s = %d but summary.draft_check documents %s",
                                      nm, w[[nm]], format(documented[[nm]] %||% "nothing")))
      }
    }
  } else if (isTRUE(walk)) {
    problems <- rbind(problems, .ig_advisor_problem_table(
      "walk_skipped", "answer paths were not enumerated because of structural problems"))
  }
  rownames(problems) <- NULL
  structure(list(valid = nrow(problems) == 0L, problems = problems, walk = w,
                 versions = rs$versions),
            class = "advisor_rules_validation")
}

#' @export
format.advisor_rules_validation <- function(x, ...) {
  head <- sprintf("<advisor_rules_validation> %s: %d problem(s)",
                  if (isTRUE(x$valid)) "valid" else "INVALID", nrow(x$problems))
  walk <- if (!is.null(x$walk)) {
    sprintf("  %d answer paths, max %d questions, %d outcome-status combinations",
            x$walk$answer_paths, x$walk$max_questions_on_a_path, x$walk$outcome_status_combinations)
  }
  probs <- if (nrow(x$problems)) sprintf("  - [%s] %s", x$problems$check, x$problems$message)
  c(head, walk, utils::head(probs, 20L),
    if (length(probs) > 20L) sprintf("  ... and %d more", length(probs) - 20L))
}

#' @export
print.advisor_rules_validation <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}
