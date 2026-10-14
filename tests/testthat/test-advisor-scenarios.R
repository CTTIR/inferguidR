# The acceptance scenarios of DIGITAL_ADVISOR.md (specs/advisor/acceptance-scenarios.json,
# installed as inst/advisor/acceptance-scenarios.json), replayed through the session API.

ta_scenarios <- function() {
  path <- system.file("advisor", "acceptance-scenarios.json", package = "inferguidR")
  jsonlite::parse_json(paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = "\n"),
                       simplifyVector = FALSE)
}

# Answer the questions one by one through advisor_answer(); every answer must
# be for the question the session asks next.
ta_play <- function(answers, session = advisor_session()) {
  for (a in answers) {
    expect_identical(session$state$pending_question_id, a$question_id)
    session <- advisor_answer(session, a$answer)
  }
  session
}

ta_ids <- function(records, key = "analysis_variant_id") vapply(records, function(r) r[[key]], "")

# Check one `expected` block against a recommendation; every key must be known.
ta_check_expected <- function(exp, rec, rs, label) {
  known <- c("outcome_id", "status", "analysis_variant_id", "planning_compatibility",
             "planning_availability", "r_template_id", "assumptions_made_from", "estimand_en_contains",
             "alternatives_include", "assumption_ids_include", "status_modifiers_applied",
             "rejected_routes_include", "pending_question_id", "destination", "next_actions_include")
  expect_true(all(names(exp) %in% known), info = paste(label, "unknown expected key(s):",
                                                       paste(setdiff(names(exp), known), collapse = ", ")))
  if (!is.null(exp$outcome_id)) expect_identical(rec$outcome_id, exp$outcome_id, info = label)
  if (!is.null(exp$status)) expect_identical(rec$status, exp$status, info = label)
  if (!is.null(exp$analysis_variant_id)) {
    expect_identical(rec$preferred$analysis_variant_id, exp$analysis_variant_id, info = label)
  }
  if (!is.null(exp$planning_compatibility)) {
    expect_identical(rec$planning$planning_compatibility, exp$planning_compatibility, info = label)
  }
  if (!is.null(exp$planning_availability)) {
    expect_identical(rec$planning$planning_availability, exp$planning_availability, info = label)
  }
  if (!is.null(exp$r_template_id)) expect_identical(rec$r_implementation$template_id, exp$r_template_id, info = label)
  if ("assumptions_made_from" %in% names(exp)) {
    made <- Filter(function(a) identical(a$source, "dont_know"), rec$assumptions)
    expect_setequal(vapply(made, `[[`, "", "question_id"), as.character(unlist(exp$assumptions_made_from)))
  }
  if (!is.null(exp$estimand_en_contains)) {
    expect_true(grepl(exp$estimand_en_contains, rec$estimand$en, fixed = TRUE), info = label)
  }
  if (!is.null(exp$alternatives_include)) {
    expect_true(all(unlist(exp$alternatives_include) %in% ta_ids(rec$alternatives)), info = label)
  }
  if (!is.null(exp$assumption_ids_include)) {
    expect_true(all(unlist(exp$assumption_ids_include) %in% ta_ids(rec$assumptions, "assumption_id")), info = label)
  }
  if (!is.null(exp$status_modifiers_applied)) {
    expect_setequal(as.character(unlist(rec$rules$status_modifiers_applied)),
                    as.character(unlist(exp$status_modifiers_applied)))
    applied <- Filter(function(a) identical(a$source, "status_modifier"), rec$assumptions)
    expect_setequal(ta_ids(applied, "assumption_id"), as.character(unlist(exp$status_modifiers_applied)))
  }
  if (!is.null(exp$rejected_routes_include)) {
    expect_true(all(unlist(exp$rejected_routes_include) %in% ta_ids(rec$preferred$rejected_routes)), info = label)
  }
  if (!is.null(exp$pending_question_id)) {
    expect_identical(rec$design$pending_question_id, exp$pending_question_id, info = label)
    expect_identical(rec$missing_information, exp$pending_question_id, info = label)
    last <- rec$rule_trace[[length(rec$rule_trace)]]
    expect_identical(last$transition_id, "engine.incomplete_policy", info = label)
  }
  if (!is.null(exp$destination)) expect_identical(rec$design$destination, exp$destination, info = label)
  if (!is.null(exp$next_actions_include)) {
    expect_true(all(unlist(exp$next_actions_include) %in% rec$next_actions), info = label)
  }
}

# Operations on a played scenario.
ta_run_operation <- function(op, session, rec, rs, label) {
  switch(op$op,
         change_answer = {
           s2 <- advisor_answer(session, op$answer, question_id = op$question_id)
           expect_setequal(s2$invalidated$question_ids, as.character(unlist(op$expect_invalidated)))
           expect_identical(s2$state$pending_question_id, op$then_answers[[1]]$question_id, info = label)
           ta_play(op$then_answers, s2)
         },
         override = {
           res <- tryCatch(advisor_override(session, op$analysis_variant_id),
                           inferguidR_override_refused = function(e) e)
           if (identical(op$expect, "allowed")) {
             expect_s3_class(res, "advisor_session")
             r2 <- advisor_recommendation(res)
             expect_identical(r2$preferred$analysis_variant_id, op$analysis_variant_id, info = label)
             expect_true(r2$status %in% c("conditional", "needs_information", "outside_supported_scope"),
                         info = label) # never better than conditional
             if (!is.null(op$expect_status)) expect_identical(r2$status, op$expect_status, info = label)
             expect_true(isTRUE(r2$preferred$overridden))
             expect_true(paste0("override.", op$analysis_variant_id) %in% ta_ids(r2$assumptions, "assumption_id"))
             expect_identical(r2$r_implementation$template_id, paste0("tpl.", op$analysis_variant_id, ".v1"))
           } else {
             expect_s3_class(res, "inferguidR_override_refused")
             expect_s3_class(res, "inferguidR_invalid_input")
             want <- c(refused_rejected_route = "rejected_route", refused_not_listed = "not_listed")[[op$expect]]
             expect_identical(res$status$details$refusal, want, info = label)
             expect_true(nzchar(res$status$details$reason_en))
           }
           session
         },
         language = {
           de <- advisor_set_language(session, op$language)
           expect_identical(advisor_recommendation(de), rec) # identical decision and record
           for (st in de$state$steps) {
             q <- advisor_question(de, st$question_id)
             expect_identical(q$language, "de")
             expect_true(nzchar(q$text) && nzchar(q$why_it_matters), info = st$question_id)
             expect_true(all(vapply(q$options, function(o) nzchar(o$label), logical(1))), info = st$question_id)
             expect_identical(q$text, rs$questions[[st$question_id]]$text_de)
           }
           expect_true(nzchar(rec$rationale$de) && nzchar(rec$estimand$de))
           expect_true(nzchar(rec$preferred$outcome_title_de))
           session
         },
         save_restore = {
           js <- advisor_save(session)
           if (identical(op$rules_version, "same")) {
             back <- advisor_restore(js)
             expect_identical(back, session, info = label)
             expect_identical(advisor_recommendation(back), rec)
           } else {
             doc <- jsonlite::parse_json(js, simplifyVector = FALSE)
             doc$rules$rules_version <- op$rules_version
             doc$state_hash <- NULL
             doc$state_hash <- ig_hash(doc) # a consistent state saved under another rules version
             err <- tryCatch(advisor_restore(ig_to_json(doc)), inferguidR_error = function(e) e)
             expect_s3_class(err, "inferguidR_migration_required")
             expect_identical(err$status$details$reason_code, op$expect)
             expect_identical(err$status$details$saved_rules_version, op$rules_version)
           }
           session
         },
         code_export = {
           expect_identical(rec$r_implementation$template_id, op$expect_r_template_id, info = label)
           # The bundled code must use the exact selected variant.
           expect_identical(rec$r_implementation$template_status, "available")
           sketch <- advisor_r_code(session)
           expect_true(any(grepl("stats::t.test(", sketch, fixed = TRUE) &
                             grepl("var.equal = FALSE", sketch, fixed = TRUE)), info = op$expect_template_calls)
           session
         },
         deep_link = {
           res <- tryCatch(advisor_replay(unlist(op$answers_map)), inferguidR_invalid_input = function(e) e)
           if (identical(op$expect, "invalid_input")) {
             expect_s3_class(res, "inferguidR_invalid_input")
           } else {
             expect_s3_class(res, "advisor_session")
             r2 <- advisor_recommendation(res)
             expect_identical(r2$outcome_id, rec$outcome_id, info = label)
             expect_identical(r2$status, rec$status, info = label)
             expect_setequal(res$ignored_answers, as.character(unlist(op$expect_ignored)))
             expect_identical(as.character(unlist(r2$design$ignored_answers)),
                              as.character(unlist(op$expect_ignored)))
           }
           session
         },
         planning_action = {
           # No manufactured sample size: no specification, no evidence, an explanation instead.
           expect_identical(rec$planning$planning_availability, "simulation_only_if_implemented")
           expect_length(rec$planning$spec_refs, 0L)
           expect_length(rec$planning$evidence_ids, 0L)
           expect_true(nzchar(rec$planning$note_en))
           expect_true(nzchar(rec$preferred$planning_note_de))
           expect_false(any(grepl("sample_size|n_total", names(unlist(rec$planning)))))
           session
         },
         stop("unknown scenario operation: ", op$op))
}

test_that("all 45 acceptance scenarios pass", {
  doc <- ta_scenarios()
  rs <- advisor_rules()
  expect_identical(doc$rules_version, rs$versions$rules_version)
  expect_length(doc$scenarios, 45L)
  n_engine <- 0L
  for (sc in doc$scenarios) {
    label <- sc$scenario_id
    if (identical(sc$test_level, "static")) next
    n_engine <- n_engine + 1L
    session <- ta_play(sc$answers)
    if (!is.null(sc$language)) session <- advisor_set_language(session, sc$language)
    rec <- advisor_recommendation(session)
    final <- session
    for (op in sc$operations) {
      out <- ta_run_operation(op, session, rec, rs, label)
      if (identical(op$op, "change_answer")) final <- out
    }
    ta_check_expected(sc$expected, advisor_recommendation(final), rs, label)
    # Every recommendation is explainable node by node
    trace <- advisor_recommendation(final)$rule_trace
    expect_true(length(trace) >= 1L, info = label)
    expect_identical(vapply(trace[seq_along(final$state$steps)], `[[`, "", "question_id"),
                     vapply(final$state$steps, `[[`, "", "question_id"))
  }
  expect_identical(n_engine, 44L)
})

test_that("AS21: no normality-pre-test switching rule anywhere in the graph", {
  doc <- ta_scenarios()
  sc <- Filter(function(s) s$scenario_id == "AS21", doc$scenarios)[[1]]
  exp <- sc$expected
  rs <- advisor_rules()
  expect_true(all(vapply(rs$questions, function(q) identical(q$input_kind, exp$all_questions_input_kind),
                         logical(1))))
  expect_true(isTRUE(exp$forbidden_routing_inputs_declared))
  expect_true(length(rs$docs$rules$engine$forbidden_routing_inputs) > 0L)
  terms <- tolower(unlist(exp$forbidden_option_terms))
  expect_setequal(terms, tolower(.ig_advisor_forbidden_terms))
  for (q in rs$questions) {
    for (o in q$options) {
      txt <- tolower(paste(o$option_id, o$label_en, o$label_de))
      for (t in terms) expect_false(grepl(t, txt, fixed = TRUE), info = paste(q$question_id, o$option_id, t))
    }
  }
  # Conditions reference only question and option ids
  qids <- names(rs$questions)
  oids <- unique(unlist(lapply(rs$questions, function(q) vapply(q$options, `[[`, "", "option_id"))))
  conds <- c(lapply(rs$questions, `[[`, "applies_when"), lapply(rs$docs$rules$transitions, `[[`, "when"),
             unlist(lapply(rs$outcomes, function(o) lapply(o$status_modifiers, `[[`, "when")), recursive = FALSE))
  refs <- unique(unlist(lapply(conds, .ig_advisor_condition_refs)))
  expect_true(all(refs %in% c(qids, oids)))
  for (t in terms) expect_false(any(grepl(t, tolower(refs), fixed = TRUE)))
})
