# Advisor sessions: answering, backtracking, overrides, language, deep links,
# save/restore and the allowlisted condition interpreter.

ta_welch_map <- c(q.aim = "compare", q.stage = "analysis", q.outcome_type = "continuous",
                  q.unit = "individuals", q.structure = "independent", q.independence_check = "no",
                  q.n_factors = "one", q.n_groups = "two", q.covariates = "no",
                  q.estimand_continuous = "mean", q.variance = "not_justified",
                  q.direction = "two_sided", q.missingness = "none")

ta_err <- function(expr) tryCatch({
  force(expr)
  NULL
}, inferguidR_error = function(e) e)

test_that("a new session starts at the first question and asks only what applies", {
  s <- advisor_session()
  expect_s3_class(s, "advisor_session")
  expect_identical(s$state$pending_question_id, "q.aim")
  expect_null(s$state$outcome_id)
  q <- advisor_question(s)
  expect_s3_class(q, "advisor_question")
  expect_true(q$is_current)
  expect_true("dont_know" %in% vapply(q$options, `[[`, "", "option_id"))
  dk <- Filter(function(o) o$dont_know, q$options)[[1]]
  expect_true(nzchar(dk$explanation))
  # A question that does not apply yet cannot be answered
  e <- ta_err(advisor_answer(s, "yes", question_id = "q.causal_randomised"))
  expect_s3_class(e, "inferguidR_invalid_input")
  e <- ta_err(advisor_answer(s, "maybe"))
  expect_s3_class(e, "inferguidR_invalid_input")
  e <- ta_err(advisor_answer(s, "compare", question_id = "q.nonexistent"))
  expect_s3_class(e, "inferguidR_invalid_input")
  # Sessions are immutable
  expect_identical(th_condition_code(s$language <- "de"), "invalid_input")
  expect_identical(th_condition_code(advisor_session(language = "fr")), "invalid_input")
})

test_that("raw and effective answers are both kept in the trace", {
  m <- ta_welch_map
  m[["q.variance"]] <- "dont_know"
  s <- advisor_replay(m)
  st <- Filter(function(x) x$question_id == "q.variance", s$state$steps)[[1]]
  expect_identical(st$raw_answer, "dont_know")
  expect_identical(st$effective_answer, "not_justified")
  expect_identical(st$dont_know_policy, "route_without_assumption")
  rec <- advisor_recommendation(s)
  a <- Filter(function(x) x$question_id == "q.variance", rec$answers)[[1]]
  expect_identical(a$option_ids, "dont_know")
  expect_identical(a$effective_option_id, "not_justified")
  expect_identical(rec$design$answers[["q.variance"]], "not_justified")
  expect_identical(rec$status, "recommended") # Welch is valid either way: no downgrade
})

test_that("'I do not know' never invents a default", {
  # needs_information policy: stop, say what to find out
  s <- advisor_answer(advisor_session(), "dont_know")
  rec <- advisor_recommendation(s)
  expect_identical(rec$outcome_id, "out.needs_information")
  expect_identical(rec$status, "needs_information")
  expect_identical(rec$missing_information, "q.aim")
  expect_null(rec$preferred)
  # route_with_assumption: routed, assumption recorded, downgraded
  m <- ta_welch_map
  m[["q.estimand_continuous"]] <- "dont_know"
  rec <- advisor_recommendation(advisor_replay(m))
  expect_identical(rec$outcome_id, "out.two_group.welch")
  expect_identical(rec$status, "conditional")
  made <- Filter(function(a) identical(a$source, "dont_know"), rec$assumptions)
  expect_identical(vapply(made, `[[`, "", "question_id"), "q.estimand_continuous")
  expect_true(nzchar(made[[1]]$text_en) && nzchar(made[[1]]$text_de))
  expect_true("sr.dont_know_downgrade" %in% unlist(rec$rules$status_rules_applied))
  expect_identical(rec$missing_information, "q.estimand_continuous")
})

test_that("changing an earlier answer invalidates every later answer", {
  s <- advisor_replay(ta_welch_map)
  expect_identical(s$state$outcome_id, "out.two_group.welch")
  s2 <- advisor_answer(s, "paired", question_id = "q.structure")
  expect_identical(s2$invalidated$question_ids,
                   c("q.independence_check", "q.n_factors", "q.n_groups", "q.covariates",
                     "q.estimand_continuous", "q.variance", "q.direction", "q.missingness"))
  expect_identical(s2$state$pending_question_id, "q.pairing_id")
  expect_identical(vapply(s2$answers, `[[`, "", "question_id"),
                   c("q.aim", "q.stage", "q.outcome_type", "q.unit", "q.structure"))
  expect_identical(s2$answers[[5]]$option_id, "paired")
  # The invalidated questions can no longer be answered out of order
  expect_s3_class(ta_err(advisor_answer(s2, "no", question_id = "q.independence_check")),
                  "inferguidR_invalid_input")
  # The same answer again is not a change
  expect_identical(advisor_answer(s, "independent", question_id = "q.structure"), s)
  # Changing the aim invalidates everything after it
  s3 <- advisor_answer(s, "describe", question_id = "q.aim")
  expect_identical(s3$invalidated$question_ids, names(ta_welch_map)[-1])
  expect_identical(s3$state$outcome_id, "out.describe")
  # A completed walkthrough needs question_id to change an answer
  expect_s3_class(ta_err(advisor_answer(s, "none")), "inferguidR_invalid_input")
})

test_that("advisor_back retracts answers", {
  s <- advisor_replay(ta_welch_map)
  b1 <- advisor_back(s)
  expect_identical(b1$state$pending_question_id, "q.missingness")
  expect_identical(b1$invalidated$question_ids, "q.missingness")
  b2 <- advisor_back(s, to = "q.n_groups")
  expect_identical(b2$state$pending_question_id, "q.n_groups")
  expect_identical(b2$invalidated$question_ids, names(ta_welch_map)[8:13])
  expect_s3_class(ta_err(advisor_back(advisor_session())), "inferguidR_invalid_input")
  expect_s3_class(ta_err(advisor_back(b2, to = "q.variance")), "inferguidR_invalid_input")
})

test_that("overrides downgrade the status and are recorded with the why-not text", {
  s <- advisor_replay(ta_welch_map)
  rec0 <- advisor_recommendation(s)
  expect_identical(rec0$status, "recommended")
  o <- advisor_override(s, "rank_sum_mann_whitney")
  rec <- advisor_recommendation(o)
  expect_identical(rec$status, "conditional")
  expect_identical(rec$outcome_id, "out.two_group.welch")
  expect_identical(rec$preferred$analysis_variant_id, "rank_sum_mann_whitney")
  expect_identical(rec$preferred$outcome_variant_id, "t_means_two_independent_welch")
  expect_true(rec$preferred$overridden)
  expect_identical(rec$preferred$override$status_cap, "conditional")
  why <- Filter(function(a) a$analysis_variant_id == "rank_sum_mann_whitney",
                advisor_rules()$outcomes[["out.two_group.welch"]]$alternatives)[[1]]$why_not_preferred_en
  expect_identical(rec$preferred$override$why_not_preferred_en, why)
  ov <- Filter(function(a) identical(a$source, "override"), rec$assumptions)
  expect_identical(ov[[1]]$text_en, why)
  expect_true("override_policy" %in% unlist(rec$rules$status_rules_applied))
  # The overriding variant brings its own assumptions; the replaced outcome's
  # assumptions stay listed with the variant they were stated for
  own <- Filter(function(a) identical(a$source, "override_variant"), rec$assumptions)
  expect_true("a.rank_null_identical" %in% vapply(own, `[[`, "", "assumption_id"))
  expect_true(all(vapply(own, `[[`, "", "from_outcome_id") == "out.two_group.wmw"))
  stated <- Filter(function(a) identical(a$source, "outcome"), rec$assumptions)
  expect_true(all(vapply(stated, `[[`, "", "stated_for") == "t_means_two_independent_welch"))
  ids <- vapply(rec$assumptions, `[[`, "", "assumption_id")
  expect_false(anyDuplicated(ids) > 0L)
  # A variant that is no outcome's primary variant: the gap is stated explicitly
  bos_map <- c(ta_welch_map[c("q.aim", "q.stage", "q.unit", "q.structure", "q.independence_check",
                              "q.n_factors", "q.n_groups", "q.covariates", "q.direction", "q.missingness")],
               q.outcome_type = "binary", q.binary_margins = "group_sizes_fixed")
  ml <- advisor_recommendation(advisor_override(advisor_replay(bos_map), "boschloo_two_sided_minlike"))
  gap <- Filter(function(a) identical(a$source, "override_variant"), ml$assumptions)
  expect_identical(vapply(gap, `[[`, "", "assumption_id"),
                   "override.boschloo_two_sided_minlike.assumptions_not_stated")
  expect_match(gap[[1]]$text_en, "state no method-specific assumptions")
  expect_true(nzchar(gap[[1]]$text_de))
  expect_true("a.ordering_ties" %in% vapply(ml$assumptions, `[[`, "", "assumption_id"))
  # The replaced variant is listed as an alternative; the override is not
  alt_ids <- vapply(rec$alternatives, `[[`, "", "analysis_variant_id")
  expect_true("t_means_two_independent_welch" %in% alt_ids)
  expect_false("rank_sum_mann_whitney" %in% alt_ids)
  # Code export follows the override; planning is not borrowed from another design
  expect_identical(rec$r_implementation$template_id, "tpl.rank_sum_mann_whitney.v1")
  expect_identical(rec$planning$planning_compatibility, "not_yet_validated")
  expect_identical(rec$planning$planning_availability, "unavailable")
  expect_match(rec$estimand$en, "P(X > Y)", fixed = TRUE)
  expect_match(rec$rationale$en, "User override")
  # A conditional outcome stays conditional
  m <- ta_welch_map
  m[["q.variance"]] <- "justified"
  stud <- advisor_replay(m)
  expect_identical(advisor_recommendation(stud)$status, "conditional")
  expect_identical(advisor_recommendation(advisor_override(stud, "t_means_two_independent_welch"))$status,
                   "conditional")
  # Overriding with the outcome's own variant keeps no override
  same <- advisor_override(s, "t_means_two_independent_welch")
  expect_null(same$override)
  # Clearing restores the recommendation
  cl <- advisor_clear_override(o)
  expect_null(cl$override)
  expect_identical(advisor_recommendation(cl), rec0)
  expect_identical(cl$invalidated$override, "rank_sum_mann_whitney")
})

test_that("overrides are refused for rejected, unlisted, disallowed or missing variants", {
  s <- advisor_replay(ta_welch_map)
  refusal <- function(expr) {
    e <- ta_err(expr)
    expect_s3_class(e, "inferguidR_override_refused")
    expect_s3_class(e, "inferguidR_invalid_input")
    e$status$details$refusal
  }
  expect_identical(refusal(advisor_override(s, "boschloo_two_sided_central")), "rejected_route")
  expect_identical(refusal(advisor_override(s, "anova_one_way")), "not_listed")
  expect_identical(refusal(advisor_override(s, "ancova_baseline_adjusted")), "override_not_allowed")
  expect_identical(refusal(advisor_override(s, "no_such_variant")), "unknown_variant")
  expect_identical(refusal(advisor_override(advisor_session(), "t_means_two_independent_welch")), "no_method")
  out <- advisor_replay(c(q.aim = "causal", q.causal_randomised = "no"))
  expect_identical(refusal(advisor_override(out, "regression_linear_multiple")), "no_method")
  # Paired binary data: the independent Boschloo route is refused with its reason
  mc <- advisor_replay(c(q.aim = "compare", q.stage = "analysis", q.outcome_type = "binary",
                         q.unit = "individuals", q.structure = "paired", q.pairing_id = "yes",
                         q.direction = "two_sided", q.missingness = "none"))
  e <- ta_err(advisor_override(mc, "boschloo_one_sided"))
  expect_identical(e$status$details$refusal, "rejected_route")
  expect_match(e$status$details$reason_en, "independent")
})

test_that("changing an answer invalidates an override", {
  s <- advisor_override(advisor_replay(ta_welch_map), "t_means_two_independent_student")
  s2 <- advisor_answer(s, "some", question_id = "q.missingness")
  expect_null(s2$override)
  expect_identical(s2$invalidated$override, "t_means_two_independent_student")
  expect_identical(s2$invalidated$question_ids, character())
  rec <- advisor_recommendation(s2)
  expect_identical(rec$preferred$analysis_variant_id, "t_means_two_independent_welch")
  expect_identical(rec$status, "conditional") # missingness modifier
  expect_identical(unlist(rec$rules$status_modifiers_applied), "sm.missing_complete_case")
})

test_that("the language changes texts only", {
  s <- advisor_replay(ta_welch_map[1:3])
  de <- advisor_set_language(s, "de")
  expect_identical(de$state, s$state)
  expect_identical(advisor_question(de)$text, advisor_rules()$questions[["q.unit"]]$text_de)
  expect_identical(advisor_question(s)$text, advisor_rules()$questions[["q.unit"]]$text_en)
  expect_identical(advisor_question(s, language = "de")$text, advisor_question(de)$text)
  expect_identical(advisor_recommendation(de), advisor_recommendation(s))
})

test_that("deep links replay in transition order, ignore unreached answers and reject unknown ids", {
  m <- c(ta_welch_map, q.pairing_id = "yes", q.binary_margins = "both_margins_fixed")
  s <- advisor_replay(m)
  expect_identical(s$state$outcome_id, "out.two_group.welch")
  expect_setequal(s$ignored_answers, c("q.pairing_id", "q.binary_margins"))
  # order in the map does not matter
  s_rev <- advisor_replay(rev(ta_welch_map))
  expect_identical(s_rev$state, advisor_replay(ta_welch_map)$state)
  # a list works too; an incomplete map leaves the next question pending
  p <- advisor_replay(as.list(ta_welch_map[1:6]))
  expect_identical(p$state$pending_question_id, "q.n_factors")
  rec <- advisor_recommendation(p)
  expect_identical(rec$status, "needs_information")
  expect_identical(rec$design$pending_question_id, "q.n_factors")
  for (bad in list(c(q.aim = "compare", q.structure = "crossover"),
                   c(q.aim = "compare", q.unknown = "x"),
                   c(q.aim = "compare", q.aim = "describe"),
                   c("compare"))) {
    expect_s3_class(ta_err(advisor_replay(bad)), "inferguidR_invalid_input")
  }
  expect_identical(advisor_replay(character())$state$pending_question_id, "q.aim")
})

test_that("save and restore reproduce the session exactly", {
  sessions <- list(
    advisor_session(),
    advisor_replay(ta_welch_map[1:5]),
    advisor_replay(ta_welch_map),
    advisor_override(advisor_replay(ta_welch_map), "t_means_two_independent_student"),
    advisor_set_language(advisor_replay(c(ta_welch_map, q.binary_margins = "both_margins_fixed")), "de"),
    advisor_answer(advisor_replay(ta_welch_map), "paired", question_id = "q.structure")
  )
  for (s in sessions) {
    js <- advisor_save(s)
    expect_s3_class(js, "json")
    back <- advisor_restore(js)
    expect_identical(back, s)
    expect_identical(advisor_recommendation(back), advisor_recommendation(s))
  }
  # Through a file
  f <- tempfile(fileext = ".json")
  on.exit(unlink(f), add = TRUE)
  advisor_save(sessions[[4]], file = f)
  expect_identical(advisor_restore(file = f), sessions[[4]])
})

test_that("restore refuses tampered, foreign, future and other-version states", {
  s <- advisor_replay(ta_welch_map)
  js <- advisor_save(s)
  doc <- jsonlite::parse_json(js, simplifyVector = FALSE)
  reseal <- function(d) {
    d$state_hash <- NULL
    d$state_hash <- ig_hash(d)
    ig_to_json(d)
  }
  # tampering without resealing
  t1 <- doc
  t1$answers[[13]]$option_id <- "some"
  expect_match(conditionMessage(ta_err(advisor_restore(ig_to_json(t1)))), "state_hash")
  # resealed, but the replay no longer gives the saved outcome and status
  expect_match(conditionMessage(ta_err(advisor_restore(reseal(t1)))), "replaying the saved answers")
  # answers that do not follow the graph
  t2 <- doc
  t2$answers[[2]]$question_id <- "q.variance"
  expect_s3_class(ta_err(advisor_restore(reseal(t2))), "inferguidR_invalid_input")
  # future schema version, wrong object type, unknown fields
  t3 <- doc
  t3$schema_version <- "inferguidR.advisor_session/9.0.0"
  expect_match(conditionMessage(ta_err(advisor_restore(reseal(t3)))), "schema_version")
  t4 <- doc
  t4$object_type <- "power_spec"
  expect_s3_class(ta_err(advisor_restore(reseal(t4))), "inferguidR_invalid_input")
  t5 <- doc
  t5$extra <- "x"
  expect_match(conditionMessage(ta_err(advisor_restore(reseal(t5)))), "unknown field")
  # other rules version: explicit migration required, with the change log
  t6 <- doc
  t6$rules$rules_version <- "2026-10-03.1"
  e <- ta_err(advisor_restore(reseal(t6)))
  expect_s3_class(e, "inferguidR_migration_required")
  expect_s3_class(e, "inferguidR_invalid_input")
  expect_identical(e$status$details$current_rules_version, "2026-10-26.1")
  expect_true(any(vapply(e$status$details$change_log, function(x) nzchar(x$migration_en), logical(1))))
  # same version label, different content
  t7 <- doc
  t7$rules$rules_sha256 <- strrep("0", 64)
  expect_match(conditionMessage(ta_err(advisor_restore(reseal(t7)))), "content hash")
  # duplicate keys and non-JSON input are refused; remote paths are refused
  expect_s3_class(ta_err(advisor_restore('{"a":1,"a":2}')), "inferguidR_invalid_input")
  expect_s3_class(ta_err(advisor_restore("not json")), "inferguidR_invalid_input")
  expect_s3_class(ta_err(advisor_restore(file = "https://example.org/state.json")), "inferguidR_invalid_input")
  expect_s3_class(ta_err(advisor_restore()), "inferguidR_invalid_input")
})

test_that("the condition interpreter accepts only the allowlisted grammar", {
  qs <- advisor_rules()$questions
  ok <- list(all = list(list(q = "q.aim", "in" = list("compare")),
                        list(not = list(q = "q.stage", "in" = list("planning"))),
                        list(any = list(list(always = TRUE)))))
  expect_length(.ig_advisor_condition_problems(ok, qs, "t"), 0L)
  bad <- list(
    "q.aim == 'compare'",
    list(eval = "1 + 1"),
    list(q = "q.aim"),
    list(q = "q.aim", "in" = "compare"),
    list(q = "q.aim", "in" = list()),
    list(q = "q.nope", "in" = list("x")),
    list(q = "q.aim", "in" = list("crossover")),
    list(always = FALSE),
    list(all = list()),
    list(all = list(list(q = "q.aim", "in" = list("compare"))), any = list(list(always = TRUE))),
    list(not = "TRUE"),
    list(q = "q.aim", "in" = list("compare"), extra = 1)
  )
  for (b in bad) expect_true(length(.ig_advisor_condition_problems(b, qs, "t")) > 0L, info = deparse(b))
  # A string condition is rejected without ever being evaluated
  marker <- tempfile()
  sneaky <- sprintf("file.create('%s')", marker)
  expect_true(length(.ig_advisor_condition_problems(sneaky, qs, "t")) > 0L)
  expect_error(.ig_advisor_compile_condition(sneaky))
  expect_false(file.exists(marker))
  # Truth table of the compiled forms
  f <- .ig_advisor_compile_condition(ok)
  expect_true(f(c(q.aim = "compare")))
  expect_true(f(c(q.aim = "compare", q.stage = "analysis")))
  expect_false(f(c(q.aim = "compare", q.stage = "planning")))
  expect_false(f(c(q.stage = "analysis")))
  g <- .ig_advisor_compile_condition(list(q = "q.aim", "in" = list("compare", "describe")))
  expect_false(g(stats::setNames(character(), character()))) # unanswered is false
  expect_true(g(c(q.aim = "describe")))
  expect_true(.ig_advisor_compile_condition(list(always = TRUE))(character()))
})

test_that("malformed ids and unusable paths give typed conditions, never untyped errors", {
  s <- advisor_replay(ta_welch_map)
  # advisor_back(to =) takes one question id
  for (to in list(c("q.aim", "q.stage"), 3, NA_character_, "")) {
    expect_s3_class(ta_err(advisor_back(s, to = to)), "inferguidR_invalid_input")
  }
  # advisor_question() never selects a question by position or recursively
  for (qid in list(2, 2L, c("q.aim", "q.stage"), NA_character_, TRUE)) {
    expect_s3_class(ta_err(advisor_question(s, question_id = qid)), "inferguidR_invalid_input")
  }
  expect_identical(advisor_question(s, question_id = "q.stage")$question_id, "q.stage")
  # directories, missing folders and remote paths are refused as files
  dir <- tempfile("advisor-dir-")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  typed_quietly <- function(expr) {
    e <- tryCatch(force(expr), warning = function(w) w, error = function(e) e)
    expect_s3_class(e, "inferguidR_invalid_input")
  }
  typed_quietly(advisor_restore(file = dir))
  typed_quietly(advisor_save(s, file = dir))
  typed_quietly(advisor_save(s, file = file.path(dir, "no-such-folder", "state.json")))
  typed_quietly(advisor_save(s, file = "https://example.org/state.json"))
  typed_quietly(advisor_rules(rules = dir))
  v <- tryCatch(advisor_validate_rules(rules = dir, walk = FALSE), warning = function(w) w)
  expect_s3_class(v, "advisor_rules_validation")
  expect_false(v$valid)
  expect_identical(unique(v$problems$check), "schema")
  # an empty file is a typed parse failure
  empty <- file.path(dir, "empty.json")
  file.create(empty)
  typed_quietly(advisor_restore(file = empty))
})

test_that("printing and formatting work", {
  s <- advisor_replay(ta_welch_map)
  expect_output(print(s), "out.two_group.welch")
  expect_output(print(advisor_question(advisor_session())), "q.aim")
  expect_output(print(advisor_rules()), "2026-10-26.1")
  expect_output(print(advisor_recommendation(s)), "t_means_two_independent_welch")
  expect_null(advisor_question(s))
})
