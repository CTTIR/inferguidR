# These tests describe the foundation behaviour for methods without an
# engine; they run against the bare registry (see th_registry_only()).
th_registry_only()

th_recommendation <- function(...) {
  args <- th_override(list(
    status = "recommended",
    outcome_id = "out.two_independent_numeric",
    rules = list(rules_id = "advisor_rules", rules_version = "2026-10-03.2",
                 uzh_tree_version = "2026-10-03.1", descrtab2_tree_version = NULL,
                 analysis_variants_version = "2026-10-03.2"),
    answers = list(list(question_id = "q.aim", option_ids = "compare"),
                   list(question_id = "q.stage", option_ids = "planning")),
    rule_trace = list(list(step = 1L, question_id = "q.aim", transition_id = "t.aim.04", to = "q.stage"),
                      list(step = 2L, question_id = "q.stage", transition_id = "t.stage.01",
                           to = "out.two_independent_numeric")),
    rationale = list(en = "Two independent groups with a numeric outcome.", de = NULL),
    estimand = list(en = "Difference in means.", de = "Mittelwertsdifferenz."),
    preferred = list(analysis_variant_id = "t_welch_two_sample"),
    assumptions = list(list(assumption_id = "a.independence", state = "met",
                            text_en = "Observations are independent.")),
    alternatives = list(list(analysis_variant_id = "t_student_two_sample",
                             why_not_preferred_en = "Assumes equal variances.")),
    planning = list(planning_compatibility = "gpower_comparator_different_method",
                    planning_availability = "planned_independent_validation",
                    gpower_method_id = "t_means_two_independent")
  ), list(...))
  do.call(.ig_new_recommendation, args)
}

test_that("a complete recommendation is accepted and hashed", {
  rec <- th_recommendation()
  expect_s3_class(rec, "recommendation")
  expect_match(ig_hash(rec), "^[0-9a-f]{64}$")
  expect_output(print(rec), "t_welch_two_sample")
})

test_that("recommendations reference specifications by hash, never results", {
  spec <- th_minimal_spec("t_means_two_independent", "a_priori")
  rec <- th_recommendation(planning = list(planning_compatibility = "not_yet_validated",
                                           planning_availability = "planned_independent_validation",
                                           specs = list(spec)))
  expect_identical(rec$planning$spec_refs, spec$spec_hash)
  res <- calculate_power(spec)
  expect_identical(th_condition_code(th_recommendation(design = list(result = res))), "invalid_input")
  expect_identical(th_condition_code(th_recommendation(
    planning = list(planning_compatibility = "not_yet_validated",
                    planning_availability = "unavailable", specs = list(res)))), "invalid_input")
})

test_that("status and planning claims must be consistent and truthful", {
  expect_identical(th_condition_code(th_recommendation(preferred = NULL)), "invalid_input")
  expect_identical(th_condition_code(th_recommendation(status = "needs_information")), "invalid_input")
  needs <- th_recommendation(status = "needs_information", preferred = NULL,
                             missing_information = "q.structure",
                             planning = list(planning_compatibility = "not_applicable",
                                             planning_availability = "not_applicable"))
  expect_identical(needs$status, "needs_information")
  expect_identical(th_condition_code(th_recommendation(
    status = "outside_supported_scope", preferred = NULL)), "invalid_input")
  expect_identical(th_condition_code(th_recommendation(
    planning = list(planning_compatibility = "exact_gpower_match",
                    planning_availability = "planned_with_gpower_reference_candidate"))),
    "invalid_input")
  expect_identical(th_condition_code(th_recommendation(rule_trace = list())), "invalid_input")
  expect_identical(th_condition_code(th_recommendation(status = "maybe")), "invalid_input")
  expect_identical(th_condition_code(th_recommendation(
    planning = list(planning_compatibility = "not_yet_validated",
                    planning_availability = "unavailable", gpower_method_id = "t_welch"))),
    "invalid_input")
  expect_identical(th_condition_code(th_recommendation(
    assumptions = list(list(assumption_id = "a.x", state = "probably", text_en = "x")))),
    "invalid_input")
  # Extensions are registry methods but not G*Power procedures.
  for (id in th_extension_ids) {
    cond <- tryCatch(th_recommendation(
      planning = list(planning_compatibility = "independently_validated_extension",
                      planning_availability = "planned_independent_validation",
                      gpower_method_id = id)),
      inferguidR_invalid_input = function(e) e)
    expect_match(unlist(cond$status$details$problems), "not a G\\*Power procedure", all = FALSE)
  }
  expect_identical(th_condition_code(th_recommendation(
    planning = list(planning_compatibility = "not_yet_validated",
                    planning_availability = "unavailable", note_de = c("a", "b")))),
    "invalid_input")
})

test_that("planning notes are kept in both languages", {
  rec <- th_recommendation(planning = list(
    planning_compatibility = "not_applicable", planning_availability = "not_applicable",
    note_en = "No planning is defined for this outcome.",
    note_de = "F\u00fcr dieses Ergebnis ist keine Planung festgelegt."),
    status = "needs_information", preferred = NULL, missing_information = "q.structure")
  expect_identical(rec$planning$note_de, "F\u00fcr dieses Ergebnis ist keine Planung festgelegt.")
  back <- ig_from_json(ig_to_json(rec))
  expect_identical(back, rec)
  expect_identical(back$planning$note_de, rec$planning$note_de)
  lst <- .ig_recommendation_to_list(rec)
  expect_true("note_de" %in% names(lst$planning))
  lst$planning$note_de <- "Ge\u00e4ndert."
  expect_identical(th_condition_code(ig_from_json(ig_to_json(lst))), "invalid_input")
  if (requireNamespace("jsonvalidate", quietly = TRUE)) {
    expect_true(th_schema_valid(rec, "recommendation"))
    expect_true(th_schema_valid(th_recommendation(), "recommendation"))
  }
})

test_that("recommendations round-trip through JSON and tampering is detected", {
  rec <- th_recommendation()
  expect_identical(ig_from_json(ig_to_json(rec)), rec)
  lst <- .ig_recommendation_to_list(rec)
  lst$rationale$en <- "Edited."
  expect_identical(th_condition_code(ig_from_json(ig_to_json(lst))), "invalid_input")
  expect_identical(th_condition_code(rec$status <- "conditional"), "invalid_input")
})
