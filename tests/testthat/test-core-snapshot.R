# These tests describe the foundation behaviour for methods without an
# engine; they run against the bare registry (see th_registry_only()).
th_registry_only()

th_envelope <- function(id = "t_means_one_sample", mode = "post_hoc") {
  spec <- th_minimal_spec(id, mode)
  analysis_envelope(spec, calculate_power(spec))
}

test_that("an incomplete-settings snapshot freezes envelopes with identity hashes", {
  env <- th_envelope()
  snap <- .ig_new_report_snapshot(env, report_kind = "incomplete_settings",
                                  created_at = "2026-10-03T12:00:00.000Z",
                                  display_timezone = "Europe/Berlin", title = "Settings")
  expect_s3_class(snap, "report_snapshot")
  expect_identical(snap$completeness, "incomplete")
  expect_match(snap$snapshot_id, "^snap-[0-9a-f]{24}$")
  expect_match(snap$analysis_hash, "^[0-9a-f]{64}$")
  expect_match(snap$snapshot_hash, "^[0-9a-f]{64}$")
  expect_identical(snap$envelopes[[1]], env)
  expect_output(print(snap), "incomplete_settings")
})

test_that("final-result reports are blocked while results are missing or not ok", {
  env <- th_envelope()
  cond <- tryCatch(.ig_new_report_snapshot(env), inferguidR_invalid_input = function(e) e)
  expect_match(unlist(cond$status$details$problems), "incomplete_settings", all = FALSE)
  pending <- analysis_envelope(th_minimal_spec("t_means_one_sample", "post_hoc"))
  expect_identical(th_condition_code(.ig_new_report_snapshot(pending, report_kind = "session")),
                   "invalid_input")
  expect_identical(th_condition_code(.ig_new_report_snapshot(list(env, env),
                                                             report_kind = "incomplete_settings",
                                                             created_at = "2026-10-03T12:00:00Z",
                                                             title = NULL)$report_kind),
                   NA_character_)
  expect_identical(th_condition_code(.ig_new_report_snapshot(list())), "invalid_input")
})

test_that("snapshot identity: analysis hash ignores time, snapshot id does not", {
  env <- th_envelope()
  a <- .ig_new_report_snapshot(env, report_kind = "incomplete_settings",
                               created_at = "2026-10-03T12:00:00.000Z")
  b <- .ig_new_report_snapshot(env, report_kind = "incomplete_settings",
                               created_at = "2026-10-03T12:00:01.000Z")
  expect_identical(a$analysis_hash, b$analysis_hash)
  expect_false(identical(a$snapshot_id, b$snapshot_id))
  a2 <- .ig_new_report_snapshot(env, report_kind = "incomplete_settings",
                                created_at = "2026-10-03T12:00:00.000Z")
  expect_identical(a2$snapshot_id, a$snapshot_id)
  other <- .ig_new_report_snapshot(th_envelope("t_generic"), report_kind = "incomplete_settings",
                                   created_at = "2026-10-03T12:00:00.000Z")
  expect_false(identical(other$analysis_hash, a$analysis_hash))
})

test_that("snapshot arguments are validated", {
  env <- th_envelope()
  bad <- list(list(locale = "fr"), list(created_at = "2026-10-03 12:00:00"),
              list(report_kind = "pdf"), list(texts = list("unnamed")),
              list(plot_data = list(a = list())), list(recommendation = list()))
  for (b in bad) {
    args <- c(list(envelopes = env), th_override(list(report_kind = "incomplete_settings"), b))
    expect_identical(th_condition_code(do.call(.ig_new_report_snapshot, args)), "invalid_input",
                     info = names(b))
  }
  expect_match(.ig_utc_now(), "^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}\\.\\d{3}Z$")
})

test_that("snapshots round-trip through JSON, are immutable and detect tampering", {
  env <- th_envelope()
  rec_spec <- th_minimal_spec("t_means_two_independent", "a_priori")
  rec <- .ig_new_recommendation(
    status = "conditional", outcome_id = "out.x",
    rules = list(rules_id = "advisor_rules", rules_version = "2026-10-03.2"),
    answers = list(list(question_id = "q.aim", option_ids = "compare")),
    rule_trace = list(list(question_id = "q.aim", transition_id = "t.aim.04", to = "out.x")),
    rationale = list(en = "Example."), preferred = list(analysis_variant_id = "t_welch_two_sample"),
    planning = list(planning_compatibility = "not_yet_validated",
                    planning_availability = "planned_independent_validation", specs = rec_spec))
  snap <- .ig_new_report_snapshot(env, recommendation = rec, report_kind = "incomplete_settings",
                                  created_at = "2026-10-03T12:00:00.000Z",
                                  texts = list(summary = "No engine yet."),
                                  plot_data = list(list(plot_type = "none", note = "placeholder-free")))
  expect_identical(ig_from_json(ig_to_json(snap)), snap)
  expect_identical(th_condition_code(snap$envelopes <- list()), "invalid_input")
  lst <- .ig_snapshot_to_list(snap)
  lst$texts$summary <- "Edited."
  expect_identical(th_condition_code(ig_from_json(ig_to_json(lst))), "invalid_input")
  lst <- .ig_snapshot_to_list(snap)
  lst$snapshot_id <- "snap-000000000000000000000000"
  expect_identical(th_condition_code(ig_from_json(ig_to_json(lst))), "invalid_input")
})
