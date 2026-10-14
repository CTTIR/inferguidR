# Exercise the foundation contracts independently of later method providers.
# Restore the live registry and solvers when this file finishes.
th_registry_only()

test_that("bundled schemas are draft 2020-12 documents with unique ids", {
  dir <- system.file("schemas", package = "inferguidR")
  files <- sort(list.files(dir, pattern = "\\.schema\\.json$"))
  expect_setequal(files, paste0(.ig_schema_names(), ".schema.json"))
  ids <- character()
  for (f in files) {
    s <- jsonlite::fromJSON(file.path(dir, f), simplifyVector = FALSE)
    expect_identical(s[["$schema"]], "https://json-schema.org/draft/2020-12/schema", info = f)
    ids <- c(ids, s[["$id"]])
  }
  expect_false(anyDuplicated(ids) > 0)
})

test_that("serialised objects validate against their schemas", {
  skip_if_not_installed("jsonvalidate")
  spec <- th_minimal_spec("exact_proportions_fisher", "sensitivity")
  res <- calculate_power(spec)
  env <- analysis_envelope(spec, res)
  expect_true(th_schema_valid(spec, "power_spec"))
  expect_true(th_schema_valid(res, "power_result"))
  expect_true(th_schema_valid(env, "analysis_envelope"))
  comp <- th_minimal_spec("t_means_one_sample", "compromise")
  expect_true(th_schema_valid(comp, "power_spec"))
  expect_true(th_schema_valid(calculate_power(comp), "power_result"))
  snap <- .ig_new_report_snapshot(env, report_kind = "incomplete_settings",
                                  created_at = "2026-10-03T12:00:00.000Z")
  expect_true(th_schema_valid(snap, "report_snapshot"))
  rec <- .ig_new_recommendation(
    status = "needs_information", outcome_id = "out.needs_information",
    rules = list(rules_id = "advisor_rules", rules_version = "2026-10-03.2"),
    answers = list(list(question_id = "q.aim", option_ids = "dont_know")),
    rule_trace = list(list(step = 1L, question_id = "q.aim", transition_id = "t.aim.03",
                           to = "out.needs_information")),
    rationale = list(en = "The aim was not answered."), missing_information = "q.aim")
  expect_true(th_schema_valid(rec, "recommendation"))
})

test_that("schema validation of definition-based specs and ok results", {
  skip_if_not_installed("jsonvalidate")
  th_local_definition()
  spec <- analysis_spec("t_means_two_independent", "post_hoc",
                        inputs = list(alpha = 0.05, n_total = 40, d = 0.5),
                        options = list(tails = "two"))
  expect_true(th_schema_valid(spec, "power_spec"))
  th_local_solvers(list(post_hoc = function(method, spec, entry) {
    .ig_power_result(spec, ig_status("ok"), solved = list(power = 0.5),
                     group_sizes = c(n1 = 20L, n2 = 20L),
                     quantities = list(critical = c(-2, 2), ncp = Inf),
                     provenance = list(solver = "test_double"))
  }))
  expect_true(th_schema_valid(calculate_power(spec), "power_result"))
})

test_that("schemas reject documents that break the contracts", {
  skip_if_not_installed("jsonvalidate")
  spec <- th_minimal_spec("t_means_one_sample", "post_hoc")
  res <- calculate_power(spec)
  # A number in a result without an ok status.
  lst <- .ig_result_to_list(res)
  lst$solved <- list(power = 0.8)
  expect_false(.ig_validate_schema(ig_to_json(lst), "power_result")$valid)
  # A G*Power validation claim without fixtures.
  lst <- .ig_result_to_list(res)
  lst$validation$gpower_reference_status <- "validated"
  expect_false(.ig_validate_schema(ig_to_json(lst), "power_result")$valid)
  # A non-compromise spec with two targets; an unknown top-level field.
  sl <- .ig_spec_to_list(spec)
  sl$target <- I(c("alpha", "power"))
  expect_false(.ig_validate_schema(ig_to_json(sl), "power_spec")$valid)
  sl <- .ig_spec_to_list(spec)
  sl$extra <- "x"
  expect_false(.ig_validate_schema(ig_to_json(sl), "power_spec")$valid)
  # Non-finite input values are not admissible inputs.
  sl <- .ig_spec_to_list(spec)
  sl$inputs$alpha <- Inf
  expect_false(.ig_validate_schema(ig_to_json(sl), "power_spec")$valid)
  # ig_from_json reports schema failures as invalid_input.
  sl <- .ig_spec_to_list(spec)
  sl$mode <- "simulation"
  cond <- tryCatch(ig_from_json(ig_to_json(sl)), inferguidR_invalid_input = function(e) e)
  expect_match(cond$status$reason, "JSON Schema")
})

test_that("schema validation degrades to a recorded skip without jsonvalidate", {
  res <- .ig_validate_schema("{}", "power_spec")
  if (requireNamespace("jsonvalidate", quietly = TRUE)) {
    expect_true(res$checked)
    expect_false(res$valid)
  } else {
    expect_false(res$checked)
    expect_match(res$note, "not installed")
  }
})
