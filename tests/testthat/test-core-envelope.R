# Exercise the foundation contracts independently of later method providers.
# Restore the live registry and solvers when this file finishes.
th_registry_only()

th_classification_spec <- function(...) {
  args <- th_override(list(
    input_kind = "counts_2x2",
    models = list(list(model_id = "m1", name = "Model 1", version = "1.0",
                       score_meaning = "predicted risk", score_scale = "probability",
                       positive_label = "event", direction = "higher_is_positive",
                       comparator = ">=", tie_rule = "scores equal to the cutoff alert",
                       evaluation_split = "held_out_test", missing_value_policy = "exclude")),
    operating_rule = list(score_cutoff = 0.5, selection = "manual",
                          selection_data = "predeclared_rule"),
    counts = list(tp = 40, fn = 10, fp = 20, tn = 130)
  ), list(...))
  do.call(.ig_new_classification_spec, args)
}

test_that("power envelopes pair a spec with its own result", {
  spec <- th_minimal_spec("t_means_one_sample", "post_hoc")
  res <- calculate_power(spec)
  env <- analysis_envelope(spec, res)
  expect_s3_class(env, "analysis_envelope")
  expect_identical(env$analysis_kind, "power")
  expect_identical(env$result_hash, res$result_hash)
  expect_match(env$analysis_hash, "^[0-9a-f]{64}$")
  pending <- analysis_envelope(spec)
  expect_null(pending$result)
  expect_null(pending$result_hash)
  expect_false(identical(pending$analysis_hash, env$analysis_hash))
  other <- calculate_power(th_minimal_spec("t_means_one_sample", "criterion"))
  expect_identical(th_condition_code(analysis_envelope(spec, other)), "invalid_input")
  expect_identical(th_condition_code(analysis_envelope(spec, res, kind = "classification")), "invalid_input")
  expect_identical(th_condition_code(analysis_envelope(list(a = 1))), "invalid_input")
  expect_output(print(env), "unsupported")
})

test_that("the analysis hash does not depend on provenance", {
  spec <- th_minimal_spec("t_means_one_sample", "post_hoc")
  env <- analysis_envelope(spec, calculate_power(spec))
  lst <- .ig_envelope_to_list(env)
  lst$provenance$r_version <- "9.9.9"
  back <- ig_from_json(ig_to_json(lst))
  expect_identical(back$analysis_hash, env$analysis_hash)
})

test_that("envelopes round-trip through JSON and tampering is detected", {
  spec <- th_minimal_spec("chi2_generic", "compromise")
  env <- analysis_envelope(spec, calculate_power(spec))
  expect_identical(ig_from_json(ig_to_json(env)), env)
  lst <- .ig_envelope_to_list(env)
  lst$analysis_hash <- paste(rep("0", 64), collapse = "")
  expect_identical(th_condition_code(ig_from_json(ig_to_json(lst))), "invalid_input")
})

test_that("classification specifications have a validated shape", {
  cs <- th_classification_spec()
  expect_s3_class(cs, "classification_spec")
  env <- analysis_envelope(cs)
  expect_identical(env$analysis_kind, "classification")
  expect_identical(ig_from_json(ig_to_json(env)), env)
  expect_identical(th_condition_code(th_classification_spec(input_kind = "pictures")), "invalid_input")
  expect_identical(th_condition_code(th_classification_spec(counts = list(tp = -1, fn = 1, fp = 1, tn = 1))),
                   "invalid_input")
  bad_model <- th_classification_spec()$models
  bad_model[[1]]$comparator <- "=="
  expect_identical(th_condition_code(th_classification_spec(models = bad_model)), "invalid_input")
  expect_identical(th_condition_code(th_classification_spec(
    operating_rule = list(selection = "manual", selection_data = "same_data"))), "invalid_input")
  expect_identical(th_condition_code(th_classification_spec(
    input_kind = "theoretical_binormal", counts = NULL, theory = list(auc = 1.2, sigma = 1))),
    "invalid_input")
  binormal <- th_classification_spec(input_kind = "theoretical_binormal", counts = NULL,
                                     theory = list(auc = 0.8, sigma = 1.2),
                                     prevalence = list(kind = "scenario", value = 0.1, source = "assumed"))
  expect_identical(binormal$theory$auc, 0.8)
  expect_identical(th_condition_code(th_classification_spec(
    prevalence = list(kind = "scenario", value = NULL, source = NULL))), "invalid_input")
  expect_identical(th_condition_code(th_classification_spec(
    input_kind = "scored_data", counts = NULL)), "invalid_input")
})

test_that("observed-test specifications are carried with tamper detection", {
  os <- .ig_new_observed_test_spec("boschloo_two_sided_minlike",
                                   inputs = list(x1 = 3, n1 = 12, x2 = 9, n2 = 12))
  env <- analysis_envelope(os)
  expect_identical(env$analysis_kind, "observed_test")
  expect_identical(ig_from_json(ig_to_json(env)), env)
  forged <- unclass(os)
  forged$inputs$x1 <- 4
  class(forged) <- class(os)
  expect_identical(th_condition_code(analysis_envelope(forged)), "invalid_input")
  expect_identical(th_condition_code(.ig_new_observed_test_spec("Boschloo")), "invalid_input")
})
