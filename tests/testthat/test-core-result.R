# These tests describe the foundation behaviour for methods without an
# engine; they run against the bare registry (see th_registry_only()).
th_registry_only()

test_that("every method without an engine returns a typed unsupported result and no number", {
  expect_setequal(.ig_list_method_ids(), c(th_inventory_ids, th_extension_ids))
  for (id in .ig_list_method_ids()) {
    for (mode in .ig_get_metadata(id)$modes) {
      spec <- th_minimal_spec(id, mode)
      res <- calculate_power(spec)
      expect_s3_class(res, "power_result")
      expect_identical(res$status$code, "unsupported")
      expect_identical(res$status$reason, sprintf("engine for %s not implemented yet", id))
      # No finite number anywhere outside the embedded input specification.
      expect_identical(th_finite_numbers(res, exclude = "spec"), character(), info = paste(id, mode))
      expect_identical(res$spec_hash, spec$spec_hash)
    }
  }
})

test_that("solve_design answers a priori specs and refuses other modes with a typed result", {
  ap <- th_minimal_spec("t_means_one_sample", "a_priori")
  res <- solve_design(ap)
  expect_identical(res$status$code, "unsupported")
  expect_identical(res$status$reason, "engine for t_means_one_sample not implemented yet")
  ph <- solve_design(th_minimal_spec("t_means_one_sample", "post_hoc"))
  expect_identical(ph$status$code, "invalid_input")
  expect_match(ph$status$reason, "calculate_power")
  expect_identical(th_finite_numbers(ph, exclude = "spec"), character())
})

test_that("printing a result without an ok status shows no number", {
  res <- calculate_power(th_minimal_spec("t_means_two_independent", "post_hoc"))
  out <- capture.output(print(res))
  expect_true(any(grepl("no numerical result", out)))
  body <- out[!grepl("spec hash|^<power_result>|status:", out)]
  expect_false(any(grepl("[0-9]\\.[0-9]", body)))
})

test_that("non-specifications and tampered specifications are typed invalid_input", {
  expect_identical(th_condition_code(calculate_power(list(method_id = "t_generic"))), "invalid_input")
  spec <- th_minimal_spec("t_means_one_sample", "post_hoc")
  forged <- unclass(spec)
  forged$inputs$alpha <- 0.5
  class(forged) <- class(spec)
  expect_identical(th_condition_code(calculate_power(forged)), "invalid_input")
})

test_that("a registry-only spec is refused once a definition exists, and vice versa", {
  spec <- th_minimal_spec("t_means_two_independent", "post_hoc")
  th_local_definition()
  cond <- tryCatch(calculate_power(spec), inferguidR_invalid_input = function(e) e)
  expect_match(cond$status$reason, "re-create it")
  dspec <- analysis_spec("t_means_two_independent", "post_hoc",
                         inputs = list(alpha = 0.05, n_total = 40, d = 0.5), options = list(tails = "two"))
  res <- calculate_power(dspec)
  expect_identical(res$status$code, "unsupported")
  expect_identical(res$status$reason, "solver for mode 'post_hoc' not implemented yet")
  expect_identical(th_finite_numbers(res, exclude = "spec"), character())
  .ig_unregister_method("t_means_two_independent")
  expect_identical(th_condition_code(calculate_power(dspec)), "invalid_input")
})

test_that("solver failures and malformed solver output become numerical_failure", {
  th_local_definition()
  spec <- analysis_spec("t_means_two_independent", "post_hoc",
                        inputs = list(alpha = 0.05, n_total = 40, d = 0.5), options = list(tails = "two"))
  th_local_solvers(list(post_hoc = function(method, spec, entry) stop("bracket not found")))
  res <- calculate_power(spec)
  expect_identical(res$status$code, "numerical_failure")
  expect_match(res$status$reason, "bracket not found")
  expect_identical(th_finite_numbers(res, exclude = "spec"), character())

  th_local_solvers(list(post_hoc = function(method, spec, entry) 0.8))
  res <- calculate_power(spec)
  expect_identical(res$status$code, "numerical_failure")
  expect_identical(th_finite_numbers(res, exclude = "spec"), character())

  other <- update(spec, inputs = list(d = 0.6))
  th_local_solvers(list(post_hoc = function(method, spec, entry) {
    .ig_unsupported_result(other, "wrong spec", entry)
  }))
  expect_identical(calculate_power(spec)$status$code, "numerical_failure")

  # A solver whose result breaks the no-number guard is caught as a failure.
  th_local_solvers(list(post_hoc = function(method, spec, entry) {
    .ig_power_result(spec, ig_status("infeasible", "test double"), solved = list(power = 0.42))
  }))
  res <- calculate_power(spec)
  expect_identical(res$status$code, "numerical_failure")
  expect_match(res$status$reason, "must not report numbers")
})

test_that("a matching solver result is passed through and warnings are recorded", {
  th_local_definition()
  spec <- analysis_spec("t_means_two_independent", "post_hoc",
                        inputs = list(alpha = 0.05, n_total = 40, d = 0.5), options = list(tails = "two"))
  th_local_solvers(list(post_hoc = function(method, spec, entry) {
    warning("test double warning")
    .ig_power_result(spec, ig_status("ok"), solved = list(power = 0.5),
                     quantities = list(df = 38L, ncp = NA_real_),
                     undefined = list(ncp = "test double: not computed"),
                     diagnostics = list(note = "test double, not an engine"),
                     provenance = list(entry = entry, solver = "test_double"))
  }))
  expect_identical(list_methods()$implementation_status[list_methods()$id == "t_means_two_independent"],
                   "partial")
  res <- calculate_power(spec)
  expect_identical(res$status$code, "ok")
  expect_identical(res$warnings, "test double warning")
  expect_identical(res$quantities$df, 38)
  expect_output(print(res), "solved:")
  expect_identical(ig_from_json(ig_to_json(res)), res)
})

test_that("the result constructor enforces truthful content", {
  spec <- th_minimal_spec("t_means_one_sample", "post_hoc")
  expect_error(.ig_power_result(spec, ig_status("unsupported", "x"), achieved = list(power = 0.9)),
               "must not report numbers")
  expect_error(.ig_power_result(spec, ig_status("numerical_failure", "x"), group_sizes = c(n = 10)),
               "must not report numbers")
  expect_error(.ig_power_result(spec, ig_status("unsupported", "x"),
                                method_outputs = list(moments = list(p1 = 0.3))),
               "must not report numbers")
  expect_error(.ig_power_result(spec, ig_status("ok")), "finite solved value for target 'power'")
  expect_error(.ig_power_result(spec, ig_status("ok"), solved = list(power = 0.7),
                                quantities = list(ncp = NA_real_)),
               "needs a reason")
  expect_error(.ig_power_result(spec, ig_status("ok"), solved = list(power = "0.7")), "numeric")
  expect_error(.ig_power_result(spec, ig_status("ok"), solved = list(power = 0.7),
                                group_sizes = c(n = 10.5)), "whole numbers")
  # Non-finite markers are not counted as reported numbers.
  res <- .ig_power_result(spec, ig_status("infeasible", "target power not reachable"),
                          diagnostics = list(max_power = NA_real_))
  expect_identical(res$status$code, "infeasible")
})

test_that("results round-trip through JSON and tampering is detected", {
  res <- calculate_power(th_minimal_spec("z_proportions_two_independent", "sensitivity"))
  expect_identical(ig_from_json(ig_to_json(res)), res)
  lst <- .ig_result_to_list(res)
  lst$status$reason <- "edited"
  cond <- tryCatch(ig_from_json(ig_to_json(lst), validate = FALSE),
                   inferguidR_invalid_input = function(e) e)
  expect_match(cond$status$reason, "result_hash does not match")
  lst <- .ig_result_to_list(res)
  lst$solved <- list(effect_size = 0.3)
  expect_identical(th_condition_code(ig_from_json(ig_to_json(lst), validate = FALSE)), "invalid_input")
})

test_that("results are immutable and hash independently of the runtime", {
  res <- calculate_power(th_minimal_spec("t_generic", "post_hoc"))
  expect_identical(th_condition_code(res$status <- ig_status("ok")), "invalid_input")
  lst <- .ig_result_to_list(res)
  lst$runtime$r_version <- "0.0.1"
  expect_identical(ig_from_json(ig_to_json(lst))$result_hash, res$result_hash)
  expect_match(res$r_call, "^inferguidR::calculate_power\\(inferguidR::analysis_spec\\(")
})
