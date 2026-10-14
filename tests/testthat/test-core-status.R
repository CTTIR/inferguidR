test_that("the six status codes are defined", {
  expect_identical(ig_status_codes(), c("ok", "infeasible", "unsupported", "invalid_input",
                                        "numerical_failure", "reference_pending"))
})

test_that("ig_status validates code, reason and details", {
  s <- ig_status("unsupported", "engine for t_generic not implemented yet")
  expect_true(is_ig_status(s))
  expect_identical(s$code, "unsupported")
  expect_identical(format(s), "unsupported: engine for t_generic not implemented yet")
  expect_identical(format(ig_status("ok")), "ok")
  expect_error(ig_status("fine"), "must be one of")
  expect_error(ig_status("infeasible"), "non-empty `reason`")
  expect_error(ig_status("infeasible", "   "), "non-empty `reason`")
  expect_error(ig_status("ok", NA_character_), "single string")
  expect_error(ig_status("ok", details = list(1)), "named list")
})

test_that("status details are kept in canonical order", {
  s <- ig_status("invalid_input", "x", list(b = 1, a = 2))
  expect_identical(names(s$details), c("a", "b"))
  expect_identical(names(ig_status("ok")$details), character(0))
})

test_that("typed conditions carry the status object", {
  cond <- tryCatch(.ig_abort("unsupported", "not yet"), error = function(e) e)
  expect_s3_class(cond, "inferguidR_unsupported")
  expect_s3_class(cond, "inferguidR_error")
  expect_identical(cond$status$code, "unsupported")
  cond <- tryCatch(.ig_abort_invalid(c("a is wrong", "b is wrong"), "context"),
                   error = function(e) e)
  expect_s3_class(cond, "inferguidR_invalid_input")
  expect_identical(unlist(cond$status$details$problems), c("a is wrong", "b is wrong"))
  expect_match(conditionMessage(cond), "2 problems")
})

test_that("printing a status lists its problems", {
  s <- ig_status("invalid_input", "bad", list(problems = list("p1", "p2")))
  out <- capture.output(print(s))
  expect_match(out[1], "invalid_input: bad")
  expect_true(any(grepl("p2", out)))
})
