test_that("portable report content preserves null values in nested lists", {
  content <- list(absent = NULL, nested = list(NULL, list(absent = NULL)), empty = list())
  expect_identical(.ig_report_plain(content), content)
  expect_identical(.ig_canonical_json(.ig_report_plain(content)),
                   '{"absent":null,"empty":[],"nested":[null,{"absent":null}]}')
})

test_that("power reports retain absent comparisons through portable exports", {
  spec <- analysis_spec("t_means_two_independent", "a_priori",
    inputs = list(alpha = .05, power = .8, d = .5), options = list(tails = "two"))
  report <- create_report_snapshot(analysis_envelope(spec, solve_design(spec)))
  expect_identical(report$comparisons, list(NULL))
  expect_null(report$results[[1]]$comparisons)
  expect_match(as.character(ig_to_json(report)), '"comparisons":[null]', fixed = TRUE)
  expect_identical(ig_hash(ig_from_json(ig_to_json(report))), ig_hash(report))
  for (format in c("json", "txt")) {
    destination <- tempfile(fileext = paste0(".", format))
    on.exit(unlink(destination), add = TRUE)
    expect_no_error(export_report(report, destination, format))
    expect_gt(file.info(destination)$size, 0)
  }
})
