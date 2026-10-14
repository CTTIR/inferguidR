test_that("localized decimal entry preserves values and rejects ambiguous grouping", {
  parse <- .ig_workbench_parse_number
  expect_equal(parse("  -1,25e2 ", "de"), -125)
  expect_equal(parse("+.125", "en"), 1 / 8)
  expect_equal(parse(",125", "de"), 1 / 8)
  expect_equal(parse("1e-3", "de"), .001)
  expect_equal(parse("1234", "de"), 1234)
  for (x in c("1.234", "1,234.5", "1 234", "NaN", "Inf", "1e309", "", "--1"))
    expect_error(parse(x, "de"), class = "inferguidR_invalid_input")
  for (x in c("1,234", "1.234,5", "1 234", "NaN", "Inf", "1e309", "", "--1"))
    expect_error(parse(x, "en"), class = "inferguidR_invalid_input")
})
