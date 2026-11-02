# Numerical regression blocks retained from the historical guidance suite.
test_that("the generic set supremum reproduces region and Fisher sizes", {
  reg <- boschloo_region(9, 13, 0.05, "central")
  member <- logical(10 * 14)
  member[reg$rejected$x1 + 1 + reg$rejected$x2 * 10] <- TRUE
  s <- .bl_set_sup(9, 13, member)
  expect_equal(s$value, reg$size$value, tolerance = 1e-12)
  expect_lte(s$upper, 0.05)
  cmp <- boschloo_compare_fisher(9, 13, alpha = 0.05, variant = "central")
  fs <- .bl_set_sup(9, 13, cmp$tables$fisher_reject)
  expect_equal(fs$value, cmp$size$fisher$value, tolerance = 1e-12)
  expect_identical(.bl_set_sup(9, 13, logical(140))$value, 0)
  expect_identical(.bl_set_sup(9, 13, rep(TRUE, 140))$value, 1)
  expect_error(.bl_set_sup(9, 13, TRUE), "one element per table")
})
