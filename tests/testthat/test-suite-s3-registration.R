test_that("planning results register the shared infercalcR result methods", {
  for (generic in c("as_inference", "result_data", "result_metadata")) {
    method <- utils::getS3method(generic, "power_result", optional = TRUE,
      envir = asNamespace("infercalcR"))
    expect_true(is.function(method), info = generic)
  }
  expect_false(any(c("as_inference.power_result", "result_data.power_result") %in%
    getNamespaceExports("inferguidR")))
  expect_true("planning_methods" %in% getNamespaceExports("inferguidR"))
})
