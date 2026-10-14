test_that("unrepresentable Poisson grids fail before allocation", {
  expect_identical(.ig_glm_expect(
    list(type = "poisson", lambda = 1), function(x) 1000 * x,
    function(x) rep(1, length(x)), tilt = 1000), NA_real_)
})

test_that("finite but oversized Poisson grids fail before allocation", {
  expect_identical(.ig_glm_expect(
    list(type = "poisson", lambda = 1), function(x) 30 * x,
    function(x) rep(1, length(x)), tilt = 30), NA_real_)
})

test_that("ordinary Poisson expectations retain their analytic value", {
  value <- .ig_glm_expect(list(type = "poisson", lambda = 2),
                         function(x) 0.3 * x, function(x) rep(1, length(x)), tilt = 0.3)
  expect_equal(value, exp(2 * (exp(0.3) - 1)), tolerance = 1e-12)
})
