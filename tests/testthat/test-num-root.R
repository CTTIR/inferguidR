# Bracketing, root finding and integer design search (R/num-root.R) on
# synthetic functions with known answers.

ref_reset_registry()

test_that("a smooth sign change is refined to a crossing with a small residual", {
  f <- function(x) x^3 - 2
  r <- .ig_refine_root(f, 0, 3, f(0), f(3))
  expect_identical(r$kind, "crossing")
  expect_equal(r$x, 2^(1 / 3), tolerance = 1e-12)
  expect_lt(abs(r$residual), 1e-9)
  expect_identical(r$bracket, c(0, 3))
  expect_identical(r$slope, 1)
  expect_gt(r$iterations, 0)
})

test_that("a discontinuity is classified as a jump, an unconverged root as unresolved", {
  step <- function(x) if (x < 1) -0.4 else 0.6
  r <- .ig_refine_root(step, 0, 2, -0.4, 0.6)
  expect_identical(r$kind, "jump")
  expect_equal(r$x, 1, tolerance = 1e-10)
  # A function whose value at the returned root does not match what the
  # root finder saw (and has no sign change around it) is not accepted as a
  # root. The function answers differently outside uniroot() to reach this
  # branch deterministically.
  inside_uniroot <- function() {
    any(vapply(sys.calls(), function(cl) {
      fn <- cl[[1]]
      identical(fn, quote(stats::uniroot)) || identical(fn, quote(uniroot))
    }, logical(1)))
  }
  odd <- function(x) if (inside_uniroot()) x - 1 else 0.5
  r <- .ig_refine_root(odd, 0, 2, -1, 1)
  expect_identical(r$kind, "unresolved")
  expect_identical(r$residual, 0.5)
})

test_that("a continuous but very steep crossing is accepted and told apart from a jump", {
  # Slope 1e12: the x resolution of the root finder (1e-13) leaves a
  # residual near 0.1, far above the tolerance, although f is continuous.
  steep <- function(x) 1e12 * (x - 1 / 3)
  r <- .ig_refine_root(steep, 0, 1, steep(0), steep(1))
  expect_identical(r$kind, "steep_crossing")
  expect_equal(r$x, 1 / 3, tolerance = 1e-11)
  expect_gt(abs(r$residual), .ig_solver_settings$residual_tol)
  # A jump of the same height is still a jump, also when it is irregular
  # (rounding noise around a level, no trend).
  step <- function(x) if (x < 1 / 3) -0.2 else 0.2
  expect_identical(.ig_refine_root(step, 0, 1, -0.2, 0.2)$kind, "jump")
  noisy <- function(x) {
    if (abs(x - 1 / 3) < 1e-9) return(if (sin(1e15 * x) > 0) 0.1 else -0.1)
    if (x < 1 / 3) -0.1 else 0.1
  }
  expect_identical(.ig_refine_root(noisy, 0, 1, -0.1, 0.1)$kind, "jump")
  # Grid scans use the same classification (a root finder that lands on the
  # exact double may also meet the residual tolerance: a plain crossing).
  scan <- .ig_scan_roots(steep, seq(0, 1, by = 0.25))
  expect_true(scan$roots[[1]]$kind %in% c("crossing", "steep_crossing"))
  expect_equal(scan$roots[[1]]$x, 1 / 3, tolerance = 1e-11)
  scan <- .ig_scan_roots(step, seq(0, 1, by = 0.25))
  expect_identical(scan$roots[[1]]$kind, "jump")
})

test_that("a grid scan reports every root, its slope and flat regions", {
  scan <- .ig_scan_roots(sin, seq(0.5, 10, by = 0.25))
  x <- vapply(scan$roots, `[[`, numeric(1), "x")
  expect_equal(x, pi * 1:3, tolerance = 1e-10)
  expect_identical(vapply(scan$roots, function(r) r$slope, numeric(1)), c(-1, 1, -1))
  expect_true(all(vapply(scan$roots, `[[`, character(1), "kind") == "crossing"))
  expect_identical(scan$failed, 0L)

  flat <- function(x) if (abs(x) <= 1) 0 else x
  scan <- .ig_scan_roots(flat, seq(-3, 3, by = 0.5))
  expect_length(scan$roots, 1L)
  expect_identical(scan$roots[[1]]$kind, "grid_zero")
  expect_identical(scan$flat, list(c(-1, 1)))

  # Non-finite values are counted and never treated as roots.
  hole <- function(x) if (x > 0 && x < 1) NaN else x - 0.5
  scan <- .ig_scan_roots(hole, seq(-1, 2, by = 0.25))
  expect_length(scan$roots, 0L)
  expect_identical(scan$failed, 3L)
})

test_that("monotone bracketing expands in the right direction or explains why it cannot", {
  inc <- function(x) x - 1000
  br <- .ig_bracket_monotone(inc, 0, -Inf, Inf, decreasing = FALSE)
  expect_true(br$found)
  expect_true(br$a <= 1000 && br$b >= 1000)
  expect_true(br$fa <= 0 && br$fb >= 0)
  dec <- function(x) 5 - x
  br <- .ig_bracket_monotone(dec, 100, -Inf, Inf, decreasing = TRUE)
  expect_true(br$found && br$a <= 5 && br$b >= 5)
  br <- .ig_bracket_monotone(function(x) x - 3, 3, -Inf, Inf, decreasing = FALSE)
  expect_true(br$found && br$exact)
  capped <- .ig_bracket_monotone(function(x) -1 / (1 + x^2), 0, -10, 10, decreasing = FALSE)
  expect_false(capped$found)
  expect_match(capped$reason, "no sign change before the end")
  expect_lt(capped$value_at_limit, 0)
  bad <- .ig_bracket_monotone(function(x) NaN, 0, -1, 1, decreasing = FALSE)
  expect_false(bad$found)
  expect_match(bad$reason, "not finite")
})

test_that("scan grids respect the domain bounds", {
  g <- .ig_domain_grid(ig_domain_real())
  expect_identical(min(g), -1e4)
  expect_identical(max(g), 1e4)
  expect_true(0 %in% g)
  g <- .ig_domain_grid(ig_domain_real(0, Inf, lower_inclusive = FALSE))
  expect_true(all(g > 0))
  g <- .ig_domain_grid(ig_domain_real(0, 1, lower_inclusive = FALSE, upper_inclusive = FALSE))
  expect_true(all(g > 0 & g < 1))
  expect_false(is.unsorted(g))
  g <- .ig_domain_grid(ig_domain_real(-2, 3))
  expect_identical(range(g), c(-2, 3))
})

# Integer design search ----------------------------------------------------------

even <- function(N) N %% 2 == 0
z_power <- function(N, effect) pnorm(effect * sqrt(N / 4) - qnorm(0.975))

test_that("monotone power: the smallest admissible N, proven by enumeration in short ranges", {
  res <- .ig_integer_search(function(N) z_power(N, 0.5), even, n_lo = 4, n_max = 1e6,
                            target = 0.8, discrete = FALSE)
  brute <- seq(4, 400, by = 2)
  expected <- brute[which(z_power(brute, 0.5) >= 0.8)[1]]
  expect_identical(res$status, "ok")
  expect_identical(res$n, expected)
  expect_identical(res$diagnostics$predecessor_n_total, expected - 2)
  expect_lt(res$diagnostics$predecessor_power, 0.8)
  expect_match(res$diagnostics$minimality, "^proven")
  expect_identical(res$diagnostics$enumerated_range, c(4, expected))
  expect_false(res$diagnostics$nonmonotone_detected)
})

test_that("monotone power with a large solution: bisection, predecessor and window checks", {
  res <- .ig_integer_search(function(N) z_power(N, 0.02), even, n_lo = 4, n_max = 1e9,
                            target = 0.9, discrete = FALSE)
  # Closed form: smallest even N with 0.02 sqrt(N/4) >= z_0.975 + z_0.9.
  n_cont <- 4 * ((qnorm(0.975) + qnorm(0.9)) / 0.02)^2
  expected <- 2 * ceiling(n_cont / 2)
  expect_identical(res$status, "ok")
  expect_identical(res$n, expected)
  expect_identical(res$diagnostics$search, "bisection")
  expect_identical(res$diagnostics$predecessor_n_total, expected - 2)
  expect_lt(res$diagnostics$predecessor_power, 0.9)
  expect_match(res$diagnostics$minimality, "50 admissible designs below it also fail")
  expect_lt(res$diagnostics$evaluations, 200)
})

# The power sequence of the manual's binomial example (manual 4.3, DOC-034,
# engine-decisions.md 5.3) for N = 10..20, extended upwards: power is not
# monotone in N, and the lowest N reaching 0.35 is 16, not 20.
staircase <- c(0.107, 0.086, 0.275, 0.234, 0.198, 0.167, 0.352, 0.310, 0.271, 0.237, 0.411)
step_power <- function(N) if (N <= 20) staircase[N - 9] else 0.411 + (N - 20) * 0.001

test_that("discrete staircase: enumeration finds the lowest N whose power is not below the target", {
  res <- .ig_integer_search(step_power, function(N) TRUE, n_lo = 10, n_max = 1000, target = 0.35,
                            discrete = TRUE)
  expect_identical(res$status, "ok")
  expect_identical(res$n, 16)
  expect_identical(res$power, 0.352)
  expect_identical(res$diagnostics$search, "enumeration")
  expect_identical(res$diagnostics$enumerated_range, c(10, 16))
  expect_identical(res$diagnostics$predecessor_n_total, 15)
  expect_identical(res$diagnostics$predecessor_power, 0.167)
  expect_match(res$diagnostics$monotonicity, "not assumed")
  expect_identical(res$diagnostics$power_sequence$power, staircase[1:7])
})

test_that("non-monotone power is detected even when the method claims continuity", {
  res <- .ig_integer_search(step_power, function(N) TRUE, n_lo = 10, n_max = 1000, target = 0.35,
                            discrete = FALSE)
  expect_identical(res$status, "ok")
  expect_identical(res$n, 16)
  expect_true(res$diagnostics$nonmonotone_detected)
  expect_match(res$diagnostics$monotonicity, "violated")

  # Beyond the full-scan span a design below the bisection result that
  # reaches the target triggers a full enumeration.
  n_smooth <- 2 * ceiling(4 * ((qnorm(0.975) + qnorm(0.9)) / 0.06)^2 / 2)
  expect_gt(n_smooth, 4000 + 4)
  n_spike <- n_smooth - 40
  spike <- function(N) if (N == n_spike) 0.95 else z_power(N, 0.06)
  res <- .ig_integer_search(spike, even, n_lo = 4, n_max = 1e6, target = 0.9, discrete = FALSE)
  expect_identical(res$n, n_spike)
  expect_true(res$diagnostics$nonmonotone_detected)
  expect_identical(res$diagnostics$search, "enumeration")
})

test_that("unreachable targets are infeasible with the reason, flat power is identified", {
  plateau <- function(N) 0.6 - 0.5 / N
  res <- .ig_integer_search(plateau, function(N) TRUE, n_lo = 2, n_max = 5000, target = 0.8,
                            discrete = FALSE)
  expect_identical(res$status, "infeasible")
  expect_match(res$reason, "not reached by any admissible design with N <= 5000")
  expect_false(res$diagnostics$flat)
  expect_null(res$n)
  flat <- .ig_integer_search(function(N) 0.05, function(N) TRUE, n_lo = 2, n_max = 5000,
                             target = 0.8, discrete = FALSE)
  expect_identical(flat$status, "infeasible")
  expect_true(flat$diagnostics$flat)
  expect_match(flat$reason, "power does not change with N")
  disc <- .ig_integer_search(function(N) 0.05, function(N) TRUE, n_lo = 2, n_max = 50,
                             target = 0.8, discrete = TRUE)
  expect_identical(disc$status, "infeasible")
  expect_match(disc$reason, "not reached")
})

test_that("admissible-set edge cases are handled explicitly", {
  # The smallest admissible design already reaches the target.
  res <- .ig_integer_search(function(N) 0.99, even, n_lo = 3, n_max = 100, target = 0.8,
                            discrete = FALSE)
  expect_identical(res$n, 4)
  expect_match(res$diagnostics$minimality, "smallest admissible design")
  # No admissible design at all in a short range.
  none <- .ig_integer_search(function(N) 0.99, function(N) FALSE, n_lo = 3, n_max = 100,
                             target = 0.8, discrete = FALSE)
  expect_identical(none$status, "infeasible")
  expect_match(none$reason, "no admissible design")
  # A gap wider than the stepping limit is an error, not a silent skip.
  expect_error(.ig_integer_search(function(N) 0.99, function(N) FALSE, n_lo = 3, n_max = 1e6,
                                  target = 0.8, discrete = FALSE), "no admissible design within")
  # A forward model that returns no power is an error.
  expect_error(.ig_integer_search(function(N) NA_real_, function(N) TRUE, n_lo = 3, n_max = 100,
                                  target = 0.8, discrete = FALSE), "returned no power")
})

test_that("per-domain scans reach large effects without expanding other domains", {
  domain <- ig_domain_real(0, Inf, scan = c(1e-12, 1e12))
  scan <- .ig_scan_roots(function(x) log1p(x) - log1p(26065.34), .ig_domain_grid(domain))
  expect_length(scan$roots, 1)
  expect_equal(scan$roots[[1]]$x, 26065.34, tolerance = 1e-10)
  default <- .ig_scan_roots(function(x) log1p(x) - log1p(26065.34),
                            .ig_domain_grid(ig_domain_real(0, Inf)))
  expect_length(default$roots, 0)
  for (dom in list(ig_domain_real(5, Inf, scan = c(0.01, 1e6)),
                  ig_domain_real(-Inf, -5, scan = c(0.01, 1e6)))) {
    g <- .ig_domain_grid(dom)
    expect_true(all(g >= dom$lower & g <= dom$upper))
    expect_true(all(diff(g) > 0))
  }
})
