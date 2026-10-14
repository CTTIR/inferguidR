# Generic mode solvers (R/num-modes.R) on synthetic method definitions with
# closed-form answers. The synthetic definitions exist only in these tests;
# they are registered under inventory ids for the duration of a test and the
# shipped definition (if any) is restored afterwards.

ref_reset_registry()

# Paths of finite numbers in a result outside the named sections.
nm_numbers <- function(res, exclude) {
  x <- unclass(res)
  .ig_finite_paths(x[setdiff(names(x), exclude)], "")
}

nm_local_method <- function(def, env = parent.frame()) {
  old <- .ig_get_method(def$id)
  .ig_register_method(def, replace = TRUE)
  withr::defer({
    if (is.null(old)) .ig_unregister_method(def$id) else .ig_register_method(old, replace = TRUE)
  }, envir = env)
  invisible(def)
}

# One-sample one-sided z test (synthetic): power = Phi(d sqrt(N) - z_{1-alpha}),
# beta computed directly as Phi(z_{1-alpha} - d sqrt(N)).
nm_z_power <- function(p) {
  z <- qnorm(p$alpha, lower.tail = FALSE)
  m <- p$d * sqrt(p$n_total)
  list(power = pnorm(m - z), beta = pnorm(z - m), alpha_attained = p$alpha, critical = z,
       ncp = m, n_total = p$n_total)
}

nm_z_method <- function(power = nm_z_power, effect_direction = "increasing", discrete = FALSE,
                        design = NULL) {
  ig_method(
    id = "t_means_one_sample", family = "t",
    label = list(en = "Synthetic one-sample z test (tests only)"),
    params = ig_params(
      ig_param("d", role = "effect", domain = ig_domain_real()),
      ig_param_alpha(), ig_param_power(), ig_param_error_ratio(), ig_param_n_total(lower = 2)
    ),
    modes = .ig_modes, power = power, design = design, discrete = discrete,
    effect_direction = effect_direction,
    engine = list(package = "stats", fun = "pnorm (synthetic test method)", version = NA_character_)
  )
}

nm_spec <- function(mode, ...) {
  analysis_spec("t_means_one_sample", mode, inputs = list(...))
}

test_that("forward warnings survive every public solver", {
  warned <- function(p) {
    out <- nm_z_power(p)
    out$warnings <- "Synthetic model caveat"
    out
  }
  nm_local_method(nm_z_method(power = warned))
  specs <- list(
    nm_spec("post_hoc", alpha = 0.05, n_total = 50, d = 0.3),
    nm_spec("a_priori", alpha = 0.05, power = 0.8, d = 0.3),
    nm_spec("sensitivity", alpha = 0.05, power = 0.8, n_total = 50),
    nm_spec("criterion", power = 0.8, n_total = 50, d = 0.3),
    nm_spec("compromise", error_ratio = 4, n_total = 50, d = 0.3)
  )
  for (spec in specs) {
    res <- calculate_power(spec)
    expect_identical(res$status$code, "ok")
    expect_identical(res$warnings, "Synthetic model caveat")
  }
})

test_that("a priori hook warnings are retained and validated", {
  def <- nm_z_method()
  def$a_priori_search <- function(p, n_lo, n_max) {
    list(status = "ok", n = 100, diagnostics = list(), warnings = "Adjusted design")
  }
  nm_local_method(def)
  spec <- nm_spec("a_priori", alpha = 0.05, power = 0.8, d = 0.3)
  res <- calculate_power(spec)
  expect_identical(res$status$code, "ok")
  expect_identical(res$warnings, "Adjusted design")
  expect_gte(res$achieved$power, 0.8)

  for (warning in list(1, NA_character_)) {
    def$a_priori_search <- local({
      value <- warning
      function(p, n_lo, n_max) list(status = "ok", n = 100, warnings = value)
    })
    nm_local_method(def)
    res <- calculate_power(spec)
    expect_identical(res$status$code, "numerical_failure")
    expect_match(res$status$reason, "warnings")
  }
})

test_that("malformed forward warnings produce a typed failure", {
  for (warning in list(1, NA_character_)) {
    model <- local({
      value <- warning
      function(p) { out <- nm_z_power(p); out$warnings <- value; out }
    })
    nm_local_method(nm_z_method(power = model))
    res <- calculate_power(nm_spec("post_hoc", alpha = 0.05, n_total = 50, d = 0.3))
    expect_identical(res$status$code, "numerical_failure")
    expect_match(res$status$reason, "warnings")
  }
})

test_that("all five modes reproduce the closed-form z-test answers", {
  nm_local_method(nm_z_method())
  za <- qnorm(0.95)
  zb <- qnorm(0.8)

  ap <- solve_design(nm_spec("a_priori", alpha = 0.05, power = 0.8, d = 0.3))
  expect_identical(ap$status$code, "ok")
  expect_identical(ap$solved$n_total, as.integer(ceiling(((za + zb) / 0.3)^2)))
  expect_identical(ap$diagnostics$allocation_rule, "single sample: N is the design")
  expect_gte(ap$achieved$power, 0.8)
  expect_lt(nm_z_power(list(alpha = 0.05, d = 0.3, n_total = ap$solved$n_total - 1))$power, 0.8)

  ph <- calculate_power(nm_spec("post_hoc", alpha = 0.05, n_total = 50, d = 0.3))
  expect_equal(ph$solved$power, pnorm(0.3 * sqrt(50) - za), tolerance = 1e-15)
  expect_equal(ph$achieved$beta, pnorm(za - 0.3 * sqrt(50)), tolerance = 1e-14)

  se <- calculate_power(nm_spec("sensitivity", alpha = 0.05, power = 0.8, n_total = 50))
  expect_equal(se$solved$d, (za + zb) / sqrt(50), tolerance = 1e-10)
  expect_lt(abs(se$diagnostics$residual), 1e-9)
  expect_identical(se$diagnostics$root_kind, "crossing")
  expect_length(se$diagnostics$roots_found, 1L)

  cr <- calculate_power(nm_spec("criterion", power = 0.8, n_total = 50, d = 0.3))
  expect_equal(cr$solved$alpha, pnorm(zb - 0.3 * sqrt(50)), tolerance = 1e-10)
  expect_equal(cr$quantities$critical, qnorm(cr$solved$alpha, lower.tail = FALSE), tolerance = 1e-12)
  expect_lt(abs(cr$diagnostics$residual), 1e-9)

  co <- calculate_power(nm_spec("compromise", error_ratio = 4, n_total = 50, d = 0.3))
  expect_identical(co$status$code, "ok")
  expect_equal(co$achieved$beta, 4 * co$solved$alpha, tolerance = 1e-10)
  expect_equal(co$solved$power, 1 - co$achieved$beta, tolerance = 1e-14)
  expect_equal(co$achieved$error_ratio, 4, tolerance = 1e-10)
  # An independent root of Phi(z_{1-a} - d sqrt(N)) = 4 a.
  a_ref <- uniroot(function(a) pnorm(qnorm(a, lower.tail = FALSE) - 0.3 * sqrt(50)) - 4 * a,
                   c(1e-8, 0.2), tol = 1e-14)$root
  expect_equal(co$solved$alpha, a_ref, tolerance = 1e-9)
})

test_that("the effect direction decides between sensitivity roots", {
  # Two roots (power symmetric in d).
  sym <- function(p) {
    z <- qnorm(p$alpha, lower.tail = FALSE)
    m <- abs(p$d) * sqrt(p$n_total)
    list(power = pnorm(m - z), beta = pnorm(z - m), alpha_attained = p$alpha)
  }
  nm_local_method(nm_z_method(power = sym))
  res <- calculate_power(nm_spec("sensitivity", alpha = 0.05, power = 0.8, n_total = 50))
  expect_identical(res$status$code, "ok")
  expect_gt(res$solved$d, 0)
  expect_equal(res$diagnostics$other_branch_roots, -res$solved$d, tolerance = 1e-9)
  expect_match(res$diagnostics$branch_rule, "smallest root at which power is increasing")

  nm_local_method(nm_z_method(power = sym, effect_direction = "unknown"))
  res <- calculate_power(nm_spec("sensitivity", alpha = 0.05, power = 0.8, n_total = 50))
  expect_identical(res$status$code, "numerical_failure")
  expect_match(res$status$reason, "2 roots found")
  expect_identical(nm_numbers(res, exclude = c("spec", "diagnostics")), character())

  # Power decreasing in the effect.
  dec <- function(p) nm_z_power(utils::modifyList(p, list(d = -p$d)))
  nm_local_method(nm_z_method(power = dec, effect_direction = "decreasing"))
  res <- calculate_power(nm_spec("sensitivity", alpha = 0.05, power = 0.8, n_total = 50))
  expect_equal(res$solved$d, -(qnorm(0.95) + qnorm(0.8)) / sqrt(50), tolerance = 1e-10)
  expect_match(res$diagnostics$branch_rule, "largest root at which power is decreasing")
})

test_that("unreachable targets are infeasible with an explanation and no number", {
  capped <- function(p) {
    fw <- nm_z_power(p)
    list(power = 0.5 * fw$power, alpha_attained = p$alpha)
  }
  nm_local_method(nm_z_method(power = capped))
  checks <- list(
    calculate_power(nm_spec("sensitivity", alpha = 0.05, power = 0.8, n_total = 50)),
    calculate_power(nm_spec("criterion", power = 0.8, n_total = 50, d = 0.3)),
    solve_design(nm_spec("a_priori", alpha = 0.05, power = 0.8, d = 0.3))
  )
  for (res in checks) {
    expect_identical(res$status$code, "infeasible")
    expect_match(res$status$reason, "not reached")
    expect_identical(nm_numbers(res, exclude = c("spec", "diagnostics")), character())
  }
  # A target below the power at the smallest effect is not reachable either.
  floor_power <- function(p) list(power = 0.9 + 0.05 * pnorm(p$d), alpha_attained = p$alpha)
  nm_local_method(nm_z_method(power = floor_power))
  res <- calculate_power(nm_spec("sensitivity", alpha = 0.05, power = 0.8, n_total = 50))
  expect_identical(res$status$code, "infeasible")
  expect_match(res$status$reason, "exceeds the target")
  res <- calculate_power(nm_spec("criterion", power = 0.8, n_total = 50, d = 0.3))
  expect_identical(res$status$code, "infeasible")
  expect_match(res$status$reason, "exceeds the target 0.8 even for alpha = 1e-300")
})

test_that("a jump at the solution of a continuous method is a numerical failure", {
  jumpy <- function(p) list(power = if (p$d < 0.3) 0.1 else 0.95, alpha_attained = p$alpha)
  nm_local_method(nm_z_method(power = jumpy))
  res <- calculate_power(nm_spec("sensitivity", alpha = 0.05, power = 0.8, n_total = 50))
  expect_identical(res$status$code, "numerical_failure")
  expect_match(res$status$reason, "without a continuous trend")
  expect_identical(nm_numbers(res, exclude = c("spec", "diagnostics")), character())
})

test_that("alpha is solved on the logit scale, accurately at both ends of (0, 1)", {
  nm_local_method(nm_z_method())
  m <- 0.05 * sqrt(50)
  # Criterion with a target power close to one: alpha close to one.
  # Closed form: Phi(m - z_{1-a}) = target  =>  a = Phi(qnorm(target) - m).
  for (target in c(0.9999999, 1e-6)) {
    cr <- calculate_power(nm_spec("criterion", power = target, n_total = 50, d = 0.05))
    expect_identical(cr$status$code, "ok", info = target)
    a_ref <- pnorm(qnorm(target) - m)
    expect_equal(cr$solved$alpha, a_ref, tolerance = 1e-10, info = target)
    expect_equal(1 - cr$solved$alpha, 1 - a_ref, tolerance = 1e-6, info = target)
    expect_identical(cr$diagnostics$scale, "logit(alpha)")
    expect_identical(cr$diagnostics$search_range, c(1e-300, 1 - 1e-12))
  }
  # Compromise with extreme ratios q = beta/alpha. Reference: root of
  # Phi(z_{1-a} - m) - q a on the logit scale, written out independently.
  m <- 0.3 * sqrt(50)
  for (q in c(1e-6, 1e-3, 1, 1e3, 1e6)) {
    co <- calculate_power(nm_spec("compromise", error_ratio = q, n_total = 50, d = 0.3))
    expect_identical(co$status$code, "ok", info = q)
    g <- function(u) {
      a <- plogis(u)
      log(pnorm(qnorm(a, lower.tail = FALSE) - m)) - log(q * a)
    }
    a_ref <- plogis(uniroot(g, c(-700, 30), tol = 1e-14)$root)
    expect_equal(co$solved$alpha, a_ref, tolerance = 1e-9, info = q)
    expect_equal(co$achieved$beta, q * co$solved$alpha, tolerance = 1e-9, info = q)
    expect_lt(abs(co$diagnostics$residual), 1e-9)
  }
  # A ratio so small that alpha would have to exceed 1 - 1e-12.
  tiny <- calculate_power(nm_spec("compromise", error_ratio = 1e-300, n_total = 50, d = 0.3))
  expect_identical(tiny$status$code, "infeasible")
  expect_match(tiny$status$reason, "search range \\[1e-300, 1 - 1e-12\\]")
  expect_identical(nm_numbers(tiny, exclude = c("spec", "diagnostics")), character())
})

test_that("a target equal to the minimum of the power curve is touched, not solved", {
  # Two-sided z test: power is symmetric in d with its minimum alpha at d = 0.
  two_sided <- function(p) {
    z <- qnorm(p$alpha / 2, lower.tail = FALSE)
    m <- p$d * sqrt(p$n_total)
    pw <- pnorm(m - z) + pnorm(-m - z)
    list(power = pw, beta = pnorm(z - m) - pnorm(-z - m), alpha_attained = p$alpha)
  }
  nm_local_method(nm_z_method(power = two_sided))
  res <- calculate_power(nm_spec("sensitivity", alpha = 0.05, power = 0.05, n_total = 50))
  expect_identical(res$status$code, "infeasible")
  expect_match(res$status$reason, "equals the smallest power on the scan .* touched there but not crossed")
  expect_identical(nm_numbers(res, exclude = c("spec", "diagnostics")), character())
  # Slightly above the minimum there is a genuine solution.
  res <- calculate_power(nm_spec("sensitivity", alpha = 0.05, power = 0.06, n_total = 50))
  expect_identical(res$status$code, "ok")
  expect_gt(res$solved$d, 0)
  expect_equal(res$achieved$power, 0.06, tolerance = 1e-10)
  # A target measurably above the minimum but within the residual tolerance
  # (1e-11 above alpha) is crossed, so it is not infeasible; its solution
  # cannot be separated from the minimum, which is a numerical limit. The
  # reason shows the target with all its digits.
  res <- calculate_power(nm_spec("sensitivity", alpha = 0.05, power = 0.05 + 1e-11, n_total = 50))
  expect_identical(res$status$code, "numerical_failure")
  expect_match(res$status$reason, "0.05000000001 exceeds the smallest power on the scan", fixed = TRUE)
  expect_match(res$status$reason, "cannot be separated from it numerically", fixed = TRUE)
  expect_identical(nm_numbers(res, exclude = c("spec", "diagnostics")), character())
})

test_that("forward-model failures become typed results without numbers", {
  pending <- function(p) .ig_abort("reference_pending", "branch semantics not captured yet")
  nm_local_method(nm_z_method(power = pending))
  for (mode in c("post_hoc", "sensitivity", "criterion")) {
    spec <- switch(mode,
                   post_hoc = nm_spec(mode, alpha = 0.05, n_total = 50, d = 0.3),
                   sensitivity = nm_spec(mode, alpha = 0.05, power = 0.8, n_total = 50),
                   criterion = nm_spec(mode, power = 0.8, n_total = 50, d = 0.3))
    res <- calculate_power(spec)
    expect_identical(res$status$code, "reference_pending", info = mode)
    expect_identical(res$status$reason, "branch semantics not captured yet")
    expect_identical(nm_numbers(res, exclude = "spec"), character())
  }

  broken <- function(p) stop("division by zero in the test model")
  nm_local_method(nm_z_method(power = broken))
  res <- calculate_power(nm_spec("post_hoc", alpha = 0.05, n_total = 50, d = 0.3))
  expect_identical(res$status$code, "numerical_failure")
  expect_match(res$status$reason, "division by zero")

  outside <- function(p) list(power = 1.5, alpha_attained = p$alpha)
  nm_local_method(nm_z_method(power = outside))
  res <- calculate_power(nm_spec("post_hoc", alpha = 0.05, n_total = 50, d = 0.3))
  expect_identical(res$status$code, "numerical_failure")
  expect_match(res$status$reason, "outside \\[0, 1\\]")

  silent_na <- function(p) list(power = 0.5, alpha_attained = p$alpha, ncp = NA_real_)
  nm_local_method(nm_z_method(power = silent_na))
  res <- calculate_power(nm_spec("post_hoc", alpha = 0.05, n_total = 50, d = 0.3))
  expect_identical(res$status$code, "numerical_failure")
  expect_match(res$status$reason, "undefined without a reason")

  na_power <- function(p) list(power = NA_real_, alpha_attained = p$alpha,
                               notes = list(power = "the test model has no power here"))
  nm_local_method(nm_z_method(power = na_power))
  res <- calculate_power(nm_spec("post_hoc", alpha = 0.05, n_total = 50, d = 0.3))
  expect_identical(res$status$code, "numerical_failure")
  expect_match(res$status$reason, "the test model has no power here")

  # An undefined secondary quantity with a reason is reported as such.
  na_ncp <- function(p) c(nm_z_power(p)[c("power", "beta", "alpha_attained")],
                          list(ncp = NA_real_, notes = list(ncp = "not defined for this test model")))
  nm_local_method(nm_z_method(power = na_ncp))
  res <- calculate_power(nm_spec("post_hoc", alpha = 0.05, n_total = 50, d = 0.3))
  expect_identical(res$status$code, "ok")
  expect_true(is.na(res$quantities$ncp))
  expect_identical(res$undefined$ncp, "not defined for this test model")
})

test_that("a design rule that returns non-integer group sizes is rejected", {
  bad_design <- function(p, n_total, continuous = FALSE) c(n1 = n_total / 3, n2 = 2 * n_total / 3)
  nm_local_method(nm_z_method(design = bad_design))
  res <- solve_design(nm_spec("a_priori", alpha = 0.05, power = 0.8, d = 0.3))
  expect_identical(res$status$code, "numerical_failure")
  expect_match(res$status$reason, "not whole numbers")
})

# Discrete synthetic test: exact one-sided binomial test of p0 = 0.5 against
# p1 > 0.5 with rejection region X >= k, k the smallest value whose size is
# at most the nominal alpha.
nm_binom_power <- function(p) {
  n <- p$n_total
  k <- 0:(n + 1)
  size <- pbinom(k - 1, n, 0.5, lower.tail = FALSE)
  k <- k[which(size <= p$alpha)[1]]
  list(power = pbinom(k - 1, n, p$p1, lower.tail = FALSE), beta = pbinom(k - 1, n, p$p1),
       alpha_attained = pbinom(k - 1, n, 0.5, lower.tail = FALSE), critical = k, n_total = n)
}

nm_binom_method <- function() {
  ig_method(
    id = "exact_proportion_one_sample_binomial", family = "exact",
    label = list(en = "Synthetic exact binomial test (tests only)"),
    params = ig_params(
      ig_param("p1", role = "effect",
               domain = ig_domain_real(0, 1, lower_inclusive = FALSE, upper_inclusive = FALSE)),
      ig_param_alpha(), ig_param_power(), ig_param_error_ratio(), ig_param_n_total(lower = 1)
    ),
    modes = .ig_modes, power = nm_binom_power, discrete = TRUE, effect_direction = "increasing",
    engine = list(package = "stats", fun = "pbinom (synthetic test method)", version = NA_character_)
  )
}

nm_bspec <- function(mode, ...) analysis_spec("exact_proportion_one_sample_binomial", mode, inputs = list(...))

test_that("discrete a priori: lowest N reaching the target over a non-monotone power sequence", {
  nm_local_method(nm_binom_method())
  res <- solve_design(nm_bspec("a_priori", alpha = 0.05, power = 0.8, p1 = 0.75))
  pw <- vapply(1:200, function(N) nm_binom_power(list(n_total = N, alpha = 0.05, p1 = 0.75))$power,
               numeric(1))
  expected <- which(pw >= 0.8)[1]
  expect_true(any(diff(pw[1:expected]) < 0)) # the sequence really is non-monotone
  expect_identical(res$status$code, "ok")
  expect_identical(res$solved$n_total, expected)
  expect_identical(res$diagnostics$search, "enumeration")
  expect_identical(res$diagnostics$predecessor_n_total, expected - 1)
  expect_lt(res$diagnostics$predecessor_power, 0.8)
  expect_lt(res$achieved$alpha_attained, 0.05)
})

test_that("discrete criterion and compromise keep nominal and attained alpha apart", {
  nm_local_method(nm_binom_method())
  cr <- calculate_power(nm_bspec("criterion", power = 0.6, n_total = 20, p1 = 0.75))
  expect_identical(cr$status$code, "ok")
  expect_identical(cr$diagnostics$root_kind, "jump")
  expect_match(cr$diagnostics$discrete_note, "step function")
  a <- cr$solved$alpha
  expect_gte(cr$achieved$power, 0.6)
  expect_lte(cr$achieved$alpha_attained, a)
  # Just below the reported nominal alpha the rule changes and the target is missed.
  below <- nm_binom_power(list(n_total = 20, alpha = a * (1 - 1e-9), p1 = 0.75))
  expect_lt(below$power, 0.6)
  # Brute force over the attainable sizes: the smallest size whose rule reaches 0.6.
  sizes <- pbinom(0:20 - 1, 20, 0.5, lower.tail = FALSE)
  pows <- pbinom(0:20 - 1, 20, 0.75, lower.tail = FALSE)
  expect_equal(a, min(sizes[pows >= 0.6]), tolerance = 1e-12)

  co <- calculate_power(nm_bspec("compromise", error_ratio = 1, n_total = 20, p1 = 0.75))
  expect_identical(co$status$code, "ok")
  expect_match(co$diagnostics$discrete_note, "nominal alpha")
  expect_lte(co$achieved$beta, co$solved$alpha * (1 + 1e-9))
  expect_lte(co$achieved$alpha_attained, co$solved$alpha)
  expect_equal(co$achieved$error_ratio, co$achieved$beta / co$achieved$alpha_attained)
})

test_that("discrete sensitivity: power is continuous in the effect for a fixed rule", {
  nm_local_method(nm_binom_method())
  res <- calculate_power(nm_bspec("sensitivity", alpha = 0.05, power = 0.8, n_total = 20))
  expect_identical(res$status$code, "ok")
  expect_identical(res$diagnostics$root_kind, "crossing")
  expect_lt(abs(res$achieved$power - 0.8), 1e-9)
  k <- nm_binom_power(list(n_total = 20, alpha = 0.05, p1 = 0.6))$critical
  p_ref <- uniroot(function(x) pbinom(k - 1, 20, x, lower.tail = FALSE) - 0.8, c(0.5, 0.999),
                   tol = 1e-14)$root
  expect_equal(res$solved$p1, p_ref, tolerance = 1e-9)
})

test_that("method sensitivity hooks are checked and exposed through calculate_power", {
  hook <- function(p, target_name) {
    list(status = "ok", value = (qnorm(1 - p$alpha) + qnorm(p$power)) / sqrt(p$n_total),
         diagnostics = list(target_reached = TRUE), warnings = "Method sensitivity search")
  }
  args <- unclass(nm_z_method()); args$sensitivity_search <- hook
  nm_local_method(do.call(ig_method, args))
  spec <- nm_spec("sensitivity", alpha = 0.05, power = 0.8, n_total = 50)
  result <- calculate_power(spec)
  expect_identical(result$status$code, "ok")
  expect_equal(result$achieved$power, 0.8, tolerance = 1e-10)
  expect_identical(result$warnings, "Method sensitivity search")
  expect_identical(result$diagnostics$solver, "sensitivity_method_search")

  for (answer in list(
    list(status = "ok", value = Inf),
    list(status = "ok", value = 0),
    list(status = "ok", value = 0, diagnostics = list(target_reached = FALSE)),
    list(status = "ok", value = 0, warnings = NA_character_),
    list(status = "unsupported"))) {
    args$sensitivity_search <- local({value <- answer; function(p, target_name) value})
    nm_local_method(do.call(ig_method, args))
    expect_identical(calculate_power(spec)$status$code, "numerical_failure")
  }
})

test_that("an unattained sensitivity target needs an explicit checked warning", {
  args <- unclass(nm_z_method(power = function(p) list(power = p$alpha, alpha_attained = p$alpha)))
  args$sensitivity_search <- function(p, target_name) {
    list(status = "ok", value = 0, diagnostics = list(target_reached = FALSE, max_power = p$alpha),
         warnings = "Target not reached in the searched domain")
  }
  nm_local_method(do.call(ig_method, args))
  r <- calculate_power(nm_spec("sensitivity", alpha = 0.05, power = 0.8, n_total = 50))
  expect_identical(r$status$code, "ok")
  expect_false(r$diagnostics$target_reached)
  expect_equal(r$achieved$power, 0.05, tolerance = 1e-12)
  expect_match(r$warnings, "Target not reached")
})

test_that("a sensitivity direction callback selects the option-defined branch", {
  sym <- function(p) nm_z_power(utils::modifyList(p, list(d = abs(p$d))))
  def <- nm_z_method(power = sym, effect_direction = function(p) {
    if (p$tails == "one") "decreasing" else "increasing"
  })
  def$options <- ig_options(ig_option("tails", values = c("one", "two"), default = "one"))
  nm_local_method(def)
  for (tails in c("one", "two")) {
    res <- calculate_power(analysis_spec("t_means_one_sample", "sensitivity",
      inputs = list(alpha = 0.05, power = 0.8, n_total = 50), options = list(tails = tails)))
    want <- (qnorm(0.95) + qnorm(0.8)) / sqrt(50) * if (tails == "one") -1 else 1
    expect_identical(res$status$code, "ok")
    expect_equal(res$solved$d, want, tolerance = 1e-9)
    expect_equal(res$diagnostics$other_branch_roots, -want, tolerance = 1e-9)
    expect_equal(res$diagnostics$failed_evaluations, 0)
  }
})

test_that("malformed sensitivity direction callbacks give typed failures", {
  for (bad in list(NA_character_, c("increasing", "decreasing"), "sideways", 1)) {
    nm_local_method(nm_z_method(effect_direction = function(p) bad))
    res <- calculate_power(nm_spec("sensitivity", alpha = 0.05, power = 0.8, n_total = 50))
    expect_identical(res$status$code, "numerical_failure")
    expect_match(res$status$reason, "effect_direction")
    expect_identical(nm_numbers(res, exclude = c("spec", "diagnostics")), character())
  }
  expect_error(nm_z_method(effect_direction = function() "increasing"), "function")
})

test_that("constraints receive declared mode and target context", {
  def <- nm_z_method()
  seen <- NULL
  def$constraints <- function(p, mode, target) {
    seen <<- list(mode = mode, target = target)
    if (mode == "post_hoc") "post hoc excluded for this fixture" else character()
  }
  nm_local_method(def)
  expect_error(nm_spec("post_hoc", d = 0.5, alpha = 0.05, n_total = 50), "post hoc excluded")
  expect_identical(seen, list(mode = "post_hoc", target = "power"))
  expect_s3_class(nm_spec("a_priori", d = 0.5, alpha = 0.05, power = 0.8), "power_spec")
  expect_identical(seen, list(mode = "a_priori", target = "n_total"))
})

test_that("sensitivity domain hooks restrict evaluations and reject invalid domains", {
  def <- nm_z_method()
  visited <- numeric()
  def$power <- function(p) { visited <<- c(visited, p$d); nm_z_power(p) }
  def$sensitivity_domain <- function(p) ig_domain_real(0.1, 0.8)
  nm_local_method(def)
  r <- calculate_power(nm_spec("sensitivity", alpha = 0.05, power = 0.8, n_total = 50))
  expect_identical(r$status$code, "ok")
  expect_true(all(visited >= 0.1 & visited <= 0.8))
  expect_equal(r$solved$d, (qnorm(0.95) + qnorm(0.8))/sqrt(50), tolerance = 1e-8)
  for (bad in list(NULL, c(0, 1), ig_domain_integer(0, 1), ig_domain_real(0, 1, length = 2))) {
    def$sensitivity_domain <- function(p) bad
    nm_local_method(def)
    r <- calculate_power(nm_spec("sensitivity", alpha = 0.05, power = 0.8, n_total = 50))
    expect_identical(r$status$code, "numerical_failure")
  }
  def$params$d$domain <- ig_domain_real(0, 1, lower_inclusive = FALSE)
  for (bad in list(ig_domain_real(-1, 1), ig_domain_real(0, 1))) {
    def$sensitivity_domain <- function(p) bad
    nm_local_method(def)
    r <- calculate_power(nm_spec("sensitivity", alpha = 0.05, power = 0.8, n_total = 50))
    expect_identical(r$status$code, "numerical_failure")
  }
})
