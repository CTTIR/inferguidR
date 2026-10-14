test_that("domains check type, length and range", {
  d <- ig_domain_real(0, 1, lower_inclusive = FALSE, upper_inclusive = FALSE)
  expect_null(.ig_domain_check(d, 0.5, "alpha"))
  expect_match(.ig_domain_check(d, 0, "alpha"), "must lie in \\(0, 1\\)")
  expect_match(.ig_domain_check(d, 1, "alpha"), "must lie in")
  expect_match(.ig_domain_check(d, NA_real_, "alpha"), "finite")
  expect_match(.ig_domain_check(d, Inf, "alpha"), "finite")
  expect_match(.ig_domain_check(d, "0.5", "alpha"), "numeric")
  expect_match(.ig_domain_check(d, c(0.1, 0.2), "alpha"), "length 1")
  n <- ig_domain_integer(lower = 2)
  expect_null(.ig_domain_check(n, 10, "n"))
  expect_match(.ig_domain_check(n, 10.5, "n"), "whole number")
  expect_match(.ig_domain_check(n, 1, "n"), "must lie in \\[2, Inf\\)")
  v <- ig_domain_real(length = NA)
  expect_null(.ig_domain_check(v, c(1, 2, 3), "means"))
  ch <- ig_domain_choice(c("normal", "laplace"))
  expect_null(.ig_domain_check(ch, "normal", "parent"))
  expect_match(.ig_domain_check(ch, "uniform", "parent"), "one of")
  expect_null(.ig_domain_check(ig_domain_logical(), TRUE, "flag"))
  expect_error(ig_domain_real(1, 0), "lower <= upper")
  expect_error(ig_domain_choice(c("A", "b")), "lower-case")
})

test_that("parameters validate their declaration and default", {
  p <- ig_param("d", role = "effect", domain = ig_domain_real(), default = 0.5)
  expect_s3_class(p, "ig_param")
  expect_identical(p$default, 0.5)
  expect_identical(ig_param_n_total()$domain$type, "integer")
  expect_identical(ig_param("n", "design", ig_domain_integer(), default = 5)$default, 5L)
  expect_error(ig_param("D", "effect", ig_domain_real()), "must match")
  expect_error(ig_param("d", "size", ig_domain_real()), "role")
  expect_error(ig_param("d", "effect", ig_domain_real(0, 1), default = 2), "invalid default")
  expect_error(ig_param("d", "effect", ig_domain_real(), aliases = "d"), "aliases")
  expect_error(ig_param("d", "effect", ig_domain_real(), modes = "everything"), "subset")
  expect_error(ig_param("d", "effect", ig_domain_real(), applies_when = list(when = 1)),
               "unknown condition")
  expect_error(ig_params(ig_param_alpha(), ig_param_alpha()), "duplicated")
  expect_error(ig_params(ig_param("a", "effect", ig_domain_real(), aliases = "b"),
                         ig_param("b", "effect", ig_domain_real())), "collide")
})

test_that("options validate codes, defaults and conditions", {
  o <- ig_option("tails", c("one", "two"))
  expect_null(o$default)
  expect_error(ig_option("tails", c("One", "Two")), "lower-case codes")
  expect_error(ig_option("tails", c("one", "two"), default = "both"), "one of the values")
  expect_error(ig_options(o, o), "duplicated")
  expect_length(ig_options(), 0L)
})

test_that("applicability conditions follow the shared grammar", {
  cond <- list(all = list(list(option = "tails", "in" = "two"), list(mode = c("a_priori", "post_hoc"))))
  expect_true(.ig_eval_condition(cond, "a_priori", list(tails = "two")))
  expect_false(.ig_eval_condition(cond, "criterion", list(tails = "two")))
  expect_false(.ig_eval_condition(cond, "a_priori", list(tails = "one")))
  # A condition on an option that does not apply (absent) is false.
  expect_false(.ig_eval_condition(list(option = "tails", "in" = "two"), "a_priori", list()))
  expect_true(.ig_eval_condition(list(any = list(list(option = "a", "in" = "x"),
                                                 list(option = "b", "in" = "y"))),
                                 "post_hoc", list(b = "y")))
  expect_true(.ig_eval_condition(NULL, "post_hoc", list()))
})

test_that("a complete method definition is accepted", {
  def <- th_definition()
  expect_s3_class(def, "ig_method")
  expect_identical(def$targets$a_priori, "n_total")
  expect_identical(def$targets$sensitivity, "d")
  expect_identical(def$targets$compromise, c("alpha", "power"))
  expect_output(print(def), "t_means_two_independent")
})

test_that("method definitions reject incomplete or untruthful declarations", {
  base <- list(
    id = "t_means_one_sample", family = "t", label = list(en = "x"),
    params = ig_params(ig_param("d", "effect", ig_domain_real()), ig_param_alpha(),
                       ig_param_power(), ig_param_n_total(), ig_param_error_ratio()),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = function(p) stop("no model")
  )
  expect_s3_class(do.call(ig_method, base), "ig_method")
  no_power <- base
  no_power$params <- ig_params(ig_param("d", "effect", ig_domain_real()), ig_param_alpha(),
                               ig_param_n_total())
  no_power$modes <- c("a_priori", "post_hoc")
  expect_error(do.call(ig_method, no_power), "parameter 'power' must be declared")
  no_ratio <- base
  no_ratio$params <- ig_params(ig_param("d", "effect", ig_domain_real()), ig_param_alpha(),
                               ig_param_power(), ig_param_n_total())
  expect_error(do.call(ig_method, no_ratio), "error_ratio")
  ambiguous <- base
  ambiguous$params <- ig_params(ig_param("p1", "effect", ig_domain_real()),
                                ig_param("p2", "effect", ig_domain_real()), ig_param_alpha(),
                                ig_param_power(), ig_param_n_total(), ig_param_error_ratio())
  expect_error(do.call(ig_method, ambiguous), "sensitivity target is ambiguous")
  ambiguous$targets <- list(sensitivity = c("p2", "p1"))
  expect_identical(do.call(ig_method, ambiguous)$targets$sensitivity, c("p2", "p1"))
  checked <- base
  checked$status <- "independently_checked"
  expect_error(do.call(ig_method, checked), "requires an independent `checker`")
  checked$checker <- "validation/independent/t/one_sample.R"
  expect_s3_class(do.call(ig_method, checked), "ig_method")
  validated <- base
  validated$status <- "gpower_validated"
  expect_error(do.call(ig_method, validated), "fixture_ids")
  bad_mode <- base
  bad_mode$modes <- c("a_priori", "simulation")
  expect_error(do.call(ig_method, bad_mode), "`modes`")
  bad_power <- base
  bad_power$power <- 0.8
  expect_error(do.call(ig_method, bad_power), "`power` must be a function")
  clash <- base
  clash$options <- ig_options(ig_option("d", c("a", "b")))
  expect_error(do.call(ig_method, clash), "both a parameter and an option")
})

# A synthetic one-sample one-sided z test with a method-specific a priori
# search (tests only): power = Phi(d sqrt(N) - z_{1-alpha}).
th_hook_method <- function(search, modes = c("a_priori", "post_hoc")) {
  ig_method(
    id = "t_means_one_sample", family = "t",
    label = list(en = "Synthetic z test with its own a priori search (tests only)"),
    params = ig_params(ig_param("d", role = "effect", domain = ig_domain_real()),
                       ig_param_alpha(), ig_param_power(), ig_param_n_total(lower = 2)),
    modes = modes,
    power = function(p) {
      z <- stats::qnorm(p$alpha, lower.tail = FALSE)
      list(power = stats::pnorm(p$d * sqrt(p$n_total) - z), alpha_attained = p$alpha,
           n_total = p$n_total)
    },
    engine = list(package = "stats", fun = "pnorm (synthetic test method)", version = NA_character_),
    a_priori_search = search
  )
}

th_hook_power <- function(p, N) stats::pnorm(p$d * sqrt(N) - stats::qnorm(p$alpha, lower.tail = FALSE))

test_that("ig_method() checks the a priori search hook", {
  scan <- function(p, n_lo, n_max) list(status = "ok", n = n_lo, power = 1)
  expect_null(th_hook_method(NULL)$a_priori_search)
  expect_identical(th_hook_method(scan)$a_priori_search, scan)
  expect_error(th_hook_method("scan"), "a_priori_search")
  expect_error(th_hook_method(function(p) NULL), "function\\(p, n_lo, n_max\\)")
  expect_error(th_hook_method(scan, modes = "post_hoc"), "a priori mode is not offered")
})

test_that("the a priori solver uses a method's own search and checks its answer", {
  spec_ap <- function() {
    analysis_spec("t_means_one_sample", "a_priori", inputs = list(alpha = 0.05, power = 0.8, d = 0.3))
  }
  calls <- 0L
  scan <- function(p, n_lo, n_max) {
    calls <<- calls + 1L
    for (N in n_lo:1000) {
      if (th_hook_power(p, N) >= p$power) {
        return(list(status = "ok", n = N, power = th_hook_power(p, N),
                    diagnostics = list(search = "linear scan (test hook)", n_lower = n_lo)))
      }
    }
    list(status = "infeasible", reason = "not reached up to N = 1000 (test hook)",
         diagnostics = list(search = "linear scan (test hook)"))
  }
  th_local_definition(th_hook_method(scan))
  res <- solve_design(spec_ap())
  expect_identical(calls, 1L)
  expect_identical(res$status$code, "ok")
  n_ref <- ceiling(((stats::qnorm(0.95) + stats::qnorm(0.8)) / 0.3)^2)
  expect_identical(res$solved$n_total, as.integer(n_ref))
  expect_identical(res$diagnostics$solver, "a_priori_method_search")
  expect_identical(res$diagnostics$search, "linear scan (test hook)")
  expect_identical(res$diagnostics$n_lower, 2)
  expect_identical(res$provenance$solver, "a_priori_method_search")
  expect_identical(res$achieved$power, th_hook_power(list(d = 0.3, alpha = 0.05), n_ref))

  # A typed answer without a design is passed on, without numbers.
  th_local_definition(th_hook_method(function(p, n_lo, n_max) {
    list(status = "infeasible", reason = "test hook: no design", diagnostics = list(search = "none"))
  }))
  res <- solve_design(spec_ap())
  expect_identical(res$status$code, "infeasible")
  expect_identical(res$status$reason, "test hook: no design")
  expect_identical(res$diagnostics$solver, "a_priori_method_search")
  expect_identical(th_finite_numbers(res, exclude = c("spec", "diagnostics")), character())

  # Malformed answers, and designs that do not reach the target when
  # re-evaluated, are numerical failures.
  bad_answers <- list(
    list(status = "fine"),
    list(status = "infeasible"),
    list(status = "ok", n = 10.5, power = 0.9),
    list(status = "ok", n = 1, power = 0.9),
    list(status = "ok", n = 10, power = 0.9)          # power at N = 10 is below 0.8
  )
  for (ans in bad_answers) {
    th_local_definition(th_hook_method(function(p, n_lo, n_max) ans))
    res <- solve_design(spec_ap())
    expect_identical(res$status$code, "numerical_failure", info = deparse(ans))
    expect_match(res$status$reason, "a priori search", info = deparse(ans))
    expect_identical(th_finite_numbers(res, exclude = c("spec", "diagnostics")), character())
  }
})

test_that("real domains validate optional unbounded search distances", {
  d <- ig_domain_real(0, Inf, scan = c(1e-12, 1e12))
  expect_null(.ig_domain_check(d, 2e12, "ratio"))
  for (scan in list(c(0, 10), c(-1, 10), c(10, 1), c(1, Inf), 10, c(1, NA_real_))) {
    expect_error(ig_domain_real(0, Inf, scan = scan), "scan")
  }
  expect_error(ig_domain_real(0, 1, scan = c(1e-12, 1e12)), "unbounded")
})

test_that("method definitions validate sensitivity domain hooks", {
  args <- unclass(th_hook_method(NULL, modes = "sensitivity"))
  args$sensitivity_domain <- function(p) ig_domain_real(0, 1)
  expect_s3_class(do.call(ig_method, args), "ig_method")
  args$sensitivity_domain <- function() ig_domain_real(0, 1)
  expect_error(do.call(ig_method, args), "sensitivity_domain")
  args$sensitivity_domain <- 1
  expect_error(do.call(ig_method, args), "sensitivity_domain")
  args$sensitivity_domain <- function(p) ig_domain_real(0, 1)
  args$modes <- "post_hoc"
  args$targets <- NULL
  expect_error(do.call(ig_method, args), "sensitivity mode is not offered")
})
