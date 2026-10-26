ap_bos <- function(direction = "two_sided", margins = "group_sizes_fixed") {
  advisor_replay(c(q.aim = "compare", q.stage = "planning", q.outcome_type = "binary",
    q.unit = "individuals", q.structure = "independent", q.independence_check = "no",
    q.n_factors = "one", q.n_groups = "two", q.covariates = "no",
    q.binary_margins = margins, q.direction = direction, q.missingness = "none"))
}

test_that("Boschloo planning preserves selected variant and agrees with elementary exact cases", {
  two <- ap_bos()
  one <- ap_bos("one_sided_prespecified")
  sessions <- list(two, advisor_override(two, "boschloo_two_sided_minlike"), one, one)
  variants <- c("central", "minlike", "less", "greater")
  alternatives <- list(NULL, "two.sided", "less", "greater")
  # With n1=n2=2 these alpha values reject just the two extreme tables
  # (two-sided) or the appropriate single extreme table (one-sided).
  p1 <- .2; p2 <- .7
  lower <- (1 - p1)^2 * p2^2
  upper <- p1^2 * (1 - p2)^2
  for (i in seq_along(sessions)) {
    sp <- advisor_power_spec(sessions[[i]],
      list(p1 = p1, p2 = p2, alpha = if (i <= 2) .13 else .07, n1 = 2, n2 = 2),
      mode = "post_hoc", alternative = alternatives[[i]])
    expect_s3_class(sp, "power_spec")
    expect_identical(sp$options$variant, variants[i])
    r <- calculate_power(sp)
    expect_identical(r$status$code, "ok")
    expect_equal(r$solved$power, c(lower + upper, lower + upper, lower, upper)[i], tolerance = 1e-9)
    expect_identical(ig_from_json(ig_to_json(sp)), sp)
    expect_identical(advisor_restore(advisor_save(sessions[[i]])), sessions[[i]])
  }
})

test_that("a priori bridge retains explicit allocation and supported nuisance controls", {
  for (allocation in c("n2_ceiling", "total_nearest")) {
    sp <- advisor_power_spec(ap_bos(), list(p1 = 0, p2 = 1, alpha = .13,
      power = .8, allocation_ratio = 1, n_max = 8), allocation = allocation,
      control = boschloo_control(grid_points = 501L, nuisance_tol = 1e-8, tie_tolerance = 1e-8))
    expect_identical(sp$options$allocation, allocation)
    expect_equal(sp$inputs$allocation_ratio, 1)
    expect_equal(sp$inputs$grid_points, 501)
    expect_equal(sp$inputs$nuisance_tol, 1e-8)
    expect_equal(sp$inputs$tie_tolerance, 1e-8)
    r <- solve_design(sp)
    expect_identical(r$status$code, "ok")
    expect_equal(r$solved$n_total, 4)
    expect_equal(unname(r$group_sizes), c(2, 2))
    expect_equal(r$achieved$power, 1)
  }
})

test_that("supported nondefault controls reach the fixed-design calculation", {
  input <- list(p1=.2,p2=.7,alpha=.13,n1=2,n2=2)
  control <- boschloo_control(tie_tolerance=1e-8,grid_points=501L,nuisance_tol=1e-8,cache=FALSE)
  sp <- advisor_power_spec(ap_bos(),input,"post_hoc",control=control)
  r <- calculate_power(sp)
  expect_identical(r$status$code,"ok")
  for (name in c("tie_tolerance","grid_points","nuisance_tol")) {
    expect_equal(sp$inputs[[name]],control[[name]])
    expect_equal(r$method_outputs$settings[[name]],control[[name]])
  }
  with_cache <- advisor_power_spec(ap_bos(),input,"post_hoc",
    control=boschloo_control(tie_tolerance=1e-8,grid_points=501L,nuisance_tol=1e-8,cache=TRUE))
  expect_identical(with_cache,sp)
  expect_equal(calculate_power(with_cache)$solved$power,r$solved$power)
})

test_that("the bridge refuses unsupported, ambiguous or silently discarded inputs", {
  post <- list(p1 = .2, p2 = .7, alpha = .13, n1 = 2, n2 = 2)
  prior <- list(p1 = .2, p2 = .7, alpha = .13, power = .8, allocation_ratio = 1)
  expect_error(advisor_power_spec(advisor_session(), post, "post_hoc"), class = "inferguidR_unsupported")
  unknown <- ap_bos(margins = "dont_know")
  expect_identical(advisor_recommendation(unknown)$planning$planning_availability, "unavailable")
  expect_error(advisor_power_spec(unknown, post, "post_hoc"), class = "inferguidR_unsupported")
  fisher <- ap_bos(margins = "both_margins_fixed")
  expect_error(advisor_power_spec(fisher, post, "post_hoc"), class = "inferguidR_unsupported")
  expect_error(advisor_power_spec(ap_bos("one_sided_prespecified"), post, "post_hoc"), class = "inferguidR_invalid_input")
  expect_error(advisor_power_spec(ap_bos(), post, "post_hoc", "less"), class = "inferguidR_invalid_input")
  expect_error(advisor_power_spec(ap_bos(), prior), class = "inferguidR_invalid_input")
  expect_error(advisor_power_spec(ap_bos(), prior[-5], allocation = "n2_ceiling"), class = "inferguidR_invalid_input")
  expect_error(advisor_power_spec(ap_bos(), post, "post_hoc", allocation = "n2_ceiling"), class = "inferguidR_invalid_input")
  expect_error(advisor_power_spec(ap_bos(), c(post, list(grid_points = 101L)), "post_hoc"), class = "inferguidR_invalid_input")
  expect_error(advisor_power_spec(ap_bos(), post, "sensitivity"), class = "inferguidR_invalid_input")
  expect_error(advisor_power_spec(ap_bos(), post, "post_hoc",
    control = boschloo_control(polish_tol = 1e-10)), class = "inferguidR_unsupported")
  broken <- unclass(boschloo_control()); broken$grid_points <- 1
  class(broken) <- "boschloo_control"
  expect_error(advisor_power_spec(ap_bos(), post, "post_hoc", control = broken), class = "inferguidR_invalid_input")
  expect_error(advisor_power_spec(advisor_override(fisher, "boschloo_two_sided_central"),
    post, "post_hoc"), class = "inferguidR_unsupported")
})

test_that("old planning vocabulary stays readable and rules changes require explicit replay", {
  current <- ap_bos()
  saved <- jsonlite::parse_json(advisor_save(current), simplifyVector = FALSE)
  saved$rules$rules_version <- "2026-10-26.1"
  saved$state_hash <- NULL
  saved$state_hash <- ig_hash(saved)
  expect_error(advisor_restore(ig_to_json(saved)), class = "inferguidR_migration_required")
  answers <- setNames(vapply(saved$answers, `[[`, "", "option_id"),
    vapply(saved$answers, `[[`, "", "question_id"))
  replayed <- advisor_replay(answers, language = saved$language)
  expect_identical(replayed$state$steps, current$state$steps)
  expect_identical(advisor_recommendation(replayed)$planning$planning_availability, "available")
  for (availability in c("unavailable", "planned_independent_validation", "available")) {
    doc <- .ig_recommendation_to_list(advisor_recommendation(current))
    doc$planning$planning_availability <- availability
    doc$rules$rules_version <- if (availability == "available") "2026-10-26.2" else "2026-10-26.1"
    doc$recommendation_hash <- .ig_hash_value(.ig_recommendation_body(doc))
    restored <- ig_from_json(ig_to_json(doc))
    expect_identical(restored$planning$planning_availability, availability)
    expect_identical(restored$rules$rules_version, doc$rules$rules_version)
  }
})

test_that("available rules cannot relabel an unrelated variant as a mapped extension", {
  doc <- jsonlite::read_json(system.file("advisor", "rules.json", package = "inferguidR"))
  idx <- which(vapply(doc$outcomes, function(x) identical(x$outcome_id, "out.two_group.boschloo_central"), logical(1)))
  doc$outcomes[[idx]]$analysis_variant_id <- "boschloo_one_sided"
  f <- tempfile(fileext = ".json")
  on.exit(unlink(f), add = TRUE)
  jsonlite::write_json(doc, f, auto_unbox = TRUE, null = "null")
  v <- advisor_validate_rules(rules = f, walk = FALSE)
  expect_true("planning" %in% v$problems$check)
})
