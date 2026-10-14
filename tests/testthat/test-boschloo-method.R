# Boschloo's test in the generic analysis framework: registry entry, method
# definition, specifications, the generic post hoc and a priori solvers, and
# agreement with the dedicated Boschloo functions (same regions, same power).

bl_id <- "boschloo_proportions_two_independent"

cond_class <- function(expr) {
  tryCatch({
    force(expr)
    "no condition"
  }, error = function(e) class(e)[1])
}

bl_spec <- function(mode, inputs, options) {
  analysis_spec(bl_id, mode, inputs = inputs, options = options)
}

test_that("the registry lists Boschloo as an implemented extension without a G*Power reference", {
  expect_identical(.ig_registry_problems(), character())
  lm <- list_methods(family = "boschloo")
  expect_identical(lm$id, bl_id)
  expect_identical(lm$modes, "a_priori, post_hoc")
  expect_identical(lm$implementation_status, "implemented")
  expect_identical(lm$validation_status, "independently_checked")
  meta <- .ig_get_metadata(bl_id)
  expect_identical(meta$reference_status, "not_applicable")
  expect_null(meta$inventory_index)
  expect_identical(meta$owner_phase, "09A")
  expect_identical(meta$coverage$row_count, 0L)
  # Registry options and definition options use the same codes and defaults.
  def <- .ig_get_method(bl_id)
  for (nm in names(def$options)) {
    expect_identical(meta$options[[nm]]$codes, def$options[[nm]]$values, info = nm)
    expect_identical(meta$options[[nm]]$default, def$options[[nm]]$default, info = nm)
  }
  expect_identical(sort(names(meta$options)), sort(names(def$options)))
  expect_identical(def$modes, meta$modes)
  expect_true(def$discrete)
  expect_identical(def$checker, "validation/independent/boschloo")
  expect_identical(def$reference$planning_compatibility, "independently_validated_extension")
  expect_identical(def$reference$gpower_relation, "gpower_comparator_different_method")
  vs <- validation_status(bl_id)
  expect_identical(vs$implementation_status, "implemented")
  expect_identical(vs$independent_check$status, "passed")
  expect_identical(vs$gpower_reference$status, "not_applicable")
  expect_true(any(grepl("not part of the G\\*Power inventory", vs$explanation)))
  expect_false(any(grepl("No G\\*Power application reference captures", vs$explanation)))
})

test_that("the registry file is valid, with the null inventory index of an extension", {
  skip_if_not_installed("jsonvalidate")
  f <- system.file("registry", "boschloo.json", package = "inferguidR")
  txt <- readChar(f, file.size(f), useBytes = TRUE)
  res <- .ig_validate_schema(txt, "registry")
  # The registry schema allows a null inventory index only for entries whose
  # reference_status is not_applicable (extensions outside the inventory).
  expect_true(res$valid, info = paste(res$errors, collapse = "; "))
  forged <- sub('"reference_status": "not_applicable"', '"reference_status": "reference_pending"', txt,
                fixed = TRUE)
  expect_false(identical(forged, txt))
  expect_false(.ig_validate_schema(forged, "registry")$valid)
})

test_that("specifications require the variant and record every setting", {
  e <- tryCatch(bl_spec("post_hoc", list(p1 = 0.2, p2 = 0.7, alpha = 0.05, n1 = 10, n2 = 10), list()),
                error = function(e) e)
  expect_s3_class(e, "inferguidR_invalid_input")
  expect_match(paste(unlist(e$status$details$problems), collapse = " "), "variant")
  s <- bl_spec("post_hoc", list(p1 = 0.2, p2 = 0.7, alpha = 0.05, n1 = 10, n2 = 12),
               list(variant = "minlike"))
  expect_identical(s$validation_level, "method_definition")
  expect_identical(s$inputs$tie_tolerance, 1e-7)
  expect_identical(s$inputs$grid_points, 1001L)
  expect_identical(s$inputs$nuisance_tol, 1e-9)
  expect_setequal(s$defaults_applied$inputs, c("tie_tolerance", "grid_points", "nuisance_tol"))
  expect_null(s$options$allocation)                      # a priori only
  a <- bl_spec("a_priori", list(p1 = 0.2, p2 = 0.7, alpha = 0.05, power = 0.8), list(variant = "central"))
  expect_identical(a$options$allocation, "n2_ceiling")
  expect_identical(a$inputs$allocation_ratio, 1)
  expect_identical(a$inputs$n_max, 500L)
  expect_identical(cond_class(bl_spec("post_hoc", list(p1 = 1.2, p2 = 0.7, alpha = 0.05, n1 = 10, n2 = 10),
                                      list(variant = "less"))), "inferguidR_invalid_input")
  expect_identical(cond_class(bl_spec("post_hoc", list(p1 = 0.2, p2 = 0.7, alpha = 0.05, n1 = 10, n2 = 10),
                                      list(variant = "less", allocation = "total_nearest"))),
                   "inferguidR_invalid_input")
  expect_identical(cond_class(bl_spec("sensitivity", list(p1 = 0.2, alpha = 0.05, power = 0.8, n1 = 10,
                                                          n2 = 10), list(variant = "less"))),
                   "inferguidR_invalid_input")
  expect_identical(cond_class(bl_spec("post_hoc", list(p1 = 0.2, p2 = 0.7, alpha = 0.05, n1 = 2500,
                                                       n2 = 2500), list(variant = "less"))),
                   "inferguidR_invalid_input")
})

test_that("post hoc power equals boschloo_power() and reports the attained size", {
  cases <- list(list(20, 20, 0.2, 0.7, "central"), list(9, 14, 0.15, 0.6, "minlike"),
                list(12, 7, 0.7, 0.1, "greater"), list(8, 8, 0.3, 0.3, "less"),
                list(10, 15, 0, 1, "less"))
  for (cs in cases) {
    s <- bl_spec("post_hoc", list(p1 = cs[[3]], p2 = cs[[4]], alpha = 0.05, n1 = cs[[1]], n2 = cs[[2]]),
                 list(variant = cs[[5]]))
    r <- calculate_power(s)
    ref <- boschloo_power(cs[[1]], cs[[2]], cs[[3]], cs[[4]], alpha = 0.05, variant = cs[[5]])
    info <- paste(unlist(cs), collapse = " ")
    expect_identical(r$status$code, "ok", info = info)
    expect_identical(r$solved$power, ref$power, info = info)
    expect_identical(r$achieved$alpha_attained, ref$size$value, info = info)
    expect_equal(r$achieved$beta, 1 - ref$power, tolerance = 1e-12)
    expect_identical(r$method_outputs$region_hash, ref$region_hash, info = info)
    expect_identical(r$method_outputs$tables_rejected, as.double(ref$n_rejected))
    expect_identical(r$method_outputs$variant_id, ref$variant$variant_id)
    expect_identical(unname(r$group_sizes), as.integer(c(cs[[1]], cs[[2]])))
    expect_identical(r$method_outputs$planning_compatibility, "independently_validated_extension")
  }
  # Non-default settings reach the engine.
  s <- bl_spec("post_hoc", list(p1 = 0.2, p2 = 0.75, alpha = 0.05, n1 = 10, n2 = 14, grid_points = 501),
               list(variant = "minlike"))
  r <- calculate_power(s)
  ref <- boschloo_power(10, 14, 0.2, 0.75, variant = "minlike", control = boschloo_control(grid_points = 501L))
  expect_identical(r$solved$power, ref$power)
  expect_identical(r$method_outputs$settings$grid_points, 501)
})

test_that("the generic a priori search finds the first design of boschloo_sample_size()", {
  cases <- list(list(0.2, 0.7, "central", "total_nearest", 1),
                list(0.6, 0.1, "greater", "n2_ceiling", 2),
                list(0.1, 0.55, "minlike", "total_nearest", 0.5),
                list(0.05, 0.6, "less", "n2_ceiling", 1))
  for (cs in cases) {
    a <- bl_spec("a_priori", list(p1 = cs[[1]], p2 = cs[[2]], alpha = 0.05, power = 0.8,
                                  allocation_ratio = cs[[5]]),
                 list(variant = cs[[3]], allocation = cs[[4]]))
    r <- solve_design(a)
    ss <- boschloo_sample_size(cs[[1]], cs[[2]], power = 0.8, alpha = 0.05, variant = cs[[3]],
                               allocation = cs[[4]], allocation_ratio = cs[[5]], stability_window = 2L)
    info <- paste(unlist(cs), collapse = " ")
    expect_identical(r$status$code, "ok", info = info)
    expect_identical(r$solved$n_total, ss$first$n_total, info = info)
    expect_identical(unname(r$group_sizes), c(ss$first$n1, ss$first$n2), info = info)
    expect_identical(r$achieved$power, ss$first$power, info = info)
    expect_identical(r$achieved$alpha_attained, ss$first$size, info = info)
    # The method's own search (Neyman-Pearson exclusion, then enumeration).
    expect_identical(r$diagnostics$solver, "a_priori_method_search", info = info)
    expect_match(r$diagnostics$search, "Neyman-Pearson")
    expect_match(r$diagnostics$minimality, "^proven")
    expect_identical(r$diagnostics$designs_excluded_by_bound,
                     as.double(ss$minimality$designs_excluded_by_bound), info = info)
    expect_identical(r$provenance$solver, "a_priori_method_search")
  }
  # The generic result reports the stable design of the dedicated search
  # (default window K = 10) and survives a JSON round trip.
  a <- bl_spec("a_priori", list(p1 = 0.2, p2 = 0.7, alpha = 0.05, power = 0.8),
               list(variant = "central", allocation = "total_nearest"))
  r <- solve_design(a)
  ss <- boschloo_sample_size(0.2, 0.7, power = 0.8, alpha = 0.05, variant = "central",
                             allocation = "total_nearest")
  expect_identical(r$solved$n_total, ss$first$n_total)
  expect_identical(r$diagnostics$stable_n_total, as.double(ss$stable$n_total))
  expect_identical(r$diagnostics$stability_window, 10)
  expect_true(r$diagnostics$stable_found)
  expect_identical(r$diagnostics$evaluations, as.double(ss$minimality$designs_evaluated))
  expect_identical(r$validation$gpower_reference_status, "not_applicable")
  back <- ig_from_json(ig_to_json(r))
  expect_identical(back$result_hash, r$result_hash)
  expect_identical(back$diagnostics$stable_n_total, r$diagnostics$stable_n_total)
})

test_that("unanswerable a priori requests return typed statuses without numbers", {
  run <- function(inputs, options) solve_design(bl_spec("a_priori", inputs, options))
  r1 <- run(list(p1 = 0.4, p2 = 0.2, alpha = 0.05, power = 0.8), list(variant = "less"))
  expect_identical(r1$status$code, "invalid_input")
  expect_match(r1$status$reason, "does not match")
  r2 <- run(list(p1 = 0.3, p2 = 0.3, alpha = 0.05, power = 0.8), list(variant = "central"))
  expect_identical(r2$status$code, "infeasible")
  r3 <- run(list(p1 = 0.3, p2 = 0.4, alpha = 0.05, power = 0.95, n_max = 40), list(variant = "central"))
  expect_identical(r3$status$code, "infeasible")
  expect_match(r3$status$reason, "Neyman-Pearson")
  r4 <- run(list(p1 = 0.1, p2 = 0.6, alpha = 0.05, power = 0.04), list(variant = "central"))
  expect_identical(r4$status$code, "invalid_input")
  # The bound of the search itself: power 0.8 needs N = 30 here (see above).
  r5 <- run(list(p1 = 0.2, p2 = 0.7, alpha = 0.05, power = 0.8, n_max = 20),
            list(variant = "central", allocation = "total_nearest"))
  expect_identical(r5$status$code, "infeasible")
  for (r in list(r1, r2, r3, r4, r5)) {
    expect_null(r$solved$n_total)
    expect_identical(length(.ig_finite_paths(r$achieved, "achieved")), 0L)
  }
})

test_that("Boschloo planning results travel in envelopes and survive JSON", {
  s <- bl_spec("post_hoc", list(p1 = 0.25, p2 = 0.65, alpha = 0.05, n1 = 11, n2 = 13),
               list(variant = "central"))
  r <- calculate_power(s)
  env <- analysis_envelope(s, r)
  back <- ig_from_json(ig_to_json(env))
  expect_identical(back$analysis_hash, env$analysis_hash)
  expect_identical(back$result$solved$power, r$solved$power)
  expect_identical(ig_from_json(ig_to_json(s))$spec_hash, s$spec_hash)
})

test_that("allocation sequences map totals back to designs", {
  for (al in c("n2_ceiling", "total_nearest")) {
    for (r in c(0.3, 1, 1.5, 4)) {
      d <- .bl_design_at(1:40, r, al)
      tot <- d$n1 + d$n2
      for (k in seq_along(tot)) expect_identical(.bl_index_of_total(tot[k], r, al), k)
      gaps <- setdiff(seq(min(tot), max(tot)), tot)
      for (N in gaps) expect_true(is.na(.bl_index_of_total(N, r, al)))
    }
  }
  # The first design may already exceed n_max.
  expect_identical(.bl_max_index(2L, 2, "n2_ceiling"), 0L)
  expect_identical(.bl_max_index(3L, 2, "n2_ceiling"), 1L)
  expect_identical(boschloo_sample_size(0.1, 0.9, variant = "central", allocation_ratio = 2,
                                        n_max = 2L)$status$code, "infeasible")
})
