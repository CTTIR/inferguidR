test_that("the loaded package registers its engines and solvers without problems", {
  # Real load state (before th_registry_only() below): the phase 05 slice and
  # the Boschloo extension have definitions, and all five mode solvers exist.
  expect_identical(.ig_registry_problems(), character())
  expect_setequal(names(.ig_state$solvers), .ig_modes)
  defs <- names(.ig_state$definitions)
  expect_true(all(defs %in% .ig_list_method_ids()))
  expect_true(all(c("t_means_two_independent", th_extension_ids) %in% defs))
  lm <- list_methods()
  expect_identical(lm$implementation_status[lm$id %in% defs],
                   vapply(lm$id[lm$id %in% defs], .ig_implementation_status, "", USE.NAMES = FALSE))
  expect_true(all(lm$implementation_status[!lm$id %in% defs] == "not_implemented"))
  expect_true(all(lm$validation_status[!lm$id %in% c(defs, th_extension_ids)] == "reference_pending"))
  expect_false(any(lm$validation_status == "gpower_validated"))
  for (id in th_extension_ids) {
    expect_identical(validation_status(id)$gpower_reference$status, "not_applicable", info = id)
  }
  for (id in setdiff(defs, th_extension_ids)) {
    expect_identical(validation_status(id)$gpower_reference$status, "reference_pending", info = id)
  }
})

# The remaining tests describe the foundation behaviour for methods without an
# engine; they run against the bare registry (see th_registry_only()).
th_registry_only()

test_that("the registry loads the 49 inventory method ids plus the extensions", {
  ids <- .ig_list_method_ids()
  expect_length(ids, 49L + length(th_extension_ids))
  expect_setequal(ids, c(th_inventory_ids, th_extension_ids))
  expect_false(anyDuplicated(ids) > 0)
  expect_identical(.ig_registry_problems(), character())
  for (id in th_extension_ids) {
    meta <- .ig_get_metadata(id)
    expect_identical(meta$reference_status, "not_applicable", info = id)
    expect_null(meta$inventory_index)
  }
  for (id in th_inventory_ids) {
    meta <- .ig_get_metadata(id)
    expect_identical(meta$reference_status, "reference_pending", info = id)
    expect_true(is.numeric(meta$inventory_index) && meta$inventory_index >= 1, info = id)
  }
})

test_that("family counts match the inventory", {
  lm <- list_methods()
  counts <- table(lm$family)
  expect_identical(as.integer(counts[c("exact", "F", "t", "z", "chi2")]),
                   c(9L, 17L, 12L, 8L, 3L))
  expect_identical(as.integer(counts[["boschloo"]]), 1L)
})

test_that("list_methods returns the documented columns in inventory order", {
  lm <- list_methods()
  expect_s3_class(lm, "data.frame")
  expect_identical(names(lm), c("id", "family", "label", "modes",
                                "implementation_status", "validation_status"))
  # Inventory procedures first, in inventory order; extensions (null
  # inventory index) last.
  expect_identical(lm$id, c(th_inventory_ids, th_extension_ids))
  expect_true(all(vapply(lm, is.character, logical(1))))
  expect_true(all(lm$implementation_status == "not_implemented"))
  inv <- lm$id %in% th_inventory_ids
  expect_true(all(lm$validation_status[inv] == "reference_pending"))
  expect_true(all(lm$validation_status[!inv] == "not_applicable"))
  expect_false(any(grepl("gpower_validated", unlist(lm))))
  expect_identical(nrow(list_methods(family = "chi2")), 3L)
  expect_identical(th_condition_code(list_methods(family = "nope")), "invalid_input")
})

test_that("generic procedures do not offer the a priori mode", {
  lm <- list_methods()
  generic <- lm[grepl("_generic$", lm$id) & lm$id != "exact_generic_binomial", ]
  expect_identical(nrow(generic), 4L)
  expect_false(any(grepl("a_priori", generic$modes)))
})

test_that("registry metadata keeps both navigation paths and machine codes", {
  meta <- .ig_get_metadata("exact_proportions_fisher")
  expect_true(nzchar(meta$navigation$windows_3197))
  expect_true(nzchar(meta$navigation$mac_3196))
  expect_false(identical(meta$navigation$mac_3196, "same"))
  tails <- meta$options$tails
  expect_identical(tails$codes, c("one", "two"))
  # Certified inventory: the live default of Tail(s) is "One".
  expect_identical(tails$default, "one")
  expect_match(tails$default_evidence, "live capture")
  balancing <- meta$options$two_sided_balancing
  expect_true("symmetric_fisher_statistic" %in% balancing$codes)
  expect_identical(meta$options$allocation_ratio$value_type, "continuous")
  expect_identical(meta$options$allocation_ratio$domain$upper, Inf)
  expect_identical(.ig_get_metadata("exact_proportion_one_sample_binomial")$options$alpha_balancing$default,
                   "half_each_side")
  expect_identical(meta$implementation_status, "not_implemented")
  expect_identical(meta$reference_status, "reference_pending")
})

test_that("certified option defaults come from the live controls, mapped to codes", {
  # Every inventory choice option has a default derived from the live
  # control state; only the extension's variant must still be supplied.
  for (id in th_inventory_ids) {
    meta <- .ig_get_metadata(id)
    expect_identical(meta$inventory_status, "certified", info = id)
    for (o in meta$options) {
      if (!identical(o$value_type, "choice")) next
      expect_true(.ig_is_string(o$default) && o$default %in% o$codes, info = paste(id, o$option_id))
      # Each value names the live control that selects it.
      expect_setequal(names(o$value_live), o$codes)
      expect_true(all(vapply(o$value_live, function(v) .ig_is_string(v$text), logical(1))),
                  info = paste(id, o$option_id))
    }
  }
  corr <- .ig_get_metadata("exact_correlation_bivariate_normal")$options$computation
  # The inventory default is the live radio text; the code is the machine code.
  expect_identical(corr$default, "exact")
  expect_identical(corr$value_live$exact$text, "use exact distribution if  N  <")
  expect_identical(corr$value_live$exact$match, "position")
  expect_identical(corr$live$secondary_controls[[1]]$default, "10000")
  # A G*Power typo in a radio text is matched by a small edit distance, not guessed.
  bal <- .ig_get_metadata("exact_generic_binomial")$options$alpha_balancing
  expect_identical(bal$value_live$half_each_side$text, "Assing α/2 to both sides")
  expect_identical(bal$value_live$half_each_side$match, "near_text")
  # The expanded Fisher balancing choice takes its default from the combo item.
  fb <- .ig_get_metadata("exact_proportions_fisher")$options$two_sided_balancing
  expect_identical(fb$default, "symmetric_fisher_statistic")
  expect_identical(fb$value_live$symmetric_lr_statistic$item, "Likelihood ratio")
  # Codes stay stable where certification replaced the label (coverage key).
  ud <- .ig_get_metadata("exact_proportions_unconditional")$options$effect_direction
  expect_identical(ud$codes, c("p_h1_le_p_h0", "p_h1_ge_p_h0"))
  expect_identical(ud$value_labels, c("p1≤p2", "p1≥p2"))
  expect_identical(ud$default, "p_h1_ge_p_h0")
  # Check boxes: fixed (disabled) and free.
  expect_identical(.ig_get_metadata("f_rm_manova_within_between")$options$use_mean_correlation$codes, "off")
  expect_identical(.ig_get_metadata("z_proportions_two_independent")$options$continuity_correction$default, "off")
  # Continuous options carry the live field default and the analysis types showing them.
  ar <- .ig_get_metadata("t_means_two_independent")$options$allocation_ratio
  expect_null(ar$default)
  expect_identical(ar$default_value, 1)
  expect_identical(.ig_flatten(ar$applies_when$condition$mode, "character"), "a_priori")
  expect_null(.ig_get_metadata("z_poisson_regression")$options$over_dispersion$default_value)
  expect_false(.ig_get_metadata("z_poisson_regression")$options$over_dispersion$live$present)
})

test_that("live input and output labels per analysis type are carried", {
  ll <- .ig_live_labels("t_means_two_independent", "a_priori")
  expect_identical(ll$inputs, c("Tail(s)", "Effect size d", "α err prob",
                                "Power (1-β err prob)", "Allocation ratio N2/N1"))
  expect_identical(ll$outputs, c("Noncentrality parameter δ", "Critical t", "Df",
                                 "Sample size group 1", "Sample size group 2",
                                 "Total sample size", "Actual power"))
  expect_identical(unname(ll$input_defaults[["Tail(s)"]]), "One")
  expect_identical(.ig_live_labels("t_means_two_independent", "compromise")$outputs[4:6],
                   c("α err prob", "β err prob", "Power (1-β err prob)"))
  for (id in th_inventory_ids) {
    meta <- .ig_get_metadata(id)
    expect_setequal(names(meta$outputs_by_mode_live), meta$modes)
    expect_setequal(names(meta$inputs_by_mode_live), meta$modes)
  }
  expect_identical(th_condition_code(.ig_live_labels("f_generic", "a_priori")), "invalid_input")
  # Extensions have no G*Power labels.
  expect_null(.ig_live_labels(th_extension_ids[1], "post_hoc")$outputs)
})

test_that("coverage branches pruned or added by the live capture are recorded", {
  cov <- .ig_get_metadata("f_rm_manova_between")$coverage
  expect_identical(cov$row_count, 80L)
  expect_length(cov$branches_pruned, 0L)
  expect_length(cov$branches_added, 8L)
  expect_setequal(.ig_flatten(cov$branch_dimensions$use_mean_correlation, "character"), c("on", "off"))
  added <- .ig_get_metadata("z_proportions_two_independent")$coverage$branches_added
  expect_true(length(added) > 0L)
  expect_true(all(grepl("critical_value_scale=normal_difference", vapply(added, `[[`, "", "branch_spec"))))
  expect_length(.ig_get_metadata("t_means_two_independent")$coverage$branches_pruned, 0L)
})

test_that("extensions report that no G*Power reference exists", {
  for (id in th_extension_ids) {
    vs <- validation_status(id)
    expect_identical(vs$validation_status, "not_applicable", info = id)
    expect_identical(vs$gpower_reference$status, "not_applicable", info = id)
    expect_true(any(grepl("not part of the G\\*Power inventory", vs$explanation)), info = id)
    expect_false(any(grepl("No G\\*Power application reference captures", vs$explanation)), info = id)
    rec <- .ig_validation_record(id)
    expect_identical(rec$gpower_reference_status, "not_applicable", info = id)
  }
  rec <- .ig_validation_record("t_means_two_independent")
  expect_identical(rec$gpower_reference_status, "reference_pending")
})

test_that("unknown method ids are typed invalid_input", {
  expect_identical(th_condition_code(.ig_get_metadata("t_welch")), "invalid_input")
  expect_identical(th_condition_code(validation_status("t_welch")), "invalid_input")
})

test_that("validation_status reports the honest state", {
  vs <- validation_status("t_means_two_independent")
  expect_s3_class(vs, "ig_validation_status")
  expect_identical(vs$implementation_status, "not_implemented")
  expect_identical(vs$validation_status, "reference_pending")
  expect_identical(vs$gpower_reference$status, "reference_pending")
  expect_identical(vs$independent_check$status, "none")
  expect_true(length(vs$documented_example_ids) > 0L)
  expect_true(any(grepl("not claimed", vs$explanation)))
  expect_output(print(vs), "reference_pending")
})

test_that("method definitions are checked against the registry when registered", {
  expect_error(.ig_register_method(th_definition(id = "t_means_two_independent", family = "z")),
               "registry says 't'")
  expect_error(.ig_register_method(th_definition(id = "t_generic")), "absent from the reference inventory")
  bad <- th_definition()
  bad_id <- unclass(bad)
  bad_id$id <- "t_welch"
  class(bad_id) <- "ig_method"
  expect_error(.ig_register_method(bad_id), "no registry metadata")
  th_local_definition()
  expect_error(.ig_register_method(th_definition()), "already registered")
  expect_identical(list_methods()$implementation_status[list_methods()$id == "t_means_two_independent"],
                   "definition_only")
  expect_identical(validation_status("t_means_two_independent")$validation_status, "development")
})

test_that("definition providers are collected at load time", {
  ns <- new.env()
  ns$.ig_define_methods_test <- function() list(th_definition())
  ns$.ig_define_methods_broken <- function() stop("provider failure")
  saved_defs <- .ig_state$definitions
  saved_problems <- .ig_state$problems
  withr::defer({
    .ig_state$definitions <- saved_defs
    .ig_state$problems <- saved_problems
  })
  .ig_state$definitions <- list()
  .ig_state$problems <- character()
  .ig_collect_definitions(ns)
  expect_false(is.null(.ig_get_method("t_means_two_independent")))
  expect_match(.ig_state$problems, "provider failure")
})

test_that("registry files validate against the registry schema", {
  skip_if_not_installed("jsonvalidate")
  files <- list.files(system.file("registry", package = "inferguidR"), full.names = TRUE)
  # Five inventory families plus the Boschloo extension.
  expect_length(files, 6L)
  expect_setequal(basename(files), paste0(c("chi2", "exact", "F", "t", "z", "boschloo"), ".json"))
  for (f in files) {
    txt <- readChar(f, file.size(f), useBytes = TRUE)
    res <- .ig_validate_schema(txt, "registry")
    expect_true(res$valid, info = paste(basename(f), paste(res$errors, collapse = "; ")))
  }
})

test_that("registry files are in sync with the inventory (development checkout only)", {
  root <- th_pkg_root()
  skip_if(is.null(root), "inventory not available (installed package)")
  tool <- file.path(root, "project", "tools", "build_registry.R")
  skip_if_not(file.exists(tool), "registry builder not available")
  env <- new.env()
  sys.source(tool, envir = env)
  reg <- withr::with_dir(root, env$build_all())
  for (f in names(reg)) {
    path <- file.path(root, "inst", "registry", paste0(f, ".json"))
    current <- readChar(path, file.size(path), useBytes = TRUE)
    Encoding(current) <- "UTF-8"
    expect_identical(current, env$to_json(reg[[f]]), info = f)
  }
  inv <- jsonlite::fromJSON(file.path(root, "specs", "gpower-capabilities.json"), simplifyVector = FALSE)
  expect_setequal(vapply(inv, `[[`, "", "method_id"), setdiff(.ig_list_method_ids(), th_extension_ids))
  # Extensions are outside the inventory and say so.
  expect_length(intersect(vapply(inv, `[[`, "", "method_id"), th_extension_ids), 0L)
  for (id in th_extension_ids) {
    expect_identical(.ig_get_metadata(id)$reference_status, "not_applicable", info = id)
  }
})

test_that("the registry builder maps live controls robustly and fails on ambiguity (development checkout only)", {
  root <- th_pkg_root()
  skip_if(is.null(root), "inventory not available (installed package)")
  tool <- file.path(root, "project", "tools", "build_registry.R")
  skip_if_not(file.exists(tool), "registry builder not available")
  env <- new.env()
  sys.source(tool, envir = env)
  # Normalised text: comparison signs, Greek letters, '*', case and spacing.
  expect_identical(env$norm_text("r ≤ ρ"), env$norm_text("r <= rho"))
  expect_identical(env$norm_text("as in GPower 3.0"), env$norm_text("as in G*Power 3.0"))
  labels <- c("Use exact distribution if N < x", "Use large sample approximation")
  idx <- env$match_live("t", c("Use large sample approximation", "Use exact distribution if  N  <"), labels)
  expect_identical(as.integer(idx), c(2L, 1L))
  expect_identical(attr(idx, "how"), c("text", "near_text"))
  # A live text that is a prefix of two allowed values is ambiguous.
  expect_error(env$match_live("t", "Use exact distribution",
                              c("Use exact distribution if N < x", "Use exact distribution always")),
               "prefix of several")
  # Unmatched texts are placed by position only when every text match agrees with the order.
  expect_error(env$match_live("t", c("Use large sample approximation", "zzz unknown control"), labels),
               "contradict the control order")
  expect_error(env$match_live("t", "zzz unknown control", labels), "counts differ")
  # A default that disagrees with the live control state stops the build.
  opt <- list(option_id = "computation", default = "Use large sample approximation",
              allowed_values = list("Use exact distribution if N < x", "Use large sample approximation"),
              live = list(present = TRUE, where = "Options dialog", controls = list(
                list(id = 1L, kind = "radio", text = "Use exact distribution if  N  <", checked_by_default = TRUE),
                list(id = 2L, kind = "radio", text = "Use large sample approximation", checked_by_default = FALSE))))
  vals <- env$option_values("exact_proportions_fisher", opt)
  expect_error(env$resolve_live_choice("exact_proportions_fisher", opt, vals), "disagrees with the live default")
  opt$default <- "Use exact distribution if  N  <"
  expect_identical(env$resolve_live_choice("exact_proportions_fisher", opt, vals)$default, "exact")
})

test_that("a broken registry directory is reported, not hidden", {
  expect_error(.ig_load_registry_files(tempfile("no-registry-")), "missing")
  dir <- withr::local_tempdir()
  writeLines('{"schema_version":"inferguidR.registry/9.9.9","family":"t","methods":[]}',
             file.path(dir, "t.json"))
  expect_error(.ig_load_registry_files(dir), "schema_version")
})

test_that("RM-MANOVA metadata retains the restored mean-correlation branches", {
  for (id in c("f_rm_manova_between", "f_rm_manova_within")) {
    meta <- .ig_get_metadata(id)
    opt <- Filter(function(o) o$option_id == "use_mean_correlation", meta$options)[[1]]
    expect_setequal(opt$codes, c("on", "off"))
    expect_identical(meta$coverage$row_count, 80L)
  }
  meta <- .ig_get_metadata("f_rm_manova_within_between")
  opt <- Filter(function(o) o$option_id == "use_mean_correlation", meta$options)[[1]]
  expect_identical(opt$codes, "off")
  expect_identical(meta$coverage$row_count, 40L)
})
