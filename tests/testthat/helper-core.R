# Helpers for the core tests (phase 04).

# The 49 procedure ids of the reference inventory (specs/gpower-capabilities.json,
# provisional 2026-10-03, certified against the live application 2026-10-05;
# the ids did not change).
th_inventory_ids <- c(
  "exact_correlation_bivariate_normal", "exact_proportion_one_sample_binomial",
  "exact_proportion_sign_test", "exact_proportions_mcnemar", "exact_proportions_fisher",
  "exact_proportions_unconditional", "exact_proportions_unconditional_offset",
  "exact_mrc_random_model", "exact_generic_binomial",
  "f_anova_oneway", "f_anova_special_main_interactions", "f_ancova_main_interactions",
  "f_rm_anova_between", "f_rm_anova_within", "f_rm_anova_within_between",
  "f_rm_manova_between", "f_rm_manova_within", "f_rm_manova_within_between",
  "f_hotelling_t2_one_group", "f_hotelling_t2_two_groups", "f_manova_global",
  "f_manova_special", "f_mrc_r2_deviation_from_zero", "f_mrc_r2_increase",
  "f_variances_two_groups", "f_generic",
  "t_means_one_sample", "t_means_matched_pairs", "t_means_two_independent",
  "t_wilcoxon_signed_rank_one_sample", "t_wilcoxon_signed_rank_matched_pairs",
  "t_wilcoxon_mann_whitney", "t_correlation_point_biserial", "t_linreg_one_group_slope",
  "t_linreg_two_groups_slopes", "t_linreg_two_groups_intercepts",
  "t_mrc_single_coefficient", "t_generic",
  "z_logistic_regression", "z_poisson_regression", "z_tetrachoric_correlation",
  "z_correlations_two_dependent_common_index", "z_correlations_two_dependent_no_common_index",
  "z_correlations_two_independent", "z_proportions_two_independent", "z_generic",
  "chi2_variance_one_sample", "chi2_goodness_of_fit_contingency", "chi2_generic"
)

# Registry entries that are inferguidR extensions, not procedures of the
# reference inventory (reference_status "not_applicable", inventory_index
# null).
th_extension_ids <- c("boschloo_proportions_two_independent")

# Reproduce the foundation state (registry metadata only: no method
# definitions and no mode solvers) for the duration of the calling test or
# test file, so that registry-only behaviour stays testable after engines
# and solvers are registered by later phases. Called at the top level of a
# test file, the state is restored when the file has run.
th_registry_only <- function(env = parent.frame()) {
  old_defs <- .ig_state$definitions
  old_solvers <- .ig_state$solvers
  .ig_state$definitions <- list()
  .ig_state$solvers <- list()
  withr::defer({
    .ig_state$definitions <- old_defs
    .ig_state$solvers <- old_solvers
  }, envir = env)
  invisible(NULL)
}

# Options a registry-only specification needs: the first code of every
# applicable result-affecting choice option without a documented default
# (same applicability rules as the validator). Since the inventory was
# certified, only options without a live default (the extension's variant)
# are chosen here; all others take their default.
th_minimal_options <- function(method_id, mode) {
  meta <- .ig_get_metadata(method_id)
  specs <- Filter(function(o) identical(o$value_type, "choice") && isTRUE(o$result_affecting),
                  meta$options)
  chosen <- list()
  for (iter in seq_len(length(specs) + 2L)) {
    eff <- chosen
    for (o in specs) {
      if (!is.null(o$default) && is.null(eff[[o$option_id]])) eff[[o$option_id]] <- o$default
    }
    new <- list()
    for (o in specs) {
      aw <- o$applies_when
      if (!is.null(aw) && is.null(aw$condition)) next # applicability unverified: optional
      if (!.ig_eval_condition(if (is.null(aw)) NULL else aw$condition, mode, eff)) next
      if (is.null(o$default)) new[[o$option_id]] <- o$codes[1]
    }
    if (identical(new, chosen)) break
    chosen <- new
  }
  chosen
}

# Standard inputs each mode needs for a registry-only specification.
th_minimal_inputs <- function(mode) {
  switch(mode,
         a_priori = list(alpha = 0.05, power = 0.8),
         post_hoc = list(alpha = 0.05, n_total = 40L),
         sensitivity = list(alpha = 0.05, power = 0.8, n_total = 40L),
         criterion = list(power = 0.8, n_total = 40L),
         compromise = list(n_total = 40L, error_ratio = 1))
}

th_minimal_spec <- function(method_id, mode) {
  target <- NULL
  if (mode == "sensitivity") {
    alts <- .ig_registry_targets(.ig_get_metadata(method_id), mode)
    if (length(alts) > 1L) target <- alts[1]
  }
  analysis_spec(method_id, mode, inputs = th_minimal_inputs(mode),
                options = th_minimal_options(method_id, mode), target = target)
}

# A complete test definition for a real inventory id. Its forward model is
# deliberately not implemented: it must never produce a number.
th_definition <- function(id = "t_means_two_independent", family = "t",
                          modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
                          status = "development") {
  ig_method(
    id = id, family = family,
    label = list(en = "Test definition (no forward model)"),
    params = ig_params(
      ig_param("d", role = "effect", domain = ig_domain_real(),
               label = list(en = "Effect size d"), aliases = "effect_size_d"),
      ig_param_alpha(), ig_param_power(), ig_param_n_total(lower = 4),
      ig_param_error_ratio(),
      ig_param("allocation_ratio", role = "design",
               domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE),
               default = 1, aliases = "ratio", label = list(en = "Allocation ratio N2/N1")),
      ig_param("exact_threshold", role = "control", domain = ig_domain_integer(lower = 1),
               default = 1000L, applies_when = list(option = "computation", "in" = "exact"))
    ),
    options = ig_options(
      ig_option("tails", values = c("one", "two")),
      ig_option("computation", values = c("exact", "approximate"), default = "approximate"),
      ig_option("balancing", values = c("half", "minor_tail"), default = "half",
                applies_when = list(option = "tails", "in" = "two"))
    ),
    modes = modes,
    power = function(p) stop("test definition: no forward model"),
    design = function(p, n_total) stop("test definition: no design rule"),
    engine = list(package = "stats", fun = "pt", version = NA_character_),
    status = status,
    constraints = function(p) {
      if (!is.null(p$n_total) && identical(p$allocation_ratio, 1) && p$n_total %% 2L == 1L) {
        "with allocation ratio 1 the total sample size must be even"
      } else {
        character()
      }
    }
  )
}

# Register a definition for the duration of the calling test.
th_local_definition <- function(def = th_definition(), env = parent.frame()) {
  old <- .ig_get_method(def$id)
  .ig_register_method(def, replace = TRUE)
  withr::defer({
    if (is.null(old)) .ig_unregister_method(def$id) else .ig_register_method(old, replace = TRUE)
  }, envir = env)
  invisible(def)
}

# Install test solvers for the duration of the calling test.
th_local_solvers <- function(solvers, env = parent.frame()) {
  old <- .ig_state$solvers
  .ig_state$solvers <- solvers
  withr::defer(.ig_state$solvers <- old, envir = env)
  invisible(solvers)
}

# Replace (or add) list elements; a NULL value is kept as an explicit NULL.
th_override <- function(base, over) {
  for (nm in names(over)) base[nm] <- list(over[[nm]])
  base
}

# Status code carried by an inferguidR condition (or NA when none).
th_condition_code <- function(expr) {
  tryCatch({
    force(expr)
    NA_character_
  }, inferguidR_error = function(e) e$status$code)
}

# Paths of finite numbers in a list, excluding named top-level sections.
th_finite_numbers <- function(x, exclude = character()) {
  x <- unclass(x)
  x <- x[setdiff(names(x), exclude)]
  .ig_finite_paths(x, "")
}

th_schema_valid <- function(obj, schema) {
  res <- .ig_validate_schema(ig_to_json(obj), schema)
  if (!isTRUE(res$valid)) message(paste(res$errors, collapse = "\n"))
  isTRUE(res$valid)
}

th_pkg_root <- function() {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
  if (file.exists(file.path(root, "specs", "gpower-capabilities.json"))) root else NULL
}
