# Advisor recommendations: JSON round trip against the bundled schema, no
# numerical results, R-function and manual provenance, learning lookup and
# immutability, over every reachable outcome-status combination and every
# allowed override.

tr_welch_map <- c(q.aim = "compare", q.stage = "analysis", q.outcome_type = "continuous",
                  q.unit = "individuals", q.structure = "independent", q.independence_check = "no",
                  q.n_factors = "one", q.n_groups = "two", q.covariates = "no",
                  q.estimand_continuous = "mean", q.variance = "not_justified",
                  q.direction = "two_sided", q.missingness = "none")

# One recommendation per reachable outcome-status combination, plus one per
# allowed override of each method outcome (computed once per file).
tr_cache <- new.env(parent = emptyenv())
tr_all_recommendations <- function() {
  if (is.null(tr_cache$recs)) {
    rs <- advisor_rules()
    w <- .ig_advisor_enumerate(rs)
    recs <- list()
    for (combo in names(w$examples)) {
      s <- advisor_replay(w$examples[[combo]])
      recs[[combo]] <- advisor_recommendation(s)
      for (a in rs$outcomes[[s$state$outcome_id]]$alternatives) {
        if (isTRUE(a$override_allowed)) {
          recs[[paste(combo, a$analysis_variant_id, sep = "|override=")]] <-
            advisor_recommendation(advisor_override(s, a$analysis_variant_id))
        }
      }
    }
    tr_cache$recs <- recs
  }
  tr_cache$recs
}

# Canonical JSON text of each recommendation (computed once per file).
tr_all_json <- function() {
  if (is.null(tr_cache$json)) tr_cache$json <- lapply(tr_all_recommendations(), ig_to_json)
  tr_cache$json
}

test_that("every recommendation round-trips through canonical JSON and matches its schema", {
  skip_if_not_installed("jsonvalidate")
  recs <- tr_all_recommendations()
  expect_length(recs, 189L) # 90 outcome-status combinations + 99 allowed overrides
  for (id in names(recs)) {
    rec <- recs[[id]]
    js <- tr_all_json()[[id]]
    sv <- .ig_validate_schema(unclass(js), "recommendation")
    expect_true(isTRUE(sv$checked) && isTRUE(sv$valid), info = paste(id, paste(sv$errors, collapse = "; ")))
    back <- ig_from_json(js, type = "recommendation", validate = FALSE)
    expect_identical(back, rec, info = id)
    expect_identical(ig_hash(back), rec$recommendation_hash, info = id)
  }
})

test_that("recommendations contain no numerical results and no claimed parity", {
  recs <- tr_all_recommendations()
  # Integers that index the trace or cite a manual page are the only numbers.
  allowed <- "\\$(step|transition_order|manual_page)$"
  for (id in names(recs)) {
    rec <- recs[[id]]
    nums <- th_finite_numbers(rec)
    expect_true(all(grepl(allowed, nums)), info = paste(id, paste(nums[!grepl(allowed, nums)], collapse = ", ")))
    expect_false(.ig_contains_result(unclass(rec)), info = id)
    expect_false(rec$planning$planning_compatibility == "exact_gpower_match", info = id)
    if (rec$planning$planning_compatibility == "independently_validated_extension")
      expect_true(rec$preferred$analysis_variant_id %in% c("boschloo_two_sided_central",
        "boschloo_two_sided_minlike", "boschloo_one_sided"), info = id)
    expect_length(rec$planning$spec_refs, 0L)
    expect_length(rec$planning$evidence_ids, 0L)
    expect_null(rec$planning$gpower_method_id)
    rel <- rec$preferred$gpower_manual_relation
    if (!is.null(rel)) expect_match(rel$parity_note_en, "no parity is claimed", info = id)
  }
})

test_that("the R implementation is the preferred variant's template id and exports no code", {
  rs <- advisor_rules()
  recs <- tr_all_recommendations()
  fields <- c("note_de", "note_en", "rules_template_status", "template_id", "template_status")
  for (id in names(recs)) {
    rec <- recs[[id]]
    ri <- rec$r_implementation
    expect_setequal(names(ri), fields)
    if (is.null(rec$preferred)) {
      expect_null(ri$template_id)
      expect_identical(ri$template_status, "not_applicable")
      next
    }
    vid <- rec$preferred$analysis_variant_id
    expect_identical(ri$template_id, paste0("tpl.", vid, ".v1"), info = id)
    available <- vid %in% c("t_means_two_independent_welch", "t_means_two_independent_student", "t_means_paired",
        "t_mean_one_sample", "binomial_exact_one_sample", "prop_one_sample_score",
        "correlation_pearson", "correlation_spearman", "correlation_kendall",
        "sign_test_one_sample", "sign_test_paired", "mcnemar_exact_conditional",
        "mcnemar_midp", "fisher_exact_2x2", "boschloo_two_sided_central",
        "boschloo_two_sided_minlike", "boschloo_one_sided", "signed_rank_wilcoxon_one_sample",
        "signed_rank_wilcoxon", "rank_sum_mann_whitney", "chisq_goodness_of_fit",
        "chisq_independence_rxc", "prop_two_independent_pearson_chisq",
        "mcnemar_asymptotic_uncorrected", "mcnemar_asymptotic_cc",
        "variance_f_test_two", "variance_brown_forsythe", "variance_fligner_killeen",
        "anova_one_way", "anova_one_way_welch", "kruskal_wallis", "descriptive_summary", "fisher_exact_rxc", "friedman", "cochran_q", "stuart_maxwell_marginal_homogeneity", "bowker_symmetry", "regression_linear_simple", "regression_linear_multiple", "pca_components", "diagnostic_2x2_fixed_cutoff", "brunner_munzel", "chisq_goodness_of_fit_mc", "cochran_armitage_trend", "classifier_binormal_theory", "classifier_cutoff_exploration", "classifier_predeclared_rule", "classifier_heldout_validation", "ancova_baseline_adjusted", "jonckheere_terpstra", "unconditional_exact_z_pooled", "anova_factorial", "anova_repeated_one_way_gg", "regression_logistic_binary", "regression_poisson", "regression_negative_binomial", "regression_ordinal_proportional_odds", "regression_multinomial", "lmm_repeated_measures", "lmm_cluster_random_intercept", "glmm_logistic_random_intercept", "gee_logistic_exchangeable", "prediction_model_linear", "prediction_model_logistic", "efa_common_factor_ml", "cluster_hierarchical_ward", "cluster_kmeans", "cluster_gower_pam")
    expect_identical(ri$template_status, if (available) "available" else "template_not_yet_available", info = id)
    expect_match(ri$note_en, if (available) "bundled with the package" else "not yet available; no substitute code", info = id)
    # The call sketches of analysis-variants.json are drafts: never exported
    sketches <- vapply(rs$variants[[vid]]$r_functions, `[[`, "", "call_sketch")
    txt <- unclass(tr_all_json()[[id]])
    for (sk in sketches) expect_false(grepl(sk, txt, fixed = TRUE), info = id)
  }
  # In particular, the Boschloo code route never points to exact2x2 (decision D12;
  # the rationale may cite it as the definition, with the tie caveat as an assumption)
  bos <- Filter(function(r) grepl("^boschloo", r$preferred$analysis_variant_id %||% ""), recs)
  expect_true(length(bos) >= 2L)
  for (r in bos) {
    expect_false(grepl("exact2x2", unclass(ig_to_json(r$r_implementation)), fixed = TRUE))
    # nor do the recorded engine settings (rules 2026-10-05.1)
    expect_false(grepl("exact2x2", unclass(ig_to_json(r$preferred$variant_settings)), fixed = TRUE))
    if (!is.null(r$preferred$variant_settings$control)) {
      expect_match(r$preferred$variant_settings$control, "boschloo_control", fixed = TRUE)
    }
    expect_true("a.ordering_ties" %in% vapply(r$assumptions, `[[`, "", "assumption_id"))
  }
})

test_that("manual relation, learning links and source nodes resolve against the bundled data", {
  rs <- advisor_rules()
  ids <- vapply(rs$resources, `[[`, "", "resource_id")
  recs <- tr_all_recommendations()
  for (id in names(recs)) {
    rec <- recs[[id]]
    expect_true(all(rec$learning_links %in% ids), info = id)
    expect_true(all(unlist(rec$rules$source_nodes$uzh) %in% rs$uzh_nodes), info = id)
    expect_true(all(unlist(rec$rules$source_nodes$descrtab2) %in% rs$descrtab2_nodes), info = id)
    expect_identical(rec$rules$rules_version, "2026-10-26.2")
    expect_identical(rec$rules$analysis_variants_version, "2026-10-26.1")
    expect_match(rec$rules$learning_resources_version, "^sha256:[0-9a-f]{64}$")
    if (is.null(rec$preferred)) next
    # The planning note exists in both languages, also for an override
    expect_true(nzchar(rec$planning$note_en), info = id)
    expect_true(is.character(rec$preferred$planning_note_de) && nzchar(rec$preferred$planning_note_de), info = id)
    if (isTRUE(rec$preferred$overridden)) {
      # an override's planning note names the overriding variant in both languages
      vid <- rec$preferred$analysis_variant_id
      expect_true(grepl(vid, rec$planning$note_en, fixed = TRUE), info = id)
      expect_true(grepl(vid, rec$preferred$planning_note_de, fixed = TRUE), info = id)
      if (vid %in% c("boschloo_two_sided_minlike", "boschloo_two_sided_central") &&
          rec$outcome_id %in% c("out.two_group.boschloo_central", "out.two_group.boschloo_one_sided"))
        expect_identical(rec$planning$planning_compatibility, "independently_validated_extension")
      else expect_true(rec$planning$planning_compatibility %in% c("not_yet_validated", "not_applicable"), info = id)
    }
    for (role in c("course", "tutorial", "overview")) {
      blk <- rec$preferred$learning[[role]]
      expect_true(blk$status %in% c("exact", "related_only", "broader_only", "gap"), info = paste(id, role))
      if (identical(blk$status, "gap")) {
        expect_length(blk$links, 0L)
        expect_match(blk$note_en, "^Gap: ")
      } else {
        expect_true(length(blk$links) > 0L, info = paste(id, role))
      }
    }
  }
  welch <- advisor_recommendation(advisor_replay(tr_welch_map))
  rel <- welch$preferred$gpower_manual_relation
  expect_identical(rel$relation, "different_method_comparator")
  expect_identical(rel$manual_page, 52L)
  # Boschloo: G*Power has no Boschloo procedure; the manual documents only Fisher's
  # (conditional) test, a comparator with a different method (phase 09A, provisional)
  bos_map <- c(tr_welch_map[c("q.aim", "q.stage", "q.unit", "q.structure", "q.independence_check",
                              "q.n_factors", "q.n_groups", "q.covariates", "q.direction", "q.missingness")],
               q.outcome_type = "binary", q.binary_margins = "group_sizes_fixed")
  bos <- advisor_recommendation(advisor_replay(bos_map))
  expect_identical(bos$preferred$analysis_variant_id, "boschloo_two_sided_central")
  expect_identical(bos$preferred$gpower_manual_relation$relation, "different_method_comparator")
  expect_identical(bos$planning$planning_availability, "available")
})

test_that("overrides list the overriding variant's assumptions or state that there are none", {
  rs <- advisor_rules()
  recs <- tr_all_recommendations()
  n_override <- 0L
  for (id in names(recs)) {
    rec <- recs[[id]]
    ids <- vapply(rec$assumptions, `[[`, "", "assumption_id")
    expect_false(anyDuplicated(ids) > 0L, info = id)
    if (!isTRUE(rec$preferred$overridden)) next
    n_override <- n_override + 1L
    vid <- rec$preferred$analysis_variant_id
    homes <- Filter(function(o) identical(o$analysis_variant_id, vid), rs$outcomes)
    own <- Filter(function(a) identical(a$source, "override_variant"), rec$assumptions)
    expect_true(all(vapply(own, `[[`, "", "stated_for") == vid), info = id)
    if (length(homes)) {
      # every assumption the rules state for the variant is listed (here or earlier)
      want <- unique(unlist(lapply(homes, function(o) vapply(o$assumptions, `[[`, "", "assumption_id"))))
      expect_true(all(want %in% ids), info = id)
    } else {
      expect_identical(vapply(own, `[[`, "", "assumption_id"),
                       paste0("override.", vid, ".assumptions_not_stated"), info = id)
    }
  }
  expect_identical(n_override, 99L)
})

test_that("recommendations and sessions are immutable", {
  s <- advisor_replay(tr_welch_map)
  rec <- advisor_recommendation(s)
  expect_identical(th_condition_code(rec$status <- "recommended"), "invalid_input")
  expect_identical(th_condition_code(rec[["preferred"]] <- NULL), "invalid_input")
  expect_identical(th_condition_code(s$override <- "rank_sum_mann_whitney"), "invalid_input")
  expect_identical(th_condition_code(s[["answers"]] <- list()), "invalid_input")
  expect_identical(rec, advisor_recommendation(s))
})

test_that("restore refuses resealed states with malformed fields", {
  js <- advisor_save(advisor_replay(tr_welch_map))
  doc <- jsonlite::parse_json(js, simplifyVector = FALSE)
  reseal <- function(d) {
    d$state_hash <- NULL
    d$state_hash <- ig_hash(d)
    ig_to_json(d)
  }
  bad <- list(
    ignored_answers = function(d) { d$ignored_answers <- list(list(1)); d },
    language = function(d) { d$language <- list("en", "de"); d },
    override = function(d) { d$override <- list(a = 1); d },
    invalidated = function(d) { d$invalidated <- "everything"; d },
    expected = function(d) { d$expected <- list(node = 1); d },
    answers = function(d) { d$answers <- list("q.aim=compare"); d },
    missing = function(d) { d$expected <- NULL; d }
  )
  for (nm in names(bad)) {
    e <- tryCatch(advisor_restore(reseal(bad[[nm]](doc))), error = function(e) e)
    expect_s3_class(e, "inferguidR_invalid_input")
  }
  e <- tryCatch(advisor_restore(reseal(bad$missing(doc))), error = function(e) e)
  expect_match(conditionMessage(e), "lacks the field")
})
