# Advisor rule set: bundled data, validation and the exhaustive path walk.

ta_rules_doc <- function() {
  jsonlite::parse_json(paste(readLines(system.file("advisor", "rules.json", package = "inferguidR"),
                                       encoding = "UTF-8", warn = FALSE), collapse = "\n"),
                       simplifyVector = FALSE)
}

# An independent re-implementation of the path walker. It interprets the
# conditions of rules.json directly (no compiled closures, no package
# helpers) and enumerates every answer path depth-first. It also applies the
# status rules as worded in rules.json engine.status_rules (its own code, not
# .ig_advisor_status()) and keeps one record per path, so that statuses can be
# compared with the package.
ta_walk_all <- function(doc) {
  questions <- doc$questions
  names(questions) <- vapply(questions, `[[`, "", "question_id")
  outcomes <- doc$outcomes
  names(outcomes) <- vapply(outcomes, `[[`, "", "outcome_id")
  outcome_ids <- names(outcomes)
  by_from <- split(doc$transitions, vapply(doc$transitions, `[[`, "", "from"))
  by_from <- lapply(by_from, function(ts) ts[order(vapply(ts, function(t) t$order, numeric(1)))])
  holds <- function(cond, eff) {
    if (!is.null(cond$always)) return(isTRUE(cond$always))
    if (!is.null(cond$all)) return(all(vapply(cond$all, holds, logical(1), eff = eff)))
    if (!is.null(cond$any)) return(any(vapply(cond$any, holds, logical(1), eff = eff)))
    if (!is.null(cond$not)) return(!holds(cond$not, eff))
    a <- eff[[cond$q]]
    !is.null(a) && a %in% unlist(cond[["in"]])
  }
  out <- new.env()
  out$paths <- 0L
  out$max_len <- 0L
  out$outcomes <- character()
  out$questions <- character()
  out$transitions <- character()
  out$guards <- 0L
  out$cycles <- 0L
  out$not_applicable <- 0L
  out$records <- list()
  status_of <- function(oid, eff, dk_assume) {
    o <- outcomes[[oid]]
    st <- o$status
    if (length(dk_assume) && st == "recommended") st <- "conditional"
    mods <- character()
    for (m in o$status_modifiers) {
      if (holds(m$when, eff)) {
        mods <- c(mods, m$modifier_id)
        if (st == "recommended") st <- "conditional"
      }
    }
    list(status = st, modifiers = mods)
  }
  visit <- function(node, eff, raw, dk_assume, trace, seen) {
    out$questions <- union(out$questions, node)
    if (!holds(questions[[node]]$applies_when, eff)) out$not_applicable <- out$not_applicable + 1L
    for (opt in questions[[node]]$options) {
      e <- opt$option_id
      dk2 <- dk_assume
      if (identical(e, "dont_know") &&
          isTRUE(opt$dont_know_policy %in% c("route_without_assumption", "route_with_assumption"))) {
        e <- opt$routes_as
        if (identical(opt$dont_know_policy, "route_with_assumption")) dk2 <- c(dk2, node)
      }
      eff2 <- eff
      eff2[[node]] <- e
      raw2 <- c(raw, stats::setNames(opt$option_id, node))
      tr <- NULL
      for (t in by_from[[node]]) if (holds(t$when, eff2)) { tr <- t; break }
      out$transitions <- union(out$transitions, tr$transition_id)
      trace2 <- c(trace, tr$transition_id)
      if (identical(tr$kind, "guard")) out$guards <- out$guards + 1L
      if (tr$to %in% outcome_ids) {
        out$paths <- out$paths + 1L
        out$max_len <- max(out$max_len, length(eff2))
        out$outcomes <- union(out$outcomes, tr$to)
        st <- status_of(tr$to, eff2, dk2)
        out$records[[out$paths]] <- list(answers = raw2, outcome = tr$to, status = st$status,
                                         modifiers = st$modifiers, dk_assume = dk2, trace = trace2)
      } else if (tr$to %in% c(seen, node)) {
        out$cycles <- out$cycles + 1L
      } else {
        visit(tr$to, eff2, raw2, dk2, trace2, c(seen, node))
      }
    }
  }
  visit(doc$engine$start_question, list(), character(), character(), character(), character())
  as.list(out)
}

# The independent walk is computed once per file (about two seconds).
ta_cache <- new.env(parent = emptyenv())
ta_independent_walk <- function() {
  if (is.null(ta_cache$walk)) ta_cache$walk <- ta_walk_all(ta_rules_doc())
  ta_cache$walk
}

test_that("the bundled rule set loads with its versions and provenance", {
  rs <- advisor_rules()
  expect_s3_class(rs, "advisor_rules")
  v <- rs$versions
  expect_identical(v$rules_id, "advisor_rules")
  expect_identical(v$rules_version, "2026-10-26.2")
  expect_identical(v$uzh_tree_version, "2026-10-03.1")
  expect_identical(v$analysis_variants_version, "2026-10-26.1")
  expect_match(v$descrtab2_tree_version, "2026-10-03")
  expect_match(v$learning_resources_version, "^sha256:[0-9a-f]{64}$")
  expect_named(v$data_sha256, c("rules", "analysis_variants", "uzh_tree", "descrtab2_tree",
                                "source_crosswalk", "learning_resources"))
  expect_true(all(grepl("^[0-9a-f]{64}$", unlist(v$data_sha256))))
  expect_identical(advisor_rules(), rs) # cached
  expect_length(rs$questions, 27L)
  expect_length(rs$outcomes, 59L)
})

test_that("the bundled rule set passes every validation check", {
  v <- advisor_validate_rules()
  expect_s3_class(v, "advisor_rules_validation")
  if (nrow(v$problems)) print(v$problems)
  expect_true(v$valid)
  expect_identical(nrow(v$problems), 0L)
  expect_false(is.null(v$walk))
})

test_that("the exhaustive walk: all paths terminate, no dead ends, everything reachable, counts as documented", {
  doc <- ta_rules_doc()
  summary <- doc$summary
  rs <- advisor_rules()
  w <- .ig_advisor_enumerate(rs)
  # Package enumerator
  expect_identical(w$answer_paths, as.integer(summary$draft_check$answer_paths_enumerated))
  expect_identical(w$answer_paths, 5076L)
  expect_identical(w$max_questions_on_a_path, as.integer(summary$draft_check$max_questions_on_a_path))
  expect_identical(w$outcome_status_combinations,
                   as.integer(summary$draft_check$outcome_status_combinations_reachable))
  expect_identical(unname(w$issue_counts), c(0L, 0L, 0L, 0L)) # guard, applies_when, cycle, design
  expect_length(w$unreachable_outcomes, 0L)
  expect_length(w$unreachable_questions, 0L)
  expect_length(w$unused_transitions, 0L)
  expect_identical(sum(w$outcome_hits), w$answer_paths)
  # Independent walker written in this test
  iw <- ta_independent_walk()
  expect_identical(iw$paths, as.integer(summary$draft_check$answer_paths_enumerated))
  expect_identical(iw$max_len, w$max_questions_on_a_path)
  expect_identical(iw$guards, 0L)
  expect_identical(iw$cycles, 0L)
  expect_identical(iw$not_applicable, 0L)
  expect_setequal(iw$outcomes, vapply(doc$outcomes, `[[`, "", "outcome_id"))
  expect_setequal(iw$questions, vapply(doc$questions, `[[`, "", "question_id"))
  non_guard <- vapply(Filter(function(t) t$kind != "guard", doc$transitions), `[[`, "", "transition_id")
  expect_setequal(iw$transitions, non_guard)
  # Documented counts of the graph itself
  expect_identical(length(doc$questions), as.integer(summary$n_questions))
  expect_identical(length(doc$transitions), as.integer(summary$n_transitions))
  expect_identical(length(doc$outcomes), as.integer(summary$n_outcomes))
  expect_identical(sum(vapply(doc$transitions, function(t) t$kind == "guard", logical(1))),
                   as.integer(summary$n_guard_transitions))
})

test_that("statuses agree with the independent walker, and a fixed sample of paths replays identically", {
  rs <- advisor_rules()
  iw <- ta_independent_walk()
  w <- .ig_advisor_enumerate(rs)
  # Same outcome-status combinations, found by status code that shares nothing
  # with .ig_advisor_status() or the package enumerator
  combos <- unique(vapply(iw$records, function(r) paste(r$outcome, r$status, sep = "|"), ""))
  expect_setequal(combos, w$combos)
  expect_length(combos, 90L)
  # Path-level agreement through the session API on every 40th path (127
  # paths, all fixed; the complete comparison of all 5076 paths lives in
  # validation/independent/advisor/). The deep link is given in reverse order.
  idx <- seq(1L, length(iw$records), by = 40L)
  for (i in idx) {
    r <- iw$records[[i]]
    s <- advisor_replay(rev(r$answers))
    rec <- advisor_recommendation(s)
    lab <- paste(names(r$answers), r$answers, sep = "=", collapse = " ")
    expect_identical(rec$outcome_id, r$outcome, info = lab)
    expect_identical(rec$status, r$status, info = lab)
    expect_identical(as.character(unlist(rec$rules$status_modifiers_applied)), r$modifiers, info = lab)
    made <- Filter(function(a) identical(a$source, "dont_know"), rec$assumptions)
    expect_identical(vapply(made, `[[`, "", "question_id"), r$dk_assume, info = lab)
    expect_identical(vapply(rec$rule_trace, `[[`, "", "transition_id"), r$trace, info = lab)
    expect_identical(vapply(rec$answers, function(a) a$option_ids, ""), unname(r$answers), info = lab)
    expect_length(s$ignored_answers, 0L)
    has_method <- r$status %in% c("recommended", "conditional")
    expect_identical(is.null(rec$preferred), !has_method, info = lab)
    if (has_method) {
      expect_identical(rec$preferred$analysis_variant_id, rs$outcomes[[r$outcome]]$analysis_variant_id, info = lab)
      expect_identical(rec$missing_information, r$dk_assume, info = lab)
    }
  }
})

test_that("every enumerated outcome-status combination gives a consistent recommendation", {
  rs <- advisor_rules()
  w <- .ig_advisor_enumerate(rs)
  expect_length(w$examples, 90L)
  for (combo in names(w$examples)) {
    parts <- strsplit(combo, "|", fixed = TRUE)[[1]]
    s <- advisor_replay(w$examples[[combo]])
    expect_identical(s$state$outcome_id, parts[1], info = combo)
    rec <- advisor_recommendation(s)
    expect_identical(rec$status, parts[2], info = combo)
    expect_false(rec$planning$planning_compatibility == "exact_gpower_match", info = combo)
    if (rec$planning$planning_compatibility == "independently_validated_extension")
      expect_true(parts[1] %in% c("out.two_group.boschloo_central", "out.two_group.boschloo_one_sided"), info = combo)
    expect_length(rec$planning$spec_refs, 0L)
    expect_true(nzchar(rec$rationale$en) && nzchar(rec$rationale$de), info = combo)
    o <- rs$outcomes[[parts[1]]]
    if (rec$status %in% c("recommended", "conditional")) {
      vid <- o$analysis_variant_id
      expect_identical(rec$preferred$analysis_variant_id, vid, info = combo)
      expect_identical(rec$r_implementation$template_id, paste0("tpl.", vid, ".v1"), info = combo)
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
      expect_identical(rec$r_implementation$template_status, if (available) "available" else "template_not_yet_available", info = combo)
    } else {
      expect_null(rec$preferred)
      expect_null(rec$r_implementation$template_id)
      expect_true(length(rec$next_actions) > 0L, info = combo)
    }
    expect_identical(as.character(unlist(rec$rules$source_nodes$uzh)), as.character(unlist(o$source_nodes$uzh)))
    expect_identical(as.character(unlist(rec$rules$source_nodes$descrtab2)),
                     as.character(unlist(o$source_nodes$descrtab2)))
  }
})

test_that("every question and outcome has English and German texts", {
  rs <- advisor_rules()
  str_ok <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(trimws(x))
  for (q in rs$questions) {
    for (f in c("text_en", "text_de", "why_it_matters_en", "why_it_matters_de")) {
      expect_true(str_ok(q[[f]]), info = paste(q$question_id, f))
    }
    for (o in q$options) {
      expect_true(str_ok(o$label_en) && str_ok(o$label_de), info = paste(q$question_id, o$option_id))
      if (o$option_id == "dont_know") {
        expect_true(str_ok(o$explanation_en) && str_ok(o$explanation_de), info = q$question_id)
      }
      if (!is.null(o$assumption_en)) expect_true(str_ok(o$assumption_de), info = q$question_id)
    }
  }
  for (o in rs$outcomes) {
    id <- o$outcome_id
    for (f in c("title_en", "title_de", "rationale_en", "rationale_de")) expect_true(str_ok(o[[f]]), info = paste(id, f))
    if (!is.null(o$estimand_en)) expect_true(str_ok(o$estimand_de), info = id)
    for (a in o$assumptions) expect_true(str_ok(a$text_en) && str_ok(a$text_de), info = paste(id, a$assumption_id))
    for (a in o$alternatives) {
      expect_true(str_ok(a$why_not_preferred_en) && str_ok(a$why_not_preferred_de), info = paste(id, a$analysis_variant_id))
    }
    for (a in o$rejected_routes) expect_true(str_ok(a$reason_en) && str_ok(a$reason_de), info = id)
    for (m in o$status_modifiers) expect_true(str_ok(m$assumption_en) && str_ok(m$assumption_de), info = id)
    expect_true(str_ok(o$planning$note_en) && str_ok(o$planning$note_de), info = id)
  }
  # The generic *_en/*_de sibling check of the validator agrees
  for (q in rs$questions) expect_length(.ig_advisor_bilingual_problems(q, q$question_id), 0L)
  for (o in rs$outcomes) expect_length(.ig_advisor_bilingual_problems(o, o$outcome_id), 0L)
})

test_that("no outcome claims exact_gpower_match or an unvalidated extension", {
  rs <- advisor_rules()
  compat <- vapply(rs$outcomes, function(o) o$planning$planning_compatibility, "")
  expect_false(any(compat == "exact_gpower_match"))
  expect_setequal(names(compat)[compat == "independently_validated_extension"],
    c("out.two_group.boschloo_central", "out.two_group.boschloo_one_sided"))
  expect_true(all(compat %in% .ig_recommendation_vocab$planning_compatibility))
  avail <- vapply(rs$outcomes, function(o) o$planning$planning_availability, "")
  expect_true(all(avail %in% .ig_recommendation_vocab$planning_availability))
  # Classifier destinations never carry test-power planning
  cls <- Filter(function(o) identical(o$destination, "classification_workspace"), rs$outcomes)
  expect_true(length(cls) >= 5L)
  for (o in cls) {
    expect_identical(o$planning$planning_compatibility, "not_applicable")
    expect_identical(o$planning$planning_availability, "not_applicable")
  }
})

test_that("installed advisor data are identical to the audited specifications", {
  root <- th_pkg_root()
  skip_if(is.null(root), "source tree with specs/advisor not available")
  files <- c("rules.json", "analysis-variants.json", "uzh-tree.json", "descrtab2-tree.json",
             "source-crosswalk.csv", "acceptance-scenarios.json")
  for (f in files) {
    spec <- file.path(root, "specs", "advisor", f)
    inst <- system.file("advisor", f, package = "inferguidR")
    expect_true(file.exists(spec) && nzchar(inst), info = f)
    expect_identical(unname(tools::md5sum(inst)), unname(tools::md5sum(spec)), info = f)
  }
})

test_that("learning links are looked up from the registry, with explicit gaps", {
  rs <- advisor_rules()
  welch <- .ig_advisor_learning(rs, "t_means_two_independent_welch")
  expect_identical(welch$tutorial$status, "exact")
  expect_true(length(welch$tutorial$links) > 0L)
  ids <- vapply(rs$resources, `[[`, "", "resource_id")
  rank <- c(exact = 1L, related = 2L, broader = 3L)
  for (role in c("course", "tutorial", "overview")) {
    for (l in welch[[role]]$links) {
      expect_true(l$resource_id %in% ids)
      expect_identical(l$caution, length(l$known_issues) > 0L)
      # every link says which method id it matches, and `match` is the best of these
      mm <- vapply(l$method_matches, `[[`, "", "match")
      expect_true(all(vapply(l$method_matches, `[[`, "", "method_id") %in% unlist(welch$learning_method_ids)))
      expect_identical(l$match, names(rank)[min(rank[mm])])
    }
  }
  # A paired-t lesson is linked for Welch only through its exact Cohen's d content
  pt <- Filter(function(l) l$resource_id == "cttir_tutorials_dt_paired_t", welch$tutorial$links)
  expect_length(pt, 1L)
  expect_identical(pt[[1]]$method_matches, list(list(method_id = "effect_size_cohens_d", match = "exact")))
  # Boschloo is only taught as a related topic: never presented as exact
  bos <- .ig_advisor_learning(rs, "boschloo_two_sided_central")
  for (role in c("course", "tutorial")) expect_false(identical(bos[[role]]$status, "exact"))
  # A variant without learning method ids gets explicit gaps, nothing substituted
  bm <- .ig_advisor_learning(rs, "brunner_munzel")
  expect_length(bm$learning_method_ids, 0L)
  for (role in c("course", "tutorial", "overview")) {
    expect_identical(bm[[role]]$status, "gap")
    expect_length(bm[[role]]$links, 0L)
    expect_match(bm[[role]]$note_en, "^Gap: ")
    expect_match(bm[[role]]$note_de, "^L\u00fccke: ")
  }
  efa <- advisor_recommendation(advisor_replay(c(q.aim = "structure", q.stage = "planning",
                                                 q.structure_goal = "variables", q.reduction_purpose = "latent")))
  expect_identical(efa$preferred$learning$course$status, "gap")
  expect_length(efa$learning_links, 0L)
})

test_that("the variant used_by lists and learning ids resolve", {
  rs <- advisor_rules()
  lids <- unique(unlist(lapply(rs$resources, function(r) vapply(r$crosswalk, `[[`, "", "method_id"))))
  for (v in rs$variants) {
    expect_true(all(unlist(v$learning_method_ids) %in% lids), info = v$analysis_variant_id)
  }
  p <- .ig_advisor_content_problems(rs)
  expect_identical(nrow(p), 0L)
})
