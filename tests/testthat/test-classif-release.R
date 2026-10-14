# Release tests of CLASSIFICATION_CUTOFF ("Release tests" section): extreme
# predictions, missing classes, ties, reversed scores, zero cells, perfect and
# reversed classification, imbalance, constant scores, weights and clusters,
# probability inputs, aggregate-only data, unequal-variance theory and the
# assumptions behind invariance claims.

rl_spec <- function(y, s, cutoff = 0, ...) {
  classification_spec(data.frame(y = y, s = s), truth = "y", positive_label = 1, scores = "s",
                      cutoff = cutoff, ...)
}
rl_metric <- function(res, id, what = "value") {
  m <- classification_table(res, "metrics")
  m[[what]][m$metric_id == id]
}
rl_summary <- function(res, id, what = "value") {
  s <- classification_table(res, "summaries")
  s[[what]][s$metric_id == id][1]
}

test_that("alert-none and alert-all rules give typed predictive values, never a fake PPV", {
  y <- c(1, 1, 0, 0, 0, 1, 0)
  s <- c(5, 4, 4, 3, 2, 2, 1)
  none <- evaluate_cutoff(rl_spec(y, s, cutoff = 6))
  op <- none$models[[1]]$operating_point
  expect_equal(c(op$tp, op$fp, op$n_alerts), c(0, 0, 0))
  expect_identical(rl_metric(none, "ppv", "status"), "undefined")
  expect_identical(rl_metric(none, "ppv", "reason_code"), "no_alerts")
  expect_true(is.na(rl_metric(none, "ppv")))
  expect_identical(rl_metric(none, "fdr", "status"), "undefined")
  expect_identical(rl_metric(none, "lr_pos", "reason_code"), "tpr_and_fpr_zero")
  expect_identical(rl_metric(none, "f1"), 0)
  expect_identical(rl_metric(none, "mcc", "status"), "undefined")
  allr <- evaluate_cutoff(rl_spec(y, s, cutoff = 1))
  op <- allr$models[[1]]$operating_point
  expect_equal(c(op$tp, op$fn, op$fp, op$tn), c(3, 0, 4, 0))
  expect_identical(rl_metric(allr, "npv", "reason_code"), "no_non_alerts")
  expect_identical(rl_metric(allr, "for_rate", "status"), "undefined")
  expect_equal(rl_metric(allr, "ppv"), 3 / 7)
  expect_equal(rl_metric(allr, "ppv"), rl_metric(allr, "prevalence"))
  expect_identical(rl_metric(allr, "lr_neg", "reason_code"), "fnr_and_tnr_zero")
  # both endpoints are enumerated decision sets
  sw <- cutoff_sweep(rl_spec(y, s), metrics = "ppv")
  expect_identical(sw$cutoff_kind[1], "alert_none")
  expect_identical(sw$ppv_status[1], "undefined")
  expect_equal(sw$alerts[nrow(sw)], 7)
})

test_that("a nominal cutoff equal to the largest score still alerts that score under >=", {
  r <- evaluate_cutoff(rl_spec(c(1, 0, 1), c(1, 0.5, 0.2), cutoff = 1, score_scale = "probability"))
  expect_equal(r$models[[1]]$operating_point$tp, 1)
  r2 <- evaluate_cutoff(rl_spec(c(1, 0, 1), c(1, 0.5, 0.2), cutoff = 1, comparator = ">", score_scale = "probability"))
  expect_equal(r2$models[[1]]$operating_point$n_alerts, 0)
})

test_that("a missing actual class makes class-dependent quantities unavailable with reasons", {
  r <- evaluate_cutoff(rl_spec(c(0, 0, 0, 0), c(0.1, 0.5, 0.7, 0.9), cutoff = 0.6))
  expect_identical(rl_metric(r, "tpr", "reason_code"), "no_events")
  expect_identical(rl_metric(r, "youden_j", "status"), "undefined")
  expect_equal(rl_metric(r, "fpr"), 0.5)
  expect_identical(rl_summary(r, "roc_auc", "status"), "undefined")
  expect_identical(rl_summary(r, "average_precision", "reason_code"), "no_events")
  iv <- classification_table(r, "intervals")
  expect_identical(iv$status[iv$metric_id == "roc_auc"], "undefined")
  expect_identical(iv$reason_code[iv$metric_id == "tpr" & iv$method == "wilson"], "no_events")
  expect_identical(r$models[[1]]$views$pr$reason_code, "no_events")
  expect_identical(th_condition_code(classification_plot_data(r, "roc")), "unsupported")
  sel <- select_cutoff(rl_spec(c(0, 0, 0, 0), c(0.1, 0.5, 0.7, 0.9)), "max_youden")
  expect_identical(sel$status$code, "infeasible")
  expect_identical(sel$models[[1]]$selection$reason_code, "objective_undefined_everywhere")
  expect_null(sel$models[[1]]$operating_point)
  # only events: specificity side undefined
  r1 <- evaluate_cutoff(rl_spec(c(1, 1, 1), c(1, 2, 3), cutoff = 2))
  expect_identical(rl_metric(r1, "tnr", "reason_code"), "no_non_events")
  expect_identical(rl_summary(r1, "roc_auc", "reason_code"), "no_non_events")
})

test_that("tied scores move together and ties count one half in the AUC", {
  y <- c(1, 0, 1, 0, 1, 0)
  s <- c(3, 3, 2, 2, 1, 1)
  spec <- rl_spec(y, s)
  sw <- cutoff_sweep(spec, metrics = character())
  expect_equal(nrow(sw), 4)                 # 3 distinct scores + alert none
  expect_equal(sw$alerts, c(0, 2, 4, 6))
  r <- evaluate_cutoff(update(spec, cutoff = 2))
  expect_equal(r$models[[1]]$operating_point$n_alerts, 4)
  r_gt <- evaluate_cutoff(rl_spec(y, s, cutoff = 2, comparator = ">"))
  expect_equal(r_gt$models[[1]]$operating_point$n_alerts, 2)
  expect_equal(rl_summary(r, "roc_auc"), 0.5)
  pd <- classification_plot_data(r, "roc")
  expect_true(all(pd$data$points$segment_from_previous[-1] == "tie_diagonal"))
})

test_that("reversing the declared direction mirrors counts and the AUC, without automatic flipping", {
  set.seed(4)
  y <- rbinom(40, 1, 0.4)
  s <- round(rnorm(40, y), 1)
  up <- evaluate_cutoff(rl_spec(y, s, cutoff = 0.3))
  down <- evaluate_cutoff(rl_spec(y, s, cutoff = 0.3, direction = "lower_is_positive"))
  neg <- evaluate_cutoff(rl_spec(y, -s, cutoff = -0.3, comparator = ">="))
  expect_equal(rl_summary(up, "roc_auc") + rl_summary(down, "roc_auc"), 1)
  # "alert if s <= 0.3" equals "alert if -s >= -0.3"
  expect_equal(unlist(down$models[[1]]$operating_point[c("tp", "fn", "fp", "tn")]),
               unlist(neg$models[[1]]$operating_point[c("tp", "fn", "fp", "tn")]))
  expect_identical(down$models[[1]]$operating_point$rule_text, "alert if score <= 0.3")
  expect_identical(down$models[[1]]$comparator, "<=")
})

test_that("zero cells give infinity and undefined statuses that survive serialisation", {
  r <- evaluate_cutoff(rl_spec(c(1, 1, 0, 0, 0), c(5, 4, 3, 2, 1), cutoff = 4))
  expect_identical(rl_metric(r, "lr_pos", "status"), "positive_infinity")
  expect_identical(rl_metric(r, "dor", "status"), "positive_infinity")
  expect_identical(rl_metric(r, "lr_neg"), 0)
  back <- ig_from_json(ig_to_json(analysis_envelope(classification_evaluated_spec(r), r)))
  m <- classification_table(back$result, "metrics")
  expect_identical(m$status[m$metric_id == "lr_pos"], "positive_infinity")
  expect_true(is.na(m$value[m$metric_id == "lr_pos"]))
  expect_false(grepl("Infinity|NaN", ig_to_json(r)))
})

test_that("perfect and reversed classification are reported honestly", {
  y <- c(1, 1, 1, 0, 0, 0)
  s <- c(6, 5, 4, 3, 2, 1)
  perfect <- evaluate_cutoff(rl_spec(y, s, cutoff = 4))
  expect_equal(rl_summary(perfect, "roc_auc"), 1)
  expect_equal(rl_metric(perfect, "mcc"), 1)
  iv <- classification_table(perfect, "intervals")
  expect_identical(iv$reason_code[iv$metric_id == "roc_auc"], "delong_degenerate")
  expect_equal(iv$upper[iv$metric_id == "tpr" & iv$method == "wilson"], 1)
  best <- select_cutoff(rl_spec(y, s), "max_youden")
  expect_equal(best$models[[1]]$selection$objective_value, 1)
  expect_equal(best$models[[1]]$operating_point$cutoff, 4)
  reversed <- evaluate_cutoff(rl_spec(y, -s, cutoff = -4))
  expect_equal(rl_summary(reversed, "roc_auc"), 0)
  # with the declared direction the best Youden J of a reversed score is 0, at both endpoints
  sel <- select_cutoff(rl_spec(y, -s), "max_youden")
  expect_equal(sel$models[[1]]$selection$objective_value, 0)
  expect_equal(sel$models[[1]]$selection$n_tied, 2)
  expect_identical(sel$models[[1]]$operating_point$cutoff_kind, "alert_none")
})

test_that("strong imbalance keeps exact counts and valid intervals", {
  y <- c(rep(1, 3), rep(0, 1997))
  s <- c(0.9, 0.8, 0.2, seq(0, 0.95, length.out = 1997))
  r <- evaluate_cutoff(rl_spec(y, s, cutoff = 0.85))
  expect_equal(rl_metric(r, "prevalence"), 0.0015)
  op <- r$models[[1]]$operating_point
  expect_equal(op$tp, 1)
  expect_equal(op$fp, sum(seq(0, 0.95, length.out = 1997) >= 0.85))
  iv <- classification_table(r, "intervals")
  w <- iv[iv$metric_id == "tpr" & iv$method == "wilson", ]
  expect_true(w$lower > 0 && w$upper < 1 && w$lower < 1 / 3 && w$upper > 1 / 3)
})

test_that("constant scores give two decision sets, a diagonal ROC and AUC 0.5", {
  r <- evaluate_cutoff(rl_spec(c(1, 0, 1, 0, 0), rep(3, 5), cutoff = 3))
  sw <- cutoff_sweep(r, metrics = c("tpr", "fpr"))
  expect_equal(nrow(sw), 2)
  expect_equal(sw$tpr, c(0, 1))
  expect_equal(sw$fpr, c(0, 1))
  expect_equal(rl_summary(r, "roc_auc"), 0.5)
  iv <- classification_table(r, "intervals")
  expect_identical(iv$reason_code[iv$metric_id == "roc_auc"], "delong_degenerate")
})

test_that("clustered units make independence-based intervals unsupported and enable the cluster bootstrap", {
  set.seed(9)
  d <- data.frame(y = rbinom(48, 1, 0.4), s = round(rnorm(48), 1), pid = rep(1:12, each = 4))
  spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "s", cutoff = 0, cluster = "pid",
                              uncertainty = classification_uncertainty(bootstrap_replicates = 200, seed = 11))
  r <- evaluate_cutoff(spec)
  iv <- classification_table(r, "intervals")
  expect_true(all(iv$status[iv$method %in% c("wilson", "clopper_pearson", "delong")] == "unsupported"))
  expect_true(all(iv$reason_code[iv$method == "delong"] == "clustered_units"))
  boot <- iv[iv$method == "bootstrap_percentile", ]
  expect_true(all(boot$status == "ok"))
  expect_match(boot$note[1], "whole clusters resampled")
  expect_true("clustered_units" %in% r$limitations$code)
  d$s2 <- d$s + round(rnorm(48, 0, 0.5), 1)
  cmp <- compare_classifiers(classification_spec(d, truth = "y", positive_label = 1, scores = c("s", "s2"),
                                                 cutoff = 0, cluster = "pid"))
  expect_identical(cmp$inference_status, "unsupported")
  expect_true(is.finite(cmp$difference))
  expect_true(is.na(cmp$p_value))
})

test_that("probability views need declared probabilities; endpoints make calibration fits undefined", {
  raw <- evaluate_cutoff(rl_spec(c(1, 0, 1, 0), c(3, 1, 2, 2), cutoff = 2))
  expect_identical(rl_summary(raw, "brier_score", "reason_code"), "not_probability_scale")
  expect_identical(raw$models[[1]]$views$calibration$reason_code, "not_probability_scale")
  expect_identical(th_condition_code(net_benefit(rl_spec(c(1, 0), c(3, 1)), 0.2)), "unsupported")
  expect_identical(th_condition_code(evaluate_cutoff(classification_spec(counts = c(tp = 1, fn = 1, fp = 1, tn = 1)),
                                                     threshold_probability = 0.2)), "invalid_input")
  ends <- evaluate_cutoff(rl_spec(c(1, 0, 1, 0, 1), c(1, 0, 0.7, 0.2, 0.5), cutoff = 0.5, score_scale = "probability"))
  expect_identical(rl_summary(ends, "calibration_slope", "reason_code"), "probability_endpoints")
  expect_equal(rl_summary(ends, "brier_score"), mean((c(1, 0, 0.7, 0.2, 0.5) - c(1, 0, 1, 0, 1))^2))
  low <- rl_spec(c(1, 0, 1), c(0.2, 0.8, 0.3), score_scale = "probability", direction = "lower_is_positive", cutoff = 0.5)
  expect_identical(th_condition_code(net_benefit(low, 0.2)), "unsupported")
  expect_identical(th_condition_code(net_benefit(rl_spec(c(1, 0), c(0.9, 0.1), score_scale = "probability"))),
                   "invalid_input")
})

test_that("aggregate-only data give one operating point and disable every curve view", {
  spec <- classification_spec(counts = c(tp = 30, fn = 10, fp = 15, tn = 145), cutoff = 0.4)
  r <- evaluate_cutoff(spec)
  expect_identical(r$status$code, "ok")
  expect_true(all(vapply(r$models[[1]]$views, function(v) identical(v$reason_code, "aggregate_counts_only"), TRUE)))
  sm <- classification_table(r, "summaries")
  expect_true(all(sm$status == "not_applicable"))
  expect_true(all(is.na(sm$value)))
  expect_identical(th_condition_code(cutoff_sweep(r)), "unsupported")
  expect_identical(th_condition_code(cutoff_sweep(spec)), "unsupported")
  expect_identical(th_condition_code(classification_auc(spec)), "unsupported")
  expect_identical(th_condition_code(pr_summary(spec)), "unsupported")
  expect_identical(th_condition_code(select_cutoff(spec, "max_youden")), "invalid_input")
  expect_identical(th_condition_code(evaluate_cutoff(spec, cutoff = 0.5)), "invalid_input")
  for (p in c("roc", "pr", "metric_vs_cutoff", "sweep_table", "score_distribution", "metric_pair")) {
    expect_identical(th_condition_code(classification_plot_data(r, p)), "unsupported", info = p)
  }
  ct <- classification_plot_data(r, "confusion_table")
  expect_equal(ct$data$cells$count, c(30, 15, 10, 145))
  expect_true("aggregate_counts_only" %in% r$limitations$code)
})

test_that("unequal-variance binormal theory: analytic rates, fixed cutoff, N scales counts only", {
  ex <- classification_examples("binormal_unequal_variance")
  r <- evaluate_cutoff(ex$spec)
  m <- binormal_model(0.8, 1.5)
  op <- r$models[[1]]$operating_point
  expect_equal(op$tp, 100 * pnorm((1 - m$d) / 1.5, lower.tail = FALSE), tolerance = 1e-12)
  expect_equal(op$fp, 900 * pnorm(1, lower.tail = FALSE), tolerance = 1e-12)
  expect_identical(op$count_kind, "expected")
  expect_match(op$cutoff_scale, "standardised to the non-event distribution")
  # changing the prevalence keeps TPR, FPR and AUC (fixed class-conditional distributions) but moves PPV
  r2 <- evaluate_cutoff(update(ex$spec, prevalence = 0.4))
  expect_equal(rl_metric(r2, "tpr"), rl_metric(r, "tpr"))
  expect_equal(rl_metric(r2, "fpr"), rl_metric(r, "fpr"))
  expect_equal(rl_summary(r2, "roc_auc"), rl_summary(r, "roc_auc"))
  expect_gt(rl_metric(r2, "ppv"), rl_metric(r, "ppv"))
  expect_false(any(classification_table(r, "intervals")$status == "ok"))
})

test_that("selecting a cutoff changes operating-point metrics but not the ROC AUC", {
  ex <- classification_examples("scored_synthetic", models = "model_a")
  a <- evaluate_cutoff(ex$spec, cutoff = 0.2)
  b <- evaluate_cutoff(ex$spec, cutoff = 0.6)
  c <- select_cutoff(ex$spec, "max_f1")
  expect_equal(rl_summary(a, "roc_auc"), rl_summary(b, "roc_auc"))
  expect_equal(rl_summary(a, "roc_auc"), rl_summary(c, "roc_auc"))
  expect_false(isTRUE(all.equal(rl_metric(a, "tpr"), rl_metric(b, "tpr"))))
})

test_that("precision need not rise with the cutoff, and the PPV constraint uses actual operating points", {
  y <- c(0, 1, 1, 0, 1, 0, 0, 0)
  s <- c(8, 7, 6, 5, 4, 3, 2, 1)
  sw <- cutoff_sweep(rl_spec(y, s), metrics = "ppv")
  p <- sw$ppv[sw$ppv_status == "ok"]
  expect_true(any(diff(p) > 0) && any(diff(p) < 0))
  sel <- select_cutoff(rl_spec(y, s), "constraint", constraints = list(ppv_min = 0.6))
  op <- sel$models[[1]]$operating_point
  expect_gte(op$tp / op$n_alerts, 0.6)
  ok <- sw$ppv_status == "ok" & sw$ppv >= 0.6
  expect_equal(op$tp, max(sw$tp[ok]))
})

test_that("prevalence invariance holds only as a declared scenario; case-control samples are flagged", {
  ex <- classification_examples("scored_synthetic", models = "model_a")
  s <- update(ex$spec, prevalence = 0.02, prevalence_source = "external", sampling = "case_control")
  r <- evaluate_cutoff(s)
  sc <- classification_table(r, "scenario")
  ob <- classification_table(r, "metrics")
  expect_equal(sc$value[sc$metric_id == "tpr"], ob$value[ob$metric_id == "tpr"])
  expect_equal(sc$value[sc$metric_id == "prevalence"], 0.02)
  expect_false(isTRUE(all.equal(sc$value[sc$metric_id == "ppv"], ob$value[ob$metric_id == "ppv"])))
  expect_identical(unique(sc$estimate_kind), "expected_scenario")
  expect_identical(r$models[[1]]$scenario$source, "external")
  expect_match(r$models[[1]]$scenario$assumption, "not a claim about")
  iv <- classification_table(r, "intervals")
  expect_identical(iv$reason_code[iv$metric_id == "ppv" & iv$method == "wilson"], "case_control_sampling")
  expect_identical(iv$status[iv$metric_id == "tpr" & iv$method == "wilson"], "ok")
  expect_match(ob$note[ob$metric_id == "ppv"], "sample-specific")
  expect_true(all(c("case_control_sampling", "prevalence_scenario") %in% r$limitations$code))
})
