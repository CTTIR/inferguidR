# Plot-ready tables of the linked classification panels: every panel agrees
# with the result it comes from, reference lines satisfy their definitions,
# nothing is clamped onto a log axis and realised rules stay separate.

pd_result <- function(prevalence = NULL) {
  ex <- classification_examples("scored_synthetic")
  spec <- if (is.null(prevalence)) ex$spec else update(ex$spec, prevalence = prevalence)
  evaluate_cutoff(spec)
}

test_that("the confusion table and mosaic agree with the operating point", {
  r <- pd_result()
  op <- r$models[[1]]$operating_point
  ct <- classification_plot_data(r, "confusion_table")
  expect_equal(ct$data$cells$count, c(op$tp, op$fp, op$fn, op$tn))
  expect_equal(ct$data$totals$value, c(op$n_alerts, op$n_non_alerts, op$n_events, op$n_non_events, op$n_total))
  expect_true(all(ct$data$cells$expected_marker == ""))
  expect_equal(ct$data$tiles$denominator_value[ct$data$tiles$metric_id == "ppv"], op$n_alerts)
  mo <- classification_plot_data(r, "mosaic")
  joint <- mo$data$rectangles[mo$data$rectangles$mode == "joint", ]
  expect_equal(sum((joint$x1 - joint$x0) * (joint$y1 - joint$y0)), 1)
  expect_equal(joint$area_share, c(op$tp, op$fn, op$fp, op$tn) / op$n_total)
  within <- mo$data$rectangles[mo$data$rectangles$mode == "within_class", ]
  expect_equal(within$y1[within$cell == "tp"], op$tp / op$n_events)
  expect_equal(unique(within$x1[within$cell %in% c("tp", "fn")]), 0.5)
  # theory: expected counts carry the marker
  th <- evaluate_cutoff(classification_examples("binormal_equal_variance")$spec)
  ct2 <- classification_plot_data(th, "confusion_table")
  expect_true(all(ct2$data$cells$expected_marker == "\u2248"))
  expect_true(any(grepl("expected counts", ct2$notes)))
  # German labels
  de <- classification_plot_data(r, "confusion_table", lang = "de")
  expect_true("Alarm" %in% de$data$cells$row_label)
})

test_that("the scenario basis uses expected counts at the scenario prevalence and says so", {
  r <- pd_result(prevalence = 0.05)
  ct <- classification_plot_data(r, "confusion_table", basis = "scenario")
  sc <- r$models[[1]]$scenario$counts
  expect_equal(ct$data$cells$count, c(sc$tp, sc$fp, sc$fn, sc$tn))
  expect_true(all(ct$data$cells$count_kind == "expected"))
  expect_true(any(grepl("not a claim about", ct$notes)))
  expect_identical(th_condition_code(classification_plot_data(pd_result(), "confusion_table", basis = "scenario")),
                   "invalid_input")
})

test_that("the score distribution uses exact scores; bins are display only", {
  r <- pd_result()
  sd <- classification_plot_data(r, "score_distribution", bins = 20)
  sc <- sd$data$scores
  op <- r$models[[1]]$operating_point
  expect_equal(sum(sc$events), op$n_events)
  expect_equal(sum(sc$non_events), op$n_non_events)
  expect_equal(sum(sc$events[sc$alerted]), op$tp)
  expect_equal(sum(sc$non_events[sc$alerted]), op$fp)
  expect_true(all(sc$score[sc$alerted] >= 0.3))
  expect_equal(sum(sd$data$display_bins$events), op$n_events)
  expect_true(any(grepl("display only", sd$notes)))
  expect_true(any(grepl("calibration is not implied", sd$notes)))
  th <- classification_plot_data(evaluate_cutoff(classification_examples("binormal_unequal_variance")$spec),
                                 "score_distribution")
  expect_true(all(c("density_events", "density_non_events") %in% names(th$data$densities)))
  expect_true(any(grepl("not a probability", th$notes)))
})

test_that("ROC points are the decision sets; iso-PPV lines satisfy their definition", {
  r <- pd_result()
  roc <- classification_plot_data(r, "roc", ppv_levels = c(0.4, 0.6))
  sw <- cutoff_sweep(r, model_id = "model_a", metrics = c("tpr", "fpr"))
  expect_equal(roc$data$points$tpr, sw$tpr)
  expect_equal(roc$data$points$fpr, sw$fpr)
  op <- r$models[[1]]$operating_point
  expect_equal(roc$data$current$tpr, op$tp / op$n_events)
  iso <- roc$data$iso_ppv
  pi <- 0.3
  ok <- iso$fpr > 0
  ppv <- iso$tpr[ok] * pi / (iso$tpr[ok] * pi + iso$fpr[ok] * (1 - pi))
  expect_equal(ppv, iso$ppv_level[ok])
  expect_true(all(iso$prevalence_kind == "observed_sample"))
  expect_equal(roc$data$log_axis$n_fpr_zero, sum(sw$fpr == 0))
  expect_true(any(grepl("cannot be drawn on a log-FPR axis", roc$notes)))
  rules <- data.frame(id = "pooled", tp = 40, fn = 20, fp = 30, tn = 110)
  roc2 <- classification_plot_data(r, "roc", realised_rules = rules)
  expect_identical(roc2$data$realised_rules$marker, "realised_rule")
  expect_true(any(grepl("need not lie on the fixed-threshold curve", roc2$notes)))
  # a scenario prevalence moves the iso-PPV references, not the curve
  r5 <- pd_result(prevalence = 0.05)
  roc5 <- classification_plot_data(r5, "roc")
  expect_equal(roc5$data$points$tpr, roc$data$points$tpr)
  expect_true(all(roc5$data$iso_ppv$prevalence == 0.05))
  expect_true(all(roc5$data$iso_ppv$prevalence_kind == "scenario_assumed"))
})

test_that("PR points are attained operating points only; iso-F1 curves satisfy their definition", {
  r <- pd_result()
  pr <- classification_plot_data(r, "pr")
  sw <- cutoff_sweep(r, model_id = "model_a", metrics = character())
  expect_equal(nrow(pr$data$points), sum(sw$alerts > 0))
  expect_equal(pr$data$excluded$n_zero_alert_points, sum(sw$alerts == 0))
  expect_false(any(is.na(pr$data$points$precision)))
  f <- pr$data$iso_f1
  expect_equal(2 * f$precision * f$recall / (f$precision + f$recall), f$f1_level)
  expect_equal(pr$data$baseline$precision, 0.3)
  expect_true(any(grepl("not joined by straight lines", pr$notes)))
  expect_true(any(grepl("different quantities", pr$notes)))
  # on the scenario basis precision is recomputed and realised rules are withheld
  r5 <- pd_result(prevalence = 0.05)
  rules <- data.frame(id = "pooled", tp = 40, fn = 20, fp = 30, tn = 110)
  pr5 <- classification_plot_data(r5, "pr", basis = "scenario", realised_rules = rules)
  expect_equal(pr5$data$baseline$precision, 0.05)
  expect_null(pr5$data$realised_rules)
  expect_equal(pr5$data$points$recall, pr$data$points$recall)
})

test_that("metric families cover the full range and count what a log axis cannot show", {
  r <- pd_result()
  fam <- classification_plot_data(r, "metric_vs_cutoff")
  v <- fam$data$values
  expect_setequal(unique(v$family), c("within_truth", "within_alerts", "summary", "burden", "likelihood_ratios"))
  expect_equal(fam$data$domains$lower[fam$data$domains$family == "summary"], -1)
  lr <- v[v$family == "likelihood_ratios", ]
  expect_equal(fam$data$log_axis$n_positive_infinity, sum(lr$status == "positive_infinity"))
  expect_gt(fam$data$log_axis$n_positive_infinity, 0)
  expect_false(any(v$on_cutoff_axis[v$cutoff_kind == "alert_none"]))
  expect_true(all(is.na(v$value[v$status != "ok"])))
})

test_that("the prevalence curve reproduces the observed predictive values at the sample prevalence", {
  r <- pd_result()
  pc <- classification_plot_data(r, "prevalence_curve", ppv_target = 0.5)
  m <- classification_table(r, "metrics", model_id = "model_a")
  at <- pc$data$curve[abs(pc$data$curve$prevalence - 0.3) < 1e-12, ]
  expect_equal(at$ppv, m$value[m$metric_id == "ppv"])
  expect_equal(at$npv, m$value[m$metric_id == "npv"])
  expect_identical(pc$data$markers$kind, "observed_sample")
  expect_identical(pc$data$target$status, "ok")
  expect_true(any(grepl("no sampling uncertainty", pc$notes)))
  expect_null(classification_plot_data(r, "prevalence_curve")$data$target)
})

test_that("the metric pair, network, sweep and current tables come from the same result", {
  r <- pd_result(prevalence = 0.05)
  pair <- classification_plot_data(r, "metric_pair", metric_x = "fp", metric_y = "tp")
  expect_true(all(pair$data$axes$prevalence_dependent))
  rules <- data.frame(id = "pooled", tp = 40, fn = 20, fp = 30, tn = 110)
  hidden <- classification_plot_data(r, "metric_pair", basis = "scenario", metric_x = "fp", metric_y = "tp",
                                     realised_rules = rules)
  expect_null(hidden$data$realised_rules)
  expect_true(any(grepl("realised rules are hidden", hidden$notes)))
  free <- classification_plot_data(r, "metric_pair", basis = "scenario", metric_x = "fpr", metric_y = "tpr",
                                   realised_rules = rules)
  expect_identical(free$data$realised_rules$marker, "realised_rule")
  expect_identical(th_condition_code(classification_plot_data(r, "metric_pair", metric_x = "roc_auc")), "invalid_input")
  net <- classification_plot_data(r, "metric_network")
  expect_setequal(net$data$nodes$id, metric_dictionary(scope = "operating_point")$id)
  bayes <- net$data$edges[net$data$edges$kind == "bayes", ]
  expect_setequal(unique(bayes$to), c("ppv", "npv"))
  sw <- classification_plot_data(r, "sweep_table")
  expect_equal(sum(sw$data$sweep$selected), 1)
  cur <- classification_plot_data(r, "current_table", lang = "de")
  expect_true(all(c("short", "formula", "denominator") %in% names(cur$data$metrics)))
  expect_match(cur$data$metrics$denominator[cur$data$metrics$metric_id == "ppv"], "Alarme")
})

test_that("panels refuse what the result cannot support", {
  r <- pd_result()
  expect_identical(th_condition_code(classification_plot_data(r, "roc", model_id = "nope")), "invalid_input")
  expect_identical(th_condition_code(classification_plot_data(r, "roc", ppv_levels = 1)), "invalid_input")
  no_curves <- evaluate_cutoff(classification_examples("scored_synthetic")$spec, include_curves = FALSE)
  expect_identical(th_condition_code(classification_plot_data(no_curves, "roc")), "unsupported")
  inf <- select_cutoff(classification_examples("scored_synthetic")$spec, "constraint",
                       constraints = list(tpr_min = 1, fpr_max = 0))
  expect_identical(th_condition_code(classification_plot_data(inf, "confusion_table")), "infeasible")
  # the ROC curve of an infeasible selection is still available (it does not depend on a cutoff)
  expect_null(classification_plot_data(inf, "roc")$data$current)
  expect_identical(th_condition_code(classification_plot_data(list(), "roc")), "invalid_input")
})

test_that("the model comparison panel overlays matched models with distinct markers", {
  r <- pd_result()
  cmp <- classification_plot_data(r, "model_comparison", model_id = "model_b")
  expect_setequal(unique(cmp$data$curves$model_id), c("model_a", "model_b"))
  op <- cmp$data$operating_points
  expect_identical(op$marker[op$model_id == "model_b"], "current_operating_point")
  expect_identical(op$marker[op$model_id == "model_a"], "comparison_operating_point")
  expect_equal(cmp$data$comparisons$difference, compare_classifiers(classification_evaluated_spec(r))$difference)
  expect_true(any(grepl("not an equivalent operating point", cmp$notes)))
  expect_true(any(grepl("same 200 units", cmp$notes)))
  one <- evaluate_cutoff(classification_examples("scored_synthetic", models = "model_a")$spec)
  expect_identical(th_condition_code(classification_plot_data(one, "model_comparison")), "unsupported")
})
