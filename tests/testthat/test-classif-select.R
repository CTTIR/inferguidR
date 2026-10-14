# Cutoff selection: data-split labels, feasibility and tie reporting, held-out
# transfer, capacity and cost objectives, the selection basis and the exact
# handling of bounds.

se_ex <- function(part = "all") classification_examples("scored_synthetic", part = part, models = "model_a")

test_that("a cutoff chosen on the evaluated data is labelled in-sample and descriptive", {
  r <- select_cutoff(se_ex()$spec, "max_youden")
  expect_identical(r$evaluation$label, "in_sample_selection")
  expect_match(r$evaluation$statement, "not validated operational performance")
  expect_true(all(c("in_sample_selection", "fixed_cutoff_intervals", "post_selection_not_computed") %in%
                    r$limitations$code))
  iv <- classification_table(r, "intervals")
  expect_true(all(grepl("selected on these data", iv$note[iv$method == "wilson" & iv$status == "ok"])))
  sel <- classification_table(r, "selection")
  expect_identical(sel$procedure, "max_youden")
  expect_identical(sel$selection_data, "same_data")
  expect_identical(sel$method, "exhaustive search over all distinct decision sets")
  # a manual cutoff defaults to "same data" and says so; a predeclared rule does not get the caveat
  m <- evaluate_cutoff(se_ex()$spec)
  expect_identical(m$evaluation$label, "same_data_manual_cutoff")
  p <- evaluate_cutoff(update(se_ex()$spec, rule = cutoff_rule(cutoff = 0.3, selection_data = "predeclared_rule")))
  expect_identical(p$evaluation$label, "predeclared_rule")
  expect_false("post_selection_not_computed" %in% p$limitations$code)
})

test_that("a rule chosen on tuning data is evaluated as held out on test data", {
  tune <- se_ex("tuning")$spec
  test <- se_ex("test")$spec
  chosen <- select_cutoff(tune, "max_youden")
  held <- evaluate_cutoff(test, cutoff = chosen)
  expect_identical(held$evaluation$label, "held_out_rule")
  expect_match(held$evaluation$statement, "independent test or external data")
  expect_identical(held$evaluation$selection_data, "held_out_tuning")
  expect_equal(held$models[[1]]$operating_point$cutoff, chosen$models[[1]]$operating_point$cutoff)
  tr <- held$evaluation$transfer
  expect_identical(tr$source_spec_hash, chosen$spec_hash)
  expect_identical(tr$source_procedure, "max_youden")
  expect_false("in_sample_selection" %in% held$limitations$code)
  # the same data cannot pose as held-out data
  expect_identical(th_condition_code(evaluate_cutoff(tune, cutoff = chosen)), "invalid_input")
  # a different comparator is a different rule
  d_test <- se_ex("test")$data
  test_gt <- classification_spec(d_test, truth = "outcome", positive_label = "event", scores = "model_a",
                                 cutoff = 0.3, comparator = ">", score_scale = "probability")
  expect_identical(th_condition_code(evaluate_cutoff(test_gt, cutoff = chosen)), "invalid_input")
  # a theoretical cutoff (z scale) does not carry over to data
  th <- select_cutoff(classification_examples("binormal_equal_variance")$spec, "max_youden")
  expect_identical(th_condition_code(evaluate_cutoff(test, cutoff = th)), "invalid_input")
  # an infeasible source has no rule to carry over
  bad <- select_cutoff(tune, "constraint", constraints = list(tpr_min = 1, fpr_max = 0))
  expect_identical(th_condition_code(evaluate_cutoff(test, cutoff = bad)), "invalid_input")
})

test_that("conflicting constraints give an infeasible result with a per-bound report and no number", {
  r <- select_cutoff(se_ex()$spec, "constraint", constraints = list(tpr_min = 0.99, fpr_max = 0.01))
  expect_identical(r$status$code, "infeasible")
  m <- r$models[[1]]
  expect_null(m$operating_point)
  expect_null(m$metrics)
  expect_identical(m$selection$reason_code, "constraints_infeasible")
  b <- classification_table(r, "bounds")
  expect_setequal(b$constraint, c("tpr_min", "fpr_max"))
  expect_true(all(b$n_satisfied > 0))      # each bound alone is attainable
  expect_true(all(is.na(b$satisfied_at_choice)))
  expect_equal(b$best_value[b$constraint == "tpr_min"], 1)
  # curve summaries and the AUC interval remain, fixed-cutoff intervals do not
  iv <- classification_table(r, "intervals")
  expect_false(any(iv$method %in% c("wilson", "clopper_pearson")))
  expect_identical(iv$method, "delong")
  expect_true(any(grepl("no operating point satisfies", format(r))))
  # and the result serialises with its status
  back <- ig_from_json(ig_to_json(analysis_envelope(classification_evaluated_spec(r), r)))
  expect_identical(back$result$status$code, "infeasible")
})

test_that("ties are reported in full and broken among actual operating points", {
  d <- data.frame(y = c(1, 0, 1, 0), s = c(4, 3, 2, 1))
  spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "s", cutoff = 1)
  # Youden J = 0.5 at alert {4} and at alert {4, 3, 2}
  r1 <- select_cutoff(spec, "max_youden")
  r2 <- select_cutoff(spec, "max_youden", tie_break = "most_alerts")
  t1 <- classification_table(r1, "tied")
  expect_equal(nrow(t1), 2)
  expect_equal(t1$alerts, c(1, 3))
  expect_equal(r1$models[[1]]$operating_point$n_alerts, 1)
  expect_equal(r2$models[[1]]$operating_point$n_alerts, 3)
  expect_equal(r1$models[[1]]$selection$n_tied, 2)
  expect_true("tie_break" %in% classification_evaluated_spec(r1)$operating_rule$defaults_applied)
})

test_that("capacity and cost objectives choose the right decision sets", {
  spec <- se_ex()$spec
  cap <- select_cutoff(spec, "capacity", capacity = list(max_alerts = 25))
  op <- cap$models[[1]]$operating_point
  sw <- cutoff_sweep(spec, metrics = character())
  expect_equal(op$n_alerts, max(sw$alerts[sw$alerts <= 25]))
  cap2 <- select_cutoff(spec, "capacity", capacity = list(max_alert_rate = 0.1))
  expect_lte(cap2$models[[1]]$operating_point$n_alerts, 20)
  cost <- select_cutoff(spec, "cost_utility", costs = c(fn = 10, fp = 1, tp = -2))
  cv <- (10 * sw$fn + sw$fp - 2 * sw$tp) / 200
  expect_equal(cost$models[[1]]$selection$objective_value, min(cv))
  expect_identical(cost$models[[1]]$selection$objective_metric, "expected_cost_per_unit")
  expect_identical(cost$models[[1]]$selection$objective_direction, "min")
})

test_that("a manual cutoff checks declared constraints without enforcing them", {
  rule <- cutoff_rule(cutoff = 0.3, constraints = list(fpr_max = 0.05))
  r <- evaluate_cutoff(update(se_ex()$spec, rule = rule))
  expect_identical(r$status$code, "ok")
  sel <- classification_table(r, "selection")
  expect_false(sel$constraints_met)
  b <- classification_table(r, "bounds")
  expect_false(b$satisfied_at_choice)
})

test_that("the selection basis decides which prevalence prevalence-dependent objectives see", {
  spec <- update(se_ex()$spec, prevalence = 0.02)
  obs <- select_cutoff(spec, "max_f1")
  scen <- select_cutoff(spec, "max_f1", prevalence_basis = "scenario")
  expect_identical(obs$models[[1]]$selection$prevalence_basis, "observed")
  expect_identical(scen$models[[1]]$selection$prevalence_basis, "scenario")
  # at 2 % prevalence F1 favours fewer alerts than at the sample prevalence of 30 %
  expect_lt(scen$models[[1]]$operating_point$n_alerts, obs$models[[1]]$operating_point$n_alerts)
  # class-conditional objectives do not depend on the basis
  y1 <- select_cutoff(spec, "max_youden")
  y2 <- select_cutoff(spec, "max_youden", prevalence_basis = "scenario")
  expect_equal(y1$models[[1]]$operating_point$set_index, y2$models[[1]]$operating_point$set_index)
})

test_that("bounds are honoured exactly, also for very small or very large bounds", {
  th <- classification_spec(theory = list(auc = 0.8, sigma = 1), cutoff = 0, prevalence = 0.1)
  r <- select_cutoff(th, "constraint", constraints = list(fpr_max = 1e-20))
  op <- r$models[[1]]$operating_point
  expect_identical(op$cutoff_kind, "real")
  expect_equal(op$cutoff, qnorm(1e-20, lower.tail = FALSE), tolerance = 1e-9)
  expect_lte(op$fp / op$n_non_events, 1e-20 * (1 + 1e-12))
  expect_gt(op$tp, 0)                       # not the empty rule "alert none"
  cap <- select_cutoff(th, "capacity", capacity = list(max_alert_rate = 1e-15))
  expect_identical(cap$models[[1]]$operating_point$cutoff_kind, "real")
  expect_lte(cap$models[[1]]$operating_point$n_alerts / 1000, 1e-15 * (1 + 1e-12))
  # small objective values do not tie with zero
  expect_equal(r$models[[1]]$selection$n_tied, 1)
  # scored data: a bound equal to an attained rate is met
  d <- data.frame(y = c(1, 1, 1, 1, 0, 0, 0, 0, 0), s = 9:1)
  spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "s", cutoff = 5)
  r2 <- select_cutoff(spec, "constraint", constraints = list(tpr_min = 0.75), optimise = "tnr")
  expect_equal(r2$models[[1]]$operating_point$tp, 3)
})

test_that("select_cutoff() without an objective runs the rule in the specification", {
  spec <- classification_spec(se_ex()$data, truth = "outcome", positive_label = "event", scores = "model_a",
                              rule = cutoff_rule(selection = "max_mcc"))
  r <- select_cutoff(spec)
  expect_identical(r$models[[1]]$selection$procedure, "max_mcc")
  expect_identical(r$spec_hash, spec$spec_hash)
  expect_identical(th_condition_code(select_cutoff(se_ex()$spec)), "invalid_input")
  expect_identical(th_condition_code(select_cutoff(spec, constraints = list(tpr_min = 0.5))), "invalid_input")
  expect_identical(th_condition_code(select_cutoff(se_ex()$spec, "max_accuracy")), "invalid_input")
})

test_that("theoretical selections are model properties, not estimates", {
  th <- classification_examples("binormal_unequal_variance")$spec
  r <- select_cutoff(th, "max_youden")
  expect_identical(r$evaluation$label, "theoretical_model")
  expect_match(r$evaluation$statement, "property of the model")
  expect_false(grepl("observed", r$evaluation$statement))
  sp <- classification_evaluated_spec(r)
  expect_identical(sp$operating_rule$selection_data, "predeclared_rule")
  expect_true("selection_data_theoretical" %in% sp$operating_rule$defaults_applied)
  expect_length(r$models[[1]]$selection$analytic_stationary_points, 2)
})
