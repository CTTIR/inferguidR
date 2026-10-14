# Input contract of classification_spec(): kinds of input, labels, direction,
# comparator, weights, clusters, missing values, metadata and immutability.

sp_data <- function() {
  data.frame(id = sprintf("u%02d", 1:10),
             y = c("yes", "no", "yes", "no", "no", "yes", "no", "no", "yes", "no"),
             s = c(0.9, 0.2, 0.7, 0.4, 0.4, 0.6, 0.1, 0.3, 0.4, 0.8),
             s2 = c(3, 1, 2, 2, 1, 3, 0, 1, 2, 2),
             stringsAsFactors = FALSE)
}

test_that("exactly one kind of input is accepted", {
  d <- sp_data()
  expect_identical(th_condition_code(classification_spec(cutoff = 1)), "invalid_input")
  expect_identical(th_condition_code(classification_spec(d, truth = "y", positive_label = "yes", scores = "s",
                                                         counts = c(tp = 1, fn = 1, fp = 1, tn = 1), cutoff = 0.5)),
                   "invalid_input")
  expect_identical(th_condition_code(classification_spec(counts = c(tp = 1, fn = 1, fp = 1, tn = 1),
                                                         theory = list(auc = 0.8, sigma = 1), prevalence = 0.1)),
                   "invalid_input")
  s <- classification_spec(d, truth = "y", positive_label = "yes", scores = "s", cutoff = 0.5)
  expect_s3_class(s, "classification_spec")
  expect_identical(s$input_kind, "scored_data")
  expect_identical(classification_spec(counts = c(tp = 1, fn = 1, fp = 1, tn = 1))$input_kind, "counts_2x2")
})

test_that("the positive class is explicit, never inferred, and labels are validated", {
  d <- sp_data()
  expect_identical(th_condition_code(classification_spec(d, truth = "y", scores = "s", cutoff = 0.5)), "invalid_input")
  expect_identical(th_condition_code(classification_spec(d, truth = "y", positive_label = "maybe", scores = "s",
                                                         cutoff = 0.5)), "invalid_input")
  three <- d
  three$y[1] <- "unclear"
  expect_identical(th_condition_code(classification_spec(three, truth = "y", positive_label = "yes", scores = "s",
                                                         cutoff = 0.5)), "invalid_input")
  expect_identical(th_condition_code(classification_spec(d, truth = "y", positive_label = "yes", negative_label = "nein",
                                                         scores = "s", cutoff = 0.5)), "invalid_input")
  expect_identical(th_condition_code(classification_spec(d, truth = "y", positive_label = "yes", negative_label = "yes",
                                                         scores = "s", cutoff = 0.5)), "invalid_input")
  # the choice of positive label changes which units are events
  a <- evaluate_cutoff(classification_spec(d, truth = "y", positive_label = "yes", scores = "s", cutoff = 0.5))
  b <- evaluate_cutoff(classification_spec(d, truth = "y", positive_label = "no", scores = "s", cutoff = 0.5,
                                           direction = "lower_is_positive"))
  expect_equal(a$models[[1]]$operating_point$n_events, 4)
  expect_equal(b$models[[1]]$operating_point$n_events, 6)
  expect_identical(classification_spec(d, truth = "y", positive_label = "yes", scores = "s", cutoff = 0.5)$models[[1]]$negative_label, "no")
  # factor and logical truth columns
  f <- d
  f$y <- factor(f$y)
  expect_identical(classification_spec(f, truth = "y", positive_label = "yes", scores = "s", cutoff = 0.5)$models[[1]]$positive_label, "yes")
  l <- data.frame(y = c(TRUE, FALSE, TRUE), s = 1:3)
  expect_identical(classification_spec(l, truth = "y", positive_label = TRUE, scores = "s", cutoff = 2)$models[[1]]$negative_label, "FALSE")
})

test_that("a single observed class is accepted but recorded", {
  d <- data.frame(y = c(0, 0, 0), s = c(1, 2, 3))
  s <- classification_spec(d, truth = "y", positive_label = 1, scores = "s", cutoff = 2)
  r <- evaluate_cutoff(s)
  expect_true("single_class" %in% r$limitations$code)
  expect_identical(r$models[[1]]$views$roc$available, FALSE)
  expect_identical(r$models[[1]]$views$roc$reason_code, "no_events")
})

test_that("the default direction is visible as a default and never chosen from the data", {
  d <- sp_data()
  s <- classification_spec(d, truth = "y", positive_label = "yes", scores = "s", cutoff = 0.5)
  m <- s$models[[1]]
  expect_identical(m$direction, "higher_is_positive")
  expect_identical(m$comparator, ">=")
  expect_true(all(c("direction", "comparator", "evaluation_split", "score_meaning") %in% m$defaults_applied))
  expect_match(paste(format(s), collapse = "\n"), "defaults applied \\(confirm\\)")
  # a score that ranks events low keeps its declared direction: AUC below 0.5 is reported as is
  rev <- d
  rev$s <- -rev$s
  r <- evaluate_cutoff(classification_spec(rev, truth = "y", positive_label = "yes", scores = "s", cutoff = -0.5))
  auc <- classification_table(r, "summaries")
  expect_lt(auc$value[auc$metric_id == "roc_auc"], 0.5)
  expect_identical(r$models[[1]]$direction, "higher_is_positive")
  # an explicit comparator implies the direction and is recorded as such
  s2 <- classification_spec(d, truth = "y", positive_label = "yes", scores = "s", cutoff = 0.5, comparator = "<")
  expect_identical(s2$models[[1]]$direction, "lower_is_positive")
  expect_true("direction_implied_by_comparator" %in% s2$models[[1]]$defaults_applied)
  expect_identical(th_condition_code(classification_spec(d, truth = "y", positive_label = "yes", scores = "s",
                                                         cutoff = 0.5, direction = "higher_is_positive",
                                                         comparator = "<=")), "invalid_input")
  expect_identical(th_condition_code(classification_spec(d, truth = "y", positive_label = "yes", scores = "s",
                                                         cutoff = 0.5, comparator = "=>")), "invalid_input")
})

test_that("per-model directions and metadata are matched by model id", {
  d <- sp_data()
  meta <- list(classifier_model("b", name = "Model B", version = "2", score_meaning = "points"),
               classifier_model("a", name = "Model A", population = "synthetic", unit = "visit",
                                source_hash = strrep("a", 64)))
  s <- classification_spec(d, truth = "y", positive_label = "yes", scores = c(a = "s", b = "s2"),
                           models = meta, cutoff = 1, direction = c(a = "higher_is_positive", b = "lower_is_positive"))
  expect_identical(vapply(s$models, `[[`, "", "model_id"), c("a", "b"))
  expect_identical(s$models[[1]]$name, "Model A")
  expect_identical(s$models[[1]]$score_column, "s")
  expect_identical(s$models[[2]]$score_meaning, "points")
  expect_identical(s$models[[2]]$direction, "lower_is_positive")
  expect_identical(s$models[[1]]$unit, "visit")
  expect_identical(th_condition_code(classification_spec(d, truth = "y", positive_label = "yes", scores = c("s", "s2"),
                                                         cutoff = 1, direction = c(s = "higher_is_positive"))),
                   "invalid_input")
  expect_identical(th_condition_code(classifier_model("bad id")), "invalid_input")
  expect_identical(th_condition_code(classifier_model("m", source_hash = "abc")), "invalid_input")
  expect_identical(th_condition_code(classification_spec(d, truth = "y", positive_label = "yes", scores = c("s", "s"),
                                                         cutoff = 1)), "invalid_input")
})

test_that("scores must be finite numbers; probability scores must lie in [0, 1]", {
  d <- sp_data()
  bad <- d
  bad$s[2] <- Inf
  expect_identical(th_condition_code(classification_spec(bad, truth = "y", positive_label = "yes", scores = "s",
                                                         cutoff = 0.5)), "invalid_input")
  chr <- d
  chr$s <- as.character(chr$s)
  expect_identical(th_condition_code(classification_spec(chr, truth = "y", positive_label = "yes", scores = "s",
                                                         cutoff = 0.5)), "invalid_input")
  expect_identical(th_condition_code(classification_spec(d, truth = "y", positive_label = "yes", scores = "s2",
                                                         cutoff = 1, score_scale = "probability")), "invalid_input")
  expect_identical(th_condition_code(classification_spec(d, truth = "y", positive_label = "yes", scores = "nope",
                                                         cutoff = 1)), "invalid_input")
})

test_that("weights are frequency weights only; other weights are refused as unsupported", {
  d <- sp_data()
  d$w <- c(1, 2, 0, 1, 3, 1, 1, 2, 1, 1)
  s <- classification_spec(d, truth = "y", positive_label = "yes", scores = "s", cutoff = 0.5, weights = "w")
  expect_identical(s$models[[1]]$weight_semantics, "frequency")
  d$w2 <- d$w + 0.5
  expect_identical(th_condition_code(classification_spec(d, truth = "y", positive_label = "yes", scores = "s",
                                                         cutoff = 0.5, weights = "w2")), "unsupported")
  expect_identical(th_condition_code(classification_spec(d, truth = "y", positive_label = "yes", scores = "s",
                                                         cutoff = 0.5, weights = "w", weight_type = "survey")),
                   "unsupported")
  d$w3 <- -d$w
  expect_identical(th_condition_code(classification_spec(d, truth = "y", positive_label = "yes", scores = "s",
                                                         cutoff = 0.5, weights = "w3")), "invalid_input")
  d$w4 <- 0
  expect_identical(th_condition_code(classification_spec(d, truth = "y", positive_label = "yes", scores = "s",
                                                         cutoff = 0.5, weights = "w4")), "invalid_input")
})

test_that("ids must be unique; repeated units are declared as clusters", {
  d <- sp_data()
  d$id[2] <- d$id[1]
  expect_identical(th_condition_code(classification_spec(d, truth = "y", positive_label = "yes", scores = "s",
                                                         cutoff = 0.5, id = "id")), "invalid_input")
  d <- sp_data()
  d$patient <- rep(1:5, each = 2)
  s <- classification_spec(d, truth = "y", positive_label = "yes", scores = "s", cutoff = 0.5, id = "id",
                           cluster = "patient")
  expect_identical(s$models[[1]]$cluster_column, "patient")
})

test_that("missing values are an error by default or excluded for all models alike", {
  d <- sp_data()
  d$s[3] <- NA
  d$s2[5] <- NA
  expect_identical(th_condition_code(classification_spec(d, truth = "y", positive_label = "yes", scores = c("s", "s2"),
                                                         cutoff = 1)), "invalid_input")
  s <- classification_spec(d, truth = "y", positive_label = "yes", scores = c("s", "s2"), cutoff = 1,
                           missing = "exclude")
  r <- evaluate_cutoff(s)
  expect_equal(r$diagnostics$n_excluded_rows, 2)
  expect_equal(r$diagnostics$n_input_rows, 10)
  # both models are evaluated on the same 8 units
  expect_equal(r$models[[1]]$operating_point$n_total, 8)
  expect_equal(r$models[[2]]$operating_point$n_total, 8)
  expect_true("missing_excluded" %in% r$limitations$code)
})

test_that("aggregate counts and theory refuse data-only settings", {
  expect_identical(th_condition_code(classification_spec(counts = c(tp = 1, fn = 1, fp = 1, tn = 1),
                                                         positive_label = 1)), "invalid_input")
  expect_identical(th_condition_code(classification_spec(counts = c(tp = 1.5, fn = 1, fp = 1, tn = 1))), "invalid_input")
  expect_identical(th_condition_code(classification_spec(counts = c(tp = 1, fn = 1, fp = 1, tn = 1),
                                                         rule = cutoff_rule(selection = "max_youden"))),
                   "invalid_input")
  # a table without a stated cutoff is a binary test result (convention recorded)
  s <- classification_spec(counts = c(tp = 4, fn = 1, fp = 2, tn = 3))
  expect_identical(s$operating_rule$cutoff_source, "binary_test_convention")
  expect_match(s$models[[1]]$score_meaning, "binary test result")
  expect_identical(th_condition_code(classification_spec(theory = list(auc = 0.8, sigma = 1, extra = 2),
                                                         cutoff = 0, prevalence = 0.1)), "invalid_input")
})

test_that("operating rules validate their settings and record defaults", {
  expect_identical(th_condition_code(cutoff_rule()), "invalid_input")
  expect_identical(th_condition_code(cutoff_rule(cutoff = 1, selection = "max_youden")), "invalid_input")
  expect_identical(th_condition_code(cutoff_rule(selection = "constraint")), "invalid_input")
  expect_identical(th_condition_code(cutoff_rule(selection = "constraint", constraints = list(ppv_min = 1.2))),
                   "invalid_input")
  expect_identical(th_condition_code(cutoff_rule(selection = "constraint", constraints = list(specificity = 0.9))),
                   "invalid_input")
  expect_identical(th_condition_code(cutoff_rule(selection = "constraint", constraints = list(tpr_min = 0.9),
                                                 optimise = "roc_auc")), "invalid_input")
  expect_identical(th_condition_code(cutoff_rule(selection = "capacity", capacity = list(max_alert_rate = 2))),
                   "invalid_input")
  expect_identical(th_condition_code(cutoff_rule(selection = "cost_utility", costs = c(fn = 1))), "invalid_input")
  expect_identical(th_condition_code(cutoff_rule(selection = "max_f1", selection_data = "held_out_tuning")),
                   "invalid_input")
  r <- cutoff_rule(selection = "constraint", constraints = list(fpr_max = 0.1))
  expect_identical(r$optimise, "tpr")
  expect_true(all(c("optimise", "selection_data", "tie_break", "prevalence_basis") %in% r$defaults_applied))
  # no policy threshold is built in
  expect_null(cutoff_rule(selection = "max_youden")$constraints)
  # a data-driven selection cannot claim to be predeclared
  d <- sp_data()
  expect_identical(th_condition_code(classification_spec(d, truth = "y", positive_label = "yes", scores = "s",
                                                         rule = cutoff_rule(selection = "max_youden",
                                                                            selection_data = "predeclared_rule"))),
                   "invalid_input")
  # a scenario basis needs a scenario prevalence
  expect_identical(th_condition_code(classification_spec(d, truth = "y", positive_label = "yes", scores = "s",
                                                         rule = cutoff_rule(selection = "max_f1",
                                                                            prevalence_basis = "scenario"))),
                   "invalid_input")
})

test_that("prevalence provenance and sampling design are recorded", {
  d <- sp_data()
  s <- classification_spec(d, truth = "y", positive_label = "yes", scores = "s", cutoff = 0.5)
  expect_identical(s$prevalence$kind, "observed")
  expect_identical(s$prevalence$source, "sample")
  expect_identical(s$prevalence$sampling, "unknown")
  s2 <- classification_spec(d, truth = "y", positive_label = "yes", scores = "s", cutoff = 0.5,
                            prevalence = 0.02, prevalence_source = "external", sampling = "case_control")
  expect_identical(s2$prevalence$kind, "scenario")
  expect_identical(s2$prevalence$source, "external")
  expect_identical(s2$prevalence$sampling, "case_control")
  s3 <- classification_spec(d, truth = "y", positive_label = "yes", scores = "s", cutoff = 0.5, prevalence = 0.02)
  expect_identical(s3$prevalence$source, "assumed")
  expect_true("source" %in% s3$prevalence$defaults_applied)
  expect_identical(th_condition_code(classification_spec(d, truth = "y", positive_label = "yes", scores = "s",
                                                         cutoff = 0.5, prevalence_source = "external")), "invalid_input")
  expect_identical(th_condition_code(classification_spec(d, truth = "y", positive_label = "yes", scores = "s",
                                                         cutoff = 0.5, prevalence = 1.1)), "invalid_input")
  expect_identical(th_condition_code(classification_spec(d, truth = "y", positive_label = "yes", scores = "s",
                                                         cutoff = 0.5, sampling = "convenience")), "invalid_input")
})

test_that("specifications are immutable, hashed and carry only a data fingerprint", {
  d <- sp_data()
  s <- classification_spec(d, truth = "y", positive_label = "yes", scores = "s", cutoff = 0.5)
  expect_match(s$spec_hash, "^[0-9a-f]{64}$")
  expect_match(s$data_fingerprint, "^[0-9a-f]{64}$")
  expect_error(s$input_kind <- "counts_2x2")
  expect_error(s[["models"]] <- list())
  # no raw rows in the serialised specification
  js <- ig_to_json(analysis_envelope(s))
  expect_false(grepl("0.9,0.2,0.7", js, fixed = TRUE))
  expect_false(grepl("\"yes\",\"no\"", js, fixed = TRUE))
  # the same data give the same fingerprint; a changed value does not
  s_same <- classification_spec(d, truth = "y", positive_label = "yes", scores = "s", cutoff = 0.5)
  expect_identical(s_same$spec_hash, s$spec_hash)
  d2 <- d
  d2$s[1] <- 0.91
  expect_false(identical(classification_spec(d2, truth = "y", positive_label = "yes", scores = "s",
                                             cutoff = 0.5)$data_fingerprint, s$data_fingerprint))
  # tampering is detected
  tampered <- unclass(s)
  tampered$operating_rule$score_cutoff <- 0.6
  tampered <- structure(tampered, class = class(s))
  expect_identical(th_condition_code(evaluate_cutoff(tampered)), "invalid_input")
})

test_that("update() builds a new validated specification and keeps the data", {
  d <- sp_data()
  s <- classification_spec(d, truth = "y", positive_label = "yes", scores = "s", cutoff = 0.5)
  s2 <- update(s, cutoff = 0.7)
  expect_identical(s2$operating_rule$score_cutoff, 0.7)
  expect_false(identical(s2$spec_hash, s$spec_hash))
  expect_identical(s2$data_fingerprint, s$data_fingerprint)
  expect_equal(evaluate_cutoff(s2)$models[[1]]$operating_point$n_alerts, 3)
  s3 <- update(s, prevalence = 0.05)
  expect_identical(s3$prevalence$kind, "scenario")
  s4 <- update(s3, prevalence = NA)
  expect_identical(s4$prevalence$kind, "observed")
  s5 <- update(s, evaluation_split = "held_out_test")
  expect_identical(s5$models[[1]]$evaluation_split, "held_out_test")
  expect_false("evaluation_split" %in% s5$models[[1]]$defaults_applied)
  expect_identical(th_condition_code(update(s, theory = list(auc = 0.9))), "invalid_input")
  expect_identical(th_condition_code(update(s, unknown = 1)), "invalid_input")
})
