# Scalars apply to every model; heterogeneous scales require explicit model keys.
.ig_cls_score_scales <- function(scale, ids) {
  if(!is.character(scale) || !length(scale) || anyNA(scale) ||
     !all(scale %in% c("raw_score","probability")))
    .ig_abort("invalid_input","score_scale must contain raw_score or probability")
  if(length(scale)==1L && is.null(names(scale))) return(stats::setNames(rep(scale,length(ids)),ids))
  if(is.null(names(scale)) || anyDuplicated(names(scale)) || !setequal(names(scale),ids))
    .ig_abort("invalid_input","per-model score_scale must have one value named for every model id")
  scale[ids]
}

# Classification specifications ----------------------------------------------------
#
# classification_spec() builds the immutable, hashed specification of a
# binary-classification analysis on top of the core shape
# (.ig_new_classification_spec in core-envelope.R). Three input kinds:
#
#   scored_data          truth labels plus one or more numeric score columns
#                        (one model per column, matched on the same units);
#   counts_2x2           aggregate TP/FN/FP/TN at one fixed operating point;
#   theoretical_binormal S0 ~ N(0, 1), S1 ~ N(d, sigma^2),
#                        d = qnorm(AUC) * sqrt(1 + sigma^2).
#
# The specification never contains raw rows. For scored data it holds a
# SHA-256 fingerprint of the analysed data; the prepared data (0/1 truth,
# scores, frequency weights, clusters) are attached to the R object as the
# attribute "ig_cls_data", which is not serialised. A specification read back
# from JSON needs attach_classification_data() before it can be evaluated.
#
# Every choice the user did not make explicitly (direction, comparator,
# evaluation split, selection data, tie break, theory N) is written out and
# listed in `defaults_applied` so that it is visible and confirmable.

.ig_cls_constraint_names <- c("tpr_min", "tnr_min", "ppv_min", "npv_min", "fpr_max", "fnr_max",
                              "fdr_max", "alert_rate_max")
.ig_cls_tie_breaks <- c("fewest_alerts", "most_alerts")
.ig_cls_prevalence_sources <- c("sample", "assumed", "external")
.ig_cls_sampling <- c("cross_sectional", "case_control", "unknown", "not_applicable")
.ig_cls_data_attr <- "ig_cls_data"

#' Model metadata for a classification specification
#'
#' Describes one already-scored model (or test). inferguidR does not train
#' models and never deserialises model objects; it evaluates the scores you
#' supply.
#'
#' @param model_id Stable identifier (letters, digits, `_`, `.`, `-`).
#' @param name Display name.
#' @param version Optional model version string.
#' @param score_meaning What the score is (for example "predicted 30-day
#'   risk"). Defaults depend on the input kind.
#' @param population Optional description of the population the scores come
#'   from.
#' @param unit Optional unit noun of an evaluation unit (for example
#'   "patient" or "order").
#' @param source_hash Optional SHA-256 (64 lower-case hex characters) of the
#'   source data set or model artefact.
#' @return A list of class `ig_classifier_model`.
#' @examples
#' classifier_model("logit_v2", name = "Logistic model", version = "2.0",
#'                  score_meaning = "predicted probability of readmission")
#' @concept classification
#' @seealso
#' * [classification_spec()]
#' * [evaluate_cutoff()]
#' * [plot_classification()]
#' @export
classifier_model <- function(model_id, name = model_id, version = NULL, score_meaning = NULL,
                             population = NULL, unit = NULL, source_hash = NULL) {
  p <- character()
  if (!.ig_is_string(model_id) || !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", model_id)) {
    p <- c(p, "model_id must be a non-empty identifier (letters, digits, '_', '.', '-')")
  }
  for (nm in c("name")) if (!.ig_is_string(get(nm))) p <- c(p, sprintf("%s must be a string", nm))
  for (nm in c("version", "score_meaning", "population", "unit")) {
    v <- get(nm)
    if (!is.null(v) && !.ig_is_string(v)) p <- c(p, sprintf("%s must be a string or NULL", nm))
  }
  if (!is.null(source_hash) && !(.ig_is_string(source_hash) && grepl("^[0-9a-f]{64}$", source_hash))) {
    p <- c(p, "source_hash must be a SHA-256 hex string (64 lower-case characters) or NULL")
  }
  if (length(p)) .ig_abort_invalid(p, "invalid model metadata")
  structure(list(model_id = model_id, name = name, version = version, score_meaning = score_meaning,
                 population = population, unit = unit, source_hash = source_hash),
            class = "ig_classifier_model")
}

#' Operating rule of a classification analysis
#'
#' Declares how the cutoff (operating point) is obtained: a manual cutoff, or
#' a selection objective with optional constraints. Policy thresholds are never
#' built in; every constraint is a user decision.
#'
#' @param cutoff Numeric score cutoff for a manual rule. It is a classifier
#'   score threshold, never a significance critical value.
#' @param selection `"manual"` (default when `cutoff` is given), `"max_youden"`,
#'   `"max_f1"`, `"max_mcc"`, `"constraint"` (best `optimise` value subject to
#'   `constraints`), `"capacity"` (most alerts within `capacity`) or
#'   `"cost_utility"` (smallest declared expected cost).
#' @param selection_data Where the cutoff comes from relative to the evaluated
#'   data: `"same_data"` (chosen on these data; the default because it claims
#'   least), `"held_out_tuning"` (chosen on separate tuning data) or
#'   `"predeclared_rule"` (fixed before any data were seen; manual rules only).
#' @param constraints Optional named list of bounds: `tpr_min`, `tnr_min`,
#'   `ppv_min`, `npv_min`, `fpr_max`, `fnr_max`, `fdr_max`, `alert_rate_max`
#'   (values in \[0, 1\]). An undefined metric never satisfies a bound. With a
#'   selection objective the bounds restrict the admissible operating points;
#'   with a manual cutoff they are only checked and reported.
#' @param optimise For `selection = "constraint"`: the metric id to optimise
#'   among feasible operating points (default `"tpr"`), in its preferred
#'   direction from [metric_dictionary()].
#' @param capacity For `selection = "capacity"`: `list(max_alerts = k)` or
#'   `list(max_alert_rate = r)`.
#' @param costs For `selection = "cost_utility"`: named numeric costs per unit
#'   for `fn` and `fp` (required) and optionally `tp` and `tn` (negative values
#'   are benefits).
#' @param tie_break Rule that chooses among operating points with the same
#'   objective value: `"fewest_alerts"` (default) or `"most_alerts"`. All tied
#'   operating points are reported.
#' @param prevalence_basis Counts on which a selection objective and the
#'   constraints are evaluated: `"observed"` (the evaluated sample; default for
#'   scored data) or `"scenario"` (expected counts at the specification's
#'   scenario prevalence under fixed class-conditional score distributions;
#'   the only basis of a theoretical model). It matters for
#'   prevalence-dependent quantities such as PPV, F1, MCC or costs.
#' @return A list of class `ig_cutoff_rule`.
#' @examples
#' cutoff_rule(cutoff = 0.3, selection_data = "predeclared_rule")
#' cutoff_rule(selection = "constraint", constraints = list(tpr_min = 0.9))
#' cutoff_rule(selection = "max_youden", constraints = list(fpr_max = 0.2))
#' cutoff_rule(selection = "cost_utility", costs = c(fn = 10, fp = 1))
#' @concept classification
#' @seealso
#' * [classification_spec()]
#' * [evaluate_cutoff()]
#' * [plot_classification()]
#' @export
cutoff_rule <- function(cutoff = NULL, selection = NULL, selection_data = NULL,
                        constraints = NULL, optimise = NULL, capacity = NULL,
                        costs = NULL, tie_break = NULL, prevalence_basis = NULL) {
  p <- character()
  defaults <- character()
  vocab <- .ig_classification_vocab
  if (!is.null(cutoff) && !(is.numeric(cutoff) && length(cutoff) == 1L && is.finite(cutoff))) {
    p <- c(p, "cutoff must be a single finite number")
  }
  if (is.null(selection)) {
    if (is.null(cutoff)) {
      .ig_abort("invalid_input",
                "no operating rule: give `cutoff` (manual rule) or a `selection` objective")
    }
    selection <- "manual"
  }
  if (!.ig_is_string(selection) || !selection %in% vocab$selection) {
    p <- c(p, sprintf("selection must be one of: %s", paste(vocab$selection, collapse = ", ")))
  }
  if (is.null(selection_data)) {
    selection_data <- "same_data"
    defaults <- c(defaults, "selection_data")
  }
  if (!.ig_is_string(selection_data) || !selection_data %in% vocab$selection_data) {
    p <- c(p, sprintf("selection_data must be one of: %s", paste(vocab$selection_data, collapse = ", ")))
  }
  if (.ig_is_string(selection)) {
    if (selection == "manual" && is.null(cutoff)) p <- c(p, "a manual rule needs `cutoff`")
    if (selection != "manual" && !is.null(cutoff)) {
      p <- c(p, sprintf("selection '%s' chooses the cutoff from the data; do not also give `cutoff`", selection))
    }
    if (selection != "manual" && identical(selection_data, "held_out_tuning")) {
      p <- c(p, paste(
        "a selection objective chooses the cutoff on the data it is applied to; to choose it on",
        "tuning data, run select_cutoff() on the tuning specification and pass that result as",
        "`cutoff` to evaluate_cutoff() on the evaluation specification"))
    }
    if (selection == "constraint" && !length(constraints)) p <- c(p, "selection 'constraint' needs `constraints`")
    if (selection == "capacity" && is.null(capacity)) p <- c(p, "selection 'capacity' needs `capacity`")
    if (selection == "cost_utility" && is.null(costs)) p <- c(p, "selection 'cost_utility' needs `costs`")
    if (selection != "capacity" && !is.null(capacity)) p <- c(p, "`capacity` is only used with selection 'capacity'")
    if (selection != "cost_utility" && !is.null(costs)) p <- c(p, "`costs` are only used with selection 'cost_utility'")
    if (selection != "constraint" && !is.null(optimise)) p <- c(p, "`optimise` is only used with selection 'constraint'")
  }
  if (length(constraints)) {
    cn <- names(constraints)
    if (!is.list(constraints) && !is.numeric(constraints)) p <- c(p, "constraints must be a named list")
    if (is.null(cn) || any(!nzchar(cn)) || anyDuplicated(cn)) {
      p <- c(p, "constraints must be a named list without duplicate names")
    } else {
      bad <- setdiff(cn, .ig_cls_constraint_names)
      if (length(bad)) {
        p <- c(p, sprintf("unknown constraint '%s'; allowed: %s", bad,
                          paste(.ig_cls_constraint_names, collapse = ", ")))
      }
      vals <- unlist(constraints)
      if (!is.numeric(vals) || length(vals) != length(cn) || any(!is.finite(vals)) ||
          any(vals < 0 | vals > 1)) {
        p <- c(p, "constraint bounds must be single numbers in [0, 1]")
      }
    }
    constraints <- as.list(constraints)[sort(names(constraints))]
    constraints <- lapply(constraints, as.double)
  } else {
    constraints <- NULL
  }
  if (identical(selection, "constraint")) {
    if (is.null(optimise)) {
      optimise <- "tpr"
      defaults <- c(defaults, "optimise")
    }
    entry <- .ig_cls_dictionary()$metrics[[optimise]]
    if (!.ig_is_string(optimise) || is.null(entry) || !identical(entry$scope, "operating_point") ||
        !entry$preferred_direction %in% c("higher", "lower")) {
      p <- c(p, "`optimise` must be an operating-point metric id with a preferred direction (see metric_dictionary())")
    }
  }
  if (!is.null(capacity)) {
    ok <- is.list(capacity) && length(capacity) == 1L &&
      names(capacity) %in% c("max_alerts", "max_alert_rate") &&
      is.numeric(capacity[[1]]) && length(capacity[[1]]) == 1L && is.finite(capacity[[1]]) &&
      capacity[[1]] >= 0 &&
      (names(capacity) != "max_alert_rate" || capacity[[1]] <= 1)
    if (!isTRUE(ok)) p <- c(p, "capacity must be list(max_alerts = k) with k >= 0 or list(max_alert_rate = r) with r in [0, 1]")
    else capacity <- lapply(capacity, as.double)
  }
  if (!is.null(costs)) {
    cn <- names(costs)
    vals <- suppressWarnings(as.double(unlist(costs)))
    ok <- !is.null(cn) && all(cn %in% c("tp", "fn", "fp", "tn")) && !anyDuplicated(cn) &&
      all(c("fn", "fp") %in% cn) && length(vals) == length(cn) && all(is.finite(vals))
    if (!ok) p <- c(p, "costs must be finite named numbers with fn and fp (and optionally tp, tn)")
    else costs <- as.list(stats::setNames(vals, cn))[sort(cn)]
  }
  if (is.null(tie_break)) {
    tie_break <- "fewest_alerts"
    if (!identical(selection, "manual")) defaults <- c(defaults, "tie_break")
  }
  if (!.ig_is_string(tie_break) || !tie_break %in% .ig_cls_tie_breaks) {
    p <- c(p, sprintf("tie_break must be one of: %s", paste(.ig_cls_tie_breaks, collapse = ", ")))
  }
  if (is.null(prevalence_basis)) {
    prevalence_basis <- "observed"
    if (!identical(selection, "manual") || length(constraints)) {
      defaults <- c(defaults, "prevalence_basis")
    }
  }
  if (!.ig_is_string(prevalence_basis) || !prevalence_basis %in% c("observed", "scenario")) {
    p <- c(p, "prevalence_basis must be 'observed' or 'scenario'")
  }
  if (length(p)) .ig_abort_invalid(p, "invalid operating rule")
  structure(list(score_cutoff = if (is.null(cutoff)) NULL else as.double(cutoff),
                 selection = selection, selection_data = selection_data,
                 constraints = constraints, optimise = optimise, capacity = capacity,
                 costs = costs, tie_break = tie_break, prevalence_basis = prevalence_basis,
                 defaults_applied = if (length(defaults)) sort(defaults) else NULL),
            class = "ig_cutoff_rule")
}

# Checks and adjustments of an operating rule that depend on the input kind
# and the prevalence settings (shared by classification_spec() and update()).
.ig_cls_finalise_rule <- function(input_kind, op, prevalence) {
  p <- character()
  defaults <- op$defaults_applied
  data_driven <- !identical(op$selection, "manual")
  if (input_kind == "counts_2x2" && data_driven) {
    p <- c(p, "aggregate 2x2 counts describe one fixed operating point; no cutoff can be selected from them")
  }
  if (input_kind == "theoretical_binormal") {
    # A cutoff set or optimised on an assumed model uses no outcome data, so it
    # is a rule declared before any data are seen.
    if ("selection_data" %in% defaults) {
      op$selection_data <- "predeclared_rule"
      defaults <- c(setdiff(defaults, "selection_data"), "selection_data_theoretical")
    } else if (!identical(op$selection_data, "predeclared_rule")) {
      p <- c(p, "a theoretical model involves no outcome data: selection_data must be 'predeclared_rule'")
    }
    if (identical(op$prevalence_basis, "observed") && !"prevalence_basis" %in% defaults &&
        (data_driven || length(op$constraints))) {
      p <- c(p, "a theoretical model has no observed sample: prevalence_basis must be 'scenario'")
    }
    # At prevalence 0 or 1 the expected counts lack a class: most objectives are
    # undefined or constant along a plateau that reaches the limiting rules, so
    # no optimum can be reported honestly.
    if (data_driven && is.numeric(prevalence$value) && (prevalence$value <= 0 || prevalence$value >= 1)) {
      p <- c(p, "selecting a cutoff on a theoretical model needs a prevalence strictly between 0 and 1")
    }
    op$prevalence_basis <- "scenario"
    defaults <- setdiff(defaults, "prevalence_basis")
  } else if (data_driven && identical(op$selection_data, "predeclared_rule")) {
    p <- c(p, "a cutoff selected from the evaluated data is not a predeclared rule; use 'same_data'")
  }
  if (input_kind != "theoretical_binormal" && identical(op$prevalence_basis, "scenario") &&
      !identical(prevalence$kind, "scenario")) {
    p <- c(p, "prevalence_basis 'scenario' needs a scenario `prevalence` in the specification")
  }
  if (length(p)) .ig_abort_invalid(p, "invalid operating rule")
  op["defaults_applied"] <- list(if (length(defaults)) sort(unique(defaults)) else NULL)
  op
}

#' Uncertainty settings of a classification analysis
#'
#' Intervals are conditional on the scores and, for operating-point metrics,
#' on the fixed cutoff. They do not include the uncertainty of a data-driven
#' cutoff selection.
#'
#' Bootstrap limits use type-7 empirical quantiles. Undefined replicates are
#' excluded only when at most 10 percent are undefined; otherwise the interval
#' is explicitly unavailable. Class stratification conditions on the observed
#' class totals, so it does not estimate population-prevalence uncertainty.
#' Frequency weights represent replicated units; declared clusters are sampled
#' as whole clusters. These are percentile intervals, without bias correction
#' or studentisation. Their nominal level is not a coverage guarantee: small
#' samples and boundary estimates can have severe undercoverage. Independent
#' arithmetic checks and a binomial coverage diagnostic are provided in
#' `tools/product/check-bootstrap.R`; bootstrap validation remains development.
#'
#' @param level Confidence level in (0, 1).
#' @param proportion_methods Binomial interval methods for proportions:
#'   any of `"wilson"` and `"clopper_pearson"`, or `character(0)` for none.
#' @param auc_method `"delong"` (asymptotic interval, own implementation) or
#'   `"none"`.
#' @param bootstrap_replicates Number of bootstrap replicates (0 = no
#'   bootstrap). The bootstrap needs a `seed`.
#' @param seed Integer seed for the bootstrap. The global random-number state
#'   is restored afterwards.
#' @param bootstrap_strata `"class"` (resample events and non-events
#'   separately, the default) or `"none"`. With declared clusters, whole
#'   clusters are resampled instead.
#' @param bootstrap_metrics Metric ids that receive percentile bootstrap
#'   intervals: operating-point metrics (evaluated at the fixed cutoff in every
#'   replicate) and the curve summaries `roc_auc`, `partial_auc`,
#'   `average_precision`, `pr_area_trapezoid` and `brier_score`.
#' @return A list of class `ig_cls_uncertainty`.
#' @examples
#' classification_uncertainty(level = 0.9, bootstrap_replicates = 500, seed = 1)
#' @concept classification
#' @seealso
#' * [classification_spec()]
#' * [evaluate_cutoff()]
#' * [plot_classification()]
#' @export
classification_uncertainty <- function(level = 0.95, proportion_methods = c("wilson", "clopper_pearson"),
                                       auc_method = c("delong", "none"), bootstrap_replicates = 0,
                                       seed = NULL, bootstrap_strata = c("class", "none"),
                                       bootstrap_metrics = c("tpr", "fpr", "ppv", "npv", "youden_j",
                                                             "f1", "mcc", "roc_auc",
                                                             "average_precision")) {
  p <- character()
  auc_method <- match.arg(auc_method)
  bootstrap_strata <- match.arg(bootstrap_strata)
  allowed_boot <- c(.ig_cls_operating_ids(), .ig_cls_bootstrap_curve_ids)
  if (!is.character(bootstrap_metrics) || anyNA(bootstrap_metrics) || anyDuplicated(bootstrap_metrics) ||
      !all(bootstrap_metrics %in% allowed_boot)) {
    p <- c(p, "bootstrap_metrics must be distinct operating-point metric ids or curve summaries (roc_auc, partial_auc, average_precision, pr_area_trapezoid, brier_score)")
  }
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    p <- c(p, "level must be a single number in (0, 1)")
  }
  if (!is.character(proportion_methods) || anyNA(proportion_methods) ||
      !all(proportion_methods %in% c("wilson", "clopper_pearson")) || anyDuplicated(proportion_methods)) {
    p <- c(p, "proportion_methods must be a subset of c('wilson', 'clopper_pearson')")
  }
  if (!is.numeric(bootstrap_replicates) || length(bootstrap_replicates) != 1L ||
      !is.finite(bootstrap_replicates) || bootstrap_replicates < 0 ||
      bootstrap_replicates != round(bootstrap_replicates)) {
    p <- c(p, "bootstrap_replicates must be a whole number >= 0")
  } else if (bootstrap_replicates > 0 && bootstrap_replicates < 100) {
    p <- c(p, "use at least 100 bootstrap replicates (or 0 for none)")
  }
  if (!is.null(seed) && !(is.numeric(seed) && length(seed) == 1L && is.finite(seed) && seed == round(seed))) {
    p <- c(p, "seed must be a single whole number or NULL")
  }
  if (is.numeric(bootstrap_replicates) && length(bootstrap_replicates) == 1L &&
      isTRUE(bootstrap_replicates > 0) && is.null(seed)) {
    p <- c(p, "a bootstrap needs a `seed` so that it can be reproduced")
  }
  if (length(p)) .ig_abort_invalid(p, "invalid uncertainty settings")
  boot_on <- bootstrap_replicates > 0
  structure(list(level = as.double(level),
                 proportion_methods = if (length(proportion_methods)) proportion_methods else NULL,
                 auc_method = auc_method,
                 bootstrap_replicates = as.double(bootstrap_replicates),
                 seed = if (is.null(seed)) NULL else as.double(seed),
                 bootstrap_strata = bootstrap_strata,
                 bootstrap_interval = "percentile",
                 bootstrap_metrics = if (boot_on && length(bootstrap_metrics)) bootstrap_metrics else NULL,
                 bootstrap_rng = if (boot_on) "Mersenne-Twister/Inversion/Rejection" else NULL),
            class = "ig_cls_uncertainty")
}

.ig_cls_bootstrap_curve_ids <- c("roc_auc", "partial_auc", "average_precision", "pr_area_trapezoid",
                                 "brier_score")

#' Specify a binary classification or cutoff analysis
#'
#' Builds an immutable, hashed `classification_spec`. Give exactly one kind of
#' input:
#'
#' * **Scored data:** `data` with a `truth` column, an explicit
#'   `positive_label` and one or more numeric `scores` columns (one model per
#'   column; all models are evaluated on the same units).
#' * **Aggregate counts:** `counts = c(tp =, fn =, fp =, tn =)` for one fixed
#'   operating point. No ROC curve or cutoff sweep can be computed from these.
#' * **Theory:** `theory = list(auc =, sigma =, n_total =)` for the binormal
#'   model \eqn{S_0 \sim N(0,1)}, \eqn{S_1 \sim N(d, \sigma^2)} with
#'   \eqn{d = \Phi^{-1}(AUC)\sqrt{1+\sigma^2}}; the cutoff is on the score
#'   scale standardised to the non-event distribution and needs a prevalence.
#'
#' The default direction is "higher score means positive" and is recorded as a
#' default in the model record (`defaults_applied`) so that it can be
#' confirmed; the direction is never chosen from the data. The comparator
#' (`>=`, `>`, `<=`, `<`) fixes which side of the cutoff alerts and what
#' happens to scores equal to the cutoff: under `>=` a cutoff of 1 alerts
#' units with score 1.
#'
#' Weights are frequency weights only (whole numbers meaning replicated
#' units); survey or importance weights are refused as `unsupported`. With a
#' `cluster` column, binomial and DeLong intervals are not reported and the
#' bootstrap resamples whole clusters.
#'
#' @param data A data frame for scored data.
#' @param truth Name of the truth column in `data`.
#' @param positive_label The value of `truth` that marks an event. Required:
#'   the positive class is never inferred.
#' @param scores Character vector of score column names; names (if given) are
#'   used as model ids.
#' @param counts Named counts `tp`, `fn`, `fp`, `tn` (whole numbers).
#' @param theory List with `auc` in (0, 1), `sigma` > 0 and optionally
#'   `n_total` (expected-count scale, default 1000, at most 1e100).
#' @param rule An operating rule from [cutoff_rule()].
#' @param cutoff Shorthand for `rule = cutoff_rule(cutoff = cutoff)`.
#' @param models Optional model metadata: a [classifier_model()] or a list of
#'   them, matched to the score columns by `model_id` or position.
#' @param direction `"higher_is_positive"` (default) or
#'   `"lower_is_positive"`; a single value or one per model (named by model
#'   id).
#' @param comparator `">="`/`">"` (higher is positive) or `"<="`/`"<"`
#'   (lower is positive); default `">="` or `"<="`. An explicit comparator
#'   implies the direction.
#' @param score_scale `"raw_score"` (default) or `"probability"`, shared by all
#'   models, or a character vector with one scale per model named by model id.
#'   Probability scores are
#'   predicted probabilities of the positive class; values must lie in
#'   \[0, 1\]. Brier score, calibration and net benefit need probabilities
#'   and the direction `"higher_is_positive"`; otherwise they are reported as
#'   not applicable.
#' @param evaluation_split How the evaluated data relate to model development:
#'   `"same_data"` (default, claims least), `"held_out_tuning"`,
#'   `"held_out_test"`, `"external"` or `"resampling"`.
#' @param negative_label Optional value of `truth` that marks a non-event;
#'   when given, every other value is an error.
#' @param prevalence Optional prevalence scenario in \[0, 1\] (required for
#'   theory; strictly between 0 and 1 when a cutoff is selected on a
#'   theoretical model). Observed results are kept separately.
#' @param prevalence_source Provenance of a scenario prevalence: `"assumed"`
#'   (default) or `"external"`.
#' @param sampling How the evaluation sample was drawn: `"cross_sectional"`,
#'   `"case_control"` (outcome-selected; the sample prevalence is not a
#'   population prevalence) or `"unknown"` (default).
#' @param weights Optional name of a frequency-weight column.
#' @param weight_type Only `"frequency"` is supported.
#' @param cluster Optional name of a cluster column (for example patient id
#'   when patients contribute several units).
#' @param id Optional name of a unit-id column; duplicates are an error.
#' @param missing `"error"` (default) or `"exclude"`: rows with a missing
#'   truth, score, weight or cluster value are dropped for all models alike
#'   and counted.
#' @param uncertainty Settings from [classification_uncertainty()].
#' @return An immutable object of class `classification_spec`.
#' @examples
#' # Aggregate counts at one fixed operating point
#' classification_spec(counts = c(tp = 40, fn = 10, fp = 20, tn = 130), cutoff = 0.5)
#'
#' # Scored data (synthetic)
#' d <- data.frame(y = c(1, 1, 1, 0, 0, 0, 0, 0), risk = c(.9, .8, .4, .7, .3, .2, .2, .1))
#' spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "risk",
#'                             cutoff = 0.5, score_scale = "probability")
#' spec
#'
#' # Binormal theory: AUC 0.8, unequal variance, 10 % prevalence
#' classification_spec(theory = list(auc = 0.8, sigma = 1.5, n_total = 1000),
#'                     cutoff = 1, prevalence = 0.1)
#' @seealso [evaluate_cutoff()], [select_cutoff()], [cutoff_sweep()]
#' @concept classification
#' @seealso
#' * [evaluate_cutoff()]
#' * [plot_classification()]
#' @export
classification_spec <- function(data = NULL, truth = NULL, positive_label = NULL, scores = NULL,
                                counts = NULL, theory = NULL, rule = NULL, cutoff = NULL,
                                models = NULL, direction = NULL, comparator = NULL,
                                score_scale = c("raw_score", "probability"),
                                evaluation_split = NULL, negative_label = NULL,
                                prevalence = NULL, prevalence_source = NULL, sampling = NULL,
                                weights = NULL, weight_type = "frequency", cluster = NULL, id = NULL,
                                missing = c("error", "exclude"),
                                uncertainty = classification_uncertainty()) {
  if(missing(score_scale)) score_scale <- "raw_score"
  missing <- match.arg(missing)
  given <- c(scored_data = !is.null(data) || !is.null(scores),
             counts_2x2 = !is.null(counts),
             theoretical_binormal = !is.null(theory))
  if (sum(given) != 1L) {
    .ig_abort("invalid_input",
              "give exactly one input: `data` with `scores`, `counts`, or `theory`")
  }
  input_kind <- names(given)[given]
  if (!is.null(rule) && !is.null(cutoff)) .ig_abort("invalid_input", "give either `rule` or `cutoff`, not both")
  binary_convention <- FALSE
  if (is.null(rule)) {
    if (is.null(cutoff) && input_kind == "counts_2x2") {
      # A 2x2 table without a stated cutoff is the result of a binary decision:
      # represented as a binary test result (1 = alert) with cutoff 1 under >=.
      cutoff <- 1
      binary_convention <- TRUE
    }
    rule <- cutoff_rule(cutoff = cutoff)
  }
  if (!inherits(rule, "ig_cutoff_rule")) .ig_abort("invalid_input", "`rule` must come from cutoff_rule()")
  if (!inherits(uncertainty, "ig_cls_uncertainty")) {
    .ig_abort("invalid_input", "`uncertainty` must come from classification_uncertainty()")
  }
  if (!identical(weight_type, "frequency")) {
    .ig_abort("unsupported", sprintf(
      "weight_type '%s' is not supported: only frequency weights (whole numbers meaning replicated units) have a defined interval method here; survey or importance weights are not treated as independent replicated units",
      paste(format(weight_type), collapse = "")))
  }

  prep <- NULL
  model_ids <- NULL
  if (input_kind == "scored_data") {
    if (!is.data.frame(data)) .ig_abort("invalid_input", "`data` must be a data frame")
    if (!is.character(scores) || !length(scores) || anyNA(scores)) {
      .ig_abort("invalid_input", "`scores` must name one or more score columns")
    }
    model_ids <- if (!is.null(names(scores)) && all(nzchar(names(scores)))) names(scores) else unname(scores)
    if (anyDuplicated(model_ids)) .ig_abort("invalid_input", "model ids (score column names) must be unique")
    prep <- .ig_cls_prepare_data(data, truth, positive_label, negative_label,
                                 stats::setNames(unname(scores), model_ids),
                                 weights, cluster, id, missing, score_scale)
  } else if (input_kind == "counts_2x2") {
    model_ids <- "test"
    counts <- .ig_cls_check_counts(counts, "observed")
    if (!is.null(positive_label) || !is.null(truth) || !is.null(weights) || !is.null(cluster) || !is.null(id)) {
      .ig_abort("invalid_input", "truth, positive_label, weights, cluster and id apply to scored data only")
    }
  } else {
    model_ids <- "binormal"
    if (!is.null(positive_label) || !is.null(truth) || !is.null(weights) || !is.null(cluster)) {
      .ig_abort("invalid_input", "truth, positive_label, weights and cluster apply to scored data only")
    }
  }

  score_scale <- .ig_cls_score_scales(score_scale,model_ids)
  model_records <- .ig_cls_model_records(input_kind, model_ids, models, direction, comparator,
                                         score_scale, evaluation_split, positive_label,
                                         prep, missing, truth, scores, weights, cluster, id,
                                         binary_convention)
  op_rule <- unclass(rule)
  if (binary_convention) op_rule$cutoff_source <- "binary_test_convention"
  prev <- .ig_cls_prevalence_record(input_kind, prevalence, prevalence_source, sampling)
  th <- NULL
  if (input_kind == "theoretical_binormal") th <- .ig_cls_theory_record(theory)
  .ig_cls_build_spec(input_kind, model_records, op_rule, prev,
                     counts = if (is.null(counts)) NULL else as.list(counts),
                     theory = th, prep = prep, uncertainty = unclass(uncertainty))
}

# Assemble and validate the core shape, then attach the private data.
.ig_cls_build_spec <- function(input_kind, models, operating_rule, prevalence, counts, theory,
                               prep, uncertainty) {
  operating_rule <- .ig_cls_finalise_rule(input_kind, operating_rule, prevalence)
  if (input_kind == "theoretical_binormal") {
    operating_rule$cutoff_scale <- "z_standardised_to_non_events"
  }
  if (input_kind != "scored_data" && !is.null(uncertainty$bootstrap_metrics)) {
    .ig_abort("invalid_input", "the bootstrap resamples scored data; it is not available for 2x2 counts or a theoretical model")
  }
  fingerprint <- if (is.null(prep)) NULL else .ig_cls_fingerprint(prep)
  spec <- .ig_new_classification_spec(
    input_kind = input_kind, models = unname(models), operating_rule = operating_rule,
    prevalence = prevalence, counts = counts, theory = theory,
    data_fingerprint = fingerprint, uncertainty = uncertainty)
  if (!is.null(prep)) attr(spec, .ig_cls_data_attr) <- prep
  spec
}

.ig_cls_model_records <- function(input_kind, model_ids, models, direction, comparator,
                                  score_scale, evaluation_split, positive_label, prep, missing,
                                  truth, scores, weights, cluster, id, binary_convention) {
  k <- length(model_ids)
  p <- character()
  # Metadata: a single classifier_model, or a list matched by id or position.
  if (inherits(models, "ig_classifier_model")) models <- list(models)
  if (!is.null(models)) {
    if (!is.list(models) || !all(vapply(models, inherits, logical(1), "ig_classifier_model"))) {
      .ig_abort("invalid_input", "`models` must be a classifier_model() or a list of them")
    }
    ids <- vapply(models, `[[`, "", "model_id")
    if (all(ids %in% model_ids) && !anyDuplicated(ids)) {
      meta <- lapply(model_ids, function(m) if (m %in% ids) models[[match(m, ids)]] else NULL)
    } else if (length(models) == k) {
      meta <- models
    } else {
      .ig_abort("invalid_input", "`models` must match the score columns by model_id or by position")
    }
  } else {
    meta <- vector("list", k)
  }
  per_model <- function(x, name, allowed) {
    if (is.null(x)) return(rep(list(NULL), k))
    if (length(x) == 1L && is.null(names(x))) return(rep(list(unname(x)), k))
    if (is.null(names(x)) || !setequal(names(x), model_ids)) {
      p <<- c(p, sprintf("`%s` must be a single value or one value per model, named by model id", name))
      return(rep(list(NULL), k))
    }
    lapply(model_ids, function(m) unname(x[[m]]))
  }
  dirs <- per_model(direction, "direction")
  comps <- per_model(comparator, "comparator")
  scales <- .ig_cls_score_scales(score_scale,model_ids)
  if (length(p)) .ig_abort_invalid(p, "invalid classification specification")
  if (input_kind == "theoretical_binormal" && !is.null(evaluation_split) &&
      !identical(evaluation_split, "not_applicable")) {
    .ig_abort("invalid_input", "a theoretical model has no evaluation split")
  }
  splits <- c("same_data", "held_out_tuning", "held_out_test", "external", "resampling")
  if (!is.null(evaluation_split) && input_kind != "theoretical_binormal" &&
      !(.ig_is_string(evaluation_split) && evaluation_split %in% splits)) {
    .ig_abort("invalid_input", sprintf("evaluation_split must be one of: %s", paste(splits, collapse = ", ")))
  }
  records <- vector("list", k)
  for (i in seq_len(k)) {
    defaults <- character()
    d <- dirs[[i]]
    cmp <- comps[[i]]
    if (!is.null(cmp) && !(.ig_is_string(cmp) && cmp %in% c(">=", ">", "<=", "<"))) {
      p <- c(p, "comparator must be one of: >=, >, <=, <")
      next
    }
    if (!is.null(d) && !(.ig_is_string(d) && d %in% .ig_classification_vocab$direction)) {
      p <- c(p, "direction must be 'higher_is_positive' or 'lower_is_positive'")
      next
    }
    if (is.null(d)) {
      if (is.null(cmp)) {
        d <- "higher_is_positive"
        defaults <- c(defaults, "direction")
      } else {
        d <- if (cmp %in% c(">=", ">")) "higher_is_positive" else "lower_is_positive"
        defaults <- c(defaults, "direction_implied_by_comparator")
      }
    }
    if (is.null(cmp)) {
      cmp <- if (d == "higher_is_positive") ">=" else "<="
      defaults <- c(defaults, "comparator")
    }
    if ((d == "higher_is_positive") != (cmp %in% c(">=", ">"))) {
      p <- c(p, sprintf("model '%s': comparator '%s' contradicts direction '%s'", model_ids[i], cmp, d))
      next
    }
    if (input_kind == "theoretical_binormal" && d != "higher_is_positive") {
      p <- c(p, "the binormal model is defined with higher scores meaning positive")
      next
    }
    split <- evaluation_split
    if (input_kind == "theoretical_binormal") {
      split <- "not_applicable"
    } else if (is.null(split)) {
      split <- "same_data"
      defaults <- c(defaults, "evaluation_split")
    }
    m <- meta[[i]]
    default_meaning <- switch(input_kind,
      scored_data = if (scales[[i]] == "probability") "predicted probability of the positive class" else
        "model score (raw scale)",
      counts_2x2 = if (binary_convention) "binary test result (1 = alert); no numeric cutoff was supplied" else
        "score of the test that produced the counts",
      theoretical_binormal = "binormal score standardised to the non-event distribution (z)")
    if (is.null(m) || is.null(m$score_meaning)) defaults <- c(defaults, "score_meaning")
    tie_rule <- if (cmp %in% c(">=", "<=")) "scores_equal_to_cutoff_alert" else
      "scores_equal_to_cutoff_do_not_alert"
    rec <- list(
      model_id = model_ids[i],
      name = if (is.null(m)) model_ids[i] else m$name,
      version = if (is.null(m)) NULL else m$version,
      score_meaning = if (is.null(m) || is.null(m$score_meaning)) default_meaning else m$score_meaning,
      score_scale = if (input_kind == "scored_data") unname(scales[[i]]) else "raw_score",
      positive_label = switch(input_kind,
                              scored_data = as.character(positive_label),
                              counts_2x2 = "event", theoretical_binormal = "event"),
      negative_label = if (input_kind == "scored_data") prep$negative_label else "non_event",
      direction = d, comparator = cmp, tie_rule = tie_rule,
      evaluation_split = split,
      missing_value_policy = if (input_kind == "scored_data") missing else "not_applicable",
      population = if (is.null(m)) NULL else m$population,
      unit = if (is.null(m)) NULL else m$unit,
      source_hash = if (is.null(m)) NULL else m$source_hash,
      score_column = if (input_kind == "scored_data") unname(scores[[i]]) else NULL,
      truth_column = if (input_kind == "scored_data") truth else NULL,
      weight_column = weights, cluster_column = cluster, id_column = id,
      weight_semantics = if (is.null(weights)) "none" else "frequency",
      defaults_applied = if (length(defaults)) sort(defaults) else NULL
    )
    records[[i]] <- rec
  }
  if (length(p)) .ig_abort_invalid(unique(p), "invalid classification specification")
  records
}

.ig_cls_prevalence_record <- function(input_kind, prevalence, prevalence_source, sampling) {
  p <- character()
  if (!is.null(prevalence) && !(is.numeric(prevalence) && length(prevalence) == 1L &&
                                is.finite(prevalence) && prevalence >= 0 && prevalence <= 1)) {
    p <- c(p, "prevalence must be a single probability in [0, 1]")
  }
  if (input_kind == "theoretical_binormal" && is.null(prevalence)) {
    p <- c(p, "a theoretical model needs an assumed `prevalence`")
  }
  if (!is.null(sampling) && !(.ig_is_string(sampling) && sampling %in% .ig_cls_sampling[1:3])) {
    p <- c(p, "sampling must be 'cross_sectional', 'case_control' or 'unknown'")
  }
  if (input_kind == "theoretical_binormal" && !is.null(sampling)) {
    p <- c(p, "a theoretical model has no sampling design")
  }
  if (is.null(prevalence)) {
    if (!is.null(prevalence_source) && !identical(prevalence_source, "sample")) {
      p <- c(p, "prevalence_source 'assumed' or 'external' needs a `prevalence` value")
    }
  } else if (!is.null(prevalence_source) &&
             !(.ig_is_string(prevalence_source) && prevalence_source %in% c("assumed", "external"))) {
    p <- c(p, "the source of a scenario prevalence must be 'assumed' or 'external'")
  }
  if (length(p)) .ig_abort_invalid(p, "invalid prevalence settings")
  samp <- if (input_kind == "theoretical_binormal") "not_applicable" else sampling %||% "unknown"
  if (is.null(prevalence)) {
    list(kind = "observed", value = NULL, source = "sample", sampling = samp)
  } else {
    list(kind = "scenario", value = as.double(prevalence),
         source = prevalence_source %||% "assumed", sampling = samp,
         defaults_applied = if (is.null(prevalence_source)) "source" else NULL)
  }
}

.ig_cls_theory_record <- function(theory) {
  if (!is.list(theory)) .ig_abort("invalid_input", "`theory` must be a list with auc, sigma and optionally n_total")
  bad <- setdiff(names(theory), c("auc", "sigma", "n_total"))
  if (length(bad) || is.null(names(theory))) {
    .ig_abort("invalid_input", "`theory` accepts only auc, sigma and n_total")
  }
  defaults <- NULL
  n_total <- theory$n_total
  if (is.null(n_total)) {
    n_total <- 1000
    defaults <- "n_total"
  }
  p <- character()
  num1 <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!num1(theory$auc) || theory$auc <= 0 || theory$auc >= 1) p <- c(p, "theory$auc must be in (0, 1)")
  if (!num1(theory$sigma) || theory$sigma <= 0) p <- c(p, "theory$sigma must be > 0")
  if (!num1(n_total) || n_total <= 0 || n_total > .ig_cls_max_total) {
    p <- c(p, "theory$n_total must be > 0 and at most 1e100")
  }
  if (length(p)) .ig_abort_invalid(p, "invalid theoretical model")
  list(auc = as.double(theory$auc), sigma = as.double(theory$sigma), n_total = as.double(n_total),
       model = "binormal_S0_N(0,1)_S1_N(d,sigma^2)", defaults_applied = defaults)
}

# Data preparation --------------------------------------------------------------------

.ig_cls_prepare_data <- function(data, truth, positive_label, negative_label, score_cols,
                                 weights, cluster, id, missing, score_scale) {
  p <- character()
  col_ok <- function(x) .ig_is_string(x) && x %in% names(data)
  if (!col_ok(truth)) p <- c(p, "`truth` must name a column of `data`")
  if (is.null(positive_label) || length(positive_label) != 1L || is.na(positive_label)) {
    p <- c(p, "`positive_label` is required (the positive class is never inferred)")
  }
  for (s in score_cols) if (!col_ok(s)) p <- c(p, sprintf("score column '%s' is not in `data`", s))
  for (nm in c("weights", "cluster", "id")) {
    v <- get(nm)
    if (!is.null(v) && !col_ok(v)) p <- c(p, sprintf("`%s` must name a column of `data`", nm))
  }
  if (!is.null(negative_label) && (length(negative_label) != 1L || is.na(negative_label))) {
    p <- c(p, "`negative_label` must be a single value")
  }
  if (length(p)) .ig_abort_invalid(p, "invalid scored data")
  scales <- .ig_cls_score_scales(score_scale,names(score_cols))
  n_input <- nrow(data)
  if (n_input == 0L) .ig_abort("invalid_input", "`data` has no rows")
  y_raw <- data[[truth]]
  if (!is.atomic(y_raw) || is.matrix(y_raw)) .ig_abort("invalid_input", "the truth column must be an atomic vector")
  y_chr <- as.character(y_raw)
  sc <- lapply(score_cols, function(s) data[[s]])
  for (j in seq_along(sc)) {
    if (!is.numeric(sc[[j]]) || !is.null(dim(sc[[j]]))) p <- c(p, sprintf("score column '%s' must be a numeric vector", score_cols[j]))
  }
  w <- if (is.null(weights)) NULL else data[[weights]]
  if (!is.null(w) && !is.numeric(w)) p <- c(p, "the weight column must be numeric")
  cl <- if (is.null(cluster)) NULL else data[[cluster]]
  ids <- if (is.null(id)) NULL else data[[id]]
  if (length(p)) .ig_abort_invalid(p, "invalid scored data")

  # Missing values, matched across models.
  miss <- is.na(y_chr)
  miss_by <- list(truth = sum(is.na(y_chr)))
  for (j in seq_along(sc)) {
    mj <- is.na(sc[[j]])
    miss_by[[paste0("score:", names(score_cols)[j])]] <- sum(mj)
    miss <- miss | mj
  }
  if (!is.null(w)) { miss_by$weights <- sum(is.na(w)); miss <- miss | is.na(w) }
  if (!is.null(cl)) { miss_by$cluster <- sum(is.na(cl)); miss <- miss | is.na(cl) }
  if (any(miss) && missing == "error") {
    where <- names(miss_by)[unlist(miss_by) > 0]
    .ig_abort_invalid(sprintf("%d missing value(s) in %s", unlist(miss_by)[where], where),
                      "missing values found (use missing = 'exclude' to drop those rows for all models)")
  }
  keep <- !miss
  n_excluded <- sum(miss)
  if (!any(keep)) .ig_abort("invalid_input", "no complete rows remain after excluding missing values")
  y_chr <- y_chr[keep]
  sc <- lapply(sc, function(s) as.double(s[keep]))
  if (!is.null(w)) w <- as.double(w[keep])
  if (!is.null(cl)) cl <- as.character(cl[keep])
  if (!is.null(ids)) {
    if (anyNA(ids)) p <- c(p, "the id column has missing values")
    else if (anyDuplicated(ids)) {
      p <- c(p, sprintf("duplicate unit ids (%d); declare repeated units with `cluster` and give each row its own id",
                        sum(duplicated(ids))))
    }
    ids <- as.character(ids[keep])
  }

  # Truth labels: binary, positive label explicit.
  pos <- as.character(positive_label)
  levels_seen <- sort(unique(y_chr))
  neg <- if (is.null(negative_label)) NULL else as.character(negative_label)
  if (!is.null(neg)) {
    if (identical(neg, pos)) p <- c(p, "negative_label must differ from positive_label")
    other <- setdiff(levels_seen, c(pos, neg))
    if (length(other)) p <- c(p, sprintf("truth contains values other than the declared labels: %s",
                                         paste(utils::head(other, 5), collapse = ", ")))
  } else {
    other <- setdiff(levels_seen, pos)
    if (length(other) > 1L) {
      p <- c(p, sprintf("truth must be binary; found %d values (%s)", length(levels_seen),
                        paste(utils::head(levels_seen, 6), collapse = ", ")))
    } else if (length(other) == 1L) {
      neg <- other
    }
    if (!pos %in% levels_seen && length(levels_seen) > 1L) {
      p <- c(p, sprintf("positive_label '%s' does not occur in truth (values: %s)", pos,
                        paste(levels_seen, collapse = ", ")))
    }
  }
  # Scores.
  for (j in seq_along(sc)) {
    s <- sc[[j]]
    if (any(!is.finite(s))) p <- c(p, sprintf("score column '%s' has non-finite values", score_cols[j]))
    else if (scales[[j]] == "probability" && any(s < 0 | s > 1)) {
      p <- c(p, sprintf("score column '%s' is declared as probabilities but has values outside [0, 1]",
                        score_cols[j]))
    }
  }
  # Frequency weights.
  if (!is.null(w)) {
    if (any(!is.finite(w)) || any(w < 0)) p <- c(p, "frequency weights must be finite and non-negative")
    else if (any(w != round(w))) {
      .ig_abort("unsupported", paste(
        "non-integer weights are not frequency weights; survey, importance or probability weights",
        "are not supported because the binomial and DeLong intervals would treat them as",
        "independent replicated units"))
    } else if (sum(w) <= 0) {
      p <- c(p, "the frequency weights sum to zero")
    } else if (sum(w) > .ig_cls_max_total) {
      p <- c(p, "the frequency weights must not exceed 1e100 in total")
    }
  }
  if (length(p)) .ig_abort_invalid(p, "invalid scored data")
  y <- as.integer(y_chr == pos)
  wt <- if (is.null(w)) rep(1, length(y)) else w
  n_pos <- sum(wt[y == 1L])
  n_neg <- sum(wt[y == 0L])
  list(y = y, scores = sc, w = w, cluster = cl, id = ids,
       positive_label = pos, negative_label = neg,
       n_input = n_input, n_excluded = n_excluded,
       excluded_by = miss_by[unlist(miss_by) > 0],
       n_events = n_pos, n_non_events = n_neg,
       single_class = (n_pos == 0 || n_neg == 0))
}

# SHA-256 fingerprint of the analysed data (truth, scores, weights, clusters,
# ids after exclusions). Deterministic text, built without canonical JSON for
# speed; strings are length-prefixed so that values cannot run together.
.ig_cls_fingerprint <- function(prep) {
  num <- function(x) paste(.ig_format_double(as.double(x)), collapse = ",")
  str <- function(x) paste0(nchar(enc2utf8(x), type = "bytes"), ":", enc2utf8(x), collapse = ",")
  parts <- c("inferguidR.classification_data/1",
             paste0("truth=", paste(prep$y, collapse = ",")),
             vapply(names(prep$scores), function(m) paste0("score[", m, "]=", num(prep$scores[[m]])), ""),
             paste0("weights=", if (is.null(prep$w)) "none" else num(prep$w)),
             paste0("cluster=", if (is.null(prep$cluster)) "none" else str(prep$cluster)),
             paste0("id=", if (is.null(prep$id)) "none" else str(prep$id)))
  .ig_sha256(paste(parts, collapse = "\n"))
}

# Prepared data of a specification (error if not attached).
.ig_cls_spec_data <- function(spec) {
  prep <- attr(spec, .ig_cls_data_attr, exact = TRUE)
  if (is.null(prep)) {
    .ig_abort("invalid_input", paste(
      "the scored data are not attached to this specification (specifications read from JSON",
      "carry only a fingerprint); reattach them with attach_classification_data()"))
  }
  prep
}

#' Reattach scored data to a classification specification
#'
#' Specifications serialised to JSON carry a SHA-256 fingerprint of the
#' analysed data, never the rows. To evaluate such a specification again,
#' reattach the original data: the columns named in the model records are read
#' and the fingerprint must match exactly.
#'
#' @param spec A `classification_spec` with input kind `scored_data`.
#' @param data The original data frame.
#' @return The specification with the data attached.
#' @examples
#' d <- data.frame(y = c(1, 1, 0, 0, 0), s = c(3, 2, 2, 1, 0))
#' spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "s", cutoff = 2)
#' back <- ig_from_json(ig_to_json(analysis_envelope(spec)))$spec
#' spec2 <- attach_classification_data(back, d)
#' identical(spec2$spec_hash, spec$spec_hash)
#' @concept classification
#' @seealso
#' * [classification_spec()]
#' * [evaluate_cutoff()]
#' * [plot_classification()]
#' @export
attach_classification_data <- function(spec, data) {
  .ig_cls_check_spec(spec)
  if (!identical(spec$input_kind, "scored_data")) {
    .ig_abort("invalid_input", "only scored-data specifications carry data")
  }
  m1 <- spec$models[[1]]
  score_cols <- stats::setNames(vapply(spec$models, function(m) m$score_column, ""),
                                vapply(spec$models, function(m) m$model_id, ""))
  prep <- .ig_cls_prepare_data(data, m1$truth_column, m1$positive_label,
                               negative_label = m1$negative_label, score_cols,
                               m1$weight_column, m1$cluster_column, m1$id_column,
                               m1$missing_value_policy, stats::setNames(vapply(spec$models,function(m) m$score_scale,""),names(score_cols)))
  if (!identical(.ig_cls_fingerprint(prep), spec$data_fingerprint)) {
    .ig_abort("invalid_input", "the data do not match the specification's data fingerprint")
  }
  attr(spec, .ig_cls_data_attr) <- prep
  spec
}

# Accept a classification_spec (possibly rebuilt from JSON as a plain list
# with the right class) and verify its integrity.
.ig_cls_check_spec <- function(spec) {
  if (!inherits(spec, "classification_spec")) {
    .ig_abort("invalid_input", "`spec` must be a classification_spec (see classification_spec())")
  }
  .ig_check_aux_spec(spec, "classification")
  invisible(spec)
}

#' Update a classification specification
#'
#' Specifications are immutable; `update()` builds a new, fully validated
#' specification. For a theoretical model, changing `auc` or `sigma` keeps the
#' numerical cutoff unchanged.
#'
#' @param object A `classification_spec`.
#' @param ... Not used.
#' @param rule,cutoff A new operating rule (see [cutoff_rule()]) or a new
#'   manual cutoff (which keeps the declared selection data).
#' @param prevalence,prevalence_source,sampling New prevalence settings; use
#'   `prevalence = NA` to return to the observed prevalence.
#' @param theory For theoretical models: a list with any of `auc`, `sigma`,
#'   `n_total`.
#' @param evaluation_split New evaluation split for all models.
#' @param uncertainty New settings from [classification_uncertainty()].
#' @return A new `classification_spec`.
#' @examples
#' th <- classification_spec(theory = list(auc = 0.8, sigma = 1), cutoff = 1.2, prevalence = 0.1)
#' th2 <- update(th, theory = list(auc = 0.9))
#' th2$operating_rule$score_cutoff
#' @concept classification
#' @seealso
#' * [classification_spec()]
#' * [evaluate_cutoff()]
#' * [plot_classification()]
#' @export
update.classification_spec <- function(object, ..., rule = NULL, cutoff = NULL, prevalence = NULL,
                                       prevalence_source = NULL, sampling = NULL, theory = NULL,
                                       evaluation_split = NULL, uncertainty = NULL) {
  .ig_cls_check_spec(object)
  if (length(list(...))) .ig_abort("invalid_input", "unknown arguments to update()")
  spec <- object
  prep <- attr(spec, .ig_cls_data_attr, exact = TRUE)
  if (spec$input_kind == "scored_data" && is.null(prep)) .ig_cls_spec_data(spec)
  op <- spec$operating_rule
  if (!is.null(rule) && !is.null(cutoff)) .ig_abort("invalid_input", "give either `rule` or `cutoff`, not both")
  if (!is.null(cutoff)) {
    rule <- cutoff_rule(cutoff = cutoff, selection_data = op$selection_data)
  }
  if (!is.null(rule)) {
    if (!inherits(rule, "ig_cutoff_rule")) .ig_abort("invalid_input", "`rule` must come from cutoff_rule()")
    op <- unclass(rule)
  } else {
    op$cutoff_scale <- NULL
  }
  prev <- spec$prevalence
  if (!is.null(prevalence) || !is.null(prevalence_source) || !is.null(sampling)) {
    pv <- if (is.null(prevalence)) prev$value else if (length(prevalence) == 1L && is.na(prevalence)) NULL else prevalence
    src <- if (!is.null(prevalence_source)) prevalence_source else if (is.null(prevalence)) {
      if (identical(prev$kind, "scenario") && !"source" %in% prev$defaults_applied) prev$source else NULL
    } else NULL
    smp <- sampling %||% (if (prev$sampling %in% c("unknown", "not_applicable")) NULL else prev$sampling)
    prev <- .ig_cls_prevalence_record(spec$input_kind, pv, src, smp)
  }
  th <- spec$theory
  if (!is.null(theory)) {
    if (spec$input_kind != "theoretical_binormal") .ig_abort("invalid_input", "`theory` applies to theoretical models only")
    keep_n <- if ("n_total" %in% th$defaults_applied) NULL else th$n_total
    th <- .ig_cls_theory_record(list(auc = theory$auc %||% th$auc, sigma = theory$sigma %||% th$sigma,
                                     n_total = theory$n_total %||% keep_n))
  }
  models <- spec$models
  if (!is.null(evaluation_split)) {
    if (spec$input_kind == "theoretical_binormal") .ig_abort("invalid_input", "a theoretical model has no evaluation split")
    splits <- c("same_data", "held_out_tuning", "held_out_test", "external", "resampling")
    if (!(.ig_is_string(evaluation_split) && evaluation_split %in% splits)) {
      .ig_abort("invalid_input", sprintf("evaluation_split must be one of: %s", paste(splits, collapse = ", ")))
    }
    models <- lapply(models, function(m) {
      m$evaluation_split <- evaluation_split
      d <- setdiff(m$defaults_applied, "evaluation_split")
      m["defaults_applied"] <- list(if (length(d)) d else NULL)
      m
    })
  }
  unc <- spec$uncertainty
  if (!is.null(uncertainty)) {
    if (!inherits(uncertainty, "ig_cls_uncertainty")) {
      .ig_abort("invalid_input", "`uncertainty` must come from classification_uncertainty()")
    }
    unc <- unclass(uncertainty)
  }
  .ig_cls_build_spec(spec$input_kind, models, op, prev, spec$counts, th, prep, unc)
}

#' @export
format.classification_spec <- function(x, ...) {
  kind <- x$input_kind
  op <- x$operating_rule
  m <- x$models
  rule_txt <- if (identical(op$selection, "manual")) {
    sprintf("manual cutoff %s%s", format(op$score_cutoff, digits = 7),
            if (identical(op$cutoff_source, "binary_test_convention")) " (binary test convention)" else "")
  } else {
    sprintf("selection %s", op$selection)
  }
  model_txt <- vapply(m, function(r) {
    sprintf("%s: alert if score %s cutoff (%s; %s)", r$model_id, r$comparator, r$direction,
            r$score_scale)
  }, "")
  defaults <- unique(unlist(lapply(m, `[[`, "defaults_applied")))
  prev <- x$prevalence
  prev_txt <- if (identical(prev$kind, "scenario")) sprintf("scenario %s (%s)", format(prev$value), prev$source) else
    sprintf("observed in the sample (sampling: %s)", prev$sampling)
  c(sprintf("<classification_spec> %s", kind),
    paste0("  model    ", model_txt),
    sprintf("  rule:     %s; selection data: %s", rule_txt, op$selection_data),
    sprintf("  prevalence: %s", prev_txt),
    if (kind == "theoretical_binormal") sprintf("  theory:   AUC %s, sigma %s, N %s",
                                               format(x$theory$auc), format(x$theory$sigma),
                                               format(x$theory$n_total)),
    if (kind == "counts_2x2") sprintf("  counts:   TP %s, FN %s, FP %s, TN %s", x$counts$tp, x$counts$fn,
                                     x$counts$fp, x$counts$tn),
    if (kind == "scored_data") sprintf("  data:     fingerprint %s%s", substr(x$data_fingerprint, 1, 16),
                                      if (is.null(attr(x, .ig_cls_data_attr, exact = TRUE))) " (data not attached)" else ""),
    if (length(defaults)) sprintf("  defaults applied (confirm): %s", paste(defaults, collapse = ", ")),
    sprintf("  spec hash: %s", x$spec_hash))
}

#' @export
print.classification_spec <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}
