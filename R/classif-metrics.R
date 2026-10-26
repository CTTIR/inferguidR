# Classification metrics: dictionary and confusion-count engine ---------------------
#
# The metric dictionary (inst/metadata/metric-dictionary.json) documents every
# metric: stable id, EN/DE labels, the formula actually computed, denominator,
# domain, prevalence dependence, preferred direction and the typed statuses it
# takes when undefined or infinite. This file implements each operating-point
# metric directly from the four cells, vectorised over decision sets.
#
# Status vocabulary of a metric value:
#   ok                finite value
#   undefined         mathematically undefined (value NA, reason code given)
#   positive_infinity mathematically +Inf (value NA, never a silent number)
#   not_applicable    not applicable to the input or design
#   not_computed      not requested or no method implemented
# No division convention is silent: a zero denominator always yields a typed
# status with a reason code from the dictionary.

.ig_cls_metric_statuses <- c("ok", "undefined", "positive_infinity", "not_applicable",
                             "not_computed")

# Cache for the parsed dictionary (package-private, filled on first use).
.ig_cls_cache <- new.env(parent = emptyenv())

.ig_cls_dictionary_path <- function() {
  path <- system.file("metadata", "metric-dictionary.json", package = "inferguidR")
  if (!nzchar(path)) stop("the metric dictionary (inst/metadata/metric-dictionary.json) is missing", call. = FALSE)
  path
}

# Parsed dictionary document (plain lists).
.ig_cls_dictionary <- function() {
  doc <- .ig_cls_cache$dictionary
  if (is.null(doc)) {
    path <- .ig_cls_dictionary_path()
    txt <- readChar(path, file.size(path), useBytes = TRUE)
    Encoding(txt) <- "UTF-8"
    doc <- .ig_json_parse(txt)
    if (!identical(doc$schema_version, "inferguidR.metric_dictionary/0.1.0")) {
      stop("unsupported metric dictionary schema_version", call. = FALSE)
    }
    names(doc$metrics) <- vapply(doc$metrics, `[[`, "", "id")
    .ig_cls_cache$dictionary <- doc
  }
  doc
}

.ig_cls_metric_entry <- function(id) {
  entry <- .ig_cls_dictionary()$metrics[[id]]
  if (is.null(entry)) stop("unknown metric id: ", id, call. = FALSE)
  entry
}

# Reason text for a reason code (English or German).
.ig_cls_reason_text <- function(code, lang = "en") {
  codes <- .ig_cls_dictionary()$reason_codes
  vapply(code, function(cd) {
    if (is.na(cd)) return(NA_character_)
    rec <- codes[[cd]]
    if (is.null(rec)) cd else rec[[lang]] %||% rec$en
  }, character(1), USE.NAMES = FALSE)
}

# Operating-point metrics in dictionary order (the ones computed from counts).
.ig_cls_operating_ids <- function() {
  m <- .ig_cls_dictionary()$metrics
  ids <- names(m)[vapply(m, function(e) identical(e$scope, "operating_point"), logical(1))]
  ids
}

#' Metric dictionary for binary classification
#'
#' Returns the documented definition of every classification metric that
#' inferguidR computes: stable id, label, the formula actually computed, its
#' denominator, the valid domain, whether the metric depends on prevalence,
#' the preferred direction and the typed statuses it can take when it is
#' undefined or infinite. The dictionary is bundled as
#' `inst/metadata/metric-dictionary.json`; formulas are documentation and are
#' never evaluated.
#'
#' A metric that is "invariant given fixed class-conditional distributions"
#' does not change when only the prevalence changes while the score
#' distributions within the true classes stay fixed. That is an assumption of
#' prevalence scenarios, not a claim about real population shift.
#'
#' @param lang Label language, `"en"` or `"de"`.
#' @param scope Optional filter: `"operating_point"`, `"curve"`,
#'   `"probability"` or `"decision"`.
#' @return A data frame with one row per metric and the columns `id`,
#'   `short`, `name`, `scope`, `group`, `family`, `formula`, `denominator`,
#'   `unit`, `lower`, `upper`, `prevalence_dependence`,
#'   `prevalence_dependent`, `preferred_direction`, `complement`, `inputs`,
#'   `undefined_statuses`, `interval_methods`, `source_feature_ids` and
#'   `correction`. List-like fields are collapsed with `"; "`.
#' @examples
#' d <- metric_dictionary()
#' d[d$id %in% c("ppv", "lr_pos"), c("id", "formula", "undefined_statuses")]
#' @concept classification
#' @seealso
#' * [classification_spec()]
#' * [evaluate_cutoff()]
#' * [plot_classification()]
#' @export
metric_dictionary <- function(lang = c("en", "de"), scope = NULL) {
  lang <- match.arg(lang)
  m <- .ig_cls_dictionary()$metrics
  if (!is.null(scope)) {
    if (!.ig_is_string(scope) || !scope %in% c("operating_point", "curve", "probability", "decision")) {
      .ig_abort("invalid_input", "`scope` must be one of: operating_point, curve, probability, decision")
    }
    m <- m[vapply(m, function(e) identical(e$scope, scope), logical(1))]
  }
  str_or_na <- function(x) if (is.null(x)) NA_character_ else as.character(x)
  collapse <- function(x) if (!length(x)) "" else paste(unlist(x), collapse = "; ")
  rows <- lapply(m, function(e) {
    lab <- e$label[[lang]] %||% e$label$en
    den <- e$denominator
    data.frame(
      id = e$id,
      short = lab$short,
      name = lab$name,
      scope = e$scope,
      group = e$group,
      family = str_or_na(e$family),
      formula = e$formula,
      denominator = if (is.null(den)) NA_character_ else
        sprintf("%s (%s)", den$expression, den[[lang]] %||% den$en),
      unit = e$unit,
      lower = if (is.null(e$domain$lower)) NA_real_ else as.double(e$domain$lower),
      upper = if (is.null(e$domain$upper)) NA_real_ else as.double(e$domain$upper),
      prevalence_dependence = e$prevalence_dependence,
      prevalence_dependent = !e$prevalence_dependence %in%
        "invariant_given_fixed_class_conditional_distributions",
      preferred_direction = e$preferred_direction,
      complement = str_or_na(e$complement),
      inputs = collapse(e$inputs),
      undefined_statuses = collapse(vapply(e$undefined, function(u) {
        paste0(u$status, ":", u$reason_code)
      }, "")),
      interval_methods = collapse(e$interval_methods),
      source_feature_ids = collapse(e$source_feature_ids),
      correction = str_or_na(e$correction),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Vectorised engine --------------------------------------------------------------

# Apply ordered status conditions to a vector of raw values. Each condition is
# list(mask, status, reason_code); the first condition that matches wins.
# Values whose status is not "ok" become NA.
#
# `range_mask` marks rows with a non-zero cell below .ig_cls_tiny_cell (only
# far tails of a theoretical model or user-supplied expected counts produce
# them). There a mathematically finite ratio can leave the double range
# (overflow, or underflow of both of its terms); such a value is reported as
# "not_computed" with reason "outside_double_range", never as a number or as
# infinity. Any other non-finite value that no condition explains is a defect
# of this engine and stops with an error: with all non-zero cells at least
# .ig_cls_tiny_cell and totals at most .ig_cls_max_total no formula below can
# leave the double range.
.ig_cls_finish <- function(value, conds = list(), range_mask = FALSE) {
  k <- length(value)
  status <- rep("ok", k)
  reason <- rep(NA_character_, k)
  for (cd in conds) {
    m <- rep_len(cd[[1]], k)
    m[is.na(m)] <- FALSE
    m <- m & status == "ok"
    status[m] <- cd[[2]]
    reason[m] <- cd[[3]]
  }
  out_of_range <- status == "ok" & !is.finite(value) & rep_len(range_mask, k)
  status[out_of_range] <- "not_computed"
  reason[out_of_range] <- "outside_double_range"
  value[status != "ok"] <- NA_real_
  if (any(status == "ok" & !is.finite(value))) {
    stop("internal error: unexplained non-finite metric value", call. = FALSE)
  }
  list(value = as.double(value), status = status, reason = reason)
}

# All operating-point metrics for vectors of cell counts (observed, weighted
# or expected; expected counts may be fractional). Returns a named list (by
# metric id) of list(value, status, reason).
#
# A row with a missing cell holds undefined expected counts (a prevalence
# scenario of a sample without events or without non-events transports
# nothing); every metric of such a row is undefined with reason
# "counts_undefined" instead of reaching the engine below.
.ig_cls_metrics_vec <- function(tp, fn, fp, tn) {
  tp <- as.double(tp); fn <- as.double(fn); fp <- as.double(fp); tn <- as.double(tn)
  missing_cells <- is.na(tp) | is.na(fn) | is.na(fp) | is.na(tn)
  if (any(missing_cells)) {
    k <- length(tp)
    keep <- !missing_cells
    part <- .ig_cls_metrics_vec(tp[keep], fn[keep], fp[keep], tn[keep])
    return(lapply(part, function(m) {
      value <- rep(NA_real_, k)
      status <- rep("undefined", k)
      reason <- rep("counts_undefined", k)
      value[keep] <- m$value
      status[keep] <- m$status
      reason[keep] <- m$reason
      list(value = value, status = status, reason = reason)
    }))
  }
  P <- tp + fn
  N <- fp + tn
  A <- tp + fp
  Q <- fn + tn
  n <- P + N
  no_ev <- list(P == 0, "undefined", "no_events")
  no_ne <- list(N == 0, "undefined", "no_non_events")
  no_al <- list(A == 0, "undefined", "no_alerts")
  no_na <- list(Q == 0, "undefined", "no_non_alerts")
  no_n <- list(n == 0, "undefined", "no_units")
  # rows with a tiny non-zero cell, where a finite ratio can leave the double range
  tiny <- function(x) x > 0 & x < .ig_cls_tiny_cell
  range_mask <- tiny(tp) | tiny(fn) | tiny(fp) | tiny(tn)
  fin <- function(v, conds = list()) .ig_cls_finish(v, conds, range_mask)
  ok <- function(v) fin(v)
  f_beta <- function(b2) {
    fin(((1 + b2) * tp) / ((1 + b2) * tp + b2 * fn + fp),
        list(list(tp + fn + fp == 0, "undefined", "no_tp_fp_fn")))
  }
  tpr <- tp / P
  fpr <- fp / N
  out <- list(
    tp = ok(tp), fn = ok(fn), fp = ok(fp), tn = ok(tn),
    n_events = ok(P), n_non_events = ok(N), n_alerts = ok(A), n_non_alerts = ok(Q),
    n_total = ok(n),
    prevalence = fin(P / n, list(no_n)),
    alert_rate = fin(A / n, list(no_n)),
    tpr = fin(tpr, list(no_ev)),
    fnr = fin(fn / P, list(no_ev)),
    tnr = fin(tn / N, list(no_ne)),
    fpr = fin(fpr, list(no_ne)),
    ppv = fin(tp / A, list(no_al)),
    fdr = fin(fp / A, list(no_al)),
    npv = fin(tn / Q, list(no_na)),
    for_rate = fin(fn / Q, list(no_na)),
    accuracy = fin((tp + tn) / n, list(no_n)),
    balanced_accuracy = fin((tp / P + tn / N) / 2, list(no_ev, no_ne)),
    youden_j = fin(tp / P - fp / N, list(no_ev, no_ne)),
    markedness = fin(tp / A - fn / Q, list(no_al, no_na)),
    f1 = f_beta(1),
    f05 = f_beta(0.25),
    f2 = f_beta(4),
    fowlkes_mallows = fin(tp / (sqrt(A) * sqrt(P)), list(no_ev, no_al)),
    jaccard = fin(tp / (tp + fn + fp), list(list(tp + fn + fp == 0, "undefined", "no_tp_fp_fn"))),
    # Square roots taken separately so that large counts cannot overflow.
    mcc = fin((tp * tn - fp * fn) / (sqrt(A) * sqrt(P) * sqrt(N) * sqrt(Q)),
              list(no_ev, no_ne, no_al, no_na)),
    # Exact count form of (accuracy - pe) / (1 - pe). The chance agreement pe
    # is 1 exactly when (A = 0 or N = 0) and (P = 0 or Q = 0); the condition
    # is tested on the cells so that an underflowing product cannot trigger it.
    cohen_kappa = fin(2 * (tp * tn - fp * fn) / (A * N + P * Q),
                      list(list((A == 0 | N == 0) & (P == 0 | Q == 0), "undefined", "chance_agreement_one"))),
    lr_pos = fin((tp * N) / (fp * P),
                 list(no_ev, no_ne,
                      list(fp == 0 & tp > 0, "positive_infinity", "fpr_zero_tpr_positive"),
                      list(fp == 0 & tp == 0, "undefined", "tpr_and_fpr_zero"))),
    lr_neg = fin((fn * N) / (tn * P),
                 list(no_ev, no_ne,
                      list(tn == 0 & fn > 0, "positive_infinity", "tnr_zero_fnr_positive"),
                      list(tn == 0 & fn == 0, "undefined", "fnr_and_tnr_zero"))),
    # FP x FN = 0 and TP x TN > 0 (resp. = 0), tested on the cells.
    dor = fin((tp * tn) / (fp * fn),
              list(no_ev, no_ne,
                   list((fp == 0 | fn == 0) & tp > 0 & tn > 0, "positive_infinity", "fp_fn_product_zero"),
                   list((fp == 0 | fn == 0) & (tp == 0 | tn == 0), "undefined", "zero_over_zero_products"))),
    alerts_per_tp = fin(A / tp,
                        list(no_al, list(tp == 0 & A > 0, "positive_infinity", "tp_zero_alerts_positive"))),
    fp_per_tp = fin(fp / tp,
                    list(no_al, list(tp == 0 & fp > 0, "positive_infinity", "tp_zero_fp_positive"))),
    fp_per_1000_non_events = fin(1000 * fp / N, list(no_ne)),
    missed_per_1000_total = fin(1000 * fn / n, list(no_n)),
    tp_share = fin(tp / n, list(no_n)),
    fp_share = fin(fp / n, list(no_n))
  )
  pg <- .ig_cls_ppv_gap(tp / P, fp / N, P == 0, N == 0)
  out$ppv_gap_max_prevalence <- pg
  out[.ig_cls_operating_ids()]
}

# Prevalence maximising PPV(pi) - pi for fixed TPR = a and FPR = b
# (classification-source-audit section 4): pi* = sqrt(b) / (sqrt(a) + sqrt(b))
# for a > b > 0, with typed statuses for every other case.
.ig_cls_ppv_gap <- function(a, b, no_events = rep(FALSE, length(a)),
                            no_non_events = rep(FALSE, length(a))) {
  v <- sqrt(b) / (sqrt(a) + sqrt(b))
  .ig_cls_finish(v, list(
    list(no_events | is.na(a), "undefined", "no_events"),
    list(no_non_events | is.na(b), "undefined", "no_non_events"),
    list(a == 0 & b == 0, "undefined", "tpr_and_fpr_zero"),
    list(a == b, "undefined", "equal_tpr_fpr"),
    list(a < b, "not_applicable", "tpr_below_fpr"),
    list(b == 0, "undefined", "fpr_zero_supremum")
  ))
}

# Largest total (counts, frequency weights, model or scenario N) accepted.
# Metric formulas multiply two cells (for example TP x TN in kappa, the DOR
# and LR+); totals up to 1e100 keep every product far from double overflow.
.ig_cls_max_total <- 1e100

# Non-zero cells below this size can make a finite ratio leave the double range
# (see .ig_cls_finish()); with every non-zero cell at least this large and totals
# at most .ig_cls_max_total, products of two cells and their ratios stay within it.
.ig_cls_tiny_cell <- 1e-50

# Validate a set of four cell counts. `kind` is the count kind: observed and
# frequency-weighted counts must be whole numbers.
.ig_cls_check_counts <- function(x, kind = "observed") {
  if (is.list(x)) x <- unlist(x)
  if (!is.numeric(x) || is.null(names(x)) || !all(c("tp", "fn", "fp", "tn") %in% names(x))) {
    .ig_abort("invalid_input", "counts must be a named numeric vector or list with tp, fn, fp and tn")
  }
  x <- x[c("tp", "fn", "fp", "tn")]
  problems <- character()
  if (any(!is.finite(x)) || any(x < 0)) problems <- c(problems, "counts must be finite and non-negative")
  if (kind != "expected" && all(is.finite(x)) && any(x != round(x))) {
    problems <- c(problems, "observed counts must be whole numbers (fractional values are expected counts)")
  }
  if (all(is.finite(x)) && sum(x) <= 0) problems <- c(problems, "the counts must contain at least one unit")
  if (all(is.finite(x)) && sum(x) > .ig_cls_max_total) {
    problems <- c(problems, "the counts must not exceed 1e100 in total (larger values would overflow the metric formulas)")
  }
  if (length(problems)) .ig_abort_invalid(problems, "invalid 2x2 counts")
  stats::setNames(as.double(x), c("tp", "fn", "fp", "tn"))
}

.ig_cls_count_kinds <- c("observed", "frequency_weighted", "expected")

# Numerator/denominator pairs of the proportion metrics that get binomial
# intervals (function of the four cells).
.ig_cls_proportion_parts <- function(tp, fn, fp, tn) {
  P <- tp + fn; N <- fp + tn; A <- tp + fp; Q <- fn + tn; n <- P + N
  list(prevalence = c(P, n), alert_rate = c(A, n), tpr = c(tp, P), fnr = c(fn, P),
       tnr = c(tn, N), fpr = c(fp, N), ppv = c(tp, A), fdr = c(fp, A),
       npv = c(tn, Q), for_rate = c(fn, Q), accuracy = c(tp + tn, n),
       tp_share = c(tp, n), fp_share = c(fp, n))
}

# One metric table (data frame) for a single set of counts. `estimate_kind`
# labels the evidence (observed, frequency_weighted, expected_theoretical,
# expected_scenario, sample_descriptive).
.ig_cls_metric_table <- function(counts, estimate_kind, lang = "en") {
  vals <- .ig_cls_metrics_vec(counts[["tp"]], counts[["fn"]], counts[["fp"]], counts[["tn"]])
  ids <- names(vals)
  dict <- .ig_cls_dictionary()$metrics
  parts <- .ig_cls_proportion_parts(counts[["tp"]], counts[["fn"]], counts[["fp"]], counts[["tn"]])
  # Value of the single denominator named in the dictionary (NA for metrics
  # with two denominators, such as Youden's J, or none, such as the cells).
  tp <- counts[["tp"]]; fn <- counts[["fn"]]; fp <- counts[["fp"]]; tn <- counts[["tn"]]
  P <- tp + fn; N <- fp + tn; A <- tp + fp; Q <- fn + tn
  den_value <- vapply(ids, function(id) {
    if (!is.null(parts[[id]])) return(parts[[id]][2])
    v <- switch(id,
                f1 = 2 * tp + fp + fn,
                f05 = 1.25 * tp + 0.25 * fn + fp,
                f2 = 5 * tp + 4 * fn + fp,
                fowlkes_mallows = sqrt(A) * sqrt(P),
                jaccard = tp + fn + fp,
                mcc = sqrt(A) * sqrt(P) * sqrt(N) * sqrt(Q),
                cohen_kappa = A * N + P * Q,
                lr_pos = fp / N,
                lr_neg = tn / N,
                dor = fp * fn,
                fp_per_1000_non_events = N,
                missed_per_1000_total = P + N,
                alerts_per_tp = tp, fp_per_tp = tp,
                ppv_gap_max_prevalence = sqrt(tp / P) + sqrt(fp / N),
                NA_real_)
    if (is.finite(v)) v else NA_real_
  }, double(1))
  den_label <- vapply(ids, function(id) {
    d <- dict[[id]]$denominator
    if (is.null(d)) NA_character_ else sprintf("%s (%s)", d$expression, d[[lang]] %||% d$en)
  }, character(1))
  data.frame(
    metric_id = ids,
    value = vapply(vals, function(v) v$value, double(1)),
    status = vapply(vals, function(v) v$status, character(1)),
    reason_code = vapply(vals, function(v) v$reason, character(1)),
    reason = .ig_cls_reason_text(vapply(vals, function(v) v$reason, character(1)), lang),
    denominator = unname(den_label),
    denominator_value = unname(den_value),
    unit = vapply(ids, function(id) dict[[id]]$unit, character(1)),
    estimate_kind = estimate_kind,
    stringsAsFactors = FALSE, row.names = NULL
  )
}

#' Classification metrics from 2x2 counts
#'
#' Computes every operating-point metric of the [metric_dictionary()] from the
#' four cells TP, FN, FP and TN. Undefined values are `NA` with a typed status
#' and a reason; mathematical infinity has the status `"positive_infinity"`
#' (value `NA`), never a silent number. PPV at zero alerts is undefined, never
#' 1.
#'
#' @param x A named numeric vector or list with `tp`, `fn`, `fp` and `tn`, or
#'   a `classification_result` (its observed operating-point metrics are
#'   returned, one block per model).
#' @param count_kind For counts: `"observed"` (whole numbers, the default),
#'   `"frequency_weighted"` or `"expected"` (fractional expected counts, for
#'   example from [prevalence_scenario()]).
#' @param lang Language of the reason texts and denominator labels.
#' @return A data frame with the columns `metric_id`, `value`, `status`,
#'   `reason_code`, `reason`, `denominator`, `denominator_value`, `unit` and
#'   `estimate_kind` (plus `model_id` for results).
#' @examples
#' m <- classification_metrics(c(tp = 40, fn = 10, fp = 20, tn = 130))
#' m[m$metric_id %in% c("tpr", "ppv", "lr_pos", "mcc"), c("metric_id", "value", "status")]
#' # Zero alerts: PPV is undefined, not 1
#' z <- classification_metrics(c(tp = 0, fn = 5, fp = 0, tn = 95))
#' z[z$metric_id == "ppv", c("value", "status", "reason")]
#' @concept classification
#' @seealso
#' * [classification_spec()]
#' * [evaluate_cutoff()]
#' * [plot_classification()]
#' @export
classification_metrics <- function(x, count_kind = c("observed", "frequency_weighted", "expected"),
                                   lang = c("en", "de")) {
  lang <- match.arg(lang)
  if (inherits(x, "classification_result")) {
    if (!any(vapply(x$models, function(m) !is.null(m$metrics), logical(1)))) {
      .ig_abort("invalid_input", "the result holds no operating-point metrics (status is not ok)")
    }
    return(classification_table(x, "metrics", lang = lang))
  }
  count_kind <- match.arg(count_kind)
  counts <- .ig_cls_check_counts(x, count_kind)
  .ig_cls_metric_table(counts, count_kind, lang)
}

# Prevalence relations ------------------------------------------------------------

.ig_cls_check_rate <- function(x, name) {
  if (!is.numeric(x) || !length(x) || any(!is.na(x) & (x < 0 | x > 1 | !is.finite(x)))) {
    .ig_abort("invalid_input", sprintf("`%s` must be numeric values in [0, 1] (NA allowed)", name))
  }
  as.double(x)
}

#' Prevalence that maximises PPV minus prevalence
#'
#' For a fixed operating point with true-positive rate \eqn{a} and
#' false-positive rate \eqn{b}, PPV as a function of prevalence is
#' \eqn{PPV(\pi) = a\pi / (b + (a - b)\pi)}. It is concave, so it rises
#' fastest at \eqn{\pi \to 0} (slope LR+), not at an interior prevalence. The
#' expression \eqn{(\sqrt{ab} - b)/(a - b) = \sqrt{b}/(\sqrt{a} + \sqrt{b})},
#' shown by some tools as the prevalence "where PPV rises fastest", is the
#' prevalence that maximises the gap \eqn{PPV(\pi) - \pi}: the point where
#' the PPV curve has slope 1 and where \eqn{PPV = 1 - \pi}. This function
#' returns that value under its correct name, using the numerically stable
#' form.
#'
#' Domain: a value only for \eqn{a > b > 0}. \eqn{b = 0 < a}: the supremum 1 of
#' the gap is approached as \eqn{\pi \to 0} but not attained (status
#' `undefined`, reason `fpr_zero_supremum`). \eqn{a = b}: the gap is 0 for
#' every prevalence (`equal_tpr_fpr`). \eqn{a < b}: the expression would be a
#' minimiser (`not_applicable`, `tpr_below_fpr`).
#'
#' @param tpr,fpr True- and false-positive rates in \[0, 1\] (vectors of equal
#'   length, or length 1). `NA` gives an undefined result.
#' @return A data frame with `tpr`, `fpr`, `value`, `status`, `reason_code`
#'   and `reason`.
#' @examples
#' ppv_gap_max_prevalence(tpr = c(0.8, 0.5, 0.3, 0.6), fpr = c(0.1, 0.5, 0.4, 0))
#' @concept classification
#' @seealso
#' * [classification_spec()]
#' * [evaluate_cutoff()]
#' * [plot_classification()]
#' @export
ppv_gap_max_prevalence <- function(tpr, fpr) {
  tpr <- .ig_cls_check_rate(tpr, "tpr")
  fpr <- .ig_cls_check_rate(fpr, "fpr")
  k <- max(length(tpr), length(fpr))
  tpr <- rep_len(tpr, k)
  fpr <- rep_len(fpr, k)
  r <- .ig_cls_ppv_gap(tpr, fpr, is.na(tpr), is.na(fpr))
  data.frame(tpr = tpr, fpr = fpr, value = r$value, status = r$status,
             reason_code = r$reason, reason = .ig_cls_reason_text(r$reason),
             stringsAsFactors = FALSE)
}

#' Predictive values at a given prevalence (Bayes transform)
#'
#' Transports a fixed operating point to other prevalences: with TPR and FPR
#' held fixed (fixed class-conditional score distributions),
#' \deqn{PPV(\pi) = \frac{TPR\,\pi}{TPR\,\pi + FPR\,(1-\pi)}, \quad
#' NPV(\pi) = \frac{(1-FPR)(1-\pi)}{(1-FPR)(1-\pi) + (1-TPR)\pi}.}
#' The result is a scenario, not an observed measurement; real population
#' shifts can also change TPR and FPR.
#'
#' @param tpr,fpr Rates of one operating point (length 1).
#' @param prevalence Vector of prevalences in \[0, 1\].
#' @return A data frame with `prevalence`, `ppv`, `ppv_status`,
#'   `ppv_reason_code`, `npv`, `npv_status` and `npv_reason_code`.
#' @examples
#' # Worked example: TPR 20 %, FPR 1 %
#' predictive_values_at_prevalence(0.2, 0.01, c(0.1, 0.2, 0.31))
#' @concept classification
#' @seealso
#' * [classification_spec()]
#' * [evaluate_cutoff()]
#' * [plot_classification()]
#' @export
predictive_values_at_prevalence <- function(tpr, fpr, prevalence) {
  tpr <- .ig_cls_check_rate(tpr, "tpr")
  fpr <- .ig_cls_check_rate(fpr, "fpr")
  if (length(tpr) != 1L || length(fpr) != 1L) .ig_abort("invalid_input", "`tpr` and `fpr` must be single values")
  prevalence <- .ig_cls_check_rate(prevalence, "prevalence")
  if (anyNA(prevalence)) .ig_abort("invalid_input", "`prevalence` must not contain NA")
  pi <- prevalence
  rate_na <- list(list(is.na(tpr), "undefined", "no_events"),
                  list(is.na(fpr), "undefined", "no_non_events"))
  ppv_den <- tpr * pi + fpr * (1 - pi)
  ppv <- .ig_cls_finish(tpr * pi / ppv_den,
                        c(rate_na, list(list(ppv_den == 0, "undefined", "no_alerts"))))
  npv_den <- (1 - fpr) * (1 - pi) + (1 - tpr) * pi
  npv <- .ig_cls_finish((1 - fpr) * (1 - pi) / npv_den,
                        c(rate_na, list(list(npv_den == 0, "undefined", "no_non_alerts"))))
  data.frame(prevalence = pi, ppv = ppv$value, ppv_status = ppv$status,
             ppv_reason_code = ppv$reason, npv = npv$value, npv_status = npv$status,
             npv_reason_code = npv$reason, stringsAsFactors = FALSE)
}

#' Prevalence needed to reach a PPV target
#'
#' At fixed TPR and FPR, PPV increases with prevalence (for TPR > 0). The
#' smallest prevalence at which PPV reaches a target \eqn{t} is
#' \deqn{\pi_t = \frac{t\,FPR}{(1-t)\,TPR + t\,FPR}.}
#' PPV is at least \eqn{t} for every prevalence in \eqn{[\pi_t, 1]}. The target
#' is never reachable when TPR = 0; when FPR = 0 and TPR > 0, PPV is 1 at
#' every prevalence above 0, so no positive minimum exists.
#'
#' @param tpr,fpr Rates of one operating point (vectors, recycled).
#' @param target PPV target in (0, 1). No default: targets are user decisions.
#' @return A data frame with `tpr`, `fpr`, `target`, `prevalence_needed`,
#'   `status`, `reason_code` and `reason`.
#' @examples
#' prevalence_for_ppv_target(tpr = c(0.2, 0.04, 0), fpr = c(0.01, 0.10, 0.1), target = 0.9)
#' @concept classification
#' @seealso
#' * [classification_spec()]
#' * [evaluate_cutoff()]
#' * [plot_classification()]
#' @export
prevalence_for_ppv_target <- function(tpr, fpr, target) {
  tpr <- .ig_cls_check_rate(tpr, "tpr")
  fpr <- .ig_cls_check_rate(fpr, "fpr")
  if (!is.numeric(target) || length(target) != 1L || !is.finite(target) || target <= 0 || target >= 1) {
    .ig_abort("invalid_input", "`target` must be a single PPV in (0, 1)")
  }
  k <- max(length(tpr), length(fpr))
  tpr <- rep_len(tpr, k)
  fpr <- rep_len(fpr, k)
  need <- target * fpr / ((1 - target) * tpr + target * fpr)
  r <- .ig_cls_finish(need, list(
    list(is.na(tpr), "undefined", "no_events"),
    list(is.na(fpr), "undefined", "no_non_events"),
    list(tpr == 0, "undefined", "tpr_zero_unreachable"),
    list(fpr == 0, "undefined", "fpr_zero_any_positive_prevalence")
  ))
  data.frame(tpr = tpr, fpr = fpr, target = target, prevalence_needed = r$value,
             status = r$status, reason_code = r$reason, reason = .ig_cls_reason_text(r$reason),
             stringsAsFactors = FALSE)
}
