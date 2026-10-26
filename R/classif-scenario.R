# Prevalence scenarios ----------------------------------------------------------------------
#
# A prevalence scenario transports one operating point to another prevalence
# under the explicit assumption that the class-conditional score distributions
# stay fixed: TPR and FPR are kept, and for a scenario prevalence pi and total
# N the expected cells are
#
#   TP' = TPR pi N,  FN' = (1 - TPR) pi N,  FP' = FPR (1 - pi) N,  TN' = (1 - FPR) (1 - pi) N.
#
# Expected counts are generally fractional; they are not observed units and
# carry no sampling interval. Without events or without non-events TPR or FPR
# is undefined and so is the scenario (no counts are invented). This is not a
# claim about real population shifts, which can also change TPR and FPR
# (spectrum or case-mix effects). For an outcome-selected (case-control)
# evaluation sample, the scenario at the population prevalence gives the
# population PPV and NPV that the sample itself cannot estimate.

.ig_cls_scenario_assumption <- paste(
  "Expected counts under fixed class-conditional score distributions (TPR and FPR of the",
  "operating point kept); a scenario, not an observed measurement and not a claim about",
  "real population shift. Fractional expected counts are not observed units.")

# Transport cells (vectors) to prevalence pi and total n_total. Returns a data
# frame of expected cells (NA where TPR or FPR is undefined) and the status.
.ig_cls_transport <- function(tp, fn, fp, tn, pi, n_total) {
  P <- tp + fn
  N <- fp + tn
  tpr <- tp / P
  fpr <- fp / N
  Pn <- pi * n_total
  Nn <- (1 - pi) * n_total
  status <- ifelse(P == 0, "undefined", ifelse(N == 0, "undefined", "ok"))
  reason <- ifelse(P == 0, "no_events", ifelse(N == 0, "no_non_events", NA_character_))
  out <- data.frame(tp = tpr * Pn, fn = (1 - tpr) * Pn, fp = fpr * Nn, tn = (1 - fpr) * Nn)
  out[status != "ok", ] <- NA_real_
  list(counts = out, status = status, reason = reason)
}

# Scenario record of one operating point (observed or weighted counts).
.ig_cls_scenario_record <- function(counts, pi, n_total, source, metric_note = NULL) {
  tr <- .ig_cls_transport(counts[["tp"]], counts[["fn"]], counts[["fp"]], counts[["tn"]], pi, n_total)
  rec <- list(prevalence = as.double(pi), n_total = as.double(n_total), source = source,
              assumption = .ig_cls_scenario_assumption, count_kind = "expected")
  if (tr$status != "ok") {
    rec$status <- "undefined"
    rec$reason_code <- tr$reason
    rec$counts <- NULL
    rec$metrics <- NULL
    return(rec)
  }
  ec <- stats::setNames(as.double(tr$counts[1, c("tp", "fn", "fp", "tn")]), c("tp", "fn", "fp", "tn"))
  rec$status <- "ok"
  rec$reason_code <- NA_character_
  rec$counts <- as.list(ec)
  rec$metrics <- .ig_cls_result_metric_cols(ec, "expected_scenario", context = list())
  rec
}

#' Prevalence scenarios for a fixed operating point
#'
#' Transports operating points to assumed prevalences under fixed
#' class-conditional score distributions (TPR and FPR kept) and recomputes all
#' operating-point metrics from the expected counts. The values are scenarios,
#' not new observed measurements: they assume that only the prevalence
#' changes, which real population shifts need not satisfy. Expected counts are
#' fractional, are not observed units and have no sampling interval. Without
#' events or without non-events at the operating point, the scenario is
#' undefined rather than computed with an invented rate.
#'
#' @param x Named counts `tp`, `fn`, `fp`, `tn` of one operating point, or a
#'   `classification_result` (its observed operating point per model).
#' @param prevalence Numeric vector of scenario prevalences in \[0, 1\].
#' @param n_total Total of the scenario population (default: the total of the
#'   counts); it scales the expected counts, not the rates.
#' @param source Provenance of the prevalence: `"assumed"` (default) or
#'   `"external"` (for example a registry estimate).
#' @param lang Language of reason texts and denominator labels.
#' @return A data frame with `model_id`, `scenario_prevalence`,
#'   `scenario_n_total`, `scenario_source`, `scenario_status` and the metric
#'   columns of [classification_metrics()] (estimate kind
#'   `"expected_scenario"`).
#' @examples
#' # Case-control style sample (50 % events) transported to 2 % prevalence
#' prevalence_scenario(c(tp = 80, fn = 20, fp = 10, tn = 90), prevalence = c(0.02, 0.5))
#' @concept classification
#' @seealso
#' * [classification_spec()]
#' * [evaluate_cutoff()]
#' * [plot_classification()]
#' @export
prevalence_scenario <- function(x, prevalence, n_total = NULL, source = c("assumed", "external"),
                                lang = c("en", "de")) {
  source <- match.arg(source)
  lang <- match.arg(lang)
  if (!is.numeric(prevalence) || !length(prevalence) || any(!is.finite(prevalence)) ||
      any(prevalence < 0 | prevalence > 1)) {
    .ig_abort("invalid_input", "`prevalence` must be numeric values in [0, 1]")
  }
  if (!is.null(n_total) && !(is.numeric(n_total) && length(n_total) == 1L && is.finite(n_total) && n_total > 0 &&
                             n_total <= .ig_cls_max_total)) {
    .ig_abort("invalid_input", "`n_total` must be a single number > 0 and at most 1e100")
  }
  points <- list()
  if (inherits(x, "classification_result")) {
    for (m in x$models) {
      op <- m$operating_point
      if (is.null(op) || !identical(.ig_cls_status_code(m$status), "ok")) next
      if (identical(op$count_kind, "expected")) {
        .ig_abort("invalid_input", "the operating point already holds expected counts of a theoretical model; change the prevalence in the specification instead")
      }
      points[[m$model_id]] <- stats::setNames(c(op$tp, op$fn, op$fp, op$tn), c("tp", "fn", "fp", "tn"))
    }
    if (!length(points)) .ig_abort("invalid_input", "the result holds no observed operating point")
  } else {
    points[["counts"]] <- .ig_cls_check_counts(x, "frequency_weighted")
  }
  blocks <- list()
  for (mid in names(points)) {
    cnt <- points[[mid]]
    nt <- n_total %||% sum(cnt)
    for (pi in prevalence) {
      rec <- .ig_cls_scenario_record(cnt, pi, nt, source)
      tab <- if (identical(rec$status, "ok")) .ig_cls_frame(rec$metrics, "metrics", lang) else {
        data.frame(metric_id = NA_character_, value = NA_real_, status = "undefined",
                   reason_code = rec$reason_code, reason = .ig_cls_reason_text(rec$reason_code, lang),
                   stringsAsFactors = FALSE)
      }
      blocks[[length(blocks) + 1L]] <- cbind(
        model_id = mid, scenario_prevalence = pi, scenario_n_total = nt, scenario_source = source,
        scenario_status = rec$status, tab, stringsAsFactors = FALSE)
    }
  }
  all_names <- unique(unlist(lapply(blocks, names)))
  out <- do.call(rbind, lapply(blocks, function(b) {
    for (nm in setdiff(all_names, names(b))) b[[nm]] <- NA
    b[all_names]
  }))
  rownames(out) <- NULL
  attr(out, "assumption") <- .ig_cls_scenario_assumption
  out
}
