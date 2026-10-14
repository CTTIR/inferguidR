# Exact decision-set enumeration -----------------------------------------------------
#
# A cutoff rule on original scores partitions the units into alerts and
# non-alerts. Units with equal scores always fall on the same side, so the
# distinct decision sets of a score column with m distinct values are the m + 1
# sets "alert the k highest-risk distinct values", k = 0 (alert none) ... m
# (alert all). The comparator decides which cutoff represents a set and what
# happens to scores equal to the cutoff:
#
#   >=  alert if score >= c     >  alert if score > c
#   <=  alert if score <= c     <  alert if score < c
#
# Internally scores are oriented (o = s for higher_is_positive, o = -s for
# lower_is_positive) so that alerts are always the high oriented values;
# reported cutoffs are on the original scale and are always original score
# values. A set that no original score can represent under the comparator
# (alert none under >=, alert all under >) has cutoff NA and cutoff_kind
# "alert_none" / "alert_all". Histogram bins are never used.

.ig_cls_orientation <- function(direction) if (identical(direction, "lower_is_positive")) -1 else 1
.ig_cls_inclusive <- function(comparator) comparator %in% c(">=", "<=")

# Per distinct oriented score, descending: event weight, non-event weight.
# Units with frequency weight 0 are absent and are dropped.
.ig_cls_score_groups <- function(score, y, w, direction) {
  if (is.null(w)) w <- rep(1, length(y))
  keep <- w > 0
  o <- .ig_cls_orientation(direction) * score[keep]
  y <- y[keep]
  w <- w[keep]
  u <- sort(unique(o), decreasing = TRUE)
  idx <- match(o, u)
  pos <- as.double(.ig_cls_tabulate_weighted(idx[y == 1L], w[y == 1L], length(u)))
  neg <- as.double(.ig_cls_tabulate_weighted(idx[y == 0L], w[y == 0L], length(u)))
  list(u = u, pos = pos, neg = neg, n_pos = sum(pos), n_neg = sum(neg))
}

# Weighted tabulation of integer bins 1..nbins.
.ig_cls_tabulate_weighted <- function(bin, weight, nbins) {
  out <- numeric(nbins)
  if (length(bin)) {
    s <- rowsum(weight, bin, reorder = TRUE)
    out[as.integer(rownames(s))] <- s[, 1]
  }
  out
}

# All distinct decision sets for one model (data frame ordered by alerts).
.ig_cls_decision_sets <- function(score, y, w, direction, comparator) {
  g <- .ig_cls_score_groups(score, y, w, direction)
  m <- length(g$u)
  sgn <- .ig_cls_orientation(direction)
  tp <- c(0, cumsum(g$pos))
  fp <- c(0, cumsum(g$neg))
  k <- 0:m
  if (.ig_cls_inclusive(comparator)) {
    # set k (k >= 1) = {o >= u_k}; set 0 needs a cutoff above every score
    cutoff <- c(NA_real_, sgn * g$u)
    kind <- c("alert_none", rep("score", m))
  } else {
    # set k (k < m) = {o > u_(k+1)}; set m needs a cutoff below every score
    cutoff <- c(sgn * g$u, NA_real_)
    kind <- c(rep("score", m), "alert_all")
  }
  # lowest-risk score inside each alert set (original scale), for display
  boundary <- c(NA_real_, sgn * g$u)
  data.frame(set_index = as.double(k), cutoff = cutoff, cutoff_kind = kind,
             boundary_score = boundary, alerts = tp + fp,
             tp = tp, fn = g$n_pos - tp, fp = fp, tn = g$n_neg - fp,
             stringsAsFactors = FALSE)
}

# Counts at one cutoff on the original scale, computed directly from the
# scores, plus the index of the decision set it selects.
.ig_cls_counts_at <- function(score, y, w, direction, comparator, cutoff) {
  if (is.null(w)) w <- rep(1, length(y))
  sgn <- .ig_cls_orientation(direction)
  o <- sgn * score
  thr <- sgn * cutoff
  alert <- if (.ig_cls_inclusive(comparator)) o >= thr else o > thr
  present <- w > 0
  u <- unique(o[present])
  set_index <- if (.ig_cls_inclusive(comparator)) sum(u >= thr) else sum(u > thr)
  c(tp = sum(w[alert & y == 1L]), fn = sum(w[!alert & y == 1L]),
    fp = sum(w[alert & y == 0L]), tn = sum(w[!alert & y == 0L]),
    set_index = set_index)
}

.ig_cls_default_sweep_metrics <- c("tpr", "fpr", "ppv", "npv", "f1", "mcc", "youden_j", "lr_pos")

# Add metric value/status columns to a table with tp, fn, fp, tn columns.
.ig_cls_add_metric_columns <- function(tab, metrics) {
  if (!length(metrics)) return(tab)
  vals <- .ig_cls_metrics_vec(tab$tp, tab$fn, tab$fp, tab$tn)
  for (id in metrics) {
    tab[[id]] <- vals[[id]]$value
    tab[[paste0(id, "_status")]] <- vals[[id]]$status
  }
  tab
}

#' Sweep all distinct decision sets of a scored model
#'
#' Enumerates every distinct empirical decision set (alert none, alert all and
#' every set in between) from the original scores, with exact counts and the
#' requested metrics. The cutoff of each set is an original score value under
#' the declared comparator; a set that no score can represent (alert none
#' under `>=`, alert all under `>`) has cutoff `NA` and a `cutoff_kind` that
#' says so. For a theoretical model the sweep is a display grid of real-valued
#' cutoffs with expected counts. Aggregate 2x2 counts cannot be swept.
#'
#' @param x A `classification_spec` (scored data with attached data, or
#'   theory) or a `classification_result` that holds decision sets.
#' @param metrics Metric ids (see [metric_dictionary()]) to add as columns;
#'   each value column has a companion `<id>_status` column.
#' @param model_id Optional model id (default: all models).
#' @param grid For theoretical models: number of display cutoffs (default
#'   401) or a numeric vector of cutoffs.
#' @param basis For results: `"observed"` (default) or `"scenario"` (expected
#'   counts of every decision set at the result's scenario prevalence, under
#'   fixed class-conditional score distributions).
#' @return A data frame with `model_id`, `set_index`, `cutoff`, `cutoff_kind`,
#'   `alerts`, `tp`, `fn`, `fp`, `tn`, `count_kind` and the metric columns;
#'   for results also `selected` (the reported operating point, always
#'   included).
#' @examples
#' d <- data.frame(y = c(1, 1, 1, 0, 0, 0), s = c(5, 4, 4, 4, 2, 1))
#' spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "s", cutoff = 4)
#' cutoff_sweep(spec, metrics = c("tpr", "fpr", "ppv"))
#' @family classification
#' @export
cutoff_sweep <- function(x, metrics = .ig_cls_default_sweep_metrics, model_id = NULL, grid = 401,
                         basis = c("observed", "scenario")) {
  basis <- match.arg(basis)
  bad <- setdiff(metrics, .ig_cls_operating_ids())
  if (length(bad)) .ig_abort("invalid_input", sprintf("unknown operating-point metric id(s): %s",
                                                      paste(bad, collapse = ", ")))
  if (inherits(x, "classification_result")) {
    if (identical(x$input_kind, "counts_2x2")) .ig_abort("unsupported", .ig_cls_reason_text("aggregate_counts_only"))
    if (basis == "scenario" && (!identical(x$prevalence$kind, "scenario") ||
                                identical(x$input_kind, "theoretical_binormal"))) {
      .ig_abort("invalid_input", "basis 'scenario' needs a scored-data result with a scenario prevalence")
    }
    blocks <- list()
    for (m in x$models) {
      if (!is.null(model_id) && !m$model_id %in% model_id) next
      if (is.null(m$curve)) next
      tab <- .ig_cls_frame(m$curve, "curve")
      count_kind <- if (identical(x$input_kind, "theoretical_binormal")) "expected" else m$count_kind
      if (basis == "scenario") {
        n_all <- tab$tp[1] + tab$fn[1] + tab$fp[1] + tab$tn[1]
        tr <- .ig_cls_transport(tab$tp, tab$fn, tab$fp, tab$tn, x$prevalence$value, n_all)$counts
        tab[, c("tp", "fn", "fp", "tn")] <- tr
        tab$alerts <- tab$tp + tab$fp
        count_kind <- "expected"
      }
      tab <- cbind(model_id = m$model_id, tab, count_kind = count_kind, stringsAsFactors = FALSE)
      op <- m$operating_point
      tab$selected <- if (is.null(op)) FALSE else if (identical(x$input_kind, "theoretical_binormal")) {
        tab$cutoff_kind == "operating_point" | (op$cutoff_kind %in% c("alert_none", "alert_all") &
                                                  tab$cutoff_kind == op$cutoff_kind)
      } else !is.na(tab$set_index) & tab$set_index == op$set_index
      blocks[[length(blocks) + 1L]] <- .ig_cls_add_metric_columns(tab, metrics)
    }
    if (!length(blocks)) {
      .ig_abort("unsupported", "the result holds no decision sets (computed with include_curves = FALSE, or no model matches `model_id`)")
    }
    out <- do.call(rbind, blocks)
    rownames(out) <- NULL
    return(out)
  }
  if (basis != "observed") .ig_abort("invalid_input", "basis 'scenario' is available for results only")
  .ig_cls_check_spec(x)
  if (identical(x$input_kind, "counts_2x2")) {
    .ig_abort("unsupported", .ig_cls_reason_text("aggregate_counts_only"))
  }
  if (identical(x$input_kind, "theoretical_binormal")) {
    tab <- .ig_cls_theory_grid(x, grid)
    tab <- cbind(model_id = "binormal", tab, stringsAsFactors = FALSE)
    return(.ig_cls_add_metric_columns(tab, metrics))
  }
  prep <- .ig_cls_spec_data(x)
  blocks <- list()
  for (i in seq_along(x$models)) {
    m <- x$models[[i]]
    if (!is.null(model_id) && !m$model_id %in% model_id) next
    tab <- .ig_cls_decision_sets(prep$scores[[i]], prep$y, prep$w, m$direction, m$comparator)
    tab <- cbind(model_id = m$model_id, tab, comparator = m$comparator,
                 count_kind = if (is.null(prep$w)) "observed" else "frequency_weighted",
                 stringsAsFactors = FALSE)
    blocks[[length(blocks) + 1L]] <- .ig_cls_add_metric_columns(tab, metrics)
  }
  if (!length(blocks)) .ig_abort("invalid_input", "no model matches `model_id`")
  out <- do.call(rbind, blocks)
  rownames(out) <- NULL
  out
}
