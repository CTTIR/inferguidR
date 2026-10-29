# Sampling uncertainty of classification metrics ------------------------------------------
#
# * Wilson score and Clopper-Pearson intervals for proportions with an integer
#   numerator and denominator from independent units. They are conditional on
#   the fixed cutoff (and on the class totals for class-conditional rates).
# * Restrictions are explicit, never silent: expected (theoretical or
#   scenario) counts have no sampling interval; clustered units make binomial
#   and DeLong intervals unsupported (use the cluster bootstrap); under
#   outcome-selected (case-control) sampling the prevalence-dependent
#   proportions are sample-specific, so no population interval is reported.
#   Class-conditional rates (TPR, FNR, TNR, FPR) stay valid there.
# * Bootstrap: percentile intervals from a seeded resampling of units. Units
#   with frequency weight w count as w replicated units; a replicate draws
#   multinomial counts of the replicated units (stratified by class by
#   default) or, with clusters, multinomial counts of whole clusters. The rule
#   (cutoff and comparator) stays fixed in every replicate, so the interval
#   does not include the uncertainty of a data-driven cutoff selection.
# * The global random-number state is restored after the bootstrap.

.ig_cls_interval_statuses <- c("ok", "undefined", "positive_infinity", "not_applicable",
                               "not_computed", "unsupported")

# Wilson score interval (vectorised over x and n; n > 0).
.ig_cls_wilson <- function(x, n, level) {
  z <- stats::qnorm((1 - level) / 2, lower.tail = FALSE)
  p <- x / n
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  lower <- pmax(0, centre - half)
  upper <- pmin(1, centre + half)
  lower[x == 0] <- 0
  upper[x == n] <- 1
  list(lower = lower, upper = upper)
}

# Clopper-Pearson (exact binomial) interval (vectorised; n > 0).
.ig_cls_clopper_pearson <- function(x, n, level) {
  a <- 1 - level
  lower <- ifelse(x == 0, 0, stats::qbeta(a / 2, pmax(x, 1e-300), n - x + 1))
  upper <- ifelse(x == n, 1, stats::qbeta(a / 2, x + 1, pmax(n - x, 1e-300), lower.tail = FALSE))
  list(lower = lower, upper = upper)
}

.ig_cls_interval_row <- function(metric_id, method, level, lower = NA_real_, upper = NA_real_,
                                 status = "ok", reason_code = NA_character_,
                                 conditioning = "fixed_cutoff", note = NA_character_,
                                 variant = NA_character_) {
  if (status != "ok") {
    lower <- NA_real_
    upper <- NA_real_
  }
  list(metric_id = metric_id, variant = variant, method = method, level = as.double(level),
       lower = as.double(lower), upper = as.double(upper), status = status,
       reason_code = reason_code, conditioning = conditioning, note = note)
}

# Binomial intervals of every proportion metric of the dictionary for one
# set of counts. `context` holds count_kind, clustered, case_control and an
# optional note (for example the post-selection caveat).
.ig_cls_proportion_intervals <- function(counts, methods, level, context) {
  if (!length(methods)) return(list())
  dict <- .ig_cls_dictionary()$metrics
  parts <- .ig_cls_proportion_parts(counts[["tp"]], counts[["fn"]], counts[["fp"]], counts[["tn"]])
  vals <- .ig_cls_metrics_vec(counts[["tp"]], counts[["fn"]], counts[["fp"]], counts[["tn"]])
  rows <- list()
  for (id in names(parts)) {
    entry <- dict[[id]]
    for (meth in intersect(methods, unlist(entry$interval_methods))) {
      x <- parts[[id]][1]
      n <- parts[[id]][2]
      prev_dep <- !identical(entry$prevalence_dependence,
                             "invariant_given_fixed_class_conditional_distributions")
      row <- if (context$count_kind == "expected") {
        .ig_cls_interval_row(id, meth, level, status = "not_applicable",
                             reason_code = context$expected_reason %||% "expected_counts_no_interval")
      } else if (isTRUE(context$clustered)) {
        .ig_cls_interval_row(id, meth, level, status = "unsupported", reason_code = "clustered_units")
      } else if (isTRUE(context$case_control) && prev_dep) {
        .ig_cls_interval_row(id, meth, level, status = "not_applicable", reason_code = "case_control_sampling")
      } else if (vals[[id]]$status != "ok") {
        .ig_cls_interval_row(id, meth, level, status = vals[[id]]$status, reason_code = vals[[id]]$reason)
      } else {
        ci <- if (meth == "wilson") .ig_cls_wilson(x, n, level) else .ig_cls_clopper_pearson(x, n, level)
        .ig_cls_interval_row(id, meth, level, ci$lower, ci$upper, note = context$note %||% NA_character_)
      }
      rows[[length(rows) + 1L]] <- row
    }
  }
  rows
}

# Run `expr` with a fixed seed and restore the caller's random-number state
# (including the generator kinds) afterwards.
.ig_cls_with_seed <- function(seed, expr) {
  genv <- globalenv()
  had_seed <- exists(".Random.seed", envir = genv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = genv, inherits = FALSE) else NULL
  old_kind <- RNGkind()
  on.exit({
    RNGkind(old_kind[1], old_kind[2], old_kind[3])
    if (had_seed) assign(".Random.seed", old_seed, envir = genv) else
      if (exists(".Random.seed", envir = genv, inherits = FALSE)) rm(".Random.seed", envir = genv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  force(expr)
}

# Alert indicator of a fixed rule on original scores.
.ig_cls_alert_vec <- function(score, direction, comparator, cutoff, cutoff_kind) {
  if (identical(cutoff_kind, "alert_none")) return(rep(FALSE, length(score)))
  if (identical(cutoff_kind, "alert_all")) return(rep(TRUE, length(score)))
  o <- .ig_cls_orientation(direction) * score
  thr <- .ig_cls_orientation(direction) * cutoff
  if (.ig_cls_inclusive(comparator)) o >= thr else o > thr
}

# Draw the replicate weights of one bootstrap replicate: the number of times
# each unit (row) is drawn. Units with frequency weight w count as w
# replicated units. rmultinom() gets the unnormalised weights: it normalises
# them itself, and where long double has no extra precision (macOS arm64) a
# vector normalised beforehand can draw differently from the same seed.
.ig_cls_bootstrap_draw <- function(y, w, cluster_index, n_clusters, strata) {
  n <- length(y)
  if (!is.null(cluster_index)) {
    m <- stats::rmultinom(1L, n_clusters, rep(1, n_clusters))[, 1]
    return(m[cluster_index] * w)
  }
  wb <- numeric(n)
  groups <- if (identical(strata, "class")) list(which(y == 1L), which(y == 0L)) else list(seq_len(n))
  for (g in groups) {
    g <- g[w[g] > 0]
    if (!length(g)) next
    size <- sum(w[g])
    wb[g] <- stats::rmultinom(1L, size, w[g])[, 1]
  }
  wb
}

# Percentile bootstrap intervals for one scored model at a fixed rule.
.ig_cls_bootstrap <- function(score, y, w, cluster, direction, comparator, cutoff, cutoff_kind,
                              unc, partial, score_scale, note) {
  R <- as.integer(unc$bootstrap_replicates)
  level <- unc$level
  ids <- unc$bootstrap_metrics
  strata <- if (is.null(cluster)) unc$bootstrap_strata else "cluster"
  if (is.null(w)) w <- rep(1, length(y))
  cl_index <- NULL
  n_clusters <- 0L
  if (!is.null(cluster)) {
    f <- factor(cluster)
    cl_index <- as.integer(f)
    n_clusters <- nlevels(f)
  }
  alert <- .ig_cls_alert_vec(score, direction, comparator, cutoff, cutoff_kind)
  ev <- y == 1L
  op_ids <- intersect(ids, .ig_cls_operating_ids())
  curve_ids <- intersect(ids, .ig_cls_bootstrap_curve_ids)
  if (!is.null(.ig_cls_probability_unavailable(score_scale, direction))) {
    curve_ids <- setdiff(curve_ids, "brier_score")
  }
  if (is.null(partial)) curve_ids <- setdiff(curve_ids, "partial_auc")
  need_sets <- any(c("partial_auc", "average_precision", "pr_area_trapezoid") %in% curve_ids)
  cells <- matrix(NA_real_, R, 4L, dimnames = list(NULL, c("tp", "fn", "fp", "tn")))
  cv <- matrix(NA_real_, R, length(curve_ids), dimnames = list(NULL, curve_ids))
  .ig_cls_with_seed(unc$seed, {
    for (b in seq_len(R)) {
      wb <- .ig_cls_bootstrap_draw(y, w, cl_index, n_clusters, unc$bootstrap_strata)
      cells[b, ] <- c(sum(wb[alert & ev]), sum(wb[!alert & ev]), sum(wb[alert & !ev]),
                      sum(wb[!alert & !ev]))
      if (!length(curve_ids)) next
      if ("roc_auc" %in% curve_ids) cv[b, "roc_auc"] <- .ig_cls_auc_stats(score, y, wb, direction)$auc
      if (need_sets) {
        ds <- .ig_cls_decision_sets(score, y, wb, direction, comparator)
        if ("partial_auc" %in% curve_ids && ds$tp[1] + ds$fn[1] > 0 && ds$fp[1] + ds$tn[1] > 0) {
          pts <- .ig_cls_roc_points(ds)
          cv[b, "partial_auc"] <- .ig_cls_partial_values(pts$fpr, pts$tpr, partial)$raw
        }
        if (any(c("average_precision", "pr_area_trapezoid") %in% curve_ids)) {
          pr <- .ig_cls_pr_values(ds)
          if ("average_precision" %in% curve_ids) cv[b, "average_precision"] <- pr$ap
          if ("pr_area_trapezoid" %in% curve_ids) cv[b, "pr_area_trapezoid"] <- pr$trap
        }
      }
      if ("brier_score" %in% curve_ids) cv[b, "brier_score"] <- sum(wb * (score - y)^2) / sum(wb)
    }
  })
  boot <- list()
  if (length(op_ids)) {
    mv <- .ig_cls_metrics_vec(cells[, "tp"], cells[, "fn"], cells[, "fp"], cells[, "tn"])
    for (id in op_ids) boot[[id]] <- list(value = mv[[id]]$value, status = mv[[id]]$status)
  }
  for (id in curve_ids) {
    v <- cv[, id]
    boot[[id]] <- list(value = v, status = ifelse(is.na(v), "undefined", "ok"))
  }
  design <- switch(strata, class = "stratified by class", none = "unstratified",
                   cluster = "whole clusters resampled")
  rows <- list()
  for (id in names(boot)) {
    v <- boot[[id]]$value
    st <- boot[[id]]$status
    v[st == "positive_infinity"] <- Inf
    unusable <- !(st %in% c("ok", "positive_infinity"))
    n_bad <- sum(unusable)
    txt <- sprintf("%d replicates (%s), seed %s, percentile; %d replicate(s) with an undefined value excluded",
                   R, design, format(unc$seed), n_bad)
    if (!is.na(note)) txt <- paste0(txt, "; ", note)
    variant <- if (id == "partial_auc") "raw" else NA_character_
    conditioning <- if (id %in% .ig_cls_bootstrap_curve_ids) "fixed_scores" else "fixed_cutoff"
    if (n_bad > 0.1 * R) {
      rows[[length(rows) + 1L]] <- .ig_cls_interval_row(id, "bootstrap_percentile", level, status = "undefined",
                                                        reason_code = "bootstrap_too_many_undefined",
                                                        conditioning = conditioning, note = txt,
                                                        variant = variant)
      next
    }
    q <- stats::quantile(v[!unusable], c((1 - level) / 2, 1 - (1 - level) / 2), type = 7,
                         names = FALSE)
    if (any(!is.finite(q))) {
      rows[[length(rows) + 1L]] <- .ig_cls_interval_row(id, "bootstrap_percentile", level,
                                                        status = "positive_infinity",
                                                        reason_code = "bootstrap_limit_infinite",
                                                        conditioning = conditioning, note = txt,
                                                        variant = variant)
    } else {
      rows[[length(rows) + 1L]] <- .ig_cls_interval_row(id, "bootstrap_percentile", level, q[1], q[2],
                                                        conditioning = conditioning, note = txt,
                                                        variant = variant)
    }
  }
  rows
}
