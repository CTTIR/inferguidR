# Cutoff selection ----------------------------------------------------------------------
#
# A selection objective (max Youden J, max F1, max MCC, an operating-point
# metric under constraints, capacity, declared costs) is optimised over the
# admissible operating points:
#
# * scored data: every distinct empirical decision set (exhaustive search, no
#   grid, no rounding);
# * theoretical binormal model: real-valued cutoffs plus the two limiting
#   rules "alert none" and "alert all". A dense scan brackets the optimum;
#   constraint boundaries are refined by bisection and interior optima by
#   optimize(). For Youden's J without constraints the analytic stationary
#   points are added as candidates.
#
# Constraints and capacity limits are bounds on metrics; an undefined metric
# never satisfies a bound. A metric meets a bound when it lies on the right
# side up to floating-point rounding (.ig_cls_bound_tol relative to the
# larger of the two magnitudes), so that tiny bounds such as FPR <= 1e-9 are
# honoured exactly. Values within .ig_cls_tie_tol (relative) of the best value
# are ties; every tied operating point is reported and the declared tie break
# chooses among actual operating points, never among rounded displayed values.

.ig_cls_tie_tol <- 1e-12                    # relative tolerance for equal objective values
.ig_cls_bound_tol <- 4 * .Machine$double.eps # relative rounding tolerance against a bound

.ig_cls_selection_codes <- c("max_youden", "max_f1", "max_mcc", "constraint", "capacity",
                             "cost_utility")

.ig_cls_bound_defs <- list(
  tpr_min = c(metric = "tpr", sense = "at_least"),
  tnr_min = c(metric = "tnr", sense = "at_least"),
  ppv_min = c(metric = "ppv", sense = "at_least"),
  npv_min = c(metric = "npv", sense = "at_least"),
  fpr_max = c(metric = "fpr", sense = "at_most"),
  fnr_max = c(metric = "fnr", sense = "at_most"),
  fdr_max = c(metric = "fdr", sense = "at_most"),
  alert_rate_max = c(metric = "alert_rate", sense = "at_most"),
  max_alerts = c(metric = "n_alerts", sense = "at_most"),
  max_alert_rate = c(metric = "alert_rate", sense = "at_most")
)

# Active bounds of a rule: declared constraints plus a capacity limit.
.ig_cls_bounds <- function(rule) {
  out <- list()
  add <- function(name, source, bound) {
    def <- .ig_cls_bound_defs[[name]]
    out[[length(out) + 1L]] <<- list(name = name, source = source, metric = def[["metric"]],
                                     sense = def[["sense"]], bound = as.double(bound))
  }
  for (nm in names(rule$constraints)) add(nm, "constraint", rule$constraints[[nm]])
  if (identical(rule$selection, "capacity")) {
    for (nm in names(rule$capacity)) add(nm, "capacity", rule$capacity[[nm]])
  }
  out
}

# Feasibility of cell counts (vectors) under a list of bounds. Returns
# feasible (logical), slack (minimum signed distance to the bounds; >= 0 when
# feasible; -1 for an undefined metric), the per-bound satisfaction matrix
# and the per-bound metric values and statuses.
.ig_cls_constraint_fun <- function(bounds) {
  function(tp, fn, fp, tn) {
    k <- length(tp)
    if (!length(bounds)) {
      return(list(feasible = rep(TRUE, k), slack = rep(1, k),
                  satisfied = matrix(TRUE, k, 0L), values = matrix(NA_real_, k, 0L),
                  statuses = matrix(NA_character_, k, 0L)))
    }
    vals <- .ig_cls_metrics_vec(tp, fn, fp, tn)
    slack <- matrix(NA_real_, k, length(bounds))
    sat <- matrix(FALSE, k, length(bounds))
    values <- matrix(NA_real_, k, length(bounds))
    statuses <- matrix(NA_character_, k, length(bounds))
    for (j in seq_along(bounds)) {
      b <- bounds[[j]]
      v <- vals[[b$metric]]
      s <- if (b$sense == "at_least") v$value - b$bound else b$bound - v$value
      tol <- .ig_cls_bound_tol * pmax(abs(v$value), abs(b$bound))
      ok <- v$status == "ok"
      sat[, j] <- ok & s >= -tol
      s[!ok] <- -1
      slack[, j] <- s
      values[, j] <- v$value
      statuses[, j] <- v$status
    }
    list(feasible = rowSums(!sat) == 0L, slack = apply(slack, 1L, min), satisfied = sat,
         values = values, statuses = statuses)
  }
}

# Per-bound summary over the evaluated operating points: how many satisfy the
# bound, the best value attainable and whether the chosen point satisfies it.
.ig_cls_bound_summary <- function(bounds, cf, chosen = NA_integer_) {
  if (!length(bounds)) return(NULL)
  rows <- lapply(seq_along(bounds), function(j) {
    b <- bounds[[j]]
    ok <- cf$statuses[, j] == "ok"
    best <- if (any(ok)) {
      if (b$sense == "at_least") max(cf$values[ok, j]) else min(cf$values[ok, j])
    } else NA_real_
    list(constraint = b$name, source = b$source, metric_id = b$metric, sense = b$sense,
         bound = b$bound, n_satisfied = as.double(sum(cf$satisfied[, j])),
         n_points = as.double(nrow(cf$satisfied)),
         best_value = best, best_status = if (any(ok)) "ok" else "undefined",
         satisfied_at_choice = if (is.na(chosen)) NA else cf$satisfied[chosen, j])
  })
  .ig_cls_rows_to_cols(rows)
}

# Objective of a selection rule: metric id, direction and a vectorised value
# function of the four cells (NA where undefined, Inf where +infinity).
.ig_cls_objective <- function(rule) {
  sel <- rule$selection
  if (identical(sel, "manual")) return(NULL)
  if (identical(sel, "cost_utility")) {
    cc <- c(tp = 0, fn = 0, fp = 0, tn = 0)
    cc[names(rule$costs)] <- unlist(rule$costs)
    return(list(metric = "expected_cost_per_unit", direction = "min", costs = cc,
                values = function(tp, fn, fp, tn) {
                  n <- tp + fn + fp + tn
                  v <- (cc[["tp"]] * tp + cc[["fn"]] * fn + cc[["fp"]] * fp + cc[["tn"]] * tn) / n
                  v[n == 0] <- NA_real_
                  v
                }))
  }
  metric <- switch(sel, max_youden = "youden_j", max_f1 = "f1", max_mcc = "mcc",
                   capacity = "n_alerts", constraint = rule$optimise)
  entry <- .ig_cls_metric_entry(metric)
  direction <- if (identical(sel, "capacity") || !identical(entry$preferred_direction, "lower")) "max" else "min"
  list(metric = metric, direction = direction,
       values = function(tp, fn, fp, tn) {
         v <- .ig_cls_metrics_vec(tp, fn, fp, tn)[[metric]]
         out <- v$value
         out[v$status == "positive_infinity"] <- Inf
         out
       })
}

# Purely relative: an objective value of 1e-15 does not tie with 0 (that would
# let a tie break discard the only operating point that alerts anyone).
.ig_cls_is_tied <- function(v, best) {
  if (is.infinite(best)) return(!is.na(v) & v == best)
  !is.na(v) & (v == best | abs(v - best) <= .ig_cls_tie_tol * abs(best))
}

.ig_cls_selection_reason <- function(code) .ig_cls_reason_text(code)

# Common head of a selection record.
.ig_cls_selection_head <- function(rule, obj, method) {
  list(procedure = rule$selection,
       objective_metric = if (is.null(obj)) NA_character_ else obj$metric,
       objective_direction = if (is.null(obj)) NA_character_ else obj$direction,
       prevalence_basis = rule$prevalence_basis %||% "observed",
       tie_break = rule$tie_break, selection_data = rule$selection_data,
       constraints = rule$constraints, capacity = rule$capacity, costs = rule$costs,
       method = method)
}

# Manual cutoff: no search; declared constraints are checked at that point.
.ig_cls_manual_selection <- function(rule, counts_row) {
  bounds <- .ig_cls_bounds(rule)
  rec <- .ig_cls_selection_head(rule, NULL, "manual cutoff (no search)")
  rec$status <- "ok"
  rec$reason_code <- NA_character_
  if (length(bounds)) {
    cf <- .ig_cls_constraint_fun(bounds)(counts_row$tp, counts_row$fn, counts_row$fp, counts_row$tn)
    rec$bounds <- .ig_cls_bound_summary(bounds, cf, 1L)
    rec$constraints_met <- isTRUE(cf$feasible[1])
  }
  rec
}

# Exhaustive selection over the decision sets of one scored model. `sets`
# holds set_index, cutoff, cutoff_kind and alerts (observed); `counts` holds
# the cells on the selection basis (observed or scenario expected counts).
.ig_cls_select_sets <- function(sets, counts, rule) {
  obj <- .ig_cls_objective(rule)
  bounds <- .ig_cls_bounds(rule)
  cf <- .ig_cls_constraint_fun(bounds)(counts$tp, counts$fn, counts$fp, counts$tn)
  ov <- obj$values(counts$tp, counts$fn, counts$fp, counts$tn)
  rec <- .ig_cls_selection_head(rule, obj, "exhaustive search over all distinct decision sets")
  defined <- !is.na(ov)
  cand <- defined & cf$feasible
  rec$n_candidates <- as.double(nrow(sets))
  rec$n_defined <- as.double(sum(defined))
  rec$n_feasible <- as.double(sum(cf$feasible))
  rec$n_admissible <- as.double(sum(cand))
  if (!any(cand)) {
    code <- if (!any(defined)) "objective_undefined_everywhere" else
      if (!any(cf$feasible)) "constraints_infeasible" else "no_admissible_point"
    rec$status <- "infeasible"
    rec$reason_code <- code
    rec$reason <- .ig_cls_selection_reason(code)
    rec$bounds <- .ig_cls_bound_summary(bounds, cf)
    rec$row <- NA_integer_
    return(rec)
  }
  sv <- if (obj$direction == "max") ov else -ov
  best <- max(sv[cand])
  idx <- which(cand & .ig_cls_is_tied(sv, best))
  # Decision sets are ordered by alerts (set_index), so the fewest alerts is
  # the first tied set.
  pick <- if (identical(rule$tie_break, "most_alerts")) idx[length(idx)] else idx[1L]
  rec$status <- "ok"
  rec$reason_code <- NA_character_
  rec$row <- pick
  rec$objective_value <- ov[pick]
  rec$n_tied <- as.double(length(idx))
  rec$tied <- list(set_index = sets$set_index[idx], cutoff = sets$cutoff[idx],
                   cutoff_kind = sets$cutoff_kind[idx], alerts = sets$alerts[idx],
                   objective_value = ov[idx])
  rec$bounds <- .ig_cls_bound_summary(bounds, cf, pick)
  rec
}

# Cutoffs at which a bound of the theoretical model is met with equality:
# class-conditional rates from the normal quantiles (lower-tail quantiles
# where they are more accurate), the alert rate (and the expected number of
# alerts) by root finding on the log scale, and PPV, NPV and FDR (not
# monotone in general) at every sign change of their logit gap on a wide grid.
# The scan range is extended to cover all of them.
.ig_cls_theory_bound_cutoffs <- function(th, pi, bounds) {
  d <- .ig_cls_theory_d(th)
  s <- th$sigma
  log_alert_rate <- function(z) {
    a <- log(pi) + stats::pnorm((z - d) / s, lower.tail = FALSE, log.p = TRUE)
    b <- log1p(-pi) + stats::pnorm(z, lower.tail = FALSE, log.p = TRUE)
    m <- pmax(a, b)
    m + log(exp(a - m) + exp(b - m))
  }
  rate_cut <- function(v) {
    if (!is.finite(v) || v <= 0 || v >= 1) return(NA_real_)
    lo <- min(-40, d - 40 * s)
    hi <- max(40, d + 40 * s)
    stats::uniroot(function(z) log_alert_rate(z) - log(v), c(lo, hi), tol = 1e-13)$root
  }
  cuts <- vapply(bounds, function(b) {
    v <- b$bound
    switch(b$metric,
           tpr = d + s * stats::qnorm(v, lower.tail = FALSE),
           fnr = d + s * stats::qnorm(v),
           fpr = stats::qnorm(v, lower.tail = FALSE),
           tnr = stats::qnorm(v),
           alert_rate = rate_cut(v),
           n_alerts = rate_cut(v / th$n_total),
           NA_real_)
  }, double(1))
  # Bounds on PPV, NPV and FDR need not be monotone in the cutoff. Where they
  # hold with equality is found on a wide grid of the logit gap (both tails in
  # log space) and refined by uniroot. The grid spans every cutoff at which
  # TPR or FPR is representable (|z| <= 40 on the scale of either class), so
  # a feasible region far in a tail (PPV >= 0.5 at prevalence 1e-4 with
  # AUC 0.6 needs z near 26) is reached by the scan.
  logit_gap <- function(z, metric, target) {
    if (metric == "ppv") {
      log(pi) + stats::pnorm((z - d) / s, lower.tail = FALSE, log.p = TRUE) -
        log1p(-pi) - stats::pnorm(z, lower.tail = FALSE, log.p = TRUE) - stats::qlogis(target)
    } else {
      log1p(-pi) + stats::pnorm(z, log.p = TRUE) - log(pi) - stats::pnorm((z - d) / s, log.p = TRUE) -
        stats::qlogis(target)
    }
  }
  wide <- seq(min(-40, d - 40 * s), max(40, d + 40 * s), length.out = 8001L)
  extra <- numeric()
  for (b in bounds) {
    if (!b$metric %in% c("ppv", "npv", "fdr") || pi <= 0 || pi >= 1) next
    metric <- if (b$metric == "fdr") "ppv" else b$metric
    target <- if (b$metric == "fdr") 1 - b$bound else b$bound
    if (target <= 0 || target >= 1) next
    g <- logit_gap(wide, metric, target)
    idx <- which(is.finite(g[-1L]) & is.finite(g[-length(g)]) & g[-1L] * g[-length(g)] <= 0)
    for (i in idx) {
      extra <- c(extra, stats::uniroot(function(z) logit_gap(z, metric, target), wide[c(i, i + 1L)],
                                       tol = 1e-13)$root)
    }
  }
  cuts <- c(cuts, extra)
  cuts[is.finite(cuts)]
}

# Theoretical model: optimise over real cutoffs and the two limiting rules.
# Returns a selection record with `cutoff` (Inf / -Inf for alert none / all).
.ig_cls_theory_select <- function(th, pi, rule, grid_n = 4001L) {
  obj <- .ig_cls_objective(rule)
  bounds <- .ig_cls_bounds(rule)
  cfun <- .ig_cls_constraint_fun(bounds)
  counts_at <- function(z) .ig_cls_theory_counts(z, th, pi)
  fval <- function(z) {
    ec <- counts_at(z)
    v <- obj$values(ec$tp, ec$fn, ec$fp, ec$tn)
    if (obj$direction == "min") -v else v
  }
  gval <- function(z) {
    ec <- counts_at(z)
    cfun(ec$tp, ec$fn, ec$fp, ec$tn)
  }
  feasible_at <- function(z) gval(z)$feasible & !is.na(fval(z))
  r <- .ig_cls_theory_range(th)
  r <- c(r[1] - 4, r[2] + 4)
  # The scan must reach every boundary of a monotone bound, however extreme
  # (for example FPR <= 1e-20 lies near z = 9.3).
  bound_cuts <- .ig_cls_theory_bound_cutoffs(th, pi, bounds)
  if (length(bound_cuts)) r <- c(min(r[1], bound_cuts - 1), max(r[2], bound_cuts + 1))
  cuts <- seq(r[1], r[2], length.out = grid_n)
  step <- cuts[2] - cuts[1]
  fv <- fval(cuts)
  gv <- gval(cuts)
  feas <- gv$feasible & !is.na(fv)
  rec <- .ig_cls_selection_head(rule, obj, sprintf(
    "scan of %d cutoffs on [%s, %s] plus the limiting rules alert all / alert none; constraint boundaries refined by bisection, interior optima by optimize()",
    grid_n, format(signif(r[1], 6)), format(signif(r[2], 6))))
  rec$scan_range <- r
  rec$grid_n <- as.double(grid_n)

  # Bisection between an infeasible point a and a feasible point b; the
  # returned point is feasible.
  boundary <- function(a, b) {
    for (it in seq_len(200L)) {
      mid <- (a + b) / 2
      if (mid == a || mid == b) break
      if (feasible_at(mid)) b <- mid else a <- mid
    }
    b
  }
  cand <- numeric()
  runs <- rle(feas)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  # optimize() needs finite values (it warns and substitutes otherwise): an
  # undefined objective is the worst value, +/-infinity the largest magnitude.
  safe_f <- function(z) {
    v <- fval(z)
    big <- .Machine$double.xmax
    if (is.na(v)) -big else max(-big, min(big, v))
  }
  for (j in which(runs$values)) {
    i0 <- starts[j]
    i1 <- ends[j]
    lo <- if (i0 > 1L) boundary(cuts[i0 - 1L], cuts[i0]) else cuts[i0]
    hi <- if (i1 < grid_n) boundary(cuts[i1 + 1L], cuts[i1]) else cuts[i1]
    cand <- c(cand, lo, hi)
    if (i1 - i0 >= 2L) {
      seg <- i0:i1
      v <- fv[seg]
      inner <- seq(2L, length(seg) - 1L)
      loc <- inner[v[inner] >= v[inner - 1L] & v[inner] >= v[inner + 1L]]
      # a plateau gives a run of local maxima; refine the first of each run
      if (length(loc) > 1L) loc <- loc[c(TRUE, diff(loc) > 1L)]
      for (q in loc) {
        a <- max(lo, cuts[seg[q - 1L]])
        b <- min(hi, cuts[seg[q + 1L]])
        if (b > a) {
          o <- stats::optimize(safe_f, c(a, b), maximum = TRUE, tol = 1e-11)
          cand <- c(cand, o$maximum)
        }
      }
    }
  }
  if (identical(obj$metric, "youden_j") && !length(bounds)) {
    an <- .ig_cls_theory_youden(th)
    cand <- c(cand, an$roots)
    rec$analytic_stationary_points <- an$roots
  }
  cand <- c(cand, bound_cuts, -Inf, Inf)
  cand <- unique(cand)
  ok <- vapply(cand, feasible_at, logical(1))
  # Count feasibility of the limiting rules and the scan for the diagnostics.
  rec$n_candidates <- as.double(length(cand))
  rec$n_feasible_scan <- as.double(sum(feas))
  if (!any(ok)) {
    code <- if (all(is.na(fv))) "objective_undefined_everywhere" else
      if (!any(gv$feasible)) "constraints_infeasible" else "no_admissible_point"
    rec$status <- "infeasible"
    rec$reason_code <- code
    rec$reason <- .ig_cls_selection_reason(code)
    rec$bounds <- .ig_cls_bound_summary(bounds, gv)
    rec$cutoff <- NA_real_
    return(rec)
  }
  cand <- cand[ok]
  cv <- vapply(cand, fval, double(1))
  best <- max(cv)
  tied <- cand[.ig_cls_is_tied(cv, best)]
  # the same optimum reached from neighbouring scan points counts once; an
  # analytic stationary point represents its cluster, since on which side of
  # it optimize() stops is rounding noise that differs between platforms
  tied <- sort(tied)
  if (length(tied) > 1L) {
    cluster <- cumsum(c(TRUE, diff(tied) > 1e-7 | is.infinite(diff(tied))))
    tied <- vapply(split(tied, cluster), function(z) {
      exact <- z[z %in% rec$analytic_stationary_points]
      if (length(exact)) exact[1L] else z[1L]
    }, double(1), USE.NAMES = FALSE)
  }
  # fewer alerts means a higher cutoff
  chosen <- if (identical(rule$tie_break, "most_alerts")) min(tied) else max(tied)
  # An optimum within two grid steps of the scan edge needs a look beyond the
  # edge (the grid step can be large when sigma is large). Probes farther out
  # on that side decide:
  # * something there is better: the objective keeps improving, so the
  #   supremum is not attained;
  # * the limiting rule of that side is tied with the optimum: a plateau that
  #   reaches alert all / alert none, so the limiting rule attains the optimum
  #   and is reported;
  # * the probes are tied and the tie break prefers that side: the preferred
  #   operating point lies beyond every finite cutoff and is not attained
  #   (for example PPV rounding to 1 as the cutoff grows);
  # * otherwise the point is a genuine optimum close to the edge.
  if (is.finite(chosen) && (chosen - r[1] < 2 * step || r[2] - chosen < 2 * step)) {
    left <- chosen - r[1] < 2 * step
    width <- r[2] - r[1]
    probes <- if (left) r[1] - width * c(0.25, 1, 4, 16) else r[2] + width * c(0.25, 1, 4, 16)
    pv <- vapply(probes, function(z) if (feasible_at(z)) fval(z) else NA_real_, double(1))
    better <- !is.na(pv) & pv > best & !.ig_cls_is_tied(pv, best)
    level <- .ig_cls_is_tied(pv, best)
    limit <- if (left) -Inf else Inf
    toward <- left == identical(rule$tie_break, "most_alerts")
    if (!any(better) && limit %in% tied) {
      chosen <- limit
    } else if (any(better) || (any(level) && toward)) {
      rec$status <- "infeasible"
      rec$reason_code <- "supremum_not_attained"
      rec$reason <- .ig_cls_selection_reason("supremum_not_attained")
      rec$bounds <- .ig_cls_bound_summary(bounds, gv)
      rec$cutoff <- NA_real_
      return(rec)
    }
  }
  ec <- counts_at(chosen)
  rec$status <- "ok"
  rec$reason_code <- NA_character_
  rec$cutoff <- chosen
  rec$objective_value <- if (obj$direction == "min") -best else best
  rec$n_tied <- as.double(length(tied))
  rec$tied <- list(cutoff = ifelse(is.finite(tied), tied, NA_real_),
                   cutoff_kind = ifelse(is.finite(tied), "real", ifelse(tied > 0, "alert_none", "alert_all")),
                   objective_value = rep(rec$objective_value, length(tied)))
  rec$bounds <- .ig_cls_bound_summary(bounds, cfun(ec$tp, ec$fn, ec$fp, ec$tn), 1L)
  if (length(bounds)) {
    # report attainability over the scan, but satisfaction at the choice
    scan_sum <- .ig_cls_bound_summary(bounds, gv)
    rec$bounds$n_satisfied <- scan_sum$n_satisfied
    rec$bounds$n_points <- scan_sum$n_points
    rec$bounds$best_value <- scan_sum$best_value
    rec$bounds$best_status <- scan_sum$best_status
  }
  rec
}

#' Select a cutoff by an objective or constraints
#'
#' Chooses the operating point of each model by a declared rule and evaluates
#' it with [evaluate_cutoff()]. For scored data the search is exhaustive over
#' every distinct decision set of the original scores (including alert none and
#' alert all); for a theoretical binormal model the optimum over real-valued
#' cutoffs is found numerically (and analytically for Youden's J). Ties are
#' reported in full and broken by the declared `tie_break` among actual
#' operating points.
#'
#' No policy threshold is built in: every constraint is a user decision.
#' Conflicting constraints, an objective that is undefined everywhere or a
#' supremum that no cutoff attains give a result with status `"infeasible"`
#' and a per-constraint feasibility report instead of a number.
#'
#' A cutoff chosen on the evaluated data is descriptive: the result is labelled
#' in-sample and its intervals are conditional on the selected cutoff. To
#' evaluate a rule chosen on tuning data, pass the tuning result as `cutoff` to
#' [evaluate_cutoff()] on the evaluation specification.
#'
#' @param spec A `classification_spec` (scored data or theory).
#' @param objective One of `"max_youden"`, `"max_f1"`, `"max_mcc"`,
#'   `"constraint"`, `"capacity"`, `"cost_utility"`; `NULL` uses the selection
#'   rule already in `spec`.
#' @param constraints,optimise,capacity,costs,tie_break,prevalence_basis See
#'   [cutoff_rule()].
#' @param ... Evaluation settings passed to [evaluate_cutoff()]
#'   (`partial_auc`, `threshold_probability`, `calibration_bins`,
#'   `include_curves`).
#' @return A `classification_result`. The specification it was evaluated
#'   with (the input with the selection rule) is returned by
#'   [classification_evaluated_spec()]; pair the two in an
#'   [analysis_envelope()].
#' @examples
#' ex <- classification_examples("scored_synthetic")
#' res <- select_cutoff(ex$spec, "max_youden")
#' classification_table(res, "selection")
#'
#' # Highest TPR among cutoffs with FPR at most 10 %
#' res2 <- select_cutoff(ex$spec, "constraint", constraints = list(fpr_max = 0.1))
#' res2$models[[1]]$operating_point$cutoff
#'
#' # Theory: analytic Youden optimum of an unequal-variance binormal model
#' th <- classification_spec(theory = list(auc = 0.8, sigma = 1.5), cutoff = 0, prevalence = 0.1)
#' select_cutoff(th, "max_youden")$models[[1]]$operating_point$cutoff
#' @seealso [evaluate_cutoff()], [cutoff_rule()]
#' @concept classification
#' @seealso
#' * [classification_spec()]
#' * [evaluate_cutoff()]
#' * [plot_classification()]
#' @export
select_cutoff <- function(spec, objective = NULL, constraints = NULL, optimise = NULL,
                          capacity = NULL, costs = NULL, tie_break = NULL,
                          prevalence_basis = NULL, ...) {
  .ig_cls_check_spec(spec)
  if (is.null(objective)) {
    if (identical(spec$operating_rule$selection, "manual")) {
      .ig_abort("invalid_input", "give a selection `objective` (the specification holds a manual cutoff)")
    }
    if (!is.null(constraints) || !is.null(optimise) || !is.null(capacity) || !is.null(costs) ||
        !is.null(tie_break) || !is.null(prevalence_basis)) {
      .ig_abort("invalid_input", "rule settings need an `objective`; or omit them to use the rule in `spec`")
    }
    return(evaluate_cutoff(spec, ...))
  }
  if (!.ig_is_string(objective) || !objective %in% .ig_cls_selection_codes) {
    .ig_abort("invalid_input", sprintf("`objective` must be one of: %s",
                                       paste(.ig_cls_selection_codes, collapse = ", ")))
  }
  sel_data <- if (identical(spec$input_kind, "theoretical_binormal")) "predeclared_rule" else NULL
  rule <- cutoff_rule(selection = objective, selection_data = sel_data, constraints = constraints,
                      optimise = optimise, capacity = capacity, costs = costs, tie_break = tie_break,
                      prevalence_basis = prevalence_basis)
  if (!is.null(sel_data)) {
    # recorded as the theoretical default, not as a user statement
    r <- unclass(rule)
    r$defaults_applied <- sort(unique(c(r$defaults_applied, "selection_data")))
    rule <- structure(r, class = "ig_cutoff_rule")
  }
  evaluate_cutoff(update(spec, rule = rule), ...)
}
