# Plot-ready tables for the linked classification panels ---------------------------------
#
# classification_plot_data() turns a classification_result into the data of
# one panel: plain data frames plus notes, no rendering. Everything comes from
# the result itself (decision sets, operating point, metric tables), so every
# panel, table and report built from the same result agrees. Rules carried
# over from the source audit:
#
# * operating points are the attained decision sets; PR points are never
#   joined by straight segments that imply unattainable operating points;
# * display bins of the score distribution are display only;
# * no value is clamped onto a log axis: points that a log axis cannot show
#   (zero, infinite or undefined values) are counted and reported;
# * prevalence-dependent panels state the prevalence used and whether it is
#   observed or a scenario; realised varying-cutoff rules are separate
#   markers and are not placed on a fixed-threshold curve.

.ig_cls_panels <- c("confusion_table", "mosaic", "score_distribution", "roc", "pr",
                    "metric_vs_cutoff", "prevalence_curve", "metric_pair", "metric_network",
                    "model_comparison", "sweep_table", "current_table", "calibration", "net_benefit")

.ig_cls_families <- list(
  within_truth = c("tpr", "fpr"),
  within_alerts = c("ppv", "npv"),
  summary = c("f1", "mcc", "youden_j"),
  burden = c("alert_rate", "tp_share", "fp_share"),
  likelihood_ratios = c("lr_pos", "lr_neg")
)

.ig_cls_short_label <- function(ids, lang) {
  dict <- .ig_cls_dictionary()$metrics
  vapply(ids, function(id) {
    e <- dict[[id]]
    if (is.null(e)) id else (e$label[[lang]] %||% e$label$en)$short
  }, character(1), USE.NAMES = FALSE)
}

# Curve table of one model on the requested basis (observed or scenario).
.ig_cls_model_curve <- function(x, m, basis) {
  if (is.null(m$curve)) return(NULL)
  tab <- .ig_cls_frame(m$curve, "curve")
  if (basis == "scenario") {
    n_all <- tab$tp[1] + tab$fn[1] + tab$fp[1] + tab$tn[1]
    tr <- .ig_cls_transport(tab$tp, tab$fn, tab$fp, tab$tn, x$prevalence$value, n_all)$counts
    tab[, c("tp", "fn", "fp", "tn")] <- tr
    tab$alerts <- tab$tp + tab$fp
  }
  tab
}

# Cells of the operating point on the requested basis.
.ig_cls_model_cells <- function(m, basis) {
  if (basis == "scenario") {
    if (is.null(m$scenario) || !identical(m$scenario$status, "ok")) return(NULL)
    c <- m$scenario$counts
    return(list(counts = c(tp = c$tp, fn = c$fn, fp = c$fp, tn = c$tn), count_kind = "expected",
                metrics = .ig_cls_frame(m$scenario$metrics, "metrics")))
  }
  op <- m$operating_point
  if (is.null(op)) return(NULL)
  list(counts = c(tp = op$tp, fn = op$fn, fp = op$fp, tn = op$tn), count_kind = op$count_kind,
       metrics = .ig_cls_frame(m$metrics, "metrics"))
}

.ig_cls_metric_value <- function(mt, id) {
  r <- mt[mt$metric_id == id, , drop = FALSE]
  list(value = r$value[1], status = r$status[1], reason_code = r$reason_code[1],
       denominator = r$denominator[1], denominator_value = r$denominator_value[1])
}

.ig_cls_check_levels <- function(x, name, open = TRUE) {
  ok <- is.numeric(x) && length(x) >= 1L && all(is.finite(x)) && all(x > 0 & x < 1)
  if (!ok) .ig_abort("invalid_input", sprintf("`%s` must be numeric values in (0, 1)", name))
  sort(unique(as.double(x)))
}

.ig_cls_check_realised <- function(rules) {
  if (is.null(rules)) return(NULL)
  if (!is.data.frame(rules) || !all(c("id", "tp", "fn", "fp", "tn") %in% names(rules))) {
    .ig_abort("invalid_input", "`realised_rules` must be a data frame with id, tp, fn, fp, tn (and optionally label)")
  }
  for (i in seq_len(nrow(rules))) {
    .ig_cls_check_counts(unlist(rules[i, c("tp", "fn", "fp", "tn")]), "frequency_weighted")
  }
  if (is.null(rules$label)) rules$label <- rules$id
  rules
}

#' Plot-ready data for the linked classification panels
#'
#' Returns the data of one classification panel for a `classification_result`
#' as plain data frames with notes; no plot is drawn (static renderers and the
#' application draw from these tables). All values come from the result, so
#' panels, tables and reports of the same result agree.
#'
#' Panels: `"confusion_table"` (cells, totals and headline rates with
#' denominators), `"mosaic"` (rectangles for joint and within-class widths),
#' `"score_distribution"` (exact counts per distinct score plus display-only
#' bins; densities for a theoretical model), `"roc"` (attained operating
#' points, current point, chance diagonal, iso-PPV references at a stated
#' prevalence, log-axis boundary counts, realised rules), `"pr"` (attained
#' points without interpolation, prevalence baseline, iso-F1 references),
#' `"metric_vs_cutoff"` (five metric families along the decision sets),
#' `"prevalence_curve"` (PPV and NPV versus prevalence at fixed TPR/FPR),
#' `"metric_pair"` (any two operating-point metrics), `"metric_network"`
#' (metric dependency graph from the dictionary), `"model_comparison"` (ROC
#' curves and operating points of all models evaluated on the same units,
#' with the paired DeLong comparisons; `model_id` is the focal model),
#' `"sweep_table"` and `"current_table"`. Probability-score panels
#' `"calibration"` (retained grouped observed rates and diagnostic summaries) and
#' `"net_benefit"` (explicit decision probabilities with alert-all and alert-none
#' references) use observed data, never prevalence scenarios.
#'
#' @param x A `classification_result`.
#' @param panel Panel name (see Details).
#' @param model_id Model id (default: the first model).
#' @param basis `"observed"` (default) or `"scenario"` (expected counts at the
#'   result's scenario prevalence) for prevalence-dependent panels.
#' @param lang Label language, `"en"` or `"de"`.
#' @param ppv_levels Iso-PPV reference levels of the ROC panel.
#' @param f1_levels Iso-F1 reference levels of the PR panel.
#' @param ppv_target Optional PPV target for the prevalence panel (no
#'   default: targets are user decisions).
#' @param metric_x,metric_y Metric ids of the metric-pair panel.
#' @param realised_rules Optional data frame of realised varying-cutoff rules
#'   (`id`, `label`, `tp`, `fn`, `fp`, `tn`) for the ROC, PR and metric-pair
#'   panels; they are drawn as separate markers.
#' @param bins Number of display bins of the score distribution.
#' @return A list with `panel`, `model_id`, `basis`, `data` (named list of
#'   data frames) and `notes` (character).
#' @examples
#' ex <- classification_examples("scored_synthetic")
#' res <- evaluate_cutoff(ex$spec)
#' roc <- classification_plot_data(res, "roc")
#' head(roc$data$points)
#' roc$notes
#' @concept classification
#' @seealso
#' * [classification_spec()]
#' * [evaluate_cutoff()]
#' * [plot_classification()]
#' @export
classification_plot_data <- function(x, panel = .ig_cls_panels, model_id = NULL,
                                     basis = c("observed", "scenario"), lang = c("en", "de"),
                                     ppv_levels = c(0.25, 0.5, 0.75), f1_levels = c(0.2, 0.4, 0.6, 0.8),
                                     ppv_target = NULL, metric_x = "fpr", metric_y = "tpr",
                                     realised_rules = NULL, bins = 40) {
  if (!inherits(x, "classification_result")) {
    .ig_abort("invalid_input", "`x` must be a classification_result (see evaluate_cutoff())")
  }
  panel <- match.arg(panel)
  basis <- match.arg(basis)
  lang <- match.arg(lang)
  ids <- vapply(x$models, function(m) m$model_id, "")
  if (is.null(model_id)) model_id <- ids[1]
  j <- match(model_id, ids)
  if (!.ig_is_string(model_id) || is.na(j)) .ig_abort("invalid_input", "`model_id` is not a model of the result")
  m <- x$models[[j]]
  if (basis == "scenario" && (!identical(x$prevalence$kind, "scenario") ||
                              identical(x$input_kind, "theoretical_binormal"))) {
    .ig_abort("invalid_input", "basis 'scenario' needs a result of scored data or counts with a scenario prevalence")
  }
  ppv_levels <- .ig_cls_check_levels(ppv_levels, "ppv_levels")
  f1_levels <- .ig_cls_check_levels(f1_levels, "f1_levels")
  realised_rules <- .ig_cls_check_realised(realised_rules)
  out <- switch(panel,
    confusion_table = .ig_cls_pd_confusion(x, m, basis, lang),
    mosaic = .ig_cls_pd_mosaic(x, m, basis, lang),
    score_distribution = .ig_cls_pd_distribution(x, m, bins),
    roc = .ig_cls_pd_roc(x, m, ppv_levels, realised_rules, basis),
    pr = .ig_cls_pd_pr(x, m, basis, f1_levels, realised_rules),
    metric_vs_cutoff = .ig_cls_pd_families(x, m, basis, lang),
    prevalence_curve = .ig_cls_pd_prevalence(x, m, ppv_target),
    metric_pair = .ig_cls_pd_pair(x, m, basis, metric_x, metric_y, realised_rules, lang),
    metric_network = .ig_cls_pd_network(m, basis, lang),
    model_comparison = .ig_cls_pd_comparison(x, m),
    sweep_table = .ig_cls_pd_sweep(x, m, basis),
    current_table = .ig_cls_pd_current(x, m, basis, lang),
    calibration = .ig_cls_pd_probability(x, m, panel, basis, lang),
    net_benefit = .ig_cls_pd_probability(x, m, panel, basis, lang))
  c(list(panel = panel, model_id = model_id, basis = basis), out)
}

.ig_cls_pd_need_cells <- function(m, basis) {
  cells <- .ig_cls_model_cells(m, basis)
  if (is.null(cells)) {
    why <- if (basis == "scenario" && !is.null(m$scenario) && !identical(m$scenario$status, "ok")) {
      sprintf("scenario undefined: %s", .ig_cls_reason_text(m$scenario$reason_code))
    } else sprintf("status: %s", .ig_cls_status_code(m$status))
    .ig_abort("infeasible", sprintf("model '%s' has no operating point on the %s basis (%s)",
                                    m$model_id, basis, why))
  }
  cells
}

.ig_cls_pd_need_curve <- function(x, m, basis, view) {
  v <- m$views[[view]]
  if (!is.null(v) && !isTRUE(v$available)) .ig_abort("unsupported", .ig_cls_reason_text(v$reason_code))
  tab <- .ig_cls_model_curve(x, m, basis)
  if (is.null(tab)) .ig_abort("unsupported", "the result holds no decision sets (computed with include_curves = FALSE)")
  tab
}

.ig_cls_pd_prevalence_used <- function(x, cells, basis) {
  P <- cells$counts[["tp"]] + cells$counts[["fn"]]
  n <- sum(cells$counts)
  kind <- if (basis == "scenario" || identical(x$input_kind, "theoretical_binormal")) {
    if (identical(x$input_kind, "theoretical_binormal")) "theoretical_assumed" else
      paste0("scenario_", x$prevalence$source)
  } else "observed_sample"
  list(value = if (n > 0) P / n else NA_real_, kind = kind)
}

.ig_cls_pd_confusion <- function(x, m, basis, lang) {
  cells <- .ig_cls_pd_need_cells(m, basis)
  c <- cells$counts
  expected <- cells$count_kind == "expected"
  lab <- if (lang == "de") c(alert = "Alarm", no_alert = "Kein Alarm", event = "Ereignis", non_event = "Kein Ereignis") else
    c(alert = "Alert", no_alert = "No alert", event = "Event", non_event = "No event")
  cell_df <- data.frame(
    cell = c("tp", "fp", "fn", "tn"), row = c("alert", "alert", "no_alert", "no_alert"),
    column = c("event", "non_event", "event", "non_event"),
    row_label = lab[c("alert", "alert", "no_alert", "no_alert")],
    column_label = lab[c("event", "non_event", "event", "non_event")],
    count = unname(c[c("tp", "fp", "fn", "tn")]), count_kind = cells$count_kind,
    expected_marker = if (expected) "\u2248" else "", stringsAsFactors = FALSE, row.names = NULL)
  totals <- data.frame(
    margin = c("alerts", "non_alerts", "events", "non_events", "total"),
    value = c(c[["tp"]] + c[["fp"]], c[["fn"]] + c[["tn"]], c[["tp"]] + c[["fn"]], c[["fp"]] + c[["tn"]], sum(c)),
    stringsAsFactors = FALSE)
  tiles_ids <- c("tpr", "fpr", "ppv", "alerts_per_tp")
  mt <- cells$metrics
  tiles <- do.call(rbind, lapply(tiles_ids, function(id) {
    v <- .ig_cls_metric_value(mt, id)
    data.frame(metric_id = id, short = .ig_cls_short_label(id, lang), value = v$value, status = v$status,
               reason = .ig_cls_reason_text(v$reason_code, lang), denominator = v$denominator,
               denominator_value = v$denominator_value, stringsAsFactors = FALSE)
  }))
  op <- m$operating_point
  notes <- c(sprintf("rule: %s (%s)", op$rule_text, m$direction),
             if (expected) "\u2248 marks expected counts, not observed units",
             if (basis == "scenario") .ig_cls_scenario_assumption)
  list(data = list(cells = cell_df, totals = totals, tiles = tiles), notes = notes)
}

.ig_cls_pd_mosaic <- function(x, m, basis, lang) {
  cells <- .ig_cls_pd_need_cells(m, basis)
  c <- cells$counts
  P <- c[["tp"]] + c[["fn"]]
  N <- c[["fp"]] + c[["tn"]]
  n <- P + N
  rects <- list()
  add <- function(mode, cell, x0, x1, y0, y1, share) {
    rects[[length(rects) + 1L]] <<- data.frame(mode = mode, cell = cell, x0 = x0, x1 = x1, y0 = y0, y1 = y1,
                                              area_share = share, stringsAsFactors = FALSE)
  }
  notes <- character()
  for (mode in c("joint", "within_class")) {
    w_ev <- if (mode == "joint") P / n else if (P > 0 && N > 0) 0.5 else as.numeric(P > 0)
    if (P > 0) {
      tpr <- c[["tp"]] / P
      add(mode, "tp", 0, w_ev, 0, tpr, if (mode == "joint") c[["tp"]] / n else NA_real_)
      add(mode, "fn", 0, w_ev, tpr, 1, if (mode == "joint") c[["fn"]] / n else NA_real_)
    }
    if (N > 0) {
      fpr <- c[["fp"]] / N
      add(mode, "fp", w_ev, 1, 0, fpr, if (mode == "joint") c[["fp"]] / n else NA_real_)
      add(mode, "tn", w_ev, 1, fpr, 1, if (mode == "joint") c[["tn"]] / n else NA_real_)
    }
  }
  if (P == 0) notes <- c(notes, .ig_cls_reason_text("no_events", lang))
  if (N == 0) notes <- c(notes, .ig_cls_reason_text("no_non_events", lang))
  ppv <- .ig_cls_metric_value(cells$metrics, "ppv")
  caption <- data.frame(quantity = "ppv", value = ppv$value, status = ppv$status,
                        reason = .ig_cls_reason_text(ppv$reason_code, lang),
                        formula = "TP area / (TP area + FP area) in the joint view", stringsAsFactors = FALSE)
  notes <- c(notes, "joint: column widths are the class shares (areas are joint shares)",
             "within class: equal widths; heights are class-conditional rates, which stay fixed under a prevalence change only while the class-conditional score distributions stay fixed")
  list(data = list(rectangles = do.call(rbind, rects), caption = caption), notes = notes)
}

.ig_cls_pd_distribution <- function(x, m, bins) {
  if (!is.numeric(bins) || length(bins) != 1L || !is.finite(bins) || bins < 2 || bins != round(bins)) {
    .ig_abort("invalid_input", "`bins` must be a whole number >= 2")
  }
  op <- m$operating_point
  cut <- data.frame(cutoff = if (is.null(op)) NA_real_ else op$cutoff %||% NA_real_,
                    cutoff_kind = if (is.null(op)) NA_character_ else op$cutoff_kind,
                    comparator = m$comparator, stringsAsFactors = FALSE)
  if (identical(x$input_kind, "theoretical_binormal")) {
    th <- x$theory
    pi <- x$prevalence$value
    lo <- min(-4, th$d - 4 * th$sigma)
    hi <- max(4, th$d + 4 * th$sigma)
    z <- seq(lo, hi, length.out = 401)
    dens <- data.frame(z = z, density_events = stats::dnorm(z, th$d, th$sigma),
                       density_non_events = stats::dnorm(z))
    dens$joint_events <- pi * dens$density_events
    dens$joint_non_events <- (1 - pi) * dens$density_non_events
    dens$alerted <- if (is.null(op)) NA else if (op$cutoff_kind == "alert_all") TRUE else
      if (op$cutoff_kind == "alert_none") FALSE else z >= op$cutoff
    return(list(data = list(densities = dens, cutoff = cut),
                notes = c("class densities of the binormal model (S0 ~ N(0, 1), S1 ~ N(d, sigma^2)); joint_* columns are weighted by the assumed prevalence",
                          "the x axis is the score standardised to the non-event distribution (z), not a probability",
                          "counts and metrics use analytic tails; the grid is for display only")))
  }
  tab <- .ig_cls_pd_need_curve(x, m, "observed", "score_distribution")
  sc <- tab[tab$set_index > 0, , drop = FALSE]
  prev <- tab[match(sc$set_index - 1, tab$set_index), ]
  ev <- sc$tp - prev$tp
  ne <- sc$fp - prev$fp
  P <- tab$tp[1] + tab$fn[1]
  N <- tab$fp[1] + tab$tn[1]
  alerted <- if (is.null(op)) NA else sc$set_index <= op$set_index
  scores <- data.frame(score = sc$boundary_score, events = ev, non_events = ne,
                       share_of_events = if (P > 0) ev / P else NA_real_,
                       share_of_non_events = if (N > 0) ne / N else NA_real_,
                       alerted = alerted, stringsAsFactors = FALSE)
  rng <- range(scores$score)
  brk <- if (rng[1] == rng[2]) c(rng[1] - 0.5, rng[1] + 0.5) else seq(rng[1], rng[2], length.out = bins + 1L)
  bin <- findInterval(scores$score, brk, rightmost.closed = TRUE, all.inside = TRUE)
  nb <- length(brk) - 1L
  agg <- function(v) vapply(seq_len(nb), function(b) sum(v[bin == b]), double(1))
  display <- data.frame(bin_lower = brk[-length(brk)], bin_upper = brk[-1], events = agg(ev),
                        non_events = agg(ne), stringsAsFactors = FALSE)
  display$share_of_events <- if (P > 0) display$events / P else NA_real_
  display$share_of_non_events <- if (N > 0) display$non_events / N else NA_real_
  list(data = list(scores = scores, display_bins = display, cutoff = cut),
       notes = c(sprintf("display bins of width %s are display only; counts and metrics use the original scores",
                         format(signif(diff(brk[1:2]), 4))),
                 if (identical(m$score_scale, "probability")) "x axis: declared predicted probability of the event (calibration is not implied)" else
                   "x axis: raw model score (no probability or calibration is claimed)"))
}

.ig_cls_pd_roc <- function(x, m, ppv_levels, realised, basis) {
  tab <- .ig_cls_pd_need_curve(x, m, "observed", "roc")
  P <- tab$tp[1] + tab$fn[1]
  N <- tab$fp[1] + tab$tn[1]
  pts <- data.frame(set_index = tab$set_index, cutoff = tab$cutoff, cutoff_kind = tab$cutoff_kind,
                    fpr = tab$fp / N, tpr = tab$tp / P, stringsAsFactors = FALSE)
  if (!identical(x$input_kind, "theoretical_binormal")) {
    # a tied score group moves both coordinates at once: a diagonal segment
    pts$segment_from_previous <- c(NA, ifelse(diff(tab$tp) > 0 & diff(tab$fp) > 0, "tie_diagonal",
                                              ifelse(diff(tab$tp) > 0, "vertical", "horizontal")))
  }
  op <- m$operating_point
  current <- if (is.null(op) || op$n_events == 0 || op$n_non_events == 0) NULL else
    data.frame(fpr = op$fp / op$n_non_events, tpr = op$tp / op$n_events, marker = "current_operating_point")
  pi_ref <- if (identical(x$prevalence$kind, "scenario")) x$prevalence$value else P / (P + N)
  pi_kind <- if (identical(x$prevalence$kind, "scenario")) paste0("scenario_", x$prevalence$source) else if (identical(x$input_kind,"theoretical_binormal")) "theoretical_assumed" else if (identical(x$prevalence$sampling,"case_control")) "sample_prevalence_fixed_by_design" else "observed_sample"
  iso <- NULL
  if (pi_ref > 0 && pi_ref < 1) {
    iso <- do.call(rbind, lapply(ppv_levels, function(L) {
      slope <- (1 - pi_ref) / pi_ref * L / (1 - L)
      f <- seq(0, min(1, 1 / slope), length.out = 201)
      t <- f * slope
      keep <- t <= 1
      data.frame(ppv_level = L, prevalence = pi_ref, prevalence_kind = pi_kind, fpr = f[keep], tpr = t[keep])
    }))
  }
  rr <- NULL
  if (!is.null(realised)) {
    rr <- data.frame(id = realised$id, label = realised$label,
                     fpr = realised$fp / (realised$fp + realised$tn),
                     tpr = realised$tp / (realised$tp + realised$fn), marker = "realised_rule",
                     stringsAsFactors = FALSE)
  }
  log_info <- data.frame(n_points = nrow(pts), n_fpr_zero = sum(pts$fpr == 0),
                         min_positive_fpr = if (any(pts$fpr > 0)) min(pts$fpr[pts$fpr > 0]) else NA_real_)
  list(data = list(points = pts, current = current, diagonal = data.frame(fpr = c(0, 1), tpr = c(0, 1)),
                   iso_ppv = iso, realised_rules = rr, log_axis = log_info),
       notes = c("the ROC curve and its AUC do not depend on the selected cutoff",
                 if (!is.null(iso)) sprintf("iso-PPV references at prevalence %s (%s); straight lines through the origin on linear axes, curves on a log-FPR axis",
                                            format(signif(pi_ref, 4)), pi_kind),
                 sprintf("%d operating point(s) have FPR = 0 and cannot be drawn on a log-FPR axis; show them in a boundary strip", log_info$n_fpr_zero),
                 if (!is.null(rr)) "realised varying-cutoff rules are separate markers and need not lie on the fixed-threshold curve",
                 if (basis == "scenario") "the ROC curve is the same under a prevalence scenario (fixed class-conditional distributions)"))
}

.ig_cls_pd_pr <- function(x, m, basis, f1_levels, realised) {
  tab <- .ig_cls_pd_need_curve(x, m, basis, "pr")
  P <- tab$tp[1] + tab$fn[1]
  n <- P + tab$fp[1] + tab$tn[1]
  has_alert <- tab$alerts > 0
  pts <- data.frame(set_index = tab$set_index, cutoff = tab$cutoff, cutoff_kind = tab$cutoff_kind,
                    recall = tab$tp / P, precision = ifelse(has_alert, tab$tp / tab$alerts, NA_real_),
                    stringsAsFactors = FALSE)[has_alert, , drop = FALSE]
  op <- m$operating_point
  cells <- .ig_cls_model_cells(m, basis)
  current <- NULL
  if (!is.null(cells)) {
    c <- cells$counts
    if (c[["tp"]] + c[["fp"]] > 0 && c[["tp"]] + c[["fn"]] > 0) {
      current <- data.frame(recall = c[["tp"]] / (c[["tp"]] + c[["fn"]]),
                            precision = c[["tp"]] / (c[["tp"]] + c[["fp"]]), marker = "current_operating_point")
    }
  }
  pi_used <- P / n
  iso <- do.call(rbind, lapply(f1_levels, function(f) {
    r <- seq(f / (2 - f), 1, length.out = 201)
    p <- f * r / (2 * r - f)
    p[1] <- 1 # Exact analytic intersection; avoid dropping it due to rounding.
    keep <- p <= 1
    data.frame(f1_level = f, recall = r[keep], precision = p[keep])
  }))
  rr <- NULL
  if (!is.null(realised) && basis == "observed") {
    rr <- data.frame(id = realised$id, label = realised$label, recall = realised$tp / (realised$tp + realised$fn),
                     precision = ifelse(realised$tp + realised$fp > 0, realised$tp / (realised$tp + realised$fp), NA_real_),
                     marker = "realised_rule", stringsAsFactors = FALSE)
  }
  list(data = list(points = pts, current = current,
                   baseline = data.frame(precision = pi_used,
                                         prevalence_kind = if (basis == "scenario") paste0("scenario_", x$prevalence$source) else
                                           if (identical(x$input_kind, "theoretical_binormal")) "theoretical_assumed" else "observed_sample"),
                   iso_f1 = iso, realised_rules = rr,
                   excluded = data.frame(n_zero_alert_points = sum(!has_alert), reason_code = "no_alerts")),
       notes = c("points are attained operating points; they are not joined by straight lines (linear PR interpolation implies operating points that no cutoff attains)",
                 "average precision (step sum) and the trapezoidal PR area are different quantities",
                 sprintf("precision depends on prevalence; here %s (%s)", format(signif(pi_used, 4)),
                         if (basis == "scenario") "scenario" else "as evaluated"),
                 if (!is.null(realised) && basis == "scenario") "realised rules are not shown on the scenario basis: their observed precision belongs to the observed prevalence"))
}

.ig_cls_pd_families <- function(x, m, basis, lang) {
  tab <- .ig_cls_pd_need_curve(x, m, basis, "sweep")
  ids <- unlist(.ig_cls_families, use.names = FALSE)
  vals <- .ig_cls_metrics_vec(tab$tp, tab$fn, tab$fp, tab$tn)
  long <- do.call(rbind, lapply(names(.ig_cls_families), function(fam) {
    do.call(rbind, lapply(.ig_cls_families[[fam]], function(id) {
      data.frame(family = fam, metric_id = id, short = .ig_cls_short_label(id, lang),
                 set_index = tab$set_index, cutoff = tab$cutoff, cutoff_kind = tab$cutoff_kind,
                 on_cutoff_axis = !is.na(tab$cutoff), value = vals[[id]]$value, status = vals[[id]]$status,
                 reason_code = vals[[id]]$reason, stringsAsFactors = FALSE)
    }))
  }))
  lr <- long[long$family == "likelihood_ratios", ]
  domains <- data.frame(family = names(.ig_cls_families),
                        lower = c(0, 0, -1, 0, NA), upper = c(1, 1, 1, NA, NA),
                        scale = c("linear", "linear", "linear", "linear", "log"),
                        stringsAsFactors = FALSE)
  list(data = list(values = long, domains = domains,
                   log_axis = data.frame(family = "likelihood_ratios",
                                         n_positive_infinity = sum(lr$status == "positive_infinity"),
                                         n_zero = sum(lr$status == "ok" & lr$value == 0),
                                         n_undefined = sum(lr$status == "undefined"))),
       notes = c("summary family uses the full range [-1, 1]; nothing is clipped",
                 "on the log scale of the likelihood-ratio family, +infinity, zero and undefined values are counted in log_axis and not drawn at an axis bound",
                 "alert none / alert all without a numeric cutoff have on_cutoff_axis = FALSE",
                 if (basis == "observed") "TPR and FPR stay fixed under a prevalence change only while the class-conditional score distributions stay fixed"))
}

.ig_cls_pd_prevalence <- function(x, m, ppv_target) {
  cells <- .ig_cls_pd_need_cells(m, "observed")
  c <- cells$counts
  P <- c[["tp"]] + c[["fn"]]
  N <- c[["fp"]] + c[["tn"]]
  if (P == 0 || N == 0) .ig_abort("unsupported", .ig_cls_reason_text(if (P == 0) "no_events" else "no_non_events"))
  tpr <- c[["tp"]] / P
  fpr <- c[["fp"]] / N
  sample_pi <- P / (P + N)
  grid <- sort(unique(c(10^seq(-3, log10(0.999), length.out = 200),
                        if (!identical(x$input_kind, "theoretical_binormal")) sample_pi,
                        if (identical(x$prevalence$kind, "scenario")) x$prevalence$value)))
  pv <- predictive_values_at_prevalence(tpr, fpr, grid)
  markers <- data.frame(kind = character(), prevalence = double(), stringsAsFactors = FALSE)
  if (!identical(x$input_kind, "theoretical_binormal")) {
    markers <- rbind(markers, data.frame(kind = if (identical(x$prevalence$sampling, "case_control"))
      "sample_prevalence_fixed_by_design" else "observed_sample", prevalence = sample_pi))
  }
  if (identical(x$prevalence$kind, "scenario")) {
    markers <- rbind(markers, data.frame(kind = paste0("scenario_", x$prevalence$source), prevalence = x$prevalence$value))
  }
  target <- NULL
  if (!is.null(ppv_target)) target <- prevalence_for_ppv_target(tpr, fpr, ppv_target)
  list(data = list(curve = pv, markers = markers, target = target,
                   rates = data.frame(tpr = tpr, fpr = fpr)),
       notes = c("PPV and NPV at other prevalences with TPR and FPR of the current operating point held fixed (fixed class-conditional score distributions); a scenario, not a measurement",
                 "no sampling uncertainty of TPR and FPR is shown",
                 "a log-prevalence axis starts at 0.001; PPV is undefined at prevalence 0"))
}

.ig_cls_pd_pair <- function(x, m, basis, mx, my, realised, lang) {
  ok_ids <- .ig_cls_operating_ids()
  if (!.ig_is_string(mx) || !mx %in% ok_ids || !.ig_is_string(my) || !my %in% ok_ids) {
    .ig_abort("invalid_input", "`metric_x` and `metric_y` must be operating-point metric ids")
  }
  tab <- .ig_cls_pd_need_curve(x, m, basis, "sweep")
  vals <- .ig_cls_metrics_vec(tab$tp, tab$fn, tab$fp, tab$tn)
  pts <- data.frame(set_index = tab$set_index, cutoff = tab$cutoff, cutoff_kind = tab$cutoff_kind,
                    x = vals[[mx]]$value, x_status = vals[[mx]]$status,
                    y = vals[[my]]$value, y_status = vals[[my]]$status, stringsAsFactors = FALSE)
  dep <- metric_dictionary()
  pdep <- function(id) dep$prevalence_dependent[dep$id == id]
  cells <- .ig_cls_model_cells(m, basis)
  current <- NULL
  if (!is.null(cells)) {
    v <- .ig_cls_metrics_vec(cells$counts[["tp"]], cells$counts[["fn"]], cells$counts[["fp"]], cells$counts[["tn"]])
    current <- data.frame(x = v[[mx]]$value, x_status = v[[mx]]$status, y = v[[my]]$value,
                          y_status = v[[my]]$status, marker = "current_operating_point")
  }
  rr <- NULL
  suppress <- basis == "scenario" && (pdep(mx) || pdep(my))
  if (!is.null(realised) && !suppress) {
    rv <- .ig_cls_metrics_vec(realised$tp, realised$fn, realised$fp, realised$tn)
    rr <- data.frame(id = realised$id, label = realised$label, x = rv[[mx]]$value, x_status = rv[[mx]]$status,
                     y = rv[[my]]$value, y_status = rv[[my]]$status, marker = "realised_rule",
                     stringsAsFactors = FALSE)
  }
  list(data = list(points = pts, current = current, realised_rules = rr,
                   axes = data.frame(axis = c("x", "y"), metric_id = c(mx, my),
                                     short = .ig_cls_short_label(c(mx, my), lang),
                                     prevalence_dependent = c(pdep(mx), pdep(my)),
                                     n_not_drawable = c(sum(pts$x_status != "ok"), sum(pts$y_status != "ok")),
                                     stringsAsFactors = FALSE)),
       notes = c("points with an undefined or infinite value on an axis are counted in axes$n_not_drawable and not drawn",
                 if (suppress && !is.null(realised)) "realised rules are hidden: a prevalence-dependent axis on the scenario basis would mix observed and expected counts"))
}

.ig_cls_pd_network <- function(m, basis, lang) {
  dict <- .ig_cls_dictionary()$metrics
  ids <- .ig_cls_operating_ids()
  cells <- .ig_cls_model_cells(m, basis)
  mt <- if (is.null(cells)) NULL else cells$metrics
  nodes <- do.call(rbind, lapply(ids, function(id) {
    e <- dict[[id]]
    v <- if (is.null(mt)) list(value = NA_real_, status = "not_computed", reason_code = "not_requested") else
      .ig_cls_metric_value(mt, id)
    data.frame(id = id, short = .ig_cls_short_label(id, lang), group = e$group,
               prevalence_dependent = !identical(e$prevalence_dependence,
                                                 "invariant_given_fixed_class_conditional_distributions"),
               value = v$value, status = v$status, stringsAsFactors = FALSE)
  }))
  edges <- do.call(rbind, lapply(ids, function(id) {
    e <- dict[[id]]
    a <- unlist(e$inputs)
    b <- unlist(e$bayes_inputs)
    rbind(if (length(a)) data.frame(from = a, to = id, kind = "input", stringsAsFactors = FALSE),
          if (length(b)) data.frame(from = b, to = id, kind = "bayes", stringsAsFactors = FALSE))
  }))
  list(data = list(nodes = nodes, edges = edges),
       notes = c("edges come from the metric dictionary (inputs; Bayes edges link TPR, FPR and prevalence to PPV, and TNR, FNR and prevalence to NPV)",
                 "'prevalence-free' (class-conditional) metrics are invariant only under fixed class-conditional score distributions"))
}

# ROC curves and operating points of every model of a scored result. The
# focal model's operating point is the current one; the others are comparison
# markers. Operating points of different models at the same numerical cutoff
# are not equivalent, so each model keeps its own rule.
.ig_cls_pd_comparison <- function(x, m) {
  if (length(x$models) < 2L) {
    .ig_abort("unsupported", "a model comparison needs a result with at least two models")
  }
  curves <- list()
  points <- list()
  aucs <- list()
  for (mm in x$models) {
    focal <- identical(mm$model_id, m$model_id)
    if (isTRUE(mm$views$roc$available) && !is.null(mm$curve)) {
      tab <- .ig_cls_frame(mm$curve, "curve")
      P <- tab$tp[1] + tab$fn[1]
      N <- tab$fp[1] + tab$tn[1]
      curves[[length(curves) + 1L]] <- data.frame(model_id = mm$model_id, focal = focal,
                                                  set_index = tab$set_index, cutoff = tab$cutoff,
                                                  cutoff_kind = tab$cutoff_kind, fpr = tab$fp / N,
                                                  tpr = tab$tp / P, stringsAsFactors = FALSE)
    }
    op <- mm$operating_point
    if (!is.null(op) && op$n_events > 0 && op$n_non_events > 0) {
      points[[length(points) + 1L]] <- data.frame(
        model_id = mm$model_id, focal = focal, cutoff = op$cutoff %||% NA_real_, rule_text = op$rule_text,
        fpr = op$fp / op$n_non_events, tpr = op$tp / op$n_events,
        marker = if (focal) "current_operating_point" else "comparison_operating_point",
        stringsAsFactors = FALSE)
    }
    sm <- .ig_cls_frame(mm$summaries, "summaries")
    a <- sm[sm$metric_id == "roc_auc", , drop = FALSE]
    if (nrow(a)) {
      aucs[[length(aucs) + 1L]] <- data.frame(model_id = mm$model_id, roc_auc = a$value[1],
                                              status = a$status[1], reason_code = a$reason_code[1],
                                              stringsAsFactors = FALSE)
    }
  }
  bind <- function(l) if (length(l)) do.call(rbind, l) else NULL
  dg <- x$diagnostics
  n_excl <- dg$n_excluded_rows %||% 0
  list(data = list(curves = bind(curves), operating_points = bind(points), auc = bind(aucs),
                   comparisons = .ig_cls_frame(x$comparisons, "comparisons")),
       notes = c(sprintf("all models are evaluated on the same %s units%s",
                         format(sum(unlist(x$models[[1]]$operating_point[c("n_events", "n_non_events")]))),
                         if (n_excl > 0) sprintf(" (%s row(s) with missing values excluded for every model)", format(n_excl)) else ""),
                 "each model keeps its own rule: the same numerical cutoff is not an equivalent operating point across models",
                 "comparison markers are distinct from the current operating point and from realised varying-cutoff rules",
                 "paired DeLong comparisons are conditional on the fixed score sets"))
}

.ig_cls_pd_sweep <- function(x, m, basis) {
  .ig_cls_pd_need_curve(x, m, basis, "sweep")
  tab <- cutoff_sweep(x, model_id = m$model_id, basis = basis)
  list(data = list(sweep = tab),
       notes = c("one row per distinct decision set (or theoretical display cutoff) at full precision; the selected operating point is always included (column selected)"))
}

.ig_cls_pd_current <- function(x, m, basis, lang) {
  cells <- .ig_cls_pd_need_cells(m, basis)
  mt <- cells$metrics
  if (lang == "de") mt <- .ig_cls_translate_frame(mt, "metrics")
  dict <- metric_dictionary(lang)
  mt$short <- dict$short[match(mt$metric_id, dict$id)]
  mt$formula <- dict$formula[match(mt$metric_id, dict$id)]
  list(data = list(metrics = mt), notes = c(sprintf("count kind: %s", cells$count_kind)))
}
