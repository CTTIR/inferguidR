#' Draw a classification panel
#'
#' Static views of [classification_plot_data()] for classification results. ROC
#' paths retain the ordering of distinct decision sets, including tied scores.
#' Precision-recall points are not interpolated: a straight segment need not
#' correspond to any attainable cutoff. The selected operating point is a ring.
#' The source plot data are retained in the `inferguidR_plot_data` attribute.
#'
#' @param x A `classification_result` from [evaluate_cutoff()].
#' @param panel ROC, precision-recall, mosaic, score distribution, prevalence
#'   curve, or model comparison (see the supported values below).
#' @param model_id Model identifier; defaults to the first evaluated model.
#' @param basis `"observed"` or `"scenario"`; see [classification_plot_data()].
#' @param mosaic_mode Joint population shares or within-class conditional rates.
#' @param bins Number of display bins for scored-data distributions; no metric
#'   or cutoff is computed from these bins.
#' @param metric_family For metric-versus-cutoff views: `within_truth`,
#'   `within_alerts`, `summary`, `burden`, or `likelihood_ratios`. Likelihood
#'   ratios use a logarithmic axis; other families keep their declared domains.
#' @param metric_x,metric_y Operating-point metric identifiers for a metric-pair view.
#' @param log_x,log_y Logarithmic axes for the metric-pair view. Non-positive,
#'   undefined and infinite coordinates are counted and omitted, never clamped.
#' @param metric_focus Operating-point metric whose direct dictionary dependencies
#'   are shown in a metric-network view. This is a calculation graph, not a causal graph.
#' @param interval_method Optional marginal interval bars at the current ROC or
#'   PR point: wilson, clopper_pearson or bootstrap_percentile. The method must
#'   have been computed in the result. These are separate fixed-cutoff intervals,
#'   not a joint confidence region or a curve band. Scenario intervals are not drawn.
#' @param ppv_levels,f1_levels Reference levels strictly between zero and one for ROC iso-PPV and PR iso-F1 curves.
#' @param realised_rules Optional independently supplied counts for varying-cutoff rules; see [classification_plot_data()]. Drawn as separate triangles, never as fixed-threshold curve points.
#' @param lang Label language, `"en"` or `"de"`.
#' @return A `ggplot` object. Requires the suggested package 'ggplot2'.
#'   Fixed counts cannot provide ROC, PR or score-distribution curves. Score
#'   distributions and model comparisons always use the observed/model basis.
#'   Prevalence curves hold the current class-conditional rates fixed.
#' @concept classification
#' @seealso
#' * [classification_spec()]
#' * [evaluate_cutoff()]
#' @export
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   r <- evaluate_cutoff(classification_examples("scored_synthetic")$spec)
#'   plot_classification(r, "roc", "model_b")
#' }
plot_classification <- function(x, panel = c("roc", "pr", "mosaic", "score_distribution",
                                                            "prevalence_curve", "model_comparison", "metric_vs_cutoff", "metric_pair", "metric_network", "calibration", "net_benefit"), model_id = NULL,
                                basis = c("observed", "scenario"), lang = c("en", "de"),
                                mosaic_mode = c("joint", "within_class"), bins = 40,
                                metric_family = "within_truth", metric_x = "fpr", metric_y = "tpr",
                                log_x = FALSE, log_y = FALSE, metric_focus = "ppv", interval_method = NULL,
                                ppv_levels = c(.25,.5,.75), f1_levels = c(.2,.4,.6,.8),
                                realised_rules = NULL) {
  panel <- match.arg(panel)
  if(!is.null(interval_method) && !panel %in% c("roc","pr"))
    .ig_abort("invalid_input","interval bars are available only for ROC and PR panels")
  if(!is.null(realised_rules) && !panel %in% c("roc","pr","metric_pair"))
    .ig_abort("invalid_input","realised rules are available only for ROC, PR and metric-pair panels")
  basis <- match.arg(basis)
  lang <- match.arg(lang)
  mosaic_mode <- match.arg(mosaic_mode)
  if (panel %in% c("calibration", "net_benefit")) return(.ig_plot_classification_probability(
    x, panel, model_id, basis, lang))
  if (panel == "metric_network") return(.ig_plot_classification_network(
    x, model_id, basis, lang, metric_focus))
  if (panel %in% c("metric_vs_cutoff", "metric_pair")) return(.ig_plot_classification_metrics(
    x, panel, model_id, basis, lang, metric_family, metric_x, metric_y, log_x, log_y, realised_rules))
  if (!panel %in% c("roc", "pr")) return(.ig_plot_classification_additional(
    x, panel, model_id, basis, lang, mosaic_mode, bins))
  pd <- classification_plot_data(x, panel, model_id, basis, lang, ppv_levels=ppv_levels, f1_levels=f1_levels, realised_rules=realised_rules)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    .ig_abort("unsupported", "plot_classification() requires the suggested package 'ggplot2'.")
  fpr <- tpr <- recall <- precision <- reference <- NULL
  de <- lang == "de"
  if (panel == "roc") {
    mapping <- ggplot2::aes(x = fpr, y = tpr)
    p <- ggplot2::ggplot(pd$data$points, mapping) +
      ggplot2::geom_path(colour = "#007C78", linewidth = .8) +
      ggplot2::geom_point(colour = "#007C78", size = 1.3) +
      ggplot2::geom_path(data = pd$data$diagonal, colour = "#697984", linetype = "dashed") +
      ggplot2::labs(x = if (de) "Falsch-Positiv-Rate" else "False positive rate",
                    y = if (de) "Richtig-Positiv-Rate" else "True positive rate",
                    caption = if (de) "Punkte: Entscheidungsregeln; Ring: gewaehlter Schwellenwert. Bei festen klassenbedingten Score-Verteilungen bleibt die ROC-Kurve bei geaenderter Praevalenz gleich." else
                      "Points: decision rules; ring: selected cutoff. The ROC curve is unchanged by prevalence when class-conditional score distributions stay fixed.")
  } else {
    mapping <- ggplot2::aes(x = recall, y = precision)
    prevalence <- formatC(pd$data$baseline$precision, digits = 4, format = "g",
                          decimal.mark = if (de) "," else ".")
    p <- ggplot2::ggplot(pd$data$points, mapping) +
      ggplot2::geom_point(colour = "#007C78", size = 1.7) +
      ggplot2::geom_hline(yintercept = pd$data$baseline$precision, colour = "#697984", linetype = "dashed") +
      ggplot2::labs(x = if (de) "Sensitivitaet (Recall)" else "Recall (sensitivity)",
                    y = if (de) "Positiver Vorhersagewert (Precision)" else "Precision (positive predictive value)",
                    caption = paste(if (de) "Keine lineare Interpolation zwischen Entscheidungsregeln. Ring: gewaehlter Schwellenwert.\nGestrichelt: Praevalenz" else
                      "No linear interpolation between decision rules. Ring: selected cutoff.\nDashed: prevalence", prevalence,
                      paste0("(", pd$data$baseline$prevalence_kind, ").")))
  }
  iso <- if(panel=="roc") pd$data$iso_ppv else pd$data$iso_f1
  if(!is.null(iso) && nrow(iso)) {
    level <- if(panel=="roc") iso$ppv_level else iso$f1_level
    iso$reference <- factor(format(level,trim=TRUE),levels=unique(format(level,trim=TRUE)))
    iso_mapping <- if(panel=="roc") ggplot2::aes(fpr,tpr,group=reference,linetype=reference) else
      ggplot2::aes(recall,precision,group=reference,linetype=reference)
    p <- p + ggplot2::geom_path(data=iso,mapping=iso_mapping,inherit.aes=FALSE,
      colour="#8b6b31",linewidth=.45) +
      ggplot2::scale_linetype_manual(values=stats::setNames(rep(c("dashed","dotted","dotdash","longdash","twodash","solid"),length.out=nlevels(iso$reference)),levels(iso$reference))) +
      ggplot2::labs(linetype=if(panel=="roc") "PPV" else "F1")
    note <- if(panel=="roc") paste(if(de) "Referenzlinien fuer PPV bei Praevalenz" else "PPV reference lines at prevalence",
      format(unique(iso$prevalence),digits=4),paste0(" (",unique(iso$prevalence_kind),")."),
      if(de) "Am Ursprung ist PPV nicht definiert." else "PPV is undefined at the origin.") else
      if(de) "F1-Referenzlinien verbinden keine beobachteten Entscheidungsregeln." else "F1 reference lines do not connect observed decision rules."
    p <- p + ggplot2::labs(caption=paste(p$labels$caption,note))
  }
  rr <- pd$data$realised_rules
  if(!is.null(rr) && nrow(rr)) {
    axes <- if(panel=="roc") c("fpr","tpr") else c("recall","precision")
    drawable <- is.finite(rr[[axes[1]]]) & is.finite(rr[[axes[2]]])
    if(any(drawable)) p <- p + ggplot2::geom_point(data=rr[drawable,,drop=FALSE],mapping=mapping,
      shape=17,colour="#8b3a62",size=3,inherit.aes=FALSE)
    p <- p + ggplot2::labs(caption=paste(p$labels$caption,
      if(de) "Dreiecke: separate realisierte Regeln; keine Punkte der festen Schwellenwertkurve." else
        "Triangles: separate realised rules, not fixed-threshold curve points.",
      if(de) paste(sum(!drawable),"Regeln ohne definierte Koordinaten.") else paste(sum(!drawable),"rules have undefined coordinates.")))
  } else if(!is.null(realised_rules) && basis=="scenario") {
    p <- p + ggplot2::labs(caption=paste(p$labels$caption,if(de)
      "Beobachtete realisierte Regeln werden auf der Szenariobasis nicht gezeigt." else
      "Observed realised rules are not shown on the scenario basis."))
  }
  if (!is.null(pd$data$current) && nrow(pd$data$current))
    p <- p + ggplot2::geom_point(data = pd$data$current, mapping = mapping,
                                shape = 21, fill = NA, colour = "#152330", size = 4, stroke = 1.1)
  p <- p + ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, .2)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, .2)) +
    ggplot2::labs(title = paste(if (panel == "roc") "ROC" else "Precision-recall", pd$model_id, sep = ": "),
                  subtitle = paste(if (de) "Bezugsbasis:" else "Basis:", pd$basis)) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, size = 10),
                    panel.grid.minor = ggplot2::element_blank())
  if(!is.null(interval_method)) {
    overlay <- .ig_cls_interval_overlay(x,pd,interval_method,lang)
    x <- xend <- y <- yend <- NULL
    if(nrow(overlay$drawable)) p <- p + ggplot2::geom_segment(data=overlay$drawable,
      mapping=ggplot2::aes(x=x,xend=xend,y=y,yend=yend),inherit.aes=FALSE,
      colour="#956000",linewidth=.8)
    p <- p + ggplot2::labs(caption=paste(strwrap(paste(p$labels$caption,overlay$caption),width=85),collapse="\n"))
    pd$data$intervals <- overlay$data
  }
  p <- p + ggplot2::labs(caption=paste(strwrap(p$labels$caption,width=85),collapse="\n")) +
    ggplot2::theme(legend.position="bottom",plot.caption.position="plot")
  attr(p, "inferguidR_plot_data") <- pd
  p
}
