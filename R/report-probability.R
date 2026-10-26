# Explicitly opted-in aggregates; score ranges and individual scores are omitted.
.ig_report_probability <- function(envelopes) {
  entries <- list()
  rows <- function(tab, keep) {
    if (is.null(tab) || !length(tab[[keep[1L]]])) return(list())
    lapply(seq_along(tab[[keep[1L]]]), function(i) lapply(tab[keep], `[[`, i))
  }
  for (i in seq_along(envelopes)) {
    e <- envelopes[[i]]
    if (e$analysis_kind != "classification" || is.null(e$result)) next
    for (m in e$result$models) {
      if (!identical(m$score_scale, "probability") || !identical(m$direction, "higher_is_positive") ||
          !identical(m$status$code, "ok")) next
      calibration <- rows(m$probability$calibration, c("group", "n_units", "events", "mean_predicted", "observed_rate"))
      net_benefit <- rows(m$probability$net_benefit, c("threshold_probability", "harm_weight", "tp", "fp", "n_total",
        "net_benefit", "net_benefit_alert_all", "net_benefit_alert_none"))
      if (!length(calibration) && !length(net_benefit)) next
      entries[[length(entries)+1L]] <- list(analysis = i, model_id = m$model_id,
        calibration = calibration, net_benefit = net_benefit)
    }
  }
  entries
}

.ig_report_probability_figures <- function(entries, locale) {
  out <- list(); de <- locale == "de"
  for (e in entries) for (panel in c("calibration", "net_benefit")) {
    if (!length(e[[panel]])) next
    cal <- panel == "calibration"
    out[[length(out)+1L]] <- list(analysis = e$analysis, model_id = e$model_id,
      plot_type = paste0("classifier_", panel), basis = "observed", status = "ok",
      axes = if (cal) list(x = if (de) "Mittlere vorhergesagte Wahrscheinlichkeit" else "Mean predicted probability",
        y = if (de) "Beobachtete Ereignisrate" else "Observed event rate") else
        list(x = if (de) "Entscheidungsschwelle" else "Decision threshold probability",
          y = if (de) "Nettonutzen pro Einheit" else "Net benefit per unit"),
      data = e[[panel]], caption = if (cal) {
        if (de) "Explizit freigegebene Gruppenaggregate; kleine Gruppen koennen sensible Informationen offenlegen. Keine Intervalle oder interpolierte Kalibrierungskurve." else
          "Explicitly included grouped aggregates; small groups may disclose sensitive information. No intervals or interpolated calibration curve."
      } else if (de) "Nur explizit angeforderte Wahrscheinlichkeitsschwellen. Modell, alle und niemand als Referenzen; keine Interpolation oder Konfidenzintervalle. Nutzen setzt passende Kosten und Wahrscheinlichkeiten voraus." else
        "Explicitly requested probability thresholds only. Model, alert-all and alert-none references; no interpolation or confidence intervals. Benefit presumes suitable costs and probabilities.")
  }
  out
}

.ig_report_validate_probability <- function(x) {
  bad <- function(message) .ig_abort("invalid_input", paste("Invalid probability report:", message))
  fields <- function(v, names) {
    if (!is.list(v) || anyDuplicated(names(v)) || !setequal(names(v), names)) bad("unexpected aggregate fields")
  }
  scalar <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v)
  near <- function(a,b) scalar(a) && scalar(b) && abs(a-b) <= 1e-10*max(1,abs(a),abs(b))
  entries <- x$probability_summaries %||% list()
  if(length(entries)&&!isTRUE(x$privacy$probability_aggregates_included)) bad("probability aggregates require explicit inclusion")
  if (length(entries) && !identical(x$schema_version, "inferguidR.report_document/0.4.0")) bad("aggregate content requires schema 0.4")
  if (!is.list(entries) || (length(entries) && !is.null(names(entries)))) bad("expected array")
  keys <- character()
  for (e in entries) {
    fields(e, c("analysis", "model_id", "calibration", "net_benefit"))
    if (!scalar(e$analysis) || e$analysis != floor(e$analysis) || e$analysis < 1 || e$analysis > length(x$results) ||
        !.ig_is_string(e$model_id)) bad("invalid model reference")
    key <- paste(e$analysis,e$model_id); if (key %in% keys) bad("duplicate model"); keys <- c(keys,key)
    if (x$specifications[[e$analysis]]$analysis_kind != "classification") bad("classifier required")
    m <- Filter(function(m) identical(m$model_id,e$model_id),x$results[[e$analysis]]$models)
    if (length(m) != 1L || !identical(m[[1]]$status$code,"ok") ||
        !identical(m[[1]]$score_scale,"probability") || !identical(m[[1]]$direction,"higher_is_positive")) bad("unavailable probability model")
    m <- m[[1]]; op <- m$operating_point
    if (is.null(op)) bad("missing operating margins")
    settings <- x$results[[e$analysis]]$settings
    for (kind in c("calibration","net_benefit"))
      if (!is.list(e[[kind]]) || (length(e[[kind]]) && !is.null(names(e[[kind]])))) bad("expected aggregate rows")
    if (!length(e$calibration)) bad("missing calibration aggregates")
    if (length(e$calibration)) {
      if (length(e$calibration) > settings$calibration_bins) bad("too many bins")
      for (j in seq_along(e$calibration)) {
        r <- e$calibration[[j]]
        fields(r,c("group","n_units","events","mean_predicted","observed_rate"))
        if (!all(vapply(r,scalar,TRUE)) || r$group != j || r$n_units <= 0 || r$n_units != floor(r$n_units) ||
            r$events < 0 || r$events > r$n_units || r$events != floor(r$events) ||
            r$mean_predicted < 0 || r$mean_predicted > 1 || !near(r$observed_rate,r$events/r$n_units)) bad("invalid calibration arithmetic")
      }
      total <- sum(vapply(e$calibration, `[[`, 0, "n_units")); events <- sum(vapply(e$calibration, `[[`, 0, "events"))
      if (any(diff(vapply(e$calibration,`[[`,0,"mean_predicted")) < 0)) bad("calibration bins are out of score order")
      if (!near(total,op$n_total) || !near(events,op$n_events)) bad("calibration margins disagree")
      meanp <- sum(vapply(e$calibration,function(r) r$mean_predicted*r$n_units,0))/total
      ids <- unlist(m$summaries$metric_id); j <- match("calibration_in_the_large",ids)
      if (is.na(j) || !near(meanp-events/total,m$summaries$value[[j]])) bad("calibration summary disagrees")
    }
    thresholds <- numeric()
    for (r in e$net_benefit) {
      fields(r,c("threshold_probability","harm_weight","tp","fp","n_total","net_benefit","net_benefit_alert_all","net_benefit_alert_none"))
      if (!all(vapply(r,scalar,TRUE)) || r$threshold_probability <= 0 || r$threshold_probability >= 1 ||
          r$tp < 0 || r$tp > op$n_events || r$fp < 0 || r$fp > op$n_non_events ||
          r$tp != floor(r$tp) || r$fp != floor(r$fp) || !near(r$n_total,op$n_total) ||
          !near(r$harm_weight,r$threshold_probability/(1-r$threshold_probability)) ||
          !near(r$net_benefit,(r$tp-r$fp*r$harm_weight)/r$n_total) ||
          !near(r$net_benefit_alert_all,(op$n_events-op$n_non_events*r$harm_weight)/r$n_total) ||
          r$net_benefit_alert_none != 0) bad("invalid net-benefit arithmetic")
      thresholds <- c(thresholds,r$threshold_probability)
    }
    if (!identical(thresholds,sort(unique(thresholds))) ||
        !identical(thresholds,as.double(unlist(settings$threshold_probability)))) bad("missing, unrequested or duplicated thresholds")
    if (length(e$net_benefit) && any(vapply(c("tp","fp"),function(cell)
        any(diff(vapply(e$net_benefit,`[[`,0,cell)) > 0),TRUE))) bad("alert counts increase with threshold")
  }
  if (identical(x$schema_version,"inferguidR.report_document/0.4.0")&&isTRUE(x$privacy$probability_aggregates_included)) {
    expected_keys <- character()
    for (i in seq_along(x$results)) {
      if (x$specifications[[i]]$analysis_kind != "classification") next
      for (m in x$results[[i]]$models) if (identical(m$status$code,"ok") &&
          identical(m$score_scale,"probability") && identical(m$direction,"higher_is_positive") &&
          identical(x$results[[i]]$input_kind,"scored_data")) expected_keys <- c(expected_keys,paste(i,m$model_id))
    }
    if (!setequal(keys,expected_keys)) bad("probability model aggregates are incomplete")
  }
  actual <- Filter(function(v) v$plot_type %in% c("classifier_calibration","classifier_net_benefit"),x$visualisations)
  expected <- .ig_report_probability_figures(entries,x$locale)
  if (!identical(.ig_canonical_json(actual),.ig_canonical_json(.ig_report_plain(expected)))) bad("figures disagree with aggregates")
  invisible(TRUE)
}

.ig_report_draw_probability <- function(v, locale) {
  de <- locale == "de"
  if (v$plot_type == "classifier_calibration") {
    graphics::plot(c(0,1),c(0,1),type="n",xlab=v$axes$x,ylab=v$axes$y,
      main=if(de) "Kalibrierungsdiagnostik" else "Calibration diagnostic")
    graphics::abline(0,1,col="grey60",lty=2)
    graphics::points(vapply(v$data,`[[`,0,"mean_predicted"),vapply(v$data,`[[`,0,"observed_rate"),pch=16,col="#007C78")
  } else {
    ids <- c("net_benefit","net_benefit_alert_all","net_benefit_alert_none")
    values <- unlist(lapply(v$data,function(r) unlist(r[ids])),use.names=FALSE)
    graphics::plot(range(vapply(v$data,`[[`,0,"threshold_probability")),range(values),type="n",
      xlab=v$axes$x,ylab=v$axes$y,main=if(de) "Entscheidungsnutzen" else "Decision net benefit")
    for (i in seq_along(ids)) graphics::points(vapply(v$data,`[[`,0,"threshold_probability"),
      vapply(v$data,`[[`,0,ids[i]),pch=c(16,2,4)[i],col=c("#007C78","#956000","#697984")[i])
    graphics::legend("bottomleft",legend=if(de) c("Modell","Alle","Niemand") else c("Model","Alert all","Alert none"),
      pch=c(16,2,4),col=c("#007C78","#956000","#697984"),bty="n",cex=.8)
  }
  graphics::mtext(.ig_report_classifier_subtitle(v$model_id,if(de) "Ausgangsdaten" else "Original data"),side=3,line=.6,cex=.8)
}
