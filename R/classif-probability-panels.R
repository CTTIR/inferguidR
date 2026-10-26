# Probability diagnostics use only the frozen evaluation's retained tables.
.ig_cls_pd_probability <- function(x, m, panel, basis, lang) {
  if (basis != "observed") .ig_abort("unsupported",
    if (lang == "de") "Kalibrierung und Nettonutzen beziehen sich auf beobachtete Daten, nicht auf ein Praevalenzszenario." else
      "Calibration and net benefit describe observed data, not a prevalence scenario.")
  view <- m$views[[panel]]
  if (is.null(view) || !isTRUE(view$available)) .ig_abort("unsupported",
    .ig_cls_reason_text(view$reason_code %||% "not_probability_scale", lang))
  table <- m$probability[[panel]]
  if (is.null(table)) .ig_abort("unsupported", "This result does not retain the probability diagnostic table.")
  table <- .ig_cls_frame(table, panel)
  if (!nrow(table)) .ig_abort("unsupported", "The probability diagnostic table is empty.")
  summaries <- .ig_cls_frame(m$summaries, "summaries")
  summaries <- summaries[summaries$metric_id %in% c("brier_score","calibration_in_the_large","calibration_intercept","calibration_slope"),,drop=FALSE]
  notes <- if (panel == "calibration") {
    if (lang == "de") "Gruppierte beobachtete Ereignisraten; Gruppierung nur zur Darstellung. Keine Konfidenzintervalle. Diskriminierung belegt keine Kalibrierung." else
      "Grouped observed event rates; grouping is for display only. No confidence intervals. Discrimination does not establish calibration."
  } else if (lang == "de") "Deklarierte Entscheidungsschwellen p; relatives Schadensgewicht p/(1-p). Werte nur an berechneten Schwellen, keine Interpolation. Gueltigkeit haengt von Kalibrierung und Bewertungspopulation ab." else
    "Declared decision probabilities p; relative harm weight p/(1-p). Values only at evaluated thresholds, without interpolation. Validity depends on calibration and the evaluation population."
  list(data=list(points=table,summaries=summaries),notes=notes,
    estimate_kind=m$count_kind,evaluation=x$evaluation)
}

.ig_plot_classification_probability <- function(x, panel, model_id, basis, lang) {
  pd <- classification_plot_data(x,panel,model_id,basis,lang)
  if(!requireNamespace("ggplot2",quietly=TRUE)) .ig_abort("unsupported", "plot_classification() requires the suggested package 'ggplot2'.")
  mean_predicted <- observed_rate <- n_units <- threshold_probability <- value <- rule <- NULL
  de <- lang=="de";d<-pd$data$points
  if(panel=="calibration") {
    p<-ggplot2::ggplot(d,ggplot2::aes(mean_predicted,observed_rate))+
      ggplot2::geom_abline(slope=1,intercept=0,linetype="dashed",colour="#697984")+
      ggplot2::geom_point(ggplot2::aes(size=n_units),colour="#007C78",alpha=.85)+
      ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
      ggplot2::labs(x=if(de) "Mittlere vorhergesagte Wahrscheinlichkeit" else "Mean predicted probability",
        y=if(de) "Beobachtete Ereignisrate" else "Observed event rate",
        size=if(de) "Gruppengewicht" else "Group weight",
        title=if(de) "Kalibrierungsdiagnostik" else "Calibration diagnostic")
  } else {
    ids<-c("net_benefit","net_benefit_alert_all","net_benefit_alert_none")
    points<-do.call(rbind,lapply(ids,function(id)data.frame(threshold_probability=d$threshold_probability,value=d[[id]],rule=id)))
    labels<-if(de)c("Modellregel","Alle alarmieren","Niemanden alarmieren") else c("Model rule","Alert all","Alert none")
    p<-ggplot2::ggplot(points,ggplot2::aes(threshold_probability,value,colour=rule,shape=rule))+
      ggplot2::geom_point(size=2.5)+
      ggplot2::scale_colour_manual(values=stats::setNames(c("#007C78","#956000","#697984"),ids),breaks=ids,labels=labels)+
      ggplot2::scale_shape_manual(values=stats::setNames(c(16,2,4),ids),breaks=ids,labels=labels)+
      ggplot2::labs(x=if(de) "Entscheidungsschwelle (Wahrscheinlichkeit)" else "Decision threshold probability",
        y=if(de) "Nettonutzen pro Auswertungseinheit" else "Net benefit per evaluation unit",colour=NULL,shape=NULL,
        title=if(de) "Entscheidungsnutzen" else "Decision net benefit")
  }
  caption<-paste(strwrap(pd$notes,width=100),collapse="\n")
  p<-p+ggplot2::labs(subtitle=pd$model_id,caption=caption)+ggplot2::theme_minimal(base_size=12)+
    ggplot2::theme(legend.position="bottom",plot.caption=ggplot2::element_text(hjust=0))
  attr(p,"inferguidR_plot_data")<-pd;p
}
