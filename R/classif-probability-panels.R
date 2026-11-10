
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
