# Selected ROC/PR coordinates from frozen aggregate counts; no score rows or thresholds.
.ig_report_classifier_operating_figures <- function(m, analysis, locale) {
 if(is.null(m$operating_point)) return(list())
 de<-locale=="de";original_basis<-if(m$count_kind=="expected") "theoretical" else "observed"
 coordinate<-function(cells,basis,kind) {
  p<-cells$tp+cells$fn;n<-cells$fp+cells$tn;alerts<-cells$tp+cells$fp
  x<-if(kind=="roc") {if(n>0) cells$fp/n else NULL} else {if(p>0) cells$tp/p else NULL}
  y<-if(kind=="roc") {if(p>0) cells$tp/p else NULL} else {if(alerts>0) cells$tp/alerts else NULL}
  reason<-if(!is.null(x)&&!is.null(y)) NULL else if(kind=="roc") {
   if(de) "ROC-Punkt benoetigt beide Klassen." else "ROC point requires both classes."
  } else if(de) "PR-Punkt benoetigt Ereignisse und mindestens einen Alarm." else "PR point requires events and at least one alert."
  list(basis=basis,x=x,y=y,status=if(is.null(reason)) "ok" else "undefined",reason=reason)
 }
 lapply(c("roc","pr"),function(kind) {
  points<-list(coordinate(m$operating_point,original_basis,kind))
  if(kind=="pr"&&!is.null(m$scenario)&&identical(m$scenario$status,"ok"))
   points[[2L]]<-coordinate(m$scenario$counts,"scenario",kind)
  available<-any(vapply(points,function(z) identical(z$status,"ok"),FALSE))
  list(analysis=analysis,model_id=m$model_id,plot_type=paste0("classifier_",kind,"_point"),basis=original_basis,
    status=if(available) "ok" else "not_applicable",
    reason=if(available) NULL else points[[1L]]$reason,
    axes=list(x=if(kind=="roc") {if(de) "Falsch-positiv-Rate" else "False-positive rate"} else {if(de) "Sensitivitaet (Recall)" else "Recall (sensitivity)"},
      y=if(kind=="roc") {if(de) "Sensitivitaet" else "Sensitivity"} else {if(de) "Positiver Vorhersagewert (Precision)" else "Precision (positive predictive value)"}),
    data=list(points=points,x_limits=c(0,1),y_limits=c(0,1)),
    caption=if(de) paste("Punkte am gewaehlten Schwellenwert aus den gespeicherten Zellen; keine vollstaendige Kurve.",
      if(kind=="roc") "AUC ist eine Kenngroesse der vollstaendigen Kurve und nicht dieses einzelnen Punkts." else "Praevalenzszenarien halten die klassenbedingten Raten fest; Szenariopunkte sind keine neuen Beobachtungen.",
      "Undefinierte Koordinaten werden nicht als null gezeichnet. Unsicherheitsmethoden und Verfuegbarkeit stehen separat im Bericht.") else
      paste("Points at the selected cutoff from retained cells; not a full curve.",
      if(kind=="roc") "AUC describes the full curve, not this individual point." else "Prevalence scenarios hold class-conditional rates fixed; scenario points are not new observations.",
      "Undefined coordinates are not drawn as zero. Uncertainty methods and availability are reported separately."))
 })
}
.ig_report_draw_classifier_operating <- function(v,locale) {
 de<-locale=="de";roc<-v$plot_type=="classifier_roc_point"
 graphics::plot(c(0,1),c(0,1),type="n",xaxs="r",yaxs="r",xlab=v$axes$x,ylab=v$axes$y,
   main=if(roc) {if(de) "ROC: gewaehlter Schwellenwert" else "ROC: selected cutoff"} else {if(de) "Precision-Recall: gewaehlter Schwellenwert" else "Precision-recall: selected cutoff"})
 valid<-Filter(function(z) z$status=="ok",v$data$points)
 for(j in seq_along(valid)) graphics::points(valid[[j]]$x,valid[[j]]$y,
   pch=if(valid[[j]]$basis=="scenario") 17 else 1,cex=if(valid[[j]]$basis=="scenario") 1.1 else 1.7,lwd=1.5,col=if(valid[[j]]$basis=="scenario") "#A65E00" else "#007C78")
 label<-function(basis) switch(basis,observed=if(de) "Beobachtet" else "Observed",
   theoretical=if(de) "Theoretisch" else "Theoretical",scenario=if(de) "Praevalenzszenario" else "Prevalence scenario")
 if(length(valid)) graphics::mtext(paste(vapply(valid,function(z) paste0(label(z$basis),
   if(z$basis=="scenario") {if(de) " (Dreieck)" else " (triangle)"} else {if(de) " (Kreis)" else " (circle)"}),""),collapse=" | "),
   side=1,line=3.8,cex=.8)
 graphics::mtext(.ig_report_classifier_subtitle(v$model_id,if(de) "Fester Schwellenwert" else "Fixed cutoff"),side=3,line=.6,cex=.8)
 invisible(NULL)
}
