# Retained paired AUC summaries only; no row-level scores or new inference.
.ig_report_classifier_comparison_figures <- function(r, analysis, locale) {
 x<-r$comparisons
 if(is.null(x) || !length(x$model_a)) return(list())
 de<-locale=="de";n<-length(x$model_a)
 fields<-c("model_a","model_b","difference","lower","upper","level","inference_status","reason_code")
 rows<-lapply(seq_len(n),function(j) {
  out<-lapply(x[fields],function(v) {
   value<-if(is.list(v)) v[[j]] else v[j]
   if(is.null(value)||is.na(value)) NULL else unname(value)
  })
  c(list(pair_index=j),out)
 })
 bounds<-unlist(lapply(rows,function(r) c(r$difference,r$lower,r$upper)),use.names=FALSE)
 bounds<-bounds[is.finite(bounds)]
 span<-if(length(bounds)) range(c(0,bounds)) else c(-1,1)
 if(diff(span)<.02) span<-span+c(-.01,.01)
 span<-span+c(-1,1)*diff(span)*.08
 groups<-split(seq_len(n),ceiling(seq_len(n)/8L))
 lapply(seq_along(groups),function(j) list(analysis=analysis,
  plot_type="classifier_comparison",basis="observed",status="ok",
  axes=list(x=if(de) "Differenz der ROC-AUC (A minus B)" else "ROC AUC difference (A minus B)",
   y=if(de) "Modellpaar" else "Model pair"),
  data=list(panel_index=j,panel_count=length(groups),x_limits=span,pairs=rows[groups[[j]]]),
  caption=if(de) paste("Gepaarte ROC-AUC (A minus B), unabhaengig vom Schwellenwert.",
   "Marginale Normalintervalle nach DeLong; keine Mehrfachkorrektur oder Begrenzung auf [-1,1].",
   "Nicht verfuegbare Intervalle bleiben offen. Bezeichnungen und Gruende: Vergleichstabelle.") else
   paste("Paired ROC AUC (A minus B), independent of the selected cutoff.",
    "Marginal normal DeLong intervals; no multiplicity adjustment or clipping to [-1,1].",
    "Unavailable intervals stay absent. Full labels and reasons: comparison table.")))
}

.ig_report_draw_classifier_comparison <- function(v,locale) {
 de<-locale=="de";rows<-v$data$pairs;n<-length(rows)
 span<-unlist(v$data$x_limits,use.names=FALSE)
 graphics::par(mar=c(9,12,5,2))
 graphics::plot(span,c(.5,n+.5),type="n",yaxt="n",xlab=v$axes$x,ylab="",
  main=if(de) "Gepaarte ROC-AUC-Vergleiche" else "Paired ROC AUC comparisons")
 graphics::abline(v=0,lty=2,col="#777777")
 cap<-.055*diff(graphics::par("usr")[3:4])/graphics::par("pin")[2]
 for(i in seq_along(rows)) {
  r<-rows[[i]];y<-n-i+1
  label<-paste0("[",r$pair_index,"] ",r$model_a," - ",r$model_b)
  label<-.ig_report_plot_label(label,width_in=1.8,size=9)
  graphics::mtext(label,side=2,at=y,las=1,line=.3,cex=.75)
  if(!is.null(r$difference)) graphics::points(r$difference,y,pch=16,col="#007C78",cex=1.1)
  if(identical(r$inference_status,"ok")) {
   graphics::segments(r$lower,y,r$upper,y,col="#007C78",lwd=2)
   graphics::segments(c(r$lower,r$upper),y-cap,c(r$lower,r$upper),y+cap,col="#007C78")
  } else graphics::text(span[2],y+2.5*cap,labels=if(de) "Intervall n.v." else "CI unavailable",adj=c(1,0),cex=.7,col="#555555")
 }
 levels<-unique(vapply(rows,function(r) r$level,0))
 graphics::mtext(paste(paste0(format(100*levels,trim=TRUE),"%"),
  if(de) "marginale Konfidenzintervalle" else "marginal confidence intervals",
  paste0(" | ",v$data$panel_index,"/",v$data$panel_count)),side=3,line=.6,cex=.8)
 invisible(NULL)
}
