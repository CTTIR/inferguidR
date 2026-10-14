.ig_cls_interval_overlay <- function(x, pd, method, lang) {
  if(!.ig_is_string(method) || !method %in% c("wilson","clopper_pearson","bootstrap_percentile"))
    .ig_abort("invalid_input","interval_method must be wilson, clopper_pearson or bootstrap_percentile")
  de <- lang=="de"
  ids <- if(pd$panel=="roc") c("fpr","tpr") else c("tpr","ppv")
  table <- classification_table(x,"intervals",model_id=pd$model_id,lang=lang)
  point <- pd$data$current
  xy <- if(pd$panel=="roc") c("fpr","tpr") else c("recall","precision")
  rows <- lapply(seq_along(ids),function(i) {
    row <- table[table$metric_id==ids[i] & table$method==method,,drop=FALSE]
    status <- "not_computed";reason <- if(de) "Intervallmethode nicht berechnet." else "Interval method was not computed."
    lower <- upper <- level <- NA_real_
    if(nrow(row)==1L) {
      status <- row$status;reason <- row$reason;lower <- row$lower;upper <- row$upper;level <- row$level
    }
    if(pd$basis!="observed") {
      status <- "not_applicable";reason <- if(de) "Szenariowerte sind erwartet, nicht beobachtet." else "Scenario values are expected, not observed."
      lower <- upper <- NA_real_
    }
    px <- py <- NA_real_
    if(!is.null(point) && nrow(point)==1L) {px <- point[[xy[1]]];py <- point[[xy[2]]]}
    if(!is.finite(px) || !is.finite(py)) {
      status <- "not_applicable";reason <- if(de) "Kein darstellbarer aktueller Betriebspunkt." else "No drawable current operating point."
      lower <- upper <- NA_real_
    }
    data.frame(metric_id=ids[i],axis=if(i==1L) "x" else "y",method=method,level=level,
      lower=lower,upper=upper,status=status,reason=reason,
      x=if(i==1L) lower else px,xend=if(i==1L) upper else px,
      y=if(i==2L) lower else py,yend=if(i==2L) upper else py,stringsAsFactors=FALSE)
  })
  data <- do.call(rbind,rows)
  available <- data$status=="ok" & is.finite(data$x) & is.finite(data$xend) & is.finite(data$y) & is.finite(data$yend)
  note <- if(de) "Balken: getrennte marginale Intervalle am festen Betriebspunkt; kein gemeinsamer Konfidenzbereich und kein Kurvenband. Auswahlunsicherheit ist nicht enthalten." else
    "Bars: separate marginal intervals at the fixed operating point; neither a joint confidence region nor a curve band. Cutoff-selection uncertainty is not included."
  if(any(!available)) note <- paste(note,if(de) "Nicht verfuegbare Intervalle stehen mit Grund in der Datentabelle der Abbildung." else "Unavailable intervals retain their reasons in the plot data table.")
  list(data=data,drawable=data[available,,drop=FALSE],caption=note)
}

.ig_classifier_interval_table <- function(result, model_id, basis, lang) {
  d <- classification_table(result,"intervals",model_id=model_id,lang=lang)
  if(!nrow(d)) d <- data.frame(metric_id="all",method="",level=NA_real_,lower=NA_real_,upper=NA_real_,
    status="not_computed",reason=if(lang=="de") "Keine Intervalle fuer dieses Ergebnis berechnet." else
      "No intervals were computed for this result.",conditioning="",note="",stringsAsFactors=FALSE)
  if(basis=="scenario" && nrow(d)) {
    d$lower <- d$upper <- NA_real_;d$status <- "not_applicable"
    d$reason <- if(lang=="de") "Szenariowerte sind erwartet; keine Stichprobenintervalle." else "Scenario values are expected; no sampling intervals."
  }
  d[,intersect(c("metric_id","method","level","lower","upper","status","reason","conditioning","note"),names(d)),drop=FALSE]
}

.ig_classifier_comparison_table <- function(result,lang) {
  d <- classification_table(result,"comparisons",lang=lang)
  if(!nrow(d)) d <- data.frame(inference_status="not_applicable",
    reason=if(lang=="de") "Mindestens zwei Modelle mit gepaarten Score-Daten erforderlich." else
      "At least two models with paired scored data are required.",stringsAsFactors=FALSE)
  d[,intersect(c("model_a","model_b","auc_a","auc_b","difference","lower","upper","level","p_value","method","inference_status","reason"),names(d)),drop=FALSE]
}

.ig_classifier_display_table <- function(d,lang) {
  for(nm in names(d)) if(is.numeric(d[[nm]])) d[[nm]] <- format(d[[nm]],digits=6,
    decimal.mark=if(lang=="de") "," else ".",trim=TRUE)
  labels <- if(lang=="de") c(metric_id="Kennzahl",method="Methode",level="Konfidenzniveau",
    lower="Untergrenze",upper="Obergrenze",status="Status",reason="Grund",conditioning="Bedingung",note="Hinweis",
    model_a="Modell A",model_b="Modell B",auc_a="AUC A",auc_b="AUC B",difference="AUC A - B",
    p_value="p (unkorrigiert)",inference_status="Inferenzstatus") else
    c(metric_id="Metric",method="Method",level="Confidence level",lower="Lower",upper="Upper",
      status="Status",reason="Reason",conditioning="Conditioning",note="Note",model_a="Model A",model_b="Model B",
      auc_a="AUC A",auc_b="AUC B",difference="AUC A - B",p_value="p (unadjusted)",inference_status="Inference status")
  n <- names(d);found <- n %in% names(labels);n[found] <- unname(labels[n[found]]);names(d) <- n
  d
}
