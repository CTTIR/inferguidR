# Human-readable summaries use only the immutable portable report values.
.ig_report_summary_lines <- function(x) {
  de <- identical(x$locale,"de")
  tr <- function(en,de_text) if(de) de_text else en
  val <- function(v) {
    if(is.null(v)) return(tr("not available","nicht verf\u00fcgbar"))
    if(is.list(v)) {
      if(identical(v$status,"positive_infinity")) return("+Inf")
      if(identical(v$status,"negative_infinity")) return("-Inf")
      if(!is.null(v$status)) return(as.character(v$status))
      return(paste(vapply(v,val,""),collapse=", "))
    }
    if(is.numeric(v)) return(paste(vapply(v,function(z) format(signif(z,6),digits=6,scientific=NA,trim=TRUE,decimal.mark="."),""),collapse=", "))
    paste(.ig_report_literal(v),collapse=", ")
  }
  row <- function(label,value) paste0(label,": ",val(value))
  lines <- c(.ig_report_literal(x$title),
    tr("RESULT SUMMARY","ERGEBNIS\u00dcBERSICHT"),
    row(tr("Scope","Umfang"),switch(x$report_kind,single_analysis=tr("Single analysis","Einzelanalyse"),model_cutoff=tr("Model and cutoff","Modell und Schwellenwert"),comparison=tr("Comparison","Vergleich"),session=tr("Session","Sitzung"),incomplete_settings=tr("Incomplete settings","Unvollst\u00e4ndige Einstellungen"))),row(tr("Result state","Ergebnisstatus"),if(x$locale=="de") switch(x$completeness,complete="vollstaendig",incomplete="unvollstaendig",x$completeness) else x$completeness),
    row("Snapshot",x$snapshot_id),row(tr("Created (UTC)","Erstellt (UTC)"),x$created_at),
    row(tr("Display timezone","Anzeigezeitzone"),x$display_timezone),
    tr("Summary values use a decimal point and six significant digits. Decisions use stored full precision, retained in the appendix.",
       "Zusammenfassungswerte verwenden einen Dezimalpunkt und sechs signifikante Stellen. Entscheidungen beruhen auf der im Anhang bewahrten vollen Praezision."))
  if(x$report_kind=="comparison") lines<-c(lines,if(length(x$specifications)==1L) tr(
    "Models within one frozen analysis. Paired inference is shown only where retained and available; no multiplicity adjustment is implied.",
    "Modelle innerhalb einer eingefrorenen Analyse. Gepaarte Inferenz erscheint nur, wenn sie vorliegt; keine Mehrfachkorrektur wird unterstellt.") else tr(
    "Comparison of separate frozen analyses; no automatic paired test or test of differences across analyses is inferred.",
    "Vergleich getrennt eingefrorener Analysen; keine Paarung oder Differenztests zwischen Analysen werden abgeleitet."))
  for(i in seq_along(x$specifications)) {
    e<-x$specifications[[i]];s<-e$specification;r<-x$results[[i]]
    lines<-c(lines,"",paste(tr("ANALYSIS","ANALYSE"),i),row(tr("Computation status","Berechnungsstatus"),r$status$code))
    if(!identical(r$status$code,"ok")) lines<-c(lines,row(tr("Reason","Grund"),r$status$reason))
    if(e$analysis_kind=="power") {
      lines<-c(lines,row(tr("Method","Methode"),s$method_id),row(tr("Calculation mode","Berechnungsmodus"),s$mode),
        tr("Chosen inputs:","Gew\u00e4hlte Eingaben:"))
      for(key in names(s$inputs)) lines<-c(lines,row(paste("  ",key),s$inputs[[key]]))
      if(length(r$solved)) {
        lines<-c(lines,tr("Solved quantities:","Berechnete Zielgr\u00f6ssen:"))
        for(key in names(r$solved)) lines<-c(lines,row(paste("  ",key),r$solved[[key]]))
      }
      if(length(r$achieved)) {
        lines<-c(lines,tr("Achieved values:","Erreichte Werte:"))
        for(key in names(r$achieved)) lines<-c(lines,row(paste("  ",key),r$achieved[[key]]))
      }
      if(length(r$group_sizes)) for(key in names(r$group_sizes)) lines<-c(lines,row(paste(tr("Group size","Gruppengr\u00f6sse"),key),r$group_sizes[[key]]))
      lines<-c(lines,tr("Prospective design calculation. Inputs do not have sampling confidence intervals.",
        "Prospektive Planungsrechnung. F\u00fcr Eingaben werden keine Stichproben-Konfidenzintervalle berechnet."))
    } else if(e$analysis_kind=="observed_test") {
      d<-s$inputs
      lines<-c(lines,row(tr("Test convention","Testkonvention"),s$options$variant),
        row(tr("Event definition","Ereignisdefinition"),d$event_label),
        paste(val(d$group_labels[[1]]),": ",val(d$x1)," / ",val(d$n1),tr(" events"," Ereignisse"),sep=""),
        paste(val(d$group_labels[[2]]),": ",val(d$x2)," / ",val(d$n2),tr(" events"," Ereignisse"),sep=""),
        row(tr("Boschloo p-value","Boschloo-p-Wert"),r$p_value),
        row(tr("Numerical upper bound","Numerische Obergrenze"),r$p_upper),
        row(tr("Significance level","Signifikanzniveau"),r$alpha),row(tr("Decision","Entscheidung"),if(is.null(r$decision)) NULL else switch(r$decision,reject=tr("Reject the null hypothesis","Nullhypothese verwerfen"),do_not_reject=tr("Do not reject the null hypothesis","Nullhypothese nicht verwerfen"),unresolved=tr("Unresolved at numerical precision","Bei numerischer Pr\u00e4zision unentschieden"))),
        tr("The numerical bound measures optimisation error, not sampling uncertainty. No sampling confidence interval is computed.",
           "Die numerische Grenze beschreibt den Optimierungsfehler, nicht die Stichprobenunsicherheit. Es wird kein Stichproben-Konfidenzintervall berechnet."),
        tr("Two independent binomial groups; these counts are not classifier confusion cells.",
           "Zwei unabh\u00e4ngige binomiale Gruppen; diese Anzahlen sind keine Klassifikator-Konfusionszellen."))
    } else {
      lines<-c(lines,row(tr("Input type","Eingabetyp"),s$input_kind))
      evaluation <- r$evaluation
      label <- evaluation$label
      if (identical(label,"held_out_rule")) lines<-c(lines,tr(
        "Declared held-out evaluation: independence and absence of training leakage are not established by this declaration.",
        "Deklarierte Validierung auf getrennten Daten: Die Angabe belegt weder Unabhaengigkeit noch den Ausschluss von Informationsuebertragung aus dem Modelltraining."))
      if (identical(label,"predeclared_rule")) lines<-c(lines,tr(
        "The cutoff is declared to precede inspection of these outcomes; the declaration is not independently verified.",
        "Der Schwellenwert wurde laut Angabe vor Sichtung dieser Ergebnisse festgelegt; diese Angabe ist nicht unabhaengig geprueft."))
      if (isTRUE(label %in% c("in_sample_selection","same_data_manual_cutoff"))) lines<-c(lines,tr(
        "Descriptive evaluation on data that may have informed the cutoff; this is not validated operational performance.",
        "Deskriptive Auswertung von Daten, die den Schwellenwert beeinflusst haben koennen; keine validierte Leistung im Einsatz."))
      if (!is.null(evaluation$transfer$overlap_check)) {
        check <- evaluation$transfer$overlap_check
        labels <- c(unit_identity=tr("Unit identifiers","Einheitenkennungen"),cluster_identity=tr("Cluster identifiers","Clusterkennungen"))
        for (key in names(labels)) {
          value <- if(identical(check[[key]],"no_overlap_detected")) tr("no overlap detected in provided IDs","keine Ueberschneidung der angegebenen Kennungen gefunden") else tr("not checked: private identities unavailable","nicht geprueft: private Kennungen nicht verfuegbar")
          lines <- c(lines,row(labels[[key]],value))
        }
      }

      for(m in r$models) {
        lines<-c(lines,"",row(tr("Model","Modell"),m$name),row("Model ID",m$model_id),
          row(tr("Score direction","Score-Richtung"),m$direction),row(tr("Comparator","Vergleichsoperator"),m$comparator),
          row(tr("Evaluation population","Auswertungspopulation"),m$evaluation_split),
          row(tr("Count type","Anzahltyp"),m$count_kind))
        op<-m$operating_point
        if(!is.null(op)) for(key in c("cutoff","tp","fn","fp","tn")) lines<-c(lines,row(toupper(key),op[[key]]))
        for(table in list(m$metrics,m$summaries)) if(length(table$metric_id)) {
          selected<-which(unlist(table$metric_id) %in% c("tpr","tnr","ppv","npv","accuracy","roc_auc","average_precision","pr_area_trapezoid","brier_score"))
          for(j in selected) lines<-c(lines,paste0(toupper(table$metric_id[[j]]),": ",val(table$value[[j]]),
            " [",table$status[[j]],"]",if(!identical(table$status[[j]],"ok")) paste0(" ",val(table$reason[[j]])) else ""))
        }
        if(!is.null(m$scenario)) {
          scenario<-m$scenario
          lines<-c(lines,tr("PREVALENCE SCENARIO (expected counts)","PR\u00c4VALENZSZENARIO (erwartete Anzahlen)"),
            row(tr("Scenario prevalence","Szenariopr\u00e4valenz"),scenario$prevalence),
            row(tr("Scenario status","Szenariostatus"),scenario$status))
          for(key in c("tp","fn","fp","tn")) lines<-c(lines,row(paste(tr("Expected","Erwartet"),toupper(key)),scenario$counts[[key]]))
          table<-scenario$metrics
          for(j in which(unlist(table$metric_id) %in% c("ppv","npv"))) lines<-c(lines,row(paste("Scenario",toupper(table$metric_id[[j]])),table$value[[j]]))
          lines<-c(lines,tr("TPR and FPR are held fixed. Fractional expected counts are not observed units and have no sampling interval.",
            "TPR und FPR bleiben konstant. Bruchteilige erwartete Anzahlen sind keine beobachteten Einheiten und haben kein Stichprobenintervall."))
        }
      }
      lines<-c(lines,tr("AUC describes the full curve. Target-prevalence results are scenarios, not new observed measurements.",
        "AUC beschreibt die gesamte Kurve. Ergebnisse f\u00fcr eine Zielpr\u00e4valenz sind Szenarien und keine neuen Beobachtungen."),
        tr("Uncertainty methods, conditional-cutoff assumptions and unavailable quantities are recorded in the appendix.",
        "Unsicherheitsmethoden, Annahmen zum festen Schwellenwert und nicht verf\u00fcgbare Gr\u00f6ssen stehen im Anhang."))
    }
  }
  unavailable <- Filter(function(v) identical(v$status,"unavailable"),x$visualisations)
  if(length(unavailable)) {
    lines <- c(lines,"",tr("UNAVAILABLE FIGURES","NICHT VERFUEGBARE ABBILDUNGEN"))
    for(v in unavailable) lines <- c(lines,paste0(v$plot_type,": ",v$reason %||% v$caption))
  }
  if(length(x$sweeps)) {
    lines<-c(lines,"",tr("DESIGN SWEEPS","ENTWURFSVARIANTEN"))
    for(k in seq_along(x$sweeps)) {
      sw<-x$sweeps[[k]];ok<-vapply(sw$rows,function(row) identical(row$status,"ok"),TRUE)
      lines<-c(lines,paste0(tr("Sweep ","Vergleich "),k,": ",sw$free,"; ",sum(ok),tr(" successful, "," erfolgreich, "),sum(!ok),tr(" failed rows."," fehlgeschlagene Zeilen.")),
        tr("Evaluated points only; all requested values and failure reasons are retained in the appendix.","Nur berechnete Punkte; alle angeforderten Werte und Fehlergruende stehen im Anhang."))
    }
  }
  if (identical(x$resources$status, "attached")) {
    lines <- c(lines, "", tr("TEACHING RESOURCES", "LEHRMATERIALIEN"), x$resources$evidence_note)
    for (resource in x$resources$entries) {
      lines <- c(lines, paste0(tr("Analysis ", "Analyse "), resource$analysis, ": ", resource$title),
        resource$scope, resource$access,
        vapply(resource$references, function(reference) paste0(reference$title, ": ", reference$url), ""))
    }
  }
  c(lines,"",tr("INTERPRETATION AND LIMITS","INTERPRETATION UND GRENZEN"),
    tr("Check the chosen assumptions and evaluation population. A numerical result alone does not establish scientific validity.",
       "Pr\u00fcfen Sie die gew\u00e4hlten Annahmen und Auswertungspopulation. Ein Zahlenwert allein belegt keine wissenschaftliche G\u00fcltigkeit."),
    tr("This development package has not completed numerical release certification. See the validation record in the appendix.",
       "Die numerische Freigabezertifizierung dieser Entwicklungsversion ist noch nicht abgeschlossen. Siehe Validierungsnachweis im Anhang."),
    tr("Titles, model metadata and selected cutoffs are included. Review these before sharing.",
       "Titel, Modellmetadaten und ausgew\u00e4hlte Schwellenwerte sind enthalten. Pr\u00fcfen Sie diese vor der Weitergabe."))
}

.ig_report_wrap_pdf_lines <- function(lines,width_in=7.1,size=11,family="sans") {
  gp <- grid::gpar(fontsize=size,fontfamily=family)
  fits <- function(value) grid::convertWidth(grid::grobWidth(grid::textGrob(value,gp=gp)),"inches",valueOnly=TRUE)<=width_in
  unlist(lapply(lines,function(line) {
    text <- .ig_report_literal(line)
    if(!nzchar(text)) return("")
    out <- character()
    while(nzchar(text)) {
      if(fits(text)) {out<-c(out,text);break}
      low<-1L;high<-nchar(text)
      while(low<high) {
        middle<-ceiling((low+high)/2)
        if(fits(substr(text,1L,middle))) low<-middle else high<-middle-1L
      }
      piece<-substr(text,1L,low);spaces<-gregexpr(" ",piece,fixed=TRUE)[[1]]
      cut<-if(any(spaces>1L)) max(spaces) else low
      out<-c(out,trimws(substr(text,1L,cut),which="right"))
      text<-sub("^ +","",substring(text,cut+1L))
    }
    out
  }),use.names=FALSE)
}

.ig_report_plot_label <- function(text,width_in,size,max_lines=3L) {
  lines<-.ig_report_wrap_pdf_lines(text,width_in=width_in,size=size,family="sans")
  if(length(lines)>max_lines) lines<-c(utils::head(lines,max_lines-1L),"[...]")
  paste(lines,collapse="\n")
}
