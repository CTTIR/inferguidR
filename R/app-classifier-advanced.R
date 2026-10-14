# Input keys follow column identity when the user reorders model columns.
.ig_classifier_model_key <- function(column, field) {
  paste0("model_", digest::digest(enc2utf8(column), algo="sha256", serialize=FALSE), "_", field)
}

.ig_classifier_model_ui <- function(ns, input, columns) {
  shiny::tagList(lapply(columns, function(column) {
    key <- function(field) .ig_classifier_model_key(column, field)
    value <- function(field, default="") shiny::isolate(input[[key(field)]]) %||% default
    shiny::tags$details(shiny::tags$summary(paste("Model / Modell:",column)),
      shiny::textInput(ns(key("name")),"Name",value("name",column)),
      shiny::textInput(ns(key("version")),"Version",value("version")),
      shiny::textInput(ns(key("meaning")),"Score meaning / Score-Bedeutung",value("meaning")),
      shiny::textInput(ns(key("population")),"Population",value("population")),
      shiny::textInput(ns(key("unit")),"Evaluation unit / Auswertungseinheit",value("unit")),
      shiny::selectInput(ns(key("scale")),"Score scale / Score-Skala",
        c("Use shared scale / Gemeinsame Skala"="shared","Raw score / Roher Score"="raw_score",
          "Probability / Wahrscheinlichkeit"="probability"),value("scale","shared")),
      shiny::selectInput(ns(key("direction")),"Direction / Richtung",
        c("Use shared direction / Gemeinsame Richtung"="shared","Higher is positive / Hoeher ist positiv"="higher_is_positive",
          "Lower is positive / Niedriger ist positiv"="lower_is_positive"),value("direction","shared")),
      shiny::selectInput(ns(key("comparator")),"Comparator / Vergleichsoperator",
        c("Inclusive in declared direction / Inklusiv in angegebener Richtung"="auto",">=",">","<=","<"),value("comparator","auto")))
  }))
}

.ig_classifier_model_inputs <- function(input, columns) {
  ids <- paste0("model_",seq_along(columns))
  records <- lapply(seq_along(columns),function(i) {
    value <- function(field,default=NULL) {
      x <- input[[.ig_classifier_model_key(columns[i],field)]]
      if(is.null(x) || !nzchar(trimws(x))) default else x
    }
    direction <- value("direction","shared")
    if(direction=="shared") direction <- input$direction %||% "higher_is_positive"
    scale <- value("scale","shared")
    if(scale=="shared") scale <- input$score_scale %||% "raw_score"
    comparator <- value("comparator","auto")
    if(comparator=="auto") comparator <- if(direction=="higher_is_positive") ">=" else "<="
    list(model=classifier_model(ids[i],name=value("name",columns[i]),version=value("version"),
      score_meaning=value("meaning"),population=value("population"),unit=value("unit")),
      direction=direction,comparator=comparator,scale=scale)
  })
  list(models=lapply(records,`[[`,"model"),
    direction=stats::setNames(vapply(records,`[[`,"","direction"),ids),
    comparator=stats::setNames(vapply(records,`[[`,"","comparator"),ids),
    scale=stats::setNames(vapply(records,`[[`,"","scale"),ids))
}

.ig_classifier_uncertainty_ui <- function(ns) {
  shiny::tags$details(shiny::tags$summary("Uncertainty / Unsicherheit"),
    shiny::textInput(ns("confidence_level"),"Confidence level / Konfidenzniveau","95e-2"),
    shiny::checkboxGroupInput(ns("proportion_methods"),"Proportion intervals / Anteilsintervalle",
      c("Wilson"="wilson","Clopper-Pearson"="clopper_pearson"),selected=c("wilson","clopper_pearson")),
    shiny::selectInput(ns("auc_interval"),"AUC interval / AUC-Intervall",c("DeLong"="delong","None / Keines"="none")),
    shiny::textInput(ns("bootstrap_replicates"),"Bootstrap replicates (0 disables) / Wiederholungen (0 deaktiviert)","0"),
    shiny::textInput(ns("bootstrap_seed"),"Bootstrap seed / Startwert",""),
    shiny::selectInput(ns("bootstrap_strata"),"Resampling / Stichprobenziehung",c("Within class / Innerhalb der Klasse"="class","Unstratified / Unstratifiziert"="none")),
    shiny::selectInput(ns("bootstrap_metrics"),"Bootstrap metrics / Bootstrap-Kennzahlen",
      c(.ig_cls_operating_ids(),.ig_cls_bootstrap_curve_ids),multiple=TRUE,
      selected=c("tpr","fpr","ppv","npv","youden_j","f1","mcc","roc_auc","average_precision")),
    shiny::tags$p("Intervals condition on a fixed cutoff and omit cutoff-selection uncertainty. Declared clusters are resampled as whole clusters. Expected scenario counts have no sampling interval. / Intervalle gelten bei festem Schwellenwert ohne Auswahlunsicherheit. Angegebene Cluster werden als Ganzes gezogen. Erwartete Szenario-Anzahlen haben kein Stichprobenintervall."))
}

.ig_classifier_uncertainty_input <- function(input,lang) {
  num <- function(name,default) .ig_workbench_parse_number(input[[name]] %||% default,lang)
  seed <- input$bootstrap_seed %||% ""
  classification_uncertainty(level=num("confidence_level",if(lang=="de") "0,95" else "0.95"),
    proportion_methods=if(is.null(input$confidence_level)) c("wilson","clopper_pearson") else input$proportion_methods %||% character(),
    auc_method=input$auc_interval %||% "delong",bootstrap_replicates=num("bootstrap_replicates","0"),
    seed=if(nzchar(trimws(seed))) .ig_workbench_parse_number(seed,lang) else NULL,
    bootstrap_strata=input$bootstrap_strata %||% "class",
    bootstrap_metrics=if(is.null(input$confidence_level)) c("tpr","fpr","ppv","npv","youden_j","f1","mcc","roc_auc","average_precision") else input$bootstrap_metrics %||% character())
}

.ig_classifier_counts_rule <- function(spec,input,lang) {
  if(!identical(input$selection %||% "manual","manual"))
    .ig_abort("invalid_input","Fixed counts do not support cutoff selection")
  rule <- .ig_classifier_input_rule(input,lang,spec$operating_rule$score_cutoff)
  # Keep the existing binary-decision provenance; no measured score is invented.
  rule$cutoff_source <- spec$operating_rule$cutoff_source
  update(spec,rule=rule)
}
