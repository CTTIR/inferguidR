.ig_classifier_probability_ui <- function(ns) {
  shiny::tags$details(shiny::tags$summary("Probability diagnostics / Wahrscheinlichkeitsdiagnostik"),
    shiny::textInput(ns("calibration_bins"),"Calibration groups / Kalibrierungsgruppen","10"),
    shiny::textInput(ns("threshold_probability"),"Decision probabilities, separated by ; / Entscheidungsschwellen, getrennt mit ;",""),
    shiny::tags$p("Leave blank to omit net benefit. These probabilities are distinct from a raw score cutoff. / Leer lassen, um Nettonutzen auszulassen. Diese Wahrscheinlichkeiten sind von einem Rohscore-Schwellenwert verschieden."),
    shiny::selectInput(ns("partial_auc_focus"),"Partial ROC area / Partielle ROC-Flaeche",
      c("None / Keine"="none","False positive rate / Falsch-Positiv-Rate"="fpr","True positive rate / Richtig-Positiv-Rate"="tpr")),
    shiny::textInput(ns("partial_auc_low"),"Lower limit / Untere Grenze","0"),
    shiny::textInput(ns("partial_auc_high"),"Upper limit / Obere Grenze","1"))
}

.ig_classifier_probability_thresholds <- function(text, lang) {
  if(is.null(text) || !nzchar(trimws(text))) return(NULL)
  if(grepl("(^|;)\\s*(;|$)",text)) .ig_abort("invalid_input","Decision probabilities must be separated by semicolons without empty entries.")
  values<-vapply(strsplit(text,";",fixed=TRUE)[[1]],.ig_workbench_parse_number,0,locale=lang)
  .ig_cls_check_thresholds(values)
}

.ig_classifier_evaluation_settings <- function(input, lang) {
  focus<-input$partial_auc_focus %||% "none"
  partial<-if(identical(focus,"none")) NULL else list(focus=focus,range=c(
    .ig_workbench_parse_number(input$partial_auc_low %||% "0",lang),
    .ig_workbench_parse_number(input$partial_auc_high %||% "1",lang)))
  list(partial_auc=.ig_cls_check_partial(partial),
    threshold_probability=.ig_classifier_probability_thresholds(input$threshold_probability,lang),
    calibration_bins=.ig_cls_check_bins(.ig_workbench_parse_number(input$calibration_bins %||% "10",lang)),
    include_curves=TRUE)
}

.ig_classifier_evaluation_current <- function(spec, result, settings) {
  !is.null(spec) && !is.null(result) && identical(spec$spec_hash,result$spec_hash) &&
    identical(.ig_canonical_json(settings), .ig_canonical_json(result$settings))
}
