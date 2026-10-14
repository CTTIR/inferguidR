.ig_classifier_custom_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::conditionalPanel(sprintf("input['%s'] == 'uploaded'", ns("example")),
      shiny::fileInput(ns("csv"), "CSV data / CSV-Daten", accept = ".csv"),
      shiny::selectInput(ns("separator"), "Delimiter / Trennzeichen", c("Comma / Komma" = ",", "Semicolon / Semikolon" = ";", "Tab" = "\t")),
      shiny::selectInput(ns("decimal"), "Decimal mark / Dezimalzeichen", c("." = ".", "," = ",")),
      shiny::uiOutput(ns("columns")),
      shiny::textInput(ns("positive"), "Positive outcome label / Positives Ergebnis", "1"),
      shiny::selectInput(ns("direction"), "Score direction / Score-Richtung", c("Higher is positive / Hoeher ist positiv" = "higher_is_positive", "Lower is positive / Niedriger ist positiv" = "lower_is_positive")),
      shiny::selectInput(ns("score_scale"), "Shared score scale / Gemeinsame Score-Skala", c("Raw score / Roher Score" = "raw_score", "Probability / Wahrscheinlichkeit" = "probability")),
      shiny::uiOutput(ns("model_settings")),
      shiny::selectInput(ns("sampling"), "Sampling design / Stichprobendesign", c("Unknown / Unbekannt"="unknown", "Cross-sectional / Querschnitt"="cross_sectional", "Case-control / Fall-Kontroll"="case_control")),
      shiny::selectInput(ns("missing"), "Missing values / Fehlende Werte", c("Reject / Ablehnen" = "error", "Exclude and count / Ausschliessen und zaehlen" = "exclude")),
      shiny::selectInput(ns("evaluation_split"), "Evaluation population / Auswertungsdaten", c("Same data / Gleiche Daten" = "same_data", "Held-out test / Getrennte Testdaten" = "held_out_test", "External / Extern" = "external")),
      shiny::tags$p("Scores are evaluated, not trained. Upload only data permitted for this app session. / Scores werden ausgewertet, nicht trainiert. Nur fuer diese Sitzung freigegebene Daten hochladen.")),
    shiny::conditionalPanel(sprintf("input['%s'] == 'custom_counts'", ns("example")),
      shiny::textInput(ns("tp"), "True positives / Richtig positiv", "40"),
      shiny::textInput(ns("fn"), "False negatives / Falsch negativ", "10"),
      shiny::textInput(ns("fp"), "False positives / Falsch positiv", "20"),
      shiny::textInput(ns("tn"), "True negatives / Richtig negativ", "130")),
    shiny::conditionalPanel(sprintf("input['%s'] == 'custom_theory'", ns("example")),
      shiny::textInput(ns("auc"), "AUC", ""),
      shiny::textInput(ns("sigma"), "SD ratio / SD-Verhaeltnis", "1"),
      shiny::textInput(ns("theory_n"), "Expected total N / Erwartetes Gesamt-N", "1000"),
      shiny::textInput(ns("theory_prevalence"), "Assumed prevalence / Angenommene Praevalenz", "")),
    shiny::tags$details(shiny::tags$summary("Selection and constraints / Auswahl und Grenzen"),
      shiny::selectInput(ns("selection"), "Selection rule / Auswahlregel", c("Manual cutoff / Manueller Schwellenwert" = "manual",
        "Maximum Youden" = "max_youden", "Maximum F1" = "max_f1", "Maximum MCC" = "max_mcc", "Constraints / Grenzen" = "constraint",
        "Capacity / Kapazitaet" = "capacity", "Declared costs / Angegebene Kosten" = "cost_utility")),
      shiny::selectInput(ns("selection_data"), "Cutoff origin / Herkunft des Schwellenwerts", c("Same data / Gleiche Daten" = "same_data", "Held-out tuning / Separate Abstimmungsdaten" = "held_out_tuning", "Predeclared manual rule / Vorab festgelegte manuelle Regel" = "predeclared_rule")),
      shiny::textInput(ns("tpr_min"), "Minimum sensitivity (optional) / Minimale Sensitivitaet (optional)", ""),
      shiny::textInput(ns("tnr_min"), "Minimum specificity / Minimale Spezifitaet", ""),
      shiny::textInput(ns("npv_min"), "Minimum NPV / Minimaler NPV", ""),
      shiny::textInput(ns("fnr_max"), "Maximum false-negative rate / Maximale Falsch-Negativ-Rate", ""),
      shiny::textInput(ns("fdr_max"), "Maximum false-discovery rate / Maximale Falschentdeckungsrate", ""),
      shiny::textInput(ns("fpr_max"), "Maximum false-positive rate (optional) / Maximale Falsch-Positiv-Rate (optional)", ""),
      shiny::textInput(ns("ppv_min"), "Minimum PPV (optional) / Minimaler PPV (optional)", ""),
      shiny::textInput(ns("alert_rate_max"), "Maximum alert rate (optional) / Maximale Alarmrate (optional)", ""),
      shiny::selectInput(ns("optimise"), "Optimise within constraints / Optimieren innerhalb der Grenzen", c("Sensitivity"="tpr","F1"="f1","MCC"="mcc","PPV"="ppv")),
      shiny::selectInput(ns("tie_break"), "Tied objectives / Gleiche Zielwerte", c("Fewest alerts / Wenigste Alarme"="fewest_alerts", "Most alerts / Meiste Alarme"="most_alerts")),
      shiny::selectInput(ns("capacity_kind"), "Capacity unit / Kapazitaetseinheit", c("Number of alerts / Anzahl Alarme"="max_alerts", "Alert fraction / Alarmanteil"="max_alert_rate")),
      shiny::textInput(ns("max_alert_rate"), "Capacity: maximum alert fraction / Maximaler Alarmanteil", ""),
      shiny::textInput(ns("max_alerts"), "Capacity: maximum alerts / Kapazitaet: maximale Alarme", "20"),
      shiny::textInput(ns("cost_fn"), "Cost per false negative / Kosten je falsch negativem Ergebnis", "1"),
      shiny::textInput(ns("cost_fp"), "Cost per false positive / Kosten je falsch positivem Ergebnis", "1"),
      shiny::textInput(ns("cost_tp"), "Cost per true positive (optional) / Kosten je richtig positivem Ergebnis (optional)", ""),
      shiny::textInput(ns("cost_tn"), "Cost per true negative (optional) / Kosten je richtig negativem Ergebnis (optional)", ""),
      shiny::tags$p("Constraints and costs are your decisions. Selection on evaluation data is optimistic. / Grenzen und Kosten sind Ihre Festlegungen. Auswahl an Auswertungsdaten ist optimistisch.")))
}

.ig_classifier_read_csv <- function(path, separator = ",", decimal = ".") {
  if (!.ig_is_string(path) || !file.exists(path) || file.size(path) > 25 * 1024^2)
    .ig_abort("invalid_input", "CSV input must be an existing file of at most 25 MiB")
  if (!separator %in% c(",", ";", "\t") || !decimal %in% c(".", ",") || separator == decimal)
    .ig_abort("invalid_input", "choose distinct valid delimiter and decimal marks")
  d <- utils::read.table(path, header = TRUE, sep = separator, dec = decimal, quote = "\"",
    comment.char = "", stringsAsFactors = FALSE, check.names = FALSE, nrows = 100001L,
    fileEncoding = "UTF-8", na.strings = c("", "NA"))
  if (!nrow(d) || nrow(d) > 100000L || !ncol(d) || ncol(d) > 100L ||
      anyDuplicated(names(d)) || any(!nzchar(names(d))))
    .ig_abort("invalid_input", "CSV needs 1-100000 rows, 1-100 uniquely named columns")
  d
}

.ig_classifier_input_rule <- function(input, lang, cutoff) {
  selection <- input$selection %||% "manual"
  constraints <- list()
  for (nm in .ig_cls_constraint_names) {
    value <- input[[nm]] %||% ""
    if (nzchar(trimws(value))) constraints[[nm]] <- .ig_workbench_parse_number(value, lang)
  }
  costs <- NULL
  if(selection=="cost_utility") {
    costs <- list(fn=.ig_workbench_parse_number(input$cost_fn,lang),fp=.ig_workbench_parse_number(input$cost_fp,lang))
    for(nm in c("tp","tn")) {
      value <- input[[paste0("cost_",nm)]] %||% ""
      if(nzchar(trimws(value))) costs[[nm]] <- .ig_workbench_parse_number(value,lang)
    }
  }
  capacity_name <- input$capacity_kind %||% "max_alerts"
  if(!capacity_name %in% c("max_alerts","max_alert_rate")) .ig_abort("invalid_input","Unknown capacity unit")
  cutoff_rule(cutoff = if (selection == "manual") cutoff else NULL,
    selection = selection, selection_data = if (input$example %in% c("custom_theory", "binormal_equal_variance", "binormal_unequal_variance")) "predeclared_rule" else input$selection_data %||% "same_data",
    constraints = if (length(constraints)) constraints else NULL,
    optimise = if (selection == "constraint") input$optimise %||% "tpr" else NULL,
    capacity = if (selection == "capacity") stats::setNames(list(.ig_workbench_parse_number(input[[capacity_name]],lang)),capacity_name) else NULL,
    costs=costs,tie_break=input$tie_break %||% "fewest_alerts")
}

.ig_classifier_custom_spec <- function(input, lang, data = NULL) {
  num <- function(name) .ig_workbench_parse_number(input[[name]], lang)
  kind <- input$example
  if (kind == "custom_counts") {
    if (!identical(input$selection %||% "manual", "manual"))
      .ig_abort("invalid_input", "Fixed counts do not support cutoff selection")
    return(.ig_classifier_counts_rule(classification_spec(counts = list(tp=num("tp"),fn=num("fn"),fp=num("fp"),tn=num("tn"))),input,lang))
  }
  rule <- .ig_classifier_input_rule(input, lang, if (identical(input$selection %||% "manual", "manual")) num("cutoff") else NULL)
  if (kind == "custom_theory") return(classification_spec(theory = list(auc=num("auc"),sigma=num("sigma"),n_total=num("theory_n")),
      rule=rule,prevalence=num("theory_prevalence")))
  if (kind != "uploaded" || is.null(data)) .ig_abort("invalid_input", "Upload a CSV file and choose its columns first")
  if (!.ig_is_string(input$truth_column) || !length(input$score_columns) || !all(input$score_columns %in% names(data)) ||
      !input$truth_column %in% names(data)) .ig_abort("invalid_input", "Choose outcome and score columns from the uploaded data")
  if (input$truth_column %in% input$score_columns || (input$id_column %||% "") %in% input$score_columns)
    .ig_abort("invalid_input", "Outcome and identifier columns must not be selected as model scores")
  optional_column <- function(name) {
    value <- input[[name]] %||% ""
    if(!nzchar(value)) NULL else value
  }
  auxiliary <- c(optional_column("weight_column"),optional_column("cluster_column"))
  if(any(auxiliary %in% c(input$truth_column,input$score_columns)))
    .ig_abort("invalid_input","Weight and cluster columns must be distinct from outcome and model scores")
  model_inputs <- .ig_classifier_model_inputs(input,input$score_columns)
  classification_spec(data=data,truth=input$truth_column,positive_label=input$positive,scores=stats::setNames(input$score_columns, paste0("model_",seq_along(input$score_columns))),
    id=if (nzchar(input$id_column %||% "")) input$id_column else NULL,
    models=model_inputs$models,rule=rule,direction=model_inputs$direction,comparator=model_inputs$comparator,score_scale=model_inputs$scale,
    weights=optional_column("weight_column"),cluster=optional_column("cluster_column"),sampling=input$sampling %||% "unknown",
    evaluation_split=input$evaluation_split %||% "same_data",missing=input$missing %||% "error")
}
