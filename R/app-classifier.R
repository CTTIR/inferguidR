.ig_classifier_workspace_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tags$div(class = "ig-study",
    shiny::tags$aside(class = "ig-controls",
      shiny::selectInput(ns("example"), "Input source / Datenquelle", c(
        "Scored models / Modelle mit Scores" = "scored_synthetic",
        "Binormal, equal variance / gleiche Varianz" = "binormal_equal_variance",
        "Binormal, unequal variance / ungleiche Varianz" = "binormal_unequal_variance",
        "Fixed counts / Feste Anzahlen" = "counts_fixed",
        "Upload CSV / CSV hochladen" = "uploaded", "Enter counts / Anzahlen eingeben" = "custom_counts",
        "Enter theory / Modellparameter eingeben" = "custom_theory")),
      .ig_classifier_custom_ui(id),
      .ig_classifier_uncertainty_ui(ns),
      shiny::uiOutput(ns("cutoff_field")),
      shiny::uiOutput(ns("scenario_controls")),
      shiny::actionButton(ns("evaluate"), "Evaluate / Auswerten")),
    shiny::tags$main(class = "ig-figure",
      shiny::tags$div(role = "status", `aria-live` = "polite", shiny::textOutput(ns("message"))),
      shiny::uiOutput(ns("model_selector")),
      shiny::uiOutput(ns("basis_selector")),
      shiny::tableOutput(ns("cells")), shiny::tableOutput(ns("metrics")),
      shiny::tags$details(shiny::tags$summary("Intervals / Intervalle"),
        shiny::tags$p("Fixed-cutoff intervals omit cutoff-selection uncertainty. Scenario counts have no sampling interval. / Intervalle gelten bei festem Schwellenwert ohne Auswahlunsicherheit. Szenario-Anzahlen haben kein Stichprobenintervall."),
        shiny::tags$div(style="max-height:32rem;overflow:auto;",shiny::tableOutput(ns("intervals")))),
      shiny::tags$details(shiny::tags$summary("Paired AUC comparisons / Gepaarte AUC-Vergleiche"),
        shiny::tags$p("Same evaluation units; unadjusted pairwise comparisons. AUC does not measure calibration. / Gleiche Auswertungseinheiten; paarweise Vergleiche ohne Multiplizitaetskorrektur. AUC misst keine Kalibrierung."),
        shiny::tags$div(style="max-height:32rem;overflow:auto;",shiny::tableOutput(ns("comparisons")))),
      .ig_report_ui(ns("report")),
      shiny::selectInput(ns("curve_panel"), "Panel / Ansicht", c("ROC" = "roc", "Precision-recall" = "pr",
        "Mosaic / Mosaik" = "mosaic", "Score distributions / Score-Verteilungen" = "score_distribution",
        "Prevalence / Praevalenz" = "prevalence_curve", "Model comparison / Modellvergleich" = "model_comparison",
        "Metrics by cutoff / Kennzahlen nach Schwellenwert" = "metric_vs_cutoff",
        "Metric pair / Kennzahlenpaar" = "metric_pair", "Dependencies / Abhaengigkeiten" = "metric_network", "All metrics / Alle Kennzahlen" = "current_table",
        "Cutoff sweep / Schwellendurchlauf" = "sweep_table")),
      shiny::conditionalPanel(sprintf("['roc','pr'].includes(input['%s'])",ns("curve_panel")),
        shiny::selectInput(ns("interval_overlay"),"Marginal point intervals / Marginale Punktintervalle",
          c("None / Keine"="none","Wilson"="wilson","Clopper-Pearson"="clopper_pearson","Percentile bootstrap / Perzentil-Bootstrap"="bootstrap_percentile"))),
      shiny::uiOutput(ns("metric_controls")),
      shiny::conditionalPanel(sprintf("input['%s'] == 'mosaic'", ns("curve_panel")),
        shiny::selectInput(ns("mosaic_mode"), "Mosaic / Mosaik", c("Joint / Gemeinsam" = "joint", "Within class / Innerhalb der Klasse" = "within_class"))),
      shiny::conditionalPanel(sprintf("!['current_table', 'sweep_table'].includes(input['%s'])", ns("curve_panel")),
        shiny::uiOutput(ns("curve_container"))),
      shiny::conditionalPanel(sprintf("['current_table', 'sweep_table'].includes(input['%s'])", ns("curve_panel")),
        shiny::tags$div(style = "max-height: 32rem; overflow: auto;", shiny::tableOutput(ns("primary_table")))),
      shiny::tags$details(shiny::tags$summary("Panel data / Daten der Ansicht"),
        shiny::uiOutput(ns("data_selector")),
        shiny::tags$div(style = "max-height: 24rem; overflow: auto;", shiny::tableOutput(ns("curve_data"))))),
    shiny::tags$aside(class = "ig-inspector", shiny::textOutput(ns("description"))))
}

.ig_classifier_workspace_server <- function(id, language) {
  shiny::moduleServer(id, function(input, output, session) {
    result <- shiny::reactiveVal(NULL)
    description <- shiny::reactiveVal("")
    message <- shiny::reactiveVal("")
    last_language <- "en"
    prevalence_language <- "en"
    output$cutoff_field <- shiny::renderUI({
      if (input$example %in% c("counts_fixed", "custom_counts")) return(shiny::tags$p(
        if (language() == "de") "Feste Anzahlen beschreiben einen Betriebspunkt. Andere Schwellenwerte erfordern Score-Daten." else
          "Fixed counts describe one operating point. Other cutoffs require scored data."))
      lang <- language()
      current <- shiny::isolate(input$cutoff)
      value <- if (is.null(current)) .3 else
        tryCatch(.ig_workbench_parse_number(current, last_language), error = function(e) NULL)
      shown <- if (is.null(value)) current else .ig_format_double(value)
      if (!is.null(value) && lang == "de") shown <- sub(".", ",", shown, fixed = TRUE)
      last_language <<- lang
      shiny::textInput(session$ns("cutoff"), "Cutoff / Schwellenwert", shown)
    })
    output$scenario_controls <- shiny::renderUI({
      if (is.null(input$example) || !input$example %in% c("scored_synthetic", "counts_fixed", "uploaded", "custom_counts")) return(NULL)
      lang <- language()
      current <- shiny::isolate(input$prevalence) %||% ""
      value <- tryCatch(.ig_workbench_parse_number(current, prevalence_language), error = function(e) NULL)
      shown <- if (is.null(value)) current else .ig_format_double(value)
      if (!is.null(value) && lang == "de") shown <- sub(".", ",", shown, fixed = TRUE)
      prevalence_language <<- lang
      shiny::tagList(
        shiny::checkboxInput(session$ns("scenario"), "Assumed prevalence scenario / Angenommene Praevalenz", shiny::isolate(isTRUE(input$scenario))),
        shiny::textInput(session$ns("prevalence"), "Scenario prevalence / Szenario-Praevalenz", shown))
    })
    custom_input_language <- shiny::isolate(language())
    shiny::observeEvent(language(), {
      lang <- language()
      for (nm in c("tp", "fn", "fp", "tn", "auc", "sigma", "theory_n", "theory_prevalence",
                   "tpr_min", "tnr_min", "npv_min", "fnr_max", "fdr_max", "fpr_max", "ppv_min", "alert_rate_max", "max_alerts", "max_alert_rate", "cost_fn", "cost_fp", "cost_tp", "cost_tn", "confidence_level", "bootstrap_replicates", "bootstrap_seed")) {
        value <- input[[nm]]
        if (is.null(value) || !nzchar(trimws(value))) next
        numeric <- tryCatch(.ig_workbench_parse_number(value, custom_input_language), error = function(e) NULL)
        if (is.null(numeric)) next
        shown <- .ig_format_double(numeric)
        if (lang == "de") shown <- sub(".", ",", shown, fixed = TRUE)
        shiny::updateTextInput(session, nm, value = shown)
      }
      custom_input_language <<- lang
    }, ignoreInit = TRUE)
    uploaded_data <- shiny::reactive({
      shiny::req(input$csv)
      .ig_classifier_read_csv(input$csv$datapath, input$separator %||% ",", input$decimal %||% ".")
    })
    output$columns <- shiny::renderUI({
      d <- uploaded_data()
      choices <- names(d)
      shiny::tagList(
        shiny::selectInput(session$ns("truth_column"), "Outcome / Ergebnis", choices),
        shiny::selectInput(session$ns("score_columns"), "Score columns / Score-Spalten", choices, multiple = TRUE),
        shiny::selectInput(session$ns("id_column"), "Unique ID (optional) / Eindeutige ID (optional)", c("None / Keine"="", choices)),
        shiny::selectInput(session$ns("weight_column"), "Frequency-weight column / Haeufigkeitsgewichte", c("None / Keine"="", choices)),
        shiny::selectInput(session$ns("cluster_column"), "Cluster column / Cluster-Spalte", c("None / Keine"="", choices)))
    })
    output$model_settings <- shiny::renderUI(.ig_classifier_model_ui(session$ns,input,input$score_columns %||% character()))
    custom <- shiny::reactive(input$example %in% c("uploaded", "custom_counts", "custom_theory"))
    candidate <- shiny::reactive({
      if (custom()) {
        spec <- .ig_classifier_custom_spec(input, language(), if (input$example == "uploaded") uploaded_data() else NULL)
      } else {
        spec <- classification_examples(input$example)$spec
        if (input$example != "counts_fixed") {
          cutoff <- if (identical(input$selection %||% "manual", "manual")) .ig_workbench_parse_number(input$cutoff, language()) else NULL
          spec <- update(spec, rule = .ig_classifier_input_rule(input, language(), cutoff))
        } else if (!identical(input$selection %||% "manual", "manual"))
          .ig_abort("invalid_input", "Fixed counts do not support cutoff selection")
      }
      if(input$example=="counts_fixed") spec <- .ig_classifier_counts_rule(spec,input,language())
      if (input$example %in% c("scored_synthetic", "counts_fixed", "uploaded", "custom_counts") && isTRUE(input$scenario))
        spec <- update(spec, prevalence = .ig_workbench_parse_number(input$prevalence, language()))
      update(spec,uncertainty=.ig_classifier_uncertainty_input(input,language()))
    })
    shiny::observeEvent(input$evaluate, {
      tryCatch({
        example <- if (custom()) list(description = "User-provided inputs / Eigene Eingaben") else classification_examples(input$example)
        spec <- candidate()
        value <- evaluate_cutoff(spec)
        result(value)
        description(example$description)
        message(if (input$example %in% c("counts_fixed", "custom_counts")) {
          if (language() == "de") "Feste Anzahlen: Der Schwellenwert wird hier nicht ver\u00e4ndert." else
            "Fixed counts: no cutoff selection is performed."
        } else if (language() == "de") "Modelle ausgewertet." else "Models evaluated.")
      }, error = function(e) message(paste(conditionMessage(e), if (!is.null(result()))
        if (language() == "de") "Angezeigt wird die vorherige Auswertung." else "The previous evaluation remains displayed.")))
    }, ignoreInit = TRUE)
    output$model_selector <- shiny::renderUI({
      shiny::req(result())
      ids <- vapply(result()$models, function(m) m$model_id, "")
      selected <- shiny::isolate(input$model)
      if (is.null(selected) || !selected %in% ids) selected <- ids[1]
      shiny::selectInput(session$ns("model"), if (language() == "de") "Modell" else "Model", ids, selected)
    })
    selected_model <- shiny::reactive({
      shiny::req(result())
      ids <- vapply(result()$models, function(m) m$model_id, "")
      if (is.null(input$model) || !input$model %in% ids) ids[1] else input$model
    })
    output$basis_selector <- shiny::renderUI({
      shiny::req(result())
      choices <- stats::setNames("observed", if (language() == "de") "Beobachtet / Modellbasis" else "Observed / model basis")
      if (identical(result()$prevalence$kind, "scenario") && result()$input_kind != "theoretical_binormal")
        choices <- c(choices, stats::setNames("scenario", if (language() == "de") "Angenommene Praevalenz" else "Assumed prevalence"))
      selected <- shiny::isolate(input$basis)
      if (is.null(selected) || !selected %in% choices) selected <- "observed"
      shiny::selectInput(session$ns("basis"), if (language() == "de") "Bezugsbasis" else "Basis", choices, selected)
    })
    selected_basis <- shiny::reactive({
      shiny::req(result())
      if (identical(input$basis, "scenario") && identical(result()$prevalence$kind, "scenario") &&
          result()$input_kind != "theoretical_binormal") "scenario" else "observed"
    })
    output$metric_controls <- shiny::renderUI({
      panel <- input$curve_panel %||% "roc"
      if (panel == "metric_vs_cutoff") {
        labels <- if (language() == "de") c("Nach wahrer Klasse", "Nach Entscheidung", "Zusammenfassung", "Belastung", "Likelihood-Verhaeltnisse") else
          c("Within truth", "Within decisions", "Summary", "Burden", "Likelihood ratios")
        return(shiny::selectInput(session$ns("metric_family"), if (language() == "de") "Kennzahlengruppe" else "Metric family",
          stats::setNames(names(.ig_cls_families), labels), shiny::isolate(input$metric_family) %||% "within_truth"))
      }
      if (!panel %in% c("metric_pair", "metric_network")) return(NULL)
      ids <- .ig_cls_operating_ids()
      labels <- .ig_cls_short_label(ids, language())
      choices <- stats::setNames(ids, paste(labels, ids, sep = " / "))
      if (panel == "metric_network") return(shiny::selectInput(session$ns("metric_focus"),
        if (language() == "de") "Abhaengigkeiten von" else "Dependencies of", choices,
        shiny::isolate(input$metric_focus) %||% "ppv"))
      shiny::tagList(
        shiny::selectInput(session$ns("metric_x"), "X", choices, shiny::isolate(input$metric_x) %||% "fpr"),
        shiny::selectInput(session$ns("metric_y"), "Y", choices, shiny::isolate(input$metric_y) %||% "tpr"),
        shiny::checkboxInput(session$ns("log_x"), "Log X", shiny::isolate(isTRUE(input$log_x))),
        shiny::checkboxInput(session$ns("log_y"), "Log Y", shiny::isolate(isTRUE(input$log_y))))
    })
    curve_data <- shiny::reactive({
      shiny::req(result())
      shiny::validate(shiny::need(!(identical(result()$input_kind, "counts_2x2") && (input$curve_panel %||% "roc") %in% c("roc", "pr", "score_distribution")),
        if (language() == "de") "Feste Anzahlen bestimmen keine Kurve; die Tabellen zeigen den einzelnen Betriebspunkt." else
          "Fixed counts do not define a curve; the tables show the single operating point."))
      panel_id <- input$curve_panel %||% "roc"
      basis <- if (panel_id %in% c("score_distribution", "prevalence_curve", "model_comparison")) "observed" else selected_basis()
      pd <- tryCatch(classification_plot_data(result(), panel_id, selected_model(), basis = basis, lang = language(),
        metric_x = input$metric_x %||% "fpr", metric_y = input$metric_y %||% "tpr"),
        error = function(e) shiny::validate(shiny::need(FALSE, conditionMessage(e))))
      if(panel_id %in% c("roc","pr") && !is.null(input$interval_overlay) && input$interval_overlay!="none")
        pd$data$intervals <- .ig_cls_interval_overlay(result(),pd,input$interval_overlay,language())$data
      pd
    })
    output$curve_container <- shiny::renderUI({
      height <- 420
      if (identical(input$curve_panel, "metric_network")) {
        pd <- curve_data()
        focus <- input$metric_focus %||% "ppv"
        n <- length(unique(pd$data$edges$from[pd$data$edges$to == focus]))
        height <- max(height, 160 + 90 * n)
      }
      shiny::plotOutput(session$ns("curve"), height = paste0(height, "px"))
    })
    output$curve <- shiny::renderPlot({
      pd <- curve_data()
      if (pd$panel %in% c("current_table", "sweep_table")) return(NULL)
      plot_classification(result(), pd$panel, selected_model(), basis = pd$basis, lang = language(),
                          mosaic_mode = input$mosaic_mode %||% "joint",
                          metric_family = input$metric_family %||% "within_truth",
                          metric_x = input$metric_x %||% "fpr", metric_y = input$metric_y %||% "tpr",
                          log_x = isTRUE(input$log_x), log_y = isTRUE(input$log_y), metric_focus = input$metric_focus %||% "ppv",
                          interval_method=if(pd$panel %in% c("roc","pr") && !is.null(input$interval_overlay) && input$interval_overlay!="none") input$interval_overlay else NULL)
    }, res = 110)
    output$data_selector <- shiny::renderUI({
      pd <- curve_data()
      names <- names(pd$data)[vapply(pd$data, is.data.frame, logical(1))]
      selected <- shiny::isolate(input$data_table)
      if (is.null(selected) || !selected %in% names) selected <- names[1]
      shiny::selectInput(session$ns("data_table"), if (language() == "de") "Datentabelle" else "Data table", names, selected)
    })
    displayed_table <- shiny::reactive({
      pd <- curve_data()
      names <- names(pd$data)[vapply(pd$data, is.data.frame, logical(1))]
      key <- if (!is.null(input$data_table) && input$data_table %in% names) input$data_table else names[1]
      tab <- pd$data[[key]]
      for (nm in names(tab)) if (is.numeric(tab[[nm]])) {
        tab[[nm]] <- .ig_format_double(tab[[nm]])
        if (language() == "de") tab[[nm]] <- sub(".", ",", tab[[nm]], fixed = TRUE)
      }
      tab
    })
    output$curve_data <- shiny::renderTable(displayed_table())
    output$primary_table <- shiny::renderTable(displayed_table())
    output$message <- shiny::renderText({
      if (!is.null(result())) {
        hash <- tryCatch(candidate()$spec_hash, error = function(e) NULL)
        if (!identical(hash, result()$spec_hash)) return(paste(message(), if (language() == "de")
          "Entwurf ge\u00e4ndert; angezeigt wird die vorherige Auswertung." else
          "Draft changed; the previous evaluation remains displayed."))
      }
      message()
    })
    output$description <- shiny::renderText(description())
    panel <- shiny::reactive({
      shiny::req(result())
      classification_plot_data(result(), "confusion_table", selected_model(), basis = selected_basis(), lang = language())
    })
    output$cells <- shiny::renderTable({
      tab <- panel()$data$cells[, c("row_label", "column_label", "count", "count_kind")]
      tab$count <- format(tab$count, digits = 6, decimal.mark = if (language() == "de") "," else ".", trim = TRUE)
      names(tab) <- if (language() == "de") c("Entscheidung", "Wahrer Zustand", "Anzahl", "Art") else
        c("Decision", "True state", "Count", "Kind")
      tab
    }, digits = 3)
    output$metrics <- shiny::renderTable({
      tab <- panel()$data$tiles[, c("short", "value", "status", "denominator_value")]
      tab$value <- format(tab$value, digits = 6, decimal.mark = if (language() == "de") "," else ".", trim = TRUE)
      tab$denominator_value <- format(tab$denominator_value, digits = 6, decimal.mark = if (language() == "de") "," else ".", trim = TRUE)
      names(tab) <- if (language() == "de") c("Kennzahl", "Wert", "Status", "Nenner") else
        c("Metric", "Value", "Status", "Denominator")
      tab
    }, digits = 6)
    output$intervals <- shiny::renderTable({
      shiny::req(result())
      .ig_classifier_display_table(.ig_classifier_interval_table(result(),selected_model(),selected_basis(),language()),language())
    },digits=6)
    output$comparisons <- shiny::renderTable({
      shiny::req(result())
      .ig_classifier_display_table(.ig_classifier_comparison_table(result(),language()),language())
    },digits=6)
    report_module <- .ig_report_server("report", shiny::reactive({
      analysis_envelope(candidate(), result())
    }), language, eligible = shiny::reactive({
      s <- tryCatch(candidate(), error = function(e) NULL)
      !is.null(s) && !is.null(result()) && identical(s$spec_hash, result()$spec_hash)
    }))
    list(result = shiny::reactive(result()), report = report_module)
  })
}
