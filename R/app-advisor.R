# Guided-question module; all routing remains in the tested advisor API.
.ig_advisor_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::uiOutput(ns("question")),
    shiny::uiOutput(ns("navigation")),
    shiny::tags$div(role = "status", `aria-live` = "polite", shiny::textOutput(ns("message"))),
    shiny::uiOutput(ns("recommendation")),
    shiny::uiOutput(ns("override_controls")),
    shiny::uiOutput(ns("sources")),
    shiny::uiOutput(ns("r_code_controls")),
    shiny::uiOutput(ns("worked_lesson")),
    shiny::tags$details(shiny::tags$summary("Session / Sitzung"),
      shiny::downloadButton(ns("save"), "Save / Speichern"),
      shiny::fileInput(ns("restore"), "Restore / Wiederherstellen", accept = ".json")),
    shiny::uiOutput(ns("path")))
}

.ig_advisor_server <- function(id, language) {
  shiny::moduleServer(id, function(input, output, session) {
    state <- shiny::reactiveVal(advisor_session(language = shiny::isolate(language())))
    message <- shiny::reactiveVal("")
    text <- function(en, de) if (language() == "de") de else en
    attempt <- function(expr) tryCatch({
      value <- force(expr)
      state(value)
      message("")
    }, error = function(e) message(conditionMessage(e)))
    shiny::observeEvent(language(), {
      state(advisor_set_language(state(), language()))
    }, ignoreInit = TRUE)
    output$question <- shiny::renderUI({
      q <- advisor_question(state())
      if (is.null(q)) return(NULL)
      shiny::tags$fieldset(
        shiny::tags$legend(q$text),
        shiny::tags$p(q$why_it_matters),
        shiny::radioButtons(session$ns(paste0("answer_", q$question_id)),
          label = text("Choose an answer", "Antwort ausw\u00e4hlen"),
          choices = stats::setNames(vapply(q$options, `[[`, "", "option_id"),
                                    vapply(q$options, `[[`, "", "label")), selected = character()),
        shiny::tags$details(shiny::tags$summary(text("What the answers mean", "Bedeutung der Antworten")),
          lapply(q$options, function(o) shiny::tags$p(shiny::tags$strong(o$label), ": ", o$explanation))))
    })
    output$navigation <- shiny::renderUI({
      shiny::tagList(
        if (!is.null(advisor_question(state())))
          shiny::actionButton(session$ns("next_answer"), text("Continue", "Weiter")),
        if (length(state()$answers))
          shiny::actionButton(session$ns("back"), text("Back", "Zur\u00fcck")),
        shiny::actionButton(session$ns("restart"), text("Start again", "Neu beginnen")))
    })
    shiny::observeEvent(input$next_answer, {
      q <- advisor_question(state())
      if (is.null(q)) return()
      answer <- input[[paste0("answer_", q$question_id)]]
      if (is.null(answer) || !nzchar(answer)) {
        message(text("Choose an answer before continuing.", "Bitte zuerst eine Antwort ausw\u00e4hlen."))
        return()
      }
      attempt(advisor_answer(state(), answer, question_id = q$question_id))
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$back, {
      if (length(state()$answers)) attempt(advisor_back(state()))
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$restart, attempt(advisor_session(language = language())), ignoreInit = TRUE)
    output$message <- shiny::renderText(message())
    output$recommendation <- shiny::renderUI({
      r <- advisor_recommendation(state())
      if (!is.null(advisor_question(state()))) return(NULL)
      labels <- stats::setNames(lapply(r$missing_information, function(id) advisor_question(state(), question_id = id)$text), r$missing_information)
      .ig_advisor_card(r, language(), advisor_rules(), labels)
    })
    code_result <- shiny::reactive({
      tryCatch(list(code = advisor_r_code(state(),
        alternative = if (identical(advisor_recommendation(state())$design$answers[["q.direction"]], "one_sided_prespecified")) {
          if (identical(input$code_alternative, "")) NULL else input$code_alternative
        } else NULL)),
        error = function(e) list(error = conditionMessage(e)))
    })
    output$r_code_controls <- shiny::renderUI({
      if (!is.null(advisor_question(state()))) return(NULL)
      r <- advisor_recommendation(state())
      if (!identical(r$r_implementation$template_status, "available")) return(NULL)
      one_sided <- identical(r$design$answers[["q.direction"]], "one_sided_prespecified")
      shiny::tags$details(class = "ig-r-example",
        shiny::tags$summary(text("R analysis example", "R-Analysebeispiel")),
        shiny::tags$p(text("Executable synthetic example for the selected method, using a 95% confidence setting. The example states which intervals are available. It does not analyse your data. Review the recorded assumptions and data requirements before adapting it.",
          "Ausfuehrbares synthetisches Beispiel fuer die gewaehlte Methode mit einem Konfidenzniveau von 95%. Das Beispiel gibt an, welche Intervalle verfuegbar sind. Es analysiert nicht Ihre Daten. Pruefen Sie die erfassten Annahmen und Datenanforderungen vor einer Anpassung.")),
        if (one_sided) shiny::selectInput(session$ns("code_alternative"),
          text("Pre-specified alternative (meaning stated in the example)", "Vorab festgelegte Alternative (Bedeutung im Beispiel)"),
          choices = stats::setNames(c("", "greater", "less"), c(text("Choose the pre-specified sign", "Vorab festgelegtes Vorzeichen waehlen"),
            text("Greater", "Groesser"), text("Less", "Kleiner"))), selected = ""),
        shiny::uiOutput(session$ns("r_code_download")),
        shiny::verbatimTextOutput(session$ns("r_code_text")))
    })
    output$r_code_text <- shiny::renderText(code_result()$code %||% code_result()$error)
    output$r_code_download <- shiny::renderUI({
      if (is.null(code_result()$code)) return(NULL)
      shiny::downloadButton(session$ns("save_code"), text("Download R example", "R-Beispiel herunterladen"))
    })
    output$save_code <- shiny::downloadHandler(
      filename = function() paste0("inferguidR-", advisor_recommendation(state())$preferred$analysis_variant_id, ".R"),
      content = function(file) {
        result <- shiny::isolate(code_result())
        if (is.null(result$code)) stop(result$error, call. = FALSE)
        writeLines(result$code, file, useBytes = TRUE)
      }, contentType = "text/plain")
    selected_lesson <- shiny::reactive({
      if (!is.null(advisor_question(state()))) return(NULL)
      .ig_advisor_lesson(advisor_recommendation(state()), language())
    })
    output$worked_lesson <- shiny::renderUI({
      lesson <- selected_lesson()
      if (is.null(lesson)) return(NULL)
      shiny::tags$details(class = "ig-recommendation",
        shiny::tags$summary(text("Original worked lesson (available offline)",
          "Eigenstaendige Lektion (offline verfuegbar)")),
        shiny::downloadButton(session$ns("save_lesson"), text("Download lesson", "Lektion herunterladen")),
        shiny::downloadButton(session$ns("save_lesson_code"), text("Download lesson R script", "R-Skript der Lektion herunterladen")),
        shiny::includeHTML(lesson$html))
    })
    output$save_lesson <- shiny::downloadHandler(
      filename = function() paste0("boschloo-", shiny::isolate(language()), ".md"),
      content = function(file) {
        lesson <- shiny::isolate(selected_lesson())
        if (is.null(lesson)) stop("No lesson for the selected variant.", call. = FALSE)
        if (!file.copy(lesson$markdown, file, overwrite = TRUE)) stop("Cannot write lesson.", call. = FALSE)
      }, contentType = "text/markdown; charset=utf-8")
    output$save_lesson_code <- shiny::downloadHandler(
      filename = function() "boschloo.R",
      content = function(file) {
        lesson <- shiny::isolate(selected_lesson())
        if (is.null(lesson)) stop("No lesson for the selected variant.", call. = FALSE)
        if (!file.copy(lesson$script, file, overwrite = TRUE)) stop("Cannot write lesson script.", call. = FALSE)
      }, contentType = "text/plain; charset=utf-8")
    output$sources <- shiny::renderUI({
      .ig_advisor_source_view(advisor_recommendation(state()), language(), advisor_rules())
    })
    output$override_controls <- shiny::renderUI({
      if (!is.null(advisor_question(state()))) return(NULL)
      r <- advisor_recommendation(state())
      allowed <- Filter(function(a) isTRUE(a$override_allowed), r$alternatives)
      if (!length(allowed)) return(NULL)
      ids <- vapply(allowed, `[[`, "", "analysis_variant_id")
      names <- vapply(allowed, function(a) a[[paste0("name_", language())]], "")
      shiny::tags$details(shiny::tags$summary(text("Choose a listed alternative", "Aufgefuehrte Alternative waehlen")),
        shiny::tags$p(text("Review its assumptions and tradeoffs above. A changed method is recorded as a conditional recommendation.",
          "Pruefen Sie die Annahmen und Abwaegungen oben. Ein Methodenwechsel wird als bedingte Empfehlung erfasst.")),
        shiny::selectInput(session$ns("override_variant"), text("Alternative", "Alternative"),
          c(stats::setNames("", text("Choose", "Bitte waehlen")), stats::setNames(ids, names)), selected = ""),
        shiny::actionButton(session$ns("apply_override"), text("Apply alternative", "Alternative uebernehmen")))
    })
    shiny::observeEvent(input$apply_override, {
      if (is.null(input$override_variant) || !nzchar(input$override_variant)) return()
      attempt(advisor_override(state(), input$override_variant))
    }, ignoreInit = TRUE)
    output$path <- shiny::renderUI({
      s <- state()
      if (!length(s$answers)) return(NULL)
      shiny::tags$section(shiny::tags$h3(text("Your reasoning path", "Ihr Entscheidungsweg")),
        shiny::tags$ol(lapply(s$answers, function(a) {
          q <- advisor_question(s, question_id = a$question_id)
          option <- Filter(function(o) identical(o$option_id, a$option_id), q$options)[[1]]
          shiny::tags$li(shiny::tags$strong(q$text), ": ", option$label)
        })))
    })
    output$save <- shiny::downloadHandler(
      filename = function() "inferguidR-advisor-session.json",
      content = function(file) advisor_save(shiny::isolate(state()), file = file),
      contentType = "application/json")
    shiny::observeEvent(input$restore, {
      f <- input$restore
      if (is.null(f)) return()
      if (f$size > 5 * 1024^2) {
        message(text("Session file exceeds 5 MB.", "Sitzungsdatei ist gr\u00f6\u00dfer als 5 MB."))
        return()
      }
      attempt(advisor_set_language(advisor_restore(file = f$datapath), language()))
    }, ignoreInit = TRUE)
    list(state = shiny::reactive(state()), recommendation = shiny::reactive(advisor_recommendation(state())))
  })
}
