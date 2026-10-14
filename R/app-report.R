.ig_report_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tags$details(shiny::tags$summary("Report downloads / Berichte herunterladen"),
    shiny::checkboxInput(ns("include_resources"), "Include verified teaching references when available / Gepruefte Lehrquellen aufnehmen, soweit verfuegbar", value = TRUE),
    shiny::actionButton(ns("freeze"), "Freeze report / Bericht festhalten"),
    shiny::textOutput(ns("status")), shiny::uiOutput(ns("downloads")))
}

.ig_report_server <- function(id, envelope, language, eligible = shiny::reactive(TRUE), sweeps = shiny::reactive(list())) {
  shiny::moduleServer(id, function(input, output, session) {
    snapshot <- shiny::reactiveVal(NULL)
    message <- shiny::reactiveVal("")
    shiny::observeEvent(input$freeze, {
      tryCatch({
        if (!isTRUE(eligible())) .ig_abort("invalid_input", "Calculate the current settings before freezing a final report.")
        e <- envelope()
        entries <- sweeps()
        incomplete_sweep <- any(vapply(entries,function(entry) any(entry$sweep$status!="ok"),TRUE))
        snapshot(create_report_snapshot(e, locale = language(), sweeps = entries,
          resources = if(isTRUE(input$include_resources)) .ig_report_envelope_resources(e) else list(),
          report_kind = if(incomplete_sweep) "incomplete_settings" else "single_analysis"))
        message("")
      }, error = function(e) message(conditionMessage(e)))
    }, ignoreInit = TRUE)
    output$status <- shiny::renderText({
      s <- snapshot()
      if (nzchar(message())) return(message())
      if (is.null(s)) return(if (language() == "de") "Noch kein Bericht festgehalten." else "No report frozen yet.")
      paste(if (language() == "de") "Downloads verwenden diesen festgehaltenen Stand:" else
        "Downloads use this frozen snapshot:", s$snapshot_id, s$created_at,
        if (language() == "de") "Neue Eingaben aendern ihn nicht. Bei Bedarf erneut festhalten." else
          "New inputs do not change it. Freeze again when needed.")
    })
    output$downloads <- shiny::renderUI({
      if (is.null(snapshot())) return(NULL)
      shiny::tagList(shiny::downloadButton(session$ns("json"), "JSON"),
        shiny::downloadButton(session$ns("txt"), "TXT"),
        if (report_pdf_capabilities()$available) shiny::downloadButton(session$ns("pdf"), "PDF") else
          shiny::tags$p(report_pdf_capabilities()$reason))
    })
    for (fmt in c("json", "txt", "pdf")) local({
      format <- fmt
      output[[format]] <- shiny::downloadHandler(
        filename = function() { shiny::req(snapshot()); paste0(snapshot()$snapshot_id, ".", format) },
        contentType = switch(format, json = "application/json", txt = "text/plain; charset=utf-8", pdf = "application/pdf"),
        content = function(file) {
          s <- snapshot(); shiny::req(s)
          export_report(s, file, format, overwrite = TRUE)
        })
    })
    list(snapshot = shiny::reactive(snapshot()))
  })
}
