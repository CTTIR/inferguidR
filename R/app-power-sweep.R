.ig_power_sweep_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tags$details(shiny::tags$summary("Compare input values / Eingabewerte vergleichen"),
    shiny::uiOutput(ns("controls")),
    shiny::tags$div(role="status", `aria-live`="polite", shiny::textOutput(ns("message"))),
    shiny::uiOutput(ns("curve_panel")),
    shiny::tableOutput(ns("rows")),
    shiny::uiOutput(ns("selection")),
    shiny::uiOutput(ns("preview_label")))
}

.ig_power_sweep_values <- function(value, lang) {
  if (!.ig_is_string(value) || !nzchar(trimws(value)))
    .ig_abort("invalid_input", "Enter values separated by semicolons.")
  parts <- strsplit(value, ";", fixed=TRUE)[[1]]
  if (length(parts) < 2L || length(parts) > 21L || grepl(";[[:space:]]*$",value))
    .ig_abort("invalid_input", "Enter between 2 and 21 values separated by semicolons.")
  vapply(parts, .ig_workbench_parse_number, numeric(1), locale=lang, USE.NAMES=FALSE)
}

.ig_power_sweep_server <- function(id, state, language, research, apply_value) {
  shiny::moduleServer(id, function(input, output, session) {
    sweep <- shiny::reactiveVal(NULL)
    preview_row <- shiny::reactiveVal(NULL)
    message <- shiny::reactiveVal("")
    values_language <- NULL
    text <- function(en,de) if(language()=="de") de else en
    ready <- function() {
      s <- state()
      !is.null(s) && !isTRUE(s$dirty) && identical(s$committed$result$status$code,"ok") && isTRUE(research())
    }
    current_sweep <- shiny::reactive({
      sw <- sweep();s <- state()
      if(is.null(sw) || is.null(s) || !identical(attr(sw,"base_spec_hash"),s$committed$spec_hash)) return(NULL)
      sw
    })
    output$controls <- shiny::renderUI({
      if(!ready()) return(shiny::tags$p(text("Calculate a valid design and resolve draft edits before comparing values.",
        "Berechnen Sie ein gueltiges Design und klaeren Sie Entwurfsaenderungen vor dem Wertevergleich.")))
      spec <- state()$committed$spec;def <- .ig_get_method(spec$method_id)
      choices <- names(spec$inputs)[vapply(spec$inputs,function(v) is.numeric(v) && length(v)==1L,logical(1))]
      labels <- vapply(choices,function(nm) .ig_workbench_label(def$params[[nm]]$label,language(),nm),"")
      selected <- shiny::isolate(input$free)
      if(is.null(selected) || !selected %in% choices) selected <- choices[1]
      shown_values <- shiny::isolate(input$values) %||% ""
      parsed <- tryCatch(.ig_power_sweep_values(shown_values,values_language %||% language()),error=function(e) NULL)
      if(!is.null(parsed)) {
        shown_values <- paste(vapply(parsed,.ig_format_double,""),collapse="; ")
        if(language()=="de") shown_values <- gsub(".",",",shown_values,fixed=TRUE)
      }
      values_language <<- language()
      shiny::tagList(
        shiny::tags$p(text("Uses the calculated design above. Every row is evaluated separately; failed rows retain their reason. Enter 2 to 21 values separated by semicolons. A preview changes the design only when you apply it.",
          "Verwendet das oben berechnete Design. Jede Zeile wird einzeln berechnet; fehlgeschlagene Zeilen behalten ihre Begruendung. Geben Sie 2 bis 21 Werte getrennt durch Semikolons ein. Eine Vorschau aendert das Design erst nach Uebernahme.")),
        shiny::selectInput(session$ns("free"),text("Input to vary","Zu variierende Eingabe"),stats::setNames(choices,labels),selected=selected),
        shiny::textInput(session$ns("values"),text("Values (semicolon separated)","Werte (durch Semikolon getrennt)"),value=shown_values),
        shiny::actionButton(session$ns("run"),text("Calculate comparison","Vergleich berechnen")))
    })
    shiny::observeEvent(input$run,{
      tryCatch({
        if(!ready()) .ig_abort("invalid_input","Resolve the draft and calculate before running a comparison.")
        values <- .ig_power_sweep_values(input$values,language())
        sw <- sweep_design(state()$committed$spec,input$free,values)
        sweep(sw);preview_row(NULL)
        s <- state();s$preview <- NULL;state(s);message("")
      },error=function(e) message(conditionMessage(e)))
    },ignoreInit=TRUE)
    output$message <- shiny::renderText(message())
    output$curve_panel <- shiny::renderUI({
      if(is.null(current_sweep())) return(NULL)
      shiny::plotOutput(session$ns("curve"),height="350px")
    })
    output$curve <- shiny::renderPlot({
      sw <- current_sweep();shiny::req(sw)
      solved <- grep("^solved_", names(sw), value=TRUE)
      def <- .ig_get_method(attr(sw,"method_id"))
      tab <- do.call(rbind, lapply(solved, function(nm) {
        target <- sub("^solved_", "", nm)
        data.frame(value=sw$value, solved=sw[[nm]], status=sw$status,
          target=.ig_workbench_label(def$params[[target]]$label,language(),target))
      }))
      tab <- tab[tab$status=="ok" & is.finite(tab$value) & is.finite(tab$solved),,drop=FALSE]
      shiny::validate(shiny::need(nrow(tab)>0,text("No successful finite points to plot.","Keine erfolgreichen endlichen Punkte darstellbar.")))
      value <- solved_value <- target <- NULL
      names(tab)[names(tab)=="solved"] <- "solved_value"
      ggplot2::ggplot(tab,ggplot2::aes(x=value,y=solved_value)) + ggplot2::geom_point(colour="#334E9A",size=2.5) +
        ggplot2::facet_wrap(~target,scales="free_y") +
        ggplot2::labs(x=.ig_workbench_label(def$params[[attr(sw,"free")]]$label,language(),attr(sw,"free")),
          y=text("Solved value","Berechneter Wert"),
          caption=text("Evaluated points only; no interpolation. See the table for failed rows.",
            "Nur berechnete Punkte; keine Interpolation. Fehlgeschlagene Zeilen stehen in der Tabelle.")) +
        ggplot2::theme_minimal(base_size=13)
    },res=110)
    output$rows <- shiny::renderTable({
      sw <- current_sweep();shiny::req(sw)
      columns <- c("value","status","reason",grep("^solved_",names(sw),value=TRUE),
        if(!"solved_power" %in% names(sw)) "achieved_power")
      tab <- sw[,unique(columns),drop=FALSE]
      tab <- cbind(row=seq_len(nrow(tab)),tab)
      for(nm in names(tab)) if(is.numeric(tab[[nm]])) {
        tab[[nm]] <- vapply(tab[[nm]],function(v) if(is.na(v)) "" else .ig_format_double(v),"")
        if(language()=="de") tab[[nm]] <- sub(".",",",tab[[nm]],fixed=TRUE)
      }
      def <- .ig_get_method(attr(sw,"method_id"))
      names(tab) <- vapply(names(tab),function(nm) {
        if(startsWith(nm,"solved_")) {
          target <- sub("^solved_","",nm)
          return(.ig_workbench_label(def$params[[target]]$label,language(),target))
        }
        switch(nm,row=text("Row","Zeile"),value=text("Input value","Eingabewert"),
          status="Status",reason=text("Reason","Begruendung"),achieved_power=text("Achieved power","Erreichte Teststaerke"))
      },"")
      tab
    },striped=FALSE)
    output$selection <- shiny::renderUI({
      sw <- current_sweep();if(is.null(sw) || !ready()) return(NULL)
      ok <- which(sw$status=="ok")
      if(!length(ok)) return(shiny::tags$p(text("No successful row is available for preview.","Keine erfolgreiche Zeile fuer eine Vorschau vorhanden.")))
      selected <- shiny::isolate(input$row) %||% ""
      shiny::tagList(
        shiny::selectInput(session$ns("row"),text("Row to preview","Zeile fuer Vorschau"),
          c(stats::setNames("",text("Choose a row","Zeile waehlen")),stats::setNames(as.character(ok),paste(ok,sw$value[ok],sep=": "))),selected=selected),
        shiny::actionButton(session$ns("preview"),text("Preview row","Zeile ansehen")),
        shiny::actionButton(session$ns("apply"),text("Apply preview","Vorschau uebernehmen")))
    })
    row_number <- function(sw) {
      i <- suppressWarnings(as.integer(input$row))
      if(length(i)!=1L || is.na(i) || i<1L || i>nrow(sw) || sw$status[i]!="ok")
        .ig_abort("invalid_input","Choose a successful row.")
      i
    }
    shiny::observeEvent(input$preview,{
      tryCatch({
        sw <- current_sweep()
        if(is.null(sw) || !ready()) .ig_abort("invalid_input","This comparison no longer matches the current design.")
        i <- row_number(sw)
        state(.ig_workbench_preview(state(),analysis_envelope(sw$spec[[i]],sw$result[[i]])))
        preview_row(i);message("")
      },error=function(e) message(conditionMessage(e)))
    },ignoreInit=TRUE)
    output$preview_label <- shiny::renderUI({
      sw <- current_sweep();s <- state();i <- preview_row()
      if(is.null(sw) || is.null(s$preview) || is.null(i)) return(NULL)
      shiny::tagList(shiny::tags$p(shiny::tags$strong(paste(text("Preview only \u2014 row","Nur Vorschau \u2014 Zeile"),i)),
        text(". The main result is unchanged until Apply preview.",". Das Hauptergebnis bleibt bis zur Uebernahme unveraendert.")),
        shiny::plotOutput(session$ns("preview_plot"),height="350px"))
    })
    output$preview_plot <- shiny::renderPlot({
      shiny::req(current_sweep(),preview_row(),state()$preview)
      plot_power(state()$preview$result,lang=language())
    },res=110)
    shiny::observeEvent(input$apply,{
      tryCatch({
        sw <- current_sweep();s <- state()
        if(is.null(sw) || !ready() || is.null(s$preview)) .ig_abort("invalid_input","Preview a current successful row before applying it.")
        i <- row_number(sw)
        if(!identical(i,preview_row()) || !identical(s$preview$spec_hash,sw$spec_hash[i]))
          .ig_abort("invalid_input","The selected row has changed; preview it before applying.")
        state(.ig_workbench_apply_preview(s))
        apply_value(attr(sw,"free"),sw$value[i])
        sweep(NULL);preview_row(NULL);message(text("Preview applied to the design.","Vorschau in das Design uebernommen."))
      },error=function(e) message(conditionMessage(e)))
    },ignoreInit=TRUE)
    list(sweep=shiny::reactive(current_sweep()),preview_row=shiny::reactive(preview_row()))
  })
}
