.ig_power_workspace_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tags$div(class = "ig-study",
    shiny::tags$aside(class = "ig-controls",
      shiny::uiOutput(ns("method_selector")),
      shiny::checkboxInput(ns("research"), "Enable research preview / Forschungsvorschau aktivieren", FALSE),
      shiny::selectInput(ns("mode"), "Calculate / Berechnen", choices = c(
        "Sample size / Stichprobengroesse" = "a_priori", "Power / Teststaerke" = "post_hoc",
        "Effect size / Effektstaerke" = "sensitivity", "Alpha" = "criterion", "Error ratio / Fehlerverhaeltnis" = "compromise"), selected = "post_hoc"),
      shiny::uiOutput(ns("options")),
      shiny::uiOutput(ns("target")),
      shiny::uiOutput(ns("context")),
      shiny::uiOutput(ns("fields")),
      .ig_power_effect_ui(ns("effect")),
      shiny::actionButton(ns("calculate"), "Calculate / Berechnen")),
    shiny::tags$main(class = "ig-figure",
      shiny::tags$div(role = "status", `aria-live` = "polite", shiny::textOutput(ns("status"))),
      shiny::uiOutput(ns("result_method")),
      shiny::plotOutput(ns("distribution"), height = "620px"),
      shiny::tableOutput(ns("decisions")),
      shiny::uiOutput(ns("history_controls")),
      .ig_power_sweep_ui(ns("sweep")),
      .ig_report_ui(ns("report"))),
    shiny::tags$aside(class = "ig-inspector", shiny::uiOutput(ns("assumptions"))))
}

.ig_power_workspace_server <- function(id, language) {
  shiny::moduleServer(id, function(input, output, session) {
    state <- shiny::reactiveVal(NULL)
    restoring <- shiny::reactiveVal(NULL)
    message <- shiny::reactiveVal("")
    context_notice <- shiny::reactiveVal("")
    field_languages <- new.env(parent = emptyenv())
    field_drafts <- new.env(parent = emptyenv())
    remember_field <- function(id, text, lang) {
      saved <- field_drafts[[id]]
      if (!is.null(text) && (is.null(saved) || !identical(text, saved$text))) {
        saved <- list(text = text, locale = field_languages[[id]] %||% lang)
        field_drafts[[id]] <- saved
      }
      saved
    }
    example <- list(d = .5, alpha = .05, power = .8, n_total = 30, error_ratio = 1)
    method <- shiny::reactive(input$method %||% "t_means_one_sample")
    definition <- shiny::reactive(.ig_get_method(method()))
    output$method_selector <- shiny::renderUI({
      lang <- language()
      methods <- list_methods()
      labels <- vapply(methods$id, function(id) .ig_workbench_label(.ig_get_method(id)$label, lang, id), "")
      choices <- split(stats::setNames(methods$id, labels), methods$family)
      shiny::selectInput(session$ns("method"), if (lang == "de") "Verfahren" else "Method",
                         choices, selected = restoring()$method_id %||% shiny::isolate(method()))
    })
    shiny::observeEvent(method(), {
      message("")
      modes <- definition()$modes
      labels <- c(a_priori = "Sample size / Stichprobenumfang", post_hoc = "Power / Teststaerke",
                  sensitivity = "Effect size / Effektstaerke", criterion = "Alpha", compromise = "Error balance / Fehlerbalance")
      selected <- if (!is.null(restoring())) restoring()$mode else if (!is.null(input$mode) && input$mode %in% modes) input$mode else
        if ("post_hoc" %in% modes) "post_hoc" else modes[1]
      shiny::updateSelectInput(session, "mode", choices = stats::setNames(modes, labels[modes]), selected = selected)
      context_notice(if (language() == "de") "Methodenbezogene Eingaben werden getrennt aufbewahrt. Neue Felder bitte pruefen." else
        "Inputs are kept separately for each method. Review the fields after changing method.")
    })
    input_form <- shiny::reactive({
      def <- definition()
      mode <- input$mode
      if (is.null(mode) || !mode %in% def$modes) mode <- if ("post_hoc" %in% def$modes) "post_hoc" else def$modes[1]
      options <- list()
      for (nm in names(def$options)) {
        value <- input[[.ig_workbench_field_id(method(), nm, "option")]]
        if (!is.null(value) && nzchar(value)) options[[nm]] <- value
      }
      if (method() == "t_means_one_sample" && is.null(options$tails)) options$tails <- "two"
      target <- input[[.ig_workbench_field_id(method(), mode, "target")]]
      if (identical(target, "")) target <- NULL
      .ig_workbench_form(method(), mode, options, target)
    })
    form <- shiny::reactive({
      saved <- restoring()
      if (is.null(saved)) input_form() else
        .ig_workbench_form(saved$method_id, saved$mode, saved$options, saved$target)
    })
    field_layout <- shiny::reactiveVal(NULL)
    option_layout <- shiny::reactiveVal(NULL)
    target_layout <- shiny::reactiveVal(NULL)
    restore_active <- FALSE
    restore_generation <- 0L
    shiny::observe({
      f <- form(); active <- !is.null(restoring())
      if (active && !restore_active) restore_generation <<- restore_generation + 1L
      restore_active <<- active
      layout <- list(method = f$method, parameters = f$parameters,
                     restore_generation = restore_generation)
      # Option bindings often invalidate form() without changing these fields.
      # Replacing an identical editor can erase typing before its input event
      # reaches the server. History explicitly requests one new render; ending
      # restoration does not request another replacement of the same editor.
      if (!identical(shiny::isolate(field_layout()), layout)) field_layout(layout)
      options <- list(method = f$method, choices = f$option_choices,
                      restore_generation = restore_generation)
      if (!identical(shiny::isolate(option_layout()), options)) option_layout(options)
      target <- list(method = f$method, mode = f$mode, targets = f$targets,
                     restore_generation = restore_generation)
      if (!identical(shiny::isolate(target_layout()), target)) target_layout(target)
    }, priority = 90)
    # Dynamic editors can briefly unbind inputs during method/locale changes.
    # Retain the latest draft, including invalid and deliberately empty text,
    # with the locale in which it was entered; never substitute an old result.
    shiny::observe({
      f <- input_form()
      for (p in f$parameters) {
        id <- .ig_workbench_field_id(f$method, p$name)
        remember_field(id, input[[id]], language())
      }
    }, priority = 100)
    shiny::observeEvent(form(), {
      f <- form()
      reset <- union(f$dropped_options, names(f$reset_options))
      if (length(reset) || f$reset_target) {
        labels <- vapply(reset, function(nm) .ig_workbench_label(f$definition$options[[nm]]$label, language(), nm), "")
        context_notice(paste(if (language() == "de") "Durch Modellwechsel zurueckgesetzt oder ausgeblendet:" else
          "Reset or hidden after the model change:", paste(c(labels, if (f$reset_target) "solve target"), collapse = ", ")))
      }
    }, ignoreInit = TRUE)
    output$options <- shiny::renderUI({
      shiny::req(!is.null(option_layout()))
      f <- shiny::isolate(form())
      shiny::tagList(lapply(names(f$option_choices), function(nm) {
        o <- f$definition$options[[nm]]
        codes <- f$option_choices[[nm]]
        labels <- o$labels[codes]
        if (identical(f$method, .bl_method_id) && identical(nm, "variant"))
          labels <- .ig_boschloo_variant_labels(language())[codes]
        if (is.null(labels)) labels <- codes
        labels[is.na(labels) | !nzchar(labels)] <- codes[is.na(labels) | !nzchar(labels)]
        shiny::selectInput(session$ns(.ig_workbench_field_id(f$method, nm, "option")),
          .ig_workbench_label(o$label, language(), nm), c(stats::setNames("", if (language() == "de") "Bitte waehlen" else "Choose"),
                                  stats::setNames(codes, labels)), selected = f$options[[nm]] %||% "")
      }))
    })
    output$target <- shiny::renderUI({
      shiny::req(!is.null(target_layout()))
      f <- shiny::isolate(form())
      label <- function(nm) .ig_workbench_label(f$definition$params[[nm]]$label, language(), nm)
      if (length(f$targets) > 1L && f$mode != "compromise") {
        shiny::selectInput(session$ns(.ig_workbench_field_id(f$method, f$mode, "target")),
          if (language() == "de") "Gesuchte Groesse" else "Solve for",
          stats::setNames(f$targets, vapply(f$targets, label, "")), selected = f$target)
      } else shiny::tags$p(shiny::tags$strong(if (language() == "de") "Gesucht: " else "Solving for: "),
                            paste(vapply(f$target, label, ""), collapse = ", "))
    })
    output$context <- shiny::renderUI({
      shiny::tags$div(shiny::tags$p(context_notice()),
        if (!form()$ready) shiny::tags$ul(lapply(form()$problems, shiny::tags$li)))
    })
    output$fields <- shiny::renderUI({
      lang <- language()
      f <- field_layout(); shiny::req(!is.null(f))
      restore <- shiny::isolate(restoring())
      shiny::tagList(lapply(f$parameters, function(p) {
        id <- .ig_workbench_field_id(f$method, p$name)
        current <- shiny::isolate(input[[id]])
        saved <- remember_field(id, current, lang)
        if (is.null(current) && !is.null(saved)) current <- saved$text
        value <- if (!is.null(restore)) restore$inputs[[p$name]] else if (is.null(current)) p$default %||% (if (f$method == "t_means_one_sample") example[[p$name]]) else
          tryCatch(.ig_workbench_parse_number(current, saved$locale %||% field_languages[[id]] %||% lang), error = function(e) NULL)
        shown <- if (is.null(value)) current %||% "" else .ig_format_double(value)
        if (!is.null(value) && lang == "de") shown <- sub(".", ",", shown, fixed = TRUE)
        field_languages[[id]] <- lang
        shiny::tagList(shiny::textInput(session$ns(id), .ig_workbench_label(p$label, lang, p$name), value = shown),
          shiny::tags$small(.ig_workbench_domain_hint(p, lang)))
      }))
    })
    draft <- shiny::reactive({
      f <- input_form()
      list(form = f, fields = stats::setNames(lapply(names(f$parameters), function(nm)
        input[[.ig_workbench_field_id(f$method, nm)]]), names(f$parameters)))
    })
    make_spec <- function(d) .ig_workbench_form_spec(d$form, d$fields, language())
    shiny::observeEvent(list(draft(), restoring()), {
      s <- state()
      if (is.null(s)) return()
      candidate <- tryCatch(make_spec(draft()), error = function(e) NULL)
      s$dirty <- is.null(candidate) || !identical(candidate$spec_hash, s$committed$spec_hash)
      s$draft <- if (s$dirty) draft() else list()
      state(s)
      if (!s$dirty && !is.null(restoring()) &&
          identical(input$method, restoring()$method_id) &&
          identical(input$mode, restoring()$mode)) restoring(NULL)
    }, ignoreInit = TRUE)
    restore_history <- function(direction) {
      s <- state()
      if (is.null(s) || !is.null(restoring())) return(invisible(NULL))
      next_state <- .ig_workbench_history(s, direction)
      if (identical(s$cursor, next_state$cursor)) return(invisible(NULL))
      # The result is immutable; the browser must acknowledge the matching inputs.
      next_state$dirty <- TRUE
      next_state$draft <- draft()
      state(next_state)
      restoring(next_state$committed$spec)
      message("")
      shiny::updateSelectInput(session, "mode", selected = next_state$committed$spec$mode)
    }
    shiny::observeEvent(input$history_back, restore_history(-1), ignoreInit = TRUE)
    shiny::observeEvent(input$history_forward, restore_history(1), ignoreInit = TRUE)
    output$history_controls <- shiny::renderUI({
      s <- state()
      if (is.null(s)) return(NULL)
      disabled <- !is.null(restoring())
      button <- function(id, label, unavailable) {
        b <- shiny::actionButton(session$ns(id), label)
        if (unavailable) b$attribs$disabled <- "disabled"
        b
      }
      shiny::tags$section(
        shiny::tags$p(if (language() == "de") "Verlauf" else "History", ": ",
                      s$cursor, " / ", length(s$history)),
        button("history_back", if (language() == "de") "Zurueck" else "Previous",
               disabled || s$cursor <= 1L),
        button("history_forward", if (language() == "de") "Weiter" else "Next",
               disabled || s$cursor >= length(s$history)),
        shiny::tags$p(if (language() == "de")
          "Verlaufsauswahl ersetzt den aktuellen Eingabeentwurf." else
          "Selecting history replaces the current input draft."))
    })
    shiny::observeEvent(language(), {
      if (!is.null(state())) state(.ig_workbench_locale(state(), language()))
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$calculate, {
      if (!is.null(restoring())) return()
      if (!isTRUE(input$research)) {
        message("This method is a research preview; its numerical release gate is still pending.")
        return()
      }
      tryCatch({
        spec <- make_spec(draft())
        s <- state()
        s <- if (is.null(s)) .ig_workbench_new(spec, session$token, locale = language()) else
          .ig_workbench_commit(s, spec)
        ticket <- s$pending
        state(s)
        result <- calculate_power(spec)
        state(.ig_workbench_accept(state(), ticket, result))
        message("")
        context_notice("")
      }, error = function(e) message(conditionMessage(e)))
    }, ignoreInit = TRUE)
    output$status <- shiny::renderText({
      s <- state()
      if (nzchar(message())) return(paste(message(),
        if (!is.null(s$committed)) if (language() == "de")
          "Angezeigt wird das vorherige Ergebnis." else "The previous result remains displayed."))
      if (is.null(s)) return(if (language() == "de") "Einstellungen pr\u00fcfen und Berechnen w\u00e4hlen." else "Review the settings, then choose Calculate.")
      if (!is.null(restoring())) return(if (language() == "de")
        "Verlaufseingaben werden wiederhergestellt. Ergebnis noch nicht mit Eingaben abgeglichen." else
        "Restoring history inputs. Result and inputs have not yet been matched.")
      if (s$dirty) return(if (language() == "de") "Entwurf ge\u00e4ndert. Angezeigt wird das vorherige Ergebnis." else "Draft changed. The previous result remains displayed.")
      paste(if (language() == "de") "Berechnungsstatus:" else "Calculation status:", s$committed$result$status$code)
    })
    current <- shiny::reactive({
      shiny::req(state())
      state()$committed$result
    })
    output$result_method <- shiny::renderUI({
      r <- current()
      shiny::tags$h3(.ig_workbench_label(.ig_get_method(r$method_id)$label, language(), r$method_id))
    })
    output$distribution <- shiny::renderPlot({
      r <- current()
      panel <- if (identical(r$method_id, .bl_method_id)) "lattice" else "distribution"
      pd <- power_plot_data(r, panel)
      shiny::validate(shiny::need(pd$status == "ok", pd$reason %||% "Distribution unavailable"))
      .ig_workbench_plot_labels(plot_power(r, panel = panel, lang = language()),
        session$clientData[[paste0("output_",session$ns("distribution"),"_width")]])
    }, res = 110)
    output$decisions <- shiny::renderTable({
      pd <- power_plot_data(current(), "decision_table")
      shiny::validate(shiny::need(pd$status == "ok", pd$reason %||% "Decision table unavailable"))
      tab <- pd$data[, c("hypothesis", "decision", "probability")]
      tab$probability <- formatC(tab$probability, digits = 6, format = "g",
                                decimal.mark = if (language() == "de") "," else ".")
      tab$decision <- if (language() == "de") ifelse(tab$decision == "reject", "H0 verwerfen", "H0 nicht verwerfen") else
        ifelse(tab$decision == "reject", "Reject H0", "Do not reject H0")
      names(tab) <- if (language() == "de") c("Hypothese", "Entscheidung", "Bedingte Wahrscheinlichkeit") else
        c("Hypothesis", "Decision", "Conditional probability")
      tab
    }, digits = 6, striped = FALSE)
    output$assumptions <- shiny::renderUI({
      r <- current()
      shiny::tags$section(shiny::tags$h3(if (language() == "de") "Annahmen" else "Assumptions"),
        shiny::tags$ul(lapply(.ig_power_display_assumptions(r, language()), shiny::tags$li)),
        shiny::tags$p(if (language() == "de") "Die H1-Spalte gilt nur f\u00fcr den angegebenen Effekt und Entwurf." else
          "The H1 column applies only to the specified effect and design."))
    })
    effect_module <- .ig_power_effect_server("effect", form, language,
      apply_value = function(field, value) {
        if (!is.null(restoring())) return(invisible(NULL))
        shown <- .ig_format_double(value)
        if (language() == "de") shown <- sub(".", ",", shown, fixed = TRUE)
        shiny::updateTextInput(session, .ig_workbench_field_id(form()$method, field), value = shown)
      })
    sweep_module <- .ig_power_sweep_server("sweep", state, language,
      research = shiny::reactive(isTRUE(input$research)),
      apply_value = function(free, value) {
        shown <- .ig_format_double(value)
        if (language() == "de") shown <- sub(".", ",", shown, fixed = TRUE)
        shiny::updateTextInput(session, .ig_workbench_field_id(state()$committed$spec$method_id, free), value = shown)
      })
    report_module <- .ig_report_server("report", shiny::reactive(state()$committed), language,
      eligible = shiny::reactive(!is.null(state()) && !state()$dirty && is.null(restoring())),
      sweeps = shiny::reactive({sw <- sweep_module$sweep();if(is.null(sw)) list() else list(list(analysis=1,sweep=sw))}))
    list(state = shiny::reactive(state()), result = current, sweep = sweep_module, effect = effect_module, report = report_module)
  })
}
