# Input layouts only; all calculations and domain checks belong to effect_size().
.ig_power_effect_layouts <- function(method) {
  means <- list(common_sd = c("mean_1", "mean_2", "sd"),
                two_sds = c("mean_1", "mean_2", "sd_1", "sd_2"))
  anova <- list(variances = c("var_effect", "var_error"), eta_squared = "partial_eta_squared")
  switch(method,
    t_means_one_sample = list(from_means = c("mean_h0", "mean_h1", "sd")),
    t_means_two_independent = c(means, list(point_biserial_r = c("r", "n1", "n2"))),
    t_means_matched_pairs = list(from_differences = c("mean_diff", "sd_diff"),
      from_group_parameters = c("mean_1", "mean_2", "sd_1", "sd_2", "rho")),
    t_wilcoxon_signed_rank_one_sample = list(means_sd = c("mean_h0", "mean_h1", "sd")),
    t_wilcoxon_signed_rank_matched_pairs = list(differences = c("mean_diff", "sd_diff"),
      group_parameters = c("mean_1", "mean_2", "sd_1", "sd_2", "r")),
    t_wilcoxon_mann_whitney = means,
    f_variances_two_groups = list(from_variances = c("var0", "var1")),
    chi2_variance_one_sample = list(from_variances = c("var0", "var1")),
    f_anova_oneway = c(list(means = c("means", "sizes", "sd")), anova),
    f_anova_special_main_interactions = anova,
    f_ancova_main_interactions = anova,
    list())
}

.ig_power_effect_label <- function(name, lang) {
  labels <- list(
    from_means = c("Mean under H0 and H1", "Mittelwerte unter H0 und H1"),
    common_sd = c("Two means, common SD", "Zwei Mittelwerte, gemeinsame SD"),
    two_sds = c("Two means, two SDs", "Zwei Mittelwerte, zwei SD"),
    point_biserial_r = c("Point-biserial correlation", "Punktbiseriale Korrelation"),
    from_differences = c("Paired differences", "Gepaarte Differenzen"),
    from_group_parameters = c("Paired group parameters", "Parameter gepaarter Gruppen"),
    means_sd = c("Mean under H0 and H1", "Mittelwerte unter H0 und H1"),
    differences = c("Paired differences", "Gepaarte Differenzen"),
    group_parameters = c("Paired group parameters", "Parameter gepaarter Gruppen"),
    from_variances = c("Two variances", "Zwei Varianzen"),
    variances = c("Effect and error variances", "Effekt- und Fehlervarianz"),
    eta_squared = c("Partial eta squared", "Partielles Eta-Quadrat"),
    means = c("Group means", "Gruppenmittelwerte"),
    sizes = c("Group sample sizes", "Gruppenstichprobenumfaenge"),
    mean_h0 = c("Mean under H0", "Mittelwert unter H0"),
    mean_h1 = c("Mean under H1", "Mittelwert unter H1"),
    mean_1 = c("Mean, group 1", "Mittelwert, Gruppe 1"),
    mean_2 = c("Mean, group 2", "Mittelwert, Gruppe 2"),
    sd = c("Standard deviation", "Standardabweichung"),
    sd_1 = c("SD, group 1", "SD, Gruppe 1"),
    sd_2 = c("SD, group 2", "SD, Gruppe 2"),
    mean_diff = c("Mean paired difference", "Mittlere gepaarte Differenz"),
    sd_diff = c("SD of paired differences", "SD der gepaarten Differenzen"),
    r = c("Correlation r", "Korrelation r"),
    rho = c("Within-pair correlation", "Korrelation innerhalb der Paare"),
    n1 = c("Sample size, group 1", "Stichprobenumfang, Gruppe 1"),
    n2 = c("Sample size, group 2", "Stichprobenumfang, Gruppe 2"),
    var0 = c("Denominator variance (var0)", "Varianz im Nenner (var0)"),
    var1 = c("Numerator variance (var1)", "Varianz im Zaehler (var1)"),
    var_effect = c("Effect variance", "Effektvarianz"),
    var_error = c("Error variance", "Fehlervarianz"),
    partial_eta_squared = c("Partial eta squared", "Partielles Eta-Quadrat"))
  value <- labels[[name]]
  if (is.null(value)) return(name)
  value[if (lang == "de") 2L else 1L]
}

.ig_power_effect_numbers <- function(text, lang, vector = FALSE) {
  parts <- strsplit(paste0(text, ";END"), ";", fixed = TRUE)[[1]]
  parts <- utils::head(parts, -1L)
  if ((!vector && length(parts) != 1L) || !length(parts) || length(parts) > 100L ||
      any(!nzchar(trimws(parts)))) .ig_abort("invalid_input", "Enter finite numbers; separate vector entries with semicolons (maximum 100).")
  unname(vapply(parts, .ig_workbench_parse_number, numeric(1), locale = lang))
}

.ig_power_effect_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tags$details(shiny::tags$summary("Effect-size helper / Effektstaerkenhilfe"),
    shiny::uiOutput(ns("controls")), shiny::uiOutput(ns("fields")),
    shiny::uiOutput(ns("actions")), shiny::textOutput(ns("message")),
    shiny::uiOutput(ns("result")))
}

.ig_power_effect_server <- function(id, form, language, apply_value) {
  shiny::moduleServer(id, function(input, output, session) {
    result <- shiny::reactiveVal(NULL)
    message <- shiny::reactiveVal("")
    field_languages <- new.env(parent = emptyenv())
    layouts <- shiny::reactive(.ig_power_effect_layouts(form()$method))
    mode <- shiny::reactive({
      choices <- names(layouts())
      if (isTRUE(input$mode %in% choices)) input$mode else choices[1]
    })
    fields <- shiny::reactive({
      available <- layouts(); selected <- mode()
      if (!length(available) || length(selected) != 1L || is.na(selected) ||
          !selected %in% names(available)) return(character())
      available[[selected]] %||% character()
    })
    field_id <- function(nm) paste(form()$method, mode(), nm, sep = "_")
    output$controls <- shiny::renderUI({
      if (!length(layouts())) return(shiny::tags$p(if (language() == "de")
        "Fuer dieses Verfahren ist noch keine Eingabehilfe in der App verfuegbar." else
        "This method does not yet have an input helper in the app."))
      shiny::tagList(
        shiny::tags$p(.ig_get_metadata(form()$method)$effect_size$definition),
        shiny::selectInput(session$ns("mode"), if (language() == "de") "Ausgangsgroessen" else "Calculate from",
          stats::setNames(names(layouts()), vapply(names(layouts()), .ig_power_effect_label, "", lang = language())), selected = mode()),
        shiny::tags$p(if (language() == "de")
          "Originaleinheiten verwenden. Vektoren mit Semikolon trennen. Uebernehmen aendert nur den Entwurf; danach Berechnen waehlen." else
          "Use the original measurement units. Separate vector entries with semicolons. Apply changes only the draft; then choose Calculate."))
    })
    output$fields <- shiny::renderUI({
      lang <- language()
      shiny::tagList(lapply(fields(), function(nm) {
        id <- field_id(nm)
        current <- shiny::isolate(input[[id]]) %||% ""
        value <- tryCatch(.ig_power_effect_numbers(current, field_languages[[id]] %||% lang,
                         vector = nm %in% c("means", "sizes")), error = function(e) NULL)
        if (!is.null(value)) {
          current <- paste(vapply(value, .ig_format_double, ""), collapse = ";")
          if (lang == "de") current <- gsub(".", ",", current, fixed = TRUE)
        }
        field_languages[[id]] <- lang
        shiny::textInput(session$ns(id), .ig_power_effect_label(nm, lang), current)
      }))
    })
    values <- shiny::reactive(stats::setNames(lapply(fields(), function(nm)
      .ig_power_effect_numbers(input[[field_id(nm)]] %||% "", language(),
                              vector = nm %in% c("means", "sizes"))), fields()))
    current <- shiny::reactive({
      r <- result()
      if (is.null(r) || !identical(r$method_id, form()$method) || !identical(r$input_mode, mode())) return(NULL)
      v <- tryCatch(values(), error = function(e) NULL)
      if (!identical(v, r$entered)) return(NULL)
      r$value
    })
    output$actions <- shiny::renderUI({
      if (!length(layouts())) return(NULL)
      r <- current()
      shiny::tagList(shiny::actionButton(session$ns("calculate"), if (language() == "de") "Effekt berechnen" else "Calculate effect"),
        if (!is.null(r) && identical(r$status$code, "ok") && r$effect %in% names(form()$parameters))
          shiny::actionButton(session$ns("apply"), if (language() == "de") "In Entwurf uebernehmen" else "Apply to draft"))
    })
    shiny::observeEvent(input$calculate, {
      result(NULL)
      tryCatch({
        entered <- values()
        r <- effect_size(form()$method, mode(), inputs = entered)
        result(list(method_id = form()$method, input_mode = mode(), entered = entered, value = r))
        message("")
      }, error = function(e) message(conditionMessage(e)))
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$apply, {
      r <- current()
      if (is.null(r) || !identical(r$status$code, "ok") || !r$effect %in% names(form()$parameters)) {
        message("needs_effect")
        return()
      }
      apply_value(r$effect, r$value)
      message("applied")
    }, ignoreInit = TRUE)
    output$message <- shiny::renderText({
      switch(message(),
        needs_effect = if (language() == "de") "Zuerst einen aktuellen, anwendbaren Effekt berechnen." else
          "Calculate a current, applicable effect first.",
        applied = if (language() == "de") "Entwurf aktualisiert. Berechnen waehlen." else
          "Draft updated. Choose Calculate.", message())
    })
    output$result <- shiny::renderUI({
      r <- current()
      if (is.null(r)) return(NULL)
      if (!identical(r$status$code, "ok")) return(shiny::tags$p(
        role = "status", paste(r$status$code, r$status$reason)))
      shown <- .ig_format_double(r$value)
      if (language() == "de") shown <- sub(".", ",", shown, fixed = TRUE)
      shiny::tags$section(shiny::tags$p(paste(r$effect, "=", shown)),
        shiny::tags$p(r$formula), shiny::tags$p(r$source),
        shiny::tags$ul(lapply(r$warnings, shiny::tags$li)),
        if (!r$effect %in% names(form()$parameters)) shiny::tags$p(if (language() == "de")
          "Diese Groesse ist im aktuellen Modus kein Eingabefeld." else
          "This quantity is not an input in the current calculation mode."))
    })
    list(result = current)
  })
}
