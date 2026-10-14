.ig_advisor_link <- function(label, url) {
  if (.ig_is_string(url) && grepl("^https?://[^[:space:]]+$", url))
    shiny::tags$a(label, href = url, target = "_blank", rel = "noopener noreferrer") else
    shiny::tags$span(label)
}

.ig_advisor_card <- function(r, lang, rules, question_labels = list()) {
  text <- function(en, de) if (lang == "de") de else en
  field <- function(x, prefix) x[[paste0(prefix, "_", lang)]] %||% x[[paste0(prefix, "_en")]] %||% ""
  section <- function(title, ...) shiny::tags$section(shiny::tags$h4(title), ...)
  status <- switch(r$status, recommended = text("Recommended for the recorded design", "Empfehlung fuer das erfasste Design"),
    conditional = text("Conditional recommendation", "Bedingte Empfehlung"),
    needs_information = text("More information needed", "Weitere Informationen erforderlich"),
    outside_supported_scope = text("Outside the supported scope", "Ausserhalb des unterstuetzten Bereichs"))
  overridden <- isTRUE(r$preferred$overridden)
  variant_name <- function(id) {
    if (is.null(id)) return("")
    variant <- rules$variants[[id]]
    if (is.null(variant)) id else field(variant, "name")
  }
  assumptions <- lapply(r$assumptions, function(a) {
    context <- NULL
    if (overridden) {
      context <- if (identical(a$source, "override")) {
        text("Reason for considering the selected alternative", "Begruendung der gewaehlten Alternative")
      } else if (!is.null(a$stated_for) && identical(a$stated_for, r$preferred$outcome_variant_id)) {
        paste0(text("Original recommendation: ", "Urspruengliche Empfehlung: "), variant_name(a$stated_for))
      } else if (!is.null(a$stated_for) && identical(a$stated_for, r$preferred$analysis_variant_id)) {
        paste0(text("Selected alternative: ", "Gewaehlte Alternative: "), variant_name(a$stated_for))
      } else if (!is.null(a$stated_for)) {
        paste0(text("Recorded for: ", "Erfasst fuer: "), variant_name(a$stated_for))
      } else text("Recorded design context", "Erfasster Designkontext")
    }
    shiny::tags$li(if (!is.null(context)) shiny::tags$strong(paste0(context, ": ")),
      field(a, "text"), " (", if (identical(a$state, "unknown")) text("not assessed", "nicht geprueft") else a$state, ")")
  })
  alternatives <- lapply(r$alternatives, function(a) shiny::tags$li(
    shiny::tags$strong(field(a, "name")), ": ", field(a, "why_not_preferred")))
  learning <- r$preferred$learning
  learning_ui <- lapply(c("course", "tutorial", "overview"), function(role) {
    block <- learning[[role]]
    if (is.null(block)) return(NULL)
    title <- switch(role, course = text("Course", "Kurs"), tutorial = text("Method tutorial", "Methodentutorial"),
      overview = text("Broader resources", "Weiterfuehrende Ressourcen"))
    section(title, shiny::tags$p(field(block, "note")), shiny::tags$ul(lapply(block$links, function(l) {
      topics <- vapply(l$method_matches, function(m) {
        variant <- rules$variants[[m$method_id]]
        topic <- if (!is.null(variant)) field(variant, "name") else gsub("_", " ", m$method_id)
        match <- switch(m$match, exact = text("exact topic match", "exakte Themenzuordnung"),
          related = text("related", "verwandt"), broader = text("broader", "allgemeiner"))
        paste0(topic, " (", match, ")")
      }, "")
      shiny::tags$li(.ig_advisor_link(l$title, l$published_url),
        shiny::tags$p(paste(text("Linked topics:", "Verknuepfte Themen:"), paste(topics, collapse = "; "))),
        if (isTRUE(l$caution)) shiny::tags$details(shiny::tags$summary(text("Recorded limitations", "Erfasste Einschraenkungen")),
          shiny::tags$ul(lapply(l$known_issues, shiny::tags$li))))
    })))
  })
  sources <- list()
  for (id in unlist(r$rules$source_nodes$uzh)) {
    nodes <- Filter(function(n) identical(n$node_id, id), rules$docs$uzh_tree$nodes)
    if (length(nodes)) sources[[length(sources) + 1L]] <- shiny::tags$li(
      .ig_advisor_link(paste("UZH:", field(nodes[[1]], "label")), nodes[[1]]$source_url))
  }
  for (id in unlist(r$rules$source_nodes$descrtab2)) {
    nodes <- Filter(function(n) identical(n$id, id), c(rules$docs$descrtab2_tree$nodes, rules$docs$descrtab2_tree$code_only_nodes))
    if (length(nodes)) sources[[length(sources) + 1L]] <- shiny::tags$li(
      .ig_advisor_link(paste("DescrTab2:", nodes[[1]]$source_label %||% id), rules$docs$descrtab2_tree$source$url))
  }
  name <- field(r$preferred, "name")
  shiny::tags$section(class = "ig-recommendation",
    shiny::tags$h3(if (nzchar(name)) name else status),
    shiny::tags$p(shiny::tags$strong(status)),
    shiny::tags$p(r$rationale[[lang]] %||% r$rationale$en),
    if (!is.null(r$estimand)) section(text("What is estimated", "Zielgroesse"), shiny::tags$p(r$estimand[[lang]] %||% r$estimand$en)),
    if (length(r$missing_information)) section(text("Information still needed", "Noch benoetigte Informationen"),
      shiny::tags$ul(lapply(r$missing_information, function(id) shiny::tags$li(question_labels[[id]] %||% id)))),
    section(text("Assumptions to assess", "Zu pruefende Annahmen"),
      if (length(assumptions)) shiny::tags$ul(assumptions) else shiny::tags$p(text(
        "No method-specific assumptions are recorded for this outcome.", "Fuer dieses Ergebnis sind keine methodenspezifischen Annahmen erfasst."))),
    if (length(alternatives)) section(text("Alternatives and tradeoffs", "Alternativen und Abwaegungen"), shiny::tags$ul(alternatives)),
    section(text("R analysis", "R-Analyse"), shiny::tags$p(field(r$r_implementation, "note"))),
    section(text("Sample-size planning", "Fallzahlplanung"), shiny::tags$p(field(r$planning, "note"))),
    shiny::tags$p(text("Topic matches can concern supporting concepts; they do not by themselves mean a resource teaches the selected test variant.",
      "Themenzuordnungen koennen ergaenzende Konzepte betreffen; sie bedeuten nicht automatisch, dass eine Ressource die gewaehlte Testvariante behandelt.")),
    learning_ui,
    if (!is.null(learning)) shiny::tags$p(field(learning, "evidence_note")),
    section(if (overridden) text("Sources of the original rule recommendation", "Quellen der urspruenglichen Regelempfehlung") else
      text("Methodological sources", "Methodische Quellen"),
      if (overridden) shiny::tags$p(paste0(
        text("These source links trace the original recommendation: ", "Diese Quellenlinks dokumentieren die urspruengliche Empfehlung: "),
        variant_name(r$preferred$outcome_variant_id),
        text(". They do not establish support for the selected alternative.", ". Sie belegen nicht die gewaehlte Alternative."))),
      shiny::tags$p(text("External methodological background; executable analysis uses R. Source records belong to the bundled audit version.",
        "Externer methodischer Hintergrund; ausfuehrbare Analysen verwenden R. Quellenangaben beziehen sich auf die mitgelieferte Audit-Version.")),
      shiny::tags$ul(sources),
      if (length(r$references)) shiny::tags$details(shiny::tags$summary(text("Supporting references", "Weiterfuehrende Literatur")),
        shiny::tags$ul(lapply(r$references, function(id) {
          ref <- rules$references[[id]]
          shiny::tags$li(ref$citation %||% id, shiny::tags$p(ref$verification %||% ""))
        })))),
    shiny::tags$details(shiny::tags$summary(text("Rule provenance", "Regelprovenienz")),
      shiny::tags$p(paste(text("Rule version:", "Regelversion:"), r$rules$rules_version)),
      shiny::tags$p(paste(text("Rule status:", "Regelstatus:"), r$rules$rules_status)),
      shiny::tags$p(paste("UZH:", r$rules$uzh_tree_version)),
      shiny::tags$p(paste("DescrTab2:", r$rules$descrtab2_tree_version))),
    shiny::tags$p(text("Research preview: assumptions have not been verified against your data.",
      "Forschungsvorschau: Annahmen wurden nicht anhand Ihrer Daten geprueft.")))
}
