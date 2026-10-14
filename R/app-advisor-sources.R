# Offline source views use the bundled audit, independently of advisor routing.
.ig_advisor_source_view <- function(r, lang, rules) {
  text <- function(en, de) if (lang == "de") de else en
  label <- function(n) n[[paste0("label_", lang)]] %||% n$label_en %||%
    n$source_label %||% n$code_test_name %||% n$id %||% n$node_id
  tree <- function(doc, selected, source) {
    nodes <- c(doc$nodes, doc$code_only_nodes)
    ids <- vapply(nodes, function(n) n$node_id %||% n$id, "")
    names(nodes) <- ids
    edges <- doc$edges
    linked <- selected
    repeat {
      parents <- vapply(Filter(function(e) e$to %in% linked, edges), `[[`, "", "from")
      next_linked <- union(linked, parents)
      if (identical(next_linked, linked)) break
      linked <- next_linked
    }
    node_ui <- function(id, path = character()) {
      n <- nodes[[id]]
      if (is.null(n) || id %in% path) return(NULL)
      children <- Filter(function(e) identical(e$from, id), edges)
      selected_label <- if (id %in% selected) text(" \u2014 linked source endpoint", " \u2014 verknuepfter Quellenendpunkt") else
        if (id %in% linked) text(" \u2014 source ancestor", " \u2014 uebergeordneter Quellenknoten") else ""
      content <- shiny::tagList(
        shiny::tags$p(.ig_advisor_link(text("Open source", "Quelle oeffnen"), n$source_url %||% doc$source$url %||% doc$source$root_url)),
        shiny::tags$p(n[[paste0("description_",lang)]] %||% n$description_en %||% n$description %||% n$reached_by %||% ""),
        if (length(n$source_rules)) shiny::tags$ul(lapply(n$source_rules, function(rule) shiny::tags$li(
          shiny::tags$p(shiny::tags$strong(text("Source rule: ", "Quellenregel: ")), rule$source_says_en),
          shiny::tags$p(shiny::tags$strong(text("Audit correction: ", "Audit-Korrektur: ")), rule$correction_en),
          shiny::tags$p(paste(text("Recorded status:", "Erfasster Status:"), rule$status))))),
        if (length(children)) shiny::tags$ul(lapply(children, function(e) shiny::tags$li(
          shiny::tags$p(if (lang == "de") e$condition_de %||% e$condition %||% e$diagram_label else e$condition %||% e$diagram_label),
          node_ui(e$to, c(path,id))))))
      shiny::tags$details(open = if (id %in% linked) "open" else NULL,
        `data-source-node` = id,
        shiny::tags$summary(label(n), selected_label), content)
    }
    targets <- vapply(edges, `[[`, "", "to")
    roots <- setdiff(ids, targets)
    shiny::tags$section(shiny::tags$h4(source),
      shiny::tags$p(paste(text("Bundled audit retrieved:", "Mitgeliefertes Audit abgerufen:"), paste(unlist(doc$retrieved_utc), collapse = "; "))),
      lapply(roots, node_ui))
  }
  cw <- rules$docs$crosswalk
  selected <- r$preferred$analysis_variant_id %||% ""
  original <- r$preferred$outcome_variant_id %||% ""
  keep <- (nzchar(cw$analysis_variant_id) & cw$analysis_variant_id %in% c(selected, original)) | cw$crosswalk_id %in% unlist(r$rules$crosswalk_ids)
  comparison <- cw[keep, , drop=FALSE]
  shiny::tagList(
    shiny::tags$details(class="ig-source-overview",
      shiny::tags$summary(text("Source tree overview", "Uebersicht der Quellenbaeume")),
      shiny::tags$p(text("Expand branches to inspect the bundled source audit. Highlighted ancestors lead to source endpoints linked to the original rule recommendation; they are not a replay of your answers. Unmarked branches have not been ruled out by this display.",
        "Zweige aufklappen, um das mitgelieferte Quellenaudit zu pruefen. Markierte Vorfahren fuehren zu Quellenendpunkten der urspruenglichen Regelempfehlung; sie bilden Ihre Antworten nicht nach. Unmarkierte Zweige werden durch diese Ansicht nicht ausgeschlossen.")),
      shiny::tags$p(text("UZH is a broad methods taxonomy; DescrTab2 concerns tests for summary tables. Audit annotations may be in English. Neither tree certifies numerical validity or planning compatibility.",
        "UZH ist eine breite Methodentaxonomie; DescrTab2 betrifft Tests fuer Uebersichtstabellen. Audit-Anmerkungen koennen englisch sein. Keiner der Baeume bestaetigt numerische Gueltigkeit oder Planungskompatibilitaet.")),
      tree(rules$docs$uzh_tree, unlist(r$rules$source_nodes$uzh), "UZH"),
      tree(rules$docs$descrtab2_tree, unlist(r$rules$source_nodes$descrtab2), "DescrTab2")),
    shiny::tags$details(class="ig-source-comparison",
      shiny::tags$summary(text("Compare source and application guidance", "Quellen- und Anwendungsempfehlungen vergleichen")),
      shiny::tags$p(text("Differences below are recorded audit judgments, not a new source review. Method identity does not imply identical routing or assumptions.",
        "Die folgenden Unterschiede sind dokumentierte Audit-Urteile, keine neue Quellenpruefung. Gleiche Methoden bedeuten nicht gleiche Auswahlregeln oder Annahmen.")),
      if (!nrow(comparison)) shiny::tags$p(text("No crosswalk entry is recorded for this outcome.", "Fuer dieses Ergebnis ist kein Vergleichseintrag erfasst.")) else
        lapply(seq_len(nrow(comparison)), function(i) {
          row <- comparison[i,,drop=FALSE]
          variant <- rules$variants[[row$analysis_variant_id]]
          name <- if (nzchar(row$analysis_variant_id)) variant[[paste0("name_",lang)]] %||% row$analysis_variant_id else
            paste(text("Outcome or routing context", "Ergebnis- oder Auswahlkontext"), row$crosswalk_id)
          shiny::tags$section(`data-crosswalk-id`=row$crosswalk_id,
            shiny::tags$h4(name),
            shiny::tags$p(if (!nzchar(row$analysis_variant_id)) text("No method is specified by this entry", "Dieser Eintrag benennt keine Methode") else if (identical(row$analysis_variant_id,selected)) text("Selected method", "Gewaehlte Methode") else text("Original rule recommendation", "Urspruengliche Regelempfehlung")),
            shiny::tags$p(paste("UZH:", if (nzchar(row$uzh_node_id)) row$uzh_node_id else text("No mapped endpoint", "Kein zugeordneter Endpunkt"))),
            shiny::tags$p(paste("DescrTab2:", if (nzchar(row$descrtab2_node_id)) row$descrtab2_node_id else text("No mapped endpoint", "Kein zugeordneter Endpunkt"))),
            shiny::tags$p(paste(text("Recorded relationship:", "Erfasste Beziehung:"), row$relationship)),
            shiny::tags$p(row$difference_reason), shiny::tags$p(row$notes))
        })))
}
