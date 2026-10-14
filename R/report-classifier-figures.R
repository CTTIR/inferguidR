# Static classifier views derived only from the retained aggregate cells.
.ig_report_classifier_aggregate_figures <- function(m, analysis, locale) {
  if (is.null(m$operating_point)) return(list())
  de <- locale == "de"
  cells <- m$operating_point[c("tp", "fn", "fp", "tn")]
  figures <- list()
  mosaic <- function(cells, basis, kind) {
    cts <- unlist(cells, use.names = TRUE)
    p <- cts[["tp"]] + cts[["fn"]]; n <- cts[["fp"]] + cts[["tn"]]
    share <- p / (p + n)
    rectangles <- list()
    add <- function(cell, x0, x1, y0, y1) {
      rectangles[[length(rectangles) + 1L]] <<- list(cell = cell, x0 = x0, x1 = x1,
        y0 = y0, y1 = y1, count = cts[[cell]], share = cts[[cell]] / (p + n))
    }
    if (p > 0) {rate <- cts[["tp"]] / p; add("tp", 0, share, 0, rate); add("fn", 0, share, rate, 1)}
    if (n > 0) {rate <- cts[["fp"]] / n; add("fp", share, 1, 0, rate); add("tn", share, 1, rate, 1)}
    list(analysis = analysis, model_id = m$model_id, plot_type = "classifier_mosaic", basis = basis,
      status = "ok", axes = list(x = if (de) "Ereignis links; kein Ereignis rechts" else "Event left; non-event right",
        y = if (de) "Anteil innerhalb der Klasse" else "Share within class"),
      data = list(rectangles = rectangles, count_kind = kind),
      caption = if (de) paste("Flaechen entsprechen Anteilen an allen Einheiten. TP und FP bilden gemeinsam alle Alarme.",
        if (basis == "scenario") "Erwartete Szenariowerte bei festen klassenbedingten Raten; keine neuen Beobachtungen." else "Zellen des gewaehlten Schwellenwerts.") else
        paste("Areas represent shares of all units. TP and FP together form all alerts.",
        if (basis == "scenario") "Expected scenario values at fixed class-conditional rates; not new observations." else "Cells at the selected cutoff."))
  }
  figures[[1L]] <- mosaic(cells, if (m$count_kind == "expected") "theoretical" else "observed", m$count_kind)
  if (!is.null(m$scenario) && identical(m$scenario$status, "ok"))
    figures[[length(figures) + 1L]] <- mosaic(m$scenario$counts[c("tp", "fn", "fp", "tn")], "scenario", "expected")
  cts <- unlist(cells, use.names = TRUE)
  p <- cts[["tp"]] + cts[["fn"]]; n <- cts[["fp"]] + cts[["tn"]]
  defined <- p > 0 && n > 0
  curve <- list(); selected <- list()
  if (defined) {
    tpr <- cts[["tp"]] / p; fpr <- cts[["fp"]] / n
    point <- function(pi, label = NULL) {
      alert <- pi * tpr + (1 - pi) * fpr
      quiet <- pi * (1 - tpr) + (1 - pi) * (1 - fpr)
      out <- list(prevalence = pi, ppv = if (alert > 0) pi * tpr / alert else NULL,
        npv = if (quiet > 0) (1 - pi) * (1 - fpr) / quiet else NULL)
      if (!is.null(label)) out$label <- label
      out
    }
    curve <- lapply(seq(0, 1, length.out = 101L), point)
    selected[[1L]] <- point(p / (p + n), if (m$count_kind == "expected") "assumed" else "observed")
    if (!is.null(m$scenario) && identical(m$scenario$status, "ok")) selected[[2L]] <- point(m$scenario$prevalence, "scenario")
  }
  figures[[length(figures) + 1L]] <- list(analysis = analysis, model_id = m$model_id,
    plot_type = "classifier_prevalence", basis = "scenario", status = if (defined) "ok" else "not_applicable",
    reason = if (defined) NULL else if (de) "Beide Klassen werden zur Schaetzung der bedingten Raten benoetigt." else "Both classes are needed to estimate the conditional rates.",
    axes = list(x = if (de) "Angenommene Praevalenz" else "Assumed prevalence", y = if (de) "Vorhersagewert" else "Predictive value"),
    data = list(prevalence_grid = curve, selected = selected),
    caption = if (de) "PPV und NPV bei festem Schwellenwert und festen klassenbedingten Raten. Die Kurven zeigen Szenarien, keine neuen Beobachtungen. Undefinierte Nenner bleiben offen; Unsicherheit der bedingten Raten ist nicht dargestellt." else
      "PPV and NPV at a fixed cutoff and fixed class-conditional rates. Curves describe scenarios, not new observations. Undefined denominators remain gaps; sampling uncertainty of the conditional rates is not shown.")
  figures
}

.ig_report_draw_classifier_aggregate <- function(v, locale) {
  de <- locale == "de"
  if (v$plot_type == "classifier_mosaic") {
    graphics::plot(c(0, 1), c(0, 1), type = "n", xaxs = "i", yaxs = "i",
      xlab = "", ylab = v$axes$y,
      main = if (de) "Gemeinsame Anteile der vier Zellen" else "Joint shares of the four cells")
    colours <- c(tp = "#72BCB3", fn = "#C8E8E4", fp = "#D5BA85", tn = "#ACBAE1")
    for (r in v$data$rectangles) {
      graphics::rect(r$x0, r$y0, r$x1, r$y1, col = colours[[r$cell]], border = "white")
      if ((r$x1-r$x0) > .12 && (r$y1-r$y0) > .08)
        graphics::text((r$x0+r$x1)/2, (r$y0+r$y1)/2, toupper(r$cell))
    }
    graphics::mtext(v$axes$x, side = 1, line = 2)
    graphics::legend(x = .5, y = -.075, xjust = .5, yjust = 1, xpd = NA, horiz = TRUE,
      legend = toupper(names(colours)), fill = colours, bty = "n", cex = .85)
  } else {
    graphics::plot(c(0, 1), c(0, 1), type = "n", xlab = v$axes$x, ylab = v$axes$y,
      main = if (de) "Vorhersagewerte und Praevalenz" else "Predictive values and prevalence")
    x <- vapply(v$data$prevalence_grid, `[[`, 0, "prevalence")
    for (k in seq_along(c("ppv", "npv"))) {
      metric <- c("ppv", "npv")[k]; colour <- c("#007C78", "#334E9A")[k]
      y <- vapply(v$data$prevalence_grid, function(r) if (is.null(r[[metric]])) NA_real_ else r[[metric]], 0)
      graphics::lines(x, y, col = colour, lty = k, lwd = 2)
      for (j in seq_along(v$data$selected)) {
        row <- v$data$selected[[j]]
        if (!is.null(row[[metric]])) graphics::points(row$prevalence, row[[metric]], pch = if(j==1L) 16 else 17, col = colour, cex = 1.2)
      }
    }
    marker_labels <- vapply(v$data$selected, function(row) switch(row$label,
      observed = if (de) "Beobachtete Praevalenz" else "Observed prevalence",
      assumed = if (de) "Angenommene Praevalenz" else "Assumed prevalence",
      scenario = if (de) "Szenario" else "Scenario"), "")
    graphics::legend("bottom", legend = c("PPV", "NPV", marker_labels),
      col = c("#007C78", "#334E9A", rep("#444444", length(marker_labels))),
      lty = c(1:2, rep(NA, length(marker_labels))), pch = c(NA, NA, c(16, 17)[seq_along(marker_labels)]),
      lwd = 2, bty = "n", cex = .8)
  }
  basis_label <- switch(v$basis, observed = if (de) "Ausgangsdaten" else "Original data",
    theoretical = if (de) "Theoretische Erwartungswerte" else "Theoretical expected values",
    scenario = if (de) "Szenario" else "Scenario")
  graphics::mtext(.ig_report_classifier_subtitle(v$model_id, basis_label),
    side = 3, line = .6, cex = .8)
  invisible(NULL)
}

# Keep the model/basis subtitle on its own line below the figure title.
.ig_report_classifier_subtitle <- function(model_id, basis_label) {
  lines <- .ig_report_wrap_pdf_lines(model_id, width_in = 4, size = 10, family = "sans")
  short <- paste0(lines[1L], if (length(lines) > 1L) " [...]" else "")
  paste(short, "|", basis_label)
}
