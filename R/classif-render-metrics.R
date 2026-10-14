.ig_plot_classification_metrics <- function(x, panel, model_id, basis, lang,
                                            metric_family, metric_x, metric_y, log_x, log_y, realised_rules = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    .ig_abort("unsupported", "plot_classification() requires the suggested package 'ggplot2'.")
  for (flag in list(log_x, log_y)) if (!is.logical(flag) || length(flag) != 1L || is.na(flag))
    .ig_abort("invalid_input", "Log-axis settings must be TRUE or FALSE.")
  pd <- classification_plot_data(x, panel, model_id, basis, lang, metric_x = metric_x, metric_y = metric_y, realised_rules = realised_rules)
  de <- lang == "de"
  cutoff <- value <- short <- px <- py <- NULL
  if (panel == "metric_vs_cutoff") {
    if (!.ig_is_string(metric_family) || !metric_family %in% names(.ig_cls_families))
      .ig_abort("invalid_input", "Choose a supported metric family.")
    tab <- pd$data$values[pd$data$values$family == metric_family, , drop = FALSE]
    log_scale <- metric_family == "likelihood_ratios"
    drawable <- tab$on_cutoff_axis & is.finite(tab$cutoff) & tab$status == "ok" & is.finite(tab$value)
    if (log_scale) drawable <- drawable & tab$value > 0
    omitted <- sum(!drawable)
    data <- tab[drawable, , drop = FALSE]
    if (!nrow(data)) .ig_abort("unsupported", "No values are drawable on the selected axes; inspect the data table or change the axes.")
    p <- ggplot2::ggplot(data, ggplot2::aes(cutoff, value, colour = short, shape = short)) +
      ggplot2::geom_point(size = 1.8) +
      ggplot2::scale_colour_manual(values = c("#007C78", "#334E9A", "#956000")) +
      ggplot2::scale_shape_manual(values = c(16, 17, 15)) +
      ggplot2::labs(x = if (de) "Schwellenwert" else "Cutoff", y = if (de) "Kennzahl" else "Metric",
                    colour = NULL, shape = NULL)
    domain <- pd$data$domains[pd$data$domains$family == metric_family, , drop = FALSE]
    if (log_scale) p <- p + ggplot2::scale_y_log10() else
      if (is.finite(domain$lower) && is.finite(domain$upper))
        p <- p + ggplot2::coord_cartesian(ylim = c(domain$lower, domain$upper))
    ids <- vapply(x$models, function(m) m$model_id, "")
    op <- x$models[[match(pd$model_id, ids)]]$operating_point
    if (!is.null(op) && length(op$cutoff) == 1 && is.finite(op$cutoff))
      p <- p + ggplot2::geom_vline(xintercept = op$cutoff, linetype = "dotted", colour = "#697984")
    caption <- if (de) sprintf("Punkte sind erreichte Entscheidungsregeln, ohne Interpolation. Gepunktet: aktueller Schwellenwert. %d Werte ohne darstellbare Koordinaten stehen nur in der Datentabelle.", omitted) else
      sprintf("Points are attained decision rules, without interpolation. Dotted: current cutoff. %d values without drawable coordinates remain in the data table.", omitted)
    title <- if (de) "Kennzahlen nach Schwellenwert" else "Metrics by cutoff"
  } else {
    tab <- pd$data$points
    drawable <- tab$x_status == "ok" & tab$y_status == "ok" & is.finite(tab$x) & is.finite(tab$y)
    if (log_x) drawable <- drawable & tab$x > 0
    if (log_y) drawable <- drawable & tab$y > 0
    data <- tab[drawable, , drop = FALSE]
    if (!nrow(data)) .ig_abort("unsupported", "No values are drawable on the selected axes; inspect the data table or change the axes.")
    data$px <- data$x; data$py <- data$y
    p <- ggplot2::ggplot(data, ggplot2::aes(px, py)) + ggplot2::geom_point(colour = "#007C78", size = 1.8)
    current <- pd$data$current
    current_visible <- !is.null(current) && nrow(current) == 1 && current$x_status == "ok" && current$y_status == "ok" &&
      is.finite(current$x) && is.finite(current$y) && (!log_x || current$x > 0) && (!log_y || current$y > 0)
    if (current_visible) {
      current$px <- current$x; current$py <- current$y
      p <- p + ggplot2::geom_point(data = current, shape = 21, fill = NA, size = 4, stroke = 1.1)
    }
    if (log_x) p <- p + ggplot2::scale_x_log10()
    if (log_y) p <- p + ggplot2::scale_y_log10()
    p <- p + ggplot2::labs(x = pd$data$axes$short[1], y = pd$data$axes$short[2])
    caption <- if (de) sprintf("%d Entscheidungsregeln sind auf diesen Achsen nicht darstellbar und bleiben in der Tabelle; keine Begrenzung auf Achsenraender. Aktueller Punkt: %s.", sum(!drawable), if (current_visible) "Ring" else "nicht darstellbar") else
      sprintf("%d decision rules cannot be drawn on these axes and remain in the table; no clamping to axis bounds. Current point: %s.", sum(!drawable), if (current_visible) "ring" else "not drawable")
    rr <- pd$data$realised_rules
    if(!is.null(rr) && nrow(rr)) {
      visible <- rr$x_status=="ok" & rr$y_status=="ok" & is.finite(rr$x) & is.finite(rr$y)
      if(log_x) visible <- visible & rr$x>0
      if(log_y) visible <- visible & rr$y>0
      rr$px <- rr$x;rr$py <- rr$y
      if(any(visible)) p <- p + ggplot2::geom_point(data=rr[visible,,drop=FALSE],shape=17,colour="#8b3a62",size=3)
      caption <- paste(caption,if(de) "Dreiecke: separate realisierte Regeln." else "Triangles: separate realised rules.",
        if(de) paste(sum(!visible),"realisierte Regeln nicht darstellbar.") else paste(sum(!visible),"realised rules not drawable."))
    } else if(!is.null(realised_rules)) caption <- paste(caption,if(de)
      "Beobachtete realisierte Regeln werden auf praevalenzabhaengigen Szenarioachsen nicht gezeigt." else
      "Observed realised rules are hidden on prevalence-dependent scenario axes.")
    title <- if (de) "Kennzahlenpaar" else "Metric pair"
  }
  p <- p + ggplot2::labs(title = paste(title, pd$model_id, sep = ": "),
    subtitle = paste(if (de) "Bezugsbasis:" else "Basis:", pd$basis),
    caption = paste(strwrap(caption, width = 80), collapse = "\n")) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(legend.position = "bottom", plot.caption = ggplot2::element_text(hjust = 0, size = 10),
      plot.caption.position = "plot", panel.grid.minor = ggplot2::element_blank())
  attr(p, "inferguidR_plot_data") <- pd
  p
}
