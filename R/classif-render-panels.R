.ig_plot_classification_additional <- function(x, panel, model_id, basis, lang, mosaic_mode, bins) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    .ig_abort("unsupported", "plot_classification() requires the suggested package 'ggplot2'.")
  if (panel %in% c("score_distribution", "prevalence_curve", "model_comparison")) basis <- "observed"
  pd <- classification_plot_data(x, panel, model_id, basis, lang, bins = bins)
  de <- lang == "de"
  x0 <- x1 <- y0 <- y1 <- cell <- label <- score <- value <- class <- prevalence <- metric <-
    fpr <- tpr <- model <- marker <- segment <- NULL
  if (panel == "mosaic") {
    tab <- pd$data$rectangles[pd$data$rectangles$mode == mosaic_mode, , drop = FALSE]
    tab$label <- toupper(tab$cell)
    p <- ggplot2::ggplot(tab, ggplot2::aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = cell)) +
      ggplot2::geom_rect(colour = "white", linewidth = .8) +
      ggplot2::geom_text(ggplot2::aes(x = (x0 + x1) / 2, y = (y0 + y1) / 2, label = label),
                         data = tab[(tab$x1 - tab$x0) >= .12 & (tab$y1 - tab$y0) >= .08, , drop = FALSE],
                         inherit.aes = FALSE, size = 5) +
      ggplot2::scale_fill_manual(values = c(tp = "#72BCB3", fn = "#C8E8E4", fp = "#D5BA85", tn = "#ACBAE1"),
        breaks = c("tp", "fn", "fp", "tn"), labels = if (de)
          c("TP: richtig positiv", "FN: falsch negativ", "FP: falsch positiv", "TN: richtig negativ") else
          c("TP: true positive", "FN: false negative", "FP: false positive", "TN: true negative")) +
      ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::labs(x = if (de) "Ereignis links; kein Ereignis rechts" else "Event on left; non-event on right",
                    y = NULL, fill = NULL,
                    caption = if (mosaic_mode == "joint") {
                      if (de) "Gemeinsame Anteile: Flaechen entsprechen Anteilen an allen Einheiten." else
                        "Joint shares: areas represent proportions of all units."
                    } else if (de) "Innerhalb jeder Klasse: Hoehen sind bedingte Raten. Flaechen sind keine gemeinsamen Anteile." else
                      "Within each class: heights are conditional rates. Areas are not joint shares.")
  } else if (panel == "score_distribution") {
    if (!is.null(pd$data$densities)) {
      d <- pd$data$densities
      tab <- rbind(data.frame(score = d$z, value = d$density_events, class = "event"),
                   data.frame(score = d$z, value = d$density_non_events, class = "non_event"))
      p <- ggplot2::ggplot(tab, ggplot2::aes(score, value, colour = class, linetype = class)) +
        ggplot2::geom_line(linewidth = .8) +
        ggplot2::labs(x = if (de) "Standardisierter Score" else "Standardised score", y = if (de) "Dichte" else "Density")
    } else {
      d <- pd$data$display_bins
      tab <- rbind(data.frame(score = (d$bin_lower + d$bin_upper) / 2, value = d$share_of_events, class = "event"),
                   data.frame(score = (d$bin_lower + d$bin_upper) / 2, value = d$share_of_non_events, class = "non_event"))
      p <- ggplot2::ggplot(tab, ggplot2::aes(score, value, colour = class, linetype = class)) +
        ggplot2::geom_segment(ggplot2::aes(xend = score, yend = 0), linewidth = .7,
                              position = ggplot2::position_dodge(width = diff(d$bin_lower[1:2]) * .3)) +
        ggplot2::geom_point(size = 1.5, position = ggplot2::position_dodge(width = diff(d$bin_lower[1:2]) * .3)) +
        ggplot2::labs(x = "Score", y = if (de) "Anteil innerhalb der Klasse" else "Share within class")
    }
    cut <- pd$data$cutoff$cutoff
    if (length(cut) == 1 && is.finite(cut)) p <- p + ggplot2::geom_vline(xintercept = cut, linetype = "dotted", colour = "#956000")
    p <- p + ggplot2::scale_colour_manual(values = c(event = "#007C78", non_event = "#334E9A"),
      labels = if (de) c(event = "Ereignis", non_event = "Kein Ereignis") else c(event = "Event", non_event = "Non-event")) +
      ggplot2::scale_linetype_manual(values = c(event = "solid", non_event = "longdash"),
        labels = if (de) c(event = "Ereignis", non_event = "Kein Ereignis") else c(event = "Event", non_event = "Non-event")) +
      ggplot2::labs(colour = NULL, linetype = NULL, caption = if (de)
        "Darstellungsraster bestimmt keine Kennzahlen. Gepunktet: numerischer Schwellenwert, falls vorhanden." else
        "Display bins or grid do not determine metrics. Dotted: numeric cutoff, when available.")
  } else if (panel == "prevalence_curve") {
    d <- pd$data$curve
    tab <- rbind(data.frame(prevalence = d$prevalence, value = d$ppv, metric = "PPV"),
                 data.frame(prevalence = d$prevalence, value = d$npv, metric = "NPV"))
    # Undefined regions split the path before rows are filtered.
    tab$segment <- stats::ave(as.integer(!is.finite(tab$value)), tab$metric, FUN = cumsum)
    tab <- tab[is.finite(tab$value), , drop = FALSE]
    p <- ggplot2::ggplot(tab, ggplot2::aes(prevalence, value, colour = metric, linetype = metric,
                                         group = interaction(metric, segment))) +
      ggplot2::geom_line(linewidth = .8) +
      ggplot2::scale_colour_manual(values = c(PPV = "#007C78", NPV = "#334E9A")) +
      ggplot2::scale_linetype_manual(values = c(PPV = "solid", NPV = "longdash")) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::labs(x = if (de) "Angenommene Praevalenz" else "Assumed prevalence", y = if (de) "Vorhersagewert" else "Predictive value",
        colour = NULL, linetype = NULL, caption = if (de)
          "Szenario mit festen TPR/FPR; keine Messung, keine Unsicherheitsintervalle. Raster: 0,001 bis 0,999." else
          "Scenario with fixed TPR/FPR; not a measurement; no uncertainty intervals. Grid: 0.001 to 0.999.")
  } else {
    d <- pd$data$curves
    if (is.null(d) || !nrow(d)) .ig_abort("unsupported", "No model has an available ROC curve.")
    d$model <- d$model_id
    pts <- pd$data$operating_points
    ids <- unique(d$model)
    p <- ggplot2::ggplot(d, ggplot2::aes(fpr, tpr, colour = model, linetype = model)) +
      ggplot2::scale_colour_manual(values = stats::setNames(grDevices::hcl.colors(length(ids), "Dark 3"), ids)) +
      ggplot2::scale_linetype_manual(values = stats::setNames(rep(c("solid", "longdash", "dotted", "dotdash"), length.out = length(ids)), ids)) +
      ggplot2::geom_path(linewidth = .8) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::labs(x = if (de) "Falsch-Positiv-Rate" else "False positive rate", y = if (de) "Richtig-Positiv-Rate" else "True positive rate",
        colour = if (de) "Modell" else "Model", linetype = if (de) "Modell" else "Model",
        caption = if (de) "Jedes Modell behaelt seine eigene Regel; gleiche Zahlenwerte sind keine gleichwertigen Schwellenwerte." else
          "Each model keeps its own rule; equal numeric cutoffs are not equivalent operating points.")
    if (!is.null(pts) && nrow(pts)) {
      pts$model <- pts$model_id
      p <- p + ggplot2::geom_point(data = pts, ggplot2::aes(shape = marker), size = 3, stroke = 1) +
        ggplot2::scale_shape_manual(values = c(current_operating_point = 21, comparison_operating_point = 2),
          labels = if (de) c(current_operating_point = "Aktuelles Modell", comparison_operating_point = "Vergleichsmodell") else
            c(current_operating_point = "Current model", comparison_operating_point = "Comparison model")) + ggplot2::labs(shape = NULL)
    }
  }
  titles <- if (de) c(mosaic = "Entscheidungsmosaik", score_distribution = "Score-Verteilungen",
    prevalence_curve = "Vorhersagewerte und Praevalenz", model_comparison = "Modellvergleich") else
    c(mosaic = "Decision mosaic", score_distribution = "Score distributions",
      prevalence_curve = "Predictive values and prevalence", model_comparison = "Model comparison")
  p <- p + ggplot2::labs(title = paste(titles[[panel]], pd$model_id, sep = ": "),
                          subtitle = paste(if (de) "Bezugsbasis:" else "Basis:", pd$basis)) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(legend.position = "bottom", plot.caption = ggplot2::element_text(hjust = 0, size = 10),
                    panel.grid.minor = ggplot2::element_blank())
  p$labels$caption <- paste(strwrap(p$labels$caption, width = 80), collapse = "\n")
  p <- p + ggplot2::theme(plot.caption.position = "plot", plot.title.position = "plot")
  attr(p, "inferguidR_plot_data") <- pd
  p
}
