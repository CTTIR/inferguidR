.ig_plot_classification_network <- function(x, model_id, basis, lang, metric_focus) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    .ig_abort("unsupported", "plot_classification() requires the suggested package 'ggplot2'.")
  pd <- classification_plot_data(x, "metric_network", model_id, basis, lang)
  if (!.ig_is_string(metric_focus) || !metric_focus %in% pd$data$nodes$id)
    .ig_abort("invalid_input", "Choose an operating-point metric for the dependency view.")
  de <- lang == "de"
  edges <- pd$data$edges[pd$data$edges$to == metric_focus, , drop = FALSE]
  parents <- unique(edges$from)
  ids <- c(parents, metric_focus)
  tab <- pd$data$nodes[match(ids, pd$data$nodes$id), , drop = FALSE]
  tab$id <- ids
  tab$short <- .ig_cls_short_label(ids, lang)
  tab$label <- paste(tab$short, ifelse(is.finite(tab$value),
    formatC(tab$value, digits = 5, format = "g", decimal.mark = if (de) "," else "."), tab$status), sep = "\n")
  tab$px <- c(rep(.1, length(parents)), .9)
  tab$py <- c(if (length(parents)) seq(.1, .9, length.out = length(parents)), .5)
  edges$px <- rep(.26, nrow(edges))
  edges$py <- tab$py[match(edges$from, tab$id)]
  edges$ex <- rep(.74, nrow(edges))
  edges$ey <- rep(.5, nrow(edges))
  px <- py <- ex <- ey <- label <- kind <- NULL
  p <- ggplot2::ggplot(tab, ggplot2::aes(px, py))
  if (nrow(edges)) p <- p + ggplot2::geom_segment(data = edges,
    ggplot2::aes(x = px, y = py, xend = ex, yend = ey, linetype = kind),
    inherit.aes = FALSE, colour = "#697984", linewidth = .7)
  p <- p + ggplot2::geom_label(ggplot2::aes(label = label), fill = "#EDF5F4", colour = "#152330", size = 4) +
    ggplot2::coord_cartesian(xlim = c(-.1, 1.1), ylim = c(0, 1), clip = "off")
  if (nrow(edges)) p <- p +
    ggplot2::scale_linetype_manual(values = c(input = "solid", bayes = "longdash"),
      labels = if (de) c(input = "Direkte Eingaben", bayes = "Bayes-Beziehung") else
        c(input = "Direct inputs", bayes = "Bayes relationship"))
  dict <- metric_dictionary(lang)
  formula <- dict$formula[match(metric_focus, dict$id)]
  caption <- paste(if (de) "Links: Eingaben; rechts: ausgewaehlte Kennzahl. Gezeigt werden direkte Abhaengigkeiten, keine kausalen Beziehungen." else
    "Left: inputs; right: selected metric. Direct calculation dependencies, not causal relationships.", formula)
  p <- p + ggplot2::labs(title = paste(if (de) "Abhaengigkeiten:" else "Dependencies:", .ig_cls_short_label(metric_focus, lang)),
    subtitle = paste(if (de) "Bezugsbasis:" else "Basis:", pd$basis),
    caption = paste(strwrap(caption, width = 80), collapse = "\n"), linetype = NULL) +
    ggplot2::theme_void(base_size = 14) +
    ggplot2::theme(legend.position = "bottom", plot.caption = ggplot2::element_text(hjust = 0, size = 10),
      plot.caption.position = "plot", plot.margin = ggplot2::margin(10, 15, 10, 15))
  attr(p, "inferguidR_plot_data") <- pd
  p
}
