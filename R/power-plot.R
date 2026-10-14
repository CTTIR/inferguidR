#' Draw a power analysis panel
#'
#' Uses [power_plot_data()] for both the plotted values and availability
#' checks. The plot is a static counterpart for application and report use.
#' Distribution curves are display samples; decision probabilities come from
#' the engine. The returned plot retains its source data in the
#' `inferguidR_plot_data` attribute.
#'
#' @param result An unmodified `power_result`.
#' @param panel `"distribution"`, `"decision_table"`, `"mosaic"`, or `"lattice"`.
#' @param lang Label language, `"en"` or `"de"`.
#' @return A `ggplot` object. Requires the suggested package 'ggplot2'.
#'   An unavailable panel raises a typed unsupported condition.
#' @export
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   s <- analysis_spec("t_means_one_sample", "post_hoc",
#'     inputs = list(d = .5, alpha = .05, n_total = 30),
#'     options = list(tails = "two"))
#'   plot_power(calculate_power(s))
#' }
plot_power <- function(result, panel = c("distribution", "decision_table", "mosaic", "lattice"),
                       lang = c("en", "de")) {
  panel <- match.arg(panel)
  lang <- match.arg(lang)
  pd <- power_plot_data(result, panel)
  if (!identical(pd$status, "ok")) .ig_abort("unsupported", pd$reason)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    .ig_abort("unsupported", "plot_power() requires the suggested package 'ggplot2'.")
  if (panel == "lattice") return(.ig_draw_power_lattice(pd, lang))
  # Local bindings keep non-standard aesthetic evaluation explicit to checks.
  x <- y <- hypothesis <- decision <- probability <- xmin <- xmax <- ymin <- ymax <- label <- error <- group <- NULL
  de <- lang == "de"
  hypothesis_labels <- if (de) c("H0: Nullhypothese", "H1: Alternative") else
    c("H0: null hypothesis", "H1: alternative")
  caption <- if (panel == "distribution") {
    if (de) paste("Fehlerwahrscheinlichkeiten stehen in der Entscheidungstabelle.",
      if (pd$ordinate == "density") "Dargestellt: zentrale 99,998 % jeder Verteilung und kritische Werte." else
        "Dargestellt: vollst\u00e4ndiger diskreter Tr\u00e4ger.", sep = "\n") else
      paste("Engine error probabilities are in the decision table.",
        if (pd$ordinate == "density") "Shown: central 99.998% of each distribution and critical values." else
          "Shown: full discrete support.", sep = "\n")
  } else {
    if (de) "Jede Spalte bedingt auf eine Hypothese; keine gemeinsamen oder pr\u00e4diktiven Wahrscheinlichkeiten." else
      "Each column conditions on one hypothesis; these are not joint or predictive probabilities."
  }
  tab <- pd$data
  if (panel == "distribution") {
    bands <- .ig_power_error_bands(tab, pd$critical, pd$ordinate)
    error_labels <- .ig_power_error_labels(result$achieved$alpha_attained, result$achieved$beta, lang)
    p <- ggplot2::ggplot(tab, ggplot2::aes(x = x, y = y, colour = hypothesis, linetype = hypothesis))
    if (identical(pd$ordinate, "density") && nrow(bands)) p <- p +
      ggplot2::geom_ribbon(data = bands, ggplot2::aes(x=x,ymin=0,ymax=y,fill=error,group=group),
        inherit.aes=FALSE,alpha=.4)
    p <- p + ggplot2::scale_fill_manual(values=c(alpha="#D98920",beta="#A9C9C3"),
      breaks=c("alpha","beta"),labels=error_labels,drop=FALSE,name=NULL)
    if (identical(pd$ordinate, "probability_mass")) {
      p <- p + ggplot2::geom_segment(ggplot2::aes(xend = x, yend = 0), linewidth = .55) +
        ggplot2::geom_point(size = 1.5)
    } else p <- p + ggplot2::geom_line(linewidth = .85)
    if (identical(pd$ordinate, "probability_mass") && nrow(bands)) p <- p +
      ggplot2::geom_point(data=bands,ggplot2::aes(x=x,y=y,fill=error),inherit.aes=FALSE,
        shape=22,size=2.5,colour="#444444")
    caption <- paste(caption, if(de) "Markierte Bereiche: alpha unter H0 und beta unter H1; Werte aus dem Rechenkern." else
      "Marked regions: alpha under H0 and beta under H1; values come from the engine.",sep="\n")
    finite_critical <- pd$critical[is.finite(pd$critical)]
    if (length(finite_critical)) p <- p + ggplot2::geom_vline(
      xintercept = finite_critical, colour = "#956000", linetype = "dotted", linewidth = .6)
    p <- p + ggplot2::scale_colour_manual(values = c(H0 = "#334E9A", H1 = "#007C78"),
                                         breaks = c("H0", "H1"), labels = hypothesis_labels) +
      ggplot2::scale_linetype_manual(values = c(H0 = "solid", H1 = "longdash"),
                                    breaks = c("H0", "H1"), labels = hypothesis_labels) +
      ggplot2::labs(subtitle = paste(if (de) "Gepunktete Linien: kritische Werte" else "Dotted lines: critical values",
                                      paste(formatC(finite_critical, digits = 5, format = "g",
                                                    decimal.mark = if (de) "," else "."), collapse = "; ")),
                    x = if (de) "Teststatistik" else "Test statistic",
                    y = if (pd$ordinate == "density") {
                      if (de) "Dichte" else "Density"
                    } else if (de) "Wahrscheinlichkeitsmasse" else "Probability mass",
                    colour = NULL, linetype = NULL)
  } else {
    tab$label <- formatC(tab$probability, digits = 5L, format = "g", decimal.mark = if (de) "," else ".")
    tab$hypothesis <- factor(tab$hypothesis, levels = c("H0", "H1"), labels = hypothesis_labels)
    tab$decision <- factor(tab$decision, levels = c("reject", "retain"),
                            labels = if (de) c("H0 verwerfen", "H0 nicht verwerfen") else c("Reject H0", "Do not reject H0"))
    if (panel == "decision_table") {
      p <- ggplot2::ggplot(tab, ggplot2::aes(x = hypothesis, y = decision)) +
        ggplot2::geom_tile(fill = "#F2F5F8", colour = "white", linewidth = 2) +
        ggplot2::geom_text(ggplot2::aes(label = label), size = 5, colour = "#152330") +
        ggplot2::labs(x = NULL, y = NULL)
    } else {
      tab$x <- (tab$xmin + tab$xmax) / 2
      tab$y <- (tab$ymin + tab$ymax) / 2
      p <- ggplot2::ggplot(tab, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = decision)) +
        ggplot2::geom_rect(colour = "white", linewidth = .8) +
        ggplot2::scale_fill_manual(values = c("#334E9A", "#DCE4EA")) +
        ggplot2::scale_x_continuous(breaks = c(.5, 1.5), labels = hypothesis_labels) +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::labs(x = NULL, y = if (de) "Bedingte Wahrscheinlichkeit" else "Conditional probability", fill = NULL)
    }
  }
  p <- p + ggplot2::labs(caption = caption) +
    ggplot2::theme_minimal(base_size = 14, base_family = "sans") +
    ggplot2::theme(text = ggplot2::element_text(colour = "#152330"),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "#DCE4EA"),
                   legend.position = "bottom", plot.caption = ggplot2::element_text(hjust = 0, size = 10))
  attr(p, "inferguidR_plot_data") <- pd
  p
}


.ig_draw_power_lattice <- function(pd, lang) {
  de <- identical(lang, "de")
  x1 <- x2 <- probability <- hypothesis <- NULL
  tab <- pd$data
  tab$hypothesis <- factor(tab$hypothesis, levels = c("H0", "H1"))
  pi0 <- formatC(pd$fixed_inputs$null_probability_at, digits = 5, format = "g",
    decimal.mark = if (de) "," else ".")
  caption <- if (de) paste0("Kreuze: H0 verwerfen; keine beobachteten Daten.
",
    "H0 am Maximum der Fehlerwahrscheinlichkeit (gemeinsames p = ", pi0, ").
",
    "Getrennte bedingte Verteilungen; keine pr\u00e4diktiven Wahrscheinlichkeiten.") else
    paste0("Crosses: reject H0; these are not observed data.
",
    "H0 at the maximum rejection probability (common p = ", pi0, ").
",
    "Separate conditional distributions; not predictive probabilities.")
  p <- ggplot2::ggplot(tab, ggplot2::aes(x = x1, y = x2, fill = probability)) +
    ggplot2::geom_tile() +
    ggplot2::geom_point(data = tab[tab$reject, , drop = FALSE], shape = 4, size = 2, colour = "#111111") +
    ggplot2::facet_wrap(~hypothesis) + ggplot2::coord_equal() +
    ggplot2::scale_fill_gradient(low = "#F0F6FF", high = "#167C80",
      labels = function(z) format(z, trim = TRUE, decimal.mark = if (de) "," else "."),
      guide = ggplot2::guide_colourbar(barwidth = grid::unit(5, "cm"))) +
    ggplot2::scale_x_continuous(breaks = if (pd$fixed_inputs$n1 <= 12) 0:pd$fixed_inputs$n1 else
      unique(round(pretty(c(0, pd$fixed_inputs$n1))))) +
    ggplot2::scale_y_continuous(breaks = if (pd$fixed_inputs$n2 <= 12) 0:pd$fixed_inputs$n2 else
      unique(round(pretty(c(0, pd$fixed_inputs$n2))))) +
    ggplot2::labs(x = if (de) "Ereignisse in Gruppe 1" else "Events in group 1",
      y = if (de) "Ereignisse in Gruppe 2" else "Events in group 2",
      fill = if (de) "Wahrscheinlichkeit" else "Probability", caption = caption) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", plot.caption = ggplot2::element_text(hjust = 0))
  attr(p, "inferguidR_plot_data") <- pd
  p
}
