# Split display segments before filtering invalid points, so gaps cannot close.
.ig_xy_render_segments <- function(x) {
  if (!inherits(x, "ig_xy_plot_data"))
    .ig_abort("invalid_input", "a computed ig_xy_plot_data object is required")
  tab <- x$points
  for (i in which(tab$status == "ok")) {
    result <- tab$result[[i]]
    envelope <- analysis_envelope(result$spec, result)
    values <- .ig_xy_base_values(result, .ig_get_method(result$method_id))
    if (!identical(tab$spec_hash[i], envelope$spec_hash) ||
        !identical(tab$result_hash[i], envelope$result_hash) ||
        !identical(as.double(tab$x[i]), as.double(values[[x$axes$x$quantity]])) ||
        !identical(as.double(tab$y[i]), as.double(values[[x$axes$y$quantity]])))
      .ig_abort("invalid_input", "sweep coordinates or hashes differ from the stored point result")
  }
  tab$point_index <- seq_len(nrow(tab))
  tab$segment <- NA_character_
  for (curve in unique(tab$curve)) {
    rows <- which(tab$curve == curve)
    ok <- tab$status[rows] == "ok" & is.finite(tab$x[rows]) & is.finite(tab$y[rows])
    groups <- cumsum(!ok)
    tab$segment[rows[ok]] <- paste(curve, groups[ok], sep = ":")
  }
  tab
}

#' Draw an evaluated design sweep without crossing missing results
#'
#' @param x Data returned by [xy_plot_data()].
#' @param lang Label language, `"en"` or `"de"`.
#' @return A static ggplot object retaining all point data and result identities.
#' @export
plot_design_sweep <- function(x, lang = c("en", "de")) {
  lang <- match.arg(lang)
  tab <- .ig_xy_render_segments(x)
  def <- .ig_get_method(x$method_id)
  axis_label <- function(role) .ig_xy_local_label(def, role, lang)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    .ig_abort("unsupported", "plot_design_sweep() requires the suggested package 'ggplot2'.")
  good <- tab[!is.na(tab$segment), , drop = FALSE]
  if (!nrow(good)) .ig_abort("unsupported", "the sweep contains no finite successful points")
  good$curve_label <- if (is.null(x$curves)) x$method_id else
    paste(axis_label(x$curves$quantity), format(good$curve_value, digits = 6, decimal.mark = if (lang == "de") "," else "."))
  good$curve_label <- factor(good$curve_label, levels = unique(good$curve_label))
  sizes <- table(good$segment)
  connected <- good[good$segment %in% names(sizes)[sizes >= 2], , drop = FALSE]
  x_value <- y_value <- curve_label <- segment <- NULL
  names(good)[names(good) == "x"] <- "x_value"
  names(good)[names(good) == "y"] <- "y_value"
  names(connected)[names(connected) == "x"] <- "x_value"
  names(connected)[names(connected) == "y"] <- "y_value"
  p <- ggplot2::ggplot(good, ggplot2::aes(x = x_value, y = y_value, colour = curve_label,
                                        shape = curve_label, group = segment))
  if (nrow(connected)) {
    p <- p + if (x$render == "step") ggplot2::geom_step(data = connected, direction = "hv") else
      ggplot2::geom_line(data = connected)
  }
  de <- lang == "de"
  missing <- nrow(tab) - nrow(good)
  caption <- if (de) sprintf("%d Punkte ohne Ergebnis: keine Verbindung ueber Luecken. Jeder Punkt speichert seinen vollstaendigen Entwurf.", missing) else
    sprintf("%d points without results: no connections across gaps. Each point retains its full design.", missing)
  if (x$render == "step") caption <- paste(caption, if (de)
    "Stufen verbinden berechnete Entwuerfe; Uebergaenge zwischen Stuetzstellen sind nicht bestimmt." else
    "Steps connect evaluated designs; transition locations between samples are unresolved.", sep = "\n")
  palette <- rep(c("#334E9A", "#007C78", "#956000", "#783D8C", "#4F626E"), length.out = nlevels(good$curve_label))
  p <- p + ggplot2::scale_colour_manual(values = palette) +
    ggplot2::scale_shape_manual(values = seq_len(nlevels(good$curve_label)) - 1L) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(x = axis_label(x$axes$x$quantity), y = axis_label(x$axes$y$quantity),
                  title = paste(strwrap(def$label[[lang]], 70), collapse = "\n"),
                  subtitle = .ig_xy_local_fixed(x, def, lang),
                  colour = NULL, shape = NULL, caption = paste(unlist(lapply(strsplit(caption, "\n", fixed = TRUE)[[1]], strwrap, width = 100)), collapse = "\n")) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(legend.position = "bottom", plot.caption = ggplot2::element_text(hjust = 0, size = 10))
  if (x$axes$y$quantity %in% c("power", "alpha"))
    p <- p + ggplot2::scale_y_continuous(limits = c(0, 1))
  if (is.null(x$curves)) p <- p + ggplot2::guides(colour = "none", shape = "none")
  attr(p, "inferguidR_plot_data") <- x
  p
}

.ig_xy_local_label <- function(def, role, lang) {
  if (lang == "en") return(.ig_xy_label(def, role))
  switch(role, power = "Testst\u00e4rke (1-beta)", sample_size = "Gesamtstichprobengr\u00f6\u00dfe",
    alpha = "Fehlerwahrscheinlichkeit alpha", effect_size = {
      p <- def$params[[.ig_xy_effect_param(def)]]
      p$label$de %||% p$label$en
    })
}

.ig_xy_local_fixed <- function(x, def, lang) {
  labels <- names(x$fixed)
  values <- x$fixed
  if (lang == "de") {
    for (role in .ig_xy_roles) labels[labels == .ig_xy_label(def, role)] <- .ig_xy_local_label(def, role, lang)
    labels[labels == "Allocation ratio N2/N1"] <- "Zuteilungsverh\u00e4ltnis N2/N1"
    meta <- .ig_get_metadata(x$method_id)
    for (id in names(def$options)) {
      old <- meta$options[[id]]$ui_label %||% id
      labels[labels == old] <- def$options[[id]]$label$de %||% def$options[[id]]$label$en %||% old
    }
    values <- lapply(values, function(v) {
      if (identical(v, "One")) "Einseitig" else if (identical(v, "Two")) "Zweiseitig" else v
    })
  }
  show <- function(v) if (is.numeric(v)) format(v, digits = 6, decimal.mark = if (lang == "de") "," else ".") else as.character(v)
  paste(strwrap(paste(paste(labels, vapply(values, show, ""), sep = " = "), collapse = "; "), 85), collapse = "\n")
}
