#' Check the offline PDF report renderer
#' @return A list with available and reason. No software is installed.
#' @export
report_pdf_capabilities <- function() {
  ok <- isTRUE(capabilities("cairo"))
  list(available = ok, renderer = "R Cairo PDF with grid and base graphics",
       reason = if (ok) "" else "This R installation needs Cairo support for UTF-8 PDF reports; JSON and TXT remain available.")
}

.ig_report_literal <- function(x) gsub("[[:cntrl:]]", " ", enc2utf8(as.character(x)))

.ig_report_lines <- function(x, prefix = "") {
  if (is.null(x)) return(paste0(prefix, ": null"))
  if (is.list(x)) {
    if (!is.null(x$metric_id)) {
      n <- length(x$metric_id)
      keys <- intersect(c("metric_id", "variant", "value", "status", "reason", "denominator",
        "denominator_value", "unit", "estimate_kind", "method", "level", "lower", "upper", "conditioning"), names(x))
      return(vapply(seq_len(n), function(i) {
        cells <- vapply(keys, function(key) {
          col <- x[[key]]
          val <- if (length(col) >= i) col[[i]] else NULL
          paste0(key, "=", if (is.null(val)) "null" else if (is.numeric(val)) paste(vapply(val, .ig_format_double, ""), collapse = ",") else paste(unlist(val), collapse = ","))
        }, "")
        paste0(prefix, ": ", paste(cells, collapse = " | "))
      }, ""))
    }
    if (!length(x)) return(paste0(prefix, ": []"))
    n <- names(x)
    if (is.null(n)) n <- as.character(seq_along(x))
    return(unlist(lapply(seq_along(x), function(i) .ig_report_lines(x[[i]],
      if (!nzchar(prefix)) n[i] else paste(prefix, n[i], sep = "."))), use.names = FALSE))
  }
  value <- if (is.numeric(x)) vapply(x, .ig_format_double, "") else .ig_report_literal(x)
  paste0(prefix, ": ", paste(value, collapse = "; "))
}

.ig_report_appendix_text <- function(x, pdf = FALSE) {
  labels <- if (x$locale == "de") c("Einstellungen", "Empfehlung", "Ergebnisse", "Unsicherheit",
    "Validierung", "Provenienz", "Reproduzierbarkeit", "Ressourcen", "Abbildungen", "Entwurfsvarianten") else
    c("Settings", "Recommendation", "Results", "Uncertainty", "Validation", "Provenance", "Reproducibility", "Resources", "Visualisations", "Design sweeps")
  sections <- list(x$specifications, x$recommendation, x$results, x$uncertainty,
    x$validation, x$provenance, x$reproducibility, x$resources, x$visualisations, x$sweeps)
  if (pdf) sections[[9]] <- lapply(x$visualisations, function(v) v[setdiff(names(v), "data")])
  if (!is.null(x$probability_summaries)) {
    sections <- c(sections,list(x$probability_summaries))
    labels <- c(labels,if(x$locale=="de") "Freigegebene Wahrscheinlichkeitsaggregate" else "Included probability aggregates")
  }
  lines <- c(.ig_report_literal(x$title), paste("Snapshot:", x$snapshot_id),
    paste("Analysis hash:", x$analysis_hash), paste("Created UTC:", x$created_at),
    paste("Display timezone:", x$display_timezone), paste("Language:", x$locale),
    paste("Scope:", x$report_kind, "|", x$completeness),
    "Numeric values use a decimal point and round-trip precision.",
    .ig_report_lines(x$privacy, "Privacy"))
  for (i in seq_along(sections)) lines <- c(lines, "", labels[i], .ig_report_lines(sections[[i]]))
  enc2utf8(lines)
}

.ig_report_text <- function(x, pdf = FALSE) {
  c(.ig_report_summary_lines(x), "", if(x$locale=="de") "TECHNISCHER ANHANG" else "TECHNICAL APPENDIX",
    .ig_report_appendix_text(x,pdf=pdf))
}

# Retain every wrapped line while avoiding a very short trailing text page.
.ig_report_page_chunks <- function(lines, capacity, minimum_tail = floor(capacity / 2), paragraphs = FALSE) {
  if (!length(lines)) return(list())
  if (isTRUE(paragraphs) && any(!nzchar(lines))) {
    ends <- unique(c(which(!nzchar(lines)), length(lines)))
    starts <- c(1L, utils::head(ends, -1L) + 1L)
    pages <- list(); current <- character()
    for (i in seq_along(ends)) {
      block <- lines[seq.int(starts[i], ends[i])]
      if (length(block) > capacity) {
        if (length(current)) {
          room <- capacity - length(current)
          if (room >= 3L) {
            current <- c(current, utils::head(block, room))
            block <- utils::tail(block, -room)
          }
          pages[[length(pages) + 1L]] <- current; current <- character()
        }
        parts <- .ig_report_page_chunks(block, capacity, minimum_tail, paragraphs = FALSE)
        if (length(parts) > 1L) pages <- c(pages, utils::head(parts, -1L))
        current <- parts[[length(parts)]]
      } else {
        if (length(current) + length(block) > capacity) {
          pages[[length(pages) + 1L]] <- current; current <- character()
        }
        current <- c(current, block)
      }
    }
    if (length(current)) pages[[length(pages) + 1L]] <- current
    if (length(pages) > 1L) for (i in seq.int(length(pages), 2L)) {
      while (length(pages[[i]]) < minimum_tail) {
        previous <- pages[[i - 1L]]
        boundaries <- which(!nzchar(previous))
        boundaries <- boundaries[boundaries < length(previous)]
        if (!length(boundaries)) break
        cut <- utils::tail(boundaries, 1L)
        moved <- previous[seq.int(cut + 1L, length(previous))]
        if (cut < minimum_tail || length(moved) + length(pages[[i]]) > capacity) break
        pages[[i - 1L]] <- utils::head(previous, cut)
        pages[[i]] <- c(moved, pages[[i]])
      }
    }
    return(pages)
  }
  ends <- seq.int(capacity, length(lines) + capacity - 1L, by = capacity)
  ends[length(ends)] <- length(lines)
  if (length(ends) > 1L && ends[length(ends)] - ends[length(ends)-1L] < minimum_tail)
    ends[length(ends)-1L] <- ends[length(ends)] - minimum_tail
  starts <- c(1L, utils::head(ends,-1L)+1L)
  Map(function(a,b) lines[seq.int(a,b)], starts, ends)
}

.ig_report_pdf <- function(x, file) {
  cap <- report_pdf_capabilities()
  if (!cap$available) .ig_abort("unsupported", cap$reason)
  grDevices::cairo_pdf(file, width = 8.27, height = 11.69, onefile = TRUE)
  device <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(device), add = TRUE)
  page <- 0L
  text_pages <- function(lines, appendix = FALSE) {
    lines <- .ig_report_wrap_pdf_lines(lines, width_in=7.1,
      size=if(appendix) 9 else 11, family=if(appendix) "mono" else "sans")
    count <- if(appendix) 54L else 50L
    pages <- .ig_report_page_chunks(lines, count, paragraphs = !appendix)
    for(content in pages) {
      if(page>0L) grid::grid.newpage()
      page <<- page + 1L
      heading <- if(appendix) {
        if(x$locale=="de") "Technischer Anhang" else "Technical appendix"
      } else if(x$locale=="de") "inferguidR Analysebericht" else "inferguidR analysis report"
      grid::grid.text(heading,x=.06,y=.96,just="left",gp=grid::gpar(fontsize=14,fontface="bold"))
      grid::grid.text(content,x=.06,y=seq(.91,by=if(appendix) -.015 else -.0168,length.out=length(content)),
        just=c("left","top"),gp=grid::gpar(fontsize=if(appendix) 9 else 11,fontfamily=if(appendix) "mono" else "sans"))
      grid::grid.text(paste(x$snapshot_id,"|",page),x=.06,y=.035,just="left",gp=grid::gpar(fontsize=8,col="#555555"))
    }
  }
  text_pages(.ig_report_summary_lines(x))
  for (i in seq_along(x$visualisations)) {
    v <- x$visualisations[[i]]
    if (!identical(v$status, "ok")) next
    graphics::par(mar = c(9, 5, 5, 2), oma = c(2, 0, 1, 0))
    if (v$plot_type == "distribution" && length(v$data)) {
      .ig_report_draw_distribution(v, x$results[[v$analysis]], x$locale)
    } else if (v$plot_type == "power_sweep") {
      .ig_report_draw_sweep(v,x$locale)
    } else if (v$plot_type == "power_lattice") {
      .ig_report_draw_lattice(v,x$locale)
    } else if (v$plot_type == "two_group_counts") {
      values <- rbind(unlist(v$data$event), unlist(v$data$non_event))
      labels <- vapply(unlist(v$data$group_labels),.ig_report_plot_label,"",width_in=2.4,size=9.6)
      centres <- graphics::barplot(values, names.arg = rep("", length(labels)), cex.names = .8,
        ylab = if (x$locale == "de") "Beobachtete Anzahl" else "Observed count",
        main = if (x$locale == "de") "Unabh\u00e4ngige Gruppen" else "Independent groups",
        col = c("#304f9c", "#b7c8d9"),
        legend.text = c(.ig_report_plot_label(v$data$event_label,width_in=2.3,size=9.6,max_lines=2L),
          if (x$locale == "de") "Kein Ereignis" else "Non-event"),
        args.legend = list(x = "topright", bty = "n", cex = .8), ylim = c(0, max(colSums(values)) * 1.25))
      for (j in seq_along(labels)) {
        label_lines <- strsplit(labels[j], "\n", fixed = TRUE)[[1L]]
        for (k in seq_along(label_lines))
          graphics::mtext(label_lines[k], side = 1, at = centres[j], line = 1 + (k - 1L) * .8, cex = .8)
      }
    } else if (v$plot_type == "confusion_counts") {
      values <- unlist(v$data, use.names = TRUE)
      graphics::barplot(values, ylab = if(x$locale=="de") "Anzahl" else "Count", main = .ig_report_plot_label(paste(if(x$locale=="de") "Gewaehlter Schwellenwert:" else "Selected cutoff:",v$model_id),width_in=6,size=14), col = "#304f9c")
    } else if (v$plot_type %in% c("classifier_roc_point", "classifier_pr_point")) {
      .ig_report_draw_classifier_operating(v, x$locale)
    } else if (v$plot_type == "classifier_comparison") {
      .ig_report_draw_classifier_comparison(v, x$locale)
    } else if (v$plot_type %in% c("classifier_mosaic", "classifier_prevalence")) {
      .ig_report_draw_classifier_aggregate(v, x$locale)
    } else if (v$plot_type %in% c("classifier_calibration", "classifier_net_benefit")) {
      .ig_report_draw_probability(v, x$locale)
    } else next
    caption_lines <- .ig_report_wrap_pdf_lines(paste(.ig_report_literal(.ig_report_caption(v,x$locale)),if(x$locale=="de") "Vollstaendige Bezeichnungen stehen im Anhang." else "Full labels are retained in the appendix."),width_in=if(v$plot_type=="classifier_comparison") 4.8 else 6.5,size=8.4,family="sans")
    if(length(caption_lines)>4L) caption_lines <- c(utils::head(caption_lines,3L),if(x$locale=="de") "Vollstaendige Bildunterschrift im technischen Anhang." else "Full caption in the technical appendix.")
    for (k in seq_along(caption_lines)) {
      graphics::mtext(caption_lines[k], side = 1, line = 5 + (k - 1L) * .8, cex = .7)
    }
    page <- page + 1L
    graphics::mtext(paste(x$snapshot_id, "|", page, if(x$locale=="de") "| Abbildung" else "| Figure", i), side = 1, outer = TRUE, line = 1, cex = .7)
  }
  text_pages(.ig_report_appendix_text(x,pdf=TRUE),appendix=TRUE)
  invisible(file)
}

#' Export a frozen report as PDF, JSON or UTF-8 text
#' @param snapshot A report_document from create_report_snapshot().
#' @param file Destination file. Existing files are preserved unless overwrite is TRUE.
#' @param format One of pdf, json or txt.
#' @param overwrite Allow replacement of an existing file.
#' @return Invisibly, the destination and its SHA-256 hash with the snapshot ID.
#' @details Rendering uses a private temporary file. JSON and TXT remain
#'   available if the offline Cairo PDF capability is unavailable. Text is
#'   rendered literally; report content is never evaluated as R or markup.
#'   PDF and TXT begin with a readable summary rounded to six significant
#'   digits; stored decisions and the technical appendix retain full precision.
#'   PDF figures precede the appendix and have adjacent captions. Very long
#'   plot labels may be abbreviated with `[...]` while full labels remain in
#'   the appendix. PDF text wrapping uses measured font widths.
#' @export
export_report <- function(snapshot, file, format = c("pdf", "json", "txt"), overwrite = FALSE) {
  format <- match.arg(format)
  if (!inherits(snapshot, "report_document")) .ig_abort("invalid_input", "snapshot must be a portable report_document")
  snapshot <- .ig_report_document_from_list(unclass(snapshot))
  if (!.ig_is_string(file) || !nzchar(file) || grepl("^[A-Za-z][A-Za-z0-9+.-]*://", file) ||
      !dir.exists(dirname(file)) || dir.exists(file)) .ig_abort("invalid_input", "file must be a local destination in an existing directory")
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) .ig_abort("invalid_input", "overwrite must be TRUE or FALSE")
  if (file.exists(file) && !overwrite) .ig_abort("invalid_input", "destination exists; use overwrite = TRUE to replace it")
  temporary <- tempfile(pattern = "inferguidR-report-", fileext = paste0(".", format))
  on.exit(unlink(temporary), add = TRUE)
  switch(format,
    json = writeLines(as.character(ig_to_json(snapshot, pretty = TRUE)), temporary, useBytes = TRUE),
    txt = writeLines(.ig_report_text(snapshot), temporary, useBytes = TRUE),
    pdf = .ig_report_pdf(snapshot, temporary))
  if (!file.exists(temporary) || file.size(temporary) <= 0) .ig_abort("numerical_failure", "report renderer produced no file")
  if (!file.copy(temporary, file, overwrite = overwrite)) .ig_abort("invalid_input", "could not write report destination")
  invisible(list(file = file, format = format, snapshot_id = snapshot$snapshot_id,
                 sha256 = digest::digest(file = file, algo = "sha256")))
}

# Localise known explanatory prose while retaining canonical snapshot metadata.
.ig_report_caption <- function(v, locale) {
  value <- v$caption
  if (locale != "de") return(value)
  source <- c(
    "Display uses central 99.998% quantiles of each law plus finite critical values. Curve sampling and display extent do not determine reported probabilities.",
    "Full discrete support; boundary counts are included in rejection.",
    "Displayed densities or probability masses describe the fixed design.",
    "Observed event and non-event counts in independent groups. These are not classifier confusion cells.",
    "Selected-cutoff cells. AUC summarises the full curve; prevalence-scenario values are expected, not observed measurements.")
  translated <- c(
    "Gezeigt werden die zentralen 99,998 % jeder Verteilung sowie endliche kritische Werte. Kurvenpunkte und Darstellungsbereich bestimmen nicht die berichteten Wahrscheinlichkeiten.",
    "Vollstaendiger diskreter Wertebereich; Randwerte sind im Ablehnungsbereich enthalten.",
    "Dichten und Wahrscheinlichkeitsmassen beschreiben das feste Design.",
    "Beobachtete Ereignisse und Nicht-Ereignisse in unabhaengigen Gruppen. Dies sind keine Klassifikationszaehlungen.",
    "Zaehlungen am ausgewaehlten Schwellenwert. Die AUC beschreibt die gesamte Kurve. Praevalenzszenarien zeigen erwartete statt beobachteter Werte.")
  for (i in seq_along(source)) value <- gsub(source[i], translated[i], value, fixed = TRUE)
  value
}
