# Display formatting of classification metric values ------------------------------------
#
# Formatting is presentation only: selection, constraints and every export use
# the unrounded values. The glyph conventions of the source explorer are kept
# (an approximately-equal sign for expected counts, an infinity sign for +Inf)
# with one correction: an undefined value is never a bare dash, it always comes
# with its reason (classification-source-audit ENG-12, VIS-15).

.ig_cls_expected_kinds <- c("expected", "expected_theoretical", "expected_scenario", "theoretical")

# Number formatting with locale marks; `kind` selects the rule.
.ig_cls_format_number <- function(v, kind, lang, digits) {
  dec <- if (lang == "de") "," else "."
  big <- if (lang == "de") "." else ","
  fixed <- function(x, d) formatC(x, format = "f", digits = d, decimal.mark = dec, big.mark = big)
  sig <- function(x, d) trimws(formatC(signif(x, d), format = "g", digits = d, decimal.mark = dec))
  switch(kind,
    percent = {
      p <- 100 * v
      d <- if (!is.null(digits)) digits else if (abs(p) >= 10) 1L else if (abs(p) >= 1) 2L else 3L
      txt <- if (p != 0 && abs(p) < 10^(-d)) sig(p, 2) else fixed(p, d)
      paste0(txt, if (lang == "de") "\u00a0%" else "%")
    },
    count = fixed(v, 0L),
    expected_count = paste0("\u2248", if (v != 0 && abs(v) < 0.1) sig(v, 2) else fixed(v, 1L)),
    ratio = if (abs(v) >= 100) fixed(v, 0L) else sig(v, if (is.null(digits)) 3L else digits),
    per_1000 = fixed(v, if (is.null(digits)) 1L else digits),
    fixed(v, if (is.null(digits)) 3L else digits))
}

#' Format classification metric values for display
#'
#' Adds display strings to a table of classification metrics without changing
#' its values. Finite values are formatted by unit: proportions as
#' percentages, counts as whole numbers (expected counts with an
#' approximately-equal sign and one decimal, because they are not observed
#' units), ratios with three significant digits (whole numbers from 100 on),
#' rates per 1,000 with one decimal and all other quantities with three
#' decimals. Positive infinity is shown as an infinity sign. An undefined, not
#' applicable or not computed value is shown as a dash, and its reason is
#' always given in `display_reason`: a dash is never left unexplained. German
#' output uses a decimal comma.
#'
#' Formatting never feeds a decision: cutoff selection, constraints and every
#' export use the unrounded values.
#'
#' @param x A data frame with `metric_id`, `value` and `status`, and
#'   optionally `reason_code`, `reason`, `unit` and `estimate_kind`, as
#'   returned by [classification_metrics()], [classification_table()] (tables
#'   `"metrics"`, `"scenario"` and `"summaries"`) or [prevalence_scenario()].
#' @param lang `"en"` or `"de"`.
#' @param percent Show proportions as percentages (default `TRUE`); otherwise
#'   with three decimals.
#' @param digits Optional fixed number of decimals (significant digits for
#'   ratios) instead of the magnitude-dependent default.
#' @return `x` with the added columns `display` and `display_reason`.
#' @examples
#' m <- format_metric_values(classification_metrics(c(tp = 0, fn = 5, fp = 0, tn = 95)))
#' m[m$metric_id %in% c("tn", "tpr", "ppv", "lr_neg"), c("metric_id", "display", "display_reason")]
#' sc <- prevalence_scenario(c(tp = 8, fn = 2, fp = 9, tn = 81), prevalence = 0.02)
#' format_metric_values(sc, lang = "de")[1:4, c("metric_id", "display")]
#' @concept classification
#' @seealso
#' * [classification_spec()]
#' * [evaluate_cutoff()]
#' * [plot_classification()]
#' @export
format_metric_values <- function(x, lang = c("en", "de"), percent = TRUE, digits = NULL) {
  lang <- match.arg(lang)
  if (!is.data.frame(x) || !all(c("metric_id", "value", "status") %in% names(x))) {
    .ig_abort("invalid_input", "`x` must be a data frame with metric_id, value and status columns")
  }
  if (!.ig_is_flag(percent)) .ig_abort("invalid_input", "`percent` must be TRUE or FALSE")
  if (!is.null(digits) && !(is.numeric(digits) && length(digits) == 1L && is.finite(digits) &&
                            digits >= 0 && digits <= 10 && digits == round(digits))) {
    .ig_abort("invalid_input", "`digits` must be NULL or a whole number between 0 and 10")
  }
  if (!is.null(digits)) digits <- as.integer(digits)
  dict <- .ig_cls_dictionary()$metrics
  unit <- if ("unit" %in% names(x)) as.character(x$unit) else rep(NA_character_, nrow(x))
  miss <- is.na(unit)
  unit[miss] <- vapply(as.character(x$metric_id[miss]), function(id) {
    e <- if (is.na(id)) NULL else dict[[id]]
    if (is.null(e)) NA_character_ else e$unit
  }, character(1))
  kind_col <- if ("estimate_kind" %in% names(x)) as.character(x$estimate_kind) else rep(NA_character_, nrow(x))
  reason_code <- if ("reason_code" %in% names(x)) as.character(x$reason_code) else rep(NA_character_, nrow(x))
  reason_txt <- .ig_cls_reason_text(reason_code, lang)
  if ("reason" %in% names(x)) {
    fill <- is.na(reason_txt) & !is.na(x$reason)
    reason_txt[fill] <- as.character(x$reason[fill])
  }
  display <- character(nrow(x))
  display_reason <- rep(NA_character_, nrow(x))
  for (i in seq_len(nrow(x))) {
    st <- as.character(x$status[i])
    v <- x$value[i]
    if (identical(st, "ok") && is.finite(v)) {
      u <- if (is.na(unit[i])) "" else unit[i]
      kind <- switch(u,
        proportion = if (percent) "percent" else "fixed",
        count = if (kind_col[i] %in% .ig_cls_expected_kinds) "expected_count" else "count",
        ratio = "ratio",
        per_1000 = "per_1000",
        "fixed")
      display[i] <- .ig_cls_format_number(v, kind, lang, digits)
    } else if (identical(st, "positive_infinity")) {
      display[i] <- "\u221e"
      display_reason[i] <- if (is.na(reason_txt[i])) .ig_cls_status_text(st, lang) else reason_txt[i]
    } else {
      display[i] <- "\u2013"
      display_reason[i] <- if (is.na(reason_txt[i])) .ig_cls_status_text(st, lang) else reason_txt[i]
    }
  }
  x$display <- display
  x$display_reason <- display_reason
  x
}

# Text of a metric status from the dictionary (fallback for rows without a
# reason code).
.ig_cls_status_text <- function(status, lang = "en") {
  rec <- .ig_cls_dictionary()$statuses[[status]]
  if (is.null(rec)) status else rec[[lang]] %||% rec$en
}
