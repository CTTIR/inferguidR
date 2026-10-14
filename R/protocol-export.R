# Protocol export: a text record of analyses in the layout of the G*Power protocol --------
#
# Ownership: phase 10 utilities (group G12). Method notes:
# specs/method-notes/utilities.md, section "Protocol export".
#
# 'G*Power' 3.1.9.7 appends one record per calculation to its protocol window
# (tab 'Protocol of power analyses'; specs/gpower-utilities.json
# util_protocol, certified live). A record, as captured from the running
# program (validation/raw/emulator/cases/util_protocol/), is
#
#   [1] -- Monday, October 05, 2026 -- 13:57:49
#   t tests - Means: Difference between two independent means (two groups)
#   Options:<TAB><option texts>                       (only for tests with options)
#   Analysis:<TAB>Post hoc: Compute achieved power <-- analysis type up to "- given", trailing blank kept
#   Input:<TAB>Tail(s)<TAB>=<TAB>Two
#   <TAB>Effect size d<TAB>=<TAB>0.5
#   ...
#   Output:<TAB>Noncentrality parameter delta<TAB>=<TAB>2.5000000
#   <TAB>Critical t<TAB>=<TAB>1.9844675
#   ...
#   <empty line>
#
# with CR LF line ends. protocol_record() builds the record of a power_result
# from the live input and output labels of its method and analysis type
# (registry inputs_by_mode_live / outputs_by_mode_live); protocol_text()
# formats one or more records. Numbers are formatted as G*Power displays them
# (.ig_gpower_number(); rule derived from 39,523 real-valued displays of the
# emulator-capture corpus, method notes). An output label that inferguidR
# cannot fill is printed as "(not reported by inferguidR)" and listed in the
# record; no value is invented.
#
# A method can supply its own label mapping with a function
# `.ig_protocol_outputs_<method_id>(result)` returning a named list keyed by
# live output labels (and `.ig_protocol_options_<method_id>(spec)` for the
# text of the 'Options:' line); otherwise the generic mapping below is used.

# Live test-family texts of the G*Power test list (specs/gpower-app-capture.json,
# procedures.*.test_family_live; checked by tests/testthat/test-protocol-export.R).
.ig_protocol_family_text <- c(t = "t tests", F = "F tests", chi2 = "\u03c7\u00b2 tests", z = "z tests", exact = "Exact")

# Analysis-type texts (specs/gpower-app-capture.json analysis_types_live): the
# five standard types and the four types of the generic tests.
.ig_protocol_analysis_live <- list(
  standard = c(
    a_priori = "A priori: Compute required sample size - given \u03b1, power, and effect size",
    compromise = "Compromise: Compute implied \u03b1 & power - given \u03b2/\u03b1 ratio, sample size, and effect size",
    criterion = "Criterion: Compute required \u03b1 - given power, effect size, and sample size",
    post_hoc = "Post hoc: Compute achieved power - given \u03b1, sample size, and effect size",
    sensitivity = "Sensitivity: Compute required effect size - given \u03b1, power, and sample size"),
  generic = c(
    compromise = "Compromise: Compute \u03b1 and \u03b2 - given \u03b2/\u03b1 ratio, and noncentrality parameter",
    criterion = "Criterion: Compute alpha - given \u03b2, and noncentrality parameter",
    post_hoc = "Post hoc: Compute power - given \u03b1, and noncentrality parameter",
    sensitivity = "Sensitivity: Compute noncentrality parameter - given \u03b1, and power")
)

.ig_protocol_analysis_type <- function(method_id, mode) {
  set <- if (grepl("_generic$", method_id)) "generic" else "standard"
  unname(.ig_protocol_analysis_live[[set]][mode])
}

# The protocol writes the analysis type up to "- given" and keeps the blank
# before it ("Post hoc: Compute achieved power ").
.ig_protocol_analysis_text <- function(method_id, mode) {
  full <- .ig_protocol_analysis_type(method_id, mode)
  if (is.na(full)) return(NA_character_)
  sub("- given.*$", "", full)
}

# Output labels shown as integers (counts and integral degrees of freedom).
# The repeated-measures procedures show their degrees of freedom as reals
# ("19.0000000"; every capture of these methods, protocol corpus comparison
# of validation/compare/utilities/util_protocol), because a nonsphericity
# correction can make them fractional.
.ig_protocol_real_df_methods <- c("f_rm_anova_between", "f_rm_anova_within", "f_rm_anova_within_between",
                                  "f_rm_manova_between", "f_rm_manova_within", "f_rm_manova_within_between")
.ig_protocol_integer_label <- function(label, method_id = NULL) {
  if (!is.null(method_id) && method_id %in% .ig_protocol_real_df_methods && grepl("df$", label)) return(FALSE)
  grepl("^(Df|Numerator df|Denominator df|Total sample size|Sample size( group [0-9]+)?|Lower critical N|Upper critical N|Number of .*)$",
        label)
}

# A real number as the G*Power output fields show it: printf("%.7g") when
# |x| >= 100 or |x| < 0.01 (zero included), printf("%.7f") otherwise;
# non-finite values as the Microsoft C runtime prints them.
.ig_gpower_number <- function(x) {
  if (length(x) != 1L || is.na(x) && !is.nan(x)) return(NA_character_)
  if (is.nan(x)) return("-nan(ind)")
  if (is.infinite(x)) return(if (x > 0) "inf" else "-inf")
  a <- abs(x)
  if (a >= 100 || a < 0.01) sprintf("%.7g", x) else sprintf("%.7f", x)
}

# An input value as typed into a G*Power field: the shortest decimal form
# with up to 15 significant digits.
.ig_protocol_input_number <- function(x) {
  if (is.logical(x)) return(if (isTRUE(x)) "TRUE" else "FALSE")
  if (!is.numeric(x)) return(as.character(x))
  if (length(x) != 1L) return(paste(vapply(x, .ig_protocol_input_number, ""), collapse = ", "))
  if (is.na(x)) return("NA")
  if (is.infinite(x)) return(if (x > 0) "inf" else "-inf")
  sprintf("%.15g", x)
}

# Generic mapping of a result to live output labels.
.ig_protocol_generic_outputs <- function(res, labels, def) {
  q <- res$quantities
  a <- res$achieved
  s <- res$solved
  gs <- res$group_sizes
  ncp <- q$ncp
  crit <- q$critical
  effect_par <- NULL
  if (!is.null(def)) {
    for (p in def$params) if (identical(p$role, "effect") && .ig_is_string(p$label$en)) effect_par[[p$label$en]] <- p$name
  }
  # Method outputs whose name equals the label up to case and punctuation
  # (for example wilks_u for "Wilks U", pillai_v for "Pillai V").
  norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  mo_num <- Filter(function(v) is.numeric(v) && length(v) == 1L, res$method_outputs %||% list())
  mo_norm <- if (length(mo_num)) stats::setNames(names(mo_num), norm(names(mo_num))) else character()
  out <- list()
  for (lab in labels) {
    v <- NULL
    mo <- res$method_outputs[[lab]]
    if (is.numeric(mo) && length(mo) == 1L) {
      v <- mo
    } else if (!is.na(mo_norm[norm(lab)]) && nzchar(norm(lab))) {
      v <- mo_num[[mo_norm[[norm(lab)]]]]
    } else if (grepl("^Noncentrality parameter ", lab)) {
      v <- ncp
    } else if (grepl("^Lower critical ", lab) && length(crit) && !grepl(" N$", lab)) {
      v <- min(crit)
    } else if (grepl("^Upper critical ", lab) && length(crit) && !grepl(" N$", lab)) {
      v <- max(crit)
    } else if (grepl("^Critical ", lab) && length(crit)) {
      # One critical value is shown; for two tails it is the one on the side
      # of the noncentrality parameter (captured for the t tests, mapping of
      # validation/compare/t_means_two_independent).
      v <- if (length(crit) == 1L) crit else if (isTRUE(ncp < 0)) min(crit) else max(crit)
    } else if (identical(lab, "Df") && length(q$df) == 1L) {
      v <- q$df
    } else if (identical(lab, "Numerator df") && length(q$df) == 2L) {
      v <- q$df[1]
    } else if (identical(lab, "Denominator df") && length(q$df) == 2L) {
      v <- q$df[2]
    } else if (grepl("^Sample size group [0-9]+$", lab) && length(gs)) {
      k <- as.integer(sub("^Sample size group ", "", lab))
      nm <- paste0("n", k)
      v <- if (nm %in% names(gs)) gs[[nm]] else if (k <= length(gs)) gs[[k]] else NULL
    } else if (lab %in% c("Total sample size", "Sample size")) {
      v <- s$n_total %||% q$n_total
    } else if (lab %in% c("Actual power", "Power (1-\u03b2 err prob)")) {
      v <- a$power
    } else if (identical(lab, "\u03b1 err prob")) {
      v <- s$alpha
    } else if (identical(lab, "\u03b2 err prob")) {
      v <- a$beta
    } else if (identical(lab, "Actual \u03b1")) {
      v <- a$alpha_attained
    } else if (identical(lab, "Actual \u03b2/\u03b1 ratio")) {
      v <- a$error_ratio
    } else if (!is.null(effect_par[[lab]])) {
      v <- s[[effect_par[[lab]]]]
    }
    if (!is.null(v)) out[[lab]] <- unname(v)
  }
  out
}

# Generic 'Options:' text: the live texts of the selected values of the
# result-affecting options that are set in the Options dialog (not main-window
# choices, which are inputs), joined by ", ". Methods whose G*Power wording
# differs supply .ig_protocol_options_<method_id>().
.ig_protocol_generic_options <- function(spec, meta) {
  txt <- character()
  for (o in meta$options) {
    if (identical(o$kind, "main_window_input")) next
    code <- spec$options[[o$option_id]]
    if (is.null(code)) next
    lv <- o$value_live[[as.character(code)]]
    t <- if (!is.null(lv$text)) lv$text else if (is.numeric(code)) sprintf("%s = %s", o$ui_label, .ig_protocol_input_number(code)) else NULL
    # Dialog texts carry runs of blanks that align the controls.
    if (!is.null(t) && nzchar(t)) txt <- c(txt, trimws(gsub("[[:space:]]+", " ", t)))
  }
  if (length(txt)) paste(txt, collapse = ", ") else NULL
}

#' Protocol record of an analysis
#'
#' Builds the record that the 'G\*Power' protocol window shows for a
#' calculation: test family and procedure, analysis type, the inputs and
#' the outputs with their 'G\*Power' labels and display formatting. The
#' labels are the live labels of the method and analysis type in the method
#' registry. Use [protocol_text()] to format records.
#'
#' @param result A `power_result` (see [calculate_power()]).
#' @return An object of class `ig_protocol_record`: a list with `method_id`,
#'   `mode`, `status`, `family_line`, `options_line` (`NULL` without
#'   options), `analysis_line`, `inputs` and `outputs` (data frames with
#'   `label`, `text` and `value`), `unresolved` (labels inferguidR could not
#'   fill), `notes` (where the record may differ from G\*Power's layout),
#'   `validation` (validation status of the method), `spec_hash`,
#'   `result_hash` and `r_call`.
#' @examples
#' res <- calculate_power(analysis_spec("t_means_two_independent", "post_hoc",
#'   inputs = list(d = 0.5, alpha = 0.05, n1 = 50, n2 = 50), options = list(tails = "two")))
#' rec <- protocol_record(res)
#' rec$outputs
#' @seealso [protocol_text()]
#' @export
protocol_record <- function(result) {
  if (!inherits(result, "power_result")) .ig_abort_invalid("`result` must be a power_result")
  id <- result$method_id
  mode <- result$mode
  spec <- result$spec
  meta <- .ig_get_metadata(id)
  def <- .ig_get_method(id)
  live <- tryCatch(.ig_live_labels(id, mode), inferguidR_error = function(e) list(inputs = NULL, outputs = NULL))
  family_text <- unname(.ig_protocol_family_text[meta$family])
  if (is.na(family_text)) family_text <- meta$family
  test_label <- meta$label_raw$windows_3197_live %||% meta$label$en
  family_line <- paste0(family_text, " - ", test_label)
  opt_fun <- get0(paste0(".ig_protocol_options_", id), envir = topenv(environment(protocol_record)), inherits = FALSE)
  options_line <- if (is.function(opt_fun)) opt_fun(spec) else .ig_protocol_generic_options(spec, meta)
  generic_options <- !is.function(opt_fun) && !is.null(options_line)
  analysis <- .ig_protocol_analysis_text(id, mode)

  # Inputs: live labels in live order; a label is an option of the main
  # window (choice: its live text; continuous: its value) or a parameter.
  in_labels <- live$inputs
  if (is.null(in_labels)) {
    in_labels <- unlist(lapply(def$params, function(p) {
      if (p$name %in% names(spec$inputs)) p$label$en else NULL
    }))
  }
  unresolved <- character()
  notes <- character()
  if (generic_options) {
    # Corpus comparison (validation/compare/utilities/util_protocol): G*Power
    # words this line itself ('A.R.E. method', 'Exact distribution', ...);
    # those texts are not in the registry.
    notes <- c(notes, paste("the 'Options:' line lists the chosen options by their dialog texts; G*Power writes its own",
                            "wording there"))
  }
  inputs <- data.frame(label = character(), text = character(), value = numeric(), stringsAsFactors = FALSE)
  shown <- character()
  for (lab in in_labels) {
    txt <- NULL
    val <- NA_real_
    opt <- Filter(function(o) identical(o$ui_label, lab) || identical(o$live$label, lab), meta$options)
    par <- Filter(function(p) identical(p$label$en, lab), def$params %||% list())
    if (!length(par) && !length(opt)) {
      # A live label that shortens a parameter label ("Total sample size" for
      # "Total sample size (number of pairs)") is that parameter when exactly
      # one parameter label starts with it. The reverse (a live label with an
      # added suffix such as "(V)") is not resolved this way: the suffix names
      # an option variant.
      cand <- Filter(function(p) .ig_is_string(p$label$en) && startsWith(p$label$en, paste0(lab, " ")), def$params %||% list())
      if (length(cand) == 1L) par <- cand
    }
    if (length(par) && !(par[[1]]$name %in% names(spec$inputs))) {
      # A field of the default layout that this specification does not use
      # (another option variant, or the solve target): not written.
      next
    }
    if (length(par)) shown <- c(shown, par[[1]]$name)
    if (length(par) && par[[1]]$name %in% names(spec$inputs)) {
      v <- spec$inputs[[par[[1]]$name]]
      txt <- .ig_protocol_input_number(v)
      if (is.numeric(v) && length(v) == 1L) val <- as.double(v)
    } else if (length(opt)) {
      o <- opt[[1]]
      code <- spec$options[[o$option_id]] %||% spec$inputs[[o$option_id]]
      if (!is.null(code)) {
        lv <- o$value_live[[as.character(code)]]
        txt <- if (!is.null(lv$text)) lv$text else .ig_protocol_input_number(code)
        if (is.numeric(code) && length(code) == 1L) val <- as.double(code)
      }
    }
    if (is.null(txt)) {
      unresolved <- c(unresolved, paste0("input: ", lab))
      txt <- "(not reported by inferguidR)"
    }
    inputs[nrow(inputs) + 1L, ] <- list(lab, txt, val)
  }
  # Inputs of the specification that no live label of the default layout
  # covers (an option variant with other fields) are written with their
  # parameter labels, so that no input is hidden.
  for (nm in setdiff(names(spec$inputs), shown)) {
    p <- def$params[[nm]]
    if (is.null(p)) next
    lab <- p$label$en %||% nm
    if (lab %in% inputs$label) next
    v <- spec$inputs[[nm]]
    inputs[nrow(inputs) + 1L, ] <- list(lab, .ig_protocol_input_number(v), if (is.numeric(v) && length(v) == 1L) as.double(v) else NA_real_)
    notes <- c(notes, sprintf("input '%s' is not a field of G*Power's default layout for this analysis type", lab))
  }
  # Options that change labels: the registry lists the default layout only.
  dialog <- Filter(function(o) !identical(o$kind, "main_window_input") && isTRUE(o$result_affecting), meta$options)
  non_default <- vapply(dialog, function(o) {
    v <- spec$options[[o$option_id]]
    !is.null(v) && !is.null(o$default) && !identical(as.character(v), as.character(o$default))
  }, logical(1))
  if (any(non_default)) {
    notes <- c(notes, paste("labels follow G*Power's default layout; for the options chosen here G*Power may show other",
                            "labels or option texts (option-variant layouts are not in the registry yet)"))
  }

  # Outputs.
  out_labels <- live$outputs %||% character()
  ok <- identical(result$status$code, "ok")
  outputs <- data.frame(label = character(), text = character(), value = numeric(), stringsAsFactors = FALSE)
  if (ok) {
    hook <- get0(paste0(".ig_protocol_outputs_", id), envir = topenv(environment(protocol_record)), inherits = FALSE)
    vals <- if (is.function(hook)) hook(result) else .ig_protocol_generic_outputs(result, out_labels, def)
    for (lab in out_labels) {
      v <- vals[[lab]]
      if (is.null(v) || length(v) != 1L || (is.na(v) && !is.nan(v))) {
        unresolved <- c(unresolved, paste0("output: ", lab))
        outputs[nrow(outputs) + 1L, ] <- list(lab, "(not reported by inferguidR)", NA_real_)
        next
      }
      v <- as.double(v)
      txt <- if (.ig_protocol_integer_label(lab, id) && is.finite(v) && v == round(v)) sprintf("%.0f", v) else .ig_gpower_number(v)
      outputs[nrow(outputs) + 1L, ] <- list(lab, txt, v)
    }
  }
  structure(list(
    method_id = id, mode = mode, status = result$status,
    family_line = family_line, options_line = options_line,
    analysis_line = analysis,
    inputs = inputs, outputs = outputs, unresolved = unresolved, notes = notes,
    validation = list(method_status = result$validation$method_validation_status %||% NA_character_,
                      gpower_reference = result$validation$gpower_reference_status %||% NA_character_),
    spec_hash = result$spec_hash, result_hash = result$result_hash, r_call = result$r_call
  ), class = "ig_protocol_record")
}

.ig_protocol_weekdays <- c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday")
.ig_protocol_months <- c("January", "February", "March", "April", "May", "June", "July", "August",
                         "September", "October", "November", "December")

# "[k] -- Monday, October 05, 2026 -- 13:57:49" (English names in every locale).
.ig_protocol_header <- function(k, time) {
  lt <- as.POSIXlt(time)
  sprintf("[%d] -- %s, %s %02d, %04d -- %02d:%02d:%02d", as.integer(k), .ig_protocol_weekdays[lt$wday + 1L],
          .ig_protocol_months[lt$mon + 1L], lt$mday, lt$year + 1900L, lt$hour, lt$min, floor(lt$sec))
}

.ig_protocol_lines <- function(rec, k, time, provenance) {
  tab <- "\t"
  lines <- c(.ig_protocol_header(k, time), rec$family_line)
  if (!is.null(rec$options_line)) lines <- c(lines, paste0("Options:", tab, rec$options_line))
  lines <- c(lines, paste0("Analysis:", tab, rec$analysis_line))
  for (i in seq_len(nrow(rec$inputs))) {
    lines <- c(lines, paste0(if (i == 1L) "Input:" else "", tab, rec$inputs$label[i], tab, "=", tab, rec$inputs$text[i]))
  }
  if (identical(rec$status$code, "ok")) {
    for (i in seq_len(nrow(rec$outputs))) {
      lines <- c(lines, paste0(if (i == 1L) "Output:" else "", tab, rec$outputs$label[i], tab, "=", tab, rec$outputs$text[i]))
    }
  } else {
    lines <- c(lines, paste0("Output:", tab, "(no result: ", rec$status$code, ": ", rec$status$reason, ")"))
  }
  if (provenance) {
    lines <- c(lines,
               paste0("inferguidR:", tab, sprintf("version %s; method status %s; G*Power reference %s",
                                                  .ig_package_version(), rec$validation$method_status,
                                                  rec$validation$gpower_reference)),
               paste0(tab, "spec ", rec$spec_hash, "; result ", rec$result_hash),
               paste0(tab, rec$r_call))
    if (length(rec$unresolved)) lines <- c(lines, paste0(tab, "not reported by inferguidR: ", paste(rec$unresolved, collapse = "; ")))
    if (length(rec$notes)) lines <- c(lines, paste0(tab, "note: ", rec$notes))
  }
  c(lines, "")
}

#' Protocol of power analyses as text
#'
#' Formats one or more analyses as the protocol window of 'G\*Power' 3.1.9.7
#' does: a numbered, time-stamped record per analysis with the test, the
#' analysis type, and the inputs and outputs by their 'G\*Power' labels and
#' display formatting. With `provenance = TRUE` (default) each record ends
#' with the inferguidR version, the validation status of the method, the
#' specification and result hashes and the reproducible R call.
#'
#' @param ... `power_result` objects (or `ig_protocol_record` objects).
#' @param results A list of further results.
#' @param time Time stamp(s) of the records: one `POSIXct` value used for all
#'   records or one per record. Defaults to the current time.
#' @param start Number of the first record.
#' @param provenance Append the inferguidR provenance lines to each record.
#' @param eol Line end: `"\n"` (default) or `"\r\n"` as in the 'G\*Power'
#'   export.
#' @return A single string.
#' @examples
#' res <- calculate_power(analysis_spec("t_means_two_independent", "post_hoc",
#'   inputs = list(d = 0.5, alpha = 0.05, n1 = 50, n2 = 50), options = list(tails = "two")))
#' cat(protocol_text(res, time = as.POSIXct("2026-10-05 13:57:49", tz = "UTC")))
#' @seealso [protocol_record()]
#' @export
protocol_text <- function(..., results = list(), time = Sys.time(), start = 1L, provenance = TRUE, eol = "\n") {
  items <- c(list(...), results)
  if (!length(items)) .ig_abort_invalid("no results to write")
  recs <- lapply(items, function(x) {
    if (inherits(x, "ig_protocol_record")) return(x)
    if (inherits(x, "power_result")) return(protocol_record(x))
    .ig_abort_invalid("every item must be a power_result or an ig_protocol_record")
  })
  if (!inherits(time, "POSIXt") || !(length(time) %in% c(1L, length(recs)))) {
    .ig_abort_invalid("`time` must be one POSIXct value or one per record")
  }
  if (!(identical(eol, "\n") || identical(eol, "\r\n"))) .ig_abort_invalid("`eol` must be \"\\n\" or \"\\r\\n\"")
  if (!.ig_is_flag(provenance)) .ig_abort_invalid("`provenance` must be TRUE or FALSE")
  if (!is.numeric(start) || length(start) != 1L || start < 1 || start != round(start)) {
    .ig_abort_invalid("`start` must be a whole number >= 1")
  }
  time <- rep_len(time, length(recs))
  lines <- unlist(lapply(seq_along(recs), function(i) .ig_protocol_lines(recs[[i]], start + i - 1L, time[i], provenance)))
  paste0(paste(lines, collapse = eol), eol)
}

#' @export
format.ig_protocol_record <- function(x, ...) {
  strsplit(protocol_text(x, time = as.POSIXct("1970-01-01", tz = "UTC"), provenance = TRUE), "\n", fixed = TRUE)[[1]][-1]
}

#' @export
print.ig_protocol_record <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}
