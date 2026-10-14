# Distribution calculator: evaluation and result object ----------------------------------
#
# Ownership: phase 10 utilities (group G12). Method notes:
# specs/method-notes/utilities.md, section "Distribution calculator".
#
# distribution_calculator() parses a script (R/calc-parse.R) and evaluates
# its syntax tree with the allowlisted function table (R/calc-functions.R).
# Variables live in a plain named list; no R environment, eval() or parse()
# is involved, so pasted text can only ever compute numbers.
#
# Outcomes are typed: a syntax error, an unknown name, a wrong number of
# arguments or a parameter outside a function's domain gives status
# invalid_input with the reason and no value (G*Power shows a message box and
# a placeholder 0); a valid parameter beyond the range the calculator computes
# (R/calc-functions.R, .ig_calc_range) gives unsupported; a failed numerical
# routine gives numerical_failure.
# IEEE results of elementary operations (Inf, -Inf, NaN) are values, as in
# G*Power, and are listed in `notes`.

.ig_calc_error <- function(kind, message, position = NA_integer_) {
  stop(structure(class = c("ig_calc_eval_error", "error", "condition"),
                 list(message = message, call = NULL, kind = kind, position = position)))
}

# Evaluate one syntax-tree node. `ctx` is an environment holding the
# variables (a named list), the function table and collected notes; it is
# internal state of one calculator run, never an evaluation environment.
.ig_calc_eval_node <- function(node, ctx) {
  switch(node$type,
    num = node$value,
    group = .ig_calc_eval_node(node$arg, ctx),
    neg = -.ig_calc_eval_node(node$arg, ctx),
    var = {
      if (!is.null(ctx$table[[node$name]])) {
        .ig_calc_error("syntax", sprintf("function %s must be called with ()", node$name), node$pos)
      }
      v <- ctx$vars[[node$name]]
      if (is.null(v)) .ig_calc_error("syntax", sprintf("unknown variable %s", node$name), node$pos)
      v
    },
    bin = {
      a <- .ig_calc_eval_node(node$lhs, ctx)
      b <- .ig_calc_eval_node(node$rhs, ctx)
      switch(node$op, "+" = a + b, "-" = a - b, "*" = a * b, "/" = a / b, "^" = a^b)
    },
    call = {
      fn <- ctx$table[[node$fun]]
      if (is.null(fn)) .ig_calc_error("syntax", sprintf("unknown function %s", node$fun), node$pos)
      k <- length(fn$args)
      if (length(node$args) != k) {
        .ig_calc_error("syntax", sprintf("function %s takes %d argument%s (%s), %d given", node$fun, k,
                                         if (k == 1L) "" else "s", paste(fn$args, collapse = ", "),
                                         length(node$args)), node$pos)
      }
      args <- vapply(node$args, .ig_calc_eval_node, numeric(1), ctx = ctx)
      ctx$used <- union(ctx$used, node$fun)
      if (anyNA(args) && fn$category == "distribution") {
        .ig_calc_error("math", sprintf("%s: an argument is not a number (NaN)", node$fun), node$pos)
      }
      if (is.function(fn$domain)) {
        bad <- fn$domain(args)
        if (!is.null(bad)) {
          if (!startsWith(bad, paste0(node$fun, ":"))) bad <- paste0(node$fun, ": ", bad)
          .ig_calc_error("math", bad, node$pos)
        }
      }
      if (is.function(fn$range)) {
        bad <- fn$range(args)
        if (!is.null(bad)) .ig_calc_error("range", bad, node$pos)
      }
      v <- tryCatch(suppressWarnings(fn$impl(args)),
                    inferguidR_error = function(e) {
                      kind <- switch(e$status$code, invalid_input = "math", unsupported = "range", "numerical")
                      .ig_calc_error(kind, sprintf("%s: %s", node$fun, conditionMessage(e)), node$pos)
                    },
                    error = function(e) {
                      # A typed calculator error raised by the handler above
                      # passes through (a later handler of the same tryCatch()
                      # also sees conditions signalled by an earlier one). Any
                      # other failure of a numerical routine (quadrature,
                      # optimisation, root finding) is typed, never an R error.
                      if (inherits(e, "ig_calc_eval_error")) stop(e)
                      .ig_calc_error("numerical", sprintf("%s: %s", node$fun, conditionMessage(e)), node$pos)
                    })
      if (!is.numeric(v) || length(v) != 1L) {
        .ig_calc_error("numerical", sprintf("%s did not return one number", node$fun), node$pos)
      }
      if (is.na(v) && !is.nan(v)) v <- NaN
      if (is.nan(v) && fn$category == "distribution") {
        # A distribution function has a value for every argument that passed
        # the domain checks; NaN here means the numerical routine gave up.
        .ig_calc_error("numerical", sprintf("%s(%s) could not be computed", node$fun,
                                            paste(vapply(args, .ig_calc_display, ""), collapse = ",")), node$pos)
      }
      if (!is.finite(v) && all(is.finite(args))) {
        ctx$notes <- c(ctx$notes, sprintf("%s(%s) = %s", node$fun,
                                          paste(vapply(args, .ig_calc_display, ""), collapse = ","),
                                          .ig_calc_display(v)))
      }
      as.double(v)
    },
    .ig_calc_error("syntax", sprintf("unknown syntax node '%s'", node$type))
  )
}

#' Distribution calculator
#'
#' Evaluates a script in the syntax of the 'G\*Power' 3.1.9.7 calculator:
#' arithmetic with `+ - * / ^`, assignments, several statements separated by
#' `;` or line breaks, `#` comments, and the mathematical and distribution
#' functions listed by [calculator_functions()] (central and noncentral t,
#' F and chi-square, normal, beta, binomial, Poisson, hypergeometric, sample
#' correlation, squared multiple correlation, log-normal, Laplace,
#' exponential and uniform CDFs, densities and quantiles). The result is the
#' value of the last statement.
#'
#' The script is never evaluated as R code. It is parsed by an allowlisted
#' parser; only numbers, variables, the five operators and the listed
#' functions exist.
#'
#' Syntax conventions follow 'G\*Power' (captured from the running program):
#' `^` is left-associative (`2^3^2` is 64) and a leading minus binds tighter
#' than `^` (`-2^2` is 4, `2^-1` is 0.5); there is no unary `+`; names are
#' letters and digits and are case-sensitive; there are no predefined
#' constants; function arguments are separated by `,`. Where 'G\*Power'
#' computes something other than the documented function (for example
#' `chi2pdf`, `ncx2pdf` and `unipdf`), the documented function is computed;
#' see the `gpower_note` column of [calculator_functions()].
#'
#' Large sine, cosine and tangent arguments use a bundled high-precision
#' argument reduction before evaluating a small angle, avoiding dependence
#' on platform-specific large-argument reduction in the system math library.
#'
#' @param expression A single string: the calculator script.
#' @return An object of class `ig_calc_result`: a list with `expression`,
#'   `status` (an [ig_status()]: `ok`; `invalid_input` for syntax and math
#'   errors; `unsupported` for valid parameters beyond the range the
#'   calculator computes, for example a noncentrality parameter above 1e8;
#'   or `numerical_failure`), `value` (the result, `NA` unless the
#'   status is `ok`; `Inf`, `-Inf` and `NaN` are results of IEEE arithmetic
#'   as in 'G\*Power'), `display` (the value as the 'G\*Power' result field
#'   shows it: 15 significant digits, `inf`, `-inf`, `-nan(ind)`),
#'   `statements` (a data frame with one row per statement), `variables`,
#'   `functions_used`, `notes` (non-finite intermediate results and
#'   reference notes of the functions used) and `r_call`.
#' @examples
#' distribution_calculator("1 - normcdf(1.96, 0, 1)")
#' distribution_calculator("N = 16; tc = tinv(0.95, N - 1); 1 - nctcdf(tc, N - 1, 2)")$value
#' distribution_calculator("2^3^2")$value   # 64: '^' is left-associative
#' @seealso [calculator_functions()]
#' @export
distribution_calculator <- function(expression) {
  if (!is.character(expression) || length(expression) != 1L || is.na(expression)) {
    .ig_abort_invalid("`expression` must be a single string")
  }
  text <- enc2utf8(expression)
  r_call <- sprintf("inferguidR::distribution_calculator(%s)", .ig_r_string(text))
  parsed <- .ig_calc_parse(text)
  if (!isTRUE(parsed$ok)) {
    return(.ig_new_calc_result(text, ig_status("invalid_input", paste0("syntax error: ", parsed$message),
                                               list(kind = parsed$kind, position = parsed$position)),
                               r_call = r_call))
  }
  ctx <- new.env(parent = emptyenv())
  ctx$vars <- list()
  ctx$table <- .ig_calc_table()
  ctx$notes <- character()
  ctx$used <- character()
  stmts <- parsed$statements
  rows <- list()
  value <- NULL
  for (i in seq_along(stmts)) {
    s <- stmts[[i]]
    if (identical(s$kind, "empty")) next
    out <- tryCatch({
      if (identical(s$kind, "assign") && !is.null(ctx$table[[s$name]])) {
        .ig_calc_error("syntax", sprintf("function %s must be called with ()", s$name), s$pos)
      }
      .ig_calc_eval_node(s$tree, ctx)
    }, ig_calc_eval_error = function(e) e)
    if (inherits(out, "ig_calc_eval_error")) {
      code <- switch(out$kind, numerical = "numerical_failure", range = "unsupported", "invalid_input")
      prefix <- switch(out$kind, syntax = "syntax error: ", math = "math error: ", numerical = "numerical failure: ",
                       range = "outside the computable range: ")
      return(.ig_new_calc_result(text, ig_status(code, paste0(prefix, conditionMessage(out)),
                                                 list(kind = out$kind, position = out$position, statement = i)),
                                 statements = rows, variables = ctx$vars, used = ctx$used, notes = ctx$notes,
                                 r_call = r_call))
    }
    if (identical(s$kind, "assign")) ctx$vars[[s$name]] <- out
    rows[[length(rows) + 1L]] <- list(index = i, kind = s$kind, name = s$name,
                                      text = .ig_calc_deparse(s$tree), value = out)
    value <- out
  }
  if (is.null(value)) {
    return(.ig_new_calc_result(text, ig_status("invalid_input", paste(
      "the script contains no expression (G*Power displays 0 for such a script; here there is no result)"),
      list(kind = "syntax")), r_call = r_call))
  }
  extra <- character()
  if (grepl("[0-9.][eE][+][0-9]", text)) {
    # Captured live (SUP-G12-UT-CALC-F05, -F07, -F08, -P11): G*Power 3.1.9.7
    # reads the exponent digits after 'e+' a second time as '+digits'.
    extra <- paste("reference defect: G*Power 3.1.9.7 evaluates a number with an explicit '+' in the",
                   "exponent differently (1e+3 is 1003 and 2.5e+2 is 252 there); here 1e+3 = 1000")
  }
  .ig_new_calc_result(text, ig_status("ok"), value = value, statements = rows, variables = ctx$vars,
                      used = ctx$used, notes = ctx$notes, extra_notes = extra, r_call = r_call)
}

.ig_new_calc_result <- function(text, status, value = NA_real_, statements = list(), variables = list(),
                                used = character(), notes = character(), extra_notes = character(),
                                r_call = NA_character_) {
  tab <- .ig_calc_table()
  ok <- identical(status$code, "ok")
  ref_notes <- unlist(lapply(used, function(f) {
    g <- tab[[f]]$gpower_note
    if (is.null(g) || is.na(g)) NULL else sprintf("%s: %s", f, g)
  }))
  stmt_df <- data.frame(
    index = vapply(statements, function(r) as.integer(r$index), integer(1)),
    kind = vapply(statements, `[[`, "", "kind"),
    name = vapply(statements, function(r) as.character(r$name), ""),
    text = vapply(statements, `[[`, "", "text"),
    value = vapply(statements, `[[`, numeric(1), "value"),
    display = vapply(statements, function(r) .ig_calc_display(r$value), ""),
    stringsAsFactors = FALSE
  )
  structure(list(
    schema_version = "1.0.0",
    object_type = "calculator_result",
    expression = text,
    status = status,
    value = if (ok) as.double(value) else NA_real_,
    display = if (ok) .ig_calc_display(value) else NA_character_,
    statements = stmt_df,
    variables = variables,
    functions_used = sort(unique(used)),
    notes = c(if (length(notes)) paste("non-finite result:", unique(notes)), extra_notes, ref_notes),
    r_call = r_call
  ), class = "ig_calc_result")
}

# The value as the G*Power result field shows it: printf("%.15g") of the
# Microsoft C runtime; non-finite values as that runtime prints them (the
# default quiet NaN of an invalid operation prints as -nan(ind)).
.ig_calc_display <- function(x) {
  if (length(x) != 1L || is.na(x) && !is.nan(x)) return(NA_character_)
  if (is.nan(x)) return("-nan(ind)")
  if (is.infinite(x)) return(if (x > 0) "inf" else "-inf")
  sprintf("%.15g", x)
}

# A string literal for a reproducible R call.
.ig_r_string <- function(x) {
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub("\"", "\\\"", x, fixed = TRUE)
  x <- gsub("\n", "\\n", x, fixed = TRUE)
  x <- gsub("\r", "\\r", x, fixed = TRUE)
  x <- gsub("\t", "\\t", x, fixed = TRUE)
  paste0("\"", x, "\"")
}

#' @export
format.ig_calc_result <- function(x, ...) {
  head <- sprintf("<ig_calc_result> %s", format(x$status))
  body <- if (identical(x$status$code, "ok")) sprintf("  result: %s", x$display) else character()
  c(head, body, if (length(x$notes)) paste0("  note: ", x$notes))
}

#' @export
print.ig_calc_result <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}
