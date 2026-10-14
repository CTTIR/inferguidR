# Allowlisted condition interpreter ---------------------------------------------
#
# Rule conditions are JSON objects in a closed grammar (rules.json,
# engine.condition_grammar):
#
#   {"q": <question id>, "in": [<option ids>]}  effective answer is one of the options
#   {"all": [<conditions>]}                      every condition holds
#   {"any": [<conditions>]}                      at least one condition holds
#   {"not": <condition>}                         negation
#   {"always": true}                             always holds
#
# Conditions are checked against this grammar and turned into closures built
# from the five fixed forms above. No rule text is ever passed to parse() or
# eval(); a string, an unknown key or an extra key is a structural error.
#
# Effective answers are a named character vector (question id -> option id).
# A question that was not answered makes every {"q", "in"} test on it false.

.ig_advisor_condition_keys <- c("q", "in", "all", "any", "not", "always")

# Problems of one condition (character vector, empty if valid). `questions`
# is the named list of question records, used to resolve ids; `path` names
# the location for messages.
.ig_advisor_condition_problems <- function(cond, questions, path) {
  if (is.character(cond)) {
    return(sprintf("%s: conditions must be JSON objects; expression strings are never evaluated", path))
  }
  if (!is.list(cond) || !.ig_is_named(cond)) {
    return(sprintf("%s: a condition must be a JSON object", path))
  }
  keys <- names(cond)
  unknown <- setdiff(keys, .ig_advisor_condition_keys)
  if (length(unknown)) {
    return(sprintf("%s: condition operator '%s' is not allowlisted", path, unknown[1]))
  }
  if (setequal(keys, c("q", "in")) && length(keys) == 2L) {
    q <- cond$q
    opts <- cond[["in"]]
    if (!.ig_is_string(q)) return(sprintf("%s: 'q' must be a question id", path))
    if (!is.list(opts) || !length(opts) ||
        !all(vapply(opts, .ig_is_string, logical(1)))) {
      return(sprintf("%s: 'in' must be a non-empty array of option ids", path))
    }
    qrec <- questions[[q]]
    if (is.null(qrec)) return(sprintf("%s: unknown question id '%s'", path, q))
    known <- vapply(Filter(function(o) .ig_is_string(o$option_id), qrec$options %||% list()),
                    `[[`, "", "option_id")
    bad <- setdiff(unlist(opts), known)
    if (length(bad)) {
      return(sprintf("%s: unknown option id '%s' for question '%s'", path, bad[1], q))
    }
    return(character())
  }
  if (length(keys) != 1L) {
    return(sprintf("%s: a condition must have exactly one operator (or the pair q/in)", path))
  }
  switch(keys,
         all = ,
         any = {
           parts <- cond[[keys]]
           if (!is.list(parts) || !length(parts) || !is.null(names(parts))) {
             return(sprintf("%s: '%s' must be a non-empty array of conditions", path, keys))
           }
           unlist(lapply(seq_along(parts), function(i) {
             .ig_advisor_condition_problems(parts[[i]], questions, sprintf("%s/%s[%d]", path, keys, i))
           }))
         },
         not = .ig_advisor_condition_problems(cond$not, questions, paste0(path, "/not")),
         always = if (isTRUE(cond$always)) character() else
           sprintf("%s: 'always' must be literally true", path),
         sprintf("%s: incomplete condition (q needs in)", path))
}

# Compile a (structurally valid) condition into a closure function(eff).
.ig_advisor_compile_condition <- function(cond) {
  if (!is.list(cond) || !.ig_is_named(cond)) {
    stop("internal error: condition is not a JSON object", call. = FALSE)
  }
  keys <- names(cond)
  if (length(keys) == 2L && setequal(keys, c("q", "in"))) {
    q <- cond$q
    opts <- as.character(unlist(cond[["in"]]))
    return(function(eff) {
      a <- eff[q]
      !is.na(a) && a %in% opts
    })
  }
  switch(keys[1L],
         all = {
           fs <- lapply(cond$all, .ig_advisor_compile_condition)
           function(eff) {
             for (f in fs) if (!f(eff)) return(FALSE)
             TRUE
           }
         },
         any = {
           fs <- lapply(cond$any, .ig_advisor_compile_condition)
           function(eff) {
             for (f in fs) if (f(eff)) return(TRUE)
             FALSE
           }
         },
         not = {
           f <- .ig_advisor_compile_condition(cond$not)
           function(eff) !f(eff)
         },
         always = function(eff) TRUE,
         stop("internal error: condition operator not allowlisted", call. = FALSE))
}

# Question and option ids referenced by a condition (for static checks).
.ig_advisor_condition_refs <- function(cond) {
  if (!is.list(cond)) return(character())
  out <- character()
  if (.ig_is_string(cond$q)) out <- c(out, cond$q, unlist(cond[["in"]]))
  for (k in c("all", "any")) {
    if (is.list(cond[[k]])) out <- c(out, unlist(lapply(cond[[k]], .ig_advisor_condition_refs)))
  }
  if (is.list(cond$not)) out <- c(out, .ig_advisor_condition_refs(cond$not))
  as.character(out)
}
