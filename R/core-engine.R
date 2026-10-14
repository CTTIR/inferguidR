# Generic entry points ---------------------------------------------------------------
#
# calculate_power() and solve_design() route a validated specification to
# the generic mode solver for its mode (phase 05). Until a method has an
# engine, or a mode has a solver, the answer is a power_result with status
# "unsupported" and no number. Solver errors become "numerical_failure"
# results; a solver that returns anything but a matching power_result is
# treated as a failure too.

#' Calculate a power analysis
#'
#' Solves the specification for its declared target: sample size (a priori),
#' power (post hoc), effect size (sensitivity), alpha and critical values
#' (criterion) or alpha and power jointly (compromise).
#'
#' The result always carries a typed status. Methods without an implemented
#' engine return status `"unsupported"` with the reason
#' `"engine for <method> not implemented yet"`; such results contain no
#' numbers.
#'
#' @param spec A `power_spec` created by [analysis_spec()].
#' @return A `power_result`.
#' @examples
#' spec <- analysis_spec("t_means_one_sample", "post_hoc",
#'   inputs = list(alpha = 0.05, n_total = 20, d = 0.5),
#'   options = list(tails = "two"))
#' res <- calculate_power(spec)
#' res$status
#' @seealso [solve_design()], [analysis_spec()]
#' @export
calculate_power <- function(spec) .ig_run(spec, "calculate_power")

#' Solve for a study design
#'
#' Solves an a priori specification for the smallest admissible sample size
#' (and group sizes) that reaches the target power. Specifications in other
#' modes are answered with status `"invalid_input"`; use [calculate_power()]
#' for them.
#'
#' @inheritParams calculate_power
#' @return A `power_result`.
#' @examples
#' spec <- analysis_spec("t_means_one_sample", "a_priori",
#'   inputs = list(alpha = 0.05, power = 0.8, d = 0.5),
#'   options = list(tails = "two"))
#' solve_design(spec)$status
#' @seealso [calculate_power()]
#' @export
solve_design <- function(spec) .ig_run(spec, "solve_design")

.ig_run <- function(spec, entry) {
  .ig_check_spec_integrity(spec)
  if (entry == "solve_design" && !identical(spec$mode, "a_priori")) {
    return(.ig_power_result(spec, ig_status(
      "invalid_input",
      sprintf("solve_design() solves a priori specifications for the sample size; this specification is in mode '%s' (use calculate_power())",
              spec$mode)), provenance = list(entry = entry)))
  }
  def <- .ig_get_method(spec$method_id)
  if (is.null(def)) {
    return(.ig_unsupported_result(spec, sprintf("engine for %s not implemented yet", spec$method_id),
                                  entry))
  }
  solver <- .ig_solver_for(spec$mode)
  if (!is.function(solver)) {
    return(.ig_unsupported_result(spec, sprintf("solver for mode '%s' not implemented yet", spec$mode),
                                  entry))
  }
  caught <- character()
  res <- withCallingHandlers(
    tryCatch(solver(def, spec, entry), error = function(e) e),
    warning = function(w) {
      caught <<- c(caught, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (inherits(res, "error")) {
    return(.ig_power_result(spec, ig_status(
      "numerical_failure", paste("the solver stopped with an error:", conditionMessage(res)),
      list(condition_class = class(res)[1])), provenance = list(entry = entry),
      warnings = caught))
  }
  if (!inherits(res, "power_result") || !identical(res$spec_hash, spec$spec_hash) ||
      !identical(.ig_result_hash(res), res$result_hash)) {
    return(.ig_power_result(spec, ig_status(
      "numerical_failure", "the solver returned an object that is not a matching power_result"),
      provenance = list(entry = entry), warnings = caught))
  }
  if (length(caught)) res <- .ig_result_add_warnings(res, caught)
  res
}

# A specification is used only if it is a power_spec whose content still
# matches its hash and whose validation level matches the current registry.
.ig_check_spec_integrity <- function(spec) {
  if (!inherits(spec, "power_spec")) {
    .ig_abort("invalid_input", "`spec` must be a power_spec created by analysis_spec()")
  }
  key <- .ig_spec_key(spec$method_id, spec$mode, spec$target, spec$inputs, spec$options)
  if (!identical(.ig_hash_value(key), spec$spec_hash)) {
    .ig_abort("invalid_input", "the specification was modified after creation (hash mismatch)")
  }
  has_def <- !is.null(.ig_get_method(spec$method_id))
  level <- if (has_def) "method_definition" else "registry_metadata"
  if (!identical(spec$validation_level, level)) {
    .ig_abort("invalid_input", sprintf(
      "the specification was validated against %s, but %s; re-create it with analysis_spec()",
      gsub("_", " ", spec$validation_level),
      if (has_def) "a method definition is now registered" else "no method definition is registered"))
  }
  invisible(TRUE)
}
