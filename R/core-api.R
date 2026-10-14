# Public analysis API around the mode solvers -----------------------------------------
#
# calculate_power() and solve_design() (R/core-engine.R) dispatch a
# specification to the generic mode solver of its mode (R/num-modes.R). This
# file adds the functions built on top of them:
#
# * effect_size(): effect-size drawers, looked up as `.ig_drawer_<method_id>()`
#   in the namespace (for example R/effect-t-means.R);
# * sweep_design(): one input varied over a set of values, every row a full
#   specification with its own result;
# * compare_designs(): several specifications or results side by side.
#
# None of them computes statistics itself; every number comes from the
# method engines through calculate_power().

# Effect-size drawers -----------------------------------------------------------------

#' Compute an effect size from raw inputs
#'
#' Effect-size drawers turn descriptive quantities (for example group means
#' and standard deviations) into the effect-size parameter of a method, as
#' the "Determine" dialog of 'G*Power' does. The formula, its source and any
#' caveats are returned with the value.
#'
#' Available drawers and input modes:
#'
#' * `"t_means_two_independent"` (effect `d`):
#'   - `"common_sd"`: `mean_1`, `mean_2`, `sd` (one SD common to both groups),
#'     d = (mean_1 - mean_2) / sd;
#'   - `"two_sds"`: `mean_1`, `mean_2`, `sd_1`, `sd_2`, with the mean SD
#'     sqrt((sd_1^2 + sd_2^2) / 2) (intended for equal group sizes);
#'   - `"point_biserial_r"`: `r`, `n1`, `n2`, d = r / sqrt(1 - r^2) *
#'     (n1 + n2) / sqrt(n1 * n2) (the hand conversion of the 'G*Power' manual,
#'     section 16.3).
#'
#' * `"t_means_one_sample"`: `"from_means"` takes `mean_h0`, `mean_h1`, `sd`.
#' * `"t_means_matched_pairs"`: `"from_differences"` takes `mean_diff`,
#'   `sd_diff`; `"from_group_parameters"` takes `mean_1`, `mean_2`, `sd_1`,
#'   `sd_2`, `rho` (the within-pair correlation).
#' * `"t_wilcoxon_signed_rank_one_sample"`: `"means_sd"` takes `mean_h0`,
#'   `mean_h1`, `sd`.
#' * `"t_wilcoxon_signed_rank_matched_pairs"`: `"differences"` takes
#'   `mean_diff`, `sd_diff`; `"group_parameters"` takes `mean_1`, `mean_2`,
#'   `sd_1`, `sd_2`, `r` (the within-pair correlation).
#' * `"t_wilcoxon_mann_whitney"`: `"common_sd"` takes `mean_1`, `mean_2`,
#'   `sd`; `"two_sds"` takes `mean_1`, `mean_2`, `sd_1`, `sd_2`.
#'
#' Repeated-measures drawers:
#'
#' * `"f_rm_anova_between"`, `"f_rm_anova_within"`,
#'   `"f_rm_anova_within_between"`, `"f_rm_manova_between"` and
#'   `"f_rm_manova_within"` offer `"direct"` (`partial_eta2`),
#'   `"from_variances"` (`variance_effect`, `variance_within`) and `"convert"`.
#'   The two between-subjects methods additionally offer `"from_means"`
#'   (`means`, `sd`, optional group sizes `n`).
#' * `"convert"` takes `f`, `from`, `to`. ANOVA convention codes are
#'   `"gpower30"`, `"gpower30_implicit_rho"`, `"spss"`, `"cohen1988"`;
#'   MANOVA codes are `"on"` and `"off"`
#'   for using the mean correlation. Supply `measurements` and `rho` for
#'   `"gpower30"` or `"on"`, and `measurements`, `n_total`, `groups` for
#'   `"spss"`. This conversion preserves noncentrality and is an extension
#'   of the reference application's drawer modes.
#' * `"f_rm_manova_within_between"` offers `"from_criterion"` (`criterion`,
#'   `approximation`, `groups`, `measurements`, optional `n_total`) and
#'   `"from_means_cov"` (a group-by-measurement matrix `means`,
#'   `approximation`, and either covariance matrix `sigma` or `sd` with
#'   correlation matrix `corr` or common correlation `rho`; optional group
#'   sizes `n` and `n_total`). Approximation codes are `"pillai_v"`,
#'   `"wilks_u"`, `"hotelling_lawley_t1"` and `"hotelling_lawley_t2"`.
#'
#' Use the full drawer input names, particularly `groups` and `measurements`.
#' R can partially match names such as `m` or `i` to the formal arguments
#' `method` or `input_mode`; alternatively supply drawer values in `inputs`.
#' The reference application's mean-based t-test drawers display magnitudes,
#' whereas its matched-pairs "From differences" drawer retains the sign.
#'
#' The returned `value` is signed as defined by the method (for two
#' independent means d = (mean_1 - mean_2) / sigma); `magnitude` is its
#' absolute value.
#'
#' @param method Method id (see [list_methods()]).
#' @param input_mode Drawer input mode (see Details).
#' @param ... Named drawer inputs.
#' @param inputs Named list of drawer inputs (alternative to `...`; a name
#'   given in both places must carry the same value).
#' @return An object of class `ig_effect_size`: a list with `method_id`,
#'   `input_mode`, `status` (an [ig_status()]), `effect` (the parameter
#'   name), `value`, `magnitude`, `formula`, `source`, `inputs` and
#'   `warnings`, plus named lists `outputs` (derived quantities) and `details`
#'   (supporting calculations), empty when the drawer supplies none.
#'   A method without a drawer gives status `"unsupported"` and
#'   `value = NA`. Invalid inputs signal an `inferguidR_invalid_input`
#'   condition.
#' @examples
#' es <- effect_size("t_means_two_independent", "two_sds",
#'   mean_1 = 0, mean_2 = 1, sd_1 = 0.5, sd_2 = 0.5)
#' es$value
#' es$magnitude
#' @seealso [analysis_spec()]
#' @export
effect_size <- function(method, input_mode, ..., inputs = list()) {
  meta <- .ig_get_metadata(method)
  merged <- .ig_merge_named(list(inputs, list(...)), "input")
  if (length(merged$problems)) .ig_abort_invalid(merged$problems, "invalid effect-size input")
  drawer <- get0(paste0(".ig_drawer_", method), envir = topenv(environment(effect_size)),
                 inherits = FALSE)
  if (!is.function(drawer)) {
    return(.ig_new_effect_size(method, if (.ig_is_string(input_mode)) input_mode else NA_character_,
                               ig_status("unsupported", sprintf("no effect-size drawer for %s yet",
                                                                method)),
                               effect = NA_character_, value = NA_real_, inputs = merged$values))
  }
  out <- drawer(input_mode, merged$values)
  if (!is.numeric(out$value) || length(out$value) != 1L || !is.finite(out$value)) {
    return(.ig_new_effect_size(method, input_mode, ig_status(
      "numerical_failure", sprintf("the drawer formula %s gave no finite value", out$formula)),
      effect = out$effect, value = NA_real_, formula = out$formula, source = out$source,
      inputs = out$inputs))
  }
  .ig_new_effect_size(method, input_mode, ig_status("ok"), effect = out$effect, value = out$value,
                      formula = out$formula, source = out$source, inputs = out$inputs,
                      warnings = out$warnings, outputs = out$outputs %||% list(),
                      details = out$details %||% list())
}

.ig_new_effect_size <- function(method_id, input_mode, status, effect, value, formula = NA_character_,
                                source = NA_character_, inputs = list(), warnings = character(),
                                outputs = list(), details = list()) {
  structure(list(method_id = method_id, input_mode = input_mode, status = status, effect = effect,
                 value = as.double(value), magnitude = abs(as.double(value)),
                 formula = formula, source = source, inputs = inputs,
                 warnings = as.character(warnings), outputs = outputs, details = details),
            class = "ig_effect_size")
}

#' @export
format.ig_effect_size <- function(x, ...) {
  c(sprintf("<ig_effect_size> %s | %s", x$method_id, x$input_mode),
    sprintf("  status:  %s", format(x$status)),
    if (identical(x$status$code, "ok")) c(
      sprintf("  %s = %s (magnitude %s)", x$effect, format(x$value, digits = 8),
              format(x$magnitude, digits = 8)),
      sprintf("  formula: %s", x$formula),
      sprintf("  source:  %s", x$source)),
    if (length(x$warnings)) paste0("  warning: ", x$warnings))
}

#' @export
print.ig_effect_size <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}

# Sweeps --------------------------------------------------------------------------------

#' Sweep one input of a specification
#'
#' Varies one input of a specification over the given values and solves each
#' resulting specification with [calculate_power()]. Every row is a complete,
#' separately validated specification with its own result. Values are kept in
#' the order given; rows that are not computable (an inadmissible value, an
#' infeasible target, an unsupported branch) keep their status and reason and
#' have no numbers, so gaps and jumps are preserved rather than interpolated.
#'
#' @param spec A `power_spec` created by [analysis_spec()].
#' @param free Name of the input to vary. It must be an input of `spec`
#'   (not its solve target).
#' @param values Values of the free input.
#' @return A data frame with one row per value and columns `value`,
#'   `admissible` (the value gives a valid specification), `status`,
#'   `reason`, one `solved_<target>` column per solve target,
#'   `achieved_power`, `achieved_beta`, `alpha_attained`, `n_total`,
#'   `spec_hash`, `result_hash` and the list columns `spec` and `result`.
#'   Attributes `free`, `fixed_inputs`, `fixed_options`, `method_id`, `mode`,
#'   `target` and `base_spec_hash` describe the sweep.
#' @examples
#' spec <- analysis_spec("t_means_two_independent", "a_priori",
#'   inputs = list(alpha = 0.05, power = 0.8, d = 0.5), options = list(tails = "two"))
#' sw <- sweep_design(spec, "d", c(0.4, 0.5, 0.6))
#' sw[, c("value", "status", "solved_n_total", "achieved_power")]
#' @seealso [compare_designs()], [calculate_power()]
#' @export
sweep_design <- function(spec, free, values) {
  .ig_check_spec_integrity(spec)
  if (!.ig_is_string(free)) .ig_abort_invalid("`free` must be the name of one input")
  if (free %in% spec$target) {
    .ig_abort_invalid(sprintf("'%s' is the solve target of this specification and cannot be swept", free))
  }
  if (!(free %in% names(spec$inputs))) {
    .ig_abort_invalid(sprintf("'%s' is not an input of this specification; inputs: %s", free,
                              paste(names(spec$inputs), collapse = ", ")))
  }
  if (!is.numeric(values) || !length(values)) {
    .ig_abort_invalid("`values` must be a non-empty numeric vector")
  }
  rows <- lapply(as.double(values), function(v) {
    s <- tryCatch(stats::update(spec, inputs = stats::setNames(list(v), free)),
                  inferguidR_invalid_input = function(e) e)
    if (inherits(s, "condition")) {
      return(list(value = v, admissible = FALSE, status = "invalid_input",
                  reason = conditionMessage(s), spec = NULL, result = NULL))
    }
    r <- calculate_power(s)
    list(value = v, admissible = TRUE, status = r$status$code, reason = r$status$reason,
         spec = s, result = r)
  })
  get_num <- function(r, f) {
    if (is.null(r$result) || !identical(r$result$status$code, "ok")) return(NA_real_)
    v <- f(r$result)
    if (is.null(v) || length(v) != 1L) NA_real_ else as.double(v)
  }
  out <- data.frame(
    value = vapply(rows, `[[`, numeric(1), "value"),
    admissible = vapply(rows, `[[`, logical(1), "admissible"),
    status = vapply(rows, `[[`, character(1), "status"),
    reason = vapply(rows, `[[`, character(1), "reason"),
    stringsAsFactors = FALSE
  )
  for (tg in spec$target) {
    out[[paste0("solved_", tg)]] <- vapply(rows, get_num, numeric(1), f = function(x) x$solved[[tg]])
  }
  out$achieved_power <- vapply(rows, get_num, numeric(1), f = function(x) x$achieved$power)
  out$achieved_beta <- vapply(rows, get_num, numeric(1), f = function(x) x$achieved$beta)
  out$alpha_attained <- vapply(rows, get_num, numeric(1), f = function(x) x$achieved$alpha_attained)
  out$n_total <- vapply(rows, get_num, numeric(1), f = function(x) {
    x$solved$n_total %||% x$quantities$n_total
  })
  out$spec_hash <- vapply(rows, function(r) if (is.null(r$spec)) NA_character_ else r$spec$spec_hash,
                          character(1))
  out$result_hash <- vapply(rows, function(r) {
    if (is.null(r$result)) NA_character_ else r$result$result_hash
  }, character(1))
  out$spec <- I(lapply(rows, `[[`, "spec"))
  out$result <- I(lapply(rows, `[[`, "result"))
  attr(out, "free") <- free
  attr(out, "fixed_inputs") <- spec$inputs[setdiff(names(spec$inputs), free)]
  attr(out, "fixed_options") <- spec$options
  attr(out, "method_id") <- spec$method_id
  attr(out, "mode") <- spec$mode
  attr(out, "target") <- spec$target
  attr(out, "base_spec_hash") <- spec$spec_hash
  out
}

# Comparisons ------------------------------------------------------------------------------

#' Compare several designs
#'
#' Computes (where needed) and tabulates several analyses side by side.
#' Specifications are solved with [calculate_power()]; existing results are
#' used as they are. Results that are not `"ok"` keep their status and
#' reason and contribute no numbers.
#'
#' @param ... `power_spec` or `power_result` objects, optionally named.
#' @param designs A list of further `power_spec` or `power_result` objects.
#' @return A data frame with one row per design and columns `design`,
#'   `method_id`, `mode`, `target`, `status`, `reason`, `achieved_power`,
#'   `achieved_beta`, `alpha_nominal`, `alpha_attained`, `n_total`, `df`,
#'   `ncp`, `spec_hash`, `result_hash` and the list columns `solved`,
#'   `group_sizes`, `critical` and `result`.
#' @examples
#' base <- analysis_spec("t_means_two_independent", "post_hoc",
#'   inputs = list(alpha = 0.05, n1 = 30, n2 = 30, d = 0.5), options = list(tails = "two"))
#' compare_designs(equal = base, unequal = update(base, inputs = list(n1 = 20, n2 = 40)))[,
#'   c("design", "status", "achieved_power", "n_total")]
#' @seealso [sweep_design()]
#' @export
compare_designs <- function(..., designs = list()) {
  items <- c(list(...), designs)
  if (!length(items)) .ig_abort_invalid("no designs to compare")
  nms <- names(items)
  if (is.null(nms)) nms <- rep("", length(items))
  nms[is.na(nms) | !nzchar(nms)] <- paste0("design_", seq_along(items))[is.na(nms) | !nzchar(nms)]
  if (anyDuplicated(nms)) .ig_abort_invalid("design names must be unique")
  results <- lapply(seq_along(items), function(i) {
    x <- items[[i]]
    if (inherits(x, "power_result")) return(x)
    if (inherits(x, "power_spec")) return(calculate_power(x))
    .ig_abort_invalid(sprintf("design '%s' is neither a power_spec nor a power_result", nms[i]))
  })
  num <- function(f) vapply(results, function(r) {
    if (!identical(r$status$code, "ok")) return(NA_real_)
    v <- f(r)
    if (is.null(v) || length(v) != 1L) NA_real_ else as.double(v)
  }, numeric(1))
  ok_list <- function(f) I(lapply(results, function(r) if (identical(r$status$code, "ok")) f(r) else NULL))
  data.frame(
    design = nms,
    method_id = vapply(results, `[[`, character(1), "method_id"),
    mode = vapply(results, `[[`, character(1), "mode"),
    target = vapply(results, function(r) paste(r$target, collapse = "+"), character(1)),
    status = vapply(results, function(r) r$status$code, character(1)),
    reason = vapply(results, function(r) r$status$reason, character(1)),
    achieved_power = num(function(r) r$achieved$power),
    achieved_beta = num(function(r) r$achieved$beta),
    alpha_nominal = num(function(r) r$quantities$alpha_nominal),
    alpha_attained = num(function(r) r$achieved$alpha_attained),
    n_total = num(function(r) r$solved$n_total %||% r$quantities$n_total),
    df = num(function(r) r$quantities$df),
    ncp = num(function(r) r$quantities$ncp),
    spec_hash = vapply(results, `[[`, character(1), "spec_hash"),
    result_hash = vapply(results, `[[`, character(1), "result_hash"),
    solved = ok_list(function(r) unlist(r$solved)),
    group_sizes = ok_list(function(r) r$group_sizes),
    critical = ok_list(function(r) r$quantities$critical),
    result = I(results),
    stringsAsFactors = FALSE
  )
}
