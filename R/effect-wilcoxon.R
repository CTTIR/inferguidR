# Effect-size drawers: Wilcoxon tests (group G10) ---------------------------------------------
#
# G*Power's "Determine =>" drawers of the three Wilcoxon procedures are the
# drawers of the corresponding t tests (Mac nibs drawer_T_OneMeanWilcox,
# drawer_T_PairedMeansWilcox, drawer_T_2GroupsWilcox; registry
# effect_size$drawer). They are reached through effect_size()
# (R/core-api.R), which looks up `.ig_drawer_<method_id>()`.
#
# Sign conventions, from the emulator captures of G*Power 3.1.9.7
# (util_effect_size_drawers records GPW-UT-DRAW-UTL-U049..U053) and the
# manual examples DOC-127/DOC-128:
#
# * one sample: d = (mean_h1 - mean_h0) / sd; every capture has
#   mean_h1 > mean_h0, so the sign is the one of the t-test drawer
#   definition (not tested by a capture);
# * matched pairs, "From differences": dz = mean_diff / sd_diff, displayed
#   with its sign (U050-UE01: mean -5, SD 0.1 gives -50);
# * matched pairs, "From group parameters": dz = |mean_1 - mean_2| / sd_z
#   with sd_z = sqrt(sd_1^2 + sd_2^2 - 2 r sd_1 sd_2); G*Power displays the
#   magnitude for both orders of the means (U051-UE01, DOC-127);
# * two groups: d = |mean_1 - mean_2| / sigma, sigma the common SD ("n1 !=
#   n2") or sqrt((sd_1^2 + sd_2^2) / 2) ("n1 = n2"); G*Power displays the
#   magnitude for both orders of the means (U052-UO01, U052-UE01).
#
# `value` is signed as defined here (mean_h1 - mean_h0, mean_1 - mean_2 or
# the mean difference); `magnitude` is |value|. The Lehmann method uses |d|,
# the A.R.E. method uses the sign of d only for the side of a one-sided test.

.ig_drawer_t_wilcoxon_signed_rank_one_sample <- function(input_mode, inputs) {
  .ig_wilcoxon_drawer(input_mode, inputs, modes = list(
    means_sd = list(required = c("mean_h0", "mean_h1", "sd"), positive = "sd")),
    compute = function(mode, v) list(
      value = (v$mean_h1 - v$mean_h0) / v$sd,
      formula = "d = (mean_h1 - mean_h0) / sd",
      source = paste("manual ch. 22 (effect size as in the one-sample t test); G*Power drawer 'from means",
                     "and SD' (inputs Mean H0, Mean H1, SD sigma)"),
      warnings = character()),
    effect = "d")
}

.ig_drawer_t_wilcoxon_signed_rank_matched_pairs <- function(input_mode, inputs) {
  .ig_wilcoxon_drawer(input_mode, inputs, modes = list(
    differences = list(required = c("mean_diff", "sd_diff"), positive = "sd_diff"),
    group_parameters = list(required = c("mean_1", "mean_2", "sd_1", "sd_2", "r"),
                            positive = c("sd_1", "sd_2"))),
    compute = function(mode, v) {
      if (identical(mode, "differences")) {
        return(list(value = v$mean_diff / v$sd_diff, formula = "dz = mean_diff / sd_diff",
                    source = "manual ch. 23; G*Power drawer 'From differences' (signed, as displayed)",
                    warnings = character()))
      }
      if (abs(v$r) > 1) .ig_abort_invalid("'r' must lie in [-1, 1]", "invalid effect-size input")
      var_z <- v$sd_1^2 + v$sd_2^2 - 2 * v$r * v$sd_1 * v$sd_2
      if (!(var_z > 0)) {
        .ig_abort_invalid("the SD of the differences sqrt(sd_1^2 + sd_2^2 - 2 r sd_1 sd_2) is zero",
                          "invalid effect-size input")
      }
      list(value = (v$mean_1 - v$mean_2) / sqrt(var_z),
           formula = "dz = (mean_1 - mean_2) / sqrt(sd_1^2 + sd_2^2 - 2 r sd_1 sd_2)",
           source = paste("manual ch. 23 (Fig. 31, DOC-127); G*Power drawer 'From group parameters', which",
                          "displays the magnitude |dz|"),
           warnings = character())
    },
    effect = "dz")
}

.ig_drawer_t_wilcoxon_mann_whitney <- function(input_mode, inputs) {
  .ig_wilcoxon_drawer(input_mode, inputs, modes = list(
    common_sd = list(required = c("mean_1", "mean_2", "sd"), positive = "sd"),
    two_sds = list(required = c("mean_1", "mean_2", "sd_1", "sd_2"), positive = c("sd_1", "sd_2"))),
    compute = function(mode, v) {
      if (identical(mode, "common_sd")) {
        return(list(value = (v$mean_1 - v$mean_2) / v$sd, formula = "d = (mean_1 - mean_2) / sd",
                    source = paste("manual ch. 24 (example 24.3, DOC-128); G*Power drawer 'n1 != n2' (common",
                                   "SD within each group), which displays the magnitude |d|"),
                    warnings = character()))
      }
      list(value = (v$mean_1 - v$mean_2) / sqrt((v$sd_1^2 + v$sd_2^2) / 2),
           formula = "d = (mean_1 - mean_2) / sqrt((sd_1^2 + sd_2^2) / 2)",
           source = "G*Power drawer 'n1 = n2' (as for two independent means, manual 21.1), which displays |d|",
           warnings = paste("the mean SD sqrt((sd_1^2 + sd_2^2) / 2) is meant for equal group sizes; with",
                            "clearly unequal groups it can give misleading power (Cohen 1969, p. 42)"))
    },
    effect = "d")
}

# Shared input checks of the drawers.
.ig_wilcoxon_drawer <- function(input_mode, inputs, modes, compute, effect) {
  if (!.ig_is_string(input_mode) || !(input_mode %in% names(modes))) {
    .ig_abort_invalid(sprintf("`input_mode` must be one of: %s", paste(names(modes), collapse = ", ")),
                      "invalid effect-size input")
  }
  spec <- modes[[input_mode]]
  problems <- character()
  unknown <- setdiff(names(inputs), spec$required)
  if (length(unknown)) {
    problems <- c(problems, sprintf("input '%s' is not used in input mode '%s' (inputs: %s)", unknown,
                                    input_mode, paste(spec$required, collapse = ", ")))
  }
  vals <- list()
  for (nm in spec$required) {
    v <- inputs[[nm]]
    if (is.null(v)) {
      problems <- c(problems, sprintf("missing input '%s'", nm))
      vals[[nm]] <- NA_real_
    } else if (!is.numeric(v) || length(v) != 1L || !is.finite(v)) {
      problems <- c(problems, sprintf("'%s' must be a single finite number", nm))
      vals[[nm]] <- NA_real_
    } else {
      vals[[nm]] <- as.double(v)
    }
  }
  for (nm in spec$positive) {
    if (!is.na(vals[[nm]]) && vals[[nm]] <= 0) problems <- c(problems, sprintf("'%s' must be positive", nm))
  }
  if (length(problems)) .ig_abort_invalid(problems, "invalid effect-size input")
  out <- compute(input_mode, vals)
  list(effect = effect, value = out$value, magnitude = abs(out$value), formula = out$formula,
       source = out$source, inputs = vals, warnings = out$warnings)
}
