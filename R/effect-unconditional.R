# Effect-size conversions: unconditional exact tests of two proportions (group G9) ------------
#
# G*Power 3.1.9.7 offers no effect-size drawer for "Proportions: Inequality, two
# independent groups (unconditional)" or its offset variant (the certified
# inventory lists no Determine button for either test). Its Options dialog
# instead selects the "Effect size index" in which the effect is entered and,
# in sensitivity mode, displayed: Proportion p1, Difference (p1-p2), Risk ratio
# (p1/p2) or Odds ratio, always together with Proportion p2. The drawers below
# perform these conversions so that an effect known on one scale can be
# entered on another; they are reached through effect_size().
#
# Input modes (both methods; for the offset test apply them to the |H1 and the
# |H0 value separately, with the same p2):
#   difference, risk_ratio, odds_ratio        index of the pair (p1, p2)
#   p1_from_difference, p1_from_risk_ratio,   p1 from an index value and p2
#   p1_from_odds_ratio

.ig_drawer_exact_proportions_unconditional <- function(input_mode, inputs) {
  .ig_unc_drawer(input_mode, inputs)
}

.ig_drawer_exact_proportions_unconditional_offset <- function(input_mode, inputs) {
  .ig_unc_drawer(input_mode, inputs)
}

.ig_unc_drawer <- function(input_mode, inputs) {
  modes <- c("difference", "risk_ratio", "odds_ratio", "p1_from_difference", "p1_from_risk_ratio",
             "p1_from_odds_ratio")
  if (!.ig_is_string(input_mode) || !(input_mode %in% modes)) {
    .ig_abort_invalid(sprintf("`input_mode` must be one of: %s", paste(modes, collapse = ", ")),
                      "invalid effect-size input")
  }
  required <- switch(input_mode,
                     difference = , risk_ratio = , odds_ratio = c("p1", "p2"),
                     p1_from_difference = c("difference", "p2"),
                     p1_from_risk_ratio = c("risk_ratio", "p2"),
                     p1_from_odds_ratio = c("odds_ratio", "p2"))
  problems <- character()
  unknown <- setdiff(names(inputs), required)
  if (length(unknown)) {
    problems <- c(problems, sprintf("input '%s' is not used in input mode '%s' (inputs: %s)", unknown,
                                    input_mode, paste(required, collapse = ", ")))
  }
  vals <- list()
  for (nm in required) {
    v <- inputs[[nm]]
    if (is.null(v)) {
      problems <- c(problems, sprintf("missing input '%s'", nm))
    } else if (!is.numeric(v) || length(v) != 1L || !is.finite(v)) {
      problems <- c(problems, sprintf("'%s' must be a single finite number", nm))
    } else {
      vals[[nm]] <- as.double(v)
    }
  }
  inside <- function(nm) {
    if (!is.null(vals[[nm]]) && (vals[[nm]] <= 0 || vals[[nm]] >= 1)) {
      problems <<- c(problems, sprintf("'%s' must lie in (0, 1)", nm))
    }
  }
  for (nm in intersect(c("p1", "p2"), required)) inside(nm)
  for (nm in intersect(c("risk_ratio", "odds_ratio"), required)) {
    if (!is.null(vals[[nm]]) && vals[[nm]] <= 0) problems <- c(problems, sprintf("'%s' must be positive", nm))
  }
  if (length(problems)) .ig_abort_invalid(problems, "invalid effect-size input")
  src <- "G*Power's 'Effect size index' conversions of the unconditional tests (Options dialog; no drawer exists)"
  if (input_mode %in% c("difference", "risk_ratio", "odds_ratio")) {
    idx <- .ig_unc_indices(vals$p1, vals$p2)
    value <- idx[[input_mode]]
    formula <- switch(input_mode,
                      difference = "difference = p1 - p2",
                      risk_ratio = "risk ratio = p1 / p2",
                      odds_ratio = "odds ratio = (p1 / (1 - p1)) / (p2 / (1 - p2))")
    return(list(effect = input_mode, value = value, formula = formula, source = src, inputs = vals,
                warnings = character()))
  }
  index <- sub("^p1_from_", "", input_mode)
  p1 <- .ig_unc_p1_from_index(index, vals[[index]], vals$p2)
  if (is.na(p1) || p1 <= 0 || p1 >= 1) {
    .ig_abort_invalid(sprintf("%s = %s with p2 = %s gives p1 = %s outside (0, 1)", index, format(vals[[index]]),
                              format(vals$p2), format(p1)), "invalid effect-size input")
  }
  formula <- switch(index,
                    difference = "p1 = p2 + difference",
                    risk_ratio = "p1 = risk ratio * p2",
                    odds_ratio = "p1 = odds ratio * p2 / (1 - p2 + odds ratio * p2)")
  list(effect = "p1", value = p1, formula = formula, source = src, inputs = vals, warnings = character())
}
