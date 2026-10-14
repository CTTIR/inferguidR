# Effect-size drawers: proportions (group G8) ------------------------------------------------
#
# Reached through effect_size() (R/core-api.R), which looks up
# `.ig_drawer_<method_id>()`. Method notes: specs/method-notes/<method_id>.md.
#
# * exact_proportion_one_sample_binomial: G*Power's drawer (manual 4.1,
#   Windows RT_DIALOG 193) computes g = P2 - P1 from the constant proportion P1
#   and P2 given directly or through a difference, ratio or odds ratio.
# * exact_proportions_fisher: G*Power's drawer "Calc P1 from ..." (RT_DIALOG
#   202) computes the proportion p1 from p2 and a difference, ratio or odds
#   ratio (certified drawer capture: radios "difference P1 - P2", "ratio P1/P2",
#   "odds ratio").
# * exact_proportions_mcnemar: not a G*Power drawer. Conversion of the
#   discordant cell probabilities into the odds ratio and pi_D, with the
#   corrected form of manual 5.1 (OR = (pi_D + d) / (pi_D - d) for
#   d = pi12 - pi21; the manual prints (d + pi_D) / (d - pi_D), which is
#   negative whenever |d| < pi_D; registry note).
# * z_proportions_two_independent: not a G*Power drawer. Cohen's h from two
#   proportions (the option input = "cohens_h" uses h instead of p2).

.ig_drawer_inputs <- function(input_mode, inputs, modes, required, ranges = list()) {
  if (!.ig_is_string(input_mode) || !(input_mode %in% modes)) {
    .ig_abort_invalid(sprintf("`input_mode` must be one of: %s", paste(modes, collapse = ", ")),
                      "invalid effect-size input")
  }
  need <- required[[input_mode]]
  problems <- character()
  unknown <- setdiff(names(inputs), need)
  if (length(unknown)) {
    problems <- c(problems, sprintf("input '%s' is not used in input mode '%s' (inputs: %s)", unknown, input_mode,
                                    paste(need, collapse = ", ")))
  }
  vals <- list()
  for (nm in need) {
    v <- inputs[[nm]]
    if (is.null(v)) {
      problems <- c(problems, sprintf("missing input '%s'", nm))
      next
    }
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v)) {
      problems <- c(problems, sprintf("'%s' must be a single finite number", nm))
      next
    }
    rg <- ranges[[nm]]
    if (!is.null(rg) && !(v > rg[1] && v < rg[2])) {
      problems <- c(problems, sprintf("'%s' must lie in (%s, %s)", nm, format(rg[1]), format(rg[2])))
      next
    }
    vals[[nm]] <- as.double(v)
  }
  if (length(problems)) .ig_abort_invalid(problems, "invalid effect-size input")
  vals
}

# Proportion from a reference proportion and a relation (difference, ratio,
# odds ratio); NA when the result leaves (0, 1).
.ig_prop_from_relation <- function(ref, mode, value) {
  out <- switch(mode,
                difference = ref + value,
                ratio = ref * value,
                odds_ratio = value * ref / (1 - ref + value * ref))
  if (!is.finite(out) || out <= 0 || out >= 1) NA_real_ else out
}

.ig_drawer_exact_proportion_one_sample_binomial <- function(input_mode, inputs) {
  modes <- c("direct", "difference", "ratio", "odds_ratio")
  unit <- c(0, 1)
  v <- .ig_drawer_inputs(input_mode, inputs, modes,
                         required = list(direct = c("p1", "p2"), difference = c("p1", "difference"),
                                         ratio = c("p1", "ratio"), odds_ratio = c("p1", "odds_ratio")),
                         ranges = list(p1 = unit, p2 = unit, ratio = c(0, Inf), odds_ratio = c(0, Inf)))
  p2 <- if (input_mode == "direct") v$p2 else .ig_prop_from_relation(v$p1, input_mode, v[[input_mode]])
  if (is.na(p2)) {
    .ig_abort_invalid(sprintf("the %s %s with P1 = %s gives a proportion P2 outside (0, 1)", gsub("_", " ", input_mode),
                              format(v[[input_mode]]), format(v$p1)), "invalid effect-size input")
  }
  g <- p2 - v$p1
  list(effect = "g", value = g,
       formula = switch(input_mode,
                        direct = "g = P2 - P1",
                        difference = "P2 = P1 + difference; g = P2 - P1",
                        ratio = "P2 = P1 * ratio; g = P2 - P1",
                        odds_ratio = "P2 = OR P1 / (1 - P1 + OR P1); g = P2 - P1"),
       source = "manual 4.1 (drawer of the one-sample binomial test; P1 = constant proportion)",
       inputs = c(v, list(p2 = p2, constant_proportion = v$p1)),
       warnings = if (abs(g) <= .ig_prop_eps) "|g| <= 1e-6: G*Power refuses such an effect in a priori analyses (manual 4.1)" else character())
}

.ig_drawer_exact_proportions_fisher <- function(input_mode, inputs) {
  modes <- c("direct", "difference", "ratio", "odds_ratio")
  unit <- c(0, 1)
  v <- .ig_drawer_inputs(input_mode, inputs, modes,
                         required = list(direct = c("p1", "p2"), difference = c("p2", "difference"),
                                         ratio = c("p2", "ratio"), odds_ratio = c("p2", "odds_ratio")),
                         ranges = list(p1 = unit, p2 = unit, ratio = c(0, Inf), odds_ratio = c(0, Inf)))
  p1 <- if (input_mode == "direct") v$p1 else .ig_prop_from_relation(v$p2, input_mode, v[[input_mode]])
  if (is.na(p1)) {
    .ig_abort_invalid(sprintf("the %s %s with P2 = %s gives a proportion P1 outside (0, 1)", gsub("_", " ", input_mode),
                              format(v[[input_mode]]), format(v$p2)), "invalid effect-size input")
  }
  list(effect = "p1", value = p1,
       formula = switch(input_mode,
                        direct = "P1 as given",
                        difference = "P1 = P2 + difference (difference P1 - P2)",
                        ratio = "P1 = P2 * ratio (ratio P1/P2)",
                        odds_ratio = "P1 = OR P2 / (1 - P2 + OR P2) (odds ratio (P1/(1-P1)) / (P2/(1-P2)))"),
       source = "G*Power drawer 'Calc P1 from ...' of Fisher's exact test (certified drawer capture)",
       inputs = c(v, list(p1 = p1, difference = p1 - v$p2, ratio = p1 / v$p2,
                          odds_ratio = (p1 / (1 - p1)) / (v$p2 / (1 - v$p2)))),
       warnings = character())
}

.ig_drawer_exact_proportions_mcnemar <- function(input_mode, inputs) {
  modes <- c("cells", "difference")
  v <- .ig_drawer_inputs(input_mode, inputs, modes,
                         required = list(cells = c("p12", "p21"), difference = c("pd", "difference")),
                         ranges = list(p12 = c(0, 1), p21 = c(0, 1), pd = c(0, 1), difference = c(-1, 1)))
  if (input_mode == "cells") {
    if (v$p12 + v$p21 >= 1) .ig_abort_invalid("p12 + p21 must be below 1", "invalid effect-size input")
    pd <- v$p12 + v$p21
    d <- v$p12 - v$p21
  } else {
    pd <- v$pd
    d <- v$difference
    if (abs(d) >= pd) .ig_abort_invalid("|difference| must be smaller than pi_D (both cells positive)", "invalid effect-size input")
  }
  or <- (pd + d) / (pd - d)
  list(effect = "or", value = or,
       formula = "OR = pi12 / pi21 = (pi_D + d) / (pi_D - d) with pi_D = pi12 + pi21, d = pi12 - pi21",
       source = paste("manual 5.1 with the corrected sign (the manual prints (d + pi_D)/(d - pi_D), negative for",
                      "|d| < pi_D); not a G*Power drawer"),
       inputs = c(v, list(pd = pd, difference = d, p12 = (pd + d) / 2, p21 = (pd - d) / 2)),
       warnings = character())
}

.ig_drawer_z_proportions_two_independent <- function(input_mode, inputs) {
  v <- .ig_drawer_inputs(input_mode, inputs, "proportions", required = list(proportions = c("p1", "p2")),
                         ranges = list(p1 = c(0, 1), p2 = c(0, 1)))
  h <- 2 * asin(sqrt(v$p2)) - 2 * asin(sqrt(v$p1))
  list(effect = "h", value = h,
       formula = "h = 2 asin(sqrt(p2)) - 2 asin(sqrt(p1))",
       source = "Cohen (1988), effect size h; sign as G*Power's option 'Use Cohen's effect size index h' (p2 relative to p1)",
       inputs = v,
       warnings = character())
}
