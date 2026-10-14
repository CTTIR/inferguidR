# Effect-size drawer: logistic regression ----------------------------------------------------
#
# G*Power's drawer for "Logistic regression" (manual 29.1; drawer_Z_LogRegress
# on the Mac, RT_DIALOG 236 on Windows) converts between the two event
# probabilities and the odds ratio:
#
#   OR = [p2 / (1 - p2)] / [p1 / (1 - p1)],    beta1 = log(OR)
#
# with p1 = Pr(Y=1|X=1) under H0 (the event probability at X = 0 of the H1
# model) and p2 = Pr(Y=1|X=1) under H1. It is reached through effect_size()
# (R/core-api.R), which looks up `.ig_drawer_<method_id>()`. Poisson
# regression has no drawer in G*Power (registry: drawer null).
#
# Input modes:
#   "probabilities": p1, p2 -> odds ratio (the drawer's transfer to the main
#                    window; captures GPW-UT-DRAW-UTL-U068-*)
#   "odds_ratio":    p1, odds_ratio -> p2

.ig_drawer_z_logistic_regression <- function(input_mode, inputs) {
  modes <- c("probabilities", "odds_ratio")
  if (!.ig_is_string(input_mode) || !(input_mode %in% modes)) {
    .ig_abort_invalid(sprintf("`input_mode` must be one of: %s", paste(modes, collapse = ", ")),
                      "invalid effect-size input")
  }
  required <- switch(input_mode, probabilities = c("p1", "p2"), odds_ratio = c("p1", "odds_ratio"))
  problems <- character()
  unknown <- setdiff(names(inputs), required)
  if (length(unknown)) {
    problems <- c(problems, sprintf("input '%s' is not used in input mode '%s' (inputs: %s)", unknown,
                                    input_mode, paste(required, collapse = ", ")))
  }
  num <- function(nm) {
    v <- inputs[[nm]]
    if (is.null(v)) {
      problems <<- c(problems, sprintf("missing input '%s'", nm))
      return(NA_real_)
    }
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v)) {
      problems <<- c(problems, sprintf("'%s' must be a single finite number", nm))
      return(NA_real_)
    }
    as.double(v)
  }
  vals <- lapply(stats::setNames(required, required), num)
  for (nm in intersect(c("p1", "p2"), required)) {
    v <- vals[[nm]]
    if (!is.na(v) && !(v > 0 && v < 1)) problems <- c(problems, sprintf("'%s' must lie in (0, 1)", nm))
  }
  if (!is.null(vals$odds_ratio) && !is.na(vals$odds_ratio) && !(vals$odds_ratio > 0)) {
    problems <- c(problems, "'odds_ratio' must be positive")
  }
  if (length(problems)) .ig_abort_invalid(problems, "invalid effect-size input")

  b0 <- stats::qlogis(vals$p1)
  if (input_mode == "probabilities") {
    b1 <- stats::qlogis(vals$p2) - b0
    list(effect = "odds_ratio", value = exp(b1),
         formula = "OR = [p2 / (1 - p2)] / [p1 / (1 - p1)]",
         source = "manual 29.1; G*Power drawer 'probabilities <-> odds ratio'",
         inputs = vals, warnings = character())
  } else {
    list(effect = "p2", value = stats::plogis(b0 + log(vals$odds_ratio)),
         formula = "p2 = OR p1 / (1 - p1 + OR p1)",
         source = "manual 29.1 (inverse of the odds-ratio definition)",
         inputs = vals, warnings = character())
  }
}
