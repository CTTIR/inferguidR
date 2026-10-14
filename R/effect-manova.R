# Effect-size drawers: MANOVA global effects and special effects ----------------------
#
# Ownership: group G4 (with R/methods-manova.R). G*Power's drawers for "MANOVA:
# Global effects" (drawer_F_Global, RT_DIALOG 192/201) and "MANOVA: Special
# effects and interactions" (drawer_F_Special, RT_DIALOG 204) convert a value
# of the multivariate criterion of the selected approximation into the effect
# size f2 of that approximation (G*Power's f²(V), f²(U), f²(T1), f²(T2)).
# Reached through effect_size() (R/core-api.R) as .ig_drawer_<method_id>().
#
# Inputs (live drawer fields, certified inventory and supplement captures
# SUP-G4-MANG-DRW-*, SUP-G4-MANS-DRW-*): the criterion value, the number of
# groups k (global effects only) and the response variables p; special effects
# take the number of predictors q instead of k. With q = k - 1 (global) or the
# number of predictors (special) and s = min(p, q):
#
#   Pillai V             f2 = V / (s - V)            (0 < V < s)
#   Wilks U              f2 = U^(-1/t) - 1           (0 < U < 1; Rao's t, see R/methods-manova.R)
#   Hotelling-Lawley T1  f2 = T / s                  (T > 0)
#   Hotelling-Lawley T2  f2 = T / c, c = (df2 - 2) / (nu_e - p - 1) with McKeon's df2 at
#                        error df nu_e = N - k       (T > 0)
#
# These are the inverses of the criterion values reported by the forward
# model (R/methods-manova.R). McKeon's scale c depends on the error df, so the
# T2 conversion needs the total sample size N and the number of groups k; the
# G*Power drawer takes both from the main window (in a priori mode it uses the
# value of the hidden 'Total sample size' field and says so in a message).

.ig_drawer_f_manova_global <- function(input_mode, inputs) .ig_drawer_f_manova(input_mode, inputs, "global")

.ig_drawer_f_manova_special <- function(input_mode, inputs) .ig_drawer_f_manova(input_mode, inputs, "special")

.ig_drawer_f_manova <- function(input_mode, inputs, effect) {
  global <- identical(effect, "global")
  modes <- names(.ig_manova_criterion_labels)
  if (!.ig_is_string(input_mode) || !(input_mode %in% modes)) {
    .ig_abort_invalid(sprintf("`input_mode` must be one of: %s", paste(modes, collapse = ", ")),
                      "invalid effect-size input")
  }
  design <- if (global) c("n_groups", "n_responses") else c("n_responses", "n_predictors")
  if (identical(input_mode, "hotelling_lawley_t2")) design <- union(c("n_total", "n_groups"), design)
  required <- c(input_mode, design)
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
  minimum <- c(n_groups = 2, n_responses = 1, n_predictors = 1, n_total = 3)
  for (nm in intersect(names(minimum), required)) {
    v <- vals[[nm]]
    if (!is.na(v) && (v != round(v) || v < minimum[[nm]])) {
      problems <- c(problems, sprintf("'%s' must be a whole number >= %s", nm, minimum[[nm]]))
    }
  }
  if (length(problems)) .ig_abort_invalid(problems, "invalid effect-size input")

  p <- vals$n_responses
  q <- if (global) vals$n_groups - 1 else vals$n_predictors
  s <- min(p, q)
  x <- vals[[input_mode]]
  label <- .ig_manova_criterion_labels[[input_mode]]
  if (identical(input_mode, "pillai_v")) {
    if (!(x > 0 && x < s)) problems <- sprintf("Pillai V must lie in (0, s) = (0, %s)", .ig_format_double(s))
    value <- x / (s - x)
    formula <- "f2 = V / (s - V), s = min(p, q)"
  } else if (identical(input_mode, "wilks_u")) {
    if (!(x > 0 && x < 1)) problems <- "Wilks U must lie in (0, 1)"
    value <- x^(-1 / .ig_manova_rao_t(p, q)) - 1
    formula <- "f2 = U^(-1/t) - 1, t = sqrt((p^2 q^2 - 4) / (p^2 + q^2 - 5)) (t = 1 when p^2 + q^2 - 5 <= 0)"
  } else if (identical(input_mode, "hotelling_lawley_t1")) {
    if (!(x > 0)) problems <- "Hotelling-Lawley T1 must be positive"
    value <- x / s
    formula <- "f2 = T / s, s = min(p, q)"
  } else {
    if (!(x > 0)) problems <- c(problems, "Hotelling-Lawley T2 must be positive")
    if (!(vals$n_total > vals$n_groups)) {
      problems <- c(problems, "the total sample size must be greater than the number of groups")
    }
    tm <- .ig_manova_terms(vals$n_total, vals$n_groups, p, q, "hotelling_lawley_t2")
    if (!is.finite(tm$g) || tm$g <= 0) {
      problems <- c(problems, sprintf(paste("McKeon's scale c = (df2 - 2) / (nu_e - p - 1) is %s at error df %s;",
                                            "the T2 conversion needs a larger total sample size"),
                                      .ig_format_double(tm$g), .ig_format_double(tm$error_df)))
    }
    value <- x / tm$g
    formula <- paste("f2 = T / c, c = (df2 - 2) / (nu_e - p - 1), df2 = 4 + (p q + 2) / (b - 1),",
                     "b = (nu_e + q - p - 1)(nu_e - 1) / ((nu_e - p - 3)(nu_e - p)), nu_e = N - k")
  }
  if (length(problems)) .ig_abort_invalid(problems, "invalid effect-size input")
  list(effect = "f2", value = value, formula = formula,
       source = paste("inverse of the", label, "criterion of the", if (global) "global" else "special",
                      "effects procedure (R/methods-manova.R; G*Power drawer 'Effect size from criterion')"),
       inputs = vals,
       warnings = if (identical(input_mode, "hotelling_lawley_t2")) {
         "the Hotelling-Lawley T2 conversion depends on the total sample size N and the number of groups k"
       } else character())
}
