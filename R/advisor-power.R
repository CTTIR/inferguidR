# Explicit mappings are keyed by the selected analysis variant, including overrides.
.ig_advisor_boschloo_variant <- function(variant_id) {
  if (!.ig_is_string(variant_id)) return(NULL)
  switch(variant_id %||% "", boschloo_two_sided_central = "central",
    boschloo_two_sided_minlike = "minlike", boschloo_one_sided = "direction_required", NULL)
}

.ig_advisor_boschloo_outcome <- function(outcome_id, variant_id) {
  identical(outcome_id, "out.two_group.boschloo_central") &&
    identical(variant_id, "boschloo_two_sided_central") ||
    identical(outcome_id, "out.two_group.boschloo_one_sided") &&
    identical(variant_id, "boschloo_one_sided")
}

.ig_advisor_boschloo_planning_note <- function(lang = "en") {
  if (lang == "de") paste(
    "Eine gepruefte Zuordnung zum eigenstaendigen Boschloo-Planer ist ueber advisor_power_spec() verfuegbar.",
    "Sie verwendet die ausgewaehlte zentrale, Minlike- oder vorab gerichtete Variante; p1 und p2 bezeichnen dieselben Gruppen wie im Analysebeispiel.",
    "Ereigniswahrscheinlichkeiten, Alpha, Zielpower und Allokation muessen ausdruecklich angegeben werden.",
    "Die numerischen Kontrollen und die Allokationsregel werden in der Planungsspezifikation festgehalten.",
    "Offene Angaben muessen zuerst geklaert werden. Keine G*Power-Paritaet oder vollstaendige numerische Zertifizierung wird behauptet.") else paste(
    "A checked mapping to the standalone Boschloo planner is available through advisor_power_spec().",
    "It uses the selected central, minlike or prespecified directional variant; p1 and p2 identify the same groups as the analysis example.",
    "Event probabilities, alpha, target power and allocation must be supplied explicitly.",
    "The numerical controls and allocation policy are retained in the planning specification.",
    "Resolve missing information first. No G*Power parity or complete numerical certification is claimed.")
}

#' Create a power specification for an explicitly mapped advisor variant
#'
#' Converts a completed advisor session to a power specification without
#' running a calculation. Currently the three Boschloo analysis variants are
#' mapped to the same standalone Boschloo engine: central, minlike and an
#' explicitly prespecified less/greater direction. Other methods and
#' unrelated overrides return a typed unsupported condition. Missing answers
#' must be resolved first. Existing saved sessions require the usual explicit
#' rules-version migration/replay; this function does not relabel them.
#'
#' Group 1 is `p1` and group 2 is `p2`, as in the exported analysis example.
#' `alternative = "greater"` means p1 > p2; `"less"` means p1 < p2. The
#' chosen direction is never inferred from the supplied probabilities.
#' The mapping establishes the identity of analysis and planning procedures,
#' not exhaustive numerical certification or agreement with G*Power.
#'
#' @param session A completed `advisor_session`.
#' @param inputs Named list with explicit `p1`, `p2`, `alpha`, and either
#'   `power`, `allocation_ratio` (a priori) or `n1`, `n2` (post hoc).
#'   Optional `n_max` bounds an a priori search. Numerical controls are
#'   supplied through `control`, not duplicated in `inputs`.
#' @param mode Either `"a_priori"` or `"post_hoc"`.
#' @param alternative Required `"less"` or `"greater"` for a one-sided
#'   recommendation; NULL or `"two.sided"` for a two-sided variant.
#' @param allocation Explicit a priori policy: `"n2_ceiling"` or
#'   `"total_nearest"`. Must be NULL for a fixed post-hoc design.
#' @param control Numerical settings from [boschloo_control()]. The generic
#'   planner supports its default `polish_tol`, `max_levels` and
#'   `ambiguity_band`; changed values are refused rather than discarded.
#'   Caching does not affect the numerical specification.
#' @param rules Optional matching `advisor_rules` object.
#' @return A `power_spec` for [calculate_power()] or [solve_design()].
#'   Retain the session/recommendation alongside the specification to preserve
#'   the question trace and scientific assumptions.
#' @family advisor
#' @export
advisor_power_spec <- function(session, inputs, mode = "a_priori", alternative = NULL,
                               allocation = NULL, control = boschloo_control(), rules = NULL) {
  rec <- advisor_recommendation(session, rules = rules)
  if (!identical(rec$planning$planning_availability, "available"))
    .ig_abort("unsupported", "This recommendation has no available same-test planning mapping; resolve missing information or use the stated separate procedure.")
  variant <- .ig_advisor_boschloo_variant(rec$preferred$analysis_variant_id)
  if (is.null(variant)) .ig_abort("unsupported", "No planning mapping exists for this selected analysis variant.")
  if (!.ig_is_string(mode) || !mode %in% c("a_priori", "post_hoc"))
    .ig_abort("invalid_input", "The Boschloo advisor bridge supports a_priori and post_hoc only.")
  if (variant == "direction_required") {
    if (!.ig_is_string(alternative) || !alternative %in% c("less", "greater"))
      .ig_abort("invalid_input", "Specify the prespecified alternative: less (p1 < p2) or greater (p1 > p2).")
    variant <- alternative
  } else if (!is.null(alternative) && !identical(alternative, "two.sided"))
    .ig_abort("invalid_input", "The selected two-sided variant cannot be changed to a directional test.")
  required <- c("p1", "p2", "alpha", if (mode == "a_priori") c("power", "allocation_ratio") else c("n1", "n2"))
  allowed <- c(required, if (mode == "a_priori") "n_max")
  if (!is.list(inputs) || is.null(names(inputs)) || anyDuplicated(names(inputs)) ||
      !all(required %in% names(inputs)) || any(!names(inputs) %in% allowed))
    .ig_abort("invalid_input", paste("Supply explicit inputs:", paste(required, collapse = ", "),
      "; only n_max is additionally supported in a_priori mode."))
  if (mode == "a_priori") {
    if (!.ig_is_string(allocation) || !allocation %in% c("n2_ceiling", "total_nearest"))
      .ig_abort("invalid_input", "Select an explicit allocation policy: n2_ceiling or total_nearest.")
  } else if (!is.null(allocation)) .ig_abort("invalid_input", "A fixed design uses explicit n1 and n2; allocation must be NULL.")
  .bl_check_control(control)
  control <- do.call(boschloo_control, unclass(control))
  defaults <- boschloo_control()
  for (id in c("polish_tol", "max_levels", "ambiguity_band"))
    if (!identical(control[[id]], defaults[[id]]))
      .ig_abort("unsupported", paste("The generic Boschloo planner does not support changing", id))
  inputs <- c(inputs, control[c("tie_tolerance", "grid_points", "nuisance_tol")])
  options <- list(variant = variant)
  if (mode == "a_priori") options$allocation <- allocation
  analysis_spec("boschloo_proportions_two_independent", mode, inputs, options)
}
