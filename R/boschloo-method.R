# Boschloo's test in the generic analysis framework ---------------------------------
#
# Registers Boschloo's test as the method `boschloo_proportions_two_independent`
# (family "boschloo"; metadata in inst/registry/boschloo.json) so that a
# planning question can be written as an analysis_spec() and answered by the
# generic mode solvers of R/num-modes.R:
#
# * post_hoc: exact power of a fixed design (n1, n2);
# * a_priori: smallest design of the declared allocation sequence whose exact
#   power reaches the target. The method brings its own a priori search
#   (ig_method(a_priori_search =), .bl_method_a_priori_search()): the search
#   of boschloo_sample_size(), which skips the designs that the Neyman-Pearson
#   bound excludes for every level-alpha test, evaluates every later design
#   of the allocation sequence up to the first one reaching the target
#   (minimality by enumeration), and reports the stable design as well.
#
# The forward model evaluates the same rejection region as boschloo_region(),
# boschloo_power(), boschloo_sample_size() and the decisions of
# boschloo_test(), so power_result objects from calculate_power() and
# solve_design() agree with the dedicated functions to the last bit.
#
# Boschloo's test is not a G*Power procedure. Its registry entry therefore has
# reference_status "not_applicable", and the planning compatibility is
# "independently_validated_extension" (specs/method-notes/boschloo.md).
#
# Sensitivity, criterion and compromise modes are not offered. Power is a
# polynomial in p2 with one branch on each side of p1 for the two-sided
# variants (no unique effect root), and alpha acts only through discrete
# rejection regions; both would need their own discrete solvers.

.bl_method_id <- "boschloo_proportions_two_independent"

# Numerical settings that the specification does not expose. They are fixed
# for this method and recorded in every result (method_outputs$settings).
.bl_method_fixed_control <- list(polish_tol = 1e-12, max_levels = 60L,
                                 ambiguity_band = c(1e-12, 1e-6))

.ig_define_methods_boschloo <- function() {
  list(.bl_method_definition())
}

.bl_method_definition <- function() {
  info <- .bl_variant_info
  ig_method(
    id = .bl_method_id,
    family = "boschloo",
    label = list(en = "Proportions: Inequality, two independent groups (Boschloo's unconditional exact test)",
                 de = "Anteile: Ungleichheit, zwei unabh\u00e4ngige Gruppen (unbedingter exakter Test nach Boschloo)"),
    params = ig_params(
      ig_param("p1", role = "effect", domain = ig_domain_real(0, 1),
               label = list(en = "Proportion p1 (group 1)", de = "Anteil p1 (Gruppe 1)"),
               description = paste("Event probability of group 1 (the reference group) under the",
                                   "alternative. Baseline-plus-effect input is converted with",
                                   "boschloo_probabilities().")),
      ig_param("p2", role = "effect", domain = ig_domain_real(0, 1),
               label = list(en = "Proportion p2 (group 2)", de = "Anteil p2 (Gruppe 2)"),
               description = "Event probability of group 2 under the alternative."),
      ig_param_alpha(),
      ig_param_power(),
      ig_param("n_total", role = "design", domain = ig_domain_integer(lower = 2, upper = 3000),
               unit = "subjects",
               label = list(en = "Total sample size", de = "Gesamtstichprobenumfang"),
               description = "Solved in a priori mode (n1 + n2).",
               modes = "a_priori"),
      ig_param("n1", role = "design", domain = ig_domain_integer(lower = 1, upper = 2999),
               unit = "subjects",
               label = list(en = "Sample size group 1", de = "Stichprobenumfang Gruppe 1"),
               aliases = "sample_size_group_1", modes = "post_hoc"),
      ig_param("n2", role = "design", domain = ig_domain_integer(lower = 1, upper = 2999),
               unit = "subjects",
               label = list(en = "Sample size group 2", de = "Stichprobenumfang Gruppe 2"),
               aliases = "sample_size_group_2", modes = "post_hoc"),
      ig_param("allocation_ratio", role = "design",
               domain = ig_domain_real(0, 100, lower_inclusive = FALSE), default = 1,
               label = list(en = "Allocation ratio N2/N1", de = "Allokationsverh\u00e4ltnis N2/N1"),
               description = paste("r = n2 / n1 of the a priori allocation sequence (see the option",
                                   "`allocation`). The default 1 is an inferguidR default."),
               aliases = "allocation_ratio_n2_n1", modes = "a_priori"),
      ig_param("n_max", role = "control", domain = ig_domain_integer(lower = 2, upper = 3000),
               default = 500L, unit = "subjects",
               label = list(en = "Largest total sample size searched",
                            de = "Gr\u00f6\u00dfter durchsuchter Gesamtstichprobenumfang"),
               description = paste("Bounds the run time of the a priori search; designs with",
                                   "n1 + n2 above n_max are not evaluated."),
               modes = "a_priori"),
      ig_param("tie_tolerance", role = "control",
               domain = ig_domain_real(.bl_method_fixed_control$ambiguity_band[1],
                                       .bl_method_fixed_control$ambiguity_band[2]),
               default = 1e-7,
               label = list(en = "Relative tie tolerance of the ordering statistic",
                            de = "Relative Bindungstoleranz der Ordnungsstatistik"),
               description = paste("See boschloo_control(); it must lie within the fixed ambiguity band",
                                   "[1e-12, 1e-6] of this method.")),
      ig_param("grid_points", role = "control", domain = ig_domain_integer(lower = 3, upper = 1e5),
               default = 1001L,
               label = list(en = "Nuisance grid points", de = "Gitterpunkte des St\u00f6rparameters"),
               description = "Global grid of the nuisance supremum before branch-and-bound refinement; see boschloo_control()."),
      ig_param("nuisance_tol", role = "control",
               domain = ig_domain_real(0, 1e-4, lower_inclusive = FALSE), default = 1e-9,
               label = list(en = "Tolerance of the rigorous supremum bound",
                            de = "Toleranz der strengen Supremumsschranke"),
               description = "See boschloo_control().")
    ),
    options = ig_options(
      ig_option("variant", values = names(info),
                labels = vapply(info, `[[`, "", "label"),
                label = list(en = "Boschloo variant", de = "Boschloo-Variante")),
      ig_option("allocation", values = .bl_allocations, default = "n2_ceiling",
                labels = c(n2_ceiling = "n1 = 1, 2, ...; n2 = max(1, ceiling(r * n1))",
                           total_nearest = "N = 2, 3, ...; n1 = round(N / (1 + r)), n2 = N - n1"),
                applies_when = list(mode = "a_priori"),
                label = list(en = "Allocation policy", de = "Allokationsregel"))
    ),
    modes = c("a_priori", "post_hoc"),
    power = .bl_method_power,
    design = .bl_method_design,
    discrete = TRUE,
    effect_direction = "unknown",
    engine = list(package = "inferguidR",
                  fun = paste("Boschloo engine (R/boschloo-*.R): Fisher ordering from stats::dhyper/phyper in",
                              "log space with exact integer fallback; Bernstein-form nuisance supremum with a",
                              "rigorous branch-and-bound bound; stats::dbinom double sum over the region"),
                  version = NA_character_),
    checker = "validation/independent/boschloo",
    reference = list(
      gpower_relation = .bl_evidence$gpower_relation,
      gpower_comparators = c("exact_proportions_fisher", "exact_proportions_unconditional"),
      gpower_classification_status = .bl_evidence$gpower_classification_status,
      planning_compatibility = .bl_evidence$planning_compatibility,
      allocation_rule = paste("declared allocation sequence: n2_ceiling (n1 = 1, 2, ...; n2 = max(1, ceiling(r n1)))",
                              "or total_nearest (N = 2, 3, ...; n1 = round(N / (1 + r)), n2 = N - n1)"),
      method_notes = "specs/method-notes/boschloo.md"
    ),
    status = .bl_evidence$validation_status,
    constraints = .bl_method_constraints,
    a_priori_search = .bl_method_a_priori_search
  )
}

# Control object of a fully specified input list.
.bl_method_control <- function(p) {
  fx <- .bl_method_fixed_control
  boschloo_control(tie_tolerance = p$tie_tolerance, grid_points = p$grid_points,
                   nuisance_tol = p$nuisance_tol, polish_tol = fx$polish_tol,
                   max_levels = fx$max_levels, ambiguity_band = fx$ambiguity_band)
}

.bl_method_constraints <- function(p) {
  if (is.null(p$n1) || is.null(p$n2)) return(character())
  if ((as.double(p$n1) + 1) * (as.double(p$n2) + 1) > .bl_max_tables) {
    return(sprintf("the design n1 = %s, n2 = %s has more than %.0f tables, the enumeration limit",
                   format(p$n1), format(p$n2), .bl_max_tables))
  }
  character()
}

# Forward model: exact power and attained size of the design (p$n1, p$n2).
.bl_method_power <- function(p) {
  info <- .bl_variant(p$variant)
  control <- .bl_method_control(p)
  n1 <- as.integer(p$n1)
  n2 <- as.integer(p$n2)
  .bl_check_design(n1, n2)
  reg <- .bl_region(n1, n2, as.double(p$alpha), info, control)
  if (length(reg$problems)) {
    .ig_abort("numerical_failure", sprintf("Boschloo region for n1 = %d, n2 = %d not certified: %s",
                                           n1, n2, paste(reg$problems, collapse = "; ")))
  }
  P <- outer(stats::dbinom(0:n1, n1, p$p1), stats::dbinom(0:n2, n2, p$p2))
  pw <- .bl_power_idx(reg$idx, n1, n2, p$p1, p$p2)
  # beta is summed over the accepted tables directly, so that it keeps its
  # relative accuracy when the power is close to 1.
  beta <- if (length(reg$idx)) min(1, sum(P[-reg$idx])) else 1
  region <- .bl_region_object(reg, control)
  outputs <- list(
    variant = info$variant,
    variant_id = info$variant_id,
    variant_label = info$label,
    alternative = info$alternative,
    attained_size_upper = reg$size$upper,
    tables_rejected = length(reg$idx),
    tables_total = length(reg$member),
    region_hash = region$region_hash,
    planning_compatibility = .bl_evidence$planning_compatibility,
    gpower_relation = sprintf("%s (%s)", .bl_evidence$gpower_relation,
                              .bl_evidence$gpower_classification_status),
    settings = .bl_control_record(control)
  )
  if (is.finite(reg$size$argmax)) outputs$attained_size_argmax <- reg$size$argmax
  if (length(region$warnings)) outputs$region_warnings <- region$warnings
  list(
    power = pw,
    beta = beta,
    alpha_attained = reg$size$value,
    n_total = n1 + n2,
    group_sizes = c(n1 = n1, n2 = n2),
    outputs = outputs,
    assumptions = c(
      "two independent binomial samples X1 ~ Bin(n1, p1) (group 1) and X2 ~ Bin(n2, p2) (group 2) with fixed group sizes; H0: p1 = p2",
      sprintf("%s; tables are ordered by %s", info$label, info$ordering_statistic),
      paste("rejection region = tables whose rigorous upper bound of the Boschloo p-value is at most alpha;",
            "power = double binomial sum over the region; attained size = supremum over the common event",
            "probability of the null rejection probability (alpha_attained)"),
      "no mid-p, no randomisation, no Berger-Boos restriction, no null margin"
    )
  )
}

# Allocation rule of the a priori search (see .bl_design_at()). Total sample
# sizes outside the allocation sequence are not admissible (NULL). Requests
# that cannot be answered are reported with a typed status when the search
# reaches its first design, and the search stops with status "infeasible"
# beyond n_max.
.bl_method_design <- function(p, n_total) {
  N <- as.integer(round(n_total))
  r <- as.double(p$allocation_ratio)
  allocation <- p$allocation
  first <- .bl_design_at(1L, r, allocation)
  if (N == first$n1 + first$n2) .bl_method_a_priori_checks(p)
  if (N > p$n_max) {
    .ig_abort("infeasible", sprintf(
      "no design of the allocation sequence with n1 + n2 <= n_max = %d reaches power %s (every design up to n_max was evaluated)",
      as.integer(p$n_max), format(p$power)))
  }
  m <- .bl_index_of_total(N, r, allocation)
  if (is.na(m)) return(NULL)
  d <- .bl_design_at(m, r, allocation)
  c(n1 = d$n1, n2 = d$n2)
}

# Index m of the design with total N in the allocation sequence (NA if none).
.bl_index_of_total <- function(N, r, allocation) {
  if (allocation == "total_nearest") return(if (N >= 2L) N - 1L else NA_integer_)
  total <- function(m) {
    d <- .bl_design_at(m, r, allocation)
    d$n1 + d$n2
  }
  m <- max(1L, as.integer(floor(N / (1 + r))))
  while (m > 1L && total(m) > N) m <- m - 1L
  while (total(m) < N) m <- m + 1L
  if (total(m) == N) m else NA_integer_
}

# Checks run once per a priori search (at its first design). Each failure is
# a typed status: the request is well formed but has no answer.
.bl_method_a_priori_checks <- function(p) {
  info <- .bl_variant(p$variant)
  if (p$power <= p$alpha) {
    .ig_abort("invalid_input", sprintf("the target power (%s) must exceed alpha (%s)",
                                       format(p$power), format(p$alpha)))
  }
  if ((info$variant == "less" && p$p1 >= p$p2) || (info$variant == "greater" && p$p1 <= p$p2)) {
    .ig_abort("invalid_input", sprintf(
      "the one-sided variant '%s' (H1: %s) does not match the alternative p1 = %s, p2 = %s",
      info$variant, info$alternative, format(p$p1), format(p$p2)))
  }
  if (p$p1 == p$p2) {
    .ig_abort("infeasible", sprintf(
      "with p1 = p2 the rejection probability is a null probability, at most alpha = %s < power = %s",
      format(p$alpha), format(p$power)))
  }
  r <- as.double(p$allocation_ratio)
  m_max <- .bl_max_index(p$n_max, r, p$allocation)
  if (m_max < 1L) return(invisible(TRUE))     # the search reports n_max itself
  d <- .bl_design_at(m_max, r, p$allocation)
  bound <- .bl_np_bound(d$n1, d$n2, p$p1, p$p2, p$alpha)
  if (bound < p$power - .bl_np_margin) {
    .ig_abort("infeasible", sprintf(
      "no level-%s test can reach power %s with n1 + n2 <= %d under this allocation (Neyman-Pearson bound %s at n1 = %d, n2 = %d)",
      format(p$alpha), format(p$power), as.integer(p$n_max), format(bound, digits = 6), d$n1, d$n2))
  }
  invisible(TRUE)
}

# A priori search of the generic solver (hook ig_method(a_priori_search =)).
# It runs the dedicated search of boschloo_sample_size(): the designs before
# the first one whose Neyman-Pearson power bound reaches the target are
# excluded for every level-alpha test, every later design of the allocation
# sequence is evaluated up to the first one reaching the target, and the
# search continues to the stable design (the next K = 10 designs also reach
# the target; no claim beyond that window). The answer has the shape of
# .ig_integer_search() results; the generic solver re-evaluates the design
# with .bl_method_power(), which uses the same cached rejection region.
.bl_method_a_priori_search <- function(p, n_lo, n_max) {
  ss <- boschloo_sample_size(p$p1, p$p2, power = p$power, alpha = p$alpha, variant = p$variant,
                             allocation_ratio = p$allocation_ratio, allocation = p$allocation,
                             n_max = min(as.integer(p$n_max), as.integer(n_max)),
                             control = .bl_method_control(p))
  np <- ss$neyman_pearson
  base <- list(search = "boschloo_sample_size: Neyman-Pearson exclusion, then enumeration of the allocation sequence",
               allocation = ss$policy$rule, n_upper_limit = ss$policy$n_max)
  if (!identical(ss$status$code, "ok")) {
    diag <- base
    if (!is.null(np) && !is.null(np$first_index) && !is.na(np$first_index)) {
      diag$neyman_pearson_first_n_total <- np$n_total
    }
    return(list(status = ss$status$code, reason = ss$status$reason, diagnostics = diag))
  }
  st <- ss$staircase
  mi <- ss$minimality
  tail_rows <- utils::tail(seq_len(nrow(st)), 10L)
  list(status = "ok", n = as.double(ss$first$n_total), power = ss$first$power,
       diagnostics = c(base, list(
         minimality = paste("proven: every earlier design of the allocation sequence was evaluated and failed,",
                            "or is excluded for every level-alpha test by the Neyman-Pearson bound"),
         monotonicity = "not assumed (discrete test)",
         designs_excluded_by_bound = as.double(mi$designs_excluded_by_bound),
         neyman_pearson_first_n_total = as.double(np$n_total),
         evaluations = as.double(mi$designs_evaluated),
         stability_window = as.double(mi$stability_window),
         stable_found = isTRUE(mi$stable_found)),
         if (!is.null(ss$stable)) list(stable_n_total = as.double(ss$stable$n_total)),
         list(power_sequence = list(n_total = as.double(st$n_total[tail_rows]),
                                    power = as.double(st$power[tail_rows]))),
         if (length(ss$notes)) list(notes = paste(ss$notes, collapse = "; "))))
}
