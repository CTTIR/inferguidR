# Generic mode solvers ------------------------------------------------------------------
#
# The five G*Power analysis types are solved generically from a method's
# forward model `power(p)` and allocation rule `design(p, n_total)` (see
# ig_method()). The solvers are registered through the hook
# .ig_mode_solvers(), which R/core-registry.R collects at load time;
# calculate_power() and solve_design() dispatch to them.
#
# Forward-model contract. `power(p)` receives the specification's inputs and
# options as one named list (plus the quantities a solver sets: `n_total`
# and the design's group sizes, `alpha`, or the effect parameter) and
# returns a named list with
#   power           achieved power (required; NA with a reason in `notes`)
#   beta            non-rejection probability under H1, computed directly
#                   (optional; 1 - power is used otherwise)
#   alpha_attained  attained size of the rejection rule (default: alpha)
#   critical        critical value(s) on the test-statistic scale
#   df, ncp         degrees of freedom and noncentrality (optional)
#   n_total, group_sizes   the design evaluated (optional)
#   statistic_distribution plot-data hint for H0 and H1 (optional list)
#   outputs         further method-specific numbers (optional named list)
#   notes           reasons for NA quantities (named list of strings)
#   assumptions     character vector (optional)
#   warnings        character vector without missing values (optional)
# A forward model or design rule may signal a typed inferguidR condition
# (for example reference_pending for a branch whose reference semantics are
# unknown); the solver then returns a result with that status.
#
# Result layout shared by all modes:
#   solved      the solve target(s) of the mode
#   achieved    power, beta and alpha_attained at the reported design (and
#               the realised beta/alpha ratio in compromise mode)
#   quantities  alpha_nominal, df, ncp, critical, n_total (unless solved)
#   group_sizes integer group sizes of the design
#   diagnostics solver evidence (bracket, iterations, residual, roots,
#               minimality, ...)

.ig_mode_solvers <- function() {
  list(
    a_priori = .ig_solve_a_priori,
    post_hoc = .ig_solve_post_hoc,
    sensitivity = .ig_solve_sensitivity,
    criterion = .ig_solve_criterion,
    compromise = .ig_solve_compromise
  )
}

# Shared helpers --------------------------------------------------------------------------

.ig_solver_inputs <- function(spec) c(spec$inputs, spec$options)

# Run a solver body; typed inferguidR conditions become results with their
# status, other errors propagate (the dispatcher reports numerical_failure).
.ig_solver_guard <- function(spec, entry, solver_id, expr) {
  tryCatch(expr, inferguidR_error = function(e) {
    .ig_power_result(spec, e$status, provenance = list(entry = entry, solver = solver_id),
                     r_call = .ig_result_call(spec, entry))
  })
}

.ig_result_call <- function(spec, entry) {
  fun <- if (identical(entry, "solve_design")) "solve_design" else "calculate_power"
  paste0("inferguidR::", fun, "(", .ig_spec_call(spec), ")")
}

# Evaluate and normalise the forward model.
.ig_forward <- function(method, p) {
  out <- method$power(p)
  if (!is.list(out) || !is.numeric(out$power) || length(out$power) != 1L) {
    stop(sprintf("the forward model of %s must return a list with one numeric `power`", method$id),
         call. = FALSE)
  }
  pw <- as.double(out$power)
  if (!is.na(pw) && (pw < -1e-12 || pw > 1 + 1e-12)) {
    stop(sprintf("the forward model of %s returned power %s outside [0, 1]", method$id, format(pw)),
         call. = FALSE)
  }
  if (!is.na(pw)) pw <- min(max(pw, 0), 1)
  if (!is.null(out$warnings) && (!is.character(out$warnings) || anyNA(out$warnings))) {
    stop(sprintf("the forward model of %s returned invalid `warnings`", method$id), call. = FALSE)
  }
  beta <- if (is.null(out$beta)) 1 - pw else as.double(out$beta)
  if (!is.na(beta)) beta <- min(max(beta, 0), 1)
  out$power <- pw
  out$beta <- beta
  if (is.null(out$alpha_attained)) out$alpha_attained <- if (is.null(p$alpha)) NA_real_ else p$alpha
  out
}

# Power only (for objective functions).
.ig_power_only <- function(method, p) .ig_forward(method, p)$power

# Build an ok result from a forward evaluation.
.ig_ok_result <- function(method, spec, entry, solver_id, p, fw, solved, achieved_extra = list(),
                          diagnostics = list(), warnings = character()) {
  achieved <- c(list(power = fw$power, beta = fw$beta, alpha_attained = as.double(fw$alpha_attained)),
                achieved_extra)
  quantities <- list()
  if (!is.null(p$alpha)) quantities$alpha_nominal <- as.double(p$alpha)
  for (nm in c("df", "ncp", "critical")) {
    if (!is.null(fw[[nm]])) quantities[[nm]] <- as.double(fw[[nm]])
  }
  if (!is.null(fw$n_total) && !("n_total" %in% spec$target)) quantities$n_total <- as.double(fw$n_total)
  undefined <- list()
  notes <- fw$notes %||% list()
  for (sec in list(achieved, quantities)) {
    for (nm in names(sec)) {
      if (anyNA(sec[[nm]])) {
        reason <- notes[[nm]]
        if (!.ig_is_string(reason)) {
          stop(sprintf("the forward model of %s left '%s' undefined without a reason", method$id, nm),
               call. = FALSE)
        }
        undefined[[nm]] <- reason
      }
    }
  }
  method_outputs <- fw$outputs %||% list()
  if (!is.null(fw$statistic_distribution)) {
    method_outputs$statistic_distribution <- fw$statistic_distribution
  }
  group_sizes <- if (is.null(fw$group_sizes)) NULL else {
    gs <- fw$group_sizes
    if (any(abs(gs - round(gs)) > 0)) {
      stop("group sizes of a reported design must be whole numbers", call. = FALSE)
    }
    stats::setNames(as.double(round(gs)), names(gs))
  }
  assumptions <- unique(c(as.character(fw$assumptions %||% character())))
  warnings <- unique(c(warnings, fw$warnings %||% character()))
  .ig_power_result(
    spec, ig_status("ok"), solved = solved, achieved = achieved, group_sizes = group_sizes,
    quantities = quantities, undefined = undefined, method_outputs = method_outputs,
    diagnostics = c(list(solver = solver_id), diagnostics),
    provenance = list(entry = entry, solver = solver_id,
                      engine = list(package = as.character(method$engine$package),
                                    fun = as.character(method$engine$fun))),
    assumptions = assumptions, warnings = warnings, r_call = .ig_result_call(spec, entry)
  )
}

.ig_status_result <- function(spec, entry, solver_id, code, reason, diagnostics = list()) {
  .ig_power_result(spec, ig_status(code, reason), diagnostics = c(list(solver = solver_id), diagnostics),
                   provenance = list(entry = entry, solver = solver_id),
                   r_call = .ig_result_call(spec, entry))
}

# Post hoc: power at the given design ------------------------------------------------------

.ig_solve_post_hoc <- function(method, spec, entry) {
  .ig_solver_guard(spec, entry, "post_hoc_direct", {
    p <- .ig_solver_inputs(spec)
    fw <- .ig_forward(method, p)
    if (is.na(fw$power)) {
      return(.ig_status_result(spec, entry, "post_hoc_direct", "numerical_failure",
                               paste("the power could not be computed:",
                                     fw$notes$power %||% "the forward model returned NA")))
    }
    .ig_ok_result(method, spec, entry, "post_hoc_direct", p, fw, solved = list(power = fw$power),
                  diagnostics = list(method = "direct evaluation of the forward model"))
  })
}

# A priori: smallest admissible integer design ----------------------------------------------

.ig_solve_a_priori <- function(method, spec, entry) {
  hook <- method$a_priori_search
  solver_id <- if (is.function(hook)) "a_priori_method_search" else "a_priori_integer_search"
  .ig_solver_guard(spec, entry, solver_id, {
    p0 <- .ig_solver_inputs(spec)
    target <- p0$power
    dom <- method$params$n_total$domain
    n_lo <- max(1, ceiling(dom$lower))
    n_max <- min(.ig_solver_settings$n_max, floor(dom$upper))
    design_fun <- method$design
    design_at <- function(N, continuous = FALSE) {
      if (is.null(design_fun)) return(NULL)
      gs <- if (continuous) design_fun(p0, N, continuous = TRUE) else design_fun(p0, N)
      if (is.null(gs)) return(NULL)
      if (!is.numeric(gs) || is.null(names(gs)) || anyNA(gs)) {
        stop("the design rule must return named group sizes or NULL", call. = FALSE)
      }
      if (!continuous && (any(gs != round(gs)) || abs(sum(gs) - N) > 0)) {
        stop(sprintf("the design rule returned group sizes that are not whole numbers adding up to N = %s",
                     format(N, scientific = FALSE)), call. = FALSE)
      }
      gs
    }
    point <- function(N, continuous = FALSE) {
      p <- p0
      p$n_total <- N
      gs <- design_at(N, continuous)
      if (!is.null(gs)) for (nm in names(gs)) p[[nm]] <- gs[[nm]]
      p
    }
    admissible <- function(N) is.null(design_fun) || !is.null(design_at(N))

    # Continuous relaxation (when the design rule supports fractional group
    # sizes): used only as a starting point and reported as a diagnostic.
    relax <- list()
    n_start <- NA_real_
    if (!is.function(hook) && !is.null(design_fun) && "continuous" %in% names(formals(design_fun)) &&
        !method$discrete) {
      relax <- tryCatch({
        f <- function(u) .ig_power_only(method, point(exp(u), continuous = TRUE)) - target
        br <- .ig_bracket_monotone(f, log(max(n_lo, 2)), log(n_lo), log(n_max), decreasing = FALSE)
        if (isTRUE(br$found)) {
          r <- if (isTRUE(br$exact)) list(x = br$a, kind = "crossing") else
            .ig_refine_root(f, br$a, br$b, br$fa, br$fb)
          if (r$kind %in% c("crossing", "steep_crossing")) list(n_total = exp(r$x)) else list()
        } else list()
      }, inferguidR_error = function(e) list(), error = function(e) list())
      if (length(relax)) n_start <- relax$n_total
    }

    # A method may bring its own search (ig_method(a_priori_search =)), for
    # example a discrete test that excludes designs by a proven bound; its
    # answer is checked here and the design is re-evaluated below.
    search <- if (is.function(hook)) {
      .ig_check_a_priori_search(hook(p0, n_lo = n_lo, n_max = n_max), method, n_lo, n_max,
                                admissible)
    } else {
      .ig_integer_search(
        power_at = function(N) .ig_power_only(method, point(N)),
        admissible = admissible, n_lo = n_lo, n_max = n_max, target = target,
        discrete = method$discrete, n_start = n_start
      )
    }
    if (!identical(search$status, "ok")) {
      return(.ig_status_result(spec, entry, solver_id, search$status, search$reason, search$diagnostics))
    }
    p <- point(search$n)
    fw <- .ig_forward(method, p)
    diag <- search$diagnostics
    diag$target_power <- target
    if (length(relax)) {
      diag$continuous_relaxation <- list(
        n_total = relax$n_total,
        note = "root of the power curve with fractional group sizes; not an admissible design, used only to start the integer search"
      )
    }
    diag$allocation_rule <- method$reference$allocation_rule %||%
      (if (is.null(design_fun)) "single sample: N is the design" else "method design rule")
    if (is.function(hook) && (is.na(fw$power) || fw$power < target)) {
      stop(sprintf("the a priori search of %s returned N = %s, whose recomputed power %s does not reach the target %s",
                   method$id, format(search$n, scientific = FALSE), format(fw$power), format(target)),
           call. = FALSE)
    }
    .ig_ok_result(method, spec, entry, solver_id, p, fw, solved = list(n_total = search$n),
                  diagnostics = diag, warnings = search$warnings %||% character())
  })
}

# Check the answer of a method-specific a priori search: the shape of
# .ig_integer_search() results, a known status code, and for "ok" an
# admissible whole-number design inside [n_lo, n_max]. Violations are
# programming errors (the dispatcher reports numerical_failure).
.ig_check_a_priori_search <- function(search, method, n_lo, n_max, admissible) {
  bad <- function(what) {
    stop(sprintf("the a priori search of %s returned %s", method$id, what), call. = FALSE)
  }
  if (!is.list(search) || !.ig_is_string(search$status) || !(search$status %in% ig_status_codes())) {
    bad("no list with a valid `status`")
  }
  if (!is.null(search$diagnostics) && !is.list(search$diagnostics)) bad("diagnostics that are not a list")
  if (!is.null(search$warnings) && (!is.character(search$warnings) || anyNA(search$warnings))) {
    bad("warnings that are not a character vector without missing values")
  }
  if (!identical(search$status, "ok")) {
    if (!.ig_is_string(search$reason)) bad(sprintf("status '%s' without a `reason`", search$status))
    return(search)
  }
  n <- search$n
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n != round(n)) {
    bad("no whole-number `n`")
  }
  if (n < n_lo || n > n_max) {
    bad(sprintf("N = %s outside the admissible range [%s, %s]", format(n, scientific = FALSE),
                format(n_lo, scientific = FALSE), format(n_max, scientific = FALSE)))
  }
  if (!admissible(n)) bad(sprintf("N = %s, which the design rule does not admit", format(n, scientific = FALSE)))
  search$n <- as.double(n)
  search
}

# Sensitivity: effect size reaching the target power ---------------------------------------

.ig_solve_sensitivity <- function(method, spec, entry) {
  if (is.function(method$sensitivity_search)) return(.ig_solve_sensitivity_hook(method, spec, entry))
  solver_id <- "sensitivity_scan_root"
  .ig_solver_guard(spec, entry, solver_id, {
    p0 <- .ig_solver_inputs(spec)
    direction <- method$effect_direction
    if (is.function(direction)) direction <- direction(p0)
    if (!.ig_is_string(direction) || !(direction %in% .ig_effect_directions)) {
      stop("method effect_direction must return increasing, decreasing or unknown", call. = FALSE)
    }
    target <- p0$power
    eff <- spec$target
    dom <- method$params[[eff]]$domain
    if (is.function(method$sensitivity_domain)) {
      supplied <- method$sensitivity_domain(p0)
      if (!inherits(supplied, "ig_domain") || !identical(supplied$type, "real") ||
          !identical(supplied$length, 1L)) {
        stop("sensitivity_domain must return a scalar real domain", call. = FALSE)
      }
      narrowed <- ig_domain_real(supplied$lower, supplied$upper,
        supplied$lower_inclusive, supplied$upper_inclusive, scan = supplied$scan)
      if (narrowed$lower >= narrowed$upper || narrowed$lower < dom$lower ||
          narrowed$upper > dom$upper ||
          (narrowed$lower == dom$lower && narrowed$lower_inclusive && !dom$lower_inclusive) ||
          (narrowed$upper == dom$upper && narrowed$upper_inclusive && !dom$upper_inclusive)) {
        stop("sensitivity_domain must be a nonempty interval within the effect parameter domain", call. = FALSE)
      }
      dom <- narrowed
    }
    at <- function(e) {
      p <- p0
      p[[eff]] <- e
      p
    }
    f <- function(e) .ig_power_only(method, at(e)) - target
    grid <- .ig_domain_grid(dom)
    scan <- .ig_scan_roots(f, grid)
    roots <- scan$roots
    root_x <- vapply(roots, `[[`, numeric(1), "x")
    root_slope <- vapply(roots, function(r) as.double(r$slope), numeric(1))
    base_diag <- list(scan_grid = list(points = length(grid), from = min(grid), to = max(grid)),
                      evaluations = scan$evaluations, failed_evaluations = scan$failed)
    if (length(roots)) base_diag$roots_found <- root_x
    if (length(scan$flat)) base_diag$flat_regions <- lapply(scan$flat, as.double)

    vals <- scan$values
    ok <- !is.na(vals)
    if (!any(ok)) {
      return(.ig_status_result(spec, entry, solver_id, "numerical_failure",
                               "the power could not be evaluated anywhere on the effect-size scan",
                               base_diag))
    }
    v <- vals[ok]
    g <- grid[ok]
    pw <- v + target

    # A target within the residual tolerance of the smallest power on the
    # scan cannot be separated from that minimum: every effect size near it
    # meets the tolerance. If the target equals the minimum up to rounding
    # (for example a target power equal to alpha, the power at no effect) it
    # is touched but never crossed, and sign changes found there are rounding
    # noise. If it lies measurably above the minimum, a solution next to the
    # minimum exists but cannot be resolved, which is a numerical limit, not
    # infeasibility.
    gap <- -min(v)
    if (abs(gap) <= .ig_solver_settings$residual_tol) {
      at_min <- sprintf("%s at %s = %s", format(min(pw), digits = 15), eff,
                        format(signif(g[which.min(v)], 7)))
      if (gap <= .ig_solver_settings$power_noise) {
        return(.ig_status_result(spec, entry, solver_id, "infeasible", sprintf(paste(
          "the target power %s equals the smallest power on the scan (%s) up to rounding and is touched",
          "there but not crossed, so no effect size on a rising or falling part of the power curve has the",
          "target power; choose a target above that minimum"), format(target, digits = 15), at_min),
          base_diag))
      }
      return(.ig_status_result(spec, entry, solver_id, "numerical_failure", sprintf(paste(
        "the target power %s exceeds the smallest power on the scan (%s) by only %s, less than the",
        "residual tolerance %s, so a solution next to that minimum cannot be separated from it numerically;",
        "choose a target further above the minimum"), format(target, digits = 15), at_min,
        format(signif(gap, 3)), format(.ig_solver_settings$residual_tol)), base_diag))
    }

    if (!length(roots)) {
      if (any(v < 0) && any(v > 0)) {
        return(.ig_status_result(spec, entry, solver_id, "numerical_failure", sprintf(
          "the power crosses the target only between scan points where it could not be evaluated (%d failed evaluations)",
          scan$failed), base_diag))
      }
      reason <- if (all(v < 0)) {
        sprintf("the target power %s is not reached for any %s on the scan from %s to %s (largest power %s at %s = %s)",
                format(target), eff, format(min(grid)), format(max(grid)), format(signif(max(pw), 7)),
                eff, format(signif(g[which.max(pw)], 7)))
      } else {
        sprintf("the power exceeds the target %s for every %s on the scan (smallest power %s at %s = %s), so no effect size has exactly the target power",
                format(target), eff, format(signif(min(pw), 7)), eff, format(signif(g[which.min(pw)], 7)))
      }
      return(.ig_status_result(spec, entry, solver_id, "infeasible", reason, base_diag))
    }

    branch <- switch(direction,
                     increasing = which(root_slope > 0),
                     decreasing = which(root_slope < 0),
                     seq_along(roots))
    rule <- switch(direction,
                   increasing = "power increases with the effect: the smallest root at which power is increasing",
                   decreasing = "power decreases with the effect: the largest root at which power is decreasing",
                   "direction unknown: the unique root")
    if (!length(branch)) {
      return(.ig_status_result(spec, entry, solver_id, "infeasible", sprintf(
        "the target power is reached only where power %s with %s, not on the declared branch (%s)",
        if (direction == "increasing") "does not increase" else "does not decrease", eff, rule), base_diag))
    }
    if (direction == "unknown" && length(branch) > 1L) {
      return(.ig_status_result(spec, entry, solver_id, "numerical_failure", sprintf(
        "%d roots found and the method declares no effect direction to choose between them",
        length(branch)), base_diag))
    }
    pick <- if (direction == "decreasing") branch[which.max(root_x[branch])] else branch[which.min(root_x[branch])]
    r <- roots[[pick]]
    warn <- character()
    if (length(branch) > 1L) {
      warn <- sprintf("%d roots on the selected branch; the %s one is reported (see diagnostics$roots_found)",
                      length(branch), if (direction == "decreasing") "largest" else "smallest")
    }
    x <- r$x
    note <- NULL
    if (r$kind == "jump") {
      if (!method$discrete) {
        return(.ig_status_result(spec, entry, solver_id, "numerical_failure",
                                 .ig_unresolved_jump_reason(sprintf("the power as a function of %s", eff)),
                                 base_diag))
      }
      x <- .ig_snap_to_attaining(f, r)
      note <- "discrete test: power jumps across the target; the reported effect is the boundary value whose power is at least the target"
    } else if (r$kind == "unresolved") {
      return(.ig_status_result(spec, entry, solver_id, "numerical_failure", sprintf(
        "root refinement for %s did not reach the residual tolerance (residual %s)", eff,
        format(signif(r$residual, 3))), base_diag))
    }
    p <- at(x)
    fw <- .ig_forward(method, p)
    diag <- c(base_diag, list(
      branch_rule = rule, bracket = r$bracket, iterations = r$iterations,
      residual = fw$power - target, root_kind = r$kind, target_power = target))
    if (!is.null(note)) diag$discrete_note <- note
    if (r$kind == "steep_crossing") diag$steep_note <- .ig_steep_note(eff)
    if (direction != "unknown" && length(roots) > length(branch)) {
      diag$other_branch_roots <- root_x[-branch]
    }
    .ig_ok_result(method, spec, entry, solver_id, p, fw, solved = stats::setNames(list(x), eff),
                  diagnostics = diag, warnings = warn)
  })
}

# Validate a method's sensitivity answer and recompute its forward power.
.ig_solve_sensitivity_hook <- function(method, spec, entry) {
  solver_id <- "sensitivity_method_search"
  .ig_solver_guard(spec, entry, solver_id, {
    p <- .ig_solver_inputs(spec)
    target <- p$power
    effect <- spec$target
    answer <- method$sensitivity_search(p, effect)
    bad <- function(reason) stop(sprintf("the sensitivity search of %s returned %s", method$id, reason), call. = FALSE)
    if (!is.list(answer) || !.ig_is_string(answer$status) || !(answer$status %in% ig_status_codes())) {
      bad("no list with a valid status")
    }
    diagnostics <- answer$diagnostics %||% list()
    warnings <- answer$warnings %||% character()
    if (!is.list(diagnostics)) bad("invalid diagnostics")
    if (!is.character(warnings) || anyNA(warnings)) bad("invalid warnings")
    if (!is.null(diagnostics$target_reached) && !.ig_is_flag(diagnostics$target_reached)) {
      bad("an invalid target_reached flag")
    }
    if (!identical(answer$status, "ok")) {
      if (!.ig_is_string(answer$reason)) bad("a failure without a reason")
      return(.ig_status_result(spec, entry, solver_id, answer$status, answer$reason, diagnostics))
    }
    value <- answer$value
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        !is.null(.ig_domain_check(method$params[[effect]]$domain, value, effect))) {
      bad("an effect outside its parameter domain")
    }
    p[[effect]] <- value
    fw <- .ig_forward(method, p)
    if (!is.finite(fw$power)) bad("an effect whose power cannot be evaluated")
    residual <- fw$power - target
    if (identical(diagnostics$target_reached, FALSE)) {
      maximum <- diagnostics$max_power
      if (!length(warnings) || !any(nzchar(trimws(warnings))) ||
          !is.numeric(maximum) || length(maximum) != 1L || !is.finite(maximum) ||
          abs(maximum - fw$power) > .ig_solver_settings$residual_tol || fw$power >= target) {
        bad("an unattained target without a checked maximum and explanatory warning")
      }
    } else {
      if (abs(residual) > .ig_solver_settings$residual_tol) bad("an effect that does not reach the target")
      diagnostics$target_reached <- TRUE
    }
    diagnostics$target_power <- target
    diagnostics$residual <- residual
    .ig_ok_result(method, spec, entry, solver_id, p, fw,
                  solved = stats::setNames(list(value), effect),
                  diagnostics = diagnostics, warnings = warnings)
  })
}

# For a jump root of a discrete method, the boundary point of the bracket on
# the side where f >= 0 (power at least the target), located by bisection
# down to adjacent doubles. Whichever side of the jump that is (smallest
# attaining value for power increasing in x, largest for decreasing), the
# returned point is the attaining one.
.ig_snap_to_attaining <- function(f, r) {
  a <- r$bracket[1]
  b <- r$bracket[2]
  fa <- r$f_bracket[1]
  # Keep a on the failing side and b on the attaining side.
  if (fa >= 0) {
    tmp <- a
    a <- b
    b <- tmp
  }
  for (i in seq_len(200L)) {
    m <- (a + b) / 2
    if (m == a || m == b) break
    if (.ig_eval_finite(f, m) >= 0) b <- m else a <- m
  }
  b
}

# Alpha as a solve target ------------------------------------------------------------------
#
# Criterion and compromise solve for alpha on the logit scale
# u = log(alpha / (1 - alpha)). It behaves like log(alpha) for small alpha
# and like -log(1 - alpha) for alpha near one, so the objectives stay well
# conditioned at both ends of (0, 1) (on the log(alpha) scale alone the
# compromise objective becomes so steep near alpha = 1 that the residual
# check cannot be met). The search range is fixed and reported.
.ig_alpha_search <- list(lower = 1e-300, upper = 1 - 1e-12, text = "[1e-300, 1 - 1e-12]")

.ig_alpha_u <- function(alpha) stats::qlogis(alpha)

.ig_u_alpha <- function(u) stats::plogis(u)

.ig_alpha_range_text <- function() .ig_alpha_search$text

.ig_steep_note <- function(what) {
  sprintf(paste("the objective is steep at the solution: its residual is limited by the",
                "double-precision resolution of %s, not by the root finder (a sign change was",
                "confirmed on both sides of the reported value)"), what)
}

# A sign change without a continuous trend at the solution of a continuous
# method: a genuine discontinuity, or rounding noise where the forward model
# has lost its relative precision (for example a beta of 1e-9 formed from
# tail probabilities near 0.5). Either way no number is reported.
.ig_unresolved_jump_reason <- function(what) {
  sprintf(paste("%s changes sign at the solution without a continuous trend (a discontinuity, or",
                "rounding noise where the forward model has lost relative precision), so the solution",
                "of this continuous method cannot be resolved"), what)
}

# Criterion: alpha reaching the target power ----------------------------------------------

.ig_solve_criterion <- function(method, spec, entry) {
  solver_id <- "criterion_alpha_root"
  .ig_solver_guard(spec, entry, solver_id, {
    p0 <- .ig_solver_inputs(spec)
    target <- p0$power
    at <- function(u) {
      p <- p0
      p$alpha <- .ig_u_alpha(u)
      p
    }
    # Power is non-decreasing in alpha for nested rejection regions, so the
    # root is bracketed by expansion on the logit(alpha) scale.
    f <- function(u) .ig_power_only(method, at(u)) - target
    lower <- .ig_alpha_u(.ig_alpha_search$lower)
    upper <- .ig_alpha_u(.ig_alpha_search$upper)
    br <- .ig_bracket_monotone(f, .ig_alpha_u(0.05), lower, upper, decreasing = FALSE)
    diag <- list(scale = "logit(alpha)", search_range = c(.ig_alpha_search$lower, .ig_alpha_search$upper),
                 bracket_trail = .ig_u_alpha(br$trail), target_power = target,
                 monotonicity = "power is non-decreasing in alpha (nested rejection regions)")
    if (!isTRUE(br$found)) {
      reason <- if (!is.null(br$value_at_limit) && !is.na(br$value_at_limit) && br$value_at_limit < 0) {
        sprintf("the target power %s is not reached for any alpha in the search range %s (power at the upper end: %s)",
                format(target), .ig_alpha_range_text(), format(signif(br$value_at_limit + target, 7)))
      } else if (!is.null(br$value_at_limit) && !is.na(br$value_at_limit)) {
        sprintf("the power exceeds the target %s even for alpha = %s, the lower end of the search range %s",
                format(target), format(.ig_alpha_search$lower), .ig_alpha_range_text())
      } else {
        paste("no bracket for alpha:", br$reason)
      }
      return(.ig_status_result(spec, entry, solver_id, "infeasible", reason, diag))
    }
    r <- if (isTRUE(br$exact)) {
      list(x = br$a, kind = "crossing", bracket = c(br$a, br$a), iterations = 0L)
    } else {
      .ig_refine_root(f, br$a, br$b, br$fa, br$fb)
    }
    u <- r$x
    if (r$kind == "jump") {
      if (!method$discrete) {
        return(.ig_status_result(spec, entry, solver_id, "numerical_failure",
                                 .ig_unresolved_jump_reason("the power as a function of alpha"), diag))
      }
      u <- .ig_snap_to_attaining(f, r)
      diag$discrete_note <- paste(
        "discrete test: power is a step function of the nominal alpha; the reported alpha is the smallest",
        "nominal alpha whose rejection rule reaches the target power (attained size in achieved$alpha_attained)")
    } else if (r$kind == "unresolved") {
      return(.ig_status_result(spec, entry, solver_id, "numerical_failure",
                               "root refinement for alpha did not reach the residual tolerance", diag))
    } else if (r$kind == "steep_crossing") {
      diag$steep_note <- .ig_steep_note("alpha")
    }
    p <- at(u)
    fw <- .ig_forward(method, p)
    diag$bracket <- .ig_u_alpha(r$bracket)
    diag$iterations <- r$iterations
    diag$residual <- fw$power - target
    diag$root_kind <- r$kind
    .ig_ok_result(method, spec, entry, solver_id, p, fw, solved = list(alpha = p$alpha), diagnostics = diag)
  })
}

# Compromise: alpha and beta with beta = q * alpha ------------------------------------------

.ig_solve_compromise <- function(method, spec, entry) {
  solver_id <- "compromise_ratio_root"
  .ig_solver_guard(spec, entry, solver_id, {
    p0 <- .ig_solver_inputs(spec)
    q <- p0$error_ratio
    at <- function(u) {
      p <- p0
      p$alpha <- .ig_u_alpha(u)
      p
    }
    # beta(alpha) / (q alpha) - 1 is decreasing in alpha (beta is
    # non-increasing), finite for every alpha in (0, 1), positive as alpha ->
    # 0 and -1 at alpha = 1 (beta = 0); the root is unique.
    f <- function(u) {
      p <- at(u)
      .ig_forward(method, p)$beta / (q * p$alpha) - 1
    }
    lower <- .ig_alpha_u(.ig_alpha_search$lower)
    upper <- .ig_alpha_u(.ig_alpha_search$upper)
    start <- .ig_alpha_u(min(0.5, 1 / (1 + q)))
    br <- .ig_bracket_monotone(f, start, lower, upper, decreasing = TRUE)
    diag <- list(scale = "logit(alpha)", search_range = c(.ig_alpha_search$lower, .ig_alpha_search$upper),
                 error_ratio = q, bracket_trail = .ig_u_alpha(br$trail),
                 objective = "beta(alpha) / (q * alpha) - 1",
                 monotonicity = "beta is non-increasing in alpha (nested rejection regions)")
    if (!isTRUE(br$found)) {
      end <- if (!is.null(br$limit)) format(signif(.ig_u_alpha(br$limit), 7)) else "?"
      return(.ig_status_result(spec, entry, solver_id, "infeasible", sprintf(
        "no alpha in the search range %s satisfies beta = q * alpha with q = %s (%s; searched up to alpha = %s)",
        .ig_alpha_range_text(), format(q), br$reason, end), diag))
    }
    r <- if (isTRUE(br$exact)) {
      list(x = br$a, kind = "crossing", bracket = c(br$a, br$a), iterations = 0L)
    } else {
      .ig_refine_root(f, br$a, br$b, br$fa, br$fb)
    }
    u <- r$x
    # A small residual on the failing side can conceal a discrete jump:
    # beta/(q*alpha) may approach one immediately before a rejection rule
    # changes. Locate the attaining side instead of accepting that limit.
    if (method$discrete && !isTRUE(br$exact) && r$kind != "unresolved" && f(u) > 0) {
      u <- .ig_snap_to_attaining(function(x) -f(x),
        list(bracket = r$bracket, f_bracket = -r$f_bracket))
      if (abs(f(u)) > .ig_solver_settings$residual_tol) r$kind <- "jump"
    }
    if (r$kind == "jump") {
      if (!method$discrete) {
        return(.ig_status_result(spec, entry, solver_id, "numerical_failure",
                                 .ig_unresolved_jump_reason("beta / (q * alpha)"), diag))
      }
      # Smallest nominal alpha whose rule satisfies beta <= q * alpha.
      u <- .ig_snap_to_attaining(function(x) -f(x), list(bracket = r$bracket, f_bracket = -r$f_bracket))
      diag$discrete_note <- paste(
        "discrete test: beta is a step function of the nominal alpha and beta = q * alpha has no exact solution;",
        "the reported alpha is the smallest nominal alpha whose rejection rule has beta <= q * alpha.",
        "Nominal alpha (solved$alpha), attained size (achieved$alpha_attained) and the realised ratio",
        "beta / attained alpha (achieved$error_ratio) are reported separately.")
    } else if (r$kind == "unresolved") {
      return(.ig_status_result(spec, entry, solver_id, "numerical_failure",
                               "root refinement for the compromise alpha did not reach the residual tolerance",
                               diag))
    } else {
      if (r$kind == "steep_crossing") diag$steep_note <- .ig_steep_note("alpha")
      if (method$discrete) {
        # Between two attainable rejection rules beta is constant while the
        # nominal alpha moves, so beta = q * alpha is met on the nominal scale.
        diag$discrete_note <- paste(
          "discrete test: beta = q * alpha is solved on the nominal alpha; the rejection rule's attained size",
          "(achieved$alpha_attained) is not above it, and the realised ratio beta / attained alpha",
          "(achieved$error_ratio) is reported separately.")
      }
    }
    p <- at(u)
    fw <- .ig_forward(method, p)
    attained <- as.double(fw$alpha_attained)
    ratio <- if (is.finite(attained) && attained > 0) fw$beta / attained else NA_real_
    if (is.na(ratio)) {
      fw$notes$error_ratio <- paste("the rejection rule at the solved alpha has attained size 0 (it never rejects),",
                                    "so beta / attained alpha is undefined")
    }
    diag$bracket <- .ig_u_alpha(r$bracket)
    diag$iterations <- r$iterations
    diag$residual <- fw$beta / (q * p$alpha) - 1
    diag$root_kind <- r$kind
    .ig_ok_result(method, spec, entry, solver_id, p, fw,
                  solved = list(alpha = p$alpha, power = fw$power),
                  achieved_extra = list(error_ratio = ratio),
                  diagnostics = diag)
  })
}
