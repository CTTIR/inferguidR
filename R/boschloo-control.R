# Boschloo's test: variants, numerical controls and the result cache ------------
#
# Boschloo's test compares two independent binomial samples X1 ~ Bin(n1, p1)
# (group 1) and X2 ~ Bin(n2, p2) (group 2) under the equality null p1 = p2.
# Tables (x1, x2) are ordered by a Fisher conditional p-value, and the
# p-value is the supremum over the common nuisance probability of the null
# probability of the tables that are at least as extreme. Every output names
# the variant, the direction in terms of p1 and p2, the tie tolerance, the
# nuisance grid and the error bound (contract in specs/method-notes/boschloo.md).
#
# Variants (fixed group order: group 1 is (x1, n1), group 2 is (x2, n2); x1 and
# x2 count events):
#
# * "minlike" (variant id boschloo_minlike): two-sided, ordered by the
#   two-sided Fisher p-value that sums conditional probabilities not larger
#   than that of the observed table (as stats::fisher.test, Exact
#   tsmethod = "square", exact2x2 tsmethod = "minlike").
# * "central" (boschloo_central): two-sided, min(1, 2 * min(p_less, p_greater))
#   with each one-sided Boschloo p-value maximised separately.
# * "less" (boschloo_lower): one-sided, H1: p1 < p2.
# * "greater" (boschloo_upper): one-sided, H1: p1 > p2.

.bl_variant_info <- list(
  minlike = list(
    variant = "minlike", variant_id = "boschloo_minlike",
    advisor_variant_id = "boschloo_two_sided_minlike",
    sides = "minlike", two_sided = TRUE, alternative = "p1 != p2",
    label = "Boschloo test, two-sided, minlike ordering (two-sided Fisher p-value)",
    ordering_statistic = paste(
      "two-sided Fisher p-value: sum of the conditional null probabilities",
      "P(X1 = y | t) that do not exceed P(X1 = x1 | t) (relative tolerance applied)"),
    construction = "tail set ordered by the two-sided Fisher p-value",
    exact_call = "Exact::exact.test(..., method = \"boschloo\", alternative = \"two.sided\", tsmethod = \"square\")",
    exact2x2_call = "exact2x2::boschloo(..., alternative = \"two.sided\", tsmethod = \"minlike\")"),
  central = list(
    variant = "central", variant_id = "boschloo_central",
    advisor_variant_id = "boschloo_two_sided_central",
    sides = c("lower", "upper"), two_sided = TRUE, alternative = "p1 != p2",
    label = "Boschloo test, two-sided, central (twice the smaller one-sided p-value)",
    ordering_statistic = paste(
      "one-sided Fisher p-values P(X1 <= x1 | t) and P(X1 >= x1 | t),",
      "each maximised separately"),
    construction = "p = min(1, 2 * min(p_less, p_greater))",
    exact_call = paste(
      "Exact::exact.test(..., method = \"boschloo\", alternative = \"two.sided\", tsmethod = \"central\")",
      "(doubles the one-sided p in the observed direction; differs when the sample proportions are equal)"),
    exact2x2_call = "exact2x2::boschloo(..., alternative = \"two.sided\", tsmethod = \"central\")"),
  less = list(
    variant = "less", variant_id = "boschloo_lower",
    advisor_variant_id = "boschloo_one_sided",
    sides = "lower", two_sided = FALSE, alternative = "p1 < p2",
    label = "Boschloo test, one-sided, H1: p1 < p2",
    ordering_statistic = "one-sided Fisher p-value P(X1 <= x1 | t)",
    construction = "tail set ordered by the lower one-sided Fisher p-value",
    exact_call = "Exact::exact.test(..., method = \"boschloo\", alternative = \"less\")",
    exact2x2_call = "exact2x2::boschloo(..., alternative = \"greater\") (group 2 versus group 1)"),
  greater = list(
    variant = "greater", variant_id = "boschloo_upper",
    advisor_variant_id = "boschloo_one_sided",
    sides = "upper", two_sided = FALSE, alternative = "p1 > p2",
    label = "Boschloo test, one-sided, H1: p1 > p2",
    ordering_statistic = "one-sided Fisher p-value P(X1 >= x1 | t)",
    construction = "tail set ordered by the upper one-sided Fisher p-value",
    exact_call = "Exact::exact.test(..., method = \"boschloo\", alternative = \"greater\")",
    exact2x2_call = "exact2x2::boschloo(..., alternative = \"less\") (group 2 versus group 1)")
)

.bl_variants <- names(.bl_variant_info)

# Evidence status of the Boschloo engine (shared by all variants). The
# independent checker (validation/independent/boschloo/) shares no numerical
# routine with the engine: exact rational ordering statistics, Pascal-rule
# binomial coefficients and its own grid plus golden-section maximisation.
# G*Power has no Boschloo procedure. The classification of the related
# G*Power procedures is taken from the manual and the option resources of the
# binaries only (no application output was used) and stays provisional until
# a capture of their options dialogs confirms the runtime option lists
# (specs/method-notes/boschloo.md, section 11).
.bl_evidence <- list(
  validation_status = "independently_checked",
  checker = "validation/independent/boschloo",
  planning_compatibility = "independently_validated_extension",
  gpower_relation = "gpower_comparator_different_method",
  gpower_classification_status = "provisional_pending_capture",
  gpower_note = paste(
    "G*Power 3.1 has no Boschloo procedure. Its unconditional two-proportion test",
    "(exact_proportions_unconditional, no manual chapter) orders tables by z (pooled or",
    "unpooled, with or without continuity correction), Mantel-Haenszel, likelihood-ratio,",
    "t or score statistics according to the option resources of the binaries; none of them",
    "is the Fisher p-value ordering. Fisher's test (manual chapter 6) is the conditional test.",
    "Both are labelled comparators only; the runtime option lists are unverified until the",
    "options dialogs are captured in the application.")
)

# Largest total sample size for which all binomial coefficients C(N, k) are
# exactly representable as doubles (C(56, 28) < 2^53 < C(57, 28)); the exact
# integer fallback of the ordering is available up to this size.
.bl_exact_max_n <- 56L

# Largest number of tables (n1 + 1) * (n2 + 1) that is enumerated.
.bl_max_tables <- 4e6

#' Numerical settings for Boschloo's test
#'
#' Collects the settings that affect Boschloo p-values, rejection regions,
#' power and sample sizes. They are recorded in every result and are part of
#' every cache key, so results computed with different settings are never
#' mixed.
#'
#' The nuisance supremum `sup_pi P_pi(A)` of a tail set `A` is written in
#' Bernstein form `sum_k c_k dbinom(k, N, pi)`, where `c_k` is the conditional
#' null probability of `A` given `k` events in total. It is maximised by a
#' global grid of `grid_points` equally spaced values of `pi` (both endpoints
#' included), followed by branch-and-bound bisection: every grid interval is
#' bounded with a rigorous local bound of the second derivative of the
#' Bernstein polynomial, and intervals whose bound exceeds the best value
#' found by more than `nuisance_tol` are bisected. A final Brent step
#' ([stats::optimize()], tolerance `polish_tol`) polishes the maximiser. Each
#' supremum is reported as a pair: the refined value (`p_value`, a lower
#' bound attained at a computed point) and a rigorous upper bound
#' (`p_upper`), with `p_upper - p_value <= nuisance_tol` plus a
#' floating-point allowance. Tables are rejected only when `p_upper <= alpha`.
#'
#' Ordering statistics are compared with the relative tolerance
#' `tie_tolerance`: values whose relative difference is at most this
#' tolerance are treated as ties (both tables enter the tail set together).
#' Values above 0.5 are compared on the complementary tail, which is computed
#' directly so that values near 1 keep their relative accuracy. A pair of
#' distinct ordering values whose relative gap falls in `ambiguity_band` is
#' flagged as an ambiguous tie. For `n1 + n2 <= 56` the ordering is then
#' recomputed in exact integer arithmetic; above that size a result whose
#' tail set depends on the ambiguous pair is refused (status
#' `"numerical_failure"`).
#'
#' Gaps at or below the lower edge of `ambiguity_band` are taken as exact
#' ties, and gaps above the upper edge as exact distinctions, without a flag.
#' The tie tolerance therefore has to lie inside the band: a tolerance below
#' the lower edge would split exactly tied tables whose computed values
#' differ by rounding (for 1/5 versus 6/7 with `variant = "less"` it would
#' give p = 0.0174 instead of 0.0232), and a tolerance above the upper edge
#' would silently merge distinct values. For the same reason the lower edge
#' is at most `1e-12` (distinct values come as close as about `4e-10` at
#' `(500, 1500)`), and the upper edge is at least `1e-10` (exactly tied values
#' are computed at most about `4.6e-13` apart up to `n1 + n2 = 2000`).
#'
#' @param tie_tolerance Relative tie tolerance for ordering statistics
#'   (default `1e-7`, the tolerance used by [stats::fisher.test()]). It must
#'   lie in `[ambiguity_band[1], ambiguity_band[2]]`.
#' @param grid_points Number of points of the global nuisance grid on
#'   `[0, 1]` (default `1001`).
#' @param nuisance_tol Absolute tolerance of the rigorous supremum bound
#'   (default `1e-9`).
#' @param polish_tol Tolerance of the final Brent refinement of the
#'   maximiser (default `1e-12`).
#' @param max_levels Maximum number of bisection levels (default `60`).
#' @param ambiguity_band Relative gaps between distinct ordering values that
#'   are flagged as ambiguous (default `c(1e-12, 1e-6)`, a half-open interval
#'   `(lower, upper]` with `lower <= 1e-12` and `1e-10 <= upper <= 1e-3`).
#' @param cache Cache orderings and rejection regions in memory (keyed by
#'   design, variant, alpha and all settings above).
#' @return An object of class `boschloo_control`.
#' @examples
#' boschloo_control()
#' boschloo_control(grid_points = 2001L)
#' @family Boschloo
#' @export
boschloo_control <- function(tie_tolerance = 1e-7, grid_points = 1001L, nuisance_tol = 1e-9,
                             polish_tol = 1e-12, max_levels = 60L,
                             ambiguity_band = c(1e-12, 1e-6), cache = TRUE) {
  p <- character()
  num1 <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!num1(tie_tolerance) || tie_tolerance < 0 || tie_tolerance > 1e-4) {
    p <- c(p, "`tie_tolerance` must be a number in [0, 1e-4]")
  }
  if (!num1(grid_points) || grid_points != round(grid_points) || grid_points < 3 ||
      grid_points > 1e5) {
    p <- c(p, "`grid_points` must be a whole number between 3 and 100000")
  }
  if (!num1(nuisance_tol) || nuisance_tol <= 0 || nuisance_tol > 1e-4) {
    p <- c(p, "`nuisance_tol` must be a number in (0, 1e-4]")
  }
  if (!num1(polish_tol) || polish_tol <= 0 || polish_tol > 1e-6) {
    p <- c(p, "`polish_tol` must be a number in (0, 1e-6]")
  }
  if (!num1(max_levels) || max_levels != round(max_levels) || max_levels < 1 || max_levels > 200) {
    p <- c(p, "`max_levels` must be a whole number between 1 and 200")
  }
  band_ok <- is.numeric(ambiguity_band) && length(ambiguity_band) == 2L && !anyNA(ambiguity_band) &&
    ambiguity_band[1] >= 0 && ambiguity_band[1] <= 1e-12 &&
    ambiguity_band[2] >= 1e-10 && ambiguity_band[2] <= 1e-3
  if (!band_ok) {
    p <- c(p, "`ambiguity_band` must be c(lower, upper) with 0 <= lower <= 1e-12 and 1e-10 <= upper <= 1e-3")
  }
  # Unflagged gaps are classified by the band edges alone; the tolerance
  # decides only inside the band, where every result is checked.
  if (band_ok && num1(tie_tolerance) &&
      (tie_tolerance < ambiguity_band[1] || tie_tolerance > ambiguity_band[2])) {
    p <- c(p, sprintf(paste("`tie_tolerance` (%s) must lie within `ambiguity_band` [%s, %s]; outside it,",
                            "tables would be split or merged without an ambiguity flag"),
                      format(tie_tolerance), format(ambiguity_band[1]), format(ambiguity_band[2])))
  }
  if (!.ig_is_flag(cache)) p <- c(p, "`cache` must be TRUE or FALSE")
  if (length(p)) .ig_abort_invalid(p, "invalid Boschloo control settings")
  structure(list(tie_tolerance = as.double(tie_tolerance),
                 grid_points = as.integer(grid_points),
                 nuisance_tol = as.double(nuisance_tol),
                 polish_tol = as.double(polish_tol),
                 max_levels = as.integer(max_levels),
                 ambiguity_band = as.double(ambiguity_band),
                 cache = cache),
            class = "boschloo_control")
}

#' @export
print.boschloo_control <- function(x, ...) {
  cat("<boschloo_control>\n",
      sprintf("  tie tolerance (relative): %s\n", format(x$tie_tolerance)),
      sprintf("  nuisance grid points:     %d (plus branch-and-bound bisection)\n", x$grid_points),
      sprintf("  supremum bound tolerance: %s\n", format(x$nuisance_tol)),
      sprintf("  Brent polish tolerance:   %s\n", format(x$polish_tol)),
      sprintf("  ambiguity band:           (%s, %s]\n", format(x$ambiguity_band[1]),
              format(x$ambiguity_band[2])),
      sprintf("  cache:                    %s\n", x$cache), sep = "")
  invisible(x)
}

.bl_check_control <- function(control) {
  if (!inherits(control, "boschloo_control")) {
    .ig_abort("invalid_input", "`control` must be created with boschloo_control()")
  }
  control
}

# Settings as a plain list for specifications and results (no cache flag:
# caching never changes a number).
.bl_control_record <- function(control) {
  list(tie_tolerance = control$tie_tolerance, grid_points = control$grid_points,
       nuisance_tol = control$nuisance_tol, polish_tol = control$polish_tol,
       max_levels = control$max_levels, ambiguity_band = I(control$ambiguity_band))
}

# Rebuild a control object from a recorded settings list (results and their
# JSON round trips store whole numbers and vectors in plain form).
.bl_control_from_record <- function(rec) {
  if (!is.list(rec)) .ig_abort("invalid_input", "the object carries no Boschloo settings")
  boschloo_control(tie_tolerance = as.double(rec$tie_tolerance),
                   grid_points = as.integer(rec$grid_points),
                   nuisance_tol = as.double(rec$nuisance_tol),
                   polish_tol = as.double(rec$polish_tol),
                   max_levels = as.integer(rec$max_levels),
                   ambiguity_band = as.double(unlist(rec$ambiguity_band)))
}

.bl_control_key <- function(control) {
  paste(.ig_format_double(c(control$tie_tolerance, control$nuisance_tol, control$polish_tol,
                            control$ambiguity_band)),
        control$grid_points, control$max_levels, sep = "|", collapse = "|")
}

.bl_variant <- function(variant) {
  if (missing(variant) || is.null(variant)) {
    .ig_abort_invalid(sprintf("`variant` must be given explicitly; one of: %s",
                              paste(.bl_variants, collapse = ", ")))
  }
  if (!.ig_is_string(variant) || !(variant %in% .bl_variants)) {
    .ig_abort_invalid(sprintf("unknown Boschloo variant '%s'; one of: %s",
                              paste(format(variant), collapse = ""),
                              paste(.bl_variants, collapse = ", ")))
  }
  .bl_variant_info[[variant]]
}

# Shared input checks ---------------------------------------------------------------

.bl_is_count <- function(x, lower = 0) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x == round(x) && x >= lower
}

.bl_check_design <- function(n1, n2) {
  p <- character()
  if (!.bl_is_count(n1, 1)) p <- c(p, "`n1` must be a whole number >= 1")
  if (!.bl_is_count(n2, 1)) p <- c(p, "`n2` must be a whole number >= 1")
  if (length(p)) .ig_abort_invalid(p, "invalid Boschloo design")
  if ((n1 + 1) * (n2 + 1) > .bl_max_tables) {
    .ig_abort("unsupported", sprintf(
      "the design n1 = %d, n2 = %d has %.0f tables; at most %.0f tables are enumerated",
      as.integer(n1), as.integer(n2), (n1 + 1) * (n2 + 1), .bl_max_tables))
  }
  invisible(TRUE)
}

.bl_check_alpha <- function(alpha, name = "alpha") {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    .ig_abort_invalid(sprintf("`%s` must be a number in (0, 1)", name))
  }
  as.double(alpha)
}

.bl_check_prob <- function(p, name) {
  if (!is.numeric(p) || length(p) != 1L || !is.finite(p) || p < 0 || p > 1) {
    .ig_abort_invalid(sprintf("`%s` must be a probability in [0, 1]", name))
  }
  as.double(p)
}

# In-memory cache ---------------------------------------------------------------------
#
# Orderings and rejection regions are pure functions of their key (design,
# side or variant, alpha and every numerical setting), so caching never
# changes a result. The cache is bounded; the oldest entries are dropped.

.bl_cache <- new.env(parent = emptyenv())
.bl_cache_max <- 256L

.bl_cache_get <- function(key, control) {
  if (!isTRUE(control$cache)) return(NULL)
  entry <- .bl_cache[[key]]
  if (is.null(entry)) NULL else entry$value
}

.bl_cache_counter <- new.env(parent = emptyenv())
.bl_cache_counter$n <- 0

.bl_cache_set <- function(key, value, control) {
  if (!isTRUE(control$cache)) return(value)
  keys <- ls(.bl_cache, all.names = TRUE)
  if (length(keys) >= .bl_cache_max) {
    stamps <- vapply(keys, function(k) .bl_cache[[k]]$stamp, double(1))
    rm(list = keys[order(stamps)][seq_len(length(keys) - .bl_cache_max + 1L)], envir = .bl_cache)
  }
  .bl_cache_counter$n <- .bl_cache_counter$n + 1
  assign(key, list(value = value, stamp = .bl_cache_counter$n), envir = .bl_cache)
  value
}

.bl_cache_clear <- function() {
  rm(list = ls(.bl_cache, all.names = TRUE), envir = .bl_cache)
  invisible(NULL)
}

# Provenance shared by all Boschloo results. Package and R versions are not
# part of the result body (they would make result hashes runtime-dependent);
# analysis_envelope() records them as envelope provenance.
.bl_provenance <- function() {
  list(engine = "inferguidR Boschloo engine",
       ordering = "stats::dhyper / stats::phyper in log space; exact integer fallback for n1 + n2 <= 56",
       nuisance = "Bernstein form with stats::dbinom; grid + branch-and-bound + stats::optimize polish",
       power = "double binomial sum (stats::dbinom) over the enumerated rejection region")
}

#' Boschloo test variants and their evidence status
#'
#' Lists the four Boschloo variants with their identifiers, hypotheses,
#' ordering statistics and two-sided constructions, the matching calls of the
#' comparator packages, and the evidence status shared by all variants:
#' the independent validation status, the planning-compatibility class and
#' the provisional classification of the related G*Power procedures.
#'
#' Group 1 is `(x1, n1)` and group 2 is `(x2, n2)` throughout; `x1` and `x2`
#' count events. The null hypothesis is always `p1 = p2`. Mid-p values,
#' randomisation, Berger-Boos confidence-set restrictions and non-zero null
#' margins are not part of any variant.
#'
#' @return A data frame with one row per variant: `variant` (the value of the
#'   `variant` argument), `variant_id`, `advisor_variant_id`, `alternative`,
#'   `label`, `ordering_statistic`, `construction`, `exact_call`,
#'   `exact2x2_call`, `validation_status`, `planning_compatibility`,
#'   `gpower_relation` and `gpower_classification_status`.
#' @examples
#' boschloo_variants()[, c("variant", "variant_id", "alternative")]
#' @family Boschloo
#' @export
boschloo_variants <- function() {
  rows <- lapply(.bl_variant_info, function(v) {
    data.frame(variant = v$variant, variant_id = v$variant_id,
               advisor_variant_id = v$advisor_variant_id, alternative = v$alternative,
               label = v$label, ordering_statistic = v$ordering_statistic,
               construction = v$construction, exact_call = v$exact_call,
               exact2x2_call = v$exact2x2_call,
               validation_status = .bl_evidence$validation_status,
               planning_compatibility = .bl_evidence$planning_compatibility,
               gpower_relation = .bl_evidence$gpower_relation,
               gpower_classification_status = .bl_evidence$gpower_classification_status,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "gpower_note") <- .bl_evidence$gpower_note
  out
}
