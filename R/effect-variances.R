# Effect-size drawers: variance tests and goodness of fit (group G11) --------------------
#
# Reached through effect_size() (R/core-api.R), which looks up
# `.ig_drawer_<method_id>()`.
#
# * f_variances_two_groups, chi2_variance_one_sample: G*Power's drawer "from two
#   variances" (Windows RT_DIALOG 194; Mac drawer_F_VarEquality /
#   drawer_Chi_Variance) computes Ratio var1/var0 from the fields "Variance var0"
#   and "Variance var1" (live labels; the manual writes V0 and V1). Documented
#   examples DOC-096 and DOC-136 (V0 = 1.5, V1 = 1 gives 0.6666667) and the
#   drawer captures GPW-UT-DRAW-UTL-U043 / U074 fix the orientation:
#   ratio = var1 / var0.
# * chi2_goodness_of_fit_contingency: G*Power's drawer is a table of cell
#   probabilities p(H0), p(H1) (RT_DIALOG 157/212, MFCGridCtrl; not automatable,
#   campaign cases GPW-UT-DRAW-UTL-U075 skipped). Cohen's w =
#   sqrt(sum((p1 - p0)^2 / p0)) (Cohen 1988, eq. 7.2.1). A contingency-table
#   input mode (joint H1 probabilities; p0 = product of their margins) is an
#   inferguidR convenience, not a G*Power drawer mode.

.ig_g11_drawer_inputs <- function(input_mode, inputs, modes, required) {
  if (!.ig_is_string(input_mode) || !(input_mode %in% modes)) {
    .ig_abort_invalid(sprintf("`input_mode` must be one of: %s", paste(modes, collapse = ", ")),
                      "invalid effect-size input")
  }
  req <- required[[input_mode]]
  problems <- character()
  unknown <- setdiff(names(inputs), req)
  if (length(unknown)) {
    problems <- c(problems, sprintf("input '%s' is not used in input mode '%s' (inputs: %s)", unknown,
                                    input_mode, paste(req, collapse = ", ")))
  }
  missing <- setdiff(req, names(inputs))
  if (length(missing)) problems <- c(problems, sprintf("missing input '%s'", missing))
  problems
}

.ig_g11_variance_ratio_drawer <- function(input_mode, inputs) {
  problems <- .ig_g11_drawer_inputs(input_mode, inputs, "from_variances",
                                    list(from_variances = c("var0", "var1")))
  vals <- list()
  for (nm in intersect(c("var0", "var1"), names(inputs))) {
    v <- inputs[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v)) {
      problems <- c(problems, sprintf("'%s' must be a single finite number", nm))
    } else if (v < 0 || (nm == "var0" && v == 0)) {
      problems <- c(problems, sprintf("'%s' must be %s", nm, if (nm == "var0") "positive" else "non-negative"))
    } else {
      vals[[nm]] <- as.double(v)
    }
  }
  if (length(problems)) .ig_abort_invalid(problems, "invalid effect-size input")
  list(effect = "ratio", value = vals$var1 / vals$var0,
       formula = "ratio = var1 / var0",
       source = paste("G*Power drawer 'from two variances' (fields 'Variance var0' and 'Variance var1');",
                      "documented examples DOC-096 / DOC-136"),
       inputs = vals[c("var0", "var1")], warnings = character())
}

.ig_drawer_f_variances_two_groups <- function(input_mode, inputs) {
  out <- .ig_g11_variance_ratio_drawer(input_mode, inputs)
  out$source <- paste(out$source, "(manual 15.1: var1 is the variance of group 1, var0 that of group 2)")
  out
}

.ig_drawer_chi2_variance_one_sample <- function(input_mode, inputs) {
  out <- .ig_g11_variance_ratio_drawer(input_mode, inputs)
  out$source <- paste(out$source, "(manual 26.1: var1 is the population variance under H1, var0 the constant)")
  out
}

# Cohen's w from H0 and H1 cell probabilities. Each vector must be a
# probability distribution (sum 1 within 1e-9 relative); p0 must be positive.
.ig_g11_cohen_w <- function(p0, p1) {
  sqrt(sum((p1 - p0)^2 / p0))
}

.ig_drawer_chi2_goodness_of_fit_contingency <- function(input_mode, inputs) {
  modes <- c("cell_probabilities", "contingency_table")
  problems <- .ig_g11_drawer_inputs(input_mode, inputs, modes,
                                    list(cell_probabilities = c("p0", "p1"), contingency_table = "p1"))
  check_prob <- function(nm, v, positive) {
    if (!is.numeric(v) || !length(v) || any(!is.finite(v))) {
      return(sprintf("'%s' must be finite numbers", nm))
    }
    out <- character()
    if (positive && any(v <= 0)) out <- c(out, sprintf("every cell of '%s' must be positive", nm))
    if (!positive && any(v < 0)) out <- c(out, sprintf("no cell of '%s' may be negative", nm))
    if (any(v > 1)) out <- c(out, sprintf("no cell of '%s' may exceed 1", nm))
    if (abs(sum(v) - 1) > 1e-9) {
      out <- c(out, sprintf("the cells of '%s' must sum to 1 (sum %s)", nm, format(sum(v), digits = 15)))
    }
    out
  }
  if (!length(problems) && identical(input_mode, "cell_probabilities")) {
    p0 <- inputs$p0
    p1 <- inputs$p1
    problems <- c(problems, check_prob("p0", p0, TRUE), check_prob("p1", p1, FALSE))
    if (is.numeric(p0) && is.numeric(p1) && length(p0) != length(p1)) {
      problems <- c(problems, "'p0' and 'p1' must have the same number of cells")
    }
    if (is.numeric(p0) && length(p0) < 2L) problems <- c(problems, "at least 2 cells are needed")
    if (length(problems)) .ig_abort_invalid(problems, "invalid effect-size input")
    w <- .ig_g11_cohen_w(as.double(p0), as.double(p1))
    return(list(effect = "w", value = w,
                formula = "w = sqrt(sum((p1 - p0)^2 / p0))",
                source = paste("Cohen (1988) eq. 7.2.1; G*Power drawer 'Goodness-of-fit tests: Contingency",
                               "tables' (cell probabilities p(H0), p(H1)); Df = cells - 1"),
                inputs = list(p0 = as.double(p0), p1 = as.double(p1)),
                warnings = character()))
  }
  if (!length(problems)) {
    p1 <- inputs$p1
    if (!is.matrix(p1) || !is.numeric(p1) || nrow(p1) < 2L || ncol(p1) < 2L) {
      problems <- c(problems, "'p1' must be a numeric matrix with at least 2 rows and 2 columns")
    } else {
      problems <- c(problems, check_prob("p1", as.vector(p1), FALSE))
      if (!length(problems) && (any(rowSums(p1) <= 0) || any(colSums(p1) <= 0))) {
        problems <- c(problems, "every row and column of 'p1' needs a positive margin")
      }
    }
  }
  if (length(problems)) .ig_abort_invalid(problems, "invalid effect-size input")
  p0 <- outer(rowSums(p1), colSums(p1))
  w <- .ig_g11_cohen_w(as.vector(p0), as.vector(p1))
  list(effect = "w", value = w,
       formula = "w = sqrt(sum((p1_ij - p_i. p_.j)^2 / (p_i. p_.j))) with the margins of p1",
       source = paste("Cohen (1988) ch. 7 (independence test of a contingency table); inferguidR convenience,",
                      "not a G*Power drawer mode; Df = (rows - 1)(columns - 1)"),
       inputs = list(p1 = p1),
       warnings = sprintf("use Df = %d (rows - 1)(columns - 1) with this w", (nrow(p1) - 1L) * (ncol(p1) - 1L)))
}
