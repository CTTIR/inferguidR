# Effect-size drawers: fixed-effects ANOVA, ANCOVA and Hotelling T^2 (group G2) -------------------
#
# Reached through effect_size() (R/core-api.R), which looks up `.ig_drawer_<method_id>()`.
#
# ANOVA and ANCOVA (effect f; manual 10.1, 11.1 and the drawer captures
# GPW-UT-DRAW-UTL-U013..U018):
#
#   "means"        (one-way only) group means mu_i, group sizes n_i and the common SD sigma:
#                  f = sigma_m / sigma, sigma_m^2 = sum n_i (mu_i - mu_bar)^2 / N,
#                  mu_bar = sum n_i mu_i / N, N = sum n_i (manual 10.1; DOC-057, DOC-079);
#   "variances"    variance explained by the effect s_m^2 and the error (within-group) variance
#                  s^2: f = sqrt(s_m^2 / s^2), partial eta^2 = s_m^2 / (s_m^2 + s^2)
#                  (manual 10.1 and 11.1; DOC-058, DOC-059, DOC-072);
#   "eta_squared"  partial eta^2 directly: f = sqrt(eta^2 / (1 - eta^2)) (drawer mode "Direct").
#
# Hotelling T^2 (effect Delta, no manual chapter): the Mahalanobis distance
# Delta = sqrt(d' Sigma^-1 d) of the mean difference d (mu - mu0 for one group, mu1 - mu2 for two
# groups) with the common covariance matrix Sigma, given either as a covariance matrix
# ("covariance") or as standard deviations and a correlation matrix ("sd_correlation"). G*Power's
# drawer (RT_DIALOG 172) also offers correlation patterns and a "Multiply all means by" factor; no
# capture of that drawer exists, so only the definition of Delta is implemented (method notes).
#
# Derived quantities (partial eta^2, total N) are returned in `outputs`; effect_size() keeps the
# value of the effect only.

.ig_g2_drawer_problems <- function(inputs, required, optional = character()) {
  unknown <- setdiff(names(inputs), c(required, optional))
  missing <- setdiff(required, names(inputs))
  c(if (length(unknown)) sprintf("input '%s' is not used here (inputs: %s)", unknown,
                                 paste(c(required, optional), collapse = ", ")),
    if (length(missing)) sprintf("missing input '%s'", missing))
}

.ig_g2_num <- function(v, nm, len = 1L) {
  if (!is.numeric(v) || anyNA(v) || any(!is.finite(v)) || (!is.na(len) && length(v) != len) || !length(v)) {
    return(sprintf("'%s' must be %s", nm, if (identical(len, 1L)) "a single finite number" else
      "a vector of finite numbers"))
  }
  NULL
}

# f from variances or partial eta^2 (shared by the three ANOVA/ANCOVA drawers).
.ig_g2_drawer_f_variance <- function(input_mode, inputs, method_text) {
  if (input_mode == "variances") {
    problems <- .ig_g2_drawer_problems(inputs, c("var_effect", "var_error"))
    if (!length(problems)) {
      problems <- c(.ig_g2_num(inputs$var_effect, "var_effect"), .ig_g2_num(inputs$var_error, "var_error"))
    }
    if (!length(problems)) {
      if (inputs$var_effect < 0) problems <- c(problems, "'var_effect' must be >= 0")
      if (inputs$var_error <= 0) problems <- c(problems, "'var_error' must be positive")
    }
    if (length(problems)) .ig_abort_invalid(problems, "invalid effect-size input")
    ve <- as.double(inputs$var_effect)
    vw <- as.double(inputs$var_error)
    return(list(effect = "f", value = sqrt(ve / vw),
                formula = "f = sqrt(var_effect / var_error); partial eta^2 = var_effect / (var_effect + var_error)",
                source = paste(method_text, "drawer mode 'From variances' (manual 10.1, 11.1)"),
                inputs = list(var_effect = ve, var_error = vw),
                outputs = list(partial_eta_squared = ve / (ve + vw)), warnings = character()))
  }
  problems <- .ig_g2_drawer_problems(inputs, "partial_eta_squared")
  if (!length(problems)) problems <- .ig_g2_num(inputs$partial_eta_squared, "partial_eta_squared")
  if (!length(problems) && (inputs$partial_eta_squared < 0 || inputs$partial_eta_squared >= 1)) {
    problems <- "'partial_eta_squared' must lie in [0, 1)"
  }
  if (length(problems)) .ig_abort_invalid(problems, "invalid effect-size input")
  e2 <- as.double(inputs$partial_eta_squared)
  list(effect = "f", value = sqrt(e2 / (1 - e2)), formula = "f = sqrt(eta^2 / (1 - eta^2))",
       source = paste(method_text, "drawer mode 'Direct' (partial eta^2)"),
       inputs = list(partial_eta_squared = e2), outputs = list(partial_eta_squared = e2),
       warnings = character())
}

.ig_g2_check_mode <- function(input_mode, modes) {
  if (!.ig_is_string(input_mode) || !(input_mode %in% modes)) {
    .ig_abort_invalid(sprintf("`input_mode` must be one of: %s", paste(modes, collapse = ", ")),
                      "invalid effect-size input")
  }
}

.ig_drawer_f_anova_oneway <- function(input_mode, inputs) {
  .ig_g2_check_mode(input_mode, c("means", "variances", "eta_squared"))
  if (input_mode != "means") return(.ig_g2_drawer_f_variance(input_mode, inputs, "one-way ANOVA"))
  problems <- .ig_g2_drawer_problems(inputs, c("means", "sd"), c("sizes", "n_per_group"))
  has_sizes <- !is.null(inputs$sizes)
  has_equal <- !is.null(inputs$n_per_group)
  if (has_sizes == has_equal) problems <- c(problems, "give either 'sizes' (one per group) or 'n_per_group' (equal n)")
  if (!length(problems)) {
    problems <- c(.ig_g2_num(inputs$means, "means", NA), .ig_g2_num(inputs$sd, "sd"),
                  if (has_sizes) .ig_g2_num(inputs$sizes, "sizes", length(inputs$means)) else
                    .ig_g2_num(inputs$n_per_group, "n_per_group"))
  }
  if (!length(problems)) {
    if (length(inputs$means) < 2) problems <- c(problems, "'means' needs at least two groups")
    if (inputs$sd <= 0) problems <- c(problems, "'sd' must be positive")
    sz <- if (has_sizes) inputs$sizes else inputs$n_per_group
    if (any(sz < 1 | sz != round(sz))) problems <- c(problems, "group sizes must be whole numbers >= 1")
  }
  if (length(problems)) .ig_abort_invalid(problems, "invalid effect-size input")
  mu <- as.double(inputs$means)
  n_i <- if (has_sizes) as.double(inputs$sizes) else rep(as.double(inputs$n_per_group), length(mu))
  n <- sum(n_i)
  mu_bar <- sum(n_i * mu) / n
  sigma_m <- sqrt(sum(n_i * (mu - mu_bar)^2) / n)
  f <- sigma_m / as.double(inputs$sd)
  list(effect = "f", value = f,
       formula = "f = sigma_m / sd, sigma_m^2 = sum n_i (mean_i - mean_w)^2 / N, mean_w = sum n_i mean_i / N",
       source = "one-way ANOVA drawer mode 'Effect size from means' (manual 10.1)",
       inputs = list(means = mu, sizes = n_i, sd = as.double(inputs$sd)),
       outputs = list(n_total = n, sigma_m = sigma_m, partial_eta_squared = f^2 / (1 + f^2)),
       warnings = character())
}

.ig_drawer_f_anova_special_main_interactions <- function(input_mode, inputs) {
  .ig_g2_check_mode(input_mode, c("variances", "eta_squared"))
  .ig_g2_drawer_f_variance(input_mode, inputs, "ANOVA special effects")
}

.ig_drawer_f_ancova_main_interactions <- function(input_mode, inputs) {
  .ig_g2_check_mode(input_mode, c("variances", "eta_squared"))
  .ig_g2_drawer_f_variance(input_mode, inputs, "ANCOVA")
}

# Common covariance matrix from the drawer inputs, checked for symmetry and positive definiteness.
.ig_g2_covariance <- function(input_mode, inputs, p) {
  problems <- character()
  if (input_mode == "covariance") {
    S <- inputs$covariance
    if (!is.matrix(S) || !is.numeric(S) || any(!is.finite(S)) || nrow(S) != p || ncol(S) != p) {
      return(list(problems = sprintf("'covariance' must be a finite numeric %d x %d matrix", p, p)))
    }
  } else {
    sds <- inputs$sds
    R <- inputs$correlation
    if (!is.numeric(sds) || length(sds) != p || any(!is.finite(sds)) || any(sds <= 0)) {
      problems <- c(problems, sprintf("'sds' must be %d positive finite numbers", p))
    }
    if (!is.matrix(R) || !is.numeric(R) || any(!is.finite(R)) || nrow(R) != p || ncol(R) != p) {
      problems <- c(problems, sprintf("'correlation' must be a finite numeric %d x %d matrix", p, p))
    } else {
      if (any(abs(diag(R) - 1) > 1e-12)) problems <- c(problems, "'correlation' must have a unit diagonal")
      if (any(abs(R) > 1)) problems <- c(problems, "correlations must lie in [-1, 1]")
    }
    if (length(problems)) return(list(problems = problems))
    S <- outer(sds, sds) * R
  }
  S <- unname(S)
  if (max(abs(S - t(S))) > 1e-12 * max(1, max(abs(S)))) {
    return(list(problems = "the covariance (correlation) matrix must be symmetric"))
  }
  S <- (S + t(S)) / 2
  ch <- tryCatch(chol(S), error = function(e) NULL)
  if (is.null(ch) || min(diag(ch)) <= sqrt(.Machine$double.eps) * sqrt(max(diag(S)))) {
    return(list(problems = "the covariance (correlation) matrix must be positive definite"))
  }
  list(problems = character(), S = S, chol = ch)
}

.ig_g2_drawer_hotelling <- function(input_mode, inputs, groups) {
  .ig_g2_check_mode(input_mode, c("covariance", "sd_correlation"))
  mean_names <- if (groups == 1) c("means", "null_means") else c("means_1", "means_2")
  cov_names <- if (input_mode == "covariance") "covariance" else c("sds", "correlation")
  problems <- .ig_g2_drawer_problems(inputs, c(mean_names, cov_names))
  if (!length(problems)) {
    problems <- c(.ig_g2_num(inputs[[mean_names[1]]], mean_names[1], NA),
                  .ig_g2_num(inputs[[mean_names[2]]], mean_names[2], NA))
  }
  if (!length(problems) && length(inputs[[mean_names[1]]]) != length(inputs[[mean_names[2]]])) {
    problems <- sprintf("'%s' and '%s' must have the same length", mean_names[1], mean_names[2])
  }
  if (length(problems)) .ig_abort_invalid(problems, "invalid effect-size input")
  d <- as.double(inputs[[mean_names[1]]]) - as.double(inputs[[mean_names[2]]])
  p <- length(d)
  cv <- .ig_g2_covariance(input_mode, inputs, p)
  if (length(cv$problems)) .ig_abort_invalid(cv$problems, "invalid effect-size input")
  z <- backsolve(cv$chol, d, transpose = TRUE)
  delta <- sqrt(sum(z^2))
  list(effect = "delta", value = delta,
       formula = sprintf("Delta = sqrt(d' Sigma^-1 d), d = %s - %s", mean_names[1], mean_names[2]),
       source = paste("Mahalanobis distance of the mean", if (groups == 1) "vector from its null value" else
         "vectors of two groups", "with the common covariance matrix (Hotelling T^2; no G*Power drawer capture)"),
       inputs = inputs, outputs = list(response_variables = p, delta_squared = sum(z^2)),
       warnings = character())
}

.ig_drawer_f_hotelling_t2_one_group <- function(input_mode, inputs) {
  .ig_g2_drawer_hotelling(input_mode, inputs, 1)
}

.ig_drawer_f_hotelling_t2_two_groups <- function(input_mode, inputs) {
  .ig_g2_drawer_hotelling(input_mode, inputs, 2)
}
