# Effect-size drawers: t tests for means ---------------------------------------------------
#
# G*Power's "Determine =>" drawers of the three t tests for means (group G1).
# They are reached through effect_size() (R/core-api.R), which looks up
# `.ig_drawer_<method_id>()`. Every drawer returns the signed value as the
# method defines it and its magnitude.
#
# Sign conventions displayed by G*Power 3.1.9.7 (emulator captures):
#
# * two independent means, both input modes: the magnitude |d| (mean group 1 >
#   mean group 2 gives a positive d: SUP-G1-DRAW-MEAN2-UD01, -UD02; mean
#   group 2 > mean group 1 also gives a positive d: GPW-UT-DRAW-UTL-U047/U048);
# * one sample: the magnitude |d| (Mean H1 < Mean H0 gives 0.625:
#   SUP-G1-DRAW-MEAN1-UD01);
# * matched pairs "From group parameters": the magnitude |dz|
#   (SUP-G1-DRAW-PAIRED-UD01); "From differences": the signed
#   dz = mean of difference / SD of difference (GPW-UT-DRAW-UTL-U045-UE01
#   shows -50).
#
# The `value` returned here is signed as the method's effect parameter is
# defined (manual 19.1, 20.1, 21.1); `magnitude` is what G*Power's drawer
# displays where it drops the sign.

# Shared input checks: `required` names, finite numbers, positive SDs.
.ig_drawer_inputs_g1 <- function(input_mode, inputs, modes, required, positive = character(),
                                 open_unit = character()) {
  if (!.ig_is_string(input_mode) || !(input_mode %in% modes)) {
    .ig_abort_invalid(sprintf("`input_mode` must be one of: %s", paste(modes, collapse = ", ")),
                      "invalid effect-size input")
  }
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
  for (nm in intersect(positive, required)) {
    if (!is.na(vals[[nm]]) && vals[[nm]] <= 0) problems <- c(problems, sprintf("'%s' must be positive", nm))
  }
  for (nm in intersect(open_unit, required)) {
    if (!is.na(vals[[nm]]) && abs(vals[[nm]]) >= 1) problems <- c(problems, sprintf("'%s' must lie in (-1, 1)", nm))
  }
  list(values = vals, problems = problems)
}

# Means: Difference between two independent means (two groups) -------------------------------
#
# Manual 21.1, Fig. 28; drawer_T_2Groups on the Mac, RT_DIALOG 150 on
# Windows. Two input modes, plus the r -> d conversion the manual uses by
# hand in 16.3.

.ig_drawer_t_means_two_independent <- function(input_mode, inputs) {
  modes <- c("common_sd", "two_sds", "point_biserial_r")
  required <- if (.ig_is_string(input_mode)) switch(input_mode,
                                                    common_sd = c("mean_1", "mean_2", "sd"),
                                                    two_sds = c("mean_1", "mean_2", "sd_1", "sd_2"),
                                                    point_biserial_r = c("r", "n1", "n2"),
                                                    character()) else character()
  chk <- .ig_drawer_inputs_g1(input_mode, inputs, modes, required, positive = c("sd", "sd_1", "sd_2"))
  vals <- chk$values
  problems <- chk$problems
  if (input_mode == "point_biserial_r") {
    if (!is.na(vals$r) && abs(vals$r) >= 1) problems <- c(problems, "'r' must lie in (-1, 1)")
    for (nm in c("n1", "n2")) {
      v <- vals[[nm]]
      if (!is.na(v) && (v < 1 || v != round(v))) problems <- c(problems, sprintf("'%s' must be a whole number >= 1", nm))
    }
  }
  if (length(problems)) .ig_abort_invalid(problems, "invalid effect-size input")

  warnings <- character()
  if (input_mode == "common_sd") {
    d <- (vals$mean_1 - vals$mean_2) / vals$sd
    formula <- "d = (mean_1 - mean_2) / sd"
    source <- paste("manual 21.1; G*Power drawer input mode 'n1 != n2' (common SD within each group);",
                    "G*Power displays |d|")
  } else if (input_mode == "two_sds") {
    sigma <- sqrt((vals$sd_1^2 + vals$sd_2^2) / 2)
    d <- (vals$mean_1 - vals$mean_2) / sigma
    formula <- "d = (mean_1 - mean_2) / sqrt((sd_1^2 + sd_2^2) / 2)"
    source <- "manual 21.1 (Cohen 1969, p. 42); G*Power drawer input mode 'n1 = n2'; G*Power displays |d|"
    warnings <- paste("the mean SD sqrt((sd_1^2 + sd_2^2) / 2) is meant for equal group sizes; with",
                      "clearly unequal groups it can give misleading power (Cohen 1969, p. 42; manual 21.1)")
  } else {
    n <- vals$n1 + vals$n2
    d <- vals$r / sqrt(1 - vals$r^2) * n / sqrt(vals$n1 * vals$n2)
    formula <- "d = r / sqrt(1 - r^2) * (n1 + n2) / sqrt(n1 * n2)"
    source <- "manual 16.3 (hand conversion of a point-biserial r); not a G*Power drawer mode"
  }
  list(effect = "d", value = d, magnitude = abs(d), formula = formula, source = source,
       inputs = vals, warnings = warnings)
}

# Means: Difference from constant (one sample case) -----------------------------------------
#
# Manual 20.1; drawer_T_OneMean (Mac), RT_DIALOG 169 (Windows), one input
# mode: Mean H0, Mean H1, SD sigma.

.ig_drawer_t_means_one_sample <- function(input_mode, inputs) {
  required <- c("mean_h0", "mean_h1", "sd")
  chk <- .ig_drawer_inputs_g1(input_mode, inputs, "from_means", required, positive = "sd")
  if (length(chk$problems)) .ig_abort_invalid(chk$problems, "invalid effect-size input")
  vals <- chk$values
  d <- (vals$mean_h1 - vals$mean_h0) / vals$sd
  list(effect = "d", value = d, magnitude = abs(d), formula = "d = (mean_h1 - mean_h0) / sd",
       source = paste("manual 20.1 (d = (mu - mu0) / sigma); G*Power drawer 'from means and SD' (Mean H0, Mean H1,",
                      "SD sigma); G*Power displays |d|"),
       inputs = vals, warnings = character())
}

# Means: Difference between two dependent means (matched pairs) ------------------------------
#
# Manual 19.1; drawer_T_PairedMeans (Mac), RT_DIALOG 170 (Windows), two input
# modes: from the mean and SD of the differences, or from the group means,
# SDs and the correlation between the paired measurements,
#   sigma_z = sqrt(sd_1^2 + sd_2^2 - 2 rho sd_1 sd_2).

.ig_drawer_t_means_matched_pairs <- function(input_mode, inputs) {
  modes <- c("from_differences", "from_group_parameters")
  required <- if (.ig_is_string(input_mode)) switch(input_mode,
                                                    from_differences = c("mean_diff", "sd_diff"),
                                                    from_group_parameters = c("mean_1", "mean_2", "sd_1", "sd_2", "rho"),
                                                    character()) else character()
  chk <- .ig_drawer_inputs_g1(input_mode, inputs, modes, required, positive = c("sd_diff", "sd_1", "sd_2"),
                              open_unit = "rho")
  vals <- chk$values
  problems <- chk$problems
  if (input_mode == "from_group_parameters" && !length(problems)) {
    var_z <- vals$sd_1^2 + vals$sd_2^2 - 2 * vals$rho * vals$sd_1 * vals$sd_2
    if (!(var_z > 0)) problems <- c(problems, "the SD of the differences sqrt(sd_1^2 + sd_2^2 - 2 rho sd_1 sd_2) is 0")
  }
  if (length(problems)) .ig_abort_invalid(problems, "invalid effect-size input")
  if (input_mode == "from_differences") {
    dz <- vals$mean_diff / vals$sd_diff
    formula <- "dz = mean_diff / sd_diff"
    source <- "manual 19.1; G*Power drawer input mode 'From differences' (signed, as G*Power displays it)"
  } else {
    sd_z <- sqrt(vals$sd_1^2 + vals$sd_2^2 - 2 * vals$rho * vals$sd_1 * vals$sd_2)
    dz <- (vals$mean_1 - vals$mean_2) / sd_z
    formula <- "dz = (mean_1 - mean_2) / sqrt(sd_1^2 + sd_2^2 - 2 rho sd_1 sd_2)"
    source <- paste("manual 19.1 (dz = |mu_x - mu_y| / sigma_z); G*Power drawer input mode 'From group",
                    "parameters'; G*Power displays |dz|")
  }
  list(effect = "dz", value = dz, magnitude = abs(dz), formula = formula, source = source,
       inputs = vals, warnings = character())
}
