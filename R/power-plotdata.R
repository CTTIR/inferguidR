#' Plot-ready power distributions and conditional decision probabilities
#'
#' All decision probabilities are taken from the computed result. Display
#' ordinates are densities or probability masses, never estimates of power.
#' Unsupported distribution descriptors return an explicit unavailable panel.
#' The decision table conditions on each hypothesis separately: its columns
#' must not be interpreted as joint probabilities or predictive values.
#' The Boschloo lattice uses the solved group sizes and recorded numerical
#' settings. Its H0 panel uses a size-maximising common event probability;
#' its H1 panel uses the specified alternative. Crosses mark rejection tables,
#' not observed data. Larger lattices return an explicit unavailable status.
#'
#' @param result An unmodified `power_result`.
#' @param panel One of `"distribution"`, `"decision_table"`, `"mosaic"`, or `"lattice"` (Boschloo table probabilities).
#' @param points Sampling resolution (51 to 2001). Each continuous law uses
#'   this many quantile points and equally spaced points, plus critical values.
#'   The lattice enumerates all tables (limit 5000), independently of this resolution.
#' @return A list with status, source hashes, plain data frames, and notes.
#' @export
#' @examples
#' s <- analysis_spec("t_means_one_sample", "post_hoc",
#'   inputs = list(d = .5, alpha = .05, n_total = 30),
#'   options = list(tails = "two"))
#' power_plot_data(calculate_power(s), "decision_table")
power_plot_data <- function(result, panel = c("distribution", "decision_table", "mosaic", "lattice"),
                            points = 401L) {
  .ig_power_plot_data(result, panel, points)
}

.ig_power_plot_data <- function(result, panel = c("distribution", "decision_table", "mosaic", "lattice"),
                                points = 401L, legacy_t = FALSE) {
  panel <- match.arg(panel)
  plot_law <- if (legacy_t) .ig_plot_law_legacy_t else .ig_plot_law
  if (!inherits(result, "power_result"))
    .ig_abort("invalid_input", "`result` must be a power_result")
  envelope <- analysis_envelope(result$spec, result)
  if (!is.numeric(points) || length(points) != 1L || !is.finite(points) ||
      points != floor(points) || points < 51 || points > 2001)
    .ig_abort("invalid_input", "`points` must be an integer from 51 to 2001")
  out <- list(schema_version = "1.0.0", analysis_kind = "power",
              panel = panel, status = "unavailable", reason = NULL,
              spec_hash = envelope$spec_hash, result_hash = envelope$result_hash,
              analysis_hash = envelope$analysis_hash, data = data.frame(),
              fixed_inputs = result$spec$inputs, provenance = result$provenance,
              units = list(probability = "conditional_probability"),
              critical = result$quantities$critical, notes = character())
  if (result$status$code != "ok") {
    out$reason <- paste("Analysis status:", result$status$code)
    return(out)
  }
  if (panel == "lattice") return(.ig_power_boschloo_lattice(result, out))
  if (panel != "distribution") {
    a <- result$achieved$alpha_attained
    p <- result$achieved$power
    b <- result$achieved$beta
    if (length(a) != 1L || length(p) != 1L || length(b) != 1L ||
        any(!is.finite(c(a, p, b))) || any(c(a, p, b) < 0 | c(a, p, b) > 1)) {
      out$reason <- "Finite attained alpha, power and beta are required."
      return(out)
    }
    out$data <- data.frame(hypothesis = rep(c("H0", "H1"), each = 2L),
                           decision = rep(c("reject", "retain"), 2L),
                           probability = c(a, 1 - a, p, b),
                           quantity = c("alpha_attained", "1-alpha_attained", "power", "beta"))
    if (panel == "mosaic") {
      out$data$xmin <- rep(c(0, 1), each = 2L)
      out$data$xmax <- out$data$xmin + 1
      out$data$ymin <- c(0, a, 0, p)
      out$data$ymax <- c(a, 1, p, 1)
    }
    out$notes <- "Each hypothesis is a separate conditional column; equal column widths do not imply equal prior probabilities."
    out$status <- "ok"
    return(out)
  }
  dist <- result$method_outputs$statistic_distribution
  supported <- c("normal", "t", "F", "chisq", "scaled_F", "scaled_chisq", "binomial")
  laws <- list(H0 = dist$h0, H1 = dist$h1)
  if (!all(vapply(laws, function(law) is.list(law) && length(law$family) == 1L &&
                  law$family %in% supported, logical(1)))) {
    out$reason <- "This method needs a specialised distribution panel."
    return(out)
  }
  discrete <- vapply(laws, function(law) identical(law$family, "binomial"), logical(1))
  if (any(discrete) && !all(discrete)) {
    out$reason <- "Mixed discrete and continuous descriptors are not supported."
    return(out)
  }
  if (all(discrete) && any(vapply(laws, function(law) law$size > 2000, logical(1)))) {
    out$reason <- "The full discrete support exceeds the 2001-point display limit."
    return(out)
  }
  side <- dist$rejection_side
  crit <- if (all(discrete)) unlist(dist$rejection_region, use.names = FALSE) else out$critical
  if (length(side) != 1L || !side %in% c("lower", "upper", "both") ||
      !is.numeric(crit) || anyNA(crit) ||
      (!all(discrete) && length(crit) != if (side == "both") 2L else 1L) ||
      (all(discrete) && length(crit) != 2L)) {
    out$reason <- "A supported rejection region is not available."
    return(out)
  }
  reject <- function(x) {
    if (all(discrete)) return(x <= crit[1] | x >= crit[2])
    switch(side, lower = x <= crit[1], upper = x >= crit[1],
           both = x <= crit[1] | x >= crit[2])
  }
  area <- function(law) {
    if (all(discrete)) return(.ig_p_binom(crit[1], law$size, law$prob) +
      .ig_p_binom(crit[2] - 1, law$size, law$prob, lower_tail = FALSE))
    switch(side,
      lower = .ig_plot_cdf(law, crit[1], TRUE),
      upper = .ig_plot_cdf(law, crit[1], FALSE),
      both = .ig_plot_cdf(law, crit[1], TRUE) + .ig_plot_cdf(law, crit[2], FALSE))
  }
  expected <- c(result$achieved$alpha_attained, result$achieved$power)
  areas <- vapply(laws, area, numeric(1))
  # This is a display consistency guard, not a reference acceptance tolerance.
  if (length(expected) != 2L || any(!is.finite(expected)) ||
      any(!is.finite(areas)) || any(abs(areas - expected) > 1e-12 + 1e-6 * abs(expected))) {
    out$reason <- "Distribution descriptors do not reproduce the engine rejection probabilities; use the decision table."
    return(out)
  }
  out$region_probabilities <- data.frame(hypothesis = names(laws),
                                        engine = expected, distribution = unname(areas))
  curves <- tryCatch(lapply(names(laws), function(h) {
    law <- laws[[h]]
    x <- if (all(discrete)) seq.int(0, law$size) else {
      # Quantile-spaced points retain both narrow peaks and long tails.
      q <- plot_law(law, seq(1e-5, 1 - 1e-5, length.out = points), "q")
      if (any(!is.finite(q))) q else
        sort(unique(c(q, seq(q[1], q[length(q)], length.out = points), crit[is.finite(crit)])))
    }
    y <- plot_law(law, x, "d")
    data.frame(hypothesis = h, x = x, y = y, reject = reject(x))
  }), warning = identity, error = identity)
  if (inherits(curves, "condition")) {
    out$reason <- paste("Distribution display could not be computed:", conditionMessage(curves))
    return(out)
  }
  out$data <- do.call(rbind, curves)
  if (any(!is.finite(out$data$x)) || any(!is.finite(out$data$y))) {
    out$data <- data.frame()
    out$reason <- "Non-finite display coordinates; use the decision table."
    return(out)
  }
  out$ordinate <- if (all(discrete)) "probability_mass" else "density"
  out$units <- list(x = "test_statistic", y = out$ordinate)
  out$notes <- if (all(discrete)) "Full discrete support; boundary counts are included in rejection." else
    "Display uses central 99.998% quantiles of each law plus finite critical values. Curve sampling and display extent do not determine reported probabilities."
  out$status <- "ok"
  out
}

.ig_plot_cdf <- function(law, q, lower) {
  ncp <- law$ncp %||% 0
  scale <- law$scale %||% 1
  switch(law$family,
    normal = .ig_p_norm(q, law$mean, law$sd, lower_tail = lower),
    t = .ig_p_t(q, law$df, ncp = ncp, lower_tail = lower),
    F = .ig_p_f(q, law$df1, law$df2, ncp = ncp, lower_tail = lower),
    chisq = .ig_p_chisq(q, law$df, ncp = ncp, lower_tail = lower),
    scaled_F = .ig_p_f(q / scale, law$df1, law$df2, lower_tail = lower),
    scaled_chisq = .ig_p_chisq(q / scale, law$df, lower_tail = lower))
}

.ig_plot_law <- function(law, x, operation) {
  ncp <- law$ncp %||% 0
  scale <- law$scale %||% 1
  quantile <- identical(operation, "q")
  switch(law$family,
    normal = if (quantile) stats::qnorm(x, law$mean, law$sd) else stats::dnorm(x, law$mean, law$sd),
    t = if (quantile) .ig_q_t(x, law$df, ncp = ncp) else .ig_plot_t_density(x, law$df, ncp),
    F = if (quantile) stats::qf(x, law$df1, law$df2, ncp = ncp) else stats::df(x, law$df1, law$df2, ncp = ncp),
    chisq = if (quantile) stats::qchisq(x, law$df, ncp = ncp) else stats::dchisq(x, law$df, ncp = ncp),
    scaled_F = if (quantile) scale * stats::qf(x, law$df1, law$df2) else stats::df(x / scale, law$df1, law$df2) / scale,
    scaled_chisq = if (quantile) scale * stats::qchisq(x, law$df) else stats::dchisq(x / scale, law$df) / scale,
    binomial = stats::dbinom(x, law$size, law$prob))
}

# Use the checked mixture density: nmath can lose precision without a warning.
.ig_plot_t_density <- function(x,df,ncp) {
  if(ncp==0) return(stats::dt(x,df))
  if(is.infinite(df)) return(stats::dnorm(x,mean=ncp))
  vapply(x,function(q) .ig_calc_nct_pdf(q,df,ncp),0)
}


.ig_power_boschloo_lattice <- function(result, out) {
  if (!identical(result$method_id, .bl_method_id)) {
    out$reason <- "The table lattice is available for Boschloo designs only."
    return(out)
  }
  n <- unlist(result$group_sizes[c("n1", "n2")], use.names = FALSE)
  if (length(n) != 2L || any(!is.finite(n)) || any(n < 1) ||
      any(n != floor(n)) || prod(n + 1) > 5000) {
    out$reason <- "The complete table lattice exceeds the 5000-table display limit; use the decision table."
    return(out)
  }
  inputs <- result$spec$inputs
  control <- .bl_control_from_record(result$method_outputs$settings)
  pd <- tryCatch(boschloo_plotdata_lattice(n[1], n[2], alpha = inputs$alpha,
    variant = result$spec$options$variant, p1 = inputs$p1, p2 = inputs$p2,
    pvalues = FALSE, control = control), error = identity, warning = identity)
  if (inherits(pd, "condition") || !identical(pd$status$code, "ok")) {
    out$reason <- "The complete table lattice could not be reproduced; use the decision table."
    return(out)
  }
  tab <- pd$data
  sums <- c(sum(tab$prob_null[tab$in_region]), sum(tab$prob_alternative[tab$in_region]),
            sum(tab$prob_alternative[!tab$in_region]))
  expected <- c(result$achieved$alpha_attained, result$achieved$power, result$achieved$beta)
  if (length(expected) != 3L || any(!is.finite(c(sums, expected))) ||
      any(abs(sums - expected) > 1e-12 + 1e-6 * abs(expected))) {
    out$reason <- "Lattice probabilities do not reproduce the computed result; use the decision table."
    return(out)
  }
  out$data <- rbind(data.frame(hypothesis = "H0", x1 = tab$x1, x2 = tab$x2,
    probability = tab$prob_null, reject = tab$in_region),
    data.frame(hypothesis = "H1", x1 = tab$x1, x2 = tab$x2,
    probability = tab$prob_alternative, reject = tab$in_region))
  out$fixed_inputs <- c(pd$fixed_inputs, list(settings = result$method_outputs$settings))
  out$region_probabilities <- data.frame(quantity = c("alpha_attained", "power", "beta"),
    engine = expected, lattice = sums)
  out$units <- list(x = "group_1_event_count", y = "group_2_event_count",
    probability = "conditional_table_probability")
  out$ordinate <- "table_probability_mass"
  out$notes <- paste("Full possible-table support; crosses mark rejection, not observed data.",
    "H0 uses a size-maximising common event probability; H1 uses the specified alternative.",
    "The two panels are separate conditional distributions, not joint or predictive probabilities.")
  out$status <- "ok"
  out
}


# Schema 0.1 snapshots can contain the former R noncentral-t display grid.
# Reproduce only that known renderer when validating retained coordinates.
.ig_plot_law_legacy_t <- function(law, x, operation) {
  if (!identical(law$family, "t")) return(.ig_plot_law(law, x, operation))
  if (operation == "q") stats::qt(x, law$df, ncp = law$ncp %||% 0) else
    stats::dt(x, law$df, ncp = law$ncp %||% 0)
}
