# Logistic and Poisson regression (z tests; group G6): definitions, forward
# models against an independent reference, G*Power conventions of the
# enumeration procedure, typed outcomes, mode inversion, the sensitivity
# search and the effect-size drawer.
#
# The reference functions here share no numerical routine with the engine
# (decision D18): Poisson-regression variances from the closed-form moment
# generating functions (manual 30.5.2 with the inventory corrections), the
# logistic information by an own Gauss-Legendre rule (Golub-Welsch nodes) on
# the standardised covariate, binary covariates in closed form, and normal
# tail probabilities by integrating the density. The full independent check
# (all distributions, procedures and modes) is
# validation/independent/g6/check_g6.R. Tolerances: probabilities abs 1e-8 +
# rel 1e-6 (validation/tolerances.md section 5).

lg <- "z_logistic_regression"
ps <- "z_poisson_regression"

# Reference building blocks -------------------------------------------------------------

g6r_gl <- local({
  n <- 48L
  i <- seq_len(n - 1L)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- b
  J[cbind(i + 1L, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = 2 * e$vectors[1, ]^2)
})
g6r_int <- function(f, a, b, panels = 200L) {
  edges <- seq(a, b, length.out = panels + 1L)
  h <- diff(edges) / 2
  mid <- (edges[-1L] + edges[-length(edges)]) / 2
  x <- outer(g6r_gl$x, h) + rep(mid, each = length(g6r_gl$x))
  sum(colSums(matrix(f(as.vector(x)), nrow = length(g6r_gl$x)) * g6r_gl$w) * h)
}
g6r_phi <- function(z) exp(-z^2 / 2) / sqrt(2 * pi)
g6r_Phi <- function(x) if (x > 0) 1 - g6r_Phi(-x) else if (x < -38) 0 else g6r_int(g6r_phi, max(x - 12, -40), x, 16L)
g6r_zq <- function(a) {
  lo <- -40; hi <- 40
  for (i in 1:120) { m <- (lo + hi) / 2; if (g6r_Phi(-m) > a) lo <- m else hi <- m }
  (lo + hi) / 2
}
g6r_zpower <- function(m1, s1, alpha, tails) {
  m <- abs(m1)
  if (tails == "one") { z <- g6r_zq(alpha); g6r_Phi((m - z) / s1) } else {
    z <- g6r_zq(alpha / 2); g6r_Phi((m - z) / s1) + g6r_Phi((-m - z) / s1) }
}
# Logistic v for a normal covariate (centred), by Gauss-Legendre.
g6r_logit_v_normal <- function(mu, sigma, b0, b1) {
  w <- function(z) { e <- exp(-abs(b0 + b1 * (mu + sigma * z))); g6r_phi(z) * e / (1 + e)^2 }
  a00 <- g6r_int(w, -40, 40); a01 <- g6r_int(function(z) sigma * z * w(z), -40, 40)
  a11 <- g6r_int(function(z) (sigma * z)^2 * w(z), -40, 40)
  a00 / (a00 * a11 - a01^2)
}
g6r_mean_p_normal <- function(mu, sigma, b0, b1) g6r_int(function(z) g6r_phi(z) / (1 + exp(-(b0 + b1 * (mu + sigma * z)))), -40, 40)

g6_spec <- function(method, mode, inputs, options = list(), target = NULL) {
  args <- list(method = method, mode = mode, inputs = inputs, options = options)
  if (!is.null(target)) args$target <- target
  do.call(analysis_spec, args)
}
g6_close <- function(a, b, rel = 1e-6, abs = 1e-8) isTRUE(abs(a - b) <= abs + rel * max(abs(a), abs(b)))

# Definitions ---------------------------------------------------------------------------

test_that("both methods are registered with the certified options, defaults and all five modes", {
  for (id in c(lg, ps)) {
    def <- .ig_get_method(id)
    expect_s3_class(def, "ig_method")
    expect_identical(def$modes, c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"))
    expect_identical(def$status, "independently_checked")
    expect_identical(validation_status(id)$gpower_reference$status, "reference_pending")
    meta <- .ig_get_metadata(id)
    for (o in meta$options) {
      if (identical(o$value_type, "continuous")) next
      if (identical(o$option_id, "over_dispersion")) next
      expect_true(o$option_id %in% names(def$options), info = paste(id, o$option_id))
      expect_setequal(def$options[[o$option_id]]$values, o$codes)
      expect_identical(def$options[[o$option_id]]$default, o$default, info = paste(id, o$option_id))
    }
    expect_identical(def$params$r2_other_x$default, 0)
    expect_identical(def$params$enumeration_threshold$default, 10000L)
  }
  # Over-dispersion was never displayed by the running application (registry live.present = FALSE).
  expect_false(.ig_get_metadata(ps)$options$over_dispersion$live$present)
  expect_null(.ig_get_method(ps)$params$over_dispersion)
})

# Forward models against the independent reference -----------------------------------------

test_that("logistic Demidenko power matches the reference (normal covariate, both corrections, both tails)", {
  for (cs in list(list(mu = 0, sigma = 1, or = 1.5, p1 = 0.5, n = 200, r2 = 0),
                  list(mu = 5, sigma = 0.1, or = 10, p1 = 0.05, n = 5979346, r2 = 0.9),
                  list(mu = -1, sigma = 2, or = 0.7, p1 = 0.2, n = 150, r2 = 0.3))) {
    b0 <- qlogis(cs$p1); b1 <- log(cs$or)
    v1 <- g6r_logit_v_normal(cs$mu, cs$sigma, b0, b1)
    mu <- g6r_mean_p_normal(cs$mu, cs$sigma, b0, b1)
    v0 <- 1 / (mu * (1 - mu) * cs$sigma^2)
    for (vc in c("on", "off")) for (tails in c("one", "two")) {
      a <- if (vc == "on") 1 else 0
      ref <- g6r_zpower(sqrt(cs$n * (1 - cs$r2) / v1) * b1, sqrt((a * v0 + (1 - a) * v1) / v1), 0.05, tails)
      r <- calculate_power(g6_spec(lg, "post_hoc", list(odds_ratio = cs$or, p1 = cs$p1, alpha = 0.05, n_total = cs$n,
                                                        r2_other_x = cs$r2, x_mu = cs$mu, x_sigma = cs$sigma),
                                   list(x_distribution = "normal", variance_correction = vc, tails = tails)))
      expect_true(g6_close(r$achieved$power, ref), info = paste(cs$mu, cs$sigma, vc, tails, r$achieved$power, ref))
    }
  }
})

test_that("logistic binary covariate: closed-form variances, a = 0.85, Hsieh equation (2)", {
  p1 <- 0.05; p2 <- 0.1; B <- 0.2; n <- 2158
  v1 <- 1 / ((1 - B) * p1 * (1 - p1)) + 1 / (B * p2 * (1 - p2))
  mu <- (1 - B) * p1 + B * p2
  v0 <- 1 / (mu * (1 - mu) * B * (1 - B))
  b1 <- qlogis(p2) - qlogis(p1)
  ref <- g6r_zpower(sqrt(n / v1) * b1, sqrt((0.85 * v0 + 0.15 * v1) / v1), 0.05, "two")
  r <- calculate_power(g6_spec(lg, "post_hoc", list(p2 = p2, p1 = p1, alpha = 0.05, n_total = n, x_pi = B),
                               list(x_distribution = "binomial", tails = "two", effect_size_input_mode = "two_probabilities")))
  expect_true(g6_close(r$achieved$power, ref))
  # Hsieh et al. (1998) eq. (2): only the rejection region on the side of the effect.
  z <- g6r_zq(0.025)
  arg <- (sqrt(n * (1 - B)) * (p2 - p1) - z * sqrt(mu * (1 - mu) / B)) / sqrt(p1 * (1 - p1) + p2 * (1 - p2) * (1 - B) / B)
  r <- calculate_power(g6_spec(lg, "post_hoc", list(p2 = p2, p1 = p1, alpha = 0.05, n_total = n, x_pi = B),
                               list(x_distribution = "binomial", tails = "two", test_procedure = "hsieh",
                                    effect_size_input_mode = "two_probabilities")))
  expect_true(g6_close(r$achieved$power, g6r_Phi(arg)))
  expect_identical(r$method_outputs$procedure, "hsieh")
  # Equation (1), standard normal X; other normal covariates switch to Demidenko with correction.
  r1 <- calculate_power(g6_spec(lg, "post_hoc", list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, n_total = 317),
                                list(tails = "two", test_procedure = "hsieh")))
  expect_true(g6_close(r1$achieved$power, g6r_Phi(sqrt(317 * 0.25) * log(1.5) - g6r_zq(0.025))))
  r2 <- calculate_power(g6_spec(lg, "post_hoc", list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, n_total = 200, x_sigma = 2),
                                list(tails = "one", test_procedure = "hsieh")))
  r3 <- calculate_power(g6_spec(lg, "post_hoc", list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, n_total = 200, x_sigma = 2),
                                list(tails = "one", test_procedure = "demidenko", variance_correction = "on")))
  expect_identical(r2$achieved$power, r3$achieved$power)
  expect_true(any(grepl("Hsieh et al. \\(1998\\) covers only", r2$assumptions)))
  # G*Power shows 0.9999049 for this request (supplement capture SUP-G6-Z-LOGIT-POS-H01--r01).
  expect_true(abs(r2$achieved$power - 0.9999049) <= 5e-8)
})

test_that("Poisson regression variances equal the closed-form moment generating functions", {
  b0 <- log(0.85); b1 <- log(1.3)
  mgf_v <- list(
    binomial = list(par = list(x_pi = 0.3), v1 = function() (1 / 0.7 + 1 / (0.3 * exp(b1))) / exp(b0),
                    v0 = function() 1 / (0.3 * 0.7 * exp(b0))),
    normal = list(par = list(x_mu = 0.5, x_sigma = 2), v1 = function() exp(-(b1 * 0.5 + (b1 * 2)^2 / 2)) / 4 / exp(b0),
                  v0 = function() 1 / (4 * exp(b0))),
    exponential = list(par = list(x_lambda = 3), v1 = function() (3 - b1)^3 / 3 / exp(b0), v0 = function() 9 / exp(b0)),
    poisson = list(par = list(x_lambda = 0.5), v1 = function() exp(-b1 + 0.5 - exp(b1) * 0.5) / 0.5 / exp(b0),
                   v0 = function() 1 / (0.5 * exp(b0))),
    uniform = list(par = list(x_left = -1, x_right = 2), v1 = function() {
      A <- -1; B <- 2; D <- 3
      F0 <- function(x) exp(b1 * x) / b1; F1 <- function(x) exp(b1 * x) * (x / b1 - 1 / b1^2)
      F2 <- function(x) exp(b1 * x) * (x^2 / b1 - 2 * x / b1^2 + 2 / b1^3)
      i0 <- (F0(B) - F0(A)) / D; i1 <- (F1(B) - F1(A)) / D; i2 <- (F2(B) - F2(A)) / D
      i0 / (i0 * i2 - i1^2) / exp(b0) }, v0 = function() 12 / 9 / exp(b0)))
  for (xd in names(mgf_v)) {
    m <- mgf_v[[xd]]
    r <- calculate_power(g6_spec(ps, "post_hoc", c(list(exp_beta1 = 1.3, base_rate = 0.85, alpha = 0.05, n_total = 200), m$par),
                                 list(x_distribution = xd, test_procedure = "signorini")))
    expect_true(g6_close(r$method_outputs$v1, m$v1(), rel = 1e-10, abs = 0), info = xd)
    expect_true(g6_close(r$method_outputs$v0, m$v0(), rel = 1e-10, abs = 0), info = xd)
    # Signorini: m1 = beta1 sqrt(N t / v0), s1 = sqrt(v1 / v0)
    ref <- g6r_zpower(b1 * sqrt(200 / m$v0()), sqrt(m$v1() / m$v0()), 0.05, "one")
    expect_true(g6_close(r$achieved$power, ref), info = xd)
  }
  # Manual 30.3 example: normal X, beta1 = -0.1, beta0 = 0.5, N = 200, Demidenko, two tails: 0.444593.
  r <- calculate_power(g6_spec(ps, "post_hoc", list(exp_beta1 = exp(-0.1), base_rate = exp(0.5), alpha = 0.05, n_total = 200),
                               list(tails = "two", variance_correction = "off")))
  expect_true(abs(r$achieved$power - 0.444593) <= 5e-7)
  expect_equal(min(r$quantities$critical), -1.959964, tolerance = 5e-7)
})

test_that("variance correction: a = 0.75 for lognormal X, 1 otherwise; without correction s1 = 1", {
  base <- list(exp_beta1 = 0.9, base_rate = 1, alpha = 0.05, n_total = 300, x_mu = 0, x_sigma = 0.5)
  on <- calculate_power(g6_spec(ps, "post_hoc", base, list(x_distribution = "lognormal", variance_correction = "on")))
  off <- calculate_power(g6_spec(ps, "post_hoc", base, list(x_distribution = "lognormal", variance_correction = "off")))
  mo <- on$method_outputs
  expect_identical(mo$variance_correction_a, 0.75)
  expect_true(g6_close(mo$s1, sqrt((0.75 * mo$v0_star + 0.25 * mo$v1) / mo$v1), rel = 1e-14, abs = 0))
  expect_identical(off$method_outputs$s1, 1)
  lgon <- calculate_power(g6_spec(lg, "post_hoc", list(odds_ratio = 2, p1 = 0.3, alpha = 0.05, n_total = 100, x_pi = 0.4),
                                  list(x_distribution = "binomial")))
  expect_identical(lgon$method_outputs$variance_correction_a, 0.85)
})

test_that("power equals alpha at the null; one tail rejects on the side of beta1", {
  r <- calculate_power(g6_spec(lg, "post_hoc", list(odds_ratio = 1, p1 = 0.3, alpha = 0.05, n_total = 100)))
  expect_true(g6_close(r$achieved$power, 0.05))
  neg <- calculate_power(g6_spec(ps, "post_hoc", list(exp_beta1 = 0.8, base_rate = 1, alpha = 0.05, n_total = 100)))
  expect_lt(neg$quantities$critical, 0)
  expect_identical(neg$method_outputs$statistic_distribution$rejection_side, "lower")
})

# Enumeration procedure ----------------------------------------------------------------------

test_that("enumeration: Blom grid of round-half-up N (1 - R2) points and the chi-square(1) test", {
  expect_identical(.ig_glm_grid_size(30925L, 0.5), 15463)        # R's round() would give 15462
  expect_identical(.ig_glm_grid_size(705L, 0.9), 70)             # 705 * (1 - 0.9) = 70.49999999999999
  expect_identical(.ig_glm_grid_size(20L, 0.9), 2)
  expect_identical(.ig_glm_grid(list(type = "uniform", a = 0, b = 1), 3L), ((1:3) - 0.375) / 3.25)
  # Discrete covariates: the strict rule F(x) > q puts the middle position of an odd grid (q = 0.5) at X = 1.
  x <- .ig_glm_grid(list(type = "binomial", pi = 0.5), 615L)
  expect_identical(sum(x), 308)
  xp <- .ig_glm_grid(list(type = "poisson", lambda = 1), 200L)
  q <- ((1:200) - 0.375) / 200.25
  expect_true(all(stats::ppois(xp, 1) > q & (xp == 0 | stats::ppois(xp - 1, 1) <= q)))
  # Wald noncentrality from the explicit 2 x 2 information matrix of the grid.
  b0 <- 0; b1 <- log(1.5); ng <- 200
  xg <- stats::qnorm(((1:ng) - 0.375) / (ng + 0.25))
  pr <- 1 / (1 + exp(-(b0 + b1 * xg))); w <- pr * (1 - pr)
  I <- matrix(c(sum(w), sum(w * xg), sum(w * xg), sum(w * xg^2)), 2)
  lam_ref <- b1^2 / solve(I)[2, 2]
  r <- calculate_power(g6_spec(lg, "post_hoc", list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, n_total = 200,
                                                    enumeration_threshold = 100000L),
                               list(computation = "enumeration", test_statistic = "wald", tails = "one")))
  expect_true(g6_close(r$quantities$ncp, lam_ref, rel = 1e-12, abs = 0))
  # G*Power 3.1.9.7 shows 7.2446808 (capture GPW-Z-LOGIT-POS-B066-O01--r01).
  expect_true(abs(r$quantities$ncp - 7.2446808) <= 5e-8)
  # One tail: critical chi-square at 2 alpha; the power counts both tails.
  expect_true(g6_close(r$quantities$critical, g6r_zq(0.05)^2, rel = 1e-10, abs = 0))
  rl <- sqrt(lam_ref); z <- g6r_zq(0.05)
  expect_true(g6_close(r$achieved$power, g6r_Phi(rl - z) + g6r_Phi(-rl - z)))
  expect_identical(r$quantities$df, 1)
  # Likelihood ratio: deviance difference to the fitted intercept-only model.
  pb <- mean(pr)
  lr_ref <- 2 * (sum(pr * log(pr) + (1 - pr) * log(1 - pr)) - ng * (pb * log(pb) + (1 - pb) * log(1 - pb)))
  r <- calculate_power(g6_spec(lg, "post_hoc", list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, n_total = 200,
                                                    enumeration_threshold = 100000L),
                               list(computation = "enumeration", test_statistic = "likelihood_ratio")))
  expect_true(g6_close(r$quantities$ncp, lr_ref, rel = 1e-10, abs = 0))
  # N not below the threshold: Demidenko with variance correction (G*Power: SUP-G6-Z-LOGIT-POS-T01--r01, 0.8696578).
  r <- calculate_power(g6_spec(lg, "post_hoc", list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, n_total = 200,
                                                    enumeration_threshold = 100L),
                               list(computation = "enumeration", test_statistic = "wald")))
  expect_identical(r$method_outputs$procedure, "demidenko")
  expect_true(abs(r$achieved$power - 0.8696578) <= 5e-8)
  # Poisson: distinguish the untruncated expected-data value from G*Power's
  # finite count enumeration (SUP-G6-Z-POIS-POS-P01--r01; rule P06).
  r <- calculate_power(g6_spec(ps, "post_hoc", list(exp_beta1 = 2, base_rate = 1, alpha = 0.05, n_total = 2, x_pi = 0.5,
                                                    enumeration_threshold = 100000L),
                               list(computation = "enumeration", x_distribution = "binomial", tails = "two")))
  expect_true(g6_close(r$method_outputs$ncp_untruncated, 2 * log(2)^2 / 3, rel = 1e-14, abs = 0))
  # Independent explicit response enumeration. With two distinct X values,
  # each fitted weighted mean is T_i / S_i, so its information weight is T_i.
  # The engine obtains these sufficient statistics through Poisson CDFs and
  # Newton fitting; here sum the individual count probabilities directly.
  y <- 0:10
  expected_events <- c(sum(y * stats::dpois(y, 1)), sum(y * stats::dpois(y, 2)))
  lambda_ref <- log(2)^2 * prod(expected_events) / sum(expected_events)
  expect_equal(r$method_outputs$enumeration_count_range, c(0, 10))
  expect_true(g6_close(r$quantities$ncp, lambda_ref, rel = 1e-14, abs = 0))
  expect_lt(r$quantities$ncp, r$method_outputs$ncp_untruncated)
})

# Typed outcomes --------------------------------------------------------------------------

test_that("undefined models and inconsistent requests are refused without numbers", {
  # Lognormal covariate with Exp(beta1) > 1: infinite information integrals.
  expect_error(g6_spec(ps, "post_hoc", list(exp_beta1 = 1.3, base_rate = 1, alpha = 0.05, n_total = 100),
                       list(x_distribution = "lognormal")), class = "inferguidR_invalid_input")
  expect_error(g6_spec(ps, "post_hoc", list(exp_beta1 = exp(2), base_rate = 1, alpha = 0.05, n_total = 100, x_lambda = 1.5),
                       list(x_distribution = "exponential")), class = "inferguidR_invalid_input")
  # Manual X is not offered with the enumeration procedure.
  expect_error(g6_spec(lg, "post_hoc", list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, n_total = 100, x_var_h0 = 4, x_var_h1 = 4.2),
                       list(x_distribution = "manual", computation = "enumeration")), class = "inferguidR_invalid_input")
  # Target power below alpha (G*Power: 'Power must be greater than or equal to alpha.').
  expect_error(g6_spec(lg, "a_priori", list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, power = 0.04)),
               class = "inferguidR_invalid_input")
  expect_error(g6_spec(lg, "post_hoc", list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, n_total = 100, x_left = 2, x_right = 1),
                       list(x_distribution = "uniform")), class = "inferguidR_invalid_input")
  # Manual X in sensitivity mode: v1 depends on the effect being solved (manual 30.5.2).
  r <- calculate_power(g6_spec(ps, "sensitivity", list(base_rate = 1, alpha = 0.05, power = 0.8, n_total = 100),
                               list(x_distribution = "manual")))
  expect_identical(r$status$code, "invalid_input")
  expect_null(r$solved$exp_beta1)
  # No effect: a priori infeasible.
  r <- calculate_power(g6_spec(lg, "a_priori", list(odds_ratio = 1, p1 = 0.5, alpha = 0.05, power = 0.8)))
  expect_identical(r$status$code, "infeasible")
})

# Modes --------------------------------------------------------------------------------------

test_that("every mode inverts the forward model", {
  common <- list(odds_ratio = 0.7, p1 = 0.3, x_pi = 0.4)
  opts <- list(x_distribution = "binomial", tails = "two")
  ap <- calculate_power(g6_spec(lg, "a_priori", c(common, alpha = 0.05, power = 0.9), opts))
  N <- ap$solved$n_total
  ph <- function(n, a = 0.05) calculate_power(g6_spec(lg, "post_hoc", c(common, alpha = a, n_total = n), opts))
  expect_gte(ph(N)$achieved$power, 0.9)
  expect_lt(ph(N - 1)$achieved$power, 0.9)
  cr <- calculate_power(g6_spec(lg, "criterion", c(common, power = 0.9, n_total = 400), opts))
  expect_true(g6_close(ph(400, cr$solved$alpha)$achieved$power, 0.9, rel = 1e-9))
  co <- calculate_power(g6_spec(lg, "compromise", c(common, error_ratio = 2, n_total = 400), opts))
  expect_true(g6_close(co$achieved$beta, 2 * co$solved$alpha, rel = 1e-8))
  # Poisson enumeration a priori: smallest N on the grid-size-dependent power curve.
  pe <- list(exp_beta1 = 1.3, base_rate = 0.85, x_pi = 0.5, enumeration_threshold = 100000L)
  po <- list(computation = "enumeration", x_distribution = "binomial", tails = "one")
  ap <- calculate_power(g6_spec(ps, "a_priori", c(pe, alpha = 0.05, power = 0.8), po))
  Np <- ap$solved$n_total
  pp <- function(n) calculate_power(g6_spec(ps, "post_hoc", c(pe, alpha = 0.05, n_total = n), po))$achieved$power
  expect_gte(pp(Np), 0.8)
  expect_true(all(vapply(max(2, Np - 20):(Np - 1), pp, numeric(1)) < 0.8))
})

test_that("sensitivity search: first crossing from the null in the selected direction", {
  model <- .ig_power_z_logistic_regression
  search <- .ig_glm_sensitivity_search(model)
  spec <- g6_spec(lg, "sensitivity", list(p1 = 0.5, alpha = 0.05, power = 0.8, n_total = 200, x_pi = 0.5),
                  list(x_distribution = "binomial", tails = "one", effect_direction = "p2_le_p1", variance_correction = "off"))
  p <- c(spec$inputs, spec$options)
  s <- search(p, "odds_ratio")
  expect_identical(s$status, "ok")
  expect_true(s$diagnostics$target_reached)
  expect_lt(s$value, 1)
  # Two crossings exist for a binary covariate below the null; the one closest to the null is reported.
  expect_gte(s$diagnostics$sign_changes_on_grid, 2)
  q <- p; q$odds_ratio <- s$value
  expect_true(abs(model(q)$power - 0.8) < 1e-9)
  # G*Power 3.1.9.7: odds ratio 0.4832796 (capture GPW-Z-LOGIT-SEN-B012-O01--r01).
  expect_true(abs(s$value - 0.4832796) <= 5e-8)
  # Every effect between the null and the solution has power below the target.
  ors <- exp(seq(log(s$value) * 0.999, 0, length.out = 50))
  pw <- vapply(ors, function(o) { q$odds_ratio <- o; model(q)$power }, numeric(1))
  expect_true(all(pw < 0.8))
  # Unattainable target: the power-maximising effect, with a warning and target_reached = FALSE.
  spec2 <- g6_spec(lg, "sensitivity", list(p1 = 0.05, alpha = 0.05, power = 0.95, n_total = 20, r2_other_x = 0.9, x_pi = 0.02),
                   list(x_distribution = "binomial", tails = "one", effect_direction = "p2_le_p1", variance_correction = "off"))
  s2 <- search(c(spec2$inputs, spec2$options), "odds_ratio")
  expect_identical(s2$status, "ok")
  expect_false(s2$diagnostics$target_reached)
  expect_match(s2$warnings, "not reached")
  # G*Power: actual power 0.0535486 at the power-maximising odds ratio (GPW-Z-LOGIT-SEN-B012-E01--r01).
  expect_true(abs(s2$diagnostics$max_power - 0.0535486) <= 5e-8)
  # Poisson: lognormal covariate, direction beta1 >= 0 has no defined solution.
  sp <- g6_spec(ps, "sensitivity", list(base_rate = 1, alpha = 0.05, power = 0.8, n_total = 100),
                list(x_distribution = "lognormal", effect_direction = "beta1_ge_0"))
  expect_identical(.ig_glm_sensitivity_search(.ig_power_z_poisson_regression)(c(sp$inputs, sp$options), "exp_beta1")$status,
                   "invalid_input")
})

test_that("the public sensitivity search stays inside the selected direction", {
  r <- calculate_power(g6_spec(ps, "sensitivity", list(base_rate = 0.85, alpha = 0.05, power = 0.8, n_total = 200),
                               list(effect_direction = "beta1_le_0")))
  expect_identical(r$status$code, "ok")
  expect_identical(r$diagnostics$solver, "sensitivity_method_search")
  expect_lt(r$solved$exp_beta1, 1)
  q <- calculate_power(g6_spec(ps, "post_hoc", list(exp_beta1 = r$solved$exp_beta1, base_rate = 0.85, alpha = 0.05,
                                                     n_total = 200)))
  expect_true(abs(q$achieved$power - 0.8) < 1e-8)
})

test_that("results are deterministic (the variance cache returns what a fresh evaluation returns)", {
  spec <- g6_spec(lg, "post_hoc", list(odds_ratio = 1.3, p1 = 0.2, alpha = 0.05, n_total = 150, x_mu = 1, x_sigma = 0.5),
                  list(x_distribution = "lognormal"))
  r1 <- calculate_power(spec)
  rm(list = ls(.ig_glm_cache, all.names = TRUE), envir = .ig_glm_cache)
  r2 <- calculate_power(spec)
  expect_identical(r1$achieved, r2$achieved)
  expect_identical(r1$result_hash, r2$result_hash)
})

# Effect-size drawer ---------------------------------------------------------------------------

test_that("the logistic drawer converts between the two probabilities and the odds ratio", {
  es <- effect_size(lg, "probabilities", p1 = 0.5, p2 = 0.6)
  expect_identical(es$status$code, "ok")
  expect_equal(es$value, 1.5, tolerance = 1e-14)
  # G*Power shows 998001 for p1 = 0.001, p2 = 0.999 (capture GPW-UT-DRAW-UTL-U068-UE01--r01).
  expect_equal(effect_size(lg, "probabilities", p1 = 0.001, p2 = 0.999)$value, 998001, tolerance = 1e-12)
  back <- effect_size(lg, "odds_ratio", p1 = 0.2, odds_ratio = 1.3)
  expect_equal(back$value, 0.2 * 1.3 / (0.8 + 0.2 * 1.3), tolerance = 1e-14)
  expect_error(effect_size(lg, "probabilities", p1 = 0.5, p2 = 1), class = "inferguidR_invalid_input")
  expect_identical(effect_size(ps, "x", a = 1)$status$code, "unsupported")
})
