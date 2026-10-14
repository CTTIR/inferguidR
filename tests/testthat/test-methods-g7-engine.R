# Group G7: correlation z tests -- definitions, forward models, modes ---------------------
#
# Methods: z_correlations_two_independent, z_correlations_two_dependent_common_index,
# z_correlations_two_dependent_no_common_index, z_tetrachoric_correlation
# (R/methods-correlations-z.R, R/effect-correlations-z.R; method notes
# specs/method-notes/<method_id>.md). The references below are small
# test-local re-derivations that share no routine with the engine: Fisher-z
# variances written out, Steiger covariances from Isserlis' theorem and the
# delta method, tetrachoric cells by Gauss-Legendre quadrature and the SE from
# Pearson's closed form (manual 31.5.1, first term corrected). The exhaustive
# checks are in validation/independent/g7/.

ref_reset_registry()

g7_ids <- c("z_correlations_two_independent", "z_correlations_two_dependent_common_index",
            "z_correlations_two_dependent_no_common_index", "z_tetrachoric_correlation")

g7_calc <- function(method, mode, inputs, options = list()) {
  calculate_power(analysis_spec(method, mode, inputs = inputs, options = options))
}

g7_code <- function(expr) {
  tryCatch({
    force(expr)
    "ok"
  }, inferguidR_error = function(e) e$status$code)
}

# Gauss-Legendre nodes (Golub-Welsch), for the test-local references.
g7_gl <- function(n = 60) {
  k <- seq_len(n - 1)
  b <- k / sqrt(4 * k^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1)] <- b
  J[cbind(k + 1, k)] <- b
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = 2 * e$vectors[1, ]^2)
}
g7_nodes <- g7_gl(60)

# P(X > a, Y > b) by composite Gauss-Legendre quadrature over x in [a, a + 30].
g7_orthant <- function(a, b, rho) {
  s <- sqrt(1 - rho^2)
  f <- function(x) dnorm(x) * pnorm((b - rho * x) / s, lower.tail = FALSE)
  edges <- seq(max(a, -30), max(a, -30) + 30, length.out = 31)
  sum(vapply(seq_len(30), function(i) {
    lo <- edges[i]; hi <- edges[i + 1]
    (hi - lo) / 2 * sum(g7_nodes$w * f((hi - lo) / 2 * g7_nodes$x + (hi + lo) / 2))
  }, numeric(1)))
}

g7_cells <- function(rho, px, py) {
  h <- qnorm(1 - px); k <- qnorm(1 - py)
  c(a = g7_orthant(-h, -k, rho), b = g7_orthant(h, -k, -rho), c = g7_orthant(-h, k, -rho),
    d = g7_orthant(h, k, rho), h = h, k = k)
}

g7_se_pearson <- function(cl, rho, n) {
  a <- cl[["a"]]; b <- cl[["b"]]; c <- cl[["c"]]; d <- cl[["d"]]; h <- cl[["h"]]; k <- cl[["k"]]
  s <- sqrt(1 - rho^2)
  F1 <- pnorm((h - rho * k) / s) - 0.5
  F2 <- pnorm((k - rho * h) / s) - 0.5
  v <- (a + d) * (b + c) / 4 + (a + c) * (b + d) * F2^2 + (a + b) * (c + d) * F1^2 +
    2 * (a * d - b * c) * F1 * F2 - (a * b - c * d) * F2 - (a * c - b * d) * F1
  sqrt(v / n) / (exp(-(h^2 - 2 * rho * h * k + k^2) / (2 * s^2)) / (2 * pi * s))
}

# Steiger covariance (N - 3) Cov(z_ij, z_kl) from Isserlis' theorem and the delta method.
g7_cov_z <- function(R, ij, kl) {
  i <- ij[1]; j <- ij[2]; k <- kl[1]; l <- kl[2]
  cs <- function(a, b, c, d) R[a, c] * R[b, d] + R[a, d] * R[b, c]   # N Cov(s_ab, s_cd)
  r1 <- R[i, j]; r2 <- R[k, l]
  # r_ij = s_ij - r_ij (s_ii + s_jj)/2 to first order
  terms <- list(list(1, i, j), list(-r1 / 2, i, i), list(-r1 / 2, j, j))
  terms2 <- list(list(1, k, l), list(-r2 / 2, k, k), list(-r2 / 2, l, l))
  v <- 0
  for (t1 in terms) for (t2 in terms2) v <- v + t1[[1]] * t2[[1]] * cs(t1[[2]], t1[[3]], t2[[2]], t2[[3]])
  v / ((1 - r1^2) * (1 - r2^2))
}

g7_zpower <- function(mu, sd1, alpha, tails, side) {
  if (tails == "two") {
    zc <- qnorm(1 - alpha / 2)
    return(pnorm((zc - mu) / sd1, lower.tail = FALSE) + pnorm((-zc - mu) / sd1))
  }
  zc <- qnorm(1 - alpha)
  if (side > 0) pnorm((zc - mu) / sd1, lower.tail = FALSE) else pnorm((-zc - mu) / sd1)
}

# Registration and inventory conformity ---------------------------------------------------

test_that("the four G7 methods are registered with the inventory's modes, options and defaults", {
  for (id in g7_ids) {
    def <- .ig_get_method(id)
    expect_false(is.null(def), info = id)
    meta <- .ig_get_metadata(id)
    expect_setequal(def$modes, meta$modes)
    expect_identical(def$family, "z")
    for (o in meta$options) {
      if (identical(o$value_type, "continuous")) {
        expect_identical(def$params[[o$option_id]]$default, as.double(o$default_value), info = paste(id, o$option_id))
        next
      }
      d <- def$options[[o$option_id]]
      expect_false(is.null(d), info = paste(id, o$option_id))
      expect_setequal(d$values, o$codes)
      expect_identical(d$default, o$default, info = paste(id, o$option_id))
    }
    expect_identical(def$status, "independently_checked")
    root <- th_pkg_root()
    if (!is.null(root) && dir.exists(file.path(root, "validation"))) {
      expect_true(file.exists(file.path(root, def$checker)), info = id)   # development checkout only
    }
    vs <- validation_status(id)
    expect_identical(vs$gpower_reference$status, "reference_pending")
  }
})

test_that("the effect direction option applies in sensitivity mode only", {
  spec <- analysis_spec("z_correlations_two_dependent_common_index", "sensitivity",
                        inputs = list(alpha = 0.05, power = 0.8, n_total = 144, rho_ab = 0.4, rho_bc = -0.6))
  expect_identical(spec$options$effect_direction, "rho_ac_le_rho_ab")
  expect_true("effect_direction" %in% spec$defaults_applied$options)
  expect_identical(g7_code(analysis_spec("z_correlations_two_dependent_common_index", "post_hoc",
                                         inputs = list(rho_ac = 0.2, alpha = 0.05, n_total = 144, rho_ab = 0.4, rho_bc = 0.5),
                                         options = list(effect_direction = "rho_ac_ge_rho_ab"))), "invalid_input")
  spec <- analysis_spec("z_tetrachoric_correlation", "post_hoc",
                        inputs = list(rho_h1 = 0.3, rho_h0 = 0, marginal_x = 0.5, marginal_y = 0.5, n_total = 200, alpha = 0.05))
  expect_identical(spec$options$computation_of_correlation, "exact_brown_benedetti")
  expect_identical(spec$options$tails, "one")
})

test_that("inconsistent inputs are refused as G*Power refuses them (typed, no numbers)", {
  # target power below alpha (a priori and sensitivity)
  expect_identical(g7_code(analysis_spec("z_correlations_two_independent", "a_priori",
                                         inputs = list(q = 0.3, alpha = 0.05, power = 0.04))), "invalid_input")
  expect_identical(g7_code(analysis_spec("z_tetrachoric_correlation", "sensitivity",
                                         inputs = list(alpha = 0.05, power = 0.04, n_total = 200, rho_h0 = 0,
                                                       marginal_x = 0.5, marginal_y = 0.5))), "invalid_input")
  # invalid correlation matrix under H0 (manual 28.3.1: rho_bd 0.8 -> 0.9)
  expect_identical(g7_code(analysis_spec("z_correlations_two_dependent_no_common_index", "a_priori",
                                         inputs = list(rho_cd = 0.2, alpha = 0.05, power = 0.8, rho_ab = 0.1, rho_ac = 0.5,
                                                       rho_ad = 0.4, rho_bc = -0.4, rho_bd = 0.9))), "invalid_input")
  # invalid matrix under H1 only
  expect_identical(g7_code(analysis_spec("z_correlations_two_dependent_common_index", "post_hoc",
                                         inputs = list(rho_ac = -0.4, alpha = 0.05, n_total = 50, rho_ab = 0.4, rho_bc = 1))),
                   "invalid_input")
  # domains: group size 2, beta/alpha ratio 0, marginal 1.5
  expect_identical(g7_code(analysis_spec("z_correlations_two_independent", "post_hoc",
                                         inputs = list(q = 0.3, alpha = 0.05, n1 = 2, n2 = 50))), "invalid_input")
  expect_identical(g7_code(analysis_spec("z_correlations_two_independent", "compromise",
                                         inputs = list(q = 0.3, error_ratio = 0, n1 = 100, n2 = 100))), "invalid_input")
  expect_identical(g7_code(analysis_spec("z_tetrachoric_correlation", "post_hoc",
                                         inputs = list(rho_h1 = 0.3, rho_h0 = 0, marginal_x = 1.5, marginal_y = 0.5,
                                                       n_total = 200, alpha = 0.05))), "invalid_input")
})

# Two independent Pearson r's --------------------------------------------------------------

test_that("two independent r's: post hoc power, beta and critical z follow manual 27.5", {
  for (cs in list(list(-0.4028126, 260, 51, 0.05, "two"), list(0.3, 100, 100, 0.05, "one"), list(-0.3, 4, 4, 0.05, "one"),
                  list(0.01, 1e5, 1e5, 1e-4, "two"), list(2, 7, 300, 0.5, "one"))) {
    q <- cs[[1]]; n1 <- cs[[2]]; n2 <- cs[[3]]; a <- cs[[4]]; tl <- cs[[5]]
    r <- g7_calc("z_correlations_two_independent", "post_hoc", list(q = q, alpha = a, n1 = n1, n2 = n2), list(tails = tl))
    mu <- q / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
    ref <- g7_zpower(mu, 1, a, tl, sign(q))
    expect_identical(r$status$code, "ok")
    expect_lt(abs(r$achieved$power - ref), 1e-12)
    expect_lt(abs(r$achieved$beta - (1 - ref)), 1e-12)
    expect_equal(r$quantities$ncp, mu, tolerance = 1e-14)
    if (tl == "one") expect_equal(r$quantities$critical, sign(q) * qnorm(1 - a), tolerance = 1e-14)
  }
  # n = 3: infinite Fisher-z variance, power = alpha (G*Power accepts 3; SUP-G7-Z-CORR2-POS-A01)
  r <- g7_calc("z_correlations_two_independent", "post_hoc", list(q = 0.5, alpha = 0.05, n1 = 3, n2 = 50))
  expect_equal(r$achieved$power, 0.05, tolerance = 1e-15)
})

test_that("two independent r's: a priori allocation rule (equal groups for ratio 1, split rule otherwise)", {
  frac <- function(q, a, pw, r, tl) {
    for (N in 8:5000) {
      f1 <- N / (1 + r); f2 <- r * N / (1 + r)
      if (f1 <= 3 || f2 <= 3) next
      if (r == 1 && N %% 2) next
      if (g7_zpower(q / sqrt(1 / (f1 - 3) + 1 / (f2 - 3)), 1, a, tl, 1) >= pw) return(N)
    }
  }
  # ratio 1: fractional optimum 105 is odd, equal groups give 106 = 53/53 (SUP-G7-Z-CORR2-APR-A01)
  r <- g7_calc("z_correlations_two_independent", "a_priori", list(q = 0.5, alpha = 0.05, power = 0.8))
  expect_equal(r$solved$n_total, 106)
  expect_equal(unname(r$group_sizes), c(53, 53))
  expect_equal(r$solved$n_total, frac(0.5, 0.05, 0.8, 1, "one"))
  # ratio 1.0001: split rule, 105 = 52/53 (SUP-G7-Z-CORR2-APR-A03)
  r <- g7_calc("z_correlations_two_independent", "a_priori", list(q = 0.5, alpha = 0.05, power = 0.8, allocation_ratio = 1.0001))
  expect_equal(r$solved$n_total, 105)
  expect_equal(unname(r$group_sizes), c(52, 53))
  expect_gte(r$achieved$power, 0.8)
  # ratio 0.5, two tails: 304 = 203/101 (SUP-G7-Z-CORR2-APR-A06)
  r <- g7_calc("z_correlations_two_independent", "a_priori",
               list(q = 0.4, alpha = 0.05, power = 0.9, allocation_ratio = 0.5), list(tails = "two"))
  expect_equal(c(r$solved$n_total, unname(r$group_sizes)), c(304, 203, 101))
  # the G*Power split misses the target (GPW-Z-CORR2-APR-B001-E01: 25 = 4/21, power 0.785): the engine
  # reports the smallest larger N whose split reaches the target (27 = 5/22) and records G*Power's design
  r <- g7_calc("z_correlations_two_independent", "a_priori", list(q = 2.5, alpha = 0.05, power = 0.8, allocation_ratio = 5))
  expect_identical(r$status$code, "ok")
  expect_equal(c(r$solved$n_total, unname(r$group_sizes)), c(27, 5, 22))
  expect_gte(r$achieved$power, 0.8)
  expect_true(isTRUE(r$diagnostics$adjusted_from_gpower_rule))
  gp <- r$diagnostics$gpower_allocation_rule
  expect_equal(c(gp$n_total, gp$n1, gp$n2), c(25, 4, 21))
  expect_lt(abs(gp$power - 0.7847884), 5e-8)
  # q = 0: power alpha for every N
  r <- g7_calc("z_correlations_two_independent", "a_priori", list(q = 0, alpha = 0.05, power = 0.8))
  expect_identical(r$status$code, "infeasible")
})

test_that("two independent r's: inverse modes recover the forward model", {
  r <- g7_calc("z_correlations_two_independent", "sensitivity", list(alpha = 0.05, power = 0.8, n1 = 100, n2 = 100))
  expect_identical(r$status$code, "ok")
  expect_gt(r$solved$q, 0)
  f <- g7_calc("z_correlations_two_independent", "post_hoc", list(q = r$solved$q, alpha = 0.05, n1 = 100, n2 = 100))
  expect_lt(abs(f$achieved$power - 0.8), 1e-9)
  # one tail, alpha above 0.5 (GPW-Z-CORR2-CRI-B001-D01: alpha 0.9309604, critical z -1.4829820)
  r <- g7_calc("z_correlations_two_independent", "criterion", list(q = 0.01, power = 0.9999, n1 = 1e5, n2 = 1e5))
  expect_lt(abs(r$solved$alpha - 0.9309604), 5e-8)
  expect_equal(r$quantities$critical, qnorm(1 - r$solved$alpha), tolerance = 1e-12)
  expect_lt(abs(r$quantities$critical + 1.4829820), 5e-8)
  r <- g7_calc("z_correlations_two_independent", "compromise", list(q = 0.3, error_ratio = 1, n1 = 100, n2 = 100),
               list(tails = "two"))
  expect_lt(abs(r$achieved$beta / r$solved$alpha - 1), 1e-8)
})

# Dependent correlations ---------------------------------------------------------------------

test_that("dependent r's: Steiger covariances equal the Isserlis/delta-method covariance", {
  set.seed(7)
  for (i in 1:6) {
    A <- matrix(rnorm(16), 4)
    S <- crossprod(A) + diag(4)
    R <- diag(1 / sqrt(diag(S))) %*% S %*% diag(1 / sqrt(diag(S)))
    expect_equal(.ig_g7_c_no_common(R[1, 2], R[1, 3], R[1, 4], R[2, 3], R[2, 4], R[3, 4]),
                 g7_cov_z(R, c(1, 2), c(3, 4)), tolerance = 1e-12)
    expect_equal(.ig_g7_c_common(R[1, 2], R[1, 3], R[2, 3]), g7_cov_z(R[1:3, 1:3], c(1, 2), c(1, 3)), tolerance = 1e-12)
  }
})

test_that("dependent r's: power, sign of the critical z and documented a priori results", {
  # manual 28.3.2 / GPW-Z-DEPCI-APR-B001-O01: N = 144
  r <- g7_calc("z_correlations_two_dependent_common_index", "a_priori",
               list(rho_ac = 0.2, alpha = 0.05, power = 0.8, rho_ab = 0.4, rho_bc = 0.5))
  expect_equal(r$solved$n_total, 144)
  expect_lt(r$quantities$critical, 0)   # rho_ac < rho_ab: critical z on the lower side
  R <- matrix(c(1, 0.4, 0.2, 0.4, 1, 0.5, 0.2, 0.5, 1), 3)
  R0 <- R; R0[1, 3] <- R0[3, 1] <- 0.4
  c0 <- g7_cov_z(R0, c(1, 2), c(1, 3)); c1 <- g7_cov_z(R, c(1, 2), c(1, 3))
  s0 <- sqrt((2 - 2 * c0) / 141)
  ref <- g7_zpower((atanh(0.2) - atanh(0.4)) / s0, sqrt((1 - c1) / (1 - c0)), 0.05, "one", -1)
  expect_lt(abs(r$achieved$power - ref), 1e-12)
  # manual 28.3.1: N = 886 (no common index), 1183 with zero dependencies
  r <- g7_calc("z_correlations_two_dependent_no_common_index", "a_priori",
               list(rho_cd = 0.2, alpha = 0.05, power = 0.8, rho_ab = 0.1, rho_ac = 0.5, rho_ad = 0.4, rho_bc = -0.4, rho_bd = 0.8))
  expect_equal(r$solved$n_total, 886)
  r <- g7_calc("z_correlations_two_dependent_no_common_index", "a_priori",
               list(rho_cd = 0.2, alpha = 0.05, power = 0.8, rho_ab = 0.1, rho_ac = 0, rho_ad = 0, rho_bc = 0, rho_bd = 0))
  expect_equal(r$solved$n_total, 1183)
})

test_that("dependent r's: sensitivity follows the effect direction and the PSD range", {
  # manual 28.3.3: <= branch 0.047702; >= branch has no PSD solution (G*Power shows 2)
  r <- g7_calc("z_correlations_two_dependent_common_index", "sensitivity",
               list(alpha = 0.05, power = 0.8, n_total = 144, rho_ab = 0.4, rho_bc = -0.6))
  expect_identical(r$status$code, "ok")
  expect_lt(abs(r$solved$rho_ac - 0.0477022), 5e-8)
  expect_lt(r$quantities$critical, 0)
  r <- g7_calc("z_correlations_two_dependent_common_index", "sensitivity",
               list(alpha = 0.05, power = 0.8, n_total = 144, rho_ab = 0.4, rho_bc = -0.6),
               list(effect_direction = "rho_ac_ge_rho_ab"))
  expect_identical(r$status$code, "infeasible")
  expect_null(r$solved$rho_ac)
  # root just inside the PSD bound 0.99373 (SUP-G7-Z-DEPCI-SEN-A01: 0.9877916)
  r <- g7_calc("z_correlations_two_dependent_common_index", "sensitivity",
               list(alpha = 0.05, power = 0.8, n_total = 5, rho_ab = 0.4, rho_bc = 0.5),
               list(tails = "two", effect_direction = "rho_ac_ge_rho_ab"))
  expect_identical(r$status$code, "ok")
  expect_lt(abs(r$solved$rho_ac - 0.9877916), 5e-8)
  seg <- r$method_outputs$admissible_segment
  expect_equal(seg[2], 0.2 + sqrt(0.84 * 0.75), tolerance = 1e-10)
  expect_identical(r$method_outputs$evaluated_at, r$solved$rho_ac)
  # SUP-G7-Z-DEPNCI-SEN-A01 / A02: no solution at N = 530, a solution at N = 540
  pr <- list(alpha = 0.05, power = 0.8, rho_ab = 0.1, rho_ac = 0.5, rho_ad = 0.4, rho_bc = -0.4, rho_bd = 0.8)
  r <- g7_calc("z_correlations_two_dependent_no_common_index", "sensitivity", c(pr, n_total = 530))
  expect_identical(r$status$code, "infeasible")
  r <- g7_calc("z_correlations_two_dependent_no_common_index", "sensitivity", c(pr, n_total = 540))
  expect_lt(abs(r$solved$rho_cd - 0.2275172), 5e-8)
})

# Tetrachoric correlation --------------------------------------------------------------------

test_that("tetrachoric: cell probabilities and the exact-mode standard error", {
  for (pp in list(c(0.2399846, 0.6019313, 0.5815451), c(-0.6, 0.3, 0.7), c(0.8, 0.05, 0.95), c(0.5, 0.2, 0.7))) {
    eng <- .ig_g7_tetra_cells(pp[1], pp[2], pp[3])
    ref <- g7_cells(pp[1], pp[2], pp[3])
    expect_equal(unname(eng$p), unname(ref[c("a", "b", "c", "d")]), tolerance = 1e-9)
    if (min(ref[c("a", "b", "c", "d")]) > 1e-6) {
      expect_equal(.ig_g7_tetra_se_exact(eng$p, eng$h, eng$k, pp[1], 100), g7_se_pearson(ref, pp[1], 100), tolerance = 1e-9)
    }
  }
  # DOC-229 (manual 31.3): se0 0.074465 at N = 463; the printed first term (a + d)(b + d)/4 would give 0.083673
  cl <- .ig_g7_tetra_cells(0, 0.6019313, 0.5815451)
  expect_lt(abs(.ig_g7_tetra_se_exact(cl$p, cl$h, cl$k, 0, 463) - 0.074465), 5e-7)
})

test_that("tetrachoric: power in both computation modes, critical r and the Bonett-Price interval", {
  r <- g7_calc("z_tetrachoric_correlation", "a_priori",
               list(rho_h1 = 0.2399846, rho_h0 = 0, marginal_x = 0.6019313, marginal_y = 0.5815451, alpha = 0.05, power = 0.95))
  expect_equal(r$solved$n_total, 463)    # DOC-229, GPW-Z-TETRA-APR-B001-O01
  expect_equal(r$method_outputs$critical_r, r$method_outputs$se_h0 * qnorm(0.95), tolerance = 1e-12)
  # approximation mode at symmetric margins: rho* = rho and the same SE as exact mode
  ex <- g7_calc("z_tetrachoric_correlation", "post_hoc",
                list(rho_h1 = 0.22, rho_h0 = 0.2, marginal_x = 0.5, marginal_y = 0.5, n_total = 1e5, alpha = 1e-4), list(tails = "two"))
  ap <- g7_calc("z_tetrachoric_correlation", "post_hoc",
                list(rho_h1 = 0.22, rho_h0 = 0.2, marginal_x = 0.5, marginal_y = 0.5, n_total = 1e5, alpha = 1e-4),
                list(tails = "two", computation_of_correlation = "approx_bonett_price"))
  expect_equal(ap$method_outputs$rho_h1_used, 0.22, tolerance = 1e-12)
  expect_equal(ap$method_outputs$se_h0, ex$method_outputs$se_h0, tolerance = 1e-10)
  expect_equal(ap$achieved$power, ex$achieved$power, tolerance = 1e-9)
  # two tails, approximation mode: Wald bounds and the Bonett-Price interval differ
  wb <- ap$method_outputs$critical_r
  bp <- ap$method_outputs$bonett_price_interval_h0
  expect_length(bp, 2)
  expect_gt(abs(bp[1] - wb[1]), 1e-6)
  expect_null(ex$method_outputs$bonett_price_interval_h0)
  # asymmetric margins (SUP-G7-Z-TETRA-POS-A02): power 0.5986893, r* 0.3983425 / 0.0956036
  ap <- g7_calc("z_tetrachoric_correlation", "post_hoc",
                list(rho_h1 = 0.4, rho_h0 = 0.1, marginal_x = 0.3, marginal_y = 0.7, n_total = 150, alpha = 0.05),
                list(tails = "two", computation_of_correlation = "approx_bonett_price"))
  expect_lt(abs(ap$achieved$power - 0.5986893), 5e-8)
  expect_lt(abs(ap$method_outputs$rho_h1_used - 0.3983425), 5e-8)
})

test_that("tetrachoric: sensitivity branches choose the requested crossing", {
  r <- g7_calc("z_tetrachoric_correlation", "sensitivity",
               list(alpha = 0.05, power = 0.8, n_total = 200, rho_h0 = 0, marginal_x = 0.5, marginal_y = 0.5),
               list(effect_direction = "rho_h1_le_h0"))
  expect_lt(abs(r$solved$rho_h1 + 0.2712855), 5e-8)   # GPW-Z-TETRA-SEN-B002-O01
  expect_lt(r$quantities$critical, 0)
  r <- g7_calc("z_tetrachoric_correlation", "sensitivity",
               list(alpha = 0.05, power = 0.8, n_total = 150, rho_h0 = 0.1, marginal_x = 0.3, marginal_y = 0.7),
               list(tails = "two", effect_direction = "rho_h1_le_h0"))
  expect_lt(abs(r$solved$rho_h1 + 0.2886073), 5e-8)   # SUP-G7-Z-TETRA-SEN-A02
  # One tail, upper branch, asymmetric margins: the power reaches the target at
  # rho 0.4402 (G*Power's answer, SUP-G7-Z-TETRA-SEN-A01) and falls below it
  # again near rho 1 (the H1 standard error grows without bound); the generic
  # solver selects the first increasing crossing on the requested half-interval.
  r <- g7_calc("z_tetrachoric_correlation", "sensitivity",
               list(alpha = 0.05, power = 0.8, n_total = 150, rho_h0 = 0.1, marginal_x = 0.3, marginal_y = 0.7))
  expect_identical(r$status$code, "ok")
  expect_lt(abs(r$solved$rho_h1 - 0.4402148), 5e-7)
  expect_gte(r$diagnostics$scan_grid$from, 0.1)
  expect_true(any(abs(r$diagnostics$roots_found - 0.4402148) < 5e-7))
})

test_that("tetrachoric: criterion keeps the sign convention of the critical z for alpha above 0.5", {
  r <- g7_calc("z_tetrachoric_correlation", "criterion",
               list(rho_h1 = 0.3, rho_h0 = 0, marginal_x = 0.5, marginal_y = 0.5, n_total = 10, power = 0.95))
  expect_lt(abs(r$solved$alpha - 0.8251929), 5e-8)   # SUP-G7-Z-TETRA-CRI-A01
  expect_lt(abs(r$quantities$critical + 0.9353378), 5e-8)
  # the engine's critical r is rho0 + z_c se0 (negative); G*Power displays +0.4646098 (|z_c|)
  expect_equal(r$method_outputs$critical_r, r$quantities$critical * r$method_outputs$se_h0, tolerance = 1e-12)
})

# Drawers ------------------------------------------------------------------------------------

test_that("G7 drawers reproduce the manual's conversions", {
  es <- effect_size("z_correlations_two_independent", "two_correlations", r1 = 0.75, r2 = 0.88)
  expect_lt(abs(es$value + 0.4028126), 5e-8)   # DOC-139
  es <- effect_size("z_correlations_two_dependent_common_index", "cohens_q", rho_ab = 0.4, q = 0.3)
  expect_lt(abs(es$value - 0.1230226), 5e-8)   # GPW-UT-DRAW-UTL-U071-UO01
  es <- effect_size("z_correlations_two_dependent_no_common_index", "cohens_q", rho_ab = 0.1, q = 1)
  expect_lt(abs(es$value + 0.7161345), 5e-8)   # DOC-149
  es <- effect_size("z_tetrachoric_correlation", "cell_probabilities", p11 = 0.25, p12 = 0.34, p22 = 0.2)
  expect_lt(abs(es$value + 0.1380485), 5e-8)   # DOC-235 (p21 filled as 1 - sum)
  d <- .ig_drawer_z_tetrachoric_correlation("ci_from_frequencies",
                                            list(f11 = 25, f12 = 34, f21 = 12, f22 = 25, position = 0, confidence = 0.95))
  expect_lt(abs(d$details$std_error_r - 0.1627707), 5e-8)   # DOC-236 / GPW-UT-DRAW-UTL-U069-UO03
  expect_lt(abs(d$value + 0.1631428), 5e-8)
  expect_identical(g7_code(effect_size("z_tetrachoric_correlation", "cell_probabilities",
                                       p11 = 0.25, p12 = 0.34, p21 = 0.25, p22 = 0.2)), "invalid_input")
  expect_identical(g7_code(effect_size("z_tetrachoric_correlation", "ci_from_frequencies", f11 = 1, f12 = 2, f21 = 3,
                                       f22 = 4, position = 0, confidence = 0.95, computation = "approx_bonett_price")),
                   "reference_pending")
})
