# Independent spot checks of the repeated-measures engines (group G3) and their
# drawers. The full check is validation/independent/g3/check_rm.R; here a fast
# subset uses the definition-based references of helper-reference.R (ref_pf:
# noncentral F by integrating the convolution noncentral chi-square over the
# denominator density; no pbeta/qbeta, no inferguidR routine; decision D18) and
# effects built from explicit populations.

g3i_qf_upper <- function(alpha, df1, df2) {
  # Illinois regula falsi on the reference central F tail (ref_pf with
  # lambda = 0); stats::qf() only places the starting bracket.
  h <- function(x) log(ref_pf(x, df1, df2, 0, lower_tail = FALSE)) - log(alpha)
  q0 <- stats::qf(alpha, df1, df2, lower.tail = FALSE)
  a <- q0 * (1 - 1e-3)
  b <- q0 * (1 + 1e-3)
  ha <- h(a)
  hb <- h(b)
  while (ha < 0) { a <- a / 2; ha <- h(a) }
  while (hb > 0) { b <- b * 2; hb <- h(b) }
  side <- 0
  for (i in 1:60) {
    x <- (a * hb - b * ha) / (hb - ha)
    hx <- h(x)
    if (hx > 0) {
      a <- x; ha <- hx
      if (side == 1) hb <- hb / 2
      side <- 1
    } else {
      b <- x; hb <- hx
      if (side == -1) ha <- ha / 2
      side <- -1
    }
    if (abs(b - a) <= 1e-13 * x || hx == 0) break
  }
  x
}

g3i_cs <- function(m, rho) {
  S <- matrix(rho, m, m)
  diag(S) <- 1
  S
}

g3i_contrasts <- function(m) {
  C <- sapply(seq_len(m - 1), function(j) {
    v <- c(rep(1, j), -j, rep(0, m - j - 1))
    v / sqrt(sum(v^2))
  })
  matrix(C, m)
}

test_that("post hoc power, beta and critical F agree with the integral reference", {
  cases <- list(
    list("f_rm_anova_within", list(f = 0.2, alpha = 0.05, n_total = 45, k = 3, m = 4, rho = 0.6, epsilon = 0.75), list()),
    list("f_rm_anova_between", list(f = 0.3, alpha = 0.01, n_total = 40, k = 4, m = 3, rho = 0.2),
         list(effect_size_specification = "gpower30")),
    list("f_rm_manova_within_between", list(f = 0.3, alpha = 0.05, n_total = 80, k = 3, m = 4),
         list(approximation = "wilks_u"))
  )
  for (cs in cases) {
    r <- calculate_power(analysis_spec(cs[[1]], "post_hoc", inputs = cs[[2]], options = cs[[3]]))
    df <- unname(r$quantities$df)
    crit <- g3i_qf_upper(cs[[2]]$alpha, df[1], df[2])
    expect_equal(r$quantities$critical, crit, tolerance = 1e-8, info = cs[[1]])
    expect_equal(r$achieved$power, ref_pf(crit, df[1], df[2], r$quantities$ncp, lower_tail = FALSE),
                 tolerance = 1e-7, info = cs[[1]])
    expect_equal(r$achieved$beta, ref_pf(crit, df[1], df[2], r$quantities$ncp, lower_tail = TRUE),
                 tolerance = 1e-7, info = cs[[1]])
  }
})

test_that("noncentralities equal the quadratic forms of explicit populations", {
  N <- 60; k <- 3; m <- 4; rho <- 0.35; f <- 0.27
  S <- g3i_cs(m, rho)
  C <- g3i_contrasts(m)
  z <- function(n) { v <- seq_len(n) - (n + 1) / 2; v / sqrt(mean(v^2)) }
  # within: measurement effects with variance f^2
  mu <- f * z(m)
  lam_w <- N * as.numeric(t(t(C) %*% mu) %*% solve(t(C) %*% S %*% C, t(C) %*% mu))
  r <- calculate_power(analysis_spec("f_rm_anova_within", "post_hoc",
                                     inputs = list(f = f, alpha = 0.05, n_total = N, k = k, m = m, rho = rho)))
  expect_equal(r$quantities$ncp, lam_w, tolerance = 1e-12)
  # between: group effects with variance f^2 acting on every measurement (test of the subject sums)
  a <- f * z(k)
  lam_b <- N * mean((a * m)^2) / sum(S)
  r <- calculate_power(analysis_spec("f_rm_anova_between", "post_hoc",
                                     inputs = list(f = f, alpha = 0.05, n_total = N, k = k, m = m, rho = rho)))
  expect_equal(r$quantities$ncp, lam_b, tolerance = 1e-12)
  # MANOVA within with s = 1: exact Hotelling noncentrality N theta and df (m - 1, N - k - m + 2)
  r <- calculate_power(analysis_spec("f_rm_manova_within", "post_hoc",
                                     inputs = list(f = f, alpha = 0.05, n_total = N, k = k, m = m, rho = rho)))
  expect_equal(r$quantities$ncp, lam_w, tolerance = 1e-12)
  expect_equal(unname(r$quantities$df), c(m - 1, N - k - m + 2))
})

test_that("a priori N is the smallest admissible design (linear scan with the reference)", {
  spec <- list(f = 0.3, alpha = 0.05, power = 0.85, k = 3, m = 3, rho = 0.4, epsilon = 0.9)
  r <- calculate_power(analysis_spec("f_rm_anova_within_between", "a_priori", inputs = spec))
  pw <- function(N) {
    eps <- 0.9
    df1 <- 4 * eps
    df2 <- (N - 3) * 2 * eps
    lam <- 0.09 * N * 3 * eps / 0.6
    ref_pf(g3i_qf_upper(0.05, df1, df2), df1, df2, lam, lower_tail = FALSE)
  }
  # power increases with N (lambda and df2 grow), so the solution and its
  # admissible predecessor decide minimality
  N <- r$solved$n_total
  expect_equal(N %% 3, 0)
  expect_gte(pw(N), 0.85)
  expect_lt(pw(N - 3), 0.85)
})

test_that("inverse modes meet their targets under the reference forward model", {
  r <- calculate_power(analysis_spec("f_rm_manova_within_between", "compromise",
                                     inputs = list(f = 0.02, error_ratio = 0.01, n_total = 100000, k = 10, m = 10),
                                     options = list(approximation = "hotelling_lawley_t1")))
  df <- unname(r$quantities$df)
  crit <- g3i_qf_upper(r$solved$alpha, df[1], df[2])
  beta <- ref_pf(crit, df[1], df[2], r$quantities$ncp, lower_tail = TRUE)
  # G*Power (GPW-F-RMMWB-COM-B006-D01--r01) shows beta = 0 and alpha 1.05e-15; here beta = q alpha
  expect_equal(beta / (0.01 * r$solved$alpha), 1, tolerance = 1e-6)
  expect_gt(r$solved$alpha, 1e-12)
})

test_that("drawers: direct, from variances, from means and conversions", {
  es <- effect_size("f_rm_anova_within", "direct", partial_eta2 = 0.06)
  expect_equal(es$value, sqrt(0.06 / 0.94))
  expect_true(abs(es$value - 0.2526456) <= 5e-8) # GPW-UT-DRAW-UTL-U022-UO01--r01
  dr <- .ig_drawer_f_rm_manova_between("from_variances", list(variance_effect = 0.0625, variance_within = 1))
  expect_equal(dr$value, 0.25)
  expect_true(abs(dr$partial_eta2 - 0.05882353) <= 5e-9) # GPW-UT-DRAW-UTL-U093-UO01--r01
  es <- effect_size("f_rm_anova_between", "from_means", means = c(1, 2, 4), sd = 2)
  expect_equal(es$value, sqrt(mean((c(1, 2, 4) - 7 / 3)^2)) / 2)
  # f (GPower 3.0) -> f(V): the f(V) of GPW-F-RMBET-SEN-B002-O01 from the f of -B001-O01
  es <- effect_size("f_rm_anova_between", "convert", f = 0.3359817, from = "gpower30", to = "cohen1988",
                    measurements = 3, rho = 0.5)
  expect_true(abs(es$value - 0.4114918) <= 1e-7)
  es <- effect_size("f_rm_manova_within", "convert", f = 0.1680711, from = "on", to = "off", measurements = 3,
                    rho = 0.5)
  expect_true(abs(es$value - 0.4116885) <= 1e-7) # SUP-G3-RMMW-SEN-R05 vs -SEN-U04
  expect_error(effect_size("f_rm_anova_within", "from_means", means = 1:3, sd = 1), class = "inferguidR_invalid_input")
  expect_error(effect_size("f_rm_anova_within", "direct", partial_eta2 = 1), class = "inferguidR_invalid_input")
})

test_that("drawers: MANOVA criterion and means/covariance", {
  # live drawer capture (inventory D2): Pillai V 0.4, k = 2, m = 4 -> f(V) 0.8164966
  es <- effect_size("f_rm_manova_within_between", "from_criterion", criterion = 0.4, approximation = "pillai_v",
                    groups = 2, measurements = 4)
  expect_true(abs(es$value - 0.8164966) <= 5e-8)
  M <- rbind(c(0, 0.2, 0.5), c(0, 0.1, 0.1), c(0, -0.1, 0.2))
  S <- g3i_cs(3, 0.5) * 4
  dr <- .ig_drawer_f_rm_manova_within_between("from_means_cov", list(means = M, sigma = S, approximation = "pillai_v"))
  C <- g3i_contrasts(3)
  B <- Reduce(`+`, lapply(1:3, function(g) tcrossprod(M[g, ] - colMeans(M)) / 3))
  th <- Re(eigen(solve(t(C) %*% S %*% C, t(C) %*% B %*% C), only.values = TRUE)$values)
  expect_equal(dr$criterion, sum(th / (1 + th)), tolerance = 1e-12)
  back <- effect_size("f_rm_manova_within_between", "from_criterion", criterion = dr$criterion,
                      approximation = "pillai_v", groups = 3, measurements = 3)
  expect_equal(back$value, dr$value, tolerance = 1e-12)
  expect_error(.ig_drawer_f_rm_manova_within_between("from_means_cov",
                                                      list(means = M, sigma = matrix(1, 3, 3), approximation = "pillai_v")),
               class = "inferguidR_invalid_input")
  expect_error(effect_size("f_rm_manova_within_between", "from_criterion", criterion = 1.2, approximation = "pillai_v",
                           groups = 2, measurements = 4), class = "inferguidR_invalid_input")
})
