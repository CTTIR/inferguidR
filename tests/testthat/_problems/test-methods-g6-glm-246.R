# Extracted from test-methods-g6-glm.R:246

# prequel ----------------------------------------------------------------------
lg <- "z_logistic_regression"
ps <- "z_poisson_regression"
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

# test -------------------------------------------------------------------------
expect_identical(.ig_glm_grid_size(30925L, 0.5), 15463)
expect_identical(.ig_glm_grid_size(705L, 0.9), 70)
expect_identical(.ig_glm_grid_size(20L, 0.9), 2)
expect_identical(.ig_glm_grid(list(type = "uniform", a = 0, b = 1), 3L), ((1:3) - 0.375) / 3.25)
x <- .ig_glm_grid(list(type = "binomial", pi = 0.5), 615L)
expect_identical(sum(x), 308)
xp <- .ig_glm_grid(list(type = "poisson", lambda = 1), 200L)
q <- ((1:200) - 0.375) / 200.25
expect_true(all(stats::ppois(xp, 1) > q & (xp == 0 | stats::ppois(xp - 1, 1) <= q)))
b0 <- 0
b1 <- log(1.5)
ng <- 200
xg <- stats::qnorm(((1:ng) - 0.375) / (ng + 0.25))
pr <- 1 / (1 + exp(-(b0 + b1 * xg)))
w <- pr * (1 - pr)
I <- matrix(c(sum(w), sum(w * xg), sum(w * xg), sum(w * xg^2)), 2)
lam_ref <- b1^2 / solve(I)[2, 2]
r <- calculate_power(g6_spec(lg, "post_hoc", list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, n_total = 200,
                                                    enumeration_threshold = 100000L),
                               list(computation = "enumeration", test_statistic = "wald", tails = "one")))
expect_true(g6_close(r$quantities$ncp, lam_ref, rel = 1e-12, abs = 0))
expect_true(abs(r$quantities$ncp - 7.2446808) <= 5e-8)
expect_true(g6_close(r$quantities$critical, g6r_zq(0.05)^2, rel = 1e-10, abs = 0))
rl <- sqrt(lam_ref)
z <- g6r_zq(0.05)
expect_true(g6_close(r$achieved$power, g6r_Phi(rl - z) + g6r_Phi(-rl - z)))
expect_identical(r$quantities$df, 1)
pb <- mean(pr)
lr_ref <- 2 * (sum(pr * log(pr) + (1 - pr) * log(1 - pr)) - ng * (pb * log(pb) + (1 - pb) * log(1 - pb)))
r <- calculate_power(g6_spec(lg, "post_hoc", list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, n_total = 200,
                                                    enumeration_threshold = 100000L),
                               list(computation = "enumeration", test_statistic = "likelihood_ratio")))
expect_true(g6_close(r$quantities$ncp, lr_ref, rel = 1e-10, abs = 0))
r <- calculate_power(g6_spec(lg, "post_hoc", list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, n_total = 200,
                                                    enumeration_threshold = 100L),
                               list(computation = "enumeration", test_statistic = "wald")))
expect_identical(r$method_outputs$procedure, "demidenko")
expect_true(abs(r$achieved$power - 0.8696578) <= 5e-8)
r <- calculate_power(g6_spec(ps, "post_hoc", list(exp_beta1 = 2, base_rate = 1, alpha = 0.05, n_total = 2, x_pi = 0.5,
                                                    enumeration_threshold = 100000L),
                               list(computation = "enumeration", x_distribution = "binomial", tails = "two")))
expect_true(g6_close(r$quantities$ncp, 2 * log(2)^2 / 3, rel = 1e-14, abs = 0))
