# Numerical regression blocks retained from the historical guidance suite.
un_z <- function(level) {
  half <- function(z) integrate(dnorm, 0, z, rel.tol = 1e-13, abs.tol = 0)$value
  uniroot(function(z) 2 * half(z) - level, c(0.01, 10), tol = 1e-15)$root
}
un_wilson <- function(x, n, z) {
  ph <- x / n
  f <- function(p) (ph - p)^2 - z^2 * p * (1 - p) / n
  c(if (x == 0) 0 else uniroot(f, c(0, if (x == n) 1 - 1e-12 else ph), tol = 1e-15)$root,
    if (x == n) 1 else uniroot(f, c(if (x == 0) 1e-12 else ph, 1), tol = 1e-15)$root)
}
un_tail <- function(n, p, k) {
  v <- lchoose(n, k) + k * log(p) + (n - k) * log1p(-p)
  m <- max(v)
  exp(m + log(sum(exp(v - m))))
}
un_bisect <- function(f, lo, hi) {
  flo <- f(lo)
  for (i in 1:200) {
    mid <- (lo + hi) / 2
    if (mid == lo || mid == hi) break
    fm <- f(mid)
    if (sign(fm) == sign(flo)) { lo <- mid; flo <- fm } else hi <- mid
  }
  (lo + hi) / 2
}
un_cp <- function(x, n, level) {
  a2 <- (1 - level) / 2
  c(if (x == 0) 0 else un_bisect(function(p) un_tail(n, p, x:n) - a2, 1e-300, 1 - 1e-16),
    if (x == n) 1 else un_bisect(function(p) un_tail(n, p, 0:x) - a2, 1e-300, 1 - 1e-16))
}
test_that("Wilson and Clopper-Pearson limits equal independent references", {
  for (level in c(0.9, 0.95)) {
    z <- un_z(level)
    for (n in c(1, 2, 7, 25)) for (x in 0:n) {
      w <- infercalcR:::.ig_cls_wilson(x, n, level)
      cp <- infercalcR:::.ig_cls_clopper_pearson(x, n, level)
      expect_equal(c(w$lower, w$upper), un_wilson(x, n, z), tolerance = 1e-12, info = paste(level, x, n))
      expect_equal(c(cp$lower, cp$upper), un_cp(x, n, level), tolerance = 1e-12, info = paste(level, x, n))
    }
  }
  # the reported intervals use the specification's level and the right numerator/denominator
  r <- evaluate_cutoff(classification_spec(counts = c(tp = 17, fn = 8, fp = 12, tn = 63),
                                           uncertainty = classification_uncertainty(level = 0.9)))
  iv <- classification_table(r, "intervals")
  ppv <- iv[iv$metric_id == "ppv" & iv$method == "clopper_pearson", ]
  expect_equal(c(ppv$lower, ppv$upper), un_cp(17, 29, 0.9), tolerance = 1e-12)
  npv <- iv[iv$metric_id == "npv" & iv$method == "wilson", ]
  expect_equal(c(npv$lower, npv$upper), un_wilson(63, 71, un_z(0.9)), tolerance = 1e-12)
  expect_true(all(iv$level == 0.9))
  expect_true(all(iv$conditioning[iv$method != "delong"] == "fixed_cutoff"))
})
test_that("the cluster bootstrap draws whole clusters", {
  set.seed(5)
  cl <- rep(1:6, times = c(1, 2, 3, 4, 2, 3))
  draws <- infercalcR:::.ig_cls_with_seed(3, replicate(50, infercalcR:::.ig_cls_bootstrap_draw(
    y = rep(0:1, length.out = 15), w = rep(1, 15), cluster_index = cl, n_clusters = 6L, strata = "class")))
  for (k in seq_len(ncol(draws))) {
    per_cluster <- tapply(draws[, k], cl, function(v) length(unique(v)))
    expect_true(all(per_cluster == 1))
    expect_equal(sum(tapply(draws[, k], cl, `[`, 1)), 6)
  }
})
