# Sampling uncertainty: binomial intervals against independent references,
# explicit restrictions, and the properties of the seeded bootstrap.

# Independent references (no qnorm/qbeta/pbinom): z from integrate(dnorm),
# Wilson limits as roots of the score-test equation, Clopper-Pearson limits
# by bisection on binomial tails summed from lchoose().
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
      w <- inferguidR:::.ig_cls_wilson(x, n, level)
      cp <- inferguidR:::.ig_cls_clopper_pearson(x, n, level)
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

test_that("interval methods can be switched off or restricted", {
  r <- evaluate_cutoff(classification_spec(counts = c(tp = 5, fn = 5, fp = 5, tn = 5),
                                           uncertainty = classification_uncertainty(proportion_methods = "wilson",
                                                                                    auc_method = "none")))
  iv <- classification_table(r, "intervals")
  expect_identical(unique(iv$method), "wilson")
  r0 <- evaluate_cutoff(classification_spec(counts = c(tp = 5, fn = 5, fp = 5, tn = 5),
                                            uncertainty = classification_uncertainty(proportion_methods = character(),
                                                                                     auc_method = "none")))
  expect_identical(nrow(classification_table(r0, "intervals")), 0L)
  expect_identical(th_condition_code(classification_uncertainty(level = 1)), "invalid_input")
  expect_identical(th_condition_code(classification_uncertainty(proportion_methods = "wald")), "invalid_input")
})

test_that("expected counts have no sampling interval", {
  sc <- update(classification_examples("scored_synthetic", models = "model_a")$spec, prevalence = 0.05)
  r <- evaluate_cutoff(sc)
  expect_true(is.null(r$models[[1]]$scenario$intervals))
  th <- evaluate_cutoff(classification_examples("binormal_equal_variance")$spec)
  iv <- classification_table(th, "intervals")
  expect_true(all(iv$status == "not_applicable"))
})

test_that("the bootstrap is reproducible from its seed and restores the caller's random-number state", {
  ex <- classification_examples("scored_synthetic", models = "model_a")
  unc <- classification_uncertainty(bootstrap_replicates = 200, seed = 42,
                                    bootstrap_metrics = c("tpr", "ppv", "roc_auc", "average_precision",
                                                          "brier_score", "prevalence"))
  spec <- update(ex$spec, uncertainty = unc)
  set.seed(1)
  before <- .Random.seed
  r1 <- evaluate_cutoff(spec)
  expect_identical(.Random.seed, before)
  r2 <- evaluate_cutoff(spec)
  b1 <- classification_table(r1, "intervals")
  b2 <- classification_table(r2, "intervals")
  expect_identical(b1, b2)
  boot <- b1[b1$method == "bootstrap_percentile", ]
  expect_setequal(boot$metric_id, c("tpr", "ppv", "roc_auc", "average_precision", "brier_score", "prevalence"))
  expect_true(all(boot$status == "ok"))
  expect_match(boot$note[1], "200 replicates \\(stratified by class\\), seed 42, percentile")
  # stratified by class: every replicate keeps the class totals, so the prevalence does not vary
  prev <- boot[boot$metric_id == "prevalence", ]
  expect_equal(prev$lower, 0.3)
  expect_equal(prev$upper, 0.3)
  # a different seed gives different limits
  r3 <- evaluate_cutoff(update(spec, uncertainty = classification_uncertainty(
    bootstrap_replicates = 200, seed = 43, bootstrap_metrics = "tpr")))
  b3 <- classification_table(r3, "intervals")
  expect_false(identical(b3$lower[b3$method == "bootstrap_percentile"], boot$lower[boot$metric_id == "tpr"]))
  # without a prior seed in the session none is left behind
  had <- exists(".Random.seed", envir = globalenv())
  saved <- if (had) get(".Random.seed", envir = globalenv())
  rm(".Random.seed", envir = globalenv())
  evaluate_cutoff(spec)
  expect_false(exists(".Random.seed", envir = globalenv()))
  if (had) assign(".Random.seed", saved, envir = globalenv())
})

test_that("bootstrap settings are validated and limited to scored data", {
  expect_identical(th_condition_code(classification_uncertainty(bootstrap_replicates = 200)), "invalid_input")
  expect_identical(th_condition_code(classification_uncertainty(bootstrap_replicates = 50, seed = 1)), "invalid_input")
  expect_identical(th_condition_code(classification_uncertainty(bootstrap_replicates = 200, seed = 1,
                                                                bootstrap_metrics = "lr_plus")), "invalid_input")
  unc <- classification_uncertainty(bootstrap_replicates = 200, seed = 1)
  expect_identical(th_condition_code(classification_spec(counts = c(tp = 1, fn = 1, fp = 1, tn = 1),
                                                         uncertainty = unc)), "invalid_input")
  expect_identical(th_condition_code(classification_spec(theory = list(auc = 0.8, sigma = 1), cutoff = 0,
                                                         prevalence = 0.1, uncertainty = unc)), "invalid_input")
})

test_that("the cluster bootstrap draws whole clusters", {
  set.seed(5)
  cl <- rep(1:6, times = c(1, 2, 3, 4, 2, 3))
  draws <- inferguidR:::.ig_cls_with_seed(3, replicate(50, inferguidR:::.ig_cls_bootstrap_draw(
    y = rep(0:1, length.out = 15), w = rep(1, 15), cluster_index = cl, n_clusters = 6L, strata = "class")))
  for (k in seq_len(ncol(draws))) {
    per_cluster <- tapply(draws[, k], cl, function(v) length(unique(v)))
    expect_true(all(per_cluster == 1))
    expect_equal(sum(tapply(draws[, k], cl, `[`, 1)), 6)
  }
})

test_that("unstable bootstrap quantities are reported as undefined or infinite, not as numbers", {
  # one alerted event among 12 units: in many replicates nobody is alerted, so PPV is undefined
  d <- data.frame(y = c(1, rep(0, 5), rep(1, 3), rep(0, 3)), s = c(10, 1:11))
  spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "s", cutoff = 10.5,
                              uncertainty = classification_uncertainty(bootstrap_replicates = 200, seed = 7,
                                                                       bootstrap_metrics = c("ppv", "lr_pos")))
  iv <- classification_table(evaluate_cutoff(spec), "intervals")
  b <- iv[iv$method == "bootstrap_percentile", ]
  expect_identical(b$status[b$metric_id == "ppv"], "undefined")
  expect_identical(b$reason_code[b$metric_id == "ppv"], "bootstrap_too_many_undefined")
  expect_true(all(is.na(b$lower[b$status != "ok"])))
})
