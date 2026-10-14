# Independent checks of the classifier engine.
#
# Everything defined in this file is written from the definitions with plain
# loops and arithmetic; it calls no inferguidR function except the engine
# entry points under test. Counts come from comparing every unit with every
# candidate cutoff, AUC from all event/non-event pairs, partial AUC from a
# dense numerical integration of the ROC polygon, DeLong from the explicit
# pairwise kernel.

# Reference implementation --------------------------------------------------------

bf_alert <- function(s, comparator, cutoff) {
  switch(comparator, ">=" = s >= cutoff, ">" = s > cutoff, "<=" = s <= cutoff, "<" = s < cutoff)
}

bf_counts <- function(score, y, w, comparator, cutoff) {
  tp <- fn <- fp <- tn <- 0
  for (i in seq_along(score)) {
    a <- bf_alert(score[i], comparator, cutoff)
    if (y[i] == 1) {
      if (a) tp <- tp + w[i] else fn <- fn + w[i]
    } else {
      if (a) fp <- fp + w[i] else tn <- tn + w[i]
    }
  }
  c(tp = tp, fn = fn, fp = fp, tn = tn)
}

# Distinct decision sets: try every original score and one value beyond each
# end as the cutoff; distinct (TP, FP) pairs are the decision sets.
bf_sets <- function(score, y, w, comparator) {
  keep <- w > 0
  cand <- c(sort(unique(score[keep])), min(score) - 1, max(score) + 1)
  rows <- t(vapply(cand, function(cc) bf_counts(score, y, w, comparator, cc), numeric(4)))
  rows <- unique(rows)
  rows <- rows[order(rows[, "tp"] + rows[, "fp"]), , drop = FALSE]
  as.data.frame(rows)
}

bf_ratio <- function(num, den) {
  if (den == 0) {
    if (num > 0) list(value = NA_real_, status = "positive_infinity") else list(value = NA_real_, status = "undefined")
  } else list(value = num / den, status = "ok")
}

bf_metrics <- function(tp, fn, fp, tn) {
  P <- tp + fn; N <- fp + tn; A <- tp + fp; Q <- fn + tn; n <- P + N
  r <- list()
  r$tpr <- bf_ratio(tp, P); r$fnr <- bf_ratio(fn, P); r$tnr <- bf_ratio(tn, N); r$fpr <- bf_ratio(fp, N)
  r$ppv <- bf_ratio(tp, A); r$fdr <- bf_ratio(fp, A); r$npv <- bf_ratio(tn, Q); r$for_rate <- bf_ratio(fn, Q)
  r$prevalence <- bf_ratio(P, n); r$alert_rate <- bf_ratio(A, n); r$accuracy <- bf_ratio(tp + tn, n)
  r$tp_share <- bf_ratio(tp, n); r$fp_share <- bf_ratio(fp, n)
  r$fp_per_1000_non_events <- if (N == 0) list(value = NA_real_, status = "undefined") else list(value = 1000 * fp / N, status = "ok")
  r$missed_per_1000_total <- if (n == 0) list(value = NA_real_, status = "undefined") else list(value = 1000 * fn / n, status = "ok")
  both <- P > 0 && N > 0
  r$youden_j <- if (both) list(value = tp / P + tn / N - 1, status = "ok") else list(value = NA_real_, status = "undefined")
  r$balanced_accuracy <- if (both) list(value = (tp / P + tn / N) / 2, status = "ok") else list(value = NA_real_, status = "undefined")
  r$markedness <- if (A > 0 && Q > 0) list(value = tp / A + tn / Q - 1, status = "ok") else list(value = NA_real_, status = "undefined")
  fb <- function(b) {
    den <- (1 + b^2) * tp + b^2 * fn + fp
    if (den == 0) list(value = NA_real_, status = "undefined") else list(value = (1 + b^2) * tp / den, status = "ok")
  }
  r$f1 <- fb(1); r$f05 <- fb(0.5); r$f2 <- fb(2)
  r$jaccard <- if (tp + fn + fp == 0) list(value = NA_real_, status = "undefined") else list(value = tp / (tp + fn + fp), status = "ok")
  r$fowlkes_mallows <- if (A > 0 && P > 0) list(value = sqrt((tp / A) * (tp / P)), status = "ok") else list(value = NA_real_, status = "undefined")
  r$mcc <- if (A > 0 && P > 0 && N > 0 && Q > 0) list(value = (tp * tn - fp * fn) / sqrt(A * P * N * Q), status = "ok") else
    list(value = NA_real_, status = "undefined")
  # Cohen's kappa in its textbook form (observed and chance agreement)
  if (n > 0) {
    po <- (tp + tn) / n
    pe <- (A / n) * (P / n) + (Q / n) * (N / n)
    r$cohen_kappa <- if (pe == 1) list(value = NA_real_, status = "undefined") else list(value = (po - pe) / (1 - pe), status = "ok")
  }
  if (both) {
    r$lr_pos <- bf_ratio(tp / P, fp / N)
    r$lr_neg <- bf_ratio(fn / P, tn / N)
    r$dor <- bf_ratio(tp * tn, fp * fn)
  } else {
    r$lr_pos <- r$lr_neg <- r$dor <- list(value = NA_real_, status = "undefined")
  }
  r$alerts_per_tp <- if (A == 0) list(value = NA_real_, status = "undefined") else bf_ratio(A, tp)
  r$fp_per_tp <- if (A == 0) list(value = NA_real_, status = "undefined") else bf_ratio(fp, tp)
  r
}

bf_auc <- function(score, y, w) {
  num <- 0
  for (i in which(y == 1)) for (j in which(y == 0)) {
    k <- if (score[i] > score[j]) 1 else if (score[i] == score[j]) 0.5 else 0
    num <- num + w[i] * w[j] * k
  }
  num / (sum(w[y == 1]) * sum(w[y == 0]))
}

# AP and trapezoidal PR area from decision sets (ordered by alerts).
bf_pr <- function(sets) {
  P <- sets$tp[1] + sets$fn[1]
  ap <- 0
  trap <- 0
  prev_r <- 0
  prev_p <- 1
  for (k in 2:nrow(sets)) {
    r <- sets$tp[k] / P
    p <- sets$tp[k] / (sets$tp[k] + sets$fp[k])
    if (r > prev_r) {
      ap <- ap + (r - prev_r) * p
      trap <- trap + (r - prev_r) * (p + prev_p) / 2
    }
    prev_r <- r
    prev_p <- p
  }
  c(ap = ap, trap = trap)
}

# Partial AUC by dense midpoint integration of the ROC polygon over x in
# [lo, hi] (x = FPR, y = TPR, or x = TPR, y = 1 - FPR). Points are in curve
# order (alert none to alert all); where several points share an x value the
# curve enters at the first and leaves at the last of them. The midpoint rule
# is exact on linear pieces; each jump of size j adds at most h * j, and the
# jumps sum to at most 1, so the error is below h = (hi - lo) / n.
bf_partial_n <- 1e6
bf_partial <- function(x, y, lo, hi, n = bf_partial_n) {
  ux <- unique(x)                      # x is non-decreasing along the curve
  y_in <- vapply(ux, function(v) y[which(x == v)[1]], 0)
  y_out <- vapply(ux, function(v) { i <- which(x == v); y[i[length(i)]] }, 0)
  h <- (hi - lo) / n
  mids <- lo + (seq_len(n) - 0.5) * h
  k <- findInterval(mids, ux)
  val <- y_out[k] + (y_in[k + 1] - y_out[k]) * (mids - ux[k]) / (ux[k + 1] - ux[k])
  sum(val) * h
}

bf_delong <- function(score, y) {
  xe <- score[y == 1]; xn <- score[y == 0]
  psi <- function(a, b) if (a > b) 1 else if (a == b) 0.5 else 0
  v10 <- vapply(xe, function(a) mean(vapply(xn, function(b) psi(a, b), 0)), 0)
  v01 <- vapply(xn, function(b) mean(vapply(xe, function(a) psi(a, b), 0)), 0)
  auc <- mean(v10)
  list(auc = auc, var = var(v10) / length(xe) + var(v01) / length(xn), v10 = v10, v01 = v01)
}

# Synthetic data with ties.
bf_data <- function(seed, n = 60, weights = FALSE) {
  set.seed(seed)
  y <- rbinom(n, 1, 0.35)
  if (all(y == y[1])) y[1:2] <- c(0, 1)
  s <- round(rnorm(n, mean = 0.9 * y, sd = 1), 1)
  w <- if (weights) sample(0:3, n, replace = TRUE, prob = c(0.1, 0.4, 0.3, 0.2)) else rep(1, n)
  data.frame(y = y, s = s, w = w)
}

# Tests --------------------------------------------------------------------------------

test_that("decision sets equal an independent enumeration for every comparator and direction", {
  for (seed in 1:6) for (cmp in c(">=", ">", "<=", "<")) for (wt in c(FALSE, TRUE)) {
    d <- bf_data(seed, weights = wt)
    spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "s", cutoff = 0,
                                comparator = cmp, weights = if (wt) "w" else NULL)
    sw <- cutoff_sweep(spec, metrics = character())
    bf <- bf_sets(d$s, d$y, d$w, cmp)
    expect_equal(nrow(sw), nrow(bf), info = paste(seed, cmp, wt))
    expect_equal(sw$tp, bf$tp, info = paste(seed, cmp, wt))
    expect_equal(sw$fp, bf$fp, info = paste(seed, cmp, wt))
    expect_equal(sw$fn, bf$fn)
    expect_equal(sw$tn, bf$tn)
    # each reported cutoff is an original score that reproduces its set
    has <- !is.na(sw$cutoff)
    expect_true(all(sw$cutoff[has] %in% d$s))
    for (k in which(has)) {
      expect_equal(unname(bf_counts(d$s, d$y, d$w, cmp, sw$cutoff[k])), c(sw$tp[k], sw$fn[k], sw$fp[k], sw$tn[k]))
    }
    # the unrepresentable endpoint is labelled, never given a fake cutoff
    if (cmp %in% c(">=", "<=")) expect_identical(sw$cutoff_kind[1], "alert_none") else
      expect_identical(sw$cutoff_kind[nrow(sw)], "alert_all")
  }
})

test_that("every operating-point metric equals independent arithmetic, including undefined and infinite cases", {
  tables <- list(c(40, 10, 20, 130), c(0, 5, 0, 95), c(5, 0, 95, 0), c(7, 3, 0, 10), c(0, 4, 6, 0),
                 c(3, 0, 0, 7), c(0, 0, 4, 6), c(4, 6, 0, 0), c(2, 3, 5, 7), c(1, 0, 0, 0), c(12.5, 3.25, 0.5, 80.75))
  for (tb in tables) {
    kind <- if (any(tb != round(tb))) "expected" else "observed"
    m <- classification_metrics(c(tp = tb[1], fn = tb[2], fp = tb[3], tn = tb[4]), count_kind = kind)
    ref <- bf_metrics(tb[1], tb[2], tb[3], tb[4])
    for (id in names(ref)) {
      row <- m[m$metric_id == id, ]
      expect_identical(row$status, ref[[id]]$status, info = sprintf("%s at %s", id, paste(tb, collapse = "/")))
      if (ref[[id]]$status == "ok") {
        expect_equal(row$value, ref[[id]]$value, tolerance = 1e-12,
                     info = sprintf("%s at %s", id, paste(tb, collapse = "/")))
      } else {
        expect_true(is.na(row$value))
        expect_false(is.na(row$reason_code))
      }
    }
  }
})

test_that("ROC AUC with ties equals the pairwise Mann-Whitney definition in both directions", {
  for (seed in 1:5) for (wt in c(FALSE, TRUE)) {
    d <- bf_data(seed, weights = wt)
    for (dir in c("higher_is_positive", "lower_is_positive")) {
      spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "s", cutoff = 0,
                                  direction = dir, weights = if (wt) "w" else NULL)
      a <- classification_auc(spec)
      o <- if (dir == "higher_is_positive") d$s else -d$s
      expect_equal(a$value[a$metric_id == "roc_auc"], bf_auc(o, d$y, d$w), tolerance = 1e-12)
    }
  }
})

test_that("the AUC equals the trapezoidal area under the enumerated ROC polygon", {
  d <- bf_data(11)
  spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "s", cutoff = 0)
  bf <- bf_sets(d$s, d$y, d$w, ">=")
  fpr <- bf$fp / (bf$fp[1] + bf$tn[1]); tpr <- bf$tp / (bf$tp[1] + bf$fn[1])
  trap <- sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2)
  expect_equal(classification_auc(spec)$value[1], trap, tolerance = 1e-12)
})

test_that("partial AUC equals a dense numerical integration of the ROC polygon", {
  for (seed in c(2, 7, 9)) {
    d <- bf_data(seed)
    spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "s", cutoff = 0)
    bf <- bf_sets(d$s, d$y, d$w, ">=")
    fpr <- bf$fp / (bf$fp[1] + bf$tn[1]); tpr <- bf$tp / (bf$tp[1] + bf$fn[1])
    for (rg in list(c(0, 0.2), c(0.1, 0.35), c(0, 1))) {
      bound <- 1.01 * diff(rg) / bf_partial_n
      a <- classification_auc(spec, partial_auc = list(focus = "fpr", range = rg))
      raw <- a$value[a$metric_id == "partial_auc" & a$variant == "raw"]
      expect_lt(abs(raw - bf_partial(fpr, tpr, rg[1], rg[2])), bound)
      expect_equal(a$value[a$variant %in% "normalised"], raw / diff(rg), tolerance = 1e-12)
      a2 <- classification_auc(spec, partial_auc = list(focus = "tpr", range = rg))
      raw2 <- a2$value[a2$metric_id == "partial_auc" & a2$variant == "raw"]
      expect_lt(abs(raw2 - bf_partial(tpr, 1 - fpr, rg[1], rg[2])), bound)
    }
    full <- classification_auc(spec, partial_auc = list(focus = "fpr", range = c(0, 1)))
    expect_equal(full$value[full$variant %in% "raw"], full$value[full$metric_id == "roc_auc"], tolerance = 1e-12)
  }
})

test_that("average precision and the trapezoidal PR area are distinct and equal their definitions", {
  for (seed in 1:5) {
    d <- bf_data(seed)
    spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "s", cutoff = 0)
    pr <- pr_summary(spec)
    ref <- bf_pr(bf_sets(d$s, d$y, d$w, ">="))
    expect_equal(pr$value[pr$metric_id == "average_precision"], unname(ref["ap"]), tolerance = 1e-12)
    expect_equal(pr$value[pr$metric_id == "pr_area_trapezoid"], unname(ref["trap"]), tolerance = 1e-12)
  }
  # with ties the two areas differ
  d <- bf_data(3)
  pr <- pr_summary(classification_spec(d, truth = "y", positive_label = 1, scores = "s", cutoff = 0))
  expect_false(isTRUE(all.equal(pr$value[1], pr$value[2])))
})

test_that("the DeLong variance equals the explicit pairwise-kernel computation", {
  for (seed in 1:4) {
    d <- bf_data(seed)
    spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "s", cutoff = 0)
    a <- classification_auc(spec, level = 0.9)
    ref <- bf_delong(d$s, d$y)
    z <- qnorm(0.95)
    expect_equal(a$ci_lower[1], max(0, ref$auc - z * sqrt(ref$var)), tolerance = 1e-12)
    expect_equal(a$ci_upper[1], min(1, ref$auc + z * sqrt(ref$var)), tolerance = 1e-12)
  }
})

test_that("frequency weights act exactly as replicated units", {
  d <- bf_data(5, weights = TRUE)
  long <- d[rep(seq_len(nrow(d)), d$w), c("y", "s")]
  sw <- classification_spec(d, truth = "y", positive_label = 1, scores = "s", cutoff = 0.4, weights = "w")
  sl <- classification_spec(long, truth = "y", positive_label = 1, scores = "s", cutoff = 0.4)
  rw <- evaluate_cutoff(sw)
  rl <- evaluate_cutoff(sl)
  mw <- classification_table(rw, "metrics")
  ml <- classification_table(rl, "metrics")
  expect_equal(mw$value, ml$value)
  expect_identical(unique(mw$estimate_kind), "frequency_weighted")
  sw_auc <- classification_table(rw, "summaries")
  sl_auc <- classification_table(rl, "summaries")
  expect_equal(sw_auc$value[1:4], sl_auc$value[1:4])
  iw <- classification_table(rw, "intervals")
  il <- classification_table(rl, "intervals")
  expect_equal(iw$lower, il$lower)
  expect_equal(iw$upper, il$upper)
})

test_that("the paired DeLong comparison equals the explicit kernel covariance", {
  set.seed(42)
  n <- 50
  y <- rbinom(n, 1, 0.4)
  a <- round(rnorm(n, y), 1)
  b <- round(a + rnorm(n, 0, 0.8), 1)
  d <- data.frame(y = y, a = a, b = b)
  spec <- classification_spec(d, truth = "y", positive_label = 1, scores = c("a", "b"), cutoff = 0)
  cmp <- compare_classifiers(spec)
  ra <- bf_delong(a, y); rb <- bf_delong(b, y)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  s10 <- cov(cbind(ra$v10, rb$v10)); s01 <- cov(cbind(ra$v01, rb$v01))
  S <- s10 / n1 + s01 / n0
  v <- S[1, 1] + S[2, 2] - 2 * S[1, 2]
  expect_equal(cmp$difference, ra$auc - rb$auc, tolerance = 1e-12)
  expect_equal(cmp$se, sqrt(v), tolerance = 1e-12)
  expect_equal(cmp$p_value, 2 * pnorm(-abs((ra$auc - rb$auc) / sqrt(v))), tolerance = 1e-12)
})

test_that("exhaustive selection equals an independent search over all decision sets", {
  for (seed in 1:6) {
    d <- bf_data(seed, n = 80)
    spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "s", cutoff = 0)
    bf <- bf_sets(d$s, d$y, d$w, ">=")
    vals <- lapply(seq_len(nrow(bf)), function(k) bf_metrics(bf$tp[k], bf$fn[k], bf$fp[k], bf$tn[k]))
    for (obj in c("max_youden", "max_f1", "max_mcc")) {
      id <- c(max_youden = "youden_j", max_f1 = "f1", max_mcc = "mcc")[[obj]]
      v <- vapply(vals, function(r) if (r[[id]]$status == "ok") r[[id]]$value else NA_real_, 0)
      best <- max(v, na.rm = TRUE)
      tied <- which(!is.na(v) & abs(v - best) <= 1e-12 * abs(best))
      res <- select_cutoff(spec, obj)
      op <- res$models[[1]]$operating_point
      expect_equal(c(op$tp, op$fp), c(bf$tp[tied[1]], bf$fp[tied[1]]), info = paste(seed, obj))
      expect_equal(res$models[[1]]$selection$n_tied, length(tied))
      res2 <- select_cutoff(spec, obj, tie_break = "most_alerts")
      op2 <- res2$models[[1]]$operating_point
      expect_equal(c(op2$tp, op2$fp), c(bf$tp[tied[length(tied)]], bf$fp[tied[length(tied)]]))
    }
    # highest TPR subject to FPR <= 0.2
    fpr <- vapply(vals, function(r) r$fpr$value, 0)
    tpr <- vapply(vals, function(r) r$tpr$value, 0)
    ok <- which(fpr <= 0.2)
    best <- ok[tpr[ok] == max(tpr[ok])][1]
    res <- select_cutoff(spec, "constraint", constraints = list(fpr_max = 0.2))
    op <- res$models[[1]]$operating_point
    expect_equal(c(op$tp, op$fp), c(bf$tp[best], bf$fp[best]))
  }
})
