# Independent metric arithmetic and explicit type-7 interpolation. Only the
# documented R multinomial RNG is shared with the implementation under test.
reference_bootstrap <- function(data, design, replicates, seed) {
  ids <- c("tpr", "fpr", "ppv", "accuracy", "brier_score", "roc_auc")
  values <- matrix(NA_real_, replicates, length(ids), dimnames = list(NULL, ids))
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed)
  for (b in seq_len(replicates)) {
    w <- numeric(nrow(data))
    if (design == "cluster") {
      groups <- sort(unique(data$cluster))
      draw <- stats::rmultinom(1, length(groups), rep(1, length(groups)))[, 1]
      w <- draw[match(data$cluster, groups)] * data$w
    } else {
      groups <- if (design == "class") list(which(data$y == 1), which(data$y == 0)) else list(seq_len(nrow(data)))
      for (g in groups) w[g] <- stats::rmultinom(1, sum(data$w[g]), data$w[g])[, 1]
    }
    alert <- data$p >= .5
    tp <- sum(w * (data$y == 1 & alert)); fn <- sum(w * (data$y == 1 & !alert))
    fp <- sum(w * (data$y == 0 & alert)); tn <- sum(w * (data$y == 0 & !alert))
    pos <- which(data$y == 1); neg <- which(data$y == 0)
    numerator <- 0
    for (i in pos) for (j in neg)
      numerator <- numerator + w[i] * w[j] * ((data$p[i] > data$p[j]) + .5 * (data$p[i] == data$p[j]))
    values[b, ] <- c(tp/(tp+fn), fp/(fp+tn), tp/(tp+fp), (tp+tn)/sum(w),
      sum(w*(data$p-data$y)^2)/sum(w), numerator/(sum(w[pos])*sum(w[neg])))
  }
  percentile <- function(x, probability) {
    x <- sort(x[is.finite(x)])
    h <- 1 + (length(x) - 1) * probability
    lo <- floor(h); hi <- ceiling(h)
    x[lo] + (h - lo) * (x[hi] - x[lo])
  }
  lapply(ids, function(id) {
    x <- values[, id]; bad <- sum(!is.finite(x))
    list(metric_id = id, status = if (bad > .1 * replicates) "undefined" else "ok",
      lower = percentile(x, .025), upper = percentile(x, .975), excluded = bad)
  })
}

# Exact repeated-sampling coverage of the infinite-replicate percentile
# sensitivity interval for n independent event units. This is a diagnostic
# limit, not a claim about finite-B interpolation or all classifier metrics.
reference_binomial_coverage <- function(n, true_probability) {
  coefficients <- 1
  for (j in seq_len(n)) coefficients <- c(0, coefficients) + c(coefficients, 0)
  mass <- function(p) coefficients * p^(0:n) * (1-p)^(n:0)
  intervals <- vapply(0:n, function(x) {
    cdf <- cumsum(mass(x/n))
    c(which(cdf >= .025)[1]-1, which(cdf >= .975)[1]-1)/n
  }, numeric(2))
  sum(mass(true_probability) * (intervals[1, ] <= true_probability & intervals[2, ] >= true_probability))
}
