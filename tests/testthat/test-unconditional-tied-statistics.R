test_that("pooled criterion boundaries retain whole rational equality classes", {
  # All six tables have Z^2 = 10/3 exactly. Separately rounded proportions
  # and standard errors used to split this class into four double values.
  x1 <- c(2, 0, 6, 2, 8, 6)
  x2 <- c(0, 4, 4, 8, 8, 12)
  z <- .ig_unc_stat("z_pooled", x1, x2, 8, 12, list(scale = "difference", value = 0))
  expect_identical(abs(z), rep(sqrt(10 / 3), 6))
  mh <- .ig_unc_stat("mantel_haenszel", x1, x2, 8, 12, list(scale = "difference", value = 0))
  expect_identical(abs(mh), rep(sqrt(19 / 6), 6))
  r <- calculate_power(analysis_spec("exact_proportions_unconditional", "criterion",
    list(p1 = .8, p2 = .3, n1 = 8, n2 = 12, power = .7),
    list(effect_size_index = "proportion_p1", statistic = "z_pooled", tails = "two", computation = "exact")))
  expect_identical(r$status$code, "ok")
  g <- expand.grid(x1 = 0:8, x2 = 0:12)
  cross <- g$x1 * 12 - g$x2 * 8
  total <- g$x1 + g$x2
  # Integer cross multiplication against the exact rational threshold 10/3.
  # No engine statistic, quantile, or search helper is used for membership.
  member <- total > 0 & total < 20 & 3 * 20 * cross^2 >= 10 * 8 * 12 * total * (20 - total)
  probabilities <- dbinom(g$x1, 8, .8) * dbinom(g$x2, 12, .3)
  expect_equal(r$achieved$power, sum(probabilities[member]), tolerance = 1e-13)
  expect_gte(r$achieved$power, .7)
})

test_that("pooled, corrected pooled and Mantel statistics preserve small-grid rational ties", {
  gcd <- function(a, b) { while (b != 0) { r <- a %% b; a <- b; b <- r }; a }
  for (n1 in 1:8) for (n2 in 1:8) {
    g <- expand.grid(x1 = 0:n1, x2 = 0:n2)
    cross <- g$x1 * n2 - g$x2 * n1; total <- g$x1 + g$x2; N <- n1 + n2
    for (stat in c("z_pooled", "z_pooled_cc", "mantel_haenszel")) {
      num <- if (stat == "z_pooled_cc") N * pmax(2 * abs(cross) - N, 0)^2 else
        (if (stat == "mantel_haenszel") N - 1 else N) * cross^2
      den <- n1 * n2 * total * (N - total) * if (stat == "z_pooled_cc") 4 else 1
      keys <- vapply(seq_along(num), function(i) {
        if (num[i] == 0) return("0/1")
        d <- gcd(num[i], den[i]); paste(num[i] / d, den[i] / d, sep = "/")
      }, "")
      z <- .ig_unc_stat(stat, g$x1, g$x2, n1, n2, list(scale = "difference", value = 0))
      expect_true(all(vapply(split(abs(z), keys), function(v) length(unique(v)) == 1L, FALSE)),
        info = paste(stat, n1, n2))
      expect_identical(z[c(1, nrow(g))], c(0, 0))
    }
  }
})

test_that("exact compromise refuses a target unattainable by deterministic regions", {
  g <- expand.grid(x1 = 0:4, x2 = 0:6)
  cross <- g$x1 * 6 - g$x2 * 4
  total <- g$x1 + g$x2
  z <- sign(cross) * sqrt(ifelse(total %in% c(0, 10), 0, 10 * cross^2 / (24 * total * (10 - total))))
  p0 <- dbinom(g$x1, 4, .3) * dbinom(g$x2, 6, .3)
  p1 <- dbinom(g$x1, 4, .7) * dbinom(g$x2, 6, .3)
  residuals <- vapply(sort(unique(c(-Inf, z, Inf))), function(cut) {
    member <- z >= cut
    sum(p1[!member]) - sum(p0[member])
  }, 0.0)
  expect_true(all(abs(residuals) > 1e-12))
  for (method in c("exact_proportions_unconditional", "exact_proportions_unconditional_offset")) {
    inputs <- if (method == "exact_proportions_unconditional") list(p1 = .7, p2 = .3) else
      list(p1_h0 = .4, p1_h1 = .7, p2 = .3)
    r <- calculate_power(analysis_spec(method, "compromise", c(inputs, list(n1 = 4, n2 = 6, error_ratio = 1)),
      list(effect_size_index = "proportion_p1", computation = "exact")))
    expect_identical(r$status$code, "reference_pending")
    expect_null(r$solved$alpha)
  }
})

test_that("rational ties remain tied up to the default exact design limit", {
  x1 <- c(2, 0, 6, 2, 8, 6); x2 <- c(0, 4, 4, 8, 8, 12)
  for (k in 1:500) for (stat in c("z_pooled", "mantel_haenszel")) {
    z <- .ig_unc_stat(stat, k * x1, k * x2, 8 * k, 12 * k, list(scale = "difference", value = 0))
    expect_identical(length(unique(abs(z))), 1L, info = paste(stat, k))
  }
  # Complementary tables are exact ties even at very unequal allocation.
  for (sizes in list(c(1, 9999), c(97, 9903), c(3999, 6001))) {
    a <- floor(sizes[1] / 3); b <- floor(sizes[2] / 7)
    for (stat in c("z_pooled", "z_pooled_cc", "mantel_haenszel")) {
      z <- .ig_unc_stat(stat, c(a, sizes[1] - a), c(b, sizes[2] - b), sizes[1], sizes[2],
        list(scale = "difference", value = 0))
      expect_identical(abs(z[1]), abs(z[2]), info = paste(stat, sizes[1], sizes[2]))
    }
  }
  expect_error(.ig_unc_stat("z_pooled", 12000, 3000, 30000, 30000,
    list(scale = "difference", value = 0)), class = "inferguidR_numerical_failure")
})
