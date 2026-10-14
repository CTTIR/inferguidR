# Boschloo engine versus Exact 3.3 (a separate code base with the same
# definition; layer-2b comparator, Suggests only). Exact's rows are the
# groups (cond.row = TRUE); its "less" is H1: p1 < p2. Exact's central
# variant doubles the one-sided p-value in the observed direction, so tables
# with equal sample proportions are compared only through the one-sided
# values. Tolerance 1e-7 on p-values (the precision of Exact's refinement).

exact_p <- function(x1, n1, x2, n2, variant) {
  tab <- matrix(c(x1, n1 - x1, x2, n2 - x2), 2, byrow = TRUE)
  args <- switch(variant,
                 less = list(alternative = "less"),
                 greater = list(alternative = "greater"),
                 minlike = list(alternative = "two.sided", tsmethod = "square"),
                 central = list(alternative = "two.sided", tsmethod = "central"))
  do.call(Exact::exact.test, c(list(data = tab, method = "boschloo", to.plot = FALSE,
                                    ref.pvalue = TRUE), args))$p.value
}

named_tables <- list(c(1, 5, 6, 7), c(6, 7, 1, 5), c(10, 35, 24, 61), c(7, 15, 12, 15),
                     c(2, 7, 6, 13), c(3, 15, 9, 15), c(0, 10, 5, 10), c(0, 10, 0, 10),
                     c(10, 10, 10, 10), c(0, 5, 3, 5), c(5, 5, 5, 5), c(12, 40, 22, 40),
                     c(4, 4, 0, 6), c(1, 20, 8, 20))

test_that("p-values of the phase 02 named tables agree with Exact 3.3", {
  skip_if_not_installed("Exact")
  for (z in named_tables) {
    for (v in c("less", "greater", "minlike", "central")) {
      ours <- boschloo_test(z[1], z[2], z[3], z[4], variant = v)
      theirs <- exact_p(z[1], z[2], z[3], z[4], v)
      expect_lt(abs(ours$p_value - theirs), 1e-7)
    }
  }
})

test_that("rejection regions agree with Exact::exact.reject.region", {
  skip_if_not_installed("Exact")
  cases <- list(list(10, 12, "less", 0.025), list(12, 10, "greater", 0.025),
                list(9, 9, "minlike", 0.05), list(7, 13, "minlike", 0.05),
                list(11, 8, "central", 0.05), list(6, 15, "less", 0.05))
  for (cs in cases) {
    n1 <- cs[[1]]; n2 <- cs[[2]]; v <- cs[[3]]; a <- cs[[4]]
    args <- switch(v,
                   less = list(alternative = "less"),
                   greater = list(alternative = "greater"),
                   minlike = list(alternative = "two.sided", tsmethod = "square"),
                   central = list(alternative = "two.sided", tsmethod = "central"))
    m <- do.call(Exact::exact.reject.region, c(list(n1 = n1, n2 = n2, alpha = a, method = "boschloo"), args))
    theirs <- which(m == 1, arr.ind = TRUE)
    theirs <- data.frame(x1 = theirs[, 1] - 1, x2 = theirs[, 2] - 1)
    reg <- boschloo_region(n1, n2, alpha = a, variant = v)
    key <- function(df) sort(paste(df$x1, df$x2))
    expect_identical(key(reg$rejected), key(theirs),
                     info = sprintf("(%d, %d) %s alpha %s", n1, n2, v, a))
  }
})

test_that("power agrees with Exact::power.exact.test", {
  skip_if_not_installed("Exact")
  cases <- list(list(5, 7, "less", 0.1, 0.8), list(5, 7, "minlike", 0.5, 0.5),
                list(15, 20, "greater", 0.6, 0.2), list(20, 20, "minlike", 0.2, 0.6),
                list(14, 9, "central", 0.15, 0.7))
  for (cs in cases) {
    args <- switch(cs[[3]],
                   less = list(alternative = "less"),
                   greater = list(alternative = "greater"),
                   minlike = list(alternative = "two.sided", tsmethod = "square"),
                   central = list(alternative = "two.sided", tsmethod = "central"))
    theirs <- do.call(Exact::power.exact.test,
                      c(list(p1 = cs[[4]], p2 = cs[[5]], n1 = cs[[1]], n2 = cs[[2]], alpha = 0.05,
                             method = "boschloo"), args))$power
    ours <- boschloo_power(cs[[1]], cs[[2]], cs[[4]], cs[[5]], alpha = 0.05, variant = cs[[3]])$power
    expect_lt(abs(ours - theirs), 1e-12)
  }
})

test_that("the central convention differs from Exact only for equal sample proportions", {
  skip_if_not_installed("Exact")
  # 14/30 vs 21/45 has equal sample proportions. Both one-sided values agree;
  # inferguidR doubles the smaller one, Exact doubles the 'less' one.
  lo <- boschloo_test(14, 30, 21, 45, variant = "less")$p_value
  up <- boschloo_test(14, 30, 21, 45, variant = "greater")$p_value
  expect_lt(abs(lo - exact_p(14, 30, 21, 45, "less")), 1e-7)
  expect_lt(abs(up - exact_p(14, 30, 21, 45, "greater")), 1e-7)
  ce <- boschloo_test(14, 30, 21, 45, variant = "central")$p_value
  expect_equal(ce, min(1, 2 * min(lo, up)))
  expect_lt(abs(ce - 0.9959702), 5e-8)
  expect_equal(exact_p(14, 30, 21, 45, "central"), 1)
})
