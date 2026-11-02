# Numerical regression blocks retained from the historical guidance suite.
p_of <- function(x1, n1, x2, n2, v, ...) boschloo_test(x1, n1, x2, n2, variant = v, ...)$p_value
test_that("tied tables enter the tail together (relative tie tolerance)", {
  # n = (5, 7): (1, 6) and (2, 7) have the same lower Fisher p-value 1/22;
  # dropping (2, 7), as a floating-point equality test does, gives 0.0174.
  r <- boschloo_test(1, 5, 6, 7, variant = "less")
  expect_lt(abs(r$p_value - 0.02317922), 5e-9)
  ord <- .bl_ordering(5L, 7L, "lower", boschloo_control())
  idx <- function(x1, x2) x1 + 1L + x2 * 6L
  expect_identical(ord$level[idx(1L, 6L)], ord$level[idx(2L, 7L)])
  tail_without <- setdiff(.bl_level_tables(ord, ord$level[idx(1L, 6L)]), idx(2L, 7L))
  dropped <- .bl_sup(.bl_tail_c(ord, tail_without), boschloo_control())$value
  expect_gt(r$p_value - dropped, 0.005)
  # n = (15, 15) minlike: (3, 9), (9, 3), (6, 12) and (12, 6) share the
  # two-sided Fisher p-value.
  r2 <- boschloo_test(3, 15, 9, 15, variant = "minlike")
  expect_lt(abs(r2$p_value - 0.03006801), 5e-9)
  ord2 <- .bl_ordering(15L, 15L, "minlike", boschloo_control())
  id2 <- function(x1, x2) x1 + 1L + x2 * 16L
  lv <- ord2$level[c(id2(3L, 9L), id2(9L, 3L), id2(6L, 12L), id2(12L, 6L))]
  expect_identical(length(unique(lv)), 1L)
})
test_that("p-values are monotone along the nested tail sets", {
  ctl <- boschloo_control()
  for (v in c("less", "greater", "minlike")) {
    info <- .bl_variant(v)
    pv <- .bl_pvalues_all(9L, 13L, info, ctl)
    fisher <- .bl_fisher_p(9L, 13L, info, ctl)
    o <- order(fisher)
    # Ties in the ordering share the same tail and hence the same p-value;
    # otherwise the supremum over a larger tail is at least as large.
    expect_true(all(diff(pv$value[o]) >= -1e-9), info = v)
  }
  # The central p-value is monotone within each one-sided ordering, not along
  # the doubled Fisher value; it is twice the smaller one-sided p-value.
  cen <- .bl_pvalues_all(9L, 13L, .bl_variant("central"), ctl)$value
  lo <- .bl_pvalues_all(9L, 13L, .bl_variant("less"), ctl)$value
  up <- .bl_pvalues_all(9L, 13L, .bl_variant("greater"), ctl)$value
  expect_identical(cen, pmin(1, 2 * pmin(lo, up)))
})
test_that("ambiguous ties above n1 + n2 = 56 are refused exactly where they could change a result", {
  # Band (1e-12, 1e-4]: gaps in it are ambiguous. Alternative A separates all
  # of them (tie tolerance 1e-12), alternative B merges all of them (tie
  # tolerance 1e-4). A result that is not flagged must be the same under A,
  # under B and under the default tolerance.
  band <- boschloo_control(ambiguity_band = c(1e-12, 1e-4))
  # The alternatives are hypothetical orderings without ambiguity flags.
  # boschloo_control() refuses a tolerance outside the band, so they are built
  # from a valid object here, only to show what each tie decision would give.
  hypothetical <- function(tie) {
    cc <- boschloo_control()
    cc$tie_tolerance <- tie
    cc$ambiguity_band <- c(Inf, Inf)
    cc
  }
  alt_a <- hypothetical(1e-12)
  alt_b <- hypothetical(1e-4)
  info <- .bl_variant("less")
  ord <- .bl_ordering(30L, 45L, "lower", band)
  expect_identical(ord$arithmetic, "floating_tolerance")
  expect_gt(ord$ambiguous_gaps, 0L)
  flag <- vapply(seq_along(ord$level), function(i) .bl_ambiguity_touches(ord, ord$level[i], "tail"),
                 logical(1))
  expect_gt(sum(flag), 0L)
  pv <- .bl_pvalues_all(30L, 45L, info, band)$value
  pa <- .bl_pvalues_all(30L, 45L, info, alt_a)$value
  pb <- .bl_pvalues_all(30L, 45L, info, alt_b)$value
  expect_lt(max(abs(pv[!flag] - pa[!flag])), 1e-12)
  expect_lt(max(abs(pv[!flag] - pb[!flag])), 1e-12)
  # The public test refuses the flagged tables and answers the others.
  for (i in c(which(flag), which(!flag)[c(1, 200, 700)])) {
    r <- boschloo_test(ord$x1[i], 30, ord$x2[i], 45, variant = "less", control = band)
    if (flag[i]) {
      expect_identical(r$status$code, "numerical_failure")
      expect_true(is.na(r$p_value))
      expect_match(r$status$reason, "ambiguous tie")
    } else {
      expect_identical(r$status$code, "ok")
      expect_lt(abs(r$p_value - pv[i]), 1e-12)
    }
  }
  # Region cut at the level just before the ambiguous low gap: refused (B
  # would merge the next level into it); one level earlier: certified and
  # the same region under A and B.
  amb <- ord$ambiguity
  j <- amb$level_low[!amb$merged & amb$value_low <= 0.5][1]
  p_level <- function(k) pv[ord$perm[ord$level_end[k]]]
  at_cut <- boschloo_region(30, 45, alpha = p_level(j) + 1e-10, variant = "less", control = band)
  expect_identical(at_cut$status$code, "numerical_failure")
  expect_identical(nrow(at_cut$rejected), 0L)
  b_cut <- boschloo_region(30, 45, alpha = p_level(j) + 1e-10, variant = "less", control = alt_b)
  d_cut <- boschloo_region(30, 45, alpha = p_level(j) + 1e-10, variant = "less", control = alt_a)
  expect_false(identical(b_cut$rejected, d_cut$rejected))
  earlier <- (p_level(j - 1L) + p_level(j)) / 2
  ok <- boschloo_region(30, 45, alpha = earlier, variant = "less", control = band)
  expect_identical(ok$status$code, "ok")
  expect_identical(ok$rejected, boschloo_region(30, 45, alpha = earlier, variant = "less", control = alt_a)$rejected)
  expect_identical(ok$rejected, boschloo_region(30, 45, alpha = earlier, variant = "less", control = alt_b)$rejected)
})
