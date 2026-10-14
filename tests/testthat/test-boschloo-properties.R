# Structural properties of Boschloo's test (layer-3 checks): symmetries,
# boundary tables, ties, dominance over Fisher, nested tails, size control,
# decision/region consistency and independence of the numerical settings.

p_of <- function(x1, n1, x2, n2, v, ...) boschloo_test(x1, n1, x2, n2, variant = v, ...)$p_value

test_that("all-zero and all-success tables give p = 1 in every variant", {
  for (nn in list(c(1, 1), c(4, 9), c(10, 10), c(30, 45))) {
    for (v in c("less", "greater", "minlike", "central")) {
      expect_identical(p_of(0, nn[1], 0, nn[2], v), 1)
      expect_identical(p_of(nn[1], nn[1], nn[2], nn[2], v), 1)
    }
  }
})

test_that("swapping the groups swaps the one-sided directions and keeps two-sided values", {
  tabs <- list(c(1, 5, 6, 7), c(3, 12, 8, 9), c(0, 4, 5, 11), c(7, 20, 7, 13))
  for (z in tabs) {
    expect_equal(p_of(z[3], z[4], z[1], z[2], "greater"), p_of(z[1], z[2], z[3], z[4], "less"),
                 tolerance = 1e-10)
    expect_equal(p_of(z[3], z[4], z[1], z[2], "less"), p_of(z[1], z[2], z[3], z[4], "greater"),
                 tolerance = 1e-10)
    for (v in c("minlike", "central")) {
      expect_equal(p_of(z[3], z[4], z[1], z[2], v), p_of(z[1], z[2], z[3], z[4], v), tolerance = 1e-10)
    }
  }
})

test_that("recoding the event swaps the one-sided directions and keeps two-sided values", {
  tabs <- list(c(1, 5, 6, 7), c(3, 12, 8, 9), c(2, 15, 9, 10))
  for (z in tabs) {
    flip <- c(z[2] - z[1], z[2], z[4] - z[3], z[4])
    expect_equal(p_of(flip[1], flip[2], flip[3], flip[4], "less"),
                 p_of(z[1], z[2], z[3], z[4], "greater"), tolerance = 1e-10)
    for (v in c("minlike", "central")) {
      expect_equal(p_of(flip[1], flip[2], flip[3], flip[4], v), p_of(z[1], z[2], z[3], z[4], v),
                   tolerance = 1e-10)
    }
  }
})

test_that("the ordering statistic is Fisher's p-value and Boschloo never exceeds it", {
  tabs <- list(c(1, 5, 6, 7), c(3, 15, 9, 15), c(10, 35, 24, 61), c(4, 4, 0, 6), c(2, 7, 6, 13))
  for (z in tabs) {
    m <- matrix(c(z[1], z[2] - z[1], z[3], z[4] - z[3]), 2, byrow = TRUE)
    f_two <- stats::fisher.test(m)$p.value
    f_less <- stats::fisher.test(m, alternative = "less")$p.value
    f_greater <- stats::fisher.test(m, alternative = "greater")$p.value
    r_min <- boschloo_test(z[1], z[2], z[3], z[4], variant = "minlike")
    r_less <- boschloo_test(z[1], z[2], z[3], z[4], variant = "less")
    r_gre <- boschloo_test(z[1], z[2], z[3], z[4], variant = "greater")
    r_cen <- boschloo_test(z[1], z[2], z[3], z[4], variant = "central")
    expect_equal(r_min$fisher_p_value, f_two, tolerance = 1e-12)
    expect_equal(r_less$fisher_p_value, f_less, tolerance = 1e-12)
    expect_equal(r_gre$fisher_p_value, f_greater, tolerance = 1e-12)
    expect_equal(r_cen$fisher_p_value, min(1, 2 * min(f_less, f_greater)), tolerance = 1e-12)
    for (r in list(r_min, r_less, r_gre, r_cen)) expect_lte(r$p_value, r$fisher_p_value * (1 + 1e-9))
  }
})

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

test_that("the two two-sided variants are different tests", {
  expect_lt(abs(p_of(1, 5, 6, 7, "minlike") - 0.03857422), 5e-9)
  expect_lt(abs(p_of(1, 5, 6, 7, "central") - 0.04635844), 5e-9)
  expect_lt(abs(p_of(10, 35, 24, 61, "minlike") - 0.3380542), 5e-8)
  expect_lt(abs(p_of(10, 35, 24, 61, "central") - 0.3025764), 5e-8)
  # StatXact 11 values transcribed in exact2x2's slow tests (documented
  # example, compared at 4 decimals): 10/35 vs 24/61, one-sided p2 > p1 and
  # two-sided minlike.
  expect_identical(round(p_of(10, 35, 24, 61, "less"), 4), 0.1513)
  expect_identical(round(p_of(10, 35, 24, 61, "minlike"), 4), 0.3381)
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

test_that("the attained size never exceeds alpha and null power stays below it", {
  for (nn in list(c(4, 4), c(7, 19), c(16, 16), c(25, 10))) {
    for (v in c("less", "greater", "minlike", "central")) {
      a <- if (v %in% c("less", "greater")) 0.025 else 0.05
      reg <- boschloo_region(nn[1], nn[2], alpha = a, variant = v)
      expect_lte(reg$size$upper, a)
      expect_lte(reg$size$value, reg$size$upper)
      for (p in c(0, 0.13, 0.5, reg$size$argmax, 1)) {
        if (is.na(p)) next
        pw <- boschloo_power(nn[1], nn[2], p, p, alpha = a, variant = v)
        expect_lte(pw$power, reg$size$upper + 1e-12)
      }
    }
  }
})

test_that("observed decisions and the design rejection region coincide", {
  for (v in c("less", "greater", "minlike", "central")) {
    reg <- boschloo_region(6, 9, alpha = 0.05, variant = v)
    inreg <- paste(reg$rejected$x1, reg$rejected$x2)
    for (x1 in 0:6) for (x2 in 0:9) {
      dec <- boschloo_test(x1, 6, x2, 9, variant = v, alpha = 0.05)$decision
      expect_identical(dec == "reject", paste(x1, x2) %in% inreg)
      expect_false(dec == "unresolved")
    }
  }
})

test_that("results do not depend on the grid size (rigorous bound and refinement)", {
  ctl <- list(boschloo_control(grid_points = 11L), boschloo_control(),
              boschloo_control(grid_points = 20001L, cache = FALSE))
  for (z in list(c(3, 15, 9, 15), c(10, 35, 24, 61), c(0, 30, 6, 40))) {
    vals <- vapply(ctl, function(cc) p_of(z[1], z[2], z[3], z[4], "central", control = cc), double(1))
    expect_lt(diff(range(vals)), 1e-9)
    res <- boschloo_test(z[1], z[2], z[3], z[4], variant = "minlike", control = ctl[[1]])
    expect_true(all(vapply(res$sides, function(s) s$supremum$converged, logical(1))))
    expect_lte(res$p_upper - res$p_value, 1e-9 + 1e-12)
  }
  regs <- lapply(ctl, function(cc) boschloo_region(20, 17, 0.05, "minlike", control = cc)$rejected)
  expect_identical(regs[[1]], regs[[2]])
  expect_identical(regs[[2]], regs[[3]])
})

test_that("tiny designs and extreme probabilities behave as defined", {
  reg <- boschloo_region(1, 1, alpha = 0.05, variant = "minlike")
  expect_identical(reg$n_rejected, 0L)
  expect_identical(reg$size$value, 0)
  expect_identical(boschloo_power(1, 1, 0, 1, alpha = 0.05, variant = "minlike")$power, 0)
  # n = (3, 3), one-sided: (0, 3) is the only table in the 0.05 region, so
  # power at p1 = 0, p2 = 1 is 1 and at p1 = p2 = 0 it is 0.
  reg3 <- boschloo_region(3, 3, alpha = 0.05, variant = "less")
  expect_identical(paste(reg3$rejected$x1, reg3$rejected$x2), "0 3")
  expect_identical(boschloo_power(3, 3, 0, 1, alpha = 0.05, variant = "less")$power, 1)
  expect_identical(boschloo_power(3, 3, 0, 0, alpha = 0.05, variant = "less")$power, 0)
  expect_lt(abs(reg3$size$value - 0.015625), 1e-12)   # sup_pi (1 - pi)^3 pi^3 = 1/64
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

test_that("the exact fallback agrees with the minlike tolerance rule for every design up to N = 56", {
  # The minlike ordering compares conditional probabilities with a relative
  # tolerance of 1e-7 (as fisher.test). The exact fallback compares the exact
  # weights C(n1, x1) C(n2, x2) without tolerance; both agree if no two
  # distinct weights of a stratum are within 1e-7 of each other.
  pascal <- lapply(0:56, function(n) {
    row <- 1
    for (i in seq_len(n)) row <- c(row, 0) + c(0, row)
    row
  })
  smallest <- Inf
  for (n1 in 1:55) for (n2 in 1:(56 - n1)) {
    c1 <- pascal[[n1 + 1]]
    c2 <- pascal[[n2 + 1]]
    for (t in 0:(n1 + n2)) {
      x1 <- max(0, t - n2):min(t, n1)
      u <- sort(unique(c1[x1 + 1] * c2[t - x1 + 1]))
      if (length(u) > 1L) smallest <- min(smallest, min(diff(u) / u[-length(u)]))
    }
  }
  expect_gt(smallest, 1e-7)
})
