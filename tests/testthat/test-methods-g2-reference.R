# G2 engines against the definition-based reference of helper-reference.R (noncentral chi-square by
# convolution with central pchisq, noncentral F by integrating it over the denominator chi-square
# density; no qbeta/pbeta/dpois-weighted series and no noncentral routine). A few points only; the
# full independent check is validation/independent/g2/check_g2.R (157 checks).
#
# Tolerance: abs 1e-8 + rel 1e-6 (VALIDATION_CONTRACT starting rule; the full check observed
# differences below 1e-11).

ref_reset_registry()

g2_ref_crit <- function(alpha, df1, df2) {
  # bisection on log c with the reference central tail (lambda = 0)
  g <- function(lc) log(ref_pf(exp(lc), df1, df2, 0, lower_tail = FALSE)) - log(alpha)
  lo <- -3; hi <- 6
  for (i in 1:80) {
    m <- (lo + hi) / 2
    if (g(m) > 0) lo <- m else hi <- m
  }
  exp((lo + hi) / 2)
}

test_that("power, beta and critical F agree with the integral reference", {
  cases <- list(
    list(id = "f_anova_oneway", inputs = list(f = 0.25, n_total = 120L, n_groups = 4L, alpha = 0.05), l = 7.5, d = c(3, 116)),
    list(id = "f_ancova_main_interactions", inputs = list(f = 0.25, n_total = 120L, numerator_df = 2L, n_groups = 6L,
                                                          n_covariates = 1L, alpha = 0.05), l = 7.5, d = c(2, 113)),
    list(id = "f_hotelling_t2_two_groups", inputs = list(delta = 0.6, n1 = 25L, n2 = 35L, n_response = 3L, alpha = 0.01),
         l = 25 * 35 / 60 * 0.36, d = c(3, 56))
  )
  for (cs in cases) {
    r <- calculate_power(analysis_spec(cs$id, "post_hoc", inputs = cs$inputs))
    crit <- g2_ref_crit(cs$inputs$alpha, cs$d[1], cs$d[2])
    expect_true(ref_close(r$quantities$critical, crit, 1e-8, 1e-6), info = cs$id)
    expect_true(ref_close(r$achieved$power, ref_pf(crit, cs$d[1], cs$d[2], cs$l, lower_tail = FALSE), 1e-8, 1e-6), info = cs$id)
    expect_true(ref_close(r$achieved$beta, ref_pf(crit, cs$d[1], cs$d[2], cs$l, lower_tail = TRUE), 1e-8, 1e-6), info = cs$id)
  }
})

test_that("a small beta keeps its relative accuracy", {
  # lambda 120, df (2, 60): beta is about 1e-9; formed as 1 - power it would lose most digits.
  t1 <- .ig_g2_f_test(120, 2, 60, 0.05)
  ref <- ref_pf(t1$critical, 2, 60, 120, lower_tail = TRUE)
  expect_lt(ref, 1e-6)
  expect_true(ref_close(t1$beta, ref, 0, 1e-6))
})
