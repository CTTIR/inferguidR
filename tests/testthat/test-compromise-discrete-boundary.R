test_that("discrete compromise returns the attaining side of a boundary", {
  # Two Bernoulli trials under H0 p=.5: the upper endpoint has probability
  # exactly .25. Under H1 p=.99 it has probability .99^2. A value just below
  # alpha=.25 has an empty rejection region and fails beta <= 4*alpha.
  for (options in list(list(tails="one"),
      list(tails="two",alpha_balancing="minor_tail_half_rest_major"),
      list(tails="two",alpha_balancing="half_each_side_then_increase"))) {
    r <- calculate_power(analysis_spec("exact_generic_binomial", "compromise",
      inputs=list(p1=.5,p2=.99,n_total=2,error_ratio=4),options=options))
    expect_identical(r$status$code,"ok")
    expect_true(r$achieved$beta <= 4*r$solved$alpha)
    expect_equal(r$solved$alpha,.25,tolerance=1e-15)
    expect_equal(r$achieved$power,.99^2,tolerance=1e-14)
    expect_equal(r$achieved$alpha_attained,.25,tolerance=1e-14)
    expect_equal(r$achieved$error_ratio,(1-.99^2)/.25,tolerance=1e-14)
    expect_identical(r$diagnostics$root_kind,"jump")
  }
})

test_that("an interior nominal compromise retains its distinct attained size", {
  # The first attainable upper-tail rule needs alpha=.65^2 > .25.
  # The nominal equation has beta=1, q=4 at .25 before any rejection.
  r <- calculate_power(analysis_spec("exact_proportion_one_sample_binomial",
    "compromise",inputs=list(p0=.65,g=.3,n_total=2,error_ratio=4),
    options=list(tails="one")))
  expect_identical(r$status$code,"ok")
  expect_equal(r$solved$alpha,.25,tolerance=1e-15)
  expect_equal(r$achieved$power,0)
  expect_equal(r$achieved$alpha_attained,0)
  expect_true(is.na(r$achieved$error_ratio))
  expect_true(r$achieved$beta <= 4*r$solved$alpha)
})


test_that("boundary correction preserves sign and McNemar definitions", {
  cases <- list(
    list(method="exact_proportion_sign_test", inputs=list(g=.45),
         options=list(tails="one"), power=.95^2),
    list(method="exact_generic_binomial", inputs=list(p1=.5,p2=.01),
         options=list(tails="one"), power=.99^2),
    list(method="exact_proportions_mcnemar", inputs=list(or=20,pd=.95),
         options=list(tails="one",computation="fast_approx"), power=(20/21)^2),
    list(method="exact_proportions_mcnemar", inputs=list(or=1/20,pd=.95),
         options=list(tails="one",computation="fast_approx"), power=(20/21)^2))
  for (case in cases) {
    r <- calculate_power(analysis_spec(case$method,"compromise",
      inputs=c(case$inputs,list(n_total=2,error_ratio=4)),options=case$options))
    expect_identical(r$status$code,"ok")
    expect_true(r$achieved$beta <= 4*r$solved$alpha)
    expect_equal(r$solved$alpha,.25,tolerance=1e-15)
    expect_equal(r$achieved$power,case$power,tolerance=1e-14)
    expect_equal(r$achieved$alpha_attained,.25,tolerance=1e-14)
    expect_equal(r$achieved$error_ratio,(1-case$power)/.25,tolerance=1e-14)
  }
})
