test_that("observed Boschloo reports retain independent groups and numerical bounds", {
  for (variant in c("greater", "less", "central", "minlike")) {
    x1 <- if (variant == "less") 0L else 1L
    x2 <- 1L - x1
    r <- boschloo_test(x1, 1L, x2, 1L, variant, alpha = .05,
      groups = c("Treatment", "Control"), event = "Response")
    # For one Bernoulli observation per group the directional rejection
    # probability is p*(1-p), whose maximum is 1/4. Central doubles it.
    # Minlike: every possible table has two-sided Fisher p=1, so the
    # unconditional tail includes all tables and has probability one.
    expected <- if (variant %in% c("greater", "less")) .25 else if (variant == "central") .5 else 1
    expect_equal(r$p_value, expected, tolerance = 1e-8)
    x <- create_report_snapshot(analysis_envelope(r$spec,r))
    expect_equal(x$results[[1]]$p_value, expected, tolerance = 1e-8)
    expect_equal(x$metrics[[1]]$p_value, expected, tolerance = 1e-8)
    expect_identical(x$results[[1]]$decision,"do_not_reject")
    expect_identical(x$uncertainty[[1]]$status,"not_computed")
    expect_match(x$uncertainty[[1]]$reason,"not sampling uncertainty",fixed=TRUE)
    expect_lte(x$uncertainty[[1]]$numerical_bounds$lower, expected + 1e-8)
    expect_gte(x$uncertainty[[1]]$numerical_bounds$upper, expected - 1e-8)
    v <- x$visualisations[[1]]
    expect_identical(v$plot_type,"two_group_counts")
    expect_equal(unlist(v$data$event),c(x1,x2))
    expect_equal(unlist(v$data$non_event),c(1-x1,1-x2))
    expect_identical(unlist(v$data$group_labels),c("Treatment","Control"))
    expect_false(x$reproducibility[[1]]$private_data_required)
    replay <- eval(parse(text=x$reproducibility[[1]]$r_code))
    expect_identical(replay$spec_hash,r$spec_hash)
    expect_identical(ig_to_json(x),ig_to_json(ig_from_json(ig_to_json(x))))
    expect_match(paste(.ig_report_text(x),collapse="\n"),"Response",fixed=TRUE)
  }
})

test_that("observed reports preserve incomplete state and mixed session scope", {
  r <- boschloo_test(1L,1L,0L,1L,"greater")
  e <- analysis_envelope(r$spec,r)
  missing <- analysis_envelope(r$spec)
  expect_error(create_report_snapshot(missing),class="inferguidR_invalid_input")
  x <- create_report_snapshot(missing,report_kind="incomplete_settings")
  expect_identical(x$completeness,"incomplete")
  expect_identical(x$results[[1]]$status$code,"not_computed")
  expect_null(x$metrics[[1]]$p_value)
  expect_identical(ig_to_json(x),ig_to_json(ig_from_json(ig_to_json(x))))
  s <- analysis_spec("t_means_one_sample","post_hoc",list(d=.5,alpha=.05,n_total=30),list(tails="two"))
  x <- create_report_snapshot(list(e,analysis_envelope(s,calculate_power(s))),report_kind="session",locale="de")
  expect_length(x$results,2L)
  expect_identical(x$specifications[[1]]$analysis_kind,"observed_test")
  expect_identical(x$specifications[[2]]$analysis_kind,"power")
  expect_identical(ig_to_json(x),ig_to_json(ig_from_json(ig_to_json(x))))
  changed <- unclass(e);changed$result <- unclass(changed$result);changed$result$p_value <- .1
  class(changed$result) <- class(r);class(changed) <- class(e)
  expect_error(create_report_snapshot(changed),class="inferguidR_invalid_input")
})

test_that("minlike one-observation ordering is independently all tied", {
  for (a in 0:1) for (b in 0:1) {
    tab <- matrix(c(a,1-a,b,1-b),nrow=2,byrow=TRUE)
    expect_equal(stats::fisher.test(tab,alternative="two.sided")$p.value,1)
  }
})
