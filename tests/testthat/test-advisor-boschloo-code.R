test_that("Boschloo advisor examples preserve ordering and directional meaning", {
 ids <- c("boschloo_two_sided_central","boschloo_two_sided_minlike","boschloo_one_sided")
 for(id in ids) {
  alts <- if(id=="boschloo_one_sided") c("greater","less") else "two.sided"
  for(alt in alts) {
   e <- new.env();e$alternative <- alt;e$confidence_level <- .95
   path <- system.file("advisor","templates",paste0("tpl.",id,".v1.R"),package="inferguidR")
   expect_true(nzchar(path));invisible(capture.output(sys.source(path,e)))
   tab <- if(alt=="less") matrix(c(0,1,1,0),2,byrow=TRUE) else matrix(c(1,0,0,1),2,byrow=TRUE)
   r <- e$analyse(tab,alt)
   expected <- if(id=="boschloo_two_sided_central") .5 else if(id=="boschloo_two_sided_minlike") 1 else .25
   expect_equal(r$p_value,expected,tolerance=1e-8)
   expect_gte(r$p_upper,expected-1e-10)
   expect_identical(r$status,"ok");expect_identical(r$alternative,alt)
   expect_equal(r$event_probability_difference,if(alt=="less") -1 else 1)
   expect_null(r$confidence_interval);expect_identical(r$interval_status,"not_computed")
   expect_match(r$numerical_bound_note,"not sampling uncertainty",fixed=TRUE)
   expect_error(e$analyse(matrix(0,2,2),alt))
   expect_error(e$analyse(matrix(c(1,.5,1,0),2,byrow=TRUE),alt))
   expect_error(e$analyse(matrix(c(1,NA,1,0),2,byrow=TRUE),alt))
   expect_error(e$analyse(tab,alt,conf.level=1))
   expect_error(e$analyse(tab,if(alt=="two.sided") "greater" else "two.sided"))
   expect_equal(e$analyse(matrix(c(0,1,0,1),2,byrow=TRUE),alt)$p_value,1)
  }
 }
})

test_that("Boschloo advisor routes export the selected engine convention", {
 base <- c(q.aim="compare",q.stage="analysis",q.outcome_type="binary",q.unit="individuals",
  q.structure="independent",q.independence_check="no",q.n_factors="one",q.n_groups="two",
  q.covariates="no",q.binary_margins="group_sizes_fixed",q.direction="two_sided",q.missingness="none")
 for(id in c("boschloo_two_sided_central","boschloo_two_sided_minlike","boschloo_one_sided")) {
  ans <- base;if(id=="boschloo_one_sided") ans["q.direction"] <- "one_sided_prespecified"
  s <- advisor_replay(ans)
  if(id=="boschloo_two_sided_minlike") s <- advisor_override(s,id)
  r <- advisor_recommendation(s)
  expect_identical(r$preferred$analysis_variant_id,id)
  expect_identical(r$r_implementation$template_status,"available")
  alts <- if(id=="boschloo_one_sided") c("greater","less") else "two.sided"
  for(alt in alts) {
   code <- advisor_r_code(s,alternative=alt)
   expect_match(code,"inferguidR::boschloo_test",fixed=TRUE)
   expect_match(code,r$recommendation_hash,fixed=TRUE)
   e <- new.env();invisible(capture.output(eval(parse(text=code),e)))
   expect_identical(e$example_result$alternative,alt)
   expect_identical(e$example_result$observed_result$variant$variant,if(alt!="two.sided") alt else if(id=="boschloo_two_sided_central") "central" else "minlike")
  }
 }
})
