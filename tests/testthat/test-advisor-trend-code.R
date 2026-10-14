.trend_env<-function() {
 e<-new.env(parent=baseenv());e$alternative<-"two.sided";e$confidence_level<-.95
 invisible(capture.output(sys.source(system.file("advisor","templates","tpl.cochran_armitage_trend.v1.R",package="inferguidR"),e)));e
}
test_that("ordered-proportion trend matches independent binomial score references", {
 e<-.trend_env();refs<-jsonlite::fromJSON(test_path("fixtures","trend-independent.json"),simplifyVector=FALSE)
 for(ref in refs) for(alt in c("two.sided","greater","less")) {
  x<-unlist(ref$events);n<-unlist(ref$totals);scores<-unlist(ref$scores);names(x)<-names(n)<-names(scores)<-LETTERS[seq_along(x)]
  r<-e$analyse(x,n,scores,TRUE,alt)
  expect_identical(r$status,"ok");expect_equal(r$signed_statistic,as.numeric(ref$z),tolerance=2e-12)
  expect_equal(r$statistic,as.numeric(ref$statistic),tolerance=2e-12)
  expect_equal(r$p_value,as.numeric(ref$tails[[alt]]),tolerance=2e-12)
 }
})
test_that("trend keeps score choices explicit and rejects non-binomial data", {
 e<-.trend_env();x<-c(A=3,B=8,C=13);n<-c(A=20,B=25,C=30);s<-c(A=0,B=1,C=2)
 r<-e$analyse(x,n,s,TRUE)
 expect_equal(e$analyse(x,n,10+5*s,TRUE)$p_value,r$p_value)
 expect_equal(e$analyse(n-x,n,s,TRUE)$signed_statistic,-r$signed_statistic)
 expect_equal(e$analyse(x,n,s,TRUE,"greater")$p_value,e$analyse(n-x,n,s,TRUE,"less")$p_value)
 expect_identical(e$analyse(x*0,n,s,TRUE)$status,"unavailable");expect_null(e$analyse(n,n,s,TRUE)$p_value)
 expect_error(e$analyse(x,n,s));expect_error(e$analyse(x,n,s,FALSE))
 expect_error(e$analyse(x,n,rev(s),TRUE));expect_error(e$analyse(x,n,c(A=0,B=0,C=1),TRUE))
 expect_error(e$analyse(x,n,unname(s),TRUE));expect_error(e$analyse(unname(x),n,s,TRUE))
 expect_error(e$analyse(x[1:2],n[1:2],s[1:2],TRUE));expect_error(e$analyse(x+.5,n,s,TRUE))
 expect_error(e$analyse(x,c(A=2,B=25,C=30),s,TRUE));expect_error(e$analyse(x,n,s,TRUE,conf.level=1))
 expect_null(r$confidence_interval);expect_identical(r$interval_status,"not_computed")
 expect_equal(r$statistic,r$signed_statistic^2)
 # Unequally spaced prespecified scores generally describe a different trend contrast.
 expect_gt(abs(e$analyse(x,n,c(A=0,B=1,C=5),TRUE)$signed_statistic-r$signed_statistic),.01)
})
test_that("trend override requires recorded binary outcomes in both languages", {
 for(lang in c("en","de")) {
  a<-c(q.aim="compare",q.stage="analysis",q.outcome_type="binary",q.unit="individuals",q.structure="independent",q.independence_check="no",q.n_factors="one",q.n_groups="more_than_two",q.covariates="no",q.small_expected="no",q.missingness="none")
  s<-advisor_override(advisor_replay(a,language=lang),"cochran_armitage_trend")
  expect_identical(advisor_recommendation(s)$preferred$analysis_variant_id,"cochran_armitage_trend")
  e<-new.env();invisible(capture.output(eval(parse(text=advisor_r_code(s)),e)))
  expect_identical(e$example_result$status,"ok");expect_identical(e$example_result$alternative,"two.sided")
  expect_error(advisor_r_code(s,alternative="greater"),class="inferguidR_invalid_input")
  a["q.outcome_type"]<-"nominal"
  expect_error(advisor_override(advisor_replay(a,language=lang),"cochran_armitage_trend"),"incompatible_outcome_type")
  alts<-advisor_recommendation(advisor_replay(a,language=lang))$alternatives
  shown<-Filter(function(v) v$analysis_variant_id=="cochran_armitage_trend",alts)
  expect_length(shown,1L);expect_false(shown[[1]]$override_allowed)
 }
})
