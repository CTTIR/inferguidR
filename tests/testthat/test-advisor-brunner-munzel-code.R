.bm_env<-function() {
 e<-new.env(parent=baseenv());e$alternative<-"two.sided";e$confidence_level<-.95
 invisible(capture.output(sys.source(system.file("advisor","templates","tpl.brunner_munzel.v1.R",package="inferguidR"),e)));e
}
test_that("Brunner-Munzel matches independent pairwise placements and t tails", {
 e<-.bm_env();cases<-jsonlite::fromJSON(test_path("fixtures","bm-cases.json"),simplifyVector=FALSE)
 refs<-jsonlite::fromJSON(test_path("fixtures","bm-independent.json"),simplifyVector=FALSE)
 for(i in seq_along(cases)) for(alt in c("two.sided","greater","less")) {
  r<-e$analyse(unlist(cases[[i]]$x),unlist(cases[[i]]$y),alternative=alt);ref<-refs[[i]]
  expect_identical(r$status,"ok")
  for(key in c("estimate","std_error","df","statistic")) expect_equal(r[[key]],as.numeric(ref[[key]]),tolerance=2e-11)
  expect_equal(r$p_value,as.numeric(ref$tails[[alt]]$p_value),tolerance=2e-12)
  expect_equal(r$confidence_interval,as.numeric(unlist(ref$tails[[alt]][c("lower","upper")])),tolerance=2e-11)
 }
})
test_that("Brunner-Munzel states orientation and refuses degenerate or unsupported inference", {
 e<-.bm_env();x<-e$example_x;y<-e$example_y;r<-e$analyse(x,y)
 expect_equal(e$analyse(y,x)$estimate,1-r$estimate)
 expect_equal(e$analyse(y,x)$statistic,-r$statistic)
 expect_equal(e$analyse(x,y,alternative="greater")$p_value,e$analyse(y,x,alternative="less")$p_value)
 expect_equal(e$analyse(10+3*x,10+3*y)$p_value,r$p_value)
 expect_equal(e$analyse(rev(x),rev(y))$p_value,r$p_value)
 expect_identical(e$analyse(rep(1,10),rep(1,10))$status,"unavailable")
 expect_null(e$analyse(rep(1,10),rep(2,10))$p_value)
 expect_error(e$analyse(x[1:9],y));expect_error(e$analyse(c(x,NA),y));expect_error(e$analyse(x,c(y,Inf)))
 expect_error(e$analyse(x,y,alternative="auto"));expect_error(e$analyse(x,y,conf.level=1))
 r<-e$analyse(x,y,conf.level=1-.Machine$double.eps/2)
 expect_true(all(is.finite(r$unconstrained_interval)));expect_true(all(r$confidence_interval>=0 & r$confidence_interval<=1))
 expect_equal(r$confidence_interval,pmax(0,pmin(1,r$unconstrained_interval)))
})
test_that("Brunner-Munzel override exports recorded direction in both languages", {
 for(lang in c("en","de")) for(alt in c("two.sided","greater","less")) {
  a<-c(q.aim="compare",q.stage="analysis",q.outcome_type="continuous",q.unit="individuals",q.structure="independent",q.independence_check="no",q.n_factors="one",q.n_groups="two",q.covariates="no",q.estimand_continuous="rank",q.direction=if(alt=="two.sided") "two_sided" else "one_sided_prespecified",q.missingness="none")
  s<-advisor_override(advisor_replay(a,language=lang),"brunner_munzel")
  expect_identical(advisor_recommendation(s)$preferred$analysis_variant_id,"brunner_munzel")
  e<-new.env();invisible(capture.output(eval(parse(text=advisor_r_code(s,alternative=alt)),e)))
  expect_identical(e$example_result$alternative,alt);expect_identical(e$example_result$status,"ok")
 }
})
