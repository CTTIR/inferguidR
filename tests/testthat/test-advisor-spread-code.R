.spread_code_env <- function(id) {
 e<-new.env();e$alternative<-"two.sided";e$confidence_level<-.95
 invisible(capture.output(sys.source(system.file("advisor","templates",paste0("tpl.",id,".v1.R"),package="inferguidR"),e)));e
}
.spread_route <- function(groups="more_than_two") c(q.aim="compare",q.stage="analysis",
 q.outcome_type="continuous",q.unit="individuals",q.structure="independent",
 q.independence_check="no",q.n_factors="one",q.n_groups=groups,q.covariates="no",
 q.estimand_continuous="spread",q.missingness="none")

test_that("spread overrides honor the recorded group count and rejected routes", {
 for(groups in c("two","more_than_two")) {
  s<-advisor_replay(.spread_route(groups));r<-advisor_recommendation(s)
  expect_identical(r$preferred$analysis_variant_id,"variance_brown_forsythe")
  f<-Filter(function(a) a$analysis_variant_id=="variance_f_test_two",r$alternatives)[[1]]
  expect_identical(f$override_allowed,groups=="two")
  expect_error(advisor_override(s,"variance_chisq_one_sample"),"rejected_route")
  ids<-c("variance_brown_forsythe","variance_fligner_killeen")
  if(groups=="two") ids<-c(ids,"variance_f_test_two") else
   expect_error(advisor_override(s,"variance_f_test_two"),"incompatible_group_count")
  for(id in ids) {
   chosen<-advisor_override(s,id)
   expect_identical(advisor_restore(advisor_save(chosen)),chosen)
   code<-advisor_r_code(chosen);expect_match(code,id,fixed=TRUE)
   e<-new.env();invisible(capture.output(eval(parse(text=code),e)))
   expect_true(is.finite(e$example_result$p_value))
   expect_error(advisor_r_code(chosen,alternative="greater"),"nondirectional")
  }
 }
})

test_that("F spread example agrees with rational F(2,2) reference and rejects invalid data", {
 e<-.spread_code_env("variance_f_test_two")
 for(alt in c("two.sided","greater","less")) {
  r<-e$analyse(c(-2,0,2),c(-1,0,1),alternative=alt)
  expect_equal(r$statistic,4);expect_equal(r$variance_ratio,4)
  expect_equal(r$p_value,switch(alt,two.sided=.4,greater=.2,less=.8))
  q<-function(p) p/(1-p)
  ci<-switch(alt,two.sided=c(4/q(.975),4/q(.025)),greater=c(4/q(.95),Inf),less=c(0,4/q(.05)))
  expect_equal(r$confidence_interval,ci)
 }
 for(x in list(c(1),c(1,1),c(NA,2),c(Inf,2),matrix(1:3,1),c(-1e200,0,1e200)))
  expect_error(e$analyse(x,c(1,2,3)))
 for(null in list(0,-1,NA,Inf,c(1,2))) expect_error(e$analyse(1:3,2:4,ratio0=null))
 expect_error(e$analyse(c(-1e150,0,1e150),c(-1e-150,0,1e-150)))
})

test_that("robust spread examples preserve omnibus and degenerate behavior", {
 for(id in c("variance_brown_forsythe","variance_fligner_killeen")) {
  e<-.spread_code_env(id);x<-c(-1,0,1,-2,0,2,-3,0,3);g<-rep(c("A","B","C"),each=3)
  r<-e$analyse(x,g);expect_identical(r$status,"ok");expect_null(r$confidence_interval)
  expect_true(r$p_value>=0&&r$p_value<=1)
  expect_equal(e$analyse(x,factor(g,levels=c("A","B","C","unused")))$statistic,r$statistic)
  if(id=="variance_brown_forsythe") {
   s<-6*sum((1:3-2)^2)/sum((1:3)^2)
   expect_equal(r$statistic,s);expect_equal(r$p_value,(3/(3+s))^3)
  }
  deg<-e$analyse(rep(1,9),g);expect_identical(deg$status,"not_applicable");expect_null(deg$p_value)
  expect_error(e$analyse(x,g,"less"))
  for(bad in list(rep("A",9),1:9,c(g[-1],NA),c(g[-1],""),g[-1])) expect_error(e$analyse(x,bad))
  for(bad in list(c(NA,x[-1]),c(Inf,x[-1]),matrix(x,3))) expect_error(e$analyse(bad,g))
 }
})

test_that("restored F overrides cannot bypass the group-count restriction", {
 s<-advisor_override(advisor_replay(.spread_route("two")),"variance_f_test_two")
 doc<-jsonlite::parse_json(advisor_save(s),simplifyVector=FALSE)
 for(i in seq_along(doc$answers)) if(doc$answers[[i]]$question_id=="q.n_groups")
  doc$answers[[i]]$option_id<-"more_than_two"
 doc$state_hash<-NULL;doc$state_hash<-ig_hash(doc)
 expect_error(advisor_restore(ig_to_json(doc)),"incompatible_group_count")
})
