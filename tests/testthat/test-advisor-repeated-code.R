.repeated_env <- function(method) {
 e<-new.env(parent=baseenv());e$alternative<-"two.sided";e$confidence_level<-.95
 invisible(capture.output(sys.source(system.file("advisor","templates",paste0("tpl.",method,".v1.R"),package="inferguidR"),e)));e
}
test_that("repeated tests agree with independent block covariance references", {
 cases<-jsonlite::fromJSON(test_path("fixtures","repeated-cases.json"),simplifyVector=FALSE)
 refs<-read.csv(test_path("fixtures","repeated-independent.csv"))
 for(method in c("friedman","cochran_q")) {
  e<-.repeated_env(method)
  for(i in which(refs$method==method)) {
   x<-do.call(rbind,lapply(cases[[i]]$data,unlist));r<-e$analyse(data.frame(id=seq_len(nrow(x)),x))
   expect_identical(r$status,"ok");expect_equal(r$df,refs$df[i])
   expect_lte(abs(r$statistic-refs$statistic[i]),1e-12);expect_lte(abs(r$p_value-refs$p_value[i]),1e-12)
  }
 }
})
test_that("matching, missingness and degeneracy are explicit", {
 for(method in c("friedman","cochran_q")) {
  e<-.repeated_env(method);d<-e$example_data
  bad<-d;bad$id[2]<-bad$id[1];expect_error(e$analyse(bad))
  bad<-d;bad$id[2]<-NA;expect_error(e$analyse(bad))
  expect_error(e$analyse(d[,-1]));expect_error(e$analyse(d[,1:3]))
  bad<-d;bad$a[1]<-Inf;expect_error(e$analyse(bad,missing="complete_cases"))
  bad<-d;bad$a[1]<-NA;expect_error(e$analyse(bad))
  r<-e$analyse(bad,missing="complete_cases");direct<-e$analyse(d[-1,])
  expect_identical(r$n_excluded,1L);expect_equal(r$n_complete,5)
  expect_equal(r$p_value,direct$p_value);expect_equal(r$statistic,direct$statistic)
  expect_null(r$confidence_interval);expect_identical(r$missing_strategy,"complete_cases")
  expect_error(e$analyse(d,missing="impute"));expect_error(e$analyse(d,alternative="greater"))
  for(level in c(0,1,NA,Inf)) expect_error(e$analyse(d,conf.level=level))
  all_tied<-d;all_tied[,-1]<-0;r<-e$analyse(all_tied)
  expect_identical(r$status,"unavailable");expect_null(r$p_value)
  expect_equal(e$analyse(d[nrow(d):1,])$p_value,e$analyse(d)$p_value)
  expect_equal(e$analyse(d[,c(1,4,2,3)])$p_value,e$analyse(d)$p_value)
 }
 e<-.repeated_env("cochran_q");d<-e$example_data;d$a[1]<-2
 expect_error(e$analyse(d));d$b[1]<-NA;expect_error(e$analyse(d,missing="complete_cases"))
 d<-e$example_data;d$a<-factor(d$a);expect_error(e$analyse(d))
 e<-.repeated_env("friedman");d<-e$example_data;d[,-1]<-lapply(d[,-1],function(x) 10*x+7)
 expect_equal(e$analyse(d)$p_value,e$example_result$p_value)
})
test_that("repeated examples export through bilingual advisor routes", {
 for(type in c("ordinal","binary","continuous")) for(lang in c("en","de")) {
  a<-c(q.aim="compare",q.stage="analysis",q.outcome_type=type,q.unit="individuals",
    q.structure="repeated",q.pairing_id="yes",q.estimand_continuous="rank",q.missingness="some")
  s<-advisor_replay(a,language=lang);r<-advisor_recommendation(s)
  expect_identical(r$preferred$analysis_variant_id,if(type=="binary") "cochran_q" else "friedman")
  expect_identical(r$status,"conditional")
  code<-advisor_r_code(s);e<-new.env();invisible(capture.output(eval(parse(text=code),e)))
  expect_identical(e$example_result$status,"ok");expect_null(e$example_result$confidence_interval)
  expect_error(advisor_r_code(s,alternative="greater"),class="inferguidR_invalid_input")
 }
})
