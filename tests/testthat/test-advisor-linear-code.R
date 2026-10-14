.linear_env<-function(method) {
 e<-new.env(parent=baseenv());e$alternative<-"two.sided";e$confidence_level<-.95
 invisible(capture.output(sys.source(system.file("advisor","templates",paste0("tpl.",method,".v1.R"),package="inferguidR"),e)));e
}
test_that("linear examples agree with independent high-precision fits and tails", {
 cases<-jsonlite::fromJSON(test_path("fixtures","linear-cases.json"),simplifyVector=FALSE)
 refs<-jsonlite::fromJSON(test_path("fixtures","linear-independent.json"),simplifyVector=FALSE)
 check<-function(actual,expected) expect_lte(abs(actual-as.numeric(expected))/max(1,abs(as.numeric(expected))),2e-11)
 for(method in c("regression_linear_simple","regression_linear_multiple")) {
  e<-.linear_env(method)
  for(i in seq_along(cases)) {
   z<-cases[[i]];if((z$k==1)!=(method=="regression_linear_simple")) next
   x<-do.call(rbind,lapply(z$x,unlist));namesx<-paste0("x",seq_len(z$k));colnames(x)<-namesx;d<-data.frame(id=seq_len(z$n),y=unlist(z$y),x)
   r<-if(z$k==1) e$analyse(d,predictor=namesx) else e$analyse(d,predictors=namesx,tested_predictors=tail(namesx,1));ref<-refs[[i]]
   expect_identical(r$status,"ok")
   for(j in seq_len(z$k+1)) for(key in c("estimate","std_error","statistic","p_value","lower","upper")) check(r$coefficient_table_two_sided[[key]][j],ref$coefficients[[j]][[key]])
   check(r$r_squared,ref$r_squared)
   for(key in c("statistic","p_value")) check(r$omnibus[[key]],ref$omnibus[[key]])
   if(z$k==1) for(alternative in c("greater","less")) {
    a<-e$analyse(d,predictor=namesx,alternative=alternative);check(a$p_value,ref$coefficients[[2]][[paste0(alternative,"_p")]])
    check(a$confidence_interval[if(alternative=="greater") 1 else 2],ref$coefficients[[2]][[if(alternative=="greater") "greater_lower" else "less_upper"]])
   } else {
    check(r$adjusted_r_squared,ref$adjusted_r_squared)
    for(key in names(ref$partial)) check(r$partial_test[[key]],ref$partial[[key]])
   }
  }
 }
})
test_that("linear examples preserve exclusions, predictor meaning and rank boundaries", {
 for(method in c("regression_linear_simple","regression_linear_multiple")) {
  e<-.linear_env(method);d<-e$example_data
  a<-d;a$id[2]<-a$id[1];expect_error(e$analyse(a))
  a<-d;a$y[1]<-NA;expect_error(e$analyse(a));r<-e$analyse(a,missing="complete_cases")
  expect_identical(r$n_excluded,1L);expect_equal(r$coefficient_table_two_sided,e$analyse(d[-1,])$coefficient_table_two_sided)
  a<-d;a$y[1]<-Inf;expect_error(e$analyse(a,missing="complete_cases"))
  a<-d;a$y<-1;expect_error(e$analyse(a))
  a<-d;a$y<-2*d[[if(method=="regression_linear_simple") "x" else "x1"]]+1;expect_error(e$analyse(a))
  expect_error(e$analyse(d[1:2,]));expect_error(e$analyse(d,conf.level=1));expect_error(e$analyse(d,missing="impute"))
  expect_true(all(is.finite(as.matrix(e$analyse(d,conf.level=1-.Machine$double.eps/2)$coefficient_table_two_sided[,c("lower","upper")]))))
  a<-d;a$y<-a$y+10;r<-e$analyse(a)$coefficient_table_two_sided;ref<-e$analyse(d)$coefficient_table_two_sided
  expect_equal(r$estimate[-1],ref$estimate[-1]);expect_equal(r$estimate[1],ref$estimate[1]+10)
  expect_equal(e$analyse(d[nrow(d):1,])$coefficient_table_two_sided,ref)
 }
 e<-.linear_env("regression_linear_multiple");d<-e$example_data
 a<-d;a$x2<-a$x1;expect_error(e$analyse(a))
 expect_error(e$analyse(d,predictors=c("x1","x1")))
 expect_error(e$analyse(d,tested_predictors="missing"));expect_error(e$analyse(d,alternative="greater"))
 r<-e$analyse(d,tested_predictors=c("x1","x2"));expect_equal(r$partial_test$p_value,r$omnibus$p_value)
 expect_null(e$analyse(d)$partial_test)
 # Literal non-syntactic column names cannot become formula expressions.
 a<-d;names(a)[match(c("x1","x2"),names(a))]<-c("odd name","stop('must not execute')")
 expect_equal(e$analyse(a,predictors=c("odd name","stop('must not execute')"))$p_value,e$analyse(d)$p_value)
})
test_that("linear models export their recorded hypotheses in both languages", {
 for(np in c("one","several")) for(lang in c("en","de")) {
  a<-c(q.aim="model",q.stage="analysis",q.outcome_type="continuous",q.unit="individuals",q.records_per_individual="one",q.n_predictors=np,q.missingness="none")
  s<-advisor_replay(a,language=lang);method<-if(np=="one") "regression_linear_simple" else "regression_linear_multiple"
  expect_identical(advisor_recommendation(s)$preferred$analysis_variant_id,method)
  e<-new.env();invisible(capture.output(eval(parse(text=advisor_r_code(s)),e)));expect_identical(e$example_result$status,"ok")
  expect_error(advisor_r_code(s,alternative="greater"),class="inferguidR_invalid_input")
 }
 a<-c(q.aim="associate",q.stage="analysis",q.unit="individuals",q.records_per_individual="one",q.assoc_types="both_continuous",q.assoc_form="linear",q.direction="one_sided_prespecified",q.missingness="none")
 expect_error(advisor_override(advisor_replay(a),"regression_linear_simple"),"override_not_allowed")
})
