.lme_env<-function(kind){e<-new.env(parent=baseenv());e$alternative<-'two.sided';e$confidence_level<-.95;invisible(capture.output(sys.source(system.file("advisor","templates",paste0("tpl.",kind,".v1.R"),package="inferguidR"),e)));e}
test_that('random-intercept REML models match independent block covariance likelihoods',{
 refs<-jsonlite::read_json(test_path("fixtures","mixed-lme-oracle.json"))$cases
 for(ref in refs){
  e<-.lme_env(ref$kind);d<-do.call(rbind,lapply(ref$data,as.data.frame));d$factor<-factor(d$level);d$y[d$id%in%unlist(ref$missing_ids)]<-NA
  r<-e$analyse(d,TRUE,TRUE,missing_outcomes='available_case_MAR')
  expect_identical(r$status,'ok');expect_equal(r$factor_means$estimate,unlist(ref$means),tolerance=2e-6)
  expect_equal(r$factor_means$std_error,unlist(ref$se),tolerance=2e-6)
  expect_equal(r$residual_variance,ref$sigma2,tolerance=2e-6);expect_equal(r$random_intercept_variance,ref$tau2,tolerance=2e-6)
  expect_equal(r$statistic,ref$F,tolerance=2e-6);expect_equal(r$df_numerator,ref$df1);expect_equal(r$df_denominator,ref$df2)
  expect_equal(r$p_value,ref$p_value,tolerance=2e-6);expect_equal(r$n_excluded,length(ref$missing_ids))
 }
})
test_that('mixed-model designs and missing outcomes require explicit assumptions',{
 for(kind in c('lmm_repeated_measures','lmm_cluster_random_intercept')){
  e<-.lme_env(kind);d<-e$example_data;run<-function(x=d,...)e$analyse(x,TRUE,TRUE,...)
  expect_error(e$analyse(d));expect_error(e$analyse(d,TRUE,FALSE));expect_error(run(alternative='less'))
  bad<-d;bad$y[1]<-NA;expect_error(run(bad));r<-run(bad,missing_outcomes='available_case_MAR')
  expect_identical(r$status,'ok');expect_equal(r$n_excluded,1)
  bad<-d;bad$id[2]<-bad$id[1];expect_error(run(bad))
  bad<-d;bad$cluster[1]<-NA;expect_error(run(bad))
  bad<-d;bad$y[1]<-Inf;expect_error(run(bad))
  bad<-d;bad$factor<-factor(bad$factor,levels=c(levels(bad$factor),'unused'));expect_error(run(bad))
  bad<-d;bad$y<-1;expect_identical(run(bad)$status,'not_applicable')
  bad<-d;bad$y<-d$y*100+1e10;expect_equal(run(bad)$p_value,run()$p_value,tolerance=2e-6)
  rr<-run(d[rev(seq_len(nrow(d))),]);expect_equal(rr$p_value,run()$p_value,tolerance=2e-6)
  if(kind=='lmm_repeated_measures'){
   bad<-d;bad$factor[2]<-bad$factor[1];expect_error(run(bad))
  }else{
   bad<-d;bad$factor[2]<-levels(bad$factor)[2];expect_error(run(bad))
  }
 }
})

.binary_mixed_env<-function(kind){e<-new.env(parent=baseenv());e$alternative<-'two.sided';e$confidence_level<-.95;invisible(capture.output(sys.source(system.file("advisor","templates",paste0("tpl.",kind,".v1.R"),package="inferguidR"),e)));e}
test_that('conditional and marginal binary models match independent likelihood/sandwich references',{
 refs<-jsonlite::read_json(test_path("fixtures","mixed-binary-oracle.json"))
 for(kind in c('glmm_logistic_random_intercept','gee_logistic_exchangeable')){
  e<-.binary_mixed_env(kind);r<-e$example_result;gee<-kind=='gee_logistic_exchangeable';ref<-refs[[if(gee)'gee' else 'glmm']]
  expect_identical(r$status,'ok');expect_equal(unname(r$coefficients),unlist(ref$beta),tolerance=if(gee)1e-9 else 2e-5)
  expect_equal(unname(r$covariance),do.call(rbind,lapply(ref$covariance,unlist)),tolerance=if(gee)1e-9 else 2e-4)
  expect_equal(r$p_value,ref$p_value,tolerance=if(gee)1e-9 else 2e-4)
  if(gee)expect_equal(r$working_correlation,ref$working_correlation,tolerance=1e-9) else expect_equal(r$random_intercept_variance,ref$tau2,tolerance=2e-5)
  expect_equal(r$clusters,30);expect_equal(r$n_used,180);expect_true(all(is.finite(r$confidence_interval)))
  d<-e$example_data;d$group<-factor(d$group,levels=rev(levels(d$group)))
  reverse<-e$analyse(d,TRUE,TRUE);expect_equal(reverse$log_odds_ratio,-r$log_odds_ratio,tolerance=if(gee)1e-9 else 2e-5)
  expect_equal(reverse$p_value,r$p_value,tolerance=if(gee)1e-9 else 2e-4)
 }
})
test_that('binary cluster models retain explicit estimand and missingness boundaries',{
 for(kind in c('glmm_logistic_random_intercept','gee_logistic_exchangeable')){
  e<-.binary_mixed_env(kind);d<-e$example_data;gee<-kind=='gee_logistic_exchangeable';run<-function(x=d,...)e$analyse(x,TRUE,TRUE,...)
  expect_error(e$analyse(d));expect_error(e$analyse(d,TRUE,FALSE));expect_error(run(alternative='greater'))
  bad<-d;bad$id[2]<-bad$id[1];expect_error(run(bad))
  bad<-d;bad$group[2]<-levels(bad$group)[2];expect_error(run(bad))
  bad<-d;bad$y[1]<-.5;expect_error(run(bad))
  bad<-d;bad$y[1]<-NA;expect_error(run(bad));r<-run(bad,missing_outcomes=if(gee)'available_case_MCAR' else 'available_case_MAR')
  expect_identical(r$status,'ok');expect_equal(r$n_excluded,1)
  if(gee)expect_error(run(bad,missing_outcomes='available_case_MAR'))
  bad<-d;bad$y<-0;expect_error(run(bad))
  expect_error(run(d[d$cluster<=10,]))
  shuffled<-d[rev(seq_len(nrow(d))),];expect_equal(run(shuffled)$p_value,run()$p_value,tolerance=if(gee)1e-9 else 2e-4)
 }
})

test_that('clustered and repeated advisor routes preserve covariance and estimand choices',{
 routes<-list(
  lmm_repeated_measures=c(q.aim='compare',q.stage='analysis',q.outcome_type='continuous',q.unit='individuals',q.structure='repeated',q.pairing_id='yes',q.estimand_continuous='mean',q.missingness='some'),
  lmm_cluster_random_intercept=c(q.aim='compare',q.stage='analysis',q.outcome_type='continuous',q.unit='clusters'),
  glmm_logistic_random_intercept=c(q.aim='compare',q.stage='analysis',q.outcome_type='binary',q.unit='clusters'))
 for(lang in c('en','de'))for(kind in names(routes)){
  s<-advisor_replay(routes[[kind]],language=lang);r<-advisor_recommendation(s)
  expect_identical(r$preferred$analysis_variant_id,kind);expect_identical(r$r_implementation$template_status,'available')
  code<-advisor_r_code(s);expect_match(code,r$recommendation_hash,fixed=TRUE)
  e<-new.env(parent=baseenv());invisible(capture.output(eval(parse(text=code),e)))
  expect_identical(e$example_result$status,'ok');expect_true(is.finite(e$example_result$p_value))
  expect_error(advisor_r_code(s,alternative='greater'))
  if(kind=='glmm_logistic_random_intercept'){
   ss<-advisor_override(s,'gee_logistic_exchangeable');rr<-advisor_recommendation(ss)
   expect_identical(rr$preferred$analysis_variant_id,'gee_logistic_exchangeable');expect_identical(rr$r_implementation$template_status,'available')
   e<-new.env(parent=baseenv());invisible(capture.output(eval(parse(text=advisor_r_code(ss)),e)))
   expect_identical(e$example_result$status,'ok');expect_match(e$example_result$estimand,'marginal',fixed=TRUE)
  }
 }
})
