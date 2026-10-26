.glm_env<-function(kind){e<-new.env(parent=baseenv());e$alternative<-'two.sided';e$confidence_level<-.95;invisible(capture.output(sys.source(system.file("advisor","templates",paste0("tpl.regression_",kind,".v1.R"),package="inferguidR"),e)));e}
test_that('finite log-scale regression fits do not report overflowing or underflowing ratios',{
 for(kind in c('logistic_binary','poisson','negative_binomial')) {
  e<-.glm_env(kind);d<-e$example_data
  base<-e$analyse(d,c('x1','x2'),TRUE)
  expect_identical(base$status,'ok')
  expect_true(all(is.na(base$coefficient_table[1,c('ratio','ratio_lower','ratio_upper')])))
  expect_true(is.finite(base$coefficient_table$estimate[1]))
  for(direction in c(-1,1)) {
   tiny<-d;tiny$x1<-d$x1*direction*1e-9
   r<-e$analyse(tiny,c('x1','x2'),TRUE)
   expect_identical(r$status,'not_applicable')
   expect_match(r$reason,'outside the numerical range',fixed=TRUE)
   expect_null(r$p_value);expect_null(r$coefficient_table)
  }
  rescaled<-d;rescaled$x1<-d$x1*10
  r<-e$analyse(rescaled,c('x1','x2'),TRUE)
  expect_identical(r$status,'ok')
  expect_equal(r$coefficient_table$estimate[2],base$coefficient_table$estimate[2]/10,tolerance=1e-9)
  expect_equal(r$p_value,base$p_value,tolerance=1e-9)
 }
})
test_that('GLM coefficients and likelihood tests agree with independent 70-digit references',{
 refs<-jsonlite::read_json(test_path("fixtures","glm-oracle.json"))$cases
 for(ref in refs){
  e<-.glm_env(ref$kind);d<-as.data.frame(do.call(rbind,lapply(ref$x,unlist)));names(d)<-paste0('x',seq_len(ncol(d)));cols<-names(d)
  d$id<-seq_len(nrow(d));d$y<-unlist(ref$y);d$exposure<-unlist(ref$exposure);r<-e$analyse(d,cols,TRUE)
  expect_identical(r$status,'ok');expect_equal(r$coefficient_table$estimate,as.numeric(unlist(ref$beta)),tolerance=1e-9)
  expect_equal(r$coefficient_table$std_error,as.numeric(unlist(ref$se)),tolerance=1e-9)
  expect_equal(r$fitted,as.numeric(unlist(ref$fitted)),tolerance=1e-9)
  expect_equal(r$omnibus$statistic,as.numeric(ref$LR),tolerance=1e-9)
  expect_equal(r$p_value,as.numeric(ref$p_value),tolerance=1e-9)
  expect_equal(r$n_excluded,0);expect_match(r$interval_method,if(ref$kind=='logistic_binary')'profile-likelihood intervals for slope' else 'not profile-likelihood',fixed=TRUE)
  reverse<-e$analyse(d[rev(seq_len(nrow(d))),],cols,TRUE)
  expect_equal(reverse$coefficient_table,r$coefficient_table,tolerance=1e-8)
  expect_equal(reverse$p_value,r$p_value,tolerance=1e-9)
 }
})
test_that('GLM input, missingness, rank, boundary and coding contracts are explicit',{
 for(kind in c('logistic_binary','poisson')){
  e<-.glm_env(kind);d<-e$example_data;cols<-c('x1','x2');run<-function(data=d,...)e$analyse(data,cols,TRUE,...)
  expect_error(e$analyse(d,cols));expect_error(run(alternative='greater'));expect_error(run(conf.level=0))
  bad<-d;bad$id[2]<-bad$id[1];expect_error(run(bad))
  bad<-d;bad$y[1]<-NA;expect_error(run(bad));r<-run(bad,complete_cases=TRUE)
  expect_identical(r$status,'ok');expect_equal(r$n_excluded,1);expect_equal(r$p_value,run(d[-1,])$p_value)
  bad<-d;bad$x1[1]<-Inf;expect_error(run(bad))
  bad<-d;bad$x2<-bad$x1;expect_identical(run(bad)$status,'not_applicable')
  bad<-d;bad$x1<-1;expect_identical(run(bad)$status,'not_applicable')
  bad<-d;bad$y<-0;expect_identical(run(bad)$status,'not_applicable');expect_null(run(bad)$p_value)
  bad<-d;bad$y[1]<-.5;expect_error(run(bad))
  bad<-d;bad$x1<-factor(bad$x1);expect_error(run(bad))
  bad<-d;names(bad)[names(bad)=='x1']<-'x1 + system("invalid")'
  renamed<-e$analyse(bad,c('x1 + system("invalid")','x2'),TRUE)
  expect_identical(renamed$status,'ok');expect_equal(renamed$p_value,run()$p_value)
  shifted<-d;shifted$x1<-d$x1*100+1e6
  expect_equal(run(shifted)$p_value,run()$p_value,tolerance=1e-8)
  if(kind=='logistic_binary'){
   bad<-d;bad$y<-as.numeric(bad$x1>0);expect_identical(run(bad)$status,'not_applicable')
   bad<-d;bad$y<-1-d$y;rev<-run(bad);expect_equal(rev$coefficient_table$estimate,-run()$coefficient_table$estimate,tolerance=1e-9)
  }else{
   bad<-d;bad$exposure[1]<-0;expect_error(run(bad))
   bad<-d;bad$exposure<-bad$exposure*10;rr<-run(bad)
   expect_equal(rr$coefficient_table$estimate[-1],run()$coefficient_table$estimate[-1],tolerance=1e-9)
   expect_equal(rr$coefficient_table$estimate[1],run()$coefficient_table$estimate[1]-log(10),tolerance=1e-9)
   expect_equal(rr$fitted,run()$fitted,tolerance=1e-9)
  }
 }
})
test_that('logistic profile endpoints match independent joint likelihood roots',{
 refs<-jsonlite::read_json(test_path("fixtures","glm-oracle.json"))$cases
 limits<-jsonlite::read_json(test_path("fixtures","glm-profile-oracle.json"))$cases
 for(ref in refs){
  if(ref$kind!='logistic_binary')next
  e<-.glm_env(ref$kind);d<-as.data.frame(do.call(rbind,lapply(ref$x,unlist)));names(d)<-paste0('x',seq_len(ncol(d)));cols<-names(d)
  d$id<-seq_len(nrow(d));d$y<-unlist(ref$y);r<-e$analyse(d,cols,TRUE)
  target<-limits[[which(vapply(limits,function(z)z$name,'')==ref$name)]]
  expected<-do.call(rbind,lapply(target$limits,function(z)as.numeric(unlist(z))))
  expect_equal(r$coefficient_table$lower[-1],expected[,1],tolerance=1e-8)
  expect_equal(r$coefficient_table$upper[-1],expected[,2],tolerance=1e-8)
  expect_true(is.na(r$coefficient_table$ratio[1]))
  expect_equal(r$coefficient_table$ratio_lower[-1],exp(expected[,1]),tolerance=1e-8)
  expect_equal(r$coefficient_table$ratio_upper[-1],exp(expected[,2]),tolerance=1e-8)
 }
})
test_that('NB2 estimates re-estimate dispersion in the null and full likelihoods',{
 refs<-jsonlite::read_json(test_path("fixtures","glm-nb-oracle.json"))$cases;e<-.glm_env('negative_binomial')
 for(ref in refs){
  d<-as.data.frame(do.call(rbind,lapply(ref$x,unlist)));names(d)<-paste0('x',seq_len(ncol(d)));cols<-names(d)
  d$id<-seq_len(nrow(d));d$y<-unlist(ref$y);d$exposure<-unlist(ref$exposure);r<-e$analyse(d,cols,TRUE)
  expect_identical(r$status,'ok');expect_equal(r$coefficient_table$estimate,as.numeric(unlist(ref$beta)),tolerance=1e-6)
  expect_equal(r$coefficient_table$std_error,as.numeric(unlist(ref$se)),tolerance=1e-6)
  expect_equal(r$theta,as.numeric(ref$theta),tolerance=1e-6);expect_equal(r$omnibus$statistic,as.numeric(ref$LR),tolerance=1e-8)
  expect_equal(r$p_value,as.numeric(ref$p_value),tolerance=1e-8);expect_equal(r$fitted,as.numeric(unlist(ref$fitted)),tolerance=1e-6)
 }
 d<-e$example_data;d$y<-rep(1,nrow(d));r<-e$analyse(d,c('x1','x2'),TRUE)
 expect_identical(r$status,'not_applicable');expect_null(r$p_value)
})
test_that('ordinal and multinomial likelihood fits match independent 60-digit references',{
 refs<-jsonlite::read_json(test_path("fixtures","glm-categorical-oracle.json"))$cases
 for(ref in refs){
  e<-.glm_env(ref$kind);ordinal<-ref$kind=='ordinal_proportional_odds'
  r<-if(ordinal)e$analyse(e$example_data,'x1',TRUE,TRUE) else e$analyse(e$example_data,'x1',TRUE,reference_category='1')
  expect_identical(r$status,'ok');beta<-r$coefficient_table$estimate;se<-r$coefficient_table$std_error
  if(ordinal){beta<-c(beta,r$thresholds$estimate);se<-c(se,r$thresholds$std_error)}
  # MASS/nnet use double-precision numerical optimization/Hessians.
  expect_equal(beta,as.numeric(unlist(ref$beta)),tolerance=1e-5)
  expect_equal(se,as.numeric(unlist(ref$se)),tolerance=2e-5)
  expect_equal(r$omnibus$statistic,as.numeric(ref$LR),tolerance=1e-8)
  expect_equal(r$p_value,as.numeric(ref$p_value),tolerance=1e-8)
  target<-do.call(rbind,lapply(ref$fitted,function(z)as.numeric(unlist(z))))
  expect_equal(unname(r$fitted[c(1,21,41),]),target,tolerance=1e-6)
  run<-function(d=e$example_data,...)if(ordinal)e$analyse(d,'x1',TRUE,TRUE,...)else e$analyse(d,'x1',TRUE,reference_category='1',...)
  expect_error(e$analyse(e$example_data,'x1'));expect_error(run(alternative='less'));expect_error(run(conf.level=1))
  d<-e$example_data;d$id[2]<-d$id[1];expect_error(run(d))
  d<-e$example_data;d$x1[1]<-NA;expect_error(run(d));rr<-run(d,complete_cases=TRUE)
  expect_equal(rr$n_excluded,1);expect_equal(rr$p_value,run(e$example_data[-1,])$p_value)
  d<-e$example_data;d$x1<-1;expect_identical(run(d)$status,'not_applicable')
  d<-e$example_data;d$y<-factor(d$y,levels=c('1','2','3','unused'),ordered=ordinal);expect_error(run(d))
  d<-e$example_data;d$x1<-d$x1*100+100;rr<-run(d)
  expect_equal(rr$p_value,r$p_value,tolerance=1e-8)
  expect_equal(rr$coefficient_table$estimate[if(ordinal)1 else c(2,4)],r$coefficient_table$estimate[if(ordinal)1 else c(2,4)]/100,tolerance=1e-5)
  if(ordinal){
   expect_error(e$analyse(e$example_data,'x1',TRUE,FALSE))
   d<-e$example_data;d$y<-factor(d$y,levels=rev(levels(d$y)),ordered=TRUE)
   expect_equal(run(d)$coefficient_table$estimate,-r$coefficient_table$estimate,tolerance=1e-5)
  }else{
   rr<-e$analyse(e$example_data,'x1',TRUE,reference_category='3')
   expect_equal(rr$p_value,r$p_value,tolerance=1e-8);expect_equal(rr$reference_category,'3')
   expect_error(e$analyse(e$example_data,'x1',TRUE,reference_category='unknown'))
  }
 }
})

test_that('regression advisor exports the selected families and NB override in EN and DE',{
 variants<-c(binary='regression_logistic_binary',ordinal='regression_ordinal_proportional_odds',nominal='regression_multinomial',count='regression_poisson')
 for(lang in c('en','de'))for(type in names(variants))for(missing in c('none','some')){
  a<-c(q.aim='model',q.stage='analysis',q.outcome_type=type,q.unit='individuals',q.records_per_individual='one',q.n_predictors='several',q.missingness=missing)
  s<-advisor_replay(a,language=lang);r<-advisor_recommendation(s)
  expect_identical(r$preferred$analysis_variant_id,unname(variants[type]));expect_identical(r$r_implementation$template_status,'available')
  code<-advisor_r_code(s);expect_match(code,r$recommendation_hash,fixed=TRUE)
  e<-new.env(parent=baseenv());invisible(capture.output(eval(parse(text=code),e)));expect_identical(e$example_result$status,'ok')
  expect_true(is.finite(e$example_result$p_value));expect_error(advisor_r_code(s,alternative='less'))
  if(type=='count'){
   ss<-advisor_override(s,'regression_negative_binomial');rr<-advisor_recommendation(ss)
   expect_identical(rr$preferred$analysis_variant_id,'regression_negative_binomial');expect_identical(rr$r_implementation$template_status,'available')
   e<-new.env(parent=baseenv());invisible(capture.output(eval(parse(text=advisor_r_code(ss)),e)))
   expect_identical(e$example_result$status,'ok');expect_gt(e$example_result$theta,0)
  }
 }
})
