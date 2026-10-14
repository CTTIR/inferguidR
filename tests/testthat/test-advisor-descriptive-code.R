.desc_code_env<-function() {
 e<-new.env();e$confidence_level<-.95
 invisible(capture.output(sys.source(system.file('advisor','templates','tpl.descriptive_summary.v1.R',package='inferguidR'),e)));e
}
test_that('descriptive examples require measurement meaning and retain interval estimands', {
 e<-.desc_code_env()
 r<-e$analyse(c(-1,1),'continuous');h<-1/tan(pi*.025)
 expect_equal(r$confidence_interval,c(-h,h));expect_equal(r$sd,sqrt(2))
 expect_equal(r$quantiles,c(-.5,0,.5));expect_null(r$p_value)
 expect_identical(r$interval_estimand,'population mean')
 for(x in list(1,rep(2,4))) {
  r<-e$analyse(x,'continuous');expect_null(r$confidence_interval);expect_identical(r$interval_status,'not_applicable')
 }
 for(x in list(c(NA,1),c(Inf,1),matrix(1:4,2),'1')) expect_error(e$analyse(x,'continuous'))
 expect_error(e$analyse(1:3));expect_error(e$analyse(1:3,'unknown'))
 for(conf in list(0,1,NA,Inf,c(.9,.95))) expect_error(e$analyse(1:3,'continuous',conf))
})
test_that('binary summaries preserve absent events with declared levels', {
 e<-.desc_code_env()
 for(k in c(0,5)) {
  r<-e$analyse(factor(rep(if(k==0) 'no' else 'yes',5),levels=c('no','yes')),'binary',event='yes')
  expected<-if(k==0) c(0,1-.025^(1/5)) else c(.025^(1/5),1)
  expect_equal(r$confidence_interval,expected);expect_equal(r$events,k);expect_null(r$p_value)
 }
 expect_equal(e$analyse(c(TRUE,FALSE,TRUE),'binary',event=TRUE)$proportion,2/3)
 expect_error(e$analyse(c('yes','no'),'binary',event='yes'))
 expect_error(e$analyse(factor(c('yes','no')),'binary',event='typo'))
 expect_error(e$analyse(factor('yes'),'binary',event='yes'))
 expect_error(e$analyse(c(TRUE,FALSE),'binary'))
})
test_that('ordered summaries use observed categories and marginal category intervals', {
 e<-.desc_code_env();x<-ordered(c('low','low','mid','high'),levels=c('low','mid','high','unused'))
 r<-e$analyse(x,'ordinal');expect_identical(r$ordered_quantile_categories,c('low','low','mid'))
 expect_identical(r$quantile_type,1L);expect_equal(r$categories$count,c(2,1,1,0))
 expect_equal(r$categories$upper[4],1-.025^(1/4));expect_identical(r$interval_scope,'marginal, not simultaneous')
 expect_null(r$p_value);expect_error(e$analyse(factor(x,ordered=FALSE),'ordinal'))
 nominal<-e$analyse(c('a','a','b'),'nominal');expect_equal(nominal$categories$count,c(2,1))
 expect_null(nominal$ordered_quantile_categories)
 expect_error(e$analyse(c('a',''),'nominal'));expect_error(e$analyse(c('a',NA),'nominal'))
})
test_that('descriptive advisor export does not invent a hypothesis direction', {
 s<-advisor_replay(c(q.aim='describe'));r<-advisor_recommendation(s)
 expect_identical(r$preferred$analysis_variant_id,'descriptive_summary')
 code<-advisor_r_code(s);expect_false(grepl('alternative <-',code,fixed=TRUE))
 e<-new.env();invisible(capture.output(eval(parse(text=code),e)));expect_null(e$example_result$p_value)
 expect_error(advisor_r_code(s,alternative='two.sided'),'no hypothesis')
 expect_error(advisor_r_code(s,alternative='greater'),'no hypothesis')
})
test_that('descriptive examples retain recorded measurement meaning on no-reference routes', {
 for(type in c('continuous','binary','ordinal','nominal')) {
  s<-advisor_replay(c(q.aim='compare',q.stage='analysis',q.outcome_type=type,
   q.unit='individuals',q.structure='one_sample',q.reference_value='no'))
  expect_identical(advisor_recommendation(s)$preferred$analysis_variant_id,'descriptive_summary')
  code<-advisor_r_code(s);e<-new.env();invisible(capture.output(eval(parse(text=code),e)))
  expect_identical(e$example_result$measurement,type);expect_null(e$example_result$p_value)
 }
})
