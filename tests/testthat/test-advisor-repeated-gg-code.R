e<-new.env(parent=baseenv());e$alternative<-'two.sided';e$confidence_level<-.95
invisible(capture.output(sys.source(system.file("advisor","templates","tpl.anova_repeated_one_way_gg.v1.R",package="inferguidR"),e)))
.gg_fraction<-function(x){z<-as.numeric(strsplit(x,'/',fixed=TRUE)[[1]]);if(length(z)==1)z else z[1]/z[2]}
.gg_data<-function(ref){y<-do.call(rbind,lapply(ref$y,unlist));d<-as.data.frame(y);names(d)<-paste0('time',seq_len(ncol(y)));d$id<-seq_len(nrow(y));d}
test_that('GG corrected repeated ANOVA matches exact rational references',{
 refs<-jsonlite::read_json(test_path("fixtures","repeated-gg-rational.json"))$cases
 for(ref in refs){
  d<-.gg_data(ref);cols<-setdiff(names(d),'id');r<-e$analyse(d,cols,TRUE)
  expect_identical(r$status,'ok');expect_equal(r$statistic,.gg_fraction(ref$F),tolerance=1e-12)
  expect_equal(r$epsilon_gg,.gg_fraction(ref$epsilon),tolerance=1e-12)
  expect_equal(r$df_numerator,.gg_fraction(ref$df1),tolerance=1e-12)
  expect_equal(r$df_denominator,.gg_fraction(ref$df2),tolerance=1e-12)
  expect_equal(r$p_value,as.numeric(ref$p_value),tolerance=1e-12)
  expect_equal(r$partial_eta_squared,.gg_fraction(ref$partial_eta),tolerance=1e-12)
  expect_null(r$confidence_interval);expect_equal(r$n_excluded,0)
  rr<-e$analyse(d[rev(seq_len(nrow(d))),],rev(cols),TRUE)
  expect_equal(rr$statistic,r$statistic,tolerance=1e-12)
  expect_equal(rr$epsilon_gg,r$epsilon_gg,tolerance=1e-12)
  expect_equal(rr$p_value,r$p_value,tolerance=1e-12)
 }
})
test_that('GG inputs require identified independent subjects and explicit exclusions',{
 d<-e$example_data;cols<-c('before','during','after');run<-function(x=d,...)e$analyse(x,cols,TRUE,...)
 expect_error(e$analyse(d,cols));expect_error(run(alternative='less'));expect_error(run(conf.level=1))
 bad<-d;bad$id[2]<-bad$id[1];expect_error(run(bad))
 bad<-d;bad$id[2]<-NA;expect_error(run(bad))
 bad<-d;bad$before[1]<-Inf;expect_error(run(bad,complete_cases=TRUE))
 bad<-d;bad$before<-factor(bad$before);expect_error(run(bad))
 bad<-d;bad$before[1]<-NA;expect_error(run(bad));r<-run(bad,complete_cases=TRUE)
 expect_equal(r$n_used,5);expect_equal(r$n_excluded,1);expect_equal(r$n_total,6)
 expect_equal(r$p_value,run(d[-1,])$p_value)
 bad$after[2:5]<-NA;expect_error(run(bad,complete_cases=TRUE))
 expect_error(e$analyse(d,c('before','during'),TRUE));expect_error(e$analyse(d,c('before','before','after'),TRUE))
 exact<-d;exact[cols]<-lapply(1:3,function(i)rep(i,nrow(d)))
 expect_identical(run(exact)$status,'not_applicable');expect_null(run(exact)$p_value)
 exact[cols]<-1;expect_identical(run(exact)$status,'not_applicable')
 extreme<-d;extreme$before<-1e308;extreme$during<--1e308
 expect_identical(run(extreme)$status,'not_applicable')
 scaled<-d;scaled[cols]<-lapply(d[cols],function(x)x*1e150)
 expect_equal(run(scaled)$p_value,run()$p_value,tolerance=1e-12)
 scaled[cols]<-lapply(d[cols],function(x)x*1e-150)
 expect_equal(run(scaled)$epsilon_gg,run()$epsilon_gg,tolerance=1e-12)
 shifted<-d;shifted[cols]<-lapply(d[cols],function(x)x+seq_len(nrow(d))*1e12)
 expect_equal(run(shifted)$p_value,run()$p_value,tolerance=1e-12)
})

test_that('repeated-measures advisor exports routine GG examples in EN and DE',{
 base<-c(q.aim='compare',q.stage='analysis',q.outcome_type='continuous',q.unit='individuals',
  q.structure='repeated',q.pairing_id='yes',q.estimand_continuous='mean',q.missingness='none')
 for(lang in c('en','de')) {
  s<-advisor_replay(base,language=lang);r<-advisor_recommendation(s)
  expect_identical(r$preferred$analysis_variant_id,'anova_repeated_one_way_gg')
  expect_identical(r$r_implementation$template_status,'available');expect_identical(r$status,'recommended')
  code<-advisor_r_code(s);expect_match(code,r$recommendation_hash,fixed=TRUE)
  env<-new.env(parent=baseenv());invisible(capture.output(eval(parse(text=code),env)))
  expect_identical(env$example_result$status,'ok');expect_equal(env$example_result$statistic,405/19,tolerance=1e-12)
  expect_equal(env$example_result$n_used,6);expect_match(env$example_result$method,'routine Greenhouse-Geisser',fixed=TRUE)
  expect_error(advisor_r_code(s,alternative='greater'),'omnibus');expect_error(advisor_r_code(s,alternative='less'),'omnibus')
  missing<-base;missing['q.missingness']<-'some';rm<-advisor_recommendation(advisor_replay(missing,language=lang))
  expect_identical(rm$preferred$analysis_variant_id,'lmm_repeated_measures')
 }
})
