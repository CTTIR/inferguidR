e<-new.env(parent=baseenv());e$alternative<-'two.sided';e$confidence_level<-.95
invisible(capture.output(sys.source(system.file("advisor","templates","tpl.anova_factorial.v1.R",package="inferguidR"),e)))
.factorial_data<-function(ref) {
 m<-matrix(unlist(ref$groups),ncol=length(ref$levels),byrow=TRUE)
 d<-data.frame(id=seq_len(nrow(m)),y=as.numeric(unlist(ref$y)))
 for(j in seq_len(ncol(m))) d[[LETTERS[j]]]<-factor(m[,j],levels=seq_len(ref$levels[[j]])-1L)
 d
}
.fraction_number<-function(x) {z<-as.numeric(strsplit(x,'/',fixed=TRUE)[[1]]);if(length(z)==1L)z else z[1]/z[2]}
test_that('type-III factorial hypotheses match exact rational matrix references',{
 refs<-jsonlite::read_json(test_path("fixtures","factorial-rational.json"))$cases
 for(ref in refs) {
  d<-.factorial_data(ref);cols<-LETTERS[seq_along(ref$levels)]
  r<-e$analyse(d,cols,TRUE,TRUE)
  expect_identical(r$status,'ok');expect_equal(r$residual_df,ref$residual_df)
  expect_equal(r$residual_sum_squares,as.numeric(ref$residual_sum_squares),tolerance=1e-13)
  for(t in ref$terms) {
   a<-r$terms[r$terms$term==t$term,,drop=FALSE]
   expect_equal(nrow(a),1);expect_equal(a$df,t$df)
   expect_equal(a$sum_squares,.fraction_number(t$sum_squares),tolerance=2e-12)
   expect_equal(a$F,as.numeric(t$F),tolerance=2e-12)
   expect_equal(a$p_value,as.numeric(t$p_value),tolerance=2e-12)
   expect_equal(a$partial_eta_squared,as.numeric(t$partial_eta_squared),tolerance=2e-12)
  }
  expect_equal(r$terms$p_holm,p.adjust(r$terms$p_value,'holm'))
  expect_null(r$confidence_interval);expect_equal(r$n_excluded,0)
  reordered<-d[rev(seq_len(nrow(d))),,drop=FALSE]
  for(col in cols)reordered[[col]]<-factor(as.character(reordered[[col]]),levels=rev(levels(d[[col]])))
  rr<-e$analyse(reordered,cols,TRUE,TRUE)
  expect_equal(rr$terms$F,r$terms$F,tolerance=2e-11)
 }
})
test_that('factorial design, missingness and degeneracy are explicit',{
 d<-e$example_data;run<-function(data=d,cols=c('A','B'),...)e$analyse(data,cols,TRUE,TRUE,...)
 expect_error(e$analyse(d,c('A','B')));expect_error(e$analyse(d,c('A','B'),TRUE,FALSE))
 expect_error(run(alternative='greater'));expect_error(run(conf.level=1))
 bad<-d;bad$id[2]<-bad$id[1];expect_error(run(bad))
 bad<-d;bad$id[2]<-NA;expect_error(run(bad))
 bad<-d;bad$A<-as.character(bad$A);expect_error(run(bad))
 bad<-d;bad$A<-ordered(bad$A);expect_error(run(bad))
 bad<-d;bad$A<-factor(bad$A,levels=c(levels(bad$A),'unused'));expect_error(run(bad))
 bad<-d;bad$y[1]<-Inf;expect_error(run(bad,complete_cases=TRUE))
 bad<-d;bad$y[1]<-NA;expect_error(run(bad));r<-run(bad,complete_cases=TRUE)
 expect_equal(r$n_excluded,1);expect_equal(r$n_used,nrow(d)-1)
 # Dropping two observations leaves one declared cell with only one replicate.
 bad$y[7]<-NA;expect_error(run(bad,complete_cases=TRUE))
 exact<-d;exact$y<-2*as.integer(d$A)+3*as.integer(d$B)
 r<-run(exact);expect_identical(r$status,'not_applicable');expect_null(r$p_value)
 expect_error(run(cols='A'));expect_error(run(cols=c('A','A')));expect_error(run(cols=c('A','y')))
 changed<-d;changed$y<-3*d$y+1e12
 expect_equal(run(changed)$terms$F,run()$terms$F,tolerance=2e-11)
 # Factor names matching internal implementation names must not change labels.
 unusual<-d;names(unusual)[names(unusual)=='A']<-'.factor2'
 labels<-e$analyse(unusual,c('.factor2','B'),TRUE,TRUE)$terms$term
 expect_equal(labels,c('.factor2','B','.factor2:B'))
})

test_that('factorial advisor exports prespecified omnibus models in EN and DE',{
 base<-c(q.aim='compare',q.stage='analysis',q.outcome_type='continuous',q.unit='individuals',
  q.structure='independent',q.independence_check='no',q.n_factors='two_or_more',q.missingness='none')
 for(lang in c('en','de'))for(missing in c('none','some')) {
  ans<-base;ans['q.missingness']<-missing;s<-advisor_replay(ans,language=lang)
  r<-advisor_recommendation(s)
  expect_identical(r$preferred$analysis_variant_id,'anova_factorial')
  expect_identical(r$r_implementation$template_status,'available')
  expect_identical(r$status,if(missing=='some')'conditional' else 'recommended')
  code<-advisor_r_code(s);expect_match(code,r$recommendation_hash,fixed=TRUE)
  env<-new.env(parent=baseenv());invisible(capture.output(eval(parse(text=code),env)))
  expect_identical(env$example_result$status,'ok')
  expect_equal(env$example_result$terms$F,c(72,37.5,1.5),tolerance=1e-12)
  expect_equal(env$example_result$terms$term,c('A','B','A:B'))
  expect_match(env$example_result$method,'type-III',fixed=TRUE)
  expect_error(advisor_r_code(s,alternative='greater'),'omnibus')
  expect_error(advisor_r_code(s,alternative='less'),'omnibus')
 }
})
