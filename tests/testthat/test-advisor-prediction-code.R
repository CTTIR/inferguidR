.prediction_env<-function(kind){e<-new.env(parent=baseenv());e$alternative<-'two.sided';e$confidence_level<-.95;invisible(capture.output(sys.source(system.file("advisor","templates",paste0("tpl.prediction_model_",kind,".v1.R"),package="inferguidR"),e)));e}
test_that('prediction validation agrees with independent fixed-resample model fits',{
 refs<-jsonlite::read_json(test_path("fixtures","prediction-oracle.json"))$cases
 for(ref in refs){
  e<-.prediction_env(ref$kind);d<-data.frame(id=seq_along(ref$y),y=unlist(ref$y),x1=unlist(ref$x1),x2=unlist(ref$x2));indices<-matrix(unlist(ref$resamples),nrow=nrow(d))
  r<-e$analyse(d,c('x1','x2'),TRUE,B=20,resamples=indices)
  expect_identical(r$status,'ok');expect_equal(names(r$apparent),unlist(ref$metric_names))
  expect_equal(unname(r$apparent),unlist(ref$apparent),tolerance=1e-8)
  expect_equal(unname(r$corrected),unlist(ref$corrected),tolerance=1e-8)
  expect_equal(unname(r$monte_carlo_se),unlist(ref$mcse),tolerance=1e-8)
  expect_equal(unname(r$bootstrap_differences),do.call(rbind,lapply(ref$differences,unlist)),tolerance=1e-8)
  expect_equal(e$predict_new(r,data.frame(x1=c(-2,0,3),x2=c(0,1,0))),unlist(ref$new_predictions),tolerance=1e-8)
  expect_equal(r$failed_resamples,0);expect_null(r$confidence_interval)
  expect_equal(r$corrected,r$apparent+colMeans(r$bootstrap_differences))
 }
})
test_that('prediction bootstrap preserves RNG state and never drops failed fits',{
 for(kind in c('linear','logistic')){
  e<-.prediction_env(kind);d<-e$example_data;run<-function(x=d,...)e$analyse(x,c('x1','x2'),TRUE,B=20,...)
  set.seed(782);before<-.Random.seed;kind_before<-RNGkind();r<-run()
  expect_identical(.Random.seed,before);expect_identical(RNGkind(),kind_before)
  again<-run();expect_identical(again$resample_indices,r$resample_indices);expect_equal(again$corrected,r$corrected)
  bad_indices<-r$resample_indices;bad_indices[,1]<-1L;bad<-run(resamples=bad_indices)
  expect_identical(bad$status,'validation_incomplete');expect_null(bad$corrected);expect_gte(bad$failed_resamples,1)
  expect_identical(.Random.seed,before);expect_error(e$predict_new(bad,d))
  expect_error(e$analyse(d,c('x1','x2')));expect_error(run(alternative='greater'));expect_error(run(conf.level=1))
  dd<-d;dd$id[2]<-dd$id[1];expect_error(run(dd))
  dd<-d;dd$y[1]<-NA;expect_error(run(dd));rr<-run(dd,complete_cases=TRUE)
  expect_equal(rr$n_excluded,1);expect_equal(rr$n_used,nrow(d)-1)
  dd<-d;dd$x1[1]<-Inf;expect_error(run(dd))
  dd<-d;dd$x2<-dd$x1;expect_error(run(dd))
  expect_error(run(resamples=matrix(0,nrow(d),20)))
  expect_error(e$predict_new(r,data.frame(x1=NA_real_,x2=1)))
  if(kind=='logistic'){
   expect_error(e$measure_prediction(rep(1,10),rep(.5,10)))
   expect_equal(e$measure_prediction(c(0,1,0,1,0,1),c(.2,.2,.5,.5,.8,.8))['auc'],c(auc=.5))
  }
 }
})

test_that('prediction advisor exports internal validation without hypothesis-power substitution',{
 for(lang in c('en','de'))for(type in c('continuous','binary'))for(missing in c('none','some')){
  a<-c(q.aim='predict',q.stage='analysis',q.outcome_type=type,q.unit='individuals',q.records_per_individual='one',q.missingness=missing)
  s<-advisor_replay(a,language=lang);r<-advisor_recommendation(s)
  expect_identical(r$preferred$analysis_variant_id,if(type=='continuous')'prediction_model_linear' else 'prediction_model_logistic')
  expect_identical(r$r_implementation$template_status,'available');expect_identical(r$status,'conditional')
  code<-advisor_r_code(s);expect_match(code,r$recommendation_hash,fixed=TRUE)
  e<-new.env(parent=baseenv());invisible(capture.output(eval(parse(text=code),e)))
  expect_identical(e$example_result$status,'ok');expect_equal(e$example_result$B,200L)
  expect_true(all(is.finite(e$example_result$corrected)));expect_null(e$example_result$confidence_interval)
  expect_error(advisor_r_code(s,alternative='greater'))
  expect_error(advisor_r_code(s,alternative='two.sided'),'no hypothesis alternative')
  expect_false(grepl('alternative <-',code,fixed=TRUE))
 }
})
