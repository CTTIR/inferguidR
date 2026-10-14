.uncertainty_env<-function(kind){e<-new.env(parent=baseenv());e$alternative<-'two.sided';e$confidence_level<-.95;invisible(capture.output(sys.source(system.file("advisor","templates",paste0("tpl.prediction_model_",kind,".v1.R"),package="inferguidR"),e)));e}
test_that('prediction sampling intervals match independent references',{
 refs<-jsonlite::read_json(test_path("fixtures","prediction-uncertainty-oracle.json"))$cases
 for(ref in refs){
  e<-.uncertainty_env(ref$kind);d<-data.frame(id=seq_along(ref$y),y=unlist(ref$y),x1=unlist(ref$x1),x2=unlist(ref$x2));new<-data.frame(x1=c(-2,0,3),x2=c(0,1,0))
  indices<-if(ref$kind=='logistic')matrix(unlist(ref$resamples),nrow=nrow(d)) else NULL
  r<-e$analyse(d,c('x1','x2'),TRUE,B=200,resamples=indices)
  expect_identical(r$status,'ok');expect_length(r$bootstrap_models,200)
  if(ref$kind=='linear'){
   for(mode in c('mean','prediction')){z<-e$predict_new(r,new,interval=mode);expect_identical(z$status,'ok');expect_equal(z$lower,unlist(ref[[mode]]$lower),tolerance=1e-9);expect_equal(z$upper,unlist(ref[[mode]]$upper),tolerance=1e-9);expect_equal(z$df,nrow(d)-3)}
   expect_true(all(e$predict_new(r,new,interval='prediction')$lower<e$predict_new(r,new,interval='mean')$lower))
   expect_identical(r$seed,104729L)
   saved_variance<-r$model$residual_variance;r$model$residual_variance<-0;expect_identical(e$predict_new(r,new,interval='prediction')$status,'unavailable');r$model$residual_variance<-saved_variance
  }else{
   z<-e$predict_new(r,new,interval='risk');expect_identical(z$status,'ok');expect_equal(z$lower,unlist(ref$lower),tolerance=1e-8);expect_equal(z$upper,unlist(ref$upper),tolerance=1e-8);expect_null(r$seed)
   expect_match(r$resampling_source,'supplied indices');expect_equal(e$predict_new(r,new[1,,drop=FALSE],interval='risk')$lower,z$lower[1])
   expect_identical(e$predict_new(r,new,interval='risk',conf.level=.99)$status,'insufficient_tail_resolution')
  }
  expect_error(e$predict_new(r,new,conf.level=1));expect_error(e$predict_new(r,new,interval='unknown'))
  r$requested_confidence_level<-.9;expect_equal(e$predict_new(r,new,interval=if(ref$kind=='linear')'mean' else 'risk')$confidence_level,.9)
 }
})
