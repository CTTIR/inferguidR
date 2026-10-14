# Prediction development with prespecified predictors and bootstrap optimism correction.
# No variable selection, cutoff optimization, or external-validation claim.
fit_prediction <- function(x,y) {
 center<-colMeans(x);z<-sweep(x,2L,center,'-');scale<-apply(abs(z),2L,max)
 if(any(!is.finite(z))||any(!is.finite(scale))||any(scale<=0))stop('A predictor is constant or numerically unresolved in this resample.')
 z<-sweep(z,2L,scale,'/');X<-cbind(1,z)
 if(qr(X)$rank<ncol(X)||rcond(crossprod(X))<1e-12)stop('A resampled design is singular or ill-conditioned.')
 if(length(unique(y))!=2L)stop('A resample contains only one outcome class.')
 alerts<-character();fit<-withCallingHandlers(stats::glm.fit(X,y,family=stats::binomial(),control=stats::glm.control(epsilon=1e-11,maxit=100)),warning=function(w){alerts<<-c(alerts,conditionMessage(w));invokeRestart('muffleWarning')})
 beta<-fit$coefficients
 if(length(alerts)||!isTRUE(fit$converged)||any(!is.finite(beta))||any(fit$fitted.values<1e-8|fit$fitted.values>1-1e-8))stop('A logistic prediction fit is unresolved or near a probability boundary.')
 list(coefficients=beta,center=center,scale=scale)
}
predict_model <- function(model,x) {
 eta<-as.vector(cbind(1,sweep(sweep(x,2L,model$center,'-'),2L,model$scale,'/'))%*%model$coefficients)
 if(any(!is.finite(eta)))stop('New-data predictions are outside the numerical range.')
 stats::plogis(eta)
}
measure_prediction <- function(y,prediction) {
 if(length(unique(y))!=2L||any(prediction<=0|prediction>=1))stop('Both outcome classes and interior predicted probabilities are required.')
 lp<-stats::qlogis(prediction);alerts<-character()
 calibration<-withCallingHandlers(stats::glm.fit(cbind(1,lp),y,family=stats::binomial(),control=stats::glm.control(epsilon=1e-11,maxit=100)),warning=function(w){alerts<<-c(alerts,conditionMessage(w));invokeRestart('muffleWarning')})
 if(length(alerts)||!isTRUE(calibration$converged)||calibration$rank!=2L||any(!is.finite(calibration$coefficients))||
  any(calibration$fitted.values<1e-8|calibration$fitted.values>1-1e-8))stop('Calibration fit is unresolved; no selective validation metric is reported.')
 positive<-prediction[y==1];negative<-prediction[y==0]
 # Fixed higher-score-is-event orientation; half credit for ties.
 auc<-sum(outer(positive,negative,'>')+.5*outer(positive,negative,'=='))/(length(positive)*length(negative))
 c(brier=mean((y-prediction)^2),log_loss=-mean(y*log(prediction)+(1-y)*log1p(-prediction)),auc=auc,
  calibration_intercept=unname(calibration$coefficients[1]),calibration_slope=unname(calibration$coefficients[2]))
}
analyse <- function(data,predictors,independent_units=FALSE,alternative=NULL,conf.level=.95,
 B=200L,seed=104729L,complete_cases=FALSE,resamples=NULL) {
 stopifnot(isTRUE(independent_units),is.null(alternative),is.data.frame(data),!anyDuplicated(names(data)),
  is.character(predictors),length(predictors)>=1L,length(predictors)<=10L,!anyNA(predictors),!anyDuplicated(predictors),
  !any(predictors%in%c('id','y')),all(c('id','y',predictors)%in%names(data)),is.atomic(data$id),is.null(dim(data$id)),
  !anyNA(data$id),!anyDuplicated(data$id),all(nzchar(as.character(data$id))),
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1,
  is.numeric(B),length(B)==1L,is.finite(B),B==floor(B),B>=20L,B<=2000L,
  is.numeric(seed),length(seed)==1L,is.finite(seed),seed==floor(seed),seed>=0,seed<=.Machine$integer.max,
  is.logical(complete_cases),length(complete_cases)==1L,!is.na(complete_cases),nrow(data)<=2000L)
 values<-data[c('y',predictors)];stopifnot(all(vapply(values,function(x)is.numeric(x)&&is.null(dim(x)),logical(1))),!any(is.infinite(as.matrix(values))))
 stopifnot(all(is.na(values$y)|values$y%in%c(0,1)))
 keep<-stats::complete.cases(values)
 if(!all(keep)&&!complete_cases)stop('Incomplete records require an explicit complete-case choice; imputation is not supplied.')
 values<-values[keep,,drop=FALSE];n<-nrow(values);stopifnot(n>=20L,n>length(predictors)+2L)
 x<-as.matrix(values[predictors]);y<-values$y
 final<-fit_prediction(x,y);apparent<-measure_prediction(y,predict_model(final,x))
 supplied_resamples<-!is.null(resamples)
 if(is.null(resamples)) {
  kind<-RNGkind();had_seed<-exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)
  if(had_seed)saved_seed<-get('.Random.seed',envir=.GlobalEnv,inherits=FALSE)
  on.exit({do.call(RNGkind,as.list(kind));if(had_seed)assign('.Random.seed',saved_seed,envir=.GlobalEnv) else if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE))rm('.Random.seed',envir=.GlobalEnv)},add=TRUE)
  RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(as.integer(seed))
  resamples<-matrix(sample.int(n,n*B,replace=TRUE),nrow=n,ncol=B)
 } else stopifnot(is.matrix(resamples),is.numeric(resamples),identical(dim(resamples),c(as.integer(n),as.integer(B))),
  all(is.finite(resamples)),all(resamples==floor(resamples)),all(resamples>=1&resamples<=n))
 bootstrap_models<-vector('list',B)
 optimism<-matrix(NA_real_,nrow=B,ncol=length(apparent),dimnames=list(NULL,names(apparent)));failures<-character()
 # Every resample refits preprocessing and the complete prespecified model.
 for(b in seq_len(B)) {
  index<-resamples[,b]
  value<-tryCatch({model<-fit_prediction(x[index,,drop=FALSE],y[index]);
   train<-measure_prediction(y[index],predict_model(model,x[index,,drop=FALSE]));test<-measure_prediction(y,predict_model(model,x));
   bootstrap_models[[b]]<-model;test-train},error=function(e)e)
  if(inherits(value,'error'))failures<-c(failures,paste0('Resample ',b,': ',conditionMessage(value))) else optimism[b,]<-value
 }
 # Failed resamples are never silently discarded or replaced by easier samples.
 if(length(failures))return(list(status='validation_incomplete',method='Prespecified logistic risk model with bootstrap optimism correction',apparent=apparent,
  corrected=NULL,failed_resamples=length(failures),failures=failures,B=B,seed=if(supplied_resamples)NULL else seed,model=final,
  n_total=nrow(data),n_used=n,n_excluded=sum(!keep),reason='At least one bootstrap fit/metric failed; no selectively averaged validation estimate is reported.'))
 correction<-colMeans(optimism);corrected<-apparent+correction
 list(status='ok',method='Prespecified logistic risk model with bootstrap optimism correction',predictors=predictors,model=final,apparent=apparent,corrected=corrected,
  mean_test_minus_training=correction,monte_carlo_se=apply(optimism,2L,stats::sd)/sqrt(B),B=B,seed=if(supplied_resamples)NULL else seed,resampling_source=if(supplied_resamples)'supplied indices; ordinary independent bootstrap design must be justified' else 'generated ordinary independent bootstrap',bootstrap_models=bootstrap_models,failed_resamples=0L,
  resample_indices=resamples,bootstrap_differences=optimism,n_total=nrow(data),n_used=n,n_excluded=sum(!keep),
  confidence_interval=NULL,interval_status='not_computed',requested_confidence_level=conf.level,
  validation='Ordinary individual bootstrap optimism correction: apparent + mean(original-development performance minus resample-training performance). All preprocessing and model fitting occur inside each resample.',
  interpretation='y=1 is the event. Brier/log loss are lower-is-better; AUC is higher-is-better with fixed event orientation. Calibration fits outcome on predicted log odds; ideal intercept=0 and slope=1. Calibration correction estimates are reported without clipping or a sampling interval; predict_new(interval="risk") separately quantifies risk-estimation uncertainty.',
  limitations='Internal validation estimates, not external validation or deployment approval. Monte Carlo SE describes resampling simulation variability, not sampling uncertainty. The sample-size floor is not a prediction-model sample-size calculation. Numerical correction estimates are not clipped to ideal metric ranges. No cutoff is selected and no causal interpretation is supplied.')
}
predict_new <- function(result,newdata,interval=c('none','risk'),conf.level=result$requested_confidence_level) {
 interval<-match.arg(interval)
 stopifnot(identical(result$status,'ok'),is.data.frame(newdata),nrow(newdata)>0L,!anyDuplicated(names(newdata)),all(result$predictors%in%names(newdata)),
  all(vapply(newdata[result$predictors],function(x)is.numeric(x)&&is.null(dim(x)),logical(1))),
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1)
 x<-as.matrix(newdata[result$predictors]);stopifnot(all(is.finite(x)));estimate<-predict_model(result$model,x)
 if(interval=='none')return(estimate)
 tail_probability<-(1-conf.level)/2
 if(result$B*tail_probability<5-1e-12)return(list(status='insufficient_tail_resolution',estimate=estimate,
  required_resamples=ceiling(5/tail_probability),reason='At least five expected bootstrap draws per interval tail are required by this implementation; increase B. This floor does not guarantee interval accuracy.'))
 stopifnot(length(result$bootstrap_models)==result$B,all(vapply(result$bootstrap_models,is.list,logical(1))))
 draws<-vapply(result$bootstrap_models,function(model)predict_model(model,x),numeric(nrow(x)))
 if(is.null(dim(draws)))draws<-matrix(draws,nrow=nrow(x))
 bounds<-t(apply(draws,1L,stats::quantile,probs=c(tail_probability,1-tail_probability),type=6,names=FALSE))
 list(status='ok',estimate=estimate,lower=bounds[,1],upper=bounds[,2],interval='risk',confidence_level=conf.level,B=result$B,
  method='Pointwise ordinary-bootstrap percentile interval for the fitted event risk; R quantile type6.',
  interpretation='Uncertainty in the estimated risk from development-sample estimation, not an interval for the binary outcome, validation performance or bootstrap Monte Carlo precision.',
  assumptions='Independent representative development observations, ordinary individual bootstrap, prespecified logistic model and preprocessing, all bootstrap fits successful. No claim of exact coverage, multiplicity control, transportability or protection against model misspecification.')
}
example_data<-data.frame(id=1:200,x1=rep(-4:5,20),x2=rep(c(0,1),each=100));example_data$y<-as.numeric(((seq_len(200)*37)%%101)/101 < stats::plogis(-.4+.18*example_data$x1+.35*example_data$x2))
example_result<-analyse(example_data,c('x1','x2'),TRUE,conf.level=confidence_level,B=200L)
print(example_result[c('status','method','apparent','corrected','monte_carlo_se','B','n_used')])
