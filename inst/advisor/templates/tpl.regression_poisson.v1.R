# Prespecified numeric predictors, independent observations and explicit event/exposure coding.
# Asymptotic likelihood-ratio omnibus test and marginal Wald coefficient intervals.
# No automatic selection, robust/penalized fallback, or causal interpretation.
analyse <- function(data,predictors,independent_units=FALSE,alternative="two.sided",conf.level=.95,complete_cases=FALSE) {
 stopifnot(isTRUE(independent_units),identical(alternative,"two.sided"),is.data.frame(data),
  !anyDuplicated(names(data)),is.character(predictors),length(predictors)>=1L,length(predictors)<=10L,
  !anyNA(predictors),!anyDuplicated(predictors),all(nzchar(predictors)),
  !any(predictors%in%c('id','y','exposure')),all(c('id','y',predictors)%in%names(data)),
  is.atomic(data$id),is.null(dim(data$id)),!anyNA(data$id),!anyDuplicated(data$id),all(nzchar(as.character(data$id))),
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1,
  is.logical(complete_cases),length(complete_cases)==1L,!is.na(complete_cases),nrow(data)<=10000L)
 values<-data[c('y',predictors ,'exposure')]
 stopifnot(all(vapply(values,function(x)is.numeric(x)&&is.null(dim(x)),logical(1))),!any(is.infinite(as.matrix(values))))
 stopifnot('exposure'%in%names(data),all(is.na(values$y)|(values$y>=0&values$y==floor(values$y))),all(is.na(values$exposure)|values$exposure>0))
 keep<-stats::complete.cases(values)
 if(!all(keep)&&!complete_cases)stop('Incomplete observations require an explicit complete-case choice.')
 values<-values[keep,,drop=FALSE];n<-nrow(values);q<-length(predictors)
 stopifnot(n>q+2L)
 unavailable<-function(reason,warnings=character())list(status='not_applicable',reason=reason,warnings=warnings,p_value=NULL,
  confidence_interval=NULL,n_total=nrow(data),n_used=n,n_excluded=sum(!keep))
 y<-values$y; if(sum(y)==0)return(unavailable('At least one observed event is required.'))
 x<-as.matrix(values[predictors]);center<-colMeans(x);z<-sweep(x,2L,center,'-')
 scale<-apply(abs(z),2L,max)
 if(any(!is.finite(z))||any(!is.finite(scale))||any(scale<=0))return(unavailable('Predictor variation is absent or outside the numerical range.'))
 z<-sweep(z,2L,scale,'/');X<-cbind(1,z)
 if(qr(X)$rank<ncol(X)||rcond(crossprod(X))<1e-12)return(unavailable('The prespecified design is singular or numerically unresolved.'))
 offset<-log(values$exposure);family<-stats::poisson(link="log");warnings<-character()
 fitted<-tryCatch(withCallingHandlers({
  full<-stats::glm.fit(X,y,offset=offset,family=family,control=stats::glm.control(epsilon=1e-12,maxit=100),singular.ok=FALSE)
  null<-stats::glm.fit(matrix(1,n,1),y,offset=offset,family=family,control=stats::glm.control(epsilon=1e-12,maxit=100),singular.ok=FALSE)
  list(full=full,null=null)
 },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
 if(inherits(fitted,'error'))return(unavailable(conditionMessage(fitted),warnings))
 f<-fitted$full;null<-fitted$null
 if(length(warnings)||!isTRUE(f$converged)||!isTRUE(null$converged)||isTRUE(f$boundary)||
  any(!is.finite(f$coefficients))||any(!is.finite(f$weights))||any(f$weights<=0))
  return(unavailable('Convergence, boundary or working-weight diagnostics prevent ordinary maximum-likelihood inference.',warnings))
 if(any(f$fitted.values<1e-8))return(unavailable('Near-zero fitted means require boundary/sparsity review; ordinary Wald inference is refused.'))
 # glm.fit stores working weights from the preceding iteration. Recompute the
 # information at the final fitted means for coefficient standard errors.
 final_weights<-family$mu.eta(f$linear.predictors)^2/family$variance(f$fitted.values)
 info<-crossprod(X,X*final_weights)
 if(!is.finite(rcond(info))||rcond(info)<1e-12)return(unavailable('The fitted information matrix is numerically unresolved.'))
 # Transform standardized coefficients and their covariance back to original units.
 A<-diag(c(1,1/scale));A[1,-1]<--center/scale
 beta<-as.vector(A%*%f$coefficients);V<-A%*%solve(info)%*%t(A);se<-sqrt(diag(V))
 if(any(!is.finite(beta))||any(!is.finite(se))||any(se<=0))return(unavailable('Original-unit coefficient inference is outside the numerical range.'))
 statistic<-beta/se;critical<-stats::qnorm((1-conf.level)/2,lower.tail=FALSE)
 lo<-beta-critical*se;hi<-beta+critical*se
 tab<-data.frame(term=c('(Intercept)',predictors),estimate=beta,std_error=se,z=statistic,
  p_value=2*stats::pnorm(abs(statistic),lower.tail=FALSE),lower=lo,upper=hi,
  ratio=exp(beta),ratio_lower=exp(lo),ratio_upper=exp(hi),stringsAsFactors=FALSE)
 LR<-null$deviance-f$deviance
 if(LR<0&&abs(LR)<=64*.Machine$double.eps*max(1,null$deviance,f$deviance))LR<-0
 if(!is.finite(LR)||LR<0)return(unavailable('Nested-model likelihood difference is unresolved.'))
 p<-stats::pchisq(LR,q,lower.tail=FALSE)
 list(status='ok',method='Poisson rate regression with log exposure offset',coefficient_table=tab,interval_method='Marginal asymptotic Wald intervals on the link scale, exponentiated for ratios; not profile-likelihood intervals.',
  p_value=p,omnibus=list(statistic=LR,df=q,null='All prespecified slope coefficients equal zero',method='Asymptotic likelihood-ratio chi-square'),
  n_total=nrow(data),n_used=n,n_excluded=sum(!keep),predictors=predictors,confidence_level=conf.level,
  fitted=as.numeric(f$fitted.values),residual_deviance=f$deviance,residual_df=f$df.residual,
  pearson_dispersion=sum((y-f$fitted.values)^2/family$variance(f$fitted.values))/f$df.residual,
  alpha=1-conf.level,event_or_exposure='y is the nonnegative integer event count; exposure is positive time/opportunity in a common prespecified unit. Ratios are rates per exposure unit.',assumptions='Independent Poisson counts, log exposure offset, correct log-linear mean and conditional variance equal to the mean. Pearson dispersion is descriptive and does not trigger a model switch. Complete-case inference needs justified missingness assumptions.',
  limitations='One to ten prespecified numeric predictors, independent observations, no weights or automatic variable selection. Wald coefficient intervals are marginal and not multiplicity-adjusted. Numerical guards may refuse extreme but finite fits; they do not replace design or model checks. No effect is interpreted as causal without design support.')
}
example_data<-data.frame(id=1:18,x1=rep(c(-1,0,1),6),x2=rep(c(0,1),each=9),exposure=rep(c(1,2,3),6),y=c(1,3,8,2,4,7,0,5,9,2,6,11,1,4,13,3,7,10))
example_result<-analyse(example_data,c('x1','x2'),TRUE,alternative,confidence_level)
print(example_result)
