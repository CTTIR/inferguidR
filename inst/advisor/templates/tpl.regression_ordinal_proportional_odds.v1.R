# Prespecified numeric predictors and explicit category coding; no selection.
# Independent observations; asymptotic likelihood-ratio and marginal Wald inference.
analyse <- function(data,predictors,independent_units=FALSE,proportional_odds_justified=FALSE,alternative="two.sided",conf.level=.95,complete_cases=FALSE) {
 if(!requireNamespace('MASS',quietly=TRUE))stop('This example requires the recommended MASS package.')
 stopifnot(isTRUE(independent_units),identical(alternative,'two.sided'),is.data.frame(data),
  !anyDuplicated(names(data)),is.character(predictors),length(predictors)>=1L,length(predictors)<=10L,
  !anyNA(predictors),!anyDuplicated(predictors),all(nzchar(predictors)),!any(predictors%in%c('id','y')),
  all(c('id','y',predictors)%in%names(data)),is.factor(data$y),is.ordered(data$y),
  nlevels(data$y)>=3L,nlevels(data$y)<=8L,is.atomic(data$id),is.null(dim(data$id)),!anyNA(data$id),
  !anyDuplicated(data$id),all(nzchar(as.character(data$id))),nrow(data)<=10000L,
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1,
  is.logical(complete_cases),length(complete_cases)==1L,!is.na(complete_cases),
  all(vapply(data[predictors],function(x)is.numeric(x)&&is.null(dim(x)),logical(1))),
  !any(is.infinite(as.matrix(data[predictors]))))
 stopifnot(isTRUE(proportional_odds_justified))
 keep<-stats::complete.cases(data[c('y',predictors)])
 if(!all(keep)&&!complete_cases)stop('Incomplete observations require an explicit complete-case choice.')
 d<-data[keep,,drop=FALSE];n<-nrow(d);q<-length(predictors);K<-nlevels(d$y)
 stopifnot(all(table(d$y)>0),n>q+K-1L+2L)
 unavailable<-function(reason,warnings=character())list(status='not_applicable',reason=reason,warnings=warnings,p_value=NULL,
  confidence_interval=NULL,n_total=nrow(data),n_used=n,n_excluded=sum(!keep))
 x<-as.matrix(d[predictors]);center<-colMeans(x);z<-sweep(x,2L,center,'-');scale<-apply(abs(z),2L,max)
 if(any(!is.finite(z))||any(!is.finite(scale))||any(scale<=0))return(unavailable('Predictor variation is absent or numerically unresolved.'))
 z<-sweep(z,2L,scale,'/');X<-cbind(1,z)
 if(qr(X)$rank!=ncol(X)||rcond(crossprod(X))<1e-12)return(unavailable('The prespecified design is singular or ill-conditioned.'))
 frame<-as.data.frame(z);names(frame)<-paste0('x',seq_len(q));frame$y<-d$y
 form<-stats::reformulate(names(frame)[seq_len(q)],response='y');warnings<-character()
 fits<-tryCatch(withCallingHandlers({
  full<-MASS::polr(form,data=frame,method='logistic',Hess=TRUE,model=TRUE,control=list(maxit=1000,reltol=1e-12))
  null<-MASS::polr(y~1,data=frame,method='logistic',Hess=TRUE,model=TRUE,control=list(maxit=1000,reltol=1e-12))
  list(full=full,null=null)
 },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
 if(inherits(fits,'error'))return(unavailable(conditionMessage(fits),warnings))
 f<-fits$full;null<-fits$null
 if(length(warnings)||f$convergence!=0L||null$convergence!=0L||any(!is.finite(f$Hessian))||rcond(f$Hessian)<1e-10)
  return(unavailable('Convergence or Hessian diagnostics prevent ordinary maximum-likelihood inference.',warnings))
 if(min(eigen(f$Hessian,symmetric=TRUE,only.values=TRUE)$values)<=0)
  return(unavailable('The likelihood Hessian is not positive definite.'))
 if(any(f$fitted.values<1e-8|f$fitted.values>1-1e-8))
  return(unavailable('Near-boundary category probabilities require separation/sparsity review.'))
 critical<-stats::qnorm((1-conf.level)/2,lower.tail=FALSE)
 V<-stats::vcov(f);A<-diag(c(1/scale,rep(1,K-1L)))
 A[seq.int(q+1L,q+K-1L),seq_len(q)]<-matrix(rep(center/scale,K-1L),nrow=K-1L,byrow=TRUE)
 coef<-as.vector(A%*%c(stats::coef(f),f$zeta));covariance<-A%*%V%*%t(A);se<-sqrt(diag(covariance))
 beta<-coef[seq_len(q)];ses<-se[seq_len(q)];lower<-beta-critical*ses;upper<-beta+critical*ses
 tab<-data.frame(term=predictors,estimate=beta,std_error=ses,z=beta/ses,p_value=2*stats::pnorm(abs(beta/ses),lower.tail=FALSE),lower=lower,upper=upper,odds_ratio=exp(beta),ratio_lower=exp(lower),ratio_upper=exp(upper),stringsAsFactors=FALSE)
 thresholds<-data.frame(boundary=names(f$zeta),estimate=coef[-seq_len(q)],std_error=se[-seq_len(q)])
 LR<-2*(as.numeric(stats::logLik(f))-as.numeric(stats::logLik(null)))
 if(LR<0&&abs(LR)<=64*.Machine$double.eps*max(1,abs(as.numeric(stats::logLik(f)))))LR<-0
 if(!is.finite(LR)||LR<0||any(!is.finite(as.matrix(tab[vapply(tab,is.numeric,logical(1))]))))
  return(unavailable('Likelihood or original-unit coefficient inference is numerically unresolved.'))
 df<-q;p<-stats::pchisq(LR,df,lower.tail=FALSE)
 list(status='ok',method='Proportional-odds cumulative logistic regression',coefficient_table=tab,thresholds=thresholds,sign_convention='logit P(Y <= category k) = threshold[k] - x beta; exp(beta) is the common odds ratio for higher categories.',
  p_value=p,omnibus=list(statistic=LR,df=df,method='Asymptotic likelihood-ratio chi-square',null='All prespecified slopes equal zero'),
  n_total=nrow(data),n_used=n,n_excluded=sum(!keep),category_order=levels(d$y),fitted=f$fitted.values,
  confidence_level=conf.level,interval_method='Marginal asymptotic Wald intervals; no profile or multiplicity adjustment.',
  assumptions='Independent ordered outcomes, common slopes across cumulative logits, correct additive logit specification and finite maximum likelihood. The proportional-odds assumption is explicitly justified before using this template; no automatic pretest-driven substitution.',limitations='One to ten numeric predictors and three to eight observed categories. No automatic variable selection, category pooling, penalization or model substitution. Numerical diagnostics do not replace separation and model-assumption review. Complete-case inference needs justified missingness assumptions; no unsupported causal interpretation.')
}
example_data<-data.frame(id=1:60,x1=rep(c(-1,0,1),each=20),y=factor(c(rep(1,10),rep(2,6),rep(3,4),rep(1,6),rep(2,8),rep(3,6),rep(1,3),rep(2,6),rep(3,11)),levels=1:3,ordered=TRUE))
example_result<-analyse(example_data,'x1',TRUE,TRUE,alternative=alternative,conf.level=confidence_level)
print(example_result)
