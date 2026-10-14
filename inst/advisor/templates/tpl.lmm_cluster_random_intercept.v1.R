# Prespecified Gaussian random-intercept model; covariance choice is explicit.
# REML estimates, marginal fixed-effect F test, nlme containment degrees of freedom.
# No comparison of REML likelihoods across different fixed-effect models.
analyse <- function(data,independent_clusters=FALSE,random_intercept_justified=FALSE,
 alternative='two.sided',conf.level=.95,missing_outcomes='error') {
 if(!requireNamespace('nlme',quietly=TRUE))stop('This example requires the recommended nlme package.')
 stopifnot(isTRUE(independent_clusters),isTRUE(random_intercept_justified),identical(alternative,'two.sided'),
  is.data.frame(data),!anyDuplicated(names(data)),all(c('id','cluster','y','factor')%in%names(data)),
  is.atomic(data$id),is.null(dim(data$id)),!anyNA(data$id),!anyDuplicated(data$id),all(nzchar(as.character(data$id))),
  is.atomic(data$cluster),is.null(dim(data$cluster)),!anyNA(data$cluster),all(nzchar(as.character(data$cluster))),
  is.factor(data$factor),!is.ordered(data$factor),!anyNA(data$factor),nlevels(data$factor)>=2L,nlevels(data$factor)<=30L,
  is.numeric(data$y),is.null(dim(data$y)),!any(is.infinite(data$y)),nrow(data)<=10000L,
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1,
  length(missing_outcomes)==1L,missing_outcomes%in%c('error','available_case_MAR'))
 stopifnot(all(vapply(split(as.character(data$factor),data$cluster),function(x)length(unique(x))==1L,logical(1))))
 keep<-!is.na(data$y)
 if(!all(keep)&&missing_outcomes=='error')stop('Missing outcomes require a prespecified available-case MAR assumption.')
 d<-data[keep,,drop=FALSE];d$cluster<-factor(d$cluster);G<-nlevels(d$cluster);K<-nlevels(d$factor)
 stopifnot(G>=6L,all(table(d$factor)>0),sum(table(d$cluster)>=2L)>=3L,nrow(d)>G+K)
 unavailable<-function(reason,warnings=character())list(status='not_applicable',reason=reason,warnings=warnings,p_value=NULL,
  confidence_interval=NULL,n_total=nrow(data),n_used=nrow(d),n_excluded=sum(!keep),clusters=G)
 center<-mean(d$y);scale<-max(abs(d$y-center))
 if(!is.finite(scale)||scale<=0)return(unavailable('Outcome variation is zero or numerically unresolved.'))
 d$y<-(d$y-center)/scale;d<-d[order(d$cluster),,drop=FALSE]
 alerts<-character();fit<-tryCatch(withCallingHandlers(nlme::lme(y~factor,random=~1|cluster,data=d,
  method='REML',na.action=stats::na.fail,contrasts=list(factor=stats::contr.sum(K)),
  control=nlme::lmeControl(maxIter=100,msMaxIter=200,tolerance=1e-8,msTol=1e-9,returnObject=FALSE)),
  warning=function(w){alerts<<-c(alerts,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
 if(inherits(fit,'error'))return(unavailable(conditionMessage(fit),alerts))
 if(length(alerts)||!is.matrix(fit$apVar)||any(!is.finite(fit$apVar)))
  return(unavailable('Variance-component convergence or curvature is unresolved.',alerts))
 variances<-as.numeric(nlme::VarCorr(fit)[,'Variance'])
 if(any(!is.finite(variances))||any(variances<=0)||variances[1]/sum(variances)<1e-8)
  return(unavailable('Random-intercept variance is at or near its boundary; no automatic fixed-model fallback.'))
 table<-stats::anova(fit,type='marginal');test<-table['factor',]
 if(any(!is.finite(unlist(test)))||test$denDF<=0)return(unavailable('The fixed-effect F test is unresolved.'))
 # Estimated factor means, with explicitly marginal t intervals.
 grid<-data.frame(factor=factor(levels(d$factor),levels=levels(d$factor)))
 L<-stats::model.matrix(~factor,grid,contrasts.arg=list(factor=stats::contr.sum(K)))
 b<-nlme::fixef(fit);V<-stats::vcov(fit);means<-as.vector(L%*%b)*scale+center
 se<-sqrt(diag(L%*%V%*%t(L)))*scale
 df_mean<-min(fit$fixDF$X);critical<-stats::qt((1-conf.level)/2,df_mean,lower.tail=FALSE)
 if(any(!is.finite(c(means,se,variances*scale^2))))return(unavailable('Original-unit estimates exceed the numerical range.'))
 list(status='ok',method='Cluster-level factor Gaussian random-intercept model',covariance='Random intercept plus independent equal-variance errors: compound symmetry; not unstructured or AR(1).',
  statistic=unname(test[['F-value']]),df_numerator=unname(test$numDF),df_denominator=unname(test$denDF),p_value=unname(test[['p-value']]),
  factor_means=data.frame(level=levels(d$factor),estimate=means,std_error=se,lower=means-critical*se,upper=means+critical*se),
  mean_interval_df=df_mean,mean_interval_method='Marginal t intervals using the minimum containment df among fixed coefficients; not simultaneous or pairwise comparisons.',
  random_intercept_variance=variances[1]*scale^2,residual_variance=variances[2]*scale^2,icc=variances[1]/sum(variances),
  n_total=nrow(data),n_used=nrow(d),n_excluded=sum(!keep),clusters=G,clusters_excluded=length(unique(data$cluster))-G,
  confidence_level=conf.level,missing_outcomes=missing_outcomes,
  assumptions='Independent clusters, Gaussian random intercepts independent of Gaussian homoscedastic residuals, fixed factor effects and a justified compound-symmetric covariance. Available-case likelihood inference requires ignorable missingness (MAR given the model); missing covariates are not imputed.',
  limitations='At least six observed clusters is an implementation floor, not a guarantee of reliable small-sample inference. Denominator df use nlme containment, not Satterthwaite or Kenward-Roger. Covariance structure is prespecified; no pretest-driven switch. No variance-component test, sampling interval for ICC, or automatic post-hoc contrast is computed.')
}
example_data<-expand.grid(replicate=1:3,cluster=1:12);example_data$id<-seq_len(nrow(example_data))
example_data$factor<-factor(ifelse(example_data$cluster<=6,'control','treatment'),levels=c('control','treatment'))
example_data$y<-rep(c(-3,-2,-1,1,2,3,0,1,2,4,5,6),each=3)+rep(c(-1,0,1),12)
example_result<-analyse(example_data,TRUE,TRUE,alternative,confidence_level)
print(example_result)
