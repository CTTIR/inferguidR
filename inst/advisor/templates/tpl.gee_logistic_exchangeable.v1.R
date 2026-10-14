# Prespecified two-group clustered binary example, with explicit estimand.
# A unique id identifies each observation; cluster identifies dependent units.
analyse <- function(data,independent_clusters=FALSE,model_assumptions_justified=FALSE,
 alternative='two.sided',conf.level=.95,missing_outcomes='error') {
 if(!requireNamespace('geepack',quietly=TRUE))stop('This example requires geepack.')
 stopifnot(isTRUE(independent_clusters),isTRUE(model_assumptions_justified),identical(alternative,'two.sided'),
  is.data.frame(data),!anyDuplicated(names(data)),all(c('id','cluster','group','y')%in%names(data)),
  is.atomic(data$id),is.null(dim(data$id)),!anyNA(data$id),!anyDuplicated(data$id),all(nzchar(as.character(data$id))),
  is.atomic(data$cluster),is.null(dim(data$cluster)),!anyNA(data$cluster),all(nzchar(as.character(data$cluster))),
  is.factor(data$group),!is.ordered(data$group),nlevels(data$group)==2L,!anyNA(data$group),
  is.numeric(data$y),is.null(dim(data$y)),all(is.na(data$y)|data$y%in%c(0,1)),nrow(data)<=10000L,
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1,
  length(missing_outcomes)==1L,missing_outcomes%in%c('error','available_case_MCAR'),
  all(vapply(split(as.character(data$group),data$cluster),function(x)length(unique(x))==1L,logical(1))))
 keep<-!is.na(data$y)
 if(!all(keep)&&missing_outcomes=='error')stop('Missing outcomes require the explicitly justified missingness strategy.')
 d<-data[keep,,drop=FALSE];d$cluster<-factor(d$cluster);d$group_numeric<-as.integer(d$group)-1L
 d<-d[order(d$cluster,d$id),,drop=FALSE];G<-nlevels(d$cluster);sizes<-table(d$cluster)
 stopifnot(G>=30L,all(sizes>=2L),all(table(d$group)>0),length(unique(d$y))==2L,
  all(table(tapply(d$group_numeric,d$cluster,function(x)x[1]))>=5L))
 unavailable<-function(reason,warnings=character())list(status='not_applicable',reason=reason,warnings=warnings,p_value=NULL,
  confidence_interval=NULL,n_total=nrow(data),n_used=nrow(d),n_excluded=sum(!keep),clusters=G)
 alerts<-character()
 # geeglm evaluates an initial glm call in its caller; bind it locally.
 glm<-stats::glm;model.frame<-stats::model.frame
 fit<-tryCatch(withCallingHandlers(geepack::geeglm(y~group_numeric,id=cluster,data=d,family=stats::binomial(),
  corstr='exchangeable',std.err='san.se',scale.fix=TRUE,na.action=stats::na.fail,
  control=geepack::geese.control(epsilon=1e-10,maxit=100)),
  warning=function(w){alerts<<-c(alerts,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
 if(inherits(fit,'error'))return(unavailable(conditionMessage(fit),alerts))
 f<-fit;alpha<-unname(f$geese$alpha)
 if(length(alerts)||f$geese$error!=0L||length(alpha)!=1L||!is.finite(alpha)||
  alpha<=-1/(max(sizes)-1)+1e-8||alpha>=1-1e-8||any(f$fitted.values<1e-8|f$fitted.values>1-1e-8))
  return(unavailable('GEE convergence, correlation or mean-boundary diagnostics prevent ordinary inference.',alerts))

 beta<-as.numeric(stats::coef(f));V<-f$geese$vbeta
 if(any(!is.finite(beta))||any(!is.finite(V))||any(diag(V)<=0)||rcond(V)<1e-12)
  return(unavailable('Coefficient covariance is numerically unresolved.'))
 se<-unname(sqrt(diag(V)));z<-beta[2]/se[2];critical<-stats::qnorm((1-conf.level)/2,lower.tail=FALSE)
 limits<-beta[2]+c(-1,1)*critical*se[2]
 if(any(!is.finite(exp(c(beta[2],limits)))))return(unavailable('The odds-ratio estimate or interval exceeds the numerical range.'))
 list(status='ok',method='Logistic GEE with exchangeable working correlation',estimand='Population-averaged marginal event odds ratio; not a conditional mixed-model odds ratio.',
  group_reference=levels(data$group)[1],group_comparison=levels(data$group)[2],event='y=1',
  log_odds_ratio=beta[2],std_error=se[2],statistic=z,p_value=2*stats::pnorm(abs(z),lower.tail=FALSE),
  odds_ratio=exp(beta[2]),confidence_interval=exp(limits),confidence_level=conf.level,
  interval_method='Asymptotic robust sandwich Wald normal interval; no small-cluster correction.',working_correlation=alpha,sandwich='geepack robust empirical sandwich, uncorrected for small cluster counts',
  n_total=nrow(data),n_used=nrow(d),n_excluded=sum(!keep),clusters=G,
  coefficients=stats::setNames(beta,c('(Intercept)','comparison_vs_reference')),covariance=V,
  missing_outcomes=missing_outcomes,assumptions='Independent clusters, correctly specified marginal mean, exchangeable working correlation and enough clusters for sandwich asymptotics. Unweighted available-case GEE requires MCAR, not merely MAR; no weighting or imputation is supplied.',
  limitations='The cluster minimum is an implementation floor, not a guarantee of accurate small-sample inference. Only a prespecified binary cluster-level group contrast is demonstrated. No automatic model, estimand or covariance switch; no causal claim without supporting design. Boundary or convergence issues require specialist review.')
}
counts<-c(0,1,1,2,2,3,4,4,5,6,0,2,3,5,6,1,2,3,3,4,4,5,5,6,6,2,3,4,5,6)
example_data<-data.frame(id=1:180,cluster=rep(1:30,each=6),group=factor(rep(c('control','treatment'),each=90),levels=c('control','treatment')),y=unlist(lapply(counts,function(s)c(rep(1,s),rep(0,6-s)))))
example_result<-analyse(example_data,TRUE,TRUE,alternative,confidence_level)
print(example_result)
