# Exploratory common-factor parallel analysis with squared multiple correlations.
# Independent column permutations preserve marginal values and break association.
# The result supports a factor-count decision; it never fits or selects an ML model.
parallel_analysis <- function(data,variables,independent_units=FALSE,continuous_measurements=FALSE,
 B=199L,seed=104729L,quantile_probability=.95,complete_cases=FALSE,permutations=NULL) {
 stopifnot(isTRUE(independent_units),isTRUE(continuous_measurements),is.data.frame(data),!anyDuplicated(names(data)),
  is.character(variables),length(variables)>=3L,length(variables)<=30L,!anyNA(variables),!anyDuplicated(variables),
  !('id'%in%variables),all(c('id',variables)%in%names(data)),is.atomic(data$id),is.null(dim(data$id)),!anyNA(data$id),!anyDuplicated(data$id),
  all(nzchar(as.character(data$id))),all(vapply(data[variables],function(x)is.numeric(x)&&is.null(dim(x)),logical(1))),
  !any(is.infinite(as.matrix(data[variables]))),is.logical(complete_cases),length(complete_cases)==1L,!is.na(complete_cases),
  is.numeric(B),length(B)==1L,is.finite(B),B==floor(B),B>=100L,B<=2000L,
  is.numeric(seed),length(seed)==1L,is.finite(seed),seed==floor(seed),seed>=0,seed<=.Machine$integer.max,
  is.numeric(quantile_probability),length(quantile_probability)==1L,is.finite(quantile_probability),quantile_probability>=.5,quantile_probability<1,
  B*(1-quantile_probability)>=5-1e-12)
 keep<-stats::complete.cases(data[variables]);if(!all(keep)&&!complete_cases)stop('Incomplete cases require explicit complete-case selection.')
 x<-as.matrix(data[keep,variables,drop=FALSE]);n<-nrow(x);p<-ncol(x);stopifnot(n>p,n<=10000L,n*p*B<=5000000)
 unavailable<-function(reason)list(status='unavailable',reason=reason,suggested_leading_factors=NULL)
 centered<-sweep(x,2L,colMeans(x),'-');scale<-apply(abs(centered),2L,max)
 if(any(!is.finite(centered))||any(!is.finite(scale))||any(scale<=0))return(unavailable('Variable variation is absent or unresolved.'))
 x<-sweep(centered,2L,scale,'/')
 spectrum<-function(z){C<-stats::cor(z);if(any(!is.finite(C))||rcond(C)<1e-10)stop('Correlation matrix is singular or unresolved.')
  smc<-1-1/diag(solve(C));if(any(!is.finite(smc))||any(smc< -1e-10)||any(smc>1+1e-10))stop('Squared multiple correlations are unresolved.')
  diag(C)<-pmin(1,pmax(0,smc));eigen(C,symmetric=TRUE,only.values=TRUE)$values}
 observed<-tryCatch(spectrum(x),error=function(e)e);if(inherits(observed,'error'))return(unavailable(conditionMessage(observed)))
 supplied<-!is.null(permutations)
 if(supplied){stopifnot(is.array(permutations),identical(dim(permutations),as.integer(c(n,p,B))),is.numeric(permutations),
  all(is.finite(permutations)),all(permutations==floor(permutations)),all(permutations>=1L),all(permutations<=n))
  for(b in seq_len(B))for(j in seq_len(p))stopifnot(identical(sort(as.integer(permutations[,j,b])),seq_len(n)))
 }else{
  had_seed<-exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE);if(had_seed)saved_seed<-get('.Random.seed',envir=.GlobalEnv)
  saved_kind<-RNGkind();on.exit({do.call(RNGkind,as.list(saved_kind));if(had_seed)assign('.Random.seed',saved_seed,envir=.GlobalEnv) else if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE))rm('.Random.seed',envir=.GlobalEnv)},add=TRUE)
  RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(seed)
 }
 null<-matrix(NA_real_,B,p);failures<-list()
 for(b in seq_len(B)){
  permuted<-x
  for(j in seq_len(p)){idx<-if(supplied)permutations[,j,b] else sample.int(n);permuted[,j]<-x[idx,j]}
  v<-tryCatch(spectrum(permuted),error=function(e)e)
  if(inherits(v,'error'))failures[[length(failures)+1L]]<-list(replicate=b,reason=conditionMessage(v)) else null[b,]<-v
 }
 if(length(failures))return(list(status='simulation_incomplete',failures=failures,failed_replicates=length(failures),B=B,suggested_leading_factors=NULL,reason='No failed replicate is discarded; no factor count is reported.'))
 threshold<-apply(null,2L,stats::quantile,probs=quantile_probability,type=6,names=FALSE)
 above<-observed>threshold & observed>0;leading<-sum(cumprod(as.integer(above)))
 list(status='ok',method='Common-factor parallel analysis; SMC diagonal, independent column permutations, type-6 quantiles',
  observed_eigenvalues=observed,null_eigenvalues=null,null_thresholds=threshold,above_threshold=above,suggested_leading_factors=leading,
  n_total=nrow(data),n_used=n,n_excluded=sum(!keep),B=B,seed=if(supplied)NULL else seed,
  resampling_source=if(supplied)'supplied column permutations' else 'generated independent column permutations',quantile_probability=quantile_probability,
  limitations='Exploratory heuristic with Monte Carlo variation; not a hypothesis test, interval, guarantee of a true factor count, or fit/identifiability check. SMCs may reflect minor factors. Prespecify the null scheme and quantile; consider theory, scree and replication. Repeated or clustered observations violate the permutation design. No automatic transfer to the ML factor fit.')
}

# Prespecified number of common factors; maximum likelihood and oblique promax rotation.
# Numeric continuous measurements; no PCA substitution or Kaiser-only selection.
analyse <- function(data,variables,factors,independent_units=FALSE,continuous_measurements=FALSE,
 alternative='two.sided',conf.level=.95,complete_cases=FALSE) {
 stopifnot(isTRUE(independent_units),isTRUE(continuous_measurements),identical(alternative,'two.sided'),
  is.data.frame(data),!anyDuplicated(names(data)),is.character(variables),length(variables)>=3L,length(variables)<=30L,
  !anyNA(variables),!anyDuplicated(variables),!('id'%in%variables),all(c('id',variables)%in%names(data)),
  is.atomic(data$id),is.null(dim(data$id)),!anyNA(data$id),!anyDuplicated(data$id),all(nzchar(as.character(data$id))),
  is.numeric(factors),length(factors)==1L,is.finite(factors),factors==floor(factors),factors>=1L,
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1,
  is.logical(complete_cases),length(complete_cases)==1L,!is.na(complete_cases),nrow(data)<=10000L,
  all(vapply(data[variables],function(x)is.numeric(x)&&is.null(dim(x)),logical(1))),!any(is.infinite(as.matrix(data[variables]))))
 p<-length(variables);stopifnot(factors<p,((p-factors)^2-p-factors)/2>0)
 keep<-stats::complete.cases(data[variables])
 if(!all(keep)&&!complete_cases)stop('Incomplete cases require explicit complete-case selection; no imputation is performed.')
 x<-as.matrix(data[keep,variables,drop=FALSE]);stopifnot(nrow(x)>p)
 unavailable<-function(reason)list(status='not_applicable',reason=reason,p_value=NULL,confidence_interval=NULL,
  n_total=nrow(data),n_used=nrow(x),n_excluded=sum(!keep))
 # Normalize each column before covariance calculation to avoid squared overflow.
 center<-colMeans(x);z<-sweep(x,2L,center,'-');scale<-apply(abs(z),2L,max)
 if(any(!is.finite(z))||any(!is.finite(scale))||any(scale<=0))return(unavailable('Variable variation is absent or numerically unresolved.'))
 z<-sweep(z,2L,scale,'/');C<-stats::cor(z)
 if(any(!is.finite(C))||rcond(C)<1e-10)return(unavailable('The correlation matrix is singular or numerically unresolved.'))
 alerts<-character()
 fit<-tryCatch(withCallingHandlers(stats::factanal(covmat=C,n.obs=nrow(x),factors=factors,rotation='promax',scores='none',
  start=matrix(rep(c(.2,.5,.8),each=p),nrow=p),control=list(lower=.005,opt=list(maxit=2000,factr=1e7))),
  warning=function(w){alerts<<-c(alerts,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
 if(inherits(fit,'error'))return(unavailable(conditionMessage(fit)))
 if(length(alerts)||!isTRUE(fit$converged)||any(!is.finite(fit$uniquenesses))||any(fit$uniquenesses<=.005001))
  return(unavailable('Optimization or near-Heywood uniqueness boundary prevents accepting this factor solution.'))
 L<-unclass(fit$loadings);Phi<-attr(fit$loadings,'covariance')
 if(is.null(Phi))Phi<-if(is.null(fit$rotmat))diag(factors) else solve(crossprod(fit$rotmat))
 reproduced<-L%*%Phi%*%t(L)+diag(fit$uniquenesses)
 if(any(!is.finite(reproduced))||any(!is.finite(Phi)))return(unavailable('The rotated factor covariance is unresolved.'))
 list(status='ok',method='Maximum-likelihood common-factor analysis with oblique promax rotation',factors=factors,
  loadings=L,factor_correlation=Phi,uniquenesses=fit$uniquenesses,reproduced_correlation=reproduced,
  observed_correlation=C,residual_correlation=C-reproduced,objective=unname(fit$criteria['objective']),
  goodness_of_fit=list(statistic=unname(fit$STATISTIC),df=fit$dof,p_value=unname(fit$PVAL),
   interpretation='Approximate normal-theory goodness-of-fit test for the prespecified factor count; not a rule for selecting the factor count.'),
  n_total=nrow(data),n_used=nrow(x),n_excluded=sum(!keep),confidence_interval=NULL,interval_status='not_computed',requested_confidence_level=conf.level,
  assumptions='Independent observations of continuous measurements; common-factor model with independent uniquenesses and multivariate-normal likelihood assumptions. Oblique factors may correlate. Complete-case validity requires justified missingness assumptions.',
  limitations='Factor count must be explicitly chosen; the separate parallel_analysis helper provides exploratory guidance without selecting the fitted model. Three deterministic uniqueness starts reduce but do not eliminate local-optimum risk. Signs/order/rotation are not uniquely identified. No sampling intervals or claim of factor stability is supplied; no PCA or Kaiser-only substitution.')
}
t<-seq_len(120);f1<-sin(t*.13);f2<-.3*f1+cos(t*.21)
example_data<-data.frame(id=t)
for(j in 1:6)example_data[[paste0('v',j)]]<-(if(j<=3).8*f1+.2*f2 else .2*f1+.8*f2)+.7*sin(t*(.61+j*.1))
example_result<-analyse(example_data,paste0('v',1:6),2,TRUE,TRUE)
print(example_result)

example_parallel<-parallel_analysis(example_data,paste0('v',1:6),TRUE,TRUE,B=100L)
print(example_parallel[c('status','observed_eigenvalues','null_thresholds','suggested_leading_factors','limitations')])
