# Explicit mixed measurement types, weights and k; no silent type downgrading.
gower_data <- function(data,variables,types) {
 x<-data[variables]
 for(j in seq_along(variables)){
  value<-x[[j]];type<-types[j];stopifnot(!all(is.na(value)),is.null(dim(value)))
  if(type=='continuous')stopifnot(is.numeric(value),!any(is.infinite(value)))
  if(type=='nominal')stopifnot(is.factor(value),!is.ordered(value))
  if(type=='ordinal')stopifnot(is.ordered(value))
  if(type%in%c('symmetric_binary','asymmetric_binary')){
   stopifnot(is.numeric(value)||is.logical(value),all(is.na(value)|value%in%c(0,1)));x[[j]]<-as.numeric(value)
  }
 }
 x
}
gower_fit <- function(x,k,types,weights) {
 distance<-cluster::daisy(x,metric='gower',type=list(symm=which(types=='symmetric_binary'),asymm=which(types=='asymmetric_binary')),
  weights=weights,warnBin=FALSE,warnConst=FALSE)
 if(any(!is.finite(distance))||max(distance)<=0)stop('Gower dissimilarities have no variation or a pair has no comparable information.')
 fit<-cluster::pam(distance,k=k,diss=TRUE,keep.diss=TRUE,keep.data=FALSE)
 groups<-match(fit$clustering,unique(fit$clustering));medoids<-fit$id.med[match(seq_len(k),groups[fit$id.med])]
 list(groups=groups,medoids=medoids,distance=as.matrix(distance),mean_dissimilarity=unname(fit$objective[2]),
  silhouette=fit$silinfo$avg.width)
}
adjusted_rand <- function(a,b) {
 choose2<-function(n)n*(n-1)/2;tab<-table(a,b);total<-choose2(sum(tab));same<-sum(choose2(tab))
 left<-sum(choose2(rowSums(tab)));right<-sum(choose2(colSums(tab)));expected<-left*right/total;denominator<-(left+right)/2-expected
 if(denominator==0)1 else (same-expected)/denominator
}
analyse <- function(data,variables,types,weights,k,independent_units=FALSE,alternative='two.sided',conf.level=.95,
 missing='error',B=50L,seed=104729L) {
 if(!requireNamespace('cluster',quietly=TRUE))stop('This example requires the recommended cluster package.')
 stopifnot(isTRUE(independent_units),identical(alternative,'two.sided'),is.data.frame(data),!anyDuplicated(names(data)),
  is.character(variables),length(variables)>=2L,length(variables)<=30L,!anyNA(variables),!anyDuplicated(variables),
  !('id'%in%variables),all(c('id',variables)%in%names(data)),is.atomic(data$id),is.null(dim(data$id)),!anyNA(data$id),
  !anyDuplicated(data$id),all(nzchar(as.character(data$id))),is.character(types),!is.null(names(types)),!anyDuplicated(names(types)),
  setequal(names(types),variables),all(types%in%c('continuous','nominal','ordinal','symmetric_binary','asymmetric_binary')),
  is.numeric(weights),!is.null(names(weights)),!anyDuplicated(names(weights)),setequal(names(weights),variables),all(is.finite(weights)),all(weights>0),
  is.numeric(k),length(k)==1L,is.finite(k),k==floor(k),k>=2L,k<=20L,
  is.numeric(B),length(B)==1L,is.finite(B),B==floor(B),B>=20L,B<=200L,
  is.numeric(seed),length(seed)==1L,is.finite(seed),seed==floor(seed),seed>=0,seed<=.Machine$integer.max,
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1,
  length(missing)==1L,missing%in%c('error','pairwise'),nrow(data)<=500L)
 d<-data[order(enc2utf8(as.character(data$id)),method='radix'),,drop=FALSE];types<-types[variables];weights<-weights[variables]
 x<-gower_data(d,variables,types);if(anyNA(x)&&missing=='error')stop('Missing measurements require an explicit pairwise-comparison choice; no imputation is supplied.')
 n<-nrow(x);subn<-floor(.8*n);stopifnot(n>=5L,k<subn)
 kind<-RNGkind();had_seed<-exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)
 if(had_seed)saved_seed<-get('.Random.seed',envir=.GlobalEnv,inherits=FALSE)
 on.exit({do.call(RNGkind,as.list(kind));if(had_seed)assign('.Random.seed',saved_seed,envir=.GlobalEnv) else if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE))rm('.Random.seed',envir=.GlobalEnv)},add=TRUE)
 RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(as.integer(seed));full<-gower_fit(x,k,types,weights)
 stability<-rep(NA_real_,B);failures<-character()
 for(b in seq_len(B)){
  index<-sort(sample.int(n,subn,replace=FALSE));candidate<-tryCatch(gower_fit(x[index,,drop=FALSE],k,types,weights),error=function(e)e)
  if(inherits(candidate,'error'))failures<-c(failures,paste0('Subsample ',b,': ',conditionMessage(candidate))) else stability[b]<-adjusted_rand(full$groups[index],candidate$groups)
 }
 list(status=if(length(failures))'stability_incomplete' else 'ok',method='PAM on explicitly typed weighted Gower dissimilarities',k=k,
  membership=data.frame(id=d$id,cluster=full$groups),medoid_ids=d$id[full$medoids],dissimilarity=full$distance,
  mean_dissimilarity=full$mean_dissimilarity,mean_silhouette=full$silhouette,variable_types=types,weights=weights,missing=missing,
  stability=list(method='80% subsampling, recomputing ranges/Gower distances and PAM; adjusted Rand index versus the restricted full partition',scores=stability,
   mean=if(length(failures))NA_real_ else mean(stability),failed=length(failures),failures=failures,B=B,seed=seed),
  confidence_interval=NULL,interval_status='not_computed',requested_confidence_level=conf.level,
  interpretation='Continuous and ordinal differences use observed-range scaling; nominal mismatches and declared binary conventions enter with explicit weights. Pairwise missing values and joint asymmetric-binary zeros contribute neither distance nor weight for that variable.',
  limitations='Exploratory grouping; k and weights are prespecified. PAM is an optimization heuristic, not a general proof of global optimality. Missing-pairwise comparisons may use different variable subsets. Stability and silhouette do not establish real populations or external validity; failed subsamples are not selectively averaged.')
}
example_data<-data.frame(id=sprintf('case%02d',1:12),amount=c(1,2,1.5,2.5,1.8,2.2,8,9,8.5,9.5,8.8,9.2),kind=factor(rep(c('a','b'),each=6)),grade=ordered(rep(c('low','high'),each=6),levels=c('low','middle','high')),flag=rep(c(0,1),each=6))
example_types<-c(amount='continuous',kind='nominal',grade='ordinal',flag='asymmetric_binary')
example_result<-analyse(example_data,names(example_types),example_types,stats::setNames(rep(1,4),names(example_types)),2,TRUE)
print(example_result[c('status','method','medoid_ids','mean_silhouette','stability')])
