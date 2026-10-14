# Exploratory clustering with prespecified scaling, k and subsampling stability.
# Case identifiers establish deterministic tie ordering; labels have no intrinsic meaning.
cluster_scale <- function(x,standardize) {
 center<-colMeans(x);z<-sweep(x,2L,center,'-')
 if(any(!is.finite(z)))stop('Centered measurements exceed the numerical range.')
 if(standardize){
  a<-apply(abs(z),2L,max);if(any(a<=0))stop('A variable is constant in this sample.')
  z<-sweep(z,2L,a,'/');b<-apply(z,2L,stats::sd);if(any(b<=0)||any(!is.finite(b)))stop('Variable scaling is unresolved.')
  z<-sweep(z,2L,b,'/')
 }
 z
}
cluster_fit <- function(x,k,standardize) {
 z<-cluster_scale(x,standardize);if(nrow(unique(as.data.frame(z)))<k)stop('Fewer distinct profiles than prespecified clusters.')
 distance<-stats::dist(z);if(any(!is.finite(distance)))stop('Euclidean distances are unresolved.')
 tree<-stats::hclust(distance,method='ward.D2');groups<-stats::cutree(tree,k=k)
 groups<-match(groups,unique(groups));list(groups=groups,scaled=z,tree=tree)
}
adjusted_rand <- function(a,b) {
 choose2<-function(n)n*(n-1)/2;tab<-table(a,b);n<-sum(tab);total<-choose2(n)
 same<-sum(choose2(tab));left<-sum(choose2(rowSums(tab)));right<-sum(choose2(colSums(tab)))
 expected<-left*right/total;denominator<-(left+right)/2-expected
 if(denominator==0)return(1)
 (same-expected)/denominator
}
silhouette_values <- function(distance,groups) {
 vapply(seq_along(groups),function(i){
  same<-which(groups==groups[i]);if(length(same)==1L)return(0)
  a<-mean(distance[i,setdiff(same,i)]);others<-setdiff(unique(groups),groups[i])
  b<-min(vapply(others,function(g)mean(distance[i,groups==g]),numeric(1)))
  if(max(a,b)==0)0 else (b-a)/max(a,b)
 },numeric(1))
}
analyse <- function(data,variables,k,standardize,independent_units=FALSE,
 alternative='two.sided',conf.level=.95,complete_cases=FALSE,B=50L,seed=104729L) {
 stopifnot(isTRUE(independent_units),identical(alternative,'two.sided'),is.data.frame(data),!anyDuplicated(names(data)),
  is.character(variables),length(variables)>=1L,length(variables)<=30L,!anyNA(variables),!anyDuplicated(variables),
  !('id'%in%variables),all(c('id',variables)%in%names(data)),is.atomic(data$id),is.null(dim(data$id)),!anyNA(data$id),
  !anyDuplicated(data$id),all(nzchar(as.character(data$id))),is.logical(standardize),length(standardize)==1L,!is.na(standardize),
  is.numeric(k),length(k)==1L,is.finite(k),k==floor(k),k>=2L,k<=20L,
  is.numeric(B),length(B)==1L,is.finite(B),B==floor(B),B>=20L,B<=500L,
  is.numeric(seed),length(seed)==1L,is.finite(seed),seed==floor(seed),seed>=0,seed<=.Machine$integer.max,
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1,
  is.logical(complete_cases),length(complete_cases)==1L,!is.na(complete_cases),nrow(data)<=1000L,
  all(vapply(data[variables],function(x)is.numeric(x)&&is.null(dim(x)),logical(1))),!any(is.infinite(as.matrix(data[variables]))))
 keep<-stats::complete.cases(data[variables]);if(!all(keep)&&!complete_cases)stop('Missing measurements require explicit complete-case selection; no mean imputation.')
 d<-data[keep,,drop=FALSE];d<-d[order(enc2utf8(as.character(d$id)),method='radix'),,drop=FALSE];x<-as.matrix(d[variables]);n<-nrow(x);subn<-floor(.8*n)
 stopifnot(n>=5L,k<subn)
 kind<-RNGkind();had_seed<-exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)
 if(had_seed)saved_seed<-get('.Random.seed',envir=.GlobalEnv,inherits=FALSE)
 on.exit({do.call(RNGkind,as.list(kind));if(had_seed)assign('.Random.seed',saved_seed,envir=.GlobalEnv) else if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE))rm('.Random.seed',envir=.GlobalEnv)},add=TRUE)
 RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(as.integer(seed))
 full<-cluster_fit(x,k,standardize);distance<-as.matrix(stats::dist(full$scaled));if(any(!is.finite(distance)))stop('Euclidean distances are unresolved.')
 widths<-silhouette_values(distance,full$groups);stability<-rep(NA_real_,B);failures<-character()
 for(b in seq_len(B)){
  index<-sort(sample.int(n,subn,replace=FALSE));candidate<-tryCatch(cluster_fit(x[index,,drop=FALSE],k,standardize),error=function(e)e)
  if(inherits(candidate,'error'))failures<-c(failures,paste0('Subsample ',b,': ',conditionMessage(candidate))) else stability[b]<-adjusted_rand(full$groups[index],candidate$groups)
 }
 list(status=if(length(failures))'stability_incomplete' else 'ok',method='Hierarchical Ward.D2 clustering on Euclidean distances',k=k,standardize=standardize,
  membership=data.frame(id=d$id,cluster=full$groups,silhouette=widths),mean_silhouette=mean(widths),fit=full,
  stability=list(method='80% subsampling without replacement; preprocessing refitted; adjusted Rand index versus full partition restricted to retained cases',
   scores=stability,mean=if(length(failures))NA_real_ else mean(stability),B=B,seed=seed,failed=length(failures),failures=failures),
  n_total=nrow(data),n_used=n,n_excluded=sum(!keep),confidence_interval=NULL,interval_status='not_computed',requested_confidence_level=conf.level,
  limitations='Exploratory grouping only; k is not selected automatically. Silhouette and subsampling stability do not establish real populations or external validity. No p-values, uncertainty interval, or optimal k claim is supplied. Failed subsamples are reported and not selectively averaged. Ward.D2 squares Euclidean dissimilarities internally; no squared input distances are supplied.')
}
example_data<-data.frame(id=sprintf('case%02d',1:18),x=c(-5,-4.6,-5.2,-4.8,-5.1,-4.9,0,.4,-.2,.2,-.1,.1,5,5.4,4.8,5.2,4.9,5.1),y=c(-4,-4.2,-3.8,-4.1,-3.9,-4,1,.8,1.2,.9,1.1,1,5,4.8,5.2,4.9,5.1,5))
example_result<-analyse(example_data,c('x','y'),3,TRUE,TRUE)
print(example_result[c('status','method','k','mean_silhouette','stability')])
