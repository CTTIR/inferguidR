# Ordered-group rank trend under exchangeability of independent observations.
# Choose the group order and alternative before inspecting outcomes.
jt_twice_statistic <- function(y,g) {
 pairs<-which(outer(g,g,"<"),arr.ind=TRUE)
 sum(2*(y[pairs[,2]]>y[pairs[,1]])+(y[pairs[,2]]==y[pairs[,1]]))
}
analyse <- function(data,group_order,order_prespecified,independent_units,
                    seed,alternative="two.sided",B=1999L,conf.level=.95,
                    mc.level=.99,complete_cases=FALSE) {
 stopifnot(is.data.frame(data),all(c("id","group","y") %in% names(data)),
  nrow(data)>0L,nrow(data)<=200L,is.atomic(data$id),is.null(dim(data$id)),
  !anyNA(data$id),!anyDuplicated(data$id),all(nzchar(as.character(data$id))),
  is.character(group_order),length(group_order)>=2L,length(group_order)<=10L,
  !anyNA(group_order),all(nzchar(group_order)),!anyDuplicated(group_order),
  identical(order_prespecified,TRUE),identical(independent_units,TRUE),
  is.character(alternative),length(alternative)==1L,
  alternative %in% c("two.sided","greater","less"),
  is.numeric(seed),length(seed)==1L,is.finite(seed),seed==floor(seed),seed>=0,seed<=.Machine$integer.max,
  is.numeric(B),length(B)==1L,is.finite(B),B==floor(B),B>=999L,B<=1000000L,
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1,
  is.numeric(mc.level),length(mc.level)==1L,is.finite(mc.level),mc.level>0,mc.level<1,
  is.logical(complete_cases),length(complete_cases)==1L,!is.na(complete_cases),
  is.null(dim(data$y)),is.numeric(data$y)||is.ordered(data$y),
  is.null(dim(data$group)),is.character(data$group)||is.factor(data$group))
 y<-if(is.ordered(data$y)) as.integer(data$y) else data$y
 stopifnot(!any(is.infinite(y)),all(is.na(data$group)|as.character(data$group) %in% group_order))
 keep<-!is.na(y)&!is.na(data$group)
 if(!all(keep)&&!complete_cases) stop("Missing outcomes/groups require explicit complete_cases=TRUE.")
 n_total<-nrow(data);y<-y[keep];g<-match(as.character(data$group[keep]),group_order)
 counts<-tabulate(g,nbins=length(group_order));stopifnot(all(counts>0L))
 pairs<-which(outer(g,g,"<"),arr.ind=TRUE);pair_count<-nrow(pairs)
 if(B*pair_count>1e7) stop("Teaching limit: B times the number of cross-group pairs must not exceed ten million.")
 score<-function(v) sum(2*(v[pairs[,2]]>v[pairs[,1]])+(v[pairs[,2]]==v[pairs[,1]]))
 observed<-score(y);degenerate<-length(unique(y))==1L
 old_kind<-RNGkind();had_seed<-exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
 if(had_seed) old_seed<-get(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
 on.exit({do.call(RNGkind,as.list(old_kind));if(had_seed) assign(".Random.seed",old_seed,envir=.GlobalEnv) else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv)},add=TRUE)
 RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(as.integer(seed))
 lower<-upper<-0L
 if(!degenerate) for(i in seq_len(B)) {
  value<-score(y[sample.int(length(y))]);lower<-lower+(value<=observed);upper<-upper+(value>=observed)
 }
 tail_interval<-function(k,alpha) c(if(k==0L) 0 else stats::qbeta(alpha/2,k,B-k+1),
                                  if(k==B) 1 else stats::qbeta(alpha/2,k+1,B-k,lower.tail=FALSE))
 if(degenerate) {
  p<-1;interval<-c(1,1);tails<-c(lower=1,upper=1);used<-0L
  exceedances<-c(lower=NA_integer_,upper=NA_integer_)
 } else {
  tails<-c(lower=(lower+1)/(B+1),upper=(upper+1)/(B+1));used<-as.integer(B)
  exceedances<-c(lower=lower,upper=upper)
  p<-switch(alternative,less=tails[["lower"]],greater=tails[["upper"]],two.sided=min(1,2*min(tails)))
  if(alternative=="two.sided") {
   # Simultaneous Bonferroni tail intervals, then a monotone transformation.
   lo<-tail_interval(lower,(1-mc.level)/2);hi<-tail_interval(upper,(1-mc.level)/2)
   interval<-pmin(1,2*pmin(lo,hi))
  } else interval<-tail_interval(if(alternative=="less") lower else upper,1-mc.level)
 }
 list(status=if(degenerate) "no_rank_variation" else "ok",
  method="Jonckheere-Terpstra ordered-group Monte Carlo permutation test",
  group_order=group_order,group_sizes=stats::setNames(counts,group_order),order_prespecified=TRUE,
  n_total=n_total,n_used=length(y),n_excluded=sum(!keep),independent_units=TRUE,
  outcome_scale=if(is.ordered(data$y)) "ordered_categories" else "numeric_ranks",
  statistic=observed/2,cross_group_pairs=pair_count,descriptive_concordance=observed/(2*pair_count),
  alternative=alternative,p_value=p,seed=seed,requested_replicates=B,simulated_replicates=used,
  rng_kind=c("Mersenne-Twister","Inversion","Rejection"),simulated_exceedances=exceedances,
  plus_one_tails=tails,minimum_p_value=if(degenerate) 1 else if(alternative=="two.sided") 2/(B+1) else 1/(B+1),
  monte_carlo_interval=interval,monte_carlo_confidence_level=mc.level,
  precision_crosses_alpha=interval[1]<=1-conf.level&&interval[2]>=1-conf.level,
  alpha=1-conf.level,reject=p<1-conf.level,confidence_interval=NULL,interval_status="not_computed",
  inference_note="Greater means increasing outcomes along group_order; less means decreasing. Inclusive tails use (b+1)/(B+1). Two-sided p is min(1,2*min(lower,upper)); no normal or untied-exact fallback. All-equal outcomes have p=1 without simulation.",
  monte_carlo_note="The interval describes computational uncertainty in the conditional permutation p-value, not a population effect interval. Two-sided bounds transform simultaneous Bonferroni Clopper-Pearson tail intervals. Fix seed and B before analysis; no selective rerunning.",
  assumptions="Independent units and exchangeability under identical group distributions; group order and direction prespecified. Complete-case inference requires justified missingness. Concordance is descriptive and pair-weighted; this is not a mean, median, or causal-effect test.")
}
example_data<-data.frame(id=paste0("unit",1:9),group=rep(c("low","middle","high"),each=3),y=c(1,3,2,2,4,3,5,4,6))
example_result<-analyse(example_data,c("low","middle","high"),TRUE,TRUE,seed=104729L,
 alternative=alternative,B=1999L,conf.level=confidence_level)
print(example_result)
