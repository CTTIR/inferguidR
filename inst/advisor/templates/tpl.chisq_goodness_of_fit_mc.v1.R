# Independent multinomial counts against probabilities specified before seeing the data.
# Monte Carlo uncertainty is computational, not a sampling interval for population effects.
analyse <- function(observed,probabilities,seed,B=19999L,conf.level=.95,mc.level=.99) {
 stopifnot(is.numeric(observed),length(observed)>=2L,all(is.finite(observed)),all(observed>=0),all(observed==floor(observed)),
   is.numeric(probabilities),length(probabilities)==length(observed),all(is.finite(probabilities)),all(probabilities>0),
   !is.null(names(observed)),!anyNA(names(observed)),all(nzchar(names(observed))),!anyDuplicated(names(observed)),identical(names(observed),names(probabilities)),
   is.numeric(B),length(B)==1L,is.finite(B),B==floor(B),B>=999L,B<=1000000L,
   is.numeric(seed),length(seed)==1L,is.finite(seed),seed==floor(seed),seed>=0,seed<=.Machine$integer.max,
   is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1,
   is.numeric(mc.level),length(mc.level)==1L,is.finite(mc.level),mc.level>0,mc.level<1)
 n<-sum(observed);stopifnot(n>=1,n<=.Machine$integer.max,abs(sum(probabilities)-1)<=64*.Machine$double.eps)
 if(B*length(observed)>1e7) stop("Teaching limit: B times category count must not exceed ten million simulated cells.")
 supplied_probabilities<-probabilities;probabilities<-probabilities/sum(probabilities)
 # Explicit reproducible stream; restore both RNG kind and caller state even after failure.
 old_kind<-RNGkind();had_seed<-exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
 if(had_seed) old_seed<-get(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
 on.exit({do.call(RNGkind,as.list(old_kind));if(had_seed) assign(".Random.seed",old_seed,envir=.GlobalEnv) else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv)},add=TRUE)
 RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(as.integer(seed))
 fit<-stats::chisq.test(x=observed,p=probabilities,rescale.p=FALSE,simulate.p.value=TRUE,B=as.integer(B))
 exceed<-as.integer(round(fit$p.value*(B+1)-1));stopifnot(exceed>=0L,exceed<=B)
 tail<-(1-mc.level)/2
 interval<-c(if(exceed==0L) 0 else stats::qbeta(tail,exceed,B-exceed+1),
   if(exceed==B) 1 else stats::qbeta(tail,exceed+1,B-exceed,lower.tail=FALSE))
 list(method="Multinomial goodness of fit with Monte Carlo p-value",status="ok",observed=observed,
   supplied_probabilities=supplied_probabilities,probabilities=probabilities,expected=fit$expected,statistic=unname(fit$statistic),p_value=fit$p.value,
   n_total=n,B=B,seed=seed,rng_kind=c("Mersenne-Twister","Inversion","Rejection"),
   simulated_exceedances=exceed,minimum_p_value=1/(B+1),
   monte_carlo_standard_error=sqrt((exceed/B)*(1-exceed/B)/B),
   monte_carlo_interval=interval,monte_carlo_confidence_level=mc.level,
   monte_carlo_note="Plus-one p-value (b+1)/(B+1). The Clopper-Pearson interval describes simulation error in the tail probability; it is not an effect-size confidence interval. Plug-in Monte Carlo SE is zero at b=0 or B; the interval retains uncertainty.",
   confidence_interval=NULL,interval_status="not_computed",alpha=1-conf.level,reject=fit$p.value<1-conf.level,
   precision_crosses_alpha=interval[1]<=1-conf.level && interval[2]>=1-conf.level,
   assumptions="Independent multinomial observations; fully specified category probabilities, not estimated on these data. No automatic asymptotic replacement or selective rerunning until a desired decision. Category names and order must match.")
}
example_result<-analyse(c(A=1,B=3,C=6),c(A=.2,B=.3,C=.5),seed=104729L,B=19999L,conf.level=confidence_level)
print(example_result)
