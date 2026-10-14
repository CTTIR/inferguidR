# Synthetic teaching example for prespecified count-data inference.
# Counts must represent independent units, except for explicitly paired McNemar tables.
# No missing observations, survey weights, fitted null probabilities or automatic test switching.
check_settings <- function(alternative, conf.level) {
 stopifnot(identical(alternative,"two.sided"),is.numeric(conf.level),length(conf.level)==1L,
           is.finite(conf.level),conf.level>0,conf.level<1)
}
check_counts <- function(x) {
 stopifnot(is.numeric(x),length(x)>1L,all(is.finite(x)),all(x>=0),
           all(x==floor(x)),sum(x)>0,sum(x)<=2^31-1)
}
# Null: the prespecified category probabilities p; no parameters estimated from these counts.
# Pearson's omnibus statistic uses a right-tail chi-squared reference distribution.
analyse <- function(counts, p, alternative="two.sided", conf.level=.95) {
 check_settings(alternative,conf.level);check_counts(counts)
 stopifnot(is.null(dim(counts)),is.numeric(p),is.null(dim(p)),length(p)==length(counts),
           all(is.finite(p)),all(p>0),abs(sum(p)-1)<1e-12)
 messages<-character()
 fit<-withCallingHandlers(stats::chisq.test(counts,p=p,rescale.p=FALSE,simulate.p.value=FALSE),
   warning=function(w) {messages<<-c(messages,conditionMessage(w));invokeRestart("muffleWarning")})
 list(method=fit$method,status="ok",n=sum(counts),observed=counts,null_probabilities=p,
      observed_proportions=counts/sum(counts),expected=fit$expected,
      statistic=unname(fit$statistic),df=unname(fit$parameter),p_value=fit$p.value,
      cohen_w=sqrt(unname(fit$statistic)/sum(counts)),pearson_residuals=fit$residuals,
      confidence_interval=NULL,interval_status="not_computed",requested_confidence_level=conf.level,
      alpha=1-conf.level,reject=fit$p.value<1-conf.level,warnings=messages,
      assumptions="Independent multinomial units; prespecified probabilities; chi-squared approximation. Small expected cells can invalidate the reference approximation.")
}
example_result<-analyse(c(A=20,B=15,C=25),c(.3,.3,.4),alternative,confidence_level)
print(example_result)
