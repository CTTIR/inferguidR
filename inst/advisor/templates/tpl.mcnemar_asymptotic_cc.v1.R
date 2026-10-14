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
# Rows are time1 and columns time2, with event first at both times.
# Each participant contributes exactly one cell. The null is equal discordance probabilities.
# The descriptive marginal difference is time1 event probability minus time2 event probability.
analyse <- function(counts, alternative="two.sided", conf.level=.95) {
 check_settings(alternative,conf.level);check_counts(counts)
 stopifnot(is.matrix(counts),identical(dim(counts),c(2L,2L)))
 b<-counts[1,2];c<-counts[2,1];discordant<-b+c
 if(discordant==0) return(list(status="not_applicable",reason="No discordant pairs; the asymptotic statistic is undefined.",
   n=sum(counts),discordant_pairs=0,marginal_difference=0,p_value=NULL,confidence_interval=NULL,interval_status="not_computed"))
 fit<-stats::mcnemar.test(counts,correct=TRUE)
 list(method=fit$method,status="ok",n=sum(counts),observed=counts,discordant_pairs=discordant,
      marginal_difference=(b-c)/sum(counts),statistic=unname(fit$statistic),df=1,
      p_value=fit$p.value,continuity_correction=TRUE,
      confidence_interval=NULL,interval_status="not_computed",requested_confidence_level=conf.level,
      alpha=1-conf.level,reject=fit$p.value<1-conf.level,
      warnings=if(discordant<25) "Few discordant pairs; chi-squared approximation may be inaccurate." else character(),
      assumptions="Independent pairs; paired binary outcomes with fixed coding. No automatic conditional or mid-p substitution.")
}
example_counts<-rbind(Event=c(Event=45,Non_event=20),Non_event=c(Event=10,Non_event=25))
example_result<-analyse(example_counts,alternative,confidence_level)
print(example_result)
