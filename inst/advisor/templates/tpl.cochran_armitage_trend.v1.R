# Independent binomial groups with an order and scores fixed before observing outcomes.
# Score spacing is part of the hypothesis; arbitrary monotone recoding can change the test.
analyse <- function(events,totals,scores,scores_prespecified=FALSE,alternative="two.sided",conf.level=.95) {
 stopifnot(isTRUE(scores_prespecified),is.numeric(events),is.numeric(totals),is.numeric(scores),
   length(events)>=3L,length(events)==length(totals),length(events)==length(scores),
   all(is.finite(events)),all(is.finite(totals)),all(is.finite(scores)),
   all(totals>0),all(totals==floor(totals)),all(events>=0),all(events<=totals),all(events==floor(events)),
   all(diff(scores)>0),!is.null(names(events)),!anyNA(names(events)),all(nzchar(names(events))),!anyDuplicated(names(events)),
   identical(names(events),names(totals)),identical(names(events),names(scores)),
   is.character(alternative),length(alternative)==1L,alternative %in% c("two.sided","greater","less"),
   is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1)
 n<-sum(totals);stopifnot(n<=.Machine$integer.max)
 probability<-sum(events)/n
 if(probability==0||probability==1) return(list(method="Cochran-Armitage asymptotic trend test",status="unavailable",
   reason="No event variation: the binomial score variance is zero.",p_value=NULL,confidence_interval=NULL))
 span<-max(scores)-min(scores);stopifnot(is.finite(span),span>0)
 scaled<-(scores-min(scores))/span
 centred<-scaled-sum(totals*scaled)/n
 numerator<-sum(centred*(events-totals*probability))
 variance<-probability*(1-probability)*sum(totals*centred^2)
 stopifnot(is.finite(variance),variance>0)
 z<-numerator/sqrt(variance)
 fit<-stats::prop.trend.test(events,totals,score=scaled)
 p<-switch(alternative,two.sided=unname(fit$p.value),greater=stats::pnorm(z,lower.tail=FALSE),less=stats::pnorm(z))
 list(method="Cochran-Armitage asymptotic trend test",status="ok",events=events,totals=totals,scores=scores,
   scores_prespecified=TRUE,event_proportions=events/totals,pooled_event_probability=probability,
   signed_statistic=z,statistic=unname(fit$statistic),df=1L,p_value=p,alternative=alternative,
   direction="greater means increasing event probability across increasing prespecified scores",
   confidence_interval=NULL,interval_status="not_computed",requested_confidence_level=conf.level,
   expected_events=totals*probability,expected_non_events=totals*(1-probability),
   alpha=1-conf.level,reject=p<1-conf.level,
   assumptions="Independent binomial units with fixed group sizes and prespecified ordered scores. Asymptotic normal/chi-squared reference; sparse counts can make it inaccurate. No exact, continuity-corrected or automatic fallback. Scores are not fitted from outcomes. This is a trend score test, not an estimated common risk difference or logistic slope interval.")
}
example_result<-analyse(c(low=3,medium=8,high=13),c(low=20,medium=25,high=30),c(low=0,medium=1,high=2),scores_prespecified=TRUE,alternative=alternative,conf.level=confidence_level)
print(example_result)
