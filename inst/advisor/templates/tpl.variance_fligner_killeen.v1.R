# Synthetic teaching example for an independently prespecified spread question.
# Do not use a spread test as a preliminary rule for selecting a mean test.
check_settings <- function(alternative, conf.level, omnibus=FALSE) {
 stopifnot(is.character(alternative),length(alternative)==1L,!is.na(alternative),
  alternative%in%c("two.sided","greater","less"),is.numeric(conf.level),length(conf.level)==1L,
  is.finite(conf.level),conf.level>0,conf.level<1)
 if(omnibus) stopifnot(identical(alternative,"two.sided"))
}
check_sample <- function(x) stopifnot(is.numeric(x),is.null(dim(x)),length(x)>=2L,all(is.finite(x)))
check_groups <- function(x,group) {
 check_sample(x)
 stopifnot(is.factor(group)||is.character(group),length(group)==length(x),!anyNA(group),all(nzchar(as.character(group))))
 group<-droplevels(factor(group));stopifnot(nlevels(group)>=2L,all(table(group)>=2L))
 group
}
# Fligner-Killeen: normal scores of ranked absolute deviations from group medians.
# The approximate chi-squared test is omnibus; no variance-ratio interval is implied.
analyse <- function(x,group,alternative="two.sided",conf.level=.95) {
 check_settings(alternative,conf.level,TRUE);g<-check_groups(x,group)
 z<-abs(x-stats::ave(x,g,FUN=stats::median))
 stopifnot(all(is.finite(z)))
 if(length(unique(z))<2L) return(list(status="not_applicable",p_value=NULL,
  reason="All absolute deviations are identical; the normal-score variance is zero.",confidence_interval=NULL,interval_status="not_computed"))
 fit<-stats::fligner.test(x,g)
 stopifnot(is.finite(fit$statistic),is.finite(fit$p.value))
 list(method=fit$method,status="ok",group_sizes=table(g),group_medians=tapply(x,g,stats::median),
  statistic=unname(fit$statistic),df=unname(fit$parameter),p_value=fit$p.value,
  confidence_interval=NULL,interval_status="not_computed",requested_confidence_level=conf.level,
  alpha=1-conf.level,reject=fit$p.value<1-conf.level,
  assumptions="Independent units and groups; median centering and average ranks for ties; approximate chi-squared reference. A pure spread interpretation needs comparable distribution shapes.")
}
example_x<-c(1,2,3,4,5, 1,2,4,7,9, 2,3,3,4,4)
example_group<-rep(c("A","B","C"),each=5)
example_result<-analyse(example_x,example_group,alternative,confidence_level)
print(example_result)
