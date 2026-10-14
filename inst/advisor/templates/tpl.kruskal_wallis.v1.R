# Synthetic independent-group omnibus example; diagnostics never choose a test.
check_input <- function(x,group,alternative,conf.level) {
 stopifnot(is.numeric(x),is.null(dim(x)),length(x)>=4L,all(is.finite(x)),
  is.factor(group)||is.character(group),length(group)==length(x),!anyNA(group),
  all(nzchar(as.character(group))),identical(alternative,"two.sided"),
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1)
 g<-droplevels(factor(group));stopifnot(nlevels(g)>=2L,all(table(g)>=2L));g
}
# Rank omnibus test under identical group distributions, with average tied ranks.
# It is not a test of medians unless a common-shape location-shift model is justified.
# A rejection with unequal shapes/spreads need not establish stochastic ordering.
analyse <- function(x,group,alternative="two.sided",conf.level=.95) {
 g<-check_input(x,group,alternative,conf.level)
 if(length(unique(x))<2L) return(list(status="not_applicable",p_value=NULL,
  reason="All observations are tied; rank variation is zero.",confidence_interval=NULL,interval_status="not_computed"))
 fit<-stats::kruskal.test(x,g);stopifnot(is.finite(fit$statistic),is.finite(fit$p.value))
 list(method=fit$method,status="ok",group_sizes=table(g),group_mean_ranks=tapply(rank(x,ties.method="average"),g,mean),
  statistic=unname(fit$statistic),df=unname(fit$parameter),p_value=fit$p.value,
  confidence_interval=NULL,interval_status="not_computed",requested_confidence_level=conf.level,
  alpha=1-conf.level,reject=fit$p.value<1-conf.level,
  posthoc="No pairwise comparisons computed. Specify contrasts and multiplicity control before analysis.",
  assumptions="Independent groups and units, at least ordinal outcomes; identical distributions under the null, tie correction and approximate chi-squared reference.")
}
example_x<-c(1,2,3,4,5, 2,3,5,8,9, 4,5,6,6,7)
example_group<-rep(c("A","B","C"),each=5)
example_result<-analyse(example_x,example_group,alternative,confidence_level)
print(example_result)
