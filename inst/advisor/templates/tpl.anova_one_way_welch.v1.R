# Synthetic independent-group omnibus example; diagnostics never choose a test.
check_input <- function(x,group,alternative,conf.level) {
 stopifnot(is.numeric(x),is.null(dim(x)),length(x)>=4L,all(is.finite(x)),
  is.factor(group)||is.character(group),length(group)==length(x),!anyNA(group),
  all(nzchar(as.character(group))),identical(alternative,"two.sided"),
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1)
 g<-droplevels(factor(group));stopifnot(nlevels(g)>=2L,all(table(g)>=2L));g
}
# Estimand: equality of population means, not a signed pairwise contrast.
# Independent normal errors are required for exact pooled F inference;
# Welch uses approximate degrees of freedom and permits unequal variances.
# Equal variances for the pooled test must be justified before inspecting data.
analyse <- function(x,group,alternative="two.sided",conf.level=.95) {
 g<-check_input(x,group,alternative,conf.level)
 v<-tapply(x,g,stats::var);means<-tapply(x,g,mean)
 stopifnot(all(is.finite(v)),all(is.finite(means)))
 if(any(v<=0)) return(list(status="not_applicable",p_value=NULL,
  reason="The required within-group variance is zero or unresolved.",confidence_interval=NULL,interval_status="not_computed"))
 fit<-stats::oneway.test(x~g,var.equal=FALSE)
 if(!is.finite(fit$statistic)||!is.finite(fit$p.value)||any(!is.finite(fit$parameter)))
  return(list(status="not_applicable",p_value=NULL,reason="The F calculation is not numerically resolved.",confidence_interval=NULL,interval_status="not_computed"))
 list(method=fit$method,status="ok",group_sizes=table(g),group_means=means,group_variances=v,
  statistic=unname(fit$statistic),df=unname(fit$parameter),p_value=fit$p.value,
  confidence_interval=NULL,interval_status="not_computed",requested_confidence_level=conf.level,
  alpha=1-conf.level,reject=fit$p.value<1-conf.level,
  posthoc="No pairwise comparisons computed. Specify contrasts and multiplicity control before analysis.",
  assumptions="Independent groups and units; normal within-group errors for the F model. Welch permits unequal population variances; reference degrees of freedom are approximate.")
}
example_x<-c(1,2,3,4,5, 2,3,5,8,9, 4,5,6,6,7)
example_group<-rep(c("A","B","C"),each=5)
example_result<-analyse(example_x,example_group,alternative,confidence_level)
print(example_result)
