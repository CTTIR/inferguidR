# Synthetic teaching example for an independently prespecified spread question.
# Do not use a spread test as a preliminary rule for selecting a mean test.
check_settings <- function(alternative, conf.level, omnibus=FALSE) {
 stopifnot(is.character(alternative),length(alternative)==1L,!is.na(alternative),
  alternative%in%c("two.sided","greater","less"),is.numeric(conf.level),length(conf.level)==1L,
  is.finite(conf.level),conf.level>0,conf.level<1)
 if(omnibus) stopifnot(identical(alternative,"two.sided"))
}
check_sample <- function(x) stopifnot(is.numeric(x),is.null(dim(x)),length(x)>=2L,all(is.finite(x)))
# Estimand: population variance of x divided by population variance of y.
# Exact F inference requires independent normal samples; it is sensitive to non-normality.
analyse <- function(x,y,ratio0=1,alternative="two.sided",conf.level=.95) {
 check_settings(alternative,conf.level);check_sample(x);check_sample(y)
 stopifnot(is.numeric(ratio0),length(ratio0)==1L,is.finite(ratio0),ratio0>0,
  is.finite(stats::var(x)),is.finite(stats::var(y)),stats::var(x)>0,stats::var(y)>0)
 stopifnot(is.finite(stats::var(x)/stats::var(y)),
  is.finite((stats::var(x)/stats::var(y))/ratio0),stats::var(x)/stats::var(y)>0,
  (stats::var(x)/stats::var(y))/ratio0>0)
 fit<-stats::var.test(x,y,ratio=ratio0,alternative=alternative,conf.level=conf.level)
 list(method=fit$method,status="ok",n_x=length(x),n_y=length(y),
  variance_x=stats::var(x),variance_y=stats::var(y),variance_ratio=unname(fit$estimate),
  null_ratio=ratio0,alternative=alternative,statistic=unname(fit$statistic),
  df=unname(fit$parameter),p_value=fit$p.value,confidence_interval=as.numeric(fit$conf.int),
  confidence_level=conf.level,interval_estimand="population variance x / population variance y",
  alpha=1-conf.level,reject=fit$p.value<1-conf.level)
}
example_x<-c(1,2,4,6,7,9);example_y<-c(2,3,4,4,5,6)
example_result<-analyse(example_x,example_y,alternative=alternative,conf.level=confidence_level)
print(example_result)
