# Synthetic teaching example: fix the hypothesis and method before analysis.
# Independent evaluation units are required; no survey or fractional weights.
check_settings <- function(alternative, conf.level) {
 stopifnot(is.character(alternative), length(alternative)==1L,
  alternative %in% c("two.sided","greater","less"),
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1)
}
# Estimand: probability of a positive difference among nonzero differences.
# Interpreting this as a median test needs the corresponding continuity/design assumptions.
analyse <- function(x, mu0, alternative="two.sided", conf.level=.95) {
 check_settings(alternative,conf.level)
 stopifnot(is.numeric(x),is.null(dim(x)),length(x)>0L,all(is.finite(x)),
  is.numeric(mu0),length(mu0)==1L,is.finite(mu0))
 d <- x-mu0
 stopifnot(all(is.finite(d)))
 positive <- sum(d>0);negative <- sum(d<0);ties <- sum(d==0)
 n <- positive+negative
 if(n==0L) return(list(status="not_applicable",reason="All differences are zero; no informative signs.",
  n_total=length(d),n_nonzero=0L,n_ties=ties,p_value=NULL,confidence_interval=NULL))
 fit <- stats::binom.test(positive,n,p=.5,alternative=alternative,conf.level=conf.level)
 list(status="ok",method="Exact conditional sign test",n_total=length(d),n_nonzero=n,n_ties=ties,
  n_positive=positive,n_negative=negative,positive_probability=positive/n,
  null_probability=.5,p_value=fit$p.value,alternative=alternative,
  confidence_interval=as.numeric(fit$conf.int),confidence_level=conf.level,
  interval_method="Clopper-Pearson interval for P(positive | nonzero)",
  limitations="Zeros excluded explicitly; interval is for sign probability, not the median or mean difference.")
}
example_x <- c(5,6,7,4,3,8,7,2,5,9)
example_reference <- 4
example_result <- analyse(example_x,example_reference,alternative,confidence_level)
print(example_result)
