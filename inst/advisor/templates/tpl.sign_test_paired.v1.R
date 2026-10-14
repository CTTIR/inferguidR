# Synthetic teaching example: fix the hypothesis and method before analysis.
# Independent evaluation units are required; no survey or fractional weights.
check_settings <- function(alternative, conf.level) {
 stopifnot(is.character(alternative), length(alternative)==1L,
  alternative %in% c("two.sided","greater","less"),
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1)
}
# Estimand: probability of a positive difference among nonzero differences.
# Interpreting this as a median test needs the corresponding continuity/design assumptions.
analyse <- function(data, alternative="two.sided", conf.level=.95) {
 check_settings(alternative,conf.level)
 stopifnot(is.data.frame(data),all(c("id","x","y") %in% names(data)),nrow(data)>0L,
  is.atomic(data$id),is.null(dim(data$id)),!anyNA(data$id),!anyDuplicated(data$id),all(nzchar(as.character(data$id))),
  is.numeric(data$x),is.null(dim(data$x)),is.numeric(data$y),is.null(dim(data$y)),all(is.finite(data$x)),all(is.finite(data$y)))
 d <- data$x-data$y
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
example_data <- data.frame(id=paste0("unit",1:10),x=c(5,6,7,5,3,8,7,4,5,9),y=c(4,4,6,5,4,6,6,5,4,7))
example_result <- analyse(example_data,alternative,confidence_level)
print(example_result)
