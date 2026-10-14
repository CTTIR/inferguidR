# Synthetic teaching example. The selected method and hypothesis are fixed before analysis.
# Reject missing and infinite inputs; document an explicit missingness strategy first.
check_settings <- function(alternative, conf.level) {
 stopifnot(length(alternative)==1L, alternative %in% c("two.sided","greater","less"),
           is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1)
}
# Independent binary trials with one fixed event definition and reference probability.
# The interval estimates the event probability, not a difference from the reference.
analyse <- function(events,n,p0,alternative="two.sided",conf.level=.95) {
 check_settings(alternative,conf.level)
 stopifnot(is.numeric(events),length(events)==1L,is.finite(events),events>=0,events==floor(events),
           is.numeric(n),length(n)==1L,is.finite(n),n>0,n==floor(n),events<=n,
           is.numeric(p0),length(p0)==1L,is.finite(p0),p0>0,p0<1)
 fit <- stats::binom.test(events,n,p=p0,alternative=alternative,conf.level=conf.level)
 list(method=fit$method,events=events,n=n,proportion=events/n,reference=p0,
      confidence_interval=fit$conf.int,p_value=fit$p.value,alternative=alternative)
}
example_result <- analyse(7,20,.2,alternative,confidence_level)
print(example_result)
graphics::barplot(c(event=7,non_event=13),ylab="Number of independent trials")
