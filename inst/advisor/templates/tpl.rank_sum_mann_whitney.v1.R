# Synthetic teaching example: pre-specify the hypothesis and two independent groups.
# Rank inference tests identical distributions, not P(superiority)=.5 under
# unequal shapes/spreads. A median/mean interpretation needs stronger assumptions.
check_settings <- function(alternative,conf.level,exact,location_shift) {
 stopifnot(is.character(alternative),length(alternative)==1L,
  alternative %in% c("two.sided","greater","less"),
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1,
  is.logical(exact),length(exact)==1L,!is.na(exact),
  is.logical(location_shift),length(location_shift)==1L,!is.na(location_shift))
}
check_vector <- function(x) stopifnot(is.numeric(x),is.null(dim(x)),length(x)>0L,all(is.finite(x)))
# Ordered numeric codes may represent ordinal levels when location_shift=FALSE.
# location_shift=TRUE is an explicit interval-scale, common-shape shift assumption;
# never enable it for ordinal codes. It requests a separate shift interval, not
# an interval for the probability of superiority.
analyse <- function(x,y,alternative="two.sided",conf.level=.95,exact=TRUE,location_shift=FALSE) {
 check_settings(alternative,conf.level,exact,location_shift);check_vector(x);check_vector(y)
 ties <- anyDuplicated(c(x,y))>0L
 if(exact) stopifnot(!ties,length(x)+length(y)<=30L)
 ranks <- rank(c(x,y));u <- sum(ranks[seq_along(x)])-length(x)*(length(x)+1)/2
 superiority <- u/(length(x)*length(y))
 if(length(unique(c(x,y)))==1L) return(list(status="not_applicable",reason="All pooled values are equal; rank variance is zero.",
  p_value=NULL,probability_superiority=superiority,confidence_interval=NULL))
 # Preserve the same unrounded ranks used for descriptive superiority.
 messages <- character()
 fit <- withCallingHandlers(stats::wilcox.test(x,y,alternative=alternative,exact=exact,
   correct=TRUE,conf.int=location_shift,conf.level=conf.level,digits.rank=Inf),
  warning=function(w) {messages <<- c(messages,conditionMessage(w));invokeRestart("muffleWarning")})
 list(status=if(is.finite(fit$p.value)) "ok" else "undefined",method=fit$method,
  n_x=length(x),n_y=length(y),tied_ranks=ties,statistic=unname(fit$statistic),
  p_value=if(is.finite(fit$p.value)) fit$p.value else NULL,alternative=alternative,
  probability_superiority=superiority,rank_biserial=2*superiority-1,
  superiority_interval_status="not_computed",
  shift_estimate=if(location_shift) unname(fit$estimate) else NULL,
  confidence_interval=if(location_shift) as.numeric(fit$conf.int) else NULL,
  requested_confidence_level=conf.level,
  attained_confidence_level=if(location_shift) attr(fit$conf.int,"conf.level") else NULL,
  interval_status=if(location_shift) "shift_model_only" else "not_requested",
  interval_method="Hodges-Lehmann x-minus-y shift interval only under the explicitly assumed location-shift model",
  exact=exact,location_shift=location_shift,diagnostics=messages,
  limitations="Greater means x tends to exceed y. Probability of superiority is descriptive; the rank test assumes identical distributions under its null. No mean/median interpretation for general ordinal or unequal-shape data.")
}
example_x <- c(2,4,6,8,10)
example_y <- c(1,3,5,7,9)
example_result <- analyse(example_x,example_y,alternative,confidence_level,location_shift=FALSE)
print(example_result)
