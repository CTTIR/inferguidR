# Synthetic teaching example: pre-specify the method and event orientation.
# Independent binomial groups with fixed group sizes; the event margin is random.
# No paired observations, fractional counts, survey weights or automatic fallback.
# Greater means the event probability in row 1 exceeds that in row 2.
check_settings <- function(alternative, conf.level) {
 stopifnot(is.character(alternative),length(alternative)==1L,
  alternative %in% c("greater","less"),
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1)
}
analyse <- function(table,alternative="greater",conf.level=.95) {
 check_settings(alternative,conf.level)
 stopifnot(is.matrix(table),is.numeric(table),identical(dim(table),c(2L,2L)),
  all(is.finite(table)),all(table>=0),all(table==floor(table)),all(rowSums(table)>0))
 # Rows are group x, group y; columns are event, non-event.
 variant <- alternative
 fit <- inferguidR::boschloo_test(table[1,1],sum(table[1,]),table[2,1],sum(table[2,]),
  variant=variant,alpha=1-conf.level,groups=c("group_x","group_y"),event="event")
 list(status=fit$status$code,reason=fit$status$reason,method=fit$variant$label,
  table=table,n_total=sum(table),event_probability_difference=table[1,1]/sum(table[1,])-table[2,1]/sum(table[2,]),
  null_difference=0,p_value=fit$p_value,p_upper=fit$p_upper,alternative=alternative,
  confidence_interval=NULL,confidence_level=conf.level,interval_status="not_computed",
  interval_method="No sampling confidence interval is implemented by this example.",
  numerical_bound_note="p_upper bounds numerical optimisation error, not sampling uncertainty.",
  decision=fit$decision,observed_result=fit,
  limitations="Research preview: full numerical release certification remains pending. No substitution of a Fisher or asymptotic test.")
}
example_table <- matrix(c(4,1,1,4),nrow=2,byrow=TRUE,
 dimnames=list(c("group_x","group_y"),c("event","non_event")))
example_result <- analyse(example_table,alternative,confidence_level)
print(example_result)
