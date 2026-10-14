# Synthetic teaching example: fix the hypothesis and method before analysis.
# Independent evaluation units are required; no survey or fractional weights.
check_settings <- function(alternative, conf.level) {
 stopifnot(is.character(alternative), length(alternative)==1L,
  alternative %in% c("two.sided","greater","less"),
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1)
}
# Two independent groups with a prespecified binary event definition.
# Rows: group x, group y. Columns: event, non-event. The odds ratio is x/y.
# Exact conditional inference fixes both table margins, not a McNemar paired test.
analyse <- function(table,alternative="two.sided",conf.level=.95) {
 check_settings(alternative,conf.level)
 stopifnot(is.matrix(table),is.numeric(table),identical(dim(table),c(2L,2L)),
  all(is.finite(table)),all(table>=0),all(table==floor(table)),
  all(rowSums(table)>0),all(colSums(table)>0))
 fit <- stats::fisher.test(table,alternative=alternative,conf.level=conf.level,or=1)
 list(status="ok",method=fit$method,n_total=sum(table),table=table,
  conditional_odds_ratio=unname(fit$estimate),null_odds_ratio=1,
  p_value=fit$p.value,alternative=alternative,
  confidence_interval=as.numeric(fit$conf.int),confidence_level=conf.level,
  interval_method="Exact conditional odds-ratio interval",
  limitations="Conditional MLE differs from the sample cross-product odds ratio; two-sided p-value uses probability ordering.")
}
example_table <- matrix(c(8,2,1,9),nrow=2,byrow=TRUE,
 dimnames=list(c("group_x","group_y"),c("event","non_event")))
example_result <- analyse(example_table,alternative,confidence_level)
print(example_result)
