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
analyse <- function(table, alternative = "two.sided", conf.level = 0.95) {
    spec <- infercalcR::method_spec("fisher_exact_2x2", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, table = table)$payload$native
}
example_table <- matrix(c(8,2,1,9),nrow=2,byrow=TRUE,
 dimnames=list(c("group_x","group_y"),c("event","non_event")))
example_result <- analyse(example_table,alternative,confidence_level)
print(example_result)
