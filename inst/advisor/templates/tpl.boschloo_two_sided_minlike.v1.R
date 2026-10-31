# Synthetic teaching example: pre-specify the method and event orientation.
# Independent binomial groups with fixed group sizes; the event margin is random.
# No paired observations, fractional counts, survey weights or automatic fallback.
# Greater means the event probability in row 1 exceeds that in row 2.
check_settings <- function(alternative, conf.level) {
 stopifnot(is.character(alternative),length(alternative)==1L,
  alternative %in% "two.sided",
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1)
}
analyse <- function(table, alternative = "two.sided", conf.level = 0.95) {
    spec <- infercalcR::method_spec("boschloo_two_sided_minlike", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, table = table)$payload$native
}
example_table <- matrix(c(4,1,1,4),nrow=2,byrow=TRUE,
 dimnames=list(c("group_x","group_y"),c("event","non_event")))
example_result <- analyse(example_table,alternative,confidence_level)
print(example_result)
