# Synthetic teaching example: independent binomial groups with fixed row sizes.
# Columns are event and non-event. Greater means row 1 has higher event probability.
# The null difference is zero. Two-sided extremeness uses absolute pooled Z.
# This is not Fisher ordering, Boschloo ordering, or a doubled one-sided p-value.
# Standalone teaching implementation; Exact::exact.test is a separate reference.
check_settings <- function(alternative, conf.level) {
 stopifnot(is.character(alternative),length(alternative)==1L,
  alternative %in% c("two.sided","less","greater"),is.numeric(conf.level),
  length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1)
}
analyse <- function(table, alternative = "two.sided", conf.level = 0.95, independent_binomial_groups = FALSE,
    tolerance = 1e-08, max_subdivisions = 4096L) {
    spec <- infercalcR::method_spec("unconditional_exact_z_pooled", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, table = table, independent_binomial_groups = independent_binomial_groups,
        tolerance = tolerance, max_subdivisions = max_subdivisions)$payload$native
}
example_table<-matrix(c(4,1,1,4),2,byrow=TRUE,
 dimnames=list(c("group_x","group_y"),c("event","non_event")))
example_result<-analyse(example_table,alternative,confidence_level,independent_binomial_groups=TRUE)
print(example_result)
